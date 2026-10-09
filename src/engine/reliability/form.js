/**
 * Structural reliability: first-order (FORM) and second-order (SORM) reliability methods and importance sampling
 * around the design point.
 *
 * Limit-state convention: failure when g(x) < 0. Random variables are declared as `{ name: distSpec }`
 * (math/stats.js `dist`), correlations as described in reliability/sampling.js.
 */
import { makeRng } from '../math/rng.js';
import { normCdf, normInv, normPdf } from '../math/stats.js';
import { resolveVars, corrMatrix, makeTransform } from './sampling.js';

function setup(g, vars, opts, name) {
  if (typeof g !== 'function') throw new RangeError(`${name}: g must be a function`);
  const { names, dists } = resolveVars(vars);
  const R = corrMatrix(names, opts.corr);
  const tr = makeTransform(names, dists, R, { nataf: opts.nataf });
  const x = {};
  let evals = 0;
  const G = (u) => {
    evals++;
    const v = g(tr.toX(u, x));
    if (typeof v !== 'number' || Number.isNaN(v)) throw new Error(`${name}: limit-state function returned a non-numeric value`);
    return v;
  };
  return { names, dists, tr, G, count: () => evals };
}

/**
 * First-order reliability method.
 *
 * The variables are mapped to independent standard normal space u by the isoprobabilistic transformation
 * x_i = F_i⁻¹(Φ(z_i)), z = L u (Nataf model with Cholesky factor L of the adjusted correlation matrix when `corr`
 * is given, identity otherwise). The design point u* (point of the limit-state surface closest to the origin) is
 * found with the Hasofer–Lind–Rackwitz–Fiessler recursion
 *   u⁺ = [(∇G·u - G)/|∇G|²] ∇G
 * stabilised by the merit-function line search of the improved HL-RF algorithm; gradients are central finite
 * differences in u-space. β = α·u* with α = -∇G/|∇G| (so u* = β α; α_i > 0 for load-type variables whose
 * increase drives failure, α_i < 0 for resistance-type variables), pf = Φ(-β).
 *
 * Reference: Hasofer & Lind (1974), J. Eng. Mech. Div. ASCE 100, 111–121; Rackwitz & Fiessler (1978), Comput.
 * Struct. 9, 489–494; Zhang & Der Kiureghian (1995), "Two improved algorithms for reliability analysis";
 * Liu & Der Kiureghian (1986) for the Nataf model.
 * Validity: g continuous and differentiable near the design point; |β| <~ 8 (probabilities below 1e-15 are
 * outside double precision of the transformation); a single dominant design point.
 *
 * @param {(x: Object<string, number>) => number} g limit-state function, failure when g < 0 (argument object is
 *   reused between calls)
 * @param {Object<string, object>} vars { name: distSpec }
 * @param {{ corr?: object, tol?: number, maxIter?: number, fdStep?: number, nataf?: boolean, u0?: number[] }} [opts]
 *   corr: correlations; tol: convergence tolerance on |G/G0| and on the misalignment of u with ∇G (1e-6);
 *   maxIter (100); fdStep: finite-difference step in u-space (1e-4); u0: starting point in u-space (origin)
 * @returns {{ beta: number, pf: number, alpha: Object<string, number>, designPoint: Object<string, number>,
 *   iterations: number, converged: boolean, u: number[], importance: Object<string, number>, gradient: number[],
 *   evaluations: number }} reliability index [-], failure probability [-], unit direction cosines α [-],
 *   design point in physical units, iteration count, convergence flag, design point in u-space, α² importance
 *   factors, ∇G in u-space at the design point, number of g evaluations
 */
export function form(g, vars, opts = {}) {
  const s = setup(g, vars, opts, 'form');
  return formCore(s, opts);
}

function formCore(s, opts) {
  const { names, tr, G } = s;
  const n = names.length;
  const tol = opts.tol ?? 1e-6;
  const maxIter = opts.maxIter ?? 100;
  const h = opts.fdStep ?? 1e-4;
  if (!(tol > 0)) throw new RangeError('form: tol must be > 0');
  let u = opts.u0 ? Array.from(opts.u0) : new Array(n).fill(0);
  if (u.length !== n) throw new RangeError('form: u0 must have one entry per variable');
  const grad = new Array(n).fill(0);
  const up = new Array(n);
  const gradient = (uu) => {
    for (let i = 0; i < n; i++) up[i] = uu[i];
    for (let i = 0; i < n; i++) {
      up[i] = uu[i] + h;
      const gp = G(up);
      up[i] = uu[i] - h;
      const gm = G(up);
      up[i] = uu[i];
      grad[i] = (gp - gm) / (2 * h);
    }
  };
  let Gu = G(u);
  const G0 = Math.abs(Gu) > 0 ? Math.abs(Gu) : 1;
  let converged = false;
  let it = 0;
  let gnorm = 0;
  const trial = new Array(n);
  for (; it < maxIter; ) {
    gradient(u);
    gnorm = Math.hypot(...grad);
    if (!(gnorm > 0) || !Number.isFinite(gnorm)) {
      throw new Error('form: the limit-state function has a zero or non-finite gradient — it does not depend on the random variables at the current point');
    }
    let gu = 0;
    let unorm = 0;
    for (let i = 0; i < n; i++) {
      gu += grad[i] * u[i];
      unorm += u[i] * u[i];
    }
    unorm = Math.sqrt(unorm);
    // misalignment between u and the gradient direction
    let mis = 0;
    for (let i = 0; i < n; i++) {
      const r = u[i] - (gu / (gnorm * gnorm)) * grad[i];
      mis += r * r;
    }
    mis = Math.sqrt(mis);
    if (Math.abs(Gu) <= tol * G0 && mis <= tol * (1 + unorm)) {
      converged = true;
      break;
    }
    it++;
    const lam = (gu - Gu) / (gnorm * gnorm);
    let dn = 0;
    for (let i = 0; i < n; i++) {
      trial[i] = lam * grad[i] - u[i]; // search direction
      dn += trial[i] * trial[i];
    }
    dn = Math.sqrt(dn);
    // merit function m(u) = ½|u|² + c|G(u)|
    let c = 2 * (unorm / gnorm) + 1e-12;
    if (Math.abs(Gu) > tol * G0) {
      let un2 = 0;
      for (let i = 0; i < n; i++) un2 += (u[i] + trial[i]) * (u[i] + trial[i]);
      c = Math.max(c, (un2 / Math.abs(Gu)) * 1.0);
    }
    c = Math.max(c, 1e-8);
    const m0 = 0.5 * unorm * unorm + c * Math.abs(Gu);
    let step = 1;
    let accepted = false;
    const un = new Array(n);
    let Gn = Gu;
    for (let ls = 0; ls < 12; ls++) {
      let s2 = 0;
      for (let i = 0; i < n; i++) {
        un[i] = u[i] + step * trial[i];
        s2 += un[i] * un[i];
      }
      Gn = G(un);
      if (Number.isFinite(Gn) && 0.5 * s2 + c * Math.abs(Gn) < m0) {
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) {
      // no descent of the merit function: accept the smallest step if it is negligible, else give up
      if (step * dn <= tol * (1 + unorm) && Math.abs(Gu) <= Math.sqrt(tol) * G0) converged = true;
      break;
    }
    u = un.slice();
    Gu = Gn;
  }
  // direction cosines and reliability index at the final point
  let au = 0;
  let unorm = 0;
  const alphaArr = new Array(n);
  for (let i = 0; i < n; i++) {
    alphaArr[i] = -grad[i] / gnorm;
    au += alphaArr[i] * u[i];
    unorm += u[i] * u[i];
  }
  unorm = Math.sqrt(unorm);
  const beta = au >= 0 ? unorm : -unorm;
  const alpha = {};
  const importance = {};
  for (let i = 0; i < n; i++) {
    alpha[names[i]] = alphaArr[i];
    importance[names[i]] = alphaArr[i] * alphaArr[i];
  }
  const designPoint = { ...tr.toX(u, {}) };
  return {
    beta,
    pf: normCdf(-beta),
    alpha,
    designPoint,
    iterations: it,
    converged,
    u: u.slice(),
    importance,
    gradient: grad.slice(),
    evaluations: s.count(),
  };
}

/** Eigenvalues of a small symmetric matrix by cyclic Jacobi rotations. */
function symEigenvalues(Ain) {
  const n = Ain.length;
  const A = Ain.map((r) => r.slice());
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0;
    let diag = 0;
    for (let i = 0; i < n; i++) {
      diag += A[i][i] * A[i][i];
      for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j];
    }
    if (off <= 1e-28 * (diag + 1e-300)) break;
    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        if (A[p][q] === 0) continue;
        const theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1);
        const sn = t * c;
        for (let k = 0; k < n; k++) {
          const akp = A[k][p];
          const akq = A[k][q];
          A[k][p] = c * akp - sn * akq;
          A[k][q] = sn * akp + c * akq;
        }
        for (let k = 0; k < n; k++) {
          const apk = A[p][k];
          const aqk = A[q][k];
          A[p][k] = c * apk - sn * aqk;
          A[q][k] = sn * apk + c * aqk;
        }
      }
    }
  }
  return A.map((r, i) => r[i]).sort((a, b) => a - b);
}

/**
 * Second-order reliability method (curvature fitting, Breitung's asymptotic formula).
 *
 * After a FORM analysis the Hessian of G in u-space is evaluated by central finite differences at the design
 * point, rotated so that the last axis coincides with α, and the principal curvatures κ_i of the limit-state
 * surface are the eigenvalues of (R H Rᵀ)_{(n-1)×(n-1)} / |∇G| (κ_i > 0: surface curving away from the origin).
 *   pfBreitung = Φ(-β) Π (1 + β κ_i)^(-1/2)
 * `pfHohenbichler` applies the Hohenbichler–Rackwitz correction (β replaced by φ(β)/Φ(-β) inside the product),
 * which is more accurate for moderate β. betaGeneralised = -Φ⁻¹(pfBreitung).
 * Reference: Breitung (1984), J. Eng. Mech. 110(3), 357–366; Hohenbichler & Rackwitz (1988), J. Eng. Mech.
 * 114(12), 2195–2199; Der Kiureghian, Lin & Hwang (1987), J. Eng. Mech. 113(8), 1208–1225.
 * Validity: β > 0 and 1 + β κ_i > 0 for all i (otherwise `valid` is false and pfBreitung falls back to the FORM
 * value); asymptotically exact as β → ∞. Cost: FORM + (n² + n)·2 limit-state evaluations.
 *
 * @param {(x: Object<string, number>) => number} g limit-state function, failure when g < 0
 * @param {Object<string, object>} vars { name: distSpec }
 * @param {{ corr?: object, tol?: number, maxIter?: number, fdStep?: number, hessStep?: number }} [opts] as for
 *   {@link form}; hessStep: finite-difference step for the Hessian in u-space (5e-3)
 * @returns {{ beta: number, pfBreitung: number, betaGeneralised: number, curvatures: number[], pfForm: number,
 *   pfHohenbichler: number, valid: boolean, form: object }} FORM index, SORM probability, generalised index,
 *   principal curvatures [1/-] (ascending), FORM probability, Hohenbichler–Rackwitz probability, validity flag,
 *   the underlying FORM result
 */
export function sorm(g, vars, opts = {}) {
  const s = setup(g, vars, opts, 'sorm');
  const f = formCore(s, opts);
  const n = s.names.length;
  const u = f.u;
  const h = opts.hessStep ?? 5e-3;
  const G = s.G;
  const gnorm = Math.hypot(...f.gradient);
  const pfForm = f.pf;
  if (n === 1) {
    return { beta: f.beta, pfBreitung: pfForm, betaGeneralised: f.beta, curvatures: [], pfForm, pfHohenbichler: pfForm, valid: true, form: f };
  }
  const H = Array.from({ length: n }, () => new Array(n).fill(0));
  const w = u.slice();
  const g0 = G(w);
  for (let i = 0; i < n; i++) {
    w[i] = u[i] + h;
    const gp = G(w);
    w[i] = u[i] - h;
    const gm = G(w);
    w[i] = u[i];
    H[i][i] = (gp - 2 * g0 + gm) / (h * h);
    for (let j = 0; j < i; j++) {
      w[i] = u[i] + h;
      w[j] = u[j] + h;
      const pp = G(w);
      w[j] = u[j] - h;
      const pm = G(w);
      w[i] = u[i] - h;
      const mm = G(w);
      w[j] = u[j] + h;
      const mp = G(w);
      w[i] = u[i];
      w[j] = u[j];
      H[i][j] = (pp - pm - mp + mm) / (4 * h * h);
      H[j][i] = H[i][j];
    }
  }
  // orthonormal basis whose last vector is α (Gram–Schmidt against the coordinate axes)
  const alpha = s.names.map((k) => f.alpha[k]);
  const basis = [alpha.slice()];
  const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => Math.abs(alpha[a]) - Math.abs(alpha[b]));
  for (const k of order) {
    if (basis.length === n) break;
    const v = new Array(n).fill(0);
    v[k] = 1;
    for (let pass = 0; pass < 2; pass++) {
      for (const bvec of basis) {
        let d = 0;
        for (let i = 0; i < n; i++) d += v[i] * bvec[i];
        for (let i = 0; i < n; i++) v[i] -= d * bvec[i];
      }
    }
    const nv = Math.hypot(...v);
    if (nv < 1e-8) continue;
    basis.push(v.map((q) => q / nv));
  }
  const tang = basis.slice(1);
  const m = n - 1;
  const A = Array.from({ length: m }, () => new Array(m).fill(0));
  for (let a = 0; a < m; a++) {
    const Hv = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Hv[i] += H[i][j] * tang[a][j];
    for (let b = 0; b <= a; b++) {
      let sum = 0;
      for (let i = 0; i < n; i++) sum += tang[b][i] * Hv[i];
      A[a][b] = sum / gnorm;
      A[b][a] = A[a][b];
    }
  }
  const curvatures = symEigenvalues(A);
  const beta = f.beta;
  let valid = beta > 0;
  let prod = 1;
  let prodHr = 1;
  const psi = normPdf(beta) / normCdf(-beta);
  for (const k of curvatures) {
    const t = 1 + beta * k;
    const t2 = 1 + psi * k;
    if (!(t > 0) || !(t2 > 0)) {
      valid = false;
      break;
    }
    prod *= 1 / Math.sqrt(t);
    prodHr *= 1 / Math.sqrt(t2);
  }
  const pfBreitung = valid ? Math.min(pfForm * prod, 1) : pfForm;
  const pfHohenbichler = valid ? Math.min(pfForm * prodHr, 1) : pfForm;
  return {
    beta,
    pfBreitung,
    betaGeneralised: -normInv(pfBreitung),
    curvatures,
    pfForm,
    pfHohenbichler,
    valid,
    form: f,
  };
}

/**
 * Importance sampling with a standard normal sampling density centred on the design point in u-space:
 *   pf = (1/n) Σ I[g < 0]·φ_n(u)/φ_n(u - u*),  u ~ N(u*, I).
 * When `designPoint` is omitted a FORM analysis supplies it.
 * Reference: Melchers (1989), Struct. Saf. 6, 3–10; Melchers & Beck, "Structural Reliability Analysis and
 * Prediction", 3rd ed., §3.4. Validity: a single dominant design point; `cov` is the estimated coefficient of
 * variation of the pf estimator (aim for < 0.05–0.10).
 * @param {(x: Object<string, number>) => number} g limit-state function, failure when g < 0
 * @param {Object<string, object>} vars { name: distSpec }
 * @param {{ designPoint?: Object<string, number>, n?: number, seed?: number, corr?: object }} [opts]
 *   designPoint in physical units keyed by name (e.g. `form(...).designPoint`); n samples (2000); seed (1)
 * @returns {{ pf: number, cov: number, beta: number, nFail: number, n: number }} failure probability [-],
 *   coefficient of variation of the estimate [-], generalised reliability index [-], failures sampled
 */
export function importanceSampling(g, vars, opts = {}) {
  const s = setup(g, vars, opts, 'importanceSampling');
  const n = opts.n ?? 2000;
  if (!Number.isInteger(n) || n < 2) throw new RangeError('importanceSampling: n must be an integer >= 2');
  const d = s.names.length;
  let us;
  if (opts.designPoint) {
    for (const k of s.names) {
      if (typeof opts.designPoint[k] !== 'number') throw new RangeError(`importanceSampling: designPoint.${k} is missing`);
    }
    us = s.tr.toU(opts.designPoint);
    for (let i = 0; i < d; i++) if (!Number.isFinite(us[i])) throw new RangeError(`importanceSampling: designPoint.${s.names[i]} is outside the support of its distribution`);
  } else {
    us = formCore(s, opts).u;
  }
  const rng = makeRng(opts.seed ?? 1);
  let b2 = 0;
  for (let i = 0; i < d; i++) b2 += us[i] * us[i];
  const u = new Float64Array(d);
  let s1 = 0;
  let s2 = 0;
  let nFail = 0;
  for (let k = 0; k < n; k++) {
    let dot = 0;
    for (let i = 0; i < d; i++) {
      const e = rng.normal();
      u[i] = us[i] + e;
      dot += us[i] * e;
    }
    if (s.G(u) < 0) {
      // φ(u)/φ(u - u*) = exp(-u*·e - ½|u*|²)
      const wgt = Math.exp(-dot - 0.5 * b2);
      s1 += wgt;
      s2 += wgt * wgt;
      nFail++;
    }
  }
  const pf = s1 / n;
  const varEst = (s2 / n - pf * pf) / (n - 1);
  const cov = pf > 0 ? Math.sqrt(Math.max(varEst, 0)) / pf : Infinity;
  return { pf, cov, beta: pf > 0 && pf < 1 ? -normInv(pf) : Infinity, nFail, n };
}
