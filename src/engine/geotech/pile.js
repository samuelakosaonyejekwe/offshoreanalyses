// Single-pile analysis: beam on nonlinear Winkler foundation (lateral) and t–z / Q–z load transfer (axial).
//
// Conventions (lateral): z [m] downwards from the original mudline; the pile head is at z = −stickUp.
// y [m] lateral deflection, positive in the direction of a positive head load H.
// theta [rad] = −dy/dz, positive when the head tilts in the direction of H; a positive head moment M is the
// moment produced by H acting above the head (M = H·e), so H, M do positive work on y, theta.
// With this convention the head flexibility terms are all positive and the stiffness coupling KHM is negative.
// moment(z) and shear(z) are the internal actions in the same sense: moment(head) = M, shear(head) = H and
// d(moment)/dz = shear where no soil acts.

import { pyCurveAt, pyLinear, zeroCurve } from './py.js';
import { tzClay, tzSand, qz } from './tz.js';
import { axialCapacity, unitResistanceAt } from './axialcapacity.js';
import { propertiesAt } from './soilprofile.js';

function pos(name, v, allowZero = false) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < 0 || (!allowZero && v === 0)) {
    throw new RangeError(`${name}: must be a finite number ${allowZero ? '>= 0' : '> 0'} (got ${v})`);
  }
}

// Banded Gaussian elimination without pivoting (matrix rows stored as A[i*w + (j − i + hb)], w = 2 hb + 1).
// The pile matrices are symmetric with a dominant diagonal, so pivoting is not needed; a vanishing pivot means
// a mechanism or buckling and is reported.
function bandSolve(A, b, n, hb) {
  const w = 2 * hb + 1;
  for (let i = 0; i < n; i++) {
    const piv = A[i * w + hb];
    if (!(Math.abs(piv) > 1e-300)) throw new Error(`Singular pile stiffness matrix at degree of freedom ${i} (mechanism or buckling).`);
    const jm = Math.min(n - 1, i + hb);
    for (let r = i + 1; r <= jm; r++) {
      const f = A[r * w + (i - r + hb)] / piv;
      if (f === 0) continue;
      for (let c = i; c <= jm; c++) A[r * w + (c - r + hb)] -= f * A[i * w + (c - i + hb)];
      b[r] -= f * b[i];
    }
  }
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    const jm = Math.min(n - 1, i + hb);
    for (let c = i + 1; c <= jm; c++) s -= A[i * w + (c - i + hb)] * b[c];
    b[i] = s / A[i * w + hb];
  }
  return b;
}

function bandMul(A, x, out, n, hb) {
  const w = 2 * hb + 1;
  for (let i = 0; i < n; i++) {
    let s = 0;
    const j0 = Math.max(0, i - hb), j1 = Math.min(n - 1, i + hb);
    for (let j = j0; j <= j1; j++) s += A[i * w + (j - i + hb)] * x[j];
    out[i] = s;
  }
  return out;
}

function tubeSection(D, t) {
  const Di = D - 2 * t;
  return { A: (Math.PI * (D * D - Di * Di)) / 4, I: (Math.PI * (D ** 4 - Di ** 4)) / 64 };
}

function asCurve(c, name) {
  if (c === undefined || c === null) return null;
  if (typeof c === 'number') return pyLinear(c);
  if (typeof c.p === 'function' && typeof c.k === 'function') return c;
  throw new RangeError(`${name}: must be a stiffness number or a curve object with p(y) and k(y)`);
}

// Element stiffness in (y1, θ1, y2, θ2) with θ = −dy/dz. phi = 12 EI / (κ G A h²) is the Timoshenko shear
// parameter (0 for Euler–Bernoulli); P is the axial compression for the geometric (P–Δ) term.
function elementMatrix(EI, h, phi, P, ke) {
  const c = EI / ((1 + phi) * h * h * h);
  const a = 12 * c, b = 6 * h * c, d = (4 + phi) * h * h * c, e = (2 - phi) * h * h * c;
  ke[0] = a; ke[1] = -b; ke[2] = -a; ke[3] = -b;
  ke[4] = -b; ke[5] = d; ke[6] = b; ke[7] = e;
  ke[8] = -a; ke[9] = b; ke[10] = a; ke[11] = b;
  ke[12] = -b; ke[13] = e; ke[14] = b; ke[15] = d;
  if (P !== 0) {
    const g = P / (30 * h);
    const g36 = 36 * g, g3 = 3 * h * g, g4 = 4 * h * h * g, g1 = h * h * g;
    ke[0] -= g36; ke[1] += g3; ke[2] += g36; ke[3] += g3;
    ke[4] += g3; ke[5] -= g4; ke[6] -= g3; ke[7] += g1;
    ke[8] += g36; ke[9] -= g3; ke[10] -= g36; ke[11] -= g3;
    ke[12] += g3; ke[13] += g1; ke[14] -= g3; ke[15] -= g4;
  }
  return ke;
}

function buildLateralModel(o) {
  const { L, stickUp = 0, nElem = 60, E = 210e9, nu = 0.3, timoshenko = false, shearCoeff = 0.5, axialLoad = 0 } = o;
  pos('L', L); pos('stickUp', stickUp, true); pos('E', E);
  if (!(Number.isInteger(nElem) && nElem >= 2)) throw new RangeError('nElem: must be an integer >= 2');
  if (typeof axialLoad !== 'number' || !Number.isFinite(axialLoad)) throw new RangeError('axialLoad: must be a finite number');
  const hEmb = L / nElem;
  const nUp = stickUp > 0 ? Math.max(1, Math.ceil(stickUp / hEmb - 1e-9)) : 0;
  const ne = nUp + nElem;
  const n = ne + 1;
  const z = new Float64Array(n);
  for (let i = 0; i <= nUp; i++) z[i] = nUp ? -stickUp + (stickUp * i) / nUp : 0;
  for (let i = 1; i <= nElem; i++) z[nUp + i] = hEmb * i;

  const sectionAt = (zz) => {
    if (o.sections) {
      for (const s of o.sections) if (zz >= s.zTop - 1e-9 && zz <= s.zBot + 1e-9) return s;
      throw new RangeError(`sections: no section covers z = ${zz} m`);
    }
    return { D: o.D, t: o.t };
  };
  const elem = [];
  for (let e = 0; e < ne; e++) {
    const zm = 0.5 * (z[e] + z[e + 1]);
    const s = sectionAt(zm);
    let A, I, D = s.D;
    if (o.EI !== undefined && !o.sections) { pos('EI', o.EI); I = o.EI / E; A = o.A; } else { pos('D', s.D); pos('t', s.t); if (s.t >= s.D / 2) throw new RangeError('t: must be less than D/2'); ({ A, I } = tubeSection(s.D, s.t)); }
    let phi = 0;
    const h = z[e + 1] - z[e];
    if (timoshenko) {
      if (!(A > 0)) throw new RangeError('A: cross-section area is required for Timoshenko elements');
      phi = (12 * E * I) / (shearCoeff * (E / (2 * (1 + nu))) * A * h * h);
    }
    elem.push({ h, EI: E * I, I, A, D, t: s.t, phi });
  }

  // tributary lengths: soil acts on the embedded part only
  const trib = new Float64Array(n);
  for (let i = nUp; i < n; i++) trib[i] = (i > nUp ? 0.5 * hEmb : 0) + (i < n - 1 ? 0.5 * hEmb : 0);
  const pm = o.pMult === undefined ? () => 1 : typeof o.pMult === 'function' ? o.pMult : () => o.pMult;
  const curves = new Array(n);
  const mCurves = new Array(n).fill(null);
  for (let i = 0; i < n; i++) {
    if (i < nUp) { curves[i] = zeroCurve(); continue; }
    const D = o.sections ? sectionAt(z[i]).D : o.D;
    if (o.pyAt) curves[i] = asCurve(o.pyAt(z[i], D), 'pyAt(z)') ?? zeroCurve();
    else {
      if (!o.profile) throw new RangeError('profile: a soil profile or a pyAt(z) function is required');
      pos('D', D);
      curves[i] = pyCurveAt(o.profile, z[i], { D, loading: o.loading ?? 'static', ...(o.pyOptions ?? {}), pMult: pm(z[i]), yMult: o.yMult ?? 1 });
    }
    if (o.mThetaAt) mCurves[i] = asCurve(o.mThetaAt(z[i]), 'mThetaAt(z)');
  }
  const ndof = 2 * n, hb = 3, w = 7;
  const K = new Float64Array(ndof * w);
  const ke = new Float64Array(16);
  for (let e = 0; e < ne; e++) {
    elementMatrix(elem[e].EI, elem[e].h, elem[e].phi, axialLoad, ke);
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) { const i = 2 * e + a, j = 2 * e + b; K[i * w + (j - i + hb)] += ke[4 * a + b]; }
  }
  return {
    n, ndof, hb, w, z, elem, trib, curves, mCurves, K, nUp, axialLoad, E,
    baseShear: asCurve(o.base?.shear, 'base.shear'), baseMoment: asCurve(o.base?.moment, 'base.moment'),
  };
}

function safeTangent(curve, x, f) {
  let k = curve.k(x);
  if (!Number.isFinite(k)) k = curve.p(1e-9) / 1e-9;
  // softening branch: use the secant so the iteration matrix stays positive definite
  if (k < 0) k = x !== 0 ? f / x : 0;
  return k;
}

function secantOf(curve, x) {
  if (Math.abs(x) < 1e-12) { const k = curve.k(0); return Number.isFinite(k) ? k : curve.p(1e-9) / 1e-9; }
  return curve.p(x) / x;
}

// Internal force vector of all nonlinear springs; optionally fills the diagonal stiffness (tangent or secant).
function springForces(m, u, F, diag, mode) {
  const { n, trib, curves, mCurves } = m;
  for (let i = 0; i < n; i++) {
    let f = 0, k = 0, fm = 0, km = 0;
    if (trib[i] > 0) {
      const c = curves[i], y = u[2 * i];
      f = c.p(y) * trib[i];
      if (diag) k = (mode === 'secant' ? secantOf(c, y) : safeTangent(c, y, c.p(y))) * trib[i];
      const mc = mCurves[i];
      if (mc) { const th = u[2 * i + 1]; fm = mc.p(th) * trib[i]; if (diag) km = (mode === 'secant' ? secantOf(mc, th) : safeTangent(mc, th, mc.p(th))) * trib[i]; }
    }
    if (i === n - 1) {
      if (m.baseShear) { const y = u[2 * i]; const p = m.baseShear.p(y); f += p; if (diag) k += mode === 'secant' ? secantOf(m.baseShear, y) : safeTangent(m.baseShear, y, p); }
      if (m.baseMoment) { const th = u[2 * i + 1]; const p = m.baseMoment.p(th); fm += p; if (diag) km += mode === 'secant' ? secantOf(m.baseMoment, th) : safeTangent(m.baseMoment, th, p); }
    }
    F[2 * i] = f; F[2 * i + 1] = fm;
    if (diag) { diag[2 * i] = k; diag[2 * i + 1] = km; }
  }
}

function normaliseHead(head = {}) {
  const h = { H: head.H ?? 0, M: head.M ?? 0, y: head.y ?? null, theta: head.theta ?? null, kTrans: head.kTrans ?? 0, kRot: head.kRot ?? 0 };
  if (head.fixedRotation && h.theta === null) h.theta = 0;
  for (const key of ['H', 'M', 'kTrans', 'kRot']) if (typeof h[key] !== 'number' || !Number.isFinite(h[key])) throw new RangeError(`head.${key}: must be a finite number`);
  if (h.kTrans < 0 || h.kRot < 0) throw new RangeError('head.kTrans: spring stiffnesses must be >= 0');
  return h;
}

function solveModel(m, head, u0, tol, maxIter) {
  const { ndof, hb, w, K } = m;
  const u = new Float64Array(ndof);
  if (u0) u.set(u0);
  const fixY = head.y !== null, fixT = head.theta !== null;
  if (fixY) u[0] = head.y;
  if (fixT) u[1] = head.theta;
  const Fs = new Float64Array(ndof), diag = new Float64Array(ndof), R = new Float64Array(ndof), Ku = new Float64Array(ndof);
  const A = new Float64Array(ndof * w), trial = new Float64Array(ndof);
  let floor = 0;
  const residual = (uu, withDiag) => {
    springForces(m, uu, Fs, withDiag ? diag : null, 'tangent');
    bandMul(K, uu, Ku, ndof, hb);
    if (withDiag) {
      // round-off floor of K·u: with stiff piles the products cancel to about eps·|K||u|
      let a2 = 0;
      for (let i = 0; i < ndof; i++) {
        let a = 0;
        const j0 = Math.max(0, i - hb), j1 = Math.min(ndof - 1, i + hb);
        for (let j = j0; j <= j1; j++) a += Math.abs(K[i * w + (j - i + hb)] * uu[j]);
        a2 += a * a;
      }
      floor = 256 * Number.EPSILON * Math.sqrt(a2);
    }
    let r2 = 0, s2 = 0;
    for (let i = 0; i < ndof; i++) {
      let ext = 0, int = Ku[i] + Fs[i];
      if (i === 0) { ext = head.H; int += head.kTrans * uu[0]; }
      if (i === 1) { ext = head.M; int += head.kRot * uu[1]; }
      R[i] = int - ext;
      if ((i === 0 && fixY) || (i === 1 && fixT)) R[i] = 0;
      r2 += R[i] * R[i]; s2 += Ku[i] * Ku[i] + ext * ext;
    }
    return { r: Math.sqrt(r2), s: Math.sqrt(s2) };
  };
  let res = residual(u, true), it = 0;
  const ok = (q) => q.r <= Math.max(tol * q.s, floor, 1e-12);
  for (; it < maxIter; it++) {
    if (ok(res)) break;
    A.set(K);
    for (let i = 0; i < ndof; i++) A[i * w + hb] += diag[i];
    A[hb] += head.kTrans; A[w + hb] += head.kRot;
    const du = new Float64Array(ndof);
    for (let i = 0; i < ndof; i++) du[i] = -R[i];
    for (const d of [fixY ? 0 : -1, fixT ? 1 : -1]) {
      if (d < 0) continue;
      for (let j = Math.max(0, d - hb); j <= Math.min(ndof - 1, d + hb); j++) { A[d * w + (j - d + hb)] = 0; A[j * w + (d - j + hb)] = 0; }
      A[d * w + hb] = 1; du[d] = 0;
    }
    // a vanishing pivot means every spring along the pile is at its limit: collapse, reported as non-convergence
    try { bandSolve(A, du, ndof, hb); } catch (err) { return { u, converged: false, iterations: it, residual: res.s > 0 ? res.r / res.s : res.r }; }
    let step = 1, next = null;
    const r0 = res.r;
    for (let ls = 0; ls < 12; ls++) {
      for (let i = 0; i < ndof; i++) trial[i] = u[i] + step * du[i];
      next = residual(trial, false);
      if (next.r < r0 || ok(next)) break;
      step *= 0.5;
    }
    u.set(trial);
    res = residual(u, true);
  }
  const converged = ok(res);
  return { u, converged, iterations: it, residual: res.s > 0 ? res.r / res.s : res.r };
}

function postLateral(m, head, sol) {
  const { n, z, elem, trib, curves } = m;
  const u = sol.u;
  const y = new Float64Array(n), theta = new Float64Array(n), moment = new Float64Array(n), shear = new Float64Array(n);
  const p = new Float64Array(n), pu = new Float64Array(n), mobilisation = new Float64Array(n), sigmaBending = new Float64Array(n);
  const ke = new Float64Array(16);
  const Fs = new Float64Array(m.ndof);
  springForces(m, u, Fs, null, 'tangent');
  for (let i = 0; i < n; i++) {
    y[i] = u[2 * i]; theta[i] = u[2 * i + 1];
    if (trib[i] > 0) {
      p[i] = curves[i].p(y[i]); pu[i] = curves[i].pu;
      mobilisation[i] = Number.isFinite(pu[i]) && pu[i] > 0 ? Math.abs(p[i]) / pu[i] : 0;
    }
  }
  for (let e = 0; e < n - 1; e++) {
    // elastic end forces only: moment = EI·curvature at the nodes, shear = d(moment)/dz
    elementMatrix(elem[e].EI, elem[e].h, elem[e].phi, 0, ke);
    const f = [0, 0, 0, 0];
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) f[a] += ke[4 * a + b] * u[2 * e + b];
    shear[e] = f[0]; moment[e] = f[1];
    if (e === n - 2) { shear[n - 1] = -f[2]; moment[n - 1] = -f[3]; }
  }
  // actions just above the head-node springs equal the external head load (net of head springs)
  shear[0] += Fs[0]; moment[0] += Fs[1];
  let iM = 0, zero = null;
  for (let i = 0; i < n; i++) {
    if (Math.abs(moment[i]) > Math.abs(moment[iM])) iM = i;
    const el = elem[Math.min(i, n - 2)];
    sigmaBending[i] = el.I > 0 && el.D ? (Math.abs(moment[i]) * el.D) / (2 * el.I) : NaN;
    if (zero === null && i > 0 && z[i] > 0 && y[i - 1] * y[i] < 0) zero = z[i - 1] + ((z[i] - z[i - 1]) * y[i - 1]) / (y[i - 1] - y[i]);
  }
  return {
    z, y, theta, moment, shear, p, pu, mobilisation, sigmaBending,
    head: { y: y[0], theta: theta[0], H: shear[0], M: moment[0] },
    toe: { y: y[n - 1], theta: theta[n - 1] },
    maxMoment: Math.abs(moment[iM]), depthMaxMoment: z[iM], depthFirstZeroDeflection: zero,
    maxMobilisation: Math.max(...mobilisation),
    convergence: { converged: sol.converged, iterations: sol.iterations, residual: sol.residual, method: 'Newton–Raphson with secant fallback on softening springs and step halving' },
  };
}

function runLateral(m, head, o, u0) {
  const tol = o.tol ?? 1e-9, maxIter = o.maxIter ?? 80;
  let sol = solveModel(m, head, u0, tol, maxIter);
  // if a direct solve stalls, ramp the load so each increment starts close to its solution
  for (const steps of [4, 16, 64]) {
    if (sol.converged) break;
    let u = u0 ?? null, ok = true, iters = 0;
    for (let s = 1; s <= steps && ok; s++) {
      const f = s / steps;
      const hs = { ...head, H: head.H * f, M: head.M * f, y: head.y === null ? null : head.y * f, theta: head.theta === null ? null : head.theta * f };
      sol = solveModel(m, hs, u, tol, maxIter);
      ok = sol.converged; u = sol.u; iters += sol.iterations;
    }
    sol.iterations = iters;
  }
  return sol;
}

/**
 * Laterally loaded pile as a beam on a nonlinear Winkler foundation (finite elements).
 * Cubic Hermite Euler–Bernoulli elements (optionally Timoshenko with the shear parameter Φ = 12EI/(κGAh²)),
 * geometric stiffness for a constant axial compression (P–Δ), lumped p–y springs at the nodes with tributary
 * lengths, optional base shear / base moment springs and distributed moment (m–θ) springs supplied as curves.
 * Head conditions: free (H, M); fixed rotation (fixedRotation: true); springs (kTrans, kRot); prescribed
 * displacement and/or rotation (y, theta) — the corresponding reaction is reported in result.head.
 * Reference: Reese & Van Impe (2011) Single Piles and Pile Groups under Lateral Loading, (governing
 * equation EI y'''' + P y'' + p(y) = 0); API RP 2GEO (2011) §8.5 for the soil curves.
 * Validity: small rotations; linear-elastic pile; soil response uncoupled between depths (Winkler).
 * @param {object} o
 *   L [m] embedded length below the original mudline; stickUp [m] free length above the mudline (default 0);
 *   nElem [-] elements over the embedded length (default 60; the stick-up is meshed at the same size);
 *   D [m], t [m] outer diameter and wall thickness, or sections: [{ zTop, zBot, D, t }] (z from mudline, m), or
 *   EI [N·m²] directly (with A [m²] if Timoshenko); E [Pa] (default 210e9); nu [-] (0.3);
 *   timoshenko [bool] (false); shearCoeff κ [-] (default 0.5, thin-walled tube);
 *   axialLoad [N] compression positive, constant along the pile (default 0);
 *   profile soil profile and loading 'static'|'cyclic', pyOptions (passed to pyCurveAt), pMult [-] number or
 *   function of z, yMult [-]; or pyAt(z, D) → curve object | stiffness [N/m²] to override the soil model;
 *   base: { shear, moment } curve objects or stiffnesses [N/m], [N·m/rad]; mThetaAt(z) → curve [N·m/m vs rad];
 *   head: { H [N], M [N·m], y [m], theta [rad], kTrans [N/m], kRot [N·m/rad], fixedRotation [bool] };
 *   tol [-] relative residual tolerance (1e-9); maxIter [-] (80); throwOnFail [bool] (true)
 * @returns {object} { z [m], y [m], theta [rad], moment [N·m], shear [N], p [N/m] soil reaction (resisting,
 *   same sign as y), pu [N/m], mobilisation |p|/pu [-], sigmaBending [Pa] (NaN if no section given) — all
 *   Float64Array per node; head { y, theta, H, M } (loads net of head springs, or reactions when prescribed);
 *   toe { y, theta }; maxMoment [N·m]; depthMaxMoment [m]; depthFirstZeroDeflection [m]|null;
 *   maxMobilisation [-]; convergence { converged, iterations, residual [-], method } }
 */
export function solveLateralPile(o) {
  const m = buildLateralModel(o);
  const head = normaliseHead(o.head);
  const sol = runLateral(m, head, o, null);
  if (!sol.converged && (o.throwOnFail ?? true)) throw new Error(`Lateral pile analysis did not converge (relative residual ${sol.residual.toExponential(2)}); the load may exceed the lateral capacity.`);
  return postLateral(m, head, sol);
}

/**
 * Pile-head load–displacement and moment–rotation curves: the head loads of o.head are scaled by each factor.
 * Reference: see solveLateralPile.
 * @param {object} o solveLateralPile options (o.head holds the full load)
 * @param {{ factors?: number[] }} [c] factors [-] ascending load multipliers (default 0.1, 0.2 … 1.0)
 * @returns {{ factor: number[], H: number[], M: number[], y: number[], theta: number[], maxMoment: number[] }}
 *   H [N], M [N·m], y [m], theta [rad], maxMoment [N·m] for each factor
 */
export function pileHeadCurves(o, c = {}) {
  const factors = c.factors ?? Array.from({ length: 10 }, (_, i) => (i + 1) / 10);
  const m = buildLateralModel(o);
  const head = normaliseHead(o.head);
  const out = { factor: [], H: [], M: [], y: [], theta: [], maxMoment: [] };
  let u = null;
  for (const f of factors) {
    const hs = { ...head, H: head.H * f, M: head.M * f, y: head.y === null ? null : head.y * f, theta: head.theta === null ? null : head.theta * f };
    const sol = runLateral(m, hs, o, u);
    if (!sol.converged) throw new Error(`Lateral pile analysis did not converge at load factor ${f}.`);
    u = sol.u;
    const r = postLateral(m, hs, sol);
    out.factor.push(f); out.H.push(r.head.H); out.M.push(r.head.M); out.y.push(r.head.y); out.theta.push(r.head.theta); out.maxMoment.push(r.maxMoment);
  }
  return out;
}

function condense(m, u, mode) {
  const { ndof, hb, w } = m;
  const Fs = new Float64Array(ndof), diag = new Float64Array(ndof);
  springForces(m, u, Fs, diag, mode);
  const cols = [];
  for (const d of [0, 1]) {
    const A = new Float64Array(m.K);
    for (let i = 0; i < ndof; i++) A[i * w + hb] += diag[i];
    const b = new Float64Array(ndof); b[d] = 1;
    bandSolve(A, b, ndof, hb);
    cols.push([b[0], b[1]]);
  }
  const fHH = cols[0][0], fHM = 0.5 * (cols[0][1] + cols[1][0]), fMM = cols[1][1];
  const det = fHH * fMM - fHM * fHM;
  return { KHH: fMM / det, KHM: -fHM / det, KMM: fHH / det, fHH, fHM, fMM, asymmetry: Math.abs(cols[0][1] - cols[1][0]) / Math.abs(fHM) };
}

/**
 * Pile-head stiffness matrix at the working load of o.head, for coupling to a structural model:
 *   [H; M] = [[KHH, KHM], [KHM, KMM]] · [y; theta].
 * tangent: stiffness of the linearised system at the converged state (incremental loads).
 * secant: every soil spring replaced by its secant p/y at the converged state; this matrix is symmetric positive
 * definite and reproduces the working-point response exactly (K_secant · [y; theta] = [H; M]).
 * Head springs (kTrans, kRot) and head restraints are used to find the working state but are not part of the
 * matrices: they describe the pile and soil below the head only.
 * Reference: static condensation of the Winkler model; see solveLateralPile. Validity: as solveLateralPile.
 * @param {object} o solveLateralPile options
 * @returns {{ tangent: object, secant: object, working: object }} tangent and secant: { KHH [N/m], KHM [N]
 *   (negative), KMM [N·m/rad], fHH [m/N], fHM [1/N], fMM [rad/(N·m)], asymmetry [-] }; working = head state
 *   { y, theta, H, M }
 */
export function pileHeadStiffness(o) {
  const m = buildLateralModel(o);
  const head = normaliseHead(o.head);
  const sol = runLateral(m, head, o, null);
  if (!sol.converged) throw new Error('Lateral pile analysis did not converge at the working load.');
  return { tangent: condense(m, sol.u, 'tangent'), secant: condense(m, sol.u, 'secant'), working: postLateral(m, head, sol).head };
}

/**
 * Equivalent cantilever (apparent fixity) lengths for a pile-head stiffness matrix, from the cantilever
 * flexibilities fHH = Le³/(3EI), fHM = Le²/(2EI), fMM = Le/EI and the fixed-head stiffness KHH = 12EI/Le³.
 * A single length cannot match all terms; each value below matches one of them.
 * Reference: elementary beam theory (first principles). Validity: prismatic equivalent member of stiffness EI.
 * @param {{ KHH: number, KHM: number, KMM: number, EI: number }} o KHH [N/m]; KHM [N]; KMM [N·m/rad]; EI [N·m²]
 * @returns {{ lateral: number, rotation: number, coupling: number, fixedHead: number }} lengths [m] measured from
 *   the point where the stiffness matrix is defined: lateral (free-head deflection under H), rotation (rotation
 *   under M), coupling (cross term), fixedHead (translation without rotation)
 */
export function equivalentCantilever({ KHH, KHM, KMM, EI }) {
  pos('KHH', KHH); pos('KMM', KMM); pos('EI', EI);
  const det = KHH * KMM - KHM * KHM;
  if (!(det > 0)) throw new RangeError('KHM: stiffness matrix must be positive definite');
  const fHH = KMM / det, fHM = -KHM / det, fMM = KHH / det;
  return { lateral: Math.cbrt(3 * EI * fHH), rotation: EI * fMM, coupling: Math.sqrt(2 * EI * Math.abs(fHM)), fixedHead: Math.cbrt((12 * EI) / KHH) };
}

/**
 * Ultimate lateral capacity by load stepping: H is increased with M = H·eccentricity (or zero head rotation
 * when o.head.fixedRotation) until the analysis no longer converges (soil failure mechanism) or the head
 * deflection exceeds yLimit; the last step is refined by bisection.
 * The pile is elastic — compare maxMoment with the section capacity separately (a plastic hinge is not modelled).
 * Reference: load-controlled limit analysis of the Winkler model; see solveLateralPile. Validity: curves with a
 * finite ultimate resistance (not linear springs).
 * @param {object} o solveLateralPile options (o.head.H and o.head.M are ignored)
 * @param {{ eccentricity?: number, yLimit?: number, nSteps?: number, tol?: number }} [c] eccentricity [m] of H
 *   above the head (default 0); yLimit [m] head deflection limit (default 0.1 D, a displacement criterion chosen
 *   by the user); nSteps [-] load steps up to the summed soil resistance (default 40); tol [-] relative
 *   bisection tolerance on H (default 0.005)
 * @returns {{ Hult: number, criterion: string, yAtHult: number, maxMoment: number, curve: { H: number[],
 *   y: number[], theta: number[] } }} Hult [N]; criterion 'non-convergence' | 'displacement limit'; yAtHult [m];
 *   maxMoment [N·m] at Hult; curve of converged steps
 */
export function ultimateLateralCapacity(o, c = {}) {
  const m = buildLateralModel(o);
  const e = c.eccentricity ?? 0, nSteps = c.nSteps ?? 40, tol = c.tol ?? 0.005;
  const Dref = o.D ?? o.sections?.[0]?.D;
  const yLimit = c.yLimit ?? (Dref ? 0.1 * Dref : NaN);
  pos('yLimit', yLimit);
  let Hsum = 0;
  for (let i = 0; i < m.n; i++) Hsum += m.curves[i].pu * m.trib[i];
  if (!Number.isFinite(Hsum) || Hsum <= 0) throw new RangeError('pyAt: ultimate capacity needs soil curves with a finite ultimate resistance');
  const base = normaliseHead({ ...o.head, H: 0, M: 0 });
  const opts = { ...o, maxIter: o.maxIter ?? 60 };
  const curve = { H: [], y: [], theta: [] };
  let u = null, Hok = 0, last = null, crit = 'non-convergence';
  const attempt = (H) => {
    const hs = { ...base, H, M: H * e };
    const sol = solveModel(m, hs, u, opts.tol ?? 1e-9, opts.maxIter);
    if (!sol.converged) return null;
    return { sol, r: postLateral(m, hs, sol) };
  };
  const dH = Hsum / nSteps;
  let Hfail = null;
  for (let s = 1; s <= nSteps; s++) {
    const a = attempt(s * dH);
    if (!a) { Hfail = s * dH; break; }
    if (Math.abs(a.r.head.y) > yLimit) { Hfail = s * dH; crit = 'displacement limit'; break; }
    u = a.sol.u; Hok = s * dH; last = a.r;
    curve.H.push(Hok); curve.y.push(a.r.head.y); curve.theta.push(a.r.head.theta);
  }
  if (Hfail !== null) {
    while (Hfail - Hok > tol * Hfail) {
      const Hm = 0.5 * (Hok + Hfail);
      const a = attempt(Hm);
      if (a && Math.abs(a.r.head.y) <= yLimit) { u = a.sol.u; Hok = Hm; last = a.r; curve.H.push(Hm); curve.y.push(a.r.head.y); curve.theta.push(a.r.head.theta); } else Hfail = Hm;
    }
  }
  if (!last) throw new Error('Ultimate lateral capacity: the first load step did not converge.');
  return { Hult: Hok, criterion: crit, yAtHult: last.head.y, maxMoment: last.maxMoment, curve };
}

/**
 * Critical pile length: the shortest embedded length for which the head deflection and rotation under the
 * working load are within tol of those of the longest pile considered (further length does not stiffen the
 * response).
 * Reference: definition used for laterally loaded piles, e.g. Reese & Van Impe (2011); see solveLateralPile.
 * Validity: Lmax must be long enough to be "flexible"; check that the response has levelled off in the table.
 * @param {object} o solveLateralPile options without L
 * @param {{ Lmin: number, Lmax: number, nLengths?: number, tol?: number }} c Lmin, Lmax [m]; nLengths [-]
 *   (default 16); tol [-] relative tolerance (default 0.02)
 * @returns {{ Lcrit: number|null, L: number[], y: number[], theta: number[], toeY: number[] }} Lcrit [m] (null if
 *   only Lmax qualifies trivially); head y [m], theta [rad] and toe deflection [m] for each trial length (NaN where the pile is too
 *   short to reach equilibrium)
 */
export function criticalPileLength(o, c) {
  const { Lmin, Lmax, nLengths = 16, tol = 0.02 } = c;
  pos('Lmin', Lmin); pos('Lmax', Lmax); pos('tol', tol);
  if (!(Lmax > Lmin)) throw new RangeError('Lmax: must exceed Lmin');
  const hRef = Lmax / (o.nElem ?? 60);
  const out = { Lcrit: null, L: [], y: [], theta: [], toeY: [] };
  for (let i = 0; i < nLengths; i++) {
    const L = Lmin + ((Lmax - Lmin) * i) / (nLengths - 1);
    const r = solveLateralPile({ ...o, L, nElem: Math.max(8, Math.round(L / hRef)), throwOnFail: false });
    const ok = r.convergence.converged; // a pile too short to carry the load has no equilibrium state
    out.L.push(L); out.y.push(ok ? r.head.y : NaN); out.theta.push(ok ? r.head.theta : NaN); out.toeY.push(ok ? r.toe.y : NaN);
  }
  const yRef = out.y[nLengths - 1], tRef = out.theta[nLengths - 1];
  if (!Number.isFinite(yRef)) throw new Error('Critical pile length: the analysis did not converge for Lmax.');
  for (let i = nLengths - 2; i >= 0; i--) {
    const ok = Math.abs(out.y[i] / yRef - 1) <= tol && Math.abs(out.theta[i] / tRef - 1) <= tol;
    if (!ok) break;
    out.Lcrit = out.L[i];
  }
  return out;
}

/**
 * Elastic steel stress utilisation of a tubular pile section under axial force and bending:
 *   σ = |N|/A + |M| D / (2 I), utilisation = σ / (usageFactor · fy).
 * This is a first-yield stress check, not a member code check (no buckling or hoop-stress interaction).
 * Reference: elementary beam theory. Validity: elastic section, D/t low enough that local buckling does not govern.
 * @param {{ N: number, M: number, D: number, t: number, fy: number, usageFactor?: number }} o N [N] axial force
 *   (sign ignored); M [N·m] bending moment; D [m]; t [m]; fy [Pa] yield strength; usageFactor [-] allowable
 *   fraction of yield (default 1)
 * @returns {{ sigmaAxial: number, sigmaBending: number, sigmaMax: number, utilisation: number }} stresses [Pa]
 */
export function steelUtilisation({ N, M, D, t, fy, usageFactor = 1 }) {
  pos('D', D); pos('t', t); pos('fy', fy); pos('usageFactor', usageFactor);
  if (t >= D / 2) throw new RangeError('t: must be less than D/2');
  const { A, I } = tubeSection(D, t);
  const sigmaAxial = Math.abs(N) / A, sigmaBending = (Math.abs(M) * D) / (2 * I);
  return { sigmaAxial, sigmaBending, sigmaMax: sigmaAxial + sigmaBending, utilisation: (sigmaAxial + sigmaBending) / (usageFactor * fy) };
}

// ---------------------------------------------------------------------------------------------- axial analysis

function buildAxialModel(o) {
  const { profile, L, D, t, E = 210e9, nElem = 50, mode = 'compression', residualRatio = 0.9, zPeakSand = 0.0025, tensionFactor = 1, setupFactor = 1 } = o;
  pos('L', L); pos('D', D); pos('E', E);
  if (!(Number.isInteger(nElem) && nElem >= 1)) throw new RangeError('nElem: must be an integer >= 1');
  if (mode !== 'compression' && mode !== 'tension') throw new RangeError("mode: must be 'compression' or 'tension'");
  let EA = o.EA;
  if (EA === undefined) { pos('t', t); EA = E * tubeSection(D, t).A; }
  pos('EA', EA);
  const n = nElem + 1, h = L / nElem;
  const z = new Float64Array(n), trib = new Float64Array(n), curves = new Array(n);
  let Qs = 0;
  const fac = setupFactor * (mode === 'tension' ? tensionFactor : 1);
  for (let i = 0; i < n; i++) {
    z[i] = i * h;
    trib[i] = (i === 0 || i === n - 1 ? 0.5 : 1) * h * Math.PI * D; // shaft area served by the node [m²]
    if (o.tzAt) curves[i] = o.tzAt(z[i]);
    else {
      const r = unitResistanceAt(profile, z[i], { betaFactor: o.betaFactor ?? 1 });
      const tmax = r.f * fac;
      const type = propertiesAt(profile, z[i]).type;
      curves[i] = r.type === 'clay' ? tzClay({ tmax, D, residualRatio }) : tzSand({ tmax, zPeak: type === 'rock' ? (o.zPeakRock ?? zPeakSand) : zPeakSand });
    }
    Qs += curves[i].pu * trib[i];
  }
  let tip = null, Qp = 0;
  if (mode === 'compression') {
    if (o.tipCurve) { tip = o.tipCurve; Qp = tip.pu; } else {
      Qp = axialCapacity({ profile, L, D, t, closedEnded: o.closedEnded, plugMode: o.plugMode, internalFrictionFactor: o.internalFrictionFactor, setupFactor, betaFactor: o.betaFactor, dz: o.dz }).Qp;
      tip = qz({ Qp, D });
    }
  }
  return { n, h, z, trib, curves, tip, EA, Qs, Qp, mode };
}

function solveAxialModel(m, wHead, load, w0, tol, maxIter) {
  const { n, h, trib, curves, tip, EA } = m;
  const ke = EA / h;
  const w = new Float64Array(n);
  if (w0) w.set(w0);
  const fixed = wHead !== null;
  if (fixed) w[0] = wHead;
  const R = new Float64Array(n), kd = new Float64Array(n), A = new Float64Array(n * 3), trial = new Float64Array(n);
  let floor = 0;
  const residual = (ww, withK) => {
    let r2 = 0, s2 = 0, a2 = 0;
    for (let i = 0; i < n; i++) {
      let el = 0;
      if (i > 0) el += ke * (ww[i] - ww[i - 1]);
      if (i < n - 1) el += ke * (ww[i] - ww[i + 1]);
      const a = 2 * ke * Math.abs(ww[i]);
      a2 += a * a;
      const c = curves[i];
      let fs = c.p(ww[i]) * trib[i];
      let k = withK ? safeTangent(c, ww[i], c.p(ww[i])) * trib[i] : 0;
      if (i === n - 1 && tip) { const q = tip.p(ww[i]); fs += q; if (withK) k += safeTangent(tip, ww[i], q); }
      const ext = i === 0 && !fixed ? load : 0;
      R[i] = i === 0 && fixed ? 0 : el + fs - ext;
      if (withK) kd[i] = k;
      r2 += R[i] * R[i]; s2 += fs * fs + ext * ext;
    }
    // round-off floor of the bar forces EA/h·(w_i − w_j) for very stiff piles
    if (withK) floor = 256 * Number.EPSILON * Math.sqrt(a2);
    return { r: Math.sqrt(r2), s: Math.sqrt(s2) };
  };
  let res = residual(w, true), it = 0, singular = false;
  const ok = (q) => q.r <= Math.max(tol * q.s, floor, 1e-9);
  const done = () => ok(res);
  for (; it < maxIter && !done(); it++) {
    for (let i = 0; i < n; i++) {
      A[i * 3] = i > 0 ? -ke : 0; A[i * 3 + 2] = i < n - 1 ? -ke : 0;
      A[i * 3 + 1] = (i > 0 ? ke : 0) + (i < n - 1 ? ke : 0) + kd[i];
    }
    const dw = new Float64Array(n);
    for (let i = 0; i < n; i++) dw[i] = -R[i];
    if (fixed) { A[1] = 1; A[2] = 0; A[3] = 0; dw[0] = 0; }
    // all springs at their limit under load control: no equilibrium (load above capacity)
    try { bandSolve(A, dw, n, 1); } catch (err) { singular = true; break; }
    let step = 1;
    const r0 = res.r;
    for (let ls = 0; ls < 12; ls++) {
      for (let i = 0; i < n; i++) trial[i] = w[i] + step * dw[i];
      const nx = residual(trial, false);
      if (nx.r < r0 || ok(nx)) break;
      step *= 0.5;
    }
    w.set(trial);
    res = residual(w, true);
  }
  // axial force at each node (just below the node), positive in the loading sense
  const force = new Float64Array(n), tMob = new Float64Array(n);
  for (let i = 0; i < n - 1; i++) force[i] = ke * (w[i] - w[i + 1]);
  const tipLoad = tip ? tip.p(w[n - 1]) : 0;
  force[n - 1] = tipLoad;
  for (let i = 0; i < n; i++) tMob[i] = curves[i].p(w[i]);
  const headLoad = force[0] + tMob[0] * trib[0];
  return { w, force, tMob, headLoad, tipLoad, converged: !singular && done(), iterations: it, residual: res.s > 0 ? res.r / res.s : res.r };
}

/**
 * Axial load-transfer (t–z / Q–z) analysis of a pile: bar finite elements with nonlinear shaft springs at the
 * nodes and a tip spring. Compression: shaft t–z plus tip Q–z; tension: shaft t–z only (× tensionFactor).
 * Default curves: API clay / sand t–z with tmax from the API main-text unit friction, and the API Q–z with Qp
 * from axialCapacity (for a coring pile the internal friction is part of Qp and is therefore mobilised with the
 * tip curve). Give headDisplacement (displacement control, robust past peak) or headLoad (load control).
 * Reference: API RP 2GEO (2011) §8.3–8.4 (API RP 2A-WSD §6.6–6.7); Coyle & Reese (1966) load-transfer method.
 * Validity: static monotonic loading; pile weight not included.
 * @param {object} o profile; L [m]; D [m]; t [m]; E [Pa] (210e9) or EA [N]; nElem [-] (50); mode
 *   'compression'|'tension'; headDisplacement [m] or headLoad [N] (both positive in the loading sense);
 *   residualRatio [-] clay tres/tmax (0.9); zPeakSand [m] (0.0025); tensionFactor [-] (1); setupFactor [-] (1);
 *   plugMode, closedEnded, internalFrictionFactor, betaFactor, dz as axialCapacity; tzAt(z) → curve [Pa vs m]
 *   and tipCurve [N vs m] override the soil model; tol [-] (1e-10); maxIter [-] (200)
 * @returns {object} { z [m], w [m] displacement, force [N] axial force below each node, tMob [Pa] mobilised
 *   shaft friction — Float64Array per node; headLoad [N]; headDisplacement [m]; tipLoad [N]; tipDisplacement [m];
 *   QsUltimate [N] sum of nodal peak shaft resistance; QpUltimate [N]; convergence { converged, iterations,
 *   residual } }
 */
export function solveAxialPile(o) {
  const m = buildAxialModel(o);
  const hasW = o.headDisplacement !== undefined && o.headDisplacement !== null;
  if (hasW === (o.headLoad !== undefined && o.headLoad !== null)) throw new RangeError('headDisplacement: give exactly one of headDisplacement or headLoad');
  if (hasW) pos('headDisplacement', o.headDisplacement, true); else pos('headLoad', o.headLoad, true);
  const sol = solveAxialModel(m, hasW ? o.headDisplacement : null, hasW ? 0 : o.headLoad, null, o.tol ?? 1e-10, o.maxIter ?? 200);
  if (!sol.converged) throw new Error(hasW ? 'Axial pile analysis did not converge.' : 'Axial pile analysis did not converge: the head load probably exceeds the axial capacity.');
  return {
    z: m.z, w: sol.w, force: sol.force, tMob: sol.tMob, headLoad: sol.headLoad, headDisplacement: sol.w[0], tipLoad: sol.tipLoad,
    tipDisplacement: sol.w[m.n - 1], QsUltimate: m.Qs, QpUltimate: m.Qp,
    convergence: { converged: sol.converged, iterations: sol.iterations, residual: sol.residual },
  };
}

/**
 * Pile-head load–settlement (or load–uplift) curve by displacement control, with the axial capacity taken as
 * the peak head load and the head stiffness as the initial secant.
 * Reference: see solveAxialPile.
 * @param {object} o solveAxialPile options (headDisplacement / headLoad ignored)
 * @param {{ wMax?: number, steps?: number }} [c] wMax [m] largest head displacement (default 0.15 D, enough to
 *   mobilise the Q–z curve fully at z/D = 0.1 plus shaft compression for typical piles); steps [-] (default 60)
 * @returns {{ w: number[], Q: number[], tipLoad: number[], capacity: number, wAtCapacity: number,
 *   headStiffness: number, forceAtCapacity: Float64Array, z: Float64Array, QsUltimate: number,
 *   QpUltimate: number }} w [m], Q [N], capacity [N], headStiffness [N/m], forceAtCapacity [N] vs z [m]
 */
export function axialLoadSettlement(o, c = {}) {
  const m = buildAxialModel(o);
  const wMax = c.wMax ?? 0.15 * o.D, steps = c.steps ?? 60;
  pos('wMax', wMax);
  if (!(steps >= 2)) throw new RangeError('steps: must be >= 2');
  const out = { w: [0], Q: [0], tipLoad: [0], capacity: 0, wAtCapacity: 0, headStiffness: 0, forceAtCapacity: null, z: m.z, QsUltimate: m.Qs, QpUltimate: m.Qp };
  let prev = null;
  for (let s = 1; s <= steps; s++) {
    const wh = wMax * (s / steps) ** 2; // finer steps where the curve is steep
    const sol = solveAxialModel(m, wh, 0, prev, o.tol ?? 1e-10, o.maxIter ?? 200);
    if (!sol.converged) throw new Error(`Axial pile analysis did not converge at head displacement ${wh} m.`);
    prev = sol.w;
    out.w.push(wh); out.Q.push(sol.headLoad); out.tipLoad.push(sol.tipLoad);
    if (s === 1) out.headStiffness = sol.headLoad / wh;
    if (sol.headLoad > out.capacity) { out.capacity = sol.headLoad; out.wAtCapacity = wh; out.forceAtCapacity = sol.force; }
  }
  return out;
}

// ------------------------------------------------------------------------------------- closed-form benchmarks

/**
 * Semi-infinite beam on a linear Winkler foundation loaded at its end by H and M (Hetényi):
 *   β = (k / 4EI)^(1/4); y0 = 2Hβ/k + 2Mβ²/k; θ0 = 2Hβ²/k + 4Mβ³/k;
 *   y(z) = (2β/k) e^(−βz) [H cos βz + Mβ (cos βz − sin βz)];
 *   M(z) = e^(−βz) [(H/β) sin βz + M (cos βz + sin βz)]; under H alone Mmax = 0.3224 H/β at βz = π/4.
 * Reference: Hetényi (1946) Beams on Elastic Foundation, ch. II. Validity: βL > ~5 ("long" pile), constant k.
 * @param {{ H?: number, M?: number, k: number, EI: number, z?: number[] }} o H [N]; M [N·m]; k [N/m²] foundation
 *   modulus (force per length per deflection); EI [N·m²]; z [m] optional depths for the distributions
 * @returns {{ beta: number, y0: number, theta0: number, KHH: number, KHM: number, KMM: number, y?: number[],
 *   theta?: number[], moment?: number[] }} beta [1/m]; y0 [m]; theta0 [rad]; head stiffness KHH = 4EIβ³ [N/m],
 *   KHM = −2EIβ² [N], KMM = 2EIβ [N·m/rad]; distributions at z when requested
 */
export function hetenyiSemiInfinite({ H = 0, M = 0, k, EI, z }) {
  pos('k', k); pos('EI', EI);
  const b = Math.pow(k / (4 * EI), 0.25);
  const out = { beta: b, y0: (2 * H * b) / k + (2 * M * b * b) / k, theta0: (2 * H * b * b) / k + (4 * M * b ** 3) / k, KHH: 4 * EI * b ** 3, KHM: -2 * EI * b * b, KMM: 2 * EI * b };
  if (z) {
    out.y = []; out.theta = []; out.moment = [];
    for (const zz of z) {
      const e = Math.exp(-b * zz), c = Math.cos(b * zz), s = Math.sin(b * zz);
      out.y.push(((2 * b) / k) * e * (H * c + M * b * (c - s)));
      out.theta.push(((2 * b * b) / k) * e * (H * (c + s) + 2 * M * b * c));
      out.moment.push(e * ((H / b) * s + M * (c + s)));
    }
  }
  return out;
}

/**
 * Broms ultimate lateral resistance of a free-head pile in clay (uniform su, zero resistance over the top 1.5 D,
 * 9 su D below). Short (rigid) pile: f = Hu/(9 su D), L = 1.5 D + f + g, Mmax = Hu (e + 1.5 D + 0.5 f) =
 * 2.25 su D g². Long pile: Mmax = My.
 * Reference: Broms (1964a) J. Soil Mech. Found. Div. ASCE 90(SM2). Validity: free head, uniform cohesive soil.
 * @param {{ su: number, D: number, L: number, e?: number, My?: number }} o su [Pa]; D [m]; L [m] embedment;
 *   e [m] load eccentricity above ground (default 0); My [N·m] pile yield moment (optional)
 * @returns {{ HuShort: number, HuLong: number|null, Hu: number, mode: string, f: number, Mmax: number }} Hu [N];
 *   mode 'short' | 'long'; f [m] depth of zero shear below 1.5 D; Mmax [N·m] in the pile at Hu
 */
export function bromsClay({ su, D, L, e = 0, My }) {
  pos('su', su); pos('D', D); pos('L', L); pos('e', e, true);
  if (L <= 1.5 * D) throw new RangeError('L: must exceed 1.5 D');
  const q = 9 * su * D;
  const solve = (fn, hi) => { let a = 0, b = hi; for (let i = 0; i < 200; i++) { const m = 0.5 * (a + b); if (fn(m) > 0) b = m; else a = m; } return 0.5 * (a + b); };
  // short pile: H (e + 1.5D + 0.5 f) = 2.25 su D g² with g = L − 1.5D − f, f = H/q
  const HuShort = solve((H) => { const f = H / q, g = L - 1.5 * D - f; return H * (e + 1.5 * D + 0.5 * f) - 2.25 * su * D * g * g; }, q * (L - 1.5 * D));
  let HuLong = null;
  if (My !== undefined) { pos('My', My); HuLong = solve((H) => H * (e + 1.5 * D + (0.5 * H) / q) - My, Math.max(q * L, My / D) * 10); }
  const long = HuLong !== null && HuLong < HuShort;
  const Hu = long ? HuLong : HuShort;
  const f = Hu / q;
  return { HuShort, HuLong, Hu, mode: long ? 'long' : 'short', f, Mmax: Hu * (e + 1.5 * D + 0.5 * f) };
}

/**
 * Broms ultimate lateral resistance of a free-head pile in sand (soil resistance 3 Kp γ′ z D):
 * short pile Hu = 0.5 γ′ D L³ Kp / (e + L); long pile Hu (e + 2f/3) = My with f = 0.82 sqrt(Hu / (γ′ D Kp)).
 * Reference: Broms (1964b) J. Soil Mech. Found. Div. ASCE 90(SM3). Validity: free head, uniform cohesionless soil.
 * @param {{ gammaSub: number, phiDeg: number, D: number, L: number, e?: number, My?: number }} o gammaSub [N/m³];
 *   phiDeg [deg]; D [m]; L [m]; e [m] (default 0); My [N·m] (optional)
 * @returns {{ HuShort: number, HuLong: number|null, Hu: number, mode: string, Kp: number }} Hu [N]; Kp [-]
 */
export function bromsSand({ gammaSub, phiDeg, D, L, e = 0, My }) {
  pos('gammaSub', gammaSub); pos('phiDeg', phiDeg); pos('D', D); pos('L', L); pos('e', e, true);
  const Kp = Math.tan(Math.PI / 4 + (phiDeg * Math.PI) / 360) ** 2;
  const HuShort = (0.5 * gammaSub * D * L ** 3 * Kp) / (e + L);
  let HuLong = null;
  if (My !== undefined) {
    pos('My', My);
    let a = 0, b = Math.max(HuShort, My / D) * 10;
    const fn = (H) => H * (e + (2 / 3) * 0.82 * Math.sqrt(H / (gammaSub * D * Kp))) - My;
    for (let i = 0; i < 200; i++) { const m = 0.5 * (a + b); if (fn(m) > 0) b = m; else a = m; }
    HuLong = 0.5 * (a + b);
  }
  const long = HuLong !== null && HuLong < HuShort;
  return { HuShort, HuLong, Hu: long ? HuLong : HuShort, mode: long ? 'long' : 'short', Kp };
}

/**
 * Flexible pile in an elastic continuum under head load H and moment M, Randolph (1981):
 *   Ep = EI / (π r0⁴/4); G* = G (1 + 3ν/4); critical length lc = 2 r0 (Ep/Gc)^(2/7);
 *   y0 = [(Ep/Gc)^(1/7) / (ρc Gc)] [0.27 H (lc/2)^−1 + 0.30 M (lc/2)^−2];
 *   θ0 = [(Ep/Gc)^(1/7) / (ρc Gc)] [0.30 H (lc/2)^−2 + 0.80 √ρc M (lc/2)^−3];
 * Gc is the mean G* over lc and ρc = G*(lc/4)/Gc. Here G varies linearly, G = G0 + m z, and lc is found by
 * iteration.
 * Reference: Randolph (1981) Géotechnique 31(2) 247–259. Validity: L >= lc (flexible pile), linear elastic soil.
 * @param {{ H?: number, M?: number, EI: number, D: number, G0: number, Ggrad?: number, nu: number }} o H [N];
 *   M [N·m]; EI [N·m²]; D [m]; G0 [Pa] shear modulus at the surface; Ggrad [Pa/m] (default 0); nu [-]
 * @returns {{ y0: number, theta0: number, lc: number, Gc: number, rhoc: number, KHH: number, KHM: number,
 *   KMM: number }} y0 [m]; theta0 [rad]; lc [m]; Gc [Pa]; rhoc [-]; head stiffness terms [N/m], [N], [N·m/rad]
 */
export function randolphLateral({ H = 0, M = 0, EI, D, G0, Ggrad = 0, nu }) {
  pos('EI', EI); pos('D', D); pos('G0', G0, true); pos('Ggrad', Ggrad, true);
  if (!(nu >= 0 && nu <= 0.5)) throw new RangeError('nu: must be within 0..0.5');
  if (G0 === 0 && Ggrad === 0) throw new RangeError('G0: soil stiffness must be positive');
  const r0 = D / 2, Ep = EI / ((Math.PI * r0 ** 4) / 4), fs = 1 + 0.75 * nu;
  let lc = 10 * D, Gc = 0;
  for (let i = 0; i < 200; i++) {
    Gc = fs * (G0 + 0.5 * Ggrad * lc);
    const next = 2 * r0 * Math.pow(Ep / Gc, 2 / 7);
    if (Math.abs(next - lc) < 1e-12 * lc) { lc = next; break; }
    lc = next;
  }
  Gc = fs * (G0 + 0.5 * Ggrad * lc);
  const rhoc = (fs * (G0 + 0.25 * Ggrad * lc)) / Gc;
  const c = Math.pow(Ep / Gc, 1 / 7) / (rhoc * Gc), a = lc / 2;
  const fHH = (0.27 * c) / a, fHM = (0.3 * c) / (a * a), fMM = (0.8 * Math.sqrt(rhoc) * c) / a ** 3;
  const det = fHH * fMM - fHM * fHM;
  return { y0: fHH * H + fHM * M, theta0: fHM * H + fMM * M, lc, Gc, rhoc, KHH: fMM / det, KHM: -fHM / det, KMM: fHH / det };
}

/**
 * Axial head stiffness of a pile in an elastic soil, Randolph & Wroth (1978):
 *   P/(G_L r0 w) = [4η/((1−ν)ξ) + ρ (2π/ζ) (tanh μL / μL) (L/r0)] / [1 + (1/(πλ)) (4η/((1−ν)ξ)) (tanh μL / μL) (L/r0)]
 * with ζ = ln(rm/r0), rm = 2.5 ρ (1−ν) L, μL = (L/r0) √(2/(ζλ)), λ = Ep/G_L, ρ = G_avg/G_L, η = r_base/r0 = 1,
 * ξ = G_L/G_base.
 * Reference: Randolph & Wroth (1978) J. Geotech. Eng. Div. ASCE 104(GT12); Fleming et al., Piling Engineering.
 * Validity: linear elastic soil with shear modulus constant or increasing linearly with depth; solid-equivalent
 * pile modulus Ep = EA/(π r0²).
 * @param {{ L: number, D: number, EA: number, GL: number, rho?: number, xi?: number, nu: number }} o L [m];
 *   D [m]; EA [N]; GL [Pa] soil shear modulus at pile tip level; rho [-] G_avg/G_L (1 uniform, 0.5 Gibson);
 *   xi [-] G_L/G_base (1 for a floating pile); nu [-]
 * @returns {{ stiffness: number, zeta: number, muL: number }} stiffness [N/m]; zeta [-]; muL [-]
 */
export function randolphWrothAxial({ L, D, EA, GL, rho = 1, xi = 1, nu }) {
  pos('L', L); pos('D', D); pos('EA', EA); pos('GL', GL); pos('rho', rho); pos('xi', xi);
  if (!(nu >= 0 && nu < 0.5 + 1e-12)) throw new RangeError('nu: must be within 0..0.5');
  const r0 = D / 2, lambda = EA / (Math.PI * r0 * r0) / GL;
  const zeta = Math.log((2.5 * rho * (1 - nu) * L) / r0);
  if (!(zeta > 0)) throw new RangeError('L: pile too short for the Randolph–Wroth solution (rm <= r0)');
  const muL = (L / r0) * Math.sqrt(2 / (zeta * lambda));
  const th = (Math.tanh(muL) / muL) * (L / r0);
  const base = 4 / ((1 - nu) * xi);
  return { stiffness: (GL * r0 * (base + ((rho * 2 * Math.PI) / zeta) * th)) / (1 + (base * th) / (Math.PI * lambda)), zeta, muL };
}
