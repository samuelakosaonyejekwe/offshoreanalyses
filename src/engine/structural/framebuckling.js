/**
 * Frame stability: geometric stiffness, linear (eigenvalue) buckling, iterative second-order (P-Δ / P-δ)
 * analysis and effective-length back-calculation.
 *
 * Geometric stiffness of the beam element, N ∫ w′ᵀ w′ dx with the deflection shape functions of the element
 * itself. Without shear flexibility these are the Hermite cubics and the bending block in both planes is
 * N/(30L)·[36 3L −36 3L; 3L 4L² −3L −L²; −36 −3L 36 −3L; 3L −L² −3L 4L²] (Przemieniecki 1968 §15.3; McGuire,
 * Gallagher & Ziemian, Matrix Structural Analysis, 2nd ed., Eq. 9.18 without the axial and bending–torsion coupling
 * terms). With the shear parameter Φ = 12 E I/(G As L²) of the elastic stiffness the deflection also carries the
 * shear angle and the block becomes
 *   N/(30L(1+Φ)²)·[a 3L −a 3L; 3L c −3L −d; −a −3L a −3L; 3L −d −3L c],
 *   a = 36 + 60Φ + 30Φ², c = L²(4 + 5Φ + 2.5Φ²), d = L²(1 + 5Φ + 2.5Φ²),
 * which is the same displacement field as the Timoshenko stiffness matrix (Engesser's shear-flexible column). The
 * pair is then a Rayleigh–Ritz discretisation: buckling factors are upper bounds and fall monotonically when every
 * element is halved. The Hermite matrix on a shear-flexible element is not (it counts the shear angle twice in part
 * and gives factors that move up and down with the mesh by Pcr/(G As) times a few tenths).
 * The Wagner term N (Ip/A)/L acts on the twist. Axial force N is tension positive, so compression softens the member.
 * Because the consistent matrix carries the member-curvature (P-δ) terms, the member bow effect converges with
 * mesh refinement (two to four elements per member are normally sufficient).
 * SI units.
 */

import { eigSym } from '../math/eigen.js';
import {
  buildFrameSystem, factorFrame, solveFrame, resolveLoads, loadParts, elementToGlobal, elementMatVec,
} from './frame3d.js';

const XY = [1, 5, 7, 11];
const XZ = [2, 4, 8, 10];
const XZS = [1, -1, 1, -1];

/**
 * Local 12×12 geometric stiffness matrix of a beam (or bar) carrying axial force N.
 * @param {object} o
 * @param {number} o.N axial force [N], tension positive
 * @param {number} o.L length [m]
 * @param {number} [o.r2=0] polar radius of gyration squared Ip/A [m²] for the torsional (Wagner) term
 * @param {boolean} [o.truss=false] bar element: string stiffness N/L on the transverse translations only
 * @param {number} [o.phiY=0] shear parameter 12 E Iz/(G Asy L²) of the local x–y plane [-]; 0 = Euler–Bernoulli
 * @param {number} [o.phiZ=0] shear parameter 12 E Iy/(G Asz L²) of the local x–z plane [-]
 * @returns {Float64Array} row-major 12×12 matrix [N/m, N, N·m]
 */
export function beamGeometricStiffness({ N, L, r2 = 0, truss = false, phiY = 0, phiZ = 0 }) {
  if (!(L > 0)) throw new RangeError('L must be positive');
  if (!Number.isFinite(N)) throw new RangeError('N must be finite');
  if (!(phiY >= 0) || !(phiZ >= 0)) throw new RangeError('phiY and phiZ must be non-negative');
  const k = new Float64Array(144);
  if (truss) {
    const c = N / L;
    for (const [a, b] of [[1, 7], [2, 8]]) { k[a * 12 + a] = c; k[b * 12 + b] = c; k[a * 12 + b] = -c; k[b * 12 + a] = -c; }
    return k;
  }
  const block = (phi) => {
    const c = N / (30 * L * (1 + phi) ** 2);
    const a = 36 + 60 * phi + 30 * phi * phi;
    const s = L * L * (4 + 5 * phi + 2.5 * phi * phi);
    const d = L * L * (1 + 5 * phi + 2.5 * phi * phi);
    return [
      [c * a, c * 3 * L, -c * a, c * 3 * L],
      [c * 3 * L, c * s, -c * 3 * L, -c * d],
      [-c * a, -c * 3 * L, c * a, -c * 3 * L],
      [c * 3 * L, -c * d, -c * 3 * L, c * s],
    ];
  };
  const by = block(phiY); const bz = block(phiZ);
  for (let i = 0; i < 4; i++) {
    for (let j = 0; j < 4; j++) {
      k[XY[i] * 12 + XY[j]] = by[i][j];
      k[XZ[i] * 12 + XZ[j]] = XZS[i] * XZS[j] * bz[i][j];
    }
  }
  const t = N * r2 / L;
  k[3 * 12 + 3] = t; k[9 * 12 + 9] = t; k[3 * 12 + 9] = -t; k[9 * 12 + 3] = -t;
  return k;
}

/**
 * Geometric stiffness matrices of all members of a system for given axial forces, each with the shear parameters
 * of its own elastic stiffness.
 * @param {object} sys system from buildFrameSystem
 * @param {ArrayLike<number>} axial member axial forces [N], tension positive
 * @param {object} [o]
 * @param {boolean} [o.consistentShear=true] false: Hermite (Euler–Bernoulli) matrix also on shear-flexible elements
 * @returns {{global:Float64Array[], local:Float64Array[]}} per-element matrices in global node dofs and in local dofs
 */
export function frameGeometricStiffness(sys, axial, { consistentShear = true } = {}) {
  if (!axial || axial.length !== sys.elems.length) throw new RangeError('axial must hold one force per member');
  const global = []; const local = [];
  for (let e = 0; e < sys.elems.length; e++) {
    const el = sys.elems[e];
    const kl = beamGeometricStiffness({
      N: axial[e], L: el.L, r2: (el.sec.Ip || el.sec.Iy + el.sec.Iz) / el.sec.A, truss: el.truss,
      phiY: consistentShear ? el.phiY || 0 : 0, phiZ: consistentShear ? el.phiZ || 0 : 0,
    });
    local.push(kl);
    global.push(elementToGlobal(el, kl));
  }
  return { global, local };
}

function memberAxial(res) {
  return res.memberForces.map((mf) => 0.5 * (mf.endA.N + mf.endB.N));
}

/* Solve the reduced problem A z = μ B z with B symmetric positive semi-definite (rank-deficient directions dropped). */
function reducedPencil(A, B, q) {
  const eb = eigSym(B);
  const bmax = eb.values[q - 1];
  const keep = [];
  for (let k = 0; k < q; k++) if (eb.values[k] > 1e-12 * bmax) keep.push(k);
  const r = keep.length;
  // T = W μ^-½ (q × r)
  const T = Array.from({ length: q }, () => new Array(r));
  for (let c = 0; c < r; c++) {
    const s = 1 / Math.sqrt(eb.values[keep[c]]);
    const w = eb.vectors[keep[c]];
    for (let i = 0; i < q; i++) T[i][c] = w[i] * s;
  }
  const G = Array.from({ length: r }, () => new Array(r).fill(0));
  for (let a = 0; a < r; a++) {
    for (let b = a; b < r; b++) {
      let s = 0;
      for (let i = 0; i < q; i++) { let t = 0; for (let j = 0; j < q; j++) t += A[i][j] * T[j][b]; s += T[i][a] * t; }
      G[a][b] = s; G[b][a] = s;
    }
  }
  const eg = eigSym(G);
  // columns ordered by descending μ
  const order = eg.values.map((v, i) => i).sort((i, j) => eg.values[j] - eg.values[i]);
  const mu = order.map((i) => eg.values[i]);
  const Z = order.map((i) => {
    const z = new Array(q).fill(0);
    for (let k = 0; k < q; k++) for (let c = 0; c < r; c++) z[k] += T[k][c] * eg.vectors[i][c];
    return z;
  });
  return { mu, Z };
}

/**
 * Linear (eigenvalue) buckling of a frame: (K + λ Kg(N_ref)) φ = 0 with the member axial forces of a reference
 * load case from a first-order analysis; the lowest positive factors λ are returned.
 *
 * Method. The pencil has factors of both signs: members in tension under the reference load buckle under the
 * reversed load (negative λ), and for slender ties those lie much closer to zero than the wanted positive ones.
 * Plain inverse iteration on (−Kg) φ = (1/λ) K φ converges to the factors of smallest magnitude of either sign and
 * loses the positive ones as soon as there are more reversed-load modes below them than iteration vectors. The
 * iteration is therefore run on the shifted operator
 *     (K + σ Kg)⁻¹ (−Kg) φ = θ φ,   θ = 1/(λ − σ),
 * for which every factor in (σ, 2σ) has θ > 1/σ while every negative factor has |θ| < 1/σ. The shift is found
 * from the inertia of the factorisation: K + σ Kg has as many negative pivots as there are positive factors below
 * σ (Sylvester's law; Bathe 1996 §11.4.3, Sturm sequence check). σ is doubled from an Euler estimate of the most
 * loaded member until a factor not yet found lies below it and then bisected until the bracket [σ0, σ1] with no
 * and at least one such factor satisfies σ1 <= 2 σ0; σ = σ0 then makes the next wanted mode dominant.
 * Subspace iteration follows (Bathe 1996 §11.6), with a Rayleigh–Ritz projection of the original pencil
 * (Zᵀ K Z, Zᵀ (−Kg) Z; K is positive definite, so the Ritz factors are upper bounds); the subspace is enlarged
 * when the Ritz values show a cluster of near-equal factors reaching beyond it.
 * Modes are accepted in ascending order when the factor has stopped changing (tol) and the residual
 *     ρ = ‖K φ − λ (−Kg) φ‖ / ‖K φ‖ in the K⁻¹ (energy) norm
 * is below resTol; the error of an accepted factor is at most ρ λ and of order ρ² λ for a separated one.
 * When wanted factors lie beyond 2σ, where reversed-load modes may outnumber the iteration vectors, the accepted
 * modes are locked, the iteration vectors are kept K-orthogonal to them (deflation; the operator is self-adjoint in
 * the K inner product) and the shift is moved up by the same inertia search, counting only factors not yet found.
 * A group of exactly equal factors wider than the subspace (identical members) converges at the rate set by the
 * next distinct factor and any basis of it is returned. A last inertia count just below the highest returned
 * factor checks that none was skipped (sturmCheck).
 * Effective lengths are back-calculated per compressed member: Pcr = λ1 |N|, Le = π √(E I/Pcr), K = Le/L.
 * Validity: small pre-buckling displacements; tension-only / compression-only members are treated as ordinary bars.
 * @param {object} model frame model with load cases
 * @param {object} o
 * @param {string} o.loadCase name of the reference load case or combination
 * @param {number} [o.nModes=3] number of buckling modes
 * @param {boolean} [o.shearDeformation=true]
 * @param {boolean} [o.consistentShear=true] geometric stiffness from the shear-flexible shape functions (see the file
 *   header): upper-bound factors that fall monotonically under mesh halving, with an error of order
 *   (Pcr/(G As))·(h/L)² from the shear angle taken constant per element. false = Hermite matrix on every element:
 *   fourth-order in h but not a bound and not monotone, and it converges to a factor low by about 0.2 (Pcr/(G As))²
 * @param {number} [o.tol=1e-9] relative tolerance on the change of λ between two sweeps
 * @param {number} [o.resTol=1e-7] tolerance on the relative residual ρ of each mode
 * @param {number} [o.maxIter=200] cap on the total number of sweeps
 * @returns {{factors:number[], modes:Float64Array[], residuals:number[], axial:number[], effectiveLengths:Array<{member:number,
 *   N:number, Pcr:number, Ley:number, Lez:number, Ky:number, Kz:number}>, converged:boolean, iterations:number, shift:number,
 *   shifts:number[], sturmCheck:(boolean|null), system:object, reference:object}} buckling load factors [-] ascending, modes over
 *   all dofs (max component 1), residuals ρ [-] of the returned modes, member axial forces [N], effective lengths [m] for
 *   buckling in the local x–z plane (about y: Ley, Ky) and x–y plane (Lez, Kz); converged = nModes factors (or all the positive
 *   factors there are) meet both tolerances and the inertia count finds none skipped; shift = first σ used [-], shifts = all of
 *   them; sturmCheck = false when the inertia count shows a skipped factor, null when the count could not be taken
 */
export function bucklingAnalysis(model, {
  loadCase, nModes = 3, shearDeformation = true, consistentShear = true, tol = 1e-9, resTol = 1e-7, maxIter = 200,
} = {}) {
  if (typeof loadCase !== 'string') throw new RangeError('loadCase must name a load case or combination');
  if (!Number.isInteger(nModes) || nModes < 1) throw new RangeError('nModes must be a positive integer');
  if (!(tol > 0) || !(resTol > 0)) throw new RangeError('tol and resTol must be positive');
  const sys = buildFrameSystem(model, { shearDeformation });
  const fac0 = factorFrame(sys);
  const ref = solveFrame(sys, fac0, resolveLoads(sys, loadParts(model, loadCase)));
  const axial = memberAxial(ref);
  let nmax = 0;
  for (const v of axial) nmax = Math.max(nmax, Math.abs(v));
  if (!(nmax > 0)) throw new Error('The reference load case produces no member axial forces: nothing can buckle.');
  const NONE = 'No positive buckling factor: the reference load case does not compress any member critically.';
  // Euler estimate of the factor of the most critical compressed member (pinned ends), the scale of the shift search
  let est = Infinity;
  for (let e = 0; e < sys.elems.length; e++) {
    if (!(axial[e] < -1e-12 * nmax)) continue;
    const el = sys.elems[e];
    const I = el.truss ? 0 : Math.min(el.sec.Iy, el.sec.Iz);
    const pe = I > 0 ? Math.PI ** 2 * el.mat.E * I / (el.L * el.L) : el.mat.E * el.sec.A;
    est = Math.min(est, pe / -axial[e]);
  }
  if (!Number.isFinite(est)) throw new Error(NONE);
  const kg = frameGeometricStiffness(sys, axial, { consistentShear }).global;
  const ke = sys.elems.map((el) => el.ke);
  const n = sys.nFree;
  const p = Math.min(nModes, n);
  const dot = (a, b) => { let s = 0; for (let k = 0; k < n; k++) s += a[k] * b[k]; return s; };
  const negKg = (x) => { const y = elementMatVec(sys, kg, x); for (let i = 0; i < n; i++) y[i] = -y[i]; return y; };
  const mulK = (x) => {
    const y = elementMatVec(sys, ke, x);
    for (const s of sys.springs) {
      for (let i = 0; i < 6; i++) {
        const qi = sys.free[6 * s.node + i];
        if (qi < 0) continue;
        for (let j = 0; j < 6; j++) { const qj = sys.free[6 * s.node + j]; if (qj >= 0) y[qi] += s.k[i * 6 + j] * x[qj]; }
      }
    }
    return y;
  };
  // locked modes {lam, x, kx, res}, x K-normalised, ascending
  const found = [];
  // number of positive factors below s (negative pivots of K + s Kg); a singular matrix counts as "at a factor"
  const inertia = (s) => {
    try { const f = factorFrame(sys, { add: kg, addFactor: s }); return { count: f.negativePivots, fac: f }; } catch (e) { return { count: -1, fac: null }; }
  };
  // factors below s that are not locked yet (−1: s is at a factor)
  const missing = (t, s) => {
    if (t.count < 0) return -1;
    let have = 0;
    for (const f of found) if (f.lam < s) have++;
    return Math.max(0, t.count - have);
  };
  // shift bracket above lo: s0 with no missing factor below it and s1 <= 2 s0 with at least one
  const bracket = (lo, start, facLo) => {
    let s0 = lo; let s1 = start; let fac = facLo; let hiFound = false;
    for (let k = 0; k < 60; k++) {
      const t = inertia(s1); const u = missing(t, s1);
      if (u === 0) { s0 = s1; fac = t.fac; s1 *= 2; } else if (u > 0) { hiFound = true; break; } else s1 *= 1.0625;
    }
    if (!hiFound) return null;
    for (let k = 0; k < 60 && !(s0 > 0 && s1 <= 2 * s0); k++) {
      const mid = s0 > 0 ? Math.sqrt(s0 * s1) : 0.5 * s1;
      const t = inertia(mid); const u = missing(t, mid);
      if (u === 0) { s0 = mid; fac = t.fac; } else if (u > 0) s1 = mid; else s1 = mid * (1 - 1e-3) + s0 * 1e-3;
      if (s1 < 1e-9 * est) break;
    }
    return { sigma: s0, fac };
  };
  let br = bracket(0, 0.25 * est, fac0);
  if (!br) throw new Error(NONE);
  const shifts = [br.sigma];
  const qMax = Math.min(n, Math.max(8 * p, p + 32));
  let q = Math.min(n, Math.max(2 * p, p + 8));
  // deterministic start vectors (no randomness needed: smooth + alternating patterns)
  let nStart = 0;
  const startVector = () => {
    const j = nStart++;
    const x = new Float64Array(n);
    for (let i = 0; i < n; i++) x[i] = Math.sin((i + 1) * (j + 1) * 0.7071 + j) + (j === 0 ? 1 : 0);
    return x;
  };
  const combine = (V, z, r) => {
    const x = new Float64Array(n);
    for (let j = 0; j < r; j++) { const zj = z[j]; if (zj === 0) continue; const v = V[j]; for (let k = 0; k < n; k++) x[k] += zj * v[k]; }
    return x;
  };
  let X = [];
  for (let j = 0; j < q; j++) X.push(startVector());
  let converged = false; let exhausted = false; let it = 0;
  let prev = null; let sweeps = 0; let lastGrow = 0; let pending = [];
  while (it < maxIter && !converged && !exhausted) {
    it++; sweeps++;
    const { sigma, fac } = br;
    const need = p - found.length;
    // Z = (K + σ Kg)⁻¹ (−Kg) X, K-orthogonal to the locked modes and scaled to unit K-norm
    const Z = []; const KZ = []; const GZ = [];
    for (let j = 0; j < X.length; j++) {
      const z = fac.solve(negKg(X[j]));
      for (let pass = 0; pass < 2; pass++) {
        for (const f of found) { const c = dot(f.kx, z); for (let k = 0; k < n; k++) z[k] -= c * f.x[k]; }
      }
      const kz = mulK(z);
      const e = dot(z, kz);
      if (!(e > 0) || !Number.isFinite(e)) continue;
      const sc = 1 / Math.sqrt(e);
      for (let k = 0; k < n; k++) { z[k] *= sc; kz[k] *= sc; }
      Z.push(z); KZ.push(kz); GZ.push(negKg(z));
    }
    const r = Z.length;
    const A = Array.from({ length: r }, () => new Array(r).fill(0));
    const B = Array.from({ length: r }, () => new Array(r).fill(0));
    for (let i = 0; i < r; i++) {
      for (let j = i; j < r; j++) {
        A[i][j] = 0.5 * (dot(Z[i], GZ[j]) + dot(Z[j], GZ[i])); A[j][i] = A[i][j];
        B[i][j] = 0.5 * (dot(Z[i], KZ[j]) + dot(Z[j], KZ[i])); B[j][i] = B[i][j];
      }
    }
    // Ritz pairs of (−Kg) y = μ K y, μ = 1/λ descending: the positive factors come first, ascending
    const red = reducedPencil(A, B, r);
    const m = red.mu.length;
    const lam = [];
    for (let c = 0; c < m && red.mu[c] > 0; c++) {
      if (!(red.mu[c] > 1e-14 * red.mu[0])) break;
      lam.push(1 / red.mu[c]);
    }
    const Xn = red.Z.map((z) => combine(Z, z, r));
    // leading modes that have stopped changing and satisfy the residual test
    const lead = Math.min(need, lam.length);
    const res = [];
    for (let c = 0; c < lead; c++) {
      if (!prev || c >= prev.length || Math.abs(lam[c] - prev[c]) > tol * Math.abs(lam[c])) break;
      const kx = combine(KZ, red.Z[c], r); const gx = combine(GZ, red.Z[c], r);
      const rv = new Float64Array(n);
      for (let k = 0; k < n; k++) rv[k] = kx[k] - lam[c] * gx[k];
      const rho = Math.sqrt(Math.max(0, dot(rv, fac0.solve(rv))) / dot(Xn[c], kx));
      if (!(rho <= resTol)) break;
      res.push({ kx, rho });
    }
    const lock = (count) => {
      for (let c = 0; c < count; c++) {
        const sc = 1 / Math.sqrt(dot(Xn[c], res[c].kx));
        const x = Xn[c]; const kx = res[c].kx;
        for (let k = 0; k < n; k++) { x[k] *= sc; kx[k] *= sc; }
        found.push({ lam: lam[c], x, kx, res: res[c].rho });
      }
      Xn.splice(0, count); lam.splice(0, count);
    };
    // |θ| of the Ritz pairs: how far the wanted modes are from the edge of the subspace
    const theta = red.mu.map((mu) => Math.abs(mu / (1 - sigma * mu)));
    const want = lead > 0 ? Math.abs(1 / (lam[lead - 1] - sigma)) : Infinity;
    let ahead = 0; let edge = Infinity;
    for (const t of theta) { if (t > want * (1 + 1e-9)) ahead++; edge = Math.min(edge, t); }
    if (res.length >= need) {
      lock(need); converged = true;
    } else if (res.length > 0 && sweeps >= 3 && (lam.length < need || ahead >= m - 2 || sweeps >= 25)) {
      // the remaining modes are crowded out by reversed-load modes (or slow): lock, deflate and move the shift up
      lock(res.length);
      const top = found[found.length - 1].lam;
      const next = bracket(sigma, Math.max(2 * sigma, 1.02 * top), fac);
      if (!next) exhausted = true;
      else { br = next; shifts.push(next.sigma); }
      prev = null; sweeps = 0; lastGrow = 0;
    } else {
      prev = lam.slice();
      // a cluster of near-equal factors reaching beyond the subspace: enlarge it (see math/eigen.js subspaceIteration)
      if (q < qMax && sweeps >= 2 && sweeps - lastGrow >= 2 && lam.length >= need && m >= q && edge > 0.7 * want) {
        q += Math.min(q, qMax - q); lastGrow = sweeps;
      }
    }
    X = Xn; pending = lam;
    const room = Math.min(q, n - found.length);
    while (X.length > room) X.pop();
    while (X.length < room) X.push(startVector()); // refill dropped and locked directions
    if (room < 1) exhausted = true;
  }
  // sweep cap reached: the best estimates of the missing modes are returned with their residuals, converged = false
  for (let c = 0; !converged && !exhausted && found.length < p && c < pending.length && c < X.length; c++) {
    const x = X[c]; const kx = mulK(x); const gx = negKg(x);
    const e = dot(x, kx);
    const rv = new Float64Array(n);
    for (let k = 0; k < n; k++) rv[k] = kx[k] - pending[c] * gx[k];
    found.push({ lam: pending[c], x, kx, res: Math.sqrt(Math.max(0, dot(rv, fac0.solve(rv))) / e) });
  }
  if (found.length === 0) throw new Error(NONE);
  found.sort((a, b) => a.lam - b.lam);
  // Sturm check: the factors below the highest returned one must all be among the returned ones
  let sturmCheck = null;
  {
    const s = found[found.length - 1].lam * (1 - 1e-6);
    const t = inertia(s);
    if (t.count >= 0) { let have = 0; for (const f of found) if (f.lam < s) have++; sturmCheck = t.count === have; }
  }
  const factors = []; const modes = []; const residuals = [];
  for (const f of found) {
    factors.push(f.lam); residuals.push(f.res);
    const full = new Float64Array(sys.nDof);
    let big = 0;
    for (let k = 0; k < n; k++) if (Math.abs(f.x[k]) > Math.abs(big)) big = f.x[k];
    for (let k = 0; k < n; k++) full[sys.eqDof[k]] = f.x[k] / big;
    modes.push(full);
  }
  const effectiveLengths = [];
  for (let e = 0; e < sys.elems.length; e++) {
    if (!(axial[e] < 0)) continue;
    const el = sys.elems[e];
    const Pcr = factors[0] * -axial[e];
    const Ley = Math.PI * Math.sqrt(el.mat.E * el.sec.Iy / Pcr);
    const Lez = Math.PI * Math.sqrt(el.mat.E * el.sec.Iz / Pcr);
    effectiveLengths.push({ member: e, N: axial[e], Pcr, Ley, Lez, Ky: Ley / el.L, Kz: Lez / el.L });
  }
  return {
    factors, modes, residuals, axial, effectiveLengths, converged: (converged || exhausted) && sturmCheck !== false,
    iterations: it, shift: shifts[0], shifts, sturmCheck, system: sys, reference: ref,
  };
}

/**
 * Second-order elastic analysis by iteration on the geometric stiffness: solve (K + Kg(N)) u = F, update the
 * member axial forces and repeat until they stop changing (P-Δ from the node translations and P-δ from the
 * consistent geometric stiffness of each element). Member end forces include the geometric-stiffness
 * contribution, i.e. the amplified moments.
 * Validity: moderate displacements (the geometry itself is not updated); the load level must be below the elastic
 * critical load — a zero pivot is reported as instability, and a load beyond the critical load is detected from the
 * negative pivots of the factorisation being rejected by the equilibrium check or from non-convergence.
 * @param {object} model frame model
 * @param {object} o
 * @param {string} o.loadCase load case or combination name
 * @param {number} [o.loadFactor=1] multiplier on the whole load set
 * @param {number} [o.tol=1e-8] relative tolerance on member axial forces
 * @param {number} [o.maxIter=30]
 * @param {boolean} [o.shearDeformation=true]
 * @param {'x'|'y'|'z'} [o.verticalAxis='z'] global axis of gravity; sway is the translation normal to it
 * @returns {{result:object, firstOrder:object, iterations:number, converged:boolean, amplification:number,
 *   swayAmplification:number, sway:{firstOrder:number, secondOrder:number, node:number}, system:object}} result and
 *   firstOrder as solveFrame output; amplification = largest nodal translation of any direction, second-order /
 *   first-order — under gravity-dominated loading this is governed by vertical sag, which P-Δ hardly changes, so it
 *   understates the sway effect; swayAmplification = lateral (normal to verticalAxis) translation second-order /
 *   first-order at the node with the largest first-order lateral translation, the quantity comparable to the
 *   storey amplifier 1/(1 − 1/αcr) [-] (1 when the first-order lateral translation is negligible, below 1e-9 of the
 *   largest translation); sway = those two lateral translations [m] and the node index
 */
export function pDeltaAnalysis(model, { loadCase, loadFactor = 1, tol = 1e-8, maxIter = 30, shearDeformation = true, verticalAxis = 'z' } = {}) {
  if (typeof loadCase !== 'string') throw new RangeError('loadCase must name a load case or combination');
  const va = ['x', 'y', 'z'].indexOf(verticalAxis);
  if (va < 0) throw new RangeError("verticalAxis must be 'x', 'y' or 'z'");
  const sys = buildFrameSystem(model, { shearDeformation });
  const parts = loadParts(model, loadCase).map((p) => ({ load: p.load, factor: p.factor * loadFactor }));
  const loads = resolveLoads(sys, parts);
  const first = solveFrame(sys, factorFrame(sys), loads);
  let axial = memberAxial(first);
  let res = first; let converged = false; let it = 0;
  while (it < maxIter) {
    it++;
    const kg = frameGeometricStiffness(sys, axial);
    let fac;
    try {
      fac = factorFrame(sys, { add: kg.global });
    } catch (e) {
      throw new Error(`Second-order analysis: the structure is unstable at this load level (${e.message})`);
    }
    res = solveFrame(sys, fac, loads, { add: kg.global, addLocal: kg.local });
    const next = memberAxial(res);
    let dmax = 0; let nmax = 0;
    for (let e = 0; e < next.length; e++) { dmax = Math.max(dmax, Math.abs(next[e] - axial[e])); nmax = Math.max(nmax, Math.abs(next[e])); }
    axial = next;
    if (dmax <= tol * Math.max(nmax, 1e-30)) { converged = true; break; }
  }
  if (!converged) throw new Error('Second-order analysis did not converge: the load is probably above the elastic critical load.');
  const h1 = (va + 1) % 3; const h2 = (va + 2) % 3;
  const lateral = (r, i) => Math.hypot(r.displacements[6 * i + h1], r.displacements[6 * i + h2]);
  let node = 0; let s1 = 0;
  for (let i = 0; i < sys.nNode; i++) { const v = lateral(first, i); if (v > s1) { s1 = v; node = i; } }
  const s2 = lateral(res, node);
  return {
    result: res, firstOrder: first, iterations: it, converged,
    amplification: first.maxDisplacement > 0 ? res.maxDisplacement / first.maxDisplacement : 1,
    swayAmplification: s1 > 1e-9 * first.maxDisplacement && s1 > 0 ? s2 / s1 : 1,
    sway: { firstOrder: s1, secondOrder: s2, node }, system: sys,
  };
}

/**
 * Effective length of a member from its elastic critical force: Le = π √(E I/Pcr), K = Le/L (Euler).
 * @param {object} o
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.I second moment of area for the buckling plane [m⁴]
 * @param {number} o.Pcr elastic critical axial force of the member [N]
 * @param {number} o.L member length [m]
 * @returns {{Le:number, K:number}} effective length [m] and factor [-]
 */
export function effectiveLength({ E, I, Pcr, L }) {
  for (const [k, v] of Object.entries({ E, I, Pcr, L })) if (!(v > 0)) throw new RangeError(`${k} must be positive`);
  const Le = Math.PI * Math.sqrt(E * I / Pcr);
  return { Le, K: Le / L };
}
