/**
 * Frame stability: geometric stiffness, linear (eigenvalue) buckling, iterative second-order (P-Δ / P-δ)
 * analysis and effective-length back-calculation.
 *
 * Geometric stiffness of the cubic beam element (Przemieniecki 1968 §15.3; McGuire, Gallagher & Ziemian,
 * Matrix Structural Analysis, 2nd ed., Eq. 9.18 without the axial and bending–torsion coupling terms):
 * bending block N/(30L)·[36 3L −36 3L; 3L 4L² −3L −L²; −36 −3L 36 −3L; 3L −L² −3L 4L²] in both planes and the
 * Wagner term N (Ip/A)/L on the twist. Axial force N is tension positive, so compression softens the member.
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
 * @returns {Float64Array} row-major 12×12 matrix [N/m, N, N·m]
 */
export function beamGeometricStiffness({ N, L, r2 = 0, truss = false }) {
  if (!(L > 0)) throw new RangeError('L must be positive');
  if (!Number.isFinite(N)) throw new RangeError('N must be finite');
  const k = new Float64Array(144);
  if (truss) {
    const c = N / L;
    for (const [a, b] of [[1, 7], [2, 8]]) { k[a * 12 + a] = c; k[b * 12 + b] = c; k[a * 12 + b] = -c; k[b * 12 + a] = -c; }
    return k;
  }
  const c = N / (30 * L);
  const b4 = [
    [36, 3 * L, -36, 3 * L],
    [3 * L, 4 * L * L, -3 * L, -L * L],
    [-36, -3 * L, 36, -3 * L],
    [3 * L, -L * L, -3 * L, 4 * L * L],
  ];
  for (let i = 0; i < 4; i++) {
    for (let j = 0; j < 4; j++) {
      k[XY[i] * 12 + XY[j]] = c * b4[i][j];
      k[XZ[i] * 12 + XZ[j]] = c * XZS[i] * XZS[j] * b4[i][j];
    }
  }
  const t = N * r2 / L;
  k[3 * 12 + 3] = t; k[9 * 12 + 9] = t; k[3 * 12 + 9] = -t; k[9 * 12 + 3] = -t;
  return k;
}

/**
 * Geometric stiffness matrices of all members of a system for given axial forces.
 * @param {object} sys system from buildFrameSystem
 * @param {ArrayLike<number>} axial member axial forces [N], tension positive
 * @returns {{global:Float64Array[], local:Float64Array[]}} per-element matrices in global node dofs and in local dofs
 */
export function frameGeometricStiffness(sys, axial) {
  if (!axial || axial.length !== sys.elems.length) throw new RangeError('axial must hold one force per member');
  const global = []; const local = [];
  for (let e = 0; e < sys.elems.length; e++) {
    const el = sys.elems[e];
    const kl = beamGeometricStiffness({ N: axial[e], L: el.L, r2: (el.sec.Ip || el.sec.Iy + el.sec.Iz) / el.sec.A, truss: el.truss });
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
 * load case from a first-order analysis. The lowest positive factors λ are found by subspace iteration on the
 * inverse pencil (−Kg) φ = (1/λ) K φ (Bathe 1996 §11.6 applied with K as the definite matrix).
 * Effective lengths are back-calculated per compressed member: Pcr = λ1 |N|, Le = π √(E I/Pcr), K = Le/L.
 * Validity: small pre-buckling displacements; tension-only / compression-only members are treated as ordinary bars.
 * @param {object} model frame model with load cases
 * @param {object} o
 * @param {string} o.loadCase name of the reference load case or combination
 * @param {number} [o.nModes=3] number of buckling modes
 * @param {boolean} [o.shearDeformation=true]
 * @param {number} [o.tol=1e-9] relative tolerance on 1/λ
 * @param {number} [o.maxIter=200]
 * @returns {{factors:number[], modes:Float64Array[], axial:number[], effectiveLengths:Array<{member:number, N:number, Pcr:number, Ley:number,
 *   Lez:number, Ky:number, Kz:number}>, converged:boolean, iterations:number, system:object, reference:object}} buckling load factors
 *   [-] ascending, modes over all dofs (max component 1), member axial forces [N], effective lengths [m] for buckling in the local x–z
 *   plane (about y: Ley, Ky) and x–y plane (Lez, Kz)
 */
export function bucklingAnalysis(model, { loadCase, nModes = 3, shearDeformation = true, tol = 1e-9, maxIter = 200 } = {}) {
  if (typeof loadCase !== 'string') throw new RangeError('loadCase must name a load case or combination');
  if (!Number.isInteger(nModes) || nModes < 1) throw new RangeError('nModes must be a positive integer');
  const sys = buildFrameSystem(model, { shearDeformation });
  const fac = factorFrame(sys);
  const ref = solveFrame(sys, fac, resolveLoads(sys, loadParts(model, loadCase)));
  const axial = memberAxial(ref);
  let nmax = 0;
  for (const v of axial) nmax = Math.max(nmax, Math.abs(v));
  if (!(nmax > 0)) throw new Error('The reference load case produces no member axial forces: nothing can buckle.');
  const kg = frameGeometricStiffness(sys, axial).global;
  const n = sys.nFree;
  const p = Math.min(nModes, n);
  const q = Math.min(n, Math.max(2 * p, p + 8));
  // deterministic start vectors (no randomness needed: smooth + alternating patterns)
  let X = [];
  for (let j = 0; j < q; j++) {
    const x = new Float64Array(n);
    for (let i = 0; i < n; i++) x[i] = Math.sin((i + 1) * (j + 1) * 0.7071 + j) + (j === 0 ? 1 : 0);
    X.push(x);
  }
  const negKg = (x) => { const y = elementMatVec(sys, kg, x); for (let i = 0; i < n; i++) y[i] = -y[i]; return y; };
  let prev = null; let mu = null; let converged = false; let it = 0;
  while (it < maxIter) {
    it++;
    const Y = X.map(negKg);
    const Xb = Y.map((y) => fac.solve(y));
    const W = Xb.map(negKg);
    const A = Array.from({ length: q }, () => new Array(q).fill(0));
    const B = Array.from({ length: q }, () => new Array(q).fill(0));
    for (let i = 0; i < q; i++) {
      for (let j = i; j < q; j++) {
        let a = 0; let b = 0;
        for (let k = 0; k < n; k++) { a += Xb[i][k] * W[j][k]; b += Xb[i][k] * Y[j][k]; }
        A[i][j] = a; A[j][i] = a; B[i][j] = b; B[j][i] = b;
      }
    }
    const red = reducedPencil(A, B, q);
    mu = red.mu;
    const Xn = [];
    for (let c = 0; c < q; c++) {
      const x = new Float64Array(n);
      if (c < red.Z.length) {
        const z = red.Z[c];
        for (let j = 0; j < q; j++) { const zj = z[j]; if (zj === 0) continue; const xb = Xb[j]; for (let k = 0; k < n; k++) x[k] += zj * xb[k]; }
      } else {
        // refill dropped directions
        for (let k = 0; k < n; k++) x[k] = Math.cos((k + 1) * (c + 2) * 1.618 + it);
      }
      let nr = 0;
      for (let k = 0; k < n; k++) nr = Math.max(nr, Math.abs(x[k]));
      if (nr > 0) for (let k = 0; k < n; k++) x[k] /= nr;
      Xn.push(x);
    }
    X = Xn;
    if (prev) {
      let ok = true;
      for (let c = 0; c < p; c++) {
        if (!(mu[c] > 0)) break;
        if (Math.abs(mu[c] - prev[c]) > tol * Math.abs(mu[c])) { ok = false; break; }
      }
      if (ok) { converged = true; break; }
    }
    prev = mu.slice();
  }
  const factors = []; const modes = [];
  for (let c = 0; c < Math.min(p, mu.length); c++) {
    if (!(mu[c] > 1e-14 * Math.abs(mu[0]))) break;
    factors.push(1 / mu[c]);
    const full = new Float64Array(sys.nDof);
    let big = 0;
    for (let k = 0; k < n; k++) if (Math.abs(X[c][k]) > Math.abs(big)) big = X[c][k];
    for (let k = 0; k < n; k++) full[sys.eqDof[k]] = X[c][k] / big;
    modes.push(full);
  }
  if (factors.length === 0) throw new Error('No positive buckling factor: the reference load case does not compress any member critically.');
  const effectiveLengths = [];
  for (let e = 0; e < sys.elems.length; e++) {
    if (!(axial[e] < 0)) continue;
    const el = sys.elems[e];
    const Pcr = factors[0] * -axial[e];
    const Ley = Math.PI * Math.sqrt(el.mat.E * el.sec.Iy / Pcr);
    const Lez = Math.PI * Math.sqrt(el.mat.E * el.sec.Iz / Pcr);
    effectiveLengths.push({ member: e, N: axial[e], Pcr, Ley, Lez, Ky: Ley / el.L, Kz: Lez / el.L });
  }
  return { factors, modes, axial, effectiveLengths, converged, iterations: it, system: sys, reference: ref };
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
 * @returns {{result:object, firstOrder:object, iterations:number, converged:boolean, amplification:number, system:object}} result and
 *   firstOrder as solveFrame output; amplification = max |u| second-order / first-order
 */
export function pDeltaAnalysis(model, { loadCase, loadFactor = 1, tol = 1e-8, maxIter = 30, shearDeformation = true } = {}) {
  if (typeof loadCase !== 'string') throw new RangeError('loadCase must name a load case or combination');
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
  return {
    result: res, firstOrder: first, iterations: it, converged,
    amplification: first.maxDisplacement > 0 ? res.maxDisplacement / first.maxDisplacement : 1, system: sys,
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
