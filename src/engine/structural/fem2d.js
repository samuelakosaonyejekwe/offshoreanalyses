/**
 * Two-dimensional continuum, plate and facet-shell finite elements.
 *
 * Contents
 *  - mesh generators (mapped blocks, rectangle, annulus, tube wall, quarter plate with hole, stiffened panel)
 *  - plane stress / plane strain / axisymmetric elasticity with T3, Q4, T6, Q8 elements
 *  - steady heat conduction on the same meshes
 *  - Mindlin–Reissner plate bending (MITC4): statics, linear buckling, free vibration
 *  - flat facet shell (Q4 membrane with incompatible modes + MITC4 plate + drilling stabilisation)
 *  - nodal stress recovery, Zienkiewicz–Zhu error indicator, uniform refinement
 *  - a small generic assembled-system driver (`feSystem`) shared with the solid elements
 *
 * Mesh object used throughout this file:
 *   { type: 'T3'|'Q4'|'T6'|'Q8', nen, nn, ne, nodes: Float64Array(2 nn) [m], conn: Int32Array(nen ne),
 *     edges?: { name: number[][] }, nodeSets?: { name: Int32Array } }
 * Node order: corners counter-clockwise, then mid-side nodes (edge 0-1, 1-2, 2-3/2-0, 3-0).
 * A boundary edge is [n1, n2] (linear) or [n1, n2, nMid] (quadratic), oriented with the material on its left.
 *
 * References: Zienkiewicz, Taylor & Zhu, "The Finite Element Method" (6th ed.) vol. 1 ch. 4–6, 10, 14;
 * Bathe & Dvorkin (1985) Int. J. Numer. Meth. Engng 21, 367–383 (MITC4);
 * Taylor, Beresford & Wilson (1976) Int. J. Numer. Meth. Engng 10, 1211–1219 (incompatible modes, QM6);
 * Zienkiewicz & Zhu (1987) Int. J. Numer. Meth. Engng 24, 337–357 (error estimator).
 */
import { CooBuilder, pcg, SkylineSym, reverseCuthillMcKee } from '../math/sparse.js';
import { eigGenSym } from '../math/eigen.js';

const S3 = 1 / Math.sqrt(3);
const S6 = Math.sqrt(0.6);

/* ------------------------------------------------------------------------------------------------------------
 * Small numerical helpers
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Gauss–Legendre points and weights on [-1, 1] for n = 1..4 (exact closed-form/standard values).
 * @param {number} n number of points (1–4)
 * @returns {{x: number[], w: number[]}} abscissae [-] and weights [-]
 */
export function gauss1d(n) {
  switch (n) {
    case 1: return { x: [0], w: [2] };
    case 2: return { x: [-S3, S3], w: [1, 1] };
    case 3: return { x: [-S6, 0, S6], w: [5 / 9, 8 / 9, 5 / 9] };
    case 4: {
      const a = Math.sqrt(3 / 7 - (2 / 7) * Math.sqrt(6 / 5));
      const b = Math.sqrt(3 / 7 + (2 / 7) * Math.sqrt(6 / 5));
      const wa = (18 + Math.sqrt(30)) / 36;
      const wb = (18 - Math.sqrt(30)) / 36;
      return { x: [-b, -a, a, b], w: [wb, wa, wa, wb] };
    }
    default: throw new RangeError('gauss1d: n must be 1..4');
  }
}

/**
 * In-place inverse of a small dense matrix stored row-major (Gauss–Jordan, partial pivoting).
 * @param {Float64Array} a n×n matrix, overwritten by its inverse
 * @param {number} n order
 * @returns {Float64Array} the same array
 */
export function invertSmall(a, n) {
  const piv = new Int32Array(n);
  for (let c = 0; c < n; c++) {
    let p = c;
    let best = Math.abs(a[c * n + c]);
    for (let r = c + 1; r < n; r++) {
      const v = Math.abs(a[r * n + c]);
      if (v > best) { best = v; p = r; }
    }
    if (best === 0) throw new Error('invertSmall: singular matrix');
    piv[c] = p;
    if (p !== c) for (let k = 0; k < n; k++) { const t = a[c * n + k]; a[c * n + k] = a[p * n + k]; a[p * n + k] = t; }
    const d = 1 / a[c * n + c];
    a[c * n + c] = 1;
    for (let k = 0; k < n; k++) a[c * n + k] *= d;
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const f = a[r * n + c];
      if (f === 0) continue;
      a[r * n + c] = 0;
      for (let k = 0; k < n; k++) a[r * n + k] -= f * a[c * n + k];
    }
  }
  for (let c = n - 1; c >= 0; c--) {
    const p = piv[c];
    if (p !== c) for (let r = 0; r < n; r++) { const t = a[r * n + c]; a[r * n + c] = a[r * n + p]; a[r * n + p] = t; }
  }
  return a;
}

function quadRule(n) {
  const g = gauss1d(n);
  const xi = [];
  const eta = [];
  const w = [];
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { xi.push(g.x[i]); eta.push(g.x[j]); w.push(g.w[i] * g.w[j]); }
  return { xi, eta, w };
}

/**
 * Quadrature rules on the unit triangle (weights sum to 1/2): 1 point (exact to degree 1), 3 points (degree 2)
 * or 7 points (degree 5; Radon's rule).
 * @param {number} n number of points: 1, 3 or 7
 * @returns {{xi: number[], eta: number[], w: number[]}} natural coordinates and weights [-]
 */
export function gaussTriangle(n) {
  if (n !== 1 && n !== 3 && n !== 7) throw new RangeError('gaussTriangle: n must be 1, 3 or 7');
  return triRule(n);
}

function triRule(n) {
  if (n === 1) return { xi: [1 / 3], eta: [1 / 3], w: [0.5] };
  if (n === 3) return { xi: [1 / 6, 2 / 3, 1 / 6], eta: [1 / 6, 1 / 6, 2 / 3], w: [1 / 6, 1 / 6, 1 / 6] };
  const s15 = Math.sqrt(15);
  const a1 = (6 - s15) / 21;
  const a2 = (6 + s15) / 21;
  const w1 = (155 - s15) / 2400;
  const w2 = (155 + s15) / 2400;
  return {
    xi: [1 / 3, a1, 1 - 2 * a1, a1, a2, 1 - 2 * a2, a2],
    eta: [1 / 3, a1, a1, 1 - 2 * a1, a2, a2, 1 - 2 * a2],
    w: [9 / 80, w1, w1, w1, w2, w2, w2],
  };
}

const QUAD_XI = [-1, 1, 1, -1, 0, 1, 0, -1];
const QUAD_ETA = [-1, -1, 1, 1, -1, 0, 1, 0];

const ELEM = {
  T3: { nen: 3, tri: true, nc: 3, rule: triRule(1), mass: triRule(3) },
  Q4: { nen: 4, tri: false, nc: 4, rule: quadRule(2), mass: quadRule(2) },
  T6: { nen: 6, tri: true, nc: 3, rule: triRule(7), mass: triRule(7) },
  Q8: { nen: 8, tri: false, nc: 4, rule: quadRule(3), mass: quadRule(3) },
};

/**
 * Shape functions and natural derivatives of the 2-D elements.
 * @param {'T3'|'Q4'|'T6'|'Q8'} type element type
 * @param {number} xi natural coordinate [-]
 * @param {number} eta natural coordinate [-]
 * @param {Float64Array|number[]} N out: shape functions
 * @param {Float64Array|number[]} dxi out: dN/dxi
 * @param {Float64Array|number[]} deta out: dN/deta
 */
export function shape2d(type, xi, eta, N, dxi, deta) {
  if (type === 'Q4') {
    for (let i = 0; i < 4; i++) {
      const a = QUAD_XI[i];
      const b = QUAD_ETA[i];
      N[i] = 0.25 * (1 + a * xi) * (1 + b * eta);
      dxi[i] = 0.25 * a * (1 + b * eta);
      deta[i] = 0.25 * b * (1 + a * xi);
    }
  } else if (type === 'T3') {
    N[0] = 1 - xi - eta; N[1] = xi; N[2] = eta;
    dxi[0] = -1; dxi[1] = 1; dxi[2] = 0;
    deta[0] = -1; deta[1] = 0; deta[2] = 1;
  } else if (type === 'T6') {
    const l1 = 1 - xi - eta;
    N[0] = l1 * (2 * l1 - 1); N[1] = xi * (2 * xi - 1); N[2] = eta * (2 * eta - 1);
    N[3] = 4 * l1 * xi; N[4] = 4 * xi * eta; N[5] = 4 * eta * l1;
    dxi[0] = 1 - 4 * l1; dxi[1] = 4 * xi - 1; dxi[2] = 0;
    dxi[3] = 4 * (l1 - xi); dxi[4] = 4 * eta; dxi[5] = -4 * eta;
    deta[0] = 1 - 4 * l1; deta[1] = 0; deta[2] = 4 * eta - 1;
    deta[3] = -4 * xi; deta[4] = 4 * xi; deta[5] = 4 * (l1 - eta);
  } else if (type === 'Q8') {
    for (let i = 0; i < 4; i++) {
      const a = QUAD_XI[i];
      const b = QUAD_ETA[i];
      N[i] = 0.25 * (1 + a * xi) * (1 + b * eta) * (a * xi + b * eta - 1);
      dxi[i] = 0.25 * a * (1 + b * eta) * (2 * a * xi + b * eta);
      deta[i] = 0.25 * b * (1 + a * xi) * (a * xi + 2 * b * eta);
    }
    for (let i = 4; i < 8; i++) {
      const a = QUAD_XI[i];
      const b = QUAD_ETA[i];
      if (a === 0) {
        N[i] = 0.5 * (1 - xi * xi) * (1 + b * eta);
        dxi[i] = -xi * (1 + b * eta);
        deta[i] = 0.5 * b * (1 - xi * xi);
      } else {
        N[i] = 0.5 * (1 + a * xi) * (1 - eta * eta);
        dxi[i] = 0.5 * a * (1 - eta * eta);
        deta[i] = -eta * (1 + a * xi);
      }
    }
  } else {
    throw new RangeError(`shape2d: unknown element type ${type}`);
  }
}

function checkMesh(mesh, allowed, fn) {
  if (!mesh || !mesh.nodes || !mesh.conn || !ELEM[mesh.type]) throw new RangeError(`${fn}: mesh is not a valid 2-D mesh`);
  if (allowed && !allowed.includes(mesh.type)) throw new RangeError(`${fn}: mesh type ${mesh.type} not supported (use ${allowed.join('/')})`);
}

function positive(v, name, fn) {
  if (!(typeof v === 'number' && v > 0 && Number.isFinite(v))) throw new RangeError(`${fn}: ${name} must be a positive number`);
}

function checkNu(nu, fn) {
  if (!(typeof nu === 'number' && nu > -1 && nu < 0.5)) throw new RangeError(`${fn}: nu must lie in (-1, 0.5)`);
}

/* ------------------------------------------------------------------------------------------------------------
 * Generic assembled system
 * ---------------------------------------------------------------------------------------------------------- */

function rcmOrder(nNodes, nElem, elemNodes) {
  const seen = new Set();
  const edges = [];
  for (let e = 0; e < nElem; e++) {
    const nd = elemNodes(e);
    for (let a = 0; a < nd.length; a++) {
      for (let b = a + 1; b < nd.length; b++) {
        const i = Math.min(nd[a], nd[b]);
        const j = Math.max(nd[a], nd[b]);
        const key = i * nNodes + j;
        if (!seen.has(key)) { seen.add(key); edges.push([i, j]); }
      }
    }
  }
  const inv = new Int32Array(nNodes);
  let ok = false;
  try {
    const perm = reverseCuthillMcKee(nNodes, edges);
    if (perm && perm.length === nNodes) {
      inv.fill(-1);
      ok = true;
      for (let k = 0; k < nNodes; k++) {
        const o = perm[k];
        if (!(o >= 0 && o < nNodes) || inv[o] !== -1) { ok = false; break; }
        inv[o] = k;
      }
    }
  } catch { ok = false; }
  if (!ok) for (let k = 0; k < nNodes; k++) inv[k] = k;
  return inv;
}

/**
 * Assemble a symmetric finite-element system once and solve it for one or more right-hand sides.
 * Dirichlet conditions are imposed by elimination at assembly (unit diagonal on prescribed dofs), so the reduced
 * matrix stays positive definite and reactions are recovered from the unconstrained operator.
 *
 * @param {object} o
 * @param {number} o.nNodes number of nodes
 * @param {number} o.dpn degrees of freedom per node (dof = node*dpn + k)
 * @param {number} o.nElem number of elements (may include edge/spring "elements" of any node count)
 * @param {(e: number) => ArrayLike<number>} o.elemNodes node list of element e
 * @param {(e: number) => Float64Array} o.elemMatrix row-major m×m stiffness (m = nodes·dpn); the buffer may be reused
 * @param {Float64Array} o.fixed length nNodes·dpn; NaN = free, otherwise the prescribed value
 * @param {Float64Array} [o.springs] diagonal grounded-spring stiffness per dof
 * @param {'auto'|'skyline'|'pcg'} [o.solver] direct skyline LDLᵀ or preconditioned conjugate gradients
 * @param {number} [o.tol] PCG relative residual tolerance (default 1e-10)
 * @param {number} [o.maxIter] PCG iteration cap
 * @param {'jacobi'|'ssor'|'ic0'} [o.precond] PCG preconditioner (default 'ic0')
 * @returns {{ndof: number, solver: string, isFixed: Uint8Array, lastIterations: number,
 *   solve: (f: Float64Array) => Float64Array, solveHomogeneous: (b: Float64Array) => Float64Array,
 *   mulK: (x: Float64Array, y?: Float64Array) => Float64Array}}
 *   solve(f) honours the prescribed values; solveHomogeneous(b) solves with zero prescribed values (b on fixed
 *   dofs ignored) for eigen-iterations; mulK applies the unconstrained operator (elements + springs).
 */
export function feSystem(o) {
  const { nNodes, dpn, nElem, elemNodes, elemMatrix, fixed, springs = null, tol = 1e-10, maxIter, precond = 'ic0' } = o;
  const ndof = nNodes * dpn;
  if (fixed.length !== ndof) throw new RangeError('feSystem: fixed must have nNodes*dpn entries');
  const isFixed = new Uint8Array(ndof);
  let nonzeroBc = false;
  let nFixed = 0;
  for (let i = 0; i < ndof; i++) {
    if (!Number.isNaN(fixed[i])) { isFixed[i] = 1; nFixed++; if (fixed[i] !== 0) nonzeroBc = true; }
  }
  let solver = o.solver || 'auto';
  let order = null;
  let firstRow = null;
  const profileOf = () => {
    order = rcmOrder(nNodes, nElem, elemNodes);
    firstRow = new Int32Array(ndof);
    for (let i = 0; i < ndof; i++) firstRow[i] = i;
    for (let e = 0; e < nElem; e++) {
      const nd = elemNodes(e);
      let mn = ndof;
      for (let a = 0; a < nd.length; a++) mn = Math.min(mn, order[nd[a]] * dpn);
      for (let a = 0; a < nd.length; a++) {
        const base = order[nd[a]] * dpn;
        for (let k = 0; k < dpn; k++) if (firstRow[base + k] > mn) firstRow[base + k] = mn;
      }
    }
    let prof = 0;
    let work = 0;
    for (let i = 0; i < ndof; i++) { const hgt = i - firstRow[i]; prof += hgt + 1; work += hgt * hgt; }
    return { prof, work };
  };
  if (solver === 'auto') {
    // A direct factor is preferred while its profile stays within ~200 MB and the factorisation within a few
    // 1e9 operations (Σ column-height²); beyond that iterate.
    let direct = ndof <= 150000;
    if (direct) { const est = profileOf(); direct = est.prof <= 2.5e7 && est.work <= 4e9; }
    solver = direct ? 'skyline' : 'pcg';
  } else if (solver === 'skyline') {
    profileOf();
  } else if (solver !== 'pcg') {
    throw new RangeError("feSystem: solver must be 'auto', 'skyline' or 'pcg'");
  }

  const gd = new Int32Array(64);
  let sky = null;
  let csr = null;
  if (solver === 'skyline') {
    const heights = new Int32Array(ndof);
    for (let i = 0; i < ndof; i++) heights[i] = i - firstRow[i];
    sky = new SkylineSym(ndof, heights);
    for (let e = 0; e < nElem; e++) {
      const nd = elemNodes(e);
      const ke = elemMatrix(e);
      const m = nd.length * dpn;
      for (let a = 0; a < nd.length; a++) for (let k = 0; k < dpn; k++) gd[a * dpn + k] = nd[a] * dpn + k;
      for (let p = 0; p < m; p++) {
        const op = gd[p];
        if (isFixed[op]) continue;
        const np = order[nd[(p / dpn) | 0]] * dpn + (p % dpn);
        for (let q = 0; q < m; q++) {
          const oq = gd[q];
          if (isFixed[oq]) continue;
          const nq = order[nd[(q / dpn) | 0]] * dpn + (q % dpn);
          // Each unordered dof pair is visited as (p,q) and (q,p); only the upper-triangle visit is stored.
          if (np <= nq) sky.add(np, nq, ke[p * m + q]);
        }
      }
    }
    for (let n = 0; n < nNodes; n++) {
      for (let k = 0; k < dpn; k++) {
        const od = n * dpn + k;
        const nw = order[n] * dpn + k;
        if (isFixed[od]) sky.add(nw, nw, 1);
        else if (springs && springs[od] !== 0) sky.add(nw, nw, springs[od]);
      }
    }
    try {
      sky.factor();
    } catch (err) {
      throw new Error(`Stiffness matrix is singular: the model is insufficiently supported or contains a mechanism (${err.message}).`);
    }
  } else {
    const coo = new CooBuilder(ndof);
    for (let e = 0; e < nElem; e++) {
      const nd = elemNodes(e);
      const ke = elemMatrix(e);
      const m = nd.length * dpn;
      for (let a = 0; a < nd.length; a++) for (let k = 0; k < dpn; k++) gd[a * dpn + k] = nd[a] * dpn + k;
      for (let p = 0; p < m; p++) {
        const op = gd[p];
        if (isFixed[op]) continue;
        for (let q = 0; q < m; q++) {
          const oq = gd[q];
          if (isFixed[oq]) continue;
          const v = ke[p * m + q];
          if (v !== 0) coo.add(op, oq, v);
        }
      }
    }
    for (let i = 0; i < ndof; i++) {
      if (isFixed[i]) coo.add(i, i, 1);
      else if (springs && springs[i] !== 0) coo.add(i, i, springs[i]);
    }
    csr = coo.toCsr();
  }

  const sys = { ndof, solver, isFixed, nFixed, lastIterations: 0 };

  const rawSolve = (rhs) => {
    if (sky) {
      const b = new Float64Array(ndof);
      for (let n = 0; n < nNodes; n++) for (let k = 0; k < dpn; k++) b[order[n] * dpn + k] = rhs[n * dpn + k];
      const x = sky.solve(b);
      const u = new Float64Array(ndof);
      for (let n = 0; n < nNodes; n++) for (let k = 0; k < dpn; k++) u[n * dpn + k] = x[order[n] * dpn + k];
      return u;
    }
    const r = pcg(csr, rhs, { tol, maxIter: maxIter || Math.max(2000, 4 * ndof), precond });
    sys.lastIterations = r.iterations;
    if (!r.converged) throw new Error(`Iterative solver did not converge (${r.iterations} iterations, residual ${r.residual}). Check supports for rigid-body modes.`);
    return Float64Array.from(r.x);
  };

  sys.mulK = (x, y) => {
    const out = y || new Float64Array(ndof);
    out.fill(0);
    for (let e = 0; e < nElem; e++) {
      const nd = elemNodes(e);
      const ke = elemMatrix(e);
      const m = nd.length * dpn;
      for (let a = 0; a < nd.length; a++) for (let k = 0; k < dpn; k++) gd[a * dpn + k] = nd[a] * dpn + k;
      for (let p = 0; p < m; p++) {
        let s = 0;
        const row = p * m;
        for (let q = 0; q < m; q++) s += ke[row + q] * x[gd[q]];
        out[gd[p]] += s;
      }
    }
    if (springs) for (let i = 0; i < ndof; i++) out[i] += springs[i] * x[i];
    return out;
  };

  sys.solve = (f) => {
    const rhs = Float64Array.from(f);
    if (nonzeroBc) {
      for (let e = 0; e < nElem; e++) {
        const nd = elemNodes(e);
        let touched = false;
        for (let a = 0; a < nd.length && !touched; a++) {
          for (let k = 0; k < dpn; k++) { const d = nd[a] * dpn + k; if (isFixed[d] && fixed[d] !== 0) { touched = true; break; } }
        }
        if (!touched) continue;
        const ke = elemMatrix(e);
        const m = nd.length * dpn;
        for (let a = 0; a < nd.length; a++) for (let k = 0; k < dpn; k++) gd[a * dpn + k] = nd[a] * dpn + k;
        for (let q = 0; q < m; q++) {
          const oq = gd[q];
          if (!isFixed[oq] || fixed[oq] === 0) continue;
          const ub = fixed[oq];
          for (let p = 0; p < m; p++) if (!isFixed[gd[p]]) rhs[gd[p]] -= ke[p * m + q] * ub;
        }
      }
    }
    for (let i = 0; i < ndof; i++) if (isFixed[i]) rhs[i] = fixed[i];
    return rawSolve(rhs);
  };

  sys.solveHomogeneous = (b) => {
    const rhs = Float64Array.from(b);
    for (let i = 0; i < ndof; i++) if (isFixed[i]) rhs[i] = 0;
    return rawSolve(rhs);
  };

  return sys;
}

/**
 * Translate a list of support specifications into the `fixed` vector used by `feSystem`.
 * @param {number} nNodes number of nodes
 * @param {string[]} keys dof names in nodal order, e.g. ['ux', 'uy']
 * @param {Array<object>} fix entries { nodes: ArrayLike<number>, [key]: prescribed value } (a key that is present
 *   and numeric constrains that dof; units follow the dof: m or rad or K)
 * @param {string} fn caller name for error messages
 * @returns {Float64Array} NaN for free dofs, prescribed value otherwise
 */
export function buildFixed(nNodes, keys, fix, fn) {
  const dpn = keys.length;
  const fixed = new Float64Array(nNodes * dpn).fill(NaN);
  for (const s of fix || []) {
    const nodes = s.nodes;
    if (!nodes || typeof nodes.length !== 'number') throw new RangeError(`${fn}: every fix entry needs a nodes list`);
    for (let i = 0; i < nodes.length; i++) {
      const n = nodes[i];
      if (!(n >= 0 && n < nNodes)) throw new RangeError(`${fn}: fix node index ${n} out of range`);
      for (let k = 0; k < dpn; k++) {
        const v = s[keys[k]];
        if (typeof v === 'number') fixed[n * dpn + k] = v;
        else if (v === true) fixed[n * dpn + k] = 0;
      }
    }
  }
  return fixed;
}

/* ------------------------------------------------------------------------------------------------------------
 * Mesh generation and topology
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Boundary edges of a 2-D mesh (edges belonging to exactly one element), oriented as in their element
 * (material on the left for counter-clockwise elements).
 * @param {object} mesh 2-D mesh
 * @returns {number[][]} edges [n1, n2] or [n1, n2, nMid]
 */
export function boundaryEdges(mesh) {
  checkMesh(mesh, null, 'boundaryEdges');
  const { nen, ne, conn, nn } = mesh;
  const nc = ELEM[mesh.type].nc;
  const count = new Map();
  for (let e = 0; e < ne; e++) {
    for (let k = 0; k < nc; k++) {
      const a = conn[e * nen + k];
      const b = conn[e * nen + ((k + 1) % nc)];
      const key = Math.min(a, b) * nn + Math.max(a, b);
      count.set(key, (count.get(key) || 0) + 1);
    }
  }
  const out = [];
  for (let e = 0; e < ne; e++) {
    for (let k = 0; k < nc; k++) {
      const a = conn[e * nen + k];
      const b = conn[e * nen + ((k + 1) % nc)];
      const key = Math.min(a, b) * nn + Math.max(a, b);
      if (count.get(key) === 1) out.push(nen > nc ? [a, b, conn[e * nen + nc + k]] : [a, b]);
    }
  }
  return out;
}

/**
 * Nodes satisfying a geometric predicate.
 * @param {object} mesh 2-D mesh
 * @param {(x: number, y: number) => boolean} pred predicate on nodal coordinates [m]
 * @returns {Int32Array} node indices
 */
export function nodesWhere(mesh, pred) {
  const out = [];
  for (let i = 0; i < mesh.nn; i++) if (pred(mesh.nodes[2 * i], mesh.nodes[2 * i + 1])) out.push(i);
  return Int32Array.from(out);
}

/**
 * Boundary edges whose nodes all satisfy a geometric predicate.
 * @param {object} mesh 2-D mesh
 * @param {(x: number, y: number) => boolean} pred predicate on nodal coordinates [m]
 * @returns {number[][]} oriented boundary edges
 */
export function edgesWhere(mesh, pred) {
  return boundaryEdges(mesh).filter((ed) => ed.every((n) => pred(mesh.nodes[2 * n], mesh.nodes[2 * n + 1])));
}

function attachSets(mesh, preds) {
  const all = boundaryEdges(mesh);
  mesh.edges = {};
  mesh.nodeSets = {};
  for (const name of Object.keys(preds)) {
    const p = preds[name];
    mesh.edges[name] = all.filter((ed) => ed.every((n) => p(mesh.nodes[2 * n], mesh.nodes[2 * n + 1])));
    mesh.nodeSets[name] = nodesWhere(mesh, p);
  }
  return mesh;
}

/**
 * Structured mesh of a four-sided block defined by a mapping of the unit square.
 * @param {object} o
 * @param {number} o.ns elements along s
 * @param {number} o.nt elements along t
 * @param {(s: number, t: number) => [number, number]} o.map (s, t) in [0,1]² -> (x, y) [m]; either orientation
 * @param {'T3'|'Q4'|'T6'|'Q8'} [o.type] element type (default 'Q4'); mid-side nodes are placed on the mapped geometry
 * @returns {object} mesh (see file header) with nodeSets s0, s1, t0, t1 (sides s = 0, s = 1, t = 0, t = 1)
 *   and matching oriented boundary edge lists in `edges`
 */
export function mappedMesh({ ns, nt, map, type = 'Q4' }) {
  if (!(Number.isInteger(ns) && ns >= 1)) throw new RangeError('mappedMesh: ns must be a positive integer');
  if (!(Number.isInteger(nt) && nt >= 1)) throw new RangeError('mappedMesh: nt must be a positive integer');
  const el = ELEM[type];
  if (!el) throw new RangeError(`mappedMesh: unknown element type ${type}`);
  const quad = el.nen > el.nc;
  const m = quad ? 2 : 1;
  const gs = m * ns + 1;
  const gt = m * nt + 1;
  const id = new Int32Array(gs * gt).fill(-1);
  const coords = [];
  let nn = 0;
  const node = (i, j) => {
    const k = j * gs + i;
    if (id[k] < 0) {
      const p = map(i / (gs - 1), j / (gt - 1));
      coords.push(p[0], p[1]);
      id[k] = nn++;
    }
    return id[k];
  };
  const c0 = map(0.5, 0.5);
  const cs = map(0.5 + 1e-4, 0.5);
  const ct = map(0.5, 0.5 + 1e-4);
  const flip = (cs[0] - c0[0]) * (ct[1] - c0[1]) - (cs[1] - c0[1]) * (ct[0] - c0[0]) < 0;
  const conn = [];
  const push = (arr, perm) => { for (let k = 0; k < arr.length; k++) conn.push(arr[flip ? perm[k] : k]); };
  for (let j = 0; j < nt; j++) {
    for (let i = 0; i < ns; i++) {
      const a = m * i;
      const b = m * j;
      if (type === 'Q4') {
        push([node(a, b), node(a + 1, b), node(a + 1, b + 1), node(a, b + 1)], [0, 3, 2, 1]);
      } else if (type === 'Q8') {
        push([node(a, b), node(a + 2, b), node(a + 2, b + 2), node(a, b + 2),
          node(a + 1, b), node(a + 2, b + 1), node(a + 1, b + 2), node(a, b + 1)], [0, 3, 2, 1, 7, 6, 5, 4]);
      } else if (type === 'T3') {
        push([node(a, b), node(a + 1, b), node(a + 1, b + 1)], [0, 2, 1]);
        push([node(a, b), node(a + 1, b + 1), node(a, b + 1)], [0, 2, 1]);
      } else {
        push([node(a, b), node(a + 2, b), node(a + 2, b + 2), node(a + 1, b), node(a + 2, b + 1), node(a + 1, b + 1)],
          [0, 2, 1, 5, 4, 3]);
        push([node(a, b), node(a + 2, b + 2), node(a, b + 2), node(a + 1, b + 1), node(a + 1, b + 2), node(a, b + 1)],
          [0, 2, 1, 5, 4, 3]);
      }
    }
  }
  const mesh = { type, nen: el.nen, nn, ne: conn.length / el.nen, nodes: Float64Array.from(coords), conn: Int32Array.from(conn) };
  const side = { s0: [], s1: [], t0: [], t1: [] };
  for (let j = 0; j < gt; j++) {
    if (id[j * gs] >= 0) side.s0.push(id[j * gs]);
    if (id[j * gs + gs - 1] >= 0) side.s1.push(id[j * gs + gs - 1]);
  }
  for (let i = 0; i < gs; i++) {
    if (id[i] >= 0) side.t0.push(id[i]);
    if (id[(gt - 1) * gs + i] >= 0) side.t1.push(id[(gt - 1) * gs + i]);
  }
  const all = boundaryEdges(mesh);
  mesh.nodeSets = {};
  mesh.edges = {};
  for (const name of Object.keys(side)) {
    const set = new Set(side[name]);
    mesh.nodeSets[name] = Int32Array.from(side[name]);
    mesh.edges[name] = all.filter((ed) => ed.every((n) => set.has(n)));
  }
  return mesh;
}

function renameSets(mesh, names) {
  const edges = {};
  const nodeSets = {};
  for (const k of Object.keys(names)) { edges[names[k]] = mesh.edges[k]; nodeSets[names[k]] = mesh.nodeSets[k]; }
  mesh.edges = edges;
  mesh.nodeSets = nodeSets;
  return mesh;
}

/**
 * Structured rectangle mesh.
 * @param {object} o
 * @param {number} o.lx length in x [m]
 * @param {number} o.ly length in y [m]
 * @param {number} o.nx elements along x
 * @param {number} o.ny elements along y
 * @param {'T3'|'Q4'|'T6'|'Q8'} [o.type] element type (default 'Q4')
 * @param {number} [o.x0] origin x [m]
 * @param {number} [o.y0] origin y [m]
 * @param {(x: number, y: number) => [number, number]} [o.distort] optional smooth node relocation (for patch tests)
 * @returns {object} mesh with sets left, right, bottom, top
 */
export function rectangleMesh({ lx, ly, nx, ny, type = 'Q4', x0 = 0, y0 = 0, distort = null }) {
  positive(lx, 'lx', 'rectangleMesh');
  positive(ly, 'ly', 'rectangleMesh');
  const mesh = mappedMesh({
    ns: nx, nt: ny, type,
    map: (s, t) => {
      const x = x0 + s * lx;
      const y = y0 + t * ly;
      return distort ? distort(x, y) : [x, y];
    },
  });
  return renameSets(mesh, { s0: 'left', s1: 'right', t0: 'bottom', t1: 'top' });
}

/**
 * Structured mesh of an annular sector (plane models of rings, thick cylinders in plane strain).
 * @param {object} o
 * @param {number} o.ri inner radius [m]
 * @param {number} o.ro outer radius [m]
 * @param {number} o.nr elements through the wall
 * @param {number} o.ntheta elements around the sector
 * @param {number} [o.theta0] start angle [rad] (default 0)
 * @param {number} [o.theta1] end angle [rad] (default π/2, quarter model)
 * @param {'T3'|'Q4'|'T6'|'Q8'} [o.type] element type
 * @returns {object} mesh with sets inner, outer, start (θ = theta0), end (θ = theta1)
 */
export function annulusMesh({ ri, ro, nr, ntheta, theta0 = 0, theta1 = Math.PI / 2, type = 'Q4' }) {
  positive(ri, 'ri', 'annulusMesh');
  if (!(ro > ri)) throw new RangeError('annulusMesh: ro must exceed ri');
  if (!(theta1 > theta0) || theta1 - theta0 >= 2 * Math.PI) throw new RangeError('annulusMesh: need theta0 < theta1 < theta0 + 2π (sector)');
  const mesh = mappedMesh({
    ns: nr, nt: ntheta, type,
    map: (s, t) => {
      const r = ri + s * (ro - ri);
      const th = theta0 + t * (theta1 - theta0);
      return [r * Math.cos(th), r * Math.sin(th)];
    },
  });
  return renameSets(mesh, { s0: 'inner', s1: 'outer', t0: 'start', t1: 'end' });
}

/**
 * Axisymmetric (r, z) mesh of a tube wall.
 * @param {object} o
 * @param {number} o.ri inner radius [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.length axial length [m]
 * @param {number} o.nr elements through the wall
 * @param {number} o.nz elements along the axis
 * @param {'T3'|'Q4'|'T6'|'Q8'} [o.type] element type
 * @returns {object} mesh (x = r, y = z) with sets inner, outer, bottom (z = 0), top (z = length)
 */
export function tubeWallMesh({ ri, t, length, nr, nz, type = 'Q4' }) {
  if (!(ri >= 0)) throw new RangeError('tubeWallMesh: ri must be >= 0');
  positive(t, 't', 'tubeWallMesh');
  positive(length, 'length', 'tubeWallMesh');
  const mesh = rectangleMesh({ lx: t, ly: length, nx: nr, ny: nz, type, x0: ri });
  return renameSets(mesh, { left: 'inner', right: 'outer', bottom: 'bottom', top: 'top' });
}

/**
 * Merge meshes of the same element type, fusing coincident nodes.
 * @param {object[]} meshes meshes of identical type
 * @param {number} [tol] coincidence tolerance [m] (default 1e-9 × bounding size)
 * @returns {object} merged mesh (no edge/node sets)
 */
export function mergeMeshes(meshes, tol) {
  if (!Array.isArray(meshes) || meshes.length === 0) throw new RangeError('mergeMeshes: meshes must be a non-empty array');
  const type = meshes[0].type;
  let size = 0;
  for (const m of meshes) {
    if (m.type !== type) throw new RangeError('mergeMeshes: all meshes must share one element type');
    for (let i = 0; i < m.nodes.length; i++) size = Math.max(size, Math.abs(m.nodes[i]));
  }
  const h = tol || 1e-9 * (size || 1);
  const key = (x, y) => `${Math.round(x / h)}_${Math.round(y / h)}`;
  const seen = new Map();
  const coords = [];
  const conn = [];
  for (const m of meshes) {
    const local = new Int32Array(m.nn);
    for (let i = 0; i < m.nn; i++) {
      const x = m.nodes[2 * i];
      const y = m.nodes[2 * i + 1];
      const k = key(x, y);
      let g = seen.get(k);
      if (g === undefined) { g = coords.length / 2; coords.push(x, y); seen.set(k, g); }
      local[i] = g;
    }
    for (let i = 0; i < m.conn.length; i++) conn.push(local[m.conn[i]]);
  }
  const nen = ELEM[type].nen;
  return { type, nen, nn: coords.length / 2, ne: conn.length / nen, nodes: Float64Array.from(coords), conn: Int32Array.from(conn) };
}

/**
 * Quarter model of a rectangular plate with a central circular hole (hole centre at the origin, quarter in
 * x >= 0, y >= 0). Two mapped blocks meet on the line from the hole at 45° to the outer corner.
 * @param {object} o
 * @param {number} o.a hole radius [m]
 * @param {number} o.w half-width of the plate (x extent of the quarter) [m]
 * @param {number} o.h half-height of the plate (y extent of the quarter) [m]
 * @param {number} o.nr elements from hole to outer edge
 * @param {number} o.ntheta elements along the hole per block (two blocks)
 * @param {number} [o.grading] ratio of outermost to innermost radial element size (default 8, refines at the hole)
 * @param {'T3'|'Q4'|'T6'|'Q8'} [o.type] element type
 * @returns {object} mesh with sets hole, symX (edge on y = 0), symY (edge on x = 0), right (x = w), top (y = h)
 */
export function plateWithHoleMesh({ a, w, h, nr, ntheta, grading = 8, type = 'Q4' }) {
  positive(a, 'a', 'plateWithHoleMesh');
  if (!(w > a)) throw new RangeError('plateWithHoleMesh: w must exceed a');
  if (!(h > a)) throw new RangeError('plateWithHoleMesh: h must exceed a');
  positive(grading, 'grading', 'plateWithHoleMesh');
  // Geometric grading: element sizes grow by a constant ratio so that last/first = grading.
  const q = nr > 1 ? Math.pow(grading, 1 / (nr - 1)) : 1;
  const grade = (s) => {
    if (Math.abs(q - 1) < 1e-12) return s;
    return (Math.pow(q, s * nr) - 1) / (Math.pow(q, nr) - 1);
  };
  const block = (thA, thB, outer) => mappedMesh({
    ns: nr, nt: ntheta, type,
    map: (s, t) => {
      const th = thA + t * (thB - thA);
      const pi = [a * Math.cos(th), a * Math.sin(th)];
      const po = outer(t);
      const g = grade(s);
      return [pi[0] + g * (po[0] - pi[0]), pi[1] + g * (po[1] - pi[1])];
    },
  });
  const b1 = block(0, Math.PI / 4, (t) => [w, t * h]);
  const b2 = block(Math.PI / 4, Math.PI / 2, (t) => [w * (1 - t), h]);
  const mesh = mergeMeshes([b1, b2], 1e-9 * Math.max(w, h));
  const eps = 1e-7 * Math.max(w, h);
  return attachSets(mesh, {
    hole: (x, y) => Math.abs(Math.hypot(x, y) - a) < 1e-6 * a + eps,
    symX: (x, y) => Math.abs(y) < eps,
    symY: (x) => Math.abs(x) < eps,
    right: (x) => Math.abs(x - w) < eps,
    top: (x, y) => Math.abs(y - h) < eps,
  });
}

/**
 * Uniform refinement: every T3 is split into four T3 and every Q4 into four Q4 (edge mid-points are placed on
 * straight edges; regenerate curved geometry from its generator instead when boundary curvature matters).
 * @param {object} mesh T3 or Q4 mesh
 * @returns {object} refined mesh (element size halved; edge and node sets are not carried over)
 */
export function refineUniform(mesh) {
  checkMesh(mesh, ['T3', 'Q4'], 'refineUniform');
  const { nn, ne, nen, conn, nodes } = mesh;
  const coords = Array.from(nodes);
  const mids = new Map();
  let count = nn;
  const mid = (a, b) => {
    const key = Math.min(a, b) * nn + Math.max(a, b);
    let m = mids.get(key);
    if (m === undefined) {
      m = count++;
      coords.push(0.5 * (nodes[2 * a] + nodes[2 * b]), 0.5 * (nodes[2 * a + 1] + nodes[2 * b + 1]));
      mids.set(key, m);
    }
    return m;
  };
  const out = [];
  for (let e = 0; e < ne; e++) {
    const c = conn.subarray(e * nen, e * nen + nen);
    if (nen === 3) {
      const m01 = mid(c[0], c[1]);
      const m12 = mid(c[1], c[2]);
      const m20 = mid(c[2], c[0]);
      out.push(c[0], m01, m20, m01, c[1], m12, m20, m12, c[2], m01, m12, m20);
    } else {
      const m01 = mid(c[0], c[1]);
      const m12 = mid(c[1], c[2]);
      const m23 = mid(c[2], c[3]);
      const m30 = mid(c[3], c[0]);
      const cc = count++;
      coords.push(0.25 * (nodes[2 * c[0]] + nodes[2 * c[1]] + nodes[2 * c[2]] + nodes[2 * c[3]]),
        0.25 * (nodes[2 * c[0] + 1] + nodes[2 * c[1] + 1] + nodes[2 * c[2] + 1] + nodes[2 * c[3] + 1]));
      out.push(c[0], m01, cc, m30, m01, c[1], m12, cc, cc, m12, c[2], m23, m30, cc, m23, c[3]);
    }
  }
  return { type: mesh.type, nen, nn: count, ne: out.length / nen, nodes: Float64Array.from(coords), conn: Int32Array.from(out) };
}

/**
 * Element areas and total area of a 2-D mesh (numerically integrated Jacobian, exact for the element geometry).
 * @param {object} mesh 2-D mesh
 * @returns {{area: number, elemArea: Float64Array}} total and per-element area [m²]
 */
export function meshArea(mesh) {
  checkMesh(mesh, null, 'meshArea');
  const el = ELEM[mesh.type];
  const { nen, ne } = mesh;
  const N = new Float64Array(nen);
  const dxi = new Float64Array(nen);
  const deta = new Float64Array(nen);
  const elemArea = new Float64Array(ne);
  let area = 0;
  const rule = el.mass;
  for (let e = 0; e < ne; e++) {
    let a = 0;
    for (let g = 0; g < rule.w.length; g++) {
      shape2d(mesh.type, rule.xi[g], rule.eta[g], N, dxi, deta);
      let j11 = 0; let j12 = 0; let j21 = 0; let j22 = 0;
      for (let i = 0; i < nen; i++) {
        const n = mesh.conn[e * nen + i];
        j11 += dxi[i] * mesh.nodes[2 * n]; j12 += dxi[i] * mesh.nodes[2 * n + 1];
        j21 += deta[i] * mesh.nodes[2 * n]; j22 += deta[i] * mesh.nodes[2 * n + 1];
      }
      a += rule.w[g] * (j11 * j22 - j12 * j21);
    }
    elemArea[e] = a;
    area += a;
  }
  return { area, elemArea };
}

/* ------------------------------------------------------------------------------------------------------------
 * Plane / axisymmetric elasticity
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Isotropic constitutive matrix for the strain vector [εxx, εyy, εzz (or εθθ), γxy], row-major 4×4.
 * Plane stress has σzz = 0 (third row/column zero); plane strain keeps the third row so σzz = ν(σxx + σyy).
 * @param {number} E Young's modulus [Pa]
 * @param {number} nu Poisson's ratio [-]
 * @param {'stress'|'strain'|'axisym'} mode kinematic assumption
 * @returns {Float64Array} D [Pa]
 */
export function planeD(E, nu, mode) {
  const D = new Float64Array(16);
  const G = E / (2 * (1 + nu));
  if (mode === 'stress') {
    const c = E / (1 - nu * nu);
    D[0] = c; D[1] = c * nu; D[4] = c * nu; D[5] = c;
  } else if (mode === 'strain' || mode === 'axisym') {
    const lam = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) D[i * 4 + j] = lam + (i === j ? 2 * G : 0);
  } else {
    throw new RangeError("planeD: mode must be 'stress', 'strain' or 'axisym'");
  }
  D[15] = G;
  return D;
}

function makeScratch(nen) {
  return { N: new Float64Array(nen), dxi: new Float64Array(nen), deta: new Float64Array(nen), dx: new Float64Array(nen), dy: new Float64Array(nen), det: 0, r: 0, x: 0, y: 0, ji: new Float64Array(4) };
}

/** Shape functions, Cartesian derivatives and Jacobian determinant at a natural point. */
function geom(type, nen, xi, eta, xe, sc) {
  shape2d(type, xi, eta, sc.N, sc.dxi, sc.deta);
  let j11 = 0; let j12 = 0; let j21 = 0; let j22 = 0; let x = 0; let y = 0;
  for (let i = 0; i < nen; i++) {
    const X = xe[2 * i];
    const Y = xe[2 * i + 1];
    j11 += sc.dxi[i] * X; j12 += sc.dxi[i] * Y; j21 += sc.deta[i] * X; j22 += sc.deta[i] * Y;
    x += sc.N[i] * X; y += sc.N[i] * Y;
  }
  const det = j11 * j22 - j12 * j21;
  if (!(det > 0)) throw new Error('Element with non-positive Jacobian (inverted or badly distorted element).');
  sc.ji[0] = j22 / det; sc.ji[1] = -j12 / det; sc.ji[2] = -j21 / det; sc.ji[3] = j11 / det;
  for (let i = 0; i < nen; i++) {
    sc.dx[i] = sc.ji[0] * sc.dxi[i] + sc.ji[1] * sc.deta[i];
    sc.dy[i] = sc.ji[2] * sc.dxi[i] + sc.ji[3] * sc.deta[i];
  }
  sc.det = det; sc.r = x; sc.x = x; sc.y = y;
}

function fillB(nen, sc, axisym, B) {
  const nd = 2 * nen;
  for (let i = 0; i < nen; i++) {
    B[2 * i] = sc.dx[i];
    B[nd + 2 * i + 1] = sc.dy[i];
    B[2 * nd + 2 * i] = axisym ? sc.N[i] / sc.r : 0;
    B[3 * nd + 2 * i] = sc.dy[i];
    B[3 * nd + 2 * i + 1] = sc.dx[i];
  }
}

/** ke += f · Bᵀ D B for a 4-row B with nd columns. */
function addBtDB(B, D, nd, f, ke, DB) {
  for (let r = 0; r < 4; r++) {
    for (let q = 0; q < nd; q++) {
      DB[r * nd + q] = D[r * 4] * B[q] + D[r * 4 + 1] * B[nd + q] + D[r * 4 + 2] * B[2 * nd + q] + D[r * 4 + 3] * B[3 * nd + q];
    }
  }
  for (let p = 0; p < nd; p++) {
    const b0 = B[p]; const b1 = B[nd + p]; const b2 = B[2 * nd + p]; const b3 = B[3 * nd + p];
    for (let q = 0; q < nd; q++) {
      ke[p * nd + q] += f * (b0 * DB[q] + b1 * DB[nd + q] + b2 * DB[2 * nd + q] + b3 * DB[3 * nd + q]);
    }
  }
}

/**
 * Q4 membrane with two incompatible bending modes per displacement component (QM6). The enhanced gradients use
 * the centre Jacobian scaled by det J0 / det J, so the element passes the patch test for any shape.
 * Returns the statically condensed 8×8 stiffness and the 4×8 recovery operator α = rec · u.
 */
function q4IncompatibleKe(xe, D, thick, work) {
  const { sc, B, DB, kuu, kua, kaa, G, rec, ke } = work;
  kuu.fill(0); kua.fill(0); kaa.fill(0);
  geom('Q4', 4, 0, 0, xe, sc);
  const j0 = sc.det;
  const i00 = sc.ji[0]; const i01 = sc.ji[1]; const i10 = sc.ji[2]; const i11 = sc.ji[3];
  const rule = ELEM.Q4.rule;
  const DG = work.DG;
  for (let g = 0; g < 4; g++) {
    const xi = rule.xi[g];
    const eta = rule.eta[g];
    geom('Q4', 4, xi, eta, xe, sc);
    fillB(4, sc, false, B);
    const f = rule.w[g] * sc.det * thick;
    addBtDB(B, D, 8, f, kuu, DB);
    const s = j0 / sc.det;
    // Cartesian gradients of the modes (1 - ξ²) and (1 - η²).
    const m1x = s * i00 * (-2 * xi); const m1y = s * i10 * (-2 * xi);
    const m2x = s * i01 * (-2 * eta); const m2y = s * i11 * (-2 * eta);
    G.fill(0);
    G[0] = m1x; G[12] = m1y; G[4 + 1] = m1y; G[12 + 1] = m1x;
    G[2] = m2x; G[12 + 2] = m2y; G[4 + 3] = m2y; G[12 + 3] = m2x;
    for (let r = 0; r < 4; r++) for (let q = 0; q < 4; q++) {
      DG[r * 4 + q] = D[r * 4] * G[q] + D[r * 4 + 1] * G[4 + q] + D[r * 4 + 2] * G[8 + q] + D[r * 4 + 3] * G[12 + q];
    }
    for (let p = 0; p < 8; p++) for (let q = 0; q < 4; q++) {
      kua[p * 4 + q] += f * (B[p] * DG[q] + B[8 + p] * DG[4 + q] + B[16 + p] * DG[8 + q] + B[24 + p] * DG[12 + q]);
    }
    for (let p = 0; p < 4; p++) for (let q = 0; q < 4; q++) {
      kaa[p * 4 + q] += f * (G[p] * DG[q] + G[4 + p] * DG[4 + q] + G[8 + p] * DG[8 + q] + G[12 + p] * DG[12 + q]);
    }
  }
  invertSmall(kaa, 4);
  for (let a = 0; a < 4; a++) for (let q = 0; q < 8; q++) {
    let s = 0;
    for (let b = 0; b < 4; b++) s += kaa[a * 4 + b] * kua[q * 4 + b];
    rec[a * 8 + q] = -s;
  }
  for (let p = 0; p < 8; p++) for (let q = 0; q < 8; q++) {
    let s = kuu[p * 8 + q];
    for (let a = 0; a < 4; a++) s += kua[p * 4 + a] * rec[a * 8 + q];
    ke[p * 8 + q] = s;
  }
  return ke;
}

function makeQ4Work() {
  return { sc: makeScratch(4), B: new Float64Array(32), DB: new Float64Array(32), kuu: new Float64Array(64), kua: new Float64Array(32), kaa: new Float64Array(16), G: new Float64Array(16), DG: new Float64Array(16), rec: new Float64Array(32), ke: new Float64Array(64) };
}

/** Incompatible-mode strain contribution at a natural point (adds G·α to the 4-component strain). */
function q4EnhancedStrain(xe, xi, eta, alpha, sc, eps) {
  geom('Q4', 4, 0, 0, xe, sc);
  const j0 = sc.det;
  const i00 = sc.ji[0]; const i01 = sc.ji[1]; const i10 = sc.ji[2]; const i11 = sc.ji[3];
  geom('Q4', 4, xi, eta, xe, sc);
  const s = j0 / sc.det;
  const m1x = s * i00 * (-2 * xi); const m1y = s * i10 * (-2 * xi);
  const m2x = s * i01 * (-2 * eta); const m2y = s * i11 * (-2 * eta);
  eps[0] += m1x * alpha[0] + m2x * alpha[2];
  eps[1] += m1y * alpha[1] + m2y * alpha[3];
  eps[3] += m1y * alpha[0] + m1x * alpha[1] + m2y * alpha[2] + m2x * alpha[3];
}

function evalLoad(v, x, y) {
  return typeof v === 'function' ? v(x, y) : (v || 0);
}

/** Consistent nodal loads from tractions/pressure on boundary edges; adds into f (dofs ux, uy). */
function addEdgeLoads(mesh, edgeLoads, widthAt, f, fn) {
  const g2 = gauss1d(2);
  const g3 = gauss1d(3);
  for (const L of edgeLoads || []) {
    if (!Array.isArray(L.edges)) throw new RangeError(`${fn}: every edgeLoads entry needs an edges list`);
    for (const ed of L.edges) {
      const quad = ed.length === 3;
      const g = quad ? g3 : g2;
      for (let k = 0; k < g.x.length; k++) {
        const s = g.x[k];
        const N = quad ? [0.5 * s * (s - 1), 0.5 * s * (s + 1), 1 - s * s] : [0.5 * (1 - s), 0.5 * (1 + s)];
        const dN = quad ? [s - 0.5, s + 0.5, -2 * s] : [-0.5, 0.5];
        let x = 0; let y = 0; let tx = 0; let ty = 0;
        for (let i = 0; i < ed.length; i++) {
          x += N[i] * mesh.nodes[2 * ed[i]]; y += N[i] * mesh.nodes[2 * ed[i] + 1];
          tx += dN[i] * mesh.nodes[2 * ed[i]]; ty += dN[i] * mesh.nodes[2 * ed[i] + 1];
        }
        // Outward normal of a counter-clockwise boundary is the tangent rotated clockwise; |t| dξ = ds.
        const p = evalLoad(L.pressure, x, y);
        const qx = evalLoad(L.tx, x, y) * Math.hypot(tx, ty) - p * ty;
        const qy = evalLoad(L.ty, x, y) * Math.hypot(tx, ty) + p * tx;
        const wgt = g.w[k] * widthAt(x);
        for (let i = 0; i < ed.length; i++) {
          f[2 * ed[i]] += N[i] * qx * wgt;
          f[2 * ed[i] + 1] += N[i] * qy * wgt;
        }
      }
    }
  }
}

/** Recovery (sampling) points and the extrapolation matrix to the element nodes. */
function recoveryScheme(type) {
  const el = ELEM[type];
  if (type === 'T3') return { xi: [1 / 3], eta: [1 / 3], w: [0.5], ex: [[1], [1], [1]] };
  if (el.tri) {
    const L = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]];
    return { xi: [1 / 6, 2 / 3, 1 / 6], eta: [1 / 6, 1 / 6, 2 / 3], w: [1 / 6, 1 / 6, 1 / 6], ex: L.map((l) => l.map((v) => 2 * v - 1 / 3)) };
  }
  const r = quadRule(2);
  const ex = [];
  for (let i = 0; i < el.nen; i++) {
    const row = [];
    for (let g = 0; g < 4; g++) {
      row.push(0.25 * (1 + 3 * QUAD_XI[i] * r.xi[g]) * (1 + 3 * QUAD_ETA[i] * r.eta[g]));
    }
    ex.push(row);
  }
  return { xi: r.xi, eta: r.eta, w: r.w, ex };
}

/**
 * Linear elastic plane-stress, plane-strain or axisymmetric analysis.
 *
 * Elements: T3 (constant strain), Q4 (default: incompatible modes QM6, free of shear locking in bending;
 * `q4: 'full'` gives the plain 2×2 isoparametric element; axisymmetric Q4 always uses 2×2), T6 and Q8 (full
 * integration). Axisymmetric models use x = r, y = z and loads are totals over the full 2π circumference.
 *
 * @param {object} o
 * @param {object} o.mesh 2-D mesh (see file header)
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} [o.thickness] out-of-plane thickness [m] (plane modes; default 1)
 * @param {'stress'|'strain'|'axisym'} [o.mode] default 'stress'
 * @param {Array<{nodes: ArrayLike<number>, ux?: number, uy?: number}>} o.fix prescribed displacements [m]
 * @param {Array<{node?: number, nodes?: ArrayLike<number>, fx?: number, fy?: number}>} [o.nodalLoads] point loads
 *   [N] (axisymmetric: total ring load [N])
 * @param {Array<{edges: number[][], tx?: number|Function, ty?: number|Function, pressure?: number|Function}>} [o.edgeLoads]
 *   tractions [Pa] in global axes and/or pressure [Pa] acting against the outward normal; functions of (x, y) allowed
 * @param {[number, number]} [o.bodyForce] body force density [N/m³]
 * @param {'incompatible'|'full'} [o.q4] Q4 formulation
 * @param {'auto'|'skyline'|'pcg'} [o.solver] linear solver
 * @returns {{u: Float64Array, reactions: Float64Array, reactionSum: {fx: number, fy: number},
 *   appliedSum: {fx: number, fy: number}, stress: {sxx: Float64Array, syy: Float64Array, szz: Float64Array,
 *   sxy: Float64Array, vonMises: Float64Array, s1: Float64Array, s2: Float64Array}, gauss: {n: number,
 *   stress: Float64Array}, strainEnergy: number, maxVonMises: number, maxDisp: number, solver: string}}
 *   u [m] (ux, uy per node); reactions [N] per dof (non-zero at supports); nodal-averaged stresses [Pa]
 *   (szz is the out-of-plane or hoop stress); gauss.stress holds 4 components per recovery point per element
 *   [Pa]; strainEnergy [J]
 */
export function solvePlane(o) {
  const { mesh, E, nu, thickness = 1, mode = 'stress', fix, nodalLoads = [], edgeLoads = [], bodyForce = null, q4 = 'incompatible', solver = 'auto' } = o;
  const fn = 'solvePlane';
  checkMesh(mesh, null, fn);
  positive(E, 'E', fn);
  checkNu(nu, fn);
  positive(thickness, 'thickness', fn);
  if (q4 !== 'incompatible' && q4 !== 'full') throw new RangeError(`${fn}: q4 must be 'incompatible' or 'full'`);
  const axisym = mode === 'axisym';
  const D = planeD(E, nu, mode);
  const { type, nen, nn, ne, conn, nodes } = mesh;
  const el = ELEM[type];
  const nd = 2 * nen;
  const useIncomp = type === 'Q4' && q4 === 'incompatible' && !axisym;
  const sc = makeScratch(nen);
  const B = new Float64Array(4 * nd);
  const DB = new Float64Array(4 * nd);
  const ke = new Float64Array(nd * nd);
  const xe = new Float64Array(nd);
  const work = useIncomp ? makeQ4Work() : null;
  const TWO_PI = 2 * Math.PI;
  const loadXe = (e) => {
    for (let i = 0; i < nen; i++) { const n = conn[e * nen + i]; xe[2 * i] = nodes[2 * n]; xe[2 * i + 1] = nodes[2 * n + 1]; }
  };
  const elemNodes = (e) => conn.subarray(e * nen, e * nen + nen);
  const elemMatrix = (e) => {
    loadXe(e);
    if (useIncomp) return q4IncompatibleKe(xe, D, thickness, work);
    ke.fill(0);
    const rule = el.rule;
    for (let g = 0; g < rule.w.length; g++) {
      geom(type, nen, rule.xi[g], rule.eta[g], xe, sc);
      if (axisym && !(sc.r > 0)) throw new Error('Axisymmetric model has an integration point at r <= 0.');
      fillB(nen, sc, axisym, B);
      addBtDB(B, D, nd, rule.w[g] * sc.det * (axisym ? TWO_PI * sc.r : thickness), ke, DB);
    }
    return ke;
  };

  const f = new Float64Array(2 * nn);
  for (const L of nodalLoads) {
    const list = L.nodes || [L.node];
    for (let i = 0; i < list.length; i++) {
      const n = list[i];
      if (!(n >= 0 && n < nn)) throw new RangeError(`${fn}: nodalLoads node index out of range`);
      f[2 * n] += L.fx || 0;
      f[2 * n + 1] += L.fy || 0;
    }
  }
  addEdgeLoads(mesh, edgeLoads, (x) => (axisym ? TWO_PI * x : thickness), f, fn);
  if (bodyForce) {
    const rule = el.mass;
    for (let e = 0; e < ne; e++) {
      loadXe(e);
      for (let g = 0; g < rule.w.length; g++) {
        geom(type, nen, rule.xi[g], rule.eta[g], xe, sc);
        const w = rule.w[g] * sc.det * (axisym ? TWO_PI * sc.r : thickness);
        for (let i = 0; i < nen; i++) {
          const n = conn[e * nen + i];
          f[2 * n] += sc.N[i] * bodyForce[0] * w;
          f[2 * n + 1] += sc.N[i] * bodyForce[1] * w;
        }
      }
    }
  }
  const fixed = buildFixed(nn, ['ux', 'uy'], fix, fn);
  if (fixed.every(Number.isNaN)) throw new Error('Model has no supports: add displacement constraints.');
  const sys = feSystem({ nNodes: nn, dpn: 2, nElem: ne, elemNodes, elemMatrix, fixed, solver });
  const u = sys.solve(f);
  const ku = sys.mulK(u);
  const reactions = new Float64Array(2 * nn);
  let energy = 0;
  let rx = 0; let ry = 0; let ax = 0; let ay = 0;
  for (let i = 0; i < 2 * nn; i++) {
    energy += 0.5 * u[i] * ku[i];
    if (sys.isFixed[i]) { reactions[i] = ku[i] - f[i]; if (i % 2 === 0) rx += reactions[i]; else ry += reactions[i]; }
    if (i % 2 === 0) ax += f[i]; else ay += f[i];
  }

  // Stress recovery: sample at the recovery points, extrapolate to the element nodes, average at the nodes.
  const rs = recoveryScheme(type);
  const ng = rs.w.length;
  const gstress = new Float64Array(ne * ng * 4);
  const acc = new Float64Array(nn * 4);
  const cnt = new Uint16Array(nn);
  const eps = new Float64Array(4);
  const ue = new Float64Array(nd);
  const alpha = new Float64Array(4);
  for (let e = 0; e < ne; e++) {
    loadXe(e);
    for (let i = 0; i < nen; i++) { const n = conn[e * nen + i]; ue[2 * i] = u[2 * n]; ue[2 * i + 1] = u[2 * n + 1]; }
    if (useIncomp) {
      q4IncompatibleKe(xe, D, thickness, work);
      for (let a = 0; a < 4; a++) { let s = 0; for (let q = 0; q < 8; q++) s += work.rec[a * 8 + q] * ue[q]; alpha[a] = s; }
    }
    for (let g = 0; g < ng; g++) {
      geom(type, nen, rs.xi[g], rs.eta[g], xe, sc);
      eps.fill(0);
      for (let i = 0; i < nen; i++) {
        eps[0] += sc.dx[i] * ue[2 * i];
        eps[1] += sc.dy[i] * ue[2 * i + 1];
        if (axisym) eps[2] += (sc.N[i] / sc.r) * ue[2 * i];
        eps[3] += sc.dy[i] * ue[2 * i] + sc.dx[i] * ue[2 * i + 1];
      }
      if (useIncomp) q4EnhancedStrain(xe, rs.xi[g], rs.eta[g], alpha, sc, eps);
      const base = (e * ng + g) * 4;
      for (let r = 0; r < 4; r++) gstress[base + r] = D[r * 4] * eps[0] + D[r * 4 + 1] * eps[1] + D[r * 4 + 2] * eps[2] + D[r * 4 + 3] * eps[3];
    }
    for (let i = 0; i < nen; i++) {
      const n = conn[e * nen + i];
      cnt[n]++;
      for (let g = 0; g < ng; g++) {
        const w = rs.ex[i][g];
        const base = (e * ng + g) * 4;
        for (let r = 0; r < 4; r++) acc[n * 4 + r] += w * gstress[base + r];
      }
    }
  }
  const stress = { sxx: new Float64Array(nn), syy: new Float64Array(nn), szz: new Float64Array(nn), sxy: new Float64Array(nn), vonMises: new Float64Array(nn), s1: new Float64Array(nn), s2: new Float64Array(nn) };
  let maxVm = 0;
  for (let n = 0; n < nn; n++) {
    const c = cnt[n] || 1;
    const sx = acc[n * 4] / c; const sy = acc[n * 4 + 1] / c; const sz = acc[n * 4 + 2] / c; const txy = acc[n * 4 + 3] / c;
    stress.sxx[n] = sx; stress.syy[n] = sy; stress.szz[n] = sz; stress.sxy[n] = txy;
    const vm = Math.sqrt(Math.max(0, sx * sx + sy * sy + sz * sz - sx * sy - sy * sz - sz * sx + 3 * txy * txy));
    stress.vonMises[n] = vm;
    if (vm > maxVm) maxVm = vm;
    const c0 = 0.5 * (sx + sy);
    const rr = Math.hypot(0.5 * (sx - sy), txy);
    stress.s1[n] = c0 + rr; stress.s2[n] = c0 - rr;
  }
  let maxDisp = 0;
  for (let n = 0; n < nn; n++) maxDisp = Math.max(maxDisp, Math.hypot(u[2 * n], u[2 * n + 1]));
  return {
    u, reactions, reactionSum: { fx: rx, fy: ry }, appliedSum: { fx: ax, fy: ay }, stress,
    gauss: { n: ng, stress: gstress }, strainEnergy: energy, maxVonMises: maxVm, maxDisp, solver: sys.solver,
  };
}

/**
 * Zienkiewicz–Zhu a-posteriori error indicator in the energy norm: the recovered (nodal-averaged) stress field
 * σ* is compared with the element stresses σh, ‖e‖² = ∫ (σ* − σh)ᵀ C (σ* − σh) dΩ with C the elastic compliance.
 * @param {object} o
 * @param {object} o.mesh mesh used for the analysis
 * @param {object} o.result result of `solvePlane` on that mesh
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} [o.thickness] thickness [m] (plane modes; default 1)
 * @param {'stress'|'strain'|'axisym'} [o.mode] kinematic assumption used in the analysis
 * @returns {{elemError: Float64Array, errorNorm: number, energyNorm: number, eta: number}}
 *   elemError: element error energy norm [√J]; errorNorm: global ‖e‖ [√J]; energyNorm: ‖u‖ = √(2U) [√J];
 *   eta = ‖e‖ / √(‖u‖² + ‖e‖²) relative error [-]
 */
export function zzError({ mesh, result, E, nu, thickness = 1, mode = 'stress' }) {
  checkMesh(mesh, null, 'zzError');
  positive(E, 'E', 'zzError');
  if (!result || !result.gauss) throw new RangeError('zzError: result must come from solvePlane');
  const { type, nen, ne, conn, nodes } = mesh;
  const rs = recoveryScheme(type);
  const ng = rs.w.length;
  const sc = makeScratch(nen);
  const xe = new Float64Array(2 * nen);
  const comp = [result.stress.sxx, result.stress.syy, result.stress.szz, result.stress.sxy];
  const G = E / (2 * (1 + nu));
  const elemError = new Float64Array(ne);
  let tot = 0;
  const d = new Float64Array(4);
  for (let e = 0; e < ne; e++) {
    for (let i = 0; i < nen; i++) { const n = conn[e * nen + i]; xe[2 * i] = nodes[2 * n]; xe[2 * i + 1] = nodes[2 * n + 1]; }
    let s = 0;
    for (let g = 0; g < ng; g++) {
      geom(type, nen, rs.xi[g], rs.eta[g], xe, sc);
      const base = (e * ng + g) * 4;
      for (let r = 0; r < 4; r++) {
        let v = 0;
        for (let i = 0; i < nen; i++) v += sc.N[i] * comp[r][conn[e * nen + i]];
        d[r] = v - result.gauss.stress[base + r];
      }
      const dens = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - 2 * nu * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])) / E + (d[3] * d[3]) / G;
      s += rs.w[g] * sc.det * (mode === 'axisym' ? 2 * Math.PI * sc.r : thickness) * dens;
    }
    elemError[e] = Math.sqrt(Math.max(0, s));
    tot += Math.max(0, s);
  }
  const u2 = 2 * result.strainEnergy;
  return { elemError, errorNorm: Math.sqrt(tot), energyNorm: Math.sqrt(u2), eta: Math.sqrt(tot / (u2 + tot || 1)) };
}

/* ------------------------------------------------------------------------------------------------------------
 * Steady heat conduction
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Steady-state heat conduction ∇·(k ∇T) + Q = 0 on a 2-D mesh (plane with unit or given thickness, or
 * axisymmetric with x = r, y = z) with prescribed temperature, prescribed flux and convection boundaries.
 * @param {object} o
 * @param {object} o.mesh 2-D mesh
 * @param {number} o.k thermal conductivity [W/(m·K)]
 * @param {number} [o.thickness] thickness [m] (plane; default 1)
 * @param {boolean} [o.axisym] axisymmetric model (default false)
 * @param {Array<{nodes: ArrayLike<number>, T: number}>} [o.fixT] prescribed temperatures [K or °C, consistent]
 * @param {Array<{edges: number[][], q: number}>} [o.flux] heat flux into the body on boundary edges [W/m²]
 * @param {Array<{edges: number[][], h: number, Tinf: number}>} [o.convection] film coefficient [W/(m²·K)] and
 *   ambient temperature
 * @param {number} [o.source] volumetric heat generation [W/m³]
 * @param {'auto'|'skyline'|'pcg'} [o.solver] linear solver
 * @returns {{T: Float64Array, flux: {qx: Float64Array, qy: Float64Array}, heatIn: number, maxT: number, minT: number}}
 *   nodal temperature, element-centroid conductive flux [W/m²] and net heat input through prescribed-temperature
 *   nodes [W] (axisymmetric: over the full circumference)
 */
export function solveHeat2d(o) {
  const { mesh, k, thickness = 1, axisym = false, fixT = [], flux = [], convection = [], source = 0, solver = 'auto' } = o;
  const fn = 'solveHeat2d';
  checkMesh(mesh, null, fn);
  positive(k, 'k', fn);
  positive(thickness, 'thickness', fn);
  const { type, nen, nn, ne, conn, nodes } = mesh;
  const el = ELEM[type];
  const sc = makeScratch(nen);
  const xe = new Float64Array(2 * nen);
  const ke = new Float64Array(nen * nen);
  const TWO_PI = 2 * Math.PI;
  const width = (x) => (axisym ? TWO_PI * x : thickness);
  const loadXe = (e) => {
    for (let i = 0; i < nen; i++) { const n = conn[e * nen + i]; xe[2 * i] = nodes[2 * n]; xe[2 * i + 1] = nodes[2 * n + 1]; }
  };
  // Convection edges enter the operator as extra boundary "elements".
  const cedges = [];
  for (const c of convection) {
    if (!(c.h >= 0)) throw new RangeError(`${fn}: convection h must be >= 0`);
    for (const ed of c.edges) cedges.push({ ed, h: c.h, Tinf: c.Tinf });
  }
  const g2 = gauss1d(2);
  const g3 = gauss1d(3);
  const edgeShape = (ed, s) => {
    const quad = ed.length === 3;
    const N = quad ? [0.5 * s * (s - 1), 0.5 * s * (s + 1), 1 - s * s] : [0.5 * (1 - s), 0.5 * (1 + s)];
    const dN = quad ? [s - 0.5, s + 0.5, -2 * s] : [-0.5, 0.5];
    let x = 0; let tx = 0; let ty = 0;
    for (let i = 0; i < ed.length; i++) { x += N[i] * nodes[2 * ed[i]]; tx += dN[i] * nodes[2 * ed[i]]; ty += dN[i] * nodes[2 * ed[i] + 1]; }
    return { N, ds: Math.hypot(tx, ty) * width(x) };
  };
  const kedge = new Float64Array(9);
  const elemNodes = (e) => (e < ne ? conn.subarray(e * nen, e * nen + nen) : cedges[e - ne].ed);
  const elemMatrix = (e) => {
    if (e >= ne) {
      const c = cedges[e - ne];
      const m = c.ed.length;
      const out = m === 3 ? kedge : kedge.subarray(0, 4);
      out.fill(0);
      const g = m === 3 ? g3 : g2;
      for (let q = 0; q < g.x.length; q++) {
        const s = edgeShape(c.ed, g.x[q]);
        for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) out[i * m + j] += c.h * s.N[i] * s.N[j] * s.ds * g.w[q];
      }
      return out;
    }
    loadXe(e);
    ke.fill(0);
    const rule = el.rule;
    for (let g = 0; g < rule.w.length; g++) {
      geom(type, nen, rule.xi[g], rule.eta[g], xe, sc);
      const w = rule.w[g] * sc.det * width(sc.r) * k;
      for (let i = 0; i < nen; i++) for (let j = 0; j < nen; j++) ke[i * nen + j] += w * (sc.dx[i] * sc.dx[j] + sc.dy[i] * sc.dy[j]);
    }
    return ke;
  };
  const f = new Float64Array(nn);
  for (const c of cedges) {
    const g = c.ed.length === 3 ? g3 : g2;
    for (let q = 0; q < g.x.length; q++) {
      const s = edgeShape(c.ed, g.x[q]);
      for (let i = 0; i < c.ed.length; i++) f[c.ed[i]] += c.h * c.Tinf * s.N[i] * s.ds * g.w[q];
    }
  }
  for (const L of flux) {
    for (const ed of L.edges) {
      const g = ed.length === 3 ? g3 : g2;
      for (let q = 0; q < g.x.length; q++) {
        const s = edgeShape(ed, g.x[q]);
        for (let i = 0; i < ed.length; i++) f[ed[i]] += L.q * s.N[i] * s.ds * g.w[q];
      }
    }
  }
  if (source) {
    const rule = el.mass;
    for (let e = 0; e < ne; e++) {
      loadXe(e);
      for (let g = 0; g < rule.w.length; g++) {
        geom(type, nen, rule.xi[g], rule.eta[g], xe, sc);
        const w = rule.w[g] * sc.det * width(sc.r) * source;
        for (let i = 0; i < nen; i++) f[conn[e * nen + i]] += sc.N[i] * w;
      }
    }
  }
  const fixed = buildFixed(nn, ['T'], fixT, fn);
  if (fixed.every(Number.isNaN) && cedges.length === 0) throw new Error('Heat model needs at least one prescribed temperature or convection boundary.');
  const sys = feSystem({ nNodes: nn, dpn: 1, nElem: ne + cedges.length, elemNodes, elemMatrix, fixed, solver });
  const T = sys.solve(f);
  const kT = sys.mulK(T);
  let heatIn = 0;
  for (let i = 0; i < nn; i++) if (sys.isFixed[i]) heatIn += kT[i] - f[i];
  const qx = new Float64Array(ne);
  const qy = new Float64Array(ne);
  const cx = el.tri ? 1 / 3 : 0;
  for (let e = 0; e < ne; e++) {
    loadXe(e);
    geom(type, nen, cx, cx, xe, sc);
    let gx = 0; let gy = 0;
    for (let i = 0; i < nen; i++) { const t = T[conn[e * nen + i]]; gx += sc.dx[i] * t; gy += sc.dy[i] * t; }
    qx[e] = -k * gx; qy[e] = -k * gy;
  }
  let maxT = -Infinity; let minT = Infinity;
  for (let i = 0; i < nn; i++) { if (T[i] > maxT) maxT = T[i]; if (T[i] < minT) minT = T[i]; }
  return { T, flux: { qx, qy }, heatIn, maxT, minT };
}

/* ------------------------------------------------------------------------------------------------------------
 * Mindlin–Reissner plate (MITC4)
 * ---------------------------------------------------------------------------------------------------------- */

function makePlateWork() {
  return {
    sc: makeScratch(4), ke: new Float64Array(144), Bb: new Float64Array(36), Bs: new Float64Array(24),
    tie: new Float64Array(48), cov: new Float64Array(24),
  };
}

/**
 * MITC4 plate-bending stiffness. Nodal dofs (w, rx, ry): rx, ry are right-handed rotations about the x and y
 * axes, so the fibre displacements are u = z·ry, v = −z·rx and the shear strains γxz = w,x + ry, γyz = w,y − rx.
 * The covariant shear strains are sampled at the four edge mid-points and interpolated linearly, which removes
 * shear locking without spurious zero-energy modes.
 * @param {Float64Array} xe nodal coordinates (x, y) × 4 in the element plane [m]
 * @param {number} Db bending rigidity E t³ / (12 (1 − ν²)) [N·m]
 * @param {number} nu Poisson's ratio [-]
 * @param {number} Ds shear rigidity κ G t [N/m]
 * @param {object} wk work arrays from makePlateWork()
 * @returns {Float64Array} 12×12 stiffness (row-major, reused buffer)
 */
function mitc4Ke(xe, Db, nu, Ds, wk) {
  const { sc, ke, Bb, Bs, tie } = wk;
  ke.fill(0);
  // Tying points: A (0, 1) and C (0, -1) carry γξ; B (-1, 0) and D (1, 0) carry γη.
  const tp = [[0, 1], [0, -1], [-1, 0], [1, 0]];
  for (let k = 0; k < 4; k++) {
    shape2d('Q4', tp[k][0], tp[k][1], sc.N, sc.dxi, sc.deta);
    const d = k < 2 ? sc.dxi : sc.deta;
    let xs = 0; let ys = 0;
    for (let i = 0; i < 4; i++) { xs += d[i] * xe[2 * i]; ys += d[i] * xe[2 * i + 1]; }
    for (let i = 0; i < 4; i++) {
      tie[k * 12 + 3 * i] = d[i];
      tie[k * 12 + 3 * i + 1] = -sc.N[i] * ys;
      tie[k * 12 + 3 * i + 2] = sc.N[i] * xs;
    }
  }
  const rule = ELEM.Q4.rule;
  const d2 = 0.5 * (1 - nu) * Db;
  for (let g = 0; g < 4; g++) {
    const xi = rule.xi[g];
    const eta = rule.eta[g];
    geom('Q4', 4, xi, eta, xe, sc);
    const w = rule.w[g] * sc.det;
    for (let i = 0; i < 4; i++) {
      Bb[3 * i] = 0; Bb[3 * i + 1] = 0; Bb[3 * i + 2] = sc.dx[i];
      Bb[12 + 3 * i] = 0; Bb[12 + 3 * i + 1] = -sc.dy[i]; Bb[12 + 3 * i + 2] = 0;
      Bb[24 + 3 * i] = 0; Bb[24 + 3 * i + 1] = -sc.dx[i]; Bb[24 + 3 * i + 2] = sc.dy[i];
    }
    // Covariant shear interpolation, then to Cartesian with J⁻¹ (ji = [ξ,x  η,x; ξ,y  η,y]).
    const a = 0.5 * (1 + eta); const c = 0.5 * (1 - eta);
    const b = 0.5 * (1 - xi); const d = 0.5 * (1 + xi);
    for (let q = 0; q < 12; q++) {
      const gxi = a * tie[q] + c * tie[12 + q];
      const get = b * tie[24 + q] + d * tie[36 + q];
      Bs[q] = sc.ji[0] * gxi + sc.ji[1] * get;
      Bs[12 + q] = sc.ji[2] * gxi + sc.ji[3] * get;
    }
    for (let p = 0; p < 12; p++) {
      const m0 = Db * (Bb[p] + nu * Bb[12 + p]);
      const m1 = Db * (nu * Bb[p] + Bb[12 + p]);
      const m2 = d2 * Bb[24 + p];
      const s0 = Ds * Bs[p];
      const s1 = Ds * Bs[12 + p];
      for (let q = 0; q < 12; q++) {
        ke[p * 12 + q] += w * (m0 * Bb[q] + m1 * Bb[12 + q] + m2 * Bb[24 + q] + s0 * Bs[q] + s1 * Bs[12 + q]);
      }
    }
  }
  return ke;
}

function plateSetup(o, fn) {
  const { mesh, E, nu, t, kShear = 5 / 6 } = o;
  checkMesh(mesh, ['Q4'], fn);
  positive(E, 'E', fn);
  checkNu(nu, fn);
  positive(t, 't', fn);
  positive(kShear, 'kShear', fn);
  const Db = (E * t * t * t) / (12 * (1 - nu * nu));
  const Ds = (kShear * E * t) / (2 * (1 + nu));
  const wk = makePlateWork();
  const xe = new Float64Array(8);
  const { conn, nodes, nn, ne } = mesh;
  const loadXe = (e) => {
    for (let i = 0; i < 4; i++) { const n = conn[e * 4 + i]; xe[2 * i] = nodes[2 * n]; xe[2 * i + 1] = nodes[2 * n + 1]; }
    return xe;
  };
  const elemNodes = (e) => conn.subarray(e * 4, e * 4 + 4);
  const elemMatrix = (e) => mitc4Ke(loadXe(e), Db, nu, Ds, wk);
  return { Db, Ds, wk, xe, loadXe, elemNodes, elemMatrix, nn, ne, conn, nodes };
}

/**
 * Static bending of a moderately thick or thin plate with MITC4 Mindlin–Reissner elements.
 * @param {object} o
 * @param {object} o.mesh Q4 mesh in the plate plane
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} o.t plate thickness [m]
 * @param {number} [o.kShear] shear correction factor (default 5/6)
 * @param {number|Function} [o.pressure] lateral pressure [Pa], positive in +z; number or function of (x, y)
 * @param {Array<{node: number, fz?: number, mx?: number, my?: number}>} [o.pointLoads] nodal force [N] and
 *   moments [N·m] about x and y
 * @param {Array<{nodes: ArrayLike<number>, w?: number, rx?: number, ry?: number}>} o.fix supports: w [m],
 *   rx, ry [rad] (rotations about x and y). Simply supported edge: w = 0 (hard support adds the rotation about
 *   the edge normal); clamped: w = rx = ry = 0
 * @param {number} [o.foundation] Winkler foundation modulus [N/m³] (lumped at nodes)
 * @param {'auto'|'skyline'|'pcg'} [o.solver] linear solver (direct recommended for thin plates)
 * @returns {{w: Float64Array, rx: Float64Array, ry: Float64Array, u: Float64Array, reactions: Float64Array,
 *   reactionSum: number, appliedSum: number, moments: {mx: Float64Array, my: Float64Array, mxy: Float64Array},
 *   shear: {qx: Float64Array, qy: Float64Array}, strainEnergy: number, maxW: number}}
 *   w [m], rotations [rad]; reactions per dof [N, N·m]; nodal-averaged moments per unit width [N·m/m]
 *   (mx = D(κx + ν κy) with κx = ∂ry/∂x, positive for sagging under +z load), element-centre shear forces
 *   [N/m], strain energy [J]
 */
export function solvePlate(o) {
  const fn = 'solvePlate';
  const S = plateSetup(o, fn);
  const { mesh, nu, pressure = 0, pointLoads = [], fix, foundation = 0, solver = 'auto' } = o;
  const { nn, ne, conn, wk, Db, Ds } = S;
  const sc = wk.sc;
  const f = new Float64Array(3 * nn);
  const springs = foundation ? new Float64Array(3 * nn) : null;
  const rule = ELEM.Q4.rule;
  if (pressure || foundation) {
    for (let e = 0; e < ne; e++) {
      const xe = S.loadXe(e);
      for (let g = 0; g < 4; g++) {
        geom('Q4', 4, rule.xi[g], rule.eta[g], xe, sc);
        const p = evalLoad(pressure, sc.x, sc.y);
        for (let i = 0; i < 4; i++) {
          f[3 * conn[e * 4 + i]] += sc.N[i] * p * sc.det * rule.w[g];
          if (springs) springs[3 * conn[e * 4 + i]] += sc.N[i] * foundation * sc.det * rule.w[g];
        }
      }
    }
  }
  for (const L of pointLoads) {
    if (!(L.node >= 0 && L.node < nn)) throw new RangeError(`${fn}: pointLoads node index out of range`);
    f[3 * L.node] += L.fz || 0; f[3 * L.node + 1] += L.mx || 0; f[3 * L.node + 2] += L.my || 0;
  }
  const fixed = buildFixed(nn, ['w', 'rx', 'ry'], fix, fn);
  if (fixed.every(Number.isNaN) && !foundation) throw new Error('Plate has no supports.');
  const sys = feSystem({ nNodes: nn, dpn: 3, nElem: ne, elemNodes: S.elemNodes, elemMatrix: S.elemMatrix, fixed, springs, solver });
  const u = sys.solve(f);
  const ku = sys.mulK(u);
  const reactions = new Float64Array(3 * nn);
  let energy = 0; let rsum = 0; let asum = 0;
  for (let i = 0; i < 3 * nn; i++) {
    energy += 0.5 * u[i] * ku[i];
    if (sys.isFixed[i]) { reactions[i] = ku[i] - f[i]; if (i % 3 === 0) rsum += reactions[i]; }
    if (i % 3 === 0) asum += f[i];
  }
  const w = new Float64Array(nn); const rx = new Float64Array(nn); const ry = new Float64Array(nn);
  let maxW = 0;
  for (let n = 0; n < nn; n++) { w[n] = u[3 * n]; rx[n] = u[3 * n + 1]; ry[n] = u[3 * n + 2]; maxW = Math.max(maxW, Math.abs(w[n])); }

  const rs = recoveryScheme('Q4');
  const acc = new Float64Array(3 * nn);
  const cnt = new Uint16Array(nn);
  const gm = new Float64Array(12);
  const qx = new Float64Array(ne); const qy = new Float64Array(ne);
  for (let e = 0; e < ne; e++) {
    const xe = S.loadXe(e);
    for (let g = 0; g < 4; g++) {
      geom('Q4', 4, rs.xi[g], rs.eta[g], xe, sc);
      let kx = 0; let ky = 0; let kxy = 0;
      for (let i = 0; i < 4; i++) {
        const n = conn[e * 4 + i];
        kx += sc.dx[i] * ry[n]; ky -= sc.dy[i] * rx[n]; kxy += sc.dy[i] * ry[n] - sc.dx[i] * rx[n];
      }
      gm[3 * g] = Db * (kx + nu * ky); gm[3 * g + 1] = Db * (nu * kx + ky); gm[3 * g + 2] = 0.5 * (1 - nu) * Db * kxy;
    }
    for (let i = 0; i < 4; i++) {
      const n = conn[e * 4 + i];
      cnt[n]++;
      for (let g = 0; g < 4; g++) for (let r = 0; r < 3; r++) acc[3 * n + r] += rs.ex[i][g] * gm[3 * g + r];
    }
    // Centre shear strains from the centre-evaluated compatible field (adequate for reporting).
    geom('Q4', 4, 0, 0, xe, sc);
    let gxz = 0; let gyz = 0;
    for (let i = 0; i < 4; i++) {
      const n = conn[e * 4 + i];
      gxz += sc.dx[i] * w[n] + sc.N[i] * ry[n];
      gyz += sc.dy[i] * w[n] - sc.N[i] * rx[n];
    }
    qx[e] = Ds * gxz; qy[e] = Ds * gyz;
  }
  const mx = new Float64Array(nn); const my = new Float64Array(nn); const mxy = new Float64Array(nn);
  for (let n = 0; n < nn; n++) { const c = cnt[n] || 1; mx[n] = acc[3 * n] / c; my[n] = acc[3 * n + 1] / c; mxy[n] = acc[3 * n + 2] / c; }
  return { w, rx, ry, u, reactions, reactionSum: rsum, appliedSum: asum, moments: { mx, my, mxy }, shear: { qx, qy }, strainEnergy: energy, maxW };
}

/**
 * Largest eigenpairs of B x = μ K x by subspace iteration with a factorised K (K positive definite after
 * constraints, B symmetric and possibly indefinite or semi-definite). Used for buckling (B = geometric
 * stiffness) and vibration (B = mass), where μ = 1/λ. Reference: Bathe, "Finite Element Procedures" §11.6.
 * @param {object} o
 * @param {object} o.sys system from `feSystem` (direct solver recommended)
 * @param {(x: Float64Array) => Float64Array} o.mulB product B·x with rows of prescribed dofs zeroed
 * @param {number} o.n number of dofs
 * @param {number} o.nModes number of eigenpairs wanted
 * @param {number} [o.tol] relative eigenvalue tolerance (default 1e-9)
 * @param {number} [o.maxIter] iteration cap (default 80)
 * @returns {{mu: number[], vectors: Float64Array[], iterations: number}} eigenvalues μ (descending) and unit
 *   2-norm eigenvectors
 */
export function subspaceModes({ sys, mulB, n, nModes, tol = 1e-9, maxIter = 80 }) {
  const free = [];
  for (let i = 0; i < n; i++) if (!sys.isFixed[i]) free.push(i);
  const p = Math.min(free.length, Math.max(2 * nModes, nModes + 6));
  if (p < 1) throw new Error('No free degrees of freedom for the eigen-analysis.');
  let X = [];
  // Deterministic, linearly independent start vectors (low-order sinusoids over the free dofs).
  for (let k = 0; k < p; k++) {
    const v = new Float64Array(n);
    for (let j = 0; j < free.length; j++) v[free[j]] = Math.cos(((k + 0.5) * (j + 0.37) * Math.PI) / free.length) + (k === 0 ? 1 : 0) + 0.01 * ((j * 7919 + k * 104729) % 97) / 97;
    X.push(v);
  }
  let prev = null;
  let mu = [];
  let iterations = 0;
  for (let it = 0; it < maxIter; it++) {
    iterations = it + 1;
    const Y = X.map((x) => { const b = mulB(x); return sys.solveHomogeneous(b); });
    const KY = Y.map((y) => { const k = sys.mulK(y); for (let i = 0; i < n; i++) if (sys.isFixed[i]) k[i] = 0; return k; });
    const BY = Y.map((y) => mulB(y));
    const Kr = [];
    const Br = [];
    for (let a = 0; a < p; a++) {
      Kr.push(new Array(p).fill(0));
      Br.push(new Array(p).fill(0));
    }
    for (let a = 0; a < p; a++) for (let b = a; b < p; b++) {
      let sk = 0; let sb = 0;
      const ya = Y[a]; const kb = KY[b]; const bb = BY[b];
      for (let i = 0; i < n; i++) { sk += ya[i] * kb[i]; sb += ya[i] * bb[i]; }
      Kr[a][b] = sk; Kr[b][a] = sk; Br[a][b] = sb; Br[b][a] = sb;
    }
    // Scale so the reduced K has unit diagonal; keeps the small generalised problem well conditioned.
    const sc = Kr.map((row, a) => 1 / Math.sqrt(row[a] > 0 ? row[a] : 1));
    for (let a = 0; a < p; a++) for (let b = 0; b < p; b++) { Kr[a][b] *= sc[a] * sc[b]; Br[a][b] *= sc[a] * sc[b]; }
    const eg = eigGenSym(Br, Kr, {});
    const idx = eg.values.map((v, i) => i).sort((i, j) => eg.values[j] - eg.values[i]);
    mu = idx.map((i) => eg.values[i]);
    const Xn = [];
    for (let k = 0; k < p; k++) {
      const z = eg.vectors[idx[k]];
      const v = new Float64Array(n);
      for (let a = 0; a < p; a++) { const c = z[a] * sc[a]; const ya = Y[a]; for (let i = 0; i < n; i++) v[i] += c * ya[i]; }
      let nrm = 0;
      for (let i = 0; i < n; i++) nrm += v[i] * v[i];
      nrm = Math.sqrt(nrm) || 1;
      for (let i = 0; i < n; i++) v[i] /= nrm;
      Xn.push(v);
    }
    X = Xn;
    if (prev) {
      let conv = true;
      for (let k = 0; k < nModes; k++) if (Math.abs(mu[k] - prev[k]) > tol * Math.abs(mu[k])) { conv = false; break; }
      if (conv) break;
    }
    prev = mu.slice();
  }
  return { mu: mu.slice(0, nModes), vectors: X.slice(0, nModes), iterations };
}

/**
 * Linear (bifurcation) buckling of a plate under uniform in-plane membrane forces, MITC4 elements.
 * Solves (K − λ Kg) φ = 0 with Kg = ∫ [w,x w,y] [[Nx, Nxy], [Nxy, Ny]] [w,x; w,y] dA for the smallest positive λ.
 * @param {object} o
 * @param {object} o.mesh Q4 mesh
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} o.t thickness [m]
 * @param {number} [o.Nx] reference compressive force per unit width in x [N/m] (compression positive)
 * @param {number} [o.Ny] reference compressive force per unit width in y [N/m]
 * @param {number} [o.Nxy] reference shear flow [N/m]
 * @param {Array<object>} o.fix supports as in `solvePlate`
 * @param {number} [o.nModes] number of positive load factors wanted (default 1)
 * @returns {{factors: number[], modes: Float64Array[], iterations: number}} load factors λ [-] on the reference
 *   forces (ascending, positive) and mode vectors (w, rx, ry per node)
 */
export function bucklePlate(o) {
  const fn = 'bucklePlate';
  const S = plateSetup(o, fn);
  const { Nx = 0, Ny = 0, Nxy = 0, fix, nModes = 1 } = o;
  if (Nx === 0 && Ny === 0 && Nxy === 0) throw new RangeError(`${fn}: at least one of Nx, Ny, Nxy must be non-zero`);
  const { nn, ne, conn, wk } = S;
  const sc = wk.sc;
  const n = 3 * nn;
  const fixed = buildFixed(nn, ['w', 'rx', 'ry'], fix, fn);
  const sys = feSystem({ nNodes: nn, dpn: 3, nElem: ne, elemNodes: S.elemNodes, elemMatrix: S.elemMatrix, fixed, solver: 'skyline' });
  // Element geometric matrices act on w only; store 4×4 blocks once.
  const kg = new Float64Array(ne * 16);
  const rule = ELEM.Q4.rule;
  for (let e = 0; e < ne; e++) {
    const xe = S.loadXe(e);
    for (let g = 0; g < 4; g++) {
      geom('Q4', 4, rule.xi[g], rule.eta[g], xe, sc);
      const w = rule.w[g] * sc.det;
      for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
        kg[e * 16 + i * 4 + j] += w * (Nx * sc.dx[i] * sc.dx[j] + Ny * sc.dy[i] * sc.dy[j] + Nxy * (sc.dx[i] * sc.dy[j] + sc.dy[i] * sc.dx[j]));
      }
    }
  }
  const mulB = (x) => {
    const y = new Float64Array(n);
    for (let e = 0; e < ne; e++) {
      for (let i = 0; i < 4; i++) {
        let s = 0;
        for (let j = 0; j < 4; j++) s += kg[e * 16 + i * 4 + j] * x[3 * conn[e * 4 + j]];
        y[3 * conn[e * 4 + i]] += s;
      }
    }
    for (let i = 0; i < n; i++) if (sys.isFixed[i]) y[i] = 0;
    return y;
  };
  // Shear buckling has ± pairs; ask for extra modes so enough positive ones survive.
  const want = Nx === 0 && Ny === 0 ? 2 * nModes + 2 : nModes + 2;
  const r = subspaceModes({ sys, mulB, n, nModes: want });
  const factors = [];
  const modes = [];
  for (let k = 0; k < r.mu.length && factors.length < nModes; k++) {
    if (r.mu[k] > 1e-14) { factors.push(1 / r.mu[k]); modes.push(r.vectors[k]); }
  }
  if (factors.length === 0) throw new Error('No positive buckling factor found: the reference loading is not compressive.');
  return { factors, modes, iterations: r.iterations };
}

/**
 * Natural frequencies of a plate (MITC4, consistent mass including rotary inertia).
 * @param {object} o
 * @param {object} o.mesh Q4 mesh
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} o.t thickness [m]
 * @param {number} o.rho density [kg/m³]
 * @param {number} [o.addedMass] non-structural mass per unit area [kg/m²]
 * @param {Array<object>} o.fix supports as in `solvePlate`
 * @param {number} [o.nModes] number of modes (default 4)
 * @returns {{frequencies: number[], omega: number[], modes: Float64Array[]}} natural frequencies [Hz], circular
 *   frequencies [rad/s], mode vectors (w, rx, ry per node)
 */
export function modalPlate(o) {
  const fn = 'modalPlate';
  const S = plateSetup(o, fn);
  const { rho, addedMass = 0, fix, nModes = 4, t } = o;
  positive(rho, 'rho', fn);
  const { nn, ne, conn, wk } = S;
  const sc = wk.sc;
  const n = 3 * nn;
  const fixed = buildFixed(nn, ['w', 'rx', 'ry'], fix, fn);
  const sys = feSystem({ nNodes: nn, dpn: 3, nElem: ne, elemNodes: S.elemNodes, elemMatrix: S.elemMatrix, fixed, solver: 'skyline' });
  const mm = new Float64Array(ne * 16);
  const rule = ELEM.Q4.rule;
  for (let e = 0; e < ne; e++) {
    const xe = S.loadXe(e);
    for (let g = 0; g < 4; g++) {
      geom('Q4', 4, rule.xi[g], rule.eta[g], xe, sc);
      const w = rule.w[g] * sc.det;
      for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) mm[e * 16 + i * 4 + j] += w * sc.N[i] * sc.N[j];
    }
  }
  const mw = rho * t + addedMass;
  const mr = (rho * t * t * t) / 12;
  const mulB = (x) => {
    const y = new Float64Array(n);
    for (let e = 0; e < ne; e++) {
      for (let i = 0; i < 4; i++) {
        let s0 = 0; let s1 = 0; let s2 = 0;
        for (let j = 0; j < 4; j++) {
          const m = mm[e * 16 + i * 4 + j];
          const d = 3 * conn[e * 4 + j];
          s0 += m * x[d]; s1 += m * x[d + 1]; s2 += m * x[d + 2];
        }
        const d = 3 * conn[e * 4 + i];
        y[d] += mw * s0; y[d + 1] += mr * s1; y[d + 2] += mr * s2;
      }
    }
    for (let i = 0; i < n; i++) if (sys.isFixed[i]) y[i] = 0;
    return y;
  };
  const r = subspaceModes({ sys, mulB, n, nModes });
  const omega = r.mu.map((m) => 1 / Math.sqrt(m));
  return { frequencies: omega.map((w) => w / (2 * Math.PI)), omega, modes: r.vectors };
}

/* ------------------------------------------------------------------------------------------------------------
 * Flat facet shell
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Structured 3-D shell mesh of a plate stiffened by equally spaced flat-bar or T stiffeners running along x.
 * The plate lies in z = 0; webs extend to z = −hw.
 * @param {object} o
 * @param {number} o.a panel length along the stiffeners [m]
 * @param {number} o.b panel width across the stiffeners [m]
 * @param {number} o.nStiff number of stiffeners (spacing b / (nStiff + 1))
 * @param {number} o.hw web height [m]
 * @param {number} o.tp plate thickness [m]
 * @param {number} o.tw web thickness [m]
 * @param {number} [o.bf] flange width [m] (0 = flat bar)
 * @param {number} [o.tf] flange thickness [m]
 * @param {number} o.nx elements along the length
 * @param {number} o.nyPerBay plate elements between stiffeners
 * @param {number} [o.nWeb] elements over the web height (default 2)
 * @returns {{nodes: Float64Array, quads: Int32Array, thickness: Float64Array, nn: number, ne: number,
 *   elemSets: {plate: Int32Array, web: Int32Array, flange: Int32Array}, nodeSets: {x0: Int32Array, x1: Int32Array,
 *   y0: Int32Array, y1: Int32Array, plate: Int32Array}}} shell mesh: nodes (x, y, z) [m], quads (4 nodes each),
 *   per-element thickness [m]; node sets on the four panel edges (x-ends include web and flange nodes)
 */
export function stiffenedPanelMesh(o) {
  const { a, b, nStiff, hw, tp, tw, bf = 0, tf = 0, nx, nyPerBay, nWeb = 2 } = o;
  const fn = 'stiffenedPanelMesh';
  positive(a, 'a', fn); positive(b, 'b', fn); positive(tp, 'tp', fn);
  if (!(Number.isInteger(nStiff) && nStiff >= 0)) throw new RangeError(`${fn}: nStiff must be a non-negative integer`);
  if (nStiff > 0) { positive(hw, 'hw', fn); positive(tw, 'tw', fn); }
  if (bf > 0) positive(tf, 'tf', fn);
  const ny = nyPerBay * (nStiff + 1);
  const coords = [];
  const quads = [];
  const thick = [];
  const sets = { plate: [], web: [], flange: [] };
  const add = (x, y, z) => { coords.push(x, y, z); return coords.length / 3 - 1; };
  const pid = (i, j) => j * (nx + 1) + i;
  for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++) add((a * i) / nx, (b * j) / ny, 0);
  const nPlateNodes = coords.length / 3;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    sets.plate.push(quads.length / 4);
    quads.push(pid(i, j), pid(i + 1, j), pid(i + 1, j + 1), pid(i, j + 1));
    thick.push(tp);
  }
  const extraX0 = [];
  const extraX1 = [];
  for (let s = 1; s <= nStiff; s++) {
    const j = s * nyPerBay;
    const y = (b * j) / ny;
    let prev = [];
    for (let i = 0; i <= nx; i++) prev.push(pid(i, j));
    for (let k = 1; k <= nWeb; k++) {
      const row = [];
      for (let i = 0; i <= nx; i++) row.push(add((a * i) / nx, y, (-hw * k) / nWeb));
      extraX0.push(row[0]); extraX1.push(row[nx]);
      for (let i = 0; i < nx; i++) {
        sets.web.push(quads.length / 4);
        quads.push(prev[i], prev[i + 1], row[i + 1], row[i]);
        thick.push(tw);
      }
      prev = row;
    }
    if (bf > 0) {
      const lo = []; const hi = [];
      for (let i = 0; i <= nx; i++) { lo.push(add((a * i) / nx, y - bf / 2, -hw)); hi.push(add((a * i) / nx, y + bf / 2, -hw)); }
      extraX0.push(lo[0], hi[0]); extraX1.push(lo[nx], hi[nx]);
      for (let i = 0; i < nx; i++) {
        sets.flange.push(quads.length / 4, quads.length / 4 + 1);
        quads.push(lo[i], lo[i + 1], prev[i + 1], prev[i]);
        quads.push(prev[i], prev[i + 1], hi[i + 1], hi[i]);
        thick.push(tf, tf);
      }
    }
  }
  const x0 = []; const x1 = []; const y0 = []; const y1 = []; const plate = [];
  for (let j = 0; j <= ny; j++) { x0.push(pid(0, j)); x1.push(pid(nx, j)); }
  for (let i = 0; i <= nx; i++) { y0.push(pid(i, 0)); y1.push(pid(i, ny)); }
  for (let i = 0; i < nPlateNodes; i++) plate.push(i);
  return {
    nodes: Float64Array.from(coords), quads: Int32Array.from(quads), thickness: Float64Array.from(thick),
    nn: coords.length / 3, ne: quads.length / 4,
    elemSets: { plate: Int32Array.from(sets.plate), web: Int32Array.from(sets.web), flange: Int32Array.from(sets.flange) },
    nodeSets: { x0: Int32Array.from(x0.concat(extraX0)), x1: Int32Array.from(x1.concat(extraX1)), y0: Int32Array.from(y0), y1: Int32Array.from(y1), plate: Int32Array.from(plate) },
  };
}

/**
 * Linear static analysis of a 3-D shell structure meshed with flat four-node facets: Q4 membrane with
 * incompatible modes, MITC4 plate bending and a small drilling stiffness that leaves rigid in-plane rotation
 * free (Zienkiewicz & Taylor, vol. 2, flat-shell approach). Warped quadrilaterals are projected onto their mean
 * plane; keep warping small (curved shells converge with mesh refinement).
 * @param {object} o
 * @param {Float64Array|number[]} o.nodes nodal coordinates (x, y, z) [m]
 * @param {Int32Array|number[]} o.quads element connectivity, 4 nodes per element
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number|ArrayLike<number>} o.t thickness [m], uniform or per element
 * @param {Array<object>} o.fix supports { nodes, ux?, uy?, uz? [m], rx?, ry?, rz? [rad] } in global axes
 * @param {Array<object>} [o.nodalLoads] { node, fx, fy, fz [N], mx, my, mz [N·m] }
 * @param {number|ArrayLike<number>} [o.pressure] pressure [Pa] along each element's +normal (right-hand rule on
 *   the node order), uniform or per element
 * @param {{rho: number, g: [number, number, number]}} [o.gravity] self-weight: density [kg/m³], acceleration [m/s²]
 * @param {[number, number, number]} [o.surfaceLoad] load per unit shell area in global axes [N/m²]
 * @param {number} [o.kShear] shear correction factor (default 5/6)
 * @param {number} [o.drillingFactor] drilling stiffness as a fraction of E·t·A (default 1e-3; typical 1e-4…3e-2)
 * @param {'auto'|'skyline'|'pcg'} [o.solver] linear solver
 * @returns {{u: Float64Array, reactions: Float64Array, reactionSum: {fx: number, fy: number, fz: number},
 *   appliedSum: {fx: number, fy: number, fz: number}, elem: {nx: Float64Array, ny: Float64Array, nxy: Float64Array,
 *   mx: Float64Array, my: Float64Array, mxy: Float64Array, vmTop: Float64Array, vmBottom: Float64Array,
 *   vmMid: Float64Array}, nodalVonMises: Float64Array, maxVonMises: number, maxDisp: number, strainEnergy: number}}
 *   u: 6 dofs per node [m, rad]; element-centre membrane forces [N/m] and moments [N·m/m] in element axes;
 *   von Mises stress at top/bottom/mid surface [Pa]; nodal average of the larger surface value [Pa]
 */
export function solveShell(o) {
  const fn = 'solveShell';
  const { E, nu, fix, nodalLoads = [], pressure = 0, gravity = null, surfaceLoad = null, kShear = 5 / 6, drillingFactor = 1e-3, solver = 'auto' } = o;
  positive(E, 'E', fn);
  checkNu(nu, fn);
  positive(drillingFactor, 'drillingFactor', fn);
  const nodes = o.nodes;
  const quads = o.quads;
  if (!nodes || nodes.length % 3 !== 0) throw new RangeError(`${fn}: nodes must hold 3 coordinates per node`);
  if (!quads || quads.length % 4 !== 0 || quads.length === 0) throw new RangeError(`${fn}: quads must hold 4 nodes per element`);
  const nn = nodes.length / 3;
  const ne = quads.length / 4;
  const tOf = (e) => (typeof o.t === 'number' ? o.t : o.t[e]);
  const pOf = (e) => (typeof pressure === 'number' ? pressure : pressure[e]);
  for (let e = 0; e < ne; e++) if (!(tOf(e) > 0)) throw new RangeError(`${fn}: t must be positive for every element`);
  const D = planeD(E, nu, 'stress');
  const q4w = makeQ4Work();
  const pw = makePlateWork();
  const xl = new Float64Array(8);
  const T = new Float64Array(9);
  const kl = new Float64Array(576);
  const kg = new Float64Array(576);
  const tmp = new Float64Array(9);
  const frame = (e) => {
    const p = [0, 1, 2, 3].map((i) => { const n = quads[e * 4 + i]; return [nodes[3 * n], nodes[3 * n + 1], nodes[3 * n + 2]]; });
    const e1 = [0, 0, 0]; const v = [0, 0, 0]; const c = [0, 0, 0];
    for (let k = 0; k < 3; k++) {
      e1[k] = p[1][k] + p[2][k] - p[0][k] - p[3][k];
      v[k] = p[2][k] + p[3][k] - p[0][k] - p[1][k];
      c[k] = 0.25 * (p[0][k] + p[1][k] + p[2][k] + p[3][k]);
    }
    const n3 = [e1[1] * v[2] - e1[2] * v[1], e1[2] * v[0] - e1[0] * v[2], e1[0] * v[1] - e1[1] * v[0]];
    const l1 = Math.hypot(e1[0], e1[1], e1[2]);
    const l3 = Math.hypot(n3[0], n3[1], n3[2]);
    if (!(l1 > 0 && l3 > 0)) throw new Error(`Shell element ${e} is degenerate.`);
    for (let k = 0; k < 3; k++) { e1[k] /= l1; n3[k] /= l3; }
    const e2 = [n3[1] * e1[2] - n3[2] * e1[1], n3[2] * e1[0] - n3[0] * e1[2], n3[0] * e1[1] - n3[1] * e1[0]];
    for (let k = 0; k < 3; k++) { T[k] = e1[k]; T[3 + k] = e2[k]; T[6 + k] = n3[k]; }
    for (let i = 0; i < 4; i++) {
      const d0 = p[i][0] - c[0]; const d1 = p[i][1] - c[1]; const d2 = p[i][2] - c[2];
      xl[2 * i] = d0 * e1[0] + d1 * e1[1] + d2 * e1[2];
      xl[2 * i + 1] = d0 * e2[0] + d1 * e2[1] + d2 * e2[2];
    }
  };
  const area = () => 0.5 * ((xl[4] - xl[0]) * (xl[7] - xl[3]) - (xl[6] - xl[2]) * (xl[5] - xl[1]));
  const localK = (e) => {
    frame(e);
    const t = tOf(e);
    kl.fill(0);
    const km = q4IncompatibleKe(xl, D, t, q4w);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let a = 0; a < 2; a++) for (let b = 0; b < 2; b++) {
      kl[(6 * i + a) * 24 + 6 * j + b] = km[(2 * i + a) * 8 + 2 * j + b];
    }
    const kp = mitc4Ke(xl, (E * t * t * t) / (12 * (1 - nu * nu)), nu, (kShear * E * t) / (2 * (1 + nu)), pw);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) {
      kl[(6 * i + 2 + a) * 24 + 6 * j + 2 + b] = kp[(3 * i + a) * 12 + 3 * j + b];
    }
    const kd = drillingFactor * E * t * Math.abs(area());
    // Relative drilling stiffness (no resistance to a rigid in-plane rotation) plus 1 % absolute part, which
    // keeps coplanar regions non-singular at a negligible spurious restraint.
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) kl[(6 * i + 5) * 24 + 6 * j + 5] = i === j ? 1.01 * kd : -kd / 3;
    return kl;
  };
  const elemNodes = (e) => quads.subarray ? quads.subarray(e * 4, e * 4 + 4) : quads.slice(e * 4, e * 4 + 4);
  const elemMatrix = (e) => {
    localK(e);
    // K_global block = Tᵀ K_local block T for each 3×3 block (translations and rotations transform alike).
    for (let bi = 0; bi < 8; bi++) for (let bj = 0; bj < 8; bj++) {
      for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) {
        tmp[r * 3 + c] = kl[(3 * bi + r) * 24 + 3 * bj] * T[c] + kl[(3 * bi + r) * 24 + 3 * bj + 1] * T[3 + c] + kl[(3 * bi + r) * 24 + 3 * bj + 2] * T[6 + c];
      }
      for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) {
        kg[(3 * bi + r) * 24 + 3 * bj + c] = T[r] * tmp[c] + T[3 + r] * tmp[3 + c] + T[6 + r] * tmp[6 + c];
      }
    }
    return kg;
  };

  const f = new Float64Array(6 * nn);
  for (const L of nodalLoads) {
    if (!(L.node >= 0 && L.node < nn)) throw new RangeError(`${fn}: nodalLoads node index out of range`);
    const d = 6 * L.node;
    f[d] += L.fx || 0; f[d + 1] += L.fy || 0; f[d + 2] += L.fz || 0; f[d + 3] += L.mx || 0; f[d + 4] += L.my || 0; f[d + 5] += L.mz || 0;
  }
  const rule = ELEM.Q4.rule;
  const sc = pw.sc;
  if (pressure || gravity || surfaceLoad) {
    for (let e = 0; e < ne; e++) {
      frame(e);
      const t = tOf(e);
      const p = pOf(e) || 0;
      let qx = p * T[6]; let qy = p * T[7]; let qz = p * T[8];
      if (gravity) { qx += gravity.rho * t * gravity.g[0]; qy += gravity.rho * t * gravity.g[1]; qz += gravity.rho * t * gravity.g[2]; }
      if (surfaceLoad) { qx += surfaceLoad[0]; qy += surfaceLoad[1]; qz += surfaceLoad[2]; }
      for (let g = 0; g < 4; g++) {
        geom('Q4', 4, rule.xi[g], rule.eta[g], xl, sc);
        const w = rule.w[g] * sc.det;
        for (let i = 0; i < 4; i++) {
          const d = 6 * quads[e * 4 + i];
          f[d] += sc.N[i] * qx * w; f[d + 1] += sc.N[i] * qy * w; f[d + 2] += sc.N[i] * qz * w;
        }
      }
    }
  }
  const fixed = buildFixed(nn, ['ux', 'uy', 'uz', 'rx', 'ry', 'rz'], fix, fn);
  if (fixed.every(Number.isNaN)) throw new Error('Shell model has no supports.');
  const sys = feSystem({ nNodes: nn, dpn: 6, nElem: ne, elemNodes, elemMatrix, fixed, solver });
  const u = sys.solve(f);
  const ku = sys.mulK(u);
  const reactions = new Float64Array(6 * nn);
  const rs = { fx: 0, fy: 0, fz: 0 };
  const as = { fx: 0, fy: 0, fz: 0 };
  let energy = 0;
  const key = ['fx', 'fy', 'fz'];
  for (let i = 0; i < 6 * nn; i++) {
    energy += 0.5 * u[i] * ku[i];
    const k = i % 6;
    if (sys.isFixed[i]) { reactions[i] = ku[i] - f[i]; if (k < 3) rs[key[k]] += reactions[i]; }
    if (k < 3) as[key[k]] += f[i];
  }
  const out = {};
  for (const name of ['nx', 'ny', 'nxy', 'mx', 'my', 'mxy', 'vmTop', 'vmBottom', 'vmMid']) out[name] = new Float64Array(ne);
  const acc = new Float64Array(nn);
  const cnt = new Uint16Array(nn);
  const ul = new Float64Array(24);
  const vm = (sx, sy, txy) => Math.sqrt(Math.max(0, sx * sx + sy * sy - sx * sy + 3 * txy * txy));
  let maxVm = 0;
  for (let e = 0; e < ne; e++) {
    frame(e);
    const t = tOf(e);
    for (let i = 0; i < 4; i++) {
      const d = 6 * quads[e * 4 + i];
      for (let h = 0; h < 2; h++) for (let r = 0; r < 3; r++) {
        ul[6 * i + 3 * h + r] = T[3 * r] * u[d + 3 * h] + T[3 * r + 1] * u[d + 3 * h + 1] + T[3 * r + 2] * u[d + 3 * h + 2];
      }
    }
    geom('Q4', 4, 0, 0, xl, sc);
    let ex = 0; let ey = 0; let gxy = 0; let kx = 0; let ky = 0; let kxy = 0;
    for (let i = 0; i < 4; i++) {
      ex += sc.dx[i] * ul[6 * i]; ey += sc.dy[i] * ul[6 * i + 1]; gxy += sc.dy[i] * ul[6 * i] + sc.dx[i] * ul[6 * i + 1];
      kx += sc.dx[i] * ul[6 * i + 4]; ky -= sc.dy[i] * ul[6 * i + 3]; kxy += sc.dy[i] * ul[6 * i + 4] - sc.dx[i] * ul[6 * i + 3];
    }
    const c = (E * t) / (1 - nu * nu);
    const Db = (c * t * t) / 12;
    const Nx = c * (ex + nu * ey); const Ny = c * (nu * ex + ey); const Nxy = 0.5 * c * (1 - nu) * gxy;
    const Mx = Db * (kx + nu * ky); const My = Db * (nu * kx + ky); const Mxy = 0.5 * Db * (1 - nu) * kxy;
    out.nx[e] = Nx; out.ny[e] = Ny; out.nxy[e] = Nxy; out.mx[e] = Mx; out.my[e] = My; out.mxy[e] = Mxy;
    const s6 = 6 / (t * t);
    out.vmMid[e] = vm(Nx / t, Ny / t, Nxy / t);
    out.vmTop[e] = vm(Nx / t + s6 * Mx, Ny / t + s6 * My, Nxy / t + s6 * Mxy);
    out.vmBottom[e] = vm(Nx / t - s6 * Mx, Ny / t - s6 * My, Nxy / t - s6 * Mxy);
    const m = Math.max(out.vmTop[e], out.vmBottom[e]);
    if (m > maxVm) maxVm = m;
    for (let i = 0; i < 4; i++) { acc[quads[e * 4 + i]] += m; cnt[quads[e * 4 + i]]++; }
  }
  const nodalVonMises = new Float64Array(nn);
  let maxDisp = 0;
  for (let n = 0; n < nn; n++) {
    nodalVonMises[n] = acc[n] / (cnt[n] || 1);
    maxDisp = Math.max(maxDisp, Math.hypot(u[6 * n], u[6 * n + 1], u[6 * n + 2]));
  }
  return { u, reactions, reactionSum: rs, appliedSum: as, elem: out, nodalVonMises, maxVonMises: maxVm, maxDisp, strainEnergy: energy };
}
