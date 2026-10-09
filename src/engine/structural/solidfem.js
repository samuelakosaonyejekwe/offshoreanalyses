/**
 * Three-dimensional linear elasticity with solid finite elements, and a voxel mesher that turns any closed
 * triangulated surface into an analysable hexahedral mesh.
 *
 * Elements
 *  - hex8: trilinear hexahedron with nine incompatible bending modes (Wilson et al. 1973; Taylor, Beresford &
 *    Wilson 1976). The enhanced gradients use the centre Jacobian scaled by det J0 / det J, so the element
 *    passes the patch test for any shape and is free of shear locking in bending. `hex: 'full'` gives the plain
 *    2×2×2 isoparametric element (locks in bending — for comparison only).
 *  - tet4: constant-strain tetrahedron. tet10: quadratic tetrahedron, 4-point integration.
 *  - hex20 cells of imported models are analysed on their corner nodes as hex8.
 *
 * Solvers: skyline LDLᵀ for small models, CSR + preconditioned CG for general meshes, and a matrix-free
 * element-by-element CG for uniform voxel meshes (one 24×24 element matrix for the whole model, so memory is
 * a few vectors: about 60 000 cells on a laptop and 15 000 on a phone).
 *
 * Volume model (CONTRACT §4): { kind: 'volume', positions (x, y, z per node) [m], cells: [{ type, conn }],
 *   sets?: { name: FaceSet } } with FaceSet = { quad4?: Uint32Array, tri3?: Uint32Array, tri6?: Uint32Array,
 *   nodes: Uint32Array } (faces listed by their node ids, outward orientation).
 *
 * Strain/stress component order: xx, yy, zz, xy, yz, zx (engineering shear strains).
 *
 * References: Zienkiewicz, Taylor & Zhu, "The Finite Element Method" (6th ed.) vol. 1 ch. 4–6, 10;
 * Hughes, "The Finite Element Method" (1987) ch. 3–4; Bathe, "Finite Element Procedures" (1996) ch. 5, 11.
 */
import { feSystem, invertSmall, gauss1d, gaussTriangle, subspaceModes } from './fem2d.js';

const HEX_XI = [-1, 1, 1, -1, -1, 1, 1, -1];
const HEX_ETA = [-1, -1, 1, 1, -1, -1, 1, 1];
const HEX_ZETA = [-1, -1, -1, -1, 1, 1, 1, 1];
/** Outward-oriented local faces. */
const HEX_FACES = [[0, 4, 7, 3], [1, 2, 6, 5], [0, 1, 5, 4], [3, 7, 6, 2], [0, 3, 2, 1], [4, 5, 6, 7]];
const TET_FACES = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
const TET10_FACES = [[0, 2, 1, 6, 5, 4], [0, 1, 3, 4, 8, 7], [1, 2, 3, 5, 9, 8], [2, 0, 3, 6, 7, 9]];
const TET10_EDGES = [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]];
const NEN = { tet4: 4, tet10: 10, hex8: 8 };
const TYPE_ID = { tet4: 0, tet10: 1, hex8: 2 };
const TYPE_NAME = ['tet4', 'tet10', 'hex8'];

const TET_A = 0.5854101966249685;
const TET_B = 0.1381966011250105;

function hexRule() {
  const g = gauss1d(2);
  const r = { xi: [], eta: [], zeta: [], w: [] };
  for (let k = 0; k < 2; k++) for (let j = 0; j < 2; j++) for (let i = 0; i < 2; i++) {
    r.xi.push(g.x[i]); r.eta.push(g.x[j]); r.zeta.push(g.x[k]); r.w.push(1);
  }
  return r;
}

/** Quadrature per element type; tetrahedron weights sum to 1/6 (4-point rule, exact to degree 2). */
const RULE = {
  tet4: { xi: [0.25], eta: [0.25], zeta: [0.25], w: [1 / 6] },
  tet10: {
    xi: [TET_B, TET_A, TET_B, TET_B], eta: [TET_B, TET_B, TET_A, TET_B], zeta: [TET_B, TET_B, TET_B, TET_A],
    w: [1 / 24, 1 / 24, 1 / 24, 1 / 24],
  },
  hex8: hexRule(),
};

function positive(v, name, fn) {
  if (!(typeof v === 'number' && v > 0 && Number.isFinite(v))) throw new RangeError(`${fn}: ${name} must be a positive number`);
}

/**
 * Shape functions and natural derivatives of the solid elements.
 * @param {'tet4'|'tet10'|'hex8'} type element type
 * @param {number} xi natural coordinate
 * @param {number} eta natural coordinate
 * @param {number} zeta natural coordinate
 * @param {Float64Array} N out: shape functions
 * @param {Float64Array} d1 out: dN/dxi
 * @param {Float64Array} d2 out: dN/deta
 * @param {Float64Array} d3 out: dN/dzeta
 */
export function shape3d(type, xi, eta, zeta, N, d1, d2, d3) {
  if (type === 'hex8') {
    for (let i = 0; i < 8; i++) {
      const a = 1 + HEX_XI[i] * xi;
      const b = 1 + HEX_ETA[i] * eta;
      const c = 1 + HEX_ZETA[i] * zeta;
      N[i] = 0.125 * a * b * c;
      d1[i] = 0.125 * HEX_XI[i] * b * c;
      d2[i] = 0.125 * HEX_ETA[i] * a * c;
      d3[i] = 0.125 * HEX_ZETA[i] * a * b;
    }
  } else if (type === 'tet4') {
    N[0] = 1 - xi - eta - zeta; N[1] = xi; N[2] = eta; N[3] = zeta;
    d1[0] = -1; d1[1] = 1; d1[2] = 0; d1[3] = 0;
    d2[0] = -1; d2[1] = 0; d2[2] = 1; d2[3] = 0;
    d3[0] = -1; d3[1] = 0; d3[2] = 0; d3[3] = 1;
  } else if (type === 'tet10') {
    const L = [1 - xi - eta - zeta, xi, eta, zeta];
    // dL/d(xi, eta, zeta)
    const dL = [[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let i = 0; i < 4; i++) {
      N[i] = L[i] * (2 * L[i] - 1);
      const f = 4 * L[i] - 1;
      d1[i] = f * dL[i][0]; d2[i] = f * dL[i][1]; d3[i] = f * dL[i][2];
    }
    for (let k = 0; k < 6; k++) {
      const a = TET10_EDGES[k][0];
      const b = TET10_EDGES[k][1];
      N[4 + k] = 4 * L[a] * L[b];
      d1[4 + k] = 4 * (dL[a][0] * L[b] + L[a] * dL[b][0]);
      d2[4 + k] = 4 * (dL[a][1] * L[b] + L[a] * dL[b][1]);
      d3[4 + k] = 4 * (dL[a][2] * L[b] + L[a] * dL[b][2]);
    }
  } else {
    throw new RangeError(`shape3d: unknown element type ${type}`);
  }
}

function makeScratch() {
  const n = 10;
  return {
    N: new Float64Array(n), d1: new Float64Array(n), d2: new Float64Array(n), d3: new Float64Array(n),
    dx: new Float64Array(n), dy: new Float64Array(n), dz: new Float64Array(n), ji: new Float64Array(9), det: 0,
    x: 0, y: 0, z: 0,
  };
}

/** Cartesian derivatives, Jacobian inverse (ji maps natural to Cartesian gradients) and determinant. */
function geom3(type, nen, xi, eta, zeta, xe, sc) {
  shape3d(type, xi, eta, zeta, sc.N, sc.d1, sc.d2, sc.d3);
  let a = 0; let b = 0; let c = 0; let d = 0; let e = 0; let f = 0; let g = 0; let h = 0; let i = 0;
  let x = 0; let y = 0; let z = 0;
  for (let k = 0; k < nen; k++) {
    const X = xe[3 * k]; const Y = xe[3 * k + 1]; const Z = xe[3 * k + 2];
    a += sc.d1[k] * X; b += sc.d1[k] * Y; c += sc.d1[k] * Z;
    d += sc.d2[k] * X; e += sc.d2[k] * Y; f += sc.d2[k] * Z;
    g += sc.d3[k] * X; h += sc.d3[k] * Y; i += sc.d3[k] * Z;
    x += sc.N[k] * X; y += sc.N[k] * Y; z += sc.N[k] * Z;
  }
  const det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
  if (!(det > 0)) throw new Error('Solid element with non-positive Jacobian (inverted or badly distorted element).');
  const ji = sc.ji;
  ji[0] = (e * i - f * h) / det; ji[1] = (c * h - b * i) / det; ji[2] = (b * f - c * e) / det;
  ji[3] = (f * g - d * i) / det; ji[4] = (a * i - c * g) / det; ji[5] = (c * d - a * f) / det;
  ji[6] = (d * h - e * g) / det; ji[7] = (b * g - a * h) / det; ji[8] = (a * e - b * d) / det;
  for (let k = 0; k < nen; k++) {
    sc.dx[k] = ji[0] * sc.d1[k] + ji[1] * sc.d2[k] + ji[2] * sc.d3[k];
    sc.dy[k] = ji[3] * sc.d1[k] + ji[4] * sc.d2[k] + ji[5] * sc.d3[k];
    sc.dz[k] = ji[6] * sc.d1[k] + ji[7] * sc.d2[k] + ji[8] * sc.d3[k];
  }
  sc.det = det; sc.x = x; sc.y = y; sc.z = z;
}

/** Strain–displacement rows for n "nodes" with gradients (gx, gy, gz); B is 6 × nc with column offset. */
function fillB3(n, gx, gy, gz, B, nc, off) {
  for (let i = 0; i < n; i++) {
    const c = off + 3 * i;
    B[c] = gx[i];
    B[nc + c + 1] = gy[i];
    B[2 * nc + c + 2] = gz[i];
    B[3 * nc + c] = gy[i]; B[3 * nc + c + 1] = gx[i];
    B[4 * nc + c + 1] = gz[i]; B[4 * nc + c + 2] = gy[i];
    B[5 * nc + c] = gz[i]; B[5 * nc + c + 2] = gx[i];
  }
}

/**
 * Isotropic elasticity matrix (6×6 row-major) for strains [εxx, εyy, εzz, γxy, γyz, γzx].
 * @param {number} E Young's modulus [Pa]
 * @param {number} nu Poisson's ratio [-]
 * @returns {Float64Array} D [Pa]
 */
export function solidD(E, nu) {
  positive(E, 'E', 'solidD');
  if (!(nu > -1 && nu < 0.5)) throw new RangeError('solidD: nu must lie in (-1, 0.5)');
  const D = new Float64Array(36);
  const G = E / (2 * (1 + nu));
  const lam = (E * nu) / ((1 + nu) * (1 - 2 * nu));
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) D[i * 6 + j] = lam + (i === j ? 2 * G : 0);
  D[21] = G; D[28] = G; D[35] = G;
  return D;
}

/** out (6 × nc) = D · B (6 × nc). */
function mulDB(D, B, nc, out) {
  for (let r = 0; r < 6; r++) {
    const d0 = D[r * 6]; const d1 = D[r * 6 + 1]; const d2 = D[r * 6 + 2]; const d3 = D[r * 6 + 3]; const d4 = D[r * 6 + 4]; const d5 = D[r * 6 + 5];
    for (let q = 0; q < nc; q++) {
      out[r * nc + q] = d0 * B[q] + d1 * B[nc + q] + d2 * B[2 * nc + q] + d3 * B[3 * nc + q] + d4 * B[4 * nc + q] + d5 * B[5 * nc + q];
    }
  }
}

/** K (na × nb) += f · Aᵀ (6 × na) · C (6 × nb). */
function addAtC(A, na, C, nb, f, K) {
  for (let p = 0; p < na; p++) {
    const a0 = f * A[p]; const a1 = f * A[na + p]; const a2 = f * A[2 * na + p]; const a3 = f * A[3 * na + p]; const a4 = f * A[4 * na + p]; const a5 = f * A[5 * na + p];
    const row = p * nb;
    for (let q = 0; q < nb; q++) {
      K[row + q] += a0 * C[q] + a1 * C[nb + q] + a2 * C[2 * nb + q] + a3 * C[3 * nb + q] + a4 * C[4 * nb + q] + a5 * C[5 * nb + q];
    }
  }
}

function makeWork() {
  return {
    sc: makeScratch(), B: new Float64Array(6 * 30), DB: new Float64Array(6 * 30), ke: { tet4: new Float64Array(144), tet10: new Float64Array(900), hex8: new Float64Array(576) },
    G: new Float64Array(54), DG: new Float64Array(54), kua: new Float64Array(24 * 9), kaa: new Float64Array(81), rec: new Float64Array(9 * 24),
    mg: [new Float64Array(3), new Float64Array(3), new Float64Array(3)], j0: new Float64Array(9),
    S: new Float64Array(8 * 6 * 30),
  };
}

/** Enhanced-mode gradients at a Gauss point: fills wk.G (6 × 9). Requires wk.j0 and detJ0 from the centre. */
function hexEnhanced(xi, eta, zeta, det0, det, wk) {
  const s = det0 / det;
  const q = [-2 * xi * s, -2 * eta * s, -2 * zeta * s];
  for (let m = 0; m < 3; m++) {
    wk.mg[0][m] = wk.j0[m] * q[m];
    wk.mg[1][m] = wk.j0[3 + m] * q[m];
    wk.mg[2][m] = wk.j0[6 + m] * q[m];
  }
  fillB3(3, wk.mg[0], wk.mg[1], wk.mg[2], wk.G, 9, 0);
}

/**
 * Element stiffness. With `stressOp` the routine also fills wk.S, the operator giving the stresses at the
 * integration points from the element displacement vector (ng × 6 rows, m columns).
 */
function elementK(type, xe, D, hexMode, wk, stressOp) {
  const nen = NEN[type];
  const m = 3 * nen;
  const { sc, B, DB } = wk;
  const ke = wk.ke[type];
  ke.fill(0);
  const rule = RULE[type];
  const ng = rule.w.length;
  const incompatible = type === 'hex8' && hexMode === 'incompatible';
  let det0 = 0;
  if (incompatible) {
    geom3(type, nen, 0, 0, 0, xe, sc);
    det0 = sc.det;
    wk.j0.set(sc.ji);
    wk.kua.fill(0); wk.kaa.fill(0);
  }
  for (let g = 0; g < ng; g++) {
    geom3(type, nen, rule.xi[g], rule.eta[g], rule.zeta[g], xe, sc);
    fillB3(nen, sc.dx, sc.dy, sc.dz, B, m, 0);
    mulDB(D, B, m, DB);
    const f = rule.w[g] * sc.det;
    addAtC(B, m, DB, m, f, ke);
    if (stressOp) for (let k = 0; k < 6 * m; k++) wk.S[g * 6 * m + k] = DB[k];
    if (incompatible) {
      hexEnhanced(rule.xi[g], rule.eta[g], rule.zeta[g], det0, sc.det, wk);
      mulDB(D, wk.G, 9, wk.DG);
      addAtC(B, m, wk.DG, 9, f, wk.kua);
      addAtC(wk.G, 9, wk.DG, 9, f, wk.kaa);
    }
  }
  if (incompatible) {
    invertSmall(wk.kaa, 9);
    for (let a = 0; a < 9; a++) for (let q = 0; q < 24; q++) {
      let s = 0;
      for (let b = 0; b < 9; b++) s += wk.kaa[a * 9 + b] * wk.kua[q * 9 + b];
      wk.rec[a * 24 + q] = -s;
    }
    for (let p = 0; p < 24; p++) for (let q = 0; q < 24; q++) {
      let s = 0;
      for (let a = 0; a < 9; a++) s += wk.kua[p * 9 + a] * wk.rec[a * 24 + q];
      ke[p * 24 + q] += s;
    }
    if (stressOp) {
      for (let g = 0; g < ng; g++) {
        geom3(type, nen, rule.xi[g], rule.eta[g], rule.zeta[g], xe, sc);
        hexEnhanced(rule.xi[g], rule.eta[g], rule.zeta[g], det0, sc.det, wk);
        mulDB(D, wk.G, 9, wk.DG);
        for (let r = 0; r < 6; r++) for (let q = 0; q < 24; q++) {
          let s = 0;
          for (let a = 0; a < 9; a++) s += wk.DG[r * 9 + a] * wk.rec[a * 24 + q];
          wk.S[(g * 6 + r) * 24 + q] += s;
        }
      }
    }
  }
  return ke;
}

/** Extrapolation weights from integration points to element nodes. */
const EXTRAP = (() => {
  const hex = [];
  const r = RULE.hex8;
  for (let i = 0; i < 8; i++) {
    const row = [];
    for (let g = 0; g < 8; g++) {
      row.push(0.125 * (1 + 3 * HEX_XI[i] * r.xi[g]) * (1 + 3 * HEX_ETA[i] * r.eta[g]) * (1 + 3 * HEX_ZETA[i] * r.zeta[g]));
    }
    hex.push(row);
  }
  // Tetrahedron: linear field through the four points; rule point g has barycentric L_g = a, the others b.
  const corner = [];
  for (let i = 0; i < 4; i++) {
    const row = [];
    for (let g = 0; g < 4; g++) row.push(((i === g ? 1 : 0) - TET_B) / (TET_A - TET_B));
    corner.push(row);
  }
  const tet10 = corner.map((r0) => r0.slice());
  for (const [a, b] of TET10_EDGES) tet10.push(corner[a].map((v, g) => 0.5 * (v + corner[b][g])));
  return { hex8: hex, tet10, tet4: [[1], [1], [1], [1]] };
})();

/* ------------------------------------------------------------------------------------------------------------
 * Model preparation, topology and sets
 * ---------------------------------------------------------------------------------------------------------- */

function tetVolume6(p, a, b, c, d) {
  const ax = p[3 * b] - p[3 * a]; const ay = p[3 * b + 1] - p[3 * a + 1]; const az = p[3 * b + 2] - p[3 * a + 2];
  const bx = p[3 * c] - p[3 * a]; const by = p[3 * c + 1] - p[3 * a + 1]; const bz = p[3 * c + 2] - p[3 * a + 2];
  const cx = p[3 * d] - p[3 * a]; const cy = p[3 * d + 1] - p[3 * a + 1]; const cz = p[3 * d + 2] - p[3 * a + 2];
  return ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
}

/**
 * Flatten a volume model into element arrays, correcting inverted elements (negative volume) by renumbering.
 * Returns { nn, ne, pos, etype, eptr, econn, uniform, counts }.
 */
function prepare(model, fn) {
  if (!model || !model.positions || !Array.isArray(model.cells)) throw new RangeError(`${fn}: model must be a volume model with positions and cells`);
  const pos = model.positions;
  if (pos.length % 3 !== 0 || pos.length === 0) throw new RangeError(`${fn}: model.positions must hold 3 coordinates per node`);
  const nn = pos.length / 3;
  let ne = 0;
  let total = 0;
  for (const c of model.cells) {
    const src = c.type === 'hex20' ? 20 : NEN[c.type];
    if (!src) throw new RangeError(`${fn}: cell type ${c.type} is not supported (use tet4, tet10, hex8 or hex20)`);
    if (c.conn.length % src !== 0) throw new RangeError(`${fn}: connectivity length of a ${c.type} block is not a multiple of ${src}`);
    const k = c.conn.length / src;
    ne += k;
    total += k * (c.type === 'hex20' ? 8 : src);
  }
  if (ne === 0) throw new RangeError(`${fn}: model has no cells`);
  const etype = new Uint8Array(ne);
  const eptr = new Int32Array(ne + 1);
  const econn = new Int32Array(total);
  const counts = { tet4: 0, tet10: 0, hex8: 0 };
  let e = 0;
  let p = 0;
  for (const c of model.cells) {
    const src = c.type === 'hex20' ? 20 : NEN[c.type];
    const type = c.type === 'hex20' ? 'hex8' : c.type;
    const nen = NEN[type];
    const k = c.conn.length / src;
    for (let i = 0; i < k; i++) {
      etype[e] = TYPE_ID[type];
      eptr[e] = p;
      for (let a = 0; a < nen; a++) {
        const n = c.conn[i * src + a];
        if (!(n < nn)) throw new RangeError(`${fn}: connectivity refers to node ${n} beyond positions`);
        econn[p + a] = n;
      }
      if (type === 'hex8') {
        // Sign from the corner tetrahedron at node 0 (edges to nodes 1, 3, 4).
        if (tetVolume6(pos, econn[p], econn[p + 1], econn[p + 3], econn[p + 4]) < 0) {
          for (let a = 0; a < 4; a++) { const t = econn[p + a]; econn[p + a] = econn[p + 4 + a]; econn[p + 4 + a] = t; }
        }
      } else if (tetVolume6(pos, econn[p], econn[p + 1], econn[p + 2], econn[p + 3]) < 0) {
        let t = econn[p + 1]; econn[p + 1] = econn[p + 2]; econn[p + 2] = t;
        if (type === 'tet10') {
          t = econn[p + 4]; econn[p + 4] = econn[p + 6]; econn[p + 6] = t;
          t = econn[p + 8]; econn[p + 8] = econn[p + 9]; econn[p + 9] = t;
        }
      }
      counts[type]++;
      p += nen;
      e++;
    }
  }
  eptr[ne] = p;
  const uniform = !!(model.grid && counts.hex8 === ne && model.grid.h > 0);
  return { nn, ne, pos, etype, eptr, econn, uniform, counts };
}

function faceKey(ids) {
  const s = Array.from(ids).sort((a, b) => a - b);
  return s.join(',');
}

/**
 * Boundary faces of a volume model (faces that belong to exactly one cell), outward oriented.
 * @param {object} model volume model (tet4, tet10, hex8; hex20 by its corners)
 * @returns {{quad4: Uint32Array, tri3: Uint32Array, tri6: Uint32Array, nodes: Uint32Array}} face set
 */
export function boundaryFaces(model) {
  const P = prepare(model, 'boundaryFaces');
  const seen = new Map();
  for (let e = 0; e < P.ne; e++) {
    const t = P.etype[e];
    const c = P.eptr[e];
    const faces = t === 2 ? HEX_FACES : TET_FACES;
    for (let k = 0; k < faces.length; k++) {
      const key = faceKey(faces[k].map((a) => P.econn[c + a]));
      const prev = seen.get(key);
      if (prev === undefined) seen.set(key, e * 8 + k);
      else seen.set(key, -1);
    }
  }
  const quad4 = []; const tri3 = []; const tri6 = [];
  for (const v of seen.values()) {
    if (v < 0) continue;
    const e = v >> 3;
    const k = v & 7;
    const t = P.etype[e];
    const c = P.eptr[e];
    if (t === 2) for (const a of HEX_FACES[k]) quad4.push(P.econn[c + a]);
    else if (t === 0) for (const a of TET_FACES[k]) tri3.push(P.econn[c + a]);
    else for (const a of TET10_FACES[k]) tri6.push(P.econn[c + a]);
  }
  return finishFaceSet({ quad4: Uint32Array.from(quad4), tri3: Uint32Array.from(tri3), tri6: Uint32Array.from(tri6) });
}

function finishFaceSet(set) {
  const ids = new Set();
  for (const k of ['quad4', 'tri3', 'tri6']) {
    if (!set[k]) set[k] = new Uint32Array(0);
    for (let i = 0; i < set[k].length; i++) ids.add(set[k][i]);
  }
  set.nodes = Uint32Array.from(ids).sort();
  return set;
}

function faceGeometry(pos, ids) {
  const nc = ids.length === 4 ? 4 : 3;
  let cx = 0; let cy = 0; let cz = 0;
  for (let i = 0; i < nc; i++) { cx += pos[3 * ids[i]]; cy += pos[3 * ids[i] + 1]; cz += pos[3 * ids[i] + 2]; }
  cx /= nc; cy /= nc; cz /= nc;
  // Area vector by the cross product of diagonals (quad) or edges (triangle).
  const a = ids[0]; const b = ids[1]; const c = ids[2]; const d = nc === 4 ? ids[3] : ids[0];
  const ux = pos[3 * c] - pos[3 * a]; const uy = pos[3 * c + 1] - pos[3 * a + 1]; const uz = pos[3 * c + 2] - pos[3 * a + 2];
  const vx = nc === 4 ? pos[3 * d] - pos[3 * b] : pos[3 * b] - pos[3 * a];
  const vy = nc === 4 ? pos[3 * d + 1] - pos[3 * b + 1] : pos[3 * b + 1] - pos[3 * a + 1];
  const vz = nc === 4 ? pos[3 * d + 2] - pos[3 * b + 2] : pos[3 * b + 2] - pos[3 * a + 2];
  let nx; let ny; let nz;
  if (nc === 4) { nx = 0.5 * (uy * vz - uz * vy); ny = 0.5 * (uz * vx - ux * vz); nz = 0.5 * (ux * vy - uy * vx); } else { nx = 0.5 * (vy * uz - vz * uy); ny = 0.5 * (vz * ux - vx * uz); nz = 0.5 * (vx * uy - vy * ux); }
  const area = Math.hypot(nx, ny, nz) || 1e-300;
  return { centroid: [cx, cy, cz], normal: [nx / area, ny / area, nz / area], area };
}

/**
 * Select faces of a face set (default: the whole boundary) by a predicate on face centroid, unit outward
 * normal and area — the way to build load and support sets on imported meshes.
 * @param {object} model volume model
 * @param {(f: {centroid: number[], normal: number[], area: number}) => boolean} pred selection predicate
 *   (centroid [m], normal [-], area [m²])
 * @param {object} [from] face set to filter (default `boundaryFaces(model)`)
 * @returns {{quad4: Uint32Array, tri3: Uint32Array, tri6: Uint32Array, nodes: Uint32Array}} face set
 */
export function selectFaces(model, pred, from) {
  const src = from || boundaryFaces(model);
  const pos = model.positions;
  const out = { quad4: [], tri3: [], tri6: [] };
  for (const [k, n] of [['quad4', 4], ['tri3', 3], ['tri6', 6]]) {
    const arr = src[k] || [];
    for (let i = 0; i < arr.length; i += n) {
      const ids = arr.subarray(i, i + n);
      if (pred(faceGeometry(pos, ids))) for (let a = 0; a < n; a++) out[k].push(ids[a]);
    }
  }
  return finishFaceSet({ quad4: Uint32Array.from(out.quad4), tri3: Uint32Array.from(out.tri3), tri6: Uint32Array.from(out.tri6) });
}

/**
 * Total area and area vector of a face set.
 * @param {object} model volume model
 * @param {object} set face set
 * @returns {{area: number, areaVector: number[]}} area [m²] and ∑ n·dA [m²]
 */
export function faceSetArea(model, set) {
  const pos = model.positions;
  let area = 0;
  const v = [0, 0, 0];
  for (const [k, n] of [['quad4', 4], ['tri3', 3], ['tri6', 6]]) {
    const arr = set[k] || [];
    for (let i = 0; i < arr.length; i += n) {
      const g = faceGeometry(pos, arr.subarray(i, i + n));
      area += g.area;
      for (let c = 0; c < 3; c++) v[c] += g.normal[c] * g.area;
    }
  }
  return { area, areaVector: v };
}

/**
 * Structured mesh of a rectangular block, for verification and simple solid parts.
 * @param {object} o
 * @param {number} o.lx size in x [m]
 * @param {number} o.ly size in y [m]
 * @param {number} o.lz size in z [m]
 * @param {number} o.nx divisions in x
 * @param {number} o.ny divisions in y
 * @param {number} o.nz divisions in z
 * @param {'hex8'|'tet4'|'tet10'} [o.type] element type (tetrahedra: six per block cell, Kuhn subdivision)
 * @param {(x: number, y: number, z: number) => number[]} [o.distort] optional node relocation (patch tests)
 * @returns {object} volume model with face sets xmin, xmax, ymin, ymax, zmin, zmax (named by the undistorted
 *   block faces)
 */
export function boxMesh({ lx, ly, lz, nx, ny, nz, type = 'hex8', distort = null }) {
  const fn = 'boxMesh';
  positive(lx, 'lx', fn); positive(ly, 'ly', fn); positive(lz, 'lz', fn);
  for (const [v, name] of [[nx, 'nx'], [ny, 'ny'], [nz, 'nz']]) if (!(Number.isInteger(v) && v >= 1)) throw new RangeError(`${fn}: ${name} must be a positive integer`);
  if (!NEN[type]) throw new RangeError(`${fn}: type must be hex8, tet4 or tet10`);
  const gx = nx + 1; const gy = ny + 1; const gz = nz + 1;
  const ref = [];
  const coords = [];
  const id = (i, j, k) => (k * gy + j) * gx + i;
  for (let k = 0; k < gz; k++) for (let j = 0; j < gy; j++) for (let i = 0; i < gx; i++) {
    const x = (lx * i) / nx; const y = (ly * j) / ny; const z = (lz * k) / nz;
    ref.push(x, y, z);
    const p = distort ? distort(x, y, z) : [x, y, z];
    coords.push(p[0], p[1], p[2]);
  }
  const conn = [];
  const kuhn = [[0, 1, 2, 6], [0, 5, 1, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 4, 5, 6], [0, 7, 4, 6]];
  const mids = new Map();
  const nGrid = gx * gy * gz;
  const mid = (a, b) => {
    const key = Math.min(a, b) * nGrid + Math.max(a, b);
    let m = mids.get(key);
    if (m === undefined) {
      m = coords.length / 3;
      for (let c = 0; c < 3; c++) { coords.push(0.5 * (coords[3 * a + c] + coords[3 * b + c])); ref.push(0.5 * (ref[3 * a + c] + ref[3 * b + c])); }
      mids.set(key, m);
    }
    return m;
  };
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const h = [id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k), id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)];
    if (type === 'hex8') { conn.push(...h); continue; }
    for (const t of kuhn) {
      const c = t.map((a) => h[a]);
      if (type === 'tet4') conn.push(...c);
      else conn.push(...c, ...TET10_EDGES.map(([a, b]) => mid(c[a], c[b])));
    }
  }
  const model = { kind: 'volume', name: 'box', positions: Float64Array.from(coords), cells: [{ type, conn: Uint32Array.from(conn) }] };
  // Sets are identified on the undistorted block so that names keep their meaning after distortion.
  const refModel = { positions: Float64Array.from(ref), cells: model.cells };
  const all = boundaryFaces(refModel);
  const tol = 1e-9 * Math.max(lx, ly, lz);
  const pick = (axis, val) => selectFaces(refModel, (f) => Math.abs(f.centroid[axis] - val) < tol, all);
  model.sets = { xmin: pick(0, 0), xmax: pick(0, lx), ymin: pick(1, 0), ymax: pick(1, ly), zmin: pick(2, 0), zmax: pick(2, lz) };
  return model;
}

/**
 * Mesh statistics of a volume model.
 * @param {object} model volume model
 * @returns {{nNodes: number, nCells: number, nDof: number, cellCounts: {tet4: number, tet10: number, hex8: number},
 *   volume: number, bbox: {min: number[], max: number[]}, edgeMin: number, edgeMax: number, maxAspect: number}}
 *   volume [m³], bounding box [m], shortest/longest element edge [m], worst edge-length ratio within one element
 */
export function meshStats(model) {
  const P = prepare(model, 'meshStats');
  const wk = makeWork();
  const xe = new Float64Array(30);
  let volume = 0; let edgeMin = Infinity; let edgeMax = 0; let maxAspect = 1;
  const min = [Infinity, Infinity, Infinity];
  const max = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < P.nn; i++) for (let c = 0; c < 3; c++) { const v = P.pos[3 * i + c]; if (v < min[c]) min[c] = v; if (v > max[c]) max[c] = v; }
  const hexEdges = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]];
  for (let e = 0; e < P.ne; e++) {
    const type = TYPE_NAME[P.etype[e]];
    const nen = NEN[type];
    const c = P.eptr[e];
    for (let a = 0; a < nen; a++) for (let k = 0; k < 3; k++) xe[3 * a + k] = P.pos[3 * P.econn[c + a] + k];
    const rule = RULE[type];
    for (let g = 0; g < rule.w.length; g++) {
      geom3(type, nen, rule.xi[g], rule.eta[g], rule.zeta[g], xe, wk.sc);
      volume += rule.w[g] * wk.sc.det;
    }
    let lo = Infinity; let hi = 0;
    for (const [a, b] of type === 'hex8' ? hexEdges : TET10_EDGES) {
      const l = Math.hypot(xe[3 * a] - xe[3 * b], xe[3 * a + 1] - xe[3 * b + 1], xe[3 * a + 2] - xe[3 * b + 2]);
      if (l < lo) lo = l;
      if (l > hi) hi = l;
    }
    if (lo < edgeMin) edgeMin = lo;
    if (hi > edgeMax) edgeMax = hi;
    if (hi / lo > maxAspect) maxAspect = hi / lo;
  }
  return { nNodes: P.nn, nCells: P.ne, nDof: 3 * P.nn, cellCounts: P.counts, volume, bbox: { min, max }, edgeMin, edgeMax, maxAspect };
}

/* ------------------------------------------------------------------------------------------------------------
 * Voxel mesher
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Parity pass along one axis. For every triangle the rays through the cell centres that pierce it are found
 * in the projection plane with a half-open edge rule: a point exactly on a projected edge is attributed to the
 * triangle lying on the positive side of the canonically ordered edge, which is equivalent to an infinitesimal
 * shift of the ray and therefore counts shared edges and vertices exactly once (Simulation-of-Simplicity style
 * tie handling). Each crossing toggles a difference array; a prefix XOR along the ray gives inside/outside.
 */
function parityPass(pos, idx, d, origin, h, dims, votes) {
  const u = (d + 1) % 3;
  const v = (d + 2) % 3;
  const nd = dims[d]; const nu = dims[u]; const nv = dims[v];
  const diff = new Uint8Array(nu * nv * (nd + 1));
  const ou = origin[u] + 0.5 * h;
  const ov = origin[v] + 0.5 * h;
  const od = origin[d] + 0.5 * h;
  const nt = idx.length / 3;
  const P = [0, 0, 0];
  const edge = (a, b, c, pu, pv) => {
    // Returns +1 if (pu, pv) is attributed to the side of vertex c, -1 otherwise.
    let au = pos[3 * a + u]; let av = pos[3 * a + v]; let bu = pos[3 * b + u]; let bv = pos[3 * b + v];
    if (au > bu || (au === bu && (av > bv || (av === bv && a > b)))) { const tu = au; au = bu; bu = tu; const tv = av; av = bv; bv = tv; }
    const cp = (bu - au) * (pv - av) - (bv - av) * (pu - au);
    const cc = (bu - au) * (pos[3 * c + v] - av) - (bv - av) * (pos[3 * c + u] - au);
    if (cc === 0) return -1;
    if (cp === 0) return cc > 0 ? 1 : -1;
    return (cp > 0) === (cc > 0) ? 1 : -1;
  };
  for (let t = 0; t < nt; t++) {
    P[0] = idx[3 * t]; P[1] = idx[3 * t + 1]; P[2] = idx[3 * t + 2];
    const a = P[0]; const b = P[1]; const c = P[2];
    const au = pos[3 * a + u]; const av = pos[3 * a + v]; const ad = pos[3 * a + d];
    const bu = pos[3 * b + u]; const bv = pos[3 * b + v]; const bd = pos[3 * b + d];
    const cu = pos[3 * c + u]; const cv = pos[3 * c + v]; const cd = pos[3 * c + d];
    // Normal components: n_d is twice the signed projected area.
    const nD = (bu - au) * (cv - av) - (bv - av) * (cu - au);
    if (nD === 0) continue;
    const nU = (bv - av) * (cd - ad) - (bd - ad) * (cv - av);
    const nV = (bd - ad) * (cu - au) - (bu - au) * (cd - ad);
    const i0 = Math.max(0, Math.ceil((Math.min(au, bu, cu) - ou) / h));
    const i1 = Math.min(nu - 1, Math.floor((Math.max(au, bu, cu) - ou) / h));
    const j0 = Math.max(0, Math.ceil((Math.min(av, bv, cv) - ov) / h));
    const j1 = Math.min(nv - 1, Math.floor((Math.max(av, bv, cv) - ov) / h));
    for (let j = j0; j <= j1; j++) {
      const pv = ov + j * h;
      for (let i = i0; i <= i1; i++) {
        const pu = ou + i * h;
        if (edge(a, b, c, pu, pv) < 0 || edge(b, c, a, pu, pv) < 0 || edge(c, a, b, pu, pv) < 0) continue;
        const w = ad - (nU * (pu - au) + nV * (pv - av)) / nD;
        // First cell whose centre lies beyond the crossing.
        let k = Math.floor((w - od) / h) + 1;
        if (k < 0) k = 0;
        if (k > nd) k = nd;
        diff[(j * nu + i) * (nd + 1) + k] ^= 1;
      }
    }
  }
  const stride = [1, dims[0], dims[0] * dims[1]];
  for (let j = 0; j < nv; j++) {
    for (let i = 0; i < nu; i++) {
      const base = (j * nu + i) * (nd + 1);
      const cell0 = i * stride[u] + j * stride[v];
      let s = 0;
      for (let k = 0; k < nd; k++) {
        s ^= diff[base + k];
        if (s) votes[cell0 + k * stride[d]]++;
      }
    }
  }
}

/**
 * Enclosed volume and area of a triangulated surface (divergence theorem).
 * @param {object} surface surface model { positions, indices }
 * @returns {{volume: number, area: number, bbox: {min: number[], max: number[]}}} signed volume [m³] (positive
 *   for outward-oriented triangles), area [m²], bounding box [m]
 */
export function surfaceVolume(surface) {
  if (!surface || !surface.positions || !surface.indices) throw new RangeError('surfaceVolume: surface must have positions and indices');
  const p = surface.positions;
  const ix = surface.indices;
  if (p.length % 3 !== 0 || ix.length % 3 !== 0 || ix.length === 0) throw new RangeError('surfaceVolume: surface positions/indices must be triples');
  const min = [Infinity, Infinity, Infinity];
  const max = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < p.length; i += 3) for (let c = 0; c < 3; c++) { const v = p[i + c]; if (v < min[c]) min[c] = v; if (v > max[c]) max[c] = v; }
  // Shift to the box centre to keep the signed-volume sum well conditioned.
  const o = [0.5 * (min[0] + max[0]), 0.5 * (min[1] + max[1]), 0.5 * (min[2] + max[2])];
  let vol = 0; let area = 0;
  for (let t = 0; t < ix.length; t += 3) {
    const a = 3 * ix[t]; const b = 3 * ix[t + 1]; const c = 3 * ix[t + 2];
    const ax = p[a] - o[0]; const ay = p[a + 1] - o[1]; const az = p[a + 2] - o[2];
    const bx = p[b] - o[0]; const by = p[b + 1] - o[1]; const bz = p[b + 2] - o[2];
    const cx = p[c] - o[0]; const cy = p[c + 1] - o[1]; const cz = p[c + 2] - o[2];
    vol += ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
    area += 0.5 * Math.hypot((by - ay) * (cz - az) - (bz - az) * (cy - ay), (bz - az) * (cx - ax) - (bx - ax) * (cz - az), (bx - ax) * (cy - ay) - (by - ay) * (cx - ax));
  }
  return { volume: vol / 6, area, bbox: { min, max } };
}

/**
 * Voxel mesher: builds a hex8 mesh of the solid enclosed by a closed triangulated surface by classifying the
 * centres of a regular grid of spacing h as inside or outside. Classification is by ray-crossing parity along
 * each of the three axes with a majority vote, which tolerates small gaps and inconsistent triangle orientation;
 * exact ray/edge and ray/vertex coincidences are resolved consistently (see `parityPass`). Cells not connected
 * face-to-face to the largest body are removed so the model has no loose parts.
 *
 * The staircase boundary converges to the true geometry at first order in h: displacement-type results converge
 * at about O(h), surface stresses on inclined boundaries are indicative only.
 *
 * @param {object} o
 * @param {{positions: ArrayLike<number>, indices: ArrayLike<number>}} o.surface closed triangulated surface [m]
 * @param {number} [o.h] cell size [m]; give either h or targetCells
 * @param {number} [o.targetCells] approximate number of cells wanted (h = ∛(enclosed volume / targetCells))
 * @param {Array<{name: string, point: number[], normal: number[], tol?: number}>} [o.planes] user planes: boundary
 *   faces within tol (default 0.75 h) of the plane and facing along its normal form the set `name`
 * @param {boolean} [o.keepLargest] keep only the largest face-connected body (default true)
 * @param {number} [o.maxGridCells] guard on the background grid size (default 3e7)
 * @param {string} [o.name] model name
 * @returns {object} volume model { kind: 'volume', name, positions: Float64Array, cells: [{ type: 'hex8', conn }],
 *   sets, grid: { h, origin, nx, ny, nz }, stats } where sets holds face sets '+x', '-x', '+y', '-y', '+z', '-z'
 *   (all boundary faces by outward normal), 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax' (faces on the
 *   extreme layers), 'boundary' and the user-plane sets; stats = { nCells, nNodes, nDof, h, volume [m³],
 *   surfaceVolume [m³], volumeError [-], boundaryArea [m²], voteDisagreement [-] (fraction of candidate cells
 *   where the three ray directions disagreed; > ~0.01 suggests a leaky surface), removedCells, bbox }
 */
export function voxelMesh(o) {
  const fn = 'voxelMesh';
  const { surface, planes = [], keepLargest = true, maxGridCells = 3e7, name = 'voxel' } = o;
  const sv = surfaceVolume(surface);
  const pos = surface.positions;
  const idx = surface.indices;
  const nV = pos.length / 3;
  for (let i = 0; i < idx.length; i++) if (!(idx[i] < nV)) throw new RangeError(`${fn}: surface.indices refers to a vertex beyond positions`);
  for (let i = 0; i < pos.length; i++) if (!Number.isFinite(pos[i])) throw new RangeError(`${fn}: surface.positions contains non-finite values`);
  const ext = [sv.bbox.max[0] - sv.bbox.min[0], sv.bbox.max[1] - sv.bbox.min[1], sv.bbox.max[2] - sv.bbox.min[2]];
  if (!(ext[0] > 0 && ext[1] > 0 && ext[2] > 0)) throw new RangeError(`${fn}: surface is flat or empty (zero bounding-box extent)`);
  let h = o.h;
  if (h === undefined) {
    positive(o.targetCells, 'targetCells (or h)', fn);
    const vol = Math.abs(sv.volume) > 1e-12 * ext[0] * ext[1] * ext[2] ? Math.abs(sv.volume) : ext[0] * ext[1] * ext[2];
    h = Math.cbrt(vol / o.targetCells);
  }
  positive(h, 'h', fn);
  const dims = [0, 0, 0];
  const origin = [0, 0, 0];
  // Single-precision vertex data carries a relative error of ~6e-8 of the coordinate magnitude; without this
  // allowance a body that is an exact multiple of h would gain a spurious layer and shift by half a cell.
  let maxAbs = 0;
  for (let c = 0; c < 3; c++) maxAbs = Math.max(maxAbs, Math.abs(sv.bbox.min[c]), Math.abs(sv.bbox.max[c]));
  const slack = 1e-9 + ((pos instanceof Float32Array ? 2.4e-7 : 1e-13) * maxAbs) / h;
  for (let c = 0; c < 3; c++) {
    dims[c] = Math.max(1, Math.ceil(ext[c] / h - slack));
    origin[c] = 0.5 * (sv.bbox.min[c] + sv.bbox.max[c]) - 0.5 * dims[c] * h;
  }
  const nCellGrid = dims[0] * dims[1] * dims[2];
  if (nCellGrid > maxGridCells) throw new RangeError(`${fn}: h = ${h} gives a ${dims[0]}×${dims[1]}×${dims[2]} background grid (> ${maxGridCells} cells); increase h`);
  const votes = new Uint8Array(nCellGrid);
  for (let d = 0; d < 3; d++) parityPass(pos, idx, d, origin, h, dims, votes);
  let candidates = 0; let disagree = 0;
  for (let i = 0; i < nCellGrid; i++) {
    if (votes[i] > 0) { candidates++; if (votes[i] < 3) disagree++; }
    votes[i] = votes[i] >= 2 ? 1 : 0;
  }
  const [nx, ny, nz] = dims;
  const sx = 1; const sy = nx; const sz = nx * ny;
  let removed = 0;
  if (keepLargest) {
    // Face-connected components by flood fill; label 1 = unvisited inside, labels >= 2 = components.
    const stack = new Int32Array(nCellGrid);
    let best = 0; let bestLabel = 0; let label = 1;
    for (let s = 0; s < nCellGrid && label < 254; s++) {
      if (votes[s] !== 1) continue;
      label++;
      let top = 0; let size = 0;
      stack[top++] = s;
      votes[s] = label;
      while (top > 0) {
        const c = stack[--top];
        size++;
        const i = c % nx; const j = ((c / nx) | 0) % ny; const k = (c / sz) | 0;
        if (i > 0 && votes[c - sx] === 1) { votes[c - sx] = label; stack[top++] = c - sx; }
        if (i < nx - 1 && votes[c + sx] === 1) { votes[c + sx] = label; stack[top++] = c + sx; }
        if (j > 0 && votes[c - sy] === 1) { votes[c - sy] = label; stack[top++] = c - sy; }
        if (j < ny - 1 && votes[c + sy] === 1) { votes[c + sy] = label; stack[top++] = c + sy; }
        if (k > 0 && votes[c - sz] === 1) { votes[c - sz] = label; stack[top++] = c - sz; }
        if (k < nz - 1 && votes[c + sz] === 1) { votes[c + sz] = label; stack[top++] = c + sz; }
      }
      if (size > best) { best = size; bestLabel = label; }
    }
    for (let i = 0; i < nCellGrid; i++) {
      if (votes[i] === 0) continue;
      if (votes[i] === bestLabel) votes[i] = 1; else { votes[i] = 0; removed++; }
    }
  }
  let nCells = 0;
  for (let i = 0; i < nCellGrid; i++) nCells += votes[i];
  if (nCells === 0) throw new Error('Voxel mesher found no interior cells: the surface is not closed, or h is larger than the body.');
  const gx = nx + 1; const gy = ny + 1;
  const nodeId = new Int32Array(gx * gy * (nz + 1)).fill(-1);
  const conn = new Uint32Array(8 * nCells);
  const coords = [];
  let nn = 0;
  const node = (i, j, k) => {
    const g = (k * gy + j) * gx + i;
    if (nodeId[g] < 0) { nodeId[g] = nn++; coords.push(origin[0] + i * h, origin[1] + j * h, origin[2] + k * h); }
    return nodeId[g];
  };
  const faces = { '-x': [], '+x': [], '-y': [], '+y': [], '-z': [], '+z': [] };
  const extreme = { xmin: [], xmax: [], ymin: [], ymax: [], zmin: [], zmax: [] };
  const lo = [nx, ny, nz];
  const hi = [-1, -1, -1];
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    if (!votes[(k * ny + j) * nx + i]) continue;
    if (i < lo[0]) lo[0] = i; if (j < lo[1]) lo[1] = j; if (k < lo[2]) lo[2] = k;
    if (i > hi[0]) hi[0] = i; if (j > hi[1]) hi[1] = j; if (k > hi[2]) hi[2] = k;
  }
  let e = 0;
  const hn = new Uint32Array(8);
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const c = (k * ny + j) * nx + i;
    if (!votes[c]) continue;
    hn[0] = node(i, j, k); hn[1] = node(i + 1, j, k); hn[2] = node(i + 1, j + 1, k); hn[3] = node(i, j + 1, k);
    hn[4] = node(i, j, k + 1); hn[5] = node(i + 1, j, k + 1); hn[6] = node(i + 1, j + 1, k + 1); hn[7] = node(i, j + 1, k + 1);
    conn.set(hn, 8 * e);
    e++;
    const nb = [i > 0 && votes[c - sx], i < nx - 1 && votes[c + sx], j > 0 && votes[c - sy], j < ny - 1 && votes[c + sy], k > 0 && votes[c - sz], k < nz - 1 && votes[c + sz]];
    const names = ['-x', '+x', '-y', '+y', '-z', '+z'];
    const ext2 = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'];
    const atExt = [i === lo[0], i === hi[0], j === lo[1], j === hi[1], k === lo[2], k === hi[2]];
    for (let f = 0; f < 6; f++) {
      if (nb[f]) continue;
      const q = HEX_FACES[f];
      faces[names[f]].push(hn[q[0]], hn[q[1]], hn[q[2]], hn[q[3]]);
      if (atExt[f]) extreme[ext2[f]].push(hn[q[0]], hn[q[1]], hn[q[2]], hn[q[3]]);
    }
  }
  const positions = Float64Array.from(coords);
  const sets = {};
  const all = [];
  let nFaces = 0;
  for (const k of Object.keys(faces)) {
    sets[k] = finishFaceSet({ quad4: Uint32Array.from(faces[k]) });
    nFaces += faces[k].length / 4;
    for (let i = 0; i < faces[k].length; i++) all.push(faces[k][i]);
  }
  for (const k of Object.keys(extreme)) sets[k] = finishFaceSet({ quad4: Uint32Array.from(extreme[k]) });
  sets.boundary = finishFaceSet({ quad4: Uint32Array.from(all) });
  const model = {
    kind: 'volume', name, positions, cells: [{ type: 'hex8', conn }], sets,
    grid: { h, origin, nx, ny, nz },
  };
  for (const pl of planes) {
    if (!pl || typeof pl.name !== 'string' || !pl.point || !pl.normal) throw new RangeError(`${fn}: each plane needs name, point and normal`);
    const nl = Math.hypot(pl.normal[0], pl.normal[1], pl.normal[2]);
    if (!(nl > 0)) throw new RangeError(`${fn}: plane ${pl.name} has a zero normal`);
    const n = [pl.normal[0] / nl, pl.normal[1] / nl, pl.normal[2] / nl];
    const tol = pl.tol === undefined ? 0.75 * h : pl.tol;
    sets[pl.name] = selectFaces(model, (f) => {
      const dist = (f.centroid[0] - pl.point[0]) * n[0] + (f.centroid[1] - pl.point[1]) * n[1] + (f.centroid[2] - pl.point[2]) * n[2];
      return Math.abs(dist) <= tol && f.normal[0] * n[0] + f.normal[1] * n[1] + f.normal[2] * n[2] > 0.5;
    }, sets.boundary);
  }
  const volume = nCells * h * h * h;
  model.stats = {
    nCells, nNodes: nn, nDof: 3 * nn, h, volume, surfaceVolume: Math.abs(sv.volume),
    volumeError: Math.abs(sv.volume) > 0 ? (volume - Math.abs(sv.volume)) / Math.abs(sv.volume) : NaN,
    boundaryArea: nFaces * h * h, voteDisagreement: candidates ? disagree / candidates : 0, removedCells: removed,
    bbox: sv.bbox,
  };
  return model;
}

/* ------------------------------------------------------------------------------------------------------------
 * Matrix-free system for uniform voxel meshes
 * ---------------------------------------------------------------------------------------------------------- */

function ebeSystem({ nn, ne, conn, ke, fixed, springs, tol, maxIter }) {
  const ndof = 3 * nn;
  const isFixed = new Uint8Array(ndof);
  let nFixed = 0;
  for (let i = 0; i < ndof; i++) if (!Number.isNaN(fixed[i])) { isFixed[i] = 1; nFixed++; }
  const xe = new Float64Array(24);
  const mulK = (x, y) => {
    const out = y || new Float64Array(ndof);
    out.fill(0);
    for (let e = 0; e < ne; e++) {
      const c = 8 * e;
      for (let a = 0; a < 8; a++) { const d = 3 * conn[c + a]; xe[3 * a] = x[d]; xe[3 * a + 1] = x[d + 1]; xe[3 * a + 2] = x[d + 2]; }
      // Three rows at a time with the 3×3 node blocks written out: about twice as fast as a plain 24×24 loop.
      for (let a = 0; a < 8; a++) {
        let s0 = 0; let s1 = 0; let s2 = 0;
        const r = 72 * a;
        for (let q = 0; q < 24; q += 3) {
          const x0 = xe[q]; const x1 = xe[q + 1]; const x2 = xe[q + 2];
          s0 += ke[r + q] * x0 + ke[r + q + 1] * x1 + ke[r + q + 2] * x2;
          s1 += ke[r + 24 + q] * x0 + ke[r + 25 + q] * x1 + ke[r + 26 + q] * x2;
          s2 += ke[r + 48 + q] * x0 + ke[r + 49 + q] * x1 + ke[r + 50 + q] * x2;
        }
        const d = 3 * conn[c + a];
        out[d] += s0; out[d + 1] += s1; out[d + 2] += s2;
      }
    }
    if (springs) for (let i = 0; i < ndof; i++) out[i] += springs[i] * x[i];
    return out;
  };
  // Node-block (3×3) Jacobi preconditioner.
  const blk = new Float64Array(9 * nn);
  for (let e = 0; e < ne; e++) {
    for (let a = 0; a < 8; a++) {
      const n = conn[8 * e + a];
      for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) blk[9 * n + 3 * r + c] += ke[(3 * a + r) * 24 + 3 * a + c];
    }
  }
  const b3 = new Float64Array(9);
  for (let n = 0; n < nn; n++) {
    for (let k = 0; k < 9; k++) b3[k] = blk[9 * n + k];
    for (let r = 0; r < 3; r++) {
      if (springs) b3[4 * r] += springs[3 * n + r];
      if (isFixed[3 * n + r] || b3[4 * r] === 0) {
        for (let c = 0; c < 3; c++) { b3[3 * r + c] = 0; b3[3 * c + r] = 0; }
        b3[4 * r] = 1;
      }
    }
    invertSmall(b3, 3);
    for (let k = 0; k < 9; k++) blk[9 * n + k] = b3[k];
  }
  const precond = (r, z) => {
    for (let n = 0; n < nn; n++) {
      const a = 3 * n; const b = 9 * n;
      const r0 = r[a]; const r1 = r[a + 1]; const r2 = r[a + 2];
      z[a] = blk[b] * r0 + blk[b + 1] * r1 + blk[b + 2] * r2;
      z[a + 1] = blk[b + 3] * r0 + blk[b + 4] * r1 + blk[b + 5] * r2;
      z[a + 2] = blk[b + 6] * r0 + blk[b + 7] * r1 + blk[b + 8] * r2;
    }
  };
  const sys = { ndof, solver: 'ebe', isFixed, nFixed, lastIterations: 0, mulK };
  const cg = (x, f) => {
    const r = mulK(x);
    let fn2 = 0;
    for (let i = 0; i < ndof; i++) {
      r[i] = isFixed[i] ? 0 : f[i] - r[i];
      fn2 += r[i] * r[i];
    }
    const z = new Float64Array(ndof);
    const p = new Float64Array(ndof);
    const q = new Float64Array(ndof);
    precond(r, z);
    let rz = 0;
    for (let i = 0; i < ndof; i++) { p[i] = z[i]; rz += r[i] * z[i]; }
    const target = tol * tol * fn2;
    const cap = maxIter || Math.max(3000, 2 * ndof);
    let it = 0;
    let rr = fn2;
    while (rr > target && it < cap) {
      mulK(p, q);
      let pq = 0;
      for (let i = 0; i < ndof; i++) { if (isFixed[i]) q[i] = 0; pq += p[i] * q[i]; }
      if (!(pq > 0)) throw new Error('Stiffness operator is not positive definite: the model is insufficiently supported.');
      const alpha = rz / pq;
      rr = 0;
      for (let i = 0; i < ndof; i++) { x[i] += alpha * p[i]; r[i] -= alpha * q[i]; rr += r[i] * r[i]; }
      precond(r, z);
      let rzNew = 0;
      for (let i = 0; i < ndof; i++) rzNew += r[i] * z[i];
      const beta = rzNew / rz;
      rz = rzNew;
      for (let i = 0; i < ndof; i++) p[i] = z[i] + beta * p[i];
      it++;
    }
    sys.lastIterations = it;
    if (rr > target) throw new Error(`Iterative solver did not converge in ${it} iterations. Check supports for rigid-body modes.`);
    return x;
  };
  sys.solve = (f) => {
    const x = new Float64Array(ndof);
    for (let i = 0; i < ndof; i++) if (isFixed[i]) x[i] = fixed[i];
    return cg(x, f);
  };
  sys.solveHomogeneous = (b) => cg(new Float64Array(ndof), b);
  return sys;
}

/* ------------------------------------------------------------------------------------------------------------
 * Loads, supports and solution
 * ---------------------------------------------------------------------------------------------------------- */

function resolveFaces(model, spec, fn) {
  if (spec.set !== undefined) {
    const s = model.sets && model.sets[spec.set];
    if (!s) throw new RangeError(`${fn}: model has no set named '${spec.set}'`);
    return s;
  }
  if (spec.faces) return spec.faces;
  throw new RangeError(`${fn}: a face load/support needs 'set' or 'faces'`);
}

function resolveNodes(model, spec, fn) {
  if (spec.nodes) return spec.nodes;
  if (spec.node !== undefined) return [spec.node];
  if (spec.set !== undefined || spec.faces) return resolveFaces(model, spec, fn).nodes;
  throw new RangeError(`${fn}: entry needs 'nodes', 'node', 'set' or 'faces'`);
}

const Q2 = gauss1d(2);
const T7 = gaussTriangle(7);

/**
 * Integrate a traction over a face set: calls add(node, weightVector) with ∫ N_i t dA contributions.
 * tractionAt(x, y, z, n) returns the traction vector [Pa] given position and unit outward normal.
 */
function integrateFaces(pos, set, tractionAt, f) {
  const N = [0, 0, 0, 0, 0, 0];
  const dA = [0, 0, 0, 0, 0, 0];
  const dB = [0, 0, 0, 0, 0, 0];
  const nrm = [0, 0, 0];
  const doFace = (ids, n, pts) => {
    for (const [a, b, w] of pts) {
      if (n === 4) {
        for (let i = 0; i < 4; i++) {
          const xi = [-1, 1, 1, -1][i]; const et = [-1, -1, 1, 1][i];
          N[i] = 0.25 * (1 + xi * a) * (1 + et * b); dA[i] = 0.25 * xi * (1 + et * b); dB[i] = 0.25 * et * (1 + xi * a);
        }
      } else if (n === 3) {
        N[0] = 1 - a - b; N[1] = a; N[2] = b;
        dA[0] = -1; dA[1] = 1; dA[2] = 0; dB[0] = -1; dB[1] = 0; dB[2] = 1;
      } else {
        const l = 1 - a - b;
        N[0] = l * (2 * l - 1); N[1] = a * (2 * a - 1); N[2] = b * (2 * b - 1); N[3] = 4 * l * a; N[4] = 4 * a * b; N[5] = 4 * b * l;
        dA[0] = 1 - 4 * l; dA[1] = 4 * a - 1; dA[2] = 0; dA[3] = 4 * (l - a); dA[4] = 4 * b; dA[5] = -4 * b;
        dB[0] = 1 - 4 * l; dB[1] = 0; dB[2] = 4 * b - 1; dB[3] = -4 * a; dB[4] = 4 * a; dB[5] = 4 * (l - b);
      }
      let x = 0; let y = 0; let z = 0; let ax = 0; let ay = 0; let az = 0; let bx = 0; let by = 0; let bz = 0;
      for (let i = 0; i < n; i++) {
        const X = pos[3 * ids[i]]; const Y = pos[3 * ids[i] + 1]; const Z = pos[3 * ids[i] + 2];
        x += N[i] * X; y += N[i] * Y; z += N[i] * Z;
        ax += dA[i] * X; ay += dA[i] * Y; az += dA[i] * Z;
        bx += dB[i] * X; by += dB[i] * Y; bz += dB[i] * Z;
      }
      nrm[0] = ay * bz - az * by; nrm[1] = az * bx - ax * bz; nrm[2] = ax * by - ay * bx;
      const jac = Math.hypot(nrm[0], nrm[1], nrm[2]);
      if (jac === 0) continue;
      nrm[0] /= jac; nrm[1] /= jac; nrm[2] /= jac;
      const t = tractionAt(x, y, z, nrm);
      for (let i = 0; i < n; i++) {
        const d = 3 * ids[i];
        const s = N[i] * w * jac;
        f[d] += s * t[0]; f[d + 1] += s * t[1]; f[d + 2] += s * t[2];
      }
    }
  };
  const qp = [];
  for (let j = 0; j < 2; j++) for (let i = 0; i < 2; i++) qp.push([Q2.x[i], Q2.x[j], 1]);
  const tp = T7.xi.map((v, i) => [v, T7.eta[i], T7.w[i]]);
  for (const [k, n, pts] of [['quad4', 4, qp], ['tri3', 3, tp], ['tri6', 6, tp]]) {
    const arr = set[k] || [];
    for (let i = 0; i < arr.length; i += n) doFace(arr.subarray(i, i + n), n, pts);
  }
}

/** Principal stresses of a symmetric tensor (closed form), returned descending into out[0..2]. */
function principal(sx, sy, sz, txy, tyz, tzx, out) {
  const m = (sx + sy + sz) / 3;
  const dx = sx - m; const dy = sy - m; const dz = sz - m;
  const j2 = 0.5 * (dx * dx + dy * dy + dz * dz) + txy * txy + tyz * tyz + tzx * tzx;
  if (j2 < 1e-300) { out[0] = m; out[1] = m; out[2] = m; return; }
  const j3 = dx * (dy * dz - tyz * tyz) - txy * (txy * dz - tyz * tzx) + tzx * (txy * tyz - dy * tzx);
  const r = 2 * Math.sqrt(j2 / 3);
  let c = (1.5 * Math.sqrt(3) * j3) / Math.pow(j2, 1.5);
  if (c > 1) c = 1; else if (c < -1) c = -1;
  const th = Math.acos(c) / 3;
  out[0] = m + r * Math.cos(th);
  out[1] = m + r * Math.cos(th - (2 * Math.PI) / 3);
  out[2] = m + r * Math.cos(th + (2 * Math.PI) / 3);
}

/**
 * Principal stresses of a symmetric stress tensor.
 * @param {object} s
 * @param {number} s.sxx normal stress [Pa]
 * @param {number} s.syy normal stress [Pa]
 * @param {number} s.szz normal stress [Pa]
 * @param {number} [s.sxy] shear stress [Pa]
 * @param {number} [s.syz] shear stress [Pa]
 * @param {number} [s.szx] shear stress [Pa]
 * @returns {{s1: number, s2: number, s3: number, vonMises: number, tresca: number}} principal stresses
 *   (s1 >= s2 >= s3), von Mises equivalent and Tresca (s1 − s3) [Pa]
 */
export function principalStresses({ sxx, syy, szz, sxy = 0, syz = 0, szx = 0 }) {
  for (const [v, name] of [[sxx, 'sxx'], [syy, 'syy'], [szz, 'szz']]) if (!Number.isFinite(v)) throw new RangeError(`principalStresses: ${name} must be finite`);
  const out = [0, 0, 0];
  principal(sxx, syy, szz, sxy, syz, szx, out);
  const vm = Math.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * (sxy * sxy + syz * syz + szx * szx));
  return { s1: out[0], s2: out[1], s3: out[2], vonMises: vm, tresca: out[0] - out[2] };
}

function setupSolid(o, fn) {
  const { model, E, nu, hex = 'incompatible' } = o;
  positive(E, 'E', fn);
  const D = solidD(E, nu);
  if (hex !== 'incompatible' && hex !== 'full') throw new RangeError(`${fn}: hex must be 'incompatible' or 'full'`);
  const P = prepare(model, fn);
  const wk = makeWork();
  const xe = new Float64Array(30);
  const loadXe = (e) => {
    const c = P.eptr[e];
    const nen = P.eptr[e + 1] - c;
    for (let a = 0; a < nen; a++) { const n = P.econn[c + a]; xe[3 * a] = P.pos[3 * n]; xe[3 * a + 1] = P.pos[3 * n + 1]; xe[3 * a + 2] = P.pos[3 * n + 2]; }
    return xe;
  };
  const elemNodes = (e) => P.econn.subarray(P.eptr[e], P.eptr[e + 1]);
  let cached = null;
  let cachedS = null;
  if (P.uniform) {
    // Every voxel is the same cube: one element matrix and one stress operator serve the whole mesh.
    cached = Float64Array.from(elementK('hex8', loadXe(0), D, hex, wk, true));
    cachedS = Float64Array.from(wk.S.subarray(0, 8 * 6 * 24));
  }
  const elemMatrix = (e) => cached || elementK(TYPE_NAME[P.etype[e]], loadXe(e), D, hex, wk, false);
  const stressOp = (e) => {
    if (cachedS) return cachedS;
    elementK(TYPE_NAME[P.etype[e]], loadXe(e), D, hex, wk, true);
    return wk.S;
  };
  // Nodes not referenced by any element would leave zero pivots; they are held automatically.
  const used = new Uint8Array(P.nn);
  for (let i = 0; i < P.econn.length; i++) used[P.econn[i]] = 1;
  return { P, D, wk, xe, loadXe, elemNodes, elemMatrix, stressOp, cached, used };
}

function buildSupports(model, S, fix, springsSpec, foundation, fn) {
  const { P } = S;
  const ndof = 3 * P.nn;
  const fixed = new Float64Array(ndof).fill(NaN);
  const keys = ['ux', 'uy', 'uz'];
  for (const s of fix || []) {
    const nodes = resolveNodes(model, s, fn);
    for (let i = 0; i < nodes.length; i++) {
      const n = nodes[i];
      if (!(n >= 0 && n < P.nn)) throw new RangeError(`${fn}: fix node index ${n} out of range`);
      for (let k = 0; k < 3; k++) {
        const v = s[keys[k]];
        if (typeof v === 'number') fixed[3 * n + k] = v; else if (v === true) fixed[3 * n + k] = 0;
      }
    }
  }
  let nUser = 0;
  for (let i = 0; i < ndof; i++) if (!Number.isNaN(fixed[i])) nUser++;
  for (let n = 0; n < P.nn; n++) if (!S.used[n]) { fixed[3 * n] = 0; fixed[3 * n + 1] = 0; fixed[3 * n + 2] = 0; }
  let springs = null;
  const kk = ['kx', 'ky', 'kz'];
  for (const s of springsSpec || []) {
    springs = springs || new Float64Array(ndof);
    const nodes = resolveNodes(model, s, fn);
    for (let i = 0; i < nodes.length; i++) for (let k = 0; k < 3; k++) {
      const v = s[kk[k]] || 0;
      if (v < 0) throw new RangeError(`${fn}: spring stiffness must be >= 0`);
      springs[3 * nodes[i] + k] += v;
    }
  }
  for (const s of foundation || []) {
    springs = springs || new Float64Array(ndof);
    const set = resolveFaces(model, s, fn);
    const kv = [s.kx || 0, s.ky || 0, s.kz || 0];
    if (kv.some((v) => v < 0)) throw new RangeError(`${fn}: foundation stiffness must be >= 0`);
    // ∫ N_i dA lumps the distributed stiffness [N/m³] to nodal springs [N/m].
    integrateFaces(P.pos, set, () => kv, springs);
  }
  return { fixed, springs, nUser };
}

/**
 * Linear static analysis of a solid body.
 *
 * @param {object} o
 * @param {object} o.model volume model (tet4 / tet10 / hex8 cells; `voxelMesh` or `boxMesh` output, or imported)
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {Array<object>} o.fix supports: { nodes | node | set | faces, ux?, uy?, uz? } prescribed displacement [m]
 *   (set = name in model.sets; a key that is present constrains that direction)
 * @param {Array<object>} [o.springs] nodal springs to ground { nodes | set, kx, ky, kz } [N/m per node]
 * @param {Array<object>} [o.foundation] elastic foundation on faces { set | faces, kx, ky, kz } [N/m³]
 * @param {Array<object>} [o.nodalLoads] { node | nodes, fx, fy, fz } [N] (per node)
 * @param {Array<object>} [o.faceLoads] { set | faces, pressure?, traction? }: pressure [Pa] (number or function
 *   of (x, y, z)) acts against the outward normal; traction [tx, ty, tz] [Pa] in global axes
 * @param {Array<object>} [o.hydrostatic] { set | faces, rho [kg/m³], g [m/s²] (default 9.80665), zSurface [m] }:
 *   pressure rho·g·(zSurface − z) on faces below the free surface (z up)
 * @param {{rho: number, g: number[]}} [o.gravity] self-weight: density [kg/m³] and acceleration vector [m/s²]
 * @param {number[]} [o.bodyForce] body force density [N/m³]
 * @param {'incompatible'|'full'} [o.hex] hex8 formulation (default 'incompatible')
 * @param {'auto'|'skyline'|'pcg'|'ebe'} [o.solver] 'ebe' = matrix-free CG, uniform voxel meshes only
 * @param {number} [o.tol] iterative relative residual tolerance (default 1e-8)
 * @param {number} [o.maxIter] iteration cap
 * @param {boolean} [o.stress] recover stresses (default true)
 * @returns {{u: Float64Array, reactions: Float64Array, reactionSum: {fx: number, fy: number, fz: number},
 *   appliedSum: {fx: number, fy: number, fz: number}, equilibriumResidual: number, stress: object|null,
 *   vonMises: Float64Array|null, maxPrincipal: Float64Array|null, minPrincipal: Float64Array|null,
 *   maxVonMises: number, maxDisp: number, strainEnergy: number, solver: string, iterations: number, nDof: number}}
 *   u [m] (3 per node); reactions [N] per dof (supports and springs); equilibriumResidual = ‖K u − f‖ over free
 *   dofs divided by the larger of ‖f‖ and ‖reactions‖ [-]; stress = nodal-averaged { sxx, syy, szz, sxy, syz, szx }
 *   [Pa]; strainEnergy [J]
 */
export function solveSolid(o) {
  const fn = 'solveSolid';
  const { model, nodalLoads = [], faceLoads = [], hydrostatic = [], gravity = null, bodyForce = null, tol = 1e-8, maxIter, stress = true } = o;
  const S = setupSolid(o, fn);
  const { P, wk } = S;
  const ndof = 3 * P.nn;
  const f = new Float64Array(ndof);
  for (const L of nodalLoads) {
    const nodes = resolveNodes(model, L, fn);
    for (let i = 0; i < nodes.length; i++) {
      const n = nodes[i];
      if (!(n >= 0 && n < P.nn)) throw new RangeError(`${fn}: nodalLoads node index out of range`);
      f[3 * n] += L.fx || 0; f[3 * n + 1] += L.fy || 0; f[3 * n + 2] += L.fz || 0;
    }
  }
  const tv = [0, 0, 0];
  for (const L of faceLoads) {
    const set = resolveFaces(model, L, fn);
    const tr = L.traction || [0, 0, 0];
    integrateFaces(P.pos, set, (x, y, z, n) => {
      const p = typeof L.pressure === 'function' ? L.pressure(x, y, z) : (L.pressure || 0);
      tv[0] = tr[0] - p * n[0]; tv[1] = tr[1] - p * n[1]; tv[2] = tr[2] - p * n[2];
      return tv;
    }, f);
  }
  for (const L of hydrostatic) {
    const set = resolveFaces(model, L, fn);
    positive(L.rho, 'hydrostatic.rho', fn);
    if (!Number.isFinite(L.zSurface)) throw new RangeError(`${fn}: hydrostatic.zSurface must be a number`);
    const gg = L.g === undefined ? 9.80665 : L.g;
    integrateFaces(P.pos, set, (x, y, z, n) => {
      const p = z < L.zSurface ? L.rho * gg * (L.zSurface - z) : 0;
      tv[0] = -p * n[0]; tv[1] = -p * n[1]; tv[2] = -p * n[2];
      return tv;
    }, f);
  }
  const b = [0, 0, 0];
  if (gravity) { positive(gravity.rho, 'gravity.rho', fn); for (let c = 0; c < 3; c++) b[c] += gravity.rho * gravity.g[c]; }
  if (bodyForce) for (let c = 0; c < 3; c++) b[c] += bodyForce[c];
  if (b[0] !== 0 || b[1] !== 0 || b[2] !== 0) {
    for (let e = 0; e < P.ne; e++) {
      const type = TYPE_NAME[P.etype[e]];
      const nen = NEN[type];
      const xe = S.loadXe(e);
      const rule = RULE[type];
      const c0 = P.eptr[e];
      for (let g = 0; g < rule.w.length; g++) {
        geom3(type, nen, rule.xi[g], rule.eta[g], rule.zeta[g], xe, wk.sc);
        const w = rule.w[g] * wk.sc.det;
        for (let a = 0; a < nen; a++) {
          const d = 3 * P.econn[c0 + a];
          const s = wk.sc.N[a] * w;
          f[d] += s * b[0]; f[d + 1] += s * b[1]; f[d + 2] += s * b[2];
        }
      }
    }
  }
  const sup = buildSupports(model, S, o.fix, o.springs, o.foundation, fn);
  if (sup.nUser === 0 && !sup.springs) throw new Error('Solid model has no supports: add displacement constraints or springs.');
  let solver = o.solver || 'auto';
  if (solver === 'ebe' && !S.cached) throw new RangeError(`${fn}: solver 'ebe' needs a uniform voxel mesh (model.grid)`);
  if (solver === 'auto' && S.cached && P.ne > 1500) solver = 'ebe';
  const sys = solver === 'ebe'
    ? ebeSystem({ nn: P.nn, ne: P.ne, conn: P.econn, ke: S.cached, fixed: sup.fixed, springs: sup.springs, tol, maxIter })
    : feSystem({ nNodes: P.nn, dpn: 3, nElem: P.ne, elemNodes: S.elemNodes, elemMatrix: S.elemMatrix, fixed: sup.fixed, springs: sup.springs, solver, tol: Math.min(tol, 1e-9), maxIter });
  const u = sys.solve(f);
  const ku = sys.mulK(u);
  const reactions = new Float64Array(ndof);
  const rs = [0, 0, 0];
  const as = [0, 0, 0];
  let energy = 0; let res2 = 0; let f2 = 0; let r2 = 0;
  for (let i = 0; i < ndof; i++) {
    energy += 0.5 * u[i] * ku[i];
    const spr = sup.springs ? sup.springs[i] * u[i] : 0;
    if (sys.isFixed[i]) {
      reactions[i] = ku[i] - f[i];
    } else {
      const r = ku[i] - f[i];
      res2 += r * r;
      // Spring reaction on the structure opposes the displacement.
      reactions[i] = -spr;
    }
    rs[i % 3] += reactions[i];
    as[i % 3] += f[i];
    f2 += f[i] * f[i];
    r2 += reactions[i] * reactions[i];
  }
  let maxDisp = 0;
  for (let n = 0; n < P.nn; n++) maxDisp = Math.max(maxDisp, Math.hypot(u[3 * n], u[3 * n + 1], u[3 * n + 2]));
  const out = {
    u, reactions, reactionSum: { fx: rs[0], fy: rs[1], fz: rs[2] }, appliedSum: { fx: as[0], fy: as[1], fz: as[2] },
    equilibriumResidual: Math.sqrt(res2) / (Math.sqrt(Math.max(f2, r2)) || 1),
    stress: null, vonMises: null, maxPrincipal: null, minPrincipal: null, maxVonMises: 0, maxDisp,
    strainEnergy: energy, solver: sys.solver, iterations: sys.lastIterations, nDof: ndof,
  };
  if (stress) Object.assign(out, recoverStress(S, u));
  return out;
}

function recoverStress(S, u) {
  const { P } = S;
  const nn = P.nn;
  const acc = [0, 1, 2, 3, 4, 5].map(() => new Float64Array(nn));
  const cnt = new Uint16Array(nn);
  const ue = new Float64Array(30);
  const sg = new Float64Array(48);
  for (let e = 0; e < P.ne; e++) {
    const type = TYPE_NAME[P.etype[e]];
    const nen = NEN[type];
    const m = 3 * nen;
    const c0 = P.eptr[e];
    const Sop = S.stressOp(e);
    for (let a = 0; a < nen; a++) { const d = 3 * P.econn[c0 + a]; ue[3 * a] = u[d]; ue[3 * a + 1] = u[d + 1]; ue[3 * a + 2] = u[d + 2]; }
    const ng = RULE[type].w.length;
    for (let r = 0; r < 6 * ng; r++) {
      let s = 0;
      const row = r * m;
      for (let q = 0; q < m; q++) s += Sop[row + q] * ue[q];
      sg[r] = s;
    }
    const ex = EXTRAP[type];
    for (let a = 0; a < nen; a++) {
      const n = P.econn[c0 + a];
      cnt[n]++;
      for (let g = 0; g < ng; g++) {
        const w = ex[a][g];
        for (let r = 0; r < 6; r++) acc[r][n] += w * sg[6 * g + r];
      }
    }
  }
  const vonMises = new Float64Array(nn);
  const maxPrincipal = new Float64Array(nn);
  const minPrincipal = new Float64Array(nn);
  const pr = [0, 0, 0];
  let maxVm = 0;
  for (let n = 0; n < nn; n++) {
    const c = cnt[n] || 1;
    for (let r = 0; r < 6; r++) acc[r][n] /= c;
    const sx = acc[0][n]; const sy = acc[1][n]; const sz = acc[2][n]; const txy = acc[3][n]; const tyz = acc[4][n]; const tzx = acc[5][n];
    const vm = Math.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3 * (txy * txy + tyz * tyz + tzx * tzx));
    vonMises[n] = vm;
    if (vm > maxVm) maxVm = vm;
    principal(sx, sy, sz, txy, tyz, tzx, pr);
    maxPrincipal[n] = pr[0];
    minPrincipal[n] = pr[2];
  }
  return { stress: { sxx: acc[0], syy: acc[1], szz: acc[2], sxy: acc[3], syz: acc[4], szx: acc[5] }, vonMises, maxPrincipal, minPrincipal, maxVonMises: maxVm };
}

/**
 * Lowest natural frequencies of a supported solid body (consistent mass, subspace iteration on a direct
 * factorisation — intended for small models, up to a few tens of thousands of dofs).
 * @param {object} o
 * @param {object} o.model volume model
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} o.rho density [kg/m³]
 * @param {Array<object>} o.fix supports as in `solveSolid` (the body must be restrained; free-free bodies are
 *   not supported)
 * @param {Array<object>} [o.springs] nodal springs as in `solveSolid`
 * @param {number} [o.nModes] number of modes (default 4)
 * @param {'incompatible'|'full'} [o.hex] hex8 formulation
 * @returns {{frequencies: number[], omega: number[], modes: Float64Array[], totalMass: number}} natural
 *   frequencies [Hz], circular frequencies [rad/s], mode shapes (3 per node, unit 2-norm), model mass [kg]
 */
export function modalSolid(o) {
  const fn = 'modalSolid';
  const { model, rho, nModes = 4 } = o;
  positive(rho, 'rho', fn);
  const S = setupSolid(o, fn);
  const { P, wk } = S;
  const ndof = 3 * P.nn;
  const sup = buildSupports(model, S, o.fix, o.springs, null, fn);
  if (sup.nUser === 0 && !sup.springs) throw new Error('Modal analysis needs a supported model.');
  const sys = feSystem({ nNodes: P.nn, dpn: 3, nElem: P.ne, elemNodes: S.elemNodes, elemMatrix: S.elemMatrix, fixed: sup.fixed, springs: sup.springs, solver: 'skyline' });
  // Scalar consistent mass per element (same for the three components). The 2×2×2 and 4-point rules
  // under-integrate N·N for tet10, so tetrahedra use the exact barycentric product formula instead.
  const mptr = new Int32Array(P.ne + 1);
  for (let e = 0; e < P.ne; e++) { const n = P.eptr[e + 1] - P.eptr[e]; mptr[e + 1] = mptr[e] + n * n; }
  const mm = new Float64Array(mptr[P.ne]);
  let totalMass = 0;
  const t10 = tet10MassRef();
  for (let e = 0; e < P.ne; e++) {
    const type = TYPE_NAME[P.etype[e]];
    const nen = NEN[type];
    const xe = S.loadXe(e);
    const base = mptr[e];
    if (type === 'hex8') {
      const rule = RULE.hex8;
      for (let g = 0; g < 8; g++) {
        geom3(type, nen, rule.xi[g], rule.eta[g], rule.zeta[g], xe, wk.sc);
        const w = rho * wk.sc.det;
        totalMass += w;
        for (let a = 0; a < 8; a++) for (let b = 0; b < 8; b++) mm[base + a * 8 + b] += w * wk.sc.N[a] * wk.sc.N[b];
      }
    } else {
      geom3('tet4', 4, 0.25, 0.25, 0.25, xe, wk.sc);
      const mass = (rho * wk.sc.det) / 6;
      totalMass += mass;
      if (type === 'tet4') {
        for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) mm[base + a * 4 + b] = (mass * (a === b ? 2 : 1)) / 20;
      } else {
        for (let k = 0; k < 100; k++) mm[base + k] = mass * t10[k];
      }
    }
  }
  const mulB = (x) => {
    const y = new Float64Array(ndof);
    for (let e = 0; e < P.ne; e++) {
      const c0 = P.eptr[e];
      const nen = P.eptr[e + 1] - c0;
      const base = mptr[e];
      for (let a = 0; a < nen; a++) {
        let s0 = 0; let s1 = 0; let s2 = 0;
        for (let b = 0; b < nen; b++) {
          const m = mm[base + a * nen + b];
          const d = 3 * P.econn[c0 + b];
          s0 += m * x[d]; s1 += m * x[d + 1]; s2 += m * x[d + 2];
        }
        const d = 3 * P.econn[c0 + a];
        y[d] += s0; y[d + 1] += s1; y[d + 2] += s2;
      }
    }
    for (let i = 0; i < ndof; i++) if (sys.isFixed[i]) y[i] = 0;
    return y;
  };
  const r = subspaceModes({ sys, mulB, n: ndof, nModes });
  const omega = r.mu.map((m) => 1 / Math.sqrt(m));
  return { frequencies: omega.map((w) => w / (2 * Math.PI)), omega, modes: r.vectors, totalMass };
}

/**
 * Reference consistent mass of a straight-sided tet10 divided by the element mass (10×10), from the exact
 * integral ∫ L1^a L2^b L3^c L4^d dV = 6V a! b! c! d! / (a + b + c + d + 3)!.
 */
function tet10MassRef() {
  // Each shape function as a list of monomials { c, p: [a, b, c, d] } in barycentric coordinates.
  const mono = [];
  for (let i = 0; i < 4; i++) {
    const p2 = [0, 0, 0, 0]; p2[i] = 2;
    const p1 = [0, 0, 0, 0]; p1[i] = 1;
    mono.push([{ c: 2, p: p2 }, { c: -1, p: p1 }]);
  }
  for (const [a, b] of TET10_EDGES) { const p = [0, 0, 0, 0]; p[a] = 1; p[b] = 1; mono.push([{ c: 4, p }]); }
  const fact = [1, 1, 2, 6, 24, 120, 720, 5040];
  const out = new Float64Array(100);
  for (let i = 0; i < 10; i++) for (let j = 0; j < 10; j++) {
    let s = 0;
    for (const a of mono[i]) for (const b of mono[j]) {
      const q = [0, 1, 2, 3].map((k) => a.p[k] + b.p[k]);
      const tot = q[0] + q[1] + q[2] + q[3];
      s += a.c * b.c * (6 * fact[q[0]] * fact[q[1]] * fact[q[2]] * fact[q[3]]) / fact[tot + 3];
    }
    out[i * 10 + j] = s;
  }
  return out;
}
