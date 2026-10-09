/**
 * Three-dimensional elastic frame analysis by the direct stiffness method.
 *
 * Element: 12-dof prismatic beam with optional Timoshenko shear deformation (Przemieniecki, Theory of Matrix
 * Structural Analysis, 1968, §5.6; Friedman & Kosmatka 1993 for the interdependent shape functions used for
 * consistent loads). Member releases by static condensation, rigid end offsets, truss and tension-only /
 * compression-only members (iterative), nodal springs (diagonal or full 6×6), prescribed support displacements.
 *
 * Conventions
 *  - Global axes X, Y, Z with Z up. Node dofs: ux, uy, uz, rx, ry, rz (right-hand rule) → global dof 6·node + i.
 *  - Local x runs from end A to end B (between the offset ends). Default local z lies in the vertical plane
 *    through the member and points upward (for vertical members local y = global Y); local y = z × x. `betaDeg`
 *    rotates y, z about x; `localZ` (a global vector) replaces the default reference direction.
 *  - Iy is the second moment about local y (bending in the local x–z plane, moment My); Iz about local z.
 *  - Section forces N, Vy, Vz, T, My, Mz act on the cut face whose outward normal is +x (N > 0 is tension).
 *  - SI units: m, N, N·m, Pa, kg, rad.
 */

import { SkylineSym, reverseCuthillMcKee } from '../math/sparse.js';
import { gaussLegendre } from '../math/quad.js';
import { resolveSection } from './sections.js';

const DOF_NAMES = ['ux', 'uy', 'uz', 'rx', 'ry', 'rz'];
const G0 = 9.80665;

/** Default structural steel: E [Pa], nu [-], rho [kg/m³], alpha [1/K], fy [Pa]. */
export const DEFAULT_STEEL = Object.freeze({ E: 2.1e11, nu: 0.3, rho: 7850, alpha: 1.2e-5, fy: 355e6 });

let GAUSS8 = null;
function gauss8() {
  if (!GAUSS8) GAUSS8 = gaussLegendre(8);
  return GAUSS8;
}

/* ------------------------------------------------------------------ small dense helpers (flat 12×12) */

function mul12(A, B, out) {
  for (let i = 0; i < 12; i++) {
    for (let j = 0; j < 12; j++) {
      let s = 0;
      for (let k = 0; k < 12; k++) s += A[i * 12 + k] * B[k * 12 + j];
      out[i * 12 + j] = s;
    }
  }
  return out;
}

function tAt12(T, A, out, tmp) {
  // out = Tᵀ A T
  mul12(A, T, tmp);
  for (let i = 0; i < 12; i++) {
    for (let j = 0; j < 12; j++) {
      let s = 0;
      for (let k = 0; k < 12; k++) s += T[k * 12 + i] * tmp[k * 12 + j];
      out[i * 12 + j] = s;
    }
  }
  return out;
}

function matVec12(A, x, out) {
  for (let i = 0; i < 12; i++) {
    let s = 0;
    for (let k = 0; k < 12; k++) s += A[i * 12 + k] * x[k];
    out[i] = s;
  }
  return out;
}

function matTVec12(A, x, out) {
  for (let i = 0; i < 12; i++) {
    let s = 0;
    for (let k = 0; k < 12; k++) s += A[k * 12 + i] * x[k];
    out[i] = s;
  }
  return out;
}

/* Gaussian elimination with partial pivoting for the small condensation systems; returns the inverse or null. */
function invSmall(A, n) {
  const M = [];
  for (let i = 0; i < n; i++) {
    const row = new Float64Array(2 * n);
    for (let j = 0; j < n; j++) row[j] = A[i][j];
    row[n + i] = 1;
    M.push(row);
  }
  let scale = 0;
  for (let i = 0; i < n; i++) scale = Math.max(scale, Math.abs(A[i][i]));
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    if (Math.abs(M[p][c]) <= 1e-12 * scale || scale === 0) return null;
    if (p !== c) { const t = M[p]; M[p] = M[c]; M[c] = t; }
    const piv = M[c][c];
    for (let j = c; j < 2 * n; j++) M[c][j] /= piv;
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const f = M[r][c];
      if (f === 0) continue;
      for (let j = c; j < 2 * n; j++) M[r][j] -= f * M[c][j];
    }
  }
  return M.map((row) => Array.from(row.subarray(n)));
}

/* ------------------------------------------------------------------ model normalisation */

function nodeArray(nodes) {
  if (!nodes || nodes.length === 0) throw new RangeError('model.nodes must not be empty');
  if (typeof nodes[0] === 'number') {
    if (nodes.length % 3 !== 0) throw new RangeError('model.nodes as a flat array must hold 3 numbers per node');
    return Float64Array.from(nodes);
  }
  const out = new Float64Array(3 * nodes.length);
  for (let i = 0; i < nodes.length; i++) {
    const p = nodes[i];
    const x = Array.isArray(p) || ArrayBuffer.isView(p) ? p[0] : p.x;
    const y = Array.isArray(p) || ArrayBuffer.isView(p) ? p[1] : p.y;
    const z = Array.isArray(p) || ArrayBuffer.isView(p) ? p[2] : p.z;
    if (![x, y, z].every(Number.isFinite)) throw new RangeError(`model.nodes[${i}] has non-finite coordinates`);
    out[3 * i] = x; out[3 * i + 1] = y; out[3 * i + 2] = z;
  }
  return out;
}

function resolveMaterial(model, m) {
  let mat = m;
  if (typeof m === 'string') {
    mat = model.materials && model.materials[m];
    if (!mat) throw new RangeError(`material "${m}" is not defined in model.materials`);
  }
  if (!mat) mat = DEFAULT_STEEL;
  const E = mat.E;
  if (!(E > 0)) throw new RangeError('material.E must be positive');
  const nu = mat.nu === undefined ? 0.3 : mat.nu;
  const G = mat.G === undefined ? E / (2 * (1 + nu)) : mat.G;
  return {
    E, nu, G, rho: mat.rho === undefined ? 7850 : mat.rho, alpha: mat.alpha === undefined ? 1.2e-5 : mat.alpha,
    fy: mat.fy === undefined ? 355e6 : mat.fy,
  };
}

const SECTION_CACHE = new WeakMap();
function resolveSec(model, s) {
  let sec = s;
  if (typeof s === 'string') {
    sec = model.sections && model.sections[s];
    if (!sec) throw new RangeError(`section "${s}" is not defined in model.sections`);
  }
  if (!sec || typeof sec !== 'object') throw new RangeError('every member needs a section');
  let r = SECTION_CACHE.get(sec);
  if (!r) { r = resolveSection(sec); SECTION_CACHE.set(sec, r); }
  return r;
}

function releaseFlags(r) {
  const out = new Uint8Array(6);
  if (!r) return out;
  if (r === 'pin' || r === 'pinned') { out[4] = 1; out[5] = 1; return out; }
  if (typeof r === 'string') {
    for (const tok of r.split(/[\s,]+/)) {
      if (!tok) continue;
      const map = { ux: 0, uy: 1, uz: 2, rx: 3, ry: 4, rz: 5, n: 0, vy: 1, vz: 2, t: 3, my: 4, mz: 5 };
      const k = map[tok.toLowerCase()];
      if (k === undefined) throw new RangeError(`release "${tok}" is not recognised`);
      out[k] = 1;
    }
    return out;
  }
  for (let i = 0; i < 6; i++) out[i] = r[i] ? 1 : 0;
  return out;
}

function fixFlags(f) {
  const out = new Uint8Array(6);
  if (!f) return out;
  if (f === 'fixed') return out.fill(1);
  if (f === 'pinned') { out[0] = 1; out[1] = 1; out[2] = 1; return out; }
  if (typeof f === 'string') throw new RangeError(`support fix "${f}" is not recognised`);
  for (let i = 0; i < 6; i++) out[i] = f[i] ? 1 : 0;
  return out;
}

/**
 * Local axes of a member as a 3×3 rotation matrix (rows = local x, y, z unit vectors in global components).
 * Reference: default "vertical web" convention described in the file header.
 * @param {number[]} pa end A position [m]
 * @param {number[]} pb end B position [m]
 * @param {object} [o]
 * @param {number} [o.betaDeg=0] rotation of the y, z axes about x [deg]
 * @param {number[]} [o.localZ] global vector defining the local x–z plane
 * @returns {{R:number[][], L:number}} rotation matrix and length [m]
 */
export function memberAxes(pa, pb, { betaDeg = 0, localZ } = {}) {
  const dx = pb[0] - pa[0]; const dy = pb[1] - pa[1]; const dz = pb[2] - pa[2];
  const L = Math.hypot(dx, dy, dz);
  if (!(L > 0)) throw new RangeError('member has zero length');
  const x = [dx / L, dy / L, dz / L];
  let z;
  let y;
  const ref = localZ || [0, 0, 1];
  const d = ref[0] * x[0] + ref[1] * x[1] + ref[2] * x[2];
  const refLen = Math.hypot(ref[0], ref[1], ref[2]);
  if (!(refLen > 0)) throw new RangeError('localZ must be a non-zero vector');
  if (Math.abs(d) > 0.999999 * refLen) {
    if (localZ) throw new RangeError('localZ must not be parallel to the member axis');
    y = [0, 1, 0];
    z = [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]];
  } else {
    z = [ref[0] - d * x[0], ref[1] - d * x[1], ref[2] - d * x[2]];
    const zl = Math.hypot(z[0], z[1], z[2]);
    z = [z[0] / zl, z[1] / zl, z[2] / zl];
    y = [z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]];
  }
  if (betaDeg) {
    const b = betaDeg * Math.PI / 180;
    const c = Math.cos(b); const s = Math.sin(b);
    const y2 = [c * y[0] + s * z[0], c * y[1] + s * z[1], c * y[2] + s * z[2]];
    const z2 = [-s * y[0] + c * z[0], -s * y[1] + c * z[1], -s * y[2] + c * z[2]];
    y = y2; z = z2;
  }
  return { R: [x, y, z], L };
}

/* 4×4 bending block in (v1, θ1, v2, θ2) with shear parameter phi. */
function bend4(EI, L, phi) {
  const c = EI / (L * L * L * (1 + phi));
  const L2 = L * L;
  return [
    [12 * c, 6 * L * c, -12 * c, 6 * L * c],
    [6 * L * c, (4 + phi) * L2 * c, -6 * L * c, (2 - phi) * L2 * c],
    [-12 * c, -6 * L * c, 12 * c, -6 * L * c],
    [6 * L * c, (2 - phi) * L2 * c, -6 * L * c, (4 + phi) * L2 * c],
  ];
}

const XY_IDX = [1, 5, 7, 11];
const XZ_IDX = [2, 4, 8, 10];
const XZ_SGN = [1, -1, 1, -1];

/**
 * Local 12×12 elastic stiffness matrix of a prismatic beam (Przemieniecki 1968 Eq. 5.116 with the shear
 * parameters Φy = 12 E Iz /(G Asy L²), Φz = 12 E Iy /(G Asz L²)). Dof order: ux, uy, uz, rx, ry, rz at A then B.
 * @param {object} o
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.G shear modulus [Pa]
 * @param {number} o.A area [m²]
 * @param {number} o.Iy second moment about y [m⁴]
 * @param {number} o.Iz second moment about z [m⁴]
 * @param {number} o.J torsion constant [m⁴]
 * @param {number} o.L length [m]
 * @param {number} [o.Asy=Infinity] shear area for Vy [m²] (Infinity or 0 → Euler–Bernoulli in that plane)
 * @param {number} [o.Asz=Infinity] shear area for Vz [m²]
 * @param {boolean} [o.truss=false] axial stiffness only
 * @returns {Float64Array} row-major 12×12 matrix [N/m, N, N·m]
 */
export function beamLocalStiffness({ E, G, A, Iy, Iz, J, L, Asy = Infinity, Asz = Infinity, truss = false }) {
  if (!(L > 0)) throw new RangeError('L must be positive');
  const k = new Float64Array(144);
  const ea = E * A / L;
  k[0] = ea; k[6 * 12 + 6] = ea; k[6] = -ea; k[6 * 12] = -ea;
  if (truss) return k;
  const gj = G * J / L;
  k[3 * 12 + 3] = gj; k[9 * 12 + 9] = gj; k[3 * 12 + 9] = -gj; k[9 * 12 + 3] = -gj;
  const phiY = Asy > 0 && Number.isFinite(Asy) ? 12 * E * Iz / (G * Asy * L * L) : 0;
  const phiZ = Asz > 0 && Number.isFinite(Asz) ? 12 * E * Iy / (G * Asz * L * L) : 0;
  const bz = bend4(E * Iz, L, phiY);
  const by = bend4(E * Iy, L, phiZ);
  for (let i = 0; i < 4; i++) {
    for (let j = 0; j < 4; j++) {
      k[XY_IDX[i] * 12 + XY_IDX[j]] = bz[i][j];
      k[XZ_IDX[i] * 12 + XZ_IDX[j]] = XZ_SGN[i] * XZ_SGN[j] * by[i][j];
    }
  }
  return k;
}

/* Timoshenko shape functions for (v1, θ1, v2, θ2): out[0..3] deflection, out[4..7] section rotation. */
function bendShape(xi, L, phi, truss, out) {
  if (truss) {
    out[0] = 1 - xi; out[1] = 0; out[2] = xi; out[3] = 0;
    out[4] = -1 / L; out[5] = 0; out[6] = 1 / L; out[7] = 0;
    return out;
  }
  const d = 1 + phi;
  const x2 = xi * xi; const x3 = x2 * xi;
  out[0] = (1 - 3 * x2 + 2 * x3 + phi * (1 - xi)) / d;
  out[1] = L * (xi - 2 * x2 + x3 + phi * (xi - x2) / 2) / d;
  out[2] = (3 * x2 - 2 * x3 + phi * xi) / d;
  out[3] = L * (-x2 + x3 - phi * (xi - x2) / 2) / d;
  out[4] = 6 * (-xi + x2) / (L * d);
  out[5] = (1 - 4 * xi + 3 * x2 + phi * (1 - xi)) / d;
  out[6] = 6 * (xi - x2) / (L * d);
  out[7] = (-2 * xi + 3 * x2 + phi * xi) / d;
  return out;
}

/**
 * Build the analysis system for a frame model: local axes, element matrices (with releases condensed and
 * rigid offsets), constraint flags, reverse Cuthill–McKee equation numbering and skyline profile.
 *
 * Model format
 * ```
 * { nodes: [[x,y,z], …] | Float64Array(3n),
 *   materials?: { name: { E, G?, nu?, rho?, alpha?, fy? } }, sections?: { name: section },
 *   members: [{ a, b, section, material?, type?: 'beam'|'truss'|'tension'|'compression',
 *               releases?: { a?, b? } (6 flags, 'pin', or 'ry,rz' …), betaDeg?, localZ?, offsetA?, offsetB?,
 *               shear?: boolean, flooded?, fluidRho?, extraMass? (kg/m), Cm?, group?, name? }],
 *   supports: [{ node, fix?: 'fixed'|'pinned'|[6 flags], spring?: [6] | 6×6, prescribed?: [6] }],
 *   nodeMasses?: [{ node, m, inertia?: [Ixx, Iyy, Izz] }],
 *   loadCases?: [{ name, nodal?, member?, selfWeight?, acceleration?, settlements? }],
 *   combinations?: [{ name, factors: { caseName: factor } }] }
 * ```
 * @param {object} model frame model (SI units)
 * @param {object} [opts]
 * @param {boolean} [opts.shearDeformation=true] include Timoshenko shear flexibility where the section has shear areas
 * @param {boolean} [opts.renumber=true] apply reverse Cuthill–McKee renumbering
 * @returns {object} system: { model, nNode, nDof, xyz, elems, fixed, free, nFree, eqDof, springs, autoRestrained, warnings }
 */
export function buildFrameSystem(model, opts = {}) {
  const shearDefault = opts.shearDeformation !== false;
  if (!model || typeof model !== 'object') throw new RangeError('model must be an object');
  const xyz = nodeArray(model.nodes);
  const nNode = xyz.length / 3;
  const nDof = 6 * nNode;
  if (!Array.isArray(model.members) || model.members.length === 0) throw new RangeError('model.members must not be empty');
  const warnings = [];
  const elems = [];
  const tmp = new Float64Array(144);

  for (let m = 0; m < model.members.length; m++) {
    const mem = model.members[m];
    const ia = mem.a; const ib = mem.b;
    if (!Number.isInteger(ia) || !Number.isInteger(ib) || ia < 0 || ib < 0 || ia >= nNode || ib >= nNode || ia === ib) {
      throw new RangeError(`model.members[${m}] has invalid node indices`);
    }
    const mat = resolveMaterial(model, mem.material);
    const sec = resolveSec(model, mem.section);
    const oa = mem.offsetA || [0, 0, 0];
    const ob = mem.offsetB || [0, 0, 0];
    const pa = [xyz[3 * ia] + oa[0], xyz[3 * ia + 1] + oa[1], xyz[3 * ia + 2] + oa[2]];
    const pb = [xyz[3 * ib] + ob[0], xyz[3 * ib + 1] + ob[1], xyz[3 * ib + 2] + ob[2]];
    let axes;
    try {
      axes = memberAxes(pa, pb, { betaDeg: mem.betaDeg || 0, localZ: mem.localZ });
    } catch (e) {
      throw new RangeError(`model.members[${m}]: ${e.message}`);
    }
    const { R, L } = axes;
    const type = mem.type || 'beam';
    if (!['beam', 'truss', 'tension', 'compression'].includes(type)) {
      throw new RangeError(`model.members[${m}].type "${type}" is not recognised`);
    }
    const truss = type !== 'beam';
    const useShear = mem.shear === undefined ? shearDefault : !!mem.shear;
    const Asy = useShear && sec.Asy > 0 ? sec.Asy : Infinity;
    const Asz = useShear && sec.Asz > 0 ? sec.Asz : Infinity;
    const kl = beamLocalStiffness({ E: mat.E, G: mat.G, A: sec.A, Iy: sec.Iy, Iz: sec.Iz, J: sec.J, L, Asy, Asz, truss });
    const phiY = Number.isFinite(Asy) ? 12 * mat.E * sec.Iz / (mat.G * Asy * L * L) : 0;
    const phiZ = Number.isFinite(Asz) ? 12 * mat.E * sec.Iy / (mat.G * Asz * L * L) : 0;

    // T maps the 12 global node dofs to the local dofs of the (offset) member ends.
    const T = new Float64Array(144);
    for (let e = 0; e < 2; e++) {
      const r = e === 0 ? oa : ob;
      const o = 6 * e;
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          T[(o + i) * 12 + o + j] = R[i][j];
          T[(o + 3 + i) * 12 + o + 3 + j] = R[i][j];
        }
        // u_end = u + θ × r  →  local = R·(−[r]×)·θ
        const Ri = R[i];
        T[(o + i) * 12 + o + 3] = -Ri[1] * r[2] + Ri[2] * r[1];
        T[(o + i) * 12 + o + 4] = Ri[0] * r[2] - Ri[2] * r[0];
        T[(o + i) * 12 + o + 5] = -Ri[0] * r[1] + Ri[1] * r[0];
      }
    }

    // Static condensation of released local dofs.
    const rel = new Uint8Array(12);
    if (mem.releases) {
      rel.set(releaseFlags(mem.releases.a), 0);
      rel.set(releaseFlags(mem.releases.b), 6);
    }
    const relIdx = [];
    for (let i = 0; i < 12; i++) if (rel[i] && kl[i * 12 + i] > 0) relIdx.push(i);
    let Tc = null;
    let krrInv = null;
    if (relIdx.length) {
      const nr = relIdx.length;
      const krr = relIdx.map((i) => relIdx.map((j) => kl[i * 12 + j]));
      krrInv = invSmall(krr, nr);
      if (!krrInv) {
        const names = relIdx.map((i) => `${DOF_NAMES[i % 6]}@${i < 6 ? 'A' : 'B'}`).join(', ');
        const err = new Error(`Member ${m} is a mechanism on its own: releases ${names} leave a rigid-body mode unrestrained.`);
        err.mechanism = true;
        throw err;
      }
      Tc = new Float64Array(144);
      for (let i = 0; i < 12; i++) if (!rel[i] || !(kl[i * 12 + i] > 0)) Tc[i * 12 + i] = 1;
      for (let a = 0; a < nr; a++) {
        for (let j = 0; j < 12; j++) {
          if (relIdx.includes(j)) continue;
          let s = 0;
          for (let b = 0; b < nr; b++) s += krrInv[a][b] * kl[relIdx[b] * 12 + j];
          Tc[relIdx[a] * 12 + j] = -s;
        }
      }
    }
    // TT = Tc·T maps global node dofs to full local dofs (homogeneous part).
    const TT = Tc ? mul12(Tc, T, new Float64Array(144)) : T;
    const ke = tAt12(TT, kl, new Float64Array(144), tmp);
    for (let i = 0; i < 12; i++) for (let j = i + 1; j < 12; j++) { const s = 0.5 * (ke[i * 12 + j] + ke[j * 12 + i]); ke[i * 12 + j] = s; ke[j * 12 + i] = s; }
    const dofs = new Int32Array(12);
    for (let i = 0; i < 6; i++) { dofs[i] = 6 * ia + i; dofs[6 + i] = 6 * ib + i; }
    elems.push({
      index: m, a: ia, b: ib, pa, pb, L, R, T, Tc, TT, kl, ke, relIdx, krrInv, dofs, mat, sec, type, truss, phiY, phiZ,
      nonlinear: type === 'tension' || type === 'compression', member: mem,
    });
  }

  // Supports and springs.
  const fixed = new Uint8Array(nDof);
  const prescribed = new Float64Array(nDof);
  const springs = [];
  for (const [si, sup] of (model.supports || []).entries()) {
    const n = sup.node;
    if (!Number.isInteger(n) || n < 0 || n >= nNode) throw new RangeError(`model.supports[${si}].node is invalid`);
    const f = fixFlags(sup.fix);
    for (let i = 0; i < 6; i++) if (f[i]) fixed[6 * n + i] = 1;
    if (sup.prescribed) {
      for (let i = 0; i < 6; i++) {
        const v = sup.prescribed[i] || 0;
        if (v !== 0 && !f[i]) throw new RangeError(`model.supports[${si}].prescribed[${i}] refers to a dof that is not fixed`);
        prescribed[6 * n + i] += v;
      }
    }
    if (sup.spring) {
      const k = new Float64Array(36);
      if (Array.isArray(sup.spring[0]) || ArrayBuffer.isView(sup.spring[0])) {
        for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) k[i * 6 + j] = 0.5 * (sup.spring[i][j] + sup.spring[j][i]);
      } else {
        for (let i = 0; i < 6; i++) {
          const v = sup.spring[i] || 0;
          if (v < 0) throw new RangeError(`model.supports[${si}].spring[${i}] must not be negative`);
          k[i * 6 + i] = v;
        }
      }
      springs.push({ node: n, k });
    }
  }

  // Dofs without any stiffness (e.g. rotations at pure truss nodes) are restrained automatically.
  const diag = new Float64Array(nDof);
  for (const e of elems) for (let i = 0; i < 12; i++) diag[e.dofs[i]] += e.ke[i * 12 + i];
  for (const s of springs) for (let i = 0; i < 6; i++) diag[6 * s.node + i] += s.k[i * 6 + i];
  let dmax = 0;
  for (let i = 0; i < nDof; i++) dmax = Math.max(dmax, diag[i]);
  const autoRestrained = [];
  const auto = new Uint8Array(nDof);
  for (let i = 0; i < nDof; i++) {
    if (!fixed[i] && !(diag[i] > 1e-13 * dmax)) { auto[i] = 1; autoRestrained.push(i); }
  }

  // Equation numbering (reverse Cuthill–McKee on the node graph).
  let order = null;
  if (opts.renumber !== false && nNode > 2) {
    const edges = elems.map((e) => [e.a, e.b]);
    try { order = reverseCuthillMcKee(nNode, edges); } catch (e) { order = null; }
    if (order && order.length !== nNode) order = null;
  }
  const free = new Int32Array(nDof).fill(-1);
  let nFree = 0;
  for (let q = 0; q < nNode; q++) {
    const n = order ? order[q] : q;
    for (let i = 0; i < 6; i++) {
      const d = 6 * n + i;
      if (!fixed[d] && !auto[d]) free[d] = nFree++;
    }
  }
  const eqDof = new Int32Array(nFree);
  for (let d = 0; d < nDof; d++) if (free[d] >= 0) eqDof[free[d]] = d;

  // Skyline profile: height[j] = j − (lowest coupled equation).
  const height = new Int32Array(nFree);
  const profile = (dofList) => {
    let lo = Infinity;
    for (let i = 0; i < dofList.length; i++) { const q = free[dofList[i]]; if (q >= 0 && q < lo) lo = q; }
    for (let i = 0; i < dofList.length; i++) { const q = free[dofList[i]]; if (q >= 0 && q - lo > height[q]) height[q] = q - lo; }
  };
  for (const e of elems) profile(e.dofs);
  for (const s of springs) profile([0, 1, 2, 3, 4, 5].map((i) => 6 * s.node + i));

  if (autoRestrained.length) {
    warnings.push(`${autoRestrained.length} dof(s) without stiffness were restrained automatically (e.g. rotations at truss nodes).`);
  }
  return {
    model, nNode, nDof, xyz, elems, fixed, auto, prescribed, free, nFree, eqDof, springs, height, autoRestrained, warnings,
    shearDeformation: shearDefault,
  };
}

/**
 * Human-readable name of a global dof, e.g. "node 4 rz".
 * @param {number} d global dof index (6·node + component)
 * @returns {string}
 */
export function dofName(d) {
  return `node ${Math.floor(d / 6)} ${DOF_NAMES[d % 6]}`;
}

/**
 * Transform a local 12×12 element matrix (stiffness-like: geometric stiffness, mass) to global node dofs with
 * the element's release condensation and rigid offsets: Kg = (Tc T)ᵀ m (Tc T).
 * @param {object} elem element record from buildFrameSystem
 * @param {Float64Array} mLocal row-major 12×12 local matrix
 * @returns {Float64Array} row-major 12×12 global matrix
 */
export function elementToGlobal(elem, mLocal) {
  const out = tAt12(elem.TT, mLocal, new Float64Array(144), new Float64Array(144));
  for (let i = 0; i < 12; i++) for (let j = i + 1; j < 12; j++) { const s = 0.5 * (out[i * 12 + j] + out[j * 12 + i]); out[i * 12 + j] = s; out[j * 12 + i] = s; }
  return out;
}

/**
 * y = Σ_e A_e x restricted to free equations, for element matrices given in global node dofs.
 * @param {object} sys system from buildFrameSystem
 * @param {Array<Float64Array|null>} mats one row-major 12×12 matrix per element (null entries skipped)
 * @param {Float64Array} x vector over free equations
 * @param {Float64Array} [y] output (zeroed first); allocated when omitted
 * @returns {Float64Array} y
 */
export function elementMatVec(sys, mats, x, y) {
  const out = y || new Float64Array(sys.nFree);
  out.fill(0);
  const { free, elems } = sys;
  const xe = new Float64Array(12);
  for (let e = 0; e < elems.length; e++) {
    const A = mats[e];
    if (!A) continue;
    const dofs = elems[e].dofs;
    for (let i = 0; i < 12; i++) { const q = free[dofs[i]]; xe[i] = q >= 0 ? x[q] : 0; }
    for (let i = 0; i < 12; i++) {
      const q = free[dofs[i]];
      if (q < 0) continue;
      let s = 0;
      for (let k = 0; k < 12; k++) s += A[i * 12 + k] * xe[k];
      out[q] += s;
    }
  }
  return out;
}

/**
 * Dense matrix over the free equations assembled from element matrices (and, optionally, the nodal springs).
 * Intended for small systems handed to dense eigen / time-integration routines.
 * @param {object} sys system from buildFrameSystem
 * @param {Array<Float64Array|null>} mats element matrices in global node dofs (row-major 12×12)
 * @param {boolean} [withSprings=false] add the support spring matrices
 * @returns {number[][]} nFree × nFree matrix
 */
export function assembleDense(sys, mats, withSprings = false) {
  const n = sys.nFree;
  const K = Array.from({ length: n }, () => new Array(n).fill(0));
  const { free, elems } = sys;
  for (let e = 0; e < elems.length; e++) {
    const A = mats[e];
    if (!A) continue;
    const dofs = elems[e].dofs;
    for (let i = 0; i < 12; i++) {
      const qi = free[dofs[i]];
      if (qi < 0) continue;
      for (let j = 0; j < 12; j++) {
        const qj = free[dofs[j]];
        if (qj >= 0) K[qi][qj] += A[i * 12 + j];
      }
    }
  }
  if (withSprings) {
    for (const s of sys.springs) {
      for (let i = 0; i < 6; i++) {
        const qi = free[6 * s.node + i];
        if (qi < 0) continue;
        for (let j = 0; j < 6; j++) {
          const qj = free[6 * s.node + j];
          if (qj >= 0) K[qi][qj] += s.k[i * 6 + j];
        }
      }
    }
  }
  return K;
}

/**
 * Assemble and factorise K (+ optional additional element matrices, e.g. geometric stiffness) with the
 * skyline LDLᵀ solver. A zero or negligible pivot is reported as a mechanism naming the offending dof.
 * @param {object} sys system from buildFrameSystem
 * @param {object} [o]
 * @param {Uint8Array} [o.active] per-element activity flags (inactive elements contribute no stiffness)
 * @param {Array<Float64Array|null>} [o.add] extra element matrices in global node dofs added to K
 * @param {number} [o.addFactor=1] multiplier on `add`
 * @returns {{solve:(b:Float64Array)=>Float64Array, n:number, negativePivots:number}} factorised operator on the free equations;
 *   negativePivots > 0 signals an indefinite (unstable) tangent stiffness
 */
export function factorFrame(sys, { active, add, addFactor = 1 } = {}) {
  const n = sys.nFree;
  if (n === 0) return { n: 0, solve: () => new Float64Array(0), negativePivots: 0 };
  const sky = new SkylineSym(n, sys.height);
  const { free, elems } = sys;
  for (let e = 0; e < elems.length; e++) {
    const el = elems[e];
    const on = !active || active[e];
    const A = add ? add[e] : null;
    if (!on && !A) continue;
    for (let i = 0; i < 12; i++) {
      const qi = free[el.dofs[i]];
      if (qi < 0) continue;
      for (let j = i; j < 12; j++) {
        const qj = free[el.dofs[j]];
        if (qj < 0) continue;
        let v = on ? el.ke[i * 12 + j] : 0;
        if (A) v += addFactor * A[i * 12 + j];
        if (v === 0) continue;
        if (qi <= qj) sky.add(qi, qj, v); else sky.add(qj, qi, v);
      }
    }
  }
  for (const s of sys.springs) {
    for (let i = 0; i < 6; i++) {
      const qi = free[6 * s.node + i];
      if (qi < 0) continue;
      for (let j = i; j < 6; j++) {
        const qj = free[6 * s.node + j];
        if (qj < 0) continue;
        const v = s.k[i * 6 + j];
        if (v === 0) continue;
        if (qi <= qj) sky.add(qi, qj, v); else sky.add(qj, qi, v);
      }
    }
  }
  try {
    sky.factor();
  } catch (err) {
    let q = Number.isInteger(err && err.dof) ? err.dof : (Number.isInteger(err && err.index) ? err.index : -1);
    if (q < 0) {
      const mt = /(\d+)/.exec(String(err && err.message));
      if (mt) q = Number(mt[1]);
    }
    const where = q >= 0 && q < n ? dofName(sys.eqDof[q]) : 'unknown dof';
    const e2 = new Error(`Structure is a mechanism or unstable: zero pivot at ${where}.`);
    e2.mechanism = true; e2.dof = q >= 0 && q < n ? sys.eqDof[q] : -1;
    throw e2;
  }
  return { n, solve: (b) => sky.solve(b), negativePivots: sky.negativePivots };
}

/* ------------------------------------------------------------------ loads */

function vec3(name, v) {
  if (!v || v.length < 3 || ![v[0], v[1], v[2]].every(Number.isFinite)) throw new RangeError(`${name} must be a vector of three finite numbers`);
  return v;
}

function toLocal(R, v) {
  return [
    R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
    R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
    R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2],
  ];
}

function spanOf(el, ld) {
  let s1 = 0; let s2 = el.L;
  if (ld.s1 !== undefined) s1 = ld.s1; else if (ld.xi1 !== undefined) s1 = ld.xi1 * el.L;
  if (ld.s2 !== undefined) s2 = ld.s2; else if (ld.xi2 !== undefined) s2 = ld.xi2 * el.L;
  if (!(s1 >= -1e-9 * el.L && s2 <= el.L * (1 + 1e-9) && s2 > s1)) {
    throw new RangeError(`member load span [${s1}, ${s2}] lies outside member ${el.index} (L = ${el.L})`);
  }
  return [Math.max(0, s1), Math.min(el.L, s2)];
}

/**
 * Mass per unit length of a member used for self-weight, inertia loads and mass matrices: steel + contained
 * fluid (flooded tubes or `fluidRho`) + non-structural `extraMass` + marine growth. Hydrodynamic added mass
 * is not included here (see framemodal.js).
 * @param {object} elem element record
 * @param {object} [o]
 * @param {number} [o.rhoW=1025] sea-water density used for flooded members [kg/m³]
 * @param {number} [o.mgThickness=0] marine growth thickness [m]
 * @param {number} [o.mgRho=1325] marine growth density [kg/m³]
 * @returns {{steel:number, fluid:number, extra:number, growth:number, total:number}} kg/m
 */
export function memberMassPerLength(elem, { rhoW = 1025, mgThickness = 0, mgRho = 1325 } = {}) {
  const sec = elem.sec; const mem = elem.member;
  const steel = (sec.massPerLength !== undefined) ? sec.massPerLength : elem.mat.rho * sec.A;
  let fluid = 0;
  if (sec.type === 'tube') {
    const di = sec.D - 2 * sec.t;
    const rf = mem.fluidRho !== undefined ? mem.fluidRho : (mem.flooded ? rhoW : 0);
    fluid = rf * Math.PI * di * di / 4;
  }
  const extra = mem.extraMass || 0;
  let growth = 0;
  if (mgThickness > 0 && sec.D) growth = mgRho * Math.PI / 4 * ((sec.D + 2 * mgThickness) ** 2 - sec.D ** 2);
  return { steel, fluid, extra, growth, total: steel + fluid + extra + growth };
}

function pushItem(items, e, it) { items[e].push(it); }

/**
 * Resolve a load case (or a factored set of load cases) into nodal loads, per-member local load items and
 * prescribed displacements.
 *
 * Load case format
 * ```
 * { name,
 *   nodal: [{ node, F: [Fx, Fy, Fz, Mx, My, Mz] }],
 *   member: [
 *     { member, type: 'uniform', q: [qx,qy,qz], axes?: 'global'|'local', s1?, s2? (or xi1, xi2 fractions) },
 *     { member, type: 'linear', q1, q2, axes?, s1?, s2? },
 *     { member, type: 'point', s | xi, P?: [3], M?: [3], axes? },
 *     { member, type: 'thermal', dT? [K], gradY? [K/m], gradZ? [K/m] },
 *     { member, type: 'function', fn: (s, [x,y,z]) => [fx,fy,fz] (N/m), axes?, nSeg? } ],
 *   selfWeight: { g?: [0,0,-9.80665], factor?: 1, buoyancy?: { waterline (z, m), rhoW?: 1025 },
 *                 marineGrowth?: { thickness, rho?: 1325, zTop, zBottom } },
 *   acceleration: { a?: [3] (m/s²), alpha?: [3] (rad/s²), omega?: [3] (rad/s), origin?: [3] },
 *   settlements: [{ node, d: [6] }] }
 * ```
 * Self-weight uses the "marine" (buoyant-weight) method: buoyancy is a distributed force opposite to g equal
 * to the weight of displaced water (outer envelope for unflooded tubes, steel volume for flooded tubes);
 * pressure end-cap effects are not modelled. The acceleration field applies d'Alembert loads −m·a(r) with
 * a(r) = a + α × r + ω × (ω × r), r measured from `origin` (transport / seafastening analysis).
 * @param {object} sys system from buildFrameSystem
 * @param {Array<{load:object, factor:number}>} parts load cases with factors
 * @returns {{F:Float64Array, items:Array<Array<object>>, prescribed:Float64Array}} nodal loads over all dofs [N, N·m], load items, prescribed displacements [m, rad]
 */
export function resolveLoads(sys, parts) {
  const { elems, nDof, nNode } = sys;
  const F = new Float64Array(nDof);
  const prescribed = Float64Array.from(sys.prescribed);
  const items = elems.map(() => []);
  for (const { load, factor } of parts) {
    if (!load) throw new RangeError('load case is undefined');
    const fac = factor === undefined ? 1 : factor;
    for (const [i, nl] of (load.nodal || []).entries()) {
      if (!Number.isInteger(nl.node) || nl.node < 0 || nl.node >= nNode) throw new RangeError(`nodal load ${i}: invalid node`);
      for (let k = 0; k < 6; k++) F[6 * nl.node + k] += fac * (nl.F[k] || 0);
    }
    for (const [i, ml] of (load.member || []).entries()) {
      const e = ml.member;
      if (!Number.isInteger(e) || e < 0 || e >= elems.length) throw new RangeError(`member load ${i}: invalid member`);
      const el = elems[e];
      const local = ml.axes === 'local';
      const cv = (v) => (local ? [v[0], v[1], v[2]] : toLocal(el.R, v));
      switch (ml.type) {
        case 'uniform': {
          const q = cv(vec3('q', ml.q));
          const [s1, s2] = spanOf(el, ml);
          pushItem(items, e, { t: 'd', s1, s2, fac, f: () => q, smooth: true });
          break;
        }
        case 'linear': {
          const q1 = cv(vec3('q1', ml.q1)); const q2 = cv(vec3('q2', ml.q2));
          const [s1, s2] = spanOf(el, ml);
          const len = s2 - s1;
          pushItem(items, e, {
            t: 'd', s1, s2, fac, smooth: true,
            f: (s) => { const w = (s - s1) / len; return [q1[0] + w * (q2[0] - q1[0]), q1[1] + w * (q2[1] - q1[1]), q1[2] + w * (q2[2] - q1[2])]; },
          });
          break;
        }
        case 'point': {
          const s = ml.s !== undefined ? ml.s : (ml.xi !== undefined ? ml.xi * el.L : NaN);
          if (!(s >= 0 && s <= el.L * (1 + 1e-12))) throw new RangeError(`member load ${i}: point position outside the member`);
          pushItem(items, e, { t: 'p', s: Math.min(s, el.L), fac, P: cv(ml.P || [0, 0, 0]), M: cv(ml.M || [0, 0, 0]) });
          break;
        }
        case 'thermal':
          pushItem(items, e, { t: 'th', fac, dT: ml.dT || 0, gy: ml.gradY || 0, gz: ml.gradZ || 0 });
          break;
        case 'function': {
          if (typeof ml.fn !== 'function') throw new RangeError(`member load ${i}: fn must be a function`);
          const pa = el.pa; const ex = el.R[0];
          const f = (s) => {
            const v = ml.fn(s, [pa[0] + ex[0] * s, pa[1] + ex[1] * s, pa[2] + ex[2] * s]);
            return local ? v : toLocal(el.R, v);
          };
          const [s1, s2] = spanOf(el, ml);
          // sampled at Gauss points of nSeg sub-intervals so that kinks (e.g. at the free surface) are resolved
          const nSeg = Math.max(1, Math.floor(ml.nSeg || 4));
          const h = (s2 - s1) / nSeg;
          for (let k = 0; k < nSeg; k++) pushItem(items, e, { t: 'd', s1: s1 + k * h, s2: s1 + (k + 1) * h, fac, f, smooth: false });
          break;
        }
        default:
          throw new RangeError(`member load ${i}: type "${ml.type}" is not recognised`);
      }
    }
    if (load.selfWeight) addSelfWeight(sys, load.selfWeight, fac, F, items);
    if (load.acceleration) addAcceleration(sys, load.acceleration, fac, F, items);
    for (const [i, st] of (load.settlements || []).entries()) {
      if (!Number.isInteger(st.node) || st.node < 0 || st.node >= nNode) throw new RangeError(`settlement ${i}: invalid node`);
      for (let k = 0; k < 6; k++) {
        const v = st.d[k] || 0;
        if (v !== 0 && !sys.fixed[6 * st.node + k]) throw new RangeError(`settlement ${i}: dof ${DOF_NAMES[k]} of node ${st.node} is not restrained`);
        prescribed[6 * st.node + k] += fac * v;
      }
    }
  }
  return { F, items, prescribed };
}

function addSelfWeight(sys, sw, fac, F, items) {
  const g = sw.g || [0, 0, -G0];
  vec3('selfWeight.g', g);
  const f0 = fac * (sw.factor === undefined ? 1 : sw.factor);
  const buoy = sw.buoyancy;
  const mg = sw.marineGrowth;
  const rhoW = buoy && buoy.rhoW !== undefined ? buoy.rhoW : 1025;
  if (buoy && !Number.isFinite(buoy.waterline)) throw new RangeError('selfWeight.buoyancy.waterline must be a finite elevation');
  if (mg && !(mg.thickness >= 0)) throw new RangeError('selfWeight.marineGrowth.thickness must be non-negative');
  for (let e = 0; e < sys.elems.length; e++) {
    const el = sys.elems[e];
    const za = el.pa[2]; const zb = el.pb[2];
    const cuts = [0, el.L];
    const addCut = (zc) => {
      if (zc === undefined || zb === za) return;
      const s = (zc - za) / (zb - za) * el.L;
      if (s > 1e-9 * el.L && s < el.L * (1 - 1e-9)) cuts.push(s);
    };
    if (buoy) addCut(buoy.waterline);
    if (mg) { addCut(mg.zTop); addCut(mg.zBottom); }
    cuts.sort((p, q) => p - q);
    const sec = el.sec;
    for (let c = 0; c + 1 < cuts.length; c++) {
      const s1 = cuts[c]; const s2 = cuts[c + 1];
      const zm = za + (zb - za) * (0.5 * (s1 + s2)) / el.L;
      const tg = mg && sec.D && zm <= (mg.zTop === undefined ? Infinity : mg.zTop) && zm >= (mg.zBottom === undefined ? -Infinity : mg.zBottom)
        ? mg.thickness : 0;
      const mass = memberMassPerLength(el, { rhoW, mgThickness: tg, mgRho: mg && mg.rho !== undefined ? mg.rho : 1325 });
      const wet = buoy && zm < buoy.waterline;
      let total = mass.total;
      // a flooded member holds water only below the waterline
      if (buoy && !wet && el.member.flooded && el.member.fluidRho === undefined) total -= mass.fluid;
      let displaced = 0; // m² of displaced water per unit length (outer envelope incl. marine growth)
      if (wet) displaced = sec.type === 'tube' ? Math.PI * (sec.D + 2 * tg) ** 2 / 4 : sec.A;
      const w = total - rhoW * displaced; // net mass per length acting along g
      const q = toLocal(el.R, [w * g[0], w * g[1], w * g[2]]);
      items[e].push({ t: 'd', s1, s2, fac: f0, f: () => q, smooth: true });
    }
  }
  for (const nm of sys.model.nodeMasses || []) {
    for (let k = 0; k < 3; k++) F[6 * nm.node + k] += f0 * nm.m * g[k];
  }
}

function accelAt(acc, p) {
  const a = acc.a || [0, 0, 0];
  const al = acc.alpha || [0, 0, 0];
  const om = acc.omega || [0, 0, 0];
  const o = acc.origin || [0, 0, 0];
  const r = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
  const axr = [al[1] * r[2] - al[2] * r[1], al[2] * r[0] - al[0] * r[2], al[0] * r[1] - al[1] * r[0]];
  const oxr = [om[1] * r[2] - om[2] * r[1], om[2] * r[0] - om[0] * r[2], om[0] * r[1] - om[1] * r[0]];
  const oxoxr = [om[1] * oxr[2] - om[2] * oxr[1], om[2] * oxr[0] - om[0] * oxr[2], om[0] * oxr[1] - om[1] * oxr[0]];
  return [a[0] + axr[0] + oxoxr[0], a[1] + axr[1] + oxoxr[1], a[2] + axr[2] + oxoxr[2]];
}

function addAcceleration(sys, acc, fac, F, items) {
  for (let e = 0; e < sys.elems.length; e++) {
    const el = sys.elems[e];
    const mu = memberMassPerLength(el).total;
    const pa = el.pa; const ex = el.R[0];
    const f = (s) => {
      const a = accelAt(acc, [pa[0] + ex[0] * s, pa[1] + ex[1] * s, pa[2] + ex[2] * s]);
      return toLocal(el.R, [-mu * a[0], -mu * a[1], -mu * a[2]]);
    };
    items[e].push({ t: 'd', s1: 0, s2: el.L, fac, f, smooth: true });
  }
  for (const nm of sys.model.nodeMasses || []) {
    const n = nm.node;
    const a = accelAt(acc, [sys.xyz[3 * n], sys.xyz[3 * n + 1], sys.xyz[3 * n + 2]]);
    for (let k = 0; k < 3; k++) F[6 * n + k] -= fac * nm.m * a[k];
    if (nm.inertia && acc.alpha) for (let k = 0; k < 3; k++) F[6 * n + 3 + k] -= fac * nm.inertia[k] * acc.alpha[k];
  }
}

/* Consistent local nodal load vector (uncondensed, 12) of the load items of one element. */
function elementLoadVector(el, its, out) {
  out.fill(0);
  if (!its.length) return out;
  const L = el.L;
  const gq = gauss8();
  const sy = new Float64Array(8); const sz = new Float64Array(8);
  for (const it of its) {
    const fac = it.fac;
    if (it.t === 'd') {
      const half = 0.5 * (it.s2 - it.s1); const mid = 0.5 * (it.s2 + it.s1);
      for (let g = 0; g < gq.x.length; g++) {
        const s = mid + half * gq.x[g];
        const w = gq.w[g] * half * fac;
        const q = it.f(s);
        const xi = s / L;
        out[0] += w * q[0] * (1 - xi); out[6] += w * q[0] * xi;
        bendShape(xi, L, el.phiY, el.truss, sy);
        bendShape(xi, L, el.phiZ, el.truss, sz);
        out[1] += w * q[1] * sy[0]; out[5] += w * q[1] * sy[1]; out[7] += w * q[1] * sy[2]; out[11] += w * q[1] * sy[3];
        out[2] += w * q[2] * sz[0]; out[4] -= w * q[2] * sz[1]; out[8] += w * q[2] * sz[2]; out[10] -= w * q[2] * sz[3];
      }
    } else if (it.t === 'p') {
      const xi = it.s / L;
      bendShape(xi, L, el.phiY, el.truss, sy);
      bendShape(xi, L, el.phiZ, el.truss, sz);
      const P = it.P; const M = it.M;
      out[0] += fac * P[0] * (1 - xi); out[6] += fac * P[0] * xi;
      out[3] += fac * M[0] * (1 - xi); out[9] += fac * M[0] * xi;
      out[1] += fac * (P[1] * sy[0] + M[2] * sy[4]); out[5] += fac * (P[1] * sy[1] + M[2] * sy[5]);
      out[7] += fac * (P[1] * sy[2] + M[2] * sy[6]); out[11] += fac * (P[1] * sy[3] + M[2] * sy[7]);
      out[2] += fac * (P[2] * sz[0] - M[1] * sz[4]); out[4] += fac * (-P[2] * sz[1] + M[1] * sz[5]);
      out[8] += fac * (P[2] * sz[2] - M[1] * sz[6]); out[10] += fac * (-P[2] * sz[3] + M[1] * sz[7]);
    } else if (it.t === 'th') {
      const ea = el.mat.E * el.sec.A * el.mat.alpha * it.dT * fac;
      out[0] -= ea; out[6] += ea;
      if (!el.truss) {
        const my = el.mat.E * el.sec.Iy * el.mat.alpha * it.gz * fac;
        const mz = el.mat.E * el.sec.Iz * el.mat.alpha * it.gy * fac;
        out[4] -= my; out[10] += my;
        out[5] += mz; out[11] -= mz;
      }
    }
  }
  return out;
}

/* ------------------------------------------------------------------ solution */

/**
 * Solve one resolved load set on a factorised system and recover displacements, reactions, member end
 * forces, equilibrium residuals and strain energy.
 * @param {object} sys system from buildFrameSystem
 * @param {{solve:Function}} fac factorisation from factorFrame
 * @param {{F:Float64Array, items:Array, prescribed:Float64Array}} loads resolved loads
 * @param {object} [o]
 * @param {Uint8Array} [o.active] per-element activity flags
 * @param {Array<Float64Array|null>} [o.add] extra element matrices included in K (for residual evaluation)
 * @param {number} [o.addFactor=1] multiplier on `add` and `addLocal`
 * @param {Array<Float64Array|null>} [o.addLocal] local 12×12 matrices matching `add` (e.g. geometric stiffness) so that the member
 *   end forces include their contribution (second-order moments)
 * @returns {object} see analyzeFrame for the result fields
 */
export function solveFrame(sys, fac, loads, { active, add, addFactor = 1, addLocal } = {}) {
  const { elems, nDof, nFree, free, fixed } = sys;
  const Ftot = Float64Array.from(loads.F);
  const u = new Float64Array(nDof);
  for (let d = 0; d < nDof; d++) if (fixed[d]) u[d] = loads.prescribed[d];
  const fl = new Float64Array(12);
  const fg = new Float64Array(12);
  const flAll = new Array(elems.length);
  for (let e = 0; e < elems.length; e++) {
    const el = elems[e];
    elementLoadVector(el, loads.items[e], fl);
    flAll[e] = Float64Array.from(fl);
    matTVec12(el.TT, fl, fg);
    for (let i = 0; i < 12; i++) Ftot[el.dofs[i]] += fg[i];
  }
  for (const d of sys.autoRestrained) {
    if (Math.abs(Ftot[d]) > 0) {
      const err = new Error(`Load applied to ${dofName(d)}, which has no stiffness.`);
      err.mechanism = true; err.dof = d;
      throw err;
    }
  }
  // right-hand side with the contribution of prescribed displacements
  const b = new Float64Array(nFree);
  for (let d = 0; d < nDof; d++) if (free[d] >= 0) b[free[d]] = Ftot[d];
  const de = new Float64Array(12);
  const ke = new Float64Array(144);
  const elemK = (e) => {
    const el = elems[e];
    const on = !active || active[e];
    const A = add ? add[e] : null;
    if (!on && !A) return null;
    if (on && !A) return el.ke;
    for (let i = 0; i < 144; i++) ke[i] = (on ? el.ke[i] : 0) + (A ? addFactor * A[i] : 0);
    return ke;
  };
  let anyPres = false;
  for (let d = 0; d < nDof; d++) if (u[d] !== 0) { anyPres = true; break; }
  if (anyPres) {
    for (let e = 0; e < elems.length; e++) {
      const el = elems[e];
      let has = false;
      for (let i = 0; i < 12; i++) { de[i] = u[el.dofs[i]]; if (de[i] !== 0) has = true; }
      if (!has) continue;
      const K = elemK(e);
      if (!K) continue;
      matVec12(K, de, fg);
      for (let i = 0; i < 12; i++) { const q = free[el.dofs[i]]; if (q >= 0) b[q] -= fg[i]; }
    }
  }
  const x = fac.solve(b);
  let umax = 0;
  for (let q = 0; q < nFree; q++) {
    const v = x[q];
    if (!Number.isFinite(v)) {
      const err = new Error(`Structure is a mechanism or unstable: non-finite displacement at ${dofName(sys.eqDof[q])}.`);
      err.mechanism = true; err.dof = sys.eqDof[q];
      throw err;
    }
    u[sys.eqDof[q]] = v;
    if (Math.abs(v) > umax) umax = Math.abs(v);
  }
  // g = K u − F over all dofs → reactions at restrained dofs, residual at free dofs
  const gvec = new Float64Array(nDof);
  for (let d = 0; d < nDof; d++) gvec[d] = -Ftot[d];
  let energy = 0;
  let fscale = 0;
  for (let e = 0; e < elems.length; e++) {
    const el = elems[e];
    const K = elemK(e);
    if (!K) continue;
    for (let i = 0; i < 12; i++) de[i] = u[el.dofs[i]];
    matVec12(K, de, fg);
    for (let i = 0; i < 12; i++) {
      gvec[el.dofs[i]] += fg[i];
      energy += 0.5 * de[i] * fg[i];
      if (Math.abs(fg[i]) > fscale) fscale = Math.abs(fg[i]);
    }
  }
  const springForce = new Float64Array(nDof);
  for (const s of sys.springs) {
    for (let i = 0; i < 6; i++) {
      let f = 0;
      for (let j = 0; j < 6; j++) f += s.k[i * 6 + j] * u[6 * s.node + j];
      springForce[6 * s.node + i] -= f;
      energy += 0.5 * f * u[6 * s.node + i];
    }
  }
  let nodalResidual = 0; let worstDof = -1;
  for (let d = 0; d < nDof; d++) {
    if (free[d] < 0) continue;
    const r = Math.abs(gvec[d] - springForce[d]);
    if (r > nodalResidual) { nodalResidual = r; worstDof = d; }
    if (Math.abs(Ftot[d]) > fscale) fscale = Math.abs(Ftot[d]);
  }
  if (fscale > 0 && nodalResidual > 1e-5 * fscale) {
    const err = new Error(`Structure is a mechanism or ill-conditioned: equilibrium not satisfied at ${dofName(worstDof)} `
      + `(residual ${nodalResidual.toExponential(2)} N against force scale ${fscale.toExponential(2)} N).`);
    err.mechanism = true; err.dof = worstDof;
    throw err;
  }
  // reactions (support nodes: rigid restraint + spring force)
  const reactions = [];
  const sumR = [0, 0, 0, 0, 0, 0];
  const sumF = [0, 0, 0, 0, 0, 0];
  const seen = new Set();
  const supNodes = [];
  for (const s of sys.model.supports || []) if (!seen.has(s.node)) { seen.add(s.node); supNodes.push(s.node); }
  const xyz = sys.xyz;
  const addTo = (acc, n, f) => {
    const px = xyz[3 * n]; const py = xyz[3 * n + 1]; const pz = xyz[3 * n + 2];
    acc[0] += f[0]; acc[1] += f[1]; acc[2] += f[2];
    acc[3] += f[3] + py * f[2] - pz * f[1];
    acc[4] += f[4] + pz * f[0] - px * f[2];
    acc[5] += f[5] + px * f[1] - py * f[0];
  };
  for (const n of supNodes) {
    const Rn = [0, 0, 0, 0, 0, 0];
    for (let i = 0; i < 6; i++) {
      const d = 6 * n + i;
      Rn[i] = (fixed[d] ? gvec[d] : 0) + springForce[d];
    }
    reactions.push({ node: n, R: Rn });
    addTo(sumR, n, Rn);
  }
  for (let n = 0; n < sys.nNode; n++) {
    const f = [0, 0, 0, 0, 0, 0];
    let any = false;
    for (let i = 0; i < 6; i++) { f[i] = Ftot[6 * n + i]; if (f[i] !== 0) any = true; }
    if (any) addTo(sumF, n, f);
  }
  const forceResidual = Math.hypot(sumF[0] + sumR[0], sumF[1] + sumR[1], sumF[2] + sumR[2]);
  const momentResidual = Math.hypot(sumF[3] + sumR[3], sumF[4] + sumR[4], sumF[5] + sumR[5]);

  // member end forces in local axes: f = kl (TT d + d0) − f_eq
  const memberForces = new Array(elems.length);
  const dl = new Float64Array(12);
  const endDisp = new Array(elems.length);
  for (let e = 0; e < elems.length; e++) {
    const el = elems[e];
    const on = !active || active[e];
    for (let i = 0; i < 12; i++) de[i] = u[el.dofs[i]];
    matVec12(el.TT, de, dl);
    const feq = flAll[e];
    if (el.relIdx.length) {
      for (let a = 0; a < el.relIdx.length; a++) {
        let s = 0;
        for (let c = 0; c < el.relIdx.length; c++) s += el.krrInv[a][c] * feq[el.relIdx[c]];
        dl[el.relIdx[a]] += s;
      }
    }
    const f = new Float64Array(12);
    if (on) matVec12(el.kl, dl, f);
    if (addLocal && addLocal[e]) {
      matVec12(addLocal[e], dl, fg);
      for (let i = 0; i < 12; i++) f[i] += addFactor * fg[i];
    }
    for (let i = 0; i < 12; i++) f[i] -= feq[i];
    for (const i of el.relIdx) f[i] = 0;
    // local displacements of the member ends from the node displacements (independent of releases)
    const dn = new Float64Array(12);
    matVec12(el.T, de, dn);
    endDisp[e] = dn;
    memberForces[e] = {
      member: e, active: !!on, f,
      endA: { N: -f[0], Vy: -f[1], Vz: -f[2], T: -f[3], My: -f[4], Mz: -f[5] },
      endB: { N: f[6], Vy: f[7], Vz: f[8], T: f[9], My: f[10], Mz: f[11] },
    };
  }
  const result = {
    displacements: u, reactions, reactionSum: sumR, appliedSum: sumF, memberForces,
    equilibrium: {
      forceResidual, momentResidual, nodalResidual,
      relative: fscale > 0 ? Math.max(forceResidual, nodalResidual) / fscale : 0,
    },
    strainEnergy: energy, maxDisplacement: umax,
  };
  Object.defineProperty(result, 'internal', { value: { sys, items: loads.items, endDisp, Ftot }, enumerable: false });
  return result;
}

/**
 * Member end forces in local axes for a given displacement field without member loads (f = kl·TT·d), e.g. for
 * mode shapes or modal response vectors.
 * @param {object} sys system from buildFrameSystem
 * @param {ArrayLike<number>} u displacements over all 6·nNode dofs [m, rad]
 * @returns {Float64Array[]} one 12-vector per member: forces on the member at end A (0–5) and end B (6–11) [N, N·m]
 */
export function memberForcesFromDisplacements(sys, u) {
  if (!u || u.length !== sys.nDof) throw new RangeError('u must hold 6 values per node');
  const de = new Float64Array(12); const dl = new Float64Array(12);
  return sys.elems.map((el) => {
    for (let i = 0; i < 12; i++) de[i] = u[el.dofs[i]];
    matVec12(el.TT, de, dl);
    const f = matVec12(el.kl, dl, new Float64Array(12));
    for (const i of el.relIdx) f[i] = 0;
    return f;
  });
}

function findCase(model, name) {
  const lc = (model.loadCases || []).find((c) => c.name === name);
  if (lc) return [{ load: lc, factor: 1 }];
  const cb = (model.combinations || []).find((c) => c.name === name);
  if (!cb) throw new RangeError(`load case or combination "${name}" is not defined`);
  const parts = [];
  for (const [k, f] of Object.entries(cb.factors || {})) {
    const sub = (model.loadCases || []).find((c) => c.name === k);
    if (!sub) throw new RangeError(`combination "${name}" refers to unknown load case "${k}"`);
    parts.push({ load: sub, factor: f });
  }
  return parts;
}

/**
 * Load-case parts ({load, factor}) for a named load case or combination of the model.
 * @param {object} model frame model
 * @param {string} name load case or combination name
 * @returns {Array<{load:object, factor:number}>}
 */
export function loadParts(model, name) {
  return findCase(model, name);
}

/**
 * Linear elastic analysis of a 3-D frame for all (or selected) load cases and combinations. Combinations are
 * solved as their own load sets so that tension-only / compression-only members are treated consistently.
 * Reference: direct stiffness method (McGuire, Gallagher & Ziemian, Matrix Structural Analysis, 2nd ed.).
 * @param {object} model frame model (see buildFrameSystem and resolveLoads for the format; SI units)
 * @param {object} [opts]
 * @param {boolean} [opts.shearDeformation=true] Timoshenko shear flexibility
 * @param {string[]} [opts.cases] names to solve (default: every load case and combination)
 * @param {number} [opts.maxIter=50] iteration limit for tension/compression-only members
 * @returns {{system:object, order:string[], results:Object<string,object>, warnings:string[]}} results[name] holds
 *   displacements (Float64Array 6·nNode: m, rad), reactions [{node, R[6]}] (N, N·m), reactionSum / appliedSum (resultants about the
 *   origin), memberForces [{member, active, f (local end forces on the member, 12), endA, endB {N, Vy, Vz, T, My, Mz}}],
 *   equilibrium {forceResidual [N], momentResidual [N·m], nodalResidual [N], relative [-]}, strainEnergy [J],
 *   iterations, inactiveMembers
 */
export function analyzeFrame(model, opts = {}) {
  const sys = buildFrameSystem(model, opts);
  const names = opts.cases || [
    ...(model.loadCases || []).map((c) => c.name), ...(model.combinations || []).map((c) => c.name),
  ];
  if (names.length === 0) throw new RangeError('model.loadCases must define at least one load case');
  const hasNl = sys.elems.some((e) => e.nonlinear);
  const results = {};
  let linearFac = null;
  for (const name of names) {
    const loads = resolveLoads(sys, findCase(model, name));
    let res;
    if (!hasNl) {
      if (!linearFac) linearFac = factorFrame(sys);
      res = solveFrame(sys, linearFac, loads);
      res.iterations = 1; res.inactiveMembers = [];
    } else {
      res = solveWithInactive(sys, loads, opts.maxIter || 50);
    }
    res.name = name;
    results[name] = res;
  }
  return { system: sys, order: names, results, warnings: sys.warnings };
}

function solveWithInactive(sys, loads, maxIter) {
  const { elems } = sys;
  const active = new Uint8Array(elems.length).fill(1);
  const de = new Float64Array(12); const dl = new Float64Array(12);
  for (let it = 1; it <= maxIter; it++) {
    const fac = factorFrame(sys, { active });
    const res = solveFrame(sys, fac, loads, { active });
    let changed = false;
    for (let e = 0; e < elems.length; e++) {
      const el = elems[e];
      if (!el.nonlinear) continue;
      const want = el.type === 'tension' ? 1 : -1;
      if (active[e]) {
        const mf = res.memberForces[e];
        const N = 0.5 * (mf.endA.N + mf.endB.N);
        if (N * want < -1e-9 * Math.max(1, Math.abs(N))) { active[e] = 0; changed = true; }
      } else {
        for (let i = 0; i < 12; i++) de[i] = res.displacements[el.dofs[i]];
        matVec12(el.T, de, dl);
        const elong = dl[6] - dl[0];
        if (elong * want > 1e-12 * el.L) { active[e] = 1; changed = true; }
      }
    }
    if (!changed) {
      res.iterations = it;
      res.inactiveMembers = [];
      for (let e = 0; e < elems.length; e++) if (!active[e]) res.inactiveMembers.push(e);
      return res;
    }
  }
  throw new Error('Tension-only / compression-only member iteration did not converge.');
}

/* ------------------------------------------------------------------ diagrams */

/* Section forces at station s from end-A forces and the loads on [0, s]. out = [N, Vy, Vz, T, My, Mz]. */
function sectionForces(el, f, its, s, out) {
  let Fx = f[0]; let Fy = f[1]; let Fz = f[2]; let Mx = f[3];
  let My = f[4] + s * f[2];
  let Mz = f[5] - s * f[1];
  const gq = gauss8();
  for (const it of its) {
    const fac = it.fac;
    if (it.t === 'd') {
      if (s <= it.s1) continue;
      const hi = Math.min(s, it.s2);
      const half = 0.5 * (hi - it.s1); const mid = 0.5 * (hi + it.s1);
      for (let g = 0; g < gq.x.length; g++) {
        const t = mid + half * gq.x[g];
        const w = gq.w[g] * half * fac;
        const q = it.f(t);
        Fx += w * q[0]; Fy += w * q[1]; Fz += w * q[2];
        My += w * (s - t) * q[2];
        Mz -= w * (s - t) * q[1];
      }
    } else if (it.t === 'p') {
      if (it.s > s) continue;
      Fx += fac * it.P[0]; Fy += fac * it.P[1]; Fz += fac * it.P[2];
      Mx += fac * it.M[0];
      My += fac * (it.M[1] + (s - it.s) * it.P[2]);
      Mz += fac * (it.M[2] - (s - it.s) * it.P[1]);
    }
  }
  out[0] = -Fx; out[1] = -Fy; out[2] = -Fz; out[3] = -Mx; out[4] = -My; out[5] = -Mz;
  return out;
}

/**
 * Internal force and deflection diagrams along a member for a solved load case. Section forces follow from
 * statics of the member (end forces + loads up to the station), deflections from integrating
 * θ' = M/EI (+ thermal curvature), v' = θ + V/(G As), u' = N/EA + αΔT between stations (Simpson's rule on a
 * quadratic through station, midpoint and next station) and matching the member end displacements.
 * @param {object} result one entry of analyzeFrame(...).results
 * @param {number} member member index
 * @param {number} [nStations=21] number of equally spaced stations (>= 2)
 * @returns {{s:Float64Array, N:Float64Array, Vy:Float64Array, Vz:Float64Array, T:Float64Array, My:Float64Array, Mz:Float64Array,
 *   ux:Float64Array, uy:Float64Array, uz:Float64Array, rx:Float64Array, maxAbs:object}} stations [m], forces [N], moments [N·m],
 *   local deflections [m], twist [rad]; maxAbs holds the largest absolute value of each quantity
 */
export function memberDiagram(result, member, nStations = 21) {
  const int = result && result.internal;
  if (!int) throw new RangeError('result must come from analyzeFrame / solveFrame');
  const el = int.sys.elems[member];
  if (!el) throw new RangeError('member index is out of range');
  if (!(nStations >= 2)) throw new RangeError('nStations must be at least 2');
  const n = Math.floor(nStations);
  const its = int.items[member];
  const f = result.memberForces[member].f;
  const L = el.L;
  const s = new Float64Array(n);
  const out = { s };
  const keys = ['N', 'Vy', 'Vz', 'T', 'My', 'Mz'];
  for (const k of keys) out[k] = new Float64Array(n);
  const tmp = new Float64Array(6);
  const mids = [];
  for (let i = 0; i < n; i++) {
    s[i] = L * i / (n - 1);
    sectionForces(el, f, its, s[i], tmp);
    for (let k = 0; k < 6; k++) out[keys[k]][i] = tmp[k];
    if (i < n - 1) mids.push(Float64Array.from(sectionForces(el, f, its, L * (i + 0.5) / (n - 1), tmp)));
  }
  // a point load exactly at the end station belongs to the member: enforce S(L) = end-B forces
  out.N[n - 1] = f[6]; out.Vy[n - 1] = f[7]; out.Vz[n - 1] = f[8]; out.T[n - 1] = f[9]; out.My[n - 1] = f[10]; out.Mz[n - 1] = f[11];

  const { E, G, alpha } = el.mat;
  const sec = el.sec;
  let dT = 0; let gy = 0; let gz = 0;
  for (const it of its) if (it.t === 'th') { dT += it.fac * it.dT; gy += it.fac * it.gy; gz += it.fac * it.gz; }
  const useShear = !el.truss && int.sys.shearDeformation && el.member.shear !== false;
  const kvy = useShear && sec.Asy > 0 && el.phiY > 0 ? 1 / (G * sec.Asy) : 0;
  const kvz = useShear && sec.Asz > 0 && el.phiZ > 0 ? 1 / (G * sec.Asz) : 0;
  const h = L / (n - 1);
  // particular integrals with zero initial values
  const pu = new Float64Array(n); const pr = new Float64Array(n);
  const pv = new Float64Array(n); const ptz = new Float64Array(n);  // xy plane: v, θz
  const pw = new Float64Array(n); const pty = new Float64Array(n);  // xz plane: w, θ* = −θy
  const EA = E * sec.A; const GJ = G * sec.J; const EIy = E * sec.Iy; const EIz = E * sec.Iz;
  for (let i = 0; i < n - 1; i++) {
    const m = mids[i];
    const simp = (a, mm, c) => h * (a + 4 * mm + c) / 6;
    const dbl = (a, mm) => h * h * (a / 6 + mm / 3);
    pu[i + 1] = pu[i] + simp(out.N[i], m[0], out.N[i + 1]) / EA + alpha * dT * h;
    pr[i + 1] = pr[i] + (GJ > 0 ? simp(out.T[i], m[3], out.T[i + 1]) / GJ : 0);
    // xy plane: θz' = Mz/EIz − α gy ; v' = θz + Vy/(G Asy)
    const cz0 = out.Mz[i] / EIz - alpha * gy; const czm = m[5] / EIz - alpha * gy; const cz1 = out.Mz[i + 1] / EIz - alpha * gy;
    ptz[i + 1] = ptz[i] + simp(cz0, czm, cz1);
    pv[i + 1] = pv[i] + ptz[i] * h + dbl(cz0, czm) + kvy * simp(out.Vy[i], m[1], out.Vy[i + 1]);
    // xz plane with θ* = −θy: θ*' = −My/EIy − α gz ; w' = θ* + Vz/(G Asz)
    const cy0 = -out.My[i] / EIy - alpha * gz; const cym = -m[4] / EIy - alpha * gz; const cy1 = -out.My[i + 1] / EIy - alpha * gz;
    pty[i + 1] = pty[i] + simp(cy0, cym, cy1);
    pw[i + 1] = pw[i] + pty[i] * h + dbl(cy0, cym) + kvz * simp(out.Vz[i], m[2], out.Vz[i + 1]);
  }
  const d = int.endDisp[member];
  const rel = new Uint8Array(12);
  for (const i of el.relIdx) rel[i] = 1;
  const affine = (pA, pB, relA, relB, p) => (!relA ? pA : (!relB ? pB - p[n - 1] : pA));
  const u0 = affine(d[0], d[6], rel[0], rel[6], pu);
  const r0 = affine(d[3], d[9], rel[3], rel[9], pr);
  // bending planes: unknown v0, θ0 from the first two available of [vA, vB, θA, θB]
  const plane = (vA, vB, tA, tB, rvA, rvB, rtA, rtB, p, pt) => {
    if (el.truss) { rtA = 1; rtB = 1; }
    let v0; let t0;
    if (!rvA && !rvB) { v0 = vA; t0 = (vB - vA - p[n - 1]) / L; }
    else if (!rvA) { v0 = vA; t0 = !rtA ? tA : tB - pt[n - 1]; }
    else if (!rvB) { t0 = !rtA ? tA : tB - pt[n - 1]; v0 = vB - t0 * L - p[n - 1]; }
    else { v0 = vA; t0 = tA; }
    return [v0, t0];
  };
  const [v0, tz0] = plane(d[1], d[7], d[5], d[11], rel[1], rel[7], rel[5], rel[11], pv, ptz);
  const [w0, ty0] = plane(d[2], d[8], -d[4], -d[10], rel[2], rel[8], rel[4], rel[10], pw, pty);
  out.ux = new Float64Array(n); out.uy = new Float64Array(n); out.uz = new Float64Array(n); out.rx = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    out.ux[i] = u0 + pu[i];
    out.rx[i] = r0 + pr[i];
    out.uy[i] = v0 + tz0 * s[i] + pv[i];
    out.uz[i] = w0 + ty0 * s[i] + pw[i];
  }
  out.maxAbs = {};
  for (const k of [...keys, 'ux', 'uy', 'uz']) {
    let mx = 0;
    for (let i = 0; i < n; i++) mx = Math.max(mx, Math.abs(out[k][i]));
    out.maxAbs[k] = mx;
  }
  return out;
}

/**
 * Elastic stresses in a tube or I/box-type section from section forces. Tubes: σ = N/A ± √(My² + Mz²)/Wel,
 * τ = T D/(2J) + 2√(Vy² + Vz²)/A (thin-wall peak shear). Other sections: σ = N/A ± My/Wely ± Mz/Welz,
 * τ = Vz/Asz + Vy/Asy + T·t/J when shear areas are available. von Mises equivalent σeq = √(σ² + 3τ²) is
 * evaluated with the peak normal and peak shear stress taken together (conservative).
 * @param {object} o
 * @param {object} o.section section record (tubeSection, iSection, …)
 * @param {number} [o.N=0] axial force [N] (tension positive)
 * @param {number} [o.Vy=0] shear [N]
 * @param {number} [o.Vz=0] shear [N]
 * @param {number} [o.T=0] torque [N·m]
 * @param {number} [o.My=0] bending moment [N·m]
 * @param {number} [o.Mz=0] bending moment [N·m]
 * @returns {{sigmaAxial:number, sigmaBending:number, sigmaMax:number, sigmaMin:number, tau:number, vonMises:number}} Pa
 */
export function sectionStress({ section, N = 0, Vy = 0, Vz = 0, T = 0, My = 0, Mz = 0 }) {
  const sec = section;
  if (!sec || !(sec.A > 0)) throw new RangeError('section must be a section property record');
  const sa = N / sec.A;
  let sb; let tau;
  if (sec.type === 'tube' || sec.type === 'round') {
    sb = Math.hypot(My, Mz) / sec.Wely;
    tau = Math.abs(T) * sec.D / (2 * sec.J) + (sec.type === 'tube' ? 2 : 4 / 3) * Math.hypot(Vy, Vz) / sec.A;
  } else {
    const wy = sec.Wely || (sec.Iy / (sec.cz || NaN));
    const wz = sec.Welz || (sec.Iz / (sec.cy || NaN));
    sb = Math.abs(My) / wy + Math.abs(Mz) / wz;
    const tmax = Math.max(sec.tf || 0, sec.tw || 0, sec.t || 0);
    tau = (sec.Asz > 0 ? Math.abs(Vz) / sec.Asz : 0) + (sec.Asy > 0 ? Math.abs(Vy) / sec.Asy : 0)
      + (sec.J > 0 && tmax > 0 ? Math.abs(T) * tmax / sec.J : 0);
  }
  const smax = sa + sb; const smin = sa - sb;
  const sAbs = Math.max(Math.abs(smax), Math.abs(smin));
  return { sigmaAxial: sa, sigmaBending: sb, sigmaMax: smax, sigmaMin: smin, tau, vonMises: Math.sqrt(sAbs * sAbs + 3 * tau * tau) };
}

/**
 * Mesh refinement: subdivide every member (or the listed members) into n equal elements. Releases and rigid
 * offsets stay at the original member ends, the local axes are preserved, and member loads of every load
 * case are mapped onto the sub-members.
 * @param {object} model frame model
 * @param {number} n number of elements per member (>= 1)
 * @param {object} [o]
 * @param {number[]} [o.members] indices of members to subdivide (default all)
 * @returns {{model:object, memberMap:number[][], nodeMap:number[][]}} refined model; memberMap[i] lists the new member indices of
 *   original member i (A → B); nodeMap[i] the node indices along it (n + 1 entries)
 */
export function refineModel(model, n, { members } = {}) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError('n must be a positive integer');
  const sys = buildFrameSystem(model, { renumber: false });
  const nodes = [];
  for (let i = 0; i < sys.nNode; i++) nodes.push([sys.xyz[3 * i], sys.xyz[3 * i + 1], sys.xyz[3 * i + 2]]);
  const pick = members ? new Set(members) : null;
  const out = [];
  const memberMap = []; const nodeMap = [];
  for (let m = 0; m < model.members.length; m++) {
    const mem = model.members[m];
    const el = sys.elems[m];
    const k = !pick || pick.has(m) ? n : 1;
    const ids = [mem.a];
    for (let j = 1; j < k; j++) {
      const w = j / k;
      nodes.push([el.pa[0] + w * (el.pb[0] - el.pa[0]), el.pa[1] + w * (el.pb[1] - el.pa[1]), el.pa[2] + w * (el.pb[2] - el.pa[2])]);
      ids.push(nodes.length - 1);
    }
    ids.push(mem.b);
    const list = [];
    for (let j = 0; j < k; j++) {
      const sub = { ...mem, a: ids[j], b: ids[j + 1], localZ: el.R[2].slice(), betaDeg: 0 };
      delete sub.releases; delete sub.offsetA; delete sub.offsetB;
      const rel = {};
      if (j === 0) { if (mem.offsetA) sub.offsetA = mem.offsetA; if (mem.releases && mem.releases.a) rel.a = mem.releases.a; }
      if (j === k - 1) { if (mem.offsetB) sub.offsetB = mem.offsetB; if (mem.releases && mem.releases.b) rel.b = mem.releases.b; }
      if (rel.a || rel.b) sub.releases = rel;
      sub.parent = m;
      list.push(out.length);
      out.push(sub);
    }
    memberMap.push(list); nodeMap.push(ids);
  }
  const mapLoads = (lc) => {
    const res = { ...lc };
    if (!lc.member) return res;
    res.member = [];
    for (const ml of lc.member) {
      const el = sys.elems[ml.member];
      const list = memberMap[ml.member];
      const k = list.length;
      const h = el.L / k;
      if (ml.type === 'thermal') { for (const id of list) res.member.push({ ...ml, member: id }); continue; }
      if (ml.type === 'point') {
        const s = ml.s !== undefined ? ml.s : ml.xi * el.L;
        const j = Math.min(k - 1, Math.floor(s / h));
        const cp = { ...ml, member: list[j], s: s - j * h };
        delete cp.xi;
        res.member.push(cp);
        continue;
      }
      const [s1, s2] = spanOf(el, ml);
      for (let j = 0; j < k; j++) {
        const a = Math.max(s1, j * h); const b = Math.min(s2, (j + 1) * h);
        if (!(b > a + 1e-12 * el.L)) continue;
        const cp = { ...ml, member: list[j], s1: a - j * h, s2: b - j * h };
        delete cp.xi1; delete cp.xi2;
        if (ml.type === 'linear') {
          const at = (s) => { const w = (s - s1) / (s2 - s1); return [0, 1, 2].map((c) => ml.q1[c] + w * (ml.q2[c] - ml.q1[c])); };
          cp.q1 = at(a); cp.q2 = at(b);
        } else if (ml.type === 'function') {
          const off = j * h; const fn = ml.fn;
          cp.fn = (s, p) => fn(s + off, p);
        }
        res.member.push(cp);
      }
    }
    return res;
  };
  const refined = { ...model, nodes, members: out };
  if (model.loadCases) refined.loadCases = model.loadCases.map(mapLoads);
  return { model: refined, memberMap, nodeMap };
}
