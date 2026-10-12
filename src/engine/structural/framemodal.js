/**
 * Mass matrices and free-vibration analysis of 3-D frames.
 *
 *  - Consistent element mass (cubic Hermite translational mass, Przemieniecki 1968 §11.3; rotary inertia of the
 *    cross-section neglected except about the member axis) or diagonal lumped mass by the HRZ scheme
 *    (Hinton, Rock & Zienkiewicz 1976: ½ m L per node in translation, m L³/78 in bending rotation).
 *  - Contained fluid in flooded tubes, marine growth and hydrodynamic added mass Ca ρw π D²/4 of submerged tubes
 *    (transverse directions only; DNV-RP-C205 §6.9 with Ca = 1.0 for a circular cylinder away from boundaries).
 *  - Natural frequencies and mode shapes by subspace iteration on the skyline-factorised stiffness.
 *  - Modal participation factors Γ = φᵀ M r and effective modal masses Γ² per global direction (modes are
 *    mass-normalised), Rayleigh damping coefficients.
 * SI units. Axis conventions as frame3d.js.
 */

import { subspaceIteration, eigGenSym } from '../math/eigen.js';
import {
  buildFrameSystem, factorFrame, elementToGlobal, elementMatVec, memberMassPerLength, assembleDense,
} from './frame3d.js';

const XY = [1, 5, 7, 11];
const XZ = [2, 4, 8, 10];
const XZS = [1, -1, 1, -1];

/**
 * Local 12×12 mass matrix of a prismatic member.
 * Consistent: axial m L/6 [2 1; 1 2], torsion (ρ Ip) L/6 [2 1; 1 2], bending m L/420 [156 22L 54 −13L; … 4L² 13L −3L²; … 156 −22L; … 4L²].
 * Lumped (HRZ): m L/2 per translation, m L³/78 per bending rotation, ρ Ip L/2 in torsion.
 * @param {object} o
 * @param {number} o.L length [m]
 * @param {number} o.mAxial mass per unit length effective along the member axis [kg/m]
 * @param {number} [o.mTrans=o.mAxial] mass per unit length effective in the transverse directions [kg/m]
 * @param {number} [o.rhoIp=0] polar mass moment of inertia per unit length ρ·Ip [kg·m]
 * @param {boolean} [o.lumped=false] diagonal HRZ lumping
 * @param {boolean} [o.truss=false] bar element (linear interpolation in all three translations, no rotational mass)
 * @returns {Float64Array} row-major 12×12 matrix [kg, kg·m, kg·m²]
 */
export function beamLocalMass({ L, mAxial, mTrans, rhoIp = 0, lumped = false, truss = false }) {
  if (!(L > 0)) throw new RangeError('L must be positive');
  if (!(mAxial >= 0)) throw new RangeError('mAxial must be non-negative');
  const mt = mTrans === undefined ? mAxial : mTrans;
  if (!(mt >= 0)) throw new RangeError('mTrans must be non-negative');
  const M = new Float64Array(144);
  const set = (i, j, v) => { M[i * 12 + j] = v; M[j * 12 + i] = v; };
  if (lumped) {
    set(0, 0, mAxial * L / 2); set(6, 6, mAxial * L / 2);
    for (const i of [1, 2, 7, 8]) set(i, i, mt * L / 2);
    if (!truss) {
      set(3, 3, rhoIp * L / 2); set(9, 9, rhoIp * L / 2);
      for (const i of [4, 5, 10, 11]) set(i, i, mt * L ** 3 / 78);
    }
    return M;
  }
  set(0, 0, mAxial * L / 3); set(6, 6, mAxial * L / 3); set(0, 6, mAxial * L / 6);
  if (truss) {
    for (const [a, b] of [[1, 7], [2, 8]]) { set(a, a, mt * L / 3); set(b, b, mt * L / 3); set(a, b, mt * L / 6); }
    return M;
  }
  set(3, 3, rhoIp * L / 3); set(9, 9, rhoIp * L / 3); set(3, 9, rhoIp * L / 6);
  const c = mt * L / 420;
  const b4 = [
    [156, 22 * L, 54, -13 * L],
    [22 * L, 4 * L * L, 13 * L, -3 * L * L],
    [54, 13 * L, 156, -22 * L],
    [-13 * L, -3 * L * L, -22 * L, 4 * L * L],
  ];
  for (let i = 0; i < 4; i++) {
    for (let j = 0; j < 4; j++) {
      M[XY[i] * 12 + XY[j]] = c * b4[i][j];
      M[XZ[i] * 12 + XZ[j]] = c * XZS[i] * XZS[j] * b4[i][j];
    }
  }
  return M;
}

/**
 * Assemble the mass description of a frame system: one global 12×12 matrix per element plus the diagonal of
 * lumped node masses.
 * @param {object} sys system from buildFrameSystem
 * @param {object} [o]
 * @param {boolean} [o.lumped=false] HRZ lumped element masses instead of consistent
 * @param {number} [o.waterline] still-water elevation z [m]; members below it receive added mass and (if flooded) contained water
 * @param {number} [o.rhoW=1025] sea-water density [kg/m³]
 * @param {number} [o.Ca=1.0] added-mass coefficient of submerged tubes [-] (member.Ca overrides)
 * @param {{thickness:number, rho?:number, zTop?:number, zBottom?:number}} [o.marineGrowth] marine growth zone [m, kg/m³]
 * @returns {{mats:Float64Array[], nodeDiag:Float64Array, mulM:(x:Float64Array)=>Float64Array, structuralMass:number, addedMass:number,
 *   memberMass:Array<{axial:number, transverse:number}>}} element matrices in global dofs, node-mass diagonal over all dofs, M·x operator on
 *   free equations, total structural + contained + node mass [kg], total hydrodynamic added mass [kg], per-member kg/m
 */
export function frameMass(sys, { lumped = false, waterline, rhoW = 1025, Ca = 1.0, marineGrowth } = {}) {
  const mats = [];
  const memberMass = [];
  let structuralMass = 0; let addedMass = 0;
  for (const el of sys.elems) {
    const za = el.pa[2]; const zb = el.pb[2];
    const fracBelow = (zc) => {
      if (zc === undefined) return 0;
      const lo = Math.min(za, zb); const hi = Math.max(za, zb);
      if (hi <= zc) return 1;
      if (lo >= zc) return 0;
      return (zc - lo) / (hi - lo);
    };
    const wet = fracBelow(waterline);
    let tg = 0; let growthFrac = 0;
    if (marineGrowth && marineGrowth.thickness > 0) {
      const top = marineGrowth.zTop === undefined ? (waterline === undefined ? Infinity : waterline) : marineGrowth.zTop;
      const bot = marineGrowth.zBottom === undefined ? -Infinity : marineGrowth.zBottom;
      const fTop = top === Infinity ? 1 : fracBelow(top);
      const fBot = bot === -Infinity ? 0 : fracBelow(bot);
      growthFrac = Math.max(0, fTop - fBot);
      tg = marineGrowth.thickness;
    }
    const base = memberMassPerLength(el, { rhoW, mgThickness: tg, mgRho: marineGrowth && marineGrowth.rho !== undefined ? marineGrowth.rho : 1325 });
    // contained water of a flooded member only up to the waterline (when one is given)
    const fluid = el.member.flooded && el.member.fluidRho === undefined && waterline !== undefined ? base.fluid * wet : base.fluid;
    const mAxial = base.steel + base.extra + fluid + base.growth * growthFrac;
    let ma = 0;
    if (el.member.addedMass !== undefined) ma = el.member.addedMass * (waterline === undefined ? 1 : wet);
    else if (wet > 0 && el.sec.D) {
      const De = el.sec.D + 2 * tg * (growthFrac > 0 ? 1 : 0);
      const ca = el.member.Ca !== undefined ? el.member.Ca : Ca;
      ma = ca * rhoW * Math.PI * De * De / 4 * wet;
    }
    const mTrans = mAxial + ma;
    structuralMass += mAxial * el.L;
    addedMass += ma * el.L;
    memberMass.push({ axial: mAxial, transverse: mTrans });
    const ml = beamLocalMass({ L: el.L, mAxial, mTrans, rhoIp: el.mat.rho * (el.sec.Ip || el.sec.Iy + el.sec.Iz), lumped, truss: el.truss });
    mats.push(elementToGlobal(el, ml));
  }
  const nodeDiag = new Float64Array(sys.nDof);
  for (const [i, nm] of (sys.model.nodeMasses || []).entries()) {
    if (!Number.isInteger(nm.node) || nm.node < 0 || nm.node >= sys.nNode) throw new RangeError(`model.nodeMasses[${i}].node is invalid`);
    if (!(nm.m >= 0)) throw new RangeError(`model.nodeMasses[${i}].m must be non-negative`);
    for (let k = 0; k < 3; k++) nodeDiag[6 * nm.node + k] += nm.m;
    if (nm.inertia) for (let k = 0; k < 3; k++) nodeDiag[6 * nm.node + 3 + k] += nm.inertia[k] || 0;
    structuralMass += nm.m;
  }
  const mulM = (x) => {
    const y = elementMatVec(sys, mats, x);
    for (let q = 0; q < sys.nFree; q++) y[q] += nodeDiag[sys.eqDof[q]] * x[q];
    return y;
  };
  return { mats, nodeDiag, mulM, structuralMass, addedMass, memberMass };
}

/**
 * Mass influence vectors M·r_d for the three global translations d = x, y, z, with r_d the rigid-body unit
 * translation over ALL dofs: supports move with the ground, so the mass coupling between free and restrained
 * dofs contributes to the excitation vector and the restrained mass counts in the total (Chopra §9.4–9.5).
 * @param {object} sys system from buildFrameSystem
 * @param {object} mass result of frameMass
 * @returns {{Mr:Float64Array[], total:number[]}} Mr[d] over the free equations [kg]; total rigid-body mass per direction [kg]
 */
export function massInfluence(sys, mass) {
  const n = sys.nFree;
  const Mr = []; const total = [];
  for (let d = 0; d < 3; d++) {
    const y = new Float64Array(sys.nDof);
    for (let e = 0; e < sys.elems.length; e++) {
      const A = mass.mats[e];
      const dofs = sys.elems[e].dofs;
      for (let i = 0; i < 12; i++) {
        let s = 0;
        for (let k = 0; k < 12; k++) if (dofs[k] % 6 === d) s += A[i * 12 + k];
        y[dofs[i]] += s;
      }
    }
    for (let g = d; g < sys.nDof; g += 6) y[g] += mass.nodeDiag[g];
    const mr = new Float64Array(n);
    for (let q = 0; q < n; q++) mr[q] = y[sys.eqDof[q]];
    Mr.push(mr);
    let s = 0;
    for (let g = d; g < sys.nDof; g += 6) s += y[g];
    total.push(s);
  }
  return { Mr, total };
}

/**
 * Natural frequencies, mode shapes, participation factors and effective modal masses of a frame.
 * K φ = ω² M φ is solved for the lowest modes by subspace iteration (Bathe 1996 §11.6), or by the dense
 * generalised symmetric solver for small systems; the modes are M-orthonormal. Groups of nearly equal
 * frequencies (repeated or symmetric members) are handled by enlarging the iteration subspace until it spans the
 * group (math/eigen.js subspaceIteration); `converged` is false only if a group wider than eight times nModes
 * remains unresolved on a model too large for the dense solver. Participation factor for direction d: Γ = φᵀ M r_d with r_d the rigid-body unit translation;
 * effective modal mass Γ² (Chopra, Dynamics of Structures, §13.2).
 * @param {object} model frame model (see frame3d.buildFrameSystem)
 * @param {object} [o]
 * @param {number} [o.nModes=6] number of modes (capped at the number of free equations)
 * @param {boolean} [o.lumped=false] lumped (HRZ) instead of consistent mass
 * @param {number} [o.waterline] still-water level [m] for added mass / contained water
 * @param {number} [o.rhoW=1025] [kg/m³]
 * @param {number} [o.Ca=1.0] added mass coefficient [-]
 * @param {object} [o.marineGrowth] { thickness, rho, zTop, zBottom }
 * @param {boolean} [o.shearDeformation=true]
 * @param {number} [o.tol=1e-9] relative eigenvalue tolerance
 * @param {object} [o.system] reuse a system from buildFrameSystem
 * @returns {{frequencies:number[], omegas:number[], periods:number[], modes:Float64Array[], participation:Array<{x:number,y:number,z:number}>,
 *   effectiveMass:Array<{x:number,y:number,z:number}>, massFraction:Array<{x:number,y:number,z:number}>, totalMass:{x:number,y:number,z:number},
 *   structuralMass:number, addedMass:number, converged:boolean, iterations:number, system:object, mass:object}}
 *   frequencies [Hz], omegas [rad/s], periods [s], modes over all 6·nNode dofs (mass-normalised, [1/√kg]), participation [√kg], effective
 *   mass [kg], totalMass = r_dᵀ M r_d over all dofs (rigid-body mass per direction, incl. added mass) [kg]
 */
export function modalAnalysis(model, {
  nModes = 6, lumped = false, waterline, rhoW = 1025, Ca = 1.0, marineGrowth, shearDeformation = true, tol = 1e-9, system,
} = {}) {
  const sys = system || buildFrameSystem(model, { shearDeformation });
  const n = sys.nFree;
  if (n === 0) throw new RangeError('model has no free degrees of freedom');
  if (!Number.isInteger(nModes) || nModes < 1) throw new RangeError('nModes must be a positive integer');
  const mass = frameMass(sys, { lumped, waterline, rhoW, Ca, marineGrowth });
  const fac = factorFrame(sys);
  const p = Math.min(nModes, n);
  // Small systems, or requests for a large share of the spectrum, go to the dense solver: subspace iteration
  // cannot resolve the stiffest modes of a small model within its working precision.
  let eig = null;
  if (n <= 1200 && (n <= 60 || 3 * p > n)) {
    try {
      const K = assembleDense(sys, sys.elems.map((e) => e.ke), true);
      const M = assembleDense(sys, mass.mats, false);
      for (let q = 0; q < n; q++) M[q][q] += mass.nodeDiag[sys.eqDof[q]];
      const d = eigGenSym(K, M, { nModes: p });
      eig = { values: d.values, vectors: d.vectors, converged: true, iterations: 0 };
    } catch (e) {
      if (!(e instanceof RangeError)) throw e;
      eig = null; // singular mass matrix (massless dofs): fall through to subspace iteration
    }
  }
  if (!eig) eig = subspaceIteration({ n, nModes: p, solveK: fac.solve, mulM: mass.mulM, tol, maxIter: 200 });
  if (!eig.converged && n <= 1200) {
    // a cluster of near-equal modes wider than the iteration subspace: the dense solver resolves it directly
    try {
      const K = assembleDense(sys, sys.elems.map((e) => e.ke), true);
      const M = assembleDense(sys, mass.mats, false);
      for (let q = 0; q < n; q++) M[q][q] += mass.nodeDiag[sys.eqDof[q]];
      const d = eigGenSym(K, M, { nModes: p });
      eig = { values: d.values, vectors: d.vectors, converged: true, iterations: eig.iterations };
    } catch (e) {
      if (!(e instanceof RangeError)) throw e;
    }
  }
  const dirs = ['x', 'y', 'z'];
  const Mr = [];
  const totalMass = {};
  const infl = massInfluence(sys, mass);
  for (let d = 0; d < 3; d++) { Mr.push(infl.Mr[d]); totalMass[dirs[d]] = infl.total[d]; }
  const omegas = []; const frequencies = []; const periods = []; const modes = [];
  const participation = []; const effectiveMass = []; const massFraction = [];
  for (let k = 0; k < eig.values.length; k++) {
    const lam = eig.values[k];
    if (!(lam > 0) || !Number.isFinite(lam)) continue;
    const w = Math.sqrt(lam);
    const phi = eig.vectors[k];
    // sign convention: largest component positive
    let big = 0;
    for (let q = 0; q < n; q++) if (Math.abs(phi[q]) > Math.abs(big)) big = phi[q];
    const sg = big < 0 ? -1 : 1;
    const full = new Float64Array(sys.nDof);
    for (let q = 0; q < n; q++) full[sys.eqDof[q]] = sg * phi[q];
    const g = {}; const me = {}; const fr = {};
    for (let d = 0; d < 3; d++) {
      let s = 0;
      for (let q = 0; q < n; q++) s += sg * phi[q] * Mr[d][q];
      g[dirs[d]] = s; me[dirs[d]] = s * s; fr[dirs[d]] = totalMass[dirs[d]] > 0 ? s * s / totalMass[dirs[d]] : 0;
    }
    omegas.push(w); frequencies.push(w / (2 * Math.PI)); periods.push(2 * Math.PI / w); modes.push(full);
    participation.push(g); effectiveMass.push(me); massFraction.push(fr);
  }
  return {
    frequencies, omegas, periods, modes, participation, effectiveMass, massFraction, totalMass,
    structuralMass: mass.structuralMass, addedMass: mass.addedMass, converged: eig.converged, iterations: eig.iterations, system: sys, mass,
  };
}

/**
 * Rayleigh damping coefficients C = α M + β K matching damping ratios ζ1, ζ2 at two frequencies:
 * ζ(ω) = α/(2ω) + βω/2 (Chopra §11.4.1).
 * @param {object} o
 * @param {number} o.f1 first target frequency [Hz]
 * @param {number} o.f2 second target frequency [Hz], ≠ f1
 * @param {number} o.zeta1 damping ratio at f1 [-]
 * @param {number} [o.zeta2=o.zeta1] damping ratio at f2 [-]
 * @returns {{alpha:number, beta:number, zetaAt:(f:number)=>number}} α [1/s], β [s], and ζ(f) for f in Hz
 */
export function rayleighDamping({ f1, f2, zeta1, zeta2 }) {
  if (!(f1 > 0) || !(f2 > 0) || f1 === f2) throw new RangeError('f1 and f2 must be positive and different');
  if (!(zeta1 >= 0)) throw new RangeError('zeta1 must be non-negative');
  const z2 = zeta2 === undefined ? zeta1 : zeta2;
  if (!(z2 >= 0)) throw new RangeError('zeta2 must be non-negative');
  const w1 = 2 * Math.PI * f1; const w2 = 2 * Math.PI * f2;
  const alpha = 2 * w1 * w2 * (zeta1 * w2 - z2 * w1) / (w2 * w2 - w1 * w1);
  const beta = 2 * (z2 * w2 - zeta1 * w1) / (w2 * w2 - w1 * w1);
  return { alpha, beta, zetaAt: (f) => { const w = 2 * Math.PI * f; return alpha / (2 * w) + beta * w / 2; } };
}
