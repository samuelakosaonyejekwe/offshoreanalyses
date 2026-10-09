/**
 * Unbonded flexible pipe cross-section model (API Spec 17J / API RP 17B layer terminology).
 *
 * Axisymmetric model (derived here from first principles; same family as Feret & Bournazel 1987,
 * Witz & Tan 1992, McNamara & Harte 1992):
 *
 *  Global unknowns: axial strain eps, twist per unit length tau. Per layer i: mean radial displacement u_i and
 *  thickness change d_i. Per interface j (between layers j and j+1): contact pressure pc_j >= 0.
 *
 *  Helical metallic layer (carcass, pressure armour, tensile armour): n wires of area A, lay angle alpha,
 *      eps_w = cos^2(a) eps + lay R sin(a) cos(a) tau + sin^2(a) u/R - alphaT dT
 *      Fw = E A eps_w;  N = n Fw cos(a);  Mt = lay n Fw R sin(a);  N_theta = n Fw sin^2(a)/(2 pi R cos(a))
 *    Wires are taken as radially rigid: d = t alphaT dT.
 *  Polymer / tape layer (isotropic, thin shell with through-thickness strain):
 *      strains (eps, u/R, d/t), 3-D Hooke's law with thermal strain; N = sigma_z 2 pi R t;
 *      Mt = G J tau; N_theta = sigma_theta t; sigma_r R = -(P_in R_in + P_out R_out)/2  (work-conjugate to d).
 *  Radial equilibrium of every layer (thin-ring statics): N_theta = P_in R_in - P_out R_out, where the face
 *  pressures are contact pressures plus, on sealing layers only, the fluid pressure of the adjoining zone. Fluid
 *  surrounds the individual wires of non-sealing layers and therefore produces no net layer load.
 *  Compatibility at a closed interface: (u_{j+1} - d_{j+1}/2) - (u_j + d_j/2) = 0; an interface whose contact
 *  pressure would be tensile is released (pc_j = 0, gap >= 0) by an active-set iteration.
 *  Axial equilibrium: sum N_i = T_eff + sum over sealing layers (p_in pi r_in^2 - p_out pi r_out^2), i.e.
 *  effective tension plus end-cap forces; torque equilibrium: sum Mt_i = torque.
 *
 * Fluid zones: bore pressure acts on the inner face of the innermost sealing layer (the carcass is not
 * leak-tight); external pressure on the outer face of the outermost sealing layer; the annulus between them is
 * at pAnnulus (dry) or at external pressure (flooded).
 *
 * Limitations: no wire bending/torsion stiffness in the axisymmetric response, no interlocked-profile contact
 * mechanics (carcass and pressure armour are helical strips of equivalent area), no plasticity, small strains.
 */

import { solveDense } from '../math/linalg.js';
import { helixBendingStiffness, criticalCurvature, rectWireSection } from './helix.js';

const TWO_PI = 2 * Math.PI;
const DEG = Math.PI / 180;

function pos(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function fin(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}
function nonneg(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/**
 * Permissible utilisation factors of API Spec 17J (4th ed.) Table 6, as recalled: stress utilisation against the
 * structural capacity (yield or 0.9 x UTS as specified by the manufacturer) of the armour layers for normal
 * operation (recurrent), abnormal/extreme operation and hydrostatic test (FAT). Carcass buckling factor depends
 * on water depth, see {@link carcassBucklingFactor}. These are overridable defaults: confirm against the edition
 * of the specification governing the project before use.
 */
export const API17J_FACTORS = Object.freeze({
  tensileArmour: Object.freeze({ normal: 0.67, abnormal: 0.85, hydrotest: 0.91 }),
  pressureArmour: Object.freeze({ normal: 0.55, abnormal: 0.85, hydrotest: 0.91 }),
});

/**
 * Allowable polymer sheath bending strains of API Spec 17J §5.3.1 (as recalled; overridable defaults): PE and PA
 * 7.7 %; PVDF 7.0 % for static service and storage of dynamic pipes, 3.5 % for dynamic operation.
 */
export const SHEATH_STRAIN_LIMITS = Object.freeze({
  PE: Object.freeze({ static: 0.077, dynamic: 0.077 }),
  PA: Object.freeze({ static: 0.077, dynamic: 0.077 }),
  PVDF: Object.freeze({ static: 0.07, dynamic: 0.035 }),
});

/**
 * Carcass / collapse permissible utilisation as a function of maximum water depth, API Spec 17J Table 6 (as
 * recalled): 0.67 for Dmax <= 300 m, 0.85 for Dmax >= 900 m, linear in between
 * ((Dmax - 300)/600 x 0.18 + 0.67).
 * @param {number} depth maximum water depth, m (>= 0)
 * @returns {{factor:number}} permissible utilisation, -
 */
export function carcassBucklingFactor(depth) {
  nonneg(depth, 'depth');
  if (depth <= 300) return { factor: 0.67 };
  if (depth >= 900) return { factor: 0.85 };
  return { factor: (depth - 300) / 600 * 0.18 + 0.67 };
}

/**
 * Build the layer stack of an unbonded flexible pipe from the bore outwards.
 *
 * Each layer: { name, type, thickness [m], E [Pa], nu [-], rho [kg/m^3], k [W/m/K], alphaT [1/K],
 *   yield [Pa], uts [Pa], seal [bool], strainLimit / strainLimitDynamic [-], lockStrain [-],
 *   and for helical layers: alphaDeg [deg], lay (+1/-1), fill [-] (metal fraction of the layer annulus) or
 *   nWires + wireArea [m^2], wireWidth [m], wireThickness [m] }.
 * type is one of 'carcass', 'pressureSheath', 'pressureArmour', 'backupArmour', 'antiWear', 'tensileArmour',
 * 'insulation', 'outerSheath', 'tape'. A layer with alphaDeg is helical; the others are isotropic shells.
 * pressureSheath and outerSheath are sealing layers unless seal is set explicitly. The metal area of a helical
 * layer follows from n A = fill 2 pi R t cos(alpha) (area balance on the pipe cross-section).
 * Reference: API RP 17B §4 (pipe structure), API Spec 17J §3 definitions.
 * Validity: at least one sealing layer; thickness << radius for the thin-layer statics to be accurate.
 * @param {object} o
 * @param {number} o.innerDiameter bore diameter (inside of the innermost layer), m
 * @param {object[]} o.layers layer definitions, innermost first
 * @returns {{innerRadius:number, outerRadius:number, layers:object[], sealInner:number, sealOuter:number}}
 *   radii m; layers enriched with rIn, rOut, R (mean radius) [m], kind ('helical'|'shell'), nA (total wire area,
 *   m^2), n (wire count or null), alpha [rad]; indices of the innermost and outermost sealing layers
 */
export function buildFlexibleStack({ innerDiameter, layers }) {
  pos(innerDiameter, 'innerDiameter');
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  let r = innerDiameter / 2;
  const out = [];
  layers.forEach((src, i) => {
    const tag = `layers[${i}]`;
    const t = pos(src.thickness, `${tag}.thickness`);
    const L = { ...src, index: i, name: src.name ?? src.type ?? `layer ${i + 1}`, t, rIn: r, rOut: r + t, R: r + t / 2 };
    pos(L.E, `${tag}.E`);
    L.nu = src.nu ?? 0.3;
    if (!(L.nu > -1 && L.nu < 0.5)) throw new RangeError(`${tag}.nu must be in (-1, 0.5)`);
    L.rho = nonneg(src.rho ?? 0, `${tag}.rho`);
    L.alphaT = fin(src.alphaT ?? 0, `${tag}.alphaT`);
    L.seal = src.seal ?? (src.type === 'pressureSheath' || src.type === 'outerSheath');
    if (src.alphaDeg !== undefined) {
      if (!(src.alphaDeg > 0 && src.alphaDeg < 90)) throw new RangeError(`${tag}.alphaDeg must be in (0, 90)`);
      L.kind = 'helical';
      L.alpha = src.alphaDeg * DEG;
      L.lay = src.lay ?? 1;
      if (L.lay !== 1 && L.lay !== -1) throw new RangeError(`${tag}.lay must be +1 or -1`);
      const c = Math.cos(L.alpha);
      const ring = TWO_PI * L.R * t * c;
      if (src.nWires !== undefined && (src.wireArea !== undefined || (src.wireWidth && src.wireThickness))) {
        L.n = pos(src.nWires, `${tag}.nWires`);
        L.wireArea = src.wireArea ?? src.wireWidth * src.wireThickness;
        L.nA = L.n * pos(L.wireArea, `${tag}.wireArea`);
        L.fill = L.nA / ring;
      } else {
        if (!(src.fill > 0 && src.fill <= 1)) throw new RangeError(`${tag}.fill must be in (0, 1] (or give nWires and wireArea)`);
        L.fill = src.fill;
        L.nA = src.fill * ring;
        L.n = src.nWires ?? (src.wireWidth ? Math.max(1, Math.round(src.fill * TWO_PI * L.R * c / src.wireWidth)) : null);
        L.wireArea = L.n ? L.nA / L.n : null;
      }
      if (L.seal) throw new RangeError(`${tag}: a helical layer cannot be a sealing layer`);
    } else {
      L.kind = 'shell';
    }
    out.push(L);
    r += t;
  });
  const seals = out.filter((L) => L.seal).map((L) => L.index);
  if (seals.length === 0) throw new RangeError('layers must contain at least one sealing layer (pressureSheath/outerSheath)');
  return { innerRadius: innerDiameter / 2, outerRadius: r, layers: out, sealInner: seals[0], sealOuter: seals[seals.length - 1] };
}

/** Linear layer response to (eps, tau, u, d) with temperature rise dT: axial force, torque, hoop force per unit
 *  length and the thickness-equation term. */
function layerForces(L, eps, tau, u, d, dT) {
  if (L.kind === 'helical') {
    const c = Math.cos(L.alpha), s = Math.sin(L.alpha);
    const ew = c * c * eps + L.lay * L.R * s * c * tau + s * s * u / L.R - L.alphaT * dT;
    const F = L.E * L.nA * ew; // total force of all wires along their own axes
    return { N: F * c, M: L.lay * F * L.R * s, Nth: F * s * s / (TWO_PI * L.R * c), thick: L.E * (d - L.t * L.alphaT * dT), ew };
  }
  const G = L.E / (2 * (1 + L.nu));
  const lam = L.E * L.nu / ((1 + L.nu) * (1 - 2 * L.nu));
  const eth = u / L.R, er = d / L.t;
  const tr = eps + eth + er;
  const beta = (3 * lam + 2 * G) * L.alphaT * dT;
  const sz = lam * tr + 2 * G * eps - beta;
  const st = lam * tr + 2 * G * eth - beta;
  const sr = lam * tr + 2 * G * er - beta;
  const J = Math.PI * (L.rOut ** 4 - L.rIn ** 4) / 2;
  return { N: sz * TWO_PI * L.R * L.t, M: G * J * tau, Nth: st * L.t, thick: sr * L.R, sz, st, sr, shear: G * L.R * tau };
}

function fluidOnFaces(stack, pInt, pExt, pAnn) {
  const N = stack.layers.length;
  const fIn = new Array(N).fill(0), fOut = new Array(N).fill(0);
  for (const L of stack.layers) {
    if (!L.seal) continue;
    const i = L.index;
    fIn[i] = i === stack.sealInner ? pInt : pAnn;
    fOut[i] = i === stack.sealOuter ? pExt : pAnn;
  }
  return { fIn, fOut };
}

/** Solve the linear system for a given active set. ctl = { axial: {strain}|{force}, torsion: {twist}|{torque} }. */
function solveLinear(stack, ctl, fIn, fOut, dT, active, gap0) {
  const Ls = stack.layers, N = Ls.length, m = 3 * N + 1;
  const iu = (i) => 2 + i, id = (i) => 2 + N + i, ip = (j) => 2 + 2 * N + j;
  const A = [];
  for (let r = 0; r < m; r++) A.push(new Array(m).fill(0));
  const b = new Array(m).fill(0);
  const unit = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
  for (let i = 0; i < N; i++) {
    const L = Ls[i];
    const f0 = layerForces(L, 0, 0, 0, 0, dT[i]);
    const cols = [0, 1, iu(i), id(i)];
    for (let k = 0; k < 4; k++) {
      const f = layerForces(L, unit[k][0], unit[k][1], unit[k][2], unit[k][3], dT[i]);
      const col = cols[k];
      if (ctl.axial.force !== undefined) A[0][col] += f.N - f0.N;
      if (ctl.torsion.torque !== undefined) A[1][col] += f.M - f0.M;
      A[iu(i)][col] += f.Nth - f0.Nth;
      A[id(i)][col] += f.thick - f0.thick;
    }
    if (ctl.axial.force !== undefined) b[0] -= f0.N;
    if (ctl.torsion.torque !== undefined) b[1] -= f0.M;
    b[iu(i)] = fIn[i] * L.rIn - fOut[i] * L.rOut - f0.Nth;
    b[id(i)] = -f0.thick;
    if (i > 0) A[iu(i)][ip(i - 1)] -= L.rIn;
    if (i < N - 1) A[iu(i)][ip(i)] += L.rOut;
    if (L.kind === 'shell') {
      b[id(i)] -= (fIn[i] * L.rIn + fOut[i] * L.rOut) / 2;
      if (i > 0) A[id(i)][ip(i - 1)] += L.rIn / 2;
      if (i < N - 1) A[id(i)][ip(i)] += L.rOut / 2;
    }
  }
  if (ctl.axial.force !== undefined) b[0] += ctl.axial.force; else { A[0][0] = 1; b[0] = ctl.axial.strain; }
  if (ctl.torsion.torque !== undefined) b[1] += ctl.torsion.torque; else { A[1][1] = 1; b[1] = ctl.torsion.twist; }
  for (let j = 0; j < N - 1; j++) {
    const r = ip(j);
    if (active[j]) {
      A[r][iu(j + 1)] = 1; A[r][id(j + 1)] = -0.5; A[r][iu(j)] = -1; A[r][id(j)] = -0.5;
      b[r] = -gap0[j];
    } else { A[r][r] = 1; b[r] = 0; }
  }
  // equilibrate: unknowns span ~15 orders of magnitude (strain vs pressure)
  const cs = new Array(m).fill(0), rs = new Array(m).fill(0);
  for (let c = 0; c < m; c++) { for (let r = 0; r < m; r++) cs[c] = Math.max(cs[c], Math.abs(A[r][c])); if (cs[c] === 0) cs[c] = 1; }
  for (let r = 0; r < m; r++) {
    for (let c = 0; c < m; c++) { A[r][c] /= cs[c]; rs[r] = Math.max(rs[r], Math.abs(A[r][c])); }
    if (rs[r] === 0) rs[r] = 1;
    for (let c = 0; c < m; c++) A[r][c] /= rs[r];
    b[r] /= rs[r];
  }
  const y = solveDense(A, b);
  const x = new Array(m);
  for (let c = 0; c < m; c++) {
    x[c] = y[c] / cs[c];
    if (!Number.isFinite(x[c])) throw new Error('Flexible-pipe axisymmetric system is singular for this layer stack');
  }
  const u = [], d = [], pc = [], gap = [];
  for (let i = 0; i < N; i++) { u.push(x[iu(i)]); d.push(x[id(i)]); }
  for (let j = 0; j < N - 1; j++) {
    pc.push(active[j] ? x[ip(j)] : 0);
    gap.push(gap0[j] + (u[j + 1] - d[j + 1] / 2) - (u[j] + d[j] / 2));
  }
  return { eps: x[0], tau: x[1], u, d, pc, gap };
}

function normaliseLoads(stack, loads) {
  const N = stack.layers.length;
  const o = { tension: 0, tensionType: 'effective', pInt: 0, pExt: 0, annulus: 'dry', pAnnulus: 0, torque: 0, ...loads };
  fin(o.tension, 'tension'); nonneg(o.pInt, 'pInt'); nonneg(o.pExt, 'pExt'); nonneg(o.pAnnulus, 'pAnnulus'); fin(o.torque, 'torque');
  if (o.annulus !== 'dry' && o.annulus !== 'flooded') throw new RangeError("annulus must be 'dry' or 'flooded'");
  if (o.tensionType !== 'effective' && o.tensionType !== 'wall') throw new RangeError("tensionType must be 'effective' or 'wall'");
  if (o.axialStrain !== undefined) fin(o.axialStrain, 'axialStrain');
  if (o.twist !== undefined) fin(o.twist, 'twist');
  let dT = o.deltaT ?? 0;
  if (typeof dT === 'number') dT = new Array(N).fill(fin(dT, 'deltaT'));
  else if (!Array.isArray(dT) || dT.length !== N) throw new RangeError('deltaT must be a number or an array with one entry per layer');
  let gap0 = o.initialGaps ?? new Array(N - 1).fill(0);
  if (!Array.isArray(gap0) || gap0.length !== N - 1) throw new RangeError('initialGaps must have one entry per interface');
  gap0.forEach((g, j) => nonneg(g, `initialGaps[${j}]`));
  return { ...o, dT, gap0, pAnn: o.annulus === 'flooded' ? o.pExt : o.pAnnulus };
}

function sectionMatrix(stack, active) {
  const N = stack.layers.length;
  const zeros = new Array(N).fill(0), g0 = new Array(Math.max(N - 1, 0)).fill(0);
  const col = (strain, twist) => {
    const s = solveLinear(stack, { axial: { strain }, torsion: { twist } }, zeros, zeros, zeros, active, g0);
    let F = 0, M = 0;
    stack.layers.forEach((L, i) => { const f = layerForces(L, s.eps, s.tau, s.u[i], s.d[i], 0); F += f.N; M += f.M; });
    return [F, M];
  };
  const a = col(1, 0), t = col(0, 1);
  return { K11: a[0], K21: a[1], K12: t[0], K22: t[1] };
}

/**
 * Axisymmetric response of the layer stack to tension, internal/external pressure, torque, temperature rise and
 * (optionally) prescribed axial strain or twist. See the file header for the formulation.
 *
 * Reference: first-principles layer statics as in Feret & Bournazel (1987), Witz & Tan (1992); effective tension
 * per Sparks (1984) / API RP 17B.
 * Validity: linear elastic, small strain; contact unilateral (no adhesion); friction does not enter the
 * axisymmetric response.
 * @param {object} stack result of {@link buildFlexibleStack}
 * @param {object} [loads]
 * @param {number} [loads.tension=0] effective tension (or wall force, see tensionType), N
 * @param {'effective'|'wall'} [loads.tensionType='effective'] 'wall' applies the value directly as the sum of the
 *   layer forces (no end-cap terms)
 * @param {number} [loads.pInt=0] bore pressure, Pa
 * @param {number} [loads.pExt=0] external pressure, Pa
 * @param {'dry'|'flooded'} [loads.annulus='dry'] flooded: annulus at external pressure
 * @param {number} [loads.pAnnulus=0] annulus gas pressure for a dry annulus, Pa
 * @param {number} [loads.torque=0] applied torque, N m (positive tightens lay = +1 layers)
 * @param {number} [loads.axialStrain] prescribed axial strain (replaces the tension condition), -
 * @param {number} [loads.twist] prescribed twist per unit length (replaces the torque condition), rad/m; use 0 for
 *   torsionally restrained ends
 * @param {number|number[]} [loads.deltaT=0] temperature rise of all layers or per layer, K
 * @param {number[]} [loads.initialGaps] initial radial gap at each interface, m (default 0)
 * @returns {object} { axialStrain [-], twist [rad/m], wallForce [N] (sum of layer axial forces), endCapForce [N],
 *   tensionEffective [N], torque [N m], layers: [{ name, type, kind, axialForce [N], torque [N m],
 *   hoopForce [N/m], radialDisp [m], thicknessChange [m], pIn, pOut [Pa] (total face pressures), contactIn,
 *   contactOut [Pa], wireStrain [-], wireStress [Pa] (helical) or sigmaAxial, sigmaHoop, sigmaRadial, shear,
 *   vonMises [Pa] (shell) }], interfaces: [{ contactPressure [Pa], gap [m], closed [bool] }],
 *   stiffness: { K11 [N], K12, K21 [N m], K22 [N m^2] } (tangent, [force; torque] = K [strain; twist]),
 *   residuals: { axial [N], torque [N m], radialMax [N/m], gapMax [m] }, iterations }
 */
export function axisymmetricResponse(stack, loads = {}) {
  const o = normaliseLoads(stack, loads);
  const Ls = stack.layers, N = Ls.length;
  const { fIn, fOut } = fluidOnFaces(stack, o.pInt, o.pExt, o.pAnn);
  let endCap = 0;
  if (o.tensionType === 'effective') for (const L of Ls) if (L.seal) endCap += Math.PI * (fIn[L.index] * L.rIn ** 2 - fOut[L.index] * L.rOut ** 2);
  const ctl = {
    axial: o.axialStrain !== undefined ? { strain: o.axialStrain } : { force: o.tension + endCap },
    torsion: o.twist !== undefined ? { twist: o.twist } : { torque: o.torque },
  };
  const active = new Array(N - 1).fill(true);
  for (let j = 0; j < N - 1; j++) if (o.gap0[j] > 0) active[j] = false;
  let sol, it = 0, done = false;
  const pScale = Math.max(o.pInt, o.pExt, Math.abs(o.tension) / (Math.PI * stack.outerRadius ** 2), 1);
  for (; it < 200 && !done; it++) {
    sol = solveLinear(stack, ctl, fIn, fOut, o.dT, active, o.gap0);
    const flips = [];
    for (let j = 0; j < N - 1; j++) {
      if (active[j] && sol.pc[j] < -1e-10 * pScale) flips.push({ j, v: -sol.pc[j] / pScale });
      else if (!active[j] && sol.gap[j] < -1e-13 * stack.outerRadius) flips.push({ j, v: -sol.gap[j] / stack.outerRadius * 1e3 });
    }
    if (flips.length === 0) { done = true; break; }
    // change all violated constraints at once first; fall back to single changes if that cycles
    if (it < 30) for (const f of flips) active[f.j] = !active[f.j];
    else { flips.sort((a, c) => c.v - a.v); active[flips[0].j] = !active[flips[0].j]; }
  }
  if (!done) throw new Error('Flexible-pipe contact iteration did not converge');

  const layers = [];
  let wall = 0, tq = 0, radialMax = 0;
  for (let i = 0; i < N; i++) {
    const L = Ls[i];
    const f = layerForces(L, sol.eps, sol.tau, sol.u[i], sol.d[i], o.dT[i]);
    const cIn = i > 0 ? sol.pc[i - 1] : 0, cOut = i < N - 1 ? sol.pc[i] : 0;
    const pIn = cIn + fIn[i], pOut = cOut + fOut[i];
    wall += f.N; tq += f.M;
    radialMax = Math.max(radialMax, Math.abs(f.Nth - (pIn * L.rIn - pOut * L.rOut)));
    const rec = {
      name: L.name, type: L.type, kind: L.kind, axialForce: f.N, torque: f.M, hoopForce: f.Nth,
      radialDisp: sol.u[i], thicknessChange: sol.d[i], pIn, pOut, contactIn: cIn, contactOut: cOut,
    };
    if (L.kind === 'helical') { rec.wireStrain = f.ew; rec.wireStress = L.E * f.ew; } else {
      rec.sigmaAxial = f.sz; rec.sigmaHoop = f.st; rec.sigmaRadial = f.sr; rec.shear = f.shear;
      rec.vonMises = Math.sqrt(0.5 * ((f.sz - f.st) ** 2 + (f.st - f.sr) ** 2 + (f.sr - f.sz) ** 2) + 3 * f.shear ** 2);
    }
    layers.push(rec);
  }
  let gapMax = 0;
  const interfaces = sol.pc.map((p, j) => {
    if (active[j]) gapMax = Math.max(gapMax, Math.abs(sol.gap[j]));
    return { contactPressure: p, gap: active[j] ? 0 : sol.gap[j], closed: active[j] };
  });
  return {
    axialStrain: sol.eps, twist: sol.tau, wallForce: wall, endCapForce: endCap,
    tensionEffective: wall - endCap, torque: tq, layers, interfaces,
    stiffness: sectionMatrix(stack, active),
    residuals: {
      axial: ctl.axial.force !== undefined ? wall - ctl.axial.force : 0,
      torque: ctl.torsion.torque !== undefined ? tq - ctl.torsion.torque : 0,
      radialMax, gapMax,
    },
    iterations: it + 1,
  };
}

function layerCapacity(L, basis) {
  const v = basis === 'yield' ? L.yield : basis === 'uts' ? L.uts : (L.capacity ?? (L.uts !== undefined ? 0.9 * L.uts : L.yield));
  return v;
}

/**
 * Utilisation of every metallic layer against the permissible factors of API Spec 17J Table 6
 * ({@link API17J_FACTORS}). Stress utilisation = |wire stress| / (factor x structural capacity), where the
 * structural capacity of a layer is its `capacity` field, else 0.9 x uts, else yield (17J defines the structural
 * capacity as yield or 0.9 x UTS at the manufacturer's choice). Carcass and layers without strength data are
 * reported with utilisation null (carcass collapse is checked in flexcollapse.js).
 * @param {object} stack layer stack
 * @param {object} response result of {@link axisymmetricResponse}
 * @param {object} [o]
 * @param {'normal'|'abnormal'|'hydrotest'} [o.condition='normal'] load condition
 * @param {object} [o.factors] override of API17J_FACTORS (same shape)
 * @returns {{condition:string, layers:{name:string, type:string, stress:number, capacity:number, factor:number,
 *   utilisation:number|null}[], governing:{name:string, utilisation:number}|null}} stresses Pa
 */
export function layerUtilisations(stack, response, { condition = 'normal', factors = API17J_FACTORS } = {}) {
  if (!['normal', 'abnormal', 'hydrotest'].includes(condition)) throw new RangeError("condition must be 'normal', 'abnormal' or 'hydrotest'");
  const rows = [];
  let governing = null;
  stack.layers.forEach((L, i) => {
    if (L.kind !== 'helical') return;
    const group = L.type === 'tensileArmour' ? factors.tensileArmour : (L.type === 'pressureArmour' || L.type === 'backupArmour') ? factors.pressureArmour : null;
    const cap = layerCapacity(L);
    const stress = response.layers[i].wireStress;
    const factor = group ? group[condition] : null;
    const utilisation = group && cap ? Math.abs(stress) / (factor * cap) : null;
    rows.push({ name: L.name, type: L.type, stress, capacity: cap ?? null, factor, utilisation });
    if (utilisation !== null && (!governing || utilisation > governing.utilisation)) governing = { name: L.name, utilisation };
  });
  return { condition, layers: rows, governing };
}

function scaleToLimit(stack, makeLoads, types, basis) {
  // stresses are piecewise linear in the load multiplier; iterate so the contact state is that at failure
  let lam = 1, res = null;
  for (let k = 0; k < 12; k++) {
    const r = axisymmetricResponse(stack, makeLoads(lam));
    let worst = 0, name = null;
    stack.layers.forEach((L, i) => {
      if (!types.includes(L.type)) return;
      const cap = layerCapacity(L, basis);
      if (!cap) return;
      const ut = Math.abs(r.layers[i].wireStress) / cap;
      if (ut > worst) { worst = ut; name = L.name; }
    });
    if (!(worst > 0)) return null;
    const next = lam / worst;
    res = { multiplier: next, layer: name };
    if (Math.abs(next - lam) <= 1e-9 * Math.abs(next)) break;
    lam = next;
  }
  return res;
}

/**
 * Burst pressure of the pipe with closed ends (end-cap load carried by the tensile armour): the internal pressure
 * at which the pressure armour (hoop failure mode) and the tensile armour (axial failure mode) wire stresses
 * reach the layer ultimate strength; both are reported and the lower governs.
 * Reference: API RP 17B §8 (failure modes: burst by rupture of pressure or tensile armour).
 * Validity: elastic stress distribution up to failure (no plastic load redistribution between pressure-armour
 * layers), so the value is a conservative first-yield type estimate when uts is replaced by yield.
 * @param {object} stack layer stack; armour layers need `uts` (Pa)
 * @param {object} [o]
 * @param {number} [o.pExt=0] external pressure, Pa
 * @param {number} [o.tension=0] effective tension applied together with the pressure, N
 * @param {'uts'|'yield'} [o.basis='uts'] strength used
 * @returns {{pressureArmour:number|null, tensileArmour:number|null, burst:number, governing:string}} Pa
 */
export function burstPressure(stack, { pExt = 0, tension = 0, basis = 'uts' } = {}) {
  const p0 = 1e7;
  const run = (types) => scaleToLimit(stack, (lam) => ({ pInt: lam * p0, pExt, tension }), types, basis);
  const pa = run(['pressureArmour', 'backupArmour']);
  const ta = run(['tensileArmour']);
  const pPa = pa ? pa.multiplier * p0 : null, pTa = ta ? ta.multiplier * p0 : null;
  if (pPa === null && pTa === null) throw new RangeError('stack has no armour layer with strength data (uts/yield)');
  const burst = Math.min(pPa ?? Infinity, pTa ?? Infinity);
  return { pressureArmour: pPa, tensileArmour: pTa, burst, governing: burst === pPa ? 'pressureArmour' : 'tensileArmour' };
}

/**
 * Axial failure tension: effective tension at which the most loaded tensile-armour wire reaches its strength
 * (other loads held at the given values).
 * Reference: API RP 17B §8 (tensile failure of the tensile armour).
 * @param {object} stack layer stack
 * @param {object} [o]
 * @param {number} [o.pInt=0] Pa
 * @param {number} [o.pExt=0] Pa
 * @param {'uts'|'yield'} [o.basis='uts']
 * @param {number} [o.twist] prescribed twist (0 = torsionally restrained ends), rad/m; default torque-free
 * @returns {{tension:number, layer:string}} N
 */
export function axialFailureTension(stack, { pInt = 0, pExt = 0, basis = 'uts', twist } = {}) {
  const t0 = 1e6;
  const r = scaleToLimit(stack, (lam) => ({ tension: lam * t0, pInt, pExt, twist }), ['tensileArmour'], basis);
  if (!r) throw new RangeError('stack has no tensile armour layer with strength data');
  return { tension: r.multiplier * t0, layer: r.layer };
}

/**
 * Torsional limits about a baseline load state. Torque is applied in the direction that tightens the outermost
 * tensile armour (lay-tightening) and in the opposite direction (lay-opening). For each direction two limits are
 * reported: the torque at which a tensile-armour wire reaches factor x structural capacity, and the torque at
 * which an interface adjacent to a tensile-armour layer first opens (loss of radial support: onset of layer
 * separation / bird-caging tendency in the lay-opening direction). The lower of the two governs.
 * Reference: API RP 17B §8 (torsional failure modes: wire failure, bird-caging); solved with
 * {@link axisymmetricResponse}.
 * Validity: axisymmetric (no wire buckling); separation limit is a gap-onset criterion, not a post-buckling
 * capacity.
 * @param {object} stack layer stack
 * @param {object} [o]
 * @param {object} [o.loads] baseline loads (tension, pInt, pExt, annulus ...)
 * @param {number} [o.factor=0.67] permissible stress utilisation of the tensile armour
 * @returns {{tightening:object, opening:object}} each { stressLimit [N m]|null, separation [N m]|null,
 *   governing [N m], mode: 'stress'|'separation' }, signed torques
 */
export function torsionLimits(stack, { loads = {}, factor = 0.67 } = {}) {
  pos(factor, 'factor');
  const idx = stack.layers.filter((L) => L.type === 'tensileArmour').map((L) => L.index);
  if (idx.length === 0) throw new RangeError('stack has no tensile armour layer');
  const outer = stack.layers[idx[idx.length - 1]];
  const state = (tq) => {
    const r = axisymmetricResponse(stack, { ...loads, torque: tq, twist: undefined });
    let ut = 0, open = false;
    for (const i of idx) {
      const cap = layerCapacity(stack.layers[i]);
      if (cap) ut = Math.max(ut, Math.abs(r.layers[i].wireStress) / (factor * cap));
      if (i > 0 && !r.interfaces[i - 1].closed) open = true;
      if (i < stack.layers.length - 1 && !r.interfaces[i].closed) open = true;
    }
    return { ut, open };
  };
  const base = state(0);
  const solveDir = (sign) => {
    const find = (pred) => {
      if (pred(base)) return 0;
      let hi = 1e3;
      for (let k = 0; k < 50 && !pred(state(sign * hi)); k++) hi *= 2;
      if (!pred(state(sign * hi))) return null;
      let lo = 0;
      for (let k = 0; k < 60; k++) { const mid = 0.5 * (lo + hi); if (pred(state(sign * mid))) hi = mid; else lo = mid; }
      return sign * hi;
    };
    const stressLimit = find((s) => s.ut >= 1);
    const separation = find((s) => s.open);
    const cands = [['stress', stressLimit], ['separation', separation]].filter((c) => c[1] !== null);
    cands.sort((a, c) => Math.abs(a[1]) - Math.abs(c[1]));
    return { stressLimit, separation, governing: cands.length ? cands[0][1] : null, mode: cands.length ? cands[0][0] : null };
  };
  return { tightening: solveDir(outer.lay), opening: solveDir(-outer.lay) };
}

/**
 * Bending properties of the stack: stick and slip bending stiffness, and for every tensile-armour layer the
 * slip-onset curvature and friction moment from the interface pressures of an axisymmetric solution.
 * Shell layers contribute E pi (ro^4 - ri^4)/4; helical layers contribute per
 * helix.js `helixBendingStiffness` (wire section from wireWidth x wireThickness, default thickness = layer
 * thickness; without wireWidth the wire's own bending stiffness is neglected) and `criticalCurvature` with
 * f = mu (q_in + q_out) x tributary width.
 * Reference: Saevik (2011) Marine Structures 24, 97-111; Witz & Tan (1992).
 * Validity: see helix.js; contact pressures taken as uniform round the circumference.
 * @param {object} stack layer stack
 * @param {object} [o]
 * @param {object} [o.response] axisymmetric solution supplying interface pressures (omit for stiffness only)
 * @param {number} [o.mu=0.1] interlayer friction coefficient (typical 0.05-0.2 for steel on anti-wear tape; no
 *   normative value), -
 * @param {'loxodromic'|'geodesic'} [o.curve='loxodromic']
 * @returns {{EIstick:number, EIslip:number, frictionMoment:number, kappaCrMin:number|null, layers:object[]}}
 *   N m^2; N m (sum over armour layers); 1/m (first slip); per armour layer { name, kappaCr, kappaFullSlip
 *   [1/m], sigmaFrictionMax [Pa], frictionMoment [N m], EIstick, EIslip [N m^2] }
 */
export function bendingProperties(stack, { response, mu = 0.1, curve = 'loxodromic' } = {}) {
  nonneg(mu, 'mu');
  let EIstick = 0, EIslip = 0, Mf = 0, kMin = null;
  const rows = [];
  stack.layers.forEach((L, i) => {
    if (L.kind === 'shell') {
      const ei = L.E * Math.PI * (L.rOut ** 4 - L.rIn ** 4) / 4;
      EIstick += ei; EIslip += ei;
      return;
    }
    const n = L.n ?? 1, A = L.nA / n;
    let sec = { In: 0, It: 0, J: 0 };
    if (L.wireWidth) sec = rectWireSection({ width: L.wireWidth, thickness: L.wireThickness ?? L.t });
    const scale = L.n ? 1 : 0; // without a wire count the individual-wire stiffness is unknown
    const b = helixBendingStiffness({ n, E: L.E, A, R: L.R, alpha: L.alpha, In: sec.In * scale, It: sec.It * scale, J: sec.J * scale, curve });
    EIstick += b.stick; EIslip += b.slip;
    if (L.type !== 'tensileArmour') return;
    const row = { name: L.name, EIstick: b.stick, EIslip: b.slip };
    if (response) {
      const qIn = Math.max(response.layers[i].contactIn, 0), qOut = Math.max(response.layers[i].contactOut, 0);
      const cc = criticalCurvature({ E: L.E, A, R: L.R, alpha: L.alpha, n, mu, qInner: qIn, qOuter: qOut });
      Object.assign(row, { kappaCr: cc.kappaCr, kappaFullSlip: cc.kappaFullSlip, sigmaFrictionMax: cc.sigmaFrictionMax, frictionMoment: cc.frictionMoment, frictionPerLength: cc.frictionPerLength });
      Mf += cc.frictionMoment;
      kMin = kMin === null ? cc.kappaCr : Math.min(kMin, cc.kappaCr);
    }
    rows.push(row);
  });
  return { EIstick, EIslip, frictionMoment: Mf, kappaCrMin: kMin, layers: rows };
}

/**
 * Minimum bend radius from locking of the interlocked layers and from polymer sheath strain.
 *
 * Locking (gap-closure geometry): an interlocked layer of mean pitch Lp can shorten on the intrados until the
 * profile gaps close and lengthen on the extrados until the hooks engage; with lockStrain = (free axial play per
 * pitch)/Lp the layer locks when the bending strain at its outer fibre r_out/rho reaches lockStrain:
 *     rho_lock = r_out / lockStrain.
 * Sheath strain: bending strain r_out/rho <= strain limit ({@link SHEATH_STRAIN_LIMITS}; per-layer fields
 * strainLimit / strainLimitDynamic).
 * API Spec 17J §5.3.1 (as recalled): storage MBR >= 1.1 x locking radius; static service 1.0 x storage MBR;
 * dynamic service 1.5 x storage MBR. All three factors are arguments.
 * Validity: geometric; ignores ovalisation and layer interaction.
 * @param {object} stack layers may carry lockStrain, strainLimit, strainLimitDynamic
 * @param {object} [o]
 * @param {number} [o.storageFactor=1.1] on locking radius
 * @param {number} [o.staticFactor=1.0] on storage MBR
 * @param {number} [o.dynamicFactor=1.5] on storage MBR
 * @returns {{lockingRadius:number|null, lockingLayer:string|null, sheathStaticRadius:number|null,
 *   sheathDynamicRadius:number|null, storage:number, static:number, dynamic:number}} m
 */
export function minimumBendRadius(stack, { storageFactor = 1.1, staticFactor = 1.0, dynamicFactor = 1.5 } = {}) {
  pos(storageFactor, 'storageFactor'); pos(staticFactor, 'staticFactor'); pos(dynamicFactor, 'dynamicFactor');
  let lock = null, lockLayer = null, shS = null, shD = null;
  for (const L of stack.layers) {
    if (L.lockStrain !== undefined) {
      const r = L.rOut / pos(L.lockStrain, `${L.name}.lockStrain`);
      if (lock === null || r > lock) { lock = r; lockLayer = L.name; }
    }
    if (L.strainLimit !== undefined) shS = Math.max(shS ?? 0, L.rOut / pos(L.strainLimit, `${L.name}.strainLimit`));
    const dyn = L.strainLimitDynamic ?? L.strainLimit;
    if (dyn !== undefined) shD = Math.max(shD ?? 0, L.rOut / pos(dyn, `${L.name}.strainLimitDynamic`));
  }
  if (lock === null && shS === null) throw new RangeError('no layer carries lockStrain or strainLimit');
  const storage = Math.max(lock !== null ? storageFactor * lock : 0, shS ?? 0);
  return {
    lockingRadius: lock, lockingLayer: lockLayer, sheathStaticRadius: shS, sheathDynamicRadius: shD,
    storage, static: staticFactor * storage, dynamic: Math.max(dynamicFactor * storage, shD ?? 0),
  };
}

/**
 * Linear mass and weights. Solid areas: shell 2 pi R t; helical layer nA/cos(alpha) (wire area projected on the
 * cross-section). The bore fluid fills the bore and the voids of non-sealing layers inside the inner seal; a
 * flooded annulus adds seawater in the voids of layers between the seals.
 * Reference: Archimedes; API RP 17B weight definitions (empty / seawater-filled, in air / submerged).
 * @param {object} stack layer stack (layers need rho)
 * @param {object} [o]
 * @param {number} [o.rhoContent=0] bore fluid density in operation, kg/m^3
 * @param {number} [o.rhoSeawater=1025] kg/m^3
 * @param {number} [o.g=9.80665] m/s^2
 * @param {boolean} [o.annulusFlooded=false] include seawater in annulus voids
 * @returns {{massEmpty:number, massSeawaterFilled:number, massOperating:number, annulusWaterMass:number,
 *   contentArea:number, annulusVoidArea:number, displacedMass:number, weightAirEmpty:number,
 *   submergedEmpty:number, submergedFlooded:number, submergedOperating:number}} kg/m, m^2, N/m (positive down)
 */
export function massAndWeight(stack, { rhoContent = 0, rhoSeawater = 1025, g = 9.80665, annulusFlooded = false } = {}) {
  nonneg(rhoContent, 'rhoContent'); pos(rhoSeawater, 'rhoSeawater'); pos(g, 'g');
  let m = 0, solidInside = 0, annVoid = 0;
  for (const L of stack.layers) {
    const ring = Math.PI * (L.rOut ** 2 - L.rIn ** 2);
    const solid = L.kind === 'helical' ? L.nA / Math.cos(L.alpha) : ring;
    m += L.rho * solid;
    if (L.index < stack.sealInner) solidInside += solid;
    if (L.index > stack.sealInner && L.index < stack.sealOuter) annVoid += ring - solid;
  }
  const contentArea = Math.PI * stack.layers[stack.sealInner].rIn ** 2 - solidInside;
  const annulusWaterMass = annulusFlooded ? rhoSeawater * annVoid : 0;
  const displaced = rhoSeawater * Math.PI * stack.outerRadius ** 2;
  const mSw = m + rhoSeawater * contentArea + annulusWaterMass;
  const mOp = m + rhoContent * contentArea + annulusWaterMass;
  return {
    massEmpty: m, massSeawaterFilled: mSw, massOperating: mOp, annulusWaterMass, contentArea, annulusVoidArea: annVoid,
    displacedMass: displaced, weightAirEmpty: m * g,
    submergedEmpty: (m + annulusWaterMass - displaced) * g, submergedFlooded: (mSw - displaced) * g, submergedOperating: (mOp - displaced) * g,
  };
}

/**
 * Equivalent beam properties for global (riser / flowline) analysis: axial stiffness with free and restrained
 * twist, torsional stiffness, coupling, stick and slip bending stiffness, mass and weights.
 * EA_free = K11 - K12 K21/K22 (torque-free ends), EA_fixed = K11; GJ likewise.
 * @param {object} stack layer stack
 * @param {object} [o]
 * @param {object} [o.loads] load state about which the tangent stiffness is taken (contact state dependent)
 * @param {number} [o.mu=0.1] friction coefficient for the slip-onset curvature
 * @param {object} [o.weight] options forwarded to {@link massAndWeight}
 * @returns {object} { EA, EAfixed [N], GJ, GJfixed [N m^2], coupling [N m], EIstick, EIslip [N m^2],
 *   kappaCr [1/m], frictionMoment [N m], outerDiameter, innerDiameter [m], mass: massAndWeight result }
 */
export function equivalentBeamProperties(stack, { loads = {}, mu = 0.1, weight = {} } = {}) {
  const r = axisymmetricResponse(stack, loads);
  const K = r.stiffness;
  const bend = bendingProperties(stack, { response: r, mu });
  return {
    EA: K.K11 - K.K12 * K.K21 / K.K22, EAfixed: K.K11,
    GJ: K.K22 - K.K12 * K.K21 / K.K11, GJfixed: K.K22,
    coupling: 0.5 * (K.K12 + K.K21),
    EIstick: bend.EIstick, EIslip: bend.EIslip, kappaCr: bend.kappaCrMin, frictionMoment: bend.frictionMoment,
    outerDiameter: 2 * stack.outerRadius, innerDiameter: 2 * stack.innerRadius,
    mass: massAndWeight(stack, weight),
  };
}

/**
 * Steady radial heat conduction through the layers with internal and external film coefficients.
 * Per unit length: R_layer = ln(ro/ri)/(2 pi k), R_film = 1/(2 pi r h); q = (Tin - Tout)/sum R; overall
 * coefficient referenced to the bore area U_i = 1/(2 pi r_i sum R).
 * Reference: Fourier conduction in a composite cylinder (Incropera et al., Fundamentals of Heat and Mass
 * Transfer, §3.3).
 * Validity: steady state, constant conductivities; helical layers need an effective radial conductivity `k`.
 * @param {object} stack layer stack (layers need k, W/m/K)
 * @param {object} o
 * @param {number} o.TinC bulk fluid temperature, degC
 * @param {number} o.ToutC ambient temperature, degC
 * @param {number} [o.hIn=Infinity] internal film coefficient, W/m^2/K
 * @param {number} [o.hOut=Infinity] external film coefficient, W/m^2/K
 * @returns {{heatLoss:number, UvalueInner:number, UvalueOuter:number, resistance:number,
 *   interfaceTempsC:number[], layerMeanTempsC:number[]}} W/m; W/m^2/K (bore / outer area); m K/W; temperatures at
 *   the N+1 layer boundaries (bore wall first) and layer means, degC
 */
export function thermalProfile(stack, { TinC, ToutC, hIn = Infinity, hOut = Infinity }) {
  fin(TinC, 'TinC'); fin(ToutC, 'ToutC');
  if (!(hIn > 0)) throw new RangeError('hIn must be positive'); if (!(hOut > 0)) throw new RangeError('hOut must be positive');
  const Rs = stack.layers.map((L) => Math.log(L.rOut / L.rIn) / (TWO_PI * pos(L.k, `${L.name}.k`)));
  const rFilmIn = 1 / (TWO_PI * stack.innerRadius * hIn), rFilmOut = 1 / (TWO_PI * stack.outerRadius * hOut);
  const total = rFilmIn + rFilmOut + Rs.reduce((a, c) => a + c, 0);
  const q = (TinC - ToutC) / total;
  const temps = [TinC - q * rFilmIn];
  for (const R of Rs) temps.push(temps[temps.length - 1] - q * R);
  return {
    heatLoss: q, resistance: total,
    UvalueInner: 1 / (TWO_PI * stack.innerRadius * total), UvalueOuter: 1 / (TWO_PI * stack.outerRadius * total),
    interfaceTempsC: temps, layerMeanTempsC: Rs.map((_, i) => 0.5 * (temps[i] + temps[i + 1])),
  };
}

/**
 * Pressure- and temperature-induced elongation and twist of a free pipe (zero effective tension, zero torque,
 * closed ends), from {@link axisymmetricResponse}.
 * @param {object} stack layer stack (alphaT per layer for thermal strain)
 * @param {object} [o]
 * @param {number} [o.pInt=0] Pa
 * @param {number} [o.pExt=0] Pa
 * @param {number|number[]} [o.deltaT=0] temperature rise, K (per layer or uniform)
 * @param {'dry'|'flooded'} [o.annulus='dry']
 * @returns {{strainPressure:number, strainThermal:number, strainTotal:number, twistPressure:number,
 *   twistThermal:number}} strains - (m/m), twists rad/m; pressure and thermal parts solved separately
 */
export function pressureThermalElongation(stack, { pInt = 0, pExt = 0, deltaT = 0, annulus = 'dry' } = {}) {
  const p = axisymmetricResponse(stack, { pInt, pExt, annulus });
  const t = axisymmetricResponse(stack, { deltaT });
  const all = axisymmetricResponse(stack, { pInt, pExt, annulus, deltaT });
  return { strainPressure: p.axialStrain, strainThermal: t.axialStrain, strainTotal: all.axialStrain, twistPressure: p.twist, twistThermal: t.twist };
}

const STEEL = { E: 207e9, nu: 0.3, rho: 7850, k: 45, alphaT: 1.17e-5 };
const SS316 = { E: 193e9, nu: 0.3, rho: 7980, k: 15, alphaT: 1.6e-5 };
const PA11 = { E: 350e6, nu: 0.42, rho: 1040, k: 0.25, alphaT: 1.1e-4 };
const HDPE = { E: 800e6, nu: 0.42, rho: 950, k: 0.4, alphaT: 1.8e-4 };
const TAPE = { E: 350e6, nu: 0.42, rho: 1040, k: 0.25, alphaT: 1.1e-4 };

/**
 * Generic demonstration designs. These are illustrative layer stacks with round-number dimensions and textbook
 * material constants; they are NOT manufacturer designs and must not be used for engineering decisions.
 *  - 'riser-6in'    generic 6-inch (152.4 mm bore) rough-bore dynamic production riser
 *  - 'flowline-10in' generic 10-inch (254 mm bore) rough-bore static flowline
 * @param {'riser-6in'|'flowline-10in'} id
 * @returns {{id:string, label:string, generic:true, innerDiameter:number, layers:object[]}} definition accepted
 *   by {@link buildFlexibleStack} (innerDiameter m)
 */
export function exampleFlexibleDesign(id) {
  const armour = (name, alphaDeg, lay, t, w) => ({ name, type: 'tensileArmour', thickness: t, alphaDeg, lay, fill: 0.9, wireWidth: w, wireThickness: t, ...STEEL, yield: 1100e6, uts: 1400e6 });
  const tape = (name) => ({ name, type: 'antiWear', thickness: 0.0015, ...TAPE });
  if (id === 'riser-6in') {
    return {
      id, label: 'Generic 6-inch dynamic production riser (demonstration only)', generic: true, innerDiameter: 0.1524,
      layers: [
        { name: 'Carcass', type: 'carcass', thickness: 0.007, alphaDeg: 87.5, lay: 1, fill: 0.55, ...SS316, yield: 300e6, uts: 600e6, lockStrain: 0.09 },
        { name: 'Pressure sheath', type: 'pressureSheath', thickness: 0.008, ...PA11, strainLimit: 0.077, strainLimitDynamic: 0.077 },
        { name: 'Pressure armour', type: 'pressureArmour', thickness: 0.008, alphaDeg: 87, lay: 1, fill: 0.9, ...STEEL, yield: 800e6, uts: 1000e6, lockStrain: 0.07 },
        tape('Anti-wear tape 1'),
        armour('Tensile armour 1', 35, 1, 0.004, 0.012),
        tape('Anti-wear tape 2'),
        armour('Tensile armour 2', 35, -1, 0.004, 0.012),
        { name: 'Outer sheath', type: 'outerSheath', thickness: 0.008, ...HDPE, strainLimit: 0.077, strainLimitDynamic: 0.077 },
      ],
    };
  }
  if (id === 'flowline-10in') {
    return {
      id, label: 'Generic 10-inch static flowline (demonstration only)', generic: true, innerDiameter: 0.254,
      layers: [
        { name: 'Carcass', type: 'carcass', thickness: 0.010, alphaDeg: 87.5, lay: 1, fill: 0.55, ...SS316, yield: 300e6, uts: 600e6, lockStrain: 0.09 },
        { name: 'Pressure sheath', type: 'pressureSheath', thickness: 0.010, ...HDPE, strainLimit: 0.077 },
        { name: 'Pressure armour', type: 'pressureArmour', thickness: 0.010, alphaDeg: 87, lay: 1, fill: 0.9, ...STEEL, yield: 800e6, uts: 1000e6, lockStrain: 0.07 },
        tape('Anti-wear tape 1'),
        armour('Tensile armour 1', 30, 1, 0.005, 0.012),
        tape('Anti-wear tape 2'),
        armour('Tensile armour 2', 30, -1, 0.005, 0.012),
        { name: 'Insulation', type: 'insulation', thickness: 0.012, E: 200e6, nu: 0.4, rho: 700, k: 0.15, alphaT: 1.5e-4 },
        { name: 'Outer sheath', type: 'outerSheath', thickness: 0.010, ...HDPE, strainLimit: 0.077 },
      ],
    };
  }
  throw new RangeError("id must be 'riser-6in' or 'flowline-10in'");
}
