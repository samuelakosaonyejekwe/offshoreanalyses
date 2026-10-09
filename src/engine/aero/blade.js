/**
 * Blade (and tubular tower section) structural models: rotating cantilever beam finite elements, Campbell data,
 * static deflection and tip clearance, classical laminate theory with first-ply failure, thin-walled box-spar
 * stiffness, composite fatigue (Goodman constant-life diagram) and simple scaling/fatigue-count helpers.
 */

import { eigGenSym } from '../math/eigen.js';
import { solveDense } from '../math/linalg.js';

const TWO_PI = 2 * Math.PI;
const D2R = Math.PI / 180;
const G0 = 9.80665;

function need(cond, name, msg) {
  if (!cond) throw new RangeError(`${name}: ${msg}`);
}

function interp(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let i = 0;
  while (xs[i + 1] < x) i++;
  return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i]);
}

// Assemble K, Kg (per unit omega^2) and M of a cantilever discretised with nElem equal Hermite elements.
function assemble(o) {
  const x = o.x;
  need(Array.isArray(x) || ArrayBuffer.isView(x), 'x', 'array of span stations required');
  need(x.length >= 2 && x[0] === 0, 'x', 'at least two stations starting at 0 (blade root)');
  need(o.EI && o.EI.length === x.length && o.massPerLength && o.massPerLength.length === x.length, 'EI/massPerLength', 'must match x');
  const L = x[x.length - 1];
  const ne = o.nElem ?? Math.max(24, x.length - 1);
  const le = L / ne;
  const nd = 2 * (ne + 1);
  const K = []; const Kg = []; const M = [];
  for (let i = 0; i < nd; i++) { K.push(new Array(nd).fill(0)); Kg.push(new Array(nd).fill(0)); M.push(new Array(nd).fill(0)); }
  const Rh = o.hubRadius ?? 0;
  // centrifugal tension per unit omega^2: N(s) = integral_s^L m(xi) (Rh + xi) dxi (+ tip mass), by fine trapezoid
  const nf = 8 * ne;
  const Ncum = new Float64Array(nf + 1);
  Ncum[nf] = (o.tipMass ?? 0) * (Rh + L);
  for (let i = nf - 1; i >= 0; i--) {
    const a = L * i / nf; const b = L * (i + 1) / nf;
    Ncum[i] = Ncum[i + 1] + 0.5 * (interp(x, o.massPerLength, a) * (Rh + a) + interp(x, o.massPerLength, b) * (Rh + b)) * (b - a);
  }
  const Nat = (s) => { const u = s / L * nf; const i = Math.min(nf - 1, Math.floor(u)); return Ncum[i] + (u - i) * (Ncum[i + 1] - Ncum[i]); };
  const gp = [-0.7745966692414834, 0, 0.7745966692414834];
  const gw = [5 / 9, 8 / 9, 5 / 9];
  for (let e = 0; e < ne; e++) {
    const xa = e * le; const xm = xa + 0.5 * le;
    const EI = interp(x, o.EI, xm);
    const m = interp(x, o.massPerLength, xm);
    const phi = o.GAs ? 12 * EI / (interp(x, o.GAs, xm) * le * le) : 0;
    const k = EI / ((1 + phi) * le ** 3);
    const ke = [[12 * k, 6 * le * k, -12 * k, 6 * le * k], [6 * le * k, (4 + phi) * le * le * k, -6 * le * k, (2 - phi) * le * le * k],
      [-12 * k, -6 * le * k, 12 * k, -6 * le * k], [6 * le * k, (2 - phi) * le * le * k, -6 * le * k, (4 + phi) * le * le * k]];
    const c = m * le / 420;
    const me = [[156 * c, 22 * le * c, 54 * c, -13 * le * c], [22 * le * c, 4 * le * le * c, 13 * le * c, -3 * le * le * c],
      [54 * c, 13 * le * c, 156 * c, -22 * le * c], [-13 * le * c, -3 * le * le * c, -22 * le * c, 4 * le * le * c]];
    const kg = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
    for (let q = 0; q < 3; q++) {
      const s = 0.5 * (1 + gp[q]);
      const dN = [(-6 * s + 6 * s * s) / le, 1 - 4 * s + 3 * s * s, (6 * s - 6 * s * s) / le, -2 * s + 3 * s * s];
      const w = gw[q] * 0.5 * le * Nat(xa + s * le);
      for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) kg[a][b] += w * dN[a] * dN[b];
    }
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) { K[2 * e + a][2 * e + b] += ke[a][b]; M[2 * e + a][2 * e + b] += me[a][b]; Kg[2 * e + a][2 * e + b] += kg[a][b]; }
  }
  if (o.tipMass) M[nd - 2][nd - 2] += o.tipMass;
  return { K, Kg, M, ne, le, nd, L };
}

/**
 * Natural frequencies of a rotating cantilever blade in one bending direction.
 * Euler-Bernoulli (default) or shear-flexible Timoshenko-type elements (pass `GAs`; rotary inertia neglected),
 * clamped at the root, with centrifugal (geometric) stiffness from the axial tension
 * N(x) = Omega^2 integral_x^L m (Rhub + xi) dxi. Flapwise (out-of-plane): (K + Kg) phi = w^2 M phi. Edgewise
 * (in-plane): (K + Kg - Omega^2 M) phi = w^2 M phi (centrifugal softening of lead-lag motion).
 * Reference: Bramwell's Helicopter Dynamics, ch. 7; Hansen, Aerodynamics of Wind Turbines, ch. 11. Limits: no
 * flap-edge-torsion coupling, no structural twist, no pitch-bearing or hub flexibility.
 * @param {object} o
 * @param {number[]} o.x span stations from the blade root, x[0] = 0 [m]
 * @param {number[]} o.EI bending stiffness in the chosen direction at the stations [N m2]
 * @param {number[]} o.massPerLength [kg/m]
 * @param {number} [o.omega=0] rotor speed [rad/s]
 * @param {number} [o.hubRadius=0] radius of the blade root [m]
 * @param {'flap'|'edge'} [o.direction='flap']
 * @param {number[]} [o.GAs] shear stiffness k G A at the stations [N] (Timoshenko option)
 * @param {number} [o.tipMass=0] [kg]
 * @param {number} [o.nModes=3] @param {number} [o.nElem] default max(24, stations - 1)
 * @returns {{frequencies:number[], omegaN:number[], shapes:Float64Array[], x:Float64Array}} natural frequencies [Hz]
 *   and [rad/s]; mode shapes (deflection at element nodes, unit tip deflection); node positions [m]
 */
export function bladeFrequencies(o) {
  const A = assemble(o);
  const om = o.omega ?? 0;
  const dir = o.direction ?? 'flap';
  need(dir === 'flap' || dir === 'edge', 'direction', "must be 'flap' or 'edge'");
  const nr = A.nd - 2;
  const Kt = []; const Mt = [];
  for (let i = 0; i < nr; i++) {
    const kr = new Array(nr); const mr = new Array(nr);
    for (let j = 0; j < nr; j++) {
      mr[j] = A.M[i + 2][j + 2];
      kr[j] = A.K[i + 2][j + 2] + om * om * A.Kg[i + 2][j + 2] - (dir === 'edge' ? om * om * mr[j] : 0);
    }
    Kt.push(kr); Mt.push(mr);
  }
  const nModes = Math.min(o.nModes ?? 3, nr);
  const eig = eigGenSym(Kt, Mt, { nModes });
  const xs = new Float64Array(A.ne + 1);
  for (let i = 0; i <= A.ne; i++) xs[i] = i * A.le;
  const shapes = [];
  const omegaN = [];
  for (let q = 0; q < nModes; q++) {
    omegaN.push(Math.sqrt(Math.max(0, eig.values[q])));
    const v = eig.vectors[q];
    const tip = v[nr - 2] || 1;
    const sh = new Float64Array(A.ne + 1);
    for (let i = 1; i <= A.ne; i++) sh[i] = v[2 * (i - 1)] / tip;
    shapes.push(sh);
  }
  return { frequencies: omegaN.map((w) => w / TWO_PI), omegaN, shapes, x: xs };
}

/**
 * Campbell-diagram data: flapwise and edgewise frequencies versus rotor speed with the excitation lines n P, and
 * the Southwell coefficients K fitted between standstill and the highest speed, w^2 = w0^2 + K Omega^2.
 * @param {object} o
 * @param {number[]} o.x @param {number[]} o.EIflap @param {number[]} o.EIedge @param {number[]} o.massPerLength
 * @param {number[]} o.omegas rotor speeds [rad/s], ascending
 * @param {number} [o.hubRadius=0] @param {number} [o.nModes=2] @param {number[]} [o.harmonics=[1, 3, 6]]
 * @returns {{omega:number[], flap:number[][], edge:number[][], excitation:Array<{n:number, f:number[]}>,
 *   southwellFlap:number[], southwellEdge:number[], crossings:Array<{mode:string, harmonic:number, omega:number, f:number}>}}
 *   frequencies [Hz] indexed [mode][speed]; excitation lines [Hz]; Southwell coefficients [-]; resonance crossings
 *   (linear interpolation between speeds)
 */
export function campbellData(o) {
  need(Array.isArray(o.omegas) && o.omegas.length >= 2, 'omegas', 'at least two rotor speeds required');
  const nModes = o.nModes ?? 2;
  const flap = Array.from({ length: nModes }, () => []);
  const edge = Array.from({ length: nModes }, () => []);
  for (const om of o.omegas) {
    const f = bladeFrequencies({ x: o.x, EI: o.EIflap, massPerLength: o.massPerLength, omega: om, hubRadius: o.hubRadius, nModes, direction: 'flap' });
    const e = bladeFrequencies({ x: o.x, EI: o.EIedge, massPerLength: o.massPerLength, omega: om, hubRadius: o.hubRadius, nModes, direction: 'edge' });
    for (let q = 0; q < nModes; q++) { flap[q].push(f.frequencies[q]); edge[q].push(e.frequencies[q]); }
  }
  const harmonics = o.harmonics ?? [1, 3, 6];
  const excitation = harmonics.map((n) => ({ n, f: o.omegas.map((om) => n * om / TWO_PI) }));
  const last = o.omegas.length - 1;
  const sw = (arr) => arr.map((f) => ((TWO_PI * f[last]) ** 2 - (TWO_PI * f[0]) ** 2) / (o.omegas[last] ** 2 - o.omegas[0] ** 2));
  const crossings = [];
  const scan = (arr, label) => arr.forEach((f, q) => excitation.forEach((ex) => {
    for (let i = 0; i < last; i++) {
      const d0 = f[i] - ex.f[i]; const d1 = f[i + 1] - ex.f[i + 1];
      if (d0 * d1 < 0) {
        const w = d0 / (d0 - d1);
        crossings.push({ mode: `${label} ${q + 1}`, harmonic: ex.n, omega: o.omegas[i] + w * (o.omegas[i + 1] - o.omegas[i]), f: f[i] + w * (f[i + 1] - f[i]) });
      }
    }
  }));
  scan(flap, 'flap'); scan(edge, 'edge');
  return { omega: o.omegas.slice(), flap, edge, excitation, southwellFlap: sw(flap), southwellEdge: sw(edge), crossings };
}

/**
 * Static deflection of the cantilever blade under a distributed transverse load (e.g. BEM normal force per unit
 * length), optionally with centrifugal stiffening, by the beam finite elements; bending moment and shear along the
 * span by direct integration of the load (statically determinate, independent of stiffness and stiffening relief).
 * @param {object} o
 * @param {number[]} o.x span stations from the root [m] @param {number[]} o.EI [N m2] @param {number[]} o.massPerLength [kg/m]
 * @param {number[]} o.load transverse load per unit length at the stations [N/m]
 * @param {number} [o.omega=0] rotor speed for centrifugal stiffening [rad/s] @param {number} [o.hubRadius=0]
 * @param {number} [o.nElem]
 * @returns {{x:Float64Array, deflection:Float64Array, slope:Float64Array, tipDeflection:number, moment:Float64Array,
 *   shear:Float64Array, rootMoment:number, rootShear:number}} deflection [m], slope [rad], bending moment [N m] and
 *   shear [N] at the element nodes
 */
export function bladeStaticDeflection(o) {
  need(o.load && o.load.length === o.x.length, 'load', 'must match x');
  const A = assemble(o);
  const om = o.omega ?? 0;
  const nr = A.nd - 2;
  const f = new Array(nr).fill(0);
  const le = A.le;
  for (let e = 0; e < A.ne; e++) {
    const q1 = interp(o.x, o.load, e * le); const q2 = interp(o.x, o.load, (e + 1) * le);
    // consistent nodal loads of a linearly varying load
    const fe = [le * (7 * q1 + 3 * q2) / 20, le * le * (3 * q1 + 2 * q2) / 60, le * (3 * q1 + 7 * q2) / 20, -le * le * (2 * q1 + 3 * q2) / 60];
    for (let a = 0; a < 4; a++) { const d = 2 * e + a - 2; if (d >= 0) f[d] += fe[a]; }
  }
  const Kt = [];
  for (let i = 0; i < nr; i++) { const r = new Array(nr); for (let j = 0; j < nr; j++) r[j] = A.K[i + 2][j + 2] + om * om * A.Kg[i + 2][j + 2]; Kt.push(r); }
  const u = solveDense(Kt, f);
  const n1 = A.ne + 1;
  const xs = new Float64Array(n1); const w = new Float64Array(n1); const th = new Float64Array(n1);
  const Mo = new Float64Array(n1); const V = new Float64Array(n1);
  for (let i = 0; i < n1; i++) { xs[i] = i * le; if (i > 0) { w[i] = u[2 * (i - 1)]; th[i] = u[2 * (i - 1) + 1]; } }
  for (let i = n1 - 2; i >= 0; i--) {
    const q1 = interp(o.x, o.load, xs[i]); const q2 = interp(o.x, o.load, xs[i + 1]);
    V[i] = V[i + 1] + 0.5 * (q1 + q2) * le;
    // moment of the trapezoidal load on the element about its inboard node, plus the outboard shear and moment
    Mo[i] = Mo[i + 1] + V[i + 1] * le + le * le * (q1 + 2 * q2) / 6;
  }
  return { x: xs, deflection: w, slope: th, tipDeflection: w[n1 - 1], moment: Mo, shear: V, rootMoment: Mo[0], rootShear: V[0] };
}

/**
 * Blade-tip to tower clearance for an upwind rotor with the blade pointing down.
 * Undeflected horizontal distance from the tower axis to the tip = overhang cos(tilt) + L sin(tilt + precone) +
 * prebend (L = tip radius, small hub-centre offset neglected); clearance = that distance - tower radius at tip
 * height - downwind tip deflection. Design practice requires a residual clearance under factored loads (e.g. IEC
 * 61400-1 7.6.5: no contact with the deflection multiplied by the combined partial safety factors).
 * @param {object} o
 * @param {number} o.tipDeflection downwind tip deflection [m]
 * @param {number} o.overhang [m] @param {number} o.tiltDeg [deg] @param {number} o.preconeDeg [deg] @param {number} o.tipRadius [m]
 * @param {number} o.towerRadius tower outer radius at the height of the passing tip [m]
 * @param {number} [o.prebend=0] upwind tip prebend [m]
 * @returns {{undeflected:number, clearance:number, usedFraction:number, ok:boolean}} distances [m]; fraction of the
 *   undeflected clearance consumed [-]
 */
export function tipClearance(o) {
  need(o && o.tipRadius > 0 && o.towerRadius >= 0, 'tipRadius/towerRadius', 'must be valid');
  const und = o.overhang * Math.cos(o.tiltDeg * D2R) + o.tipRadius * Math.sin((o.tiltDeg + o.preconeDeg) * D2R) + (o.prebend ?? 0) - o.towerRadius;
  const cl = und - o.tipDeflection;
  return { undeflected: und, clearance: cl, usedFraction: o.tipDeflection / und, ok: cl > 0 };
}

/**
 * Bending strain at distance c from the neutral axis: eps = M c / EI (Euler-Bernoulli).
 * @param {number} M bending moment [N m] @param {number} EI bending stiffness [N m2] @param {number} c distance to the fibre [m]
 * @returns {number} strain [-]
 */
export function bendingStrain(M, EI, c) {
  need(EI > 0, 'EI', 'must be > 0');
  return M * c / EI;
}

/**
 * GENERIC span-wise mass and stiffness distributions for a blade of given length, mass and first natural
 * frequencies (not the structural data of any real blade): mass per length falls linearly to `tipMassRatio` of the
 * root value; stiffness falls log-linearly to `tipStiffnessRatio` of the root value; both are scaled so that the
 * integrated mass and the non-rotating first flapwise and edgewise frequencies match the targets.
 * @param {object} o
 * @param {number} o.length blade length [m] @param {number} o.mass blade mass [kg]
 * @param {number} o.flapFrequency first flapwise frequency at standstill [Hz] @param {number} o.edgeFrequency [Hz]
 * @param {number} [o.tipMassRatio=0.1] @param {number} [o.tipStiffnessRatio=0.002] @param {number} [o.n=25] stations
 * @returns {{x:number[], massPerLength:number[], EIflap:number[], EIedge:number[], label:string}}
 */
export function genericBladeDistributions(o) {
  need(o && o.length > 0 && o.mass > 0, 'length/mass', 'must be > 0');
  need(o.flapFrequency > 0 && o.edgeFrequency > 0, 'flapFrequency/edgeFrequency', 'must be > 0');
  const n = o.n ?? 25;
  const tm = o.tipMassRatio ?? 0.1; const tk = o.tipStiffnessRatio ?? 0.002;
  const x = []; const m = []; const sh = [];
  const m0 = o.mass / (o.length * 0.5 * (1 + tm));
  for (let i = 0; i < n; i++) {
    const s = i / (n - 1);
    x.push(s * o.length); m.push(m0 * (1 - (1 - tm) * s)); sh.push(Math.pow(tk, s));
  }
  const f1 = bladeFrequencies({ x, EI: sh, massPerLength: m, nModes: 1 }).frequencies[0];
  const kf = (o.flapFrequency / f1) ** 2; const ke = (o.edgeFrequency / f1) ** 2;
  return { x, massPerLength: m, EIflap: sh.map((v) => v * kf), EIedge: sh.map((v) => v * ke),
    label: 'generic blade mass/stiffness distributions scaled to target mass and first frequencies — not measured blade data' };
}

/**
 * GENERIC unidirectional ply properties of typical order of magnitude (editable defaults, not qualified material
 * data — use coupon-test values for design). Units: moduli and strengths [Pa], density [kg/m3], thickness [m].
 * Strengths: Xt/Xc fibre direction tension/compression, Yt/Yc transverse, S in-plane shear (all positive).
 */
export const GENERIC_PLIES = {
  glassEpoxyUD: { label: 'generic E-glass/epoxy UD ply', E1: 40e9, E2: 10e9, G12: 4e9, nu12: 0.28, Xt: 1000e6, Xc: 600e6, Yt: 40e6, Yc: 140e6, S: 60e6, density: 1900, thickness: 0.9e-3 },
  carbonEpoxyUD: { label: 'generic carbon/epoxy UD ply', E1: 135e9, E2: 9e9, G12: 5e9, nu12: 0.30, Xt: 1800e6, Xc: 1200e6, Yt: 50e6, Yc: 200e6, S: 80e6, density: 1550, thickness: 0.6e-3 },
};

/**
 * Generic biaxial (+-45) and triaxial (0/+45/-45) lay-ups assembled from a UD ply, as lists for laminateABD.
 * @param {'biax'|'triax'|'ud'} kind @param {object} [ply=GENERIC_PLIES.glassEpoxyUD] @param {number} [repeats=1]
 * @param {boolean} [symmetric=true] mirror the stack about the mid-plane
 * @returns {Array<{ply:object, angleDeg:number, thickness:number}>} lay-up, bottom to top
 */
export function genericLayup(kind, ply = GENERIC_PLIES.glassEpoxyUD, repeats = 1, symmetric = true) {
  const seq = kind === 'biax' ? [45, -45] : kind === 'triax' ? [0, 45, -45] : kind === 'ud' ? [0] : null;
  need(seq !== null, 'kind', "must be 'biax', 'triax' or 'ud'");
  let stack = [];
  for (let i = 0; i < repeats; i++) stack = stack.concat(seq);
  if (symmetric) stack = stack.concat(stack.slice().reverse());
  return stack.map((a) => ({ ply, angleDeg: a, thickness: ply.thickness }));
}

function qbar(ply, angleDeg) {
  const nu21 = ply.nu12 * ply.E2 / ply.E1;
  const d = 1 - ply.nu12 * nu21;
  const Q11 = ply.E1 / d; const Q22 = ply.E2 / d; const Q12 = ply.nu12 * ply.E2 / d; const Q66 = ply.G12;
  const c = Math.cos(angleDeg * D2R); const s = Math.sin(angleDeg * D2R);
  const c2 = c * c; const s2 = s * s; const cs = c * s;
  return {
    Q: [[Q11, Q12, 0], [Q12, Q22, 0], [0, 0, Q66]],
    Qb: [
      [Q11 * c2 * c2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * s2 * s2, (Q11 + Q22 - 4 * Q66) * s2 * c2 + Q12 * (s2 * s2 + c2 * c2), (Q11 - Q12 - 2 * Q66) * cs * c2 + (Q12 - Q22 + 2 * Q66) * cs * s2],
      [(Q11 + Q22 - 4 * Q66) * s2 * c2 + Q12 * (s2 * s2 + c2 * c2), Q11 * s2 * s2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * c2 * c2, (Q11 - Q12 - 2 * Q66) * cs * s2 + (Q12 - Q22 + 2 * Q66) * cs * c2],
      [(Q11 - Q12 - 2 * Q66) * cs * c2 + (Q12 - Q22 + 2 * Q66) * cs * s2, (Q11 - Q12 - 2 * Q66) * cs * s2 + (Q12 - Q22 + 2 * Q66) * cs * c2, (Q11 + Q22 - 2 * Q12 - 2 * Q66) * s2 * c2 + Q66 * (s2 * s2 + c2 * c2)],
    ],
    c, s,
  };
}

function inv(m) {
  const n = m.length;
  const I = [];
  for (let i = 0; i < n; i++) { const r = new Array(n).fill(0); r[i] = 1; I.push(r); }
  const X = solveDense(m.map((r) => r.slice()), I);
  return X;
}

/**
 * Classical laminate theory: extensional A [N/m], coupling B [N] and bending D [N m] stiffness matrices
 *   A = sum Qbar_k (z_k - z_{k-1}),  B = (1/2) sum Qbar_k (z_k^2 - z_{k-1}^2),  D = (1/3) sum Qbar_k (z_k^3 - z_{k-1}^3)
 * with z measured from the mid-plane, plies listed from the bottom (z = -h/2) upward, and the effective in-plane
 * engineering constants from the compliance a = (A - B D^-1 B)^-1: Ex = 1/(h a11), Ey = 1/(h a22),
 * Gxy = 1/(h a66), nu_xy = -a12/a11. Reference: Jones, R.M., Mechanics of Composite Materials, 2nd ed., ch. 4.
 * Plane stress, perfectly bonded thin plies, linear elastic.
 * @param {Array<{ply:object, angleDeg:number, thickness?:number}>} layup plies bottom to top; ply = { E1, E2, G12,
 *   nu12 [Pa, Pa, Pa, -], thickness [m] }; angle of the fibre from the laminate x axis [deg]
 * @returns {{A:number[][], B:number[][], D:number[][], thickness:number, Ex:number, Ey:number, Gxy:number, nuxy:number,
 *   arealMass:number, z:number[]}} stiffness matrices, total thickness [m], engineering constants [Pa], areal mass
 *   [kg/m2], interface coordinates [m]
 */
export function laminateABD(layup) {
  need(Array.isArray(layup) && layup.length >= 1, 'layup', 'at least one ply required');
  let h = 0; let am = 0;
  for (const l of layup) {
    const t = l.thickness ?? l.ply.thickness;
    need(t > 0, 'layup', 'ply thickness must be > 0');
    need(l.ply.E1 > 0 && l.ply.E2 > 0 && l.ply.G12 > 0, 'layup', 'ply moduli must be > 0');
    h += t; am += t * (l.ply.density ?? 0);
  }
  const z = [-0.5 * h];
  const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; const B = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; const D = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const l of layup) {
    const t = l.thickness ?? l.ply.thickness;
    const z0 = z[z.length - 1]; const z1 = z0 + t;
    z.push(z1);
    const { Qb } = qbar(l.ply, l.angleDeg);
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
      A[i][j] += Qb[i][j] * (z1 - z0); B[i][j] += Qb[i][j] * (z1 * z1 - z0 * z0) / 2; D[i][j] += Qb[i][j] * (z1 ** 3 - z0 ** 3) / 3;
    }
  }
  const Di = inv(D);
  const red = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    let s = A[i][j];
    for (let k = 0; k < 3; k++) for (let l = 0; l < 3; l++) s -= B[i][k] * Di[k][l] * B[l][j];
    red[i][j] = s;
  }
  const a = inv(red);
  return { A, B, D, thickness: h, Ex: 1 / (h * a[0][0]), Ey: 1 / (h * a[1][1]), Gxy: 1 / (h * a[2][2]), nuxy: -a[0][1] / a[0][0], arealMass: am, z };
}

/**
 * First-ply failure of a laminate under in-plane force and moment resultants. Mid-plane strains and curvatures from
 * the inverse ABD matrix; ply stresses in material axes at the bottom and top of each ply; strength ratio R (load
 * multiplier to failure) per ply by
 *  - 'max-stress': R = min over components of strength / |stress| (sign-dependent strengths);
 *  - 'tsai-wu': F1 s1 + F2 s2 + F11 s1^2 + F22 s2^2 + F66 t12^2 + 2 F12 s1 s2 = 1 with F1 = 1/Xt - 1/Xc,
 *    F11 = 1/(Xt Xc), F2 = 1/Yt - 1/Yc, F22 = 1/(Yt Yc), F66 = 1/S^2, F12 = -0.5 sqrt(F11 F22); R is the positive
 *    root of (quadratic) R^2 + (linear) R - 1 = 0.
 * Reference: Jones (1999), ch. 2.9 and 4; Tsai and Wu (1971). First-ply failure is conservative for laminates that
 * can redistribute after matrix cracking.
 * @param {Array<object>} layup as in laminateABD; plies need Xt, Xc, Yt, Yc, S [Pa]
 * @param {object} loads { N: [Nx, Ny, Nxy] [N/m], M: [Mx, My, Mxy] [N m/m] } (M optional)
 * @param {object} [o] @param {'tsai-wu'|'max-stress'} [o.criterion='tsai-wu']
 * @returns {{strengthRatio:number, failureIndex:number, criticalPly:number, criticalAngleDeg:number, mode:string,
 *   strain0:number[], curvature:number[], plies:Array<{angleDeg:number, z:number, sigma1:number, sigma2:number, tau12:number, strengthRatio:number}>}}
 *   R [-] (failure when < 1), 1/R, index of the critical ply (0 = bottom), governing mode (max-stress), mid-plane
 *   strains [-] and curvatures [1/m], ply stresses [Pa]
 */
export function firstPlyFailure(layup, loads, o = {}) {
  const crit = o.criterion ?? 'tsai-wu';
  need(crit === 'tsai-wu' || crit === 'max-stress', 'criterion', "must be 'tsai-wu' or 'max-stress'");
  const lam = laminateABD(layup);
  const N = loads.N ?? [0, 0, 0]; const Mo = loads.M ?? [0, 0, 0];
  const big = [];
  for (let i = 0; i < 3; i++) big.push([...lam.A[i], ...lam.B[i]]);
  for (let i = 0; i < 3; i++) big.push([...lam.B[i], ...lam.D[i]]);
  const sol = solveDense(big, [...N, ...Mo]);
  const e0 = sol.slice(0, 3); const kap = sol.slice(3, 6);
  const plies = [];
  let best = { R: Infinity, ply: -1, mode: '' };
  layup.forEach((l, k) => {
    const { Q, c, s } = qbar(l.ply, l.angleDeg);
    for (const zz of [lam.z[k], lam.z[k + 1]]) {
      const ex = e0[0] + zz * kap[0]; const ey = e0[1] + zz * kap[1]; const gxy = e0[2] + zz * kap[2];
      const e1 = c * c * ex + s * s * ey + c * s * gxy;
      const e2 = s * s * ex + c * c * ey - c * s * gxy;
      const g12 = -2 * c * s * ex + 2 * c * s * ey + (c * c - s * s) * gxy;
      const s1 = Q[0][0] * e1 + Q[0][1] * e2; const s2 = Q[0][1] * e1 + Q[1][1] * e2; const t12 = Q[2][2] * g12;
      const p = l.ply;
      need(p.Xt > 0 && p.Xc > 0 && p.Yt > 0 && p.Yc > 0 && p.S > 0, 'layup', 'ply strengths Xt, Xc, Yt, Yc, S required');
      let R; let mode;
      if (crit === 'max-stress') {
        const cand = [[s1 > 0 ? p.Xt / s1 : (s1 < 0 ? p.Xc / -s1 : Infinity), s1 >= 0 ? 'fibre tension' : 'fibre compression'],
          [s2 > 0 ? p.Yt / s2 : (s2 < 0 ? p.Yc / -s2 : Infinity), s2 >= 0 ? 'matrix tension' : 'matrix compression'],
          [t12 !== 0 ? p.S / Math.abs(t12) : Infinity, 'in-plane shear']];
        cand.sort((u, v) => u[0] - v[0]);
        R = cand[0][0]; mode = cand[0][1];
      } else {
        const F1 = 1 / p.Xt - 1 / p.Xc; const F2 = 1 / p.Yt - 1 / p.Yc; const F11 = 1 / (p.Xt * p.Xc); const F22 = 1 / (p.Yt * p.Yc); const F66 = 1 / (p.S * p.S);
        const F12 = -0.5 * Math.sqrt(F11 * F22);
        const a = F11 * s1 * s1 + F22 * s2 * s2 + F66 * t12 * t12 + 2 * F12 * s1 * s2;
        const b = F1 * s1 + F2 * s2;
        R = a > 0 ? (-b + Math.sqrt(b * b + 4 * a)) / (2 * a) : (b > 0 ? 1 / b : Infinity);
        mode = 'tsai-wu';
      }
      plies.push({ angleDeg: l.angleDeg, z: zz, sigma1: s1, sigma2: s2, tau12: t12, strengthRatio: R });
      if (R < best.R) best = { R, ply: k, mode };
    }
  });
  return { strengthRatio: best.R, failureIndex: 1 / best.R, criticalPly: best.ply, criticalAngleDeg: best.ply >= 0 ? layup[best.ply].angleDeg : NaN,
    mode: best.mode, strain0: e0, curvature: kap, plies };
}

/**
 * Thin-walled rectangular box-spar section stiffness from spar-cap and shear-web laminates (walls on the box
 * mid-lines; caps at +-h/2 carry flapwise bending, webs at +-b/2 carry shear):
 *   EI_flap = 2 Ec tc b (h/2)^2 + 2 Ew tw h^3/12,   EI_edge = 2 Ew tw h (b/2)^2 + 2 Ec tc b^3/12,
 *   GJ = 4 (b h)^2 / ( 2 b/(Gc tc) + 2 h/(Gw tw) )   (Bredt-Batho single cell),   EA = 2 Ec tc b + 2 Ew tw h.
 * Own bending stiffness of the cap walls about their mid-plane is neglected (thin-wall assumption, tc << h).
 * Reference: Megson, Aircraft Structures, ch. 16-18. The aerodynamic shells are not included.
 * @param {object} o
 * @param {number} o.width b, distance between web mid-lines [m] @param {number} o.height h, distance between cap mid-lines [m]
 * @param {{E:number, G:number, thickness:number, density?:number}} o.cap effective axial modulus, shear modulus [Pa], thickness [m]
 * @param {{E:number, G:number, thickness:number, density?:number}} o.web
 * @returns {{EIflap:number, EIedge:number, GJ:number, EA:number, massPerLength:number, cFlap:number}} [N m2], [N m2],
 *   [N m2], [N], [kg/m]; extreme-fibre distance for flapwise strain [m]
 */
export function boxSparSection(o) {
  need(o && o.width > 0 && o.height > 0, 'width/height', 'must be > 0');
  const c = o.cap; const w = o.web;
  need(c && c.E > 0 && c.G > 0 && c.thickness > 0 && w && w.E > 0 && w.G > 0 && w.thickness > 0, 'cap/web', 'E, G, thickness must be > 0');
  const b = o.width; const h = o.height;
  return {
    EIflap: 2 * c.E * c.thickness * b * (h / 2) ** 2 + 2 * w.E * w.thickness * h ** 3 / 12,
    EIedge: 2 * w.E * w.thickness * h * (b / 2) ** 2 + 2 * c.E * c.thickness * b ** 3 / 12,
    GJ: 4 * (b * h) ** 2 / (2 * b / (c.G * c.thickness) + 2 * h / (w.G * w.thickness)),
    EA: 2 * c.E * c.thickness * b + 2 * w.E * w.thickness * h,
    massPerLength: 2 * (c.density ?? 0) * c.thickness * b + 2 * (w.density ?? 0) * w.thickness * h,
    cFlap: 0.5 * (h + c.thickness),
  };
}

/**
 * Stress ratio R = sigma_min / sigma_max of a cycle.
 * @param {number} amplitude stress amplitude [Pa] (>= 0) @param {number} mean mean stress [Pa]
 * @returns {number} R [-] (R = -1 fully reversed, 0.1 tension-tension, 10 compression-compression)
 */
export function stressRatio(amplitude, mean) {
  const mx = mean + amplitude;
  return mx === 0 ? -Infinity : (mean - amplitude) / mx;
}

/**
 * Cycles to failure of a composite from the linear (shifted) Goodman constant-life diagram: the static envelope
 * bounds the amplitude at mean stress sm by min(Rt - sm, |Rc| + sm) = (Rt + |Rc| - |2 sm - Rt + |Rc||)/2 (tensile
 * failure on one side, compressive on the other), and constant-life lines are obtained by scaling the envelope:
 *   N = [ (Rt + |Rc| - |2 sm - Rt + |Rc||) / (2 sa) ]^m.
 * This is the form used in the Germanischer Lloyd wind-turbine guideline (2010, 5.5.4) without its material
 * partial factors; m about 10 for glass/epoxy and about 14 for carbon/epoxy are commonly quoted and must be
 * confirmed by test. Covers any R-ratio through (sa, sm).
 * @param {object} o
 * @param {number} o.amplitude stress (or strain) amplitude [Pa or -] (> 0) @param {number} o.mean mean value (tension positive)
 * @param {number} o.Rt static tensile strength [same unit] @param {number} o.Rc static compressive strength (magnitude)
 * @param {number} [o.m=10] S-N slope [-]
 * @returns {{N:number, allowableAmplitude:number, R:number}} cycles to failure (0 if the cycle exceeds the static
 *   envelope), envelope amplitude at this mean, stress ratio
 */
export function goodmanCycles(o) {
  need(o && o.amplitude > 0, 'amplitude', 'must be > 0');
  need(o.Rt > 0 && o.Rc > 0, 'Rt/Rc', 'must be > 0 (magnitudes)');
  const m = o.m ?? 10;
  const env = 0.5 * (o.Rt + o.Rc - Math.abs(2 * o.mean - o.Rt + o.Rc));
  const R = stressRatio(o.amplitude, o.mean);
  if (env <= 0 || o.amplitude >= env) return { N: env > 0 && o.amplitude === env ? 1 : 0, allowableAmplitude: Math.max(0, env), R };
  return { N: Math.pow(env / o.amplitude, m), allowableAmplitude: env, R };
}

/**
 * Miner damage of a load history on the Goodman diagram, using rainflow cycles from a caller-supplied counter.
 * @param {object} o
 * @param {{ranges:number[], means:number[], counts:number[]}} o.cycles rainflow result (ranges and means in stress
 *   or strain units; e.g. coupled.js rainflowCount applied to a strain series)
 * @param {number} o.Rt @param {number} o.Rc @param {number} [o.m=10]
 * @param {number} [o.scale=1] multiplier on the counts (lifetime / record length)
 * @param {number} [o.gamma=1] partial factor applied to amplitudes and means [-]
 * @returns {{damage:number, life:number}} Miner sum [-] and life in units of the scaled record (1/damage)
 */
export function goodmanDamage(o) {
  need(o && o.cycles && o.cycles.ranges, 'cycles', 'rainflow result required');
  const g = o.gamma ?? 1;
  let d = 0;
  for (let i = 0; i < o.cycles.ranges.length; i++) {
    const sa = 0.5 * o.cycles.ranges[i] * g;
    if (sa <= 0) continue;
    const r = goodmanCycles({ amplitude: sa, mean: o.cycles.means[i] * g, Rt: o.Rt, Rc: o.Rc, m: o.m });
    d += r.N > 0 ? o.cycles.counts[i] * (o.scale ?? 1) / r.N : Infinity;
  }
  return { damage: d, life: d > 0 ? 1 / d : Infinity };
}

/**
 * Blade mass scaling law m = mRef (R/RRef)^n. n = 3 is geometric similarity; published trend fits of built blades
 * give lower exponents (about 2.0 - 2.5) because of improved materials and design — the exponent is an argument.
 * @param {object} o @param {number} o.massRef [kg] @param {number} o.radiusRef [m] @param {number} o.radius [m] @param {number} [o.exponent=3]
 * @returns {number} mass [kg]
 */
export function bladeMassScaling(o) {
  need(o && o.massRef > 0 && o.radiusRef > 0 && o.radius > 0, 'massRef/radiusRef/radius', 'must be > 0');
  return o.massRef * Math.pow(o.radius / o.radiusRef, o.exponent ?? 3);
}

/**
 * Tip speed and comparison with a caller-supplied limit. Leading-edge rain erosion grows steeply with impact speed
 * (incubation life falls roughly with a high power of tip speed), which is why offshore rotors are commonly limited
 * to about 80 - 100 m/s unless leading-edge protection is qualified; the limit is a project decision, not a constant
 * of this library.
 * @param {object} o @param {number} o.omega rotor speed [rad/s] @param {number} o.radius tip radius [m] @param {number} [o.limit] allowable tip speed [m/s]
 * @returns {{tipSpeed:number, ok:(boolean|null)}} [m/s]; null when no limit is given
 */
export function tipSpeedCheck(o) {
  const v = o.omega * o.radius;
  return { tipSpeed: v, ok: o.limit === undefined ? null : v <= o.limit };
}

/**
 * Gravity-driven edgewise fatigue loading of a blade: each revolution reverses the root edgewise moment of the blade
 * weight, range = 2 m g r_cg (r_cg measured from the root section considered), and the number of cycles equals the
 * number of rotor revolutions.
 * @param {object} o
 * @param {number} o.bladeMass [kg] @param {number} o.cgRadius distance of the blade centre of mass from the section [m]
 * @param {number} o.meanRotorSpeed mean operating rotor speed [rad/s] @param {number} [o.years=25] @param {number} [o.availability=0.95]
 * @returns {{momentRange:number, cycles:number}} [N m], [-]
 */
export function gravityEdgewiseCycles(o) {
  need(o && o.bladeMass > 0 && o.cgRadius > 0 && o.meanRotorSpeed >= 0, 'bladeMass/cgRadius/meanRotorSpeed', 'must be valid');
  return { momentRange: 2 * o.bladeMass * G0 * o.cgRadius, cycles: o.meanRotorSpeed / TWO_PI * (o.years ?? 25) * 365.25 * 86400 * (o.availability ?? 0.95) };
}

/**
 * Stresses in a circular tubular section (tower, monopile) under axial force, bending moment, shear and torsion
 * (thin-to-moderately-thick tube, elastic): sigma = N/A +- M/W, tau = 2 V/A (peak, at the neutral axis) and
 * T D/(2 J) (uniform, at the outer fibre); von Mises at the extreme bending fibre (bending + axial + torsion) and
 * at the neutral axis (axial + shear + torsion).
 * @param {object} o
 * @param {number} o.D outer diameter [m] @param {number} o.t wall thickness [m]
 * @param {number} [o.N=0] axial force, compression positive [N] @param {number} [o.M=0] bending moment [N m]
 * @param {number} [o.V=0] shear force [N] @param {number} [o.T=0] torque [N m]
 * @returns {{area:number, W:number, sigmaAxial:number, sigmaBending:number, sigmaCompMax:number, tauShear:number,
 *   tauTorsion:number, vonMises:number}} [m2], section modulus [m3], stresses [Pa] (compression positive)
 */
export function tubularSectionStress(o) {
  need(o && o.D > 0 && o.t > 0 && 2 * o.t < o.D, 'D/t', 'must satisfy 0 < 2t < D');
  const Di = o.D - 2 * o.t;
  const A = Math.PI / 4 * (o.D ** 2 - Di ** 2);
  const I = Math.PI / 64 * (o.D ** 4 - Di ** 4);
  const W = 2 * I / o.D;
  const sa = (o.N ?? 0) / A; const sb = Math.abs(o.M ?? 0) / W;
  const tv = 2 * Math.abs(o.V ?? 0) / A; const tt = Math.abs(o.T ?? 0) * o.D / (4 * I);
  const vm1 = Math.sqrt((Math.abs(sa) + sb) ** 2 + 3 * tt * tt);
  const vm2 = Math.sqrt(sa * sa + 3 * (tv + tt) ** 2);
  return { area: A, W, sigmaAxial: sa, sigmaBending: sb, sigmaCompMax: sa + sb, tauShear: tv, tauTorsion: tt, vonMises: Math.max(vm1, vm2) };
}

/**
 * Shell-buckling check of a tubular section with a caller-supplied capacity model (hook), e.g. EN 1993-1-6 or
 * DNV-RP-C202 implemented elsewhere. Without a hook the classical elastic critical meridional stress of a perfect
 * axially compressed cylinder is used, sigma_cr = E t / (R sqrt(3 (1 - nu^2))) (about 0.605 E t/R for steel;
 * Timoshenko and Gere, Theory of Elastic Stability), reduced by an imperfection knock-down factor that MUST be
 * chosen by the user (measured cylinders buckle at 0.2 - 0.6 of the classical value) and capped at yield.
 * @param {object} o
 * @param {number} o.D [m] @param {number} o.t [m] @param {number} o.sigmaComp design compressive meridional stress [Pa]
 * @param {number} [o.E=210e9] @param {number} [o.nu=0.3] @param {number} [o.fy] yield strength cap [Pa]
 * @param {number} [o.knockdown=1] imperfection factor on sigma_cr [-]
 * @param {function(object):number} [o.capacity] hook ({ D, t, E, nu, fy }) -> characteristic buckling stress [Pa]
 * @returns {{sigmaCr:number, capacity:number, utilisation:number}} classical critical stress [Pa], capacity used [Pa], demand/capacity
 */
export function shellBucklingCheck(o) {
  need(o && o.D > 0 && o.t > 0, 'D/t', 'must be > 0');
  const E = o.E ?? 210e9; const nu = o.nu ?? 0.3;
  const sigmaCr = E * o.t / (0.5 * (o.D - o.t) * Math.sqrt(3 * (1 - nu * nu)));
  let cap = o.capacity ? o.capacity({ D: o.D, t: o.t, E, nu, fy: o.fy }) : sigmaCr * (o.knockdown ?? 1);
  if (!o.capacity && o.fy) cap = Math.min(cap, o.fy);
  need(cap > 0, 'capacity', 'must return a positive stress');
  return { sigmaCr, capacity: cap, utilisation: o.sigmaComp / cap };
}
