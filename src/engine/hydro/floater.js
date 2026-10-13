// Frequency-domain rigid-body response of floaters.
//  (a) Slender-member floaters (semi-submersibles, spars, TLP hulls, floating wind platforms) assembled from
//      cylinders, rectangular pontoons and heave plates: strip added mass, long-wave Froude–Krylov + inertia
//      excitation, radiation damping from the Haskind relation, stochastically linearised Morison drag,
//      hydrostatic and external stiffness.
//  (b) Box-shaped ships and barges: heave, pitch and roll from the panel solution of the box (radiation.js), or
//      the closed-form heave/pitch of Jensen, Mansour & Olsen (2004) with a one-degree-of-freedom roll model
//      (Froude–Krylov moment, relative-motion diffraction, Haskind radiation damping).
// Axes: origin on the still-water plane, x forward, y to port, z up. Motions (surge, sway, heave, roll, pitch,
// yaw) are about the origin. Heading is the direction the waves travel towards, measured from +x towards +y
// (0° following sea, 90° waves travelling to port, 180° head sea). Complex convention exp(+i omega t), unit
// wave amplitude, phases are leads relative to the crest at the origin (see _wavehelp.js).
//
// Validity of (a): members slender compared with the wavelength (D/lambda < 0.2), so diffraction is neglected
// and the pressure gradient is taken as uniform over each cross-section; added mass is the frequency-independent
// strip value, with an end effect only where a cylinder end is marked free (without it the pitch added inertia
// of a single deep column is about 15 % above a panel solution, which moves its pitch natural period by 3 %);
// member-to-member hydrodynamic interaction is neglected. Validity of (b): see boxBargeRAO.

import { eigGenSym } from '../math/eigen.js';
import { G, dispersion, jonswap, spectrumMoments, csolve, needPos, needVec3 } from './_wavehelp.js';
import { rectangleAddedMassCoefficient, rectangleRollInertiaCoefficient, finiteCylinderAddedMassCoefficient } from './morison.js';
import { boxHydrodynamics, boxCoefficients, boxExcitation } from './radiation.js';

const DEG = Math.PI / 180;
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const mat6 = () => Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);

function addDir(strip, n, ma, kd) {
  const rn = cross(strip.r, n);
  strip.dirs.push({ n, q: [n[0], n[1], n[2], rn[0], rn[1], rn[2]], ma, kd });
}

/**
 * Assembles a slender-member floater model: strips with added mass and drag properties, rigid-body mass
 * matrix, strip added-mass matrix and hydrostatic stiffness about the origin (on the still-water plane).
 *
 * Member types (all coordinates in the floating position, m):
 *  - { type: 'cylinder', p1, p2, D, Ca = 1.0, Cd = 1.0, end1, end2, free1 = false, free2 = false, CdAxial = 0 }
 *      transverse added mass rho Ca pi D²/4 per unit length. end1/end2 give the axial added mass of a wet end:
 *      'none' (default), 'half-disc' = ½·(8/3) rho a³ (one face of a disc, the usual estimate for a column
 *      base), 'hemisphere' = (2/3) pi rho a³ (the water in a hemisphere on the end face, the reference volume
 *      of member ends in strip-theory codes; within 3 % of a panel solution for the free end of a column five
 *      diameters deep), 'disc' = (8/3) rho a³, or a number in kg. CdAxial is the drag coefficient on the end area.
 *      free1/free2 mark a wet end that the water can flow around (the keel of a spar or of a column without a
 *      pontoon or footing). The transverse strip added mass is then reduced at that end by the loss of a finite
 *      cylinder, (1 - C_A(2 l/D)) rho Ca pi D²/4 · l with l the wetted length and C_A from DNV-RP-C205 Table D-2
 *      (finiteCylinderAddedMassCoefficient; the member and its image in the still-water plane form the finite
 *      cylinder). It is applied as a negative point mass at the end, which is where the loss occurs; without
 *      it the pitch added inertia of a column five diameters deep is 15 % above a panel solution, with it 2 %.
 *  - { type: 'rect', p1, p2, b, h, Ca?, Cd = 2.0, endMass1 = 0, endMass2 = 0 }
 *      rectangular section, h measured along the transverse direction closest to vertical (for a vertical
 *      member: along x), b along the other. Added mass per DNV-RP-C205 Appendix D Table D-1 unless Ca = [Ca_b,
 *      Ca_h] is given (motion normal to the b and h faces respectively, referred to pi (b/2)² and pi (h/2)²).
 *  - { type: 'plate', center, normal = [0,0,1], D, Cd = 0, addedMass? }
 *      thin circular heave plate: added mass (8/3) rho a³ along the normal (both faces wet) unless given in kg.
 * Typical drag coefficients (inputs): circular 0.6–1.2, sharp-edged rectangular 1.5–2.5; heave-plate values
 * depend strongly on KC and should come from model tests. See DNV-RP-C205 §6.7.
 *
 * Reference: Faltinsen (1990) Sea Loads on Ships and Offshore Structures ch. 3 (long-wave approximation for
 * semi-submersibles); DNV-RP-C205 (2019) §6.2, §7 and Appendix D (added mass); Jonkman, Robertson & Hayman,
 * HydroDyn User's Guide and Theory Manual, NREL (hemisphere reference volume of member ends).
 * @param {object} o
 * @param {object[]} o.members see above
 * @param {number} [o.mass] kg (default: displaced mass rho·V)
 * @param {number[]} o.cog centre of gravity [x, y, z], m
 * @param {number[]} o.gyradii radii of gyration about the centre of gravity [rxx, ryy, rzz], m
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @param {number} [o.depth=Infinity] water depth, m
 * @param {number[][]} [o.Kext] external (mooring/tendon/riser) stiffness about the origin, 6×6, N/m, N, N·m
 * @param {number} [o.nStrips=12] strips per member
 * @returns {object} model: { rho, g, depth, strips, waterplanes, M, A, C, Kext (number[6][6], SI), mass kg,
 *   volume m³, cob [x,y,z] m, cog, waterplaneArea m², maxDimension m, extent m (largest horizontal distance of
 *   a strip from the origin), validOmegaMax rad/s (deep-water frequency at which the largest section reaches
 *   D/lambda = 0.2) }
 */
export function buildFloater({ members, mass, cog, gyradii, rho = 1025, g = G, depth = Infinity, Kext, nStrips = 12 }) {
  if (!Array.isArray(members) || members.length === 0) throw new RangeError('members must be a non-empty array');
  needVec3(cog, 'cog'); needVec3(gyradii, 'gyradii'); needPos(rho, 'rho');
  if (!(depth > 0)) throw new RangeError('depth must be > 0');
  if (!(nStrips >= 1)) throw new RangeError('nStrips must be >= 1');
  const strips = []; const waterplanes = [];
  let maxDim = 0;
  for (const mb of members) {
    if (mb.type === 'plate') {
      needVec3(mb.center, 'plate center'); needPos(mb.D, 'plate D');
      if (mb.center[2] > 0) continue;
      const nrm = mb.normal || [0, 0, 1];
      const nl = Math.hypot(nrm[0], nrm[1], nrm[2]);
      const a = mb.D / 2;
      const s = { r: mb.center.slice(), dV: 0, dirs: [] };
      addDir(s, [nrm[0] / nl, nrm[1] / nl, nrm[2] / nl], mb.addedMass === undefined ? (8 / 3) * rho * a ** 3 : mb.addedMass,
        0.5 * rho * (mb.Cd || 0) * Math.PI * a * a);
      strips.push(s);
      continue;
    }
    if (mb.type !== 'cylinder' && mb.type !== 'rect') throw new RangeError("member type must be 'cylinder', 'rect' or 'plate'");
    needVec3(mb.p1, 'member p1'); needVec3(mb.p2, 'member p2');
    const d = [mb.p2[0] - mb.p1[0], mb.p2[1] - mb.p1[1], mb.p2[2] - mb.p1[2]];
    const L = Math.hypot(d[0], d[1], d[2]);
    if (!(L > 0)) throw new RangeError('member end points must differ');
    const t = [d[0] / L, d[1] / L, d[2] / L];
    // section axes: n2 is the transverse direction closest to vertical (x for a vertical member)
    const ref = Math.abs(t[2]) > 0.999 ? [1, 0, 0] : [0, 0, 1];
    const rt = ref[0] * t[0] + ref[1] * t[1] + ref[2] * t[2];
    let n2 = [ref[0] - rt * t[0], ref[1] - rt * t[1], ref[2] - rt * t[2]];
    const l2 = Math.hypot(n2[0], n2[1], n2[2]);
    n2 = [n2[0] / l2, n2[1] / l2, n2[2] / l2];
    const n1 = cross(t, n2);
    let area; let ma1; let ma2; let w1; let w2; let cd1; let cd2; let iOwn1; let iOwn2;
    if (mb.type === 'cylinder') {
      needPos(mb.D, 'cylinder D');
      const ca = mb.Ca === undefined ? 1 : mb.Ca;
      area = (Math.PI * mb.D * mb.D) / 4;
      ma1 = ma2 = rho * ca * area;
      w1 = w2 = mb.D;
      cd1 = cd2 = mb.Cd === undefined ? 1 : mb.Cd;
      iOwn1 = iOwn2 = (Math.PI * mb.D ** 4) / 64;
      maxDim = Math.max(maxDim, mb.D);
    } else {
      needPos(mb.b, 'rect b'); needPos(mb.h, 'rect h');
      area = mb.b * mb.h;
      // the tabulated range ends at an aspect ratio of 0.1; flatter sections use the 0.1 value
      const caB = Array.isArray(mb.Ca) ? mb.Ca[0] : rectangleAddedMassCoefficient(Math.max(mb.b / mb.h, 0.1));
      const caH = Array.isArray(mb.Ca) ? mb.Ca[1] : rectangleAddedMassCoefficient(Math.max(mb.h / mb.b, 0.1));
      // motion along n2 (the h direction) presents the face of width b; along n1 the face of width h
      ma2 = rho * caB * Math.PI * (mb.b / 2) ** 2;
      ma1 = rho * caH * Math.PI * (mb.h / 2) ** 2;
      w2 = mb.b; w1 = mb.h;
      cd1 = cd2 = mb.Cd === undefined ? 2 : mb.Cd;
      iOwn1 = (mb.h * mb.b ** 3) / 12; // about the n2 axis: governs rotation that moves points along n1
      iOwn2 = (mb.b * mb.h ** 3) / 12;
      maxDim = Math.max(maxDim, mb.b, mb.h);
    }
    const z1 = mb.p1[2]; const z2 = mb.p2[2];
    let s0 = 0; let s1 = 1;
    if (z1 > 0 && z2 > 0) continue;
    if (z1 > 0) s0 = z1 / (z1 - z2); else if (z2 > 0) s1 = z1 / (z1 - z2);
    const at = (s) => [mb.p1[0] + s * d[0], mb.p1[1] + s * d[1], mb.p1[2] + s * d[2]];
    const dl = ((s1 - s0) * L) / nStrips;
    for (let i = 0; i < nStrips; i++) {
      // t and len let the wave kinematics be averaged exactly over the strip (see stripKinematics)
      const s = { r: at(s0 + ((i + 0.5) * (s1 - s0)) / nStrips), dV: area * dl, dirs: [], t, len: dl };
      addDir(s, n1, ma1 * dl, 0.5 * rho * cd1 * w1 * dl);
      addDir(s, n2, ma2 * dl, 0.5 * rho * cd2 * w2 * dl);
      strips.push(s);
    }
    const ends = [[0, mb.end1 !== undefined ? mb.end1 : mb.endMass1], [1, mb.end2 !== undefined ? mb.end2 : mb.endMass2]];
    for (const [sEnd, spec] of ends) {
      if (spec === undefined || spec === 'none' || spec === 0) continue;
      if (sEnd < s0 || sEnd > s1) continue;
      let me;
      if (typeof spec === 'number') me = spec;
      else if (mb.type === 'cylinder' && spec === 'half-disc') me = (4 / 3) * rho * (mb.D / 2) ** 3;
      else if (mb.type === 'cylinder' && spec === 'hemisphere') me = (2 / 3) * Math.PI * rho * (mb.D / 2) ** 3;
      else if (mb.type === 'cylinder' && spec === 'disc') me = (8 / 3) * rho * (mb.D / 2) ** 3;
      else throw new RangeError("member end added mass must be 'none', 'half-disc', 'hemisphere', 'disc' (cylinders) or a number in kg");
      const s = { r: at(sEnd), dV: 0, dirs: [] };
      addDir(s, t, me, 0.5 * rho * (mb.CdAxial || 0) * area);
      strips.push(s);
    }
    if (mb.type === 'cylinder') {
      for (const [sEnd, free] of [[0, mb.free1], [1, mb.free2]]) {
        if (!free || sEnd < s0 || sEnd > s1) continue;
        const lw = (s1 - s0) * L;
        const loss = (1 - finiteCylinderAddedMassCoefficient(Math.max((2 * lw) / mb.D, 1.2))) * ma1 * lw;
        const s = { r: at(sEnd), dV: 0, dirs: [] };
        addDir(s, n1, -loss, 0); addDir(s, n2, -loss, 0);
        strips.push(s);
      }
    }
    if ((z1 < 0) !== (z2 < 0) && Math.abs(t[2]) > 1e-6) {
      const p = at(z1 / (z1 - z2));
      // own second moments are taken as for a vertical member (obliquity of the cut neglected)
      const vertical = Math.abs(t[2]) > 0.999;
      waterplanes.push({
        x: p[0], y: p[1], A: area / Math.abs(t[2]),
        // Jxx = ∫y² dA: extent along y. For a vertical rect member b lies along y and h along x.
        ixx: mb.type === 'rect' && vertical ? iOwn1 : Math.min(iOwn1, iOwn2),
        iyy: mb.type === 'rect' && vertical ? iOwn2 : Math.min(iOwn1, iOwn2),
      });
    }
  }
  let V = 0; let Mx = 0; let My = 0; let Mz = 0; let extent = 0;
  const A = mat6();
  for (const s of strips) {
    V += s.dV; Mx += s.dV * s.r[0]; My += s.dV * s.r[1]; Mz += s.dV * s.r[2];
    extent = Math.max(extent, Math.hypot(s.r[0], s.r[1]) + 0.5 * (s.len || 0));
    for (const dd of s.dirs) {
      if (dd.ma === 0) continue;
      for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) A[i][j] += dd.ma * dd.q[i] * dd.q[j];
    }
  }
  if (!(V > 0)) throw new RangeError('members give no submerged volume');
  const m = mass === undefined ? rho * V : mass;
  needPos(m, 'mass');
  let Aw = 0; let Sx = 0; let Sy = 0; let Jxx = 0; let Jyy = 0; let Jxy = 0;
  for (const wp of waterplanes) {
    Aw += wp.A; Sx += wp.A * wp.x; Sy += wp.A * wp.y;
    Jxx += wp.ixx + wp.A * wp.y * wp.y; Jyy += wp.iyy + wp.A * wp.x * wp.x; Jxy += wp.A * wp.x * wp.y;
  }
  const C = mat6();
  const rg = rho * g;
  C[2][2] = rg * Aw;
  C[2][3] = C[3][2] = rg * Sy;
  C[2][4] = C[4][2] = -rg * Sx;
  C[3][3] = rg * Jxx + rg * Mz - m * g * cog[2];
  C[4][4] = rg * Jyy + rg * Mz - m * g * cog[2];
  C[3][4] = C[4][3] = -rg * Jxy;
  const M = mat6();
  const [xg, yg, zg] = cog;
  for (let i = 0; i < 3; i++) M[i][i] = m;
  const off = [[0, zg, -yg], [-zg, 0, xg], [yg, -xg, 0]];
  const r2 = xg * xg + yg * yg + zg * zg;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      M[i][3 + j] = m * off[i][j];
      M[3 + j][i] = m * off[i][j];
      M[3 + i][3 + j] = m * ((i === j ? r2 + gyradii[i] * gyradii[i] : 0) - cog[i] * cog[j]);
    }
  }
  const K = mat6();
  if (Kext) {
    if (Kext.length !== 6) throw new RangeError('Kext must be 6×6');
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) K[i][j] = Kext[i][j];
  }
  return {
    rho, g, depth, strips, waterplanes, M, A, C, Kext: K, mass: m, volume: V, cob: [Mx / V, My / V, Mz / V], cog: cog.slice(),
    waterplaneArea: Aw, maxDimension: maxDim, extent, validOmegaMax: maxDim > 0 ? Math.sqrt((g * 2 * Math.PI) / (5 * maxDim)) : Infinity,
  };
}

/**
 * Natural periods of a floater model: uncoupled, T_i = 2 pi sqrt((M+A)_ii / (C+Kext)_ii), and coupled, from
 * the generalised eigenproblem (C + Kext) phi = omega² (M + A) phi. Degrees of freedom without restoring have
 * no natural period (null / omitted).
 * Reference: undamped free vibration; Faltinsen (1990) ch. 3 (natural periods of semi-submersibles and TLPs).
 * Added mass is the strip value, taken as frequency independent.
 * @param {object} model from buildFloater
 * @returns {{uncoupled:(number|null)[], coupled:number[]|null}} periods in s (coupled sorted descending)
 */
export function naturalPeriods(model) {
  const uncoupled = [];
  const Mt = mat6(); const K = mat6();
  for (let i = 0; i < 6; i++) {
    for (let j = 0; j < 6; j++) { Mt[i][j] = model.M[i][j] + model.A[i][j]; K[i][j] = model.C[i][j] + model.Kext[i][j]; }
    uncoupled.push(K[i][i] > 0 ? 2 * Math.PI * Math.sqrt(Mt[i][i] / K[i][i]) : null);
  }
  let coupled = null;
  try {
    const Ks = mat6();
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) Ks[i][j] = 0.5 * (K[i][j] + K[j][i]);
    const { values } = eigGenSym(Ks, Mt);
    let vmax = 0;
    for (const v of values) vmax = Math.max(vmax, Math.abs(v));
    coupled = Array.from(values).filter((v) => v > 1e-9 * vmax).map((v) => (2 * Math.PI) / Math.sqrt(v)).sort((a, b) => b - a);
  } catch (err) { coupled = null; }
  return { uncoupled, coupled };
}

/**
 * JONSWAP spectrum sampled on a frequency grid (DNV-RP-C205 (2019) §3.5.5).
 * @param {ArrayLike<number>} omega rad/s
 * @param {{Hs:number, Tp:number, gamma?:number}} sea m, s, -
 * @returns {Float64Array} S(omega), m²·s/rad
 */
export function waveSpectrum(omega, sea) {
  const S = new Float64Array(omega.length);
  for (let i = 0; i < omega.length; i++) S[i] = jonswap(omega[i], sea);
  return S;
}

// sinh(z)/z for complex z = a + i b: the mean of exp(z s) over -1 <= s <= 1
function csinhc(a, b, out) {
  const m2 = a * a + b * b;
  if (m2 < 1e-6) { out.re = 1 + (a * a - b * b) / 6; out.im = (a * b) / 3; return out; }
  const sr = Math.sinh(a) * Math.cos(b); const si = Math.cosh(a) * Math.sin(b);
  out.re = (sr * a + si * b) / m2; out.im = (si * a - sr * b) / m2;
  return out;
}

const _sp = { re: 1, im: 0 }; const _sm = { re: 1, im: 0 }; const _kin = { chR: 0, chI: 0, shR: 0, shI: 0 };

// Complex depth-and-phase factors of the undisturbed wave averaged over the length of a strip:
// ch = <cosh k(z+d)/sinh kd · exp(-i k·x)>, sh = <sinh k(z+d)/sinh kd · exp(-i k·x)> (both exp(kz - i k·x) in deep
// water). Along a straight strip the integrand is a sum of exponentials, so the mean is the centre value times
// sinh(q)/q, q = ½ len (± k t_z - i k t_h). Sampling at the centre instead leaves an error of (k len)²/24, which
// for a deep vertical column is of the size of the heave force itself (that force is the small difference
// rho g A_wp (1 - k ∫ e^{kz} dz) = rho g A_wp e^{-kT}).
function stripKinematics(s, k, cb, sb, depth, out) {
  const ph = k * (s.r[0] * cb + s.r[1] * sb);
  const c = Math.cos(ph); const sn = Math.sin(ph);
  const half = 0.5 * (s.len || 0);
  let pr = 1; let pi = 0; let mr = 1; let mi = 0;
  if (half > 0) {
    const bh = -k * (s.t[0] * cb + s.t[1] * sb) * half; const az = k * s.t[2] * half;
    csinhc(az, bh, _sp); pr = _sp.re; pi = _sp.im;
    if (Number.isFinite(depth)) { csinhc(-az, bh, _sm); mr = _sm.re; mi = _sm.im; }
  }
  const e1 = Math.exp(k * s.r[2]);
  let chR; let chI; let shR; let shI;
  if (!Number.isFinite(depth)) { chR = shR = e1 * pr; chI = shI = e1 * pi; } else {
    const e2 = Math.exp(-k * (s.r[2] + 2 * depth)); const d = 1 - Math.exp(-2 * k * depth);
    chR = (e1 * pr + e2 * mr) / d; chI = (e1 * pi + e2 * mi) / d;
    shR = (e1 * pr - e2 * mr) / d; shI = (e1 * pi - e2 * mi) / d;
  }
  // times exp(-i ph)
  out.chR = chR * c + chI * sn; out.chI = chI * c - chR * sn;
  out.shR = shR * c + shI * sn; out.shI = shI * c - shR * sn;
  return out;
}

// Potential-flow excitation of the strip model per unit wave amplitude, added to F (re in 0..5, im in 6..11).
// When u is given it receives the strip velocities (re, im; 3 per strip) for the drag linearisation.
function waveExcitation(model, w, k, cb, sb, F, u) {
  const { rho, g, depth, strips } = model;
  for (let is = 0; is < strips.length; is++) {
    const s = strips[is];
    stripKinematics(s, k, cb, sb, depth, _kin);
    // velocity u = (w ch cb, w ch sb, i w sh); acceleration a = i w u
    const uhR = w * _kin.chR; const uhI = w * _kin.chI;
    const uzR = -w * _kin.shI; const uzI = w * _kin.shR;
    if (u) {
      u.re[3 * is] = uhR * cb; u.re[3 * is + 1] = uhR * sb; u.re[3 * is + 2] = uzR;
      u.im[3 * is] = uhI * cb; u.im[3 * is + 1] = uhI * sb; u.im[3 * is + 2] = uzI;
    }
    const aR = [-w * uhI * cb, -w * uhI * sb, -w * uzI];
    const aI = [w * uhR * cb, w * uhR * sb, w * uzR];
    if (s.dV > 0) {
      const m0 = rho * s.dV;
      const cR = cross(s.r, aR); const cI = cross(s.r, aI);
      for (let i = 0; i < 3; i++) {
        F[i] += m0 * aR[i]; F[6 + i] += m0 * aI[i];
        F[3 + i] += m0 * cR[i]; F[9 + i] += m0 * cI[i];
      }
    }
    for (const d of s.dirs) {
      if (d.ma === 0) continue;
      const pR = d.ma * (aR[0] * d.n[0] + aR[1] * d.n[1] + aR[2] * d.n[2]);
      const pI = d.ma * (aI[0] * d.n[0] + aI[1] * d.n[1] + aI[2] * d.n[2]);
      for (let i = 0; i < 6; i++) { F[i] += pR * d.q[i]; F[6 + i] += pI * d.q[i]; }
    }
  }
  for (const wp of model.waterplanes) {
    const ph = k * (wp.x * cb + wp.y * sb);
    const fR = rho * g * wp.A * Math.cos(ph); const fI = -rho * g * wp.A * Math.sin(ph);
    F[2] += fR; F[8] += fI;
    F[3] += wp.y * fR; F[9] += wp.y * fI;
    F[4] -= wp.x * fR; F[10] -= wp.x * fI;
  }
  return F;
}

const _radCache = new WeakMap();

/**
 * Wave-radiation damping of a slender-member floater from the Haskind relation: the damping follows from the
 * excitation of the same body in waves from all directions,
 *   B_ij(omega) = k / (8 pi rho g c_g) ∫₀^{2π} Re{ X_i(beta) X_j*(beta) } d beta,
 * with X the excitation per unit wave amplitude, k the wave number and c_g = (omega / 2k)(1 + 2kd / sinh 2kd) the
 * group velocity (deep water: B_ij = omega³ / (4 pi rho g³) ∫ …). X is the model's own long-wave excitation
 * (floaterRAO without drag), so the damping is as accurate as that excitation: valid for D/lambda < 0.2. The
 * matrix is symmetric and positive semi-definite by construction. The heading integral is a periodic
 * trapezoidal sum with enough points for the phase differences across the model (2 k extent + 16).
 * Checks: a vertical column in heave, B_33 = omega³ |X_3|² / (2 rho g³); in surge, B_11 = omega³ |X_1|² / (4 rho g³).
 * Reference: Newman (1962) The exciting forces on fixed bodies in waves, J. Ship Research 6(3), 10–17; Newman
 * (1977) Marine Hydrodynamics ch. 6 (Haskind relations); Falnes (2002) Ocean Waves and Oscillating Systems ch. 5.
 * @param {object} o
 * @param {object} o.model from buildFloater
 * @param {ArrayLike<number>} o.omega rad/s
 * @param {number} [o.nHeadings] number of wave directions (default: chosen per frequency as above)
 * @returns {number[][][]} per frequency a 6×6 matrix, N·s/m, N·s, N·m·s
 */
export function radiationDamping({ model, omega, nHeadings }) {
  if (!model || !model.strips) throw new RangeError('model must come from buildFloater');
  if (!omega || omega.length < 1) throw new RangeError('omega must be a non-empty array');
  const { rho, g, depth } = model;
  let cache = _radCache.get(model);
  if (!cache) { cache = new Map(); _radCache.set(model, cache); }
  let extent = model.extent;
  if (extent === undefined) { extent = 0; for (const s of model.strips) extent = Math.max(extent, Math.hypot(s.r[0], s.r[1]) + 0.5 * (s.len || 0)); }
  const F = new Float64Array(12);
  const out = [];
  for (let iw = 0; iw < omega.length; iw++) {
    const w = omega[iw];
    if (!(w > 0)) throw new RangeError('omega must be > 0');
    const key = nHeadings === undefined ? w : `${w}/${nHeadings}`;
    let B = cache.get(key);
    if (!B) {
      const k = dispersion(w, depth, g);
      const kd = k * depth;
      const cg = (w / (2 * k)) * (Number.isFinite(depth) && kd < 350 ? 1 + (2 * kd) / Math.sinh(2 * kd) : 1);
      const n = nHeadings === undefined ? Math.min(720, Math.ceil(2 * k * extent) + 16) : nHeadings;
      if (!(n >= 4)) throw new RangeError('nHeadings must be >= 4');
      B = mat6();
      for (let ih = 0; ih < n; ih++) {
        const beta = (2 * Math.PI * ih) / n;
        F.fill(0);
        waveExcitation(model, w, k, Math.cos(beta), Math.sin(beta), F);
        for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) B[i][j] += F[i] * F[j] + F[6 + i] * F[6 + j];
      }
      const c = (k / (8 * Math.PI * rho * g * cg)) * ((2 * Math.PI) / n);
      for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) B[i][j] *= c;
      cache.set(key, B);
    }
    out.push(B.map((r) => r.slice()));
  }
  return out;
}

/**
 * Six-degree-of-freedom response amplitude operators of a slender-member floater.
 * Excitation per unit wave amplitude: F = sum over strips of (rho dV a + m_a,n (a·n) n) with the undisturbed
 * fluid acceleration a averaged exactly over the length of each strip (phase by position), plus rho g A_wp eta
 * at each waterplane crossing (the Gauss-theorem closure of the Froude–Krylov pressure integral); moments by
 * r × F with the lever of the strip centre.
 * Radiation damping: from the Haskind relation with this excitation (radiationDamping above), a full 6×6 matrix
 * per frequency; it is the damping a potential-flow solution has. `extraDampingRatio` adds a fraction of critical
 * damping per degree of freedom, B_ii = 2 zeta_i sqrt((M+A)_ii (C+Kext)_ii), for damping that is not modelled
 * (viscous damping when drag is off, mooring lines, risers).
 * Viscous drag is linearised stochastically for the given sea state (Borgman): f = ½ rho Cd D sqrt(8/pi)
 * sigma_vr v_r per unit length, with sigma_vr the standard deviation of the relative normal velocity; it
 * contributes damping and excitation and is iterated to convergence.
 * Reference: Faltinsen (1990) ch. 3 and ch. 7 (viscous damping); Newman (1962) (Haskind relation); Borgman (1967)
 * spectral analysis of ocean wave forces on piling; DNV-RP-C205 (2019) §6.2, §7.
 * Valid for D/lambda < 0.2 (omega < model.validOmegaMax); no diffraction, no frequency-dependent added mass, no
 * second-order effects. Against a panel solution (Capytaine) of a single column D = 8 m, draught 40 m, hemisphere
 * end mass and free keel end (D/lambda = 0.2 at 5.1 s; src/spec/cross-solver.json): heave force within 5 % from
 * 5 s, surge force from 6 s (+6 % at 5 s, +47 % at 4 s); heave motion within 5 % from 5 s, the resonance included;
 * surge and pitch motion within 5 % outside 0.9–1.4 times the pitch natural period and within 25 % inside, where
 * the peak has radiation damping only and follows the natural period (1.2 % above the panel solution).
 * @param {object} o
 * @param {object} o.model from buildFloater
 * @param {ArrayLike<number>} o.omega rad/s, ascending
 * @param {number} [o.headingDeg=0] wave propagation heading, deg
 * @param {{Hs:number, Tp:number, gamma?:number} | {S:ArrayLike<number>}} [o.seaState] sea state for the drag
 *   linearisation (JONSWAP parameters or a spectrum on `omega`, m²·s/rad); without it drag is omitted
 * @param {string} [o.radiation='haskind'] 'haskind' (computed radiation damping) or 'none'
 * @param {number|number[]} [o.extraDampingRatio=0] additional damping, fraction of critical per DOF (typical
 *   0.01–0.05 when drag is not modelled; calibrate against model tests)
 * @param {number|number[]} [o.radiationDampingRatio] former name of extraDampingRatio, still accepted
 * @param {number[][]} [o.extraDamping] additional linear damping 6×6, N·s/m …
 * @param {boolean} [o.viscous=true]
 * @param {number} [o.maxIter=60]
 * @param {number} [o.tol=1e-4] relative convergence of the velocity standard deviations
 * @returns {{omega:Float64Array, period:Float64Array, headingDeg:number, depth:number, g:number, re:Float64Array[],
 *   im:Float64Array[], amp:Float64Array[], phase:Float64Array[], force:{re:Float64Array[], im:Float64Array[]},
 *   radiationDamping:number[][][]|null, viscousDamping:number[][], iterations:number, converged:boolean,
 *   validOmegaMax:number}} per DOF arrays over omega: m/m (translations) and rad/m (rotations); phase rad; force
 *   N/m and N·m/m (incl. linearised drag excitation); radiationDamping per frequency and viscousDamping N·s/m …
 */
export function floaterRAO({ model, omega, headingDeg = 0, seaState, radiation = 'haskind', extraDampingRatio, radiationDampingRatio, extraDamping, viscous = true, maxIter = 60, tol = 1e-4 }) {
  if (!model || !model.strips) throw new RangeError('model must come from buildFloater');
  if (!omega || omega.length < 1) throw new RangeError('omega must be a non-empty array');
  if (radiation !== 'haskind' && radiation !== 'none') throw new RangeError("radiation must be 'haskind' or 'none'");
  if (extraDampingRatio !== undefined && radiationDampingRatio !== undefined) throw new RangeError('give extraDampingRatio or its former name radiationDampingRatio, not both');
  const zeta = extraDampingRatio !== undefined ? extraDampingRatio : radiationDampingRatio !== undefined ? radiationDampingRatio : 0;
  const nw = omega.length;
  const { depth, g } = model;
  const cb = Math.cos(headingDeg * DEG); const sb = Math.sin(headingDeg * DEG);
  const Mt = mat6(); const K = mat6(); const B0 = mat6();
  for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
    Mt[i][j] = model.M[i][j] + model.A[i][j];
    K[i][j] = model.C[i][j] + model.Kext[i][j];
    if (extraDamping) B0[i][j] = extraDamping[i][j];
  }
  for (let i = 0; i < 6; i++) {
    const z = Array.isArray(zeta) ? zeta[i] : zeta;
    if (!(z >= 0)) throw new RangeError('extraDampingRatio must be >= 0');
    if (K[i][i] > 0) B0[i][i] += 2 * z * Math.sqrt(Mt[i][i] * K[i][i]);
  }
  const useDrag = viscous && !!seaState;
  let S = null;
  if (useDrag) {
    if (seaState.S) {
      if (seaState.S.length !== nw) throw new RangeError('seaState.S must have the length of omega');
      S = seaState.S;
    } else S = waveSpectrum(omega, seaState);
  }
  const dirs = [];
  if (useDrag) for (const s of model.strips) for (const d of s.dirs) if (d.kd > 0) dirs.push({ d, s });
  const nd = dirs.length;
  const unR = new Float64Array(nd * nw); const unI = new Float64Array(nd * nw);
  const Fr = Array.from({ length: 6 }, () => new Float64Array(nw));
  const Fi = Array.from({ length: 6 }, () => new Float64Array(nw));
  const stripIndex = new Map();
  model.strips.forEach((s, i) => stripIndex.set(s, i));
  const ns = model.strips.length;
  const u = { re: new Float64Array(3 * ns), im: new Float64Array(3 * ns) };
  const uR = u.re; const uI = u.im;
  const F = new Float64Array(12);
  for (let iw = 0; iw < nw; iw++) {
    const w = omega[iw];
    if (!(w > 0)) throw new RangeError('omega must be > 0');
    if (iw > 0 && !(w > omega[iw - 1])) throw new RangeError('omega must be ascending');
    const k = dispersion(w, depth, g);
    F.fill(0);
    waveExcitation(model, w, k, cb, sb, F, u);
    for (let i = 0; i < 6; i++) { Fr[i][iw] = F[i]; Fi[i][iw] = F[6 + i]; }
    for (let id = 0; id < nd; id++) {
      const is = stripIndex.get(dirs[id].s);
      const n = dirs[id].d.n;
      unR[id * nw + iw] = uR[3 * is] * n[0] + uR[3 * is + 1] * n[1] + uR[3 * is + 2] * n[2];
      unI[id * nw + iw] = uI[3 * is] * n[0] + uI[3 * is + 1] * n[1] + uI[3 * is + 2] * n[2];
    }
  }
  const Brad = radiation === 'haskind' ? radiationDamping({ model, omega }) : null;
  const xr = Array.from({ length: 6 }, () => new Float64Array(nw));
  const xi = Array.from({ length: 6 }, () => new Float64Array(nw));
  const FtR = Array.from({ length: 6 }, () => new Float64Array(nw));
  const FtI = Array.from({ length: 6 }, () => new Float64Array(nw));
  const sigma = new Float64Array(nd);
  const lin = Math.sqrt(8 / Math.PI);
  const integ = new Float64Array(nw);
  const trap = () => { let s = 0; for (let i = 1; i < nw; i++) s += 0.5 * (integ[i] + integ[i - 1]) * (omega[i] - omega[i - 1]); return s; };
  for (let id = 0; id < nd; id++) {
    for (let iw = 0; iw < nw; iw++) integ[iw] = S[iw] * (unR[id * nw + iw] ** 2 + unI[id * nw + iw] ** 2);
    sigma[id] = Math.sqrt(Math.max(trap(), 0));
  }
  const Bv = mat6();
  const bd = new Float64Array(nd);
  const solveAll = () => {
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) Bv[i][j] = 0;
    for (let id = 0; id < nd; id++) {
      bd[id] = dirs[id].d.kd * lin * sigma[id];
      const q = dirs[id].d.q;
      for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) Bv[i][j] += bd[id] * q[i] * q[j];
    }
    for (let iw = 0; iw < nw; iw++) {
      const w = omega[iw];
      const Ar = mat6(); const Ai = mat6();
      const br = new Array(6); const bi = new Array(6);
      for (let i = 0; i < 6; i++) {
        br[i] = Fr[i][iw]; bi[i] = Fi[i][iw];
        for (let j = 0; j < 6; j++) {
          Ar[i][j] = K[i][j] - w * w * Mt[i][j];
          Ai[i][j] = w * (B0[i][j] + Bv[i][j] + (Brad ? Brad[iw][i][j] : 0));
        }
      }
      for (let id = 0; id < nd; id++) {
        const q = dirs[id].d.q;
        const fr = bd[id] * unR[id * nw + iw]; const fi = bd[id] * unI[id * nw + iw];
        for (let i = 0; i < 6; i++) { br[i] += fr * q[i]; bi[i] += fi * q[i]; }
      }
      for (let i = 0; i < 6; i++) { FtR[i][iw] = br[i]; FtI[i][iw] = bi[i]; }
      const x = csolve(Ar, Ai, br, bi);
      for (let i = 0; i < 6; i++) { xr[i][iw] = x.re[i]; xi[i][iw] = x.im[i]; }
    }
  };
  let iterations = 0; let converged = !useDrag || nd === 0;
  solveAll();
  if (useDrag && nd > 0) {
    for (iterations = 1; iterations <= maxIter; iterations++) {
      let dmax = 0; let smax = 0;
      for (let id = 0; id < nd; id++) {
        const q = dirs[id].d.q;
        for (let iw = 0; iw < nw; iw++) {
          let vR = 0; let vI = 0;
          for (let i = 0; i < 6; i++) { vR += q[i] * xr[i][iw]; vI += q[i] * xi[i][iw]; }
          // body velocity = i w xi
          const rR = unR[id * nw + iw] + omega[iw] * vI;
          const rI = unI[id * nw + iw] - omega[iw] * vR;
          integ[iw] = S[iw] * (rR * rR + rI * rI);
        }
        const sNew = Math.sqrt(Math.max(trap(), 0));
        dmax = Math.max(dmax, Math.abs(sNew - sigma[id]));
        sigma[id] = 0.5 * (sigma[id] + sNew);
        smax = Math.max(smax, sigma[id]);
      }
      if (dmax <= tol * smax || smax === 0) { converged = true; break; }
      solveAll();
    }
  }
  const amp = []; const phase = []; const period = new Float64Array(nw);
  for (let i = 0; i < 6; i++) {
    amp.push(new Float64Array(nw)); phase.push(new Float64Array(nw));
    for (let iw = 0; iw < nw; iw++) { amp[i][iw] = Math.hypot(xr[i][iw], xi[i][iw]); phase[i][iw] = Math.atan2(xi[i][iw], xr[i][iw]); }
  }
  for (let iw = 0; iw < nw; iw++) period[iw] = (2 * Math.PI) / omega[iw];
  return {
    omega: Float64Array.from(omega), period, headingDeg, depth, g, re: xr, im: xi, amp, phase, force: { re: FtR, im: FtI },
    radiationDamping: Brad, viscousDamping: Bv.map((r) => r.slice()), iterations, converged, validOmegaMax: model.validOmegaMax,
  };
}

/**
 * Motion, velocity and acceleration transfer functions at a point rigidly attached to the floater:
 * s = xi_t + theta × r (small rotations), v = i omega s, a = -omega² s.
 * Reference: rigid-body kinematics; DNV-RP-C205 (2019) §7. Gravity components due to roll/pitch are not
 * added to the accelerations.
 * @param {{rao:{omega:ArrayLike<number>, re:ArrayLike<number>[], im:ArrayLike<number>[]}, point:number[]}} o
 *   RAO set from floaterRAO or hydrodb raoSolve; point [x, y, z] relative to the motion origin, m
 * @returns {{omega:Float64Array, re:Float64Array[], im:Float64Array[], amp:Float64Array[], velAmp:Float64Array[],
 *   accAmp:Float64Array[]}} arrays [x, y, z] over omega: displacement m/m, velocity (m/s)/m, acceleration (m/s²)/m
 */
export function pointMotion({ rao, point }) {
  needVec3(point, 'point');
  const nw = rao.omega.length;
  const [x, y, z] = point;
  const mk = () => [new Float64Array(nw), new Float64Array(nw), new Float64Array(nw)];
  const re = mk(); const im = mk(); const amp = mk(); const velAmp = mk(); const accAmp = mk();
  for (const [part, src] of [[re, rao.re], [im, rao.im]]) {
    for (let i = 0; i < nw; i++) {
      part[0][i] = src[0][i] + src[4][i] * z - src[5][i] * y;
      part[1][i] = src[1][i] + src[5][i] * x - src[3][i] * z;
      part[2][i] = src[2][i] + src[3][i] * y - src[4][i] * x;
    }
  }
  for (let c = 0; c < 3; c++) for (let i = 0; i < nw; i++) {
    amp[c][i] = Math.hypot(re[c][i], im[c][i]);
    velAmp[c][i] = rao.omega[i] * amp[c][i];
    accAmp[c][i] = rao.omega[i] ** 2 * amp[c][i];
  }
  return { omega: Float64Array.from(rao.omega), re, im, amp, velAmp, accAmp };
}

/**
 * Short-term statistics of a linear response in a wave spectrum: S_R = |RAO|² S, moments m0, m2, standard
 * deviation, significant amplitude 2 sqrt(m0), mean zero-crossing period 2 pi sqrt(m0/m2) and the most probable
 * maximum amplitude sqrt(2 m0 ln N), N = duration / Tz (narrow-band Rayleigh maxima).
 * Reference: DNV-RP-C205 (2019) §3.5 and §7; Ochi (1973) On prediction of extreme values, J. Ship Research 17.
 * Linear, Gaussian, stationary response.
 * @param {{omega:ArrayLike<number>, amp:ArrayLike<number>, S:ArrayLike<number>, duration?:number}} o rad/s, response
 *   per unit wave amplitude, wave spectrum m²·s/rad, duration s (default 10800)
 * @returns {{m0:number, m2:number, sigma:number, significant:number, Tz:number, nCycles:number, mpm:number,
 *   spectrum:Float64Array}} response units; Tz s; response spectrum unit²·s/rad
 */
export function responseStatistics({ omega, amp, S, duration = 10800 }) {
  if (omega.length !== amp.length || omega.length !== S.length) throw new RangeError('omega, amp and S must have equal lengths');
  needPos(duration, 'duration');
  const sr = new Float64Array(omega.length);
  for (let i = 0; i < sr.length; i++) sr[i] = amp[i] * amp[i] * S[i];
  const { m0, m2 } = spectrumMoments(omega, sr);
  const Tz = m2 > 0 ? 2 * Math.PI * Math.sqrt(m0 / m2) : Infinity;
  const nCycles = Number.isFinite(Tz) ? duration / Tz : 0;
  return {
    m0, m2, sigma: Math.sqrt(m0), significant: 2 * Math.sqrt(m0), Tz, nCycles,
    mpm: nCycles > 1 ? Math.sqrt(2 * m0 * Math.log(nCycles)) : Math.sqrt(m0), spectrum: sr,
  };
}

/**
 * Relative vertical motion between the (undisturbed) wave surface and a point on the floater:
 * r = alpha·eta(x, y) - z_p, with eta = exp(-i k (x cos beta + y sin beta)) per unit amplitude.
 * `waveFactor` alpha scales the incident elevation to allow for wave asymmetry / local amplification (input;
 * 1.0 = undisturbed linear wave).
 * Reference: DNV-RP-C205 (2019) §8.2 (air gap). Diffraction, run-up and radiated waves are not modelled.
 * @param {{rao:object, point:number[], waveFactor?:number}} o rao from floaterRAO (carries heading, depth, g)
 * @returns {{omega:Float64Array, re:Float64Array, im:Float64Array, amp:Float64Array}} m/m
 */
export function relativeMotion({ rao, point, waveFactor = 1 }) {
  const pm = pointMotion({ rao, point });
  const nw = rao.omega.length;
  const re = new Float64Array(nw); const im = new Float64Array(nw); const amp = new Float64Array(nw);
  const cb = Math.cos(rao.headingDeg * DEG); const sb = Math.sin(rao.headingDeg * DEG);
  for (let i = 0; i < nw; i++) {
    const k = dispersion(rao.omega[i], rao.depth === undefined ? Infinity : rao.depth, rao.g || G);
    const ph = k * (point[0] * cb + point[1] * sb);
    re[i] = waveFactor * Math.cos(ph) - pm.re[2][i];
    im[i] = -waveFactor * Math.sin(ph) - pm.im[2][i];
    amp[i] = Math.hypot(re[i], im[i]);
  }
  return { omega: pm.omega, re, im, amp };
}

/**
 * Air gap at a deck point in a sea state: minimum air gap = static air gap - most probable maximum relative
 * motion, and the Rayleigh probability per cycle that the relative motion exceeds the static air gap.
 * Reference: DNV-RP-C205 (2019) §8.2. Linear statistics of the undisturbed wave (scaled by `waveFactor`).
 * @param {{rao:object, point:number[], staticAirGap:number, S:ArrayLike<number>, duration?:number,
 *   waveFactor?:number}} o point m, still-water air gap m, wave spectrum m²·s/rad on rao.omega, duration s
 * @returns {{sigma:number, mpmRelative:number, minAirGap:number, exceedancePerCycle:number, Tz:number}} m, m, m, -, s
 */
export function airGap({ rao, point, staticAirGap, S, duration = 10800, waveFactor = 1 }) {
  if (!(staticAirGap >= 0)) throw new RangeError('staticAirGap must be >= 0');
  const rel = relativeMotion({ rao, point, waveFactor });
  const st = responseStatistics({ omega: rao.omega, amp: rel.amp, S, duration });
  return {
    sigma: st.sigma, mpmRelative: st.mpm, minAirGap: staticAirGap - st.mpm, Tz: st.Tz,
    exceedancePerCycle: st.m0 > 0 ? Math.exp((-staticAirGap * staticAirGap) / (2 * st.m0)) : 0,
  };
}

/**
 * Operability from a scatter diagram: a sea state is workable when every response returned by `evaluate` is
 * within its limit; operability is the summed probability of workable sea states.
 * Reference: DNV-ST-N001 (operational limiting criteria); NORDFORSK (1987) seakeeping criteria.
 * @param {object} o
 * @param {{Hs:number, Tp:number, prob:number}[]} o.seaStates probabilities (any normalisation)
 * @param {function(number, number): Object<string, number>} o.evaluate (Hs, Tp) -> named responses
 * @param {Object<string, number>} o.limits maximum allowed value per response name
 * @returns {{operability:number, downtime:number, perState:{Hs:number, Tp:number, prob:number, ok:boolean,
 *   utilisation:number, governing:string|null}[], governingCounts:Object<string, number>}} fractions 0..1;
 *   governingCounts = probability mass of the unworkable states per governing criterion
 */
export function operability({ seaStates, evaluate, limits }) {
  if (!Array.isArray(seaStates) || seaStates.length === 0) throw new RangeError('seaStates must be a non-empty array');
  if (typeof evaluate !== 'function') throw new RangeError('evaluate must be a function');
  let tot = 0; let okP = 0;
  const governingCounts = {};
  const perState = seaStates.map((ss) => {
    const r = evaluate(ss.Hs, ss.Tp);
    let util = 0; let gov = null;
    for (const key of Object.keys(limits)) {
      if (!(key in r)) throw new RangeError(`evaluate did not return the response '${key}'`);
      const u = Math.abs(r[key]) / limits[key];
      if (u > util) { util = u; gov = key; }
    }
    const ok = util <= 1;
    tot += ss.prob;
    if (ok) okP += ss.prob; else governingCounts[gov] = (governingCounts[gov] || 0) + ss.prob;
    return { Hs: ss.Hs, Tp: ss.Tp, prob: ss.prob, ok, utilisation: util, governing: gov };
  });
  if (!(tot > 0)) throw new RangeError('seaStates probabilities must sum to > 0');
  for (const key of Object.keys(governingCounts)) governingCounts[key] /= tot;
  return { operability: okP / tot, downtime: 1 - okP / tot, perState, governingCounts };
}

/**
 * Limiting significant wave height versus peak period for a linear response: Hs_lim(Tp) = limit / R(1 m, Tp),
 * where R is the significant amplitude or the most probable maximum of the response in a JONSWAP sea of unit Hs.
 * Reference: linear spectral scaling (response proportional to Hs); DNV-RP-C205 (2019) §7.
 * Not valid when the RAO itself depends on Hs (drag-linearised RAOs): use `operability` with re-evaluation.
 * @param {{omega:ArrayLike<number>, amp:ArrayLike<number>, limit:number, Tp:number[], gamma?:number,
 *   statistic?:string, duration?:number}} o statistic 'significant' (default) or 'mpm'
 * @returns {{Tp:number[], HsLimit:number[]}} s, m
 */
export function limitingHs({ omega, amp, limit, Tp, gamma = 3.3, statistic = 'significant', duration = 10800 }) {
  needPos(limit, 'limit');
  if (statistic !== 'significant' && statistic !== 'mpm') throw new RangeError("statistic must be 'significant' or 'mpm'");
  const HsLimit = Tp.map((tp) => {
    const st = responseStatistics({ omega, amp, S: waveSpectrum(omega, { Hs: 1, Tp: tp, gamma }), duration });
    const r = statistic === 'mpm' ? st.mpm : st.significant;
    return r > 0 ? limit / r : Infinity;
  });
  return { Tp: Tp.slice(), HsLimit };
}

// sin(x)/x
const sinc = (x) => (Math.abs(x) < 1e-4 ? 1 - (x * x) / 6 : Math.sin(x) / x);

/**
 * Wave-induced heave, pitch and roll of a box-shaped ship or barge L × B × T in deep water.
 *
 * method 'panel' (default at zero speed): frequency-dependent added mass, wave damping and exciting forces of
 * the box from a three-dimensional panel solution with the free-surface Green function (radiation.js:
 * boxHydrodynamics, 14 frequencies up to the wavelength equal to the beam, interpolated, cached per geometry),
 * and the equations of motion about the waterline centre in three uncoupled groups,
 *   heave:        (m + A33) z'' + B33 z' + rho g L B z = X3,
 *   surge–pitch:  2 × 2 with the mass coupling m zG, I55 = m (kyy² + zG²), C55 = rho g (B L³/12 - V T/2) - m g zG,
 *   sway–roll:    2 × 2 with the mass coupling -m zG, C44 = Δ g GMt and the added roll damping below.
 * The centre of gravity follows from GMt for the box, zG = B²/(12 T) - T/2 - GMt (on the waterline when GMt is
 * not given, which then leaves out roll). For a block coefficient below 1 heave and pitch use the box of breadth
 * B·Cb, as Jensen, Mansour & Olsen (2004) propose (same length, draught and displacement), and roll uses the
 * box of breadth B with the mass Δ = rho L B T Cb.
 * Roll inertia: with `kxxDry` the mass inertia about the centre of gravity is Δ kxxDry² and the computed added
 * inertia is added at every frequency. With `kxx` (including added inertia; default 0.4 B) the roll natural
 * period is the input, 2 pi kxx / sqrt(g GMt): the mass inertia is Δ kxx² less the added inertia of the free
 * sway–roll system at that period, A44' - A24'²/(m + A22) about the centre of gravity (not less than 5 % of
 * Δ kxx²). `rollAddedInertia` returns that added inertia at the roll natural period.
 * Roll damping: wave radiation from the panel solution (`rollRadiation: false` removes the radiation damping
 * of the sway–roll pair) plus `rollDampingRatio`, a fraction of the critical damping 2 sqrt(I44 C44) with I44
 * the total inertia at the natural period, for what potential flow does not give (bilge keels, eddies, skin
 * friction; typical 0.02–0.10 for barges with sharp bilges or bilge keels). Typical kxx 0.35–0.45 B, kyy 0.25 L.
 *
 * method 'closed-form' (default with forward speed): expressions of Jensen, Mansour & Olsen (2004) for a
 * homogeneously loaded box (deep water, optional forward speed):
 *   alpha = 1 - Fn sqrt(kL) cos(beta), k_e = |k cos(beta)|, k_y = k sin(beta),
 *   A = 2 sin(½ k B alpha²) exp(-k T alpha²)  (sectional damping),
 *   f = sqrt((1 - kT)² + (A²/(k B alpha³))²),
 *   F = kappa f s_B (2/(k_e L)) sin(k_e L/2),
 *   G = kappa f s_B (24/((k_e L)² L)) [sin(k_e L/2) - (k_e L/2) cos(k_e L/2)],
 *   eta = 1/sqrt((1 - 2 k T alpha²)² + (A²/(k B alpha²))²),  heave = eta F [m/m], pitch = eta G [rad/m],
 * with B replaced by B·Cb for other block coefficients, as in the paper. These give amplitudes only.
 * Two factors differ from the paper, which approximates the Smith factor by exp(-k_e T) and has no breadth
 * factor: kappa = exp(-k T) and s_B = sin(½ k_y B)/(½ k_y B). Both follow from integrating the undisturbed
 * pressure rho g exp(-k T) exp(-i (k_x x + k_y y)) over the flat bottom L × B, which gives rho g L B exp(-k T)
 * sinc(½ k_x L) sinc(½ k_y B): the depth decay goes with the full wave number whatever the heading. In head and
 * following seas the result is that of the paper; in beam seas the paper's form tends to 1 + kT in long waves
 * (8 % at kT = 0.08) where the heave of a floating body tends to 1 + O((kT)²).
 * Closed-form roll: NOT taken from the paper. One degree of freedom about the centre of gravity,
 *   I44 phi'' + (B_rad + 2 zeta sqrt(I44 C44)) phi' + C44 phi = M_FK + M_D,   C44 = Δ g GMt,
 *   I44 = Δ kxx² (kxx including added inertia) or Δ kxxDry² + A44,
 * - M_FK: Froude–Krylov roll moment integrated exactly over the box in the undisturbed deep-water wave;
 * - A44 = ½ beta rho pi a⁴ L: strip roll added inertia of the rectangular section about the waterline centre
 *   (a = the larger of B/2 and T, beta from rectangleRollInertiaCoefficient: DNV-RP-C205 Table D-1 for the
 *   double body, i.e. the high-frequency limit);
 * - M_D = -omega² A44 alpha_b: diffraction by the relative-motion hypothesis, the added inertia reacting to
 *   the roll relative to the wave slope. alpha_b is the wave slope carried by the pressure on the bottom,
 *   (bottom-pressure roll moment)/(rho g L B³/12) = k_y exp(-kT) sinc(½ k_x L) 3 (sin x - x cos x)/x³, x = ½ k_y B;
 * - B_rad: wave-radiation damping from the Haskind relation with this moment in waves from all directions,
 *   B_rad = omega³/(4 pi rho g³) ∫₀^{2π} |M_FK + M_D|² d beta (deep water);
 * - no sway–roll coupling (A44 is not transferred from the waterline to the centre of gravity).
 *
 * Validity, from comparisons with a panel solution (Capytaine; src/spec/cross-solver.json, docs/engine/hydro.md).
 * 'panel', three boxes (20 × 8 × 2 m, 100 × 20 × 5 m, 80 × 24 × 3 m; head, quartering and beam seas): heave and
 * pitch within 5 % (5.3 % at one period at 135°) at every compared period from 0.72 times the heave natural period
 * upward and within 1.5 % from 1.15 times; -8 % where the wavelength equals the beam of the widest box. Roll
 * within 5 % from 0.75–0.85 times the roll natural period upward and within 11 % below, where it is small; the
 * resonant roll follows the gyradius and the added damping that are entered. Wavelengths shorter than the beam
 * are outside the solved range (the coefficients of the last solved frequency are held): `validPeriodMin` is the
 * period of the wave whose length equals the beam. Not represented: hull forms other than a box, finite depth,
 * forward speed.
 * 'closed-form': heave within 5 % for wave periods above twice 2 pi sqrt(2T/g) in head and beam seas (-15 % at
 * 1.5 times, -30 to -40 % near it: added mass taken equal to the displacement, two-dimensional sectional
 * damping); pitch within 5 % above 1.5 times that period; roll within 5 % above 1.3 times the roll natural
 * period and within 15 % below. `validPeriodMin` returns the limits of the method used.
 * Reference: Newman (1977) Marine Hydrodynamics ch. 6 and the references of radiation.js (panel solution);
 * Jensen, Mansour & Olsen (2004) Estimation of ship motions using closed-form expressions, Ocean Engineering 31,
 * 61–85 (closed-form heave, pitch); Korvin-Kroukovsky & Jacobs (1957) Pitching and heaving motions of a ship in
 * regular waves, Trans. SNAME 65 (relative-motion hypothesis); Newman (1962) J. Ship Research 6(3) (Haskind
 * relation); DNV-RP-C205 (2019) Appendix D. Deep water, monohull box-like forms; 'closed-form' for Fn < ~0.3.
 * @param {object} o
 * @param {number} o.L length, m
 * @param {number} o.B beam, m
 * @param {number} o.T draft, m
 * @param {ArrayLike<number>} o.omega wave frequencies, rad/s
 * @param {number} [o.headingDeg=180] propagation heading, deg (180 = head sea, 90 = beam sea)
 * @param {number} [o.speed=0] forward speed, m/s ('closed-form' only)
 * @param {number} [o.Cb=1] block coefficient
 * @param {number} [o.GMt] transverse metacentric height, m (roll omitted when absent)
 * @param {number} [o.kxx] roll radius of gyration including added inertia, m (default 0.4 B unless kxxDry is given)
 * @param {number} [o.kxxDry] roll radius of gyration of the mass alone, m; the added inertia is added
 * @param {number} [o.kyy] pitch radius of gyration of the mass about the centre of gravity, m (default 0.25 L;
 *   'panel' only, the closed form has no mass input)
 * @param {number} [o.rollDampingRatio=0.05] damping in addition to wave radiation, fraction of critical
 * @param {boolean} [o.rollRadiation=true] include the wave-radiation damping of roll
 * @param {'panel'|'closed-form'} [o.method] default 'panel' at zero speed, 'closed-form' with forward speed
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{omega:Float64Array, heave:Float64Array, pitch:Float64Array, roll:Float64Array|null,
 *   rollPhase:Float64Array|null, heavePhase:Float64Array|null, pitchPhase:Float64Array|null,
 *   heaveNaturalPeriod:number, rollNaturalPeriod:number|null, rollAddedInertia:number|null,
 *   rollRadiationDamping:Float64Array|null, validPeriodMin:{heave:number, roll:number|null}, method:string}}
 *   m/m, rad/m, rad/m; phases rad (leads on the crest at midship; heave and pitch phases for 'panel' only);
 *   natural periods s ('panel': with the added mass at that period; 'closed-form' heave: 2 pi sqrt(2T/g));
 *   kg·m²; N·m·s (about the centre of gravity); validPeriodMin s: shortest wave period of the checked range
 */
export function boxBargeRAO({ L, B, T, omega, headingDeg = 180, speed = 0, Cb = 1, GMt, kxx, kxxDry, kyy, rollDampingRatio = 0.05, rollRadiation = true, method, rho = 1025, g = G }) {
  needPos(L, 'L'); needPos(B, 'B'); needPos(T, 'T');
  if (!(Cb > 0 && Cb <= 1)) throw new RangeError('Cb must be within (0, 1]');
  if (!(speed >= 0)) throw new RangeError('speed must be >= 0');
  if (kxx !== undefined && kxxDry !== undefined) throw new RangeError('give kxx (including added inertia) or kxxDry, not both');
  if (!(rollDampingRatio >= 0)) throw new RangeError('rollDampingRatio must be >= 0');
  const mth = method === undefined ? (speed > 0 ? 'closed-form' : 'panel') : method;
  if (mth !== 'panel' && mth !== 'closed-form') throw new RangeError("method must be 'panel' or 'closed-form'");
  if (GMt !== undefined) needPos(GMt, 'GMt');
  const args = { L, B, T, omega, headingDeg, speed, Cb, GMt, kxx, kxxDry, kyy, rollDampingRatio, rollRadiation, rho, g };
  if (mth === 'closed-form') return boxBargeClosedForm(args);
  if (speed > 0) throw new RangeError("method 'panel' is a zero-speed solution: use method 'closed-form' with forward speed");
  return boxBargePanel(args);
}

// complex division (ar + i ai)/(br + i bi)
function cdiv(ar, ai, br, bi, out) {
  const d = br * br + bi * bi;
  out.re = (ar * br + ai * bi) / d; out.im = (ai * br - ar * bi) / d;
  return out;
}

// solves the symmetric complex system [[a, b], [b, d]] x = [f1, f2]; returns x2 in out and x1 in out1
function csolve2(ar, ai, br, bi, dr, di, f1r, f1i, f2r, f2i, out1, out) {
  const detR = ar * dr - ai * di - (br * br - bi * bi); const detI = ar * di + ai * dr - 2 * br * bi;
  cdiv(f1r * dr - f1i * di - (br * f2r - bi * f2i), f1r * di + f1i * dr - (br * f2i + bi * f2r), detR, detI, out1);
  cdiv(ar * f2r - ai * f2i - (br * f1r - bi * f1i), ar * f2i + ai * f2r - (br * f1i + bi * f1r), detR, detI, out);
  return out;
}

function boxBargePanel({ L, B, T, omega, headingDeg, Cb, GMt, kxx, kxxDry, kyy, rollDampingRatio, rollRadiation, rho, g }) {
  const nw = omega.length;
  const hasRoll = GMt !== undefined;
  const B0 = B * Cb;
  const disp = rho * L * B * T * Cb;
  const zG = hasRoll ? (B * B) / (12 * T) - T / 2 - GMt : 0;
  const ky = kyy === undefined ? 0.25 * L : kyy;
  needPos(ky, 'kyy');
  // one model with all modes when the two boxes coincide; otherwise each box solves only the modes taken from it
  const hp = boxHydrodynamics({ L, B: B0, T, rho, g, modes: hasRoll && Cb === 1 ? [0, 1, 2, 3, 4] : [0, 2, 4] });
  const hr = hasRoll ? (Cb === 1 ? hp : boxHydrodynamics({ L, B, T, rho, g, modes: [1, 3] })) : null;
  const C33 = rho * g * L * B0;
  const I55 = disp * (ky * ky + zG * zG);
  const C55 = rho * g * ((B0 * L ** 3) / 12 - L * B0 * T * 0.5 * T) - disp * g * zG;
  if (!(C55 > 0)) throw new RangeError('the box has no positive pitch stiffness with this GMt (centre of gravity too high)');
  // natural frequency of a one-degree-of-freedom oscillator whose added inertia depends on frequency
  const natural = (stiff, inertia, added) => {
    let w = Math.sqrt(stiff / (inertia + added(Math.sqrt(stiff / (2 * inertia)))));
    for (let it = 0; it < 60; it++) { const wn = Math.sqrt(stiff / (inertia + added(w))); const dw = wn - w; w += 0.5 * dw; if (Math.abs(dw) < 1e-10 * w) break; }
    return w;
  };
  const heaveNaturalPeriod = (2 * Math.PI) / natural(C33, disp, (w) => boxCoefficients(hp, w).A[2][2]);
  // added inertia of the free sway–roll pair about the centre of gravity (roll about G = roll + zG × sway)
  const rollAdded = (w) => {
    const { A } = boxCoefficients(hr, w);
    const a44 = A[3][3] + 2 * zG * A[1][3] + zG * zG * A[1][1]; const a24 = A[1][3] + zG * A[1][1];
    return a44 - (a24 * a24) / (disp + A[1][1]);
  };
  let IG = 0; let I44 = 0; let C44 = 0; let wRoll = 0; let B44 = 0;
  if (hasRoll) {
    C44 = disp * g * GMt;
    if (kxxDry !== undefined) { needPos(kxxDry, 'kxxDry'); IG = disp * kxxDry * kxxDry; } else {
      const kx = kxx === undefined ? 0.4 * B : kxx;
      needPos(kx, 'kxx');
      IG = Math.max(disp * kx * kx - rollAdded(Math.sqrt((g * GMt) / (kx * kx))), 0.05 * disp * kx * kx);
    }
    wRoll = natural(C44, IG, rollAdded);
    I44 = IG + rollAdded(wRoll);
    B44 = 2 * rollDampingRatio * Math.sqrt(I44 * C44);
  }
  const heave = new Float64Array(nw); const pitch = new Float64Array(nw);
  const heavePhase = new Float64Array(nw); const pitchPhase = new Float64Array(nw);
  const roll = hasRoll ? new Float64Array(nw) : null; const rollPhase = hasRoll ? new Float64Array(nw) : null;
  const rollRad = hasRoll ? new Float64Array(nw) : null;
  const x1 = { re: 0, im: 0 }; const x2 = { re: 0, im: 0 };
  for (let i = 0; i < nw; i++) {
    const w = omega[i];
    if (!(w > 0)) throw new RangeError('omega must be > 0');
    const w2 = w * w;
    const { A, B: D } = boxCoefficients(hp, w);
    const X = boxExcitation(hp, w, headingDeg);
    cdiv(X.re[2], X.im[2], C33 - w2 * (disp + A[2][2]), w * D[2][2], x2);
    heave[i] = Math.hypot(x2.re, x2.im); heavePhase[i] = Math.atan2(x2.im, x2.re);
    csolve2(-w2 * (disp + A[0][0]), w * D[0][0], -w2 * (disp * zG + A[0][4]), w * D[0][4], C55 - w2 * (I55 + A[4][4]), w * D[4][4],
      X.re[0], X.im[0], X.re[4], X.im[4], x1, x2);
    pitch[i] = Math.hypot(x2.re, x2.im); pitchPhase[i] = Math.atan2(x2.im, x2.re);
    if (hasRoll) {
      const cr = hr === hp ? { A, B: D } : boxCoefficients(hr, w);
      const Xr = hr === hp ? X : boxExcitation(hr, w, headingDeg);
      const Ar = cr.A; const Dr = cr.B; const q = rollRadiation ? w : 0;
      rollRad[i] = rollRadiation ? Dr[3][3] + 2 * zG * Dr[1][3] + zG * zG * Dr[1][1] : 0;
      csolve2(-w2 * (disp + Ar[1][1]), q * Dr[1][1], -w2 * (-disp * zG + Ar[1][3]), q * Dr[1][3], C44 - w2 * (IG + disp * zG * zG + Ar[3][3]), q * Dr[3][3] + w * B44,
        Xr.re[1], Xr.im[1], Xr.re[3], Xr.im[3], x1, x2);
      roll[i] = Math.hypot(x2.re, x2.im); rollPhase[i] = Math.atan2(x2.im, x2.re);
    }
  }
  return {
    omega: Float64Array.from(omega), heave, pitch, roll, rollPhase, heavePhase, pitchPhase, heaveNaturalPeriod,
    rollNaturalPeriod: hasRoll ? (2 * Math.PI) / wRoll : null, rollAddedInertia: hasRoll ? I44 - IG : null, rollRadiationDamping: rollRad,
    validPeriodMin: { heave: (2 * Math.PI) / hp.omegaMax, roll: hasRoll ? (2 * Math.PI) / hr.omegaMax : null }, method: 'panel',
  };
}

function boxBargeClosedForm({ L, B, T, omega, headingDeg, speed, Cb, GMt, kxx, kxxDry, rollDampingRatio, rollRadiation, rho, g }) {
  const nw = omega.length;
  const beta = headingDeg * DEG;
  const cbt = Math.cos(beta); const sbt = Math.sin(beta);
  const B0 = B * Cb;
  const Fn = speed / Math.sqrt(g * L);
  const heave = new Float64Array(nw); const pitch = new Float64Array(nw);
  const hasRoll = GMt !== undefined;
  const roll = hasRoll ? new Float64Array(nw) : null;
  const rollPhase = hasRoll ? new Float64Array(nw) : null;
  const rollRad = hasRoll ? new Float64Array(nw) : null;
  const disp = rho * L * B * T * Cb;
  const b = B / 2;
  const A44 = 0.5 * rectangleRollInertiaCoefficient(Math.max(b / T, T / b)) * rho * Math.PI * Math.max(b, T) ** 4 * L;
  let I44;
  if (kxxDry !== undefined) { needPos(kxxDry, 'kxxDry'); I44 = disp * kxxDry * kxxDry + A44; } else {
    const kx = kxx === undefined ? 0.4 * B : kxx;
    needPos(kx, 'kxx');
    I44 = disp * kx * kx;
  }
  const C44 = hasRoll ? disp * g * GMt : 0;
  const B44 = 2 * rollDampingRatio * Math.sqrt(I44 * C44);
  // centre of gravity height relative to the waterline for a box: KG = T/2 + B²/(12 T) - GMt
  const zG = hasRoll ? T / 2 + (B * B) / (12 * T) - GMt - T : 0;
  // Imaginary part of the roll moment about G per unit wave amplitude (the moment is i·M) for the wave
  // exp(kz) exp(-i(kx x + ky y)) travelling towards (c, s):
  // Froude–Krylov  M_O = rho (∫∫∫ (y a_z - z a_y) dV) + rho g ∫∫ y eta dA,  M_G = M_O + z_G F_y,
  // plus the relative-motion diffraction moment -omega² A44 alpha_b.
  const rollMoment = (w, k, c, s) => {
    const kxw = k * c; const ky = k * s;
    const Lx = L * sinc(0.5 * kxw * L);
    let Iy0; let Iy1i; // Iy0 = ∫ exp(-i ky y) dy (real); ∫ y exp(-i ky y) dy = i·Iy1i (purely imaginary)
    if (Math.abs(ky) * b < 1e-6) { Iy0 = 2 * b; Iy1i = (-2 * ky * b ** 3) / 3; } else {
      Iy0 = (2 * Math.sin(ky * b)) / ky;
      Iy1i = (2 * (ky * b * Math.cos(ky * b) - Math.sin(ky * b))) / (ky * ky);
    }
    const e = Math.exp(-k * T);
    const Iz0 = (1 - e) / k;
    const Iz1 = (-1 + e + k * T * e) / (k * k);
    // a_z = -w² e^{kz} E, a_y = i w² sin(beta) e^{kz} E
    const mIm = rho * Lx * (-w * w * Iy1i * Iz0 - w * w * s * Iy0 * Iz1) + rho * g * Lx * Iy1i;
    const fyIm = rho * w * w * s * Lx * Iy0 * Iz0;
    // bottom pressure rho g e^{-kT} E gives the moment i rho g e^{-kT} Lx Iy1i = -i rho g (L B³/12) alpha_b
    const slopeIm = (e * Lx * Iy1i) / ((L * B ** 3) / 12);
    return mIm + zG * fyIm - w * w * A44 * slopeIm;
  };
  for (let i = 0; i < nw; i++) {
    const w = omega[i];
    if (!(w > 0)) throw new RangeError('omega must be > 0');
    const k = (w * w) / g;
    const alpha = 1 - Fn * Math.sqrt(k * L) * cbt;
    const ke = Math.abs(k * cbt);
    const kappa = Math.exp(-k * T) * sinc(0.5 * k * sbt * B0);
    const Asec = 2 * Math.sin(0.5 * k * B0 * alpha * alpha) * Math.exp(-k * T * alpha * alpha);
    const f = Math.hypot(1 - k * T, (Asec * Asec) / (k * B0 * alpha ** 3));
    const x = 0.5 * ke * L;
    const F = x < 1e-6 ? kappa * f : (kappa * f * Math.sin(x)) / x;
    const Gp = x < 1e-4 ? (kappa * f * 2 * x) / L : (kappa * f * 24 * (Math.sin(x) - x * Math.cos(x))) / ((ke * L) ** 2 * L);
    const eta = 1 / Math.hypot(1 - 2 * k * T * alpha * alpha, (Asec * Asec) / (k * B0 * alpha * alpha));
    // amplitudes only: the sin(x)/x force factors change sign beyond their first zero, which is a phase flip, not a negative amplitude
    heave[i] = Math.abs(eta * F); pitch[i] = Math.abs(eta * Gp);
    if (hasRoll) {
      const Mim = rollMoment(w, k, cbt, sbt);
      let brad = 0;
      if (rollRadiation) {
        // |M|² is even in both wave-number components: midpoint rule over one quadrant, enough points for the
        // phase variation along the longer side
        const nq = Math.min(2000, Math.max(16, Math.ceil(0.5 * k * Math.max(L, B)) + 8));
        let sum = 0;
        for (let j = 0; j < nq; j++) {
          const th = ((j + 0.5) * Math.PI) / (2 * nq);
          sum += rollMoment(w, k, Math.cos(th), Math.sin(th)) ** 2;
        }
        brad = ((w ** 3) / (4 * Math.PI * rho * g ** 3)) * ((4 * sum * Math.PI) / (2 * nq));
        rollRad[i] = brad;
      }
      const dr = C44 - w * w * I44; const di = w * (B44 + brad);
      const den = dr * dr + di * di;
      const re = (Mim * di) / den; const im = (Mim * dr) / den;
      roll[i] = Math.hypot(re, im); rollPhase[i] = Math.atan2(im, re);
    }
  }
  const heaveNaturalPeriod = 2 * Math.PI * Math.sqrt((2 * T) / g);
  const rollNaturalPeriod = hasRoll ? 2 * Math.PI * Math.sqrt(I44 / C44) : null;
  return {
    omega: Float64Array.from(omega), heave, pitch, roll, rollPhase, heavePhase: null, pitchPhase: null, heaveNaturalPeriod, rollNaturalPeriod,
    rollAddedInertia: hasRoll ? A44 : null, rollRadiationDamping: rollRad,
    validPeriodMin: { heave: 2 * heaveNaturalPeriod, roll: hasRoll ? 1.3 * rollNaturalPeriod : null }, method: 'closed-form',
  };
}
