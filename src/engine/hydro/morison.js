// Morison-type loads on slender members, coefficient models, MacCamy–Fuchs diffraction, added-mass data,
// member buoyancy, slam and vortex-shedding helpers.
// Axes: z up, z = 0 at still-water level, seabed at z = -depth. See _wavehelp.js for the wave convention.

import { gaussLegendre } from '../math/quad.js';
import { besselJ0, besselJ1, besselY0, besselY1 } from '../math/special.js';
import { G, prepWave, airyPoint, airyEta, dispersion, needPos, needVec3 } from './_wavehelp.js';

/**
 * Morison force per unit length on a cylinder of arbitrary orientation (normal-component formulation,
 * relative-velocity form when body motion is given):
 *   f = rho A a_n + rho Ca A (a_n - ab_n) + ½ rho Cd D |v_n - vb_n| (v_n - vb_n),   Cm = 1 + Ca,
 * where the subscript n denotes the component normal to the member axis.
 * Reference: DNV-RP-C205 (2019) §6.2. Valid for slender members, D/lambda < 0.2.
 * @param {object} o
 * @param {number} [o.rho=1025] fluid density, kg/m³
 * @param {number} o.D diameter, m
 * @param {number} o.Cd drag coefficient, - (typical 0.65 smooth … 1.05 rough, see dragCoefficientSteady)
 * @param {number} o.Cm inertia coefficient, - (typical 1.6 smooth … 2.0)
 * @param {number} [o.Ca=Cm-1] added-mass coefficient, -
 * @param {number[]} o.axis member axis direction (any length), [x, y, z]
 * @param {number[]} o.u fluid velocity, m/s
 * @param {number[]} o.a fluid acceleration, m/s²
 * @param {number[]} [o.vBody] member velocity, m/s
 * @param {number[]} [o.aBody] member acceleration, m/s²
 * @returns {{f:number[], fDrag:number[], fInertia:number[]}} force per unit length, N/m (total, drag, inertia)
 */
export function morisonForcePerLength({ rho = 1025, D, Cd, Cm, Ca, axis, u, a, vBody, aBody }) {
  needPos(rho, 'rho'); needPos(D, 'D');
  if (!(Cd >= 0)) throw new RangeError('Cd must be >= 0');
  if (!(Cm >= 0)) throw new RangeError('Cm must be >= 0');
  needVec3(axis, 'axis'); needVec3(u, 'u'); needVec3(a, 'a');
  const ca = Ca === undefined ? Cm - 1 : Ca;
  const L = Math.hypot(axis[0], axis[1], axis[2]);
  if (!(L > 0)) throw new RangeError('axis must be non-zero');
  const t = [axis[0] / L, axis[1] / L, axis[2] / L];
  const vb = vBody || [0, 0, 0];
  const ab = aBody || [0, 0, 0];
  const normal = (v) => {
    const d = v[0] * t[0] + v[1] * t[1] + v[2] * t[2];
    return [v[0] - d * t[0], v[1] - d * t[1], v[2] - d * t[2]];
  };
  const vr = normal([u[0] - vb[0], u[1] - vb[1], u[2] - vb[2]]);
  const an = normal(a);
  const abn = normal(ab);
  const A = (Math.PI * D * D) / 4;
  const vm = Math.hypot(vr[0], vr[1], vr[2]);
  const kd = 0.5 * rho * Cd * D * vm;
  const fDrag = [kd * vr[0], kd * vr[1], kd * vr[2]];
  const fInertia = [0, 1, 2].map((i) => rho * A * an[i] + rho * ca * A * (an[i] - abn[i]));
  return { f: [fDrag[0] + fInertia[0], fDrag[1] + fInertia[1], fDrag[2] + fInertia[2]], fDrag, fInertia };
}

const _kin = { ux: 0, uy: 0, uz: 0, ax: 0, ay: 0, az: 0, eta: 0, pd: 0 };
const STRETCH = ['none', 'wheeler', 'vertical'];

function memberSetup(o) {
  const { p1, p2, D, Cd, Cm, rho = 1025, g = G, wave, current, stretching = 'none', nGauss = 8, nSeg = 4, about } = o;
  needVec3(p1, 'p1'); needVec3(p2, 'p2'); needPos(D, 'D'); needPos(rho, 'rho');
  if (!(Cd >= 0)) throw new RangeError('Cd must be >= 0');
  if (!(Cm >= 0)) throw new RangeError('Cm must be >= 0');
  if (!STRETCH.includes(stretching)) throw new RangeError("stretching must be 'none', 'wheeler' or 'vertical'");
  if (!(nGauss >= 1 && nGauss <= 64)) throw new RangeError('nGauss must be 1..64');
  if (!(nSeg >= 1)) throw new RangeError('nSeg must be >= 1');
  const w = prepWave(wave, g);
  const d = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
  const L = Math.hypot(d[0], d[1], d[2]);
  if (!(L > 0)) throw new RangeError('p1 and p2 must differ');
  let cur = null;
  if (current) {
    const ch = ((current.headingDeg || 0) * Math.PI) / 180;
    const prof = typeof current.speed === 'function' ? current.speed : () => current.speed;
    cur = { prof, c: Math.cos(ch), s: Math.sin(ch) };
  }
  return {
    p1, d, L, t: [d[0] / L, d[1] / L, d[2] / L], D, Cd, Cm, rho, w, cur, stretching,
    gl: gaussLegendre(nGauss), nSeg, about: about || [0, 0, 0],
  };
}

// Wetted parameter range [sa, sb] of the member (s = 0 at p1, 1 at p2) at time t. A single free-surface
// crossing along the member is assumed (true for members short compared with the wavelength).
function wetRange(m, t, out) {
  const { p1, d, w } = m;
  const surf = m.stretching === 'none'
    ? (s) => p1[2] + s * d[2]
    : (s) => p1[2] + s * d[2] - airyEta(w, p1[0] + s * d[0], p1[1] + s * d[1], t);
  const g0 = surf(0); const g1 = surf(1);
  if (g0 <= 0 && g1 <= 0) { out[0] = 0; out[1] = 1; return; }
  if (g0 > 0 && g1 > 0) { out[0] = 0; out[1] = 0; return; }
  let lo = 0; let hi = 1;
  const wetAtLo = g0 <= 0;
  for (let i = 0; i < 50; i++) {
    const mid = 0.5 * (lo + hi);
    if ((surf(mid) <= 0) === wetAtLo) lo = mid; else hi = mid;
  }
  const sc = 0.5 * (lo + hi);
  if (wetAtLo) { out[0] = 0; out[1] = sc; } else { out[0] = sc; out[1] = 1; }
}

const _rng = [0, 0];

function memberForceAt(m, t, out) {
  const { p1, d, L, w, cur, gl } = m;
  const tx = m.t[0]; const ty = m.t[1]; const tz = m.t[2];
  const A = (Math.PI * m.D * m.D) / 4;
  const kI = m.rho * m.Cm * A;
  const kD = 0.5 * m.rho * m.Cd * m.D;
  wetRange(m, t, _rng);
  let Fx = 0; let Fy = 0; let Fz = 0; let Mx = 0; let My = 0; let Mz = 0;
  // the 'vertical' stretching has a slope discontinuity at z = 0: integrate each side separately
  let sBreak = -1;
  if (m.stretching === 'vertical' && d[2] !== 0) {
    const sz = -p1[2] / d[2];
    if (sz > _rng[0] + 1e-12 && sz < _rng[1] - 1e-12) sBreak = sz;
  }
  const nPart = sBreak < 0 ? 1 : 2;
  for (let part = 0; part < nPart; part++) {
    const sa = part === 0 ? _rng[0] : sBreak;
    const sb = nPart === 1 || part === 1 ? _rng[1] : sBreak;
    const ds = (sb - sa) / m.nSeg;
    if (!(ds > 0)) continue;
    for (let e = 0; e < m.nSeg; e++) {
      const s0 = sa + e * ds;
      for (let q = 0; q < gl.x.length; q++) {
        const s = s0 + 0.5 * ds * (gl.x[q] + 1);
        const wt = 0.5 * ds * gl.w[q] * L;
        const x = p1[0] + s * d[0]; const y = p1[1] + s * d[1]; const z = p1[2] + s * d[2];
        let zEff = z;
        if (m.stretching === 'wheeler') {
          const eta = airyEta(w, x, y, t);
          // Wheeler (1970): the still-water profile is stretched to the instantaneous surface.
          zEff = Number.isFinite(w.depth) ? (w.depth * (z - eta)) / (w.depth + eta) : z - eta;
        } else if (z > 0) zEff = 0;
        airyPoint(w, x, y, zEff, t, _kin);
        let ux = _kin.ux; let uy = _kin.uy; const uz = _kin.uz;
        if (cur) { const U = cur.prof(z); ux += U * cur.c; uy += U * cur.s; }
        const ut = ux * tx + uy * ty + uz * tz;
        const vx = ux - ut * tx; const vy = uy - ut * ty; const vz = uz - ut * tz;
        const at = _kin.ax * tx + _kin.ay * ty + _kin.az * tz;
        const vm = Math.sqrt(vx * vx + vy * vy + vz * vz);
        const fx = kD * vm * vx + kI * (_kin.ax - at * tx);
        const fy = kD * vm * vy + kI * (_kin.ay - at * ty);
        const fz = kD * vm * vz + kI * (_kin.az - at * tz);
        const rx = x - m.about[0]; const ry = y - m.about[1]; const rz = z - m.about[2];
        Fx += wt * fx; Fy += wt * fy; Fz += wt * fz;
        Mx += wt * (ry * fz - rz * fy); My += wt * (rz * fx - rx * fz); Mz += wt * (rx * fy - ry * fx);
      }
    }
  }
  out[0] = Fx; out[1] = Fy; out[2] = Fz; out[3] = Mx; out[4] = My; out[5] = Mz;
}

/**
 * Total Morison force and moment on a straight cylindrical member between two points at one instant, in Airy
 * wave kinematics with optional steady current, integrated along the wetted length with Gauss–Legendre points.
 * Stretching: 'none' integrates to z = 0; 'wheeler' integrates to the instantaneous surface with Wheeler
 * stretching; 'vertical' integrates to the surface with the z = 0 kinematics held constant above still water.
 * Reference: DNV-RP-C205 (2019) §6.2 and §3 (kinematics, stretching); Wheeler (1970). Valid for D/lambda < 0.2 and non-breaking waves;
 * axial (end) forces and hydrostatic buoyancy are not included.
 * @param {object} o
 * @param {number[]} o.p1 end 1, m
 * @param {number[]} o.p2 end 2, m
 * @param {number} o.D diameter, m
 * @param {number} o.Cd drag coefficient, -
 * @param {number} o.Cm inertia coefficient, -
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @param {{H:number,T:number,depth?:number,headingDeg?:number,phaseDeg?:number}} o.wave m, s, m, deg, deg
 * @param {{speed:number|function(number):number, headingDeg?:number}} [o.current] speed m/s (or profile of z)
 * @param {number} [o.t=0] time, s (crest at the origin at t = 0)
 * @param {string} [o.stretching='none']
 * @param {number} [o.nGauss=8] Gauss points per segment
 * @param {number} [o.nSeg=4] integration segments along the wetted length
 * @param {number[]} [o.about=[0,0,0]] moment reference point, m
 * @returns {{F:number[], M:number[]}} force N and moment N·m about `about`
 */
export function memberWaveForce(o) {
  const m = memberSetup(o);
  const r = new Array(6);
  memberForceAt(m, o.t || 0, r);
  return { F: [r[0], r[1], r[2]], M: [r[3], r[4], r[5]] };
}

/**
 * Time history of the Morison force and moment on a member over one wave cycle, with maxima.
 * `Fh` is the horizontal force resolved along the wave heading and `Mh` the moment about the horizontal axis
 * normal to the wave heading (overturning moment about `about`).
 * Reference and validity as memberWaveForce.
 * @param {object} o as memberWaveForce, plus
 * @param {number} [o.nSteps=72] time steps per wave period
 * @returns {{t:Float64Array, Fx:Float64Array, Fy:Float64Array, Fz:Float64Array, Mx:Float64Array, My:Float64Array,
 *   Mz:Float64Array, Fh:Float64Array, Mh:Float64Array, FhMax:number, tFhMax:number, MhMax:number, tMhMax:number,
 *   FzMax:number}} times s, forces N, moments N·m; maxima are of absolute value
 */
export function memberForceHistory(o) {
  const m = memberSetup(o);
  const n = o.nSteps || 72;
  if (!(n >= 4)) throw new RangeError('nSteps must be >= 4');
  const T = (2 * Math.PI) / m.w.omega;
  const res = { t: new Float64Array(n) };
  for (const key of ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz', 'Fh', 'Mh']) res[key] = new Float64Array(n);
  const r = new Array(6);
  let FhMax = 0; let tF = 0; let MhMax = 0; let tM = 0; let FzMax = 0;
  for (let i = 0; i < n; i++) {
    const t = (i * T) / n;
    memberForceAt(m, t, r);
    res.t[i] = t; res.Fx[i] = r[0]; res.Fy[i] = r[1]; res.Fz[i] = r[2];
    res.Mx[i] = r[3]; res.My[i] = r[4]; res.Mz[i] = r[5];
    const fh = r[0] * m.w.cb + r[1] * m.w.sb;
    const mh = -r[3] * m.w.sb + r[4] * m.w.cb;
    res.Fh[i] = fh; res.Mh[i] = mh;
    if (Math.abs(fh) > FhMax) { FhMax = Math.abs(fh); tF = t; }
    if (Math.abs(mh) > MhMax) { MhMax = Math.abs(mh); tM = t; }
    FzMax = Math.max(FzMax, Math.abs(r[2]));
  }
  return { ...res, FhMax, tFhMax: tF, MhMax, tMhMax: tM, FzMax };
}

/**
 * Closed-form Morison base shear and mudline overturning moment on a vertical surface-piercing pile in a
 * linear wave (kinematics integrated from the seabed to z = 0):
 *   F_I = rho Cm (pi D²/4) g (H/2) tanh(kd),
 *   F_D = ½ rho Cd D (omega H/2)² [sinh(2kd)/(4k) + d/2] / sinh²(kd),
 *   F(theta) = F_D cos(theta)|cos(theta)| + F_I sin(theta); F_max = F_I if F_I >= 2 F_D, else F_D + F_I²/(4 F_D).
 * Moments use the first moments of the same depth functions about the mudline.
 * Reference: Morison et al. (1950); e.g. Dean & Dalrymple (1991) Water Wave Mechanics for Engineers and Scientists. Valid for D/lambda < 0.2, no current.
 * @param {{H:number,T:number,depth:number,D:number,Cd:number,Cm:number,rho?:number,g?:number}} o m, s, m, m, -, -, kg/m³, m/s²
 * @returns {{FD:number, FI:number, Fmax:number, MD:number, MI:number, Mmax:number, k:number, phaseFmaxDeg:number}}
 *   drag and inertia force amplitudes N, maximum total force N, moment amplitudes about the mudline N·m,
 *   wave number rad/m, phase omega·t of the maximum force after the crest (0 = crest, 90 = down-crossing), deg
 */
export function pileForceClosedForm({ H, T, depth, D, Cd, Cm, rho = 1025, g = G }) {
  needPos(H, 'H'); needPos(T, 'T'); needPos(depth, 'depth'); needPos(D, 'D');
  const w = (2 * Math.PI) / T;
  const k = dispersion(w, depth, g);
  const kd = k * depth;
  const A = (Math.PI * D * D) / 4;
  const sh = Math.sinh(kd);
  const FI = rho * Cm * A * g * (H / 2) * Math.tanh(kd);
  const cD = (0.5 * rho * Cd * D * (w * H / 2) ** 2) / (sh * sh);
  const FD = cD * (Math.sinh(2 * kd) / (4 * k) + depth / 2);
  const MI = ((rho * Cm * A * w * w * (H / 2)) / sh) * ((depth * sh) / k - (Math.cosh(kd) - 1) / (k * k));
  const MD = cD * ((depth * depth) / 4 + (depth * Math.sinh(2 * kd)) / (4 * k) - (Math.cosh(2 * kd) - 1) / (8 * k * k));
  const combine = (d, i) => (d === 0 || i >= 2 * d ? i : d + (i * i) / (4 * d));
  const sinTh = FD === 0 || FI >= 2 * FD ? 1 : FI / (2 * FD);
  return { FD, FI, Fmax: combine(FD, FI), MD, MI, Mmax: combine(MD, MI), k, phaseFmaxDeg: (Math.asin(sinTh) * 180) / Math.PI };
}

/**
 * Steady-flow drag coefficient of a circular cylinder at high Reynolds number as a function of relative
 * roughness Delta = k/D: 0.65 for Delta < 1e-4 (smooth); (29 + 4 log10 Delta)/20 for 1e-4 <= Delta <= 1e-2;
 * 1.05 for Delta > 1e-2 (rough).
 * Reference: DNV-RP-C205 (2019) §6.7.1.5. Valid for Re > 1e6 and large KC.
 * @param {number} roughnessRatio Delta = k/D, -
 * @returns {number} C_DS, -
 */
export function dragCoefficientSteady(roughnessRatio) {
  if (!(roughnessRatio >= 0)) throw new RangeError('roughnessRatio must be >= 0');
  if (roughnessRatio < 1e-4) return 0.65;
  if (roughnessRatio > 1e-2) return 1.05;
  return (29 + 4 * Math.log10(roughnessRatio)) / 20;
}

/**
 * Wake amplification factor psi(KC) for the drag coefficient in oscillatory flow, C_D = C_DS psi(KC), for KC < 12:
 *   C_pi = 1.50 - 0.024 (12/C_DS - 10);
 *   psi = C_pi + 0.10 (KC - 12) for 2 <= KC < 12; C_pi - 1.00 for 0.75 <= KC < 2;
 *   C_pi - 1.00 - 2.00 (KC - 0.75) for KC < 0.75.
 * Reference: DNV-RP-C205 (2019) §6.7.2.2. For KC >= 12 the standard gives psi only as a curve (tending to 1 at
 * large KC); this function throws there and the caller supplies the value.
 * @param {number} KC Keulegan–Carpenter number, -
 * @param {number} CDS steady drag coefficient, -
 * @returns {number} psi, -
 */
export function wakeAmplification(KC, CDS) {
  if (!(KC >= 0)) throw new RangeError('KC must be >= 0');
  if (KC >= 12) throw new RangeError('KC must be < 12 (for KC >= 12 psi is given graphically in DNV-RP-C205 Fig. 6-5)');
  needPos(CDS, 'CDS');
  const cpi = 1.5 - 0.024 * (12 / CDS - 10);
  if (KC >= 2) return cpi + 0.1 * (KC - 12);
  if (KC >= 0.75) return cpi - 1;
  return cpi - 1 - 2 * (KC - 0.75);
}

/**
 * Added-mass and inertia coefficients of a circular cylinder versus KC:
 *   C_A = 1.0 for KC < 3; C_A = max(1.0 - 0.044 (KC - 3), 0.6 - (C_DS - 0.65)) for KC >= 3; C_M = 1 + C_A.
 * Reference: DNV-RP-C205 (2019) §6.9.1. Asymptotes: C_A = 0.6 smooth, 0.2 rough.
 * @param {number} KC Keulegan–Carpenter number, -
 * @param {number} CDS steady drag coefficient (0.65 … 1.05), -
 * @returns {{Ca:number, Cm:number}}
 */
export function inertiaCoefficient(KC, CDS) {
  if (!(KC >= 0)) throw new RangeError('KC must be >= 0');
  needPos(CDS, 'CDS');
  const Ca = KC < 3 ? 1 : Math.max(1 - 0.044 * (KC - 3), 0.6 - (CDS - 0.65));
  return { Ca, Cm: 1 + Ca };
}

/**
 * Keulegan–Carpenter number KC = Um T / D.
 * Reference: DNV-RP-C205 (2019) §6.
 * @param {number} Um velocity amplitude, m/s
 * @param {number} T period, s
 * @param {number} D diameter, m
 * @returns {number} KC, -
 */
export function keuleganCarpenter(Um, T, D) {
  needPos(T, 'T'); needPos(D, 'D');
  return (Math.abs(Um) * T) / D;
}

/**
 * Effective diameter and added mass per unit length of a member with marine growth of uniform thickness.
 * First-principles geometry. The growth density is an input; 1325 kg/m³ is the value commonly used with
 * NORSOK N-003 thickness tables (replace with site data).
 * @param {{D:number, thickness:number, rhoGrowth?:number, rhoWater?:number, g?:number}} o m, m, kg/m³, kg/m³, m/s²
 * @returns {{De:number, massPerLength:number, submergedWeightPerLength:number}}
 *   effective diameter m, growth mass kg/m, growth submerged weight N/m
 */
export function marineGrowth({ D, thickness, rhoGrowth = 1325, rhoWater = 1025, g = G }) {
  needPos(D, 'D');
  if (!(thickness >= 0)) throw new RangeError('thickness must be >= 0');
  const De = D + 2 * thickness;
  const area = (Math.PI / 4) * (De * De - D * D);
  return { De, massPerLength: rhoGrowth * area, submergedWeightPerLength: (rhoGrowth - rhoWater) * area * g };
}

/**
 * MacCamy–Fuchs linear diffraction force on a bottom-mounted, surface-piercing vertical circular cylinder:
 *   F = (2 rho g H / k²) A(ka) tanh(kd) cos(omega t - alpha),
 *   A(ka) = [J1'(ka)² + Y1'(ka)²]^(-1/2),  tan(alpha) = J1'(ka)/Y1'(ka),
 * the equivalent Morison inertia coefficient C_M,eff = 4 A(ka) / (pi (ka)²) (-> 2 as ka -> 0) and the moment
 * about the mudline M = F d [kd sinh(kd) - cosh(kd) + 1]/(kd sinh(kd)).
 * Reference: MacCamy & Fuchs (1954), Beach Erosion Board Tech. Memo 69; DNV-RP-C205 (2019) §7.
 * Valid for linear waves, inviscid flow; use when D/lambda > 0.2.
 * @param {{H:number,T:number,depth:number,D:number,rho?:number,g?:number}} o m, s, m, m, kg/m³, m/s²
 * @returns {{F:number, M:number, CmEff:number, phaseLagDeg:number, ka:number, A:number, k:number,
 *   forcePerLength:function(number):number}} force amplitude N, mudline moment amplitude N·m, effective C_M,
 *   phase alpha deg, ka, A(ka), k rad/m, and the force amplitude per unit length N/m at elevation z
 */
export function maccamyFuchs({ H, T, depth, D, rho = 1025, g = G }) {
  needPos(H, 'H'); needPos(T, 'T'); needPos(depth, 'depth'); needPos(D, 'D');
  if (!Number.isFinite(depth)) throw new RangeError('depth must be finite');
  const k = dispersion((2 * Math.PI) / T, depth, g);
  const ka = (k * D) / 2;
  const dJ = besselJ0(ka) - besselJ1(ka) / ka;
  const dY = besselY0(ka) - besselY1(ka) / ka;
  const A = 1 / Math.hypot(dJ, dY);
  const kd = k * depth;
  const F = ((2 * rho * g * H) / (k * k)) * A * Math.tanh(kd);
  const arm = kd > 300 ? depth * (1 - 1 / kd) : (depth * (kd * Math.sinh(kd) - Math.cosh(kd) + 1)) / (kd * Math.sinh(kd));
  return {
    F, M: F * arm, CmEff: (4 * A) / (Math.PI * ka * ka), phaseLagDeg: (Math.atan2(dJ, dY) * 180) / Math.PI, ka, A, k,
    forcePerLength: (z) => ((2 * rho * g * H) / k) * A * (Math.cosh(k * (z + depth)) / Math.cosh(kd)),
  };
}

// DNV-RP-C205 Appendix D, Table D-1: rectangle of width 2a normal to the motion and height 2b along it,
// C_A referred to the area pi a² (values after Wendel).
const RECT_RATIO = [0.1, 0.2, 0.5, 1, 2, 5, 10];
const RECT_CA = [2.23, 1.98, 1.70, 1.51, 1.36, 1.21, 1.14];

/**
 * Added-mass coefficient of a rectangular section moving normal to its side of width 2a (height 2b along the
 * motion), referred to the area pi a². Tabulated values a/b = 0.1, 0.2, 0.5, 1, 2, 5, 10 → 2.23, 1.98, 1.70,
 * 1.51, 1.36, 1.21, 1.14 and 1.0 for a/b → ∞ (flat plate); intermediate ratios are interpolated linearly in
 * log10(a/b) (interpolation rule is this library's choice), between a/b = 10 and ∞ linearly in b/a.
 * Reference: DNV-RP-C205 (2019) Appendix D, Table D-1. Valid for a/b >= 0.1.
 * @param {number} aOverB width normal to motion / height along motion, -
 * @returns {number} C_A, -
 */
export function rectangleAddedMassCoefficient(aOverB) {
  if (!(aOverB >= 0.1)) throw new RangeError('aOverB must be >= 0.1');
  if (aOverB >= 10) return 1 + (1.14 - 1) * (10 / aOverB);
  const x = Math.log10(aOverB);
  for (let i = 1; i < RECT_RATIO.length; i++) {
    if (aOverB <= RECT_RATIO[i]) {
      const x0 = Math.log10(RECT_RATIO[i - 1]); const x1 = Math.log10(RECT_RATIO[i]);
      return RECT_CA[i - 1] + ((RECT_CA[i] - RECT_CA[i - 1]) * (x - x0)) / (x1 - x0);
    }
  }
  return 1.14;
}

/**
 * Added mass per unit length of standard two-dimensional sections, m_a = rho C_A A_R:
 *   'circle' (radius a): C_A = 1.0, A_R = pi a²;  'ellipse' (semi-axis a normal to motion): 1.0, pi a²;
 *   'plate' (flat plate of width 2a normal to motion): 1.0, pi a²;  'square' (side 2a): 1.51, pi a²;
 *   'rectangle' (width 2a normal to motion, height 2b): table of rectangleAddedMassCoefficient, pi a².
 * Reference: DNV-RP-C205 (2019) Appendix D, Table D-1. Infinite fluid, no free-surface or end effects.
 * @param {{shape:string, a:number, b?:number, rho?:number}} o half-width a normal to motion m, half-height b m
 * @returns {{Ca:number, referenceArea:number, addedMassPerLength:number}} -, m², kg/m
 */
export function addedMass2D({ shape, a, b, rho = 1025 }) {
  needPos(a, 'a');
  let Ca;
  if (shape === 'circle' || shape === 'ellipse' || shape === 'plate') Ca = 1.0;
  else if (shape === 'square') Ca = 1.51;
  else if (shape === 'rectangle') { needPos(b, 'b'); Ca = rectangleAddedMassCoefficient(a / b); } else throw new RangeError("shape must be 'circle', 'ellipse', 'plate', 'square' or 'rectangle'");
  const referenceArea = Math.PI * a * a;
  return { Ca, referenceArea, addedMassPerLength: rho * Ca * referenceArea };
}

/**
 * Added mass of standard three-dimensional bodies in infinite fluid:
 *   'sphere' (radius a): m_a = ½ rho (4/3) pi a³ (C_A = 0.5);
 *   'disc' (thin circular disc of radius a moving normal to its plane): m_a = (8/3) rho a³.
 * Reference: DNV-RP-C205 (2019) Appendix D, Table D-2; Lamb (1932) Hydrodynamics.
 * @param {{shape:string, a:number, rho?:number}} o radius m, density kg/m³
 * @returns {{addedMass:number, Ca:number, referenceVolume:number}} kg, coefficient referred to (4/3) pi a³, m³
 */
export function addedMass3D({ shape, a, rho = 1025 }) {
  needPos(a, 'a');
  const referenceVolume = (4 / 3) * Math.PI * a ** 3;
  if (shape === 'sphere') return { addedMass: 0.5 * rho * referenceVolume, Ca: 0.5, referenceVolume };
  if (shape === 'disc') return { addedMass: (8 / 3) * rho * a ** 3, Ca: 2 / Math.PI, referenceVolume };
  throw new RangeError("shape must be 'sphere' or 'disc'");
}

/**
 * Buoyancy of a straight cylindrical member cut by the still-water plane z = waterline.
 * The submerged volume is the section area times the submerged centreline length, which is exact when the
 * waterplane does not intersect an end face. A flooded member displaces only its steel (pass the wall thickness).
 * Reference: Archimedes' principle.
 * @param {object} o
 * @param {number[]} o.p1 end 1, m
 * @param {number[]} o.p2 end 2, m
 * @param {number} o.D outer diameter, m
 * @param {number} [o.t] wall thickness, m (required when flooded)
 * @param {boolean} [o.flooded=false]
 * @param {number} [o.waterline=0] still-water elevation, m
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{volume:number, force:number, centre:number[]|null, submergedLength:number}} m³, N (upward), m, m
 */
export function memberBuoyancy({ p1, p2, D, t, flooded = false, waterline = 0, rho = 1025, g = G }) {
  needVec3(p1, 'p1'); needVec3(p2, 'p2'); needPos(D, 'D');
  let A = (Math.PI * D * D) / 4;
  if (flooded) {
    if (!(t > 0 && t < D / 2)) throw new RangeError('t must satisfy 0 < t < D/2 for a flooded member');
    A -= (Math.PI * (D - 2 * t) ** 2) / 4;
  }
  const L = Math.hypot(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
  const z1 = p1[2] - waterline; const z2 = p2[2] - waterline;
  let s0 = 0; let s1 = 1;
  if (z1 > 0 && z2 > 0) s1 = 0;
  else if (z1 > 0) s0 = z1 / (z1 - z2);
  else if (z2 > 0) s1 = z1 / (z1 - z2);
  const Ls = (s1 - s0) * L;
  const sm = 0.5 * (s0 + s1);
  const volume = A * Ls;
  return {
    volume, force: rho * g * volume, submergedLength: Ls,
    centre: Ls > 0 ? [p1[0] + sm * (p2[0] - p1[0]), p1[1] + sm * (p2[1] - p1[1]), p1[2] + sm * (p2[2] - p1[2])] : null,
  };
}

/**
 * Wave slamming force per unit length on a horizontal cylinder, F_s = ½ rho C_s D v².
 * C_s is an input: pi is the von Kármán (momentum) value at first contact, 5.15 the Campbell & Weynberg (1980)
 * initial value; DNV-RP-C205 (2019) §8 requires C_s >= 3.0 for a smooth circular cylinder.
 * See impact.js `cylinderSlamCoefficient` for the variation with submergence.
 * @param {{D:number, v:number, Cs:number, rho?:number}} o diameter m, relative normal velocity m/s, -, kg/m³
 * @returns {number} force per unit length, N/m
 */
export function slamForcePerLength({ D, v, Cs, rho = 1025 }) {
  needPos(D, 'D'); needPos(Cs, 'Cs');
  return 0.5 * rho * Cs * D * v * v;
}

/**
 * Vortex-shedding frequency f_s = St U / D.
 * Reference: DNV-RP-C205 (2019) §9. St ≈ 0.2 for a circular cylinder in sub-critical flow (input).
 * @param {{U:number, D:number, St?:number}} o m/s, m, -
 * @returns {number} shedding frequency, Hz
 */
export function vortexSheddingFrequency({ U, D, St = 0.2 }) {
  needPos(D, 'D'); needPos(St, 'St');
  return (St * Math.abs(U)) / D;
}
