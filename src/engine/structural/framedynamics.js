/**
 * Linear dynamic response of frames: modal superposition and direct (Newmark / HHT-α) time history under nodal
 * force histories and three-component base acceleration, response-spectrum analysis with SRSS / CQC modal
 * combination and directional combination, design spectra (EN 1998-1, ISO 19901-2), spectrum-compatible synthetic
 * accelerograms, harmonic response / frequency response functions, dynamic amplification and a Campbell-diagram
 * helper for rotating-machinery excitation.
 * SI units; frequencies in Hz unless a name says omega (rad/s); accelerations in m/s² (not g).
 */

import { newmark } from '../math/ode.js';
import { makeRng } from '../math/rng.js';
import { assembleDense, memberForcesFromDisplacements, buildFrameSystem } from './frame3d.js';
import { frameMass, massInfluence } from './framemodal.js';

const DIRS = ['x', 'y', 'z'];

/* ------------------------------------------------------------------ single-degree-of-freedom kernels */

/**
 * Coefficients of the exact recurrence for a damped SDOF oscillator ü + 2ζω u̇ + ω² u = a(t) with the excitation
 * interpolated linearly over each time step (Nigam & Jennings 1969; Chopra, Dynamics of Structures, Table 5.2.1,
 * written for unit mass): u₊ = A u + B v + C a + D a₊, v₊ = A' u + B' v + C' a + D' a₊.
 * Validity: 0 ≤ ζ < 1, ω > 0.
 * @param {number} omega natural circular frequency [rad/s]
 * @param {number} zeta damping ratio [-]
 * @param {number} dt time step [s]
 * @returns {{A:number,B:number,C:number,D:number,A1:number,B1:number,C1:number,D1:number}}
 */
export function sdofRecurrence(omega, zeta, dt) {
  if (!(omega > 0)) throw new RangeError('omega must be positive');
  if (!(zeta >= 0 && zeta < 1)) throw new RangeError('zeta must lie in [0, 1)');
  if (!(dt > 0)) throw new RangeError('dt must be positive');
  const k = omega * omega;
  const r = Math.sqrt(1 - zeta * zeta);
  const wd = omega * r;
  const e = Math.exp(-zeta * omega * dt);
  const s = Math.sin(wd * dt); const c = Math.cos(wd * dt);
  const wt = omega * dt;
  return {
    A: e * (zeta / r * s + c),
    B: e * s / wd,
    C: (2 * zeta / wt + e * (((1 - 2 * zeta * zeta) / (wd * dt) - zeta / r) * s - (1 + 2 * zeta / wt) * c)) / k,
    D: (1 - 2 * zeta / wt + e * ((2 * zeta * zeta - 1) / (wd * dt) * s + 2 * zeta / wt * c)) / k,
    A1: -e * omega / r * s,
    B1: e * (c - zeta / r * s),
    C1: (-1 / dt + e * ((omega / r + zeta / (dt * r)) * s + c / dt)) / k,
    D1: (1 - e * (zeta / r * s + c)) / (k * dt),
  };
}

/**
 * Response history of a unit-mass SDOF oscillator to a sampled excitation a(t) (force per unit mass; for base
 * excitation pass a = −üg and read relative displacement), exact for piecewise-linear excitation.
 * @param {object} o
 * @param {number} o.omega [rad/s]
 * @param {number} o.zeta [-]
 * @param {number} o.dt [s]
 * @param {ArrayLike<number>} o.a excitation samples [m/s²]
 * @param {number} [o.u0=0] initial displacement [m]
 * @param {number} [o.v0=0] initial velocity [m/s]
 * @returns {{u:Float64Array, v:Float64Array, maxAbsU:number}} displacement [m], velocity [m/s]
 */
export function sdofResponse({ omega, zeta, dt, a, u0 = 0, v0 = 0 }) {
  if (!a || !(a.length >= 1)) throw new RangeError('a must be a non-empty array');
  const k = sdofRecurrence(omega, zeta, dt);
  const n = a.length;
  const u = new Float64Array(n); const v = new Float64Array(n);
  u[0] = u0; v[0] = v0;
  let mx = Math.abs(u0);
  for (let i = 0; i < n - 1; i++) {
    u[i + 1] = k.A * u[i] + k.B * v[i] + k.C * a[i] + k.D * a[i + 1];
    v[i + 1] = k.A1 * u[i] + k.B1 * v[i] + k.C1 * a[i] + k.D1 * a[i + 1];
    if (Math.abs(u[i + 1]) > mx) mx = Math.abs(u[i + 1]);
  }
  return { u, v, maxAbsU: mx };
}

/**
 * Elastic response spectrum of an accelerogram: Sd = max|u|, pseudo-velocity Sv = ω Sd, pseudo-acceleration
 * Sa = ω² Sd, by the exact piecewise-linear SDOF recurrence. The record is zero-padded by one period of free
 * vibration so that the peak after the end of the excitation is captured.
 * @param {object} o
 * @param {ArrayLike<number>} o.acc ground acceleration samples [m/s²]
 * @param {number} o.dt sampling interval [s]
 * @param {ArrayLike<number>} o.periods oscillator periods [s]
 * @param {number} [o.zeta=0.05] damping ratio [-]
 * @returns {{periods:Float64Array, Sd:Float64Array, Sv:Float64Array, Sa:Float64Array}} [s, m, m/s, m/s²]
 */
export function responseSpectrum({ acc, dt, periods, zeta = 0.05 }) {
  if (!acc || !(acc.length >= 2)) throw new RangeError('acc must hold at least two samples');
  if (!periods || !(periods.length >= 1)) throw new RangeError('periods must be a non-empty array');
  const np = periods.length;
  const Sd = new Float64Array(np); const Sv = new Float64Array(np); const Sa = new Float64Array(np);
  for (let p = 0; p < np; p++) {
    const T = periods[p];
    if (!(T > 0)) throw new RangeError('periods must be positive');
    const w = 2 * Math.PI / T;
    const k = sdofRecurrence(w, zeta, dt);
    let u = 0; let v = 0; let mx = 0;
    const n = acc.length;
    const tail = Math.ceil(T / dt);
    for (let i = 0; i < n - 1 + tail; i++) {
      const a0 = i < n ? -acc[i] : 0;
      const a1 = i + 1 < n ? -acc[i + 1] : 0;
      const un = k.A * u + k.B * v + k.C * a0 + k.D * a1;
      v = k.A1 * u + k.B1 * v + k.C1 * a0 + k.D1 * a1;
      u = un;
      if (Math.abs(u) > mx) mx = Math.abs(u);
    }
    Sd[p] = mx; Sv[p] = w * mx; Sa[p] = w * w * mx;
  }
  return { periods: Float64Array.from(periods), Sd, Sv, Sa };
}

/**
 * Steady-state dynamic amplification of a viscously damped SDOF system under harmonic force:
 * DAF = 1/√[(1 − r²)² + (2ζr)²], phase lag φ = atan2(2ζr, 1 − r²), force transmissibility
 * TR = √[1 + (2ζr)²]·DAF, with r = f/fn (Chopra §3.2, §3.5).
 * @param {object} o
 * @param {number} o.ratio frequency ratio r = forcing frequency / natural frequency [-]
 * @param {number} o.zeta damping ratio [-]
 * @returns {{daf:number, phase:number, transmissibility:number}} DAF [-], phase lag [rad], TR [-]
 */
export function dynamicAmplification({ ratio, zeta }) {
  if (!(ratio >= 0)) throw new RangeError('ratio must be non-negative');
  if (!(zeta >= 0)) throw new RangeError('zeta must be non-negative');
  const d = Math.hypot(1 - ratio * ratio, 2 * zeta * ratio);
  if (d === 0) return { daf: Infinity, phase: Math.PI / 2, transmissibility: Infinity };
  return { daf: 1 / d, phase: Math.atan2(2 * zeta * ratio, 1 - ratio * ratio), transmissibility: Math.hypot(1, 2 * zeta * ratio) / d };
}

/* ------------------------------------------------------------------ design spectra */

/** EN 1998-1:2004 Table 3.2 (Type 1) and Table 3.3 (Type 2): S, TB, TC, TD for ground types A–E. */
export const EC8_SPECTRUM_PARAMETERS = Object.freeze({
  1: Object.freeze({
    A: { S: 1.0, TB: 0.15, TC: 0.4, TD: 2.0 }, B: { S: 1.2, TB: 0.15, TC: 0.5, TD: 2.0 }, C: { S: 1.15, TB: 0.20, TC: 0.6, TD: 2.0 },
    D: { S: 1.35, TB: 0.20, TC: 0.8, TD: 2.0 }, E: { S: 1.4, TB: 0.15, TC: 0.5, TD: 2.0 },
  }),
  2: Object.freeze({
    A: { S: 1.0, TB: 0.05, TC: 0.25, TD: 1.2 }, B: { S: 1.35, TB: 0.05, TC: 0.25, TD: 1.2 }, C: { S: 1.5, TB: 0.10, TC: 0.25, TD: 1.2 },
    D: { S: 1.8, TB: 0.10, TC: 0.30, TD: 1.2 }, E: { S: 1.6, TB: 0.05, TC: 0.25, TD: 1.2 },
  }),
});

/**
 * Horizontal elastic response spectrum of EN 1998-1:2004 §3.2.2.2 (Eq. 3.2–3.5), or the design spectrum for
 * elastic analysis of §3.2.2.5 (Eq. 3.13–3.16) when a behaviour factor q is given.
 * Elastic: Se = ag S[1 + T/TB (2.5η − 1)] (0 ≤ T ≤ TB); 2.5 ag S η (TB–TC); 2.5 ag S η TC/T (TC–TD); 2.5 ag S η TC TD/T² (T ≥ TD);
 * η = √(10/(5 + ξ)) ≥ 0.55 with ξ in percent.
 * Design: Sd = ag S[2/3 + T/TB (2.5/q − 2/3)]; 2.5 ag S/q; max(2.5 ag S/q · TC/T, β ag); max(2.5 ag S/q · TC TD/T², β ag), β = 0.2.
 * @param {object} o
 * @param {number} o.ag design ground acceleration on type A ground [m/s²]
 * @param {'A'|'B'|'C'|'D'|'E'} [o.groundType='A']
 * @param {1|2} [o.type=1] spectrum type
 * @param {number} [o.zeta=0.05] viscous damping ratio [-] (elastic spectrum only)
 * @param {number} [o.q] behaviour factor [-]; when given the design spectrum is returned
 * @param {number} [o.beta=0.2] lower bound factor of the design spectrum [-]
 * @returns {{Sa:(T:number)=>number, S:number, TB:number, TC:number, TD:number, eta:number}} Sa(T) [m/s²] for T [s] ≥ 0
 */
export function ec8Spectrum({ ag, groundType = 'A', type = 1, zeta = 0.05, q, beta = 0.2 }) {
  if (!(ag > 0)) throw new RangeError('ag must be positive');
  const tbl = EC8_SPECTRUM_PARAMETERS[type];
  if (!tbl) throw new RangeError('type must be 1 or 2');
  const p = tbl[groundType];
  if (!p) throw new RangeError('groundType must be one of A, B, C, D, E');
  if (!(zeta >= 0)) throw new RangeError('zeta must be non-negative');
  const eta = Math.max(0.55, Math.sqrt(10 / (5 + 100 * zeta)));
  const { S, TB, TC, TD } = p;
  let Sa;
  if (q === undefined) {
    Sa = (T) => {
      if (!(T >= 0)) throw new RangeError('T must be non-negative');
      if (T <= TB) return ag * S * (1 + T / TB * (2.5 * eta - 1));
      if (T <= TC) return 2.5 * ag * S * eta;
      if (T <= TD) return 2.5 * ag * S * eta * TC / T;
      return 2.5 * ag * S * eta * TC * TD / (T * T);
    };
  } else {
    if (!(q > 0)) throw new RangeError('q must be positive');
    Sa = (T) => {
      if (!(T >= 0)) throw new RangeError('T must be non-negative');
      if (T <= TB) return ag * S * (2 / 3 + T / TB * (2.5 / q - 2 / 3));
      if (T <= TC) return 2.5 * ag * S / q;
      if (T <= TD) return Math.max(2.5 * ag * S / q * TC / T, beta * ag);
      return Math.max(2.5 * ag * S / q * TC * TD / (T * T), beta * ag);
    };
  }
  return { Sa, S, TB, TC, TD, eta };
}

/**
 * Site response spectrum shape of ISO 19901-2:2004 §7.1 (also API RP 2EQ) for 5 % damping, built from the mapped
 * rock spectral accelerations at 0.2 s and 1.0 s and the site coefficients Ca, Cv:
 *   Sa,site(T) = (3T + 0.4) Ca Sa,map(0.2)   for T ≤ 0.2 s
 *   Sa,site(T) = Cv Sa,map(1.0)/T ≤ Ca Sa,map(0.2)   for T > 0.2 s
 *   Sa,site(T) = 4 Cv Sa,map(1.0)/T²   for T > 4.0 s
 * For damping η ≠ 5 % the ordinates are divided by D = ln(100/η)/ln(20) (η in percent).
 * The site coefficients depend on the site class and on the mapped accelerations (ISO 19901-2 Tables 6 and 7) and
 * are supplied by the caller; Ca = Cv = 1.0 corresponds to site class A/B (rock).
 * @param {object} o
 * @param {number} o.Sa02 mapped spectral acceleration at 0.2 s [m/s²]
 * @param {number} o.Sa10 mapped spectral acceleration at 1.0 s [m/s²]
 * @param {number} [o.Ca=1] short-period site coefficient [-]
 * @param {number} [o.Cv=1] long-period site coefficient [-]
 * @param {number} [o.zeta=0.05] damping ratio [-]
 * @param {number} [o.scale=1] scale factor applied to all ordinates [-] (e.g. NALE or seismic reserve capacity factor)
 * @returns {{Sa:(T:number)=>number, plateau:number, cornerPeriod:number, dampingFactor:number}} Sa(T) [m/s²]; plateau value; period where
 *   the 1/T branch meets the plateau [s]
 */
export function isoSeismicSpectrum({ Sa02, Sa10, Ca = 1, Cv = 1, zeta = 0.05, scale = 1 }) {
  if (!(Sa02 > 0) || !(Sa10 > 0)) throw new RangeError('Sa02 and Sa10 must be positive');
  if (!(Ca > 0) || !(Cv > 0)) throw new RangeError('Ca and Cv must be positive');
  if (!(zeta > 0 && zeta < 1)) throw new RangeError('zeta must lie in (0, 1)');
  const D = Math.log(100 / (100 * zeta)) / Math.log(20);
  const plateau = Ca * Sa02;
  const k = scale / D;
  const Sa = (T) => {
    if (!(T >= 0)) throw new RangeError('T must be non-negative');
    let v;
    if (T <= 0.2) v = (3 * T + 0.4) * Ca * Sa02;
    else if (T <= 4) v = Math.min(Cv * Sa10 / T, plateau);
    else v = Math.min(4 * Cv * Sa10 / (T * T), plateau);
    return k * v;
  };
  return { Sa, plateau: k * plateau, cornerPeriod: Cv * Sa10 / plateau, dampingFactor: D };
}

/* ------------------------------------------------------------------ modal combination */

/**
 * CQC correlation coefficient between modes i and j (Der Kiureghian 1981; Wilson, Der Kiureghian & Bayo 1981):
 * ρ = 8 √(ζiζj) (ζi + r ζj) r^1.5 / [(1 − r²)² + 4 ζiζj r (1 + r²) + 4 (ζi² + ζj²) r²], r = ωj/ωi.
 * @param {number} wi circular frequency of mode i [rad/s]
 * @param {number} wj circular frequency of mode j [rad/s]
 * @param {number} zi damping ratio of mode i [-]
 * @param {number} [zj=zi] damping ratio of mode j [-]
 * @returns {number} ρij [-] (1 for i = j)
 */
export function cqcCoefficient(wi, wj, zi, zj = zi) {
  if (!(wi > 0) || !(wj > 0)) throw new RangeError('frequencies must be positive');
  if (!(zi >= 0) || !(zj >= 0)) throw new RangeError('damping ratios must be non-negative');
  if (wi === wj) return 1;
  if (zi === 0 && zj === 0) return 0;
  const r = wj / wi;
  return 8 * Math.sqrt(zi * zj) * (zi + r * zj) * r ** 1.5
    / ((1 - r * r) ** 2 + 4 * zi * zj * r * (1 + r * r) + 4 * (zi * zi + zj * zj) * r * r);
}

/**
 * Combine modal peak responses: SRSS √(Σ rᵢ²) or CQC √(Σᵢ Σⱼ ρᵢⱼ rᵢ rⱼ).
 * @param {object} o
 * @param {ArrayLike<number>} o.peaks signed modal peak values of one response quantity
 * @param {ArrayLike<number>} [o.omegas] modal circular frequencies [rad/s] (CQC)
 * @param {number|ArrayLike<number>} [o.zeta=0.05] modal damping ratio(s) [-] (CQC)
 * @param {'SRSS'|'CQC'} [o.method='CQC']
 * @returns {number} combined peak (non-negative)
 */
export function combineModal({ peaks, omegas, zeta = 0.05, method = 'CQC' }) {
  const n = peaks.length;
  let s = 0;
  if (method === 'SRSS') {
    for (let i = 0; i < n; i++) s += peaks[i] * peaks[i];
    return Math.sqrt(s);
  }
  if (method !== 'CQC') throw new RangeError('method must be "SRSS" or "CQC"');
  if (!omegas || omegas.length !== n) throw new RangeError('omegas must match peaks for CQC');
  const z = (i) => (typeof zeta === 'number' ? zeta : zeta[i]);
  for (let i = 0; i < n; i++) {
    s += peaks[i] * peaks[i];
    for (let j = i + 1; j < n; j++) s += 2 * cqcCoefficient(omegas[i], omegas[j], z(i), z(j)) * peaks[i] * peaks[j];
  }
  return Math.sqrt(Math.max(s, 0));
}

/**
 * Directional combination of the responses to the individual earthquake components: SRSS, or the
 * 100 % / 30 % / 30 % rule max over k of (R_k + 0.3 Σ_{j≠k} R_j) (EN 1998-1 §4.3.3.5.1; ISO 19901-2 §8.2.3 uses
 * 100/100/50 or SRSS — pass `secondary` to change the 0.3).
 * @param {ArrayLike<number>} r component responses (non-negative)
 * @param {'SRSS'|'100-30-30'} [method='SRSS']
 * @param {number} [secondary=0.3] factor on the non-governing components
 * @returns {number}
 */
export function combineDirections(r, method = 'SRSS', secondary = 0.3) {
  let s = 0;
  if (method === 'SRSS') { for (const v of r) s += v * v; return Math.sqrt(s); }
  if (method !== '100-30-30') throw new RangeError('method must be "SRSS" or "100-30-30"');
  let tot = 0;
  for (const v of r) tot += Math.abs(v);
  let best = 0;
  for (const v of r) best = Math.max(best, Math.abs(v) + secondary * (tot - Math.abs(v)));
  return best;
}

/* ------------------------------------------------------------------ response spectrum analysis */

function spectrumFn(spectrum, d) {
  const s = typeof spectrum === 'function' ? spectrum : spectrum && spectrum[DIRS[d]];
  return typeof s === 'function' ? s : null;
}

/**
 * Response-spectrum analysis of a frame from its modal results. For each excitation direction d and mode n the
 * peak modal response is u_n = φ_n Γ_n,d Sa(T_n)/ω_n² (Chopra §13.7); nodal displacements, member end forces and
 * base shears are combined over the modes by CQC or SRSS and over the directions by SRSS or the 100/30/30 rule.
 * @param {object} modal result of modalAnalysis
 * @param {object} o
 * @param {((T:number)=>number)|{x?:Function, y?:Function, z?:Function}} o.spectrum pseudo-acceleration spectrum Sa(T) [m/s²] for all
 *   listed directions, or one function per direction
 * @param {Array<'x'|'y'|'z'>} [o.directions=['x','y']] excitation directions
 * @param {number|number[]} [o.zeta=0.05] modal damping ratio(s) used in the CQC coefficients [-]
 * @param {'CQC'|'SRSS'} [o.combination='CQC']
 * @param {'SRSS'|'100-30-30'} [o.directional='SRSS']
 * @returns {{displacements:Float64Array, memberForces:Float64Array[], baseShear:{x:number,y:number,z:number}, perDirection:Object<string,object>,
 *   modal:Array<{mode:number, period:number, Sa:Object<string,number>}>, massParticipation:{x:number,y:number,z:number}}}
 *   combined peak displacements over all dofs [m, rad], member end forces (12 per member, local, non-negative) [N, N·m], base shear
 *   per global direction [N], per-direction results, and the cumulative effective mass fraction of the modes used
 */
export function responseSpectrumAnalysis(modal, { spectrum, directions = ['x', 'y'], zeta = 0.05, combination = 'CQC', directional = 'SRSS' } = {}) {
  if (!modal || !modal.system || !modal.modes) throw new RangeError('modal must be a result of modalAnalysis');
  const sys = modal.system;
  const nm = modal.modes.length;
  const omegas = modal.omegas;
  const modeForces = modal.modes.map((phi) => memberForcesFromDisplacements(sys, phi));
  const perDirection = {};
  const modalInfo = modal.periods.map((T, k) => ({ mode: k + 1, period: T, Sa: {} }));
  const peaks = new Float64Array(nm);
  const comb = (get) => { for (let k = 0; k < nm; k++) peaks[k] = get(k); return combineModal({ peaks, omegas, zeta, method: combination }); };
  for (const dn of directions) {
    const d = DIRS.indexOf(dn);
    if (d < 0) throw new RangeError(`direction "${dn}" is not x, y or z`);
    const fn = spectrumFn(spectrum, d);
    if (!fn) throw new RangeError(`spectrum for direction ${dn} is missing`);
    // modal amplitude: q_n = Γ_n Sa/ω²
    const amp = new Float64Array(nm);
    for (let k = 0; k < nm; k++) {
      const Sa = fn(modal.periods[k]);
      if (!(Sa >= 0)) throw new RangeError('spectrum must return non-negative accelerations');
      modalInfo[k].Sa[dn] = Sa;
      amp[k] = modal.participation[k][dn] * Sa / (omegas[k] * omegas[k]);
    }
    const disp = new Float64Array(sys.nDof);
    for (let g = 0; g < sys.nDof; g++) disp[g] = comb((k) => amp[k] * modal.modes[k][g]);
    const mf = sys.elems.map((el, e) => {
      const f = new Float64Array(12);
      for (let i = 0; i < 12; i++) f[i] = comb((k) => amp[k] * modeForces[k][e][i]);
      return f;
    });
    // base shear in global direction e: V_n = Γ_n,e · (Γ_n,d Sa)
    const baseShear = {};
    for (const en of DIRS) baseShear[en] = comb((k) => modal.participation[k][en] * amp[k] * omegas[k] * omegas[k]);
    perDirection[dn] = { displacements: disp, memberForces: mf, baseShear };
  }
  const displacements = new Float64Array(sys.nDof);
  const vals = new Array(directions.length);
  const dirComb = (get) => { for (let i = 0; i < directions.length; i++) vals[i] = get(perDirection[directions[i]]); return combineDirections(vals, directional); };
  for (let g = 0; g < sys.nDof; g++) displacements[g] = dirComb((p) => p.displacements[g]);
  const memberForces = sys.elems.map((el, e) => {
    const f = new Float64Array(12);
    for (let i = 0; i < 12; i++) f[i] = dirComb((p) => p.memberForces[e][i]);
    return f;
  });
  const baseShear = {}; const massParticipation = {};
  for (const en of DIRS) {
    baseShear[en] = dirComb((p) => p.baseShear[en]);
    massParticipation[en] = modal.massFraction.reduce((s, m) => s + m[en], 0);
  }
  return { displacements, memberForces, baseShear, perDirection, modal: modalInfo, massParticipation };
}

/* ------------------------------------------------------------------ time history */

function sampler(h, dt, nSteps) {
  if (typeof h === 'function') {
    const out = new Float64Array(nSteps + 1);
    for (let i = 0; i <= nSteps; i++) out[i] = h(i * dt);
    return out;
  }
  if (!h || typeof h.length !== 'number') throw new RangeError('history must be an array of samples at dt or a function of time');
  const out = new Float64Array(nSteps + 1);
  for (let i = 0; i <= nSteps && i < h.length; i++) out[i] = h[i];
  return out;
}

const DOF_INDEX = { ux: 0, uy: 1, uz: 2, rx: 3, ry: 4, rz: 5 };
function dofOf(sys, o) {
  const k = typeof o.dof === 'string' ? DOF_INDEX[o.dof] : o.dof;
  if (!Number.isInteger(o.node) || o.node < 0 || o.node >= sys.nNode || !(k >= 0 && k < 6)) throw new RangeError('output / load needs a valid node and dof (0–5 or ux…rz)');
  return 6 * o.node + k;
}

/**
 * Modal-superposition time history. Each modal equation q̈ + 2ζω q̇ + ω² q = φᵀF(t) − Γ·üg(t) (mass-normalised
 * modes) is integrated with the exact piecewise-linear recurrence, so the only approximations are the modal
 * truncation and the linear interpolation of the excitation between samples.
 * @param {object} modal result of modalAnalysis
 * @param {object} o
 * @param {number} o.dt time step [s]
 * @param {number} o.nSteps number of steps
 * @param {number|number[]} [o.zeta=0.05] modal damping ratio(s) [-]
 * @param {Array<{node:number, dof:(number|string), history:(ArrayLike<number>|Function)}>} [o.nodalForces] force histories [N or N·m]
 *   sampled at dt or as functions of time
 * @param {{x?:(ArrayLike<number>|Function), y?:(ArrayLike<number>|Function), z?:(ArrayLike<number>|Function)}} [o.baseAcc] ground
 *   acceleration components [m/s²]
 * @param {Array<{node:number, dof:(number|string)}>} [o.outputs] dofs whose (relative) displacement histories are returned
 * @returns {{t:Float64Array, q:Float64Array[], outputs:Float64Array[], peak:Float64Array, peakTime:Float64Array, displacementsAt:(step:number)=>Float64Array}}
 *   modal coordinates q[n][step], output histories [m, rad], peak |displacement| over all dofs and the time it occurs
 */
export function modalTimeHistory(modal, { dt, nSteps, zeta = 0.05, nodalForces = [], baseAcc = {}, outputs = [] } = {}) {
  if (!modal || !modal.system) throw new RangeError('modal must be a result of modalAnalysis');
  if (!(dt > 0)) throw new RangeError('dt must be positive');
  if (!Number.isInteger(nSteps) || nSteps < 1) throw new RangeError('nSteps must be a positive integer');
  const sys = modal.system;
  const nm = modal.modes.length;
  const forces = nodalForces.map((f) => ({ g: dofOf(sys, f), h: sampler(f.history, dt, nSteps) }));
  const acc = DIRS.map((d) => (baseAcc[d] ? sampler(baseAcc[d], dt, nSteps) : null));
  const t = new Float64Array(nSteps + 1);
  for (let i = 0; i <= nSteps; i++) t[i] = i * dt;
  const q = [];
  const p = new Float64Array(nSteps + 1);
  for (let k = 0; k < nm; k++) {
    p.fill(0);
    const phi = modal.modes[k];
    for (const f of forces) { const c = phi[f.g]; if (c !== 0) for (let i = 0; i <= nSteps; i++) p[i] += c * f.h[i]; }
    for (let d = 0; d < 3; d++) {
      if (!acc[d]) continue;
      const g = modal.participation[k][DIRS[d]];
      for (let i = 0; i <= nSteps; i++) p[i] -= g * acc[d][i];
    }
    const z = typeof zeta === 'number' ? zeta : zeta[k];
    q.push(sdofResponse({ omega: modal.omegas[k], zeta: z, dt, a: p }).u);
  }
  const outIdx = outputs.map((o) => dofOf(sys, o));
  const out = outIdx.map((g) => {
    const h = new Float64Array(nSteps + 1);
    for (let k = 0; k < nm; k++) { const c = modal.modes[k][g]; if (c !== 0) for (let i = 0; i <= nSteps; i++) h[i] += c * q[k][i]; }
    return h;
  });
  const displacementsAt = (step) => {
    const u = new Float64Array(sys.nDof);
    for (let k = 0; k < nm; k++) { const c = q[k][step]; const phi = modal.modes[k]; for (let g = 0; g < sys.nDof; g++) u[g] += c * phi[g]; }
    return u;
  };
  const peak = new Float64Array(sys.nDof); const peakTime = new Float64Array(sys.nDof);
  const u = new Float64Array(sys.nDof);
  for (let i = 0; i <= nSteps; i++) {
    u.fill(0);
    for (let k = 0; k < nm; k++) { const c = q[k][i]; if (c === 0) continue; const phi = modal.modes[k]; for (let g = 0; g < sys.nDof; g++) u[g] += c * phi[g]; }
    for (let g = 0; g < sys.nDof; g++) if (Math.abs(u[g]) > peak[g]) { peak[g] = Math.abs(u[g]); peakTime[g] = t[i]; }
  }
  return { t, q, outputs: out, peak, peakTime, displacementsAt };
}

/**
 * Direct time integration of M ü + C u̇ + K u = F(t) − M r üg(t) by Newmark-β (average acceleration by default) or
 * HHT-α, with Rayleigh damping C = α M + β K. Dense matrices are used, so the method is intended for models up to
 * roughly a thousand free dofs; use modal superposition for larger models.
 * @param {object} model frame model
 * @param {object} o
 * @param {number} o.dt time step [s]
 * @param {number} o.nSteps number of steps
 * @param {{alpha:number, beta:number}} [o.rayleigh] Rayleigh coefficients [1/s, s] (default undamped)
 * @param {number} [o.alphaHHT=0] HHT-α parameter in [−1/3, 0]
 * @param {Array<{node:number, dof:(number|string), history:(ArrayLike<number>|Function)}>} [o.nodalForces]
 * @param {{x?:any, y?:any, z?:any}} [o.baseAcc] ground acceleration components [m/s²]
 * @param {Array<{node:number, dof:(number|string)}>} [o.outputs] dofs to record
 * @param {object} [o.massOptions] options of framemodal.frameMass (lumped, waterline, Ca, …)
 * @param {object} [o.system] reuse a system from buildFrameSystem
 * @returns {{t:Float64Array, outputs:Float64Array[], velocities:Float64Array[], accelerations:Float64Array[], peak:Float64Array}} histories
 *   of the output dofs (relative to the base) and peak |displacement| over all dofs
 */
export function directTimeHistory(model, {
  dt, nSteps, rayleigh, alphaHHT = 0, nodalForces = [], baseAcc = {}, outputs = [], massOptions = {}, system,
} = {}) {
  if (!(dt > 0)) throw new RangeError('dt must be positive');
  if (!Number.isInteger(nSteps) || nSteps < 1) throw new RangeError('nSteps must be a positive integer');
  const sys = system || buildFrameSystem(model);
  const n = sys.nFree;
  if (n > 1500) throw new RangeError('directTimeHistory uses dense matrices and is limited to 1500 free dofs; use modalTimeHistory');
  const mass = frameMass(sys, massOptions);
  const K = assembleDense(sys, sys.elems.map((e) => e.ke), true);
  const M = assembleDense(sys, mass.mats, false);
  for (let q = 0; q < n; q++) M[q][q] += mass.nodeDiag[sys.eqDof[q]];
  let C;
  if (rayleigh) {
    C = M.map((row, i) => row.map((v, j) => rayleigh.alpha * v + rayleigh.beta * K[i][j]));
  }
  const infl = massInfluence(sys, mass);
  const forces = nodalForces.map((f) => ({ q: sys.free[dofOf(sys, f)], h: sampler(f.history, dt, nSteps) }));
  const acc = DIRS.map((d) => (baseAcc[d] ? sampler(baseAcc[d], dt, nSteps) : null));
  const fvec = new Float64Array(n);
  const f = (time) => {
    // samples are addressed by step; newmark calls f at multiples of dt
    const i = Math.min(nSteps, Math.max(0, Math.round(time / dt)));
    fvec.fill(0);
    for (const ff of forces) if (ff.q >= 0) fvec[ff.q] += ff.h[i];
    for (let d = 0; d < 3; d++) {
      if (!acc[d]) continue;
      const a = acc[d][i];
      if (a === 0) continue;
      const mr = infl.Mr[d];
      for (let q = 0; q < n; q++) fvec[q] -= mr[q] * a;
    }
    return fvec;
  };
  const args = { M, K, f, u0: new Float64Array(n), v0: new Float64Array(n), dt, nSteps, alphaHHT };
  if (C) args.C = C;
  const res = newmark(args);
  const outQ = outputs.map((o) => sys.free[dofOf(sys, o)]);
  const pick = (arr) => outQ.map((q) => {
    const h = new Float64Array(nSteps + 1);
    if (q >= 0) for (let i = 0; i <= nSteps; i++) h[i] = arr[i][q];
    return h;
  });
  const peak = new Float64Array(sys.nDof);
  for (let i = 0; i <= nSteps; i++) {
    const u = res.u[i];
    for (let q = 0; q < n; q++) { const a = Math.abs(u[q]); if (a > peak[sys.eqDof[q]]) peak[sys.eqDof[q]] = a; }
  }
  return { t: res.t, outputs: pick(res.u), velocities: pick(res.v), accelerations: pick(res.a), peak };
}

/* ------------------------------------------------------------------ harmonic response */

/**
 * Steady-state harmonic response by modal superposition: U(ω) = Σ_n φ_n (φ_nᵀ F)/(ω_n² − ω² + 2 i ζ_n ω_n ω) for a set
 * of in-phase nodal force amplitudes. With a single unit force the result is the receptance frequency response
 * function between the loaded dof and each output dof.
 * @param {object} modal result of modalAnalysis
 * @param {object} o
 * @param {ArrayLike<number>} o.frequencies forcing frequencies [Hz]
 * @param {Array<{node:number, dof:(number|string), amplitude:number}>} o.forces force amplitudes [N or N·m]
 * @param {Array<{node:number, dof:(number|string)}>} o.outputs response dofs
 * @param {number|number[]} [o.zeta=0.02] modal damping ratio(s) [-]
 * @returns {{frequencies:Float64Array, amplitude:Float64Array[], phase:Float64Array[], re:Float64Array[], im:Float64Array[]}} per output:
 *   amplitude [m or rad], phase lag [rad] (positive = response lags the force), real and imaginary parts
 */
export function harmonicResponse(modal, { frequencies, forces, outputs, zeta = 0.02 } = {}) {
  if (!modal || !modal.system) throw new RangeError('modal must be a result of modalAnalysis');
  if (!frequencies || !(frequencies.length >= 1)) throw new RangeError('frequencies must be a non-empty array');
  if (!Array.isArray(forces) || forces.length === 0) throw new RangeError('forces must be a non-empty array');
  if (!Array.isArray(outputs) || outputs.length === 0) throw new RangeError('outputs must be a non-empty array');
  const sys = modal.system;
  const nm = modal.modes.length;
  const nf = frequencies.length;
  const pf = new Float64Array(nm);
  for (const f of forces) { const g = dofOf(sys, f); for (let k = 0; k < nm; k++) pf[k] += modal.modes[k][g] * f.amplitude; }
  const outIdx = outputs.map((o) => dofOf(sys, o));
  const re = outIdx.map(() => new Float64Array(nf)); const im = outIdx.map(() => new Float64Array(nf));
  for (let i = 0; i < nf; i++) {
    const w = 2 * Math.PI * frequencies[i];
    if (!(w >= 0)) throw new RangeError('frequencies must be non-negative');
    for (let k = 0; k < nm; k++) {
      const wn = modal.omegas[k];
      const z = typeof zeta === 'number' ? zeta : zeta[k];
      const a = wn * wn - w * w; const b = 2 * z * wn * w;
      const den = a * a + b * b;
      const hr = pf[k] * a / den; const hi = -pf[k] * b / den;
      for (let o = 0; o < outIdx.length; o++) { const c = modal.modes[k][outIdx[o]]; re[o][i] += c * hr; im[o][i] += c * hi; }
    }
  }
  const amplitude = re.map((r, o) => r.map((v, i) => Math.hypot(v, im[o][i])));
  const phase = re.map((r, o) => r.map((v, i) => Math.atan2(-im[o][i], v)));
  return { frequencies: Float64Array.from(frequencies), amplitude, phase, re, im };
}

/* ------------------------------------------------------------------ synthetic accelerogram */

/**
 * Spectrum-compatible synthetic accelerogram: a sum of sinusoids with seeded random phases, shaped by a
 * rise – plateau – exponential-decay intensity envelope (Jennings, Housner & Tsai 1968), whose amplitudes are
 * corrected iteratively by the ratio target / computed response spectrum at each component frequency
 * (Gasparini & Vanmarcke 1976, SIMQKE procedure). The mean is removed and the record tapers to zero at both ends.
 * @param {object} o
 * @param {(T:number)=>number} o.spectrum target pseudo-acceleration spectrum Sa(T) [m/s²]
 * @param {number} o.duration record length [s]
 * @param {number} [o.dt=0.01] time step [s]
 * @param {number} [o.zeta=0.05] damping ratio of the target spectrum [-]
 * @param {number} [o.fMin=0.2] lowest matched frequency [Hz]
 * @param {number} [o.fMax=25] highest matched frequency [Hz] (must be below 1/(2 dt))
 * @param {number} [o.nFreq=120] number of sinusoids (log-spaced)
 * @param {number} [o.iterations=10] spectral-matching iterations
 * @param {number} [o.seed=1] random seed
 * @param {{rise?:number, plateauEnd?:number, decay?:number}} [o.envelope] rise time and end of the strong phase as fractions of the
 *   duration (defaults 0.15, 0.6) and the decay exponent c [1/s] of exp(−c (t − t2)) (default 3/(duration − t2))
 * @returns {{t:Float64Array, acc:Float64Array, periods:Float64Array, target:Float64Array, achieved:Float64Array, maxError:number,
 *   meanError:number, pga:number}} time [s], acceleration [m/s²], matched periods, target and achieved Sa, relative errors
 */
export function spectrumCompatibleAccelerogram({
  spectrum, duration, dt = 0.01, zeta = 0.05, fMin = 0.2, fMax = 25, nFreq = 120, iterations = 10, seed = 1, envelope = {},
}) {
  if (typeof spectrum !== 'function') throw new RangeError('spectrum must be a function of period');
  if (!(duration > 0) || !(dt > 0)) throw new RangeError('duration and dt must be positive');
  if (!(fMin > 0 && fMax > fMin)) throw new RangeError('fMin and fMax must satisfy 0 < fMin < fMax');
  if (fMax >= 0.5 / dt) throw new RangeError('fMax must be below the Nyquist frequency 1/(2 dt)');
  const n = Math.floor(duration / dt) + 1;
  const rng = makeRng(seed);
  const t1 = (envelope.rise === undefined ? 0.15 : envelope.rise) * duration;
  const t2 = (envelope.plateauEnd === undefined ? 0.6 : envelope.plateauEnd) * duration;
  const c = envelope.decay === undefined ? 3 / Math.max(duration - t2, dt) : envelope.decay;
  const env = new Float64Array(n);
  const t = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const ti = i * dt;
    t[i] = ti;
    env[i] = ti < t1 ? (ti / t1) ** 2 : (ti <= t2 ? 1 : Math.exp(-c * (ti - t2)));
  }
  env[n - 1] = 0;
  const freqs = new Float64Array(nFreq); const periods = new Float64Array(nFreq); const target = new Float64Array(nFreq);
  const phases = new Float64Array(nFreq); const amp = new Float64Array(nFreq);
  for (let k = 0; k < nFreq; k++) {
    freqs[k] = fMin * (fMax / fMin) ** (nFreq === 1 ? 0 : k / (nFreq - 1));
    periods[k] = 1 / freqs[k];
    target[k] = spectrum(periods[k]);
    if (!(target[k] > 0)) throw new RangeError('spectrum must be positive over the matched period range');
    phases[k] = 2 * Math.PI * rng.next();
    amp[k] = target[k] / Math.sqrt(nFreq);
  }
  const acc = new Float64Array(n);
  const synth = () => {
    acc.fill(0);
    for (let k = 0; k < nFreq; k++) {
      const w = 2 * Math.PI * freqs[k]; const a = amp[k]; const ph = phases[k];
      for (let i = 0; i < n; i++) acc[i] += a * Math.sin(w * t[i] + ph);
    }
    let mean = 0;
    for (let i = 0; i < n; i++) { acc[i] *= env[i]; mean += acc[i]; }
    mean /= n;
    for (let i = 0; i < n; i++) acc[i] -= mean * env[i];
  };
  let achieved = null;
  for (let it = 0; it <= iterations; it++) {
    synth();
    achieved = responseSpectrum({ acc, dt, periods, zeta }).Sa;
    if (it === iterations) break;
    for (let k = 0; k < nFreq; k++) amp[k] *= target[k] / achieved[k];
  }
  let maxError = 0; let meanError = 0; let pga = 0;
  for (let k = 0; k < nFreq; k++) { const e = Math.abs(achieved[k] / target[k] - 1); maxError = Math.max(maxError, e); meanError += e / nFreq; }
  for (let i = 0; i < n; i++) pga = Math.max(pga, Math.abs(acc[i]));
  return { t, acc, periods, target, achieved, maxError, meanError, pga };
}

/* ------------------------------------------------------------------ Campbell diagram */

/**
 * Campbell (interference) diagram data for a structure supporting rotating machinery: excitation lines
 * f = order × speed/60 against natural frequencies with a ±margin exclusion band, the critical speeds where they
 * cross, and the interferences inside the operating speed range (separation-margin practice of API 610 / API 617
 * and DNV rules, typically ±10 % to ±20 %).
 * @param {object} o
 * @param {ArrayLike<number>} o.naturalFrequencies natural frequencies [Hz]
 * @param {number[]} [o.orders=[1, 2]] excitation orders per revolution (1× unbalance, 2× misalignment, blade-pass = number of blades, …)
 * @param {number} o.speedMin lower end of the operating speed range [rpm]
 * @param {number} o.speedMax upper end of the operating speed range [rpm]
 * @param {number} [o.margin=0.1] separation margin as a fraction of the natural frequency [-]
 * @param {number} [o.nPoints=50] points per excitation line
 * @returns {{speeds:Float64Array, lines:Array<{order:number, f:Float64Array}>, bands:Array<{mode:number, f:number, low:number, high:number}>,
 *   criticalSpeeds:Array<{mode:number, order:number, speed:number, inRange:boolean}>, interferences:Array<{mode:number, order:number,
 *   separation:number}>, ok:boolean}} speeds [rpm], frequencies [Hz]; separation = min |f_exc − f_n|/f_n over the operating range
 */
export function campbellDiagram({ naturalFrequencies, orders = [1, 2], speedMin, speedMax, margin = 0.1, nPoints = 50 }) {
  if (!naturalFrequencies || !(naturalFrequencies.length >= 1)) throw new RangeError('naturalFrequencies must be a non-empty array');
  if (!(speedMin >= 0) || !(speedMax > speedMin)) throw new RangeError('speedMin and speedMax must satisfy 0 <= speedMin < speedMax');
  if (!(margin >= 0 && margin < 1)) throw new RangeError('margin must lie in [0, 1)');
  const top = speedMax * 1.2;
  const speeds = new Float64Array(nPoints);
  for (let i = 0; i < nPoints; i++) speeds[i] = top * i / (nPoints - 1);
  const lines = orders.map((order) => {
    if (!(order > 0)) throw new RangeError('orders must be positive');
    return { order, f: speeds.map((s) => order * s / 60) };
  });
  const bands = []; const criticalSpeeds = []; const interferences = [];
  for (let k = 0; k < naturalFrequencies.length; k++) {
    const fn = naturalFrequencies[k];
    if (!(fn > 0)) throw new RangeError('naturalFrequencies must be positive');
    bands.push({ mode: k + 1, f: fn, low: fn * (1 - margin), high: fn * (1 + margin) });
    for (const order of orders) {
      const speed = 60 * fn / order;
      criticalSpeeds.push({ mode: k + 1, order, speed, inRange: speed >= speedMin && speed <= speedMax });
      const fLo = order * speedMin / 60; const fHi = order * speedMax / 60;
      const sep = fn < fLo ? (fLo - fn) / fn : (fn > fHi ? (fn - fHi) / fn : 0);
      if (sep < margin) interferences.push({ mode: k + 1, order, separation: sep });
    }
  }
  return { speeds, lines, bands, criticalSpeeds, interferences, ok: interferences.length === 0 };
}
