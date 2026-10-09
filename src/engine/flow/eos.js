// Cubic equations of state (Peng–Robinson, Soave–Redlich–Kwong) for mixtures with van der Waals one-fluid mixing
// rules: volumetric and caloric properties, PT flash with stability analysis, saturation points, phase envelope,
// and transport-property estimators. SI units: T in K, P absolute in Pa, molar volume m³/mol, molar energies
// J/mol, densities kg/m³, viscosity Pa·s.

import { R, getComponent, idealGasCp, idealGasEnthalpyEntropy, defaultKijMatrix } from './components.js';
import { solveDense } from '../math/linalg.js';
import { waterVapourPressure } from './water.js';
import { gasViscosityLee } from './blackoil.js';

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}

// Generic two-parameter cubic P = RT/(v − b) − a/((v + δ1 b)(v + δ2 b)).
const EOS_DEF = {
  PR: { omegaA: 0.45724, omegaB: 0.07780, d1: 1 + Math.SQRT2, d2: 1 - Math.SQRT2, vcOverB: 0.307401 / 0.077796 },
  SRK: { omegaA: 0.42748, omegaB: 0.08664, d1: 1, d2: 0, vcOverB: (1 / 3) / 0.08664 },
};

function mOfOmega(eos, w) {
  if (eos === 'SRK') return 0.480 + 1.574 * w - 0.176 * w * w;
  // Peng–Robinson (1976); the 1978 expression is used for heavier components (ω > 0.49)
  if (w <= 0.49) return 0.37464 + 1.54226 * w - 0.26992 * w * w;
  return 0.379642 + 1.48503 * w - 0.164423 * w * w + 0.016666 * w ** 3;
}

/**
 * Build an equation-of-state mixture model.
 *
 * Peng–Robinson: a = 0.45724 R²Tc²/Pc · α, b = 0.07780 R Tc/Pc, α = [1 + m(1 − √Tr)]²,
 *   m = 0.37464 + 1.54226 ω − 0.26992 ω² (ω ≤ 0.49), m = 0.379642 + 1.48503 ω − 0.164423 ω² + 0.016666 ω³ (ω > 0.49).
 * Soave–Redlich–Kwong: a = 0.42748 R²Tc²/Pc · α, b = 0.08664 R Tc/Pc, m = 0.480 + 1.574 ω − 0.176 ω².
 * Mixing: a = ΣΣ x_i x_j √(a_i a_j)(1 − k_ij), b = Σ x_i b_i.
 * Volume translation (Péneloux et al. 1982): v = v_eos − Σ x_i c_i. With `volumeShift: true` the PR shifts are the
 * tabulated ratios s_i = c_i/b_i of the component records (Jhaveri & Youngren 1988; 0 where none is tabulated) and
 * the SRK shifts are c_i = 0.40768 (R Tc/Pc)(0.29441 − Z_RA) with Z_RA = 0.29056 − 0.08775 ω (Yamada–Gunn).
 * An array gives the ratios s_i explicitly. The shift changes densities only, not phase equilibrium.
 * References: Peng & Robinson, Ind. Eng. Chem. Fundam. 15 (1976) 59; Robinson & Peng, GPA RR-28 (1978);
 * Soave, Chem. Eng. Sci. 27 (1972) 1197; Péneloux, Rauzy & Fréze, Fluid Phase Equilib. 8 (1982) 7.
 * Validity: non-polar and weakly polar mixtures; aqueous, methanol and glycol phases are poorly represented by
 * these mixing rules.
 * @param {(string|object)[]} components database ids or component records
 * @param {object} [o]
 * @param {'PR'|'SRK'} [o.eos='PR']
 * @param {number[][]} [o.kij] symmetric binary interaction matrix (default: `defaultKijMatrix`)
 * @param {boolean|number[]} [o.volumeShift=false]
 * @returns {object} mixture model { n, eos, comps, ids, M[], Tc[], Pc[], omega[], m[], ac[], b[], c[], kij }
 *   (kg/mol, K, Pa, –, –, Pa·m⁶/mol², m³/mol, m³/mol, –)
 */
export function makeMixture(components, { eos = 'PR', kij, volumeShift = false } = {}) {
  const def = EOS_DEF[eos];
  if (!def) throw new RangeError(`eos = ${eos} is not one of PR, SRK`);
  if (!Array.isArray(components) || components.length === 0) throw new RangeError('components must be a non-empty array');
  const comps = components.map((c) => (typeof c === 'string' ? getComponent(c) : c));
  const n = comps.length;
  const k = kij ?? defaultKijMatrix(comps);
  if (k.length !== n || k.some((row) => row.length !== n)) throw new RangeError('kij must be an n×n matrix');
  const mix = {
    n, eos, def, comps, ids: comps.map((c) => c.id), kij: k,
    M: comps.map((c) => c.M), Tc: comps.map((c) => c.Tc), Pc: comps.map((c) => c.Pc),
    omega: comps.map((c) => c.omega),
    m: comps.map((c) => mOfOmega(eos, c.omega)),
    ac: comps.map((c) => def.omegaA * R * R * c.Tc * c.Tc / c.Pc),
    b: comps.map((c) => def.omegaB * R * c.Tc / c.Pc),
  };
  mix.sqrtAc = mix.ac.map(Math.sqrt);
  if (Array.isArray(volumeShift)) {
    if (volumeShift.length !== n) throw new RangeError('volumeShift must have one entry per component');
    mix.c = volumeShift.map((s, i) => s * mix.b[i]);
  } else if (volumeShift === true) {
    mix.c = comps.map((c, i) => (eos === 'PR'
      ? (c.shiftPR ?? 0) * mix.b[i]
      : 0.40768 * (R * c.Tc / c.Pc) * (0.29441 - (0.29056 - 0.08775 * c.omega))));
  } else {
    mix.c = new Array(n).fill(0);
  }
  return mix;
}

function checkState(mix, T, P, z) {
  need('T', T, 10, 3000);
  need('P', P, 1e-3, 1e10);
  if (!z || z.length !== mix.n) throw new RangeError('composition must have one mole fraction per component');
  let s = 0;
  for (let i = 0; i < mix.n; i++) {
    if (!(z[i] >= 0)) throw new RangeError(`mole fraction z[${i}] = ${z[i]} must be ≥ 0`);
    s += z[i];
  }
  if (Math.abs(s - 1) > 1e-6) throw new RangeError(`mole fractions sum to ${s}, expected 1`);
}

function mixParams(mix, T, x) {
  const n = mix.n;
  const s = new Array(n);
  for (let i = 0; i < n; i++) s[i] = mix.sqrtAc[i] * (1 + mix.m[i] * (1 - Math.sqrt(T / mix.Tc[i])));
  const ai = new Array(n);
  let a = 0;
  let b = 0;
  for (let i = 0; i < n; i++) {
    let sum = 0;
    const ki = mix.kij[i];
    for (let j = 0; j < n; j++) sum += x[j] * s[j] * (1 - ki[j]);
    ai[i] = s[i] * sum;
    a += x[i] * ai[i];
    b += x[i] * mix.b[i];
  }
  return { a, b, ai, s };
}

function cubicRoots(p2, p1, p0) {
  // real roots of x³ + p2 x² + p1 x + p0 = 0 (trigonometric / Cardano), polished by Newton
  const Q = (p2 * p2 - 3 * p1) / 9;
  const Rr = (2 * p2 ** 3 - 9 * p2 * p1 + 27 * p0) / 54;
  const roots = [];
  if (Rr * Rr < Q * Q * Q) {
    const th = Math.acos(Math.max(-1, Math.min(1, Rr / Math.sqrt(Q * Q * Q))));
    const sq = -2 * Math.sqrt(Q);
    roots.push(sq * Math.cos(th / 3) - p2 / 3, sq * Math.cos((th + 2 * Math.PI) / 3) - p2 / 3,
      sq * Math.cos((th - 2 * Math.PI) / 3) - p2 / 3);
  } else {
    const A = -Math.sign(Rr) * Math.cbrt(Math.abs(Rr) + Math.sqrt(Math.max(0, Rr * Rr - Q * Q * Q)));
    const B = A !== 0 ? Q / A : 0;
    roots.push(A + B - p2 / 3);
  }
  for (let k = 0; k < roots.length; k++) {
    let x = roots[k];
    for (let it = 0; it < 3; it++) {
      const f = ((x + p2) * x + p1) * x + p0;
      const df = (3 * x + 2 * p2) * x + p1;
      if (df === 0) break;
      const dx = f / df;
      if (!Number.isFinite(dx) || Math.abs(dx) > 0.1 * Math.max(Math.abs(x), 1e-3)) break;
      x -= dx;
    }
    roots[k] = x;
  }
  return roots.sort((u, v) => u - v);
}

function zRootsAB(def, A, B) {
  const u = def.d1 + def.d2;
  const w = def.d1 * def.d2;
  const r = cubicRoots(u * B - B - 1, A + w * B * B - u * B - u * B * B, -(A * B + w * B * B + w * B ** 3));
  const out = [];
  for (const z of r) if (z > B * (1 + 1e-12)) out.push(z);
  return out;
}

function gibbsOfZ(def, A, B, Z) {
  return (Z - 1) - Math.log(Z - B) - A / ((def.d1 - def.d2) * B) * Math.log((Z + def.d1 * B) / (Z + def.d2 * B));
}

// Core evaluation: Z root (by selection rule) and ln φ_i for composition x.
function state(mix, T, P, x, sel = 'stable') {
  const def = mix.def;
  const { a, b, ai } = mixParams(mix, T, x);
  const RT = R * T;
  const A = a * P / (RT * RT);
  const B = b * P / RT;
  const roots = zRootsAB(def, A, B);
  if (roots.length === 0) throw new Error(`equation of state has no physical root at T = ${T} K, P = ${P} Pa`);
  let Z;
  if (roots.length === 1) {
    Z = roots[0];
  } else if (sel === 'liquid') {
    Z = roots[0];
  } else if (sel === 'vapour') {
    Z = roots[roots.length - 1];
  } else {
    const zl = roots[0];
    const zv = roots[roots.length - 1];
    Z = gibbsOfZ(def, A, B, zl) < gibbsOfZ(def, A, B, zv) ? zl : zv;
  }
  const L = Math.log((Z + def.d1 * B) / (Z + def.d2 * B));
  const c = A / ((def.d1 - def.d2) * B);
  const n = mix.n;
  const lnPhi = new Array(n);
  const lnZB = Math.log(Z - B);
  for (let i = 0; i < n; i++) {
    const bb = mix.b[i] / b;
    lnPhi[i] = bb * (Z - 1) - lnZB - c * (2 * ai[i] / a - bb) * L;
  }
  return { Z, lnPhi, A, B, a, b, nRoots: roots.length, roots };
}

function phaseLabel(mix, st) {
  // dense (liquid-like) when the molar volume is below the EOS critical volume of the same b
  return st.Z / st.B < mix.def.vcOverB ? 'liquid' : 'vapour';
}

/**
 * Compressibility-factor roots of the cubic EOS at (T, P, x).
 * Reference: analytic solution of the cubic (Cardano / trigonometric form); roots Z ≤ B are discarded.
 * @param {object} o
 * @param {object} o.mix mixture model from `makeMixture`
 * @param {number} o.T K @param {number} o.P Pa
 * @param {number[]} o.z mole fractions, –
 * @returns {{ roots: number[], stable: number, A: number, B: number }} ascending real roots (1 or 3), the root of
 *   lowest Gibbs energy, and the dimensionless A = aP/(RT)², B = bP/(RT)
 */
export function zRoots({ mix, T, P, z }) {
  checkState(mix, T, P, z);
  const st = state(mix, T, P, z, 'stable');
  return { roots: st.roots, stable: st.Z, A: st.A, B: st.B };
}

/**
 * Fugacity coefficients of all components in a phase of composition z:
 * ln φ_i = (b_i/b)(Z − 1) − ln(Z − B) − A/((δ1 − δ2)B) (2 Σ_j x_j a_ij/a − b_i/b) ln[(Z + δ1 B)/(Z + δ2 B)].
 * Reference: Peng & Robinson (1976); Soave (1972); generic form in Michelsen & Mollerup,
 * "Thermodynamic Models" (2007) ch. 3.
 * @param {object} o
 * @param {object} o.mix @param {number} o.T K @param {number} o.P Pa @param {number[]} o.z –
 * @param {'stable'|'liquid'|'vapour'} [o.phase='stable'] root selection: lowest Gibbs energy, smallest or largest Z
 * @returns {{ phi: number[], lnPhi: number[], fugacity: number[], Z: number, phase: string }} fugacities in Pa
 */
export function fugacityCoefficients({ mix, T, P, z, phase = 'stable' }) {
  checkState(mix, T, P, z);
  const st = state(mix, T, P, z, phase);
  const phi = st.lnPhi.map(Math.exp);
  return { phi, lnPhi: st.lnPhi, fugacity: phi.map((f, i) => f * z[i] * P), Z: st.Z, phase: phaseLabel(mix, st) };
}

/**
 * Single-phase thermodynamic properties from the cubic EOS and ideal-gas heat capacities.
 *   H − H_ig = RT(Z − 1) + (T a' − a)/((δ1 − δ2) b) · ln[(v + δ1 b)/(v + δ2 b)]
 *   S − S_ig(T, P) = R ln(Z − B) + a'/((δ1 − δ2) b) · ln[(v + δ1 b)/(v + δ2 b)]
 *   Cv − Cv_ig = T a''/((δ1 − δ2) b) · ln[(v + δ1 b)/(v + δ2 b)],  Cp = Cv − T (∂P/∂T)_v²/(∂P/∂v)_T
 *   w² = −v² (Cp/Cv)(∂P/∂v)_T / M,  μ_JT = [T(∂v/∂T)_P − v]/Cp
 * with a' = da/dT, a'' = d²a/dT² of the mixture. Enthalpy and entropy are relative to the ideal gas at
 * T_ref = 298.15 K, P_ref = 101 325 Pa (pure components, so the entropy includes the ideal mixing term).
 * Caloric derivatives use the untranslated volume; `rho` and `v` include the volume shift.
 * Reference: standard residual-property relations, e.g. Michelsen & Mollerup (2007) ch. 2–3.
 * Validity: single homogeneous phase; accuracy of liquid Cp and sound speed from two-parameter cubics is modest
 * (typically 5–15 %).
 * @param {object} o
 * @param {object} o.mix @param {number} o.T K @param {number} o.P Pa @param {number[]} o.z –
 * @param {'stable'|'liquid'|'vapour'} [o.phase='stable']
 * @returns {{ Z: number, v: number, rho: number, rhoMolar: number, M: number, phase: string,
 *   hDep: number, sDep: number, h: number, s: number, cp: number, cv: number, cpIdeal: number,
 *   cpMass: number, soundSpeed: number, muJT: number, dPdT: number, dPdv: number }}
 *   v m³/mol, rho kg/m³, rhoMolar mol/m³, M kg/mol, h J/mol, s J/(mol·K), cp, cv J/(mol·K), cpMass J/(kg·K),
 *   soundSpeed m/s, muJT K/Pa, dPdT Pa/K, dPdv Pa·mol/m³
 */
export function properties({ mix, T, P, z, phase = 'stable' }) {
  checkState(mix, T, P, z);
  const def = mix.def;
  const n = mix.n;
  const st = state(mix, T, P, z, phase);
  const { a, b, Z } = st;
  const { s } = mixParams(mix, T, z);
  let da = 0;
  let d2a = 0;
  const ds = new Array(n);
  const d2s = new Array(n);
  for (let i = 0; i < n; i++) {
    const rt = Math.sqrt(T * mix.Tc[i]);
    ds[i] = -mix.sqrtAc[i] * mix.m[i] / (2 * rt);
    d2s[i] = mix.sqrtAc[i] * mix.m[i] / (4 * T * rt);
  }
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const w = z[i] * z[j] * (1 - mix.kij[i][j]);
      da += w * 2 * ds[i] * s[j];
      d2a += w * 2 * (d2s[i] * s[j] + ds[i] * ds[j]);
    }
  }
  const v = Z * R * T / P;
  const dd = def.d1 - def.d2;
  const L = Math.log((v + def.d1 * b) / (v + def.d2 * b));
  const hDep = R * T * (Z - 1) + (T * da - a) / (dd * b) * L;
  const sDep = R * Math.log(Z - st.B) + da / (dd * b) * L;
  const cvDep = T * d2a / (dd * b) * L;
  const den = (v + def.d1 * b) * (v + def.d2 * b);
  const dPdT = R / (v - b) - da / den;
  const dPdv = -R * T / ((v - b) * (v - b)) + a * (2 * v + (def.d1 + def.d2) * b) / (den * den);
  let cpIdeal = 0;
  let hIg = 0;
  let sIg = 0;
  let M = 0;
  let cShift = 0;
  for (let i = 0; i < n; i++) {
    M += z[i] * mix.M[i];
    cShift += z[i] * mix.c[i];
    if (z[i] <= 0) continue;
    cpIdeal += z[i] * idealGasCp(mix.comps[i], T);
    const hs = idealGasEnthalpyEntropy(mix.comps[i], T);
    hIg += z[i] * hs.h;
    sIg += z[i] * (hs.s - R * Math.log(z[i]));
  }
  sIg -= R * Math.log(P / 101325);
  const cv = cpIdeal - R + cvDep;
  const cp = cv - T * dPdT * dPdT / dPdv;
  const vT = v - cShift;
  return {
    Z, v: vT, rho: M / vT, rhoMolar: 1 / vT, M, phase: phaseLabel(mix, st),
    hDep, sDep, h: hIg + hDep, s: sIg + sDep, cp, cv, cpIdeal, cpMass: cp / M,
    soundSpeed: Math.sqrt(Math.max(0, -v * v * (cp / cv) * dPdv / M)),
    muJT: -(T * dPdT / dPdv + v) / cp, dPdT, dPdv,
  };
}

/**
 * Wilson (1968) K-value estimate, K_i = (Pc_i/P) exp[5.373 (1 + ω_i)(1 − Tc_i/T)].
 * Reference: Wilson, G.M. (1968), "A modified Redlich–Kwong equation of state", 65th AIChE meeting, Cleveland.
 * Validity: initial estimate only, reasonable below about 3.5 MPa.
 * @param {object} mix mixture model
 * @param {number} T K
 * @param {number} P Pa
 * @returns {number[]} K-values, –
 */
export function wilsonK(mix, T, P) {
  need('T', T, 1, 5000);
  need('P', P, 1e-6, 1e10);
  return mix.Tc.map((Tc, i) => (mix.Pc[i] / P) * Math.exp(5.373 * (1 + mix.omega[i]) * (1 - Tc / T)));
}

/**
 * Solve the Rachford–Rice equation Σ z_i (K_i − 1)/(1 + β (K_i − 1)) = 0 for the vapour mole fraction β.
 * The root is sought on the full physical window (1/(1 − Kmax), 1/(1 − Kmin)), so β may fall outside [0, 1]
 * ("negative flash"). If all K < 1 the function returns 0, if all K > 1 it returns 1.
 * Reference: Rachford & Rice, JPT 4(10) (1952); Whitson & Michelsen, Fluid Phase Equilib. 53 (1989) 51.
 * @param {number[]} z feed mole fractions, –
 * @param {number[]} K equilibrium ratios y/x, –
 * @returns {number} vapour fraction β, – (mole basis)
 */
export function rachfordRice(z, K) {
  if (z.length !== K.length) throw new RangeError('z and K must have the same length');
  let kmin = Infinity;
  let kmax = 0;
  for (let i = 0; i < K.length; i++) {
    if (!(K[i] > 0) || !Number.isFinite(K[i])) throw new RangeError(`K[${i}] = ${K[i]} must be positive and finite`);
    if (z[i] > 0) { kmin = Math.min(kmin, K[i]); kmax = Math.max(kmax, K[i]); }
  }
  if (kmax <= 1) return 0;
  if (kmin >= 1) return 1;
  let lo = 1 / (1 - kmax);
  let hi = 1 / (1 - kmin);
  const span = hi - lo;
  lo += 1e-12 * span;
  hi -= 1e-12 * span;
  const f = (bt) => {
    let s = 0;
    let d = 0;
    for (let i = 0; i < z.length; i++) {
      const k1 = K[i] - 1;
      const t = k1 / (1 + bt * k1);
      s += z[i] * t;
      d -= z[i] * t * t;
    }
    return [s, d];
  };
  let bt = Math.min(Math.max(0.5, lo + 0.25 * span), hi - 0.25 * span);
  for (let it = 0; it < 200; it++) {
    const [s, d] = f(bt);
    if (s > 0) lo = bt; else hi = bt;
    let nb = bt - s / d;
    if (!(nb > lo && nb < hi)) nb = 0.5 * (lo + hi);
    if (Math.abs(nb - bt) < 1e-15 * Math.max(1, Math.abs(bt)) || hi - lo < 1e-16) { bt = nb; break; }
    bt = nb;
  }
  return bt;
}

function normalised(v) {
  let s = 0;
  for (let i = 0; i < v.length; i++) s += v[i];
  return v.map((x) => x / s);
}

function stabilityTrial(mix, T, P, z, d, W0, maxIter) {
  const n = mix.n;
  let lnW = W0.map((w) => Math.log(Math.max(w, 1e-300)));
  let prev = null;
  let sumW = 1;
  let trivial = false;
  let converged = false;
  for (let it = 0; it < maxIter; it++) {
    const W = lnW.map(Math.exp);
    sumW = 0;
    for (let i = 0; i < n; i++) sumW += W[i];
    const w = W.map((x) => x / sumW);
    const st = state(mix, T, P, w, 'stable');
    const g = new Array(n);
    let err = 0;
    let triv = 0;
    for (let i = 0; i < n; i++) {
      const target = z[i] > 0 ? d[i] - st.lnPhi[i] : -700;
      g[i] = target - lnW[i];
      if (z[i] > 0) {
        err = Math.max(err, Math.abs(g[i]));
        const r = target - Math.log(z[i]);
        triv += r * r;
      }
    }
    if (triv < 1e-8) { trivial = true; break; }
    if (err < 1e-11) { converged = true; break; }
    // dominant-eigenvalue acceleration every fifth iteration
    let acc = 1;
    if (prev && it % 5 === 4) {
      let num = 0;
      let den = 0;
      for (let i = 0; i < n; i++) { num += g[i] * prev[i]; den += prev[i] * prev[i]; }
      const lam = den > 0 ? num / den : 0;
      if (lam > 0.2 && lam < 0.98) acc = 1 / (1 - lam);
    }
    for (let i = 0; i < n; i++) lnW[i] += acc * g[i];
    prev = g;
  }
  const W = lnW.map(Math.exp);
  sumW = W.reduce((s, x) => s + x, 0);
  return { W, sumW, trivial, converged };
}

/**
 * Michelsen tangent-plane stability test of a single-phase mixture at (T, P).
 * Two trial phases (vapour-like W = zK and liquid-like W = z/K with Wilson K) are converged by successive
 * substitution on ln W_i = ln z_i + ln φ_i(z) − ln φ_i(W) with dominant-eigenvalue acceleration; the mixture is
 * unstable when a non-trivial stationary point has Σ W_i > 1 (tangent-plane distance tm = 1 − ΣW < 0).
 * Reference: Michelsen, Fluid Phase Equilib. 9 (1982) 1–19.
 * Validity: two-phase vapour–liquid instability of cubic-EOS mixtures; a third (aqueous) phase is not searched.
 * @param {object} o
 * @param {object} o.mix @param {number} o.T K @param {number} o.P Pa @param {number[]} o.z –
 * @param {number} [o.maxIter=500]
 * @returns {{ stable: boolean, tm: number, K: number[]|null }} tm = minimum modified tangent-plane distance (–);
 *   K = K-value estimate for a subsequent flash when unstable
 */
export function stabilityTest({ mix, T, P, z, maxIter = 500 }) {
  checkState(mix, T, P, z);
  const n = mix.n;
  if (n === 1) return { stable: true, tm: 0, K: null };
  const st = state(mix, T, P, z, 'stable');
  const d = z.map((zi, i) => (zi > 0 ? Math.log(zi) + st.lnPhi[i] : 0));
  const Kw = wilsonK(mix, T, P);
  const tv = stabilityTrial(mix, T, P, z, d, z.map((zi, i) => zi * Kw[i]), maxIter);
  const tl = stabilityTrial(mix, T, P, z, d, z.map((zi, i) => zi / Kw[i]), maxIter);
  const tol = 1e-9;
  const uv = !tv.trivial && tv.sumW > 1 + tol;
  const ul = !tl.trivial && tl.sumW > 1 + tol;
  let tm = 0;
  if (!tv.trivial) tm = Math.min(tm, 1 - tv.sumW);
  if (!tl.trivial) tm = Math.min(tm, 1 - tl.sumW);
  if (!uv && !ul) return { stable: true, tm, K: null };
  let K;
  if (uv && ul) K = z.map((zi, i) => (zi > 0 ? tv.W[i] / tl.W[i] : Kw[i]));
  else if (uv) K = z.map((zi, i) => (zi > 0 ? tv.W[i] / zi : Kw[i]));
  else K = z.map((zi, i) => (zi > 0 ? zi / tl.W[i] : Kw[i]));
  return { stable: false, tm, K };
}

function singlePhaseResult(mix, T, P, z, extra = {}) {
  const st = state(mix, T, P, z, 'stable');
  const ph = phaseLabel(mix, st);
  return {
    phases: 1, phase: ph, beta: ph === 'vapour' ? 1 : 0, x: z.slice(), y: z.slice(), K: null,
    ZL: st.Z, ZV: st.Z, iterations: 0, converged: true, ...extra,
  };
}

/**
 * Isothermal two-phase (vapour–liquid) PT flash.
 * Algorithm: Michelsen stability test of the feed (skipped when `K0` is given); if unstable, successive
 * substitution on ln K_i = ln φ_i^L − ln φ_i^V with Rachford–Rice inner solve and dominant-eigenvalue
 * (first-order GDEM) acceleration every fifth iteration. Each phase uses the Z root of lowest Gibbs energy.
 * Convergence: max |Δ ln K| < tol. Convergence slows close to the mixture critical point; `converged` reports it.
 * References: Michelsen, Fluid Phase Equilib. 9 (1982) 1–19 and 21–40; Crowe & Nishio, AIChE J. 21 (1975) 528.
 * Validity: two fluid phases; no aqueous or solid phase.
 * @param {object} o
 * @param {object} o.mix mixture model
 * @param {number} o.T K @param {number} o.P Pa
 * @param {number[]} o.z feed mole fractions, –
 * @param {number[]} [o.K0] initial K-values (skips the stability test)
 * @param {number} [o.tol=1e-11] @param {number} [o.maxIter=5000]
 * @returns {{ phases: 1|2, phase: string, beta: number, x: number[], y: number[], K: number[]|null,
 *   ZL: number, ZV: number, iterations: number, converged: boolean }}
 *   beta = vapour mole fraction; x, y = liquid and vapour mole fractions; phase = 'vapour' | 'liquid' |
 *   'two-phase'
 */
export function flashPT({ mix, T, P, z, K0, tol = 1e-11, maxIter = 5000 }) {
  checkState(mix, T, P, z);
  const n = mix.n;
  const zn = normalised(z);
  if (n === 1) return singlePhaseResult(mix, T, P, zn);
  let K = K0 ? K0.slice() : null;
  if (!K) {
    const st = stabilityTest({ mix, T, P, z: zn });
    if (st.stable) return singlePhaseResult(mix, T, P, zn);
    K = st.K;
  }
  let lnK = K.map(Math.log);
  let prev = null;
  let beta = 0.5;
  let x = zn;
  let y = zn;
  let sl = null;
  let sv = null;
  let converged = false;
  let it = 0;
  for (; it < maxIter; it++) {
    K = lnK.map(Math.exp);
    beta = rachfordRice(zn, K);
    x = normalised(zn.map((zi, i) => zi / (1 + beta * (K[i] - 1))));
    y = normalised(zn.map((zi, i) => K[i] * zi / (1 + beta * (K[i] - 1))));
    sl = state(mix, T, P, x, 'stable');
    sv = state(mix, T, P, y, 'stable');
    const g = new Array(n);
    let err = 0;
    let triv = 0;
    for (let i = 0; i < n; i++) {
      g[i] = sl.lnPhi[i] - sv.lnPhi[i] - lnK[i];
      if (zn[i] > 0) err = Math.max(err, Math.abs(g[i]));
      triv += lnK[i] * lnK[i];
    }
    if (triv < 1e-10) return singlePhaseResult(mix, T, P, zn, { iterations: it });
    if (err < tol) { converged = true; break; }
    let acc = 1;
    if (prev && it % 5 === 4) {
      let num = 0;
      let den = 0;
      for (let i = 0; i < n; i++) { num += g[i] * prev[i]; den += prev[i] * prev[i]; }
      const lam = den > 0 ? num / den : 0;
      if (lam > 0.2 && lam < 0.98) acc = 1 / (1 - lam);
    }
    for (let i = 0; i < n; i++) lnK[i] += acc * g[i];
    prev = g;
  }
  if (beta <= 0 || beta >= 1) return singlePhaseResult(mix, T, P, zn, { iterations: it });
  return { phases: 2, phase: 'two-phase', beta, x, y, K, ZL: sl.Z, ZV: sv.Z, iterations: it, converged };
}

// ---------------------------------------------------------------------------------------------------------------
// Saturation points

function wilsonPsat(mix, T) {
  return mix.Tc.map((Tc, i) => mix.Pc[i] * Math.exp(5.373 * (1 + mix.omega[i]) * (1 - Tc / T)));
}

function wilsonSaturationT(mix, z, P, beta) {
  // T at which Σ zK = 1 (bubble, beta = 0) or Σ z/K = 1 (dew, beta = 1) with Wilson K
  const f = (T) => {
    const K = wilsonK(mix, T, P);
    let s = 0;
    for (let i = 0; i < mix.n; i++) if (z[i] > 0) s += beta === 0 ? z[i] * K[i] : z[i] / K[i];
    return beta === 0 ? s - 1 : 1 - s;
  };
  let lo = 20;
  let hi = 2000;
  for (let it = 0; it < 80; it++) {
    const m = 0.5 * (lo + hi);
    if (f(m) > 0) hi = m; else lo = m;
  }
  return 0.5 * (lo + hi);
}

function satResidual(mix, z, beta, X, F) {
  const n = mix.n;
  const T = Math.exp(X[n]);
  const P = Math.exp(X[n + 1]);
  const x = new Array(n);
  const y = new Array(n);
  let sx = 0;
  let sy = 0;
  for (let i = 0; i < n; i++) {
    const K = Math.exp(X[i]);
    x[i] = z[i] / (1 + beta * (K - 1));
    y[i] = K * x[i];
    sx += x[i];
    sy += y[i];
  }
  for (let i = 0; i < n; i++) { x[i] /= sx; y[i] /= sy; }
  // lowest-Gibbs root for each phase, so that the same equations hold on both sides of the critical point
  const L = state(mix, T, P, x, 'stable');
  const V = state(mix, T, P, y, 'stable');
  for (let i = 0; i < n; i++) F[i] = z[i] > 0 ? X[i] + V.lnPhi[i] - L.lnPhi[i] : 0;
  F[n] = sy - sx;
  return { x, y };
}

// Newton solution of the n + 2 saturation equations with one specified variable X[spec] = S.
function satNewton(mix, z, beta, X0, spec, S, maxIter = 40) {
  const n = mix.n;
  const N = n + 2;
  const X = X0.slice();
  const F = new Array(N).fill(0);
  const Fp = new Array(N).fill(0);
  let J = null;
  let comp = null;
  for (let it = 0; it < maxIter; it++) {
    try {
      comp = satResidual(mix, z, beta, X, F);
    } catch (e) { return { converged: false }; }
    F[N - 1] = X[spec] - S;
    let err = 0;
    for (let i = 0; i < N; i++) err = Math.max(err, Math.abs(F[i]));
    if (!Number.isFinite(err)) return { converged: false };
    J = [];
    for (let i = 0; i < N; i++) J.push(new Array(N).fill(0));
    for (let j = 0; j < N; j++) {
      const h = 1e-6;
      const save = X[j];
      X[j] = save + h;
      try { satResidual(mix, z, beta, X, Fp); } catch (e) { X[j] = save; return { converged: false }; }
      X[j] = save;
      for (let i = 0; i < N - 1; i++) J[i][j] = (Fp[i] - F[i]) / h;
    }
    J[N - 1][spec] = 1;
    if (err < 1e-9) return { converged: true, X, J, iterations: it, x: comp.x, y: comp.y };
    let dX;
    try { dX = solveDense(J, F.map((v) => -v)); } catch (e) { return { converged: false }; }
    let big = 0;
    for (let i = 0; i < N; i++) {
      if (!Number.isFinite(dX[i])) return { converged: false };
      big = Math.max(big, Math.abs(dX[i]));
    }
    const damp = big > 1 ? 1 / big : 1;
    for (let i = 0; i < N; i++) X[i] += damp * dX[i];
  }
  return { converged: false };
}

function maxAbsLnK(X, z, n) {
  let m = 0;
  for (let i = 0; i < n; i++) if (z[i] > 0) m = Math.max(m, Math.abs(X[i]));
  return m;
}

// Successive substitution for a saturation pressure at fixed T (beta = 0 bubble, 1 dew).
function satPressureSS(mix, T, z, beta, P0, K0, maxIter = 200) {
  const n = mix.n;
  let P = P0;
  let K = K0.slice();
  let w = z;
  for (let it = 0; it < maxIter; it++) {
    let s = 0;
    const t = new Array(n);
    for (let i = 0; i < n; i++) { t[i] = beta === 0 ? z[i] * K[i] : z[i] / K[i]; s += t[i]; }
    w = t.map((v) => v / s);
    let feed;
    let inc;
    try {
      feed = state(mix, T, P, z, beta === 0 ? 'liquid' : 'vapour');
      inc = state(mix, T, P, w, beta === 0 ? 'vapour' : 'liquid');
    } catch (e) { return { converged: false }; }
    let triv = 0;
    for (let i = 0; i < n; i++) {
      const lnK = beta === 0 ? feed.lnPhi[i] - inc.lnPhi[i] : inc.lnPhi[i] - feed.lnPhi[i];
      K[i] = Math.exp(lnK);
      if (z[i] > 0) triv = Math.max(triv, Math.abs(lnK));
    }
    if (n > 1 && triv < 2e-3) return { converged: false, trivial: true };
    s = 0;
    for (let i = 0; i < n; i++) s += beta === 0 ? z[i] * K[i] : z[i] / K[i];
    const Pn = beta === 0 ? P * s : P / s;
    if (!Number.isFinite(Pn) || Pn <= 0) return { converged: false };
    if (Math.abs(s - 1) < 1e-11) return { converged: true, P: Pn, K, w, iterations: it };
    P = Pn;
  }
  return { converged: false, P, K, w };
}

function saturationPressure(mix, T, z, beta, P0) {
  const n = mix.n;
  const ps = wilsonPsat(mix, T);
  let Pinit = P0;
  if (Pinit == null) {
    let s = 0;
    for (let i = 0; i < n; i++) if (z[i] > 0) s += beta === 0 ? z[i] * ps[i] : z[i] / ps[i];
    Pinit = beta === 0 ? s : 1 / s;
  }
  const K0 = ps.map((p) => p / Pinit);
  const ss = satPressureSS(mix, T, z, beta, Pinit, K0, n === 1 ? 400 : 60);
  if (ss.converged) return { P: ss.P, K: ss.K, incipient: ss.w, converged: true, method: 'successive-substitution' };
  if (n > 1 && ss.P && ss.K && !ss.trivial) {
    const X = [...ss.K.map(Math.log), Math.log(T), Math.log(ss.P)];
    const nw = satNewton(mix, z, beta, X, n, Math.log(T));
    if (nw.converged && maxAbsLnK(nw.X, z, n) > 1e-4) {
      return { P: Math.exp(nw.X[n + 1]), K: nw.X.slice(0, n).map(Math.exp), incipient: beta === 0 ? nw.y : nw.x,
        converged: true, method: 'newton' };
    }
  }
  return { converged: false };
}

// A converged saturation point is accepted only if the feed is two-phase just inside and stable just outside
// (guards against the trivial solution and against solutions on a metastable root). side = +1: single phase
// lies above P (bubble / upper dew); side = −1: single phase lies below P (lower dew).
function isBoundary(mix, T, P, z, side) {
  const e = 2e-3;
  const inside = stabilityTest({ mix, T, P: P * (1 - side * e), z });
  if (inside.stable) return false;
  return stabilityTest({ mix, T, P: P * (1 + side * e), z }).stable;
}

/**
 * Upper phase-boundary pressure at temperature T located by bisection on the stability test: the highest
 * pressure at which the mixture splits into two phases. The boundary is a bubble point when the vapour fraction
 * just inside is below 0.5, otherwise an (upper, retrograde) dew point.
 * Reference: Michelsen (1982) stability analysis used as a phase-boundary detector.
 * Validity: T below the cricondentherm; throws if no two-phase pressure is found between `Pmin` and `Pmax`.
 * @param {object} o
 * @param {object} o.mix @param {number} o.T K @param {number[]} o.z –
 * @param {number} [o.Pmin=1e4] Pa @param {number} [o.Pmax=1.5e8] Pa
 * @returns {{ P: number, type: 'bubble'|'dew', beta: number }} boundary pressure Pa, type, vapour fraction just
 *   inside the envelope
 */
export function upperSaturationPressure({ mix, T, z, Pmin = 1e4, Pmax = 1.5e8 }) {
  checkState(mix, T, Pmin, z);
  const unstable = (P) => !stabilityTest({ mix, T, P, z }).stable;
  // scan downwards on a geometric grid for the first two-phase pressure
  const nScan = 60;
  const ratio = (Pmin / Pmax) ** (1 / nScan);
  let hi = Pmax;
  let lo = null;
  if (unstable(Pmax)) throw new RangeError(`Pmax = ${Pmax} Pa is still inside the two-phase region`);
  for (let i = 1; i <= nScan; i++) {
    const P = Pmax * ratio ** i;
    if (unstable(P)) { lo = P; break; }
    hi = P;
  }
  if (lo == null) throw new RangeError(`no two-phase region found at T = ${T} K between Pmin and Pmax`);
  for (let it = 0; it < 60 && (hi - lo) > 1e-9 * hi; it++) {
    const m = 0.5 * (lo + hi);
    if (unstable(m)) lo = m; else hi = m;
  }
  const fl = flashPT({ mix, T, P: lo * (1 - 1e-6), z });
  const beta = fl.phases === 2 ? fl.beta : (fl.phase === 'vapour' ? 1 : 0);
  return { P: 0.5 * (lo + hi), type: beta < 0.5 ? 'bubble' : 'dew', beta };
}

/**
 * Bubble-point pressure at temperature T: Σ z_i K_i = 1 with K_i = φ_i^L(z)/φ_i^V(y).
 * Method: successive substitution from a Wilson estimate, Newton on (ln K, ln P) if that stalls, and bisection on
 * the stability test as a last resort (then `method` = 'stability-bisection').
 * Reference: Michelsen & Mollerup (2007) ch. 12. Validity: T below the mixture critical temperature; for a pure
 * component this is the vapour pressure.
 * @param {object} o
 * @param {object} o.mix @param {number} o.T K @param {number[]} o.z –
 * @param {number} [o.P0] initial estimate, Pa
 * @returns {{ P: number, y: number[], K: number[]|null, converged: boolean, method: string }} P in Pa, incipient
 *   vapour composition
 */
export function bubblePointPressure({ mix, T, z, P0 }) {
  checkState(mix, T, 1e5, z);
  const zn = normalised(z);
  const r = saturationPressure(mix, T, zn, 0, P0);
  if (r.converged && (mix.n === 1 || isBoundary(mix, T, r.P, zn, 1))) {
    return { P: r.P, y: r.incipient, K: r.K, converged: true, method: r.method };
  }
  if (mix.n === 1) throw new RangeError(`T = ${T} K is above the critical temperature; no vapour pressure exists`);
  const u = upperSaturationPressure({ mix, T, z: zn });
  if (u.type !== 'bubble') throw new RangeError(`no bubble point at T = ${T} K (upper boundary is a dew point)`);
  const fl = flashPT({ mix, T, P: u.P * (1 - 1e-6), z: zn });
  return { P: u.P, y: fl.y, K: fl.K, converged: true, method: 'stability-bisection' };
}

/**
 * Dew-point pressure at temperature T: Σ z_i/K_i = 1. Returns the lower (normal) dew point; the upper retrograde
 * dew point, where one exists, is given by `upperSaturationPressure`.
 * Method: successive substitution from a Wilson estimate, then Newton on (ln K, ln P).
 * Reference: Michelsen & Mollerup (2007) ch. 12. Validity: T below the cricondentherm.
 * @param {object} o
 * @param {object} o.mix @param {number} o.T K @param {number[]} o.z –
 * @param {number} [o.P0] initial estimate, Pa
 * @returns {{ P: number, x: number[], K: number[], converged: boolean, method: string }} P in Pa, incipient liquid
 */
export function dewPointPressure({ mix, T, z, P0 }) {
  checkState(mix, T, 1e5, z);
  const zn = normalised(z);
  const r = saturationPressure(mix, T, zn, 1, P0);
  if (!r.converged || (mix.n > 1 && !isBoundary(mix, T, r.P, zn, -1))) throw new RangeError(`dew-point pressure did not converge at T = ${T} K (above the cricondentherm?)`);
  return { P: r.P, x: r.incipient, K: r.K, converged: true, method: r.method };
}

function saturationTemperature(mix, P, z, beta, T0) {
  const n = mix.n;
  const zn = normalised(z);
  const Tw = T0 ?? wilsonSaturationT(mix, zn, P, beta);
  const lnP = Math.log(P);
  if (n > 1) {
    const X = [...wilsonK(mix, Tw, P).map(Math.log), Math.log(Tw), lnP];
    const nw = satNewton(mix, zn, beta, X, n + 1, lnP, 60);
    if (nw.converged && maxAbsLnK(nw.X, zn, n) > 1e-4) {
      return { T: Math.exp(nw.X[n]), K: nw.X.slice(0, n).map(Math.exp), incipient: beta === 0 ? nw.y : nw.x,
        converged: true, method: 'newton' };
    }
  }
  // secant on ln Psat(T) − ln P
  const f = (T) => {
    const r = saturationPressure(mix, T, zn, beta);
    return r.converged ? { v: Math.log(r.P) - lnP, r } : null;
  };
  let t0 = Tw;
  let f0 = f(t0);
  for (let k = 0; k < 40 && !f0; k++) { t0 *= 0.97; f0 = f(t0); }
  if (!f0) return { converged: false };
  let t1 = t0 * (f0.v > 0 ? 0.98 : 1.02);
  let f1 = f(t1);
  if (!f1) { t1 = t0 * 0.99; f1 = f(t1); }
  if (!f1) return { converged: false };
  for (let it = 0; it < 60; it++) {
    if (Math.abs(f1.v) < 1e-10) {
      return { T: t1, K: f1.r.K, incipient: f1.r.incipient, converged: true, method: 'secant' };
    }
    const slope = (f1.v - f0.v) / (t1 - t0);
    if (!(slope !== 0) || !Number.isFinite(slope)) break;
    let t2 = t1 - f1.v / slope;
    if (t2 > 1.1 * t1) t2 = 1.1 * t1;
    if (t2 < 0.9 * t1) t2 = 0.9 * t1;
    let f2 = f(t2);
    for (let k = 0; k < 30 && !f2; k++) { t2 = 0.5 * (t2 + t1); f2 = f(t2); }
    if (!f2) break;
    t0 = t1; f0 = f1; t1 = t2; f1 = f2;
  }
  return { converged: false };
}

/**
 * Bubble-point temperature at pressure P (Newton on the saturation equations from a Wilson estimate, secant on
 * ln P_bub(T) as fallback). Reference: Michelsen & Mollerup (2007) ch. 12.
 * Validity: P below the mixture critical pressure (bubble branch).
 * @param {object} o
 * @param {object} o.mix @param {number} o.P Pa @param {number[]} o.z – @param {number} [o.T0] initial estimate, K
 * @returns {{ T: number, y: number[], K: number[], converged: boolean, method: string }} T in K
 */
export function bubblePointTemperature({ mix, P, z, T0 }) {
  checkState(mix, 300, P, z);
  const r = saturationTemperature(mix, P, z, 0, T0);
  if (!r.converged) throw new RangeError(`bubble-point temperature did not converge at P = ${P} Pa`);
  return { T: r.T, y: r.incipient, K: r.K, converged: true, method: r.method };
}

/**
 * Dew-point temperature at pressure P (normal dew point on the branch connected to low pressure).
 * Method and reference as `bubblePointTemperature`. Validity: P below the cricondenbar.
 * @param {object} o
 * @param {object} o.mix @param {number} o.P Pa @param {number[]} o.z – @param {number} [o.T0] initial estimate, K
 * @returns {{ T: number, x: number[], K: number[], converged: boolean, method: string }} T in K
 */
export function dewPointTemperature({ mix, P, z, T0 }) {
  checkState(mix, 300, P, z);
  const r = saturationTemperature(mix, P, z, 1, T0);
  if (!r.converged) throw new RangeError(`dew-point temperature did not converge at P = ${P} Pa`);
  return { T: r.T, x: r.incipient, K: r.K, converged: true, method: r.method };
}

/**
 * Vapour pressure of a pure component from the EOS (equal fugacity of the liquid and vapour roots).
 * Reference: Maxwell equal-area / isofugacity condition. Validity: triple point < T < Tc.
 * @param {object} o
 * @param {string|object} o.component database id or component record
 * @param {number} o.T K
 * @param {'PR'|'SRK'} [o.eos='PR']
 * @returns {number} vapour pressure, Pa
 */
export function vapourPressure({ component, T, eos = 'PR' }) {
  const mix = makeMixture([component], { eos });
  need('T', T, 10, mix.Tc[0] * (1 - 1e-6));
  return bubblePointPressure({ mix, T, z: [1] }).P;
}

function traceBranch(mix, z, beta, { Pstart, Pmax, Tmin, maxPoints }) {
  const n = mix.n;
  const N = n + 2;
  const out = { pre: [], post: [], critical: null, closed: false };
  const T0 = wilsonSaturationT(mix, z, Pstart, beta);
  let X = [...wilsonK(mix, T0, Pstart).map(Math.log), Math.log(T0), Math.log(Pstart)];
  let sol = satNewton(mix, z, beta, X, n + 1, X[n + 1], 60);
  if (!sol.converged || maxAbsLnK(sol.X, z, n) < 1e-3) {
    // fall back to the fixed-temperature solver for the start point
    const sp = saturationPressure(mix, T0, z, beta);
    if (!sp.converged) return out;
    X = [...sp.K.map(Math.log), Math.log(T0), Math.log(sp.P)];
    sol = satNewton(mix, z, beta, X, n, X[n], 60);
    if (!sol.converged) return out;
  }
  X = sol.X;
  // reference component: the one with the largest |ln K| at the start; its sign flips at the critical point
  let ref = 0;
  for (let i = 1; i < n; i++) if (z[i] > 0 && Math.abs(X[i]) > Math.abs(X[ref])) ref = i;
  const sign0 = Math.sign(X[ref]);
  let crossed = false;
  const point = () => ({ T: Math.exp(X[n]), P: Math.exp(X[n + 1]), k: X[ref] * sign0 });
  out.pre.push(point());
  let dir = null; // previous step vector, to keep marching the same way along the curve
  let frac = 0.5;
  for (let count = 0; count < maxPoints; count++) {
    const rhs = new Array(N).fill(0);
    rhs[N - 1] = 1;
    let sens;
    try { sens = solveDense(sol.J, rhs); } catch (e) { break; }
    if (sens.some((v) => !Number.isFinite(v))) break;
    // allowed change per step: ln K limited more tightly near the critical point so that it is bracketed closely
    const kLim = Math.max(0.06, Math.min(0.5, 0.4 * maxAbsLnK(X, z, n)));
    const lim = (i) => (i < n ? kLim : i === n ? 0.05 : 0.35);
    let scale = 0;
    let ns = n + 1;
    for (let i = 0; i < N; i++) {
      if (i < n && !(z[i] > 0)) continue;
      const w = Math.abs(sens[i]) / lim(i);
      if (w > scale) { scale = w; ns = i; }
    }
    if (!(scale > 0)) break;
    const tan = sens.map((v) => v / scale);
    let sign = 1;
    if (dir) {
      let d = 0;
      for (let i = 0; i < N; i++) d += tan[i] * dir[i];
      sign = d >= 0 ? 1 : -1;
    } else {
      sign = tan[n + 1] >= 0 ? 1 : -1; // start by increasing pressure
    }
    let ok = false;
    let next = null;
    for (let tries = 0; tries < 14; tries++) {
      const Xp = X.map((v, i) => v + sign * frac * tan[i]);
      next = satNewton(mix, z, beta, Xp, ns, Xp[ns], 25);
      if (next.converged && maxAbsLnK(next.X, z, n) > 1e-3
        && Math.abs(next.X[n] - X[n]) < 0.15 && Math.abs(next.X[n + 1] - X[n + 1]) < 1.0) { ok = true; break; }
      frac *= 0.5;
      if (frac < 1e-4) break;
    }
    if (!ok) break;
    const before = point();
    dir = next.X.map((v, i) => v - X[i]);
    X = next.X;
    sol = next;
    const now = point();
    if (!crossed && now.k < 0) {
      // ln K of the reference component is linear along the curve through the critical point
      const w = before.k / (before.k - now.k);
      out.critical = { T: before.T + w * (now.T - before.T), P: before.P + w * (now.P - before.P) };
      crossed = true;
    }
    (crossed ? out.post : out.pre).push(now);
    if (now.P > Pmax || now.T < Tmin) break;
    if (now.P < Pstart) { out.closed = true; break; }
    if (next.iterations <= 3) frac = Math.min(1, frac * 1.5);
    else if (next.iterations >= 6) frac *= 0.6;
  }
  return out;
}

function refineExtremum(pts, key) {
  // parabolic refinement of the maximum of pts[i][key] through its neighbours (parametrised by index)
  let k = 0;
  for (let i = 1; i < pts.length; i++) if (pts[i][key] > pts[k][key]) k = i;
  const other = key === 'P' ? 'T' : 'P';
  if (k === 0 || k === pts.length - 1) return { T: pts[k].T, P: pts[k].P };
  const y0 = pts[k - 1][key];
  const y1 = pts[k][key];
  const y2 = pts[k + 1][key];
  const den = y0 - 2 * y1 + y2;
  const s = den !== 0 ? 0.5 * (y0 - y2) / den : 0;
  const quad = (a, b, c) => b + 0.5 * s * (c - a) + 0.5 * s * s * (a - 2 * b + c);
  const out = {};
  out[key] = quad(y0, y1, y2);
  out[other] = quad(pts[k - 1][other], pts[k][other], pts[k + 1][other]);
  return out;
}

/**
 * Trace the vapour–liquid phase envelope of a mixture in the P–T plane.
 * Method: Newton continuation on the n + 2 saturation equations in (ln K_i, ln T, ln P) at fixed vapour fraction,
 * switching the specified variable to the fastest-changing one at every step (Michelsen 1980). The curve is
 * started at `Pstart` on the bubble line (vapour fraction 0) and on the dew line (vapour fraction 1); each trace
 * runs up to and through the critical point, where all ln K change sign, and continues down the other line. The
 * critical point is located by linear interpolation of ln K to zero between the two bracketing points and the
 * estimates of the two traces are averaged. The reported bubble and dew lines use the part of each trace before its
 * critical crossing, completed by the continuation of the other trace when one of them could not be finished.
 * Cricondenbar and cricondentherm are the pressure and temperature maxima of the traced points, refined by a
 * parabola through the neighbouring points.
 * Accuracy: saturation points are converged to 1e-9 in the residuals; the critical point is approximate (usually
 * within about 1 K and 1 bar of the EOS critical point) because it is interpolated, not solved for.
 * `complete` is false when neither trace could bracket the critical point.
 * Reference: Michelsen, Fluid Phase Equilib. 4 (1980) 1–10.
 * Validity: mixtures (n ≥ 2) without a separate aqueous phase.
 * @param {object} o
 * @param {object} o.mix @param {number[]} o.z –
 * @param {number} [o.Pstart=1e5] starting pressure of both traces, Pa
 * @param {number} [o.Pmax=2e8] Pa @param {number} [o.Tmin=40] K
 * @param {number} [o.maxPoints=600] per trace
 * @returns {{ bubble: {T: number, P: number}[], dew: {T: number, P: number}[],
 *   points: {T: number, P: number, type: string}[],
 *   critical: {T: number, P: number, approximate: true}|null,
 *   cricondenbar: {T: number, P: number}, cricondentherm: {T: number, P: number}, complete: boolean }}
 *   `points` runs along the envelope: bubble line from low pressure up to the critical point, then the dew line
 *   back down to low pressure. T in K, P in Pa.
 */
export function phaseEnvelope({ mix, z, Pstart = 1e5, Pmax = 2e8, Tmin = 40, maxPoints = 600 }) {
  checkState(mix, 300, Pstart, z);
  if (mix.n < 2) throw new RangeError('phaseEnvelope needs a mixture of at least two components');
  const zn = normalised(z);
  const opts = { Pstart, Pmax, Tmin, maxPoints };
  const tb = traceBranch(mix, zn, 0, opts);
  const td = traceBranch(mix, zn, 1, opts);
  const crits = [tb.critical, td.critical].filter(Boolean);
  const critical = crits.length ? {
    T: crits.reduce((s, c) => s + c.T, 0) / crits.length, P: crits.reduce((s, c) => s + c.P, 0) / crits.length,
    approximate: true } : null;
  // bubble line: low pressure → critical; dew line: low pressure → critical
  let bub = tb.pre;
  let dew = td.pre;
  if (!tb.critical && td.critical && td.post.length > 0) {
    const tail = td.post.slice().reverse();
    const lastT = bub.length ? bub[bub.length - 1].T : -Infinity;
    bub = td.closed ? tail : bub.concat(tail.filter((p) => p.T > lastT));
  }
  if (!td.critical && tb.critical && tb.post.length > 0 && tb.closed) dew = tb.post.slice().reverse();
  const strip = (p) => ({ T: p.T, P: p.P });
  const points = [];
  for (const p of bub) points.push({ T: p.T, P: p.P, type: 'bubble' });
  if (critical) points.push({ T: critical.T, P: critical.P, type: 'critical' });
  for (let i = dew.length - 1; i >= 0; i--) points.push({ T: dew[i].T, P: dew[i].P, type: 'dew' });
  if (points.length === 0) throw new Error('phase envelope could not be started; try a different Pstart');
  return {
    bubble: bub.map(strip), dew: dew.map(strip), points, critical,
    cricondenbar: refineExtremum(points, 'P'), cricondentherm: refineExtremum(points, 'T'),
    complete: critical != null,
  };
}

// ---------------------------------------------------------------------------------------------------------------
// Transport properties and auxiliary correlations

/**
 * Gas viscosity by Lee, Gonzalez & Eakin (1966); see `gasViscosityLee` in blackoil.js for the equation, reference
 * and validity.
 * @param {object} o
 * @param {number} o.T K @param {number} o.rho gas density, kg/m³ @param {number} o.M molar mass, kg/mol
 * @returns {number} Pa·s
 */
export function gasViscosityLGE({ T, rho, M }) {
  return gasViscosityLee({ T, rho, M });
}

/**
 * Dense-fluid viscosity by Lohrenz, Bray & Clark (1964):
 *   [(μ − μ*) ξ + 1e-4]^(1/4) = 0.1023 + 0.023364 ρr + 0.058533 ρr² − 0.040758 ρr³ + 0.0093324 ρr⁴,
 *   ξ = (Σ x_i Tc_i)^(1/6) / [(Σ x_i M_i)^(1/2) (Σ x_i Pc_i)^(2/3)]  (K, g/mol, atm; μ in cP),
 *   ρr = ρ_molar Σ x_i Vc_i,
 * with the dilute-gas mixture viscosity μ* by Herning & Zipperer, μ* = Σ x_i μ_i √M_i / Σ x_i √M_i, and the
 * component values by Stiel & Thodos: μ_i ξ_i = 34e-5 Tr^0.94 (Tr ≤ 1.5), 17.78e-5 (4.58 Tr − 1.67)^(5/8) (Tr > 1.5).
 * Applies to gas and liquid; the result is very sensitive to density (fourth power of a quartic in ρr), so use a
 * volume-translated density for liquids. Typical liquid accuracy without tuning Vc of the plus fraction: ±30 %.
 * References: Lohrenz, Bray & Clark, JPT 16(10) (1964) 1171–1176; Stiel & Thodos, AIChE J. 7 (1961) 611;
 * Herning & Zipperer (1936).
 * @param {object} o
 * @param {object} o.mix mixture model (all components need a critical volume)
 * @param {number} o.T K
 * @param {number[]} o.z phase composition, –
 * @param {number} o.rhoMolar molar density of the phase, mol/m³
 * @returns {number} viscosity, Pa·s
 */
export function viscosityLBC({ mix, T, z, rhoMolar }) {
  checkState(mix, T, 1e5, z);
  need('rhoMolar', rhoMolar, 0, 1e5);
  const n = mix.n;
  let Tpc = 0;
  let Ppc = 0;
  let M = 0;
  let Vpc = 0;
  let num = 0;
  let den = 0;
  for (let i = 0; i < n; i++) {
    if (z[i] <= 0) continue;
    const c = mix.comps[i];
    if (c.Vc == null) throw new RangeError(`component ${c.id} has no critical volume (needed by LBC)`);
    const Mg = c.M * 1e3;
    const Pa = c.Pc / 101325;
    Tpc += z[i] * c.Tc;
    Ppc += z[i] * Pa;
    M += z[i] * Mg;
    Vpc += z[i] * c.Vc;
    const xi = c.Tc ** (1 / 6) / (Math.sqrt(Mg) * Pa ** (2 / 3));
    const Tr = T / c.Tc;
    const mui = (Tr <= 1.5 ? 34e-5 * Tr ** 0.94 : 17.78e-5 * (4.58 * Tr - 1.67) ** 0.625) / xi;
    num += z[i] * mui * Math.sqrt(Mg);
    den += z[i] * Math.sqrt(Mg);
  }
  const muStar = num / den;
  const xi = Tpc ** (1 / 6) / (Math.sqrt(M) * Ppc ** (2 / 3));
  const r = rhoMolar * Vpc;
  const poly = 0.1023 + r * (0.023364 + r * (0.058533 + r * (-0.040758 + r * 0.0093324)));
  return (muStar + (poly ** 4 - 1e-4) / xi) * 1e-3;
}

/**
 * Thermal conductivity of a low-pressure gas from viscosity by the modified Eucken relation,
 * k M/(μ Cv) = 1.32 + 1.77 R/Cv, i.e. k = (μ/M)(1.32 Cv + 1.77 R).
 * Reference: Poling, Prausnitz & O'Connell (2001) ch. 10 (modified Eucken). Validity: dilute polyatomic gases,
 * typically within 10–15 %; dense-gas enhancement is not included (underestimates above ~10 MPa).
 * @param {object} o
 * @param {number} o.mu gas viscosity, Pa·s
 * @param {number} o.cv molar isochoric heat capacity, J/(mol·K)
 * @param {number} o.M molar mass, kg/mol
 * @returns {number} W/(m·K)
 */
export function gasThermalConductivity({ mu, cv, M }) {
  need('mu', mu, 1e-7, 1e-2);
  need('cv', cv, 1, 2000);
  need('M', M, 1e-3, 1);
  return (mu / M) * (1.32 * cv + 1.77 * R);
}

/**
 * Thermal conductivity of liquid petroleum by Cragoe (1929): k = (0.813/γ)[1 − 0.0003 (T°F − 32)]
 * Btu·in/(h·ft²·°F), i.e. k = (0.1172/γ)(1 − 0.00054 T°C) W/(m·K).
 * Reference: Cragoe, C.S., "Thermal Properties of Petroleum Products", NBS Misc. Publ. 97 (1929).
 * Validity: γ 0.78–0.95, 0–200 °C, about ±10 %.
 * @param {object} o
 * @param {number} o.T K @param {number} o.sg oil specific gravity at 60 °F, –
 * @returns {number} W/(m·K)
 */
export function oilThermalConductivity({ T, sg }) {
  need('T', T, 250, 500);
  need('sg', sg, 0.6, 1.1);
  return 0.813 * 0.144227889 / sg * (1 - 0.0003 * ((T * 1.8 - 459.67) - 32));
}

/**
 * Vapour–liquid interfacial tension by the Macleod–Sugden parachor method as extended to mixtures by Weinaug &
 * Katz: σ^(1/4) = Σ Pch_i (x_i ρ_L − y_i ρ_V), molar densities in mol/cm³, σ in dyn/cm.
 * Reference: Macleod (1923); Sugden (1924); Weinaug & Katz, Ind. Eng. Chem. 35 (1943) 239.
 * Validity: hydrocarbon systems (with N2, CO2, H2S); components with zero parachor data must be absent or given
 * through `parachors`. Accuracy ±10–20 %, poorer near the critical point.
 * @param {object} o
 * @param {object} o.mix mixture model
 * @param {number[]} o.x liquid mole fractions @param {number[]} o.y vapour mole fractions
 * @param {number} o.rhoMolarL liquid molar density, mol/m³ @param {number} o.rhoMolarV vapour molar density, mol/m³
 * @param {number[]} [o.parachors] override of component parachors
 * @returns {number} interfacial tension, N/m
 */
export function surfaceTension({ mix, x, y, rhoMolarL, rhoMolarV, parachors }) {
  need('rhoMolarL', rhoMolarL, 0, 1e5);
  need('rhoMolarV', rhoMolarV, 0, 1e5);
  let s = 0;
  for (let i = 0; i < mix.n; i++) {
    if (!(x[i] > 0) && !(y[i] > 0)) continue;
    const p = parachors ? parachors[i] : mix.comps[i].parachor;
    if (p == null) throw new RangeError(`component ${mix.ids[i]} has no parachor; supply parachors[]`);
    s += p * (x[i] * rhoMolarL - y[i] * rhoMolarV) * 1e-6;
  }
  return Math.max(s, 0) ** 4 * 1e-3;
}

/**
 * Equilibrium water content of sweet natural gas by Bukacek (1955):
 * W = 47 484 p_sat/P + B lb/MMscf, log10 B = −3083.87/T(°R) + 6.69449, with p_sat the vapour pressure of water
 * (IAPWS-IF97 here) and standard conditions 60 °F, 14.7 psia.
 * Reference: Bukacek, R.F., "Equilibrium moisture content of natural gases", IGT Research Bulletin 8 (1955).
 * Validity: 288–511 K (60–460 °F), 0.1–69 MPa, sweet gas over fresh water, above the hydrate temperature;
 * accuracy about ±5 %. No salinity, gas-gravity or acid-gas correction.
 * @param {object} o
 * @param {number} o.T K @param {number} o.P Pa
 * @returns {{ lbPerMMscf: number, kgPerSm3: number, y: number }} water content per standard volume of gas
 *   (60 °F, 101.325 kPa) and water mole fraction in the gas
 */
export function waterContentBukacek({ T, P }) {
  need('T', T, 273.15, 520);
  need('P', P, 1e4, 1e8);
  const psat = waterVapourPressure({ T }) / 6894.757293168;
  const p = P / 6894.757293168;
  const B = 10 ** (-3083.87 / (T * 1.8) + 6.69449);
  const W = 47484 * psat / p + B;
  return { lbPerMMscf: W, kgPerSm3: W * 0.45359237 / (1e6 * 0.028316846592), y: W / 47484 };
}

/** Standard reference conditions: ISO 13443 metric (15 °C, 101.325 kPa) and oilfield (60 °F, 14.696 psia). */
export const STANDARD = Object.freeze({
  metric: Object.freeze({ T: 288.15, P: 101325 }),
  field: Object.freeze({ T: 288.70555555555555, P: 101325 }),
  normal: Object.freeze({ T: 273.15, P: 101325 }),
});

/**
 * Ideal-gas molar volume at reference conditions, v = R T/P.
 * Reference: ideal-gas law (23.645 m³/kmol at 15 °C and 101.325 kPa; 379.5 scf/lbmol at 60 °F, 14.696 psia).
 * @param {{ T: number, P: number }} [ref=STANDARD.metric] reference T (K) and P (Pa)
 * @returns {number} m³/mol
 */
export function standardMolarVolume(ref = STANDARD.metric) {
  need('ref.T', ref.T, 200, 400);
  need('ref.P', ref.P, 5e4, 2e5);
  return R * ref.T / ref.P;
}

/**
 * Convert a gas rate at standard conditions to molar and mass rate (ideal gas at the reference state).
 * Reference: ideal-gas law. @param {object} o
 * @param {number} o.qStd standard volumetric rate, Sm³/s
 * @param {number} o.M molar mass, kg/mol
 * @param {{ T: number, P: number }} [o.ref=STANDARD.metric]
 * @returns {{ molar: number, mass: number }} mol/s and kg/s
 */
export function standardRateToMass({ qStd, M, ref = STANDARD.metric }) {
  need('qStd', qStd, 0, 1e9);
  need('M', M, 1e-3, 1);
  const molar = qStd / standardMolarVolume(ref);
  return { molar, mass: molar * M };
}

/**
 * Convert a mass rate to a gas rate at standard conditions (inverse of `standardRateToMass`).
 * @param {object} o
 * @param {number} o.mass kg/s @param {number} o.M kg/mol @param {{ T: number, P: number }} [o.ref=STANDARD.metric]
 * @returns {number} Sm³/s
 */
export function massToStandardRate({ mass, M, ref = STANDARD.metric }) {
  need('mass', mass, 0, 1e9);
  need('M', M, 1e-3, 1);
  return mass / M * standardMolarVolume(ref);
}

/**
 * Convert million standard cubic feet per day to m³/s at the same (60 °F, 14.696 psia) reference state.
 * Reference: exact unit conversion (1 ft³ = 0.028316846592 m³).
 * @param {number} mmscfd @returns {number} Sm³/s
 */
export function mmscfdToSm3s(mmscfd) {
  need('mmscfd', mmscfd, 0, 1e9);
  return mmscfd * 1e6 * 0.028316846592 / 86400;
}
