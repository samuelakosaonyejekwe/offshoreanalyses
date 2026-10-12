/**
 * Nonlinear material models for point (integration-point) use.
 *
 *  - Rate-independent J2 (von Mises) plasticity with the radial return algorithm in one and three dimensions:
 *    linear and Voce (saturation) isotropic hardening, linear kinematic hardening (Prager), consistent tangent.
 *    Reference: Simo & Hughes, "Computational Inelasticity" (1998), Box 1.5 (1-D) and Box 3.1/3.2 (3-D).
 *  - Ramberg–Osgood stress–strain curve and its inverse.
 *  - Johnson–Cook flow stress, Cowper–Symonds rate factor.
 *  - True ↔ engineering stress and strain.
 *  - Incompressible Neo-Hookean and Mooney–Rivlin solids in uniaxial tension / compression.
 *  - Mander model for confined concrete.
 *
 * SI units (Pa, K). 3-D stresses and strains are Voigt vectors [xx, yy, zz, xy, yz, zx] with engineering shear
 * strains (γ = 2ε) in total strains; plastic strain and back stress are stored as tensor components.
 * Small strains are assumed in the plasticity models.
 */
import { brent } from '../math/roots.js';

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}
function fin(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}

/* ------------------------------------------------------------------ hardening law */

/**
 * Isotropic hardening law σy(α) = σy0 + H α + (σ∞ − σy0)(1 − e^(−δ α)) — linear plus Voce saturation term
 * (Voce 1948; the saturation law used by Simo & Hughes ch. 3) — and its slope dσy/dα. α is the accumulated equivalent plastic strain.
 * @param {object} o
 * @param {number} o.alpha equivalent plastic strain [-] (>= 0)
 * @param {number} o.sy0 initial yield stress [Pa]
 * @param {number} [o.H=0] linear isotropic hardening modulus [Pa]
 * @param {{sInf:number, delta:number}} [o.voce] saturation stress σ∞ [Pa] and exponent δ [-]
 * @returns {{sy:number, slope:number}} yield stress [Pa] and dσy/dα [Pa]
 */
export function hardeningStress({ alpha, sy0, H = 0, voce }) {
  pos('sy0', sy0);
  if (!(alpha >= 0)) throw new RangeError('alpha must be non-negative');
  let sy = sy0 + H * alpha; let slope = H;
  if (voce) {
    const e = Math.exp(-voce.delta * alpha);
    sy += (voce.sInf - sy0) * (1 - e); slope += (voce.sInf - sy0) * voce.delta * e;
  }
  return { sy, slope };
}

/* Solve r(Δγ) = trial − c·Δγ − s·σy(α + s·Δγ) = 0 for Δγ ≥ 0 (c: elastic + kinematic modulus term, s: 1 or √(2/3)). */
function solveGamma(trial, c, s, alpha, p) {
  if (!p.voce) return (trial - s * (p.sy0 + p.H * alpha)) / (c + s * s * p.H);
  let dg = 0;
  for (let it = 0; it < 60; it++) {
    const h = hardeningStress({ alpha: alpha + s * dg, sy0: p.sy0, H: p.H, voce: p.voce });
    const r = trial - c * dg - s * h.sy;
    const step = r / (c + s * s * h.slope);
    dg += step;
    if (Math.abs(step) <= 1e-14 * Math.max(1e-12, Math.abs(dg))) return dg;
  }
  // softening laws with a non-monotone residual: bracketed fallback
  return brent((x) => trial - c * x - s * hardeningStress({ alpha: alpha + s * x, sy0: p.sy0, H: p.H, voce: p.voce }).sy, 0, trial / c, { tol: 1e-15 });
}

function checkParams(p) {
  pos('E', p.E); pos('sy0', p.sy0);
  const H = p.H === undefined ? 0 : fin('H', p.H); const Hkin = p.Hkin === undefined ? 0 : nonNeg('Hkin', p.Hkin);
  if (p.voce && !(p.voce.delta > 0 && Number.isFinite(p.voce.sInf))) throw new RangeError('voce must be { sInf, delta } with delta > 0');
  return { ...p, H, Hkin };
}

/* ------------------------------------------------------------------ 1-D plasticity */

/**
 * One-dimensional rate-independent plasticity with combined isotropic and kinematic hardening: return mapping for a
 * total strain (Simo & Hughes Box 1.5). Yield function f = |σ − q| − σy(α); flow ε̇p = γ sign(σ − q), α̇ = γ,
 * q̇ = Hkin ε̇p. Consistent tangent E (σy' + Hkin)/(E + σy' + Hkin) when yielding.
 * @param {object} o
 * @param {number} o.strain total strain at the end of the step [-]
 * @param {{epsP:number, alpha:number, back:number}} [o.state] state at the start of the step (plastic strain, accumulated plastic strain,
 *   back stress [Pa]); default virgin material
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.sy0 initial yield stress [Pa]
 * @param {number} [o.H=0] isotropic hardening modulus [Pa]
 * @param {number} [o.Hkin=0] kinematic hardening modulus [Pa]
 * @param {{sInf:number, delta:number}} [o.voce] Voce saturation hardening
 * @returns {{stress:number, state:{epsP:number, alpha:number, back:number}, tangent:number, yielded:boolean, dGamma:number}} stress [Pa], updated
 *   state, algorithmic tangent dσ/dε [Pa], plastic multiplier increment [-]
 */
export function plasticity1d(o) {
  const p = checkParams(o);
  fin('strain', o.strain);
  const st = o.state || { epsP: 0, alpha: 0, back: 0 };
  const sTr = p.E * (o.strain - st.epsP);
  const xi = sTr - st.back;
  const sy = hardeningStress({ alpha: st.alpha, sy0: p.sy0, H: p.H, voce: p.voce }).sy;
  if (Math.abs(xi) - sy <= 0) return { stress: sTr, state: { epsP: st.epsP, alpha: st.alpha, back: st.back }, tangent: p.E, yielded: false, dGamma: 0 };
  const dg = solveGamma(Math.abs(xi), p.E + p.Hkin, 1, st.alpha, p);
  const sg = Math.sign(xi);
  const alpha = st.alpha + dg;
  const slope = hardeningStress({ alpha, sy0: p.sy0, H: p.H, voce: p.voce }).slope;
  return {
    stress: sTr - p.E * dg * sg, state: { epsP: st.epsP + dg * sg, alpha, back: st.back + p.Hkin * dg * sg },
    tangent: p.E * (slope + p.Hkin) / (p.E + slope + p.Hkin), yielded: true, dGamma: dg,
  };
}

/**
 * Stress history of the 1-D plasticity model for a prescribed strain history (monotonic or cyclic).
 * @param {object} o material parameters of plasticity1d (E, sy0, H, Hkin, voce) plus
 * @param {ArrayLike<number>} o.strains total strain history [-]
 * @param {object} [o.state] initial state
 * @returns {{stress:Float64Array, epsP:Float64Array, tangent:Float64Array, state:object, dissipation:number}} stresses [Pa], plastic strains,
 *   tangents [Pa], final state, plastic work ∫σ dεp per unit volume [J/m³]
 */
export function stressStrainHistory1d(o) {
  const { strains } = o;
  if (!strains || !(strains.length > 0)) throw new RangeError('strains must be a non-empty array');
  const n = strains.length;
  const stress = new Float64Array(n); const epsP = new Float64Array(n); const tangent = new Float64Array(n);
  let state = o.state || { epsP: 0, alpha: 0, back: 0 };
  let work = 0; let sPrev = 0;
  for (let i = 0; i < n; i++) {
    const r = plasticity1d({ ...o, strain: strains[i], state });
    work += 0.5 * (r.stress + sPrev) * (r.state.epsP - state.epsP);
    stress[i] = r.stress; epsP[i] = r.state.epsP; tangent[i] = r.tangent; state = r.state; sPrev = r.stress;
  }
  return { stress, epsP, tangent, state, dissipation: work };
}

/* ------------------------------------------------------------------ 3-D J2 plasticity */

/**
 * von Mises equivalent stress of a Voigt stress vector [xx, yy, zz, xy, yz, zx], optionally relative to a back stress.
 * @param {ArrayLike<number>} s stress [Pa]
 * @param {ArrayLike<number>} [back] deviatoric back stress [Pa]
 * @returns {number} σeq = √(3 J2) [Pa]
 */
export function vonMisesStress(s, back) {
  if (!s || s.length !== 6) throw new RangeError('s must have 6 components');
  const m = (s[0] + s[1] + s[2]) / 3;
  let j2 = 0;
  for (let i = 0; i < 6; i++) {
    const d = s[i] - (i < 3 ? m : 0) - (back ? back[i] : 0);
    j2 += (i < 3 ? 0.5 : 1) * d * d;
  }
  return Math.sqrt(3 * j2);
}

/**
 * Three-dimensional J2 plasticity with isotropic (linear + Voce) and linear kinematic hardening, integrated by the
 * radial return algorithm (Wilkins 1964; Simo & Hughes Box 3.1) with the consistent elastoplastic tangent (Box 3.2):
 *   ξ = s_trial − β, f = ‖ξ‖ − √(2/3) σy(α);  ‖ξ‖ − (2G + ⅔Hkin)Δγ − √(2/3) σy(α + √(2/3)Δγ) = 0;
 *   s = s_trial − 2GΔγ n, e_p += Δγ n, β += ⅔ Hkin Δγ n, α += √(2/3) Δγ, n = ξ/‖ξ‖;
 *   C = K 1⊗1 + 2Gθ (I − ⅓ 1⊗1) − 2Gθ̄ n⊗n, θ = 1 − 2GΔγ/‖ξ‖, θ̄ = 1/(1 + (σy' + Hkin)/(3G)) − (1 − θ).
 * @param {object} o
 * @param {ArrayLike<number>} o.strain total strain [εxx, εyy, εzz, γxy, γyz, γzx] (engineering shear) [-]
 * @param {{epsP:number[], alpha:number, back:number[]}} [o.state] state at the start of the step: plastic strain tensor components
 *   [xx, yy, zz, xy, yz, zx] (tensor shear), equivalent plastic strain, back stress [Pa]; default virgin material
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} o.sy0 initial yield stress [Pa]
 * @param {number} [o.H=0] isotropic hardening modulus dσy/dα [Pa]
 * @param {number} [o.Hkin=0] kinematic hardening modulus [Pa] (uniaxial plastic modulus contribution)
 * @param {{sInf:number, delta:number}} [o.voce] Voce saturation hardening
 * @param {boolean} [o.tangent=true] also return the 6×6 consistent tangent
 * @returns {{stress:number[], state:{epsP:number[], alpha:number, back:number[]}, tangent:(number[][]|undefined), yielded:boolean, dGamma:number,
 *   vonMises:number}} stress [Pa] (Voigt), updated state, tangent ∂σ/∂ε for engineering shear strains [Pa], plastic multiplier, σeq [Pa]
 */
export function plasticity3d(o) {
  const p = checkParams(o);
  if (!(p.nu > -1 && p.nu < 0.5)) throw new RangeError('nu must lie in (-1, 0.5)');
  const eps = o.strain;
  if (!eps || eps.length !== 6) throw new RangeError('strain must have 6 components');
  const st = o.state || { epsP: [0, 0, 0, 0, 0, 0], alpha: 0, back: [0, 0, 0, 0, 0, 0] };
  const G = p.E / (2 * (1 + p.nu)); const K = p.E / (3 * (1 - 2 * p.nu));
  const tr = eps[0] + eps[1] + eps[2];
  const xi = new Array(6); const sTr = new Array(6);
  let norm2 = 0;
  for (let i = 0; i < 6; i++) {
    const e = i < 3 ? eps[i] - tr / 3 : 0.5 * eps[i];
    sTr[i] = 2 * G * (e - st.epsP[i]);
    xi[i] = sTr[i] - st.back[i];
    norm2 += (i < 3 ? 1 : 2) * xi[i] * xi[i];
  }
  const norm = Math.sqrt(norm2);
  const r23 = Math.sqrt(2 / 3);
  const sy = hardeningStress({ alpha: st.alpha, sy0: p.sy0, H: p.H, voce: p.voce }).sy;
  const wantC = o.tangent !== false;
  const elasticC = () => {
    const C = [];
    for (let i = 0; i < 6; i++) { C.push([0, 0, 0, 0, 0, 0]); }
    for (let i = 0; i < 3; i++) { for (let j = 0; j < 3; j++) C[i][j] = K - 2 * G / 3 + (i === j ? 2 * G : 0); C[i + 3][i + 3] = G; }
    return C;
  };
  if (norm - r23 * sy <= 0) {
    const stress = sTr.map((s, i) => s + (i < 3 ? K * tr : 0));
    return {
      stress, state: { epsP: st.epsP.slice(), alpha: st.alpha, back: st.back.slice() }, tangent: wantC ? elasticC() : undefined, yielded: false,
      dGamma: 0, vonMises: vonMisesStress(stress),
    };
  }
  const dg = solveGamma(norm, 2 * G + 2 / 3 * p.Hkin, r23, st.alpha, p);
  const alpha = st.alpha + r23 * dg;
  const n = xi.map((v) => v / norm);
  const stress = new Array(6); const epsP = new Array(6); const back = new Array(6);
  for (let i = 0; i < 6; i++) {
    stress[i] = sTr[i] - 2 * G * dg * n[i] + (i < 3 ? K * tr : 0);
    epsP[i] = st.epsP[i] + dg * n[i];
    back[i] = st.back[i] + 2 / 3 * p.Hkin * dg * n[i];
  }
  let C;
  if (wantC) {
    const slope = hardeningStress({ alpha, sy0: p.sy0, H: p.H, voce: p.voce }).slope;
    const theta = 1 - 2 * G * dg / norm;
    const thetaBar = 1 / (1 + (slope + p.Hkin) / (3 * G)) - (1 - theta);
    C = [];
    for (let i = 0; i < 6; i++) {
      const row = [0, 0, 0, 0, 0, 0];
      for (let j = 0; j < 6; j++) {
        let v = -2 * G * thetaBar * n[i] * n[j];
        if (i < 3 && j < 3) v += K - 2 * G * theta / 3 + (i === j ? 2 * G * theta : 0);
        else if (i === j) v += G * theta;
        row[j] = v;
      }
      C.push(row);
    }
  }
  return { stress, state: { epsP, alpha, back }, tangent: C, yielded: true, dGamma: dg, vonMises: vonMisesStress(stress) };
}

/**
 * Uniaxial-stress response of the 3-D J2 model: for each axial strain the lateral strains are iterated (Newton with the
 * consistent tangent) until the lateral stresses vanish. Used to verify the 3-D return map against the 1-D model
 * and to generate stress–strain curves with the 3-D state.
 * @param {object} o material parameters of plasticity3d plus
 * @param {ArrayLike<number>} o.strains axial strain history [-]
 * @returns {{stress:Float64Array, lateralStrain:Float64Array, alpha:Float64Array, state:object}} axial stress [Pa], lateral strain [-],
 *   equivalent plastic strain [-] at each step; final state
 */
export function uniaxialResponse3d(o) {
  const { strains } = o;
  if (!strains || !(strains.length > 0)) throw new RangeError('strains must be a non-empty array');
  const n = strains.length;
  const stress = new Float64Array(n); const lat = new Float64Array(n); const alpha = new Float64Array(n);
  let state; let eLat = 0;
  for (let i = 0; i < n; i++) {
    let r;
    for (let it = 0; it < 50; it++) {
      r = plasticity3d({ ...o, strain: [strains[i], eLat, eLat, 0, 0, 0], state });
      // lateral stress σyy = 0 with εyy = εzz: dσyy/dεlat = C[1][1] + C[1][2]
      const res = r.stress[1];
      if (Math.abs(res) <= 1e-12 * o.sy0) break;
      eLat -= res / (r.tangent[1][1] + r.tangent[1][2]);
    }
    state = r.state;
    stress[i] = r.stress[0]; lat[i] = eLat; alpha[i] = state.alpha;
  }
  return { stress, lateralStrain: lat, alpha, state };
}

/* ------------------------------------------------------------------ Ramberg–Osgood */

/**
 * Ramberg–Osgood stress–strain curve in the offset-yield form (Ramberg & Osgood 1943, NACA TN 902):
 * ε = σ/E + ε0 (σ/σy)^n, where σy is the proof stress at the plastic offset ε0 (0.2 % by default).
 * @param {object} o
 * @param {number} o.stress stress σ [Pa] (sign preserved)
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.sy proof (yield) stress at the offset [Pa]
 * @param {number} o.n hardening exponent [-] (>= 1)
 * @param {number} [o.offset=0.002] plastic strain at σy [-]
 * @returns {{strain:number, plasticStrain:number, tangent:number}} total strain, plastic strain [-], tangent modulus dσ/dε [Pa]
 */
export function rambergOsgoodStrain({ stress, E, sy, n, offset = 0.002 }) {
  fin('stress', stress); pos('E', E); pos('sy', sy); pos('offset', offset);
  if (!(n >= 1)) throw new RangeError('n must be >= 1');
  const a = Math.abs(stress);
  const ep = offset * (a / sy) ** n;
  const sg = stress < 0 ? -1 : 1;
  return { strain: sg * (a / E + ep), plasticStrain: sg * ep, tangent: 1 / (1 / E + (a > 0 ? n * ep / a : (n === 1 ? offset / sy : 0))) };
}

/**
 * Inverse of the Ramberg–Osgood curve: stress for a given total strain (Newton iteration, monotone and convex).
 * @param {object} o
 * @param {number} o.strain total strain [-]
 * @param {number} o.E [Pa]
 * @param {number} o.sy proof stress [Pa]
 * @param {number} o.n exponent [-]
 * @param {number} [o.offset=0.002]
 * @returns {{stress:number, plasticStrain:number, tangent:number}} stress [Pa], plastic strain [-], tangent modulus [Pa]
 */
export function rambergOsgoodStress({ strain, E, sy, n, offset = 0.002 }) {
  fin('strain', strain); pos('E', E); pos('sy', sy); pos('offset', offset);
  if (!(n >= 1)) throw new RangeError('n must be >= 1');
  const e = Math.abs(strain);
  if (e === 0) return { stress: 0, plasticStrain: 0, tangent: rambergOsgoodStrain({ stress: 0, E, sy, n, offset }).tangent };
  // start from the smaller of the elastic and the pure power-law estimate (both above the root → Newton descends monotonically)
  let s = Math.min(E * e, sy * (e / offset) ** (1 / n));
  for (let it = 0; it < 100; it++) {
    const f = s / E + offset * (s / sy) ** n - e;
    const df = 1 / E + n * offset * (s / sy) ** n / s;
    const ds = f / df;
    s -= ds;
    if (Math.abs(ds) <= 1e-14 * s) break;
  }
  const r = rambergOsgoodStrain({ stress: s, E, sy, n, offset });
  const sg = strain < 0 ? -1 : 1;
  return { stress: sg * s, plasticStrain: sg * Math.abs(r.plasticStrain), tangent: r.tangent };
}

/**
 * Ramberg–Osgood exponent from two points of the curve: the proof stress σy at plastic strain ε0 and the ultimate
 * stress σu at plastic strain εpu: n = ln(εpu/ε0)/ln(σu/σy).
 * @param {object} o
 * @param {number} o.sy proof stress [Pa]
 * @param {number} o.su ultimate (second point) stress [Pa]
 * @param {number} o.plasticStrainU plastic strain at σu [-]
 * @param {number} [o.offset=0.002]
 * @returns {number} n [-]
 */
export function rambergOsgoodExponent({ sy, su, plasticStrainU, offset = 0.002 }) {
  pos('sy', sy); pos('su', su); pos('plasticStrainU', plasticStrainU); pos('offset', offset);
  if (!(su > sy) || !(plasticStrainU > offset)) throw new RangeError('su must exceed sy and plasticStrainU must exceed offset');
  return Math.log(plasticStrainU / offset) / Math.log(su / sy);
}

/* ------------------------------------------------------------------ rate and temperature dependence */

/**
 * Johnson–Cook flow stress (Johnson & Cook 1983, 7th Int. Symp. on Ballistics):
 * σ = (A + B εp^n)(1 + C ln ε̇*)(1 − T*^m), ε̇* = ε̇/ε̇0, T* = (T − Tr)/(Tm − Tr).
 * The constants A, B, n, C, m, ε̇0, Tr, Tm are material data and must be supplied. For ε̇ < ε̇0 the rate factor is taken
 * as 1 and for T ≤ Tr the thermal factor as 1 (common implementation convention); the stress is zero at and above Tm.
 * @param {object} o
 * @param {number} o.epsP equivalent plastic strain [-]
 * @param {number} [o.strainRate] equivalent plastic strain rate [1/s] (omit for the quasi-static curve)
 * @param {number} [o.T] temperature [K] (omit for the reference temperature)
 * @param {number} o.A yield stress constant [Pa]
 * @param {number} o.B hardening constant [Pa]
 * @param {number} o.n hardening exponent [-]
 * @param {number} [o.C=0] strain-rate constant [-]
 * @param {number} [o.m=1] thermal softening exponent [-]
 * @param {number} [o.refRate=1] reference strain rate ε̇0 [1/s]
 * @param {number} [o.Tref=293.15] reference temperature [K]
 * @param {number} [o.Tmelt] melting temperature [K] (required when T is given)
 * @returns {{stress:number, hardening:number, rateFactor:number, thermalFactor:number}} flow stress [Pa] and its three factors
 */
export function johnsonCook({ epsP, strainRate, T, A, B, n, C = 0, m = 1, refRate = 1, Tref = 293.15, Tmelt }) {
  nonNeg('epsP', epsP); pos('A', A); nonNeg('B', B); pos('n', n); nonNeg('C', C); pos('m', m); pos('refRate', refRate);
  const hard = A + B * epsP ** n;
  let rate = 1;
  if (strainRate !== undefined) { nonNeg('strainRate', strainRate); if (strainRate > refRate) rate = 1 + C * Math.log(strainRate / refRate); }
  let th = 1;
  if (T !== undefined) {
    fin('T', T); pos('Tmelt', Tmelt);
    if (!(Tmelt > Tref)) throw new RangeError('Tmelt must exceed Tref');
    if (T >= Tmelt) th = 0; else if (T > Tref) th = 1 - ((T - Tref) / (Tmelt - Tref)) ** m;
  }
  return { stress: hard * rate * th, hardening: hard, rateFactor: rate, thermalFactor: th };
}

/**
 * Cowper–Symonds strain-rate factor on the flow stress: σd/σs = 1 + (ε̇/D)^(1/q) (Cowper & Symonds 1957). D and q are
 * material constants to be supplied (see blast.dynamicIncreaseFactor for commonly quoted values).
 * @param {object} o
 * @param {number} o.strainRate [1/s]
 * @param {number} o.D [1/s]
 * @param {number} o.q [-]
 * @returns {number} factor [-]
 */
export function cowperSymonds({ strainRate, D, q }) {
  nonNeg('strainRate', strainRate); pos('D', D); pos('q', q);
  return 1 + (strainRate / D) ** (1 / q);
}

/* ------------------------------------------------------------------ stress and strain measures */

/**
 * True (Cauchy) stress and logarithmic strain from engineering values for uniform deformation at constant volume:
 * σt = σe (1 + εe), εt = ln(1 + εe). Valid up to the onset of necking (uniform elongation).
 * @param {object} o
 * @param {number} o.stress engineering stress [Pa]
 * @param {number} o.strain engineering strain [-] (> −1)
 * @returns {{stress:number, strain:number}} true stress [Pa], true strain [-]
 */
export function trueFromEngineering({ stress, strain }) {
  fin('stress', stress);
  if (!(strain > -1)) throw new RangeError('strain must exceed -1');
  return { stress: stress * (1 + strain), strain: Math.log(1 + strain) };
}

/**
 * Engineering stress and strain from true values: εe = e^εt − 1, σe = σt e^(−εt) (inverse of trueFromEngineering).
 * @param {object} o
 * @param {number} o.stress true stress [Pa]
 * @param {number} o.strain true (logarithmic) strain [-]
 * @returns {{stress:number, strain:number}} engineering stress [Pa], engineering strain [-]
 */
export function engineeringFromTrue({ stress, strain }) {
  fin('stress', stress); fin('strain', strain);
  return { stress: stress * Math.exp(-strain), strain: Math.expm1(strain) };
}

/* ------------------------------------------------------------------ hyperelasticity */

/**
 * Incompressible Mooney–Rivlin solid in uniaxial tension or compression, W = C10 (I1 − 3) + C01 (I2 − 3)
 * (Rivlin 1948; Treloar, The Physics of Rubber Elasticity): Cauchy stress σ = 2 (λ² − 1/λ)(C10 + C01/λ), nominal
 * stress P = σ/λ. The small-strain Young's modulus is 6 (C10 + C01).
 * @param {object} o
 * @param {number} o.stretch λ = L/L0 [-] (> 0)
 * @param {number} o.C10 [Pa]
 * @param {number} [o.C01=0] [Pa] (0 gives the Neo-Hookean solid)
 * @returns {{cauchy:number, nominal:number, tangent:number, energy:number}} Cauchy stress [Pa], nominal (engineering) stress [Pa], dP/dλ [Pa],
 *   strain energy density [J/m³]
 */
export function mooneyRivlinUniaxial({ stretch, C10, C01 = 0 }) {
  pos('stretch', stretch); fin('C10', C10); fin('C01', C01);
  const l = stretch;
  const P = 2 * (l - 1 / (l * l)) * (C10 + C01 / l);
  const dP = 2 * (1 + 2 / l ** 3) * (C10 + C01 / l) - 2 * (l - 1 / (l * l)) * C01 / (l * l);
  const I1 = l * l + 2 / l; const I2 = 2 * l + 1 / (l * l);
  return { cauchy: P * l, nominal: P, tangent: dP, energy: C10 * (I1 - 3) + C01 * (I2 - 3) };
}

/**
 * Incompressible Neo-Hookean solid in uniaxial tension or compression, W = (μ/2)(I1 − 3): σ = μ (λ² − 1/λ),
 * P = μ (λ − 1/λ²); small-strain Young's modulus 3μ.
 * @param {object} o
 * @param {number} o.stretch λ [-]
 * @param {number} o.mu shear modulus [Pa]
 * @returns {{cauchy:number, nominal:number, tangent:number, energy:number}} as mooneyRivlinUniaxial
 */
export function neoHookeanUniaxial({ stretch, mu }) {
  pos('mu', mu);
  return mooneyRivlinUniaxial({ stretch, C10: mu / 2, C01: 0 });
}

/* ------------------------------------------------------------------ confined concrete */

/**
 * Effective lateral confining pressure of circular hoops or spirals (Mander, Priestley & Park 1988, J. Struct. Eng.
 * 114(8)): f'l = ½ ke ρs fyh, ρs = 4 Asp/(ds s), ke = (1 − s'/(2 ds))²/(1 − ρcc) for hoops and
 * (1 − s'/(2 ds))/(1 − ρcc) for spirals.
 * @param {object} o
 * @param {number} o.ds diameter of the hoop / spiral centreline [m]
 * @param {number} o.spacing centre-to-centre spacing or pitch s [m]
 * @param {number} [o.clearSpacing] clear spacing s' [m] (default s minus the bar diameter derived from Asp)
 * @param {number} o.Asp area of the hoop / spiral bar [m²]
 * @param {number} o.fyh yield strength of the transverse steel [Pa]
 * @param {number} [o.rhoCc=0] ratio of longitudinal steel area to core area [-]
 * @param {'hoops'|'spiral'} [o.type='hoops']
 * @returns {{fl:number, ke:number, rhoS:number}} effective confining pressure [Pa], confinement effectiveness [-], volumetric ratio [-]
 */
export function manderConfinement({ ds, spacing, clearSpacing, Asp, fyh, rhoCc = 0, type = 'hoops' }) {
  pos('ds', ds); pos('spacing', spacing); pos('Asp', Asp); pos('fyh', fyh);
  if (!(rhoCc >= 0 && rhoCc < 1)) throw new RangeError('rhoCc must lie in [0, 1)');
  if (type !== 'hoops' && type !== 'spiral') throw new RangeError('type must be "hoops" or "spiral"');
  const sc = clearSpacing === undefined ? spacing - Math.sqrt(4 * Asp / Math.PI) : nonNeg('clearSpacing', clearSpacing);
  const a = Math.max(0, 1 - sc / (2 * ds));
  const ke = (type === 'hoops' ? a * a : a) / (1 - rhoCc);
  const rhoS = 4 * Asp / (ds * spacing);
  return { fl: 0.5 * ke * rhoS * fyh, ke, rhoS };
}

/**
 * Stress–strain curve of confined concrete, Mander, Priestley & Park (1988):
 *   f'cc = f'co [−1.254 + 2.254 √(1 + 7.94 f'l/f'co) − 2 f'l/f'co],   εcc = εco [1 + 5 (f'cc/f'co − 1)],
 *   fc = f'cc x r/(r − 1 + x^r),  x = εc/εcc,  r = Ec/(Ec − Esec),  Esec = f'cc/εcc,  Ec = 5000 √f'co (MPa) by default.
 * f'l = 0 gives the unconfined curve (Popovics form).
 * @param {object} o
 * @param {number} o.strain compressive strain εc [-] (>= 0)
 * @param {number} o.fco unconfined compressive strength f'co [Pa]
 * @param {number} [o.fl=0] effective lateral confining pressure f'l [Pa] (see manderConfinement)
 * @param {number} [o.epsCo=0.002] strain at the unconfined strength [-]
 * @param {number} [o.Ec] tangent modulus of elasticity [Pa]
 * @returns {{stress:number, fcc:number, epsCc:number, Esec:number, r:number, Ec:number}} stress [Pa], confined strength [Pa], strain at the
 *   confined strength [-], secant modulus at the peak [Pa], curve parameter [-], elastic modulus [Pa]
 */
export function manderConfinedConcrete({ strain, fco, fl = 0, epsCo = 0.002, Ec }) {
  nonNeg('strain', strain); pos('fco', fco); nonNeg('fl', fl); pos('epsCo', epsCo);
  const E = Ec === undefined ? 5000e6 * Math.sqrt(fco / 1e6) : pos('Ec', Ec);
  const fcc = fco * (-1.254 + 2.254 * Math.sqrt(1 + 7.94 * fl / fco) - 2 * fl / fco);
  const epsCc = epsCo * (1 + 5 * (fcc / fco - 1));
  const Esec = fcc / epsCc;
  if (!(E > Esec)) throw new RangeError('Ec must exceed the secant modulus f\'cc/εcc');
  const r = E / (E - Esec);
  const x = strain / epsCc;
  return { stress: x === 0 ? 0 : fcc * x * r / (r - 1 + x ** r), fcc, epsCc, Esec, r, Ec: E };
}
