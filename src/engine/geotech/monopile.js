// Large-diameter monopile helpers: foundation stiffness, tower natural frequency, serviceability, sizing.
// Foundation stiffness terms follow pile.js: [H; M] = [[KLL, KLR], [KLR, KRR]] · [y; theta] with KLR negative.

import { eigGenSym } from '../math/eigen.js';
import { pileHeadStiffness, solveLateralPile } from './pile.js';

function pos(name, v, allowZero = false) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < 0 || (!allowZero && v === 0)) {
    throw new RangeError(`${name}: must be a finite number ${allowZero ? '>= 0' : '> 0'} (got ${v})`);
  }
}

/**
 * Rigid pile on Winkler springs with modulus k(z) = kRef (z/zRef)^n [N/m²] (n = 0 constant, 1 linear,
 * 0.5 parabolic), by direct integration over the embedded length:
 *   KLL = ∫k dz, KLR = −∫k z dz, KRR = ∫k z² dz.
 * For n = 0 (k = kh D) and n = 1 (k = nh z) these are the rigid-pile expressions tabulated after Poulos & Davis
 * (1980): kh D L, −kh D L²/2, kh D L³/3 and nh L²/2, −nh L³/3, nh L⁴/4.
 * Reference: first principles; Poulos & Davis (1980) Pile Foundation Analysis and Design, ch. 8.
 * Validity: pile much stiffer than the soil (rigid-body rotation); base shear and base moment neglected.
 * @param {{ L: number, kRef: number, n?: number, zRef?: number }} o L [m] embedded length; kRef [N/m²] spring
 *   modulus at depth zRef; n [-] exponent >= 0 (default 0); zRef [m] (default L)
 * @returns {{ KLL: number, KLR: number, KRR: number }} [N/m], [N], [N·m/rad]
 */
export function rigidPileStiffness({ L, kRef, n = 0, zRef }) {
  pos('L', L); pos('kRef', kRef); pos('n', n, true);
  const zr = zRef ?? L;
  pos('zRef', zr);
  const c = kRef / Math.pow(zr, n);
  return { KLL: (c * L ** (n + 1)) / (n + 1), KLR: (-c * L ** (n + 2)) / (n + 2), KRR: (c * L ** (n + 3)) / (n + 3) };
}

/**
 * Long flexible pile on Winkler springs.
 * 'constant' modulus k [N/m²]: inverse of the Hetényi end-loaded semi-infinite beam, KLL = 4EIβ³ = k/β,
 *   KLR = −2EIβ² = −k/(2β²), KRR = 2EIβ = k/(2β³), β = (k/4EI)^(1/4) (exact).
 * 'linear' modulus k = nh z: KLL = 1.074 nh^(3/5) EI^(2/5), KLR = −0.99 nh^(2/5) EI^(3/5),
 *   KRR = 1.48 nh^(1/5) EI^(4/5) (coefficients as published; inverting the Reese & Matlock (1956) deflection
 *   and slope coefficients 2.435, 1.623, 1.750 gives 1.076, −0.998, 1.497).
 * Reference: Hetényi (1946); Poulos & Davis (1980) following Barber (1953), as tabulated by Arany et al. (2017)
 * Soil Dyn. Earthq. Eng. 92. Validity: slender (flexible) piles, L greater than about 4 (EI/nh)^(1/5)
 * or βL > ~5.
 * @param {{ EI: number, modulus: string, k?: number, nh?: number }} o EI [N·m²]; modulus 'constant'|'linear';
 *   k [N/m²] for 'constant'; nh [N/m³] for 'linear'
 * @returns {{ KLL: number, KLR: number, KRR: number }} [N/m], [N], [N·m/rad]
 */
export function flexiblePileStiffness({ EI, modulus, k, nh }) {
  pos('EI', EI);
  if (modulus === 'constant') {
    pos('k', k);
    const b = Math.pow(k / (4 * EI), 0.25);
    return { KLL: 4 * EI * b ** 3, KLR: -2 * EI * b * b, KRR: 2 * EI * b };
  }
  if (modulus === 'linear') {
    pos('nh', nh);
    return { KLL: 1.074 * nh ** 0.6 * EI ** 0.4, KLR: -0.99 * nh ** 0.4 * EI ** 0.6, KRR: 1.48 * nh ** 0.2 * EI ** 0.8 };
  }
  throw new RangeError("modulus: must be 'constant' or 'linear'");
}

const EC8 = {
  homogeneous: [[1.08, 0.21], [0.22, 0.5], [0.16, 0.75]],
  linear: [[0.6, 0.35], [0.17, 0.6], [0.14, 0.8]],
  parabolic: [[0.79, 0.28], [0.24, 0.53], [0.15, 0.77]],
};

/**
 * Static head stiffness of a flexible pile in an elastic continuum for three soil stiffness profiles
 * (E = Es homogeneous; E = Es z/D linear; E = Es √(z/D) parabolic; Es is the soil Young's modulus at depth D):
 *   KLL = a D Es (Ep/Es)^p, KLR = −b D² Es (Ep/Es)^q, KRR = c D³ Es (Ep/Es)^r with
 *   homogeneous 1.08/0.21, 0.22/0.50, 0.16/0.75; linear 0.60/0.35, 0.17/0.60, 0.14/0.80;
 *   parabolic 0.79/0.28, 0.24/0.53, 0.15/0.77.
 * Ep is the Young's modulus of the equivalent solid pile, Ep = EI / (π D⁴/64).
 * Reference: EN 1998-5:2004 Annex C, Table C.1; Gazetas (1991) in Fang, Foundation Engineering Handbook.
 * Validity: flexible (long) piles, linear elastic soil.
 * @param {{ D: number, EI: number, Es: number, profile?: string }} o D [m]; EI [N·m²] pile bending stiffness;
 *   Es [Pa] soil Young's modulus at depth D; profile 'homogeneous'|'linear'|'parabolic'
 * @returns {{ KLL: number, KLR: number, KRR: number, Ep: number }} [N/m], [N], [N·m/rad], [Pa]
 */
export function flexiblePileStiffnessContinuum({ D, EI, Es, profile = 'homogeneous' }) {
  pos('D', D); pos('EI', EI); pos('Es', Es);
  const c = EC8[profile];
  if (!c) throw new RangeError("profile: must be 'homogeneous', 'linear' or 'parabolic'");
  const Ep = EI / ((Math.PI * D ** 4) / 64), r = Ep / Es;
  return { KLL: c[0][0] * D * Es * r ** c[0][1], KLR: -c[1][0] * D * D * Es * r ** c[1][1], KRR: c[2][0] * D ** 3 * Es * r ** c[2][1], Ep };
}

/**
 * Three-spring foundation stiffness (KLL, KLR, KRR) of a monopile computed numerically with the nonlinear
 * Winkler model of pile.js at the working mudline loads (tangent and secant).
 * Reference: see pileHeadStiffness (pile.js). Validity: as the p–y curves used; API curves are known to
 * misrepresent the stiffness of large-diameter rigid monopiles — supply site-calibrated curves through pyAt,
 * base and mThetaAt where available.
 * @param {object} o solveLateralPile options; o.head = { H [N], M [N·m] } working loads at the pile head
 * @returns {{ tangent: object, secant: object, working: object }} each matrix as { KLL [N/m], KLR [N],
 *   KRR [N·m/rad] }; working = { y [m], theta [rad], H [N], M [N·m] }
 */
export function monopileStiffness(o) {
  const r = pileHeadStiffness(o);
  const map = (k) => ({ KLL: k.KHH, KLR: k.KHM, KRR: k.KMM });
  return { tangent: map(r.tangent), secant: map(r.secant), working: r.working };
}

/**
 * Natural frequencies of a tower (and any substructure above the mudline) carrying a top mass, on a flexible
 * foundation represented by KLL, KLR, KRR at its base. Euler–Bernoulli beam finite elements with consistent
 * mass; each segment may taper linearly in diameter (properties evaluated at element mid-height).
 * Reference: standard beam finite elements (e.g. Cook et al., Concepts and Applications of FEA, ch. 11); the
 * three-spring foundation idealisation of Arany et al. (2017) Soil Dyn. Earthq. Eng. 92. The closed-form
 * foundation flexibility coefficients of Arany et al. are not used: the eigenproblem is solved directly.
 * Validity: linear, undamped, planar bending; no axial-load (P–Δ) softening; water added mass only through
 * addedMass.
 * @param {{ segments: object[], topMass?: number, topInertia?: number, KLL?: number, KLR?: number,
 *   KRR?: number, nModes?: number, elementsPerSegment?: number }} o
 *   segments from the base upwards: [{ L [m]; Dbot, Dtop [m] outer diameters (or D); t [m] wall thickness;
 *   E [Pa] (default 210e9); rho [kg/m³] (default 7850); addedMass [kg/m] (default 0); or EI [N·m²] and
 *   massPerLength [kg/m] directly }];
 *   topMass [kg] rotor–nacelle mass (default 0); topInertia [kg·m²] its rotary inertia (default 0);
 *   KLL [N/m], KLR [N] (<= 0), KRR [N·m/rad] foundation stiffness at the base (omit all three for a fixed base);
 *   nModes [-] (default 3); elementsPerSegment [-] (default 12)
 * @returns {{ frequencies: number[], omega: number[], modeShapes: number[][], x: number[] }} frequencies [Hz];
 *   omega [rad/s]; modeShapes[k] lateral displacement at the nodes (normalised to unit maximum); x [m] node
 *   heights above the base
 */
export function towerNaturalFrequency(o) {
  const { segments, topMass = 0, topInertia = 0, nModes = 3, elementsPerSegment = 12 } = o;
  if (!Array.isArray(segments) || segments.length === 0) throw new RangeError('segments: must be a non-empty array');
  pos('topMass', topMass, true); pos('topInertia', topInertia, true);
  const fixed = o.KLL === undefined && o.KLR === undefined && o.KRR === undefined;
  if (!fixed) {
    pos('KLL', o.KLL); pos('KRR', o.KRR);
    const klr = o.KLR ?? 0;
    if (klr > 0) throw new RangeError('KLR: must be <= 0');
    if (o.KLL * o.KRR - klr * klr <= 0) throw new RangeError('KLR: foundation stiffness matrix must be positive definite');
  }
  const els = [], x = [0];
  for (const [si, s] of segments.entries()) {
    pos(`segments[${si}].L`, s.L);
    const h = s.L / elementsPerSegment;
    for (let e = 0; e < elementsPerSegment; e++) {
      let EI, m;
      if (s.EI !== undefined) { pos('EI', s.EI); pos('massPerLength', s.massPerLength); EI = s.EI; m = s.massPerLength; } else {
        const Db = s.Dbot ?? s.D, Dt = s.Dtop ?? s.D;
        pos('D', Db); pos('D', Dt); pos('t', s.t);
        const D = Db + ((Dt - Db) * (e + 0.5)) / elementsPerSegment, Di = D - 2 * s.t;
        if (Di <= 0) throw new RangeError(`segments[${si}].t: must be less than D/2`);
        EI = ((s.E ?? 210e9) * Math.PI * (D ** 4 - Di ** 4)) / 64;
        m = ((s.rho ?? 7850) * Math.PI * (D * D - Di * Di)) / 4 + (s.addedMass ?? 0);
      }
      els.push({ h, EI, m });
      x.push(x[x.length - 1] + h);
    }
  }
  const nn = els.length + 1, nd = 2 * nn;
  const K = Array.from({ length: nd }, () => new Array(nd).fill(0));
  const M = Array.from({ length: nd }, () => new Array(nd).fill(0));
  els.forEach((el, e) => {
    const { h, EI, m } = el, c = EI / h ** 3, g = (m * h) / 420;
    const ke = [[12, 6 * h, -12, 6 * h], [6 * h, 4 * h * h, -6 * h, 2 * h * h], [-12, -6 * h, 12, -6 * h], [6 * h, 2 * h * h, -6 * h, 4 * h * h]];
    const me = [[156, 22 * h, 54, -13 * h], [22 * h, 4 * h * h, 13 * h, -3 * h * h], [54, 13 * h, 156, -22 * h], [-13 * h, -3 * h * h, -22 * h, 4 * h * h]];
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) { K[2 * e + a][2 * e + b] += c * ke[a][b]; M[2 * e + a][2 * e + b] += g * me[a][b]; }
  });
  M[nd - 2][nd - 2] += topMass; M[nd - 1][nd - 1] += topInertia;
  let Kr = K, Mr = M, off = 0;
  if (fixed) { off = 2; Kr = K.slice(2).map((r) => r.slice(2)); Mr = M.slice(2).map((r) => r.slice(2)); } else {
    // with x upwards the rotation dy/dx equals theta of pile.js, so the foundation matrix is used as is
    K[0][0] += o.KLL; K[0][1] += o.KLR ?? 0; K[1][0] += o.KLR ?? 0; K[1][1] += o.KRR;
  }
  const nm = Math.min(nModes, Kr.length);
  // Solve the inverse problem M φ = μ K φ (μ = 1/ω²): the lowest frequencies are then the largest eigenvalues
  // and keep full relative accuracy even when stiff foundation springs or a heavy top mass spread the spectrum.
  const eig = eigGenSym(Mr, Kr);
  const nAll = eig.values.length;
  const omega = [], modeShapes = [];
  for (let k = 0; k < nm; k++) {
    const mu = eig.values[nAll - 1 - k];
    omega.push(mu > 0 ? 1 / Math.sqrt(mu) : Infinity);
    const v = eig.vectors[nAll - 1 - k], shape = new Array(nn).fill(0);
    let mx = 0;
    for (let i = 0; i < nn; i++) { const d = 2 * i - off; shape[i] = d >= 0 ? v[d] : 0; if (Math.abs(shape[i]) > Math.abs(mx)) mx = shape[i]; }
    modeShapes.push(shape.map((s) => s / (mx || 1)));
  }
  return { frequencies: omega.map((w) => w / (2 * Math.PI)), omega, modeShapes, x };
}

/**
 * Serviceability check of mudline rotation: total rotation and accumulated (permanent) rotation against
 * configurable limits. Defaults 0.5° total and 0.25° accumulated, the values commonly specified for offshore
 * wind monopiles (0.5° total tilt including a 0.25° installation tolerance).
 * Reference: DNV-ST-0126 guidance on permanent mudline rotation (limits are project-specific and normally set
 * by the turbine supplier). Validity: limits are inputs — confirm against the project design basis.
 * @param {{ thetaTotal: number, thetaAccumulated?: number, limitTotalDeg?: number,
 *   limitAccumulatedDeg?: number }} o thetaTotal [rad] total mudline rotation; thetaAccumulated [rad] permanent
 *   rotation accumulated under cyclic loading (default 0); limits [deg]
 * @returns {{ totalDeg: number, accumulatedDeg: number, ucTotal: number, ucAccumulated: number, ok: boolean }}
 */
export function mudlineRotationCheck({ thetaTotal, thetaAccumulated = 0, limitTotalDeg = 0.5, limitAccumulatedDeg = 0.25 }) {
  pos('limitTotalDeg', limitTotalDeg); pos('limitAccumulatedDeg', limitAccumulatedDeg);
  if (!Number.isFinite(thetaTotal) || !Number.isFinite(thetaAccumulated)) throw new RangeError('thetaTotal: rotations must be finite');
  const totalDeg = (Math.abs(thetaTotal) * 180) / Math.PI, accumulatedDeg = (Math.abs(thetaAccumulated) * 180) / Math.PI;
  const ucTotal = totalDeg / limitTotalDeg, ucAccumulated = accumulatedDeg / limitAccumulatedDeg;
  return { totalDeg, accumulatedDeg, ucTotal, ucAccumulated, ok: ucTotal <= 1 && ucAccumulated <= 1 };
}

/**
 * Accumulated rotation of a stiff pile in sand under N load cycles, LeBlanc, Houlsby & Byrne (2010):
 *   Δθ(N) / θs = Tb(ζb, Rd) · Tc(ζc) · N^0.31,
 * θs being the rotation in a static test at the maximum cyclic load, ζb = Mmax/MR and ζc = Mmin/Mmax.
 * Tb and Tc are given by the authors as charts (Tb increasing with ζb and depending on relative density;
 * Tc = 1 for one-way loading ζc = 0, 0 for ζc = 1, with a maximum of about 4 near ζc = −0.6); they are
 * arguments here and are not reproduced.
 * Reference: LeBlanc, Houlsby & Byrne (2010) Géotechnique 60(2) 79–90. Validity: drained sand, laboratory 1 g
 * model tests on rigid piles, N up to ~10⁴–10⁵ in the tests.
 * @param {{ thetaStatic: number, N: number, Tb: number, Tc: number, exponent?: number }} o thetaStatic [rad]
 *   θs; N [-] number of cycles (>= 1); Tb [-]; Tc [-]; exponent [-] (default 0.31)
 * @returns {{ deltaTheta: number, ratio: number }} deltaTheta [rad] accumulated rotation; ratio = Δθ/θs [-]
 */
export function cyclicRotationAccumulation({ thetaStatic, N, Tb, Tc, exponent = 0.31 }) {
  pos('thetaStatic', thetaStatic, true); pos('Tb', Tb, true); pos('Tc', Tc, true); pos('exponent', exponent);
  if (!(N >= 1)) throw new RangeError('N: must be >= 1');
  const ratio = Tb * Tc * Math.pow(N, exponent);
  return { deltaTheta: ratio * thetaStatic, ratio };
}

/**
 * Embedded-length criteria for a monopile by scanning trial lengths with the Winkler model:
 *   zero toe-kick — shortest length at which the toe deflection changes sign (the toe no longer moves backwards);
 *   vertical tangent — shortest length at which the toe rotation changes sign (deflected shape vertical at the toe);
 *   critical length — shortest length whose head deflection and rotation are within tol of the longest trial pile.
 * Lengths are interpolated linearly between the scanned values.
 * Reference: embedded-length criteria commonly applied in monopile design (zero toe-kick and vertical-tangent
 * criteria of the Germanischer Lloyd offshore wind guideline; critical pile length); see solveLateralPile.
 * Validity: as the p–y model; the criteria are length heuristics, not capacity checks.
 * @param {object} o solveLateralPile options without L (o.head carries the design load)
 * @param {{ Lmin: number, Lmax: number, nLengths?: number, tol?: number }} c Lmin, Lmax [m]; nLengths [-]
 *   (default 25); tol [-] for the critical length (default 0.02)
 * @returns {{ L: number[], headY: number[], headTheta: number[], toeY: number[], toeTheta: number[],
 *   zeroToeKick: number|null, verticalTangent: number|null, critical: number|null }} lengths [m], deflections
 *   [m], rotations [rad]; null where a criterion is not met within the range
 */
export function embeddedLengthCriteria(o, c) {
  const { Lmin, Lmax, nLengths = 25, tol = 0.02 } = c;
  pos('Lmin', Lmin); pos('Lmax', Lmax); pos('tol', tol);
  if (!(Lmax > Lmin)) throw new RangeError('Lmax: must exceed Lmin');
  if (!(nLengths >= 3)) throw new RangeError('nLengths: must be >= 3');
  const hRef = Lmax / (o.nElem ?? 60);
  const out = { L: [], headY: [], headTheta: [], toeY: [], toeTheta: [], zeroToeKick: null, verticalTangent: null, critical: null };
  for (let i = 0; i < nLengths; i++) {
    const L = Lmin + ((Lmax - Lmin) * i) / (nLengths - 1);
    const r = solveLateralPile({ ...o, L, nElem: Math.max(8, Math.round(L / hRef)), throwOnFail: false });
    const ok = r.convergence.converged;
    out.L.push(L); out.headY.push(ok ? r.head.y : NaN); out.headTheta.push(ok ? r.head.theta : NaN);
    out.toeY.push(ok ? r.toe.y : NaN); out.toeTheta.push(ok ? r.toe.theta : NaN);
  }
  const crossing = (arr) => {
    for (let i = 1; i < nLengths; i++) {
      const a = arr[i - 1], b = arr[i];
      if (Number.isFinite(a) && Number.isFinite(b) && a * b <= 0 && a !== b) return out.L[i - 1] + ((out.L[i] - out.L[i - 1]) * a) / (a - b);
    }
    return null;
  };
  out.zeroToeKick = crossing(out.toeY);
  out.verticalTangent = crossing(out.toeTheta);
  const yRef = out.headY[nLengths - 1], tRef = out.headTheta[nLengths - 1];
  if (Number.isFinite(yRef)) {
    for (let i = nLengths - 2; i >= 0; i--) {
      if (!(Math.abs(out.headY[i] / yRef - 1) <= tol && Math.abs(out.headTheta[i] / tRef - 1) <= tol)) break;
      out.critical = out.L[i];
    }
  }
  return out;
}

/**
 * Minimum pile wall thickness for driving: t = 6.35 + D/100 (t and D in mm).
 * Reference: API RP 2A-WSD 21st ed. §6.10.6 (Eq. 6.10.6-1). Validity: driven steel pipe piles; larger
 * thicknesses may be needed for hard driving, buckling or fatigue.
 * @param {number} D outer diameter [m]
 * @returns {number} minimum wall thickness [m]
 */
export function apiMinimumWallThickness(D) {
  pos('D', D);
  return 0.00635 + D / 100;
}

/**
 * Preliminary diameter/thickness screening: API minimum wall, D/t ratio and a slenderness flag against a
 * user-supplied D/t limit (hook for a shell-buckling check elsewhere).
 * Reference: API RP 2A-WSD §6.10.6 for the minimum wall; D/t limit is an input (API RP 2A-WSD requires a local
 * buckling check for D/t > 60). Validity: screening only.
 * @param {{ D: number, t: number, dOverTLimit?: number }} o D [m]; t [m]; dOverTLimit [-] (default 60)
 * @returns {{ tMin: number, wallOk: boolean, dOverT: number, bucklingCheckRequired: boolean, area: number,
 *   inertia: number, massPerLength: number }} tMin [m]; area [m²]; inertia [m⁴]; massPerLength [kg/m] for steel
 *   of 7850 kg/m³
 */
export function monopileSizingScreen({ D, t, dOverTLimit = 60 }) {
  pos('D', D); pos('t', t); pos('dOverTLimit', dOverTLimit);
  if (t >= D / 2) throw new RangeError('t: must be less than D/2');
  const Di = D - 2 * t, area = (Math.PI * (D * D - Di * Di)) / 4;
  const tMin = apiMinimumWallThickness(D);
  return { tMin, wallOk: t >= tMin, dOverT: D / t, bucklingCheckRequired: D / t > dOverTLimit, area, inertia: (Math.PI * (D ** 4 - Di ** 4)) / 64, massPerLength: 7850 * area };
}
