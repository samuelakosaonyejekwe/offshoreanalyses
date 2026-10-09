// Lateral soil resistance (p–y) curves for piles.
// A curve object is { type, pu [N/m], p(y) [N/m], k(y) [N/m²], meta } with p odd in y [m] and k = dp/dy.
// The same interface is used for user-defined base-shear, base-moment and distributed-moment springs.

import { propertiesAt, scourBaseDepth, verticalEffectiveStress } from './soilprofile.js';

const RAD = Math.PI / 180;

function pos(name, v, allowZero = false) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < 0 || (!allowZero && v === 0)) {
    throw new RangeError(`${name}: must be a finite number ${allowZero ? '>= 0' : '> 0'} (got ${v})`);
  }
}

function makeCurve(type, pu, f, df, meta) {
  return {
    type, pu, meta,
    p: (y) => (y < 0 ? -f(-y) : f(y)),
    k: (y) => df(Math.abs(y)),
  };
}

/**
 * A curve with no resistance (soil removed by scour, or stick-up above the mudline).
 * @returns {{ type: string, pu: number, p: Function, k: Function, meta: object }} zero curve; pu = 0 N/m
 */
export function zeroCurve() {
  return makeCurve('none', 0, () => 0, () => 0, {});
}

/**
 * Ultimate lateral resistance of clay, Matlock (1970): pu = min[(3 su + σ′v) D + J su z, 9 su D].
 * σ′v replaces γ′X of the standard so that layered profiles and scour can be handled; for a uniform soil the two
 * are identical.
 * Reference: API RP 2GEO (2011) §8.5 (API RP 2A-WSD §6.8.2); Matlock (1970) OTC 1204.
 * Validity: soft to firm clay under static or cyclic loading; J = 0.5 soft Gulf of Mexico clay, 0.25 stiffer clay.
 * @param {{ su: number, sigmaV: number, z: number, D: number, J?: number }} o su [Pa] undrained shear strength at
 *   the depth; sigmaV [Pa] vertical effective stress; z [m] depth below the soil surface; D [m] pile diameter;
 *   J [-] empirical constant 0.25–0.5 (default 0.5)
 * @returns {{ pu: number, Np: number, shallow: boolean, zOverZR: number }} pu [N/m]; Np = pu/(su D) [-]; shallow
 *   true when the wedge expression governs; zOverZR [-] = (Np,wedge − 3)/6, equal to X/XR for constant su and γ′
 */
export function matlockPu({ su, sigmaV, z, D, J = 0.5 }) {
  pos('su', su); pos('sigmaV', sigmaV, true); pos('z', z, true); pos('D', D); pos('J', J, true);
  const NpS = 3 + sigmaV / su + (J * z) / D;
  const Np = Math.min(NpS, 9);
  return { pu: Np * su * D, Np, shallow: NpS < 9, zOverZR: (NpS - 3) / 6 };
}

/**
 * Depth of the reduced-resistance zone for constant su and γ′: XR = 6 D / (γ′ D / su + J).
 * Reference: API RP 2GEO (2011) §8.5 (API RP 2A-WSD Eq. 6.8.2-3). Validity: strength constant with depth.
 * @param {{ su: number, gammaSub: number, D: number, J?: number }} o su [Pa]; gammaSub [N/m³]; D [m]; J [-]
 * @returns {number} XR [m]
 */
export function matlockTransitionDepth({ su, gammaSub, D, J = 0.5 }) {
  pos('su', su); pos('gammaSub', gammaSub, true); pos('D', D); pos('J', J, true);
  return (6 * D) / ((gammaSub * D) / su + J);
}

/**
 * Soft-clay p–y curve, Matlock (1970). yc = 2.5 ε50 D.
 * Static: p/pu = 0.5 (y/yc)^(1/3) for y <= 8 yc, p = pu beyond.
 * Cyclic: same up to y = 3 yc (p/pu = 0.5·3^(1/3) = 0.721, tabulated as 0.72); for z >= XR constant thereafter;
 * for z < XR falling linearly to 0.72 pu z/XR at y = 15 yc and constant beyond.
 * The power law has infinite slope at the origin; for |y| < yLinearRatio·yc the curve is replaced by its secant,
 * which is the first segment of the tabulated API curve (0.23 pu at 0.1 yc). Set yLinearRatio = 0 for the pure law.
 * Reference: API RP 2GEO (2011) §8.5 (API RP 2A-WSD §6.8.3); Matlock (1970) OTC 1204.
 * Validity: soft clay; ε50 from laboratory tests (typically 0.005–0.02).
 * @param {{ su: number, sigmaV: number, z: number, D: number, eps50: number, J?: number, loading?: string,
 *   zOverZR?: number, yLinearRatio?: number }} o su [Pa]; sigmaV [Pa]; z [m] below soil surface; D [m]; eps50 [-];
 *   J [-] (default 0.5); loading 'static'|'cyclic'; zOverZR [-] ratio X/XR for the cyclic curve (default from
 *   matlockPu); yLinearRatio [-] (default 0.1)
 * @returns {object} curve { type: 'soft-clay', pu [N/m], p(y), k(y), meta: { yc [m], Np, zOverZR, loading } }
 */
export function pySoftClay({ su, sigmaV, z, D, eps50, J = 0.5, loading = 'static', zOverZR, yLinearRatio = 0.1 }) {
  pos('eps50', eps50); pos('yLinearRatio', yLinearRatio, true);
  if (loading !== 'static' && loading !== 'cyclic') throw new RangeError("loading: must be 'static' or 'cyclic'");
  const u = matlockPu({ su, sigmaV, z, D, J });
  const pu = u.pu;
  const yc = 2.5 * eps50 * D;
  const r = Math.min(zOverZR ?? u.zOverZR, 1);
  const a0 = yLinearRatio;
  const s0 = a0 > 0 ? 0.5 * Math.cbrt(a0) / a0 : 0;
  const aPeak = loading === 'static' ? 8 : 3;
  const fPeak = 0.5 * Math.cbrt(aPeak);
  const f = (y) => {
    const a = y / yc;
    if (a < a0) return pu * s0 * a;
    if (a <= aPeak) return pu * 0.5 * Math.cbrt(a);
    if (loading === 'static' || r >= 1) return pu * fPeak;
    if (a < 15) return pu * fPeak * (1 - ((1 - r) * (a - 3)) / 12);
    return pu * fPeak * r;
  };
  const df = (y) => {
    const a = y / yc;
    if (a < a0) return (pu * s0) / yc;
    if (a <= aPeak) return a === 0 ? Infinity : (pu / (6 * yc)) * Math.pow(a, -2 / 3);
    if (loading === 'static' || r >= 1 || a >= 15) return 0;
    return (-pu * fPeak * (1 - r)) / (12 * yc);
  };
  return makeCurve('soft-clay', pu, f, df, { yc, Np: u.Np, zOverZR: u.zOverZR, loading });
}

/**
 * Stiff clay without free water, Welch & Reese (1972) / Reese & Welch (1975):
 * static p/pu = 0.5 (y/y50)^(1/4) <= 1 (reached at 16 y50), y50 = 2.5 ε50 D, pu as Matlock with J = 0.5;
 * cyclic yc = ys + y50 C log10(N) with C = 9.6 (p/pu)^4, which gives p/pu = [y / (y50 (16 + 9.6 log10 N))]^(1/4).
 * This is NOT the Reese, Cox & Koop (1975) curve for stiff clay with free water, whose As/Ac coefficients are
 * given only as figures and are not reproduced here.
 * Reference: Reese & Welch (1975) J. Geotech. Eng. Div. ASCE 101(GT7). Validity: stiff clay above the water
 * table; offered for comparison only in submerged conditions.
 * @param {{ su: number, sigmaV: number, z: number, D: number, eps50: number, J?: number, loading?: string,
 *   nCycles?: number, yLinearRatio?: number }} o su [Pa]; sigmaV [Pa]; z [m]; D [m]; eps50 [-]; J [-] (0.5);
 *   loading 'static'|'cyclic'; nCycles [-] number of load cycles for cyclic (default 100); yLinearRatio [-]
 *   secant regularisation below this multiple of y50 (default 0.01)
 * @returns {object} curve { type: 'stiff-clay-no-free-water', pu [N/m], p, k, meta: { y50 [m], nCycles } }
 */
export function pyStiffClayNoFreeWater({ su, sigmaV, z, D, eps50, J = 0.5, loading = 'static', nCycles = 100, yLinearRatio = 0.01 }) {
  pos('eps50', eps50); pos('yLinearRatio', yLinearRatio, true);
  if (loading === 'cyclic' && !(nCycles >= 1)) throw new RangeError('nCycles: must be >= 1');
  const { pu } = matlockPu({ su, sigmaV, z, D, J });
  const y50 = 2.5 * eps50 * D;
  // p/pu = (y / yRef)^(1/4) with yRef = 16 y50 static, y50 (16 + 9.6 log10 N) cyclic
  const yRef = y50 * (16 + (loading === 'cyclic' ? 9.6 * Math.log10(nCycles) : 0));
  const y0 = yLinearRatio * y50;
  const s0 = y0 > 0 ? (pu * Math.pow(y0 / yRef, 0.25)) / y0 : 0;
  const f = (y) => (y < y0 ? s0 * y : y >= yRef ? pu : pu * Math.pow(y / yRef, 0.25));
  const df = (y) => (y < y0 ? s0 : y >= yRef ? 0 : y === 0 ? Infinity : (0.25 * pu * Math.pow(y / yRef, 0.25)) / y);
  return makeCurve('stiff-clay-no-free-water', pu, f, df, { y50, nCycles: loading === 'cyclic' ? nCycles : 0, loading });
}

/**
 * Coefficients C1, C2, C3 of the API sand ultimate resistance, evaluated from the Reese, Cox & Koop (1974) wedge
 * and flow-around expressions with α = φ′/2, β = 45° + φ′/2, K0 = 0.4, Ka = tan²(45° − φ′/2):
 *   C1 = K0 tanφ sinβ / (tan(β−φ) cosα) + tan²β tanα / tan(β−φ) + K0 tanβ (tanφ sinβ − tanα)
 *   C2 = tanβ / tan(β−φ) − Ka
 *   C3 = Ka (tan⁸β − 1) + K0 tanφ tan⁴β
 * Reference: API RP 2GEO (2011) §8.5, coefficient figure (curves reproduced by these expressions); Reese et al.
 * (1974) OTC 2080. Validity: 20° <= φ′ <= 45°.
 * @param {number} phiDeg effective friction angle [deg]
 * @returns {{ C1: number, C2: number, C3: number }} dimensionless coefficients
 */
export function sandCoefficients(phiDeg) {
  if (!(phiDeg >= 20 && phiDeg <= 45)) throw new RangeError(`phiDeg: must be within 20..45 deg (got ${phiDeg})`);
  const phi = phiDeg * RAD;
  const al = phi / 2;
  const be = Math.PI / 4 + phi / 2;
  const K0 = 0.4;
  const Ka = Math.tan(Math.PI / 4 - phi / 2) ** 2;
  const tp = Math.tan(phi), tb = Math.tan(be), ta = Math.tan(al), tbp = Math.tan(be - phi);
  const C1 = (K0 * tp * Math.sin(be)) / (tbp * Math.cos(al)) + (tb * tb * ta) / tbp + K0 * tb * (tp * Math.sin(be) - ta);
  const C2 = tb / tbp - Ka;
  const C3 = Ka * (tb ** 8 - 1) + K0 * tp * tb ** 4;
  return { C1, C2, C3 };
}

/**
 * Initial modulus of subgrade reaction k(φ′) for sand below the water table, closed-form approximation of the API
 * figure: k = (0.008085 φ′^2.45 − 26.09)·10⁶ N/m³ (φ′ in degrees).
 * This fit is quoted in the monopile literature as an approximation to the API RP 2GEO k–φ′ figure (below water
 * table). It returns 16.4 and 33.9 MN/m³ at 33° and 38°, i.e. the medium and dense submerged values recommended
 * by Reese et al. (1974) (16.3 and 34 MN/m³), and 4.9 MN/m³ at 29°. The deviation from the printed API curve
 * has not been quantified here; supply k directly when the figure value is required.
 * Reference: API RP 2GEO (2011) §8.5, k–φ′ figure. Validity: 29° <= φ′ <= 45°, submerged sand.
 * @param {number} phiDeg effective friction angle [deg]
 * @returns {number} k [N/m³]
 */
export function sandInitialModulus(phiDeg) {
  if (!(phiDeg >= 29 && phiDeg <= 45)) throw new RangeError(`phiDeg: k(φ′) fit valid for 29..45 deg (got ${phiDeg}); supply k explicitly`);
  return (0.008085 * Math.pow(phiDeg, 2.45) - 26.09) * 1e6;
}

/**
 * Ultimate lateral resistance of sand: pu = min[(C1 z + C2 D) σ′v, C3 D σ′v] (σ′v = γ′ z for uniform sand).
 * Reference: API RP 2GEO (2011) §8.5 (API RP 2A-WSD Eq. 6.8.6-1/-2). Validity: 20° <= φ′ <= 45°.
 * @param {{ phiDeg: number, sigmaV: number, z: number, D: number }} o phiDeg [deg]; sigmaV [Pa]; z [m] depth
 *   below soil surface; D [m]
 * @returns {{ pu: number, pus: number, pud: number, shallow: boolean }} pu, shallow (pus) and deep (pud) [N/m]
 */
export function sandPu({ phiDeg, sigmaV, z, D }) {
  pos('sigmaV', sigmaV, true); pos('z', z, true); pos('D', D);
  const { C1, C2, C3 } = sandCoefficients(phiDeg);
  const pus = (C1 * z + C2 * D) * sigmaV;
  const pud = C3 * D * sigmaV;
  return { pu: Math.min(pus, pud), pus, pud, shallow: pus <= pud };
}

/**
 * Depth at which the shallow and deep sand expressions are equal: z = (C3 − C2) D / C1.
 * Reference: follows from API RP 2GEO (2011) §8.5 equations. Validity: as sandCoefficients.
 * @param {number} phiDeg effective friction angle [deg]
 * @param {number} D pile diameter [m]
 * @returns {number} transition depth below soil surface [m]
 */
export function sandTransitionDepth(phiDeg, D) {
  pos('D', D);
  const { C1, C2, C3 } = sandCoefficients(phiDeg);
  return ((C3 - C2) * D) / C1;
}

/**
 * Sand p–y curve (O'Neill & Murchison 1983): p = A pu tanh[k z y / (A pu)], A = max(3 − 0.8 z/D, 0.9) static,
 * A = 0.9 cyclic. Initial slope dp/dy(0) = k z.
 * Reference: API RP 2GEO (2011) §8.5 (API RP 2A-WSD §6.8.7). Validity: 20° <= φ′ <= 45°; calibrated on
 * slender piles (D < ~2 m) — unconservative initial stiffness trends are reported for large monopiles.
 * @param {{ phiDeg: number, sigmaV: number, z: number, D: number, k?: number, loading?: string }} o phiDeg [deg];
 *   sigmaV [Pa]; z [m] depth below soil surface; D [m]; k [N/m³] initial modulus of subgrade reaction (default
 *   sandInitialModulus(phiDeg)); loading 'static'|'cyclic'
 * @returns {object} curve { type: 'sand', pu [N/m] (= A·pu,API, the asymptote), p, k, meta: { A, puApi [N/m],
 *   kz [N/m²], shallow } }
 */
export function pySand({ phiDeg, sigmaV, z, D, k, loading = 'static' }) {
  if (loading !== 'static' && loading !== 'cyclic') throw new RangeError("loading: must be 'static' or 'cyclic'");
  const u = sandPu({ phiDeg, sigmaV, z, D });
  const kk = k ?? sandInitialModulus(phiDeg);
  pos('k', kk);
  const A = loading === 'cyclic' ? 0.9 : Math.max(3 - (0.8 * z) / D, 0.9);
  const pa = A * u.pu;
  const kz = kk * z;
  if (!(pa > 0)) return makeCurve('sand', 0, () => 0, () => 0, { A, puApi: 0, kz, shallow: true });
  const f = (y) => pa * Math.tanh((kz * y) / pa);
  const df = (y) => { const c = Math.cosh((kz * y) / pa); return kz / (c * c); };
  return makeCurve('sand', pa, f, df, { A, puApi: u.pu, kz, shallow: u.shallow });
}

/**
 * Weak-rock p–y curve, Reese (1997):
 *   pur = αr qur b (1 + 1.4 xr/b) for xr <= 3b, 5.2 αr qur b below;
 *   Kir = kir Eir, kir = 100 + 400 xr/(3b) for xr <= 3b, 500 below; yrm = krm b;
 *   p = Kir y for y <= yA; p = (pur/2)(y/yrm)^0.25 <= pur beyond; yA = [pur / (2 yrm^0.25 Kir)]^(4/3).
 * Reference: Reese (1997) J. Geotech. Geoenviron. Eng. 123(11) 1010–1017. Validity: weak rock, qur roughly
 * 0.5–5 MPa; "interim" method based on two load tests.
 * @param {{ qur: number, Eir: number, z: number, D: number, alphaR?: number, krm?: number }} o qur [Pa] uniaxial
 *   compressive strength; Eir [Pa] initial rock-mass modulus; z [m] depth below rock surface; D [m] pile diameter;
 *   alphaR [-] strength reduction factor 1/3..1 (default 1); krm [-] strain constant 0.00005–0.0005 (required)
 * @returns {object} curve { type: 'weak-rock', pu [N/m], p, k, meta: { Kir [Pa], yrm [m], yA [m] } }
 */
export function pyWeakRock({ qur, Eir, z, D, alphaR = 1, krm }) {
  pos('qur', qur); pos('Eir', Eir); pos('z', z, true); pos('D', D); pos('krm', krm);
  if (!(alphaR > 0 && alphaR <= 1)) throw new RangeError('alphaR: must be in (0, 1]');
  const shallow = z <= 3 * D;
  const pur = shallow ? alphaR * qur * D * (1 + (1.4 * z) / D) : 5.2 * alphaR * qur * D;
  const Kir = (shallow ? 100 + (400 * z) / (3 * D) : 500) * Eir;
  const yrm = krm * D;
  const yA = Math.pow(pur / (2 * Math.pow(yrm, 0.25) * Kir), 4 / 3);
  const yU = 16 * yrm;
  const f = (y) => (y <= yA ? Kir * y : y >= yU ? pur : 0.5 * pur * Math.pow(y / yrm, 0.25));
  const df = (y) => (y <= yA ? Kir : y >= yU ? 0 : (0.125 * pur * Math.pow(y / yrm, 0.25)) / y);
  return makeCurve('weak-rock', pur, (y) => Math.min(f(y), pur), df, { Kir, yrm, yA });
}

/**
 * Clay p–y curve of Jeanjean (2009): p = pu tanh[(Gmax / (100 su)) (y/D)^0.5], pu = Np su D,
 * Np = 12 − 4 exp(−ξ z/D), ξ = 0.25 + 0.05 λ for λ < 6 and 0.55 for λ >= 6, λ = su0 / (su1 D) with
 * su = su0 + su1 z.
 * The square-root law has infinite initial slope; below yLinear the secant is used.
 * Reference: Jeanjean (2009) OTC 20158. Validity: soft to firm clay, monotonic loading, conductors/piles in
 * normally to lightly overconsolidated clay.
 * @param {{ su: number, su0: number, su1: number, Gmax: number, z: number, D: number, yLinear?: number }} o
 *   su [Pa] at the depth; su0 [Pa] mudline intercept; su1 [Pa/m] gradient (0 allowed → ξ = 0.55); Gmax [Pa];
 *   z [m]; D [m]; yLinear [m] secant regularisation limit (default 1e-4 D)
 * @returns {object} curve { type: 'jeanjean', pu [N/m], p, k, meta: { Np, xi } }
 */
export function pyJeanjean({ su, su0, su1, Gmax, z, D, yLinear }) {
  pos('su', su); pos('su0', su0, true); pos('su1', su1, true); pos('Gmax', Gmax); pos('z', z, true); pos('D', D);
  const lambda = su1 > 0 ? su0 / (su1 * D) : Infinity;
  const xi = lambda < 6 ? 0.25 + 0.05 * lambda : 0.55;
  const Np = 12 - 4 * Math.exp((-xi * z) / D);
  const pu = Np * su * D;
  const c = Gmax / (100 * su);
  const y0 = yLinear ?? 1e-4 * D;
  pos('yLinear', y0);
  const g = (y) => pu * Math.tanh(c * Math.sqrt(y / D));
  const s0 = g(y0) / y0;
  const f = (y) => (y < y0 ? s0 * y : g(y));
  const df = (y) => {
    if (y < y0) return s0;
    const r = Math.sqrt(y / D);
    const ch = Math.cosh(c * r);
    return (pu * c) / (2 * D * r * ch * ch);
  };
  return makeCurve('jeanjean', pu, f, df, { Np, xi });
}

/**
 * User-defined curve by piecewise-linear interpolation through tabulated points, odd-symmetric, constant beyond
 * the last point. Also used for base-shear (N vs m), base-moment (N·m vs rad) and distributed-moment springs.
 * Reference: none (generic). Validity: y strictly ascending and starting at 0 with p = 0.
 * @param {{ y: number[], p: number[] }} o y [m or rad] abscissae, ascending from 0; p [N/m, N or N·m] ordinates
 * @returns {object} curve { type: 'user', pu (largest ordinate), p, k, meta: {} }
 */
export function pyUser({ y, p }) {
  if (!Array.isArray(y) || !Array.isArray(p) || y.length !== p.length || y.length < 2) throw new RangeError('y: y and p must be arrays of equal length >= 2');
  if (y[0] !== 0 || p[0] !== 0) throw new RangeError('y: the first point must be (0, 0)');
  for (let i = 1; i < y.length; i++) if (!(y[i] > y[i - 1])) throw new RangeError('y: must be strictly ascending');
  const n = y.length;
  const seg = (x) => { let i = 1; while (i < n - 1 && x > y[i]) i++; return i; };
  const f = (x) => { if (x >= y[n - 1]) return p[n - 1]; const i = seg(x); return p[i - 1] + ((p[i] - p[i - 1]) * (x - y[i - 1])) / (y[i] - y[i - 1]); };
  const df = (x) => { if (x >= y[n - 1]) return 0; const i = seg(x); return (p[i] - p[i - 1]) / (y[i] - y[i - 1]); };
  return makeCurve('user', Math.max(...p), f, df, {});
}

/**
 * Linear spring as a curve object: p = stiffness · y.
 * @param {number} stiffness spring stiffness [N/m² for p–y; N/m or N·m/rad for lumped springs], >= 0
 * @returns {object} curve { type: 'linear', pu: Infinity, p, k, meta: { stiffness } }
 */
export function pyLinear(stiffness) {
  pos('stiffness', stiffness, true);
  return makeCurve('linear', Infinity, (y) => stiffness * y, () => stiffness, { stiffness });
}

/**
 * Scale a curve: p′(y) = pMult · p(y / yMult). pMult < 1 models group shadowing or cyclic degradation of
 * resistance; yMult > 1 softens the response (group y-multiplier).
 * Reference: p-multiplier concept of Brown, Morrison & Reese (1988); y-multiplier per Focht & Koch (1973).
 * Validity: pMult >= 0, yMult > 0.
 * @param {object} curve curve object
 * @param {{ pMult?: number, yMult?: number }} [o] pMult [-] (default 1); yMult [-] (default 1)
 * @returns {object} scaled curve with the same interface
 */
export function scalePy(curve, { pMult = 1, yMult = 1 } = {}) {
  pos('pMult', pMult, true); pos('yMult', yMult);
  if (pMult === 1 && yMult === 1) return curve;
  return {
    type: curve.type, pu: curve.pu * pMult, meta: { ...curve.meta, pMult, yMult },
    p: (y) => pMult * curve.p(y / yMult),
    k: (y) => (pMult / yMult) * curve.k(y / yMult),
  };
}

/**
 * Sample a curve for plotting or export to another program.
 * @param {object} curve curve object
 * @param {{ yMax: number, n?: number }} o yMax [m] largest displacement; n [-] number of points (default 41)
 * @returns {{ y: number[], p: number[] }} y [m], p [N/m]; points are quadratically clustered towards y = 0
 */
export function samplePy(curve, { yMax, n = 41 }) {
  pos('yMax', yMax);
  if (!(n >= 2)) throw new RangeError('n: must be >= 2');
  const y = [], p = [];
  for (let i = 0; i < n; i++) { const s = i / (n - 1); const yy = yMax * s * s; y.push(yy); p.push(curve.p(yy)); }
  return { y, p };
}

/**
 * Initial (small-displacement) stiffness of a curve as the secant p(y)/y at a reference displacement.
 * @param {object} curve curve object
 * @param {number} [yRef] reference displacement [m] (default 1e-6)
 * @returns {number} secant stiffness [N/m²]
 */
export function initialStiffness(curve, yRef = 1e-6) {
  pos('yRef', yRef);
  return curve.p(yRef) / yRef;
}

/**
 * p–y curve at a depth of a layered profile, selected by soil type:
 *   clay (and silt with su): Matlock soft clay. For su above stiffClaySu the same curve shape is used and, under
 *     cyclic loading, the resistance is multiplied by stiffClayCyclicFactor — the API guidance for stiff clay
 *     ("more brittle … rapid deterioration under cyclic loading") gives no equation, so the factor is an explicit
 *     engineering input (default 1 = no additional degradation);
 *   sand (and silt without su): API sand;
 *   rock: Reese (1997) weak rock using layer.qur, layer.Eir, layer.krm (and layer.alphaR).
 * Layer fields J, k override the options for that layer. Above the scoured surface a zero curve is returned.
 * Depth in the curve equations is measured from the scoured soil surface; σ′v includes the scour overburden
 * reduction. Layering corrections (e.g. Georgiadis equivalent depth) are not applied.
 * Reference: API RP 2GEO (2011) §8.5; see the individual curve functions.
 * @param {object} profile soil profile (soilprofile.js)
 * @param {number} z depth below original mudline [m]
 * @param {{ D: number, loading?: string, J?: number, stiffClaySu?: number, stiffClayCyclicFactor?: number,
 *   yLinearRatio?: number, pMult?: number, yMult?: number }} o D [m]; loading 'static'|'cyclic'; J [-] (0.5);
 *   stiffClaySu [Pa] threshold (default 96e3, API); stiffClayCyclicFactor [-] in (0, 1]; yLinearRatio [-];
 *   pMult, yMult [-] multipliers applied to the finished curve
 * @returns {object} curve object
 */
export function pyCurveAt(profile, z, o) {
  const { D, loading = 'static', J = 0.5, stiffClaySu = 96e3, stiffClayCyclicFactor = 1, yLinearRatio = 0.1, pMult = 1, yMult = 1 } = o;
  pos('D', D);
  if (!(stiffClayCyclicFactor > 0 && stiffClayCyclicFactor <= 1)) throw new RangeError('stiffClayCyclicFactor: must be in (0, 1]');
  const base = scourBaseDepth(profile);
  if (z < base - 1e-12) return zeroCurve();
  const pr = propertiesAt(profile, z);
  const zs = Math.max(z - base, 0);
  const sigmaV = verticalEffectiveStress(profile, z);
  const l = pr.layer;
  let curve;
  const cohesive = pr.type === 'clay' || (pr.type === 'silt' && pr.su !== undefined);
  if (pr.type === 'rock') {
    if (l.qur === undefined || l.Eir === undefined || l.krm === undefined) throw new RangeError(`profile: rock layer ${pr.layerIndex} needs qur, Eir and krm`);
    curve = pyWeakRock({ qur: l.qur, Eir: l.Eir, z: zs, D, alphaR: l.alphaR ?? 1, krm: l.krm });
  } else if (cohesive) {
    if (pr.eps50 === undefined) throw new RangeError(`profile: layer ${pr.layerIndex} needs eps50 for a clay p–y curve`);
    if (!(pr.su > 0)) curve = zeroCurve();
    else {
      curve = pySoftClay({ su: pr.su, sigmaV, z: zs, D, eps50: pr.eps50, J: l.J ?? J, loading, yLinearRatio });
      if (loading === 'cyclic' && pr.su > stiffClaySu) curve = scalePy(curve, { pMult: stiffClayCyclicFactor });
    }
  } else {
    if (pr.phiDeg === undefined) throw new RangeError(`profile: layer ${pr.layerIndex} needs phiDeg for a sand p–y curve`);
    curve = pySand({ phiDeg: pr.phiDeg, sigmaV, z: zs, D, k: pr.k, loading });
  }
  return scalePy(curve, { pMult, yMult });
}
