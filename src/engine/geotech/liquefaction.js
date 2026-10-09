/**
 * Liquefaction triggering and consequences: simplified (stress-based) procedure for earthquakes and
 * screening for wave-induced liquefaction of the seabed.
 *
 * SI units: depths m, stresses Pa, accelerations m/s^2, velocities m/s. Fines content FC in percent.
 * Atmospheric pressure Pa = 101 325 Pa unless stated.
 *
 * References
 *   Seed H.B. & Idriss I.M. (1971) Simplified procedure for evaluating soil liquefaction potential.
 *   Liao S.S.C. & Whitman R.V. (1986) Overburden correction factors for SPT in sand; and rd.
 *   Youd T.L. et al. (2001) Liquefaction resistance of soils: summary report from the 1996 NCEER and 1998
 *     NCEER/NSF workshops, JGGE 127(10) – equation numbers below refer to this paper.
 *   Robertson P.K. & Wride C.E. (1998) Can. Geotech. J. 35(3).
 *   Andrus R.D. & Stokoe K.H. (2000) JGGE 126(11).
 *   Iwasaki T. et al. (1982) Microzonation for soil liquefaction potential.
 *   Ishihara K. & Yoshimine M. (1992) Soils Found. 32(1); Yoshimine M. et al. (2006); Idriss I.M. &
 *     Boulanger R.W. (2008) Soil Liquefaction During Earthquakes, EERI MNO-12.
 *   Seed H.B. & Rahman M.S. (1978) Marine Geotechnology 3(2); Zen K. & Yamazaki H. (1990) Soils Found. 30(4).
 */

const PATM = 101325;
const G0 = 9.80665;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

/**
 * Stress-reduction coefficient rd.
 *   'liao-whitman' (Youd et al. 2001 eq. 2): rd = 1 − 0.00765 z (z ≤ 9.15 m); 1.174 − 0.0267 z
 *     (9.15 < z ≤ 23 m); extended by 0.744 − 0.008 z (23 < z ≤ 30 m) and 0.5 below 30 m.
 *   'blake' (Youd et al. 2001 eq. 3): rd = (1 − 0.4113 z^0.5 + 0.04052 z + 0.001753 z^1.5) /
 *     (1 − 0.4177 z^0.5 + 0.05729 z − 0.006205 z^1.5 + 0.001210 z²).
 * @param {number} z depth [m], ≥ 0
 * @param {'liao-whitman'|'blake'} [form='liao-whitman']
 * @returns {number} rd [-]
 */
export function stressReduction(z, form = 'liao-whitman') {
  need('z', z, 0);
  if (form === 'blake') {
    const r = Math.sqrt(z);
    return (1 - 0.4113 * r + 0.04052 * z + 0.001753 * z * r) / (1 - 0.4177 * r + 0.05729 * z - 0.006205 * z * r + 0.00121 * z * z);
  }
  if (form !== 'liao-whitman') throw new RangeError(`form must be 'liao-whitman' or 'blake' (got ${form})`);
  if (z <= 9.15) return 1 - 0.00765 * z;
  if (z <= 23) return 1.174 - 0.0267 * z;
  if (z <= 30) return 0.744 - 0.008 * z;
  return 0.5;
}

/**
 * Earthquake-induced cyclic stress ratio (Seed & Idriss 1971): CSR = 0.65 (amax/g)(σv/σ′v) rd.
 * @param {object} o
 * @param {number} o.amax peak horizontal ground-surface acceleration [m/s^2]
 * @param {number} o.sigmaV total vertical stress [Pa]
 * @param {number} o.sigmaVeff effective vertical stress [Pa]
 * @param {number} o.z depth [m]
 * @param {number} [o.rd] stress-reduction coefficient [-] (default Liao & Whitman)
 * @returns {number} CSR [-]
 */
export function cyclicStressRatio({ amax, sigmaV, sigmaVeff, z, rd }) {
  need('amax', amax, 0);
  need('sigmaV', sigmaV, 0);
  need('sigmaVeff', sigmaVeff, 1e-9);
  const r = rd !== undefined ? rd : stressReduction(z);
  return 0.65 * (amax / G0) * (sigmaV / sigmaVeff) * r;
}

/**
 * Overburden correction CN = (Pa/σ′v)^0.5 ≤ 1.7 (Liao & Whitman 1986; Youd et al. 2001 eq. 9).
 * @param {number} sigmaVeff effective vertical stress [Pa]
 * @returns {number} CN [-]
 */
export function overburdenCN(sigmaVeff) {
  need('sigmaVeff', sigmaVeff, 1e-9);
  return Math.min(Math.sqrt(PATM / sigmaVeff), 1.7);
}

/**
 * Clean-sand equivalent blow count (Youd et al. 2001 eqs 5–7): (N1)60cs = α + β (N1)60,
 *   FC ≤ 5 %: α = 0, β = 1;  5 < FC < 35 %: α = exp(1.76 − 190/FC²), β = 0.99 + FC^1.5/1000;
 *   FC ≥ 35 %: α = 5, β = 1.2.
 * @param {number} N160 (N1)60 [blows/0.3 m]
 * @param {number} FC fines content [%]
 * @returns {{N160cs:number,alpha:number,beta:number}}
 */
export function finesCorrectionSpt(N160, FC) {
  need('N160', N160, 0);
  need('FC', FC, 0, 100);
  let alpha = 0;
  let beta = 1;
  if (FC >= 35) { alpha = 5; beta = 1.2; } else if (FC > 5) { alpha = Math.exp(1.76 - 190 / (FC * FC)); beta = 0.99 + FC ** 1.5 / 1000; }
  return { N160cs: alpha + beta * N160, alpha, beta };
}

/**
 * Cyclic resistance ratio for M = 7.5 from the clean-sand SPT base curve (Rauch's fit, Youd et al. 2001
 * eq. 4):  CRR7.5 = 1/(34 − N) + N/135 + 50/(10 N + 45)² − 1/200,  valid for N = (N1)60cs < 30.
 * For N ≥ 30 clean granular soils are classed as non-liquefiable (Infinity returned).
 * @param {number} N160cs (N1)60cs [blows/0.3 m]
 * @returns {number} CRR7.5 [-]
 */
export function crrSpt(N160cs) {
  need('N160cs', N160cs, 0);
  if (N160cs >= 30) return Infinity;
  const N = N160cs;
  return 1 / (34 - N) + N / 135 + 50 / (10 * N + 45) ** 2 - 1 / 200;
}

/**
 * CPT soil behaviour type index and clean-sand equivalent tip resistance (Robertson & Wride 1998; Youd et
 * al. 2001 eqs 12–19):
 *   Q = [(qc − σv0)/Pa] (Pa/σ′v0)^n,  F = 100 fs/(qc − σv0) [%],
 *   Ic = √[(3.47 − log Q)² + (1.22 + log F)²];
 *   n = 1 first; if Ic > 2.6 the soil is clay-like (not liquefiable by this criterion); otherwise n = 0.5 is
 *   used and, if that gives Ic > 2.6, n = 0.75;
 *   qc1N = CQ qc/Pa, CQ = (Pa/σ′v0)^n ≤ 1.7;
 *   Kc = 1 (Ic ≤ 1.64), else −0.403 Ic⁴ + 5.581 Ic³ − 21.63 Ic² + 33.75 Ic − 17.88;  qc1Ncs = Kc qc1N.
 * @param {object} o
 * @param {number} o.qc cone tip resistance [Pa]
 * @param {number} o.fs sleeve friction [Pa]
 * @param {number} o.sigmaV total vertical stress [Pa]
 * @param {number} o.sigmaVeff effective vertical stress [Pa]
 * @param {number} [o.pa=101325] atmospheric pressure [Pa]
 * @returns {{Ic:number,n:number,Q:number,F:number,qc1N:number,Kc:number,qc1Ncs:number,clayLike:boolean}}
 */
export function cptNormalise({ qc, fs, sigmaV, sigmaVeff, pa = PATM }) {
  need('qc', qc, 1e-9);
  need('fs', fs, 1e-12);
  need('sigmaV', sigmaV, 0);
  need('sigmaVeff', sigmaVeff, 1e-9);
  if (qc <= sigmaV) throw new RangeError('qc must exceed the total vertical stress');
  const F = (100 * fs) / (qc - sigmaV);
  const icOf = (n) => {
    const Q = ((qc - sigmaV) / pa) * (pa / sigmaVeff) ** n;
    return { Q, Ic: Math.sqrt((3.47 - Math.log10(Q)) ** 2 + (1.22 + Math.log10(F)) ** 2) };
  };
  let n = 1;
  let r = icOf(1);
  let clayLike = false;
  if (r.Ic > 2.6) clayLike = true;
  else {
    n = 0.5;
    r = icOf(0.5);
    if (r.Ic > 2.6) { n = 0.75; r = icOf(0.75); }
  }
  const qc1N = Math.min((pa / sigmaVeff) ** n, 1.7) * (qc / pa);
  const Ic = r.Ic;
  const Kc = Ic <= 1.64 ? 1 : -0.403 * Ic ** 4 + 5.581 * Ic ** 3 - 21.63 * Ic ** 2 + 33.75 * Ic - 17.88;
  return { Ic, n, Q: r.Q, F, qc1N, Kc, qc1Ncs: Kc * qc1N, clayLike };
}

/**
 * CRR7.5 from the clean-sand CPT curve (Robertson & Wride 1998; Youd et al. 2001 eq. 11):
 *   qc1Ncs < 50:        CRR = 0.833 (qc1Ncs/1000) + 0.05;
 *   50 ≤ qc1Ncs < 160:  CRR = 93 (qc1Ncs/1000)³ + 0.08.
 * Above 160 the soil is classed as non-liquefiable (Infinity returned).
 * @param {number} qc1Ncs [-]
 * @returns {number} CRR7.5 [-]
 */
export function crrCpt(qc1Ncs) {
  need('qc1Ncs', qc1Ncs, 0);
  if (qc1Ncs >= 160) return Infinity;
  const x = qc1Ncs / 1000;
  return qc1Ncs < 50 ? 0.833 * x + 0.05 : 93 * x ** 3 + 0.08;
}

/**
 * CRR7.5 from overburden-corrected shear-wave velocity (Andrus & Stokoe 2000; Youd et al. 2001 eq. 22):
 *   Vs1 = Vs (Pa/σ′v)^0.25,
 *   CRR = a (Vs1/100)² + b [1/(V*s1 − Vs1) − 1/V*s1],  a = 0.022, b = 2.8,
 *   V*s1 = 215 m/s (FC ≤ 5 %), 215 − 0.5 (FC − 5) (5 < FC < 35 %), 200 m/s (FC ≥ 35 %).
 * Vs1 ≥ V*s1 is non-liquefiable (Infinity returned).
 * @param {object} o
 * @param {number} o.Vs shear-wave velocity [m/s]
 * @param {number} o.sigmaVeff effective vertical stress [Pa]
 * @param {number} [o.FC=0] fines content [%]
 * @returns {{CRR:number,Vs1:number,Vs1Limit:number}}
 */
export function crrVs({ Vs, sigmaVeff, FC = 0 }) {
  need('Vs', Vs, 1e-9);
  need('sigmaVeff', sigmaVeff, 1e-9);
  need('FC', FC, 0, 100);
  const Vs1 = Vs * (PATM / sigmaVeff) ** 0.25;
  const lim = FC <= 5 ? 215 : FC >= 35 ? 200 : 215 - 0.5 * (FC - 5);
  const CRR = Vs1 >= lim ? Infinity : 0.022 * (Vs1 / 100) ** 2 + 2.8 * (1 / (lim - Vs1) - 1 / lim);
  return { CRR, Vs1, Vs1Limit: lim };
}

/**
 * Magnitude scaling factor, revised Idriss form (Youd et al. 2001 eq. 24): MSF = 10^2.24 / Mw^2.56.
 * @param {number} Mw moment magnitude [-], 5 ≤ Mw ≤ 9
 * @returns {number} MSF [-]
 */
export function magnitudeScalingFactor(Mw) {
  need('Mw', Mw, 4, 9.5);
  return 10 ** 2.24 / Mw ** 2.56;
}

/**
 * Overburden correction Kσ = (σ′v/Pa)^(f − 1) (Hynes & Olsen 1999; Youd et al. 2001 eq. 31), with
 * f = 0.7–0.8 for Dr = 40–60 % and 0.6–0.7 for Dr = 60–80 %. Applied for σ′v > Pa only (Kσ ≤ 1).
 * @param {number} sigmaVeff [Pa]
 * @param {number} [f=0.7] exponent [-]
 * @returns {number} Kσ [-]
 */
export function kSigma(sigmaVeff, f = 0.7) {
  need('sigmaVeff', sigmaVeff, 1e-9);
  need('f', f, 0.5, 1);
  return Math.min((sigmaVeff / PATM) ** (f - 1), 1);
}

/**
 * Factor of safety against liquefaction: FS = (CRR7.5 / CSR) · MSF · Kσ (Youd et al. 2001 eqs 23, 30).
 * @param {object} o CRR75 [-], CSR [-], Mw [-], Ksigma [-] (default 1)
 * @returns {number} FS [-]
 */
export function liquefactionFs({ CRR75, CSR, Mw = 7.5, Ksigma = 1 }) {
  need('CSR', CSR, 1e-12);
  return (CRR75 * magnitudeScalingFactor(Mw) * Ksigma) / CSR;
}

/**
 * Factor-of-safety profile from SPT data through a layered column. For a submerged seabed set
 * waterTableDepth = 0: free water above the mudline carries no shear, so σv in the CSR is the total stress
 * of the soil column only (measured from the mudline, hydrostatic pressure of the water column excluded).
 * @param {object} o
 * @param {{z:number,N160:number,FC?:number,gamma:number}[]} o.points test depths [m] in increasing order,
 *   (N1)60, fines [%], and total unit weight of the soil between the previous depth and this one [N/m^3]
 * @param {number} o.waterTableDepth [m] (0 for a submerged seabed)
 * @param {number} o.amax [m/s^2]
 * @param {number} o.Mw [-]
 * @param {number} [o.f=0.7] Kσ exponent
 * @param {number} [o.gammaW=9810] [N/m^3]
 * @returns {{z:number,sigmaV:number,sigmaVeff:number,rd:number,CSR:number,N160cs:number,CRR75:number,FS:number}[]}
 */
export function sptProfile({ points, waterTableDepth, amax, Mw, f = 0.7, gammaW = 9810 }) {
  if (!Array.isArray(points) || points.length === 0) throw new RangeError('points must be a non-empty array');
  need('waterTableDepth', waterTableDepth, 0);
  let zPrev = 0;
  let sv = 0;
  const out = [];
  for (const p of points) {
    need('point.z', p.z, zPrev);
    need('point.gamma', p.gamma, 1e-9);
    sv += p.gamma * (p.z - zPrev);
    zPrev = p.z;
    const u = Math.max(p.z - waterTableDepth, 0) * gammaW;
    const sve = sv - u;
    const rd = stressReduction(p.z);
    const CSR = cyclicStressRatio({ amax, sigmaV: sv, sigmaVeff: sve, z: p.z, rd });
    const { N160cs } = finesCorrectionSpt(p.N160, p.FC || 0);
    const CRR75 = crrSpt(N160cs);
    out.push({ z: p.z, sigmaV: sv, sigmaVeff: sve, rd, CSR, N160cs, CRR75, FS: liquefactionFs({ CRR75, CSR, Mw, Ksigma: kSigma(sve, f) }) });
  }
  return out;
}

/**
 * Liquefaction potential index (Iwasaki et al. 1982): LPI = ∫0^20 F(z) w(z) dz, F = 1 − FS for FS < 1
 * (else 0), w = 10 − 0.5 z. Trapezoidal integration between the given depths; the profile is held constant
 * above the first and below the last point. LPI > 15: severe; 5–15: high; 0–5: low.
 * @param {number[]} z depths [m], increasing
 * @param {number[]} fs factors of safety [-]
 * @returns {number} LPI [-] (0–100)
 */
export function liquefactionPotentialIndex(z, fs) {
  if (!z || !fs || z.length !== fs.length || z.length === 0) throw new RangeError('z and fs must be arrays of equal non-zero length');
  const zz = [0, ...z, 20];
  const ff = [fs[0], ...fs, fs[fs.length - 1]];
  const g = (i) => (Math.min(ff[i], 1e9) < 1 ? 1 - Math.max(ff[i], 0) : 0);
  let lpi = 0;
  for (let i = 0; i < zz.length - 1; i++) {
    const a = Math.min(zz[i], 20);
    const b = Math.min(zz[i + 1], 20);
    if (b <= a) continue;
    // F is interpolated linearly, w is linear: integrate the product exactly (quadratic → Simpson)
    const fa = g(i);
    const fb = g(i + 1);
    const wa = 10 - 0.5 * a;
    const wb = 10 - 0.5 * b;
    lpi += ((b - a) / 6) * (fa * wa + 4 * (0.5 * (fa + fb)) * (0.5 * (wa + wb)) + fb * wb);
  }
  return lpi;
}

/**
 * Maximum shear strain and post-liquefaction volumetric strain (Idriss & Boulanger 2008 eqs 86–96, the
 * closed-form representation by Yoshimine et al. 2006 of the Ishihara & Yoshimine 1992 charts):
 *   γlim = 1.859 (1.1 − √(N/46))³ ≥ 0,   Fα = 0.032 + 0.69 √N − 0.13 N   (N = (N1)60cs ≥ 7),
 *   γmax = 0 (FS ≥ 2);  min{γlim, 0.035 (2 − FS)(1 − Fα)/(FS − Fα)} (Fα < FS < 2);  γlim (FS ≤ Fα),
 *   εv = 1.5 exp(−0.369 √N) · min(0.08, γmax).
 * Strains are returned as fractions.
 * @param {number} N160cs clean-sand blow count [-]
 * @param {number} FS factor of safety against triggering [-]
 * @returns {{gammaMax:number,epsV:number,gammaLim:number,Falpha:number}} [-]
 */
export function postLiquefactionStrain(N160cs, FS) {
  need('N160cs', N160cs, 0);
  need('FS', FS, 0);
  const N = Math.max(N160cs, 7);
  const r = Math.sqrt(N);
  const gammaLim = Math.max(1.859 * (1.1 - Math.sqrt(N / 46)) ** 3, 0);
  const Falpha = 0.032 + 0.69 * r - 0.13 * N;
  let gammaMax;
  if (FS >= 2) gammaMax = 0;
  else if (FS > Falpha) gammaMax = Math.min(gammaLim, (0.035 * (2 - FS) * (1 - Falpha)) / (FS - Falpha));
  else gammaMax = gammaLim;
  const epsV = 1.5 * Math.exp(-0.369 * r) * Math.min(0.08, gammaMax);
  return { gammaMax, epsV, gammaLim, Falpha };
}

/**
 * Post-liquefaction free-field settlement and lateral displacement index by integrating the strains of
 * postLiquefactionStrain over the profile (each value represents a layer of the given thickness):
 *   S = Σ εv Δz,   LDI = Σ γmax Δz.
 * @param {{thickness:number,N160cs:number,FS:number}[]} layers thickness [m]
 * @returns {{settlement:number,LDI:number}} [m, m]
 */
export function postLiquefactionSettlement(layers) {
  if (!Array.isArray(layers)) throw new RangeError('layers must be an array');
  let s = 0;
  let ldi = 0;
  for (const ly of layers) {
    need('layer.thickness', ly.thickness, 0);
    const r = postLiquefactionStrain(ly.N160cs, ly.FS);
    s += r.epsV * ly.thickness;
    ldi += r.gammaMax * ly.thickness;
  }
  return { settlement: s, LDI: ldi };
}

/**
 * Residual (liquefied) shear strength ratio, Idriss & Boulanger (2008) eqs 99–100, N = (N1)60cs-Sr:
 *   void redistribution significant:    Sr/σ′v = exp[N/16 + ((N − 16)/21.2)³ − 3.0] ≤ tan φ′;
 *   void redistribution negligible:     the above × [1 + exp(N/2.4 − 6.6)] ≤ tan φ′.
 * The fines adjustment ΔN (Seed 1987: 1, 2, 4, 5 for FC = 10, 25, 50, 75 %) is added by the caller.
 * @param {object} o
 * @param {number} o.N160csSr [-]
 * @param {number} o.phi effective friction angle for the cap [rad]
 * @param {boolean} [o.voidRedistribution=true]
 * @returns {number} Sr/σ′v0 [-]
 */
export function residualStrengthRatio({ N160csSr, phi, voidRedistribution = true }) {
  need('N160csSr', N160csSr, 0);
  need('phi', phi, 0, Math.PI / 2 - 1e-9);
  const N = N160csSr;
  let r = Math.exp(N / 16 + ((N - 16) / 21.2) ** 3 - 3.0);
  if (!voidRedistribution) r *= 1 + Math.exp(N / 2.4 - 6.6);
  return Math.min(r, Math.tan(phi));
}

/**
 * Wave-induced cyclic shear stress in a deep, homogeneous, elastic seabed (plane strain, infinite
 * thickness): amplitude of the bottom pressure p0 = γw H / (2 cosh kd) and
 *   τxz(z) = p0 k z e^(−kz),   cyclic stress ratio CSR = τxz / σ′v0 = p0 k e^(−kz) / γ′.
 * (Yamamoto 1978 / Madsen 1978 limit; used by Seed & Rahman 1978 and Ishihara & Yamazaki 1984.)
 * @param {object} o
 * @param {number} o.H wave height [m]
 * @param {number} o.k wave number 2π/L [rad/m]
 * @param {number} o.d water depth [m]
 * @param {number} o.z depth below the seabed [m]
 * @param {number} o.gammaSub submerged unit weight [N/m^3]
 * @param {number} [o.gammaW=10055] unit weight of sea water [N/m^3]
 * @returns {{p0:number,tau:number,CSR:number}} [Pa, Pa, -]
 */
export function waveCyclicStress({ H, k, d, z, gammaSub, gammaW = 10055 }) {
  need('H', H, 0);
  need('k', k, 1e-12);
  need('d', d, 1e-9);
  need('z', z, 0);
  need('gammaSub', gammaSub, 1e-9);
  const p0 = (gammaW * H) / (2 * Math.cosh(k * d));
  const e = Math.exp(-k * z);
  return { p0, tau: p0 * k * z * e, CSR: (p0 * k * e) / gammaSub };
}

/**
 * Undrained residual pore-pressure build-up under uniform cyclic loading (Seed & Rahman 1978, after
 * Seed, Martin & Lysmer 1975):  ru = u/σ′v0 = (2/π) asin[(N/Nl)^(1/(2θ))],  θ ≈ 0.7.
 * The number of cycles to liquefaction Nl follows from the laboratory cyclic-strength curve
 * τ/σ′v0 = a Nl^(−b), i.e. Nl = (a/CSR)^(1/b) – a and b are soil-specific arguments.
 * Drainage during the storm is neglected (conservative).
 * @param {object} o
 * @param {number} o.CSR applied cyclic stress ratio [-]
 * @param {number} o.N number of wave cycles [-]
 * @param {number} o.a cyclic-strength coefficient [-]
 * @param {number} o.b cyclic-strength exponent [-], > 0
 * @param {number} [o.theta=0.7] shape parameter [-]
 * @returns {{ru:number,Nl:number,liquefied:boolean}}
 */
export function seedRahmanPorePressure({ CSR, N, a, b, theta = 0.7 }) {
  need('CSR', CSR, 0);
  need('N', N, 0);
  need('a', a, 1e-12);
  need('b', b, 1e-6);
  need('theta', theta, 0.1, 2);
  if (CSR === 0) return { ru: 0, Nl: Infinity, liquefied: false };
  const Nl = (a / CSR) ** (1 / b);
  const ratio = Math.min(N / Nl, 1);
  return { ru: (2 / Math.PI) * Math.asin(ratio ** (1 / (2 * theta))), Nl, liquefied: N >= Nl };
}

/**
 * Momentary (instantaneous) liquefaction under a wave trough (Zen & Yamazaki 1990): liquefaction to depth z
 * when the upward seepage force exceeds the submerged weight,
 *     p(z) − pb ≥ γ′ z,
 * where pb is the wave pressure on the seabed surface (−p0 under the trough) and p(z) the pore pressure at
 * depth z at the same instant. For the incompressible infinite-thickness solution p(z) = pb e^(−kz), giving
 * an excess of p0 (1 − e^(−kz)); for partly saturated beds supply the actual amplitude ratio p(z)/pb as
 * `attenuation` (e.g. from seepage.js).
 * @param {object} o
 * @param {number} o.p0 amplitude of the wave pressure on the seabed [Pa]
 * @param {number} o.z depth [m]
 * @param {number} o.gammaSub [N/m^3]
 * @param {number} [o.k] wave number [rad/m] (for the default attenuation)
 * @param {number} [o.attenuation] p(z)/pb at depth z [-]
 * @returns {{excess:number,resisting:number,liquefied:boolean,ratio:number}} excess = p(z) − pb [Pa];
 *   resisting = γ′z [Pa]; ratio = excess/resisting [-]
 */
export function momentaryLiquefaction({ p0, z, gammaSub, k, attenuation }) {
  need('p0', p0, 0);
  need('z', z, 0);
  need('gammaSub', gammaSub, 1e-9);
  const att = attenuation !== undefined ? need('attenuation', attenuation, -1, 1) : Math.exp(-need('k', k, 1e-12) * z);
  const excess = p0 * (1 - att);
  const resisting = gammaSub * z;
  return { excess, resisting, liquefied: z > 0 && excess >= resisting, ratio: resisting > 0 ? excess / resisting : 0 };
}
