/**
 * Shallow-foundation bearing capacity, sliding and overturning.
 *
 * SI units throughout: lengths m, forces N, stresses and strengths Pa, unit weights N/m^3, angles rad.
 * A strip footing is described by L = Infinity; forces are then per metre run (N/m).
 *
 * References
 *   Terzaghi K. (1943) Theoretical Soil Mechanics, Wiley.
 *   Meyerhof G.G. (1963) Some recent research on the bearing capacity of foundations, Can. Geotech. J. 1(1).
 *   Brinch Hansen J. (1970) A revised and extended formula for bearing capacity, DGI Bulletin 28.
 *   Vesic A.S. (1973) Analysis of ultimate loads of shallow foundations, JSMFD ASCE 99(SM1);
 *   Vesic A.S. (1975) Bearing capacity of shallow foundations, Foundation Engineering Handbook, ch. 3.
 *   Skempton A.W. (1951) The bearing capacity of clays, Building Research Congress.
 *   Davis E.H. & Booker J.R. (1973) Geotechnique 23(4).
 *   Meyerhof G.G. & Hanna A.M. (1978) Can. Geotech. J. 15(4); Brown J.D. & Meyerhof G.G. (1969).
 *   API RP 2A-WSD (21st ed.) 6.13 and commentary; API RP 2GEO (2011) 7; ISO 19901-4:2016 7; DNV-RP-C212.
 */

const PI = Math.PI;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}
const PHI_MAX = (50 * PI) / 180;

/**
 * Prandtl–Reissner surcharge factor Nq = e^(π tan φ) tan²(π/4 + φ/2).
 * @param {number} phi friction angle [rad], 0 ≤ φ ≤ 50°
 * @returns {number} Nq [-]
 */
export function bearingNq(phi) {
  need('phi', phi, 0, PHI_MAX);
  const t = Math.tan(PI / 4 + phi / 2);
  return Math.exp(PI * Math.tan(phi)) * t * t;
}

/**
 * Prandtl cohesion factor Nc = (Nq − 1) cot φ, with the φ → 0 limit 2 + π.
 * @param {number} phi friction angle [rad], 0 ≤ φ ≤ 50°
 * @returns {number} Nc [-]
 */
export function bearingNc(phi) {
  need('phi', phi, 0, PHI_MAX);
  if (phi < 1e-9) return 2 + PI;
  return (bearingNq(phi) - 1) / Math.tan(phi);
}

/**
 * Self-weight factor Nγ by the named author.
 *   'meyerhof' (1963): (Nq − 1) tan(1.4 φ)
 *   'hansen'   (1970): 1.5 (Nq − 1) tan φ
 *   'vesic'    (1973): 2 (Nq + 1) tan φ
 * @param {number} phi friction angle [rad], 0 ≤ φ ≤ 50°
 * @param {'meyerhof'|'hansen'|'vesic'} [form='vesic']
 * @returns {number} Nγ [-]
 */
export function bearingNgamma(phi, form = 'vesic') {
  need('phi', phi, 0, PHI_MAX);
  if (phi < 1e-9) return 0;
  const Nq = bearingNq(phi);
  if (form === 'meyerhof') return (Nq - 1) * Math.tan(1.4 * phi);
  if (form === 'hansen') return 1.5 * (Nq - 1) * Math.tan(phi);
  if (form === 'vesic') return 2 * (Nq + 1) * Math.tan(phi);
  throw new RangeError(`form must be 'meyerhof', 'hansen' or 'vesic' (got ${form})`);
}

/**
 * Terzaghi (1943) factors for a rough base (general shear).
 *   Nq = e^(2(3π/4 − φ/2) tan φ) / (2 cos²(π/4 + φ/2)),  Nc = (Nq − 1) cot φ  (3π/2 + 1 = 5.71 at φ = 0).
 * Terzaghi's Nγ = ½ tan φ (Kpγ / cos²φ − 1) needs the passive coefficient Kpγ, which he gave only as a
 * chart. Nγ is therefore returned from the closed-form approximation of Coduto (2001, Foundation Design,
 * eq. 6.11), Nγ ≈ 2 (Nq + 1) tan φ / (1 + 0.4 sin 4φ), which follows the tabulated values of
 * Kumbhojkar (1993) to within about 10 % for 0–40° (e.g. 20.1 against 19.1 at 30°). Pass `Kpgamma` to use
 * Terzaghi's exact expression with a chart value instead.
 * @param {number} phi friction angle [rad], 0 ≤ φ ≤ 45°
 * @param {{Kpgamma?: number}} [opts] Kpγ [-] passive coefficient from Terzaghi's chart
 * @returns {{Nc:number,Nq:number,Ngamma:number}} factors [-]
 */
export function terzaghiFactors(phi, opts = {}) {
  need('phi', phi, 0, PI / 4);
  const c = Math.cos(PI / 4 + phi / 2);
  const Nq = Math.exp(2 * (0.75 * PI - phi / 2) * Math.tan(phi)) / (2 * c * c);
  const Nc = phi < 1e-9 ? 1.5 * PI + 1 : (Nq - 1) / Math.tan(phi);
  let Ngamma;
  if (opts.Kpgamma !== undefined) {
    need('Kpgamma', opts.Kpgamma, 0);
    Ngamma = 0.5 * Math.tan(phi) * (opts.Kpgamma / Math.cos(phi) ** 2 - 1);
  } else {
    Ngamma = phi < 1e-9 ? 0 : (2 * (Nq + 1) * Math.tan(phi)) / (1 + 0.4 * Math.sin(4 * phi));
  }
  return { Nc, Nq, Ngamma };
}

/**
 * Terzaghi (1943) ultimate bearing pressure for vertical central load:
 *   strip    qult = c Nc + q Nq + 0.5 γ B Nγ
 *   square   qult = 1.3 c Nc + q Nq + 0.4 γ B Nγ
 *   circular qult = 1.3 c Nc + q Nq + 0.3 γ B Nγ   (B = diameter)
 * Valid for D ≤ B, general shear failure.
 * @param {object} o
 * @param {'strip'|'square'|'circular'} o.shape
 * @param {number} o.B width or diameter [m]
 * @param {number} [o.c=0] cohesion or undrained strength [Pa]
 * @param {number} [o.phi=0] friction angle [rad]
 * @param {number} o.gamma effective unit weight below the base [N/m^3]
 * @param {number} [o.q=0] effective overburden at base level [Pa]
 * @param {number} [o.Kpgamma] optional Terzaghi chart value, see terzaghiFactors
 * @returns {{qult:number,Nc:number,Nq:number,Ngamma:number}} qult gross ultimate pressure [Pa]
 */
export function terzaghiBearing({ shape, B, c = 0, phi = 0, gamma, q = 0, Kpgamma }) {
  need('B', B, 1e-9);
  need('gamma', gamma, 0);
  need('c', c, 0);
  need('q', q, 0);
  const f = terzaghiFactors(phi, { Kpgamma });
  const k = { strip: [1, 0.5], square: [1.3, 0.4], circular: [1.3, 0.3] }[shape];
  if (!k) throw new RangeError(`shape must be 'strip', 'square' or 'circular' (got ${shape})`);
  return { qult: k[0] * c * f.Nc + q * f.Nq + k[1] * gamma * B * f.Ngamma, ...f };
}

/**
 * Skempton (1951) Nc for undrained clay with embedment and shape.
 * Skempton's chart runs from 5.14 (strip) and 6.2 (square/circle) at the surface to 7.5 and 9.0 at depth,
 * with Nc(rectangle) = (0.84 + 0.16 B/L) Nc(square). Two documented fits are offered:
 *   'chart'  (default): Nc = min[5.14 (1 + 0.2 D/B), 7.5] (1 + 0.2 B/L) – reproduces the four end points of
 *            the chart (5.14, 6.17, 7.5, 9.0); intermediate depths within about 5 % of the curve.
 *   'simple': Nc = 5 (1 + 0.2 D/B)(1 + 0.2 B/L), D/B ≤ 2.5 (the usual textbook approximation).
 * @param {object} o
 * @param {number} o.B width [m]
 * @param {number} [o.L=Infinity] length [m] (≥ B; Infinity = strip)
 * @param {number} [o.D=0] embedment depth [m]
 * @param {'chart'|'simple'} [o.form='chart']
 * @returns {number} Nc [-]
 */
export function skemptonNc({ B, L = Infinity, D = 0, form = 'chart' }) {
  need('B', B, 1e-9);
  need('L', L, B);
  need('D', D, 0);
  const bl = B / L;
  if (form === 'simple') return 5 * (1 + 0.2 * Math.min(D / B, 2.5)) * (1 + 0.2 * bl);
  if (form !== 'chart') throw new RangeError(`form must be 'chart' or 'simple' (got ${form})`);
  return Math.min((2 + PI) * (1 + (0.2 * D) / B), 7.5) * (1 + 0.2 * bl);
}

/**
 * Meyerhof effective area of an eccentrically loaded rectangular footing: B' = B − 2 eB, L' = L − 2 eL.
 * The returned B' is always the shorter effective side.
 * @param {object} o
 * @param {number} o.B width [m]
 * @param {number} o.L length [m] (Infinity for a strip; area is then per metre)
 * @param {number} [o.eB=0] eccentricity in the B direction, M/V [m], |eB| < B/2
 * @param {number} [o.eL=0] eccentricity in the L direction [m], |eL| < L/2
 * @returns {{Beff:number,Leff:number,Aeff:number,swapped:boolean}} [m, m, m^2 (m^2/m for a strip)];
 *   swapped is true when the effective short side lies along the original L
 */
export function effectiveAreaRect({ B, L, eB = 0, eL = 0 }) {
  need('B', B, 1e-9);
  need('L', L, 1e-9);
  need('eB', Math.abs(eB), 0, B / 2);
  if (!Number.isFinite(L)) {
    const b = B - 2 * Math.abs(eB);
    return { Beff: b, Leff: Infinity, Aeff: b, swapped: false };
  }
  need('eL', Math.abs(eL), 0, L / 2);
  const b = B - 2 * Math.abs(eB);
  const l = L - 2 * Math.abs(eL);
  const swapped = b > l;
  return { Beff: Math.min(b, l), Leff: Math.max(b, l), Aeff: b * l, swapped };
}

/**
 * Effective area of an eccentrically loaded circular footing and its equivalent rectangle
 * (API RP 2A-WSD C6.13 / API RP 2GEO 7.3; DNV-RP-C212 / DNV-ST-0126 App. G):
 *   A' = 2 [R² acos(e/R) − e √(R² − e²)],  L'/B' = √((R + e)/(R − e)),  A' = B' L'.
 * @param {object} o
 * @param {number} o.R radius [m]
 * @param {number} o.e load eccentricity M/V [m], 0 ≤ e < R
 * @returns {{Aeff:number,Beff:number,Leff:number,be:number,le:number}} [m^2, m, m, m, m]; be and le are the
 *   minor and major axes of the lens-shaped effective area
 */
export function effectiveAreaCircle({ R, e }) {
  need('R', R, 1e-9);
  need('e', e, 0, R * (1 - 1e-12));
  const Aeff = 2 * (R * R * Math.acos(e / R) - e * Math.sqrt(R * R - e * e));
  const be = 2 * (R - e);
  const le = 2 * R * Math.sqrt(1 - (1 - be / (2 * R)) ** 2);
  const Leff = Math.sqrt((Aeff * le) / be);
  const Beff = (Leff * be) / le;
  return { Aeff, Beff, Leff, be, le };
}

function meyerhofSet(phi, bl, db, H, V) {
  const Kp = Math.tan(PI / 4 + phi / 2) ** 2;
  // Meyerhof states sq = sγ = dq = dγ = 1 at φ = 0 and the Kp expressions for φ > 10°; a linear blend is
  // used in between (assumption).
  const w = Math.min(phi / ((10 * PI) / 180), 1);
  const sc = 1 + 0.2 * Kp * bl;
  const sq = 1 + w * 0.1 * Kp * bl;
  const dc = 1 + 0.2 * Math.sqrt(Kp) * db;
  const dq = 1 + w * 0.1 * Math.sqrt(Kp) * db;
  const th = Math.atan2(H, V);
  const ic = (1 - th / (PI / 2)) ** 2;
  const ig = phi < 1e-9 ? (H > 0 ? 0 : 1) : th >= phi ? 0 : (1 - th / phi) ** 2;
  return { sc, sq, sg: sq, dc, dq, dg: dq, ic, iq: ic, ig, bc: 1, bq: 1, bg: 1, gc: 1, gq: 1, gg: 1 };
}

/**
 * General bearing-capacity equation
 *   qult = c Nc sc dc ic bc gc + q Nq sq dq iq bq gq + 0.5 γ B' Nγ sγ dγ iγ bγ gγ
 * on the Meyerhof effective area, with the factor sets of
 *   'meyerhof' (1963): shape, depth and inclination factors (base/ground inclination not covered);
 *   'hansen'   Brinch Hansen (1970): shape, depth, inclination (exponents 5), base tilt η and ground slope β;
 *              for φ = 0 the additive form qult = (2 + π) su (1 + s'c + d'c − i'c − b'c − g'c) + q;
 *   'vesic'    Vesic (1973, 1975): idem with exponent m = (2 + B/L)/(1 + B/L) (load along B),
 *              (2 + L/B)/(1 + L/B) (load along L), combined as mL cos²θ + mB sin²θ. For φ = 0 on sloping
 *              ground the self-weight term uses Nγ = −2 sin β, and bc = 1 − 2η/(π + 2), gc = 1 − 2β/(π + 2).
 * Depth factors use k = D/B for D ≤ B and atan(D/B) beyond (Hansen). Set `useDepthFactors: false` when the
 * strength of the overburden cannot be relied on.
 *
 * Validity: homogeneous soil, 0 ≤ φ ≤ 50°, η + β ≤ 90°, β ≤ φ for φ > 0, H less than the sliding resistance.
 * @param {object} o
 * @param {'meyerhof'|'hansen'|'vesic'} [o.method='vesic']
 * @param {number} o.B footing width [m]
 * @param {number} [o.L=Infinity] footing length [m]
 * @param {number} [o.D=0] embedment depth of the base [m]
 * @param {number} [o.c=0] cohesion c' or undrained strength su [Pa]
 * @param {number} [o.phi=0] friction angle [rad]
 * @param {number} o.gamma effective unit weight below the base [N/m^3]
 * @param {number} [o.q] effective overburden at base level [Pa] (default γ D)
 * @param {number} [o.V] vertical load [N] (needed for Hansen/Vesic inclination when φ > 0)
 * @param {number} [o.HB=0] horizontal load along B [N]
 * @param {number} [o.HL=0] horizontal load along L [N]
 * @param {number} [o.eB=0] eccentricity along B [m]
 * @param {number} [o.eL=0] eccentricity along L [m]
 * @param {number} [o.baseTilt=0] base inclination η from horizontal [rad]
 * @param {number} [o.groundSlope=0] ground slope β [rad]
 * @param {number} [o.adhesion] base adhesion ca [Pa] (default c)
 * @param {boolean} [o.useDepthFactors=true]
 * @param {number} [o.alpha1=5] Hansen exponent on iq
 * @param {number} [o.alpha2=5] Hansen exponent on iγ
 * @returns {{qult:number,Qult:number,Aeff:number,Beff:number,Leff:number,Nc:number,Nq:number,Ngamma:number,
 *   factors:object,terms:{cohesion:number,surcharge:number,selfWeight:number}}}
 *   qult gross ultimate pressure on the effective area [Pa]; Qult = qult·A' [N] ([N/m] for a strip)
 */
export function bearingCapacity(o) {
  const {
    method = 'vesic', B, L = Infinity, D = 0, c = 0, phi = 0, gamma, V, HB = 0, HL = 0, eB = 0, eL = 0,
    baseTilt = 0, groundSlope = 0, useDepthFactors = true, alpha1 = 5, alpha2 = 5,
  } = o;
  need('B', B, 1e-9);
  need('L', L, 1e-9);
  need('D', D, 0);
  need('c', c, 0);
  need('phi', phi, 0, PHI_MAX);
  need('gamma', gamma, 0);
  need('baseTilt', baseTilt, 0, PI / 4);
  need('groundSlope', groundSlope, 0, PI / 4);
  const q = o.q !== undefined ? need('q', o.q, 0) : gamma * D;
  const ca = o.adhesion !== undefined ? need('adhesion', o.adhesion, 0) : c;
  const ea = effectiveAreaRect({ B, L, eB, eL });
  const { Beff, Leff, Aeff } = ea;
  const bl = Number.isFinite(Leff) ? Beff / Leff : 0;
  const H = Math.hypot(HB, HL);
  const und = phi < 1e-9;
  if (H > 0 && !und && !(V > 0)) throw new RangeError('V must be > 0 when a horizontal load is given');
  const db = D / Beff;
  const k = useDepthFactors ? (db <= 1 ? db : Math.atan(db)) : 0;
  const Nq = bearingNq(phi);
  const Nc = bearingNc(phi);
  const tp = Math.tan(phi);
  const eta = baseTilt;
  const beta = groundSlope;
  let Ng;
  let f;
  let qult;
  let terms;

  if (method === 'meyerhof') {
    if (eta > 0 || beta > 0) throw new RangeError('baseTilt/groundSlope are not covered by the Meyerhof factors');
    Ng = bearingNgamma(phi, 'meyerhof');
    f = meyerhofSet(phi, bl, useDepthFactors ? db : 0, H, V > 0 ? V : 1);
  } else if (method === 'hansen') {
    Ng = bearingNgamma(phi, 'hansen');
    const etaDeg = (eta * 180) / PI;
    const betaDeg = (beta * 180) / PI;
    if (und) {
      if (H > Aeff * ca) throw new RangeError('H exceeds the base adhesion resistance A′·ca');
      const sc = 0.2 * bl;
      const dc = 0.4 * k;
      const ic = H > 0 ? 0.5 - 0.5 * Math.sqrt(1 - H / (Aeff * ca)) : 0;
      const bc = etaDeg / 147;
      const gc = betaDeg / 147;
      f = { sc, dc, ic, bc, gc, additive: true };
      const coh = Nc * c * (1 + sc + dc - ic - bc - gc);
      terms = { cohesion: coh, surcharge: q, selfWeight: 0 };
      qult = coh + q;
    } else {
      const den = H > 0 ? V + (Aeff * ca) / tp : 1;
      const r1 = H > 0 ? 1 - (0.5 * H) / den : 1;
      const r2 = H > 0 ? 1 - ((0.7 - etaDeg / 450) * H) / den : 1;
      if (r1 <= 0 || r2 <= 0) throw new RangeError('H is too large for the Hansen inclination factors');
      const iq = r1 ** alpha1;
      const ig = r2 ** alpha2;
      const ic = iq - (1 - iq) / (Nq - 1);
      const gq = (1 - 0.5 * Math.tan(beta)) ** 5;
      f = {
        sc: 1 + (Nq / Nc) * bl, sq: 1 + bl * Math.sin(phi), sg: Math.max(1 - 0.4 * bl, 0.6),
        dc: 1 + 0.4 * k, dq: 1 + 2 * tp * (1 - Math.sin(phi)) ** 2 * k, dg: 1,
        ic, iq, ig,
        bc: 1 - etaDeg / 147, bq: Math.exp(-2 * eta * tp), bg: Math.exp(-2.7 * eta * tp),
        gc: 1 - betaDeg / 147, gq, gg: gq,
      };
    }
  } else if (method === 'vesic') {
    Ng = bearingNgamma(phi, 'vesic');
    let m = 1;
    if (H > 0) {
      const mB = (2 + bl) / (1 + bl);
      const mL = Number.isFinite(Leff) ? (2 + 1 / bl) / (1 + 1 / bl) : 1;
      // θ is measured from the long (L') axis; the effective sides may have swapped relative to B, L.
      const hShort = ea.swapped ? HL : HB;
      const hLong = ea.swapped ? HB : HL;
      const th = Math.atan2(Math.abs(hShort), Math.abs(hLong));
      m = mL * Math.cos(th) ** 2 + mB * Math.sin(th) ** 2;
    }
    const sq = 1 + bl * tp;
    const sg = Math.max(1 - 0.4 * bl, 0.6);
    const dq = 1 + 2 * tp * (1 - Math.sin(phi)) ** 2 * k;
    const gq = (1 - Math.tan(beta)) ** 2;
    if (und) {
      if (H > 0 && m * H >= Aeff * ca * Nc) throw new RangeError('H is too large for the Vesic inclination factor');
      f = {
        sc: 1 + bl / Nc, sq, sg, dc: 1 + 0.4 * k, dq, dg: 1,
        ic: H > 0 ? 1 - (m * H) / (Aeff * ca * Nc) : 1, iq: 1, ig: 1,
        bc: 1 - (2 * eta) / (PI + 2), bq: 1, bg: 1,
        gc: 1 - (2 * beta) / (PI + 2), gq, gg: gq, m,
      };
      if (beta > 0) Ng = -2 * Math.sin(beta);
    } else {
      const r = H > 0 ? 1 - H / (V + (Aeff * ca) / tp) : 1;
      if (r <= 0) throw new RangeError('H is too large for the Vesic inclination factors');
      const iq = r ** m;
      const bq = (1 - eta * tp) ** 2;
      f = {
        sc: 1 + (Nq / Nc) * bl, sq, sg, dc: 1 + 0.4 * k, dq, dg: 1,
        ic: iq - (1 - iq) / (Nc * tp), iq, ig: r ** (m + 1),
        bc: bq - (1 - bq) / (Nc * tp), bq, bg: bq,
        gc: gq - (1 - gq) / (Nc * tp), gq, gg: gq, m,
      };
    }
  } else {
    throw new RangeError(`method must be 'meyerhof', 'hansen' or 'vesic' (got ${method})`);
  }

  if (qult === undefined) {
    const Bg = Number.isFinite(Beff) ? Beff : B;
    terms = {
      cohesion: c * Nc * f.sc * f.dc * f.ic * f.bc * f.gc,
      surcharge: q * Nq * f.sq * f.dq * f.iq * f.bq * f.gq,
      selfWeight: 0.5 * gamma * Bg * Ng * f.sg * f.dg * f.ig * f.bg * f.gg,
    };
    qult = terms.cohesion + terms.surcharge + terms.selfWeight;
  }
  return { qult, Qult: qult * Aeff, Aeff, Beff, Leff, Nc, Nq, Ngamma: Ng, factors: f, terms };
}

/**
 * Undrained bearing capacity for strength increasing linearly with depth, su = su0 + κ z
 * (Davis & Booker 1973):  qult = F [(2 + π) su0 + κ B / 4]  for a strip under vertical load.
 * The correction factor F is a function of κB/su0 and of base roughness that Davis & Booker gave as a chart
 * (F = 1 at κ = 0; roughly 1.0–1.5 for a rough base). It is not reproduced here: supply the chart value.
 * @param {object} o
 * @param {number} o.su0 strength at foundation level [Pa]
 * @param {number} o.kappa strength gradient [Pa/m]
 * @param {number} o.B width [m]
 * @param {number} [o.F=1] Davis & Booker correction factor [-]
 * @returns {{qult:number,kappaBoverSu0:number}} qult net ultimate pressure [Pa]; κB/su0 [-]
 */
export function davisBookerStrip({ su0, kappa, B, F = 1 }) {
  need('su0', su0, 0);
  need('kappa', kappa, 0);
  need('B', B, 1e-9);
  need('F', F, 0.5, 3);
  return { qult: F * ((2 + PI) * su0 + (kappa * B) / 4), kappaBoverSu0: su0 > 0 ? (kappa * B) / su0 : Infinity };
}

/**
 * Shape-factor coefficient scv for vertical loading on clay with strength increasing with depth, as
 * tabulated in ISO 19901-4:2016 / API RP 2GEO (after Salençon & Matar 1982):
 *   κB′/su0 = 0, 2, 4, 6, 8, 10  →  scv = 0.18, 0.00, −0.05, −0.07, −0.09, −0.10 (linear interpolation).
 * The table is quoted from memory of the standard and should be checked against the edition in force; the
 * capacity function accepts an explicit `scv` to override it.
 * @param {number} kappaBoverSu0 κB′/su0 [-], ≥ 0 (clamped at 10)
 * @returns {number} scv [-]
 */
export function isoScv(kappaBoverSu0) {
  need('kappaBoverSu0', kappaBoverSu0, 0);
  const xs = [0, 2, 4, 6, 8, 10];
  const ys = [0.18, 0.0, -0.05, -0.07, -0.09, -0.1];
  const x = Math.min(kappaBoverSu0, 10);
  const i = Math.min(Math.floor(x / 2), 4);
  return ys[i] + ((ys[i + 1] - ys[i]) * (x - xs[i])) / 2;
}

/**
 * Undrained capacity of a shallow foundation to ISO 19901-4:2016 7.2 / API RP 2GEO 7.3.1:
 *   Qd = F (Nc su0 + κ B′/4) Kc A′ / γm,   Nc = 2 + π,
 *   Kc = 1 + sc + dc − ic − bc − gc,
 *   sc = scv (1 − 2 ic) B′/L′,  dc = 0.3 (su,ave/su,2) atan(D/B′),  ic = 0.5 − 0.5 √(1 − H/(A′ su0)),
 *   bc = 2ν/(π + 2),  gc = 2β/(π + 2).
 * F is the Davis & Booker factor (argument, see davisBookerStrip). The overburden term is not included
 * (net capacity).
 * @param {object} o
 * @param {number} o.Beff effective width B′ [m]
 * @param {number} o.Leff effective length L′ [m] (Infinity for strip → result per metre)
 * @param {number} o.su0 undrained strength at base level [Pa]
 * @param {number} [o.kappa=0] strength gradient [Pa/m]
 * @param {number} [o.D=0] embedment [m]
 * @param {number} [o.H=0] horizontal load [N]
 * @param {number} [o.F=1] Davis & Booker factor [-]
 * @param {number} [o.scv] shape coefficient [-] (default isoScv(κB′/su0))
 * @param {number} [o.suRatio=1] su,ave/su,2: average strength above base over strength at base [-]
 * @param {number} [o.baseTilt=0] ν [rad]
 * @param {number} [o.groundSlope=0] β [rad]
 * @param {number} [o.materialFactor=1] γm [-]
 * @returns {{Qd:number,qd:number,Kc:number,sc:number,dc:number,ic:number,bc:number,gc:number,Aeff:number}}
 *   Qd design capacity [N]; qd [Pa]
 */
export function isoUndrainedCapacity(o) {
  const { Beff, Leff, su0, kappa = 0, D = 0, H = 0, F = 1, suRatio = 1, baseTilt = 0, groundSlope = 0, materialFactor = 1 } = o;
  need('Beff', Beff, 1e-9);
  need('Leff', Leff, Beff);
  need('su0', su0, 1e-9);
  need('kappa', kappa, 0);
  need('D', D, 0);
  need('H', H, 0);
  need('materialFactor', materialFactor, 1e-9);
  const strip = !Number.isFinite(Leff);
  const Aeff = strip ? Beff : Beff * Leff;
  if (H > Aeff * su0) throw new RangeError('H exceeds A′·su0');
  const scv = o.scv !== undefined ? o.scv : isoScv((kappa * Beff) / su0);
  const ic = 0.5 - 0.5 * Math.sqrt(1 - H / (Aeff * su0));
  const sc = strip ? 0 : (scv * (1 - 2 * ic) * Beff) / Leff;
  const dc = 0.3 * suRatio * Math.atan(D / Beff);
  const bc = (2 * baseTilt) / (PI + 2);
  const gc = (2 * groundSlope) / (PI + 2);
  const Kc = 1 + sc + dc - ic - bc - gc;
  const qd = (F * ((2 + PI) * su0 + (kappa * Beff) / 4) * Kc) / materialFactor;
  return { Qd: qd * Aeff, qd, Kc, sc, dc, ic, bc, gc, Aeff };
}

/**
 * Drained (or general c–φ) capacity with the correction factors of API RP 2A-WSD 6.13.2 and commentary:
 *   Q′ = (c′ Nc Kc + q Nq Kq + 0.5 γ′ B′ Nγ Kγ) A′,   Kx = ix · sx · dx · bx · gx,
 *   iq = [1 − H/(Q + B′L′c cot φ)]^m, iγ = [...]^(m+1), ic = iq − (1 − iq)/(Nc tan φ);
 *   sc = 1 + (B′/L′)(Nq/Nc), sq = 1 + (B′/L′) tan φ, sγ = 1 − 0.4 B′/L′;
 *   dq = 1 + 2 tan φ (1 − sin φ)² D/B′, dγ = 1, dc = dq − (1 − dq)/(Nc tan φ);
 *   bq = bγ = (1 − ν tan φ)², bc = bq − (1 − bq)/(Nc tan φ);
 *   gq = gγ = (1 − tan β)², gc = gq − (1 − gq)/(Nc tan φ).
 * Nγ: RP 2A-WSD uses 2 (Nq + 1) tan φ ('vesic'); RP 2GEO / ISO 19901-4 use 1.5 (Nq − 1) tan φ ('hansen').
 * @param {object} o
 * @param {number} o.Beff B′ [m]
 * @param {number} o.Leff L′ [m] (Infinity for strip)
 * @param {number} o.phi effective friction angle [rad], > 0
 * @param {number} [o.c=0] effective cohesion [Pa]
 * @param {number} o.gamma effective unit weight [N/m^3]
 * @param {number} [o.q=0] effective overburden at base level [Pa]
 * @param {number} [o.D=0] embedment [m]
 * @param {number} [o.V] vertical load Q [N] (required when H > 0)
 * @param {number} [o.H=0] horizontal load [N]
 * @param {number} [o.thetaH=Math.PI/2] angle between the long axis and H [rad]
 * @param {number} [o.baseTilt=0] ν [rad]
 * @param {number} [o.groundSlope=0] β [rad]
 * @param {'vesic'|'hansen'} [o.nGammaForm='vesic']
 * @returns {{Q:number,q:number,Kc:number,Kq:number,Kgamma:number,Nc:number,Nq:number,Ngamma:number,Aeff:number}}
 *   Q ultimate vertical capacity [N]; q [Pa]
 */
export function apiDrainedCapacity(o) {
  const { Beff, Leff, phi, c = 0, gamma, q = 0, D = 0, V, H = 0, thetaH = PI / 2, baseTilt = 0, groundSlope = 0, nGammaForm = 'vesic' } = o;
  need('Beff', Beff, 1e-9);
  need('Leff', Leff, Beff);
  need('phi', phi, 1e-6, PHI_MAX);
  need('gamma', gamma, 0);
  const strip = !Number.isFinite(Leff);
  const bl = strip ? 0 : Beff / Leff;
  const Aeff = strip ? Beff : Beff * Leff;
  const tp = Math.tan(phi);
  const Nq = bearingNq(phi);
  const Nc = bearingNc(phi);
  const Ng = bearingNgamma(phi, nGammaForm);
  let iq = 1;
  let ig = 1;
  if (H > 0) {
    need('V', V, 1e-9);
    const mB = (2 + bl) / (1 + bl);
    const mL = strip ? 1 : (2 + 1 / bl) / (1 + 1 / bl);
    const m = mL * Math.cos(thetaH) ** 2 + mB * Math.sin(thetaH) ** 2;
    const r = 1 - H / (V + (Aeff * c) / tp);
    if (r <= 0) throw new RangeError('H is too large for the inclination factors');
    iq = r ** m;
    ig = r ** (m + 1);
  }
  const cf = (x) => x - (1 - x) / (Nc * tp);
  const dq = 1 + (2 * tp * (1 - Math.sin(phi)) ** 2 * D) / Beff;
  const bq = (1 - baseTilt * tp) ** 2;
  const gq = (1 - Math.tan(groundSlope)) ** 2;
  const Kq = iq * (1 + bl * tp) * dq * bq * gq;
  const Kgamma = ig * (1 - 0.4 * bl) * bq * gq;
  const Kc = cf(iq) * (1 + (bl * Nq) / Nc) * cf(dq) * cf(bq) * cf(gq);
  const qu = c * Nc * Kc + q * Nq * Kq + 0.5 * gamma * Beff * Ng * Kgamma;
  return { Q: qu * Aeff, q: qu, Kc, Kq, Kgamma, Nc, Nq, Ngamma: Ng, Aeff };
}

/**
 * Punch-through of a footing on a strong layer over weak clay or sand (Meyerhof & Hanna 1978):
 *   qult = qb + (γ1 H² / B)(1 + 2 D/H) Ks tan φ1 · λ + (2 ca H / B) · λ − γ1 H   ≤ qt
 * where λ = 1 for a strip and (1 + B/L) for a rectangle (2 for square/circular), qb is the capacity of the
 * lower layer at depth D + H, qt that of an infinitely thick upper layer. The punching-shear coefficient Ks
 * (a chart function of φ1 and the strength ratio) and adhesion ca are arguments.
 * @param {object} o
 * @param {number} o.B footing width or diameter [m]
 * @param {number} [o.L=Infinity] length [m]
 * @param {number} o.H thickness of the strong layer below the base [m]
 * @param {number} [o.D=0] embedment of the base [m]
 * @param {number} o.gamma1 effective unit weight of the upper layer [N/m^3]
 * @param {number} o.phi1 friction angle of the upper layer [rad]
 * @param {number} o.Ks punching shear coefficient [-] (Meyerhof & Hanna chart; typ. 1–10)
 * @param {number} [o.ca=0] adhesion on the punching surface [Pa]
 * @param {number} o.qb ultimate bearing pressure of the lower layer at its top [Pa]
 * @param {number} [o.qt=Infinity] ultimate bearing pressure of the upper layer alone [Pa]
 * @returns {{qult:number,qPunch:number,governs:'punch-through'|'upper-layer'}} [Pa]
 */
export function punchThroughMeyerhofHanna({ B, L = Infinity, H, D = 0, gamma1, phi1, Ks, ca = 0, qb, qt = Infinity }) {
  need('B', B, 1e-9);
  need('L', L, B);
  need('H', H, 1e-9);
  need('D', D, 0);
  need('gamma1', gamma1, 0);
  need('phi1', phi1, 0, PHI_MAX);
  need('Ks', Ks, 0);
  need('ca', ca, 0);
  need('qb', qb, 0);
  const lam = 1 + (Number.isFinite(L) ? B / L : 0);
  const qPunch = qb + lam * (((gamma1 * H * H) / B) * (1 + (2 * D) / H) * Ks * Math.tan(phi1) + (2 * ca * H) / B) - gamma1 * H;
  return qPunch <= qt ? { qult: qPunch, qPunch, governs: 'punch-through' } : { qult: qt, qPunch, governs: 'upper-layer' };
}

/**
 * Load-spread (projected-area) punch-through check for a circular footing or spudcan on sand over clay
 * (ISO 19905-1 A.9.3.2.6.4 / SNAME T&R 5-5A form): the load spreads at 1 horizontal : ns vertical through
 * the sand and bears on the clay over an enlarged circle,
 *   Qv = (π/4)(B + 2H/ns)² (Nc sc dc su + p′0) − W_sand_frustum (optional).
 * ns is an argument (3 to 5 is the range usually quoted).
 * @param {object} o
 * @param {number} o.B diameter [m]
 * @param {number} o.H sand thickness below the footing [m]
 * @param {number} o.ns load-spread factor (1 horizontal : ns vertical) [-]
 * @param {number} o.su clay undrained strength at the layer top [Pa]
 * @param {number} [o.NcScDc=6.17] product Nc·sc·dc for the fictitious footing [-] ((2+π)·1.2 by default)
 * @param {number} [o.p0=0] effective overburden at the clay surface [Pa]
 * @param {number} [o.sandWeight=0] submerged weight of sand between footing and clay to deduct [N]
 * @returns {{Qv:number,Bequiv:number}} Qv ultimate vertical load [N]; equivalent diameter at clay level [m]
 */
export function punchThroughLoadSpread({ B, H, ns, su, NcScDc = 1.2 * (2 + PI), p0 = 0, sandWeight = 0 }) {
  need('B', B, 1e-9);
  need('H', H, 0);
  need('ns', ns, 0.5, 20);
  need('su', su, 0);
  const Be = B + (2 * H) / ns;
  return { Qv: (PI / 4) * Be * Be * (NcScDc * su + p0) - sandWeight, Bequiv: Be };
}

/**
 * Squeezing of a soft clay layer of thickness T between the footing and a strong stratum
 * (Meyerhof & Chaplin 1953 as used in SNAME T&R 5-5A / ISO 19905-1 A.9.3.2.6.5):
 *   qult = (a + b B/T + 1.2 D/B) su + p′0   when B ≥ 3.45 T (1 + 1.1 D/B),
 * never less than the general-shear value Nc sc dc su + p′0. Coefficients a = 5, b = 0.33 are the published
 * defaults and are exposed as arguments.
 * @param {object} o
 * @param {number} o.B diameter or width [m]
 * @param {number} o.T thickness of soft clay under the base [m]
 * @param {number} [o.D=0] embedment [m]
 * @param {number} o.su undrained strength [Pa]
 * @param {number} [o.p0=0] effective overburden at base level [Pa]
 * @param {number} [o.a=5]
 * @param {number} [o.b=0.33]
 * @param {number} [o.qGeneral] general-shear ultimate pressure [Pa] (default 6 (1 + 0.2 D/B) su + p′0)
 * @returns {{qult:number,squeezing:boolean,qSqueeze:number,qGeneral:number}} [Pa]
 */
export function squeezingCapacity({ B, T, D = 0, su, p0 = 0, a = 5, b = 0.33, qGeneral }) {
  need('B', B, 1e-9);
  need('T', T, 1e-9);
  need('D', D, 0);
  need('su', su, 0);
  const qg = qGeneral !== undefined ? qGeneral : 6 * (1 + (0.2 * Math.min(D / B, 2.5))) * su + p0;
  const qs = (a + (b * B) / T + (1.2 * D) / B) * su + p0;
  const applies = B >= 3.45 * T * (1 + (1.1 * D) / B);
  return { qult: applies ? Math.max(qs, qg) : qg, squeezing: applies && qs > qg, qSqueeze: qs, qGeneral: qg };
}

/**
 * Punch-through of strong clay over weak clay for a circular footing (Brown & Meyerhof 1969 as adopted in
 * SNAME T&R 5-5A / ISO 19905-1):
 *   qult = kt (H/B) su,top + Nc sc (1 + 0.2 (D + H)/B) su,bottom + p′0  ≤ general shear in the top layer.
 * kt = 3 and Nc·sc = 6 are the published defaults, exposed as arguments.
 * @param {object} o
 * @param {number} o.B diameter [m]
 * @param {number} o.H thickness of strong clay below the base [m]
 * @param {number} [o.D=0] embedment [m]
 * @param {number} o.suTop [Pa]
 * @param {number} o.suBottom [Pa]
 * @param {number} [o.p0=0] effective overburden at base level [Pa]
 * @param {number} [o.kt=3]
 * @param {number} [o.NcSc=6]
 * @returns {{qult:number,qPunch:number,qTop:number,governs:string}} [Pa]
 */
export function punchThroughClayOverClay({ B, H, D = 0, suTop, suBottom, p0 = 0, kt = 3, NcSc = 6 }) {
  need('B', B, 1e-9);
  need('H', H, 0);
  need('suTop', suTop, 0);
  need('suBottom', suBottom, 0);
  const qPunch = (kt * H * suTop) / B + NcSc * (1 + (0.2 * Math.min((D + H) / B, 2.5))) * suBottom + p0;
  const qTop = NcSc * (1 + 0.2 * Math.min(D / B, 2.5)) * suTop + p0;
  return { qult: Math.min(qPunch, qTop), qPunch, qTop, governs: qPunch < qTop ? 'punch-through' : 'top-layer' };
}

/**
 * Sliding resistance of a shallow foundation with optional skirts.
 *   Undrained: Hbase = su,base · A.  Drained: Hbase = c′ A + V tan δ.
 * Skirt (or embedded-face) earth resistance of height h and width w normal to the load, from Rankine limit
 * states integrated over h:
 *   undrained: passive (2 su,skirt + γ′h/2) h w, active max(γ′h/2 − 2 su,skirt, 0) h w (gap assumed when
 *   the active stress is tensile); drained: ½ γ′ h² w (Kp − Ka), Kp = tan²(π/4 + φ/2) = 1/Ka.
 * The net skirt term is scaled by `skirtMobilisation` because passive resistance needs larger displacement
 * than base shear.
 * @param {object} o
 * @param {'undrained'|'drained'} o.mode
 * @param {number} o.A base (effective) area [m^2]
 * @param {number} [o.suBase] undrained strength at base/skirt-tip level [Pa]
 * @param {number} [o.V] effective vertical load [N] (drained)
 * @param {number} [o.delta] interface friction angle [rad] (drained)
 * @param {number} [o.c=0] effective cohesion/adhesion [Pa] (drained)
 * @param {number} [o.skirtDepth=0] h [m]
 * @param {number} [o.skirtWidth=0] w [m]
 * @param {number} [o.suSkirt] average undrained strength over the skirt depth [Pa]
 * @param {number} [o.phi] friction angle for skirt earth pressure [rad] (drained)
 * @param {number} [o.gamma=0] effective unit weight [N/m^3]
 * @param {number} [o.skirtMobilisation=1] fraction of net skirt resistance counted [-]
 * @param {number} [o.materialFactor=1] γm applied to the total [-]
 * @returns {{H:number,Hbase:number,Hpassive:number,Hactive:number}} resistances [N]
 */
export function slidingResistance(o) {
  const { mode, A, V, delta, c = 0, skirtDepth = 0, skirtWidth = 0, gamma = 0, skirtMobilisation = 1, materialFactor = 1 } = o;
  need('A', A, 0);
  need('skirtDepth', skirtDepth, 0);
  need('skirtWidth', skirtWidth, 0);
  need('skirtMobilisation', skirtMobilisation, 0, 1);
  need('materialFactor', materialFactor, 1e-9);
  const h = skirtDepth;
  let Hbase;
  let Hp = 0;
  let Ha = 0;
  if (mode === 'undrained') {
    need('suBase', o.suBase, 0);
    Hbase = o.suBase * A;
    if (h > 0) {
      const sus = o.suSkirt !== undefined ? need('suSkirt', o.suSkirt, 0) : o.suBase;
      Hp = (2 * sus + 0.5 * gamma * h) * h * skirtWidth;
      Ha = Math.max(0.5 * gamma * h - 2 * sus, 0) * h * skirtWidth;
    }
  } else if (mode === 'drained') {
    need('V', V, 0);
    need('delta', delta, 0, PI / 2 - 1e-6);
    Hbase = c * A + V * Math.tan(delta);
    if (h > 0) {
      need('phi', o.phi, 0, PHI_MAX);
      const Kp = Math.tan(PI / 4 + o.phi / 2) ** 2;
      Hp = 0.5 * gamma * h * h * skirtWidth * Kp;
      Ha = (0.5 * gamma * h * h * skirtWidth) / Kp;
    }
  } else {
    throw new RangeError(`mode must be 'undrained' or 'drained' (got ${mode})`);
  }
  const H = (Hbase + skirtMobilisation * (Hp - Ha)) / materialFactor;
  return { H, Hbase, Hpassive: Hp, Hactive: Ha };
}

/**
 * Overturning and uplift check of a rigid footing from statics.
 *   e = M/V; FS_overturning = V (B/2) / M (rotation about the toe); no uplift while e ≤ B/6 (rectangle) or
 *   e ≤ D/8 (circle); contact over B − 2e in the Meyerhof idealisation.
 * @param {object} o
 * @param {number} o.V vertical load [N], > 0
 * @param {number} o.M overturning moment at base level [N·m]
 * @param {number} o.B width in the plane of the moment, or diameter [m]
 * @param {'rect'|'circle'} [o.shape='rect']
 * @returns {{e:number,fsOverturning:number,kernLimit:number,fullContact:boolean,contactWidth:number}}
 *   e [m]; fsOverturning [-]; kernLimit [m]; contactWidth = B − 2e [m]
 */
export function overturningCheck({ V, M, B, shape = 'rect' }) {
  need('V', V, 1e-12);
  need('M', M, 0);
  need('B', B, 1e-9);
  const e = M / V;
  const kern = shape === 'circle' ? B / 8 : B / 6;
  return { e, fsOverturning: M > 0 ? (V * B) / (2 * M) : Infinity, kernLimit: kern, fullContact: e <= kern, contactWidth: Math.max(B - 2 * e, 0) };
}

/**
 * Design verification in either working-stress or partial-factor (LRFD) format:
 *   WSD:  utilisation = load · FS / resistance;   LRFD: utilisation = γf · load / (resistance / γm).
 * Typical values (to be confirmed against the governing edition): API RP 2A-WSD bearing FS 2.0, sliding 1.5;
 * ISO 19901-4 material factors 1.25 (drained tan φ′) and 1.5 (undrained); DNV-ST-0126 γm 1.15 / 1.25.
 * @param {object} o
 * @param {number} o.load characteristic load effect [N or Pa]
 * @param {number} o.resistance characteristic resistance [same unit]
 * @param {number} [o.safetyFactor=1] global factor [-]
 * @param {number} [o.loadFactor=1] γf [-]
 * @param {number} [o.materialFactor=1] γm [-]
 * @returns {{uc:number,designLoad:number,designResistance:number,ok:boolean}} uc = demand/capacity [-]
 */
export function designCheck({ load, resistance, safetyFactor = 1, loadFactor = 1, materialFactor = 1 }) {
  need('load', load, 0);
  need('resistance', resistance, 1e-300);
  need('safetyFactor', safetyFactor, 1e-9);
  need('loadFactor', loadFactor, 1e-9);
  need('materialFactor', materialFactor, 1e-9);
  const designLoad = load * loadFactor * safetyFactor;
  const designResistance = resistance / materialFactor;
  const uc = designLoad / designResistance;
  return { uc, designLoad, designResistance, ok: uc <= 1 };
}
