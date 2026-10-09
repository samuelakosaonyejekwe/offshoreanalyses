// Vertical breakwaters and gravity walls: Goda wave pressures (with Takahashi impulsive coefficient),
// Sainflou standing-wave pressures, sliding / overturning / bearing, hydrostatic stability and flotation.
// SI units; all forces and moments are per metre run of wall (N/m, N·m/m).

const G = 9.80665;
const TWO_PI = 2 * Math.PI;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function pos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

// Linear dispersion, Newton iteration on kh.
function waveLength(T, h, g) {
  const k0h = ((TWO_PI / T) ** 2 * h) / g;
  let kh = k0h / Math.sqrt(Math.tanh(k0h));
  for (let i = 0; i < 50; i++) {
    const th = Math.tanh(kh);
    const d = (kh * th - k0h) / (th + kh * (1 - th * th));
    kh -= d;
    if (Math.abs(d) < 1e-14 * kh) break;
  }
  return (TWO_PI * h) / kh;
}

/**
 * Takahashi et al. (1994) impulsive pressure coefficient αI = αI0 αI1 for caissons on high or wide mounds:
 *   αI0 = H/d (H <= 2d), 2 (H > 2d)
 *   αI1 = cos δ2 / cosh δ1 (δ2 <= 0);  1/(cosh δ1 · cosh^(1/2) δ2) (δ2 > 0)
 *   δ1 = 20 δ11 (δ11 <= 0), 15 δ11 (δ11 > 0);  δ2 = 4.9 δ22 (δ22 <= 0), 3 δ22 (δ22 > 0)
 *   δ11 = 0.93 (BM/L − 0.12) + 0.36 ((h − d)/h − 0.6)
 *   δ22 = −0.36 (BM/L − 0.12) + 0.93 ((h − d)/h − 0.6)
 * Reference: Takahashi, Tanimoto & Shimosako (1994) Proc. Hydro-Port '94; Goda (2000) §4.2.4.
 * Validity: rubble-mound foundation with berm width BM; used as α* = max(α2, αI) in the Goda formula.
 * @param {object} o
 * @param {number} o.H design wave height (Hmax) [m]
 * @param {number} o.d depth above the mound armour [m]
 * @param {number} o.h depth in front of the breakwater [m]
 * @param {number} o.BM berm width of the mound in front of the caisson [m]
 * @param {number} o.L wavelength at depth h [m]
 * @returns {{alphaI:number, alphaI0:number, alphaI1:number}} [-]
 */
export function takahashiImpulsive({ H, d, h, BM, L }) {
  pos(H, 'H');
  pos(d, 'd');
  pos(h, 'h');
  pos(L, 'L');
  need(Number.isFinite(BM) && BM >= 0, 'BM must be >= 0');
  const alphaI0 = H <= 2 * d ? H / d : 2;
  const a = BM / L - 0.12;
  const b = (h - d) / h - 0.6;
  const d11 = 0.93 * a + 0.36 * b;
  const d22 = -0.36 * a + 0.93 * b;
  const d1 = d11 <= 0 ? 20 * d11 : 15 * d11;
  const d2 = d22 <= 0 ? 4.9 * d22 : 3 * d22;
  const alphaI1 = d2 <= 0 ? Math.cos(d2) / Math.cosh(d1) : 1 / (Math.cosh(d1) * Math.sqrt(Math.cosh(d2)));
  return { alphaI: alphaI0 * alphaI1, alphaI0, alphaI1 };
}

/**
 * Goda wave pressure on an upright (caisson) breakwater:
 *   η* = 0.75 (1 + cos β) λ1 Hmax
 *   p1 = ½ (1 + cos β)(λ1 α1 + λ2 α* cos² β) ρ g Hmax,  p2 = p1 / cosh(2πh/L),  p3 = α3 p1
 *   p4 = p1 (1 − hc/η*) for η* > hc, else 0;  pu = ½ (1 + cos β) λ3 α1 α3 ρ g Hmax
 *   α1 = 0.6 + ½ [(4πh/L)/sinh(4πh/L)]²
 *   α2 = min{ ((hb − d)/(3 hb)) (Hmax/d)², 2d/Hmax },  α3 = 1 − (h'/h)[1 − 1/cosh(2πh/L)]
 *   α* = α2, or max(α2, αI) when a berm width BM is given (Takahashi et al. 1994)
 * and the resultants
 *   P = ½ (p1 + p3) h' + ½ (p1 + p4) hc*,  hc* = min(η*, hc)
 *   Mp = (1/6)(2p1 + p3) h'² + ½ (p1 + p4) h' hc* + (1/6)(p1 + 2p4) hc*²   (about the caisson base)
 *   U = ½ pu B,  Mu = (2/3) U B   (about the heel = harbour-side base edge)
 * Reference: Goda (1974) Proc. 14th ICCE; Goda (2000) "Random Seas and Design of Maritime Structures"
 * §4.2.
 * Validity: vertical-face composite or upright breakwaters; breaking and non-breaking waves; Hmax is the
 * highest wave at the site (1.8 H1/3 seaward of the surf zone, or the surf-zone value at a distance 5 H1/3
 * seaward of the wall); β is the angle between wave direction and the wall normal (Goda recommends
 * rotating the principal direction up to 15° toward the normal). L is evaluated at depth h.
 * @param {object} o
 * @param {number} o.Hmax design wave height [m]
 * @param {number} o.T design wave period (T1/3) [s]
 * @param {number} o.h water depth in front of the breakwater [m]
 * @param {number} o.d depth above the armour layer of the rubble mound [m]
 * @param {number} o.hPrime depth of the caisson base below still-water level h' [m]
 * @param {number} o.hc crest elevation above still-water level [m]
 * @param {number} o.hb depth at 5 H1/3 seaward of the wall [m]
 * @param {number} o.B caisson width at the base [m]
 * @param {number} [o.betaDeg] wave angle from the wall normal [deg], default 0
 * @param {number} [o.BM] mound berm width; when given the impulsive coefficient αI is included [m]
 * @param {number} [o.lambda1] modification factor for slowly varying pressure, default 1
 * @param {number} [o.lambda2] modification factor for breaking pressure, default 1
 * @param {number} [o.lambda3] modification factor for uplift, default 1
 * @param {number} [o.rho] sea-water density [kg/m³], default 1030
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{etaStar:number,L:number,alpha1:number,alpha2:number,alpha3:number,alphaI:number,alphaStar:number,
 *   p1:number,p2:number,p3:number,p4:number,pu:number,P:number,Mp:number,armP:number,U:number,Mu:number,
 *   hcStar:number}}
 *   η* [m], L [m], coefficients [-], pressures [Pa], P [N/m], Mp [N·m/m], armP = Mp/P [m], U [N/m],
 *   Mu [N·m/m], hc* [m]
 */
export function godaPressure({
  Hmax, T, h, d, hPrime, hc, hb, B, betaDeg = 0, BM, lambda1 = 1, lambda2 = 1, lambda3 = 1, rho = 1030, g = G,
}) {
  pos(Hmax, 'Hmax');
  pos(T, 'T');
  pos(h, 'h');
  pos(d, 'd');
  pos(hPrime, 'hPrime');
  pos(hb, 'hb');
  pos(B, 'B');
  need(Number.isFinite(hc) && hc >= 0, 'hc must be >= 0');
  need(d <= h + 1e-12 && hPrime <= h + 1e-12, 'd and hPrime must not exceed h');
  need(Math.abs(betaDeg) < 90, 'betaDeg must lie in (-90, 90)');
  const L = waveLength(T, h, g);
  const kh = (TWO_PI * h) / L;
  const cb = Math.cos((betaDeg * Math.PI) / 180);
  const x = 2 * kh;
  const alpha1 = 0.6 + 0.5 * (x / Math.sinh(x)) ** 2;
  const alpha2 = Math.min((Math.max(hb - d, 0) / (3 * hb)) * (Hmax / d) ** 2, (2 * d) / Hmax);
  const alpha3 = 1 - (hPrime / h) * (1 - 1 / Math.cosh(kh));
  const alphaI = BM !== undefined ? takahashiImpulsive({ H: Hmax, d, h, BM, L }).alphaI : 0;
  const alphaStar = Math.max(alpha2, alphaI);
  const w = rho * g * Hmax;
  const etaStar = 0.75 * (1 + cb) * lambda1 * Hmax;
  const p1 = 0.5 * (1 + cb) * (lambda1 * alpha1 + lambda2 * alphaStar * cb * cb) * w;
  const p2 = p1 / Math.cosh(kh);
  const p3 = alpha3 * p1;
  const p4 = etaStar > hc ? p1 * (1 - hc / etaStar) : 0;
  const pu = 0.5 * (1 + cb) * lambda3 * alpha1 * alpha3 * w;
  const hcStar = Math.min(etaStar, hc);
  const P = 0.5 * (p1 + p3) * hPrime + 0.5 * (p1 + p4) * hcStar;
  const Mp = ((2 * p1 + p3) * hPrime * hPrime) / 6 + 0.5 * (p1 + p4) * hPrime * hcStar + ((p1 + 2 * p4) * hcStar * hcStar) / 6;
  const U = 0.5 * pu * B;
  const Mu = (2 / 3) * U * B;
  return { etaStar, L, alpha1, alpha2, alpha3, alphaI, alphaStar, p1, p2, p3, p4, pu, P, Mp, armP: Mp / P, U, Mu, hcStar };
}

// Contact pressure under a rigid strip footing of width B carrying a vertical load N whose line of action
// is at distance t from the "near" edge (statics, linear distribution, no tension).
function contactPressure(N, t, B) {
  const e = 0.5 * B - t;
  if (!(N > 0) || !(t > 0) || !(t < B)) {
    return { pNear: Number.POSITIVE_INFINITY, pFar: 0, contactWidth: 0, eccentricity: e, fullContact: false };
  }
  if (t < B / 3) return { pNear: (2 * N) / (3 * t), pFar: 0, contactWidth: 3 * t, eccentricity: e, fullContact: false };
  if (t > (2 * B) / 3) {
    return { pNear: 0, pFar: (2 * N) / (3 * (B - t)), contactWidth: 3 * (B - t), eccentricity: e, fullContact: false };
  }
  return {
    pNear: ((2 * N) / B) * (2 - (3 * t) / B), pFar: ((2 * N) / B) * ((3 * t) / B - 1), contactWidth: B, eccentricity: e,
    fullContact: true,
  };
}

/**
 * Sliding, overturning and heel-pressure check of an upright section under wave thrust and uplift:
 *   SF_sliding = μ (W − U)/P,  SF_overturning = (W t − Mu)/Mp
 *   te = (W t − Mu − Mp)/(W − U);  pe = 2We/(3te) for te <= B/3,  pe = (2We/B)(2 − 3te/B) for te > B/3
 * Reference: Goda (2000) §4.3 (required safety factors 1.2; heel pressure limit
 * 400–500 kPa on a rubble mound).
 * Validity: rigid body on a rubble mound; W is the weight in still water (buoyancy deducted).
 * @param {object} o
 * @param {number} o.P horizontal wave force [N/m]
 * @param {number} o.Mp overturning moment of P about the base [N·m/m]
 * @param {number} o.U uplift force [N/m]
 * @param {number} o.Mu moment of uplift about the heel [N·m/m]
 * @param {number} o.W weight of the upright section per metre in still water [N/m]
 * @param {number} o.t horizontal distance from the heel to the line of action of W [m]
 * @param {number} o.B base width [m]
 * @param {number} [o.mu] friction coefficient caisson/mound [-], default 0.6
 * @returns {{sfSliding:number,sfOverturning:number,We:number,te:number,eccentricity:number,pHeel:number,
 *   pToe:number,contactWidth:number,fullContact:boolean}}
 *   safety factors [-], net vertical force We [N/m], te [m], eccentricity from mid-base [m], heel (maximum)
 *   and seaward-edge (minimum) bearing pressures [Pa], contact width [m]
 */
export function caissonStability({ P, Mp, U, Mu, W, t, B, mu = 0.6 }) {
  pos(P, 'P');
  pos(Mp, 'Mp');
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  need(Number.isFinite(Mu) && Mu >= 0, 'Mu must be >= 0');
  pos(W, 'W');
  pos(t, 't');
  pos(B, 'B');
  pos(mu, 'mu');
  const We = W - U;
  const te = (W * t - Mu - Mp) / We;
  const c = contactPressure(We, te, B);
  return {
    sfSliding: (mu * We) / P, sfOverturning: (W * t - Mu) / Mp, We, te, eccentricity: c.eccentricity,
    pHeel: c.pNear, pToe: c.pFar, contactWidth: c.contactWidth, fullContact: c.fullContact,
  };
}

/**
 * Sainflou (1928) standing-wave (clapotis) pressure on a vertical wall, simplified linear distribution:
 *   h0 = (π H²/L) coth(2πh/L)   (rise of the mean level of the clapotis)
 *   crest:  p = 0 at h0 + H above SWL, p1 = (p2 + ρ g h)(H + h0)/(h + H + h0) at SWL, p2 = ρ g H/cosh(2πh/L)
 *           at the bed (in excess of still-water hydrostatic pressure)
 *   trough: water surface at H − h0 below SWL; hydrostatic deficit ρ g (H − h0) at that level and p2 at
 *           the bed (net force directed seaward when the lee side is at still water)
 * Reference: Sainflou (1928) Annales des Ponts et Chaussées 98(4); SPM (1984) §7.III.2 (Miche–Rundgren
 * presentation); Goda (2000) §4.1.
 * Validity: fully reflected non-breaking waves, wall on a flat bed; overestimates for steep waves. If a
 * crest level hc is given, the pressure triangle is truncated there.
 * @param {object} o
 * @param {number} o.H incident wave height [m]
 * @param {number} o.T period [s]
 * @param {number} o.h depth at the wall [m]
 * @param {number} [o.hc] wall crest above still-water level [m] (default: not overtopped)
 * @param {number} [o.rho] density [kg/m³], default 1030
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{L:number,h0:number,p1:number,p2:number,crest:{force:number,moment:number,arm:number,pTop:number},
 *   trough:{force:number,moment:number,arm:number,pSurface:number}}}
 *   L [m], h0 [m], pressures [Pa]; crest: landward wave force in excess of still-water hydrostatic [N/m],
 *   moment about the bed [N·m/m], arm [m]; trough: seaward net force [N/m], moment about the bed, arm
 */
export function sainflou({ H, T, h, hc, rho = 1030, g = G }) {
  pos(H, 'H');
  pos(T, 'T');
  pos(h, 'h');
  const L = waveLength(T, h, g);
  const kh = (TWO_PI * h) / L;
  const w = rho * g;
  const h0 = ((Math.PI * H * H) / L) / Math.tanh(kh);
  const p2 = (w * H) / Math.cosh(kh);
  const top = H + h0;
  const p1 = ((p2 + w * h) * top) / (h + top);
  // Crest: trapezoid below SWL (p1 at z = h, p2 at z = 0, z measured up from the bed) plus the triangle above.
  const zc = hc !== undefined ? Math.min(Math.max(hc, 0), top) : top;
  const pTop = p1 * (1 - zc / top);
  const fBelow = 0.5 * (p1 + p2) * h;
  const mBelow = (h * h * (2 * p1 + p2)) / 6;
  const fAbove = 0.5 * (p1 + pTop) * zc;
  const mAbove = fAbove * h + (zc * zc * (p1 + 2 * pTop)) / 6;
  const cf = fBelow + fAbove;
  const cm = mBelow + mAbove;
  // Trough: deficit grows hydrostatically from SWL to the trough surface, then varies linearly to p2 at the bed.
  const dt = Math.min(Math.max(H - h0, 0), h);
  const ps = w * dt;
  const f1 = 0.5 * ps * dt;
  const m1 = f1 * (h - (2 * dt) / 3);
  const lower = h - dt;
  const f2 = 0.5 * (ps + p2) * lower;
  const m2 = (lower * lower * (2 * ps + p2)) / 6;
  const tf = f1 + f2;
  const tm = m1 + m2;
  return {
    L, h0, p1, p2,
    crest: { force: cf, moment: cm, arm: cm / cf, pTop },
    trough: { force: tf, moment: tm, arm: tm / tf, pSurface: ps },
  };
}

/**
 * Hydrostatic stability of a gravity wall (seawall, quay, intake or outfall headwall) with different water
 * levels on its two faces. The base is horizontal at z = 0; the high-water face is at x = 0, the
 * low-water face at x = B. Loads per metre run:
 *   Fh = ½ ρ g hHigh² at hHigh/3,  Fl = ½ ρ g hLow² at hLow/3 (resisting)
 *   uplift on the base varying linearly from ρ g hHigh to ρ g hLow (includes buoyancy); it can be scaled by
 *   upliftFactor for drained bases
 *   sliding FS = μ (W − U)/(Fh − Fl + Hextra)
 *   overturning about the low-side toe: FS = [W (B − xW) + Fl hLow/3] / [Fh hHigh/3 + U aU + Mextra]
 *   base pressure from the resultant position (no tension), effective stresses
 * Reference: statics — hydrostatic pressure, Archimedes' principle; arrangement per USACE EM 1110-2-2100
 * (2005) ch. 3–4. Validity: rigid wall, impermeable base contact, no cohesion or passive resistance.
 * @param {object} o
 * @param {number} o.W total weight of the wall in air per metre (including fill on it) [N/m]
 * @param {number} o.xW distance of the weight line of action from the high-water face [m]
 * @param {number} o.B base width [m]
 * @param {number} o.hHigh water depth above base on the high side [m]
 * @param {number} o.hLow water depth above base on the low side [m]
 * @param {number} [o.mu] base friction coefficient [-], default 0.5
 * @param {number} [o.upliftFactor] multiplier on the uplift diagram [-], default 1
 * @param {number} [o.Hextra] additional horizontal load toward the low side (waves, ice, earth) [N/m]
 * @param {number} [o.Mextra] moment of the additional load about the base [N·m/m]
 * @param {number} [o.rho] water density [kg/m³], default 1025
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Fh:number,Fl:number,U:number,armU:number,N:number,H:number,sfSliding:number,sfOverturning:number,
 *   sfFlotation:number,resultantFromToe:number,eccentricity:number,pMax:number,pMin:number,contactWidth:number,
 *   middleThird:boolean}}
 *   forces [N/m], armU = lever arm of uplift about the low-side toe [m], N = W − U [N/m], net horizontal load
 *   H [N/m], safety factors [-] (Infinity when there is no destabilising load), resultant position from the
 *   toe [m], eccentricity [m], effective base pressures [Pa], contact width [m]
 */
export function gravityWallStability({
  W, xW, B, hHigh, hLow, mu = 0.5, upliftFactor = 1, Hextra = 0, Mextra = 0, rho = 1025, g = G,
}) {
  pos(W, 'W');
  pos(B, 'B');
  need(Number.isFinite(xW) && xW >= 0 && xW <= B, 'xW must lie in [0, B]');
  need(Number.isFinite(hHigh) && hHigh >= 0, 'hHigh must be >= 0');
  need(Number.isFinite(hLow) && hLow >= 0 && hLow <= hHigh, 'hLow must lie in [0, hHigh]');
  need(upliftFactor >= 0, 'upliftFactor must be >= 0');
  const w = rho * g;
  const Fh = 0.5 * w * hHigh * hHigh;
  const Fl = 0.5 * w * hLow * hLow;
  const U = upliftFactor * w * 0.5 * (hHigh + hLow) * B;
  // Centroid of the trapezoidal uplift diagram measured from the low-side toe.
  const armU = hHigh + hLow > 0 ? (B * (2 * hHigh + hLow)) / (3 * (hHigh + hLow)) : 0.5 * B;
  const N = W - U;
  const H = Fh - Fl + Hextra;
  const mRes = W * (B - xW) + (Fl * hLow) / 3;
  const mOt = (Fh * hHigh) / 3 + U * armU + Mextra;
  const xr = (mRes - mOt) / N;
  const c = contactPressure(N, xr, B);
  return {
    Fh, Fl, U, armU, N, H,
    sfSliding: H > 0 ? (mu * N) / H : Number.POSITIVE_INFINITY,
    sfOverturning: mOt > 0 ? mRes / mOt : Number.POSITIVE_INFINITY,
    sfFlotation: U > 0 ? W / U : Number.POSITIVE_INFINITY,
    resultantFromToe: xr, eccentricity: 0.5 * B - xr, pMax: Math.max(c.pNear, c.pFar), pMin: Math.min(c.pNear, c.pFar), contactWidth: c.contactWidth,
    middleThird: xr >= B / 3 && xr <= (2 * B) / 3,
  };
}

/**
 * Flotation (uplift) safety of a buried or submerged structure: FS = (W + R)/(ρ g V + Uextra).
 * Reference: Archimedes' principle; format per USACE EM 1110-2-2100 (2005) §3-10 (required FS typically
 * 1.1–1.5 depending on load condition).
 * @param {object} o
 * @param {number} o.weight dead weight of structure and permanent ballast [N]
 * @param {number} o.displacedVolume volume below the water (or groundwater) surface [m³]
 * @param {number} [o.resisting] additional resisting force (anchors, side friction, surcharge) [N]
 * @param {number} [o.upliftExtra] additional uplift not captured by buoyancy (e.g. excess pore pressure) [N]
 * @param {number} [o.rho] fluid density [kg/m³], default 1025
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{buoyancy:number, sf:number, netDownward:number, floats:boolean, ballastForSf:(sf:number)=>number}}
 *   buoyancy [N], FS [-], net downward force [N], ballastForSf(sf) = extra weight needed to reach sf [N]
 */
export function flotationCheck({ weight, displacedVolume, resisting = 0, upliftExtra = 0, rho = 1025, g = G }) {
  need(Number.isFinite(weight) && weight >= 0, 'weight must be >= 0');
  pos(displacedVolume, 'displacedVolume');
  need(resisting >= 0, 'resisting must be >= 0');
  need(upliftExtra >= 0, 'upliftExtra must be >= 0');
  const buoyancy = rho * g * displacedVolume;
  const up = buoyancy + upliftExtra;
  const down = weight + resisting;
  return { buoyancy, sf: down / up, netDownward: down - up, floats: down < up, ballastForSf: (sf) => Math.max(0, sf * up - down) };
}
