/**
 * Suction caissons, suction anchors and bucket foundations: installation and in-place capacity.
 *
 * SI units: m, N, Pa, N/m^3, m/s, rad. Depths are measured from the mudline; suction s is the pressure
 * differential across the lid (positive = lower pressure inside).
 *
 * References
 *   Houlsby G.T. & Byrne B.W. (2005a) Design procedures for installation of suction caissons in clay and
 *     other materials, Proc. ICE Geotech. Eng. 158(GE2), 75–82.
 *   Houlsby G.T. & Byrne B.W. (2005b) Design procedures for installation of suction caissons in sand,
 *     Proc. ICE Geotech. Eng. 158(GE3), 135–144.
 *   Houlsby G.T., Kelly R.B. & Byrne B.W. (2005) The tensile capacity of suction caissons in sand under
 *     rapid loading, ISFOG.
 *   Murff J.D. & Hamilton J.M. (1993) P-ultimate for undrained analysis of laterally loaded piles, JGE 119(1).
 *   Supachawarote C., Randolph M.F. & Gourvenec S. (2004) Inclined pull-out capacity of suction caissons,
 *     ISOPE; Randolph M.F. & Gourvenec S. (2011) Offshore Geotechnical Engineering, ch. 9.
 *   Broms B.B. (1964) Lateral resistance of piles in cohesionless soils, JSMFD 90(SM3).
 *   DNV-RP-E303 (2021) Geotechnical design and installation of suction anchors in clay.
 */
import { brent } from '../math/roots.js';

const PI = Math.PI;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

function simpson(f, a, b, n = 200) {
  if (b <= a) return 0;
  const h = (b - a) / n;
  let s = 0;
  for (let i = 0; i <= n; i++) s += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * f(a + i * h);
  return (s * h) / 3;
}

/**
 * Builds an undrained-strength profile function from either a function, a linear law {su0, k} or a list
 * of layers [{zTop, su0, k}] (strength su0 at the layer top increasing at k below it).
 * @param {((z:number)=>number)|{su0:number,k?:number}|{zTop:number,su0:number,k?:number}[]} su
 * @returns {(z:number)=>number} su(z) [Pa], z [m]
 */
export function strengthProfile(su) {
  if (typeof su === 'function') return su;
  if (Array.isArray(su)) {
    if (su.length === 0) throw new RangeError('su layer list must not be empty');
    const ls = [...su].sort((a, b) => a.zTop - b.zTop);
    return (z) => {
      let cur = ls[0];
      for (const l of ls) if (z >= l.zTop) cur = l;
      return cur.su0 + (cur.k || 0) * Math.max(z - cur.zTop, 0);
    };
  }
  if (su && typeof su.su0 === 'number') return (z) => su.su0 + (su.k || 0) * z;
  throw new RangeError('su must be a function, {su0, k} or an array of layers');
}

function suBreaks(su, h) {
  const b = [0];
  if (Array.isArray(su)) for (const l of su) if (l.zTop > 0 && l.zTop < h) b.push(l.zTop);
  b.sort((x, y) => x - y);
  b.push(h);
  return b;
}

function suIntegral(su, f, h) {
  const b = suBreaks(su, h);
  let s = 0;
  for (let i = 0; i < b.length - 1; i++) s += simpson(f, b[i], b[i + 1] - (i < b.length - 2 ? 1e-12 : 0), 100);
  return s;
}

/**
 * Installation of a suction caisson in clay (Houlsby & Byrne 2005a), including layered profiles.
 *   resistance  R(h) = αo ∫su dz · πDo + αi ∫su dz · πDi + (γ′h + Nc su,tip) · π D t;
 *   with suction  V′ + s πDo²/4 = R(h)   (H&B eq. with the tip stress reduced by s on the annulus);
 *   plug-heave (reverse bearing) limit on suction
 *       s_allow = N*c su,tip + 4 αi ∫su dz / Di − 4 Do αo ∫su dz / ((m² − 1) Do²),
 *   where the last term (enhanced vertical stress outside, carried over a zone of diameter m·Do) is
 *   omitted for m = Infinity.
 * @param {object} o
 * @param {number} o.Do outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.h tip penetration [m]
 * @param {*} o.su strength profile (see strengthProfile) [Pa]
 * @param {number} [o.alphaO=0.5] outside adhesion factor [-]
 * @param {number} [o.alphaI=0.5] inside adhesion factor [-]
 * @param {number} [o.Nc=9] tip bearing factor (deep strip) [-]
 * @param {number} [o.NcPlug=9] reverse end-bearing factor for the plug N*c [-]
 * @param {number} [o.gamma=0] effective unit weight [N/m^3]
 * @param {number} [o.Vprime=0] effective vertical load (submerged weight + ballast) [N]
 * @param {number} [o.m=Infinity] outside stress-enhancement zone ratio Dm/Do [-]
 * @returns {{R:number,outside:number,inside:number,tip:number,s:number,sAllow:number,fsPlug:number}}
 *   forces [N]; required suction s ≥ 0 [Pa]; allowable suction [Pa]; fsPlug = sAllow/s [-]
 */
export function suctionInstallClay({ Do, t, h, su, alphaO = 0.5, alphaI = 0.5, Nc = 9, NcPlug = 9, gamma = 0, Vprime = 0, m = Infinity }) {
  need('Do', Do, 1e-9);
  need('t', t, 1e-9, Do / 2);
  need('h', h, 0);
  need('alphaO', alphaO, 0, 1.5);
  need('alphaI', alphaI, 0, 1.5);
  need('m', m, 1 + 1e-9, Infinity);
  const f = strengthProfile(su);
  const Di = Do - 2 * t;
  const D = Do - t;
  const I = suIntegral(su, f, h);
  const suTip = f(h);
  const outside = alphaO * I * PI * Do;
  const inside = alphaI * I * PI * Di;
  const tip = (gamma * h + Nc * suTip) * PI * D * t;
  const R = outside + inside + tip;
  const s = Math.max((R - Vprime) / ((PI * Do * Do) / 4), 0);
  const sAllow = NcPlug * suTip + (4 * alphaI * I) / Di - (Number.isFinite(m) ? (4 * Do * alphaO * I) / ((m * m - 1) * Do * Do) : 0);
  return { R, outside, inside, tip, s, sAllow, fsPlug: s > 0 ? sAllow / s : Infinity };
}

/**
 * Self-weight penetration in clay: depth at which the resistance of suctionInstallClay equals V′.
 * @param {object} o arguments of suctionInstallClay without h, plus
 * @param {number} o.L skirt length (upper bound of the search) [m]
 * @returns {number} self-weight penetration [m] (L if the caisson penetrates fully under V′)
 */
export function selfWeightPenetrationClay(o) {
  need('L', o.L, 1e-9);
  need('Vprime', o.Vprime, 0);
  const g = (h) => suctionInstallClay({ ...o, h }).R - o.Vprime;
  if (g(0) >= 0) return 0;
  if (g(o.L) <= 0) return o.L;
  return brent(g, 0, o.L, { tol: 1e-10 });
}

/**
 * Pore-pressure factor a for a caisson in sand under suction (Houlsby & Byrne 2005b): the excess pore
 * pressure at the skirt tip is −a·s.  For uniform permeability
 *     a1 = c0 − c1 [1 − exp(−(h/D)/c2)],  c0 = 0.45, c1 = 0.36, c2 = 0.48,
 * and with the plug loosened by seepage (kf = ki/ko > 1)
 *     a = a1 kf / [(1 − a1) + a1 kf].
 * @param {number} hOverD penetration ratio h/D [-]
 * @param {number} [kf=1] permeability ratio inside/outside [-]
 * @returns {number} a [-]
 */
export function seepageFactorA(hOverD, kf = 1) {
  need('hOverD', hOverD, 0);
  need('kf', kf, 1e-9);
  const a1 = 0.45 - 0.36 * (1 - Math.exp(-hOverD / 0.48));
  return (a1 * kf) / (1 - a1 + a1 * kf);
}

/**
 * Installation of a suction caisson in sand (Houlsby & Byrne 2005b).
 * Vertical effective stress next to the wall is enhanced by wall friction:
 *     σ′v = γ′ Z [exp(z/Z) − 1],   Zi = Di / (4 (K tan δ)i),   Zo = Do (m² − 1) / (4 (K tan δ)o),
 * giving the self-weight resistance
 *     R0 = γ′ Zo² [exp(h/Zo) − 1 − h/Zo] (K tan δ)o πDo + γ′ Zi² [exp(h/Zi) − 1 − h/Zi] (K tan δ)i πDi
 *          + (γ′ Zi [exp(h/Zi) − 1] Nq + ½ γ′ t Nγ) π D t.
 * Under suction the seepage gradients change the effective unit weight to γ′ + a s/h outside and
 * γ′ − (1 − a) s/h inside (which also governs the tip term), so that
 *     V′ + s πDi²/4 = (γ′ + a s/h) Fo + (γ′ − (1 − a) s/h)(Fi + Fq) + ½ γ′ t Nγ π D t,
 * which is linear in s. The critical suction for piping (zero effective stress inside) is
 *     s_crit = γ′ h / (1 − a).
 * The self-weight tip term uses ½ γ′ t Nγ (strip of width t).
 * @param {object} o
 * @param {number} o.Do outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.h tip penetration [m], > 0
 * @param {number} o.gamma effective unit weight γ′ [N/m^3]
 * @param {number} o.KtanDeltaO (K tan δ) outside [-]
 * @param {number} o.KtanDeltaI (K tan δ) inside [-]
 * @param {number} o.Nq tip bearing factor [-]
 * @param {number} o.Ngamma tip bearing factor [-]
 * @param {number} [o.m=1.5] outside load-spread ratio Dm/Do [-], > 1
 * @param {number} [o.Vprime=0] effective vertical load [N]
 * @param {number} [o.kf=1] permeability ratio ki/ko [-]
 * @param {number} [o.a] pore-pressure factor [-] (default seepageFactorA(h/Do, kf))
 * @returns {{R0:number,outside:number,inside:number,tip:number,s:number,sCrit:number,a:number,
 *   piping:boolean,Zo:number,Zi:number}} R0 self-weight resistance and its parts [N]; s required suction
 *   [Pa]; sCrit [Pa]
 */
export function suctionInstallSand({ Do, t, h, gamma, KtanDeltaO, KtanDeltaI, Nq, Ngamma, m = 1.5, Vprime = 0, kf = 1, a }) {
  need('Do', Do, 1e-9);
  need('t', t, 1e-9, Do / 2);
  need('h', h, 1e-9);
  need('gamma', gamma, 1e-9);
  need('KtanDeltaO', KtanDeltaO, 1e-9);
  need('KtanDeltaI', KtanDeltaI, 1e-9);
  need('Nq', Nq, 0);
  need('Ngamma', Ngamma, 0);
  need('m', m, 1 + 1e-9);
  const Di = Do - 2 * t;
  const D = Do - t;
  const Zo = (Do * (m * m - 1)) / (4 * KtanDeltaO);
  const Zi = Di / (4 * KtanDeltaI);
  // per unit effective weight
  const Fo = Zo * Zo * (Math.expm1(h / Zo) - h / Zo) * KtanDeltaO * PI * Do;
  const Fi = Zi * Zi * (Math.expm1(h / Zi) - h / Zi) * KtanDeltaI * PI * Di;
  const Fq = Zi * Math.expm1(h / Zi) * Nq * PI * D * t;
  const Fg = 0.5 * gamma * t * Ngamma * PI * D * t;
  const R0 = gamma * (Fo + Fi + Fq) + Fg;
  const aa = a !== undefined ? need('a', a, 0, 1 - 1e-9) : seepageFactorA(h / Do, kf);
  const den = (PI * Di * Di) / 4 - (aa * Fo) / h + ((1 - aa) * (Fi + Fq)) / h;
  const s = Math.max((R0 - Vprime) / den, 0);
  const sCrit = (gamma * h) / (1 - aa);
  return { R0, outside: gamma * Fo, inside: gamma * Fi, tip: gamma * Fq + Fg, s, sCrit, a: aa, piping: s > sCrit, Zo, Zi };
}

/**
 * Self-weight penetration in sand: depth at which R0 of suctionInstallSand equals V′.
 * @param {object} o arguments of suctionInstallSand without h, plus L skirt length [m]
 * @returns {number} [m]
 */
export function selfWeightPenetrationSand(o) {
  need('L', o.L, 1e-9);
  need('Vprime', o.Vprime, 0);
  const g = (h) => suctionInstallSand({ ...o, h }).R0 - o.Vprime;
  if (g(o.L) <= 0) return o.L;
  if (g(1e-6) >= 0) return 0;
  return brent(g, 1e-6, o.L, { tol: 1e-10 });
}

/**
 * Suction–penetration curve.
 * @param {object} o
 * @param {'clay'|'sand'} o.soil
 * @param {number} o.L skirt length [m]
 * @param {number} [o.n=50] number of depth points
 * @param {object} o.params arguments of suctionInstallClay / suctionInstallSand (without h)
 * @returns {{h:number[],s:number[],sLimit:number[],selfWeight:number}} s required and limiting suction
 *   (plug heave in clay, piping in sand) [Pa] at depths h [m]; self-weight penetration [m]
 */
export function suctionPenetrationCurve({ soil, L, n = 50, params }) {
  need('L', L, 1e-9);
  const h = [];
  const s = [];
  const sLimit = [];
  for (let i = 1; i <= n; i++) {
    const z = (L * i) / n;
    const r = soil === 'clay' ? suctionInstallClay({ ...params, h: z }) : suctionInstallSand({ ...params, h: z });
    h.push(z);
    s.push(r.s);
    sLimit.push(soil === 'clay' ? r.sAllow : r.sCrit);
  }
  const selfWeight = soil === 'clay' ? selfWeightPenetrationClay({ ...params, L }) : selfWeightPenetrationSand({ ...params, L });
  return { h, s, sLimit, selfWeight };
}

/**
 * Pump flow needed to advance the caisson: displaced water plus (in sand) Darcy seepage through the plug,
 *   Q = (π Di²/4) [ v + k (1 − a) s / (γw h) ],
 * the second term being the upward seepage velocity for the average hydraulic gradient inside the caisson.
 * @param {object} o
 * @param {number} o.Di inner diameter [m]
 * @param {number} o.v penetration rate [m/s]
 * @param {number} [o.k=0] hydraulic conductivity of the plug [m/s] (0 for clay)
 * @param {number} [o.s=0] suction [Pa]
 * @param {number} [o.h=1] penetration [m]
 * @param {number} [o.a=0] pore-pressure factor [-]
 * @param {number} [o.gammaW=10055] unit weight of sea water [N/m^3]
 * @returns {{Q:number,Qdisplacement:number,Qseepage:number,gradient:number}} flows [m^3/s]; gradient [-]
 */
export function pumpFlow({ Di, v, k = 0, s = 0, h = 1, a = 0, gammaW = 10055 }) {
  need('Di', Di, 1e-9);
  need('v', v, 0);
  need('k', k, 0);
  need('h', h, 1e-9);
  const A = (PI * Di * Di) / 4;
  const gradient = ((1 - a) * s) / (gammaW * h);
  return { Q: A * (v + k * gradient), Qdisplacement: A * v, Qseepage: A * k * gradient, gradient };
}

/**
 * Elastic buckling pressure of a long unstiffened cylinder under external (suction) pressure
 * (Bresse–Bryan):  p_cr = 2E/(1 − ν²) (t/D)³.  A lower-bound hook only: short, soil-supported or stiffened
 * shells buckle at higher pressure and need a code check (e.g. DNV-RP-C202).
 * @param {object} o E [Pa], nu [-], t [m], D mean diameter [m]
 * @returns {number} p_cr [Pa]
 */
export function ringBucklingPressure({ E, nu = 0.3, t, D }) {
  need('E', E, 1e-9);
  need('t', t, 1e-12);
  need('D', D, 1e-9);
  return ((2 * E) / (1 - nu * nu)) * (t / D) ** 3;
}

/**
 * Vertical pull-out capacity of a caisson in clay.
 *   sealed (undrained, passive suction):  V = W′ + αo ∫su dz πDo + Nc su,tip πDo²/4   (reverse end bearing,
 *       Nc ≈ 9 for L/D ≳ 1; DNV-RP-E303 notes lower values for shallow caissons);
 *   vented (top open):                    V = W′ + (αo πDo + αi πDi) ∫su dz;
 *   plug pull-out (sealed, tension failure at the base): V = W′ + αo ∫su dz πDo + W′plug.
 * The sealed capacity is the lesser of the reverse-bearing and plug mechanisms.
 * @param {object} o
 * @param {number} o.Do [m]  @param {number} o.t [m]  @param {number} o.L embedded length [m]
 * @param {*} o.su strength profile [Pa]
 * @param {number} [o.alphaO=0.65] set-up outside adhesion factor [-]
 * @param {number} [o.alphaI=0.65] inside adhesion factor [-]
 * @param {number} [o.Nc=9] reverse end-bearing factor [-]
 * @param {number} [o.W=0] submerged weight of the caisson [N]
 * @param {number} [o.gamma=0] effective unit weight of the plug [N/m^3]
 * @returns {{sealed:number,vented:number,reverseBearing:number,plugPullout:number,frictionOut:number,
 *   frictionIn:number,endBearing:number}} [N]
 */
export function caissonVerticalCapacityClay({ Do, t, L, su, alphaO = 0.65, alphaI = 0.65, Nc = 9, W = 0, gamma = 0 }) {
  need('Do', Do, 1e-9);
  need('t', t, 1e-9, Do / 2);
  need('L', L, 0);
  const f = strengthProfile(su);
  const Di = Do - 2 * t;
  const I = suIntegral(su, f, L);
  const frictionOut = alphaO * I * PI * Do;
  const frictionIn = alphaI * I * PI * Di;
  const endBearing = (Nc * f(L) * PI * Do * Do) / 4;
  const reverseBearing = W + frictionOut + endBearing;
  const plugPullout = W + frictionOut + (gamma * L * PI * Di * Di) / 4;
  return { sealed: Math.min(reverseBearing, plugPullout), vented: W + frictionOut + frictionIn, reverseBearing, plugPullout, frictionOut, frictionIn, endBearing };
}

/**
 * Lateral bearing factor Np(z) of Murff & Hamilton (1993), simplified empirical form:
 *     Np = N1 − N2 exp(−ξ z/D),  N1 = 9, N2 = 7,
 *     ξ = 0.25 + 0.05 λ (λ < 6),  0.55 (λ ≥ 6),  λ = su0 / (k D).
 * (Np = 2 at the mudline rising to N1 at depth; N1 up to 11.94 for a fully rough pile.)
 * @param {object} o
 * @param {number} o.z depth [m]  @param {number} o.D diameter [m]
 * @param {number} o.su0 mudline strength [Pa]  @param {number} [o.k=0] strength gradient [Pa/m]
 * @param {number} [o.N1=9]  @param {number} [o.N2=7]
 * @returns {number} Np [-]
 */
export function murffHamiltonNp({ z, D, su0, k = 0, N1 = 9, N2 = 7 }) {
  need('z', z, 0);
  need('D', D, 1e-9);
  need('su0', su0, 0);
  const lambda = k > 0 ? su0 / (k * D) : Infinity;
  const xi = lambda < 6 ? 0.25 + 0.05 * lambda : 0.55;
  return N1 - N2 * Math.exp((-xi * z) / D);
}

/**
 * Lateral capacity of a rigid pile or caisson from a given ultimate soil reaction p(z) [N/m]: the pile
 * rotates about a depth z0 with full reaction mobilised in opposite directions above and below. z0 follows
 * from moment equilibrium about the load point and H from horizontal equilibrium (base shear and base
 * moment neglected – conservative for stubby caissons).
 * @param {object} o
 * @param {(z:number)=>number} o.p ultimate lateral reaction per unit length [N/m]
 * @param {number} o.L embedded length [m]
 * @param {number} o.zLoad depth of the load application point [m]; negative = above the mudline
 * @param {number} [o.n=400] integration intervals
 * @returns {{H:number,z0:number,Htranslation:number,zCentroid:number}} H capacity [N]; rotation depth [m];
 *   Htranslation = ∫p dz, reached when zLoad = zCentroid [N]; centroid of the reaction [m]
 */
export function rigidPileLateralCapacity({ p, L, zLoad, n = 400 }) {
  need('L', L, 1e-9);
  need('zLoad', zLoad, -Infinity, L);
  if (typeof p !== 'function') throw new RangeError('p must be a function of depth');
  const dz = L / n;
  const P = new Float64Array(n + 1); // ∫0^z p
  const Q = new Float64Array(n + 1); // ∫0^z p (z − zLoad)
  for (let i = 1; i <= n; i++) {
    const za = (i - 1) * dz;
    const zb = i * dz;
    const zm = 0.5 * (za + zb);
    const pa = p(za);
    const pm = p(zm);
    const pb = p(zb);
    P[i] = P[i - 1] + (dz / 6) * (pa + 4 * pm + pb);
    Q[i] = Q[i - 1] + (dz / 6) * (pa * (za - zLoad) + 4 * pm * (zm - zLoad) + pb * (zb - zLoad));
  }
  const Ptot = P[n];
  const Qtot = Q[n];
  const zCentroid = Ptot > 0 ? Qtot / Ptot + zLoad : L / 2;
  const at = (arr, z) => {
    const x = Math.min(Math.max(z / dz, 0), n);
    const i = Math.min(Math.floor(x), n - 1);
    return arr[i] + (arr[i + 1] - arr[i]) * (x - i);
  };
  const f = (z) => 2 * at(Q, z) - Qtot;
  let lo;
  let hi;
  if (Qtot >= 0) { lo = Math.max(zLoad, 0); hi = L; } else { lo = 0; hi = Math.max(zLoad, 0); }
  let z0;
  if (Math.abs(Qtot) < 1e-12 * Math.abs(Ptot) * L) z0 = Qtot >= 0 ? L : 0;
  else {
    const flo = f(lo);
    for (let i = 0; i < 80; i++) {
      const mid = 0.5 * (lo + hi);
      if ((f(mid) > 0) === (flo > 0)) lo = mid; else hi = mid;
    }
    z0 = 0.5 * (lo + hi);
  }
  return { H: Math.abs(2 * at(P, z0) - Ptot), z0, Htranslation: Ptot, zCentroid };
}

/**
 * Horizontal capacity of a suction anchor / caisson in clay, Hmax = ∫ Np su D dz for pure translation
 * (load at the optimal padeye), and the reduced capacity for another load attachment depth by rigid
 * rotation (rigidPileLateralCapacity). Np from murffHamiltonNp unless a constant NpConst is given
 * (e.g. 10–12 for a quick estimate Np su L D). The optimal padeye depth is the centroid of the lateral
 * resistance (about 2L/3 for strength proportional to depth, L/2 for uniform strength and constant Np).
 * @param {object} o
 * @param {number} o.D outer diameter [m]  @param {number} o.L embedded length [m]
 * @param {number} o.su0 mudline strength [Pa]  @param {number} [o.k=0] gradient [Pa/m]
 * @param {number} [o.NpConst] constant lateral bearing factor [-]
 * @param {number} [o.N1=9]  @param {number} [o.N2=7] Murff–Hamilton parameters
 * @param {number} [o.zPadeye] load attachment depth [m] (default: optimal)
 * @returns {{Hmax:number,H:number,zOptimal:number,z0:number,NpAverage:number}} [N, N, m, m, -]
 */
export function caissonHorizontalCapacityClay({ D, L, su0, k = 0, NpConst, N1 = 9, N2 = 7, zPadeye }) {
  need('D', D, 1e-9);
  need('L', L, 1e-9);
  need('su0', su0, 0);
  need('k', k, 0);
  const p = (z) => (NpConst !== undefined ? NpConst : murffHamiltonNp({ z, D, su0, k, N1, N2 })) * (su0 + k * z) * D;
  const base = rigidPileLateralCapacity({ p, L, zLoad: 0 });
  const zOptimal = base.zCentroid;
  const r = zPadeye === undefined ? { H: base.Htranslation, z0: Infinity } : rigidPileLateralCapacity({ p, L, zLoad: zPadeye });
  const suAvg = su0 + 0.5 * k * L;
  return { Hmax: base.Htranslation, H: r.H, zOptimal, z0: r.z0, NpAverage: suAvg > 0 ? base.Htranslation / (suAvg * L * D) : 0 };
}

/**
 * Inclined-load interaction for suction caissons in clay loaded at the optimal padeye
 * (Supachawarote, Randolph & Gourvenec 2004):
 *     (H/Hmax)^a + (V/Vmax)^b = 1,   a = L/D + 0.5,   b = L/(3D) + 4.5.
 * Returns the envelope value for a load pair and the ultimate load for a given inclination.
 * @param {object} o
 * @param {number} o.Hmax [N]  @param {number} o.Vmax [N]
 * @param {number} o.L [m]  @param {number} o.D [m]
 * @param {number} [o.H] horizontal load component [N]  @param {number} [o.V] vertical component [N]
 * @param {number} [o.angle] load inclination above horizontal at the padeye [rad]
 * @param {number} [o.a]  @param {number} [o.b] exponents (override)
 * @returns {{f:number|undefined,Tult:number|undefined,a:number,b:number}} f ≤ 1 inside; Tult [N]
 */
export function inclinedInteraction({ Hmax, Vmax, L, D, H, V, angle, a, b }) {
  need('Hmax', Hmax, 1e-12);
  need('Vmax', Vmax, 1e-12);
  need('L', L, 1e-9);
  need('D', D, 1e-9);
  const ea = a !== undefined ? a : L / D + 0.5;
  const eb = b !== undefined ? b : L / (3 * D) + 4.5;
  const val = (h, v) => (Math.abs(h) / Hmax) ** ea + (Math.abs(v) / Vmax) ** eb;
  let Tult;
  if (angle !== undefined) {
    need('angle', angle, 0, PI / 2);
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    let lo = 0;
    let hi = Math.hypot(Hmax, Vmax);
    for (let i = 0; i < 80; i++) { const mid = 0.5 * (lo + hi); if (val(mid * c, mid * s) < 1) lo = mid; else hi = mid; }
    Tult = 0.5 * (lo + hi);
  }
  return { f: H !== undefined && V !== undefined ? val(H, V) : undefined, Tult, a: ea, b: eb };
}

/**
 * Moment capacity of a mono-bucket under a horizontal load H applied at height e above the mudline
 * (M = H e) – simple rigid-rotation method. The bucket and its soil plug rotate as a rigid body about a
 * depth z0; the ultimate lateral soil reaction is
 *     clay:  p = Np su(z) D            (Np constant argument, or Murff–Hamilton when omitted);
 *     sand:  p = 3 Kp γ′ z D           (Broms 1964, Kp = tan²(π/4 + φ/2)).
 * Base shear and base bearing moment are ignored, so the result is a lower estimate for L/D < 1; it is a
 * screening method, not a substitute for a V–H–M envelope analysis.
 * @param {object} o
 * @param {'clay'|'sand'} o.soil
 * @param {number} o.D [m]  @param {number} o.L skirt length [m]
 * @param {number} o.e load eccentricity M/H above the mudline [m], ≥ 0
 * @param {number} [o.su0] [Pa]  @param {number} [o.k=0] [Pa/m]  @param {number} [o.NpConst] [-] (clay)
 * @param {number} [o.phi] friction angle [rad]  @param {number} [o.gamma] γ′ [N/m^3] (sand)
 * @returns {{H:number,M:number,z0:number}} ultimate H [N], M = H·e [N·m], rotation depth [m]
 */
export function bucketMomentCapacity({ soil, D, L, e, su0, k = 0, NpConst, phi, gamma }) {
  need('D', D, 1e-9);
  need('L', L, 1e-9);
  need('e', e, 0);
  let p;
  if (soil === 'clay') {
    need('su0', su0, 0);
    p = (z) => (NpConst !== undefined ? NpConst : murffHamiltonNp({ z, D, su0, k })) * (su0 + k * z) * D;
  } else if (soil === 'sand') {
    need('phi', phi, 0, PI / 2 - 1e-6);
    need('gamma', gamma, 0);
    const Kp = Math.tan(PI / 4 + phi / 2) ** 2;
    p = (z) => 3 * Kp * gamma * z * D;
  } else throw new RangeError(`soil must be 'clay' or 'sand' (got ${soil})`);
  const r = rigidPileLateralCapacity({ p, L, zLoad: -e });
  return { H: r.H, M: r.H * e, z0: r.z0 };
}

/**
 * Drained (sustained) tensile capacity of a caisson in sand: wall friction with the vertical stress next
 * to the wall reduced by the upward shear (Houlsby, Kelly & Byrne 2005):
 *     σ′v = γ′ Z [1 − exp(−z/Z)],
 *     V = W′ + γ′ Zo² [exp(−h/Zo) − 1 + h/Zo] (K tan δ)o πDo + γ′ Zi² [exp(−h/Zi) − 1 + h/Zi] (K tan δ)i πDi.
 * Under rapid loading the capacity is augmented by passive suction under the lid, limited by cavitation:
 *     V_rapid = W′ + outside friction + min(s_cav, s_available) π Di²/4,  s_cav = p_atm + γw hw − p_vap,
 * provided the soil plug weight plus inside friction can transmit it; otherwise the supplied limit governs.
 * @param {object} o
 * @param {number} o.Do [m]  @param {number} o.t [m]  @param {number} o.h embedded length [m]
 * @param {number} o.gamma γ′ [N/m^3]
 * @param {number} o.KtanDeltaO [-]  @param {number} o.KtanDeltaI [-]
 * @param {number} [o.m=1.5] outside stress-influence ratio [-]
 * @param {number} [o.W=0] submerged weight [N]
 * @param {number} [o.waterDepth] hw [m] for the cavitation limit (rapid loading)
 * @param {number} [o.suctionLimit=Infinity] additional cap on usable suction [Pa]
 * @param {number} [o.gammaW=10055] [N/m^3]  @param {number} [o.pAtm=101325] [Pa]  @param {number} [o.pVap=2300] [Pa]
 * @returns {{drained:number,rapid:number|undefined,frictionOut:number,frictionIn:number,sCav:number|undefined}} [N]
 */
export function tensileCapacitySand({ Do, t, h, gamma, KtanDeltaO, KtanDeltaI, m = 1.5, W = 0, waterDepth, suctionLimit = Infinity, gammaW = 10055, pAtm = 101325, pVap = 2300 }) {
  need('Do', Do, 1e-9);
  need('t', t, 1e-9, Do / 2);
  need('h', h, 1e-9);
  need('gamma', gamma, 0);
  need('KtanDeltaO', KtanDeltaO, 1e-9);
  need('KtanDeltaI', KtanDeltaI, 1e-9);
  need('m', m, 1 + 1e-9);
  const Di = Do - 2 * t;
  const Zo = (Do * (m * m - 1)) / (4 * KtanDeltaO);
  const Zi = Di / (4 * KtanDeltaI);
  const fr = (Z, kt, d) => gamma * Z * Z * (Math.expm1(-h / Z) + h / Z) * kt * PI * d;
  const frictionOut = fr(Zo, KtanDeltaO, Do);
  const frictionIn = fr(Zi, KtanDeltaI, Di);
  let rapid;
  let sCav;
  if (waterDepth !== undefined) {
    sCav = pAtm + gammaW * need('waterDepth', waterDepth, 0) - pVap;
    rapid = W + frictionOut + (Math.min(sCav, suctionLimit) * PI * Di * Di) / 4;
  }
  return { drained: W + frictionOut + frictionIn, rapid, frictionOut, frictionIn, sCav };
}
