/**
 * Mooring anchors: drag-embedment anchors and embedded chain, plate anchors, anchor piles, gravity anchors,
 * dynamically installed (torpedo) anchors, and design factors.
 *
 * SI units: m, N, Pa, N/m^3, kg, s, rad. Angles of the mooring line are measured above the horizontal.
 *
 * References
 *   Neubecker S.R. & Randolph M.F. (1995) Profile and frictional capacity of embedded anchor chains,
 *     J. Geotech. Eng. 121(11), 797–803.
 *   Neubecker S.R. & Randolph M.F. (1996) The kinematic behaviour of drag anchors in clay, Can. Geotech. J.
 *     33, 584–594.
 *   Degenkamp G. & Dutta A. (1989) Soil resistances to embedded anchor chain in soft clay, JGE 115(10).
 *   Rowe R.K. & Davis E.H. (1982) The behaviour of anchor plates in clay, Geotechnique 32(1);
 *   Martin C.M. & Randolph M.F. (2001) Applications of the lower and upper bound theorems of plasticity to
 *     collapse of circular foundations, IACMAG.
 *   True D.G. (1976) Undrained vertical penetration into ocean bottom soils, PhD thesis, UC Berkeley;
 *   O'Loughlin C.D., Richardson M.D., Randolph M.F. & Gaudin C. (2013) Geotechnique 63(11).
 *   API RP 2SK (3rd ed., 2005) Tables 6–8; DNV-OS-E301 (2021); DNV-RP-E301 / E302 / E303.
 */
import { brent } from '../math/roots.js';
import { murffHamiltonNp, rigidPileLateralCapacity } from './suction.js';

const PI = Math.PI;
const G0 = 9.80665;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

/**
 * Indicative coefficients for the drag-anchor power law UHC[kN] = a · (W[t])^b.
 * These are order-of-magnitude values of the kind printed on manufacturers' design charts for modern
 * high-holding-power fluke anchors (exponent 0.92; a ≈ 161 very soft clay, 229 medium clay, 324 sand and
 * hard clay) and are NOT a substitute for the chart of the actual anchor or the curves of API RP 2SK
 * Annex. They are editable defaults only.
 */
export const DRAG_ANCHOR_INDICATIVE = Object.freeze({
  'high-efficiency / very soft clay': { a: 161, b: 0.92 },
  'high-efficiency / medium clay': { a: 229, b: 0.92 },
  'high-efficiency / sand or hard clay': { a: 324, b: 0.92 },
});

/**
 * Drag-embedment anchor ultimate holding capacity by the empirical power law used on anchor design charts
 * (NCEL 1987; API RP 2SK):  UHC[kN] = a (W[t])^b, with W the anchor mass in tonnes.
 * @param {object} o
 * @param {number} o.mass anchor mass [kg]
 * @param {number} o.a coefficient [kN per t^b] (see DRAG_ANCHOR_INDICATIVE – indicative only)
 * @param {number} o.b exponent [-]
 * @returns {{UHC:number,efficiency:number}} UHC [N]; efficiency = UHC / anchor weight in air [-]
 */
export function dragAnchorUhc({ mass, a, b }) {
  need('mass', mass, 1e-9);
  need('a', a, 0);
  need('b', b, 0.3, 1.5);
  const UHC = 1000 * a * (mass / 1000) ** b;
  return { UHC, efficiency: UHC / (mass * G0) };
}

function qIntegral(za, { d, En, Nc, su0, k }) {
  return En * d * Nc * (su0 * za + 0.5 * k * za * za);
}

/**
 * Embedded mooring chain or wire between the mudline and the anchor padeye (Neubecker & Randolph 1995):
 *     (Ta/2)(θa² − θm²) = ∫0^za Q dz = za Q̄,     Q = En d Nc su(z)   (normal soil resistance per length),
 *     T/Ta = exp[μ (θa − θ)]   →   Tm = Ta exp[μ (θa − θm)].
 * Small-angle solution; adequate for θa up to about 60° and negligible line weight.
 * Typical parameters: chain En = 2.5, wire En = 1; Nc = 7.6 (5.1 at the surface to 7.6 at depth;
 * Degenkamp & Dutta 1989); friction ratio μ = 0.4–0.6 for chain in clay, 0.1–0.3 for wire.
 * @param {object} o
 * @param {number} o.Ta line tension at the padeye [N]
 * @param {number} o.za padeye depth [m]
 * @param {number} o.d nominal chain bar or wire diameter [m]
 * @param {number} [o.En=2.5] effective width multiplier for normal resistance [-]
 * @param {number} [o.Nc=7.6] bearing factor [-]
 * @param {number} o.su0 mudline strength [Pa]
 * @param {number} [o.k=0] strength gradient [Pa/m]
 * @param {number} [o.mu=0.4] friction ratio F/Q [-]
 * @param {number} [o.thetaM=0] line angle at the mudline [rad]
 * @returns {{thetaA:number,Tm:number,Qbar:number,ratio:number}} padeye angle [rad]; mudline tension [N];
 *   average normal resistance [N/m]; Tm/Ta [-]
 */
export function embeddedChain({ Ta, za, d, En = 2.5, Nc = 7.6, su0, k = 0, mu = 0.4, thetaM = 0 }) {
  need('Ta', Ta, 1e-9);
  need('za', za, 0);
  need('d', d, 1e-9);
  need('su0', su0, 0);
  need('k', k, 0);
  need('mu', mu, 0, 2);
  need('thetaM', thetaM, 0, PI / 2);
  const zQ = qIntegral(za, { d, En, Nc, su0, k });
  const thetaA = Math.sqrt(thetaM * thetaM + (2 * zQ) / Ta);
  const ratio = Math.exp(mu * (thetaA - thetaM));
  return { thetaA, Tm: Ta * ratio, Qbar: za > 0 ? zQ / za : 0, ratio };
}

/**
 * Inverse of embeddedChain: padeye tension and angle for a given mudline tension.
 * @param {object} o as embeddedChain with Tm [N] in place of Ta
 * @returns {{Ta:number,thetaA:number,ratio:number}}
 */
export function embeddedChainFromMudline(o) {
  need('Tm', o.Tm, 1e-9);
  const g = (Ta) => embeddedChain({ ...o, Ta }).Tm - o.Tm;
  // Tm(Ta) = Ta exp(μ θa) has a minimum where θa ≈ 2/μ (far outside the small-angle range); the physical
  // root is the one on the rising branch above that minimum.
  const { d, En = 2.5, Nc = 7.6, su0, k = 0, mu = 0.4, za } = o;
  const c = 2 * qIntegral(za, { d, En, Nc, su0, k });
  const lo = Math.max((mu * mu * c) / 4, o.Tm * 1e-9);
  if (g(lo) > 0) throw new RangeError('Tm is too small to reach the padeye: the line is held by soil friction');
  if (g(o.Tm) === 0) return { Ta: o.Tm, thetaA: o.thetaM || 0, ratio: 1 };
  const Ta = brent(g, lo, o.Tm, { tol: 1e-12 });
  const r = embeddedChain({ ...o, Ta });
  return { Ta, thetaA: r.thetaA, ratio: r.ratio };
}

/**
 * Drag-anchor embedment trajectory in clay – simplified kinematic model after Neubecker & Randolph (1996).
 * Assumptions (documented simplification):
 *   • the anchor advances parallel to its fluke; the soil resistance, of magnitude Ta = Ne Af su(z), acts at
 *     a fixed angle θw to the fluke, so the fluke dips at β = θw − θa below horizontal;
 *   • θa follows from the embedded-chain equation with zero mudline angle;
 *   • hence dz/dx = tan(θw − θa(z)), and embedment stops at the ultimate depth where θa = θw.
 * Ne (capacity factor on fluke area, form factor × Nc) and θw are anchor-specific arguments
 * (N&R report f ≈ 1.4–1.6 with Nc = 9 and θw of roughly 25–45° depending on the fluke–shank angle).
 * @param {object} o
 * @param {number} o.Af fluke plan area [m^2]
 * @param {number} o.Ne anchor capacity factor Ta/(Af su) [-]
 * @param {number} o.thetaW resultant-to-fluke angle [rad]
 * @param {number} o.d chain/wire diameter [m]
 * @param {number} [o.En=2.5]  @param {number} [o.Nc=7.6]  @param {number} [o.mu=0.4]
 * @param {number} o.su0 [Pa]  @param {number} [o.k=0] [Pa/m]
 * @param {number} [o.fraction=0.95] fraction of the ultimate depth to which the trajectory is integrated
 * @param {number} [o.n=400] integration steps
 * @param {number} [o.zMax=200] search limit for the ultimate depth [m]
 * @returns {{zUlt:number,TaUlt:number,TmUlt:number,dragDistance:number,x:number[],z:number[],Ta:number[],
 *   Tm:number[],thetaA:number[]}} zUlt ultimate fluke depth [m]; padeye and mudline tension there [N];
 *   dragDistance to reach fraction·zUlt [m]; trajectory arrays
 */
export function dragAnchorTrajectory({ Af, Ne, thetaW, d, En = 2.5, Nc = 7.6, mu = 0.4, su0, k = 0, fraction = 0.95, n = 400, zMax = 200 }) {
  need('Af', Af, 1e-9);
  need('Ne', Ne, 1e-9);
  need('thetaW', thetaW, 1e-6, PI / 2 - 1e-6);
  need('d', d, 1e-9);
  need('su0', su0, 0);
  need('k', k, 0);
  need('fraction', fraction, 0.01, 0.999);
  if (su0 === 0 && k === 0) throw new RangeError('su0 and k cannot both be zero');
  const su = (z) => su0 + k * z;
  const thA = (z) => (z <= 0 ? 0 : Math.sqrt((2 * qIntegral(z, { d, En, Nc, su0, k })) / (Ne * Af * su(z))));
  const g = (z) => thA(z) - thetaW;
  if (g(zMax) < 0) throw new RangeError('ultimate depth exceeds zMax: check Ne, Af and thetaW');
  const zUlt = brent(g, 1e-9, zMax, { tol: 1e-12 });
  const zEnd = fraction * zUlt;
  const dz = zEnd / n;
  const x = [0];
  const z = [0];
  const Ta = [Ne * Af * su(0)];
  const Tm = [Ta[0]];
  const thetaA = [0];
  let xx = 0;
  for (let i = 1; i <= n; i++) {
    const zm = (i - 0.5) * dz;
    xx += dz / Math.tan(thetaW - thA(zm));
    const zi = i * dz;
    const th = thA(zi);
    const ta = Ne * Af * su(zi);
    x.push(xx); z.push(zi); Ta.push(ta); Tm.push(ta * Math.exp(mu * th)); thetaA.push(th);
  }
  const TaUlt = Ne * Af * su(zUlt);
  return { zUlt, TaUlt, TmUlt: TaUlt * Math.exp(mu * thetaW), dragDistance: xx, x, z, Ta, Tm, thetaA };
}

/**
 * Vertically loaded / plate anchor (VLA, SEPLA) capacity in clay:  T = Nc · su(z_final) · A · η,
 * with z_final = z_install − keyingLoss (loss of embedment while the plate rotates normal to the line;
 * typically 0.25–1.0 plate heights for SEPLAs – argument), and η a cyclic/set-up reduction (argument).
 * Deep-embedment bearing factors: 11.42 = 3π + 2 for a strip (exact, Rowe & Davis 1982); 12.42 (smooth) to
 * 13.11 (rough) for a circular plate (Martin & Randolph 2001). For shallow embedment supply `NcOfDepth`
 * (function of z/B) from the relevant chart; the deep value is used as a cap.
 * @param {object} o
 * @param {number} o.A plate area [m^2]
 * @param {number} o.B plate width (height) [m]
 * @param {number} o.zInstall installed depth of the plate centre [m]
 * @param {number} [o.keyingLoss=0] loss of embedment during keying [m]
 * @param {number} o.su0 [Pa]  @param {number} [o.k=0] [Pa/m]
 * @param {number} [o.Nc=11.42] deep bearing factor [-]
 * @param {(zOverB:number)=>number} [o.NcOfDepth] shallow bearing factor [-]
 * @param {number} [o.eta=1] reduction factor [-]
 * @returns {{T:number,zFinal:number,su:number,Nc:number}} capacity normal to the plate [N]; depth [m]; su [Pa]
 */
export function plateAnchorCapacity({ A, B, zInstall, keyingLoss = 0, su0, k = 0, Nc = 3 * PI + 2, NcOfDepth, eta = 1 }) {
  need('A', A, 0);
  need('B', B, 1e-9);
  need('zInstall', zInstall, 0);
  need('keyingLoss', keyingLoss, 0, zInstall);
  need('su0', su0, 0);
  need('eta', eta, 0, 2);
  const zFinal = zInstall - keyingLoss;
  const su = su0 + k * zFinal;
  const nc = NcOfDepth ? Math.min(NcOfDepth(zFinal / B), Nc) : Nc;
  return { T: nc * su * A * eta, zFinal, su, Nc: nc };
}

/**
 * Anchor pile capacity (self-contained rigid-pile solution).
 *   Clay: lateral reaction p = Np su D with Np from Murff & Hamilton (1993) or a constant; axial pull-out
 *         α ∫su dz π D + W′ (no reverse end bearing counted).
 *   Sand: lateral p = 3 Kp γ′ z D (Broms 1964); axial pull-out K tan δ ∫σ′v dz π D + W′ = ½ K tan δ γ′ L² π D + W′.
 * The lateral capacity is evaluated for the padeye depth by rigid rotation (rigidPileLateralCapacity in
 * suction.js); it applies to short stiff piles (L/D ≲ 6–10) – long flexible piles are governed by bending
 * and need a p–y analysis. Combined loading uses (H/Hult)^a + (V/Vult)^b = 1 with exponents as arguments
 * (default a = L/D + 0.5, b = L/(3D) + 4.5 after Supachawarote et al. 2004, derived for caissons in clay).
 * @param {object} o
 * @param {'clay'|'sand'} o.soil
 * @param {number} o.D [m]  @param {number} o.L embedded length [m]
 * @param {number} o.zPadeye padeye depth [m]
 * @param {number} [o.su0] [Pa]  @param {number} [o.k=0] [Pa/m]  @param {number} [o.alpha=0.65] [-]
 * @param {number} [o.NpConst] [-]
 * @param {number} [o.phi] [rad]  @param {number} [o.gamma] γ′ [N/m^3]  @param {number} [o.KtanDelta=0.4] [-]
 * @param {number} [o.W=0] submerged pile weight [N]
 * @param {number} [o.angle] line angle at the padeye above horizontal [rad]
 * @param {number} [o.a]  @param {number} [o.b] interaction exponents [-]
 * @returns {{Hult:number,Vult:number,Tult:number|undefined,z0:number,zOptimal:number}} [N, N, N, m, m]
 */
export function anchorPileCapacity({ soil, D, L, zPadeye, su0, k = 0, alpha = 0.65, NpConst, phi, gamma, KtanDelta = 0.4, W = 0, angle, a, b }) {
  need('D', D, 1e-9);
  need('L', L, 1e-9);
  need('zPadeye', zPadeye, 0, L);
  let p;
  let Vult;
  if (soil === 'clay') {
    need('su0', su0, 0);
    p = (z) => (NpConst !== undefined ? NpConst : murffHamiltonNp({ z, D, su0, k })) * (su0 + k * z) * D;
    Vult = alpha * (su0 * L + 0.5 * k * L * L) * PI * D + W;
  } else if (soil === 'sand') {
    need('phi', phi, 0, PI / 2 - 1e-6);
    need('gamma', gamma, 0);
    const Kp = Math.tan(PI / 4 + phi / 2) ** 2;
    p = (z) => 3 * Kp * gamma * z * D;
    Vult = 0.5 * KtanDelta * gamma * L * L * PI * D + W;
  } else throw new RangeError(`soil must be 'clay' or 'sand' (got ${soil})`);
  const r = rigidPileLateralCapacity({ p, L, zLoad: zPadeye });
  let Tult;
  if (angle !== undefined) {
    need('angle', angle, 0, PI / 2);
    const ea = a !== undefined ? a : L / D + 0.5;
    const eb = b !== undefined ? b : L / (3 * D) + 4.5;
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    let lo = 0;
    let hi = Math.hypot(r.H, Vult);
    for (let i = 0; i < 80; i++) {
      const mid = 0.5 * (lo + hi);
      if (((mid * c) / r.H) ** ea + ((mid * s) / Vult) ** eb < 1) lo = mid; else hi = mid;
    }
    Tult = 0.5 * (lo + hi);
  }
  return { Hult: r.H, Vult, Tult, z0: r.z0, zOptimal: r.zCentroid };
}

/**
 * Gravity (deadweight) anchor under an inclined line load T at angle θ above horizontal (statics):
 *   uplift       V = T sin θ ≤ W′;
 *   drained      H = T cos θ ≤ (W′ − V) tan δ + Hkey      →  Tult = (W′ tan δ + Hkey) / (cos θ + sin θ tan δ);
 *   undrained    H ≤ su A + Hkey                           →  Tult = min[(su A + Hkey)/cos θ, W′/sin θ].
 * Hkey is the passive resistance of shear keys/skirts (see slidingResistance in bearing.js).
 * @param {object} o
 * @param {number} o.W submerged weight [N]
 * @param {number} [o.angle=0] θ [rad]
 * @param {'drained'|'undrained'} o.mode
 * @param {number} [o.delta] interface friction angle [rad] (drained)
 * @param {number} [o.su] interface strength [Pa] and @param {number} [o.A] base area [m^2] (undrained)
 * @param {number} [o.Hkey=0] shear-key resistance [N]
 * @param {number} [o.T] applied line tension [N] for the factor of safety
 * @returns {{Tult:number,fs:number|undefined,upliftRatio:number|undefined}} Tult [N]
 */
export function gravityAnchorCapacity({ W, angle = 0, mode, delta, su, A, Hkey = 0, T }) {
  need('W', W, 0);
  need('angle', angle, 0, PI / 2);
  need('Hkey', Hkey, 0);
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  let Tult;
  if (mode === 'drained') {
    const td = Math.tan(need('delta', delta, 0, PI / 2 - 1e-6));
    Tult = (W * td + Hkey) / (c + s * td);
  } else if (mode === 'undrained') {
    need('su', su, 0);
    need('A', A, 0);
    Tult = Math.min(c > 1e-12 ? (su * A + Hkey) / c : Infinity, s > 1e-12 ? W / s : Infinity);
  } else throw new RangeError(`mode must be 'drained' or 'undrained' (got ${mode})`);
  return { Tult, fs: T !== undefined ? Tult / T : undefined, upliftRatio: T !== undefined && W > 0 ? (T * s) / W : undefined };
}

/**
 * Free-fall penetration of a torpedo / dynamically installed anchor in clay by integration of the equation
 * of motion (True 1976 framework; strain-rate treatment after O'Loughlin et al. 2013):
 *     (m + ma) dv/dt = W′ − Rf(v) [Nc su(z) Atip + α ∫ su dA_shaft] − ½ ρs Cd Atip v²,
 *     Rf = 1 + λ log10(γ̇/γ̇ref)  ('log')   or   (γ̇/γ̇ref)^β  ('power'),   γ̇ = v/D,   Rf ≥ 1.
 * The shaft area grows with penetration up to the anchor length; the hole above the anchor is assumed to
 * stay open (no resistance on the rear face). Classical RK4 with fixed step; integration stops when v ≤ 0.
 * Rate parameters λ (0.1–0.2 per log cycle) or β (0.05–0.17), α and Cd are arguments.
 * @param {object} o
 * @param {number} o.mass anchor mass [kg]
 * @param {number} [o.addedMass=0] soil added mass [kg]
 * @param {number} o.W submerged weight in soil (dry weight minus buoyancy) [N]
 * @param {number} o.D shaft diameter [m]
 * @param {number} o.Lanchor anchor length [m]
 * @param {number} [o.Atip] projected tip area [m^2] (default πD²/4)
 * @param {number} [o.perimeter] shaft (and fin) perimeter [m] (default πD)
 * @param {number} o.v0 impact velocity at the mudline [m/s]
 * @param {number} o.su0 [Pa]  @param {number} [o.k=0] [Pa/m]
 * @param {number} [o.Nc=12] tip bearing factor [-]
 * @param {number} [o.alpha=0.4] shaft adhesion factor during penetration [-]
 * @param {number} [o.Cd=0] inertial drag coefficient [-]
 * @param {number} [o.rhoS=1600] soil density [kg/m^3]
 * @param {'none'|'log'|'power'} [o.rateModel='none']
 * @param {number} [o.lambda=0.15]  @param {number} [o.beta=0.1]  @param {number} [o.gammaRef=0.1] reference strain rate [1/s]
 * @param {number} [o.dt=1e-3] time step [s]
 * @param {number} [o.tMax=30] [s]
 * @returns {{zTip:number,time:number,t:number[],z:number[],v:number[]}} final tip embedment [m]; time to rest
 *   [s]; decimated histories (≤ ~400 points)
 */
export function torpedoPenetration(o) {
  const { mass, addedMass = 0, W, D, Lanchor, v0, su0, k = 0, Nc = 12, alpha = 0.4, Cd = 0, rhoS = 1600, rateModel = 'none', lambda = 0.15, beta = 0.1, gammaRef = 0.1, dt = 1e-3, tMax = 30 } = o;
  need('mass', mass, 1e-9);
  need('W', W, 0);
  need('D', D, 1e-9);
  need('Lanchor', Lanchor, 1e-9);
  need('v0', v0, 0);
  need('su0', su0, 0);
  need('dt', dt, 1e-9);
  const Atip = o.Atip !== undefined ? need('Atip', o.Atip, 0) : (PI * D * D) / 4;
  const per = o.perimeter !== undefined ? need('perimeter', o.perimeter, 0) : PI * D;
  const M = mass + addedMass;
  const rf = (v) => {
    if (rateModel === 'none' || v <= 0) return 1;
    const r = v / D / gammaRef;
    if (rateModel === 'log') return Math.max(1 + lambda * Math.log10(Math.max(r, 1e-12)), 1);
    if (rateModel === 'power') return Math.max(r ** beta, 1);
    throw new RangeError(`rateModel must be 'none', 'log' or 'power' (got ${rateModel})`);
  };
  const resist = (z, v) => {
    const zt = Math.max(z, 0);
    const top = Math.max(zt - Lanchor, 0);
    const len = zt - top;
    const suAvg = su0 + 0.5 * k * (zt + top);
    return rf(v) * (Nc * (su0 + k * zt) * Atip + alpha * suAvg * per * len) + 0.5 * rhoS * Cd * Atip * v * v;
  };
  const acc = (z, v) => (W - resist(z, v)) / M;
  let z = 0;
  let v = v0;
  let t = 0;
  const ts = [0];
  const zs = [0];
  const vs = [v0];
  const every = Math.max(1, Math.floor(tMax / dt / 400));
  let step = 0;
  while (t < tMax) {
    const k1z = v, k1v = acc(z, v);
    const k2z = v + 0.5 * dt * k1v, k2v = acc(z + 0.5 * dt * k1z, v + 0.5 * dt * k1v);
    const k3z = v + 0.5 * dt * k2v, k3v = acc(z + 0.5 * dt * k2z, v + 0.5 * dt * k2v);
    const k4z = v + dt * k3v, k4v = acc(z + dt * k3z, v + dt * k3v);
    const zn = z + (dt / 6) * (k1z + 2 * k2z + 2 * k3z + k4z);
    const vn = v + (dt / 6) * (k1v + 2 * k2v + 2 * k3v + k4v);
    if (vn <= 0) {
      // linear interpolation to the instant of zero velocity
      const f = v / (v - vn);
      z += f * (zn - z);
      t += f * dt;
      v = 0;
      break;
    }
    z = zn; v = vn; t += dt; step++;
    if (step % every === 0) { ts.push(t); zs.push(z); vs.push(v); }
  }
  if (v > 0) throw new Error('torpedo anchor did not come to rest within tMax (weight exceeds soil resistance?)');
  ts.push(t); zs.push(z); vs.push(0);
  return { zTip: z, time: t, t: ts, z: zs, v: vs };
}

/**
 * Energy-method estimate of torpedo penetration with a constant average rate factor:
 *     ½ (m + ma) v0² + W′ z = Rf,avg ∫0^z [Nc su Atip + α su,avg(z) · perimeter · min(z, L)] dz,
 * solved for z. Agrees with torpedoPenetration when Cd = 0 and Rf is constant.
 * @param {object} o as torpedoPenetration (mass, addedMass, W, D, Lanchor, Atip, perimeter, v0, su0, k, Nc,
 *   alpha) plus
 * @param {number} [o.rateFactor=1] average strain-rate factor Rf,avg [-]
 * @returns {number} tip embedment [m]
 */
export function torpedoPenetrationEnergy(o) {
  const { mass, addedMass = 0, W, D, Lanchor, v0, su0, k = 0, Nc = 12, alpha = 0.4, rateFactor = 1 } = o;
  need('mass', mass, 1e-9);
  need('v0', v0, 0);
  need('rateFactor', rateFactor, 1e-9);
  const Atip = o.Atip !== undefined ? o.Atip : (PI * D * D) / 4;
  const per = o.perimeter !== undefined ? o.perimeter : PI * D;
  const R = (z) => {
    const top = Math.max(z - Lanchor, 0);
    return rateFactor * (Nc * (su0 + k * z) * Atip + alpha * (su0 + 0.5 * k * (z + top)) * per * (z - top));
  };
  const work = (z) => {
    const n = 400;
    const h = z / n;
    let s = 0;
    for (let i = 0; i <= n; i++) s += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * R(i * h);
    return (s * h) / 3;
  };
  const E0 = 0.5 * (mass + addedMass) * v0 * v0;
  const g = (z) => E0 + W * z - work(z);
  let hi = Math.max(Lanchor, 1);
  let guard = 0;
  while (g(hi) > 0 && guard++ < 60) hi *= 1.6;
  if (g(hi) > 0) throw new Error('energy balance has no solution: weight exceeds soil resistance');
  if (E0 === 0 && g(1e-9) <= 0) return 0;
  return brent(g, 1e-9, hi, { tol: 1e-10 });
}

/**
 * Pull-out capacity of a torpedo / dynamically installed anchor in clay after set-up (axial):
 *     T = W′ + α ∫ su dA_shaft + Nc su Aprojected   (reverse end bearing on the projected area).
 * @param {object} o
 * @param {number} o.zTip tip depth [m]  @param {number} o.Lanchor [m]  @param {number} o.D [m]
 * @param {number} [o.perimeter] [m]  @param {number} [o.Aend] projected end area [m^2]
 * @param {number} o.su0 [Pa]  @param {number} [o.k=0]
 * @param {number} [o.alpha=0.65] long-term adhesion factor [-]  @param {number} [o.Nc=9] [-]
 * @param {number} [o.W=0] submerged weight [N]
 * @returns {{T:number,friction:number,endBearing:number}} [N]
 */
export function torpedoCapacity({ zTip, Lanchor, D, perimeter, Aend, su0, k = 0, alpha = 0.65, Nc = 9, W = 0 }) {
  need('zTip', zTip, 0);
  need('Lanchor', Lanchor, 1e-9);
  need('D', D, 1e-9);
  const per = perimeter !== undefined ? perimeter : PI * D;
  const Ae = Aend !== undefined ? Aend : (PI * D * D) / 4;
  const top = Math.max(zTip - Lanchor, 0);
  const friction = alpha * (su0 + 0.5 * k * (zTip + top)) * per * (zTip - top);
  const endBearing = Nc * (su0 + k * top) * Ae;
  return { T: W + friction + endBearing, friction, endBearing };
}

/**
 * Design anchor resistance from the static capacity: R = R_static · setUp · cyclic, compared with the
 * factored line tension; and the installation proof-load check T_proof ≥ proofFactor · T_design.
 * Set-up (consolidation) and cyclic degradation factors are site- and anchor-specific arguments
 * (e.g. DNV-RP-E301: set-up 1.0–1.5 over time in soft clay; cyclic factor commonly 0.8–1.2 including
 * loading-rate effects).
 * @param {object} o
 * @param {number} o.staticCapacity [N]
 * @param {number} [o.setUpFactor=1] [-]  @param {number} [o.cyclicFactor=1] [-]
 * @param {number} o.tension characteristic or maximum line tension at the anchor [N]
 * @param {number} [o.safetyFactor=1] required factor of safety (or γ load × γm) [-]
 * @param {number} [o.proofLoad] installation proof tension [N]
 * @param {number} [o.proofFactor=1] required proof load / tension [-]
 * @returns {{resistance:number,fs:number,uc:number,ok:boolean,proofOk:boolean|undefined}}
 */
export function anchorDesignCheck({ staticCapacity, setUpFactor = 1, cyclicFactor = 1, tension, safetyFactor = 1, proofLoad, proofFactor = 1 }) {
  need('staticCapacity', staticCapacity, 0);
  need('setUpFactor', setUpFactor, 0);
  need('cyclicFactor', cyclicFactor, 0);
  need('tension', tension, 1e-12);
  need('safetyFactor', safetyFactor, 1e-9);
  const resistance = staticCapacity * setUpFactor * cyclicFactor;
  const fs = resistance / tension;
  return { resistance, fs, uc: safetyFactor / fs, ok: fs >= safetyFactor, proofOk: proofLoad === undefined ? undefined : proofLoad >= proofFactor * tension };
}

/**
 * Default anchor safety factors – editable, to be confirmed against the edition of the standard in force.
 *   API RP 2SK (3rd ed.) Tables 6–8, dynamic analysis: factor of safety on anchor holding capacity.
 *   DNV-OS-E301: partial factors on mean and dynamic tension (consequence class 1 / 2, ULS, dynamic
 *   analysis) – anchor resistance material factors are given per anchor type in DNV-RP-E301/E302/E303.
 */
export const ANCHOR_SAFETY_FACTORS = Object.freeze({
  'API RP 2SK': {
    dragAnchor: { permanent: { intact: 1.5, damaged: 1.0 }, mobile: { intact: 0.8, damaged: null } },
    pileOrSuctionLateral: { permanent: { intact: 1.6, damaged: 1.2 }, mobile: { intact: 1.2, damaged: 1.0 } },
    pileOrSuctionAxial: { permanent: { intact: 2.0, damaged: 1.5 }, mobile: { intact: 1.5, damaged: 1.2 } },
    plateAnchor: { permanent: { intact: 2.0, damaged: 1.5 }, mobile: { intact: 1.5, damaged: 1.2 } },
  },
  'DNV-OS-E301': {
    uls: { class1: { gammaMean: 1.1, gammaDyn: 1.5 }, class2: { gammaMean: 1.4, gammaDyn: 2.1 } },
    als: { class1: { gammaMean: 1.0, gammaDyn: 1.1 }, class2: { gammaMean: 1.0, gammaDyn: 1.25 } },
  },
});

/**
 * DNV-OS-E301 design tension Td = γmean · Tmean + γdyn · Tdyn.
 * @param {object} o Tmean [N], Tdyn [N], gammaMean [-], gammaDyn [-]
 * @returns {number} Td [N]
 */
export function dnvDesignTension({ Tmean, Tdyn, gammaMean, gammaDyn }) {
  need('Tmean', Tmean, 0);
  need('Tdyn', Tdyn, 0);
  need('gammaMean', gammaMean, 0);
  need('gammaDyn', gammaDyn, 0);
  return gammaMean * Tmean + gammaDyn * Tdyn;
}
