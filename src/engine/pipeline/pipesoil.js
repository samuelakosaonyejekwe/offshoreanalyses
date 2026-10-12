/**
 * Pipe–soil interaction: embedment, axial and lateral resistance, mobilisation curves, friction-factor
 * ranges, soil stiffness for free-span analysis, uplift resistance of buried pipe and cover quantities.
 *
 * Models follow the format of DNV-RP-F114 (pipe–soil interaction for submarine pipelines), DNV-RP-F105
 * (soil stiffness) and DNV-RP-F110 / Schaminée et al. (uplift). Where an empirical coefficient is not built
 * in it is an argument with its typical range stated. All quantities SI; angles in radians unless `Deg`.
 */

import { erfinv } from '../math/special.js';
import { passiveResistance } from './stability.js';

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/**
 * Geometry of a pipe embedded by z into a flat seabed: contact width, half contact angle, embedded
 * cross-section area and the wedging factor ζ = 2 sin β/(β + sin β cos β).
 * cos β = 1 − 2z/D for z ≤ D/2 (β = π/2 beyond); B = 2√(Dz − z²) (B = D beyond);
 * A_bm = (D²/4)(β − sin β cos β) for z ≤ D/2, πD²/8 + D(z − D/2) beyond.
 *
 * Reference: circle-segment geometry; wedging factor per DNV-RP-F114 (2017) §4.3 (White & Randolph 2007):
 * the integral of the normal contact stress around the arc exceeds the vertical load by ζ (1.0 at zero
 * embedment, 4/π = 1.27 at half a diameter).
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.z embedment of the invert below the seabed [m]
 * @returns {{B:number, beta:number, Abm:number, zeta:number, arc:number}}
 *   B [m]; beta [rad]; Abm [m²]; zeta [-]; arc contact arc length D·β [m]
 */
export function contactGeometry({ D, z }) {
  pos('D', D); nonNeg('z', z);
  if (z >= D / 2) {
    return { B: D, beta: Math.PI / 2, Abm: (Math.PI * D * D) / 8 + D * (z - D / 2), zeta: 4 / Math.PI, arc: (D * Math.PI) / 2 };
  }
  const beta = Math.acos(1 - (2 * z) / D);
  const s = Math.sin(beta);
  const c = Math.cos(beta);
  return { B: D * s, beta, Abm: ((D * D) / 4) * (beta - s * c), zeta: beta > 0 ? (2 * s) / (beta + s * c) : 1, arc: D * beta };
}

/**
 * Undrained vertical bearing resistance of a partially embedded pipe on clay:
 * Q_v = { min[ 6 (z/D)^0.25 ; 3.4 (10 z/D)^0.5 ] + 1.5 γ' A_bm/(D s_u) } · D · s_u, s_u taken at the invert.
 *
 * Reference: DNV-RP-F114 (2017) §4.2 (undrained penetration resistance, the SAFEBUCK-type closed form
 * including soil buoyancy with the heave factor 1.5).
 * Validity: z/D up to about 1; monotonic loading; remoulding / lay effects are handled by the caller
 * through the load (lay factor) or the strength (remoulded s_u).
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.z embedment [m]
 * @param {number} o.su0 undrained shear strength at the seabed [Pa]
 * @param {number} [o.suGrad=0] strength gradient with depth [Pa/m]
 * @param {number} [o.gammaSub=0] submerged unit weight of the soil [N/m³]
 * @returns {{Qv:number, suInvert:number, bearingFactor:number, buoyancyTerm:number}} Qv [N/m]
 */
export function verticalBearingUndrained({ D, z, su0, suGrad = 0, gammaSub = 0 }) {
  pos('D', D); nonNeg('z', z); nonNeg('su0', su0); nonNeg('suGrad', suGrad); nonNeg('gammaSub', gammaSub);
  const su = su0 + suGrad * z;
  if (!(su > 0)) {
    if (z === 0) return { Qv: 0, suInvert: 0, bearingFactor: 0, buoyancyTerm: 0 };
    throw new RangeError('su0 and suGrad cannot both be zero');
  }
  const r = z / D;
  const nc = Math.min(6 * r ** 0.25, 3.4 * Math.sqrt(10 * r));
  const { Abm } = contactGeometry({ D, z });
  const buoy = 1.5 * gammaSub * Abm;
  return { Qv: nc * D * su + buoy, suInvert: su, bearingFactor: nc, buoyancyTerm: buoy };
}

/**
 * Drained vertical bearing resistance of a partially embedded pipe on sand, as a strip footing of the
 * contact width B: Q_v = ½ γ' N_γ B² + z_0 γ' N_q d_q B,
 * N_q = e^(π tan φ) tan²(45° + φ/2); N_γ = 1.5 (N_q − 1) tan φ (Brinch Hansen) or 2 (N_q + 1) tan φ (Vesić);
 * d_q = 1 + 1.2 (z_0/B) tan φ (1 − sin φ)².
 *
 * Reference: DNV-RP-F114 (2017) §4.2 (drained penetration resistance); bearing capacity factors per
 * Brinch Hansen (1970) / Vesić (1973). The surcharge depth z_0 is the depth of the equivalent flat footing:
 * 'mean' takes A_bm/B (default), 'invert' takes z — the choice is exposed because the RP wording was not
 * available for checking.
 * Validity: fully drained response; φ between about 25° and 45°.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.z embedment [m]
 * @param {number} o.phiDeg friction angle [deg]
 * @param {number} o.gammaSub submerged unit weight [N/m³]
 * @param {'hansen'|'vesic'} [o.nGamma='hansen'] N_γ expression
 * @param {'mean'|'invert'} [o.surchargeDepth='mean']
 * @returns {{Qv:number, Nq:number, Ngamma:number, dq:number, B:number, z0:number}} Qv [N/m]
 */
export function verticalBearingDrained({ D, z, phiDeg, gammaSub, nGamma = 'hansen', surchargeDepth = 'mean' }) {
  pos('D', D); nonNeg('z', z); pos('gammaSub', gammaSub);
  if (!(phiDeg > 0 && phiDeg < 50)) throw new RangeError('phiDeg must be in (0, 50)');
  const phi = (phiDeg * Math.PI) / 180;
  const tp = Math.tan(phi);
  const Nq = Math.exp(Math.PI * tp) * Math.tan(Math.PI / 4 + phi / 2) ** 2;
  let Ngamma;
  if (nGamma === 'hansen') Ngamma = 1.5 * (Nq - 1) * tp;
  else if (nGamma === 'vesic') Ngamma = 2 * (Nq + 1) * tp;
  else throw new RangeError('nGamma must be "hansen" or "vesic"');
  const { B, Abm } = contactGeometry({ D, z });
  if (B === 0) return { Qv: 0, Nq, Ngamma, dq: 1, B: 0, z0: 0 };
  const z0 = surchargeDepth === 'invert' ? z : Abm / B;
  const dq = 1 + 1.2 * (z0 / B) * tp * (1 - Math.sin(phi)) ** 2;
  return { Qv: 0.5 * gammaSub * Ngamma * B * B + z0 * gammaSub * Nq * dq * B, Nq, Ngamma, dq, B, z0 };
}

/**
 * Pipe embedment under a vertical load: the depth at which the bearing resistance equals the load, with the
 * touchdown-zone load concentration (lay factor) and an optional dynamic embedment factor.
 * z_static solves Q_v(z) = k_lay · W; z = k_dyn · z_static.
 *
 * Reference: DNV-RP-F114 (2017) §4.2. k_lay (maximum touchdown contact force over submerged weight) follows
 * from the lay catenary — typically 1 to 3, higher for stiff pipe in shallow water; k_dyn accounts for lay
 * vessel motion and is typically 1 to 3 on clay (site-specific; argument).
 * Validity: solution searched for 0 < z ≤ zMax (default 2D). If the bearing resistance at zMax is still below
 * the load (very soft clay, heavy flooded pipe, large lay factor) the pipe is predicted to sink beyond the range
 * of the bearing model: the static embedment is then capped at zMax and `capped` is true — the pipe is at least
 * that deep (fully buried for zMax ≥ D) and the value is a lower bound, not an equilibrium. The dynamic factor is
 * not applied on top of a capped depth.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.W submerged weight in the condition considered [N/m]
 * @param {'undrained'|'drained'} o.model
 * @param {number} [o.kLay=1] lay factor [-]
 * @param {number} [o.kDyn=1] dynamic embedment factor on depth [-]
 * @param {number} [o.su0] [Pa] (undrained)
 * @param {number} [o.suGrad=0] [Pa/m]
 * @param {number} [o.gammaSub=0] [N/m³]
 * @param {number} [o.phiDeg] [deg] (drained)
 * @param {'hansen'|'vesic'} [o.nGamma='hansen']
 * @param {'mean'|'invert'} [o.surchargeDepth='mean']
 * @param {number} [o.zMax] search limit [m]
 * @returns {{z:number, zStatic:number, zOverD:number, V:number, geometry:object, capped:boolean,
 *   bearingAtCap:(number|null)}} depths [m]; V design vertical load [N/m]; capped = the load exceeds the bearing
 *   resistance at zMax and the depths are held at zMax; bearingAtCap = that resistance [N/m] (null when not capped)
 */
export function embedment({ D, W, model, kLay = 1, kDyn = 1, su0, suGrad = 0, gammaSub = 0, phiDeg, nGamma = 'hansen',
  surchargeDepth = 'mean', zMax }) {
  pos('D', D); pos('W', W); pos('kLay', kLay); pos('kDyn', kDyn);
  const V = kLay * W;
  let Q;
  if (model === 'undrained') Q = (z) => verticalBearingUndrained({ D, z, su0, suGrad, gammaSub }).Qv;
  else if (model === 'drained') Q = (z) => verticalBearingDrained({ D, z, phiDeg, gammaSub, nGamma, surchargeDepth }).Qv;
  else throw new RangeError('model must be "undrained" or "drained"');
  let hi = zMax === undefined ? 2 * D : pos('zMax', zMax);
  const qCap = Q(hi);
  if (!(qCap >= V)) {
    if (!Number.isFinite(qCap)) throw new RangeError('the bearing resistance is not finite: check the soil arguments');
    return { z: hi, zStatic: hi, zOverD: hi / D, V, geometry: contactGeometry({ D, z: hi }), capped: true, bearingAtCap: qCap };
  }
  let lo = 0;
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi);
    if (Q(mid) < V) lo = mid; else hi = mid;
  }
  const zStatic = 0.5 * (lo + hi);
  const z = kDyn * zStatic;
  return { z, zStatic, zOverD: z / D, V, geometry: contactGeometry({ D, z }), capped: false, bearingAtCap: null };
}

/**
 * Axial pipe–soil resistance (breakout and residual).
 * Undrained (total stress): F_A = α s_u · D β (interface shear over the contact arc).
 * Drained (effective stress): F_A = μ W ζ, ζ = 2 sin β/(β + sin β cos β).
 * Residual = residualRatio × breakout.
 *
 * Reference: DNV-RP-F114 (2017) §4.3. α (interface adhesion factor, about 0.2–1.0 depending on sensitivity
 * and coating roughness), μ = tan δ (interface friction, about 0.2–0.7) and the residual ratio are
 * site-specific arguments.
 *
 * @param {object} o
 * @param {'undrained'|'drained'} o.model
 * @param {number} o.D [m]
 * @param {number} o.z embedment [m]
 * @param {number} [o.W] vertical contact force [N/m] (drained)
 * @param {number} [o.mu] interface friction coefficient [-] (drained)
 * @param {number} [o.alpha] adhesion factor [-] (undrained)
 * @param {number} [o.su] undrained strength at the interface [Pa] (undrained)
 * @param {number} [o.residualRatio=1] residual / breakout [-]
 * @returns {{breakout:number, residual:number, zeta:number, equivalentFriction:(number|null)}}
 *   resistances [N/m]; equivalentFriction = breakout/W when W is given
 */
export function axialResistance({ model, D, z, W, mu, alpha, su, residualRatio = 1 }) {
  const geo = contactGeometry({ D, z });
  pos('residualRatio', residualRatio);
  let breakout;
  if (model === 'undrained') breakout = pos('alpha', alpha) * pos('su', su) * geo.arc;
  else if (model === 'drained') breakout = pos('mu', mu) * pos('W', W) * geo.zeta;
  else throw new RangeError('model must be "undrained" or "drained"');
  return { breakout, residual: breakout * residualRatio, zeta: geo.zeta, equivalentFriction: W > 0 ? breakout / W : null };
}

/**
 * Lateral pipe–soil resistance: breakout and residual.
 *
 * model 'coulomb-passive' (first-principles): F_L = μ W + passive resistance of the soil ahead of the embedded
 *   depth z treated as a smooth vertical wall — clay ½ γ' z² + 2 s_u z (Rankine, undrained), sand
 *   ½ K_p γ' z², K_p = tan²(45° + φ/2); a berm from earlier cycles is included by `bermFactor` multiplying
 *   the passive term.
 * model 'verley': Coulomb friction plus the empirical passive resistance of Verley & Sotberg (sand) /
 *   Verley & Lund (clay) as coded in DNV-RP-F109 (see stability.passiveResistance).
 * Residual: μ_res W (+ bermFactorResidual × passive term when given).
 *
 * Reference: Coulomb / Rankine earth pressure; DNV-RP-F109 (2021) §3.4.6; format of DNV-RP-F114 §4.4.
 * μ: 0.2 clay, 0.6 sand in F109; μ_res is site-specific (argument, defaults to μ).
 * Validity: monotonic breakout from the as-laid embedment; cyclic berm growth only through the factors.
 *
 * @param {object} o
 * @param {'coulomb-passive'|'verley'} [o.model='coulomb-passive']
 * @param {'sand'|'clay'} o.soil
 * @param {number} o.D [m]
 * @param {number} o.z embedment [m]
 * @param {number} o.W vertical contact force [N/m]
 * @param {number} o.mu Coulomb friction coefficient [-]
 * @param {number} [o.muResidual=o.mu] [-]
 * @param {number} [o.su] [Pa] (clay)
 * @param {number} [o.gammaSub=0] [N/m³]
 * @param {number} [o.gammaDry] dry unit weight for the Verley clay model [N/m³]
 * @param {number} [o.phiDeg] [deg] (sand)
 * @param {number} [o.bermFactor=1] multiplier on the passive term at breakout [-]
 * @param {number} [o.bermFactorResidual=0] multiplier on the passive term at residual [-]
 * @returns {{breakout:number, residual:number, passive:number, frictionBreakout:number, frictionResidual:number}}
 *   forces [N/m]; friction* = force / W [-]
 */
export function lateralResistance({ model = 'coulomb-passive', soil, D, z, W, mu, muResidual, su, gammaSub = 0, gammaDry, phiDeg,
  bermFactor = 1, bermFactorResidual = 0 }) {
  pos('D', D); nonNeg('z', z); pos('W', W); pos('mu', mu); nonNeg('bermFactor', bermFactor); nonNeg('bermFactorResidual', bermFactorResidual);
  const muR = muResidual === undefined ? mu : pos('muResidual', muResidual);
  let passive;
  if (model === 'verley') {
    passive = passiveResistance({ soil, FC: W, D, zp: z, gammaSoilSub: gammaSub, su, gammaSoilDry: gammaDry }).FR;
  } else if (model === 'coulomb-passive') {
    if (soil === 'clay') passive = 0.5 * gammaSub * z * z + 2 * pos('su', su) * z;
    else if (soil === 'sand') {
      if (!(phiDeg > 0 && phiDeg < 50)) throw new RangeError('phiDeg must be in (0, 50)');
      passive = 0.5 * Math.tan(Math.PI / 4 + (phiDeg * Math.PI) / 360) ** 2 * pos('gammaSub', gammaSub) * z * z;
    } else throw new RangeError('soil must be "sand" or "clay"');
  } else throw new RangeError('model must be "coulomb-passive" or "verley"');
  const breakout = mu * W + bermFactor * passive;
  const residual = muR * W + bermFactorResidual * passive;
  return { breakout, residual, passive, frictionBreakout: breakout / W, frictionResidual: residual / W };
}

/**
 * Piecewise-linear force–displacement (mobilisation) curve: 0 → breakout at xBreakout → residual at
 * xResidual, constant beyond. Returns the knots and an evaluator (odd in displacement).
 *
 * Reference: tri-linear idealisation used in DNV-RP-F114 (2017) §4.3–4.4 and DNV-RP-F110. Mobilisation
 * displacements are site-specific: breakout is commonly reached within a few millimetres axially and
 * 0.1–0.5 D laterally on soft clay, residual after 0.5–3 D of lateral travel — supply project values.
 *
 * @param {object} o
 * @param {number} o.breakout peak resistance [N/m]
 * @param {number} [o.residual=o.breakout] residual resistance [N/m]
 * @param {number} o.xBreakout displacement at breakout [m]
 * @param {number} [o.xResidual=o.xBreakout] displacement at which the residual is reached [m]
 * @returns {{x:number[], f:number[], force:(function(number):number), stiffness:number}}
 *   x [m]; f [N/m]; force(u) [N/m]; stiffness initial secant [N/m²]
 */
export function mobilisationCurve({ breakout, residual, xBreakout, xResidual }) {
  pos('breakout', breakout); pos('xBreakout', xBreakout);
  const res = residual === undefined ? breakout : pos('residual', residual);
  const xr = xResidual === undefined ? xBreakout : xResidual;
  if (xr < xBreakout) throw new RangeError('xResidual must not be smaller than xBreakout');
  const force = (u) => {
    const a = Math.abs(u);
    const s = Math.sign(u);
    if (a <= xBreakout) return (s * breakout * a) / xBreakout;
    if (a >= xr) return s * res;
    return s * (breakout + ((res - breakout) * (a - xBreakout)) / (xr - xBreakout));
  };
  return { x: [0, xBreakout, xr], f: [0, breakout, res], force, stiffness: breakout / xBreakout };
}

/**
 * Low / best / high estimate friction factors from a lognormal distribution fitted through two quantiles.
 * ln X ~ N(μ, σ): σ = (ln HE − ln LE)/(z_high − z_low), μ = ln LE − σ z_low; BE = median = e^μ.
 *
 * Reference: lognormal quantile matching; the LE / BE / HE convention (5 % / 50 % / 95 % by default) follows
 * DNV-RP-F114 (2017) §3.3 and DNV-RP-F110.
 *
 * @param {object} o
 * @param {number} o.LE low estimate [-]
 * @param {number} o.HE high estimate [-]
 * @param {number} [o.pLow=0.05] probability of non-exceedance of LE [-]
 * @param {number} [o.pHigh=0.95] probability of non-exceedance of HE [-]
 * @returns {{mu:number, sigma:number, BE:number, mean:number, std:number, cov:number,
 *   quantile:(function(number):number)}} mu, sigma of ln X; quantile(p) [-]
 */
export function frictionRange({ LE, HE, pLow = 0.05, pHigh = 0.95 }) {
  pos('LE', LE); pos('HE', HE);
  if (!(HE > LE)) throw new RangeError('HE must exceed LE');
  if (!(pLow > 0 && pLow < pHigh && pHigh < 1)) throw new RangeError('require 0 < pLow < pHigh < 1');
  const zq = (p) => Math.SQRT2 * erfinv(2 * p - 1);
  const zl = zq(pLow);
  const zh = zq(pHigh);
  const sigma = (Math.log(HE) - Math.log(LE)) / (zh - zl);
  const mu = Math.log(LE) - sigma * zl;
  const mean = Math.exp(mu + 0.5 * sigma * sigma);
  const std = mean * Math.sqrt(Math.exp(sigma * sigma) - 1);
  return { mu, sigma, BE: Math.exp(mu), mean, std, cov: std / mean, quantile: (p) => Math.exp(mu + sigma * zq(p)) };
}

// DNV-RP-F105 (2006) Tables 7-5 / 7-6: dynamic stiffness factors [N/m^(5/2)] and static vertical stiffness [N/m/m].
const F105_SOIL = {
  'sand-loose': { CV: 10500e3, CL: 9000e3, KVS: 250e3, nu: 0.35 },
  'sand-medium': { CV: 14500e3, CL: 12500e3, KVS: 530e3, nu: 0.35 },
  'sand-dense': { CV: 21000e3, CL: 18000e3, KVS: 1350e3, nu: 0.35 },
  'clay-very-soft': { CV: 600e3, CL: 500e3, KVS: [50e3, 100e3], nu: 0.45 },
  'clay-soft': { CV: 1400e3, CL: 1200e3, KVS: [160e3, 260e3], nu: 0.45 },
  'clay-firm': { CV: 3000e3, CL: 2600e3, KVS: [500e3, 800e3], nu: 0.45 },
  'clay-stiff': { CV: 4500e3, CL: 3900e3, KVS: [1000e3, 1600e3], nu: 0.45 },
  'clay-very-stiff': { CV: 11000e3, CL: 9500e3, KVS: [2000e3, 3000e3], nu: 0.45 },
  'clay-hard': { CV: 12000e3, CL: 10500e3, KVS: [2600e3, 4200e3], nu: 0.45 },
};

/**
 * Soil stiffness per unit length for free-span and dynamic analyses:
 * K_V = C_V/(1 − ν) · (2/3 ρ_s/ρ + 1/3) √D,  K_L = C_L (1 + ν) (2/3 ρ_s/ρ + 1/3) √D (D in metres).
 *
 * Reference: DNV-RP-F105 (2006) §7.4, Eq. 7.11–7.12 with Tables 7-5 / 7-6 (C_V, C_L; static vertical
 * stiffness K_V,S). The tabulated factors were entered from the 2006 edition and should be confirmed against
 * the edition in force; C_V, C_L, K_V,S and ν can all be overridden. Clay values apply to OCR = 1.
 * Validity: 1.2 < ρ_s/ρ < 2.0 for the mass-ratio term.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.massRatio specific mass ρ_s/ρ of the pipe (mass incl. contents, excl. added mass, over displaced mass) [-]
 * @param {string} [o.soilType] one of sand-loose, sand-medium, sand-dense, clay-very-soft, clay-soft, clay-firm,
 *   clay-stiff, clay-very-stiff, clay-hard
 * @param {number} [o.CV] override [N/m^2.5]
 * @param {number} [o.CL] override [N/m^2.5]
 * @param {number} [o.nu] Poisson's ratio of the soil [-]
 * @param {number} [o.KVS] static vertical stiffness override [N/m²]
 * @returns {{KV:number, KL:number, KVS:(number|null), KVSrange:(number[]|null), CV:number, CL:number, nu:number}}
 *   KV, KL dynamic stiffness [N/m²]; KVS static vertical stiffness [N/m²] (mid-range for clays)
 */
export function soilStiffness({ D, massRatio, soilType, CV, CL, nu, KVS }) {
  pos('D', D); pos('massRatio', massRatio);
  const rec = soilType === undefined ? {} : F105_SOIL[soilType];
  if (rec === undefined) throw new RangeError(`soilType must be one of ${Object.keys(F105_SOIL).join(', ')}`);
  const cv = CV === undefined ? rec.CV : CV;
  const cl = CL === undefined ? rec.CL : CL;
  const v = nu === undefined ? rec.nu : nu;
  pos('CV', cv); pos('CL', cl);
  if (!(v >= 0 && v < 0.5)) throw new RangeError('nu must be in [0, 0.5)');
  const m = (2 / 3) * massRatio + 1 / 3;
  let kvs = KVS === undefined ? rec.KVS : KVS;
  let range = null;
  if (Array.isArray(kvs)) { range = kvs; kvs = 0.5 * (kvs[0] + kvs[1]); }
  return { KV: (cv / (1 - v)) * m * Math.sqrt(D), KL: cl * (1 + v) * m * Math.sqrt(D), KVS: kvs === undefined ? null : kvs,
    KVSrange: range, CV: cv, CL: cl, nu: v };
}

/**
 * Uplift resistance of a buried pipe in cohesionless backfill: R = γ' H D (1 + f H/D), H = cover from the
 * top of the pipe to the surface.
 *
 * Reference: Schaminée, Zorn & Schotman (1990), OTC 6486; used in Palmer et al. (1990), OTC 6335.
 * The uplift coefficient f is about 0.1 for loose and 0.5 for dense sand (rock: about 0.6–0.7); argument.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.H cover above the pipe crown [m]
 * @param {number} o.gammaSub submerged unit weight of the backfill [N/m³]
 * @param {number} o.f uplift coefficient [-]
 * @returns {number} uplift resistance [N/m] (excluding the pipe's own submerged weight)
 */
export function upliftResistanceSand({ D, H, gammaSub, f }) {
  pos('D', D); nonNeg('H', H); pos('gammaSub', gammaSub); nonNeg('f', f);
  return gammaSub * H * D * (1 + (f * H) / D);
}

/**
 * Uplift resistance of a buried pipe in clay: the smaller of the shallow vertical-slip mechanism
 * γ' H D + 2 s_u H and the deep local (flow-around) mechanism N_c s_u D.
 *
 * Reference: Schaminée et al. (1990) for the vertical-slip form; local failure factor N_c as in
 * DNV-RP-F110 (2018) App. B (argument: about 5.14 at shallow embedment up to 9–12 for deep flow-around).
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.H cover above the crown [m]
 * @param {number} o.su undrained (remoulded backfill) shear strength [Pa]
 * @param {number} [o.gammaSub=0] [N/m³]
 * @param {number} [o.Nc=9] local failure bearing factor [-]
 * @returns {{R:number, shallow:number, deep:number, governing:'shallow'|'deep'}} [N/m]
 */
export function upliftResistanceClay({ D, H, su, gammaSub = 0, Nc = 9 }) {
  pos('D', D); nonNeg('H', H); pos('su', su); nonNeg('gammaSub', gammaSub); pos('Nc', Nc);
  const shallow = gammaSub * H * D + 2 * su * H;
  const deep = Nc * su * D;
  return { R: Math.min(shallow, deep), shallow, deep, governing: shallow <= deep ? 'shallow' : 'deep' };
}

/**
 * Cover height required to provide a given uplift resistance.
 * Sand: closed-form root of γ' f H² + γ' D H − R = 0. Clay: H = R/(γ' D + 2 s_u), provided R ≤ N_c s_u D.
 *
 * Reference: inversion of upliftResistanceSand / upliftResistanceClay.
 *
 * @param {object} o
 * @param {'sand'|'clay'} o.soil
 * @param {number} o.R required uplift resistance from the soil [N/m]
 * @param {number} o.D [m]
 * @param {number} [o.gammaSub] [N/m³]
 * @param {number} [o.f] uplift coefficient [-] (sand)
 * @param {number} [o.su] [Pa] (clay)
 * @param {number} [o.Nc=9] [-] (clay)
 * @returns {number} required cover above the crown [m] (0 when R ≤ 0)
 */
export function requiredCover({ soil, R, D, gammaSub = 0, f, su, Nc = 9 }) {
  pos('D', D);
  if (!(R > 0)) return 0;
  if (soil === 'sand') {
    pos('gammaSub', gammaSub); nonNeg('f', f);
    if (f === 0) return R / (gammaSub * D);
    return (-D + Math.sqrt(D * D + (4 * f * R) / gammaSub)) / (2 * f);
  }
  if (soil === 'clay') {
    pos('su', su);
    if (R > Nc * su * D) throw new RangeError('R exceeds the deep-failure limit Nc·su·D; cover alone cannot provide it');
    return R / (gammaSub * D + 2 * su);
  }
  throw new RangeError('soil must be "sand" or "clay"');
}

/**
 * Trench backfill or rock-berm cover over a pipe: cross-section, volume, mass and the vertical soil weight
 * directly above the pipe (prism of width D, less the upper half of the pipe).
 * Trapezoid of crest width b at height D + H above the pipe invert level with side slopes 1:n (V:H):
 * A = (b + n h) h − π D²/4, h = D + H.
 *
 * Reference: geometry. Typical rock berm: side slope 1:2.5 to 1:3, bulk density 1500–1800 kg/m³ in air
 * (porosity 0.35–0.45).
 *
 * @param {object} o
 * @param {number} o.D pipe outer diameter [m]
 * @param {number} o.H cover above the crown [m]
 * @param {number} o.crestWidth b [m]
 * @param {number} o.sideSlope n, horizontal run per unit height [-]
 * @param {number} o.length berm length [m]
 * @param {number} o.rhoSolid density of the rock / soil particles [kg/m³]
 * @param {number} [o.porosity=0.4] [-]
 * @param {number} [o.rhoWater=1025] [kg/m³]
 * @param {number} [o.g=9.80665]
 * @returns {{area:number, volume:number, mass:number, baseWidth:number, gammaSub:number, prismWeight:number}}
 *   area [m²]; volume [m³]; mass dry placed mass [kg]; baseWidth [m]; gammaSub submerged bulk unit weight
 *   [N/m³]; prismWeight submerged weight of cover directly above the pipe [N/m]
 */
export function coverQuantities({ D, H, crestWidth, sideSlope, length, rhoSolid, porosity = 0.4, rhoWater = 1025, g = 9.80665 }) {
  pos('D', D); nonNeg('H', H); nonNeg('crestWidth', crestWidth); nonNeg('sideSlope', sideSlope); pos('length', length); pos('rhoSolid', rhoSolid);
  if (!(porosity >= 0 && porosity < 1)) throw new RangeError('porosity must be in [0, 1)');
  const h = D + H;
  const area = (crestWidth + sideSlope * h) * h - (Math.PI * D * D) / 4;
  const volume = area * length;
  const gammaSub = (1 - porosity) * (rhoSolid - rhoWater) * g;
  return { area, volume, mass: volume * (1 - porosity) * rhoSolid, baseWidth: crestWidth + 2 * sideSlope * h, gammaSub,
    prismWeight: gammaSub * (H * D + (D * D) / 2 - (Math.PI * D * D) / 8) };
}
