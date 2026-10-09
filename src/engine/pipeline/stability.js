/**
 * On-bottom stability of pipelines, cables and umbilicals (DNV-RP-F109 format).
 *
 * Active hydrodynamic load model: the peak-load coefficient tables of DNV-RP-F109 (Tables 3-9 / 3-10,
 * C_Y*, C_Z* versus K* and M*) are NOT built in. Peak loads are obtained either
 *   (a) method 'morison' (default): Morison drag / inertia / lift stepped through one design oscillation
 *       and the worst phase taken, as in the simplified method of DNV RP E305 (1988), or
 *   (b) method 'coefficients': caller-supplied C_Y*, C_Z* (e.g. read from the F109 tables).
 * Every result reports which method was used in its `method` field.
 *
 * The functions apply unchanged to cables and umbilicals (use the outer diameter and submerged weight).
 * All quantities SI; angles in radians unless the name ends in `Deg`.
 */

import { pipeWeights } from './pipeprops.js';

const G = 9.80665;

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/**
 * Wave number from the linear dispersion relation ω² = g k tanh(k d) (Newton iteration).
 *
 * Reference: Airy wave theory.
 *
 * @param {number} omega angular frequency [rad/s]
 * @param {number} depth water depth [m]
 * @param {number} [g=9.80665] [m/s²]
 * @returns {number} wave number k [rad/m]
 */
export function waveNumber(omega, depth, g = G) {
  pos('omega', omega); pos('depth', depth);
  const k0 = (omega * omega) / g;
  let k = k0 / Math.sqrt(Math.tanh(k0 * depth));
  for (let i = 0; i < 50; i++) {
    const th = Math.tanh(k * depth);
    const f = g * k * th - omega * omega;
    const df = g * th + g * k * depth * (1 - th * th);
    const dk = f / df;
    k -= dk;
    if (Math.abs(dk) < 1e-14 * k) break;
  }
  return k;
}

/**
 * JONSWAP peak-enhancement factor from the sea-state steepness parameter φ = T_p/√H_s:
 * γ = 5 for φ ≤ 3.6; exp(5.75 − 1.15 φ) for 3.6 < φ < 5; 1 for φ ≥ 5.
 *
 * Reference: DNV-RP-F109 (2021) §3.4.3; DNV-RP-C205 §3.5.5.
 *
 * @param {number} Hs significant wave height [m]
 * @param {number} Tp peak period [s]
 * @returns {number} γ [-]
 */
export function jonswapGamma(Hs, Tp) {
  pos('Hs', Hs); pos('Tp', Tp);
  const phi = Tp / Math.sqrt(Hs);
  if (phi <= 3.6) return 5;
  if (phi >= 5) return 1;
  return Math.exp(5.75 - 1.15 * phi);
}

/**
 * JONSWAP surface-elevation spectrum
 * S(ω) = α g² ω⁻⁵ exp(−1.25 (ω/ω_p)⁻⁴) γ^exp(−0.5((ω − ω_p)/(σ ω_p))²),
 * α = (5/16) H_s² ω_p⁴/g² (1 − 0.287 ln γ), σ = 0.07 (ω ≤ ω_p) or 0.09.
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.5–3.7.
 * Validity: 1 ≤ γ ≤ 5 for the α normalisation.
 *
 * @param {number} omega [rad/s]
 * @param {{Hs:number, Tp:number, gamma:number}} sea Hs [m], Tp [s], γ [-]
 * @param {number} [g=9.80665] [m/s²]
 * @returns {number} S(ω) [m²·s/rad]
 */
export function jonswapSpectrum(omega, { Hs, Tp, gamma }, g = G) {
  if (!(omega > 0)) return 0;
  const wp = (2 * Math.PI) / Tp;
  const alpha = (5 / 16) * ((Hs * Hs * wp ** 4) / (g * g)) * (1 - 0.287 * Math.log(gamma));
  const sigma = omega <= wp ? 0.07 : 0.09;
  const r = Math.exp(-0.5 * ((omega - wp) / (sigma * wp)) ** 2);
  return alpha * g * g * omega ** -5 * Math.exp(-1.25 * (omega / wp) ** -4) * gamma ** r;
}

/**
 * Significant near-bed wave-induced velocity and its mean zero-up-crossing period by spectral transformation:
 * S_UU(ω) = [ω/sinh(k d)]² S_ηη(ω), U_s = 2√M0, T_u = 2π√(M0/M2).
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.8–3.12 (first-order wave theory, velocity at the seabed).
 * Validity: linear waves; the integration covers 0.2 ω_p to 8 ω_p.
 *
 * @param {object} o
 * @param {number} o.Hs significant wave height [m]
 * @param {number} o.Tp peak period [s]
 * @param {number} o.depth water depth [m]
 * @param {number} [o.gamma] peak-enhancement factor [-] (default from jonswapGamma)
 * @param {number} [o.nOmega=800] number of frequency intervals [-]
 * @param {number} [o.g=9.80665] [m/s²]
 * @returns {{Us:number, Tu:number, M0:number, M2:number, gamma:number, Tn:number}}
 *   Us [m/s]; Tu [s]; M0 [m²/s²]; M2 [m²/s⁴]; Tn = √(d/g) reference period [s]
 */
export function nearBedSpectralVelocity({ Hs, Tp, depth, gamma, nOmega = 800, g = G }) {
  pos('Hs', Hs); pos('Tp', Tp); pos('depth', depth);
  const gam = gamma === undefined ? jonswapGamma(Hs, Tp) : pos('gamma', gamma);
  const n = Math.max(50, Math.round(nOmega / 2) * 2);
  const wp = (2 * Math.PI) / Tp;
  const a = 0.2 * wp;
  const h = (8 * wp - a) / n;
  let M0 = 0;
  let M2 = 0;
  for (let i = 0; i <= n; i++) {
    const w = a + i * h;
    const kd = waveNumber(w, depth, g) * depth;
    const Gt = kd > 300 ? 0 : w / Math.sinh(kd);
    const s = Gt * Gt * jonswapSpectrum(w, { Hs, Tp, gamma: gam }, g);
    const c = i === 0 || i === n ? 1 : i % 2 ? 4 : 2;
    M0 += c * s;
    M2 += c * s * w * w;
  }
  M0 *= h / 3; M2 *= h / 3;
  return { Us: 2 * Math.sqrt(M0), Tu: M2 > 0 ? 2 * Math.PI * Math.sqrt(M0 / M2) : Tp, M0, M2, gamma: gam, Tn: Math.sqrt(depth / g) };
}

/**
 * Single design oscillation: U* = ½ U_s (√(2 ln τ) + 0.5772/√(2 ln τ)), τ = T/T_u, and its period
 * T* ÷ T_u = k_T − 5 (k_T − 1) T_n/T_u for T_n/T_u ≤ 0.2, else 1; k_T = 1.25 (γ = 1), 1.21 (γ = 3.3), 1.17 (γ = 5)
 * with linear interpolation in γ.
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.15–3.17.
 * A directionality / spreading reduction R_D (F109 §3.4.4; 1.0 for long-crested waves normal to the pipe,
 * typically 0.8–1.0 for short-crested seas depending on the spreading exponent and heading) multiplies U_s.
 *
 * @param {object} o
 * @param {number} o.Us significant near-bed velocity normal to the pipe before reduction [m/s]
 * @param {number} o.Tu zero-up-crossing period of the near-bed velocity [s]
 * @param {number} [o.durationS=10800] storm duration T [s]
 * @param {number} [o.RD=1] directionality and spreading reduction factor [-]
 * @param {number} [o.gamma=3.3] JONSWAP γ used for k_T [-]
 * @param {number} [o.Tn] reference period √(d/g) [s]; when omitted the period correction is skipped and T* = T_u
 * @param {number} [o.kT] override of k_T [-]
 * @returns {{Ustar:number, Tstar:number, tau:number, kT:number, UsReduced:number}} Ustar [m/s]; Tstar [s]
 */
export function designOscillation({ Us, Tu, durationS = 10800, RD = 1, gamma = 3.3, Tn, kT }) {
  nonNeg('Us', Us); pos('Tu', Tu); pos('durationS', durationS); pos('RD', RD);
  const tau = durationS / Tu;
  if (tau <= 1) throw new RangeError('durationS must exceed Tu');
  const s = Math.sqrt(2 * Math.log(tau));
  const UsR = Us * RD;
  const Ustar = 0.5 * UsR * (s + 0.5772 / s);
  let k = kT;
  if (k === undefined) {
    const gm = Math.min(5, Math.max(1, gamma));
    k = gm <= 3.3 ? 1.25 + ((1.21 - 1.25) * (gm - 1)) / 2.3 : 1.21 + ((1.17 - 1.21) * (gm - 3.3)) / 1.7;
  }
  let Tstar = Tu;
  if (Tn !== undefined) {
    const r = pos('Tn', Tn) / Tu;
    Tstar = Tu * (r <= 0.2 ? k - 5 * (k - 1) * r : 1);
  }
  return { Ustar, Tstar, tau, kT: k, UsReduced: UsR };
}

/**
 * Mean current velocity over the pipe diameter from a reference measurement, for a logarithmic boundary
 * layer: V_c = V(z_r) · [ (1 + z_0/D) ln(D/z_0 + 1) − 1 ] / ln(z_r/z_0 + 1) · sin θ_c.
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.3. Bottom roughness z_0 (Table 3-1): silt and clay ≈ 5e-6 m,
 * fine sand 1e-5, medium sand 4e-5, coarse sand 1e-4, gravel 3e-4, pebble 2e-3, cobble 1e-2, boulder 4e-2 m.
 *
 * @param {object} o
 * @param {number} o.Vref current speed at the reference height [m/s]
 * @param {number} o.zRef reference height above the seabed [m]
 * @param {number} o.z0 bottom roughness parameter [m]
 * @param {number} o.D outer (hydrodynamic) diameter [m]
 * @param {number} [o.angleDeg=90] angle between current and pipe axis [deg]
 * @returns {number} mean perpendicular current velocity over the pipe height [m/s]
 */
export function currentOverPipe({ Vref, zRef, z0, D, angleDeg = 90 }) {
  nonNeg('Vref', Vref); pos('zRef', zRef); pos('z0', z0); pos('D', D);
  return (Vref * ((1 + z0 / D) * Math.log(D / z0 + 1) - 1) * Math.sin((angleDeg * Math.PI) / 180)) / Math.log(zRef / z0 + 1);
}

/**
 * Load reduction factors for penetration, trenching and a permeable seabed.
 * r_pen,y = max(1 − 1.4 z_p/D, 0.3); r_pen,z = max(1 − 1.3 (z_p/D − 0.1), 0);
 * r_tr,y = 1 − 0.18 (θ − 5)^0.25, r_tr,z = 1 − 0.14 (θ − 5)^0.43 for a trench side slope 5° ≤ θ ≤ 45°;
 * r_perm,z = 0.7 on lift for a permeable seabed; r_tot = product.
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.19–3.24. The trench expressions are as recalled from the RP
 * (confirm before relying on trench shielding); they are only applied when trenchAngleDeg is given.
 *
 * @param {object} o
 * @param {number} [o.zp=0] total penetration [m]
 * @param {number} o.D outer diameter [m]
 * @param {number} [o.trenchAngleDeg] trench side-slope angle [deg], 5–45
 * @param {boolean} [o.permeable=false] permeable seabed
 * @returns {{rPenY:number, rPenZ:number, rTrY:number, rTrZ:number, rPermZ:number, rTotY:number, rTotZ:number}} [-]
 */
export function loadReductionFactors({ zp = 0, D, trenchAngleDeg, permeable = false }) {
  nonNeg('zp', zp); pos('D', D);
  const rPenY = Math.max(1 - (1.4 * zp) / D, 0.3);
  const rPenZ = Math.min(1, Math.max(1 - 1.3 * (zp / D - 0.1), 0));
  let rTrY = 1;
  let rTrZ = 1;
  if (trenchAngleDeg !== undefined) {
    if (!(trenchAngleDeg >= 5 && trenchAngleDeg <= 45)) throw new RangeError('trenchAngleDeg must be within 5–45');
    rTrY = 1 - 0.18 * (trenchAngleDeg - 5) ** 0.25;
    rTrZ = 1 - 0.14 * (trenchAngleDeg - 5) ** 0.43;
  }
  const rPermZ = permeable ? 0.7 : 1;
  return { rPenY, rPenZ, rTrY, rTrZ, rPermZ, rTotY: rPenY * rTrY, rTotZ: rPenZ * rTrZ * rPermZ };
}

function soilGroups({ soil, FC, D, gammaSoilSub, su, gammaSoilDry }) {
  if (soil === 'sand') {
    pos('gammaSoilSub', gammaSoilSub);
    return { kappa: (gammaSoilSub * D * D) / FC, Gc: null };
  }
  if (soil === 'clay') {
    pos('su', su); pos('gammaSoilDry', gammaSoilDry);
    return { kappa: (su * D) / FC, Gc: su / (D * gammaSoilDry) };
  }
  throw new RangeError('soil must be "sand" or "clay"');
}

/**
 * Passive soil resistance F_R from pipe penetration.
 * Sand: F_R/F_C = (5.0 κ_s − 0.15 κ_s²)(z_p/D)^1.25 for κ_s ≤ 26.7, κ_s (z_p/D)^1.25 above;
 *       κ_s = γ_s' D²/F_C.
 * Clay: F_R/F_C = 4.1 κ_c/G_c^0.39 (z_p/D)^1.31; κ_c = s_u D/F_C, G_c = s_u/(D γ_s).
 * F_C = w_s − F_Z is the vertical contact force.
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.25–3.29 (Verley & Sotberg 1992 for sand; Verley & Lund 1995 for clay).
 * Validity: calibrated for z_p/D up to about 0.5; F_C > 0 (returns 0 when the pipe lifts off).
 *
 * @param {object} o
 * @param {'sand'|'clay'} o.soil
 * @param {number} o.FC vertical contact force w_s − F_Z [N/m]
 * @param {number} o.D outer diameter [m]
 * @param {number} o.zp total penetration [m]
 * @param {number} [o.gammaSoilSub] submerged unit weight of sand [N/m³]
 * @param {number} [o.su] undrained shear strength of clay [Pa]
 * @param {number} [o.gammaSoilDry] dry unit weight of clay [N/m³] (≈ 18 000 N/m³)
 * @returns {{FR:number, kappa:number, Gc:(number|null)}} FR [N/m]
 */
export function passiveResistance({ soil, FC, D, zp, gammaSoilSub, su, gammaSoilDry }) {
  pos('D', D); nonNeg('zp', zp);
  if (!(FC > 0) || zp === 0) return { FR: 0, kappa: FC > 0 ? soilGroups({ soil, FC, D, gammaSoilSub, su, gammaSoilDry }).kappa : Infinity, Gc: null };
  const { kappa, Gc } = soilGroups({ soil, FC, D, gammaSoilSub, su, gammaSoilDry });
  const r = zp / D;
  if (soil === 'sand') {
    const c = kappa <= 26.7 ? 5.0 * kappa - 0.15 * kappa * kappa : kappa;
    return { FR: FC * c * r ** 1.25, kappa, Gc };
  }
  return { FR: FC * ((4.1 * kappa) / Gc ** 0.39) * r ** 1.31, kappa, Gc };
}

/**
 * Initial penetration of a pipe laid on the seabed.
 * Sand: z_pi/D = 0.037 κ_s^−0.67. Clay: z_pi/D = 0.0071 (G_c^0.3/κ_c)^3.2 + 0.062 (G_c^0.3/κ_c)^0.7.
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.26 and 3.30 (F_Z = 0 at installation unless given).
 *
 * @param {object} o
 * @param {'sand'|'clay'} o.soil
 * @param {number} o.FC vertical contact force [N/m]
 * @param {number} o.D outer diameter [m]
 * @param {number} [o.gammaSoilSub] [N/m³] (sand)
 * @param {number} [o.su] [Pa] (clay)
 * @param {number} [o.gammaSoilDry] [N/m³] (clay)
 * @returns {{zpi:number, kappa:number, Gc:(number|null)}} zpi [m]
 */
export function initialPenetration({ soil, FC, D, gammaSoilSub, su, gammaSoilDry }) {
  pos('FC', FC); pos('D', D);
  const { kappa, Gc } = soilGroups({ soil, FC, D, gammaSoilSub, su, gammaSoilDry });
  if (soil === 'sand') return { zpi: D * 0.037 * kappa ** -0.67, kappa, Gc };
  const x = Gc ** 0.3 / kappa;
  return { zpi: D * (0.0071 * x ** 3.2 + 0.062 * x ** 0.7), kappa, Gc };
}

/**
 * Vertical stability (flotation) in water: γ_w b/(w_s + b) ≤ 1.
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.1 (γ_w = 1.1 unless a sufficiently low probability of negative
 * buoyancy is otherwise documented).
 *
 * @param {object} o
 * @param {number} o.ws submerged weight [N/m]
 * @param {number} o.b buoyancy [N/m]
 * @param {number} [o.gammaW=1.1] safety factor [-]
 * @returns {{uc:number, sg:number, sgRequired:number}} uc [-]; sg = (ws + b)/b; sgRequired = γ_w
 */
export function verticalStability({ ws, b, gammaW = 1.1 }) {
  pos('b', b); pos('gammaW', gammaW);
  if (!(ws + b > 0)) throw new RangeError('ws + b must be positive');
  return { uc: (gammaW * b) / (ws + b), sg: (ws + b) / b, sgRequired: gammaW };
}

/**
 * Morison loads on a pipe on the seabed through one design oscillation.
 * U(θ) = U* cos θ + V*, dU/dt = −(2π/T*) U* sin θ;
 * F_Y = ½ ρ D C_D U|U| + ρ (π D²/4) C_M dU/dt;  F_Z = ½ ρ D C_L U².
 *
 * Reference: Morison et al. (1950); DNV RP E305 (1988) §5.3 simplified static stability, whose recommended
 * coefficients are C_D = 0.7, C_M = 3.29, C_L = 0.9 (to be used with its calibration factor F_w).
 * Validity: the coefficients are effective values for a pipe on the bed, not free-stream values.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.Ustar oscillatory velocity amplitude [m/s]
 * @param {number} o.Vstar steady current over the pipe [m/s]
 * @param {number} o.Tstar oscillation period [s]
 * @param {number} [o.rhoW=1025] [kg/m³]
 * @param {number} [o.CD=0.7] [-]
 * @param {number} [o.CM=3.29] [-]
 * @param {number} [o.CL=0.9] [-]
 * @param {number} [o.nPhase=360] phase steps per cycle [-]
 * @returns {{phaseDeg:Float64Array, FY:Float64Array, FZ:Float64Array, FYmax:number, FZmax:number, K:number, M:number}}
 *   forces [N/m]; K = U* T* /D; M = V* /U*
 */
export function morisonCycle({ D, Ustar, Vstar, Tstar, rhoW = 1025, CD = 0.7, CM = 3.29, CL = 0.9, nPhase = 360 }) {
  pos('D', D); nonNeg('Ustar', Ustar); nonNeg('Vstar', Vstar); pos('Tstar', Tstar); pos('rhoW', rhoW);
  const n = Math.max(8, Math.floor(nPhase));
  const phaseDeg = new Float64Array(n);
  const FY = new Float64Array(n);
  const FZ = new Float64Array(n);
  const om = (2 * Math.PI) / Tstar;
  let FYmax = 0;
  let FZmax = 0;
  for (let i = 0; i < n; i++) {
    const th = (2 * Math.PI * i) / n;
    const U = Ustar * Math.cos(th) + Vstar;
    const A = -om * Ustar * Math.sin(th);
    FY[i] = 0.5 * rhoW * D * CD * U * Math.abs(U) + rhoW * (Math.PI / 4) * D * D * CM * A;
    FZ[i] = 0.5 * rhoW * D * CL * U * U;
    phaseDeg[i] = (360 * i) / n;
    if (Math.abs(FY[i]) > FYmax) FYmax = Math.abs(FY[i]);
    if (FZ[i] > FZmax) FZmax = FZ[i];
  }
  return { phaseDeg, FY, FZ, FYmax, FZmax, K: (Ustar * Tstar) / D, M: Ustar > 0 ? Vstar / Ustar : Infinity };
}

/**
 * Absolute lateral static stability of a pipe on the seabed:
 *   γ_SC (F_Y* + μ F_Z*)/(μ w_s + F_R) ≤ 1   and   γ_SC F_Z* /w_s ≤ 1.
 *
 * method 'coefficients': F_Y* = r_tot,y ½ ρ D C_Y* (U* + V*)², F_Z* = r_tot,z ½ ρ D C_Z* (U* + V*)² with
 *   caller-supplied peak coefficients (DNV-RP-F109 Tables 3-9 / 3-10 are not built in).
 * method 'morison' (default): F_Y(θ), F_Z(θ) from morisonCycle, multiplied by the calibration factor F_w and
 *   by the load-reduction factors; both criteria are evaluated at every phase and the worst is reported.
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.38, 3.39 (criteria); DNV RP E305 (1988) §5.3 (Morison route).
 * Safety factor γ_SC: F109 Tables 3-5 to 3-8 tabulate it by region (North Sea, Gulf of Mexico, Southern
 * Ocean; winter storm / cyclonic) and safety class for sand/rock and clay, spanning roughly 1.0 to 1.5 —
 * the tables are not reproduced here; pass the project value. Friction μ: 0.6 sand and rock, 0.2 clay
 * (F109 §3.4.6). F_w of E305 depends on K and M (its Fig. 5.12); 1.0 leaves the Morison loads uncalibrated.
 * Validity: static criterion (no displacement allowed); passive resistance per passiveResistance.
 *
 * @param {object} o
 * @param {number} o.ws submerged weight [N/m]
 * @param {number} o.D outer diameter [m]
 * @param {number} o.mu Coulomb friction coefficient [-]
 * @param {number} o.Ustar design oscillation velocity amplitude [m/s]
 * @param {number} o.Vstar steady current over the pipe [m/s]
 * @param {number} [o.Tstar] design oscillation period [s] (required for 'morison')
 * @param {number} o.gammaSC safety factor [-] (no default: take it from the F109 table for the region, soil and safety class)
 * @param {'morison'|'coefficients'} [o.method='morison']
 * @param {number} [o.CY] peak horizontal load coefficient C_Y* [-] ('coefficients')
 * @param {number} [o.CZ] peak vertical load coefficient C_Z* [-] ('coefficients')
 * @param {number} [o.Fw=1] calibration factor on the Morison loads [-]
 * @param {number} [o.CD=0.7] [-]
 * @param {number} [o.CM=3.29] [-]
 * @param {number} [o.CL=0.9] [-]
 * @param {number} [o.rhoW=1025] [kg/m³]
 * @param {number} [o.zp=0] total penetration for load reduction and passive resistance [m]
 * @param {number} [o.trenchAngleDeg] trench slope [deg]
 * @param {boolean} [o.permeable=false]
 * @param {'sand'|'clay'|'none'} [o.soil='none'] soil model for passive resistance ('none': F_R = 0)
 * @param {number} [o.gammaSoilSub] [N/m³]
 * @param {number} [o.su] [Pa]
 * @param {number} [o.gammaSoilDry] [N/m³]
 * @param {number} [o.FR] passive resistance override [N/m]
 * @param {number} [o.nPhase=360]
 * @returns {{method:string, ucLateral:number, ucLift:number, uc:number, FY:number, FZ:number, FR:number,
 *   phaseDeg:(number|null), K:number, M:number, reduction:object}}
 *   FY, FZ, FR [N/m] at the governing lateral phase (after reduction and F_w); uc = max(ucLateral, ucLift)
 */
export function absoluteLateralStability(o) {
  const { ws, D, mu, Ustar, Vstar, Tstar, gammaSC, method = 'morison', CY, CZ, Fw = 1, CD = 0.7, CM = 3.29, CL = 0.9,
    rhoW = 1025, zp = 0, trenchAngleDeg, permeable = false, soil = 'none', gammaSoilSub, su, gammaSoilDry, FR, nPhase = 360 } = o;
  pos('ws', ws); pos('D', D); pos('mu', mu); nonNeg('Ustar', Ustar); nonNeg('Vstar', Vstar); pos('gammaSC', gammaSC); pos('Fw', Fw);
  const red = loadReductionFactors({ zp, D, trenchAngleDeg, permeable });
  const frOf = (FZ) => {
    if (FR !== undefined) return nonNeg('FR', FR);
    if (soil === 'none') return 0;
    return passiveResistance({ soil, FC: ws - FZ, D, zp, gammaSoilSub, su, gammaSoilDry }).FR;
  };
  if (method === 'coefficients') {
    pos('CY', CY); nonNeg('CZ', CZ);
    const q = 0.5 * rhoW * D * (Ustar + Vstar) ** 2;
    const FYs = red.rTotY * q * CY;
    const FZs = red.rTotZ * q * CZ;
    const fr = frOf(FZs);
    const ucLateral = (gammaSC * (FYs + mu * FZs)) / (mu * ws + fr);
    const ucLift = (gammaSC * FZs) / ws;
    return { method: 'coefficients (caller-supplied C_Y*, C_Z*)', ucLateral, ucLift, uc: Math.max(ucLateral, ucLift), FY: FYs, FZ: FZs, FR: fr,
      phaseDeg: null, K: Tstar ? (Ustar * Tstar) / D : NaN, M: Ustar > 0 ? Vstar / Ustar : Infinity, reduction: red };
  }
  if (method !== 'morison') throw new RangeError('method must be "morison" or "coefficients"');
  pos('Tstar', Tstar);
  const cyc = morisonCycle({ D, Ustar, Vstar, Tstar, rhoW, CD, CM, CL, nPhase });
  let best = { ucLateral: -Infinity };
  let ucLift = 0;
  for (let i = 0; i < cyc.FY.length; i++) {
    const fy = Math.abs(cyc.FY[i]) * Fw * red.rTotY;
    const fz = cyc.FZ[i] * Fw * red.rTotZ;
    const fr = frOf(fz);
    const u = (gammaSC * (fy + mu * fz)) / (mu * ws + fr);
    if (u > best.ucLateral) best = { ucLateral: u, FY: fy, FZ: fz, FR: fr, phaseDeg: cyc.phaseDeg[i] };
    const ul = (gammaSC * fz) / ws;
    if (ul > ucLift) ucLift = ul;
  }
  return { method: 'morison (worst phase, E305-type coefficients)', ...best, ucLift, uc: Math.max(best.ucLateral, ucLift),
    K: cyc.K, M: cyc.M, reduction: red };
}

/**
 * Submerged weight required for a target utilisation of the absolute lateral static criterion, by bisection
 * (the utilisation decreases monotonically with weight for fixed diameter).
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.38, 3.39 inverted numerically.
 *
 * @param {object} o arguments of absoluteLateralStability without `ws`
 * @param {number} [o.targetUc=1] target utilisation [-]
 * @returns {{ws:number, check:object}} required submerged weight [N/m] and the check at that weight
 */
export function requiredSubmergedWeight(o) {
  const target = o.targetUc === undefined ? 1 : pos('targetUc', o.targetUc);
  const f = (ws) => absoluteLateralStability({ ...o, ws }).uc - target;
  let lo = 1e-3;
  let hi = 100;
  let k = 0;
  while (f(hi) > 0) {
    lo = hi; hi *= 2;
    if (++k > 60) throw new RangeError('no finite submerged weight satisfies the criterion');
  }
  if (f(lo) <= 0) return { ws: lo, check: absoluteLateralStability({ ...o, ws: lo }) };
  for (let i = 0; i < 70; i++) {
    const mid = 0.5 * (lo + hi);
    if (f(mid) > 0) lo = mid; else hi = mid;
  }
  return { ws: hi, check: absoluteLateralStability({ ...o, ws: hi }) };
}

/**
 * Concrete weight-coating thickness required for vertical and absolute lateral stability. The coating changes
 * both the submerged weight and the hydrodynamic diameter, so the weight build-up and the stability check are
 * re-evaluated for each trial thickness (scan for the first passing thickness, then bisection).
 *
 * Reference: DNV-RP-F109 (2021) Eq. 3.1, 3.38, 3.39 with the weight build-up of pipeWeights.
 * Validity: checks the condition whose contents density is supplied in `pipe.rhoContents` (use 0 for the
 * empty installation case, the product density for operation).
 *
 * @param {object} o
 * @param {object} o.pipe arguments of pipeWeights (D, t, coatings without the concrete layer, rhoContents, ...)
 * @param {number} [o.concreteDensity=3040] [kg/m³]
 * @param {number} [o.concreteAbsorption=0] absorbed water per dry mass [-]
 * @param {object} o.stability arguments of absoluteLateralStability without `ws` and `D`
 * @param {number} [o.gammaW=1.1] vertical stability factor [-]
 * @param {number} [o.targetUc=1] [-]
 * @param {number} [o.tMax=0.2] largest thickness considered [m]
 * @param {number} [o.tStep=0.005] scan step [m]
 * @returns {{thickness:number, ws:number, Dtotal:number, sg:number, check:object}} thickness [m]; ws [N/m]
 */
export function requiredConcreteThickness({ pipe, concreteDensity = 3040, concreteAbsorption = 0, stability, gammaW = 1.1,
  targetUc = 1, tMax = 0.2, tStep = 0.005 }) {
  pos('concreteDensity', concreteDensity); pos('tMax', tMax); pos('tStep', tStep);
  const evalAt = (tc) => {
    const w = pipeWeights({ ...pipe, coatings: [...(pipe.coatings || []), { thickness: tc, density: concreteDensity, absorption: concreteAbsorption }] });
    const ws = w.submergedOperating;
    if (!(ws > 0)) return { ok: false, w, ws, check: null };
    const check = absoluteLateralStability({ ...stability, ws, D: w.Dhydro });
    const vert = verticalStability({ ws, b: w.buoyancy, gammaW });
    return { ok: check.uc <= targetUc && vert.uc <= 1, w, ws, check };
  };
  let lo = 0;
  let r = evalAt(0);
  if (!r.ok) {
    let hi = null;
    for (let tc = tStep; tc <= tMax + 1e-12; tc += tStep) {
      if (evalAt(tc).ok) { hi = tc; break; }
      lo = tc;
    }
    if (hi === null) throw new RangeError('no concrete thickness up to tMax satisfies the stability criteria');
    for (let i = 0; i < 40; i++) {
      const mid = 0.5 * (lo + hi);
      if (evalAt(mid).ok) hi = mid; else lo = mid;
    }
    lo = hi;
    r = evalAt(hi);
  }
  return { thickness: lo, ws: r.ws, Dtotal: r.w.Dtotal, sg: r.w.sgOperating, check: r.check };
}
