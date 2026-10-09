/**
 * Short-circuit calculations for cable systems.
 *
 * References
 *  - IEC 60949:1988 (+A1:2008) Calculation of thermally permissible short-circuit currents, taking into account
 *    non-adiabatic heating effects.
 *  - IEC 60909-0:2016 Short-circuit currents in three-phase AC systems — calculation of currents.
 *  - IEC 60865-1:2011 Short-circuit currents — calculation of effects (electromagnetic forces).
 *  - IEC 60986 / IEC 61443 (short-circuit temperature limits of cables up to 30 kV / above 30 kV).
 *
 * SI units; cross-sections in m² (converted internally to mm² where the standard constants require it).
 */

import { metalProperties } from './conductors.js';
import { cx, cadd, cdiv, cabs, toCx } from './complex.js';

const MU0 = 4e-7 * Math.PI;

function pos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a positive finite number`);
}
function nonNeg(v, name) {
  if (!(v >= 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a non-negative finite number`);
}

function constants(material, K, beta) {
  if (K !== undefined && beta !== undefined) return { K, beta };
  const m = metalProperties(material);
  const k = K ?? m.K; const b = beta ?? m.beta;
  if (!Number.isFinite(k) || !Number.isFinite(b)) throw new RangeError('material: no IEC 60949 constants for this metal; give K and beta');
  return { K: k, beta: b };
}

/**
 * Adiabatic short-circuit current of a current-carrying component:
 * I_AD²·t = K²·S²·ln((θf + β)/(θi + β)), with K in A·s^0.5/mm² and S in mm².
 * Constants (IEC 60949 Table I): copper K = 226, β = 234.5; aluminium 148, 228; lead 41, 230; steel 78, 202;
 * bronze 180, 313.
 * Reference: IEC 60949 §3.
 * @param {object} o
 * @param {string} [o.material] 'Cu' | 'Al' | 'lead' | 'steel' | 'bronze'
 * @param {number} o.area cross-section of the component, m²
 * @param {number} o.thetaInitialC temperature at the start of the short circuit, °C
 * @param {number} o.thetaFinalC permissible final temperature, °C (e.g. 250 XLPE/EPR conductor, 160 PVC, 170–200
 *        lead sheath depending on alloy and covering)
 * @param {number} o.duration short-circuit duration, s
 * @param {number} [o.K] @param {number} [o.beta] overrides of the material constants
 * @returns {{I:number, I2t:number, currentDensity:number}} adiabatic current (A, rms), I²t (A²·s) and current
 *          density (A/m²)
 */
export function adiabaticShortCircuit({ material, area, thetaInitialC, thetaFinalC, duration, K, beta }) {
  pos(area, 'area'); pos(duration, 'duration');
  const c = constants(material, K, beta);
  if (!(thetaFinalC > thetaInitialC)) throw new RangeError('thetaFinalC: must exceed thetaInitialC');
  const S = area * 1e6;
  const I2t = c.K * c.K * S * S * Math.log((thetaFinalC + c.beta) / (thetaInitialC + c.beta));
  const I = Math.sqrt(I2t / duration);
  return { I, I2t, currentDensity: I / area };
}

/**
 * Minimum cross-section for a given short-circuit duty (adiabatic): S = I·sqrt(t)/(K·sqrt(ln((θf + β)/(θi + β)))).
 * Reference: IEC 60949 §3.
 * @param {object} o
 * @param {string} [o.material] @param {number} o.I rms short-circuit current, A @param {number} o.duration s
 * @param {number} o.thetaInitialC °C @param {number} o.thetaFinalC °C @param {number} [o.K] @param {number} [o.beta]
 * @returns {number} cross-section, m²
 */
export function adiabaticMinimumArea({ material, I, duration, thetaInitialC, thetaFinalC, K, beta }) {
  pos(I, 'I'); pos(duration, 'duration');
  const c = constants(material, K, beta);
  if (!(thetaFinalC > thetaInitialC)) throw new RangeError('thetaFinalC: must exceed thetaInitialC');
  return ((I * Math.sqrt(duration)) / (c.K * Math.sqrt(Math.log((thetaFinalC + c.beta) / (thetaInitialC + c.beta))))) * 1e-6;
}

/**
 * Non-adiabatic factor ε for conductors (and spaced wire screens treated as conductors):
 * ε = sqrt(1 + X·sqrt(t/S) + Y·(t/S)), t in s, S in mm²; permissible current I = ε·I_AD.
 * X ((mm²/s)^0.5) and Y (mm²/s) depend on conductor metal and insulation (IEC 60949 Table II). The values for
 * XLPE insulation are provided (copper X = 0.41, Y = 0.12; aluminium X = 0.57, Y = 0.16); for other insulations
 * supply X and Y from the standard. ε is close to 1 for large conductors and short durations.
 * Reference: IEC 60949 §5.
 * @param {object} o
 * @param {number} o.area conductor cross-section, m²
 * @param {number} o.duration s
 * @param {'Cu'|'Al'} [o.material='Cu'] @param {'xlpe'} [o.insulation='xlpe']
 * @param {number} [o.X] @param {number} [o.Y] overrides
 * @returns {number} ε (≥ 1)
 */
export function nonAdiabaticFactorConductor({ area, duration, material = 'Cu', insulation = 'xlpe', X, Y }) {
  pos(area, 'area'); pos(duration, 'duration');
  let x = X; let y = Y;
  if (x === undefined || y === undefined) {
    if (insulation !== 'xlpe') throw new RangeError('insulation: only XLPE constants are built in; give X and Y');
    const t = { Cu: [0.41, 0.12], Al: [0.57, 0.16] }[material];
    if (!t) throw new RangeError('material: must be Cu or Al (or give X and Y)');
    x = x ?? t[0]; y = y ?? t[1];
  }
  const r = duration / (area * 1e6);
  return Math.sqrt(1 + x * Math.sqrt(r) + y * r);
}

/**
 * Non-adiabatic factor ε for sheaths, screens and armour wires:
 * ε = 1 + 0.61·M·sqrt(t) − 0.069·(M·sqrt(t))² + 0.0043·(M·sqrt(t))³,
 * M = F·(sqrt(σ2/ρ2) + sqrt(σ3/ρ3))/(2·σ1·δ·10⁻³), with δ in mm (converted internally from metres),
 * σ1 the volumetric heat capacity of the component, σ2, σ3 and ρ2, ρ3 the heat capacities (J/(K·m³)) and thermal
 * resistivities (K·m/W) of the media on either side, and F the thermal-contact factor (0.7 normally, 0.9 for
 * components completely bonded on one side).
 * Reference: IEC 60949 §6.
 * @param {object} o
 * @param {number} o.thickness thickness of the sheath/screen/armour δ, m
 * @param {number} o.duration s
 * @param {number} o.sigma1 volumetric heat capacity of the component, J/(K·m³)
 * @param {number} o.sigma2 @param {number} o.rho2 properties of the inner adjacent medium
 * @param {number} o.sigma3 @param {number} o.rho3 properties of the outer adjacent medium
 * @param {number} [o.F=0.7] thermal contact factor
 * @returns {{epsilon:number, M:number}} ε and M (s^-0.5)
 */
export function nonAdiabaticFactorSheath({ thickness, duration, sigma1, sigma2, rho2, sigma3, rho3, F = 0.7 }) {
  pos(thickness, 'thickness'); pos(duration, 'duration'); pos(sigma1, 'sigma1'); pos(sigma2, 'sigma2'); pos(rho2, 'rho2');
  pos(sigma3, 'sigma3'); pos(rho3, 'rho3'); pos(F, 'F');
  // δ·10⁻³ with δ in mm is simply the thickness in metres.
  const M = (F * (Math.sqrt(sigma2 / rho2) + Math.sqrt(sigma3 / rho3))) / (2 * sigma1 * thickness);
  const x = M * Math.sqrt(duration);
  return { epsilon: 1 + 0.61 * x - 0.069 * x * x + 0.0043 * x * x * x, M };
}

/**
 * Withstand check of a conductor, screen, sheath or armour against a fault current.
 * @param {object} o arguments of adiabaticShortCircuit plus
 * @param {number} o.Ifault rms fault current carried by the component, A
 * @param {number} [o.epsilon=1] non-adiabatic factor (from nonAdiabaticFactorConductor / nonAdiabaticFactorSheath)
 * @returns {{Ipermissible:number, utilisation:number, ok:boolean, thetaFinalAdiabaticC:number}} A, Ifault/Ipermissible,
 *          and the adiabatic final temperature actually reached for Ifault (°C)
 */
export function faultWithstand({ Ifault, epsilon = 1, material, area, thetaInitialC, thetaFinalC, duration, K, beta }) {
  nonNeg(Ifault, 'Ifault'); pos(epsilon, 'epsilon');
  const c = constants(material, K, beta);
  const ad = adiabaticShortCircuit({ material, area, thetaInitialC, thetaFinalC, duration, K: c.K, beta: c.beta });
  const Ipermissible = ad.I * epsilon;
  const S = area * 1e6;
  const thetaFinalAdiabaticC = (thetaInitialC + c.beta) * Math.exp((Ifault * Ifault * duration) / (c.K * c.K * S * S)) - c.beta;
  return { Ipermissible, utilisation: Ifault / Ipermissible, ok: Ifault <= Ipermissible, thetaFinalAdiabaticC };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* IEC 60909                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Equivalent impedance of a network feeder from its initial symmetrical short-circuit power:
 * ZQ = c·UnQ²/S″kQ, XQ = ZQ/sqrt(1 + (R/X)²), RQ = (R/X)·XQ. With R/X = 0.1 (feeders above 35 kV when unknown)
 * XQ = 0.995·ZQ.
 * Reference: IEC 60909-0 §6.2.
 * @param {object} o
 * @param {number} o.Un nominal line-to-line voltage at the connection point, V
 * @param {number} o.Sk initial symmetrical short-circuit power S″kQ, VA
 * @param {number} [o.c=1.1] voltage factor
 * @param {number} [o.rOverX=0.1]
 * @returns {{re:number, im:number}} impedance, Ω
 */
export function gridImpedance({ Un, Sk, c = 1.1, rOverX = 0.1 }) {
  pos(Un, 'Un'); pos(Sk, 'Sk'); pos(c, 'c'); nonNeg(rOverX, 'rOverX');
  const Z = (c * Un * Un) / Sk;
  const X = Z / Math.sqrt(1 + rOverX * rOverX);
  return cx(rOverX * X, X);
}

/**
 * Short-circuit impedance of a two-winding transformer with the IEC 60909 correction factor:
 * ZT = ukr·UrT²/SrT, RT = PkrT·UrT²/SrT², XT = sqrt(ZT² − RT²), KT = 0.95·cmax/(1 + 0.6·xT) with xT = XT/(UrT²/SrT);
 * the corrected impedance is KT·(RT + j·XT).
 * Reference: IEC 60909-0 §6.3.
 * @param {object} o
 * @param {number} o.Sr rated power, VA @param {number} o.Ur rated voltage of the side considered, V
 * @param {number} o.uk short-circuit voltage, per unit @param {number} [o.Pk=0] load loss, W
 * @param {number} [o.cmax=1.1] @param {boolean} [o.correction=true]
 * @returns {{Z:{re:number, im:number}, KT:number}} corrected impedance (Ω) and the factor applied
 */
export function transformerImpedance({ Sr, Ur, uk, Pk = 0, cmax = 1.1, correction = true }) {
  pos(Sr, 'Sr'); pos(Ur, 'Ur'); pos(uk, 'uk'); nonNeg(Pk, 'Pk');
  const zb = (Ur * Ur) / Sr;
  const R = (Pk / Sr) * zb;
  const Zm = uk * zb;
  if (R > Zm) throw new RangeError('Pk: load loss implies a resistance larger than the impedance');
  const X = Math.sqrt(Zm * Zm - R * R);
  const KT = correction ? (0.95 * cmax) / (1 + 0.6 * (X / zb)) : 1;
  return { Z: cx(KT * R, KT * X), KT };
}

/**
 * Peak factor κ = 1.02 + 0.98·e^(−3R/X).
 * Reference: IEC 60909-0 §8.1.1.
 * @param {number} rOverX R/X of the short-circuit impedance (≥ 0)
 * @returns {number} κ (1.02 … 2.0)
 */
export function peakFactor(rOverX) {
  nonNeg(rOverX, 'rOverX');
  return 1.02 + 0.98 * Math.exp(-3 * rOverX);
}

/**
 * Factor m for the heat effect of the DC component:
 * m = [e^(4·f·Tk·ln(κ − 1)) − 1]/(2·f·Tk·ln(κ − 1)).
 * Reference: IEC 60909-0 §10 (Figure 18).
 * @param {object} o
 * @param {number} o.kappa peak factor @param {number} o.f frequency, Hz @param {number} o.Tk fault duration, s
 * @returns {number} m
 */
export function heatFactorM({ kappa, f, Tk }) {
  pos(f, 'f'); pos(Tk, 'Tk');
  if (!(kappa > 1 && kappa <= 2)) throw new RangeError('kappa: must be in (1, 2]');
  if (kappa === 2) return 2; // limit ln(κ − 1) → 0: m → 2
  const a = 2 * f * Tk * Math.log(kappa - 1);
  return (Math.exp(2 * a) - 1) / a;
}

/**
 * Three-phase short-circuit currents at a fault location fed through a single equivalent impedance
 * (far-from-generator when `IkOverIkss` = 1):
 *   I″k = c·Un/(√3·|Zk|), ip = κ·√2·I″k, Ith = I″k·sqrt(m + n), Ib = μ·I″k.
 * n (decay of the AC component) is 1 for far-from-generator faults and must otherwise be taken from IEC 60909-0
 * Figure 19 (argument). μ is 1 far from generators; for near-to-generator faults use breakingFactorMu.
 * Reference: IEC 60909-0 §7–10.
 * @param {object} o
 * @param {number} o.Un nominal system voltage (line-to-line), V
 * @param {{re:number, im:number}} o.Zk short-circuit impedance at the fault, Ω
 * @param {number} [o.c=1.1] voltage factor (1.1 maximum for > 1 kV, 1.0 minimum)
 * @param {number} [o.f=50] Hz @param {number} [o.Tk=1] fault duration, s
 * @param {number} [o.n=1] heat factor for the AC component @param {number} [o.mu=1] breaking factor
 * @returns {{Ik:number, ip:number, kappa:number, m:number, Ith:number, Ib:number, Sk:number, rOverX:number}}
 *          I″k (A), peak (A), κ, m, thermal equivalent current (A), symmetrical breaking current (A),
 *          S″k = √3·Un·I″k (VA)
 */
export function threePhaseFault({ Un, Zk, c = 1.1, f = 50, Tk = 1, n = 1, mu = 1 }) {
  pos(Un, 'Un'); pos(c, 'c');
  const z = toCx(Zk);
  const mag = cabs(z);
  if (!(mag > 0)) throw new RangeError('Zk: must be non-zero');
  if (!(z.im > 0)) throw new RangeError('Zk: reactance must be positive');
  const Ik = (c * Un) / (Math.sqrt(3) * mag);
  const rOverX = z.re / z.im;
  const kappa = peakFactor(rOverX);
  const m = heatFactorM({ kappa, f, Tk });
  return { Ik, ip: kappa * Math.SQRT2 * Ik, kappa, m, Ith: Ik * Math.sqrt(m + n), Ib: mu * Ik,
    Sk: Math.sqrt(3) * Un * Ik, rOverX };
}

/**
 * Breaking factor μ for near-to-generator faults as a function of I″kG/IrG and the minimum time delay:
 * μ = 0.84 + 0.26·e^(−0.26·r) (tmin = 0.02 s), 0.71 + 0.51·e^(−0.30·r) (0.05 s), 0.62 + 0.72·e^(−0.32·r) (0.10 s),
 * 0.56 + 0.94·e^(−0.38·r) (≥ 0.25 s); μ = 1 when r ≤ 2. Intermediate delays are interpolated linearly.
 * Reference: IEC 60909-0 §9.1.2, Equation (71).
 * @param {object} o
 * @param {number} o.ratio I″kG/IrG @param {number} o.tmin minimum time delay, s
 * @returns {number} μ (≤ 1)
 */
export function breakingFactorMu({ ratio, tmin }) {
  nonNeg(ratio, 'ratio'); nonNeg(tmin, 'tmin');
  if (ratio <= 2) return 1;
  const ts = [0.02, 0.05, 0.10, 0.25];
  const mu = [0.84 + 0.26 * Math.exp(-0.26 * ratio), 0.71 + 0.51 * Math.exp(-0.30 * ratio),
    0.62 + 0.72 * Math.exp(-0.32 * ratio), 0.56 + 0.94 * Math.exp(-0.38 * ratio)].map((v) => Math.min(1, v));
  if (tmin <= ts[0]) return mu[0];
  if (tmin >= ts[3]) return mu[3];
  for (let i = 0; i < 3; i++) {
    if (tmin <= ts[i + 1]) return mu[i] + ((mu[i + 1] - mu[i]) * (tmin - ts[i])) / (ts[i + 1] - ts[i]);
  }
  return mu[3];
}

/**
 * Initial short-circuit current at a busbar fed by several independent sources (sum of magnitudes, which is the
 * conservative IEC 60909 treatment of non-meshed sources).
 *  - impedance sources (grid behind transformer and cable, synchronous or induction machines):
 *    I″k = c·Un/(√3·|Z|), with Z the total impedance to the fault referred to the fault voltage;
 *    for an induction machine Z = (1/(ILR/IrM))·UrM²/SrM, for a synchronous generator Z = KG·(RG + j·X″d),
 *    KG = (Un/UrG)·cmax/(1 + x″d·sin φrG);
 *  - converter-interfaced sources (full-converter wind turbines, HVDC): a current source limited to
 *    k·In with k ≈ 1.1–1.3 supplied by the manufacturer; it contributes to I″k and Ib but not to the DC component.
 * Reference: IEC 60909-0 §7.1, §6.6–6.9, §11 (2016 edition, power-station units with full-size converters).
 * @param {object} o
 * @param {number} o.Un nominal voltage at the fault, V
 * @param {number} [o.c=1.1]
 * @param {{name?:string, Z?:{re:number, im:number}, I?:number}[]} o.sources each source has either an impedance
 *        Z (Ω) or a limited current I (A, already referred to the fault voltage)
 * @returns {{Ik:number, ip:number, contributions:{name:string, Ik:number, ip:number}[]}} A
 */
export function busFault({ Un, c = 1.1, sources }) {
  pos(Un, 'Un');
  if (!Array.isArray(sources) || sources.length === 0) throw new RangeError('sources: must be a non-empty array');
  let Ik = 0; let ip = 0;
  const contributions = sources.map((s, i) => {
    let ik; let peak;
    if (s.Z !== undefined) {
      const r = threePhaseFault({ Un, Zk: s.Z, c });
      ik = r.Ik; peak = r.ip;
    } else if (s.I !== undefined) {
      nonNeg(s.I, 'sources[].I');
      ik = s.I; peak = Math.SQRT2 * s.I;
    } else throw new RangeError('sources[]: each source needs Z or I');
    Ik += ik; ip += peak;
    return { name: s.name ?? `Source ${i + 1}`, Ik: ik, ip: peak };
  });
  return { Ik, ip, contributions };
}

/**
 * Impedance of an induction machine for short-circuit calculation, ZM = UrM²/(SrM·(ILR/IrM)), split with the
 * stated R/X (IEC 60909-0 §6.10: 0.10 for MV machines ≥ 1 MW per pole pair, 0.15 below, 0.42 for LV groups).
 * @param {object} o
 * @param {number} o.Ur rated voltage, V @param {number} o.Sr rated apparent power, VA
 * @param {number} o.lockedRotorRatio ILR/IrM (typically 5–8) @param {number} [o.rOverX=0.1]
 * @returns {{re:number, im:number}} Ω
 */
export function inductionMachineImpedance({ Ur, Sr, lockedRotorRatio, rOverX = 0.1 }) {
  pos(Ur, 'Ur'); pos(Sr, 'Sr'); pos(lockedRotorRatio, 'lockedRotorRatio'); nonNeg(rOverX, 'rOverX');
  const Z = (Ur * Ur) / (Sr * lockedRotorRatio);
  const X = Z / Math.sqrt(1 + rOverX * rOverX);
  return cx(rOverX * X, X);
}

/**
 * Corrected subtransient impedance of a synchronous generator connected directly to the system:
 * ZGK = KG·(RG + j·X″d), X″d = x″d·UrG²/SrG, KG = (Un/UrG)·cmax/(1 + x″d·sin φrG).
 * Reference: IEC 60909-0 §6.6.1.
 * @param {object} o
 * @param {number} o.Ur rated voltage, V @param {number} o.Sr rated power, VA @param {number} o.xd2 x″d, per unit
 * @param {number} o.Un nominal system voltage, V @param {number} [o.powerFactor=0.85] rated power factor
 * @param {number} [o.rOverX=0.07] RG/X″d (0.05 for UrG > 1 kV and SrG ≥ 100 MVA, 0.07 for smaller HV machines,
 *        0.15 for LV machines) @param {number} [o.cmax=1.1]
 * @returns {{Z:{re:number, im:number}, KG:number}} Ω
 */
export function synchronousGeneratorImpedance({ Ur, Sr, xd2, Un, powerFactor = 0.85, rOverX = 0.07, cmax = 1.1 }) {
  pos(Ur, 'Ur'); pos(Sr, 'Sr'); pos(xd2, 'xd2'); pos(Un, 'Un');
  if (!(powerFactor > 0 && powerFactor <= 1)) throw new RangeError('powerFactor: must be in (0, 1]');
  const X = (xd2 * Ur * Ur) / Sr;
  const KG = ((Un / Ur) * cmax) / (1 + xd2 * Math.sqrt(1 - powerFactor * powerFactor));
  return { Z: cx(KG * rOverX * X, KG * X), KG };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Earth faults and forces                                                                                      */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Division of a single-phase-to-earth fault current between the metallic return (sheath/armour, bonded and
 * earthed at both ends) and the sea/earth, for a cable long compared with its end effects:
 * reduction factor r = Iearth/Ifault = 1 − zm/zs = Rs/(Rs + Re + j·X), where zs = Rs + Re + jωμ0/(2π)·ln(De/rs)
 * is the sheath self-impedance with earth return and zm = Re + jωμ0/(2π)·ln(De/rs) the conductor–sheath mutual
 * impedance (Carson first-term approximation, De = 658.87·sqrt(ρ/f), Re = π²·f·10⁻⁷ Ω/m).
 * A measured or specified reduction factor can be given instead and overrides the calculation.
 * Reference: ITU-T Directives vol. II / IEC 60909-3 (reduction factor of cables with metallic sheath).
 * @param {object} o
 * @param {number} o.Ifault fault current, A
 * @param {number} [o.Rs] resistance of the sheath (∥ armour) return path, Ω/m
 * @param {number} [o.ds] mean diameter of the return path, m
 * @param {number} [o.f=50] Hz @param {number} [o.rhoEarth=1] Ω·m (seawater ≈ 0.25)
 * @param {number} [o.reductionFactor] |r| override (0–1)
 * @returns {{reductionFactor:number, Iearth:number, Isheath:number}} |r|, current returning through sea/earth and
 *          through the metallic path (A, magnitudes; from the impedances they do not add arithmetically because
 *          of the phase shift. With a reduction-factor override the phase is unknown and Isheath = (1 − r)·Ifault
 *          is the in-phase split; for a withstand check of the sheath the full fault current is the safe value.)
 */
export function earthFaultSplit({ Ifault, Rs, ds, f = 50, rhoEarth = 1, reductionFactor }) {
  nonNeg(Ifault, 'Ifault');
  if (reductionFactor !== undefined) {
    if (!(reductionFactor >= 0 && reductionFactor <= 1)) throw new RangeError('reductionFactor: must be in [0, 1]');
    return { reductionFactor, Iearth: reductionFactor * Ifault, Isheath: (1 - reductionFactor) * Ifault };
  }
  pos(Rs, 'Rs'); pos(ds, 'ds'); pos(f, 'f'); pos(rhoEarth, 'rhoEarth');
  const De = 658.87 * Math.sqrt(rhoEarth / f);
  const Re = Math.PI * Math.PI * f * 1e-7;
  const X = 2 * Math.PI * f * (MU0 / (2 * Math.PI)) * Math.log(De / (ds / 2));
  const zs = cx(Rs + Re, X);
  const r = cdiv(cx(Rs, 0), zs);
  const sheath = cdiv(cx(Re, X), zs);
  return { reductionFactor: cabs(r), Iearth: cabs(r) * Ifault, Isheath: cabs(sheath) * Ifault };
}

/**
 * Peak electromagnetic force per unit length between parallel conductors during a short circuit:
 * F′ = k·μ0/(2π)·ip²/a. k = 1 for a line-to-line (two-conductor) fault; k = √3/2 = 0.866 for the central
 * conductor of a flat three-phase arrangement and for each conductor of a trefoil arrangement under a balanced
 * three-phase fault (force referred to the three-phase peak current).
 * Reference: IEC 60865-1 §5.2 (Equations for Fm2 and Fm3).
 * @param {object} o
 * @param {number} o.ip peak short-circuit current, A
 * @param {number} o.spacing centre-line distance between conductors a, m
 * @param {'twoPhase'|'threePhase'} [o.fault='threePhase']
 * @param {number} [o.k] configuration factor override
 * @returns {number} force per unit length, N/m
 */
export function shortCircuitForce({ ip, spacing, fault = 'threePhase', k }) {
  nonNeg(ip, 'ip'); pos(spacing, 'spacing');
  let kk = k;
  if (kk === undefined) {
    if (fault === 'threePhase') kk = Math.sqrt(3) / 2;
    else if (fault === 'twoPhase') kk = 1;
    else throw new RangeError('fault: must be twoPhase or threePhase');
  }
  return (kk * MU0 * ip * ip) / (2 * Math.PI * spacing);
}

/**
 * Series combination of impedances (grid + transformer + cable ...), all referred to the same voltage.
 * @param {({re:number, im:number}|number)[]} list impedances, Ω
 * @returns {{re:number, im:number}} Ω
 */
export function seriesImpedance(list) {
  if (!Array.isArray(list) || list.length === 0) throw new RangeError('list: must be a non-empty array');
  return list.map(toCx).reduce((a, b) => cadd(a, b));
}
