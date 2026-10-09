// Wave run-up and overtopping: Iribarren number, Hunt, EurOtop (2018) mean-value and design formulae.
// SI units. Discharges q are m³/s per metre of crest (multiply by 1000 for l/s/m).

import { gamma as gammaFn } from '../math/special.js';

const G = 9.80665;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function pos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}
function unit(v, name) {
  need(Number.isFinite(v) && v > 0 && v <= 1, `${name} must lie in (0, 1]`);
}

/**
 * Surf-similarity (Iribarren) number ξ = tan α / sqrt(H/L0), L0 = g T²/(2π).
 * With H = Hm0 and T = Tm−1,0 this is the breaker parameter ξm−1,0 of EurOtop.
 * Reference: Battjes (1974) Proc. 14th ICCE; Iribarren & Nogales (1949).
 * @param {object} o
 * @param {number} o.tanAlpha structure or beach slope tan α [-]
 * @param {number} o.H wave height at the toe [m]
 * @param {number} o.T wave period [s]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{xi:number, L0:number, steepness:number, breakerType:string}} ξ [-], L0 [m], H/L0 [-],
 *   breaker type by Battjes (1974) for breaking on the slope: spilling ξ < 0.4, plunging 0.4–2.0,
 *   collapsing/surging > 2.0
 */
export function iribarren({ tanAlpha, H, T, g = G }) {
  pos(tanAlpha, 'tanAlpha');
  pos(H, 'H');
  pos(T, 'T');
  const L0 = (g * T * T) / (2 * Math.PI);
  const steepness = H / L0;
  const xi = tanAlpha / Math.sqrt(steepness);
  const breakerType = xi < 0.4 ? 'spilling' : xi < 2.0 ? 'plunging' : 'surging';
  return { xi, L0, steepness, breakerType };
}

/**
 * Hunt (1959) run-up of regular breaking waves on a smooth plane slope: R/H = ξ.
 * Reference: Hunt (1959), J. Waterways Harbors Div. ASCE 85(WW3).
 * Validity: breaking waves, ξ below about 2.3 (the result is flagged otherwise).
 * @param {object} o
 * @param {number} o.H wave height [m]
 * @param {number} o.T period [s]
 * @param {number} o.tanAlpha slope [-]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{R:number, xi:number, valid:boolean}} run-up above still-water level [m], ξ [-], validity flag
 */
export function huntRunup({ H, T, tanAlpha, g = G }) {
  const { xi } = iribarren({ tanAlpha, H, T, g });
  return { R: H * xi, xi, valid: xi < 2.3 };
}

/**
 * Roughness influence factors γf for run-up and overtopping.
 * Reference: EurOtop (2018) chapter 6, table of roughness factors for armoured slopes (1:1.5, permeable
 * core unless noted).
 * Units: [-].
 */
export const ROUGHNESS_FACTORS = Object.freeze({
  smooth: 1.0,
  'rock-1layer-impermeable': 0.6,
  'rock-1layer-permeable': 0.45,
  'rock-2layer-impermeable': 0.55,
  'rock-2layer-permeable': 0.4,
  'cube-1layer-flat': 0.49,
  'cube-2layer-random': 0.47,
  antifer: 0.5,
  haro: 0.47,
  tetrapod: 0.38,
  dolos: 0.43,
  'accropode-I': 0.46,
  xbloc: 0.44,
  'core-loc': 0.44,
  'accropode-II': 0.44,
  'cubipod-1layer': 0.49,
  'cubipod-2layer': 0.47,
});

/**
 * Oblique-wave influence factor γβ.
 *   run-up on slopes:            γβ = 1 − 0.0022 |β|
 *   overtopping on slopes:       γβ = 1 − 0.0033 |β|
 *   rubble-mound structures:     γβ = 1 − 0.0063 |β|
 * each for 0° <= |β| <= 80°; for |β| > 80° the value at 80° applies.
 * Reference: EurOtop (2018) chapter 5 (oblique wave attack on dikes) and chapter 6 (armoured slopes).
 * Validity: short-crested seas.
 * @param {object} o
 * @param {number} o.betaDeg angle between wave direction and structure normal [deg]
 * @param {'runup'|'overtopping'|'rubble'} [o.kind] default 'overtopping'
 * @returns {number} γβ [-]
 */
export function obliqueFactor({ betaDeg, kind = 'overtopping' }) {
  need(Number.isFinite(betaDeg), 'betaDeg must be finite');
  const c = { runup: 0.0022, overtopping: 0.0033, rubble: 0.0063 }[kind];
  need(c !== undefined, 'kind must be runup, overtopping or rubble');
  return 1 - c * Math.min(Math.abs(betaDeg), 80);
}

/**
 * Berm influence factor γb = 1 − rB (1 − rdb), rB = B/Lberm,
 * rdb = 0.5 − 0.5 cos(π db/Ru2%) for a berm above still-water level (db < 0 by EurOtop's sign convention),
 * rdb = 0.5 − 0.5 cos(π db/(2 Hm0)) for a berm below it; rdb = 1 outside the influence zone;
 * 0.6 <= γb <= 1.
 * Reference: EurOtop (2018) chapter 5 (influence of berms).
 * @param {object} o
 * @param {number} o.B berm width [m]
 * @param {number} o.Lberm horizontal length of the slope between 1.0 Hm0 below and above the berm [m]
 * @param {number} o.db depth of the berm below still-water level (negative when above) [m]
 * @param {number} o.Hm0 spectral significant wave height at the toe [m]
 * @param {number} o.Ru2 2 % run-up level used for a berm above SWL [m]
 * @returns {number} γb [-]
 */
export function bermFactor({ B, Lberm, db, Hm0, Ru2 }) {
  need(Number.isFinite(B) && B >= 0, 'B must be >= 0');
  pos(Lberm, 'Lberm');
  pos(Hm0, 'Hm0');
  need(Number.isFinite(db), 'db must be finite');
  let rdb;
  if (db >= 0) rdb = db >= 2 * Hm0 ? 1 : 0.5 - 0.5 * Math.cos((Math.PI * db) / (2 * Hm0));
  else {
    pos(Ru2, 'Ru2');
    rdb = -db >= Ru2 ? 1 : 0.5 - 0.5 * Math.cos((Math.PI * db) / Ru2);
  }
  const rB = Math.min(B / Lberm, 1);
  return Math.min(1, Math.max(0.6, 1 - rB * (1 - rdb)));
}

/**
 * Roughness factor for surging waves on rough slopes: γf,surging = γf + (ξ − 1.8)(1 − γf)/8.2 for
 * 1.8 < ξ < 10, γf below and 1.0 above.
 * Reference: EurOtop (2018) §6.2 (run-up on armoured slopes).
 * @param {object} o
 * @param {number} o.gammaF roughness factor [-]
 * @param {number} o.xi breaker parameter ξm−1,0 (multiplied by γb where a berm exists) [-]
 * @returns {number} γf,surging [-]
 */
export function roughnessFactorSurging({ gammaF, xi }) {
  unit(gammaF, 'gammaF');
  pos(xi, 'xi');
  if (xi <= 1.8) return gammaF;
  if (xi >= 10) return 1;
  return gammaF + ((xi - 1.8) * (1 - gammaF)) / 8.2;
}

/**
 * 2 % wave run-up level on a sloping structure, EurOtop (2018):
 *   mean value (Eq. 5.1): Ru2%/Hm0 = 1.65 γb γf γβ ξ, maximum 1.0 γf γβ (4 − 1.5/√(γb ξ))
 *   design/assessment (Eq. 5.2): 1.75 γb γf γβ ξ, maximum 1.07 γf γβ (4 − 1.5/√(γb ξ))
 * with ξ = ξm−1,0. The roughness factor in the maximum is the surging value (roughnessFactorSurging) when
 * `surgingRoughness` is true.
 * Reference: EurOtop (2018) §5.2, §6.2.
 * Validity: 0.5 < γb ξ < 8–10, relatively gentle slopes with breaker-type run-up; an optional cap
 * `maxRatio` (Ru2%/Hm0) can be supplied for permeable-core structures.
 * @param {object} o
 * @param {number} o.Hm0 spectral significant wave height at the toe [m]
 * @param {number} o.Tm10 spectral period Tm−1,0 [s]
 * @param {number} o.tanAlpha slope [-]
 * @param {number} [o.gammaB] berm factor [-]
 * @param {number} [o.gammaF] roughness factor [-]
 * @param {number} [o.gammaBeta] obliquity factor [-]
 * @param {boolean} [o.design] use the design/assessment coefficients (mean + 1σ)
 * @param {boolean} [o.surgingRoughness] apply γf,surging in the surging maximum (default true)
 * @param {number} [o.maxRatio] optional upper limit of Ru2%/Hm0 [-]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Ru2:number, ratio:number, xi:number, regime:string}} Ru2% above SWL [m], Ru2%/Hm0 [-],
 *   ξm−1,0 [-], 'breaking' or 'surging'
 */
export function runupEurOtop({
  Hm0, Tm10, tanAlpha, gammaB = 1, gammaF = 1, gammaBeta = 1, design = false, surgingRoughness = true,
  maxRatio, g = G,
}) {
  unit(gammaB, 'gammaB');
  unit(gammaF, 'gammaF');
  unit(gammaBeta, 'gammaBeta');
  const { xi } = iribarren({ tanAlpha, H: Hm0, T: Tm10, g });
  const c1 = design ? 1.75 : 1.65;
  const c2 = design ? 1.07 : 1.0;
  const brk = c1 * gammaB * gammaF * gammaBeta * xi;
  const gfs = surgingRoughness ? roughnessFactorSurging({ gammaF, xi: gammaB * xi }) : gammaF;
  const mx = Math.max(0, c2 * gfs * gammaBeta * (4 - 1.5 / Math.sqrt(gammaB * xi)));
  let ratio = Math.min(brk, mx);
  if (maxRatio !== undefined) ratio = Math.min(ratio, maxRatio);
  return { Ru2: ratio * Hm0, ratio, xi, regime: brk <= mx ? 'breaking' : 'surging' };
}

/**
 * Mean wave overtopping discharge on a sloping structure (dike, embankment seawall), EurOtop (2018):
 *   Eq. 5.10: q/√(g Hm0³) = 0.023/√tan α · γb ξ · exp[−(2.7 Rc/(ξ Hm0 γb γf γβ γv))^1.3]
 *   Eq. 5.11 (maximum): q/√(g Hm0³) = 0.09 exp[−(1.5 Rc/(Hm0 γf γβ γ*))^1.3]
 *   design/assessment (Eq. 5.12, 5.13): 0.026 with 2.5, and 0.1035 with 1.35.
 * Reference: EurOtop (2018) §5.3.1.
 * Validity: Rc >= 0, slopes of about 1:1.3 to 1:8 (ξm−1,0 < 5); steeper slopes: see overtoppingSteepSlope.
 * @param {object} o
 * @param {number} o.Hm0 spectral significant wave height at the toe [m]
 * @param {number} o.Tm10 spectral period Tm−1,0 [s]
 * @param {number} o.tanAlpha slope [-]
 * @param {number} o.Rc crest freeboard above still-water level [m]
 * @param {number} [o.gammaB] berm factor [-]
 * @param {number} [o.gammaF] roughness factor [-]
 * @param {number} [o.gammaBeta] obliquity factor [-]
 * @param {number} [o.gammaV] factor for a vertical wall on the slope [-]
 * @param {number} [o.gammaStar] combined factor γ* for a storm wall/promenade in the maximum [-]
 * @param {boolean} [o.design] use design/assessment coefficients
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{q:number, qBreaking:number, qMax:number, qDimensionless:number, xi:number, regime:string}}
 *   q [m³/s/m], both branches [m³/s/m], q/√(g Hm0³) [-], ξm−1,0 [-], governing branch
 */
export function overtoppingSlope({
  Hm0, Tm10, tanAlpha, Rc, gammaB = 1, gammaF = 1, gammaBeta = 1, gammaV = 1, gammaStar = 1, design = false,
  g = G,
}) {
  need(Number.isFinite(Rc) && Rc >= 0, 'Rc must be >= 0');
  unit(gammaB, 'gammaB');
  unit(gammaF, 'gammaF');
  unit(gammaBeta, 'gammaBeta');
  unit(gammaV, 'gammaV');
  pos(gammaStar, 'gammaStar');
  const { xi } = iribarren({ tanAlpha, H: Hm0, T: Tm10, g });
  const a1 = design ? 0.026 : 0.023;
  const b1 = design ? 2.5 : 2.7;
  const a2 = design ? 0.1035 : 0.09;
  const b2 = design ? 1.35 : 1.5;
  const scale = Math.sqrt(g * Hm0 * Hm0 * Hm0);
  const qb = (a1 / Math.sqrt(tanAlpha)) * gammaB * xi
    * Math.exp(-Math.pow((b1 * Rc) / (xi * Hm0 * gammaB * gammaF * gammaBeta * gammaV), 1.3));
  const qm = a2 * Math.exp(-Math.pow((b2 * Rc) / (Hm0 * gammaF * gammaBeta * gammaStar), 1.3));
  const qd = Math.min(qb, qm);
  return { q: qd * scale, qBreaking: qb * scale, qMax: qm * scale, qDimensionless: qd, xi, regime: qb <= qm ? 'breaking' : 'non-breaking' };
}

/**
 * Mean overtopping discharge for very steep slopes up to vertical walls (cot α < 2), EurOtop (2018):
 * q/√(g Hm0³) = a exp[−(b Rc/(Hm0 γf γβ))^1.3],
 *   a = 0.09 − 0.01 (2 − cot α)^2.1,  b = 1.5 + 0.42 (2 − cot α)^1.5, with a maximum of b = 2.35.
 * Reference: EurOtop (2018) chapter 5 (very steep slopes up to vertical walls), mean-value approach.
 * Validity: non-breaking waves, 0 <= cot α <= 2.
 * @param {object} o
 * @param {number} o.Hm0 [m]
 * @param {number} o.cotAlpha slope cotangent [-]
 * @param {number} o.Rc crest freeboard [m]
 * @param {number} [o.gammaF] [-]
 * @param {number} [o.gammaBeta] [-]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{q:number, a:number, b:number}} q [m³/s/m], coefficients [-]
 */
export function overtoppingSteepSlope({ Hm0, cotAlpha, Rc, gammaF = 1, gammaBeta = 1, g = G }) {
  pos(Hm0, 'Hm0');
  need(Number.isFinite(cotAlpha) && cotAlpha >= 0 && cotAlpha <= 2, 'cotAlpha must lie in [0, 2]');
  need(Number.isFinite(Rc) && Rc >= 0, 'Rc must be >= 0');
  unit(gammaF, 'gammaF');
  unit(gammaBeta, 'gammaBeta');
  const a = 0.09 - 0.01 * Math.pow(2 - cotAlpha, 2.1);
  const b = Math.min(2.35, 1.5 + 0.42 * Math.pow(2 - cotAlpha, 1.5));
  const q = a * Math.exp(-Math.pow((b * Rc) / (Hm0 * gammaF * gammaBeta), 1.3)) * Math.sqrt(g * Hm0 ** 3);
  return { q, a, b };
}

/**
 * Mean overtopping discharge at a plain vertical wall, EurOtop (2018) chapter 7.
 *   No influencing foreshore:
 *     q/√(g Hm0³) = 0.047 exp[−(2.35 Rc/Hm0)^1.3]   (design: 0.054 and 2.12)
 *   With foreshore, non-impulsive, h²/(Hm0 Lm−1,0) > 0.23:
 *     q/√(g Hm0³) = 0.05 exp(−2.78 Rc/Hm0)          (design: 0.062 and 2.61)
 *   With foreshore, impulsive, h²/(Hm0 Lm−1,0) <= 0.23:
 *     q/√(g Hm0³) = 0.011 (Hm0/(h sm−1,0))^0.5 exp(−2.2 Rc/Hm0)   for 0 < Rc/Hm0 < 1.35
 *     q/√(g Hm0³) = 0.0014 (Hm0/(h sm−1,0))^0.5 (Rc/Hm0)^−3       for Rc/Hm0 >= 1.35
 *     (design: 0.0155 and 0.0020)
 * Lm−1,0 = g Tm−1,0²/(2π) and sm−1,0 = Hm0/Lm−1,0.
 * Reference: EurOtop (2018) chapter 7 (plain vertical walls, Eq. 7.1–7.10).
 * Validity: plain vertical walls without a mound; Rc > 0; Hm0 and Tm−1,0 at the toe.
 * @param {object} o
 * @param {number} o.Hm0 [m]
 * @param {number} o.Rc crest freeboard [m]
 * @param {number} [o.Tm10] spectral period [s] (required with foreshore)
 * @param {number} [o.h] water depth at the toe of the wall [m] (required with foreshore)
 * @param {boolean} [o.foreshore] true when a foreshore influences the waves (default true when h is given)
 * @param {boolean} [o.design] use design/assessment coefficients
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{q:number, qDimensionless:number, regime:string, impulsiveness:number}} q [m³/s/m],
 *   q/√(g Hm0³) [-], 'no-foreshore' | 'non-impulsive' | 'impulsive', h²/(Hm0 Lm−1,0) [-] (NaN without foreshore)
 */
export function overtoppingVertical({ Hm0, Rc, Tm10, h, foreshore, design = false, g = G }) {
  pos(Hm0, 'Hm0');
  need(Number.isFinite(Rc) && Rc >= 0, 'Rc must be >= 0');
  const fs = foreshore !== undefined ? foreshore : h !== undefined;
  const scale = Math.sqrt(g * Hm0 ** 3);
  const r = Rc / Hm0;
  if (!fs) {
    const qd = design ? 0.054 * Math.exp(-Math.pow(2.12 * r, 1.3)) : 0.047 * Math.exp(-Math.pow(2.35 * r, 1.3));
    return { q: qd * scale, qDimensionless: qd, regime: 'no-foreshore', impulsiveness: Number.NaN };
  }
  pos(Tm10, 'Tm10');
  pos(h, 'h');
  const L = (g * Tm10 * Tm10) / (2 * Math.PI);
  const imp = (h * h) / (Hm0 * L);
  if (imp > 0.23) {
    const qd = design ? 0.062 * Math.exp(-2.61 * r) : 0.05 * Math.exp(-2.78 * r);
    return { q: qd * scale, qDimensionless: qd, regime: 'non-impulsive', impulsiveness: imp };
  }
  need(Rc > 0, 'Rc must be > 0 in the impulsive regime');
  const f = Math.sqrt(Hm0 / (h * (Hm0 / L)));
  let qd;
  if (r < 1.35) qd = (design ? 0.0155 : 0.011) * f * Math.exp(-2.2 * r);
  else qd = ((design ? 0.002 : 0.0014) * f) / (r * r * r);
  return { q: qd * scale, qDimensionless: qd, regime: 'impulsive', impulsiveness: imp };
}

/**
 * Mean overtopping discharge over a rubble-mound structure: the EurOtop slope formulae (Eq. 5.10–5.13) with the armour roughness γf and the rubble-mound obliquity factor, multiplied by the crest
 * reduction Cr = 3.06 exp(−1.5 Gc/Hm0) (maximum Cr = 1) for an armoured crest berm of width Gc.
 * Reference: EurOtop (2018) §6.3 (crest reduction after Besley 1999).
 * Validity: steep armoured slopes 1:1.3 to 1:4, Rc measured to the top of the armour crest.
 * @param {object} o
 * @param {number} o.Hm0 [m]
 * @param {number} o.Tm10 [s]
 * @param {number} o.tanAlpha armour slope [-]
 * @param {number} o.Rc crest freeboard [m]
 * @param {number} o.gammaF armour roughness factor (see ROUGHNESS_FACTORS) [-]
 * @param {number} [o.gammaBeta] obliquity factor [-]
 * @param {number} [o.Gc] armoured crest width [m]
 * @param {boolean} [o.design]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{q:number, qNoCrest:number, Cr:number, xi:number, regime:string}} q [m³/s/m], q before crest
 *   reduction [m³/s/m], Cr [-], ξm−1,0 [-], governing branch
 */
export function overtoppingRubbleMound({ Hm0, Tm10, tanAlpha, Rc, gammaF, gammaBeta = 1, Gc = 0, design = false, g = G }) {
  need(Number.isFinite(Gc) && Gc >= 0, 'Gc must be >= 0');
  const r = overtoppingSlope({ Hm0, Tm10, tanAlpha, Rc, gammaF, gammaBeta, design, g });
  const Cr = Math.min(1, 3.06 * Math.exp((-1.5 * Gc) / Hm0));
  return { q: r.q * Cr, qNoCrest: r.q, Cr, xi: r.xi, regime: r.regime };
}

/**
 * Tolerable mean overtopping discharges and maximum individual volumes.
 * Each entry: { id, group, description (states the wave height Hm0 the limit refers to),
 * q (limit of mean discharge) [m³/s/m], Vmax (limit of maximum individual volume) [m³/m] }.
 * Reference: EurOtop (2018) chapter 3, tables of limits for structural design, for property behind the
 * defence and for people and vehicles (lower values where the manual gives a range). The limits depend strongly on wave height; they are guidance, not criteria,
 * and can be overridden in `classifyOvertopping`.
 */
export const TOLERABLE_OVERTOPPING = Object.freeze([
  { id: 'people-crest-Hm0-3', group: 'people', description: 'People at seawall/dike crest, clear view of the sea, Hm0 = 3 m', q: 0.3e-3, Vmax: 0.6 },
  { id: 'people-crest-Hm0-2', group: 'people', description: 'People at seawall/dike crest, clear view of the sea, Hm0 = 2 m', q: 1e-3, Vmax: 0.6 },
  { id: 'people-crest-Hm0-1', group: 'people', description: 'People at seawall/dike crest, clear view of the sea, Hm0 = 1 m', q: 10e-3, Vmax: 0.6 },
  { id: 'vehicles-crest-Hm0-3', group: 'vehicles', description: 'Cars on seawall/dike crest or railway close behind crest, Hm0 = 3 m', q: 5e-3, Vmax: 2.0 },
  { id: 'vehicles-crest-Hm0-2', group: 'vehicles', description: 'Cars on seawall/dike crest or railway close behind crest, Hm0 = 2 m', q: 10e-3, Vmax: 2.0 },
  { id: 'vehicles-crest-Hm0-1', group: 'vehicles', description: 'Cars on seawall/dike crest or railway close behind crest, Hm0 = 1 m', q: 75e-3, Vmax: 2.0 },
  { id: 'rubble-no-damage', group: 'structure', description: 'Rubble-mound breakwater, Hm0 > 5 m, no damage to crest and rear', q: 1e-3, Vmax: 2.0 },
  { id: 'rubble-rear-designed', group: 'structure', description: 'Rubble-mound breakwater, Hm0 > 5 m, rear side designed for overtopping', q: 5e-3, Vmax: 10.0 },
  { id: 'grass-maintained', group: 'structure', description: 'Grass-covered crest and landward slope, maintained closed sward, Hm0 = 1–3 m', q: 5e-3, Vmax: 2.0 },
  { id: 'grass-unmaintained', group: 'structure', description: 'Grass-covered crest and landward slope, not maintained, open spots, Hm0 = 0.5–3 m', q: 0.1e-3, Vmax: 0.5 },
  { id: 'building-elements', group: 'property', description: 'Building structure elements, Hm0 = 1–3 m', q: 1e-3, Vmax: 1.0 },
  { id: 'equipment-setback', group: 'property', description: 'Damage to equipment set back 5–10 m', q: 1e-3, Vmax: 1.0 },
]);

/**
 * Classifies a mean discharge (and optionally a maximum volume) against tolerable limits.
 * Reference: EurOtop (2018) chapter 3 (see TOLERABLE_OVERTOPPING).
 * @param {object} o
 * @param {number} o.q mean overtopping discharge [m³/s/m]
 * @param {number} [o.Vmax] maximum individual overtopping volume [m³/m]
 * @param {Array<{id:string,q:number,Vmax:number|null}>} [o.limits] table of limits, default TOLERABLE_OVERTOPPING
 * @returns {Array<{id:string,group:string,description:string,qLimit:number,VmaxLimit:number|null,
 *   utilisationQ:number,utilisationV:number,ok:boolean}>} one row per limit; utilisation = value/limit [-]
 */
export function classifyOvertopping({ q, Vmax, limits = TOLERABLE_OVERTOPPING }) {
  need(Number.isFinite(q) && q >= 0, 'q must be >= 0');
  return limits.map((l) => {
    const uq = q / l.q;
    const uv = Vmax !== undefined && l.Vmax ? Vmax / l.Vmax : Number.NaN;
    return {
      id: l.id, group: l.group, description: l.description, qLimit: l.q, VmaxLimit: l.Vmax ?? null,
      utilisationQ: uq, utilisationV: uv, ok: uq <= 1 && !(uv > 1),
    };
  });
}

/**
 * Distribution of individual overtopping volumes (two-parameter Weibull) and expected maximum volume:
 *   P(V_i <= V) = 1 − exp[−(V/a)^b],  a = q Tm /(Pov Γ(1 + 1/b)),
 *   b = 0.73 + 55 (q/(g Hm0 Tm−1,0))^0.8  (slopes),  Vmax = a (ln Now)^(1/b),
 *   Pov = Now/Nw = exp[−(√(−ln 0.02) Rc/Ru2%)²],  Nw = duration/Tm.
 * Reference: EurOtop (2018) chapter 5 (overtopping volumes: Rayleigh run-up for Pov, Weibull volumes).
 * Validity: sloping structures; b may be overridden (e.g. 0.75 of EurOtop 2007, or vertical-wall values).
 * @param {object} o
 * @param {number} o.q mean discharge [m³/s/m]
 * @param {number} o.Hm0 [m]
 * @param {number} o.Tm10 spectral period Tm−1,0 [s]
 * @param {number} o.Tm mean period [s]
 * @param {number} o.duration storm duration [s]
 * @param {number} [o.Rc] crest freeboard [m] (with Ru2 gives Pov)
 * @param {number} [o.Ru2] 2 % run-up level [m]
 * @param {number} [o.Pov] probability of overtopping per wave [-] (overrides Rc/Ru2)
 * @param {number} [o.b] Weibull shape override [-]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Pov:number,Nw:number,Now:number,a:number,b:number,Vmax:number,Vmean:number,
 *   exceedance:(V:number)=>number, volumeAt:(p:number)=>number}}
 *   a scale [m³/m], b shape [-], Vmax expected maximum volume in the duration [m³/m], Vmean mean volume
 *   per overtopping wave [m³/m]; exceedance(V) = P(V_i > V); volumeAt(p) = volume exceeded with probability p
 */
export function overtoppingVolumes({ q, Hm0, Tm10, Tm, duration, Rc, Ru2, Pov, b, g = G }) {
  pos(q, 'q');
  pos(Hm0, 'Hm0');
  pos(Tm10, 'Tm10');
  pos(Tm, 'Tm');
  pos(duration, 'duration');
  let p = Pov;
  if (p === undefined) {
    need(Number.isFinite(Rc) && Rc >= 0, 'Rc must be >= 0 when Pov is not given');
    pos(Ru2, 'Ru2');
    p = Math.exp(-Math.pow((Math.sqrt(-Math.log(0.02)) * Rc) / Ru2, 2));
  }
  need(p > 0 && p <= 1, 'Pov must lie in (0, 1]');
  const shape = b !== undefined ? b : 0.73 + 55 * Math.pow(q / (g * Hm0 * Tm10), 0.8);
  pos(shape, 'b');
  const Nw = duration / Tm;
  const Now = p * Nw;
  const Vmean = (q * Tm) / p;
  const a = Vmean / gammaFn(1 + 1 / shape);
  const Vmax = Now > 1 ? a * Math.pow(Math.log(Now), 1 / shape) : Number.NaN;
  return {
    Pov: p, Nw, Now, a, b: shape, Vmax, Vmean,
    exceedance: (V) => Math.exp(-Math.pow(Math.max(V, 0) / a, shape)),
    volumeAt: (prob) => a * Math.pow(-Math.log(prob), 1 / shape),
  };
}

/**
 * Sea-level-rise sensitivity of an overtopping (or run-up) estimate: for each rise Δ the freeboard becomes
 * Rc − Δ and the toe depth h + Δ; when `depthLimited` is true the wave height at the toe is scaled in
 * proportion to depth (Hm0 ∝ h, constant breaker index), otherwise it is unchanged.
 * Reference: first-principles bookkeeping around any overtopping formula; depth-limited scaling per the
 * constant breaker-index assumption for shallow foreshores.
 * @param {object} o
 * @param {number[]} o.rises sea-level rises Δ [m]
 * @param {number} o.Rc present crest freeboard [m]
 * @param {number} o.Hm0 present wave height at the toe [m]
 * @param {number} [o.h] present toe depth [m] (required when depthLimited)
 * @param {boolean} [o.depthLimited] scale Hm0 with depth
 * @param {(s:{Rc:number,Hm0:number,h:number})=>number} o.evaluate returns q [m³/s/m] for a state
 * @returns {Array<{rise:number,Rc:number,Hm0:number,h:number,q:number,factor:number,crestRaise:number}>}
 *   factor = q/q(Δ=0) [-]; crestRaise = crest raising that restores the present q [m] (by bisection,
 *   NaN if not bracketed within 0..10 m)
 */
export function seaLevelRiseSensitivity({ rises, Rc, Hm0, h, depthLimited = false, evaluate }) {
  need(Array.isArray(rises) || ArrayBuffer.isView(rises), 'rises must be an array');
  need(typeof evaluate === 'function', 'evaluate must be a function');
  need(Number.isFinite(Rc), 'Rc must be finite');
  pos(Hm0, 'Hm0');
  if (depthLimited) pos(h, 'h');
  const state = (d, raise) => ({
    Rc: Math.max(Rc - d + raise, 0),
    Hm0: depthLimited ? (Hm0 * (h + d)) / h : Hm0,
    h: h !== undefined ? h + d : Number.NaN,
  });
  const q0 = evaluate(state(0, 0));
  return Array.from(rises, (d) => {
    const s = state(d, 0);
    const q = evaluate(s);
    let crestRaise = 0;
    if (q > q0) {
      let lo = 0;
      let hi = 10;
      if (evaluate(state(d, hi)) > q0) crestRaise = Number.NaN;
      else {
        for (let i = 0; i < 60; i++) {
          const mid = 0.5 * (lo + hi);
          if (evaluate(state(d, mid)) > q0) lo = mid;
          else hi = mid;
        }
        crestRaise = 0.5 * (lo + hi);
      }
    }
    return { rise: d, Rc: s.Rc, Hm0: s.Hm0, h: s.h, q, factor: q / q0, crestRaise };
  });
}
