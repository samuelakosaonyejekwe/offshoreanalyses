// Gas-hydrate screening: equilibrium curves (gas-gravity correlations and a parameter-driven van der Waals–
// Platteeuw framework), thermodynamic inhibition, inhibitor injection balance, subcooling along a profile,
// kinetics (induction, growth, transport limits), slurry viscosity, plugging screening and depressurisation.
// SI units: T in K, P absolute in Pa, mass rates kg/s, standard gas rate Sm³/s.

import { brent } from '../math/roots.js';

const R = 8.314462618;
const PSI = 6894.757293168;

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}
const fToK = (F) => (F + 459.67) / 1.8;

/**
 * Hydrate equilibrium (dissociation) temperature of a sweet natural gas at pressure P from gas gravity.
 *  - 'motiee' (1991): T(°F) = −238.24469 + 78.99667 log10 P − 5.352544 (log10 P)² + 349.473877 γ − 150.854675 γ²
 *    − 27.604065 γ log10 P  (P psia).
 *  - 'towler' (Towler & Mokhatab 2005): T(°F) = 13.47 ln P + 34.27 ln γ − 1.675 ln P ln γ − 20.35  (P psia).
 *  - 'hammerschmidt' (1934): T(°F) = 8.9 P^0.285  (P psia; no gravity dependence, roughly a 0.6–0.7 gravity gas).
 * References: Motiee, Hydrocarbon Processing 70(7) (1991) 98–99; Towler & Mokhatab, Hydrocarbon Processing 84(4)
 * (2005) 61–62; Hammerschmidt, Ind. Eng. Chem. 26 (1934) 851–855. All are fits of the Katz (1945) gas-gravity
 * chart. Validity: sweet gas, γ 0.55–1.0, about 0.5–30 MPa, free water present, no inhibitor; accuracy ±1–3 K
 * (worse for gases with CO2, H2S or N2, and below 273 K). Screening only.
 * @param {object} o
 * @param {number} o.P pressure, Pa
 * @param {number} [o.gammaG=0.65] gas gravity (air = 1), –
 * @param {'motiee'|'towler'|'hammerschmidt'} [o.method='motiee']
 * @returns {number} equilibrium temperature, K
 */
export function hydrateTemperature({ P, gammaG = 0.65, method = 'motiee' }) {
  need('P', P, 1e4, 2e8);
  need('gammaG', gammaG, 0.5, 1.2);
  const p = P / PSI;
  switch (method) {
    case 'motiee': {
      const lp = Math.log10(p);
      return fToK(-238.24469 + 78.99667 * lp - 5.352544 * lp * lp + 349.473877 * gammaG - 150.854675 * gammaG * gammaG
        - 27.604065 * gammaG * lp);
    }
    case 'towler': {
      const lp = Math.log(p);
      const lg = Math.log(gammaG);
      return fToK(13.47 * lp + 34.27 * lg - 1.675 * lp * lg - 20.35);
    }
    case 'hammerschmidt': return fToK(8.9 * p ** 0.285);
    default: throw new RangeError(`method = ${method} is not one of motiee, towler, hammerschmidt`);
  }
}

/**
 * Hydrate equilibrium pressure at temperature T (inverse of `hydrateTemperature`, Brent on ln P between 0.1 and
 * 200 MPa). References and validity as `hydrateTemperature`.
 * @param {object} o
 * @param {number} o.T temperature, K @param {number} [o.gammaG=0.65] @param {string} [o.method='motiee']
 * @param {number} [o.inhibitorDT=0] equilibrium-temperature depression by inhibitor, K
 * @returns {number} equilibrium pressure, Pa (Infinity when T is above the correlation range at 200 MPa,
 *   the lower search bound 1e5 Pa when hydrates are stable even there)
 */
export function hydratePressure({ T, gammaG = 0.65, method = 'motiee', inhibitorDT = 0 }) {
  need('T', T, 200, 400);
  const f = (lnP) => hydrateTemperature({ P: Math.exp(lnP), gammaG, method }) - inhibitorDT - T;
  const lo = Math.log(1e5);
  const hi = Math.log(2e8);
  if (f(hi) < 0) return Infinity;
  if (f(lo) > 0) return 1e5;
  return Math.exp(brent(f, lo, hi, { tol: 1e-12 }));
}

/**
 * Thermodynamic inhibitor data: molar mass (kg/mol) and Hammerschmidt constant K_H in the original units
 * (°F·lb/lbmol with concentration in wt %). K_H = 2335 is Hammerschmidt's original value (methanol, ethanol);
 * for glycols values between 2335 and 4000 are quoted in the literature (GPSA uses 2335 for methanol and
 * recommends up to 4000 for glycols) — the glycol defaults here are editable mid-range values and must be confirmed
 * for design. In SI (ΔT in K, M in g/mol): K = K_H/1.8 (2335 → 1297, 2700 → 1500, 4000 → 2222).
 */
export const INHIBITORS = Object.freeze({
  methanol: Object.freeze({ M: 0.032042, KH: 2335 }),
  ethanol: Object.freeze({ M: 0.046068, KH: 2335 }),
  MEG: Object.freeze({ M: 0.062068, KH: 2700 }),
  DEG: Object.freeze({ M: 0.106120, KH: 4000 }),
});

function inhibitorData(inhibitor, M, KH) {
  const d = typeof inhibitor === 'string' ? INHIBITORS[inhibitor] : null;
  if (inhibitor != null && !d) throw new RangeError(`inhibitor = ${inhibitor} is not one of ${Object.keys(INHIBITORS).join(', ')}`);
  const m = M ?? d?.M;
  const k = KH ?? d?.KH;
  need('M', m, 0.01, 1);
  need('KH', k, 500, 10000);
  return { M: m, KH: k };
}

/**
 * Hammerschmidt (1939) hydrate-temperature depression by an inhibitor in the free-water phase:
 * ΔT(°F) = K_H w/[M (100 − w)], w = wt % inhibitor in the aqueous phase, M in g/mol; returned in kelvin (÷1.8).
 * Reference: Hammerschmidt, Gas 15(5) (1939) 30; GPSA Engineering Data Book §20.
 * Validity: methanol up to about 20–25 wt %, glycols up to about 60–70 wt %.
 * @param {object} o
 * @param {number} o.w inhibitor concentration in the aqueous phase, wt %
 * @param {string} [o.inhibitor='methanol'] key of `INHIBITORS` (or give M and KH)
 * @param {number} [o.M] molar mass, kg/mol @param {number} [o.KH] Hammerschmidt constant, °F basis
 * @returns {number} depression of the hydrate temperature, K
 */
export function hammerschmidtDepression({ w, inhibitor = 'methanol', M, KH }) {
  need('w', w, 0, 99.9);
  const d = inhibitorData(M != null && KH != null ? null : inhibitor, M, KH);
  return d.KH * w / (d.M * 1e3 * (100 - w)) / 1.8;
}

/**
 * Inhibitor concentration required for a given depression (inverse Hammerschmidt):
 * w = 100 ΔT M/(K_H + ΔT M) with ΔT in °F. Reference and validity as `hammerschmidtDepression`.
 * @param {object} o
 * @param {number} o.dT required depression, K @param {string} [o.inhibitor='methanol']
 * @param {number} [o.M] kg/mol @param {number} [o.KH]
 * @returns {number} wt % inhibitor in the aqueous phase
 */
export function hammerschmidtConcentration({ dT, inhibitor = 'methanol', M, KH }) {
  need('dT', dT, 0, 200);
  const d = inhibitorData(M != null && KH != null ? null : inhibitor, M, KH);
  const dF = dT * 1.8;
  return 100 * dF * d.M * 1e3 / (d.KH + dF * d.M * 1e3);
}

/**
 * Nielsen & Bucklin (1983) depression for concentrated methanol solutions:
 * ΔT(°C) = −72 ln(1 − x_inh) = −72 ln x_water (mole fractions in the aqueous phase; −129.6 ln x_w in °F).
 * Reference: Nielsen & Bucklin, Hydrocarbon Processing 62(4) (1983) 71–78; GPSA Engineering Data Book §20.
 * Validity: developed for methanol up to ~80 wt %; commonly applied to MEG with reduced accuracy.
 * @param {number} x mole fraction of inhibitor in the aqueous phase, –
 * @returns {number} depression, K
 */
export function nielsenBucklinDepression(x) {
  need('x', x, 0, 0.999);
  return -72 * Math.log(1 - x);
}

/**
 * Inhibitor mole fraction for a required depression (inverse Nielsen–Bucklin), x = 1 − exp(−ΔT/72).
 * @param {number} dT required depression, K @returns {number} mole fraction in the aqueous phase, –
 */
export function nielsenBucklinMoleFraction(dT) {
  need('dT', dT, 0, 300);
  return 1 - Math.exp(-dT / 72);
}

/**
 * Convert inhibitor mole fraction in water to mass fraction: w = x M/(x M + (1 − x) M_w). Reference: definition.
 * @param {number} x mole fraction, – @param {number} M inhibitor molar mass, kg/mol
 * @returns {number} mass fraction, –
 */
export function moleToMassFraction(x, M) {
  need('x', x, 0, 1);
  need('M', M, 0.001, 1);
  return x * M / (x * M + (1 - x) * 0.01801528);
}

/**
 * Convert inhibitor mass fraction in water to mole fraction (inverse of `moleToMassFraction`).
 * @param {number} w mass fraction, – @param {number} M inhibitor molar mass, kg/mol @returns {number} –
 */
export function massToMoleFraction(w, M) {
  need('w', w, 0, 1);
  need('M', M, 0.001, 1);
  return (w / M) / (w / M + (1 - w) / 0.01801528);
}

/**
 * Hydrate-temperature depression by dissolved salt, from the water-activity form of the Nielsen–Bucklin relation,
 * ΔT = −72 ln x_w, with x_w the ideal mole fraction of water counting each dissociated ion as a solute particle
 * (ν = 2 for NaCl). This is a colligative argument, not a fitted salt correlation: it reproduces measured NaCl
 * shifts for methane hydrate to within about 1 K up to roughly 10–15 wt % and under-predicts at higher salinity
 * where ion non-ideality lowers the water activity further.
 * Reference: Nielsen & Bucklin (1983) form; colligative (ideal-solution) water activity.
 * @param {object} o
 * @param {number} o.w salt concentration, wt % of the brine
 * @param {number} [o.M=0.05844] salt molar mass, kg/mol (NaCl) @param {number} [o.nu=2] ions per formula unit
 * @returns {number} depression, K
 */
export function saltDepression({ w, M = 0.05844, nu = 2 }) {
  need('w', w, 0, 30);
  need('M', M, 0.01, 1);
  need('nu', nu, 1, 5);
  const nSalt = nu * (w / 100) / M;
  const nW = (1 - w / 100) / 0.01801528;
  return -72 * Math.log(nW / (nW + nSalt));
}

/**
 * Inhibitor injection rate to protect a line to a target depression.
 * Aqueous-phase concentration w from Hammerschmidt (`method = 'hammerschmidt'`) or Nielsen–Bucklin
 * (`'nielsenBucklin'`); mass balance with an injected solution of inhibitor mass fraction `purity`:
 *   ṁ_inj = [w ṁ_w + (1 − w) ṁ_loss]/(purity − w),
 * where ṁ_w is the free (condensed + produced) water rate and ṁ_loss the inhibitor lost to the hydrocarbon phases,
 *   ṁ_loss = w (k_vap q_gas + k_cond ṁ_cond).
 * `kVap` is the vapour-phase loss in kg inhibitor per Sm³ gas per unit mass fraction in the aqueous phase and
 * `kCond` the condensate solubility in kg per kg condensate per unit mass fraction; both are arguments (default 0)
 * to be read for the outlet P and T from data such as the GPSA §20 methanol vapour-loss chart — methanol vapour losses are often
 * as large as the aqueous-phase requirement, glycol losses are small.
 * Reference: GPSA Engineering Data Book §20 (inhibitor mass balance, methanol-injection example).
 * @param {object} o
 * @param {number} o.dT required hydrate-temperature depression including margin, K
 * @param {string} [o.inhibitor='methanol'] @param {'hammerschmidt'|'nielsenBucklin'} [o.method='hammerschmidt']
 * @param {number} o.waterRate free-water mass rate at the cold end, kg/s
 * @param {number} [o.purity=1] inhibitor mass fraction of the injected solution, – (e.g. 0.9 lean MEG)
 * @param {number} [o.gasRate=0] gas rate, Sm³/s @param {number} [o.kVap=0] kg/Sm³ per mass fraction
 * @param {number} [o.condensateRate=0] kg/s @param {number} [o.kCond=0] kg/kg per mass fraction
 * @param {number} [o.M] kg/mol @param {number} [o.KH]
 * @returns {{ w: number, x: number, injectionRate: number, inhibitorRate: number, aqueousInhibitor: number,
 *   vapourLoss: number, condensateLoss: number }} w = required mass fraction in the aqueous phase (–), x = mole
 *   fraction, injected solution rate kg/s, pure inhibitor rate kg/s, inhibitor retained in the water phase kg/s,
 *   losses kg/s
 */
export function inhibitorInjection({ dT, inhibitor = 'methanol', method = 'hammerschmidt', waterRate, purity = 1,
  gasRate = 0, kVap = 0, condensateRate = 0, kCond = 0, M, KH }) {
  need('dT', dT, 0, 200);
  need('waterRate', waterRate, 0, 1e6);
  need('purity', purity, 1e-3, 1);
  need('gasRate', gasRate, 0, 1e9);
  need('kVap', kVap, 0, 10);
  need('condensateRate', condensateRate, 0, 1e6);
  need('kCond', kCond, 0, 10);
  const d = inhibitorData(M != null && KH != null ? null : inhibitor, M, KH);
  let w;
  if (method === 'hammerschmidt') w = hammerschmidtConcentration({ dT, M: d.M, KH: d.KH }) / 100;
  else if (method === 'nielsenBucklin') w = moleToMassFraction(nielsenBucklinMoleFraction(dT), d.M);
  else throw new RangeError(`method = ${method} is not one of hammerschmidt, nielsenBucklin`);
  if (w >= purity) throw new RangeError(`required aqueous concentration ${(w * 100).toFixed(1)} wt % exceeds the injected purity`);
  const vap = w * kVap * gasRate;
  const cond = w * kCond * condensateRate;
  const loss = vap + cond;
  const inj = (w * waterRate + (1 - w) * loss) / (purity - w);
  return { w, x: massToMoleFraction(w, d.M), injectionRate: inj, inhibitorRate: inj * purity,
    aqueousInhibitor: inj * purity - loss, vapourLoss: vap, condensateLoss: cond };
}

/**
 * Hydrate margin along a pressure–temperature profile (distance or time as the running coordinate):
 * subcooling ΔT_sub = T_eq(P) − ΔT_inhibitor − T (positive inside the hydrate region), the first entry into the
 * region (linear interpolation) and the maximum subcooling.
 * Reference: definition of subcooling; equilibrium by `hydrateTemperature` or a user curve.
 * @param {object} o
 * @param {number[]} o.s running coordinate (distance m, or time s), increasing
 * @param {number[]} o.P pressures, Pa @param {number[]} o.T temperatures, K
 * @param {number} [o.gammaG=0.65] @param {string} [o.method='motiee']
 * @param {number} [o.inhibitorDT=0] depression by inhibitor/salt, K
 * @param {function(number):number} [o.curve] custom equilibrium curve T_eq(P) in K (overrides the correlation)
 * @returns {{ Teq: number[], subcooling: number[], onset: number|null, maxSubcooling: number, atMax: number,
 *   inHydrateRegion: boolean, margin: number }} K; onset = coordinate of first entry (null if never);
 *   margin = −max subcooling (positive = distance from the hydrate curve), K
 */
export function subcoolingProfile({ s, P, T, gammaG = 0.65, method = 'motiee', inhibitorDT = 0, curve }) {
  if (!s || !P || !T || s.length !== P.length || s.length !== T.length || s.length === 0) {
    throw new RangeError('s, P and T must be non-empty arrays of equal length');
  }
  const Teq = [];
  const sub = [];
  let onset = null;
  let max = -Infinity;
  let atMax = s[0];
  for (let i = 0; i < s.length; i++) {
    const te = (curve ? curve(P[i]) : hydrateTemperature({ P: P[i], gammaG, method })) - inhibitorDT;
    Teq.push(te);
    sub.push(te - T[i]);
    if (sub[i] > max) { max = sub[i]; atMax = s[i]; }
    if (onset == null && sub[i] > 0) {
      onset = i === 0 ? s[0] : s[i - 1] + (s[i] - s[i - 1]) * (0 - sub[i - 1]) / (sub[i] - sub[i - 1]);
    }
  }
  return { Teq, subcooling: sub, onset, maxSubcooling: max, atMax, inHydrateRegion: max > 0, margin: -max };
}

/**
 * Nucleation induction time in the classical-nucleation-theory form t_ind = K exp(B/g), where the driving-force
 * group g is
 *   'supersaturation': g = T³ (ln S)²   (S = f/f_eq, fugacity or pressure ratio), from J ∝ exp[−16π σ³ v²/(3 k³ T³ ln²S)];
 *   'subcooling':      g = T ΔT²        (Δμ ≈ Δs ΔT at constant pressure), from ΔG_crit/(kT) ∝ σ³ v²/(k T Δs² ΔT²).
 * K (s) and B (K³ in both forms) are system-specific kinetic constants that must be fitted to induction-time
 * data; no defaults are implied. Induction is stochastic — treat
 * the result as a characteristic value.
 * Reference: Kashchiev & Firoozabadi, J. Crystal Growth 250 (2003) 499–515; Volmer–Weber/Becker–Döring classical
 * nucleation theory.
 * @param {object} o
 * @param {number} o.K pre-exponential time, s @param {number} o.B exponent constant, K³
 * @param {number} o.T temperature, K
 * @param {number} [o.S] supersaturation ratio (> 1) @param {number} [o.dT] subcooling, K (> 0)
 * @returns {number} induction time, s (Infinity without driving force)
 */
export function inductionTime({ K, B, T, S, dT }) {
  need('K', K, 0, 1e30);
  need('B', B, 0, 1e30);
  need('T', T, 100, 400);
  if ((S == null) === (dT == null)) throw new RangeError('give exactly one of S and dT');
  let g;
  if (S != null) {
    need('S', S, 0, 1e6);
    if (S <= 1) return Infinity;
    g = T ** 3 * Math.log(S) ** 2;
  } else {
    if (dT <= 0) return Infinity;
    g = T * dT * dT;
  }
  return K * Math.exp(B / g);
}

/**
 * Intrinsic hydrate growth (gas consumption) rate, Englezos et al. (1987): dn/dt = K* A_p (f − f_eq), with the
 * combined rate constant 1/K* = 1/k_r + 1/k_d (reaction and diffusion around the particle).
 * Reference: Englezos, Kalogerakis, Dholabhai & Bishnoi, Chem. Eng. Sci. 42 (1987) 2647–2658. All kinetic
 * constants are arguments (reported intrinsic values span orders of magnitude between apparatus).
 * @param {object} o
 * @param {number} o.area particle (or interfacial) surface area, m²
 * @param {number} o.f gas fugacity, Pa @param {number} o.feq three-phase equilibrium fugacity at T, Pa
 * @param {number} [o.Kstar] combined rate constant, mol/(m²·Pa·s) — or give kr and kd
 * @param {number} [o.kr] reaction rate constant, mol/(m²·Pa·s) @param {number} [o.kd] diffusion constant
 * @returns {number} gas consumption rate, mol/s (0 when f ≤ f_eq)
 */
export function growthRateEnglezos({ area, f, feq, Kstar, kr, kd }) {
  need('area', area, 0, 1e12);
  need('f', f, 0, 1e10);
  need('feq', feq, 0, 1e10);
  const k = Kstar ?? 1 / (1 / need('kr', kr, 1e-30, 1e10) + 1 / need('kd', kd, 1e-30, 1e10));
  need('Kstar', k, 0, 1e10);
  return Math.max(0, k * area * (f - feq));
}

/**
 * Semi-empirical hydrate formation rate of Vysniauskas & Bishnoi (1983):
 * r = A a_s exp(−ΔE_a/(R T)) exp(−a/ΔT^b) P^γ, with pre-exponential A, activation energy ΔE_a, the subcooling
 * constants a, b and the pressure order γ all supplied by the user (the published values are specific to the
 * authors' stirred reactor and unit system).
 * Reference: Vysniauskas & Bishnoi, Chem. Eng. Sci. 38 (1983) 1061–1072.
 * @param {object} o
 * @param {number} o.A pre-exponential constant (units fix those of the result)
 * @param {number} o.as interfacial area, m² @param {number} o.Ea activation energy, J/mol (may be negative)
 * @param {number} o.T K @param {number} o.dT subcooling, K @param {number} o.a K^b @param {number} o.b –
 * @param {number} o.P pressure, Pa @param {number} o.gamma pressure order, –
 * @returns {number} formation rate in the units implied by A (0 when ΔT ≤ 0)
 */
export function growthRateVysniauskasBishnoi({ A, as, Ea, T, dT, a, b, P, gamma }) {
  need('A', A, 0, 1e300);
  need('as', as, 0, 1e12);
  need('Ea', Ea, -1e7, 1e7);
  need('T', T, 100, 400);
  need('a', a, 0, 1e6);
  need('b', b, 0, 10);
  need('P', P, 0, 1e9);
  need('gamma', gamma, 0, 10);
  if (!(dT > 0)) return 0;
  return A * as * Math.exp(-Ea / (R * T)) * Math.exp(-a / dT ** b) * P ** gamma;
}

/**
 * Transport-limited hydrate growth by the film model: the gas consumption rate is the smallest of
 *   mass transfer:  ṅ = k_L A (c_bulk − c_eq)                  (dissolved gas to the hydrate surface),
 *   heat transfer:  ṅ = h A ΔT_sub/Δh_f                         (removal of the heat of formation),
 * and, if given, the intrinsic kinetic rate. With series resistances the actual rate is slightly below this
 * minimum; `series` gives the harmonic combination.
 * Reference: two-film theory (Lewis & Whitman 1924); Skovborg & Rasmussen, Chem. Eng. Sci. 49 (1994) 1131 (mass-
 * transfer-limited growth); heat-transfer limit after Uchida et al. / Mori (2001) shell models.
 * @param {object} o
 * @param {number} o.area interfacial area, m²
 * @param {number} [o.kL] liquid-side mass-transfer coefficient, m/s
 * @param {number} [o.cBulk] dissolved-gas concentration in the bulk, mol/m³ @param {number} [o.cEq] at equilibrium
 * @param {number} [o.h] heat-transfer coefficient, W/(m²·K) @param {number} [o.dTsub] subcooling, K
 * @param {number} [o.dHf=54000] heat of hydrate formation per mole of gas, J/mol (≈ 54 kJ/mol for methane sI)
 * @param {number} [o.intrinsic] intrinsic kinetic rate, mol/s
 * @returns {{ rate: number, series: number, massTransfer: number|null, heatTransfer: number|null,
 *   limiting: string }} mol gas/s
 */
export function transportLimitedGrowth({ area, kL, cBulk, cEq, h, dTsub, dHf = 54000, intrinsic }) {
  need('area', area, 0, 1e12);
  const mt = kL != null ? Math.max(0, need('kL', kL, 0, 10) * area * (need('cBulk', cBulk, 0, 1e6) - need('cEq', cEq, 0, 1e6))) : null;
  const ht = h != null ? Math.max(0, need('h', h, 0, 1e6) * area * need('dTsub', dTsub, -500, 500) / need('dHf', dHf, 1, 1e7)) : null;
  const cands = [['massTransfer', mt], ['heatTransfer', ht], ['intrinsic', intrinsic ?? null]].filter((c) => c[1] != null);
  if (cands.length === 0) throw new RangeError('give at least one of the mass-transfer, heat-transfer or intrinsic rates');
  let rate = Infinity;
  let limiting = '';
  let inv = 0;
  for (const [name, r] of cands) {
    if (r < rate) { rate = r; limiting = name; }
    inv += r > 0 ? 1 / r : Infinity;
  }
  return { rate, series: inv === Infinity ? 0 : 1 / inv, massTransfer: mt, heatTransfer: ht, limiting };
}

/**
 * Volume of hydrate formed per unit volume of liquid water converted, from the hydration number and densities:
 * V_h/V_w = (ρ_w/ρ_h)(n M_w + M_g)/(n M_w).
 * Reference: stoichiometry of G·nH2O (n ≈ 5.75 ideal sI, ≈ 6.0 typical for methane; ρ_h ≈ 910–940 kg/m³).
 * @param {object} [o]
 * @param {number} [o.n=6] hydration number @param {number} [o.Mg=0.016043] guest molar mass, kg/mol
 * @param {number} [o.rhoW=1000] kg/m³ @param {number} [o.rhoH=910] kg/m³
 * @returns {number} volume ratio, –
 */
export function hydrateVolumeRatio({ n = 6, Mg = 0.016043, rhoW = 1000, rhoH = 910 } = {}) {
  need('n', n, 4, 20);
  need('Mg', Mg, 0.002, 0.2);
  need('rhoW', rhoW, 500, 1500);
  need('rhoH', rhoH, 500, 1500);
  return (rhoW / rhoH) * (n * 0.01801528 + Mg) / (n * 0.01801528);
}

/**
 * Hydrate volume fraction along a pipe by one-dimensional advective bookkeeping of water conversion.
 * A liquid parcel moves at velocity v; where the subcooling exceeds `dTonset`, water converts at
 *   dX/dt = k ΔT_sub (1 − X),
 * with k (1/(s·K)) a lumped first-order kinetic/transport constant (argument — fit to flow-loop data or derive
 * from a growth model and interfacial area). Hydrate volume fraction of the liquid + hydrate slurry:
 *   φ = X α_w r_v / [1 − α_w X + X α_w r_v],  α_w = water volume fraction of the liquid, r_v = `hydrateVolumeRatio`.
 * No heat release feedback, no deposition, no slip between hydrate and liquid.
 * Reference: plug-flow material balance; rate form after the subcooling-driven models used in transient
 * flow-assurance simulators (Turner et al. 2005 kinetics in lumped form). Screening only.
 * @param {object} o
 * @param {number[]} o.x distance along the pipe, m (increasing)
 * @param {number[]} o.subcooling ΔT_sub at each point, K
 * @param {number} o.v liquid velocity, m/s (> 0)
 * @param {number} o.k rate constant, 1/(s·K)
 * @param {number} o.waterFraction water volume fraction of the liquid phase, –
 * @param {number} [o.dTonset=0] subcooling needed before formation starts, K (nucleation threshold; argument)
 * @param {number} [o.volumeRatio] hydrate/water volume ratio (default `hydrateVolumeRatio()`)
 * @returns {{ conversion: number[], phi: number[], onset: number|null }} water conversion X (–), hydrate volume
 *   fraction in the slurry (–), distance where formation starts (m)
 */
export function hydrateFractionAlongPipe({ x, subcooling, v, k, waterFraction, dTonset = 0, volumeRatio }) {
  if (!x || !subcooling || x.length !== subcooling.length || x.length < 2) throw new RangeError('x and subcooling must be arrays of equal length ≥ 2');
  need('v', v, 1e-9, 1e3);
  need('k', k, 0, 1e3);
  need('waterFraction', waterFraction, 0, 1);
  need('dTonset', dTonset, 0, 100);
  const rv = volumeRatio ?? hydrateVolumeRatio();
  const X = [0];
  let onset = null;
  let started = false;
  for (let i = 1; i < x.length; i++) {
    const dx = x[i] - x[i - 1];
    if (!(dx > 0)) throw new RangeError(`x must increase strictly (point ${i})`);
    const sub = 0.5 * (subcooling[i] + subcooling[i - 1]);
    if (!started && Math.max(subcooling[i], subcooling[i - 1]) > dTonset) { started = true; onset = x[i - 1]; }
    // exact integration of dX/dt = k ΔT (1 − X) over the residence time of the interval
    const rate = started && sub > 0 ? k * sub : 0;
    X.push(1 - (1 - X[i - 1]) * Math.exp(-rate * dx / v));
  }
  const phi = X.map((c) => c * waterFraction * rv / (1 - waterFraction * c + c * waterFraction * rv));
  return { conversion: X, phi, onset };
}

/**
 * Relative viscosity of a hard-sphere suspension, Krieger & Dougherty (1959):
 * μ_r = (1 − φ/φ_max)^(−[η] φ_max), [η] = 2.5 for spheres.
 * Reference: Krieger & Dougherty, Trans. Soc. Rheol. 3 (1959) 137–152. Validity: φ < φ_max (0.64 random close
 * packing; lower effective values for aggregating particles).
 * @param {object} o
 * @param {number} o.phi particle volume fraction, – @param {number} [o.phiMax=0.64] @param {number} [o.intrinsic=2.5]
 * @returns {number} μ_slurry/μ_liquid, – (Infinity at φ ≥ φ_max)
 */
export function slurryViscosityKriegerDougherty({ phi, phiMax = 0.64, intrinsic = 2.5 }) {
  need('phi', phi, 0, 1);
  need('phiMax', phiMax, 0.05, 1);
  need('intrinsic', intrinsic, 0.5, 20);
  if (phi >= phiMax) return Infinity;
  return (1 - phi / phiMax) ** (-intrinsic * phiMax);
}

/**
 * Hydrate-slurry relative viscosity with shear-dependent aggregation, Camargo & Palermo (2002):
 *   μ_r = (1 − φ_eff)/(1 − φ_eff/φ_max)²                                   (Mills 1985, φ_max = 4/7)
 *   φ_eff = φ (d_A/d_p)^(3 − f)                                             (fractal aggregates, f ≈ 2.5)
 *   (d_A/d_p)^(4 − f) − F_a [1 − (φ/φ_max)(d_A/d_p)^(3 − f)]² / {d_p² μ_o γ̇ [1 − φ (d_A/d_p)^(3 − f)]} = 0
 * solved for the aggregate-to-particle diameter ratio (≥ 1); F_a is the inter-particle adhesion force, an argument
 * to be fitted (it scales with particle size and interfacial tension), as are d_p and f.
 * Reference: Camargo & Palermo, Proc. 4th Int. Conf. on Gas Hydrates, Yokohama (2002) 880–885; Mills, J. Phys.
 * Lett. 46 (1985) L301.
 * Validity: oil-continuous slurries of fully converted water droplets, laminar shear; returns Infinity when the
 * effective fraction reaches φ_max (no flow).
 * @param {object} o
 * @param {number} o.phi hydrate particle volume fraction, –
 * @param {number} o.dp primary particle (droplet) diameter, m
 * @param {number} o.muOil oil viscosity, Pa·s @param {number} o.shearRate γ̇, 1/s
 * @param {number} o.Fa adhesion force between particles, N
 * @param {number} [o.fractal=2.5] fractal dimension @param {number} [o.phiMax=0.5714285714285714]
 * @returns {{ muRel: number, phiEff: number, dRatio: number }} relative viscosity, effective fraction, d_A/d_p
 */
export function slurryViscosityCamargoPalermo({ phi, dp, muOil, shearRate, Fa, fractal = 2.5, phiMax = 4 / 7 }) {
  need('phi', phi, 0, 1);
  need('dp', dp, 1e-9, 1e-1);
  need('muOil', muOil, 1e-6, 1e3);
  need('shearRate', shearRate, 1e-9, 1e7);
  need('Fa', Fa, 0, 1);
  need('fractal', fractal, 1, 3);
  need('phiMax', phiMax, 0.05, 1);
  const mills = (pe) => (pe >= phiMax ? Infinity : (1 - pe) / (1 - pe / phiMax) ** 2);
  if (phi === 0) return { muRel: 1, phiEff: 0, dRatio: 1 };
  const e3 = 3 - fractal;
  const g = (r) => {
    const pe = phi * r ** e3;
    return r ** (4 - fractal) - Fa * (1 - pe / phiMax) ** 2 / (dp * dp * muOil * shearRate * (1 - pe));
  };
  // admissible ratios: 1 ≤ r < r_max where φ_eff = φ_max
  const rMax = e3 > 0 ? (phiMax / phi) ** (1 / e3) : Infinity;
  if (phi >= phiMax) return { muRel: Infinity, phiEff: phi, dRatio: 1 };
  let r = 1;
  if (g(1) < 0) {
    const hi = Number.isFinite(rMax) ? rMax * (1 - 1e-12) : 1e6;
    // g increases with r and is positive at r_max (second term vanishes)
    r = g(hi) <= 0 ? hi : brent(g, 1, hi, { tol: 1e-12 });
  }
  const phiEff = phi * r ** e3;
  return { muRel: mills(phiEff), phiEff, dRatio: r };
}

/**
 * Hydrate plugging-risk index (0–1) — a screening heuristic, NOT a published or validated criterion. It combines
 *   thermodynamic exposure  a = clamp(ΔT_sub/ΔT_ref),
 *   hydrate loading         b = clamp(φ_eff/φ_plug),
 *   transport               c = clamp(v_ref/v)   (low velocity favours bedding and deposition),
 * as index = a · b · (0.5 + 0.5 c), with classes 'none' (no subcooling), < 0.2 'low', < 0.5 'moderate', ≥ 0.5 'high'.
 * All reference values are arguments with illustrative defaults; calibrate against field or flow-loop experience
 * before use in decisions.
 * @param {object} o
 * @param {number} o.subcooling K @param {number} o.phi hydrate (effective) volume fraction in the liquid, –
 * @param {number} [o.velocity] liquid velocity, m/s
 * @param {number} [o.dTref=3.5] subcooling at which formation is taken as certain, K
 * @param {number} [o.phiPlug=0.4] loading at which the slurry is taken as non-transportable, –
 * @param {number} [o.vRef=1] velocity below which bedding is considered likely, m/s
 * @returns {{ index: number, class: 'none'|'low'|'moderate'|'high', exposure: number, loading: number,
 *   transport: number }}
 */
export function pluggingRiskIndex({ subcooling, phi, velocity, dTref = 3.5, phiPlug = 0.4, vRef = 1 }) {
  need('subcooling', subcooling, -500, 500);
  need('phi', phi, 0, 1);
  need('dTref', dTref, 1e-3, 100);
  need('phiPlug', phiPlug, 1e-3, 1);
  need('vRef', vRef, 1e-3, 100);
  const clamp = (v) => Math.min(1, Math.max(0, v));
  const a = clamp(subcooling / dTref);
  const b = clamp(phi / phiPlug);
  const c = velocity != null ? clamp(vRef / Math.max(need('velocity', velocity, 0, 1e3), 1e-9)) : 1;
  const index = a * b * (0.5 + 0.5 * c);
  const cls = !(subcooling > 0) ? 'none' : index < 0.2 ? 'low' : index < 0.5 ? 'moderate' : 'high';
  return { index, class: cls, exposure: a, loading: b, transport: c };
}

/**
 * Hydrate dissociation by depressurisation: equilibrium temperature at the reduced pressure, the driving force
 * for melting relative to the surroundings, and an ice-formation warning. A dissociating plug cools towards the
 * equilibrium temperature at the local pressure; if that is below the ice point the released water freezes and
 * dissociation stalls (and one-sided depressurisation can launch the plug).
 * Reference: Sloan & Koh, "Clathrate Hydrates of Natural Gases", 3rd ed. (2008) ch. 8; equilibrium curve by
 * `hydrateTemperature` (note the gas-gravity correlations are least accurate below 273 K).
 * @param {object} o
 * @param {number} o.P reduced (blow-down) pressure, Pa @param {number} o.Tamb ambient temperature, K
 * @param {number} [o.gammaG=0.65] @param {string} [o.method='motiee'] @param {number} [o.inhibitorDT=0] K
 * @param {number} [o.Tice=273.15] ice point, K
 * @returns {{ Teq: number, drivingForce: number, dissociates: boolean, iceRisk: boolean, PiceLimit: number,
 *   Pdissociation: number }} Teq at P (K); T_amb − T_eq (K); whether the plug melts at ambient; whether T_eq is
 *   below the ice point; pressure at which T_eq equals the ice point (Pa) — stay above it to avoid ice;
 *   pressure below which hydrates are unstable at ambient temperature (Pa)
 */
export function depressurisation({ P, Tamb, gammaG = 0.65, method = 'motiee', inhibitorDT = 0, Tice = 273.15 }) {
  need('Tamb', Tamb, 200, 400);
  const Teq = hydrateTemperature({ P, gammaG, method }) - inhibitorDT;
  return {
    Teq, drivingForce: Tamb - Teq, dissociates: Tamb > Teq, iceRisk: Teq < Tice,
    PiceLimit: hydratePressure({ T: Tice, gammaG, method, inhibitorDT }),
    Pdissociation: hydratePressure({ T: Tamb, gammaG, method, inhibitorDT }),
  };
}

const CAVITIES = {
  // cavities per water molecule: structure I — 2 small + 6 large per 46 H2O; structure II — 16 small + 8 large per 136
  sI: [1 / 23, 3 / 23],
  sII: [2 / 17, 1 / 17],
};

/**
 * van der Waals–Platteeuw statistical-thermodynamic hydrate equilibrium (framework with user-supplied parameters).
 * Hydrate phase:  Δμ_w^H/(RT) = Σ_m ν_m ln(1 + Σ_j C_mj f_j)   (Langmuir adsorption in small and large cavities),
 * Liquid water:   Δμ_w^L/(RT) = Δμ_0/(R T_0) − ∫_{T_0}^{T} Δh/(R T²) dT + Δv P/(R T) − ln a_w,
 *                 Δh = Δh_0 + Δc_p (T − T_0)   (Holder et al. 1980 reference-property form, Δv taken constant),
 * equilibrium where the two are equal. Langmuir constants C_mj(T) = (A_mj/T) exp(B_mj/T) (Parrish–Prausnitz /
 * Munck form) or any function of T; gas fugacities f_j = φ_j y_j P with φ from `fugacity` (ideal gas by default —
 * pass an EOS callback for real gases).
 * NO parameter set is built in: the Langmuir constants and the reference properties (Δμ_0, Δh_0, Δc_p, Δv) must
 * be taken together from one published source (e.g. Parrish & Prausnitz 1972; Munck, Skjold-Jørgensen & Rasmussen
 * 1988) because they are fitted jointly.
 * References: van der Waals & Platteeuw, Adv. Chem. Phys. 2 (1959) 1; Parrish & Prausnitz, Ind. Eng. Chem.
 * Process Des. Dev. 11 (1972) 26; Holder, Corbin & Papadopoulos, Ind. Eng. Chem. Fundam. 19 (1980) 282; Munck et
 * al., Chem. Eng. Sci. 43 (1988) 2661.
 * Validity: liquid-water–hydrate–vapour equilibrium above the ice point; single hydrate structure per call.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {'sI'|'sII'} o.structure
 * @param {{ y: number, A?: number[], B?: number[], C?: function(number):number[] }[]} o.guests per guest: vapour
 *   mole fraction and Langmuir parameters for [small, large] cavities — A in K/Pa and B in K, or C(T) in 1/Pa
 * @param {{ dmu0: number, dh0: number, dcp?: number, dv: number, T0?: number }} o.ref reference properties of the
 *   empty lattice relative to liquid water: Δμ_0 J/mol, Δh_0 J/mol, Δc_p J/(mol·K), Δv m³/mol, T_0 K (273.15)
 * @param {function(number):number[]} [o.fugacity] fugacity coefficients φ_j(P) of the guests
 * @param {number} [o.waterActivity=1] activity of water (inhibitors, salts)
 * @param {number} [o.Pmin=1e4] @param {number} [o.Pmax=5e8] search bracket, Pa
 * @returns {{ P: number, occupancy: number[][], dmuRT: number }} equilibrium pressure Pa, cage occupancies
 *   θ[cavity][guest], and Δμ_w/(RT) at equilibrium
 */
export function vdwpEquilibriumPressure({ T, structure, guests, ref, fugacity, waterActivity = 1, Pmin = 1e4, Pmax = 5e8 }) {
  need('T', T, 200, 350);
  const nu = CAVITIES[structure];
  if (!nu) throw new RangeError(`structure = ${structure} is not one of sI, sII`);
  if (!Array.isArray(guests) || guests.length === 0) throw new RangeError('guests must be a non-empty array');
  if (!ref) throw new RangeError('ref (reference properties) is required');
  need('ref.dmu0', ref.dmu0, 0, 1e5);
  need('ref.dh0', ref.dh0, -1e6, 1e6);
  need('ref.dv', ref.dv, -1e-3, 1e-3);
  need('waterActivity', waterActivity, 1e-3, 1);
  const T0 = ref.T0 ?? 273.15;
  const dcp = ref.dcp ?? 0;
  const C = guests.map((g, j) => {
    need(`guests[${j}].y`, g.y, 0, 1);
    const c = g.C ? g.C(T) : [g.A[0] / T * Math.exp(g.B[0] / T), g.A[1] / T * Math.exp(g.B[1] / T)];
    if (!(c[0] >= 0) || !(c[1] >= 0)) throw new RangeError(`guests[${j}] has invalid Langmuir constants`);
    return c;
  });
  const enth = ((ref.dh0 - dcp * T0) * (1 / T0 - 1 / T) + dcp * Math.log(T / T0)) / R;
  const liquid = (P) => ref.dmu0 / (R * T0) - enth + ref.dv * P / (R * T) - Math.log(waterActivity);
  const sums = (P) => {
    const phi = fugacity ? fugacity(P) : guests.map(() => 1);
    const s = [0, 0];
    const terms = [[], []];
    guests.forEach((g, j) => {
      const f = phi[j] * g.y * P;
      for (let m = 0; m < 2; m++) { terms[m].push(C[j][m] * f); s[m] += C[j][m] * f; }
    });
    return { s, terms };
  };
  const hydrate = (P) => { const { s } = sums(P); return nu[0] * Math.log1p(s[0]) + nu[1] * Math.log1p(s[1]); };
  const F = (lnP) => hydrate(Math.exp(lnP)) - liquid(Math.exp(lnP));
  const a = Math.log(Pmin);
  const b = Math.log(Pmax);
  if (F(a) > 0) throw new RangeError('hydrate is already stable at Pmin: lower Pmin');
  if (F(b) < 0) throw new RangeError('no hydrate equilibrium below Pmax with these parameters');
  const P = Math.exp(brent(F, a, b, { tol: 1e-13 }));
  const { s, terms } = sums(P);
  return { P, occupancy: terms.map((row, m) => row.map((v) => v / (1 + s[m]))), dmuRT: liquid(P) };
}
