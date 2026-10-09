// Wax (paraffin) deposition screening: solubility curves, wax appearance temperature from an ideal-solution
// solid–liquid model, diffusion coefficients, deposition flux (molecular diffusion, shear dispersion, Matzain
// factors), deposit growth along a line with thermal feedback, pressure-drop increase, pigging interval and gel
// restart pressure.
// SI units: T in K, lengths m, mass flux kg/(m²·s), diffusivity m²/s, viscosity Pa·s.

import { brent } from '../math/roots.js';
import { nusseltInternal } from './thermalpipe.js';

const R = 8.314462618;
const CAL = 4.184;

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}

/**
 * Build a wax solubility curve C(T) (dissolved wax as mass fraction of the oil) and its temperature derivative.
 *  - `table`: measured points [{ T, C }] (e.g. from DSC or filtration), piecewise-linear
 *    interpolation, constant outside the table (dC/dT = 0 there);
 *  - `exponential`: C = a·exp(b (T − Tref)), a common fit of the precipitation curve below the WAT.
 * Above the wax appearance temperature `WAT` (if given) the oil is undersaturated and dC/dT is returned as 0.
 * Reference: representation of laboratory precipitation curves; no correlation implied.
 * @param {object} o
 * @param {{ T: number, C: number }[]} [o.table] temperature K, solubility kg/kg (T increasing)
 * @param {{ a: number, b: number, Tref: number }} [o.exponential] a kg/kg, b 1/K, Tref K
 * @param {number} [o.WAT] wax appearance temperature, K
 * @returns {{ C: function(number):number, dCdT: function(number):number, WAT: number|null }}
 */
export function makeSolubilityCurve({ table, exponential, WAT }) {
  if ((table == null) === (exponential == null)) throw new RangeError('give exactly one of table and exponential');
  let C;
  let dCdT;
  if (table) {
    if (!Array.isArray(table) || table.length < 2) throw new RangeError('table needs at least two points');
    for (let i = 1; i < table.length; i++) {
      if (!(table[i].T > table[i - 1].T)) throw new RangeError(`table T must increase strictly (row ${i})`);
    }
    const seg = (T) => {
      let i = 1;
      while (i < table.length - 1 && T > table[i].T) i++;
      return i;
    };
    C = (T) => {
      if (T <= table[0].T) return table[0].C;
      if (T >= table[table.length - 1].T) return table[table.length - 1].C;
      const i = seg(T);
      return table[i - 1].C + (table[i].C - table[i - 1].C) * (T - table[i - 1].T) / (table[i].T - table[i - 1].T);
    };
    dCdT = (T) => {
      if (T < table[0].T || T > table[table.length - 1].T) return 0;
      const i = seg(T);
      return (table[i].C - table[i - 1].C) / (table[i].T - table[i - 1].T);
    };
  } else {
    const { a, b, Tref } = exponential;
    need('exponential.a', a, 0, 1);
    need('exponential.b', b, 0, 1);
    need('exponential.Tref', Tref, 100, 500);
    C = (T) => a * Math.exp(b * (T - Tref));
    dCdT = (T) => a * b * Math.exp(b * (T - Tref));
  }
  const wat = WAT ?? null;
  return {
    C: (T) => (wat != null && T > wat ? C(wat) : C(T)),
    dCdT: (T) => (wat != null && T > wat ? 0 : dCdT(T)),
    WAT: wat,
  };
}

/**
 * Melting temperature and heat of fusion of an n-paraffin (or paraffinic pseudo-component) from molar mass by the
 * Won (1986) correlations: T_f = 374.5 + 0.02617 M − 20172/M (K), ΔH_f = 0.1426 M T_f (cal/mol), M in g/mol.
 * Reference: Won, K.W., Fluid Phase Equilib. 30 (1986) 265–279. Validity: n-alkanes of about C10–C60; ΔH_f lumps
 * the solid–solid transition and is too high for branched and cyclic species.
 * @param {number} M molar mass, kg/mol
 * @returns {{ Tf: number, dHf: number }} K and J/mol
 */
export function wonFusionProperties(M) {
  need('M', M, 0.1, 1.2);
  const Mg = M * 1e3;
  const Tf = 374.5 + 0.02617 * Mg - 20172 / Mg;
  return { Tf, dHf: 0.1426 * Mg * Tf * CAL };
}

function solidLiquidK(comp, T) {
  const p = comp.Tf != null && comp.dHf != null ? comp : wonFusionProperties(comp.M);
  // ideal solution in both phases, heat-capacity and pressure terms neglected
  return Math.exp(p.dHf / (R * T) * (1 - T / p.Tf));
}

/**
 * Wax precipitation by an ideal-solution solid–liquid equilibrium: K_i = s_i/x_i = exp[ΔH_f,i/(R T) (1 − T/T_f,i)]
 * and a Rachford–Rice balance for the solid mole fraction. Components may carry their own Tf and dHf; otherwise the
 * Won (1986) correlations are used from M.
 * Reference: Won (1986) with activity coefficients set to 1 (Pedersen, Skovborg & Rønningsen, Energy & Fuels 5
 * (1991) 924 discuss the resulting over-prediction of wax amount and WAT).
 * Validity: screening; an ideal solid solution over-predicts the WAT (typically by 5–20 K) and the wax quantity
 * for real crudes — prefer a measured WAT and precipitation curve when available.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {{ M: number, z: number, Tf?: number, dHf?: number }[]} o.components molar mass kg/mol, feed mole
 *   fraction, optional melting temperature K and heat of fusion J/mol
 * @returns {{ solidFraction: number, x: number[], s: number[], K: number[] }} solid mole fraction of the feed,
 *   liquid and solid compositions (the incipient solid composition when no solid is present), K-values
 */
export function waxPrecipitationIdeal({ T, components }) {
  need('T', T, 150, 500);
  if (!Array.isArray(components) || components.length === 0) throw new RangeError('components must be a non-empty array');
  let zs = 0;
  for (const c of components) zs += need('z', c.z, 0, 1);
  if (Math.abs(zs - 1) > 1e-6) throw new RangeError(`mole fractions sum to ${zs}, expected 1`);
  const K = components.map((c) => solidLiquidK(c, T));
  const z = components.map((c) => c.z / zs);
  const g = (S) => z.reduce((sum, zi, i) => sum + zi * (K[i] - 1) / (1 + S * (K[i] - 1)), 0);
  let S = 0;
  if (g(0) > 0) S = g(1) >= 0 ? 1 : brent(g, 0, 1, { tol: 1e-14 });
  const x = z.map((zi, i) => zi / (1 + S * (K[i] - 1)));
  const s = x.map((xi, i) => K[i] * xi);
  const ss = s.reduce((a, b) => a + b, 0);
  return { solidFraction: S, x, s: s.map((v) => v / ss), K };
}

/**
 * Wax appearance temperature from the ideal-solution model: the temperature at which Σ z_i K_i(T) = 1.
 * Reference and validity as `waxPrecipitationIdeal` (over-predicts for real crudes).
 * @param {object} o
 * @param {{ M: number, z: number, Tf?: number, dHf?: number }[]} o.components
 * @returns {number} WAT, K
 */
export function waxAppearanceIdeal({ components }) {
  const f = (T) => components.reduce((s, c) => s + c.z * solidLiquidK(c, T), 0) - 1;
  const lo = 150;
  const hi = 450;
  if (f(lo) < 0) throw new RangeError('no wax phase above 150 K for this composition');
  if (f(hi) > 0) throw new RangeError('composition is solid above 450 K: check the n-paraffin distribution');
  return brent(f, lo, hi, { tol: 1e-10 });
}

/**
 * Diffusion coefficient of wax in oil, Hayduk & Minhas (1982) correlation for normal paraffin solutions:
 * D = 13.3e-8 T^1.47 μ^(10.2/V_A − 0.791) / V_A^0.71  (cm²/s; T K, μ solvent viscosity cP, V_A solute molar volume
 * at its normal boiling point cm³/mol).
 * Reference: Hayduk & Minhas, Can. J. Chem. Eng. 60 (1982) 295–299. Validity: n-paraffin solutes C5–C32 in
 * n-paraffin solvents; average error 3.4 % on the original data.
 * @param {object} o
 * @param {number} o.T K @param {number} o.mu oil viscosity, Pa·s
 * @param {number} [o.VA=4.3e-4] wax molar volume, m³/mol (430 cm³/mol is a commonly used value for C20–C25 wax)
 * @returns {number} diffusivity, m²/s
 */
export function diffusivityHaydukMinhas({ T, mu, VA = 4.3e-4 }) {
  need('T', T, 200, 500);
  need('mu', mu, 1e-6, 100);
  need('VA', VA, 5e-5, 5e-3);
  const va = VA * 1e6;
  return 13.3e-8 * T ** 1.47 * (mu * 1e3) ** (10.2 / va - 0.791) / va ** 0.71 * 1e-4;
}

/**
 * Diffusion coefficient by Wilke & Chang (1955): D = 7.4e-8 (φ M_B)^½ T/(μ_B V_A^0.6)  (cm²/s; M_B solvent molar
 * mass g/mol, μ_B cP, V_A cm³/mol, association parameter φ = 1 for hydrocarbons).
 * Reference: Wilke & Chang, AIChE J. 1 (1955) 264–270. Validity: dilute solutions, accuracy about ±10–20 %.
 * @param {object} o
 * @param {number} o.T K @param {number} o.mu solvent viscosity, Pa·s @param {number} o.MB solvent molar mass, kg/mol
 * @param {number} [o.VA=4.3e-4] m³/mol @param {number} [o.phi=1]
 * @returns {number} diffusivity, m²/s
 */
export function diffusivityWilkeChang({ T, mu, MB, VA = 4.3e-4, phi = 1 }) {
  need('T', T, 200, 500);
  need('mu', mu, 1e-6, 100);
  need('MB', MB, 0.01, 2);
  need('VA', VA, 5e-5, 5e-3);
  need('phi', phi, 0.5, 3);
  return 7.4e-8 * Math.sqrt(phi * MB * 1e3) * T / ((mu * 1e3) * (VA * 1e6) ** 0.6) * 1e-4;
}

/**
 * Wax deposition mass flux at the wall by molecular diffusion (Fick's law with the chain rule),
 *   j = ρ_oil D_wo (dC/dT)(dT/dr),
 * C = dissolved wax mass fraction, dT/dr = radial temperature gradient in the oil at the deposit surface
 * (= h_i (T_bulk − T_surface)/k_oil), plus optionally the shear-dispersion term of Burger et al.,
 *   j_sd = k* C_w γ̇   (k* empirical, kg/m²; C_w = precipitated wax volume fraction at the wall; γ̇ wall shear rate).
 * Reference: Burger, Perkins & Striegler, JPT 33(6) (1981) 1075–1086. Validity: deposition only when the surface
 * is below the WAT and colder than the bulk. Later work (Singh et al. 2000) shows shear dispersion to be minor;
 * k* defaults to 0.
 * @param {object} o
 * @param {number} o.rhoOil kg/m³ @param {number} o.Dwo m²/s @param {number} o.dCdT 1/K @param {number} o.dTdr K/m
 * @param {number} [o.kStar=0] shear-dispersion constant, kg/m² @param {number} [o.Cw=0] –
 * @param {number} [o.shearRate=0] 1/s
 * @returns {{ flux: number, molecular: number, shearDispersion: number }} kg/(m²·s)
 */
export function depositionFlux({ rhoOil, Dwo, dCdT, dTdr, kStar = 0, Cw = 0, shearRate = 0 }) {
  need('rhoOil', rhoOil, 100, 2000);
  need('Dwo', Dwo, 0, 1e-6);
  need('dCdT', dCdT, 0, 1);
  need('dTdr', dTdr, 0, 1e7);
  need('kStar', kStar, 0, 1e3);
  need('Cw', Cw, 0, 1);
  need('shearRate', shearRate, 0, 1e7);
  const md = rhoOil * Dwo * dCdT * dTdr;
  const sd = kStar * Cw * shearRate;
  return { flux: md + sd, molecular: md, shearDispersion: sd };
}

/**
 * Matzain semi-empirical deposit growth rate:
 *   dδ/dt = Π1/(1 + Π2) · D_wo (dC/dT)(dT/dr),  Π1 = C1/(1 − C_oil/100),  Π2 = C2 N_SR^C3,
 * with C1 = 15.0, C2 = 0.055, C3 = 1.4. Π1 accounts for oil trapped in the deposit (C_oil = vol % oil) and the
 * empirical enhancement of the diffusion rate, Π2 for shear stripping. Flow-regime dependent shear number:
 *   single-phase and stratified: N_SR = ρ_o v_o δ/μ_o;  bubble and slug: N_SR = ρ_m v_o δ/μ_o;
 *   annular: N_SR = √(ρ_m ρ_o) v_o δ/μ_o     (v_o = oil velocity, ρ_m = mixture density).
 * Here C is the dissolved-wax fraction (volume basis in the original; a mass fraction is used with the ratio
 * ρ_oil/ρ_wax ≈ 1).
 * Reference: Matzain, A., PhD thesis, U. Tulsa (1999); Matzain, Apte, Zhang, Volk, Brill & Creek, J. Energy
 * Resources Technology 124 (2002) 180–186. Validity: fitted to South Pelto crude and Garden Banks condensate flow
 * loops; constants are editable.
 * @param {object} o
 * @param {number} o.Dwo m²/s @param {number} o.dCdT 1/K @param {number} o.dTdr K/m
 * @param {number} o.delta current deposit thickness, m
 * @param {number} o.rhoOil kg/m³ @param {number} o.vOil m/s @param {number} o.muOil Pa·s
 * @param {number} [o.oilFraction=0] oil trapped in the deposit, vol % (0–99)
 * @param {'single'|'stratified'|'bubble'|'slug'|'annular'} [o.regime='single'] @param {number} [o.rhoMix] kg/m³
 * @param {number} [o.C1=15.0] @param {number} [o.C2=0.055] @param {number} [o.C3=1.4]
 * @returns {{ rate: number, Pi1: number, Pi2: number, NSR: number }} growth rate m/s and factors
 */
export function matzainGrowthRate({ Dwo, dCdT, dTdr, delta, rhoOil, vOil, muOil, oilFraction = 0, regime = 'single',
  rhoMix, C1 = 15.0, C2 = 0.055, C3 = 1.4 }) {
  need('Dwo', Dwo, 0, 1e-6);
  need('dCdT', dCdT, 0, 1);
  need('dTdr', dTdr, 0, 1e7);
  need('delta', delta, 0, 1);
  need('rhoOil', rhoOil, 100, 2000);
  need('vOil', vOil, 0, 100);
  need('muOil', muOil, 1e-6, 100);
  need('oilFraction', oilFraction, 0, 99);
  let rhoEff;
  if (regime === 'single' || regime === 'stratified') rhoEff = rhoOil;
  else if (regime === 'bubble' || regime === 'slug') rhoEff = need('rhoMix', rhoMix, 1, 2000);
  else if (regime === 'annular') rhoEff = Math.sqrt(need('rhoMix', rhoMix, 1, 2000) * rhoOil);
  else throw new RangeError(`regime = ${regime} is not one of single, stratified, bubble, slug, annular`);
  const NSR = rhoEff * vOil * delta / muOil;
  const Pi1 = C1 / (1 - oilFraction / 100);
  const Pi2 = C2 * NSR ** C3;
  return { rate: Pi1 / (1 + Pi2) * Dwo * dCdT * dTdr, Pi1, Pi2, NSR };
}

function fChurchill(Re, rr) {
  const A = (2.457 * Math.log(1 / ((7 / Re) ** 0.9 + 0.27 * rr))) ** 16;
  const B = (37530 / Re) ** 16;
  return 8 * ((8 / Re) ** 12 + (A + B) ** -1.5) ** (1 / 12);
}

/**
 * Wax deposit growth along a single-phase oil line with thermal feedback, marched in time (explicit Euler on the
 * thickness) and along the pipe (segments, analytic exponential temperature decay per segment).
 * Per segment and time level:
 *   d = D − 2δ (flow diameter), v = ṁ/(ρ πd²/4), h_i from Nu(Re, Pr) (laminar 3.66 / Gnielinski),
 *   thermal resistances per unit length: film 1/(h_i π d), deposit ln(D/d)/(2π k_dep), wall + insulation +
 *   external 1/(U_wall π D); the insulating deposit lowers the heat loss and raises the surface temperature;
 *   deposit-surface temperature T_s = T_b − q' R_film; deposition when T_s < WAT and T_b > T_s;
 *   dT/dr = h_i (T_b − T_s)/k_oil; growth by
 *     'fick':    dδ/dt = ρ_oil D_wo (dC/dT)(dT/dr) / (ρ_dep F_w)   (F_w = wax mass fraction of the deposit, ageing input),
 *     'matzain': `matzainGrowthRate` (single-phase shear number).
 *   Optional shear stripping for 'fick': dδ/dt −= k_strip τ_w δ (k_strip in 1/(Pa·s), argument).
 * Frictional pressure drop uses the reduced diameter and the deposit roughness (Churchill friction factor).
 * Reference: Burger et al. (1981) molecular-diffusion mechanism; Singh, Venkatesan, Fogler & Nagarajan, AIChE J.
 * 46 (2000) 1059 (thermal feedback and ageing concept); Matzain et al. (2002).
 * Validity: single-phase, quasi-steady hydraulics and heat transfer (deposit growth is slow), constant oil
 * properties except viscosity (may be a function of T), no gelation, no pigging events inside the run.
 * @param {object} o
 * @param {number} o.D clean inner diameter, m @param {number} o.L length, m @param {number} [o.n=50] segments
 * @param {number} o.mDot oil mass rate, kg/s @param {number} o.Tin K @param {number} o.Tamb K
 * @param {number} o.Uwall conductance of wall + coatings + external film referred to πD, W/(m²·K) (no inner film)
 * @param {number} o.rhoOil kg/m³ @param {number|function(number):number} o.muOil Pa·s (or function of T)
 * @param {number} o.cpOil J/(kg·K) @param {number} o.kOil W/(m·K)
 * @param {{ dCdT: function(number):number, WAT: number|null }} o.solubility from `makeSolubilityCurve`
 * @param {number} [o.WAT] wax appearance temperature, K (default solubility.WAT)
 * @param {number|function(number, number):number} [o.Dwo] diffusivity m²/s or function (T, μ); default
 *   Hayduk–Minhas with `VA`
 * @param {number} [o.VA=4.3e-4] m³/mol
 * @param {'fick'|'matzain'} [o.model='fick']
 * @param {number} [o.waxFraction=0.3] wax mass fraction of the deposit F_w ('fick'), – (typical fresh deposits
 *   0.1–0.4, rising with age; editable)
 * @param {number} [o.rhoDeposit=880] kg/m³ @param {number} [o.kDeposit=0.25] W/(m·K) (typical; editable)
 * @param {number} [o.oilFraction=0] trapped oil, vol % ('matzain')
 * @param {number} [o.kStrip=0] shear-stripping constant, 1/(Pa·s) ('fick')
 * @param {number} [o.rough=4.5e-5] clean-pipe roughness, m @param {number} [o.roughDeposit] deposit roughness, m
 * @param {number} o.tEnd s @param {number} o.dt s
 * @param {number} [o.maxThickness] stop when the maximum thickness reaches this value, m
 * @param {number} [o.maxOutput=200]
 * @returns {{ x: number[], t: number[], thickness: number[], thicknessMax: number[], dP: number[],
 *   Tarrival: number[], waxVolume: number[], Tbulk: number[], Tsurface: number[], tStop: number }}
 *   x = segment mid-points m; thickness = final δ(x) m; histories at `t`: max thickness m, frictional ΔP Pa,
 *   arrival temperature K, deposit volume m³; final bulk and surface temperature profiles K
 */
export function waxDepositionLine({ D, L, n = 50, mDot, Tin, Tamb, Uwall, rhoOil, muOil, cpOil, kOil, solubility, WAT,
  Dwo, VA = 4.3e-4, model = 'fick', waxFraction = 0.3, rhoDeposit = 880, kDeposit = 0.25, oilFraction = 0,
  kStrip = 0, rough = 4.5e-5, roughDeposit, tEnd, dt, maxThickness, maxOutput = 200 }) {
  need('D', D, 1e-3, 5);
  need('L', L, 1e-3, 1e7);
  need('n', n, 1, 1e5);
  need('mDot', mDot, 1e-9, 1e5);
  need('Tin', Tin, 200, 600);
  need('Tamb', Tamb, 200, 600);
  need('Uwall', Uwall, 1e-6, 1e5);
  need('rhoOil', rhoOil, 100, 2000);
  need('cpOil', cpOil, 100, 1e4);
  need('kOil', kOil, 0.01, 2);
  need('waxFraction', waxFraction, 1e-3, 1);
  need('rhoDeposit', rhoDeposit, 100, 2000);
  need('kDeposit', kDeposit, 0.01, 5);
  need('kStrip', kStrip, 0, 1);
  need('tEnd', tEnd, 0, 1e10);
  need('dt', dt, 1e-6, 1e9);
  if (model !== 'fick' && model !== 'matzain') throw new RangeError(`model = ${model} is not one of fick, matzain`);
  if (!solubility || typeof solubility.dCdT !== 'function') throw new RangeError('solubility must come from makeSolubilityCurve');
  const wat = WAT ?? solubility.WAT;
  need('WAT', wat, 200, 500);
  const N = Math.round(n);
  const dx = L / N;
  const muOf = typeof muOil === 'function' ? muOil : () => muOil;
  const dwoOf = typeof Dwo === 'function' ? Dwo : Dwo != null ? () => Dwo : (T, mu) => diffusivityHaydukMinhas({ T, mu, VA });
  const eps = roughDeposit ?? rough;
  const delta = new Float64Array(N);
  const Tb = new Float64Array(N);
  const Ts = new Float64Array(N);
  const rate = new Float64Array(N);
  const Rext = 1 / (Uwall * Math.PI * D);

  const sweep = () => {
    let T = Tin;
    let dP = 0;
    for (let i = 0; i < N; i++) {
      const d = D - 2 * delta[i];
      const A = Math.PI * d * d / 4;
      const v = mDot / (rhoOil * A);
      // properties at the segment inlet temperature (deposit growth is slow; one evaluation per segment suffices)
      const mu = need('muOil(T)', muOf(T), 1e-6, 1e3);
      const Re = rhoOil * v * d / mu;
      const Pr = cpOil * mu / kOil;
      const hi = nusseltInternal({ Re, Pr }) * kOil / d;
      const Rfilm = 1 / (hi * Math.PI * d);
      const Rdep = delta[i] > 0 ? Math.log(D / d) / (2 * Math.PI * kDeposit) : 0;
      const Rtot = Rfilm + Rdep + Rext;
      const lambda = mDot * cpOil * Rtot;
      const Tout = Tamb + (T - Tamb) * Math.exp(-dx / lambda);
      const Tm = 0.5 * (T + Tout);
      const q = (Tm - Tamb) / Rtot; // W/m
      const Tsurf = Tm - q * Rfilm;
      Tb[i] = Tm;
      Ts[i] = Tsurf;
      let r = 0;
      if (Tsurf < wat && Tm > Tsurf) {
        const dTdr = hi * (Tm - Tsurf) / kOil;
        const dCdT = Math.max(0, solubility.dCdT(Tsurf));
        const Dw = dwoOf(Tsurf, need('muOil(T)', muOf(Tsurf), 1e-6, 1e3));
        if (model === 'matzain') {
          r = matzainGrowthRate({ Dwo: Dw, dCdT, dTdr, delta: delta[i], rhoOil, vOil: v, muOil: mu, oilFraction }).rate;
        } else {
          r = rhoOil * Dw * dCdT * dTdr / (rhoDeposit * waxFraction);
        }
      }
      const f = fChurchill(Math.max(Re, 1e-6), Math.min((delta[i] > 0 ? eps : rough) / d, 0.1));
      const grad = f * rhoOil * v * v / (2 * d);
      if (model === 'fick' && kStrip > 0) r -= kStrip * (f / 8 * rhoOil * v * v) * delta[i];
      rate[i] = r;
      dP += grad * dx;
      T = Tout;
    }
    return { dP, Tarr: T };
  };

  const steps = Math.max(0, Math.ceil(tEnd / dt - 1e-9));
  const stride = Math.max(1, Math.ceil(steps / maxOutput));
  const out = { t: [], thicknessMax: [], dP: [], Tarrival: [], waxVolume: [] };
  const record = (t, s) => {
    let mx = 0;
    let vol = 0;
    for (let i = 0; i < N; i++) {
      mx = Math.max(mx, delta[i]);
      vol += Math.PI * (D * delta[i] - delta[i] * delta[i]) * dx;
    }
    out.t.push(t);
    out.thicknessMax.push(mx);
    out.dP.push(s.dP);
    out.Tarrival.push(s.Tarr);
    out.waxVolume.push(vol);
    return mx;
  };
  let s = sweep();
  let tStop = tEnd;
  record(0, s);
  for (let k = 1; k <= steps; k++) {
    for (let i = 0; i < N; i++) delta[i] = Math.min(Math.max(0, delta[i] + rate[i] * dt), 0.49 * D);
    s = sweep();
    let mx = 0;
    for (let i = 0; i < N; i++) mx = Math.max(mx, delta[i]);
    const stop = maxThickness != null && mx >= maxThickness;
    if (k % stride === 0 || k === steps || stop) record(k * dt, s);
    if (stop) { tStop = k * dt; break; }
  }
  return { x: Array.from({ length: N }, (_, i) => (i + 0.5) * dx), ...out, thickness: Array.from(delta),
    Tbulk: Array.from(Tb), Tsurface: Array.from(Ts), tStop };
}

/**
 * Pigging interval from a deposition history: the first time at which the maximum deposit thickness reaches
 * `maxThickness` and/or the frictional pressure drop reaches `dPratio` times its clean value (linear interpolation
 * between stored time levels; linear extrapolation of the last growth rate when the limit is not reached).
 * Reference: operating criterion (typical limits 2–4 mm, or a 10–20 % ΔP increase); no correlation implied.
 * @param {object} o
 * @param {number[]} o.t times, s @param {number[]} o.thicknessMax m @param {number[]} [o.dP] Pa
 * @param {number} [o.maxThickness] m @param {number} [o.dPratio] – (> 1)
 * @returns {{ interval: number, governing: string, byThickness: number, byPressure: number, extrapolated: boolean }}
 *   seconds (Infinity when a criterion is absent or never approached)
 */
export function piggingInterval({ t, thicknessMax, dP, maxThickness, dPratio }) {
  if (!t || !thicknessMax || t.length !== thicknessMax.length || t.length < 2) throw new RangeError('t and thicknessMax must be arrays of equal length ≥ 2');
  if (maxThickness == null && dPratio == null) throw new RangeError('give maxThickness and/or dPratio');
  let extrapolated = false;
  const cross = (y, limit) => {
    for (let i = 1; i < y.length; i++) {
      if (y[i] >= limit) return y[i - 1] >= limit ? t[i - 1] : t[i - 1] + (t[i] - t[i - 1]) * (limit - y[i - 1]) / (y[i] - y[i - 1]);
    }
    const k = y.length - 1;
    const slope = (y[k] - y[k - 1]) / (t[k] - t[k - 1]);
    if (!(slope > 0)) return Infinity;
    extrapolated = true;
    return t[k] + (limit - y[k]) / slope;
  };
  const byT = maxThickness != null ? cross(thicknessMax, need('maxThickness', maxThickness, 1e-9, 10)) : Infinity;
  let byP = Infinity;
  if (dPratio != null) {
    if (!dP || dP.length !== t.length) throw new RangeError('dP history is required for the dPratio criterion');
    byP = cross(dP, need('dPratio', dPratio, 1, 1e3) * dP[0]);
  }
  return { interval: Math.min(byT, byP), governing: byT <= byP ? 'thickness' : 'pressure', byThickness: byT,
    byPressure: byP, extrapolated };
}

/**
 * Restart pressure to break a gelled (waxy) oil plug, from a force balance between the pressure difference and
 * the wall yield stress: ΔP = 4 τ_y L/D.
 * Reference: e.g. Davidson, Nguyen, Chang & Rønningsen, J. Non-Newtonian Fluid Mech. 123 (2004) 269 (upper-bound
 * estimate; compressibility and thixotropic breakdown reduce the pressure actually needed).
 * @param {object} o
 * @param {number} o.tauY gel yield stress at the shut-in temperature and ageing time, Pa
 * @param {number} o.L gelled length, m @param {number} o.D inner diameter, m
 * @returns {number} restart pressure difference, Pa
 */
export function gelRestartPressure({ tauY, L, D }) {
  need('tauY', tauY, 0, 1e6);
  need('L', L, 0, 1e7);
  need('D', D, 1e-4, 10);
  return 4 * tauY * L / D;
}

/**
 * Gelling screen: the line is at risk of gelling on shutdown when the lowest fluid temperature falls below the
 * pour point (plus margin). Inputs are laboratory values (ASTM D97/D5853 pour point).
 * @param {object} o
 * @param {number} o.Tmin minimum fluid temperature during shutdown, K @param {number} o.pourPoint K
 * @param {number} [o.margin=5] K
 * @returns {{ gelRisk: boolean, margin: number }} margin = Tmin − pour point, K
 */
export function gellingScreen({ Tmin, pourPoint, margin = 5 }) {
  need('Tmin', Tmin, 100, 600);
  need('pourPoint', pourPoint, 100, 600);
  need('margin', margin, 0, 100);
  return { gelRisk: Tmin < pourPoint + margin, margin: Tmin - pourPoint };
}
