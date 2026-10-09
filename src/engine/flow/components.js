// Pure-component property database, C7+ pseudo-component characterisation and default binary interaction
// parameters for cubic equations of state. SI units: molar mass kg/mol, Tc and Tb in K, Pc in Pa, Vc in m³/mol,
// heat capacity J/(mol·K).

export const R = 8.314462618;

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}

// Columns: molar mass (g/mol), Tc (K), Pc (bar), acentric factor, normal boiling point (K), Vc (cm³/mol),
// ideal-gas Cp/R = a0 + a1 T + a2 T² + a3 T³ + a4 T⁴ (T in K), parachor, PR volume-shift ratio s = c/b.
// Critical constants, acentric factors, boiling points, critical volumes and Cp polynomials: Poling, Prausnitz &
// O'Connell, "The Properties of Gases and Liquids", 5th ed. (2001), Appendix A (Cp fits valid 50–1000 K, 200–1000 K
// for C4+). Other standard tabulations (GPSA, DIPPR) differ in the last digit.
// Parachors: Weinaug & Katz (1943) as tabulated by Danesh, "PVT and Phase Behaviour of Petroleum Reservoir Fluids"
// (1998). PR volume shifts: Jhaveri & Youngren (1988) as tabulated by Whitson & Brulé, "Phase Behavior", SPE
// Monograph 20 (2000), ch. 4. Both are editable defaults (null where no value is tabulated).
const DATA = [
  // id        name                M        Tc      Pc     omega    Tb      Vc     cp a0..a4 (a1 ×1e-3, a2 ×1e-5, a3 ×1e-8, a4 ×1e-11)   Pch   sPR
  ['N2', 'Nitrogen', 28.0134, 126.20, 33.98, 0.0377, 77.35, 90.10, [3.539, -0.261, 0.007, 0.157, -0.099], 41.0, -0.1927],
  ['CO2', 'Carbon dioxide', 44.0095, 304.12, 73.74, 0.2250, 194.65, 94.07, [3.259, 1.356, 1.502, -2.374, 1.056], 78.0, -0.0817],
  ['H2S', 'Hydrogen sulphide', 34.081, 373.40, 89.63, 0.0900, 212.84, 98.00, [4.266, -3.438, 1.319, -1.331, 0.488], 80.1, -0.1288],
  ['C1', 'Methane', 16.043, 190.56, 45.99, 0.0115, 111.66, 98.60, [4.568, -8.975, 3.631, -3.407, 1.091], 77.0, -0.1595],
  ['C2', 'Ethane', 30.070, 305.32, 48.72, 0.0995, 184.55, 145.50, [4.178, -4.427, 5.660, -6.651, 2.487], 108.0, -0.1134],
  ['C3', 'Propane', 44.097, 369.83, 42.48, 0.1523, 231.02, 200.00, [3.847, 5.131, 6.011, -7.893, 3.079], 150.3, -0.0863],
  ['iC4', 'Isobutane', 58.123, 407.85, 36.40, 0.1835, 261.34, 262.70, [3.351, 17.883, 5.477, -8.099, 3.243], 181.5, -0.0844],
  ['nC4', 'n-Butane', 58.123, 425.12, 37.96, 0.2002, 272.66, 255.00, [5.547, 5.536, 8.057, -10.571, 4.134], 189.9, -0.0675],
  ['iC5', 'Isopentane', 72.150, 460.39, 33.81, 0.2275, 300.99, 306.00, [1.959, 38.191, 2.434, -5.175, 2.165], 225.0, -0.0608],
  ['nC5', 'n-Pentane', 72.150, 469.70, 33.70, 0.2515, 309.22, 311.00, [7.554, -0.368, 11.846, -14.939, 5.753], 231.5, -0.0390],
  ['nC6', 'n-Hexane', 86.177, 507.60, 30.25, 0.3013, 341.88, 368.00, [8.831, -0.166, 14.302, -18.314, 7.124], 271.0, -0.0080],
  ['nC7', 'n-Heptane', 100.204, 540.20, 27.40, 0.3495, 371.57, 428.00, [9.634, 4.156, 15.494, -20.066, 7.770], 312.5, 0.0033],
  ['nC8', 'n-Octane', 114.231, 568.70, 24.90, 0.3996, 398.82, 492.00, [10.824, 4.983, 17.751, -23.137, 8.980], 351.5, 0.0314],
  ['nC9', 'n-Nonane', 128.258, 594.60, 22.90, 0.4435, 423.97, 555.00, [12.152, 4.575, 20.416, -26.777, 10.465], 393.0, 0.0408],
  ['nC10', 'n-Decane', 142.285, 617.70, 21.10, 0.4923, 447.30, 624.00, [13.467, 4.139, 23.127, -30.477, 11.970], 433.5, 0.0655],
  ['H2O', 'Water', 18.01528, 647.14, 220.64, 0.3443, 373.15, 55.95, [4.395, -4.186, 1.405, -1.564, 0.632], null, null],
  ['methanol', 'Methanol', 32.042, 512.64, 80.97, 0.5650, 337.69, 118.00, [4.714, -6.986, 4.211, -4.443, 1.535], null, null],
  // MEG critical constants are estimates in every tabulation (decomposes below Tc); no vetted Cp polynomial is
  // carried here — supply `cp` through makeComponent() if ideal-gas enthalpies of MEG-bearing phases are needed.
  ['MEG', 'Monoethylene glycol', 62.068, 720.00, 82.00, 0.5068, 470.45, 191.00, null, null, null],
  ['H2', 'Hydrogen', 2.01588, 33.19, 13.13, -0.2160, 20.38, 64.20, [2.883, 3.681, -0.772, 0.692, -0.213], null, null],
  ['O2', 'Oxygen', 31.9988, 154.58, 50.43, 0.0222, 90.19, 73.37, [3.630, -1.794, 0.658, -0.601, 0.179], null, null],
  // Argon: monatomic ideal gas, Cp = 5R/2 exactly.
  ['Ar', 'Argon', 39.948, 150.86, 48.98, -0.0020, 87.27, 74.57, [2.5, 0, 0, 0, 0], null, null],
];

function build(row) {
  const [id, name, M, Tc, Pc, omega, Tb, Vc, cp, parachor, shiftPR] = row;
  return Object.freeze({
    id, name, M: M * 1e-3, Tc, Pc: Pc * 1e5, omega, Tb, Vc: Vc * 1e-6,
    cp: cp ? Object.freeze([cp[0], cp[1] * 1e-3, cp[2] * 1e-5, cp[3] * 1e-8, cp[4] * 1e-11]) : null,
    parachor, shiftPR, pseudo: false,
  });
}

/**
 * Pure-component database keyed by id (N2, CO2, H2S, C1, C2, C3, iC4, nC4, iC5, nC5, nC6, nC7, nC8, nC9, nC10,
 * H2O, methanol, MEG, H2, O2, Ar). Each record: { id, name, M kg/mol, Tc K, Pc Pa, omega –, Tb K, Vc m³/mol,
 * cp: [a0..a4] for Cp_ig/R = Σ a_k T^k (or null), parachor (dyn^¼·cm^11/4/mol, or null), shiftPR (–, or null) }.
 * Reference: Poling, Prausnitz & O'Connell (2001) Appendix A; see the comment above the table for the rest.
 */
export const COMPONENTS = Object.freeze(Object.fromEntries(DATA.map((r) => [r[0], build(r)])));

/** Ids of the database components in table order. */
export const COMPONENT_IDS = Object.freeze(DATA.map((r) => r[0]));

/**
 * Look up a database component.
 * @param {string} id component id (see `COMPONENT_IDS`)
 * @returns {object} frozen component record (see `COMPONENTS`)
 */
export function getComponent(id) {
  const c = COMPONENTS[id];
  if (!c) throw new RangeError(`component id = ${id} is not in the database`);
  return c;
}

/**
 * Create a component record, either from scratch or as an edited copy of a database component.
 * Reference: none (data container). Required fields when `base` is omitted: id, M, Tc, Pc, omega.
 * @param {object} o
 * @param {string} [o.base] database id to copy
 * @param {string} [o.id] @param {string} [o.name]
 * @param {number} [o.M] molar mass, kg/mol
 * @param {number} [o.Tc] critical temperature, K
 * @param {number} [o.Pc] critical pressure, Pa
 * @param {number} [o.omega] acentric factor, –
 * @param {number} [o.Tb] normal boiling point, K
 * @param {number} [o.Vc] critical volume, m³/mol
 * @param {number[]} [o.cp] ideal-gas Cp/R polynomial coefficients a0…a4 (T in K)
 * @param {number} [o.parachor] @param {number} [o.shiftPR]
 * @returns {object} component record
 */
export function makeComponent({ base, ...over }) {
  const b = base ? getComponent(base) : { name: over.id, Tb: null, Vc: null, cp: null, parachor: null, shiftPR: null, pseudo: false };
  const c = { ...b, ...over };
  if (typeof c.id !== 'string' || !c.id) throw new RangeError('id must be a non-empty string');
  need('M', c.M, 1e-3, 2);
  need('Tc', c.Tc, 5, 2000);
  need('Pc', c.Pc, 1e4, 1e8);
  need('omega', c.omega, -0.5, 2.5);
  return Object.freeze(c);
}

/**
 * Ideal-gas molar heat capacity Cp_ig = R Σ a_k T^k.
 * Reference: Poling, Prausnitz & O'Connell (2001) Appendix A polynomial. Validity: 50–1000 K (200–1000 K for C4+).
 * For C7+ pseudo-components the n-decane polynomial is scaled by molar mass (equal mass-specific Cp), a common
 * engineering approximation good to roughly ±5 % for paraffinic cuts.
 * @param {object|string} comp component record or database id
 * @param {number} T temperature, K
 * @returns {number} J/(mol·K)
 */
export function idealGasCp(comp, T) {
  const c = typeof comp === 'string' ? getComponent(comp) : comp;
  need('T', T, 20, 1500);
  if (!c.cp) throw new RangeError(`component ${c.id} has no ideal-gas Cp polynomial; supply cp via makeComponent`);
  const a = c.cp;
  return R * (a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4]))));
}

/**
 * Ideal-gas enthalpy and entropy integrals from a reference temperature:
 * h = ∫ Cp dT, s = ∫ Cp/T dT (analytic integration of the polynomial).
 * Reference: as `idealGasCp`.
 * @param {object|string} comp component record or id
 * @param {number} T temperature, K
 * @param {number} [Tref=298.15] reference temperature, K
 * @returns {{ h: number, s: number }} J/mol and J/(mol·K)
 */
export function idealGasEnthalpyEntropy(comp, T, Tref = 298.15) {
  const c = typeof comp === 'string' ? getComponent(comp) : comp;
  need('T', T, 20, 1500);
  if (!c.cp) throw new RangeError(`component ${c.id} has no ideal-gas Cp polynomial; supply cp via makeComponent`);
  const a = c.cp;
  const H = (t) => t * (a[0] + t * (a[1] / 2 + t * (a[2] / 3 + t * (a[3] / 4 + t * a[4] / 5))));
  const S = (t) => a[0] * Math.log(t) + t * (a[1] + t * (a[2] / 2 + t * (a[3] / 3 + t * a[4] / 4)));
  return { h: R * (H(T) - H(Tref)), s: R * (S(T) - S(Tref)) };
}

/**
 * Normal boiling point of a petroleum fraction from molar mass and specific gravity (Søreide 1989):
 * Tb(°R) = 1928.3 − 1.695e5 M^−0.03522 γ^3.266 exp(−4.922e-3 M − 4.7685 γ + 3.462e-3 M γ), M in g/mol.
 * Reference: Søreide, I. (1989) Dr.ing. thesis, NTH Trondheim, as given by Whitson & Brulé, "Phase Behavior", SPE Monograph 20 (2000) ch. 5.
 * Validity: C7+ fractions, M 90–300 g/mol (used with care to ~500), γ 0.70–0.95.
 * @param {number} M molar mass, kg/mol
 * @param {number} sg specific gravity (water = 1, 60 °F), –
 * @returns {number} Tb, K
 */
export function boilingPointSoreide(M, sg) {
  need('M', M, 0.07, 0.7);
  need('sg', sg, 0.6, 1.1);
  const Mg = M * 1e3;
  const TbR = 1928.3 - 1.695e5 * Mg ** -0.03522 * sg ** 3.266
    * Math.exp(-4.922e-3 * Mg - 4.7685 * sg + 3.462e-3 * Mg * sg);
  return TbR / 1.8;
}

/**
 * Kesler–Lee (1976) correlations for a petroleum fraction from boiling point and specific gravity
 * (Tb in °R, pressures psia):
 *   Tc = 341.7 + 811.1 γ + (0.4244 + 0.1174 γ) Tb + (0.4669 − 3.26238 γ)·1e5/Tb
 *   ln Pc = 8.3634 − 0.0566/γ − (0.24244 + 2.2898/γ + 0.11857/γ²)·1e-3 Tb
 *           + (1.4685 + 3.648/γ + 0.47227/γ²)·1e-7 Tb² − (0.42019 + 1.6977/γ²)·1e-10 Tb³
 *   M = −12272.6 + 9486.4 γ + (4.6523 − 3.3287 γ) Tb + (1 − 0.77084 γ − 0.02058 γ²)(1.3437 − 720.79/Tb)·1e7/Tb
 *       + (1 − 0.80882 γ + 0.02226 γ²)(1.8828 − 181.98/Tb)·1e12/Tb³
 *   ω (Tbr > 0.8) = −7.904 + 0.1352 Kw − 0.007465 Kw² + 8.359 Tbr + (1.408 − 0.01063 Kw)/Tbr
 *   ω (Tbr ≤ 0.8) = [−ln(Pc/14.7) − 5.92714 + 6.09648/Tbr + 1.28862 ln Tbr − 0.169347 Tbr⁶]
 *                   / [15.2518 − 15.6875/Tbr − 13.4721 ln Tbr + 0.43577 Tbr⁶]   (Lee–Kesler)
 * with the Watson factor Kw = Tb^(1/3)/γ.
 * Reference: Kesler & Lee, Hydrocarbon Processing 55(3) (1976) 153–158; Lee & Kesler, AIChE J. 21 (1975).
 * Validity: Tb up to about 920 K (1200 °F), M 60–650 g/mol.
 * @param {object} o
 * @param {number} o.Tb normal boiling point, K
 * @param {number} o.sg specific gravity, –
 * @returns {{ Tc: number, Pc: number, omega: number, M: number, Kw: number }} K, Pa, –, kg/mol, °R^(1/3)
 */
export function keslerLee({ Tb, sg }) {
  need('Tb', Tb, 250, 1100);
  need('sg', sg, 0.55, 1.2);
  const t = Tb * 1.8;
  const g = sg;
  const Tc = 341.7 + 811.1 * g + (0.4244 + 0.1174 * g) * t + (0.4669 - 3.26238 * g) * 1e5 / t;
  const lnPc = 8.3634 - 0.0566 / g - (0.24244 + 2.2898 / g + 0.11857 / (g * g)) * 1e-3 * t
    + (1.4685 + 3.648 / g + 0.47227 / (g * g)) * 1e-7 * t * t - (0.42019 + 1.6977 / (g * g)) * 1e-10 * t ** 3;
  const Pc = Math.exp(lnPc);
  const M = -12272.6 + 9486.4 * g + (4.6523 - 3.3287 * g) * t
    + (1 - 0.77084 * g - 0.02058 * g * g) * (1.3437 - 720.79 / t) * 1e7 / t
    + (1 - 0.80882 * g + 0.02226 * g * g) * (1.8828 - 181.98 / t) * 1e12 / t ** 3;
  const Kw = Math.cbrt(t) / g;
  const Tbr = t / Tc;
  let omega;
  if (Tbr > 0.8) {
    omega = -7.904 + 0.1352 * Kw - 0.007465 * Kw * Kw + 8.359 * Tbr + (1.408 - 0.01063 * Kw) / Tbr;
  } else {
    omega = (-Math.log(Pc / 14.7) - 5.92714 + 6.09648 / Tbr + 1.28862 * Math.log(Tbr) - 0.169347 * Tbr ** 6)
      / (15.2518 - 15.6875 / Tbr - 13.4721 * Math.log(Tbr) + 0.43577 * Tbr ** 6);
  }
  return { Tc: Tc / 1.8, Pc: Pc * 6894.757293168, omega, M: M * 1e-3, Kw };
}

/**
 * Riazi–Daubert (1980) two-parameter correlations θ = a Tb^b γ^c (Tb °R):
 *   Tc(°R) = 24.2787 Tb^0.58848 γ^0.3596, Pc(psia) = 3.12281e9 Tb^−2.3125 γ^2.3201,
 *   M = 4.5673e-5 Tb^2.1962 γ^−1.0164.
 * Reference: Riazi & Daubert, Hydrocarbon Processing 59(3) (1980) 115–116.
 * Validity: Tb 311–728 K (100–850 °F), M 70–300 g/mol.
 * @param {object} o
 * @param {number} o.Tb K @param {number} o.sg –
 * @returns {{ Tc: number, Pc: number, M: number }} K, Pa, kg/mol
 */
export function riaziDaubert({ Tb, sg }) {
  need('Tb', Tb, 250, 900);
  need('sg', sg, 0.55, 1.2);
  const t = Tb * 1.8;
  return {
    Tc: 24.2787 * t ** 0.58848 * sg ** 0.3596 / 1.8,
    Pc: 3.12281e9 * t ** -2.3125 * sg ** 2.3201 * 6894.757293168,
    M: 4.5673e-5 * t ** 2.1962 * sg ** -1.0164 * 1e-3,
  };
}

/**
 * Characterise a C7+ (or any heavy) pseudo-component from molar mass and specific gravity.
 * Boiling point by Søreide (unless given), then Tc, Pc and ω by Kesler–Lee (`method = 'keslerLee'`) or Tc, Pc by
 * Riazi–Daubert with the Lee–Kesler/Kesler–Lee ω (`method = 'riaziDaubert'`). Critical volume by the
 * Lohrenz–Bray–Clark (1964) expression Vc(ft³/lbmol) = 21.573 + 0.015122 M − 27.656 γ + 0.070615 M γ.
 * Parachor by the Weinaug–Katz linear relation for heavy fractions, Pch = 59.3 + 2.34 M (Ahmed, "Equations of
 * State and PVT Analysis"). Ideal-gas Cp: n-decane polynomial scaled by molar mass. Volume shift: none (null).
 * References: as `boilingPointSoreide`, `keslerLee`, `riaziDaubert`; Lohrenz, Bray & Clark, JPT 16(10) (1964).
 * Validity: M 0.09–0.5 kg/mol, γ 0.70–1.0.
 * @param {object} o
 * @param {number} o.M molar mass, kg/mol
 * @param {number} o.sg specific gravity, –
 * @param {number} [o.Tb] normal boiling point, K (Søreide estimate when omitted)
 * @param {string} [o.id='C7+']
 * @param {'keslerLee'|'riaziDaubert'} [o.method='keslerLee']
 * @returns {object} component record (pseudo: true) including sg and Kw
 */
export function characterizePlusFraction({ M, sg, Tb, id = 'C7+', method = 'keslerLee' }) {
  need('M', M, 0.07, 0.7);
  need('sg', sg, 0.6, 1.1);
  const tb = Tb ?? boilingPointSoreide(M, sg);
  const kl = keslerLee({ Tb: tb, sg });
  let { Tc, Pc } = kl;
  if (method === 'riaziDaubert') {
    ({ Tc, Pc } = riaziDaubert({ Tb: tb, sg }));
  } else if (method !== 'keslerLee') {
    throw new RangeError(`method = ${method} is not one of keslerLee, riaziDaubert`);
  }
  const Mg = M * 1e3;
  const VcFt3 = 21.573 + 0.015122 * Mg - 27.656 * sg + 0.070615 * Mg * sg;
  const c10 = COMPONENTS.nC10;
  const scale = M / c10.M;
  return Object.freeze({
    id, name: id, M, Tc, Pc, omega: kl.omega, Tb: tb,
    Vc: VcFt3 * 0.028316846592 / 453.59237, // ft³/lbmol → m³/mol
    cp: Object.freeze(c10.cp.map((a) => a * scale)),
    parachor: 59.3 + 2.34 * Mg, shiftPR: null, pseudo: true, sg, Kw: kl.Kw,
  });
}

// Default binary interaction parameters between non-hydrocarbons and hydrocarbons (and among non-hydrocarbons)
// for the Peng–Robinson EOS: typical values as recommended by Whitson & Brulé (2000), SPE Monograph 20, ch. 4. They are starting values to be tuned against data; all hydrocarbon–hydrocarbon pairs default to 0.
const HC_ORDER = ['C1', 'C2', 'C3', 'iC4', 'nC4', 'iC5', 'nC5', 'nC6'];
const KIJ_NONHC = {
  N2: { CO2: 0.000, H2S: 0.130, hc: [0.025, 0.010, 0.090, 0.095, 0.095, 0.100, 0.110, 0.110], heavy: 0.110 },
  CO2: { H2S: 0.135, hc: [0.105, 0.130, 0.125, 0.120, 0.115, 0.115, 0.115, 0.115], heavy: 0.115 },
  H2S: { hc: [0.070, 0.085, 0.080, 0.075, 0.075, 0.070, 0.070, 0.055], heavy: 0.050 },
};
const HEAVY = new Set(['nC7', 'nC8', 'nC9', 'nC10']);

/**
 * Default binary interaction parameter k_ij for a pair of components (symmetric; 0 when not tabulated).
 * Non-zero defaults exist only for N2, CO2 and H2S against each other and against hydrocarbons (C7 and heavier,
 * including pseudo-components, use the "C7+" column).
 * Reference: Whitson & Brulé (2000) SPE Monograph 20, ch. 4 (Peng–Robinson values). The same numbers are used
 * as a first estimate for SRK. Editable: pass your own matrix to the EOS mixture constructor.
 * @param {object|string} a component record or id
 * @param {object|string} b component record or id
 * @returns {number} k_ij, –
 */
export function defaultKij(a, b) {
  const ca = typeof a === 'string' ? getComponent(a) : a;
  const cb = typeof b === 'string' ? getComponent(b) : b;
  if (ca.id === cb.id) return 0;
  const one = (x, y) => {
    const t = KIJ_NONHC[x.id];
    if (!t) return null;
    if (t[y.id] != null) return t[y.id];
    const k = HC_ORDER.indexOf(y.id);
    if (k >= 0) return t.hc[k];
    if (HEAVY.has(y.id) || y.pseudo) return t.heavy;
    return null;
  };
  return one(ca, cb) ?? one(cb, ca) ?? 0;
}

/**
 * Full symmetric k_ij matrix of defaults for a component list.
 * Reference: as `defaultKij`.
 * @param {(object|string)[]} comps component records or ids
 * @returns {number[][]} n×n matrix, –
 */
export function defaultKijMatrix(comps) {
  const n = comps.length;
  const k = [];
  for (let i = 0; i < n; i++) {
    k.push(new Array(n).fill(0));
    for (let j = 0; j < i; j++) {
      k[i][j] = defaultKij(comps[i], comps[j]);
      k[j][i] = k[i][j];
    }
  }
  return k;
}

/**
 * Mixture molar mass Σ z_i M_i. Reference: definition.
 * @param {(object|string)[]} comps component records or ids
 * @param {number[]} z mole fractions, –
 * @returns {number} kg/mol
 */
export function mixtureMolarMass(comps, z) {
  if (comps.length !== z.length) throw new RangeError('z must have one entry per component');
  let M = 0;
  for (let i = 0; i < z.length; i++) M += z[i] * (typeof comps[i] === 'string' ? getComponent(comps[i]) : comps[i]).M;
  return M;
}
