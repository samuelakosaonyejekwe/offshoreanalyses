// Seawater and produced-water (brine) thermophysical properties.
// SI units throughout: T in K, P absolute in Pa, salinity S in g/kg (practical salinity ≈ g/kg for seawater).

const G0 = 9.80665;
const ATM = 101325;

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}

/**
 * Density of seawater by the UNESCO 1980 international equation of state (EOS-80).
 * Reference: Millero & Poisson (1981) one-atmosphere equation and Millero et al. (1980) secant bulk modulus, as
 * collected in UNESCO Technical Papers in Marine Science No. 44 (Fofonoff & Millard 1983).
 * Check values (UNESCO 44): S = 35, t = 5 °C, p = 0 → 1027.67547 kg/m³; S = 35, t = 5 °C, p = 1000 bar → 1069.48914.
 * Validity: S 0–42, t −2…40 °C, gauge pressure 0–1000 bar. Stated precision of the fit ≈ 0.009 kg/m³.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=35] practical salinity, g/kg
 * @param {number} [o.P=101325] absolute pressure, Pa
 * @returns {number} density, kg/m³
 */
export function seawaterDensity({ T, S = 35, P = ATM }) {
  need('T', T, 271.15, 313.16);
  need('S', S, 0, 42);
  need('P', P, 0, 1.002e8);
  const t = T - 273.15;
  const p = (P - ATM) / 1e5; // sea pressure, bar
  const rho0 = seawaterDensityAtm(t, S);
  if (p === 0) return rho0;
  const K = secantBulkModulusBar(t, S, p);
  return rho0 / (1 - p / K);
}

function seawaterDensityAtm(t, S) {
  const rw = 999.842594 + t * (6.793952e-2 + t * (-9.095290e-3 + t * (1.001685e-4 + t * (-1.120083e-6
    + t * 6.536332e-9))));
  const A = 8.24493e-1 + t * (-4.0899e-3 + t * (7.6438e-5 + t * (-8.2467e-7 + t * 5.3875e-9)));
  const B = -5.72466e-3 + t * (1.0227e-4 - t * 1.6546e-6);
  const C = 4.8314e-4;
  return rw + A * S + B * S * Math.sqrt(S) + C * S * S;
}

function secantBulkModulusBar(t, S, p) {
  const s15 = S * Math.sqrt(S);
  const Kw = 19652.21 + t * (148.4206 + t * (-2.327105 + t * (1.360477e-2 - t * 5.155288e-5)));
  const K0 = Kw + (54.6746 + t * (-0.603459 + t * (1.09987e-2 - t * 6.1670e-5))) * S
    + (7.944e-2 + t * (1.6483e-2 - t * 5.3009e-4)) * s15;
  const Aw = 3.239908 + t * (1.43713e-3 + t * (1.16092e-4 - t * 5.77905e-7));
  const A = Aw + (2.2838e-3 + t * (-1.0981e-5 - t * 1.6078e-6)) * S + 1.91075e-4 * s15;
  const Bw = 8.50935e-5 + t * (-6.12293e-6 + t * 5.2787e-8);
  const B = Bw + (-9.9348e-7 + t * (2.0816e-8 + t * 9.1697e-10)) * S;
  return K0 + p * (A + p * B);
}

/**
 * Secant bulk modulus of seawater K(S, t, p) of EOS-80, such that ρ(p) = ρ(0)/(1 − p/K).
 * Reference: UNESCO Technical Papers in Marine Science No. 44. Validity as `seawaterDensity`.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=35] salinity, g/kg
 * @param {number} [o.P=101325] absolute pressure, Pa
 * @returns {number} secant bulk modulus, Pa
 */
export function seawaterSecantBulkModulus({ T, S = 35, P = ATM }) {
  need('T', T, 271.15, 313.16);
  need('S', S, 0, 42);
  need('P', P, 0, 1.002e8);
  return secantBulkModulusBar(T - 273.15, S, (P - ATM) / 1e5) * 1e5;
}

/**
 * Density of produced water / oilfield brine by the McCain correlations:
 * ρ(sc) = 62.368 + 0.438603 s + 1.60074e-3 s² lb/ft³ (s = dissolved solids, wt %), and the formation volume factor
 * Bw = (1 + ΔVwp)(1 + ΔVwT) with
 * ΔVwT = −1.0001e-2 + 1.33391e-4 T + 5.50654e-7 T²,
 * ΔVwp = −1.95301e-9 p T − 1.72834e-13 p² T − 3.58922e-7 p − 2.25341e-10 p²  (T °F, p psia).
 * Reference: McCain, W.D. (1991) "Reservoir-fluid property correlations — state of the art", SPE Reservoir
 * Engineering 6(2). Validity: up to about 400 K (260 °F), 34.5 MPa (5000 psia), 0–26 wt % solids; density at
 * standard conditions within about 0.1 %, Bw within about 2 %. Dissolved gas is not included.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=0] dissolved solids, g/kg of brine
 * @param {number} [o.P=101325] absolute pressure, Pa
 * @returns {{ rho: number, rhoSc: number, Bw: number }} density kg/m³, density at standard conditions kg/m³,
 *   formation volume factor m³/Sm³
 */
export function brineDensity({ T, S = 0, P = ATM }) {
  need('T', T, 273.15, 450);
  need('S', S, 0, 300);
  need('P', P, 0, 1.5e8);
  const s = S / 10;
  const TF = T * 1.8 - 459.67;
  const p = P / 6894.757293168;
  const rhoSc = (62.368 + 0.438603 * s + 1.60074e-3 * s * s) * 16.01846337;
  const dVT = -1.0001e-2 + 1.33391e-4 * TF + 5.50654e-7 * TF * TF;
  const dVp = -1.95301e-9 * p * TF - 1.72834e-13 * p * p * TF - 3.58922e-7 * p - 2.25341e-10 * p * p;
  const Bw = (1 + dVp) * (1 + dVT);
  return { rho: rhoSc / Bw, rhoSc, Bw };
}

/**
 * Water density dispatcher: EOS-80 inside its validity range (T ≤ 40 °C, S ≤ 42 g/kg, P ≤ 100 MPa), McCain brine
 * correlation otherwise. Use `method` to force one model (the two differ by up to ~0.2 % where they overlap, so
 * a marching solver should force a single method).
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=35] salinity / dissolved solids, g/kg
 * @param {number} [o.P=101325] absolute pressure, Pa
 * @param {'auto'|'eos80'|'mccain'} [o.method='auto']
 * @returns {number} density, kg/m³
 */
export function waterDensity({ T, S = 35, P = ATM, method = 'auto' }) {
  if (method === 'eos80') return seawaterDensity({ T, S, P });
  if (method === 'mccain') return brineDensity({ T, S, P }).rho;
  if (method !== 'auto') throw new RangeError(`method = ${method} is not one of auto, eos80, mccain`);
  if (T >= 271.15 && T <= 313.15 && S <= 42 && P <= 1.0e8) return seawaterDensity({ T, S, P });
  return brineDensity({ T, S, P }).rho;
}

/**
 * Dynamic viscosity of pure water at atmospheric pressure,
 * μw = 4.2844e-5 + 1/(0.157 (t + 64.993)² − 91.296) Pa·s (t in °C).
 * Reference: Sharqawy, Lienhard & Zubair (2010), Desalination and Water Treatment 16, 354–380 (fit to IAPWS 2008).
 * Validity: 0–180 °C, stated accuracy ±0.05 %.
 * @param {number} T temperature, K
 * @returns {number} viscosity, Pa·s
 */
export function waterViscosity(T) {
  need('T', T, 273.15, 453.15);
  const t = T - 273.15;
  return 4.2844e-5 + 1 / (0.157 * (t + 64.993) ** 2 - 91.296);
}

/**
 * Dynamic viscosity of seawater, μsw = μw (1 + A s + B s²) with s in kg/kg,
 * A = 1.541 + 1.998e-2 t − 9.52e-5 t², B = 7.974 − 7.561e-2 t + 4.724e-4 t².
 * Reference: Sharqawy, Lienhard & Zubair (2010). Validity: 0–180 °C, S 0–150 g/kg, accuracy ±1.5 %.
 * Pressure dependence is neglected (below 1 % to 100 MPa at these temperatures).
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=35] salinity, g/kg
 * @returns {number} viscosity, Pa·s
 */
export function seawaterViscosity({ T, S = 35 }) {
  need('S', S, 0, 150);
  const t = T - 273.15;
  const s = S / 1000;
  const A = 1.541 + 1.998e-2 * t - 9.52e-5 * t * t;
  const B = 7.974 - 7.561e-2 * t + 4.724e-4 * t * t;
  return waterViscosity(T) * (1 + A * s + B * s * s);
}

/**
 * Viscosity of oilfield brine with pressure correction (McCain 1991):
 * μw1 = A T^B (cP, T °F), A = 109.574 − 8.40564 s + 0.313314 s² + 8.72213e-3 s³,
 * B = −1.12166 + 2.63951e-2 s − 6.79461e-4 s² − 5.47119e-5 s³ + 1.55586e-6 s⁴ (s wt %),
 * μw/μw1 = 0.9994 + 4.0295e-5 p + 3.1062e-9 p² (p psia).
 * Reference: McCain (1991), SPE Reservoir Engineering 6(2). Validity: 311–478 K (100–400 °F), 0–26 wt %,
 * pressure to 69 MPa for 303–348 K; stated accuracy about 5 % (7 % for the pressure factor above 69 MPa).
 * Below 311 K prefer `seawaterViscosity`.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=0] dissolved solids, g/kg
 * @param {number} [o.P=101325] absolute pressure, Pa
 * @returns {number} viscosity, Pa·s
 */
export function brineViscosity({ T, S = 0, P = ATM }) {
  need('T', T, 273.15, 480);
  need('S', S, 0, 260);
  need('P', P, 0, 1.5e8);
  const s = S / 10;
  const TF = T * 1.8 - 459.67;
  const p = P / 6894.757293168;
  const A = 109.574 - 8.40564 * s + 0.313314 * s * s + 8.72213e-3 * s ** 3;
  const B = -1.12166 + 2.63951e-2 * s - 6.79461e-4 * s * s - 5.47119e-5 * s ** 3 + 1.55586e-6 * s ** 4;
  const mu1 = A * TF ** B;
  return mu1 * (0.9994 + 4.0295e-5 * p + 3.1062e-9 * p * p) * 1e-3;
}

/**
 * Thermal conductivity of seawater,
 * log10 k = log10(240 + 0.0002 S) + 0.434 (2.3 − (343.5 + 0.037 S)/T)(1 − T/(647 + 0.03 S))^0.333, k in mW/(m·K).
 * Reference: Jamieson & Tudhope (1970) as recommended by Sharqawy et al. (2010).
 * Validity: 0–180 °C, S 0–160 g/kg, accuracy ±3 %.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=35] salinity, g/kg
 * @returns {number} thermal conductivity, W/(m·K)
 */
export function seawaterConductivity({ T, S = 35 }) {
  need('T', T, 273.15, 453.15);
  need('S', S, 0, 160);
  const lg = Math.log10(240 + 0.0002 * S)
    + 0.434 * (2.3 - (343.5 + 0.037 * S) / T) * (1 - T / (647 + 0.03 * S)) ** 0.333;
  return 10 ** lg * 1e-3;
}

/**
 * Specific heat capacity of seawater at constant (atmospheric) pressure, cp = A + B T + C T² + D T³ (kJ/kg/K, T K)
 * with A = 5.328 − 9.76e-2 S + 4.04e-4 S², B = −6.913e-3 + 7.351e-4 S − 3.15e-6 S²,
 * C = 9.6e-6 − 1.927e-6 S + 8.23e-9 S², D = 2.5e-9 + 1.666e-9 S − 7.125e-12 S².
 * Reference: Jamieson, Tudhope, Morris & Cartwright (1969) as given by Sharqawy et al. (2010).
 * Validity: 0–180 °C, S 0–180 g/kg, accuracy ±0.28 %.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=35] salinity, g/kg
 * @returns {number} heat capacity, J/(kg·K)
 */
export function seawaterHeatCapacity({ T, S = 35 }) {
  need('T', T, 273.15, 453.15);
  need('S', S, 0, 180);
  const A = 5.328 - 9.76e-2 * S + 4.04e-4 * S * S;
  const B = -6.913e-3 + 7.351e-4 * S - 3.15e-6 * S * S;
  const C = 9.6e-6 - 1.927e-6 * S + 8.23e-9 * S * S;
  const D = 2.5e-9 + 1.666e-9 * S - 7.125e-12 * S * S;
  return (A + T * (B + T * (C + T * D))) * 1e3;
}

/**
 * Speed of sound in seawater (nine-term equation):
 * c = 1448.96 + 4.591 t − 5.304e-2 t² + 2.374e-4 t³ + 1.340 (S − 35) + 1.630e-2 D + 1.675e-7 D²
 *     − 1.025e-2 t (S − 35) − 7.139e-13 t D³.
 * Reference: Mackenzie, K.V. (1981), J. Acoust. Soc. Am. 70(3), 807–812. Check value: t = 25 °C, S = 35,
 * D = 1000 m → 1550.744 m/s. Validity: t −2…30 °C, S 25–40 (applied here for 0–40), depth 0–8000 m; standard
 * error 0.07 m/s.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=35] salinity, g/kg
 * @param {number} [o.depth=0] depth below surface, m
 * @returns {number} sound speed, m/s
 */
export function seawaterSoundSpeed({ T, S = 35, depth = 0 }) {
  need('T', T, 271.15, 303.16);
  need('S', S, 0, 40);
  need('depth', depth, 0, 8000);
  const t = T - 273.15;
  const D = depth;
  return 1448.96 + 4.591 * t - 5.304e-2 * t * t + 2.374e-4 * t ** 3 + 1.340 * (S - 35) + 1.630e-2 * D
    + 1.675e-7 * D * D - 1.025e-2 * t * (S - 35) - 7.139e-13 * t * D ** 3;
}

const IF97_N = [0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
  -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
  0.65017534844798e3];

/**
 * Saturation (vapour) pressure of water, optionally lowered by dissolved salt.
 * `method = 'if97'`: IAPWS Industrial Formulation 1997, Region 4 saturation-pressure equation (Eq. 30), valid
 * 273.15–647.096 K; check value 300 K → 3536.58941 Pa.
 * `method = 'antoine'`: log10 p(mmHg) = 8.07131 − 1730.63/(233.426 + t °C), valid 1–100 °C (accuracy ≈ 0.3 %).
 * Salinity correction: p_sw = p_w / (1 + 0.57357 S/(1000 − S)) (Raoult-type relation, Sharqawy et al. 2010;
 * 0–240 g/kg, ±0.1 %).
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} [o.S=0] salinity, g/kg
 * @param {'if97'|'antoine'} [o.method='if97']
 * @returns {number} vapour pressure, Pa
 */
export function waterVapourPressure({ T, S = 0, method = 'if97' }) {
  need('S', S, 0, 240);
  let p;
  if (method === 'antoine') {
    need('T', T, 273.15, 373.16);
    p = 10 ** (8.07131 - 1730.63 / (233.426 + (T - 273.15))) * 133.322368;
  } else if (method === 'if97') {
    need('T', T, 273.15, 647.096);
    const n = IF97_N;
    const th = T + n[8] / (T - n[9]);
    const A = th * th + n[0] * th + n[1];
    const B = n[2] * th * th + n[3] * th + n[4];
    const C = n[5] * th * th + n[6] * th + n[7];
    p = (2 * C / (-B + Math.sqrt(B * B - 4 * A * C))) ** 4 * 1e6;
  } else {
    throw new RangeError(`method = ${method} is not one of if97, antoine`);
  }
  return p / (1 + 0.57357 * S / (1000 - S));
}

/**
 * Hydrostatic pressure at depth in a compressible water column: integrates dP/dz = ρ(T, S, P) g with the EOS-80
 * density (midpoint rule, `n` steps, second-order accurate).
 * Reference: hydrostatic equilibrium; density by UNESCO EOS-80.
 * Validity: depth ≤ 10 000 m, temperature profile inside the EOS-80 range.
 * @param {object} o
 * @param {number} o.depth depth below the free surface, m (≥ 0)
 * @param {number|function(number):number} [o.T=283.15] temperature, K, or function of depth (m) returning K
 * @param {number} [o.S=35] salinity, g/kg
 * @param {number} [o.Patm=101325] pressure at the surface, Pa
 * @param {number} [o.g=9.80665] gravity, m/s²
 * @param {number} [o.n=50] integration steps
 * @returns {{ P: number, gauge: number, rhoMean: number, rhoBottom: number }} absolute pressure Pa, gauge
 *   pressure Pa, mean and bottom density kg/m³
 */
export function hydrostaticPressure({ depth, T = 283.15, S = 35, Patm = ATM, g = G0, n = 50 }) {
  need('depth', depth, 0, 10000);
  need('n', n, 1, 1e6);
  const Tz = typeof T === 'function' ? T : () => T;
  const dz = depth / n;
  let P = Patm;
  for (let i = 0; i < n; i++) {
    const zMid = (i + 0.5) * dz;
    const Tm = Tz(zMid);
    const r1 = seawaterDensity({ T: Tm, S, P });
    const r2 = seawaterDensity({ T: Tm, S, P: P + 0.5 * r1 * g * dz });
    P += r2 * g * dz;
  }
  const rhoBottom = seawaterDensity({ T: Tz(depth), S, P });
  return { P, gauge: P - Patm, rhoMean: depth > 0 ? (P - Patm) / (g * depth) : rhoBottom, rhoBottom };
}
