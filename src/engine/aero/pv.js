// Photovoltaics: solar geometry, clear-sky and decomposition models, transposition, optical and thermal
// models, single-diode module model with datasheet parameter extraction, string/array mismatch with bypass
// diodes, inverter models and an energy-yield chain for time series or monthly climatology.
//
// Units: irradiance W/m², temperatures °C (argument names end in C or are documented), angles in degrees where
// the name ends in Deg, power W, energy Wh. Azimuths are clockwise from north (south = 180°). Longitude east
// positive. Times are milliseconds since 1970-01-01 UTC.

const DEG = Math.PI / 180;
const K_BOLTZ = 1.380649e-23; // J/K
const Q_E = 1.602176634e-19; // C
const K_EV = K_BOLTZ / Q_E; // eV/K

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function finitePos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}
const sind = (d) => Math.sin(d * DEG);
const cosd = (d) => Math.cos(d * DEG);
const mod = (a, n) => ((a % n) + n) % n;

function sunCore(time) {
  const jd = time / 86400000 + 2440587.5;
  const T = (jd - 2451545) / 36525;
  const L0 = mod(280.46646 + T * (36000.76983 + T * 0.0003032), 360);
  const M = 357.52911 + T * (35999.05029 - 0.0001537 * T);
  const e = 0.016708634 - T * (0.000042037 + 0.0000001267 * T);
  const C = sind(M) * (1.914602 - T * (0.004817 + 0.000014 * T)) + sind(2 * M) * (0.019993 - 0.000101 * T) + sind(3 * M) * 0.000289;
  const trueLong = L0 + C;
  const R = (1.000001018 * (1 - e * e)) / (1 + e * cosd(M + C));
  const omega = 125.04 - 1934.136 * T;
  const lambda = trueLong - 0.00569 - 0.00478 * sind(omega);
  const eps0 = 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - T * 0.001813))) / 60) / 60;
  const eps = eps0 + 0.00256 * cosd(omega);
  const decl = Math.asin(sind(eps) * sind(lambda)) / DEG;
  const y = Math.tan((eps / 2) * DEG) ** 2;
  const eot = (4 * (y * sind(2 * L0) - 2 * e * sind(M) + 4 * e * y * sind(M) * cosd(2 * L0) - 0.5 * y * y * sind(4 * L0) - 1.25 * e * e * sind(2 * M))) / DEG;
  return { decl, eot, R };
}

/**
 * Solar position by the NOAA Solar Calculator algorithm (Meeus, Astronomical Algorithms, low-accuracy solar
 * coordinates; accuracy about 0.01° in position for years 1900–2100).
 * Optional atmospheric refraction (NOAA piecewise correction) gives the apparent elevation.
 * @param {object} o
 * @param {number} o.time instant (ms since epoch, UTC)
 * @param {number} o.lat latitude (deg, north positive)
 * @param {number} o.lon longitude (deg, east positive)
 * @returns {{zenithDeg:number, elevationDeg:number, apparentElevationDeg:number, azimuthDeg:number,
 *   declinationDeg:number, equationOfTimeMin:number, hourAngleDeg:number, distanceAU:number}} azimuth clockwise
 *   from north; hour angle negative before solar noon
 */
export function solarPosition({ time, lat, lon }) {
  need(Number.isFinite(time), 'time must be finite');
  need(Math.abs(lat) <= 90, 'lat must be within ±90');
  need(Math.abs(lon) <= 180, 'lon must be within ±180');
  const { decl, eot, R } = sunCore(time);
  const minutes = mod(time / 60000, 1440);
  const tst = mod(minutes + eot + 4 * lon, 1440);
  const ha = tst / 4 - 180;
  const cz = Math.min(1, Math.max(-1, sind(lat) * sind(decl) + cosd(lat) * cosd(decl) * cosd(ha)));
  const zen = Math.acos(cz) / DEG;
  const sz = Math.sin(zen * DEG);
  let az;
  const den = cosd(lat) * sz;
  if (Math.abs(den) < 1e-12) az = lat > decl ? 180 : 0;
  else {
    const c = Math.min(1, Math.max(-1, (sind(lat) * cz - sind(decl)) / den));
    const a = Math.acos(c) / DEG;
    az = ha > 0 ? mod(a + 180, 360) : mod(540 - a, 360);
  }
  const el = 90 - zen;
  let refr;
  if (el > 85) refr = 0;
  else {
    const te = Math.tan(el * DEG);
    if (el > 5) refr = 58.1 / te - 0.07 / te ** 3 + 0.000086 / te ** 5;
    else if (el > -0.575) refr = 1735 + el * (-518.2 + el * (103.4 + el * (-12.79 + el * 0.711)));
    else refr = -20.772 / te;
    refr /= 3600;
  }
  return {
    zenithDeg: zen,
    elevationDeg: el,
    apparentElevationDeg: el + refr,
    azimuthDeg: az,
    declinationDeg: decl,
    equationOfTimeMin: eot,
    hourAngleDeg: ha,
    distanceAU: R,
  };
}

/**
 * Sunrise, sunset, solar noon and day length for the UTC day containing `time` (NOAA Solar Calculator):
 * hour angle at the horizon H0 = acos(cos z0 / (cos φ cos δ) − tan φ tan δ), day length = 2 H0 / 15 h.
 * z0 = 90.833° accounts for refraction and the solar semi-diameter; z0 = 90° gives the geometric day length.
 * @param {object} o
 * @param {number} o.time any instant of the day (ms since epoch, UTC)
 * @param {number} o.lat latitude (deg)
 * @param {number} o.lon longitude (deg, east positive)
 * @param {number} [o.horizonZenithDeg=90.833] zenith angle defining sunrise/sunset (deg)
 * @returns {{sunrise:number, sunset:number, solarNoon:number, dayLengthHours:number, polar:'day'|'night'|null}}
 *   instants in ms since epoch (NaN in polar day/night)
 */
export function sunriseSunset({ time, lat, lon, horizonZenithDeg = 90.833 }) {
  need(Number.isFinite(time), 'time must be finite');
  need(Math.abs(lat) <= 90, 'lat must be within ±90');
  const day0 = Math.floor(time / 86400000) * 86400000;
  const { decl, eot } = sunCore(day0 + 43200000);
  const noonMin = 720 - 4 * lon - eot;
  const c = cosd(horizonZenithDeg) / (cosd(lat) * cosd(decl)) - Math.tan(lat * DEG) * Math.tan(decl * DEG);
  const solarNoon = day0 + noonMin * 60000;
  if (c <= -1) return { sunrise: NaN, sunset: NaN, solarNoon, dayLengthHours: 24, polar: 'day' };
  if (c >= 1) return { sunrise: NaN, sunset: NaN, solarNoon, dayLengthHours: 0, polar: 'night' };
  const H0 = Math.acos(c) / DEG;
  return {
    sunrise: solarNoon - H0 * 4 * 60000,
    sunset: solarNoon + H0 * 4 * 60000,
    solarNoon,
    dayLengthHours: (2 * H0) / 15,
    polar: null,
  };
}

/**
 * Extraterrestrial normal irradiance, E0n = G_sc / R² with R the Earth–Sun distance in AU from the solar
 * position algorithm. G_sc = 1361 W/m² (Kopp & Lean 2011); 1367 W/m² is the older WMO value.
 * @param {number} time instant (ms since epoch, UTC)
 * @param {number} [solarConstant=1361] W/m²
 * @returns {number} W/m²
 */
export function extraterrestrialIrradiance(time, solarConstant = 1361) {
  const { R } = sunCore(time);
  return solarConstant / (R * R);
}

/**
 * Relative optical air mass, Kasten & Young (1989) Applied Optics 28:4735:
 *   AM = 1 / (cos z + 0.50572 (96.07995 − z)^−1.6364), z in degrees; absolute air mass = AM · p/101325.
 * @param {number} zenithDeg (apparent) zenith angle (deg), < 96.07
 * @param {number} [pressure=101325] site pressure (Pa)
 * @returns {{relative:number, absolute:number}} NaN for z beyond the formula's range
 */
export function airMass(zenithDeg, pressure = 101325) {
  if (!(zenithDeg < 96.07995)) return { relative: NaN, absolute: NaN };
  const rel = 1 / (cosd(zenithDeg) + 0.50572 * Math.pow(96.07995 - zenithDeg, -1.6364));
  return { relative: rel, absolute: (rel * pressure) / 101325 };
}

/**
 * Clear-sky global horizontal irradiance, Haurwitz (1945) as given by Reno, Hansen & Stein (2012)
 * SAND2012-2389: GHI = 1098 cos z exp(−0.059 / cos z). Zenith only; no site inputs.
 * @param {number} zenithDeg zenith angle (deg)
 * @returns {number} GHI (W/m²)
 */
export function clearSkyHaurwitz(zenithDeg) {
  const cz = cosd(zenithDeg);
  return cz > 0 ? 1098 * cz * Math.exp(-0.059 / cz) : 0;
}

/**
 * Clear-sky irradiance, Ineichen & Perez (2002) Solar Energy 73:151–157 with Linke turbidity:
 *   f_h1 = exp(−alt/8000), f_h2 = exp(−alt/1250), c_g1 = 5.09e-5 alt + 0.868, c_g2 = 3.92e-5 alt + 0.0387,
 *   GHI = c_g1 E0n cos z exp(−c_g2 AM (f_h1 + f_h2 (T_L − 1))),
 *   b = 0.664 + 0.163/f_h1,  DNI = min( b E0n exp(−0.09 AM (T_L − 1)),
 *                                       GHI [1 − (0.1 − 0.2 exp(−T_L))/(0.1 + 0.882/f_h1)] / cos z ),
 *   DHI = GHI − DNI cos z,
 * AM being the absolute (pressure-corrected) air mass; implemented as in Reno et al. (2012) SAND2012-2389.
 * @param {object} o
 * @param {number} o.zenithDeg apparent zenith angle (deg)
 * @param {number} o.linkeTurbidity Linke turbidity T_L (-), typically 2 (very clean) to 6 (polluted/humid)
 * @param {number} [o.altitude=0] site altitude (m)
 * @param {number} [o.dniExtra=1361] extraterrestrial normal irradiance (W/m²)
 * @param {number} [o.pressure] site pressure (Pa); default from altitude, 101325 exp(−alt/8434.5)
 * @returns {{ghi:number, dni:number, dhi:number}} W/m²
 */
export function clearSkyIneichen({ zenithDeg, linkeTurbidity, altitude = 0, dniExtra = 1361, pressure }) {
  finitePos(linkeTurbidity, 'linkeTurbidity');
  const cz = cosd(zenithDeg);
  if (!(cz > 0)) return { ghi: 0, dni: 0, dhi: 0 };
  const p = pressure ?? 101325 * Math.exp(-altitude / 8434.5);
  const am = airMass(zenithDeg, p).absolute;
  const tl = linkeTurbidity;
  const fh1 = Math.exp(-altitude / 8000);
  const fh2 = Math.exp(-altitude / 1250);
  const cg1 = 5.09e-5 * altitude + 0.868;
  const cg2 = 3.92e-5 * altitude + 0.0387;
  const ghi = cg1 * dniExtra * cz * Math.exp(-cg2 * am * (fh1 + fh2 * (tl - 1)));
  const b = 0.664 + 0.163 / fh1;
  const bnci = b * Math.exp(-0.09 * am * (tl - 1)) * dniExtra;
  const f2 = Math.max(0, (1 - (0.1 - 0.2 * Math.exp(-tl)) / (0.1 + 0.882 / fh1)) / cz);
  const dni = Math.min(bnci, ghi * f2);
  return { ghi, dni, dhi: ghi - dni * cz };
}

/**
 * Erbs, Klein & Duffie (1982) Solar Energy 28:293–302 diffuse-fraction correlation (hourly):
 *   k_t ≤ 0.22:        k_d = 1 − 0.09 k_t
 *   0.22 < k_t ≤ 0.80: k_d = 0.9511 − 0.1604 k_t + 4.388 k_t² − 16.638 k_t³ + 12.336 k_t⁴
 *   k_t > 0.80:        k_d = 0.165
 * with k_t = GHI / (E0n cos z). DNI = (GHI − DHI)/cos z; for zenith > 87° all irradiance is treated as diffuse.
 * @param {object} o
 * @param {number} o.ghi global horizontal irradiance (W/m²)
 * @param {number} o.zenithDeg zenith angle (deg)
 * @param {number} [o.dniExtra=1361] extraterrestrial normal irradiance (W/m²)
 * @returns {{kt:number, kd:number, dhi:number, dni:number}}
 */
export function erbsDecomposition({ ghi, zenithDeg, dniExtra = 1361 }) {
  need(Number.isFinite(ghi) && ghi >= 0, 'ghi must be >= 0');
  const cz = cosd(zenithDeg);
  if (!(cz > 0) || ghi === 0) return { kt: 0, kd: 1, dhi: ghi, dni: 0 };
  const kt = Math.min(2, ghi / (dniExtra * Math.max(cz, 0.065)));
  const kd = erbsDiffuseFraction(kt);
  if (zenithDeg > 87) return { kt, kd: 1, dhi: ghi, dni: 0 };
  const dhi = kd * ghi;
  return { kt, kd, dhi, dni: (ghi - dhi) / cz };
}

/**
 * Erbs diffuse fraction k_d as a function of the clearness index k_t (see `erbsDecomposition`).
 * @param {number} kt clearness index (-)
 * @returns {number} diffuse fraction (-)
 */
export function erbsDiffuseFraction(kt) {
  need(Number.isFinite(kt) && kt >= 0, 'kt must be >= 0');
  if (kt <= 0.22) return 1 - 0.09 * kt;
  if (kt <= 0.8) return 0.9511 - 0.1604 * kt + 4.388 * kt ** 2 - 16.638 * kt ** 3 + 12.336 * kt ** 4;
  return 0.165;
}

/**
 * Angle of incidence on a tilted plane: cos θ = cos z cos β + sin z sin β cos(γ_s − γ).
 * @param {object} o
 * @param {number} o.zenithDeg solar zenith (deg)
 * @param {number} o.azimuthDeg solar azimuth (deg, clockwise from north)
 * @param {number} o.tiltDeg surface tilt from horizontal β (deg)
 * @param {number} o.surfaceAzimuthDeg surface azimuth γ (deg, clockwise from north; 180 = facing south)
 * @returns {{cosAoi:number, aoiDeg:number}} cosAoi is not clipped (negative = sun behind the plane)
 */
export function angleOfIncidence({ zenithDeg, azimuthDeg, tiltDeg, surfaceAzimuthDeg }) {
  const c = cosd(zenithDeg) * cosd(tiltDeg) + sind(zenithDeg) * sind(tiltDeg) * cosd(azimuthDeg - surfaceAzimuthDeg);
  return { cosAoi: c, aoiDeg: Math.acos(Math.min(1, Math.max(-1, c))) / DEG };
}

/**
 * Transposition of horizontal irradiance components to the plane of array.
 *  - 'isotropic' (Liu & Jordan 1963): sky = DHI (1 + cos β)/2.
 *  - 'haydavies' (Hay & Davies 1980): sky = DHI [A_i R_b + (1 − A_i)(1 + cos β)/2], A_i = DNI/E0n,
 *    R_b = cos θ / cos z.
 * Ground (or water-surface) reflection: GHI ρ (1 − cos β)/2; calm water albedo is about 0.06–0.10.
 * @param {object} o
 * @param {number} o.dni direct normal irradiance (W/m²)
 * @param {number} o.dhi diffuse horizontal irradiance (W/m²)
 * @param {number} [o.ghi] global horizontal irradiance (W/m²); default DNI cos z + DHI
 * @param {number} o.zenithDeg solar zenith (deg)
 * @param {number} o.azimuthDeg solar azimuth (deg)
 * @param {number} o.tiltDeg surface tilt (deg)
 * @param {number} o.surfaceAzimuthDeg surface azimuth (deg)
 * @param {number} [o.albedo=0.08] surface albedo (-)
 * @param {'isotropic'|'haydavies'} [o.model='haydavies']
 * @param {number} [o.dniExtra=1361] extraterrestrial normal irradiance (W/m²), Hay–Davies only
 * @returns {{total:number, beam:number, skyDiffuse:number, groundReflected:number, aoiDeg:number,
 *   circumsolar:number}} W/m²; circumsolar is the beam-like part of the sky diffuse (Hay–Davies)
 */
export function planeOfArray({ dni, dhi, ghi, zenithDeg, azimuthDeg, tiltDeg, surfaceAzimuthDeg, albedo = 0.08, model = 'haydavies', dniExtra = 1361 }) {
  need(dni >= 0 && dhi >= 0, 'dni and dhi must be >= 0');
  need(albedo >= 0 && albedo <= 1, 'albedo must be within [0, 1]');
  need(tiltDeg >= 0 && tiltDeg <= 180, 'tiltDeg must be within [0, 180]');
  const cz = Math.max(0, cosd(zenithDeg));
  const g = ghi ?? dni * cz + dhi;
  const { cosAoi, aoiDeg } = angleOfIncidence({ zenithDeg, azimuthDeg, tiltDeg, surfaceAzimuthDeg });
  const ci = Math.max(0, cosAoi);
  const beam = cz > 0 ? dni * ci : 0;
  const iso = (1 + cosd(tiltDeg)) / 2;
  let sky;
  let circ = 0;
  if (model === 'isotropic') sky = dhi * iso;
  else {
    need(model === 'haydavies', "model must be 'isotropic' or 'haydavies'");
    const Ai = Math.min(1, dni / dniExtra);
    const Rb = ci / Math.max(cz, 0.01745);
    circ = dhi * Ai * Rb;
    sky = circ + dhi * (1 - Ai) * iso;
  }
  const ground = (g * albedo * (1 - cosd(tiltDeg))) / 2;
  return { total: beam + sky + ground, beam, skyDiffuse: sky, groundReflected: ground, aoiDeg, circumsolar: circ };
}

/**
 * Rear-side irradiance of a bifacial module by a simple isotropic view-factor model: the rear face (tilt
 * 180° − β) sees the ground with view factor (1 + cos β)/2 and the sky with (1 − cos β)/2,
 *   G_rear = ρ GHI (1 + cos β)/2 (1 − f_shade) + DHI (1 − cos β)/2,
 * where f_shade is the fraction of the ground seen by the rear that is shaded by the array (argument; row
 * geometry and height are not resolved). Bifacial gain = φ_bif G_rear / G_front.
 * @param {object} o
 * @param {number} o.ghi (W/m²)
 * @param {number} o.dhi (W/m²)
 * @param {number} o.tiltDeg module tilt (deg)
 * @param {number} o.albedo ground/water albedo (-)
 * @param {number} [o.groundShadeFactor=0] shaded fraction of the viewed ground (-)
 * @returns {number} rear irradiance (W/m²)
 */
export function rearIrradiance({ ghi, dhi, tiltDeg, albedo, groundShadeFactor = 0 }) {
  need(ghi >= 0 && dhi >= 0, 'ghi and dhi must be >= 0');
  need(albedo >= 0 && albedo <= 1, 'albedo must be within [0, 1]');
  need(groundShadeFactor >= 0 && groundShadeFactor <= 1, 'groundShadeFactor must be within [0, 1]');
  const cb = cosd(tiltDeg);
  return (albedo * ghi * (1 + cb) * (1 - groundShadeFactor)) / 2 + (dhi * (1 - cb)) / 2;
}

/**
 * Incidence-angle modifier, ASHRAE: IAM = 1 − b0 (1/cos θ − 1), b0 = 0.05 typical; zero beyond 90° or where the
 * expression turns negative. Reference: Souka & Safwat (1966); ASHRAE 93.
 * @param {number} aoiDeg angle of incidence (deg)
 * @param {number} [b0=0.05]
 * @returns {number} IAM (-)
 */
export function iamAshrae(aoiDeg, b0 = 0.05) {
  const c = cosd(aoiDeg);
  if (!(c > 0)) return 0;
  return Math.max(0, 1 - b0 * (1 / c - 1));
}

/**
 * Incidence-angle modifier from Snell and Fresnel laws with Bouguer absorption in the glazing (De Soto, Klein &
 * Beckman 2006, Solar Energy 80:78–88): θ_r = asin(sin θ / n),
 *   τ(θ) = exp(−K L / cos θ_r) [1 − ½ (sin²(θ_r − θ)/sin²(θ_r + θ) + tan²(θ_r − θ)/tan²(θ_r + θ))],
 * IAM = τ(θ)/τ(0), τ(0) = exp(−K L)[1 − ((1 − n)/(1 + n))²]. Defaults n = 1.526, K = 4 m⁻¹, L = 2 mm.
 * @param {number} aoiDeg angle of incidence (deg)
 * @param {{n?:number, K?:number, L?:number}} [o] refractive index, extinction coefficient (1/m), thickness (m)
 * @returns {number} IAM (-)
 */
export function iamPhysical(aoiDeg, { n = 1.526, K = 4, L = 0.002 } = {}) {
  if (!(aoiDeg < 90)) return 0;
  const tau0 = Math.exp(-K * L) * (1 - ((1 - n) / (1 + n)) ** 2);
  if (aoiDeg < 1e-6) return 1;
  const th = aoiDeg * DEG;
  const thr = Math.asin(Math.sin(th) / n);
  const rho = 0.5 * (Math.sin(thr - th) ** 2 / Math.sin(thr + th) ** 2 + Math.tan(thr - th) ** 2 / Math.tan(thr + th) ** 2);
  return (Math.exp((-K * L) / Math.cos(thr)) * (1 - rho)) / tau0;
}

/**
 * Cell temperature, NOCT model: T_c = T_a + (NOCT − 20) G / 800 (Ross 1980; Duffie & Beckman §23.3, without the
 * efficiency term).
 * @param {object} o
 * @param {number} o.poa plane-of-array irradiance (W/m²)
 * @param {number} o.tempAirC ambient temperature (°C)
 * @param {number} [o.noctC=45] nominal operating cell temperature (°C)
 * @returns {number} cell temperature (°C)
 */
export function cellTempNoct({ poa, tempAirC, noctC = 45 }) {
  return tempAirC + ((noctC - 20) * poa) / 800;
}

/**
 * Module temperature, Faiman (2008) Prog. Photovolt. 16:307–315: T_m = T_a + G / (U0 + U1 v).
 * Defaults U0 = 25 W/m²K, U1 = 6.84 W·s/m³K are common free-standing land values. Floating systems are
 * generally better cooled and published fits report higher heat-loss coefficients that depend strongly on the
 * float type and footprint (e.g. Liu et al. 2018; Dörenkämper et al. 2021) — pass technology-specific values. A fixed
 * `waterCoolingOffsetC` (≤ 0) can represent additional cooling by water contact or spray.
 * @param {object} o
 * @param {number} o.poa plane-of-array irradiance (W/m²)
 * @param {number} o.tempAirC ambient temperature (°C)
 * @param {number} [o.windSpeed=1] wind speed at module height (m/s)
 * @param {number} [o.u0=25] constant heat-loss coefficient (W/m²K)
 * @param {number} [o.u1=6.84] wind heat-loss coefficient (W·s/m³K)
 * @param {number} [o.waterCoolingOffsetC=0] additive temperature offset (°C)
 * @returns {number} module (≈ cell) temperature (°C)
 */
export function cellTempFaiman({ poa, tempAirC, windSpeed = 1, u0 = 25, u1 = 6.84, waterCoolingOffsetC = 0 }) {
  finitePos(u0 + u1 * windSpeed, 'u0 + u1·windSpeed');
  return tempAirC + poa / (u0 + u1 * windSpeed) + waterCoolingOffsetC;
}

/**
 * Cell temperature, Sandia model (King, Boyson & Kratochvil 2004, SAND2004-3535):
 *   T_m = G exp(a + b v) + T_a,  T_c = T_m + (G/1000) ΔT.
 * Defaults are the open-rack glass/glass set a = −3.47, b = −0.0594, ΔT = 3 °C (open-rack glass/polymer:
 * −3.56, −0.075, 3).
 * @param {object} o
 * @param {number} o.poa (W/m²)
 * @param {number} o.tempAirC (°C)
 * @param {number} [o.windSpeed=1] wind speed at 10 m (m/s)
 * @param {number} [o.a=-3.47]
 * @param {number} [o.b=-0.0594]
 * @param {number} [o.deltaT=3] cell-to-module temperature difference at 1000 W/m² (°C)
 * @returns {number} cell temperature (°C)
 */
export function cellTempSandia({ poa, tempAirC, windSpeed = 1, a = -3.47, b = -0.0594, deltaT = 3 }) {
  return poa * Math.exp(a + b * windSpeed) + tempAirC + (poa / 1000) * deltaT;
}

/**
 * Principal branch of the Lambert W function, W(x) e^W(x) = x for x ≥ −1/e (Halley iteration; Corless et al.
 * 1996).
 * @param {number} x argument, >= −1/e
 * @returns {number} W0(x)
 */
export function lambertW(x) {
  need(Number.isFinite(x) && x >= -1 / Math.E - 1e-15, 'x must be >= -1/e');
  if (x === 0) return 0;
  let w;
  if (x < -0.25) {
    const p = Math.sqrt(Math.max(0, 2 * (Math.E * x + 1)));
    w = -1 + p - (p * p) / 3 + (11 * p ** 3) / 72;
  } else if (x < 3) w = x / (1 + x) + 0.3 * Math.log1p(x) * (x > 0 ? 1 : 0);
  else w = Math.log(x) - Math.log(Math.log(x));
  for (let i = 0; i < 50; i++) {
    const ew = Math.exp(w);
    const f = w * ew - x;
    const d = ew * (w + 1) - ((w + 2) * f) / (2 * w + 2);
    if (!Number.isFinite(d) || d === 0) break;
    const dw = f / d;
    w -= dw;
    if (Math.abs(dw) <= 1e-15 * (1 + Math.abs(w))) break;
  }
  return w;
}

/** W0(exp(L)) for arguments too large to exponentiate: solves w + ln w = L by Newton. */
function lambertWofExp(L) {
  if (L < 500) return lambertW(Math.exp(L));
  let w = L - Math.log(L);
  for (let i = 0; i < 50; i++) {
    const dw = (w + Math.log(w) - L) / (1 + 1 / w);
    w -= dw;
    if (Math.abs(dw) < 1e-14 * w) break;
  }
  return w;
}

function checkParams(p) {
  need(p && Number.isFinite(p.IL) && Number.isFinite(p.I0) && p.I0 > 0, 'params need finite IL and positive I0');
  finitePos(p.a, 'params.a');
  need(p.Rs >= 0, 'params.Rs must be >= 0');
  need(p.Rsh > 0, 'params.Rsh must be > 0');
}

/**
 * Current of the single-diode model at a given voltage,
 *   I = I_L − I_0 [exp((V + I R_s)/a) − 1] − (V + I R_s)/R_sh,
 * solved explicitly with the Lambert W function (Jain & Kapoor 2004, Solar Energy Mater. Solar Cells 81:269).
 * @param {{IL:number, I0:number, a:number, Rs:number, Rsh:number}} params light current (A), saturation
 *   current (A), modified ideality a = n N_s k T/q (V), series and shunt resistance (Ω; Rsh may be Infinity)
 * @param {number} V terminal voltage (V)
 * @returns {number} current (A)
 */
export function singleDiodeCurrent(params, V) {
  checkParams(params);
  const { IL, I0, a, Rs, Rsh } = params;
  const gsh = Number.isFinite(Rsh) ? 1 / Rsh : 0;
  if (Rs === 0) return IL - I0 * Math.expm1(V / a) - V * gsh;
  const k = 1 + Rs * gsh; // (Rs + Rsh)/Rsh
  const lnTheta = Math.log((Rs * I0) / (a * k)) + (Rs * (IL + I0) + V) / (a * k);
  return (IL + I0 - V * gsh) / k - (a / Rs) * lambertWofExp(lnTheta);
}

/**
 * Voltage of the single-diode model at a given current (explicit Lambert-W form):
 *   V = (I_L + I_0 − I) R_sh − I R_s − a W( (I_0 R_sh/a) exp((I_L + I_0 − I) R_sh / a) ).
 * @param {{IL:number, I0:number, a:number, Rs:number, Rsh:number}} params
 * @param {number} I current (A)
 * @returns {number} voltage (V); −Infinity where the model has no forward solution (I ≥ I_L + I_0, R_sh = ∞)
 */
export function singleDiodeVoltage(params, I) {
  checkParams(params);
  const { IL, I0, a, Rs, Rsh } = params;
  if (!Number.isFinite(Rsh)) {
    const arg = (IL - I) / I0 + 1;
    return arg > 0 ? a * Math.log(arg) - I * Rs : -Infinity;
  }
  const z = (IL + I0 - I) * Rsh;
  return z - I * Rs - a * lambertWofExp(Math.log((I0 * Rsh) / a) + z / a);
}

/**
 * Characteristic points and sampled I–V curve of a single-diode model.
 * @param {{IL:number, I0:number, a:number, Rs:number, Rsh:number}} params
 * @param {number} [nPoints=0] number of curve samples (0 = none)
 * @returns {{isc:number, voc:number, vmp:number, imp:number, pmp:number, fillFactor:number, v:number[],
 *   i:number[]}} A, V, V, A, W, -, sampled voltages and currents
 */
export function singleDiodeCurve(params, nPoints = 0) {
  const isc = singleDiodeCurrent(params, 0);
  const voc = singleDiodeVoltage(params, 0);
  const P = (V) => V * singleDiodeCurrent(params, V);
  const gr = (Math.sqrt(5) - 1) / 2;
  let lo = 0;
  let hi = voc;
  let c = hi - gr * (hi - lo);
  let d = lo + gr * (hi - lo);
  let pc = P(c);
  let pd = P(d);
  for (let it = 0; it < 80; it++) {
    if (pc > pd) {
      hi = d;
      d = c;
      pd = pc;
      c = hi - gr * (hi - lo);
      pc = P(c);
    } else {
      lo = c;
      c = d;
      pc = pd;
      d = lo + gr * (hi - lo);
      pd = P(d);
    }
  }
  const vmp = 0.5 * (lo + hi);
  const imp = singleDiodeCurrent(params, vmp);
  const v = [];
  const i = [];
  for (let k = 0; k < nPoints; k++) {
    const vv = (voc * k) / (nPoints - 1);
    v.push(vv);
    i.push(singleDiodeCurrent(params, vv));
  }
  return { isc, voc, vmp, imp, pmp: vmp * imp, fillFactor: (vmp * imp) / (isc * voc), v, i };
}

/**
 * Translate reference single-diode parameters to operating irradiance and cell temperature (De Soto, Klein &
 * Beckman 2006; air-mass modifier omitted):
 *   a = a_ref T/T_ref,  I_L = (G/G_ref)[I_L,ref + α_Isc (T − T_ref)],
 *   I_0 = I_0,ref (T/T_ref)³ exp[(E_g,ref/T_ref − E_g/T)/k],  E_g = E_g,ref [1 − 0.0002677 (T − T_ref)],
 *   R_sh = R_sh,ref G_ref/G,  R_s = R_s,ref;  E_g,ref = 1.121 eV for silicon.
 * @param {{IL:number, I0:number, a:number, Rs:number, Rsh:number, alphaIsc:number, EgRef?:number,
 *   dEgdT?:number}} ref reference parameters at 1000 W/m², 25 °C; alphaIsc in A/K
 * @param {object} o
 * @param {number} o.G effective irradiance (W/m²)
 * @param {number} o.tempCellC cell temperature (°C)
 * @returns {{IL:number, I0:number, a:number, Rs:number, Rsh:number}}
 */
export function singleDiodeAt(ref, { G, tempCellC }) {
  checkParams(ref);
  need(G >= 0, 'G must be >= 0');
  const Tref = 298.15;
  const T = tempCellC + 273.15;
  finitePos(T, 'tempCellC + 273.15');
  const Eg0 = ref.EgRef ?? 1.121;
  const Eg = Eg0 * (1 + (ref.dEgdT ?? -0.0002677) * (T - Tref));
  const g = Math.max(G, 1e-9) / 1000;
  return {
    IL: g * (ref.IL + (ref.alphaIsc ?? 0) * (T - Tref)),
    I0: ref.I0 * (T / Tref) ** 3 * Math.exp((Eg0 / Tref - Eg / T) / K_EV),
    a: (ref.a * T) / Tref,
    Rs: ref.Rs,
    Rsh: ref.Rsh / g,
  };
}

/**
 * Five-parameter single-diode model from datasheet values at standard test conditions (De Soto, Klein &
 * Beckman 2006): the parameters satisfy
 *   (1) I(0) = I_sc, (2) I(V_oc) = 0, (3) I(V_mp) = I_mp, (4) dP/dV = 0 at the maximum-power point,
 *   (5) the model open-circuit voltage at T_ref + 10 K equals V_oc + 10 β_Voc.
 * Solution strategy: for trial (a, R_s) conditions (1)–(3) are linear in (I_L, I_0, 1/R_sh); R_s is found from
 * (4) by bracketing, and a from (5) by bracketing (or fixed from the diode ideality `n` when β_Voc is not
 * used). Throws when the datasheet values admit no physical solution (R_sh > 0, I_0 > 0).
 * @param {object} o
 * @param {number} o.Isc short-circuit current (A)
 * @param {number} o.Voc open-circuit voltage (V)
 * @param {number} o.Imp current at maximum power (A)
 * @param {number} o.Vmp voltage at maximum power (V)
 * @param {number} o.Ns number of cells in series (-)
 * @param {number} [o.alphaIsc=0] temperature coefficient of I_sc (A/K)
 * @param {number} [o.betaVoc] temperature coefficient of V_oc (V/K, negative)
 * @param {number} [o.n] diode ideality factor; used instead of condition (5) when given
 * @returns {{IL:number, I0:number, a:number, Rs:number, Rsh:number, n:number, alphaIsc:number, Ns:number,
 *   residualBetaVoc:number}} reference parameters at 1000 W/m², 25 °C; n = a q/(N_s k T_ref)
 */
export function extractSingleDiode({ Isc, Voc, Imp, Vmp, Ns, alphaIsc = 0, betaVoc, n }) {
  finitePos(Isc, 'Isc');
  finitePos(Voc, 'Voc');
  finitePos(Imp, 'Imp');
  finitePos(Vmp, 'Vmp');
  need(Number.isInteger(Ns) && Ns > 0, 'Ns must be a positive integer');
  need(Imp < Isc && Vmp < Voc, 'Imp must be below Isc and Vmp below Voc');
  need(n !== undefined || betaVoc !== undefined, 'either betaVoc or n must be given');
  const Tref = 298.15;
  const vt = Ns * K_EV * Tref;

  // (1)–(3): linear in IL, I0, Gsh for given a, Rs
  const linear = (a, Rs) => {
    const e1 = Math.expm1((Isc * Rs) / a);
    const e2 = Math.expm1(Voc / a);
    const e3 = Math.expm1((Vmp + Imp * Rs) / a);
    const v1 = Isc * Rs;
    const v3 = Vmp + Imp * Rs;
    // subtract (2) from (1) and (3): eliminates IL
    //  I0 (e2 − e1) + Gsh (Voc − v1) = Isc ;  I0 (e2 − e3) + Gsh (Voc − v3) = Imp
    const a11 = e2 - e1;
    const a12 = Voc - v1;
    const a21 = e2 - e3;
    const a22 = Voc - v3;
    const det = a11 * a22 - a12 * a21;
    const I0 = (Isc * a22 - a12 * Imp) / det;
    const Gsh = (a11 * Imp - a21 * Isc) / det;
    const IL = I0 * e2 + Voc * Gsh;
    return { IL, I0, Gsh, e3 };
  };
  const res4 = (a, Rs) => {
    const p = linear(a, Rs);
    const g = (p.I0 / a) * (p.e3 + 1) + p.Gsh;
    return { r: Imp - (Vmp * g) / (1 + Rs * g), p };
  };
  const solveRs = (a) => {
    const RsMax = (Voc - Vmp) / Imp;
    const N = 400;
    let prev = null;
    for (let i = 0; i < N; i++) {
      const Rs = (RsMax * i) / N;
      const cur = res4(a, Rs);
      const valid = cur.p.I0 > 0 && cur.p.Gsh >= 0 && Number.isFinite(cur.r);
      if (prev && valid && prev.valid && prev.r * cur.r <= 0) {
        let lo = prev.Rs;
        let hi = Rs;
        let flo = prev.r;
        for (let it = 0; it < 80; it++) {
          const mid = 0.5 * (lo + hi);
          const fm = res4(a, mid).r;
          if (flo * fm <= 0) hi = mid;
          else {
            lo = mid;
            flo = fm;
          }
        }
        const RsSol = 0.5 * (lo + hi);
        const p = linear(a, RsSol);
        if (p.I0 > 0 && p.Gsh >= 0) return { Rs: RsSol, ...p };
      }
      prev = { Rs, r: cur.r, valid };
    }
    return null;
  };
  const pack = (a, s) => ({ IL: s.IL, I0: s.I0, a, Rs: s.Rs, Rsh: s.Gsh > 0 ? 1 / s.Gsh : Infinity, alphaIsc });
  const res5 = (a) => {
    const s = solveRs(a);
    if (!s) return null;
    const ref = pack(a, s);
    const hot = singleDiodeAt(ref, { G: 1000, tempCellC: 35 });
    return { r: singleDiodeVoltage(hot, 0) - (Voc + 10 * betaVoc), s };
  };
  let aSol;
  let sSol;
  let resid = 0;
  if (n !== undefined) {
    finitePos(n, 'n');
    aSol = n * vt;
    sSol = solveRs(aSol);
    need(sSol, 'no physical single-diode solution for the given datasheet values and ideality factor');
  } else {
    let prev = null;
    let found = false;
    for (let i = 0; i <= 240 && !found; i++) {
      const a = vt * (0.5 + (2.0 * i) / 240);
      const cur = res5(a);
      if (cur && prev && prev.r * cur.r <= 0) {
        let lo = prev.a;
        let hi = a;
        let flo = prev.r;
        for (let it = 0; it < 60; it++) {
          const mid = 0.5 * (lo + hi);
          const m = res5(mid);
          if (!m) break;
          if (flo * m.r <= 0) hi = mid;
          else {
            lo = mid;
            flo = m.r;
          }
        }
        aSol = 0.5 * (lo + hi);
        const fin = res5(aSol);
        need(fin, 'single-diode extraction failed to converge');
        sSol = fin.s;
        resid = fin.r;
        found = true;
      }
      prev = cur ? { a, r: cur.r } : null;
    }
    need(found, 'no single-diode solution reproduces the datasheet values with the given betaVoc; supply n instead');
  }
  const out = pack(aSol, sSol);
  return { ...out, n: aSol / vt, Ns, residualBetaVoc: resid };
}

/**
 * PVWatts DC model: P_dc = P_dc0 (G/1000) [1 + γ (T_c − 25)] (Dobos 2014, NREL/TP-6A20-62641).
 * @param {object} o
 * @param {number} o.pdc0 array rating at STC (W)
 * @param {number} o.poa effective plane-of-array irradiance (W/m²)
 * @param {number} o.tempCellC cell temperature (°C)
 * @param {number} [o.gammaPdc=-0.0035] power temperature coefficient (1/K); −0.0026 to −0.0045 for c-Si
 * @returns {number} DC power (W)
 */
export function pvwattsDc({ pdc0, poa, tempCellC, gammaPdc = -0.0035 }) {
  need(pdc0 >= 0, 'pdc0 must be >= 0');
  need(poa >= 0, 'poa must be >= 0');
  return Math.max(0, pdc0 * (poa / 1000) * (1 + gammaPdc * (tempCellC - 25)));
}

/**
 * I–V curve of series-connected cell groups, each protected by a bypass diode, under non-uniform irradiance
 * (I–V summation at common current). A group is a fraction of a module: its voltage at current I is
 * fraction × V_module(I; G, T_c); when that falls below −V_d the bypass diode conducts and the group voltage is
 * −V_d. Reverse-bias breakdown of shaded cells before the diode conducts is not modelled.
 * @param {object} o
 * @param {{params:object, fraction?:number, G:number, tempCellC:number}[]} o.groups reference single-diode
 *   parameters of the module the group belongs to, the fraction of that module's cells in the group (default 1),
 *   irradiance (W/m²) and cell temperature (°C)
 * @param {number} [o.diodeDrop=0.5] bypass-diode forward voltage V_d (V)
 * @param {number} [o.nPoints=400] current samples
 * @returns {{i:Float64Array, v:Float64Array, pmp:number, vmp:number, imp:number, isc:number, voc:number,
 *   sumOfGroupPmp:number, mismatchLoss:number}} mismatchLoss = 1 − P_mp / Σ group P_mp
 */
export function seriesStringCurve({ groups, diodeDrop = 0.5, nPoints = 400 }) {
  need(Array.isArray(groups) && groups.length > 0, 'groups must be a non-empty array');
  need(diodeDrop >= 0, 'diodeDrop must be >= 0');
  const ps = groups.map((g) => ({ p: singleDiodeAt(g.params, { G: g.G, tempCellC: g.tempCellC }), f: g.fraction ?? 1 }));
  let iMax = 0;
  let sumP = 0;
  for (const g of ps) {
    iMax = Math.max(iMax, g.p.IL);
    sumP += g.f * singleDiodeCurve(g.p).pmp;
  }
  const vAt = (I) => {
    let v = 0;
    for (const g of ps) v += Math.max(-diodeDrop, g.f * singleDiodeVoltage(g.p, I));
    return v;
  };
  const i = new Float64Array(nPoints);
  const v = new Float64Array(nPoints);
  let best = 0;
  let bi = 0;
  for (let k = 0; k < nPoints; k++) {
    i[k] = (iMax * k) / (nPoints - 1);
    v[k] = vAt(i[k]);
    const p = i[k] * v[k];
    if (p > best) {
      best = p;
      bi = k;
    }
  }
  // local refinement around the best sample (the global curve may have several maxima)
  let lo = i[Math.max(0, bi - 1)];
  let hi = i[Math.min(nPoints - 1, bi + 1)];
  const gr = (Math.sqrt(5) - 1) / 2;
  for (let it = 0; it < 60; it++) {
    const c = hi - gr * (hi - lo);
    const d = lo + gr * (hi - lo);
    if (c * vAt(c) > d * vAt(d)) hi = d;
    else lo = c;
  }
  const imp = 0.5 * (lo + hi);
  const vmp = vAt(imp);
  // short-circuit current: V(I) = 0
  let a = 0;
  let b = iMax * 1.001 + 1e-9;
  for (let it = 0; it < 80; it++) {
    const m = 0.5 * (a + b);
    if (vAt(m) > 0) a = m;
    else b = m;
  }
  const pmp = Math.max(best, imp * vmp);
  return { i, v, pmp, vmp, imp, isc: 0.5 * (a + b), voc: v[0], sumOfGroupPmp: sumP, mismatchLoss: sumP > 0 ? 1 - pmp / sumP : 0 };
}

/**
 * Parallel connection of strings: currents add at common voltage (I–V summation). Each string curve is given as
 * samples with current ascending and voltage descending (as returned by `seriesStringCurve`).
 * @param {{i:ArrayLike<number>, v:ArrayLike<number>}[]} strings
 * @param {number} [nPoints=400] voltage samples
 * @returns {{v:Float64Array, i:Float64Array, pmp:number, vmp:number, imp:number, sumOfStringPmp:number,
 *   mismatchLoss:number}}
 */
export function parallelStrings(strings, nPoints = 400) {
  need(Array.isArray(strings) && strings.length > 0, 'strings must be a non-empty array');
  const vocMax = Math.max(...strings.map((s) => s.v[0]));
  const iAt = (s, V) => {
    // v descends with index; find the bracketing samples
    const n = s.v.length;
    if (V >= s.v[0]) return 0; // blocking diode: no reverse current into a weaker string
    if (V <= s.v[n - 1]) return s.i[n - 1];
    let lo = 0;
    let hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (s.v[mid] >= V) lo = mid;
      else hi = mid;
    }
    const t = (s.v[lo] - V) / (s.v[lo] - s.v[hi] || 1);
    return s.i[lo] + t * (s.i[hi] - s.i[lo]);
  };
  const v = new Float64Array(nPoints);
  const i = new Float64Array(nPoints);
  let pmp = 0;
  let vmp = 0;
  let imp = 0;
  for (let k = 0; k < nPoints; k++) {
    v[k] = (vocMax * k) / (nPoints - 1);
    for (const s of strings) i[k] += iAt(s, v[k]);
    if (v[k] * i[k] > pmp) {
      pmp = v[k] * i[k];
      vmp = v[k];
      imp = i[k];
    }
  }
  let sum = 0;
  for (const s of strings) {
    let m = 0;
    for (let k = 0; k < s.v.length; k++) m = Math.max(m, s.v[k] * s.i[k]);
    sum += m;
  }
  return { v, i, pmp, vmp, imp, sumOfStringPmp: sum, mismatchLoss: sum > 0 ? 1 - pmp / sum : 0 };
}

/**
 * Inverter AC output, PVWatts model (Dobos 2014, NREL/TP-6A20-62641):
 *   η = (η_nom/η_ref)(−0.0162 ζ − 0.0059/ζ + 0.9858),  ζ = P_dc/P_dc0,  P_dc0 = P_ac0/η_nom,  η_ref = 0.9637,
 * with clipping at the AC rating P_ac0.
 * @param {object} o
 * @param {number} o.pdc DC input power (W)
 * @param {number} o.pac0 AC rating (W)
 * @param {number} [o.etaNom=0.96] nominal efficiency (-)
 * @returns {{pac:number, efficiency:number, clipped:number}} AC power (W), efficiency (-), clipped power (W)
 */
export function inverterPvwatts({ pdc, pac0, etaNom = 0.96 }) {
  finitePos(pac0, 'pac0');
  need(etaNom > 0 && etaNom <= 1, 'etaNom must be within (0, 1]');
  if (!(pdc > 0)) return { pac: 0, efficiency: 0, clipped: 0 };
  const pdc0 = pac0 / etaNom;
  const z = pdc / pdc0;
  const eta = Math.max(0, (etaNom / 0.9637) * (-0.0162 * z - 0.0059 / z + 0.9858));
  const raw = eta * pdc;
  const pac = Math.min(pac0, raw);
  return { pac, efficiency: pac / pdc, clipped: raw - pac };
}

/**
 * Inverter AC output, Sandia model (King et al. 2007, SAND2007-5036):
 *   A = P_dco [1 + C1 (V_dc − V_dco)], B = P_so [1 + C2 (V_dc − V_dco)], C = C0 [1 + C3 (V_dc − V_dco)],
 *   P_ac = [P_aco/(A − B) − C (A − B)](P_dc − B) + C (P_dc − B)², limited to P_aco; −P_nt below start-up.
 * @param {object} o
 * @param {number} o.pdc DC power (W)
 * @param {number} o.vdc DC voltage (V)
 * @param {{Paco:number, Pdco:number, Vdco:number, Pso:number, C0:number, C1:number, C2:number, C3:number,
 *   Pnt?:number}} o.inverter CEC/Sandia inverter parameters
 * @returns {{pac:number, clipped:number}} W
 */
export function inverterSandia({ pdc, vdc, inverter }) {
  const { Paco, Pdco, Vdco, Pso, C0, C1, C2, C3, Pnt = 0 } = inverter ?? {};
  finitePos(Paco, 'inverter.Paco');
  finitePos(Pdco, 'inverter.Pdco');
  need([Vdco, Pso, C0, C1, C2, C3].every(Number.isFinite), 'inverter parameters Vdco, Pso, C0–C3 must be finite');
  const dv = vdc - Vdco;
  const A = Pdco * (1 + C1 * dv);
  const B = Pso * (1 + C2 * dv);
  const C = C0 * (1 + C3 * dv);
  if (!(pdc > B)) return { pac: -Math.abs(Pnt), clipped: 0 };
  const raw = (Paco / (A - B) - C * (A - B)) * (pdc - B) + C * (pdc - B) ** 2;
  const pac = Math.min(Paco, raw);
  return { pac, clipped: raw - pac };
}

/**
 * Inverter with a quadratic loss model: P_loss/P_rated = p0 + k1 p + k2 p², p = P_ac/P_rated (stand-by,
 * voltage-drop and ohmic terms; Schmidt & Sauer 1996), solved for P_ac given P_dc and clipped at the rating.
 * Also see `europeanEfficiency`.
 * @param {object} o
 * @param {number} o.pdc DC power (W)
 * @param {number} o.pacRated AC rating (W)
 * @param {number} o.p0 no-load loss fraction (-), typically 0.002–0.01
 * @param {number} o.k1 linear loss coefficient (-), typically 0.005–0.02
 * @param {number} o.k2 quadratic loss coefficient (-), typically 0.01–0.03
 * @returns {{pac:number, efficiency:number, clipped:number}}
 */
export function inverterQuadratic({ pdc, pacRated, p0, k1, k2 }) {
  finitePos(pacRated, 'pacRated');
  need(p0 >= 0 && k1 >= 0 && k2 >= 0, 'p0, k1 and k2 must be >= 0');
  const d = pdc / pacRated;
  if (!(d > p0)) return { pac: 0, efficiency: 0, clipped: 0 };
  // d = p + p0 + k1 p + k2 p²
  const p = k2 > 0 ? (-(1 + k1) + Math.sqrt((1 + k1) ** 2 + 4 * k2 * (d - p0))) / (2 * k2) : (d - p0) / (1 + k1);
  const pc = Math.min(1, p);
  return { pac: pc * pacRated, efficiency: (pc * pacRated) / pdc, clipped: (p - pc) * pacRated };
}

/**
 * European weighted efficiency (as defined for inverter rating, cf. EN 50530): η_EU = 0.03 η5 + 0.06 η10 + 0.13 η20 + 0.10 η30 + 0.48 η50 + 0.20 η100.
 * @param {(loadFraction:number)=>number} efficiencyAt inverter efficiency at a fraction of rated AC power
 * @returns {number} η_EU (-)
 */
export function europeanEfficiency(efficiencyAt) {
  need(typeof efficiencyAt === 'function', 'efficiencyAt must be a function of load fraction');
  const w = [[0.05, 0.03], [0.1, 0.06], [0.2, 0.13], [0.3, 0.1], [0.5, 0.48], [1, 0.2]];
  return w.reduce((s, [l, k]) => s + k * efficiencyAt(l), 0);
}

/**
 * Mean plane-of-array factor of a module whose tilt oscillates with the waves (floating PV): the transposed
 * irradiance is averaged over the tilt distribution β + δ, with δ sinusoidal in time (regular motion,
 * 'harmonic') or Gaussian (irregular sea, `amplitudeDeg` then being the standard deviation).
 * Rotation is about the horizontal axis normal to the surface azimuth (pitch of the tilt plane).
 * @param {object} o transposition inputs as for `planeOfArray`, plus
 * @param {number} o.amplitudeDeg tilt oscillation amplitude or standard deviation (deg)
 * @param {'harmonic'|'gaussian'} [o.distribution='harmonic']
 * @returns {{mean:number, still:number, factor:number}} mean POA (W/m²), POA at rest (W/m²), their ratio (-)
 */
export function waveTiltEffect({ amplitudeDeg, distribution = 'harmonic', ...poa }) {
  need(amplitudeDeg >= 0, 'amplitudeDeg must be >= 0');
  const still = planeOfArray(poa).total;
  if (amplitudeDeg === 0) return { mean: still, still, factor: 1 };
  const n = 24;
  let s = 0;
  let w = 0;
  for (let k = 0; k < n; k++) {
    let d;
    let wk = 1;
    if (distribution === 'harmonic') d = amplitudeDeg * Math.sin((2 * Math.PI * (k + 0.5)) / n);
    else {
      need(distribution === 'gaussian', "distribution must be 'harmonic' or 'gaussian'");
      const z = -3 + (6 * (k + 0.5)) / n;
      d = amplitudeDeg * z;
      wk = Math.exp(-0.5 * z * z);
    }
    let tilt = poa.tiltDeg + d;
    let az = poa.surfaceAzimuthDeg;
    if (tilt < 0) {
      tilt = -tilt;
      az += 180;
    }
    s += wk * planeOfArray({ ...poa, tiltDeg: Math.min(180, tilt), surfaceAzimuthDeg: az }).total;
    w += wk;
  }
  const mean = s / w;
  return { mean, still, factor: still > 0 ? mean / still : 1 };
}

/**
 * Energy yield of a PV system from a time series of global horizontal irradiance, ambient temperature and wind.
 * Chain per step: solar position → Erbs decomposition → transposition (Hay–Davies by default) with surface
 * albedo → incidence-angle modifier on the beam and circumsolar part → soiling → optional bifacial rear gain →
 * cell temperature (Faiman by default) → PVWatts DC power → DC losses → inverter (PVWatts model, clipping at
 * P_dc0 / dcAcRatio) → AC losses → availability and degradation.
 * Time stamps should mark the middle of each averaging interval.
 * @param {object} o
 * @param {number[]|Float64Array} o.times instants (ms since epoch, UTC)
 * @param {number[]|Float64Array} o.ghi global horizontal irradiance (W/m²)
 * @param {number[]|Float64Array|number} [o.tempAirC=20] ambient temperature (°C)
 * @param {number[]|Float64Array|number} [o.windSpeed=1] wind speed (m/s)
 * @param {number} o.lat latitude (deg)
 * @param {number} o.lon longitude (deg, east positive)
 * @param {number} o.tiltDeg module tilt (deg)
 * @param {number} [o.surfaceAzimuthDeg] module azimuth (deg); default 180 in the northern hemisphere, 0 south
 * @param {number} o.pdc0 array DC rating at STC (W)
 * @param {number} [o.dtHours=1] step length (h)
 * @param {number} [o.albedo=0.08] surface albedo (-); 0.06–0.10 for water, 0.2 for land
 * @param {'isotropic'|'haydavies'} [o.transposition='haydavies']
 * @param {'physical'|'ashrae'|'none'} [o.iam='physical']
 * @param {number} [o.soiling=0.02] soiling loss (-)
 * @param {number} [o.bifaciality=0] rear/front efficiency ratio (0 = monofacial)
 * @param {number} [o.groundShadeFactor=0.3] shaded fraction of the ground seen by the rear side (-)
 * @param {'faiman'|'noct'|'sandia'} [o.thermal='faiman']
 * @param {object} [o.thermalParams] coefficients for the chosen thermal model (u0, u1, noctC, a, b, deltaT,
 *   waterCoolingOffsetC)
 * @param {number} [o.gammaPdc=-0.0035] power temperature coefficient (1/K)
 * @param {number} [o.dcLoss=0.02] DC wiring/mismatch loss (-)
 * @param {number} [o.dcAcRatio=1.2] DC/AC ratio (-)
 * @param {number} [o.inverterEfficiency=0.97] nominal inverter efficiency (-)
 * @param {number} [o.acLoss=0.015] transformer + AC cable loss (-)
 * @param {number} [o.availability=0.99] (-)
 * @param {number} [o.degradationPerYear=0.005] annual degradation rate (-)
 * @param {number} [o.year=0] years since commissioning (degradation factor (1 − rate)^year)
 * @param {number} [o.tiltAmplitudeDeg=0] wave-induced tilt oscillation amplitude (deg), harmonic
 * @param {number} [o.solarConstant=1361]
 * @returns {{pac:Float64Array, pdc:Float64Array, poa:Float64Array, tempCellC:Float64Array, energyWh:number,
 *   dcEnergyWh:number, poaWhPerM2:number, ghiWhPerM2:number, clippedWh:number, specificYield:number,
 *   performanceRatio:number, capacityFactor:number, monthly:{month:string, energyWh:number,
 *   poaWhPerM2:number}[], daily:{day:string, energyWh:number}[]}} specificYield in Wh per W_dc (= kWh/kWp)
 */
export function simulatePv({
  times,
  ghi,
  tempAirC = 20,
  windSpeed = 1,
  lat,
  lon,
  tiltDeg,
  surfaceAzimuthDeg,
  pdc0,
  dtHours = 1,
  albedo = 0.08,
  transposition = 'haydavies',
  iam = 'physical',
  soiling = 0.02,
  bifaciality = 0,
  groundShadeFactor = 0.3,
  thermal = 'faiman',
  thermalParams = {},
  gammaPdc = -0.0035,
  dcLoss = 0.02,
  dcAcRatio = 1.2,
  inverterEfficiency = 0.97,
  acLoss = 0.015,
  availability = 0.99,
  degradationPerYear = 0.005,
  year = 0,
  tiltAmplitudeDeg = 0,
  solarConstant = 1361,
}) {
  need(times && ghi && times.length === ghi.length && times.length > 0, 'times and ghi must be series of equal non-zero length');
  finitePos(pdc0, 'pdc0');
  finitePos(dtHours, 'dtHours');
  finitePos(dcAcRatio, 'dcAcRatio');
  for (const [v, name] of [[soiling, 'soiling'], [dcLoss, 'dcLoss'], [acLoss, 'acLoss'], [availability, 'availability'], [degradationPerYear, 'degradationPerYear']]) {
    need(v >= 0 && v <= 1, `${name} must be within [0, 1]`);
  }
  need(['faiman', 'noct', 'sandia'].includes(thermal), "thermal must be 'faiman', 'noct' or 'sandia'");
  need(['physical', 'ashrae', 'none'].includes(iam), "iam must be 'physical', 'ashrae' or 'none'");
  const az = surfaceAzimuthDeg ?? (lat >= 0 ? 180 : 0);
  const n = times.length;
  const at = (s, i) => (typeof s === 'number' ? s : s[i]);
  const pac0 = pdc0 / dcAcRatio;
  const life = (1 - degradationPerYear) ** year;
  const pac = new Float64Array(n);
  const pdc = new Float64Array(n);
  const poaS = new Float64Array(n);
  const tc = new Float64Array(n);
  let e = 0;
  let eDc = 0;
  let hPoa = 0;
  let hGhi = 0;
  let clip = 0;
  const monthly = new Map();
  const daily = new Map();
  for (let i = 0; i < n; i++) {
    const g = Math.max(0, ghi[i]);
    const Ta = at(tempAirC, i);
    const ws = at(windSpeed, i);
    hGhi += g * dtHours;
    let pOut = 0;
    let poaTot = 0;
    let Tcell = Ta;
    if (g > 0) {
      const sp = solarPosition({ time: times[i], lat, lon });
      if (sp.zenithDeg < 90) {
        const dniExtra = solarConstant / (sp.distanceAU * sp.distanceAU);
        const dec = erbsDecomposition({ ghi: g, zenithDeg: sp.zenithDeg, dniExtra });
        const geo = { dni: dec.dni, dhi: dec.dhi, ghi: g, zenithDeg: sp.zenithDeg, azimuthDeg: sp.azimuthDeg, tiltDeg, surfaceAzimuthDeg: az, albedo, model: transposition, dniExtra };
        const pa = planeOfArray(geo);
        const f = iam === 'none' ? 1 : iam === 'ashrae' ? iamAshrae(pa.aoiDeg) : iamPhysical(pa.aoiDeg);
        let eff = (pa.beam + pa.circumsolar) * f + (pa.skyDiffuse - pa.circumsolar) + pa.groundReflected;
        poaTot = pa.total;
        if (tiltAmplitudeDeg > 0) {
          const wf = waveTiltEffect({ ...geo, amplitudeDeg: tiltAmplitudeDeg }).factor;
          eff *= wf;
          poaTot *= wf;
        }
        eff *= 1 - soiling;
        if (bifaciality > 0) eff += bifaciality * rearIrradiance({ ghi: g, dhi: dec.dhi, tiltDeg, albedo, groundShadeFactor });
        if (thermal === 'faiman') Tcell = cellTempFaiman({ poa: poaTot, tempAirC: Ta, windSpeed: ws, ...thermalParams });
        else if (thermal === 'noct') Tcell = cellTempNoct({ poa: poaTot, tempAirC: Ta, ...thermalParams });
        else Tcell = cellTempSandia({ poa: poaTot, tempAirC: Ta, windSpeed: ws, ...thermalParams });
        const dc = pvwattsDc({ pdc0, poa: eff, tempCellC: Tcell, gammaPdc }) * (1 - dcLoss) * life;
        const inv = inverterPvwatts({ pdc: dc, pac0, etaNom: inverterEfficiency });
        pdc[i] = dc;
        clip += inv.clipped * dtHours;
        pOut = inv.pac * (1 - acLoss) * availability;
      }
    }
    pac[i] = pOut;
    poaS[i] = poaTot;
    tc[i] = Tcell;
    e += pOut * dtHours;
    eDc += pdc[i] * dtHours;
    hPoa += poaTot * dtHours;
    const iso = new Date(times[i]).toISOString();
    const mk = iso.slice(0, 7);
    const mo = monthly.get(mk) ?? { month: mk, energyWh: 0, poaWhPerM2: 0 };
    mo.energyWh += pOut * dtHours;
    mo.poaWhPerM2 += poaTot * dtHours;
    monthly.set(mk, mo);
    const dk = iso.slice(0, 10);
    daily.set(dk, (daily.get(dk) ?? 0) + pOut * dtHours);
  }
  return {
    pac,
    pdc,
    poa: poaS,
    tempCellC: tc,
    energyWh: e,
    dcEnergyWh: eDc,
    poaWhPerM2: hPoa,
    ghiWhPerM2: hGhi,
    clippedWh: clip,
    specificYield: e / pdc0,
    performanceRatio: hPoa > 0 ? e / pdc0 / (hPoa / 1000) : NaN,
    capacityFactor: e / (pac0 * n * dtHours),
    monthly: [...monthly.values()],
    daily: [...daily.entries()].map(([day, energyWh]) => ({ day, energyWh })),
  };
}

/**
 * Synthetic hourly year from monthly climatology: every day of a month receives the month's mean clearness
 * index, GHI(t) = K_T · E0n(t) cos z(t) with K_T = H̄ / H̄0 (monthly mean daily global over extraterrestrial
 * horizontal insolation, Liu & Jordan 1960), so each month reproduces its mean daily insolation exactly.
 * Day-to-day variability is not represented, which slightly overstates the yield of clipping-limited systems
 * and understates diffuse fractions on clear days. Temperature: monthly mean plus an optional sinusoidal diurnal
 * swing peaking 3 h after solar noon.
 * @param {object} o
 * @param {number[]} o.monthlyGhi 12 values of mean daily global horizontal insolation (Wh/m²/day)
 * @param {number[]} [o.monthlyTempC] 12 monthly mean temperatures (°C); default 20
 * @param {number} [o.diurnalRangeC=0] daily temperature range (°C)
 * @param {number} o.lat (deg)
 * @param {number} o.lon (deg, east positive)
 * @param {number} [o.year=2023] calendar year (a non-leap year gives 8760 h)
 * @param {number} [o.solarConstant=1361]
 * @returns {{times:Float64Array, ghi:Float64Array, tempAirC:Float64Array, clearness:number[]}} hourly series
 *   with time stamps at mid-hour; clearness = K_T per month
 */
export function syntheticYear({ monthlyGhi, monthlyTempC, diurnalRangeC = 0, lat, lon, year = 2023, solarConstant = 1361 }) {
  need(Array.isArray(monthlyGhi) && monthlyGhi.length === 12, 'monthlyGhi must have 12 values');
  need(!monthlyTempC || monthlyTempC.length === 12, 'monthlyTempC must have 12 values');
  const t0 = Date.UTC(year, 0, 1);
  const t1 = Date.UTC(year + 1, 0, 1);
  const n = Math.round((t1 - t0) / 3600000);
  const times = new Float64Array(n);
  const g0 = new Float64Array(n);
  const month = new Uint8Array(n);
  const h0 = new Float64Array(12);
  const days = new Float64Array(12);
  for (let i = 0; i < n; i++) {
    times[i] = t0 + (i + 0.5) * 3600000;
    const m = new Date(times[i]).getUTCMonth();
    month[i] = m;
    const sp = solarPosition({ time: times[i], lat, lon });
    g0[i] = sp.zenithDeg < 90 ? (solarConstant / sp.distanceAU ** 2) * cosd(sp.zenithDeg) : 0;
    h0[m] += g0[i];
    days[m] += 1 / 24;
  }
  const clearness = Array.from(h0, (h, m) => {
    need(monthlyGhi[m] >= 0, 'monthlyGhi values must be >= 0');
    return h > 0 ? (monthlyGhi[m] * days[m]) / h : 0;
  });
  const ghi = new Float64Array(n);
  const tempAirC = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    ghi[i] = clearness[month[i]] * g0[i];
    const localSolarHour = mod((times[i] % 86400000) / 3600000 + lon / 15, 24);
    tempAirC[i] = (monthlyTempC ? monthlyTempC[month[i]] : 20) + 0.5 * diurnalRangeC * Math.cos(((localSolarHour - 15) * Math.PI) / 12);
  }
  return { times, ghi, tempAirC, clearness };
}

/**
 * Tilt (and optionally azimuth) maximising the annual energy of `simulatePv` by grid search with golden-section
 * refinement of the tilt.
 * @param {object} o
 * @param {object} o.system options for `simulatePv` (tiltDeg and surfaceAzimuthDeg are overridden)
 * @param {number} [o.tiltMin=0] (deg)
 * @param {number} [o.tiltMax=60] (deg)
 * @param {number} [o.tiltStep=5] coarse step (deg)
 * @param {number[]} [o.azimuthsDeg] candidate azimuths (deg); default the equator-facing azimuth only
 * @returns {{tiltDeg:number, surfaceAzimuthDeg:number, energyWh:number, table:{tiltDeg:number,
 *   surfaceAzimuthDeg:number, energyWh:number}[]}}
 */
export function optimalTilt({ system, tiltMin = 0, tiltMax = 60, tiltStep = 5, azimuthsDeg }) {
  need(system && Number.isFinite(system.lat), 'system must contain the simulatePv options');
  need(tiltMax > tiltMin && tiltStep > 0, 'tilt range and step must be positive');
  const azs = azimuthsDeg ?? [system.lat >= 0 ? 180 : 0];
  const run = (tiltDeg, surfaceAzimuthDeg) => simulatePv({ ...system, tiltDeg, surfaceAzimuthDeg }).energyWh;
  const table = [];
  let best = null;
  for (const az of azs) {
    for (let t = tiltMin; t <= tiltMax + 1e-9; t += tiltStep) {
      const energyWh = run(t, az);
      const row = { tiltDeg: t, surfaceAzimuthDeg: az, energyWh };
      table.push(row);
      if (!best || energyWh > best.energyWh) best = row;
    }
  }
  let lo = Math.max(tiltMin, best.tiltDeg - tiltStep);
  let hi = Math.min(tiltMax, best.tiltDeg + tiltStep);
  const gr = (Math.sqrt(5) - 1) / 2;
  for (let i = 0; i < 12; i++) {
    const c = hi - gr * (hi - lo);
    const d = lo + gr * (hi - lo);
    if (run(c, best.surfaceAzimuthDeg) > run(d, best.surfaceAzimuthDeg)) hi = d;
    else lo = c;
  }
  const tiltDeg = 0.5 * (lo + hi);
  const energyWh = run(tiltDeg, best.surfaceAzimuthDeg);
  return energyWh >= best.energyWh
    ? { tiltDeg, surfaceAzimuthDeg: best.surfaceAzimuthDeg, energyWh, table }
    : { ...best, table };
}
