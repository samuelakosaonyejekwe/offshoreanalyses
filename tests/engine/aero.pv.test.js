import test from 'node:test';
import assert from 'node:assert/strict';
import {
  solarPosition, sunriseSunset, extraterrestrialIrradiance, airMass, clearSkyHaurwitz, clearSkyIneichen,
  erbsDecomposition, erbsDiffuseFraction, angleOfIncidence, planeOfArray, rearIrradiance, iamAshrae,
  iamPhysical, cellTempNoct, cellTempFaiman, cellTempSandia, lambertW, singleDiodeCurrent, singleDiodeVoltage,
  singleDiodeCurve, singleDiodeAt, extractSingleDiode, pvwattsDc, seriesStringCurve, parallelStrings,
  inverterPvwatts, inverterSandia, inverterQuadratic, europeanEfficiency, waveTiltEffect, simulatePv,
  syntheticYear, optimalTilt,
} from '../../src/engine/aero/pv.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);
const absClose = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} got ${a}, expected ${b}`);

test('solar position: equinox noon on the equator is overhead; noon zenith = |latitude − declination|', () => {
  // March equinox 2024: 20 March 03:06 UTC. Local solar noon at lon 0 is 12:00 UTC − equation of time (≈ −7.4 min)
  const noonGuess = Date.UTC(2024, 2, 20, 12, 0);
  const ss = sunriseSunset({ time: noonGuess, lat: 0, lon: 0 });
  const p = solarPosition({ time: ss.solarNoon, lat: 0, lon: 0 });
  absClose(p.zenithDeg, Math.abs(p.declinationDeg), 0.01);
  assert.ok(p.zenithDeg < 0.3, `zenith ${p.zenithDeg}`);
  absClose(p.hourAngleDeg, 0, 0.01);
  absClose(p.equationOfTimeMin, -7.4, 0.3);
  const q = solarPosition({ time: ss.solarNoon, lat: 40, lon: 0 });
  absClose(q.zenithDeg, 40 - q.declinationDeg, 0.01);
  absClose(q.azimuthDeg, 180, 0.2);
  // morning sun in the east, afternoon sun in the west; east longitude advances solar time by 4 min/deg
  assert.ok(solarPosition({ time: ss.solarNoon - 3 * 3600e3, lat: 40, lon: 0 }).azimuthDeg < 180);
  assert.ok(solarPosition({ time: ss.solarNoon + 3 * 3600e3, lat: 40, lon: 0 }).azimuthDeg > 180);
  absClose(solarPosition({ time: ss.solarNoon - 3600e3, lat: 40, lon: 15 }).hourAngleDeg, 0, 0.01);
});

test('solar declination and equation of time through the year (Astronomical Almanac values)', () => {
  absClose(solarPosition({ time: Date.UTC(2023, 5, 21, 12), lat: 0, lon: 0 }).declinationDeg, 23.44, 0.01);
  absClose(solarPosition({ time: Date.UTC(2023, 11, 22, 3), lat: 0, lon: 0 }).declinationDeg, -23.44, 0.01);
  // extremes of the equation of time: about +16.4 min on 3 November and −14.2 min on 11 February
  absClose(solarPosition({ time: Date.UTC(2023, 10, 3, 12), lat: 0, lon: 0 }).equationOfTimeMin, 16.45, 0.1);
  absClose(solarPosition({ time: Date.UTC(2023, 1, 11, 12), lat: 0, lon: 0 }).equationOfTimeMin, -14.2, 0.1);
  // Earth–Sun distance: perihelion ≈ 0.9833 AU (early January), aphelion ≈ 1.0167 AU (early July)
  absClose(solarPosition({ time: Date.UTC(2023, 0, 4, 12), lat: 0, lon: 0 }).distanceAU, 0.9833, 2e-4);
  absClose(solarPosition({ time: Date.UTC(2023, 6, 6, 12), lat: 0, lon: 0 }).distanceAU, 1.0167, 2e-4);
  close(extraterrestrialIrradiance(Date.UTC(2023, 0, 4, 12)), 1361 / 0.9833 ** 2, 5e-4);
  // refraction lifts the apparent sun by about 0.5° at the horizon
  const hz = solarPosition({ time: sunriseSunset({ time: Date.UTC(2023, 2, 21), lat: 0, lon: 0, horizonZenithDeg: 90 }).sunrise, lat: 0, lon: 0 });
  absClose(hz.elevationDeg, 0, 0.05);
  absClose(hz.apparentElevationDeg - hz.elevationDeg, 0.48, 0.05);
});

test('day length at 60°N on the June solstice', () => {
  const t = Date.UTC(2023, 5, 21, 12);
  const decl = solarPosition({ time: t, lat: 60, lon: 0 }).declinationDeg;
  // geometric: (2/15) acos(−tan φ tan δ) = (2/15)·acos(−tan 60° tan 23.44°) = 18.49 h
  const geo = sunriseSunset({ time: t, lat: 60, lon: 0, horizonZenithDeg: 90 });
  close(geo.dayLengthHours, (2 / 15) * (Math.acos(-Math.tan(Math.PI / 3) * Math.tan((decl * Math.PI) / 180)) * 180) / Math.PI, 1e-9);
  absClose(geo.dayLengthHours, 18.49, 0.02);
  // with refraction and the solar semi-diameter (sunrise/sunset tables for 60°N: about 18 h 50 min)
  const std = sunriseSunset({ time: t, lat: 60, lon: 0 });
  absClose(std.dayLengthHours, 18.85, 0.08);
  // the sun is on the standard horizon at the computed sunrise
  absClose(solarPosition({ time: std.sunrise, lat: 60, lon: 0 }).zenithDeg, 90.833, 0.05);
  // equator at equinox: 12 h geometric; polar day and night
  absClose(sunriseSunset({ time: Date.UTC(2023, 2, 20, 12), lat: 0, lon: 0, horizonZenithDeg: 90 }).dayLengthHours, 12, 0.01);
  assert.equal(sunriseSunset({ time: t, lat: 80, lon: 0 }).polar, 'day');
  assert.equal(sunriseSunset({ time: Date.UTC(2023, 11, 21, 12), lat: 80, lon: 0 }).polar, 'night');
});

test('air mass (Kasten–Young) and clear-sky models by hand', () => {
  close(airMass(0).relative, 1 / (1 + 0.50572 * 96.07995 ** -1.6364), 1e-12);
  absClose(airMass(0).relative, 1.0, 5e-4);
  // z = 60°: 36.07995^−1.6364 = 0.002829; 1/(0.5 + 0.50572·0.002829) = 1.99429
  absClose(airMass(60).relative, 1.99429, 5e-5);
  absClose(airMass(90).relative, 37.92, 0.05);
  close(airMass(60, 90000).absolute, (airMass(60).relative * 90000) / 101325, 1e-12);
  // Haurwitz, overhead sun: 1098 exp(−0.059) = 1035.1 W/m²
  absClose(clearSkyHaurwitz(0), 1035.1, 0.1);
  assert.equal(clearSkyHaurwitz(95), 0);
  // Ineichen–Perez at sea level, overhead sun, TL = 3: GHI = 0.868·1361·exp(−0.0387·AM·3), AM = 0.99971
  const c = clearSkyIneichen({ zenithDeg: 0, linkeTurbidity: 3 });
  close(c.ghi, 0.868 * 1361 * Math.exp(-0.0387 * airMass(0).relative * 3), 1e-12);
  absClose(c.ghi, 1051.9, 0.5);
  // DNI: b = 0.827, 0.827·1361·exp(−0.09·AM·2) = 940.1 (first limb governs)
  absClose(c.dni, 0.827 * 1361 * Math.exp(-0.09 * airMass(0).relative * 2), 1e-9);
  close(c.dhi, c.ghi - c.dni, 1e-12);
  // more turbidity → less beam, more diffuse
  const hazy = clearSkyIneichen({ zenithDeg: 30, linkeTurbidity: 5 });
  const clean = clearSkyIneichen({ zenithDeg: 30, linkeTurbidity: 2.5 });
  assert.ok(hazy.dni < clean.dni && hazy.dhi > clean.dhi && hazy.ghi < clean.ghi);
});

test('Erbs diffuse fraction at kt = 0.1, 0.5, 0.9', () => {
  close(erbsDiffuseFraction(0.1), 0.991, 1e-12);
  // 0.9511 − 0.0802 + 1.097 − 2.07975 + 0.771 = 0.65915
  close(erbsDiffuseFraction(0.5), 0.65915, 1e-10);
  close(erbsDiffuseFraction(0.9), 0.165, 1e-12);
  // continuity at the branch points is within the correlation's published mismatch
  absClose(erbsDiffuseFraction(0.22), erbsDiffuseFraction(0.2200001), 2e-3);
  // decomposition closes: GHI = DNI cos z + DHI
  const d = erbsDecomposition({ ghi: 600, zenithDeg: 40, dniExtra: 1361 });
  close(d.kt, 600 / (1361 * Math.cos((40 * Math.PI) / 180)), 1e-12);
  close(d.dni * Math.cos((40 * Math.PI) / 180) + d.dhi, 600, 1e-12);
  close(d.dhi, erbsDiffuseFraction(d.kt) * 600, 1e-12);
  assert.equal(erbsDecomposition({ ghi: 50, zenithDeg: 89 }).dni, 0);
});

test('transposition by hand (isotropic) and Hay–Davies limits', () => {
  // sun at zenith 30° due south, module tilted 30° south: normal incidence
  const o = { dni: 800, dhi: 100, zenithDeg: 30, azimuthDeg: 180, tiltDeg: 30, surfaceAzimuthDeg: 180, albedo: 0.2 };
  const iso = planeOfArray({ ...o, model: 'isotropic' });
  absClose(iso.aoiDeg, 0, 1e-6);
  close(iso.beam, 800, 1e-12);
  // sky: 100 (1 + cos 30°)/2 = 93.301; GHI = 800 cos 30° + 100 = 792.820; ground: 792.820·0.2·0.0669873 = 10.62178
  close(iso.skyDiffuse, 93.30127, 1e-7);
  close(iso.groundReflected, 10.62178, 1e-6);
  close(iso.total, 903.92305, 1e-7);
  // Hay–Davies: A_i = 800/1361, R_b = 1/cos 30°
  const hd = planeOfArray({ ...o, model: 'haydavies' });
  const Ai = 800 / 1361;
  close(hd.skyDiffuse, 100 * (Ai / Math.cos(Math.PI / 6) + (1 - Ai) * (1 + Math.cos(Math.PI / 6)) / 2), 1e-12);
  assert.ok(hd.total > iso.total);
  // no beam → both models coincide; horizontal plane returns GHI
  close(planeOfArray({ ...o, dni: 0, model: 'haydavies' }).total, planeOfArray({ ...o, dni: 0, model: 'isotropic' }).total, 1e-12);
  close(planeOfArray({ ...o, tiltDeg: 0, model: 'isotropic' }).total, 792.8203, 1e-7);
  close(planeOfArray({ ...o, tiltDeg: 0, model: 'haydavies' }).total, 792.8203, 1e-7);
  // sun behind the module: no beam
  assert.equal(planeOfArray({ ...o, tiltDeg: 90, surfaceAzimuthDeg: 0 }).beam, 0);
  close(angleOfIncidence({ zenithDeg: 60, azimuthDeg: 180, tiltDeg: 0, surfaceAzimuthDeg: 180 }).aoiDeg, 60, 1e-12);
  // rear side over water (albedo 0.08), tilt 15°: 0.08·800·(1 + cos15°)/2·0.7 + 120·(1 − cos15°)/2 = 46.08
  close(rearIrradiance({ ghi: 800, dhi: 120, tiltDeg: 15, albedo: 0.08, groundShadeFactor: 0.3 }), 0.08 * 800 * 0.9829629 * 0.7 + 120 * 0.0170371, 1e-6);
});

test('incidence-angle modifiers', () => {
  close(iamAshrae(60), 0.95, 1e-12);
  assert.equal(iamAshrae(0), 1);
  assert.equal(iamAshrae(90), 0);
  close(iamPhysical(0), 1, 1e-12);
  // Fresnel/De Soto with n = 1.526: about 0.999 at 30°, 0.95 at 60°, 0.57 at 80° (monotonic decrease)
  absClose(iamPhysical(30), 0.9985, 2e-3);
  absClose(iamPhysical(60), 0.948, 5e-3);
  let prev = 1;
  for (let a = 5; a < 90; a += 5) {
    const v = iamPhysical(a);
    assert.ok(v <= prev + 1e-12 && v >= 0);
    prev = v;
  }
  assert.equal(iamPhysical(90), 0);
});

test('cell temperature models by hand', () => {
  close(cellTempNoct({ poa: 800, tempAirC: 20, noctC: 45 }), 45, 1e-12);
  // Faiman: 25 + 1000/(25 + 6.84·1) = 56.407
  close(cellTempFaiman({ poa: 1000, tempAirC: 25, windSpeed: 1 }), 25 + 1000 / 31.84, 1e-12);
  // better-cooled floating array (higher U-value) runs cooler; explicit water-cooling offset shifts the result
  assert.ok(cellTempFaiman({ poa: 1000, tempAirC: 25, windSpeed: 1, u0: 40, u1: 6 }) < cellTempFaiman({ poa: 1000, tempAirC: 25, windSpeed: 1 }));
  close(cellTempFaiman({ poa: 1000, tempAirC: 25, windSpeed: 1, waterCoolingOffsetC: -3 }), 22 + 1000 / 31.84, 1e-12);
  // Sandia open rack glass/glass: 1000 exp(−3.47 − 0.0594) + 25 + 3 = 57.32
  close(cellTempSandia({ poa: 1000, tempAirC: 25, windSpeed: 1 }), 1000 * Math.exp(-3.5294) + 28, 1e-12);
  absClose(cellTempSandia({ poa: 1000, tempAirC: 25, windSpeed: 1 }), 57.32, 0.02);
});

test('Lambert W: W(x) exp(W(x)) = x and special values', () => {
  for (const x of [-0.367, -0.3, -0.1, 1e-8, 0.5, 1, 10, 1e3, 1e10, 1e100]) {
    const w = lambertW(x);
    close(w * Math.exp(w), x, 1e-12, `x = ${x}`);
  }
  close(lambertW(Math.E), 1, 1e-14);
  close(lambertW(0), 0, 1e-15);
  absClose(lambertW(-1 / Math.E), -1, 1e-6);
  close(lambertW(1), 0.5671432904097838, 1e-14);
  assert.throws(() => lambertW(-1), RangeError);
});

const KC200 = { Isc: 8.21, Voc: 32.9, Imp: 7.61, Vmp: 26.3, Ns: 54, alphaIsc: 3.18e-3, betaVoc: -0.123 };

test('single-diode extraction reproduces the datasheet: Isc, Voc, MPP and the Voc temperature coefficient', () => {
  const p = extractSingleDiode(KC200);
  assert.ok(p.I0 > 0 && p.Rs > 0 && p.Rsh > 0 && p.n > 0.7 && p.n < 2);
  const c = singleDiodeCurve(p);
  close(c.isc, 8.21, 1e-9);
  close(c.voc, 32.9, 1e-9);
  close(c.pmp, 7.61 * 26.3, 1e-9);
  close(c.vmp, 26.3, 1e-5);
  close(c.imp, 7.61, 1e-5);
  // the datasheet MPP lies on the curve
  close(singleDiodeCurrent(p, 26.3), 7.61, 1e-9);
  // temperature coefficient of Voc over +10 K
  const hot = singleDiodeCurve(singleDiodeAt(p, { G: 1000, tempCellC: 35 }));
  close((hot.voc - 32.9) / 10, -0.123, 1e-6);
  close(hot.isc - c.isc, 10 * 3.18e-3, 2e-3);
  // the implicit diode equation is satisfied along the curve, and V(I) inverts I(V)
  for (const V of [0, 10, 20, 26.3, 30, 32]) {
    const I = singleDiodeCurrent(p, V);
    const res = p.IL - p.I0 * Math.expm1((V + I * p.Rs) / p.a) - (V + I * p.Rs) / p.Rsh - I;
    absClose(res, 0, 1e-9, `V = ${V}`);
    close(singleDiodeVoltage(p, I), V, 1e-8);
  }
  // irradiance translation: Isc ∝ G, Voc falls logarithmically
  const low = singleDiodeCurve(singleDiodeAt(p, { G: 200, tempCellC: 25 }));
  close(low.isc, 0.2 * 8.21, 2e-3);
  assert.ok(low.voc < 32.9 && low.voc > 29);
  // fixed ideality instead of the temperature condition
  const q = extractSingleDiode({ ...KC200, betaVoc: undefined, n: 1.3 });
  const cq = singleDiodeCurve(q);
  close(cq.isc, 8.21, 1e-9);
  close(cq.voc, 32.9, 1e-9);
  close(cq.pmp, 7.61 * 26.3, 1e-8);
  close(q.n, 1.3, 1e-12);
  assert.throws(() => extractSingleDiode({ ...KC200, Vmp: 40 }), RangeError);
});

test('PVWatts DC model and inverter models by hand', () => {
  // 5000 W · 0.8 · (1 − 0.004·20) = 3680 W
  close(pvwattsDc({ pdc0: 5000, poa: 800, tempCellC: 45, gammaPdc: -0.004 }), 3680, 1e-12);
  // PVWatts inverter at rated DC input: η = η_nom (−0.0162 − 0.0059 + 0.9858)/0.9637 = η_nom
  const inv = inverterPvwatts({ pdc: 4000 / 0.96, pac0: 4000, etaNom: 0.96 });
  close(inv.efficiency, 0.96, 1e-12);
  close(inv.pac, 4000, 1e-12);
  // half load: 0.96/0.9637 (−0.0081 − 0.0118 + 0.9858) = 0.96220
  close(inverterPvwatts({ pdc: 2000 / 0.96, pac0: 4000, etaNom: 0.96 }).efficiency, (0.96 / 0.9637) * 0.9659, 1e-12);
  const clip = inverterPvwatts({ pdc: 6000, pac0: 4000, etaNom: 0.96 });
  assert.equal(clip.pac, 4000);
  assert.ok(clip.clipped > 1500);
  // Sandia model: rated DC power at nominal voltage returns Paco exactly; below start-up power → night tare
  const sandia = { Paco: 4000, Pdco: 4150, Vdco: 400, Pso: 20, C0: -1e-5, C1: 2e-5, C2: 1e-3, C3: 5e-4, Pnt: 1 };
  close(inverterSandia({ pdc: 4150, vdc: 400, inverter: sandia }).pac, 4000, 1e-12);
  assert.equal(inverterSandia({ pdc: 10, vdc: 400, inverter: sandia }).pac, -1);
  assert.ok(inverterSandia({ pdc: 2000, vdc: 400, inverter: sandia }).pac < 2000);
  // quadratic loss model: DC = AC + losses
  const q = inverterQuadratic({ pdc: 3000, pacRated: 4000, p0: 0.005, k1: 0.01, k2: 0.02 });
  const pu = q.pac / 4000;
  close(3000 / 4000, pu + 0.005 + 0.01 * pu + 0.02 * pu * pu, 1e-12);
  close(europeanEfficiency(() => 0.97), 0.97, 1e-12);
  const eu = europeanEfficiency((l) => inverterQuadratic({ pdc: 4000 * (l + 0.005 + 0.01 * l + 0.02 * l * l), pacRated: 4000, p0: 0.005, k1: 0.01, k2: 0.02 }).efficiency);
  assert.ok(eu > 0.93 && eu < 0.98);
});

test('strings: uniform irradiance adds up; partial shading with bypass diodes; parallel strings', () => {
  const p = { ...extractSingleDiode(KC200) };
  const one = singleDiodeCurve(p);
  const groups = (shade) => Array.from({ length: 6 }, (_, i) => ({ params: p, fraction: 1 / 3, G: i === 0 ? shade : 1000, tempCellC: 25 }));
  const uni = seriesStringCurve({ groups: groups(1000), nPoints: 800 });
  close(uni.pmp, 2 * one.pmp, 1e-4);
  close(uni.voc, 2 * one.voc, 1e-9);
  close(uni.isc, one.isc, 1e-6);
  absClose(uni.mismatchLoss, 0, 1e-4);
  // one of six substrings at 20 % irradiance: the bypass diode lets the other five run near their MPP
  const sh = seriesStringCurve({ groups: groups(200), nPoints: 800 });
  assert.ok(sh.pmp > 0.75 * uni.pmp && sh.pmp < (5 / 6) * uni.pmp, `shaded pmp ratio ${sh.pmp / uni.pmp}`);
  close(sh.pmp, (5 / 6) * uni.pmp - 0.5 * sh.imp, 0.02);
  assert.ok(sh.mismatchLoss > 0.01);
  // without bypass capability (very large diode drop) the string is throttled to the shaded current
  const noBypass = seriesStringCurve({ groups: groups(200), diodeDrop: 1e4, nPoints: 800 });
  assert.ok(noBypass.pmp < sh.pmp);
  // the curve is monotonic: V falls as I rises
  for (let k = 1; k < sh.v.length; k++) assert.ok(sh.v[k] <= sh.v[k - 1] + 1e-9);
  const par = parallelStrings([uni, uni], 800);
  close(par.pmp, 2 * uni.pmp, 2e-3);
  absClose(par.mismatchLoss, 0, 2e-3);
  const mixed = parallelStrings([uni, sh], 800);
  assert.ok(mixed.pmp < uni.pmp + sh.pmp + 1e-9 && mixed.pmp > 0.8 * (uni.pmp + sh.pmp));
});

function clearDay(lat, lon, y, m, d) {
  const times = [];
  const ghi = [];
  for (let h = 0; h < 24; h++) {
    const t = Date.UTC(y, m, d, h, 30);
    times.push(t);
    ghi.push(clearSkyHaurwitz(solarPosition({ time: t, lat, lon }).zenithDeg));
  }
  return { times, ghi };
}

test('energy chain on a clear day: bookkeeping, performance ratio and loss ordering', () => {
  const { times, ghi } = clearDay(35, 0, 2023, 5, 15);
  const sys = { times, ghi, tempAirC: 25, windSpeed: 2, lat: 35, lon: 0, tiltDeg: 15, pdc0: 1e6 };
  const r = simulatePv(sys);
  close(r.energyWh, r.pac.reduce((s, x) => s + x, 0), 1e-12);
  close(r.specificYield, r.energyWh / 1e6, 1e-12);
  close(r.monthly[0].energyWh, r.energyWh, 1e-12);
  close(r.daily[0].energyWh, r.energyWh, 1e-12);
  // a clear June day at 35°N: 7.5–9 kWh/m² on the horizontal, specific yield 5.5–7.5 kWh/kWp, PR 0.75–0.9
  assert.ok(r.ghiWhPerM2 > 7500 && r.ghiWhPerM2 < 9000, `GHI ${r.ghiWhPerM2}`);
  assert.ok(r.specificYield > 5.5 && r.specificYield < 7.5, `yield ${r.specificYield}`);
  assert.ok(r.performanceRatio > 0.75 && r.performanceRatio < 0.9, `PR ${r.performanceRatio}`);
  close(r.performanceRatio, r.energyWh / 1e6 / (r.poaWhPerM2 / 1000), 1e-12);
  assert.ok(r.pac.every((p) => p <= 1e6 / 1.2 + 1e-9));
  assert.ok(Math.max(...r.tempCellC) > 40);
  // each loss mechanism lowers the yield; cooling and bifacial gain raise it
  const ideal = simulatePv({ ...sys, soiling: 0, dcLoss: 0, acLoss: 0, availability: 1, iam: 'none', dcAcRatio: 0.8, inverterEfficiency: 1 });
  assert.ok(ideal.energyWh > r.energyWh);
  assert.ok(simulatePv({ ...sys, thermalParams: { u0: 45, u1: 7 } }).energyWh > r.energyWh);
  assert.ok(simulatePv({ ...sys, bifaciality: 0.7 }).energyWh > r.energyWh);
  close(simulatePv({ ...sys, year: 10 }).dcEnergyWh, r.dcEnergyWh * 0.995 ** 10, 1e-9);
  close(simulatePv({ ...sys, availability: 0.5 }).energyWh, (r.energyWh * 0.5) / 0.99, 1e-9);
  assert.ok(simulatePv({ ...sys, dcAcRatio: 1.6 }).clippedWh > 0);
  // wave-induced tilt oscillation costs a little energy on a clear day
  const wavy = simulatePv({ ...sys, tiltAmplitudeDeg: 10 });
  assert.ok(wavy.energyWh < r.energyWh && wavy.energyWh > 0.95 * r.energyWh);
});

test('wave tilt effect: no oscillation → 1; isotropic diffuse-only sky → mean of (1 + cos β)/2', () => {
  const o = { dni: 0, dhi: 200, ghi: 200, zenithDeg: 40, azimuthDeg: 180, tiltDeg: 10, surfaceAzimuthDeg: 180, albedo: 0, model: 'isotropic' };
  assert.equal(waveTiltEffect({ ...o, amplitudeDeg: 0 }).factor, 1);
  const r = waveTiltEffect({ ...o, amplitudeDeg: 8 });
  // mean over δ = A sin φ of (1 + cos(β + δ))/2: cos β J0(A) for the cosine term
  const A = (8 * Math.PI) / 180;
  const j0 = 1 - (A * A) / 4 + A ** 4 / 64;
  close(r.mean, 200 * (1 + Math.cos((10 * Math.PI) / 180) * j0) / 2, 1e-6);
  assert.ok(r.factor < 1);
});

test('synthetic year from monthly climatology reproduces the monthly insolation; optimum tilt is plausible', () => {
  const monthlyGhi = [2200, 3100, 4400, 5600, 6600, 7200, 7300, 6500, 5000, 3600, 2400, 1900];
  const y = syntheticYear({ monthlyGhi, monthlyTempC: [8, 9, 12, 15, 19, 24, 27, 27, 23, 18, 12, 9], diurnalRangeC: 8, lat: 38, lon: -5 });
  assert.equal(y.times.length, 8760);
  const sum = new Float64Array(12);
  const days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
  for (let i = 0; i < 8760; i++) sum[new Date(y.times[i]).getUTCMonth()] += y.ghi[i];
  for (let m = 0; m < 12; m++) close(sum[m] / days[m], monthlyGhi[m], 1e-9, `month ${m}`);
  assert.ok(y.clearness.every((k) => k > 0.3 && k < 0.8));
  const system = { times: y.times, ghi: y.ghi, tempAirC: y.tempAirC, lat: 38, lon: -5, pdc0: 1e6 };
  const flat = simulatePv({ ...system, tiltDeg: 0 });
  assert.ok(flat.specificYield > 1300 && flat.specificYield < 1900, `yield ${flat.specificYield}`);
  assert.equal(flat.monthly.length, 12);
  const opt = optimalTilt({ system, tiltMin: 0, tiltMax: 60, tiltStep: 10 });
  assert.ok(opt.tiltDeg > 20 && opt.tiltDeg < 40, `tilt ${opt.tiltDeg}`);
  assert.ok(opt.energyWh > flat.energyWh);
  assert.equal(opt.surfaceAzimuthDeg, 180);
});
