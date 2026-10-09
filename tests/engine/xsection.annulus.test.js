import test from 'node:test';
import assert from 'node:assert/strict';
import {
  R_GAS, H2S_SOUR_LIMIT, permeabilityFromCm3, arrheniusPermeability, annulusPermeation, annulusFreeVolume,
  annulusPressureBuildUp, waterSaturationPressure, annulusClassification, annulusPH, annulusFlooding,
  annulusVolumeTest, rgdScreening,
} from '../../src/engine/xsection/annulus.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);

test('permeation of one component by hand: n = Pe dp 2 pi L / ln(ro/ri)', () => {
  // Pe = 1e-15 mol m/(m^2 s Pa), CO2 at 5 bar in the bore, ri = 80 mm, ro = 88 mm, L = 1000 m
  const r = annulusPermeation({ components: [{ name: 'CO2', Pe: 1e-15, fugacityBore: 5e5 }], ri: 0.08, ro: 0.088, length: 1000 });
  // 2 pi 1000 / ln(1.1) = 65923.9;  x 1e-15 x 5e5 = 3.29619e-5 mol/s
  close(r.totalMolarRate, 1e-15 * 5e5 * 2 * Math.PI * 1000 / Math.log(1.1), 1e-12);
  close(r.totalMolarRate, 3.2962e-5, 1e-4);
  // thin-wall limit: flux ~ Pe dp/t x area at the mean radius
  const thin = annulusPermeation({ components: [{ name: 'x', Pe: 1e-15, fugacityBore: 5e5 }], ri: 0.1, ro: 0.1001, length: 1 });
  close(thin.totalMolarRate, 1e-15 * 5e5 / 1e-4 * 2 * Math.PI * 0.10005, 1e-6);
  // back-pressure and shielding reduce the rate linearly
  const b = annulusPermeation({ components: [{ name: 'CO2', Pe: 1e-15, fugacityBore: 5e5, fugacityAnnulus: 1e5 }], ri: 0.08, ro: 0.088, length: 1000, areaFactor: 0.5 });
  close(b.totalMolarRate, 0.4 * r.totalMolarRate, 1e-12);
});

test('Arrhenius temperature dependence and unit conversion', () => {
  const a = arrheniusPermeability({ Pe0: 1e-9, Ea: 40e3, T: 333.15 });
  close(a.Pe, 1e-9 * Math.exp(-40e3 / (8.314462618 * 333.15)), 1e-12);
  // ratio between two temperatures depends only on Ea
  const b = arrheniusPermeability({ Pe0: 1e-9, Ea: 40e3, T: 353.15 });
  close(Math.log(b.Pe / a.Pe), 40e3 / R_GAS * (1 / 333.15 - 1 / 353.15), 1e-10);
  // 1e-8 cm3(STP) cm/(cm2 s bar) = 1e-8 * 1e-3/22414 mol m/(m2 s Pa)
  close(permeabilityFromCm3(1e-8).Pe, 4.4615e-16, 1e-4);
  const r = annulusPermeation({ components: [{ name: 'CH4', Pe0: 1e-9, Ea: 40e3, fugacityBore: 1e7 }, { name: 'H2S', Pe: 2 * a.Pe, fugacityBore: 1e5 }], ri: 0.08, ro: 0.088, length: 10, T: 333.15 });
  close(r.moleFractions.H2S, 2e5 / (1e7 + 2e5), 1e-12);
  close(r.moleFractions.CH4 + r.moleFractions.H2S, 1, 1e-12);
});

test('free volume and pressure build-up to vent opening by the ideal-gas law', () => {
  const v = annulusFreeVolume({ layers: [{ rIn: 0.09, rOut: 0.098, fill: 0.9 }, { rIn: 0.098, rOut: 0.106, fill: 0.88 }], length: 1000 });
  const a1 = Math.PI * (0.098 ** 2 - 0.09 ** 2), a2 = Math.PI * (0.106 ** 2 - 0.098 ** 2);
  close(v.volume, (0.1 * a1 + 0.12 * a2) * 1000, 1e-12);
  // n_dot = 3e-5 mol/s into 0.5 m3 at 300 K: dp/dt = 3e-5 * 8.3145 * 300/0.5 = 0.14966 Pa/s
  const b = annulusPressureBuildUp({ molarRate: 3e-5, volume: 0.5, T: 300, p0: 1e5, pVent: 3e5 });
  close(b.pressureRate, 0.149660, 1e-5);
  close(b.timeToVent, 2e5 / b.pressureRate, 1e-12); // about 15.5 days
  close(b.molesPerVent, 3e-5 * b.timeToVent, 1e-12);
  close(b.ventRateStd, 3e-5 * R_GAS * 288.15 / 101325, 1e-12);
  assert.throws(() => annulusPressureBuildUp({ molarRate: 1, volume: 1, T: 300, p0: 2e5, pVent: 1e5 }), /pVent/);
});

test('sour/sweet classification at the ISO 15156 limit of 0.3 kPa H2S', () => {
  assert.equal(H2S_SOUR_LIMIT, 300);
  // 0.3 kPa = 0.0435 psi, i.e. the 0.05 psi threshold of NACE MR0175 rounded
  close(300 / 6894.757, 0.0435, 1e-3);
  const at = (x) => annulusClassification({ moleFractions: { CH4: 1 - 0.05 - x, CO2: 0.05, H2S: x }, pAnnulus: 2e5 });
  assert.equal(at(0.0015).classification, 'sour'); // exactly 300 Pa
  assert.equal(at(0.00149).classification, 'sweet');
  close(at(0.003).h2sRatio, 2, 1e-12);
  assert.equal(annulusClassification({ moleFractions: { CH4: 1 }, pAnnulus: 2e5 }).classification, 'benign');
  // water: 5 % of 2 bar = 10 kPa; saturation at 40 degC is 7.38 kPa -> condenses; at 60 degC 19.9 kPa -> stays vapour
  close(waterSaturationPressure(100).pSat, 101325, 2e-3);
  close(waterSaturationPressure(40).pSat, 7.38e3, 5e-3);
  assert.equal(annulusClassification({ moleFractions: { H2O: 0.05 }, pAnnulus: 2e5, TC: 40 }).waterCondenses, true);
  assert.equal(annulusClassification({ moleFractions: { H2O: 0.05 }, pAnnulus: 2e5, TC: 60 }).waterCondenses, false);
});

test('pH of CO2-saturated pure water: 3.91 at 1 atm CO2 (25 degC), shifted by the confinement argument', () => {
  close(annulusPH({ pCO2: 101325 }).pH, 3.91, 1e-12);
  close(annulusPH({ pCO2: 10132.5 }).pH, 4.41, 1e-12);
  close(annulusPH({ pCO2: 101325, confinementShift: 2 }).pH, 5.91, 1e-12);
  assert.equal(annulusPH({ pCO2: 101325, pHOverride: 6.2 }).pH, 6.2);
});

test('flooding scenario and annulus volume tests', () => {
  const f = annulusFlooding({ freeVolume: 0.4, floodedLength: 500, breachDepth: 100 });
  close(f.addedMassPerLength, 1025 * 0.4 / 500); close(f.annulusPressure, 101325 + 1025 * 9.80665 * 100);
  assert.ok(f.consequences.length >= 3);
  // pressure test: 100 std litres raise the annulus by 0.25 bar at 288.15 K -> V = 1.01325e5 * 0.1/0.25e5 = 0.4053 m3
  const p = annulusVolumeTest({ method: 'pressure', p1: 1e5, p2: 1.25e5, injectedStdVolume: 0.1, designVolume: 0.8, pressureDrift: -0.01 });
  close(p.freeVolume, 101325 * 0.1 / 0.25e5, 1e-12);
  close(p.floodedFraction, 1 - p.freeVolume / 0.8, 1e-12);
  close(p.leakMolarRate, -0.01 * p.freeVolume / (R_GAS * 288.15), 1e-12);
  // vacuum test: Boyle's law mass balance holds for the returned volume
  const v = annulusVolumeTest({ method: 'vacuum', p1: 5e3, p2: 4e4, refVolume: 0.05, refPressure: 1e5 });
  close(5e3 * v.freeVolume + 1e5 * 0.05, 4e4 * (v.freeVolume + 0.05), 1e-12);
  assert.throws(() => annulusVolumeTest({ method: 'x', p1: 1, p2: 2 }), /method/);
});

test('rapid gas decompression screening against qualified rate', () => {
  const r = rgdScreening({ pressureDrop: 200e5, duration: 600, qualifiedRate: 70e5 / 60, qualifiedPressure: 300e5, initialPressure: 250e5 });
  close(r.rate, 200e5 / 600); close(r.rateUtilisation, (200e5 / 600) / (70e5 / 60));
  assert.equal(r.risk, 'within qualification');
  assert.equal(rgdScreening({ pressureDrop: 200e5, duration: 60, qualifiedRate: 70e5 / 60 }).risk, 'outside qualification');
});
