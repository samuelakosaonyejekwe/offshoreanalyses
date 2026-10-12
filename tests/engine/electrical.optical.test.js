import test from 'node:test';
import assert from 'node:assert/strict';
import * as O from '../../src/engine/electrical/optical.js';
import { besselJ0, besselJ1, besselK0, besselK1 } from '../../src/engine/math/special.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('link budget by hand: 100 km at 0.20 dB/km, 10 splices, 2 connectors → 21.5 dB; margin and maximum span', () => {
  assert.equal(O.FIBRE_ATTENUATION_DB_PER_KM['G.652'][1550], 0.20); assert.equal(O.FIBRE_ATTENUATION_DB_PER_KM['G.654'][1310], null);
  const b = O.linkBudget({ length: 100e3, attenuationDbPerKm: 0.2, launchPowerDbm: 0, sensitivityDbm: -28, nSplices: 10 });
  close(b.fibreLossDb, 20, 1e-12); close(b.spliceLossDb, 0.5, 1e-12); close(b.connectorLossDb, 1.0, 1e-12); close(b.totalLossDb, 21.5, 1e-12);
  close(b.receivedPowerDbm, -21.5, 1e-12); close(b.marginDb, 3.5, 1e-12); assert.equal(b.ok, true); assert.equal(b.nSplices, 10);
  // zero margin with the same lumped losses: 0.2 L + 1.5 = 25 → 117.5 km
  rel(b.maxLength, 117.5e3, 1e-6);
  // splice interval 5 km: ceil(100/5) − 1 + 2 = 21 splices; repair, ageing and a Raman gain enter one to one
  const s = O.linkBudget({ length: 100e3, attenuationDbPerKm: 0.2, launchPowerDbm: 0, sensitivityDbm: -28, spliceInterval: 5e3, repairMarginDb: 2, ageingDb: 0.5, amplifierGainDb: 4 });
  assert.equal(s.nSplices, 21); close(s.totalLossDb, 20 + 1.05 + 1 + 2.5, 1e-12); close(s.receivedPowerDbm, 4 - 24.55, 1e-12);
  assert.equal(O.linkBudget({ length: 150e3, attenuationDbPerKm: 0.2, launchPowerDbm: 0, sensitivityDbm: -28, nSplices: 10 }).ok, false);
});

test('spectral attenuation at 1550 nm: Rayleigh 0.9/1.55⁴ = 0.156, infrared tail 0.020, total 0.186 dB/km; OH peak at 1383 nm', () => {
  const a = O.spectralAttenuation({ wavelengthNm: 1550 });
  rel(a.rayleigh, 0.9 / 1.55 ** 4, 1e-12); rel(a.infrared, 7.81e11 * Math.exp(-48.48 / 1.55), 1e-12);
  close(a.total, 0.1559 + 0.01 + 0.0204, 2e-4); close(a.oh, 0, 1e-15);
  // Rayleigh ∝ λ⁻⁴: 1310 nm is (1550/1310)⁴ = 1.96 times higher
  rel(O.spectralAttenuation({ wavelengthNm: 1310 }).rayleigh / a.rayleigh, (1550 / 1310) ** 4, 1e-12);
  close(O.spectralAttenuation({ wavelengthNm: 1383, ohPeak: 0.5 }).oh, 0.5, 1e-12);
  rel(O.spectralAttenuation({ wavelengthNm: 1398, ohPeak: 0.5 }).oh, 0.5 * Math.exp(-0.5), 1e-9); // one σ away
});

test('chromatic dispersion of G.652 fibre: zero at 1310 nm, 17.5 ps/(nm·km) at 1550 nm; β2 = −Dλ²/(2πc)', () => {
  close(O.chromaticDispersion({ wavelengthNm: 1310 }).DPsPerNmKm, 0, 1e-12);
  const d = O.chromaticDispersion({ wavelengthNm: 1550 });
  rel(d.DPsPerNmKm, 0.092 / 4 * (1550 - 1310 ** 4 / 1550 ** 3), 1e-12); close(d.DPsPerNmKm, 17.46, 0.01);
  // D = 17.46e-6 s/m², λ² = 2.4025e-12 m²: β2 = −2.227e-26 s²/m (−22.3 ps²/km)
  rel(d.beta2, -17.46e-6 * 2.4025e-12 / (2 * Math.PI * 299792458), 1e-3);
  assert.ok(O.chromaticDispersion({ wavelengthNm: 1280 }).DPsPerNmKm < 0);
});

test('dispersion- and PMD-limited lengths, and dispersion compensation', () => {
  // broad source: L = 1/(4 B |D| σλ) = 1/(4 × 2.5e9 × 17e-12 s/km) = 5.88 km
  const b = O.dispersionLimitedLength({ bitRate: 2.5e9, DPsPerNmKm: 17, wavelengthNm: 1550, sourceWidthNm: 1 });
  rel(b.length, 1 / (4 * 2.5e9 * 17e-12) * 1000, 1e-9);
  // narrow source: L = 1/(16 B² |β2|); 10 Gbit/s at 17 ps/(nm·km) (β2 = −21.68 ps²/km) → 28.8 km; ×16 shorter at 40 Gbit/s
  const n = O.dispersionLimitedLength({ bitRate: 10e9, DPsPerNmKm: 17, wavelengthNm: 1550 });
  const beta2 = 17e-6 * (1550e-9) ** 2 / (2 * Math.PI * 299792458);
  rel(Math.abs(n.beta2), beta2, 1e-9); rel(n.length, 1 / (16 * 1e20 * beta2), 1e-9); close(n.length / 1000, 28.8, 0.1);
  rel(O.dispersionLimitedLength({ bitRate: 40e9, DPsPerNmKm: 17, wavelengthNm: 1550 }).length, n.length / 16, 1e-12);
  assert.equal(O.dispersionLimitedLength({ bitRate: 10e9, DPsPerNmKm: 0, wavelengthNm: 1310 }).length, Infinity);
  assert.notEqual(b.criterion, n.criterion);
  // PMD: 10 Gbit/s → 100 ps bit period, 10 ps allowed; 0.1 ps/√km → √L = 100 → 10 000 km; DGD of 400 km = 2 ps
  const p = O.pmdLimitedLength({ bitRate: 10e9, pmdPsPerSqrtKm: 0.1, length: 400e3 });
  rel(p.length, 1e7, 1e-9); rel(p.dgdPs, 2, 1e-12);
  // 80 km × 17 = 1360 ps/nm; DCF at −100 ps/(nm·km) → 13.6 km; 90 % compensation leaves 136 ps/nm
  const c = O.dispersionCompensation({ length: 80e3, DPsPerNmKm: 17, DdcfPsPerNmKm: -100 });
  rel(c.accumulatedPsPerNm, 1360, 1e-12); rel(c.lengthDcf, 13.6e3, 1e-12); close(c.residualPsPerNm, 0, 1e-9);
  rel(O.dispersionCompensation({ length: 80e3, DPsPerNmKm: 17, DdcfPsPerNmKm: -100, ratio: 0.9 }).residualPsPerNm, 136, 1e-9);
});

test('step-index fibre: NA, V number, cut-off wavelength and Marcuse mode-field diameter by hand', () => {
  // a = 4.1 µm, n1 = 1.4500, n2 = 1.4447: NA = √(2.1025 − 2.08716) = 0.1239; V(1550) = 2.059; λc = 1327 nm
  const f = O.stepIndexFibre({ coreRadius: 4.1e-6, n1: 1.45, n2: 1.4447, wavelengthNm: 1550 });
  const NA = Math.sqrt(1.45 ** 2 - 1.4447 ** 2);
  rel(f.NA, NA, 1e-12); close(f.NA, 0.1239, 1e-4); rel(f.V, 2 * Math.PI * 4.1e-6 * NA / 1.55e-6, 1e-12); close(f.V, 2.059, 2e-3);
  rel(f.cutoffWavelengthNm, 2 * Math.PI * 4.1e-6 * NA / 2.405 * 1e9, 1e-9); close(f.cutoffWavelengthNm, 1327, 1.5);
  assert.equal(f.singleMode, true); assert.equal(f.modes, 2); rel(f.delta, NA * NA / (2 * 1.45 ** 2), 1e-12);
  // w/a = 0.65 + 1.619 V^−1.5 + 2.879 V^−6 = 0.65 + 0.5481 + 0.0378 = 1.236 → MFD = 10.13 µm
  rel(f.modeFieldDiameter, 2 * 4.1e-6 * (0.65 + 1.619 * f.V ** -1.5 + 2.879 * f.V ** -6), 1e-12); close(f.modeFieldDiameter * 1e6, 10.13, 0.03);
  rel(f.acceptanceAngle, Math.asin(NA), 1e-12);
  // multimode 50 µm fibre at 850 nm: V = 2π × 25e-6 × 0.2/0.85e-6 = 36.96 → about 683 modes
  const mm = O.stepIndexFibre({ coreRadius: 25e-6, n1: 1.48, n2: Math.sqrt(1.48 ** 2 - 0.04), wavelengthNm: 850 });
  assert.equal(mm.singleMode, false); close(mm.modes, 683, 2); assert.ok(Number.isNaN(mm.modeFieldDiameter));
  assert.throws(() => O.stepIndexFibre({ coreRadius: 4e-6, n1: 1.44, n2: 1.45, wavelengthNm: 1550 }), RangeError);
});

test('refractive-index profile: the documented default g = Infinity is the step profile; g = 2 the parabola', () => {
  const o = { coreRadius: 25e-6, n1: 1.48, delta: 0.01 };
  // default exponent (step): n1 everywhere in the core, n1 √(1 − 2Δ) in the cladding
  close(O.refractiveIndexProfile({ ...o, r: 0 }), 1.48, 1e-15); close(O.refractiveIndexProfile({ ...o, r: 20e-6 }), 1.48, 1e-15);
  close(O.refractiveIndexProfile({ ...o, r: 25e-6 }), 1.48 * Math.sqrt(0.98), 1e-15); close(O.refractiveIndexProfile({ ...o, r: 40e-6 }), 1.48 * Math.sqrt(0.98), 1e-15);
  close(O.refractiveIndexProfile({ ...o, r: 20e-6, g: Infinity }), 1.48, 1e-15);
  // parabolic: at half radius n = n1 √(1 − 2Δ/4) = 1.48 × √0.995
  close(O.refractiveIndexProfile({ ...o, r: 12.5e-6, g: 2 }), 1.48 * Math.sqrt(0.995), 1e-15);
  // large finite exponents approach the step profile; continuity at the core edge for finite g
  close(O.refractiveIndexProfile({ ...o, r: 20e-6, g: 200 }), 1.48, 1e-12);
  close(O.refractiveIndexProfile({ ...o, r: 25e-6 * (1 - 1e-12), g: 2 }), 1.48 * Math.sqrt(0.98), 1e-9);
  assert.throws(() => O.refractiveIndexProfile({ ...o, r: 1e-6, g: 0 }), RangeError);
  assert.throws(() => O.refractiveIndexProfile({ ...o, r: 1e-6, g: -2 }), RangeError);
  assert.throws(() => O.refractiveIndexProfile({ ...o, r: 1e-6, g: NaN }), RangeError);
  assert.throws(() => O.refractiveIndexProfile({ ...o, delta: 0.6, r: 1e-6 }), RangeError);
});

test('LP01 mode satisfies the characteristic equation; Gloge approximation b ≈ (1.1428 − 0.996/V)²; bend loss trends', () => {
  for (const V of [1.5, 2.0, 2.4]) {
    const m = O.lp01Mode(V);
    close(m.u * m.u + m.w * m.w, V * V, 1e-9); close(m.b, m.w * m.w / (V * V), 1e-12);
    close(m.u * besselJ1(m.u) / besselJ0(m.u), m.w * besselK1(m.w) / besselK0(m.w), 1e-7, `V = ${V}`);
    close(m.b, (1.1428 - 0.996 / V) ** 2, 0.01);
  }
  const o = { coreRadius: 4.1e-6, n1: 1.45, n2: 1.4447, wavelengthNm: 1550 };
  const b15 = O.macroBendLoss({ ...o, bendRadius: 0.015, turns: 10 });
  const b30 = O.macroBendLoss({ ...o, bendRadius: 0.030 });
  assert.ok(b15.lossDbPerM > 0 && b30.lossDbPerM < 1e-3 * b15.lossDbPerM, `${b15.lossDbPerM} vs ${b30.lossDbPerM} dB/m`);
  rel(b15.lossDbPerTurn, b15.lossDbPerM * 2 * Math.PI * 0.015, 1e-12); rel(b15.lossDb, 10 * b15.lossDbPerTurn, 1e-12);
  rel(b15.V, 2.0586, 1e-3); assert.equal(b30.lossDb, null);
  // the weaker guidance at 1625 nm bends out more easily
  assert.ok(O.macroBendLoss({ ...o, wavelengthNm: 1625, bendRadius: 0.015 }).lossDbPerM > b15.lossDbPerM);
});

test('OSNR of an amplified chain, repeater spacing, Q factor and BER', () => {
  // 58 + 0 − 5 − 20 − 10 log10(10) = 23 dB
  close(O.osnr({ channelPowerDbm: 0, noiseFigureDb: 5, spanLossDb: 20, nSpans: 10 }), 23, 1e-12);
  close(O.osnr({ channelPowerDbm: 0, noiseFigureDb: 5, spanLossDb: 20, nSpans: 20 }), 23 - 10 * Math.log10(2), 1e-12);
  // 2000 km at 0.2 dB/km, 0 dBm, NF 5 dB, required 18 dB: N spans of 400/N dB: 53 − 400/N − 10 log N ≥ 18
  const r = O.repeaterSpacing({ length: 2000e3, attenuationDbPerKm: 0.2, channelPowerDbm: 0, noiseFigureDb: 5, requiredOsnrDb: 18 });
  const ok = (N) => 53 - 400 / N - 10 * Math.log10(N) >= 18;
  assert.ok(ok(r.nSpans) && !ok(r.nSpans - 1), `N = ${r.nSpans}`);
  assert.equal(r.nRepeaters, r.nSpans - 1); rel(r.spacing, 2000e3 / r.nSpans, 1e-12); rel(r.spanLossDb, 400 / r.nSpans, 1e-12);
  close(r.osnrDb, 53 - 400 / r.nSpans - 10 * Math.log10(r.nSpans), 1e-9);
  // Q = 6 → BER = ½ erfc(6/√2) = 9.87e-10; Q = 7 → 1.28e-12; 20 log10 6 = 15.56 dB
  const q = O.qFactorBer({ Q: 6 });
  rel(q.ber, 9.866e-10, 1e-3); close(q.QDb, 15.563, 1e-3); rel(O.qFactorBer({ Q: 7 }).ber, 1.28e-12, 5e-3);
  rel(O.qFactorBer({ ber: q.ber }).Q, 6, 1e-6); rel(O.qFactorBer({ QDb: 15.563 }).Q, 6, 1e-4);
  // hard-decision FEC limit 3.8e-3 corresponds to Q = 2.67 (8.53 dB): margin of a Q = 6 signal 7.03 dB
  const f = O.qFactorBer({ Q: 6, fecLimitBer: 3.8e-3 });
  close(f.fecLimitQDb, 8.53, 0.01); close(f.marginDb, q.QDb - f.fecLimitQDb, 1e-9);
  assert.equal(q.marginDb, null);
});

test('non-linear thresholds, power feeding, sensing range and static fatigue by hand', () => {
  // α = 0.2 dB/km = 4.605e-5 1/m; 100 km: L_eff = (1 − e^(−4.605))/α = 21.50 km
  const n = O.nonlinearThresholds({ length: 100e3, attenuationDbPerKm: 0.2, effectiveArea: 80e-12 });
  const al = 0.2 / 4.342944819 / 1000;
  rel(n.effectiveLength, (1 - Math.exp(-al * 100e3)) / al, 1e-9); close(n.effectiveLength / 1000, 21.50, 0.01);
  // P_SBS = 21 × 80e-12/(5e-11 × 21 498) = 1.563 mW (1.94 dBm); P_SRS = 16 × 80e-12/(1e-13 × 21 498) = 0.595 W
  rel(n.sbsThresholdW, 21 * 80e-12 / (5e-11 * n.effectiveLength), 1e-12); close(n.sbsThresholdDbm, 1.94, 0.01);
  rel(n.srsThresholdW, 0.5954, 1e-3); close(n.srsThresholdDbm, 10 * Math.log10(n.srsThresholdW * 1000), 1e-9);
  // short fibre: L_eff → L
  rel(O.nonlinearThresholds({ length: 100, attenuationDbPerKm: 0.2, effectiveArea: 80e-12 }).effectiveLength, 100, 3e-3);
  // 100 repeaters × 40 V + 1 A × 1 Ω/km × 6000 km = 10 kV; ±5 kV with double-end feeding; 10 kW
  const p = O.powerFeedVoltage({ nRepeaters: 100, repeaterVoltage: 40, current: 1, conductorResistance: 1e-3, length: 6e6 });
  close(p.repeaterDrop, 4000, 1e-9); close(p.cableDrop, 6000, 1e-9); close(p.totalVoltage, 10000, 1e-9); close(p.voltagePerStation, 5000, 1e-9); close(p.power, 10000, 1e-9);
  close(O.powerFeedVoltage({ nRepeaters: 100, repeaterVoltage: 40, current: 1, conductorResistance: 1e-3, length: 6e6, feed: 'single', earthPotentialMargin: 600, otherDrops: 100 }).voltagePerStation, 10700, 1e-9);
  // 20 dB of dynamic range over a two-pass 0.2 dB/km fibre: 50 km; 2 dB of one-way splice loss costs 10 km
  close(O.sensingRange({ dynamicRangeDb: 20, attenuationDbPerKm: 0.2 }), 50e3, 1e-6);
  close(O.sensingRange({ dynamicRangeDb: 20, attenuationDbPerKm: 0.2, fixedLossDb: 2 }), 40e3, 1e-6);
  close(O.sensingRange({ dynamicRangeDb: 20, attenuationDbPerKm: 0.2, passes: 1 }), 100e3, 1e-6);
  // proof 0.69 GPa for 1 s, n = 20: at one third of the proof stress t = 3^20 s = 3.49e9 s (110 years)
  const s = O.fibreStaticFatigue({ proofStress: 0.69e9, serviceStress: 0.23e9, designLife: 25 * 365.25 * 86400 });
  rel(s.lifetime, 3 ** 20, 1e-9);
  // allowable stress for 25 years: σp (1/7.889e8)^(1/20) = 0.3594 σp
  rel(s.allowableStressRatio, (1 / (25 * 365.25 * 86400)) ** (1 / 20), 1e-12); close(s.allowableStressRatio, 0.3594, 1e-3);
  rel(s.allowableStress, 0.69e9 * s.allowableStressRatio, 1e-12);
  assert.equal(O.fibreStaticFatigue({ proofStress: 0.69e9 }).lifetime, null);
});
