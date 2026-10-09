import test from 'node:test';
import assert from 'node:assert/strict';
import * as q from '../../src/engine/geotech/liquefaction.js';

const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);

test('rd (Liao & Whitman): 0.9618 at 5 m, 0.7735 at 15 m; continuous near 9.15 m; Blake fit of the mean curve within 5 %', () => {
  close(q.stressReduction(5), 1 - 0.00765 * 5, 1e-14);
  close(q.stressReduction(15), 1.174 - 0.0267 * 15, 1e-14);
  assert.ok(Math.abs(q.stressReduction(9.15) - q.stressReduction(9.1500001)) < 1e-3);
  for (const z of [2, 5, 9, 12, 18]) assert.ok(Math.abs(q.stressReduction(z, 'blake') / q.stressReduction(z) - 1) < 0.05, `z = ${z}`);
});

test('CSR by hand: amax = 0.25 g, sigma_v = 150 kPa, sigma_v\' = 90 kPa, z = 8 m', () => {
  const r = q.cyclicStressRatio({ amax: 0.25 * 9.80665, sigmaV: 150e3, sigmaVeff: 90e3, z: 8 });
  close(r, 0.65 * 0.25 * (150 / 90) * (1 - 0.00765 * 8), 1e-12);
});

test('SPT clean-sand base curve (Rauch, Youd et al. 2001 eq. 4): CRR(10) = 0.113, CRR(20) = 0.215; N >= 30 non-liquefiable', () => {
  // 1/24 + 10/135 + 50/145^2 - 1/200 = 0.04167 + 0.07407 + 0.00238 - 0.005
  close(q.crrSpt(10), 0.11312, 1e-4);
  // 1/14 + 20/135 + 50/245^2 - 1/200
  close(q.crrSpt(20), 0.21541, 1e-4);
  close(q.crrSpt(0), 1 / 34 + 50 / 2025 - 0.005, 1e-12);
  assert.equal(q.crrSpt(30), Infinity);
});

test('fines correction: FC = 15 % gives alpha = 2.498, beta = 1.0481; limits at 5 % and 35 %', () => {
  const r = q.finesCorrectionSpt(12, 15);
  close(r.alpha, Math.exp(1.76 - 190 / 225), 1e-12);
  close(r.alpha, 2.4982, 1e-4);
  close(r.beta, 0.99 + 15 ** 1.5 / 1000, 1e-12);
  close(r.N160cs, 2.4982 + 1.04809 * 12, 1e-4);
  assert.deepEqual(q.finesCorrectionSpt(12, 3), { N160cs: 12, alpha: 0, beta: 1 });
  close(q.finesCorrectionSpt(12, 50).N160cs, 5 + 1.2 * 12, 1e-12);
  // continuity at FC = 35 %: alpha -> exp(1.76 - 190/1225) = 4.98, beta -> 1.197
  close(q.finesCorrectionSpt(10, 34.999).N160cs, q.finesCorrectionSpt(10, 35).N160cs, 5e-3);
});

test('overburden correction CN = sqrt(Pa/sigma) capped at 1.7', () => {
  close(q.overburdenCN(101325), 1, 1e-12);
  close(q.overburdenCN(4 * 101325), 0.5, 1e-12);
  assert.equal(q.overburdenCN(10e3), 1.7);
});

test('CPT (Robertson & Wride 1998): CRR(100) = 0.173, CRR(30) = 0.075; Kc polynomial continuous at Ic = 1.64', () => {
  close(q.crrCpt(100), 93 * 0.001 + 0.08, 1e-12);
  close(q.crrCpt(30), 0.833 * 0.03 + 0.05, 1e-12);
  // the two branches meet near qc1Ncs = 50: 0.0917 vs 0.0916
  assert.ok(Math.abs(q.crrCpt(49.999) - q.crrCpt(50)) < 5e-4);
  const kc = (Ic) => -0.403 * Ic ** 4 + 5.581 * Ic ** 3 - 21.63 * Ic ** 2 + 33.75 * Ic - 17.88;
  assert.ok(Math.abs(kc(1.64) - 1) < 0.01);
});

test('CPT normalisation by hand: clean sand, n = 0.5 branch', () => {
  const pa = 101325;
  const o = { qc: 8e6, fs: 40e3, sigmaV: 100e3, sigmaVeff: 60e3 };
  const r = q.cptNormalise(o);
  const F = 100 * 40e3 / 7.9e6;
  const Q = (7.9e6 / pa) * (pa / 60e3) ** 0.5;
  const Ic = Math.sqrt((3.47 - Math.log10(Q)) ** 2 + (1.22 + Math.log10(F)) ** 2);
  assert.equal(r.n, 0.5);
  close(r.Ic, Ic, 1e-12);
  close(r.qc1N, Math.sqrt(pa / 60e3) * 8e6 / pa, 1e-12);
  assert.ok(Ic > 1.64 && Ic < 2.05);
  close(r.qc1Ncs, r.Kc * r.qc1N, 1e-12);
  // soft clay is flagged as clay-like (Ic > 2.6 with n = 1)
  assert.equal(q.cptNormalise({ qc: 0.6e6, fs: 20e3, sigmaV: 100e3, sigmaVeff: 60e3 }).clayLike, true);
});

test('shear-wave velocity method (Andrus & Stokoe 2000): Vs1 = 150 m/s, clean sand -> CRR = 0.0796', () => {
  const r = q.crrVs({ Vs: 150, sigmaVeff: 101325 });
  close(r.Vs1, 150, 1e-12);
  // 0.022 (1.5)^2 + 2.8 (1/65 - 1/215) = 0.0495 + 0.03005
  close(r.CRR, 0.07955, 2e-4);
  assert.equal(q.crrVs({ Vs: 150, sigmaVeff: 101325, FC: 20 }).Vs1Limit, 207.5);
  assert.equal(q.crrVs({ Vs: 220, sigmaVeff: 101325 }).CRR, Infinity);
});

test('magnitude scaling factor (Idriss): 1.00 at M 7.5, 1.77 at M 6, 0.85 at M 8 (Youd et al. 2001 Table 3)', () => {
  close(q.magnitudeScalingFactor(7.5), 1.0, 2e-3);
  close(q.magnitudeScalingFactor(6), 1.76, 6e-3);
  close(q.magnitudeScalingFactor(8), 0.84, 1e-2);
  close(q.magnitudeScalingFactor(5.5), 2.20, 6e-3);
});

test('K-sigma and factor of safety', () => {
  close(q.kSigma(4 * 101325, 0.7), 4 ** -0.3, 1e-12);
  assert.equal(q.kSigma(50e3), 1);
  close(q.liquefactionFs({ CRR75: 0.2, CSR: 0.25, Mw: 7.5 }), 0.8 * q.magnitudeScalingFactor(7.5), 1e-12);
});

test('SPT profile of a submerged seabed: stresses by hand and FS increasing with blow count', () => {
  const pts = [{ z: 2, N160: 8, gamma: 19e3 }, { z: 5, N160: 14, gamma: 19e3 }, { z: 9, N160: 26, FC: 10, gamma: 20e3 }];
  const r = q.sptProfile({ points: pts, waterTableDepth: 0, amax: 0.2 * 9.80665, Mw: 7 });
  close(r[1].sigmaV, 95e3, 1e-12);
  close(r[1].sigmaVeff, 95e3 - 5 * 9810, 1e-12);
  close(r[1].CSR, 0.65 * 0.2 * (95e3 / (95e3 - 49050)) * (1 - 0.00765 * 5), 1e-12);
  assert.ok(r[0].FS < r[1].FS && r[1].FS < r[2].FS);
});

test('liquefaction potential index (Iwasaki): uniform FS = 0.5 over 20 m gives 50; FS >= 1 gives 0', () => {
  close(q.liquefactionPotentialIndex([0, 10, 20], [0.5, 0.5, 0.5]), 50, 1e-12);
  assert.equal(q.liquefactionPotentialIndex([1, 10, 19], [1.2, 1.5, 2]), 0);
  // FS = 0 in the top 5 m only: integral of (10 - 0.5 z) from 0 to 5 = 43.75 (sharp boundary approximated)
  const z = [], fs = [];
  for (let i = 0; i <= 2000; i++) { z.push(i * 0.01); fs.push(i * 0.01 <= 5 ? 0 : 1.5); }
  close(q.liquefactionPotentialIndex(z, fs), 43.75, 2e-3);
});

test('post-liquefaction strains (Idriss & Boulanger 2008): N = 10, FS = 0.5 -> gamma_max = 47 %, eps_v = 3.7 %', () => {
  const r = q.postLiquefactionStrain(10, 0.5);
  close(r.gammaLim, 1.859 * (1.1 - Math.sqrt(10 / 46)) ** 3, 1e-12);
  close(r.gammaMax, 0.473, 2e-3);
  close(r.epsV, 1.5 * Math.exp(-0.369 * Math.sqrt(10)) * 0.08, 1e-12);
  close(r.epsV, 0.0374, 3e-3);
  // the blow-count and relative-density forms of the monograph agree to rounding: Dr = sqrt(N/46)
  const Dr = Math.sqrt(10 / 46);
  close(r.Falpha, 0.032 + 4.7 * Dr - 6.0 * Dr * Dr, 1e-2);
  assert.equal(q.postLiquefactionStrain(20, 2.1).epsV, 0);
  // dense sand: limiting strain vanishes at Dr = 1.1 (N = 55.7)
  assert.equal(q.postLiquefactionStrain(56, 0.5).gammaMax, 0);
  const s = q.postLiquefactionSettlement([{ thickness: 3, N160cs: 10, FS: 0.5 }, { thickness: 4, N160cs: 25, FS: 3 }]);
  close(s.settlement, 3 * r.epsV, 1e-12);
  close(s.LDI, 3 * r.gammaMax, 1e-12);
});

test('residual strength ratio (Idriss & Boulanger 2008): 0.091 at N = 10, 0.14 at N = 14; capped at tan(phi)', () => {
  close(q.residualStrengthRatio({ N160csSr: 10, phi: 0.6 }), Math.exp(10 / 16 + (-6 / 21.2) ** 3 - 3), 1e-12);
  close(q.residualStrengthRatio({ N160csSr: 10, phi: 0.6 }), 0.0909, 2e-3);
  close(q.residualStrengthRatio({ N160csSr: 14, phi: 0.6 }), Math.exp(0.875 - 0.00084 - 3), 1e-4);
  assert.ok(q.residualStrengthRatio({ N160csSr: 14, phi: 0.6, voidRedistribution: false }) > q.residualStrengthRatio({ N160csSr: 14, phi: 0.6 }));
  close(q.residualStrengthRatio({ N160csSr: 40, phi: 0.6 }), Math.tan(0.6), 1e-12);
});

test('wave-induced cyclic stress: CSR at the mudline p0 k / gamma\'; shear stress peaks at z = 1/k', () => {
  const k = 2 * Math.PI / 100;
  const r0 = q.waveCyclicStress({ H: 6, k, d: 20, z: 0, gammaSub: 9000 });
  close(r0.p0, 10055 * 6 / (2 * Math.cosh(k * 20)), 1e-12);
  close(r0.CSR, r0.p0 * k / 9000, 1e-12);
  const at = (z) => q.waveCyclicStress({ H: 6, k, d: 20, z, gammaSub: 9000 }).tau;
  assert.ok(at(1 / k) > at(0.9 / k) && at(1 / k) > at(1.1 / k));
  close(at(1 / k), r0.p0 / Math.E, 1e-12);
});

test('Seed & Rahman pore-pressure generation: ru = 1 at N = Nl, 0.5 at N/Nl = sin(pi/4)^(2 theta)', () => {
  const p = { CSR: 0.15, a: 0.48, b: 0.2 };
  const Nl = (0.48 / 0.15) ** 5;
  const full = q.seedRahmanPorePressure({ ...p, N: Nl });
  close(full.Nl, Nl, 1e-12);
  close(full.ru, 1, 1e-12);
  assert.equal(full.liquefied, true);
  close(q.seedRahmanPorePressure({ ...p, N: Nl * Math.sin(Math.PI / 4) ** 1.4 }).ru, 0.5, 1e-12);
});

test('momentary liquefaction criterion (Zen & Yamazaki): excess p0 (1 - e^(-kz)) against gamma\' z', () => {
  const k = 0.05;
  const r = q.momentaryLiquefaction({ p0: 30e3, z: 1, gammaSub: 9000, k });
  close(r.excess, 30e3 * (1 - Math.exp(-0.05)), 1e-12);
  assert.equal(r.liquefied, false);
  // strong attenuation in a partly saturated bed
  assert.equal(q.momentaryLiquefaction({ p0: 30e3, z: 1, gammaSub: 9000, attenuation: 0.5 }).liquefied, true);
});
