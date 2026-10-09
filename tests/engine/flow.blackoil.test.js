import test from 'node:test';
import assert from 'node:assert/strict';
import {
  apiToSg, sgToApi, gasDensitySc, gasPseudoCriticals, zFactorReduced, zFactor, gasFvf, gasViscosityLee,
  solutionGor, bubblePoint, oilFvfSaturated, oilCompressibility, deadOilViscosity, liveOilViscosity,
  undersaturatedOilViscosity, gasOilIft, gasWaterIft, oilProperties, gasProperties, fluidAt, T_SC, P_SC,
} from '../../src/engine/flow/blackoil.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const PSI = 6894.757293168;
const SCF_STB = 5.614583333333333;
const fromF = (F) => (F + 459.67) / 1.8;

test('API gravity conversions', () => {
  assert.ok(rel(apiToSg(10), 1) < 1e-12);
  assert.ok(rel(sgToApi(apiToSg(35)), 35) < 1e-12);
  // air at 60 °F and 1 atm: ρ = 28.9647e-3·101325/(8.314462618·288.7056) = 1.2226 kg/m³
  assert.ok(rel(gasDensitySc(1), 1.2226) < 1e-3);
});

test('Standing bubble point — hand calculation (Rs = 350 scf/STB, γg = 0.75, 30 °API, 200 °F)', () => {
  // pb = 18.2[(350/0.75)^0.83 · 10^(0.00091·200 − 0.0125·30) − 1.4]
  //    = 18.2[164.17 · 10^(−0.193) − 1.4] = 18.2[164.17·0.64121 − 1.4] = 18.2 · 103.87 = 1890 psia
  const pb = bubblePoint({ Rs: 350 / SCF_STB, T: fromF(200), api: 30, gammaG: 0.75 });
  assert.ok(rel(pb / PSI, 1890.4) < 5e-4, `pb = ${pb / PSI} psia`);
  // Bo = 0.9759 + 1.2e-4 [350 √(0.75/0.8762) + 1.25·200]^1.2 = 0.9759 + 1.2e-4 · 573.82^1.2 = 1.2212
  const Bo = oilFvfSaturated({ Rs: 350 / SCF_STB, T: fromF(200), api: 30, gammaG: 0.75 });
  assert.ok(rel(Bo, 1.2212) < 2e-4, `Bo = ${Bo}`);
});

test('solution GOR and bubble point are mutual inverses for Standing, Vazquez–Beggs and Glasø', () => {
  for (const method of ['standing', 'vazquezBeggs', 'glaso']) {
    for (const api of [22, 38]) {
      const args = { T: 350, api, gammaG: 0.8, method };
      const pb = bubblePoint({ Rs: 90, ...args });
      assert.ok(rel(solutionGor({ P: pb, ...args }), 90) < 1e-3, `${method} ${api}`);
    }
  }
});

test('Vazquez–Beggs Rs and Bo — hand calculation (2000 psia, 180 °F, 35 °API, γgs = 0.8)', () => {
  // API > 30: Rs = 0.0178·0.8·2000^1.187·exp(23.931·35/640)
  const rs = 0.0178 * 0.8 * 2000 ** 1.187 * Math.exp(23.931 * 35 / (180 + 460));
  const got = solutionGor({ P: 2000 * PSI, T: fromF(180), api: 35, gammaG: 0.8, method: 'vazquezBeggs' }) * SCF_STB;
  assert.ok(rel(got, rs) < 1e-10);
  const bo = 1 + 4.670e-4 * rs + (180 - 60) * (35 / 0.8) * (1.100e-5 + 1.337e-9 * rs);
  assert.ok(rel(oilFvfSaturated({ Rs: rs / SCF_STB, T: fromF(180), api: 35, gammaG: 0.8, method: 'vazquezBeggs' }), bo) < 1e-10);
  // separator correction: γgs = γg[1 + 5.912e-5·API·Tsep·log10(psep/114.7)] = 1 at a 114.7 psia separator
  const a = solutionGor({ P: 2000 * PSI, T: fromF(180), api: 35, gammaG: 0.8, method: 'vazquezBeggs', Tsep: fromF(80), Psep: 114.7 * PSI });
  assert.ok(rel(a * SCF_STB, rs) < 1e-10);
});

test('Sutton pseudo-criticals and Wichert–Aziz correction — hand calculation', () => {
  // γ = 0.7: Ppc = 756.8 − 91.7 − 1.764 = 663.34 psia, Tpc = 169.2 + 244.65 − 36.26 = 377.59 °R
  const s = gasPseudoCriticals({ gammaG: 0.7 });
  assert.ok(rel(s.Ppc / PSI, 663.336) < 1e-6);
  assert.ok(rel(s.Tpc * 1.8, 377.59) < 1e-6);
  // 10 % CO2 + 20 % H2S: ε = 120(0.3^0.9 − 0.3^1.6) + 15(0.2^0.5 − 0.2^4) °R
  const eps = 120 * (0.3 ** 0.9 - 0.3 ** 1.6) + 15 * (Math.sqrt(0.2) - 0.2 ** 4);
  const w = gasPseudoCriticals({ gammaG: 0.7, yCO2: 0.1, yH2S: 0.2 });
  assert.ok(rel(w.epsilon * 1.8, eps) < 1e-12);
  assert.ok(eps > 29 && eps < 31); // Wichert–Aziz chart value ≈ 30 °R
  assert.ok(rel(w.Tpc * 1.8, 377.59 - eps) < 1e-9);
  assert.ok(rel(w.Ppc / PSI, 663.336 * (377.59 - eps) / (377.59 + 0.2 * 0.8 * eps)) < 1e-9);
});

test('Z-factor against the Standing–Katz chart (Tpr = 1.5, Ppr = 2.0 → Z ≈ 0.82) and method agreement', () => {
  const dak = zFactorReduced({ Tpr: 1.5, Ppr: 2.0 });
  const hy = zFactorReduced({ Tpr: 1.5, Ppr: 2.0, method: 'hallYarborough' });
  assert.ok(Math.abs(dak - 0.82) < 0.01, `DAK ${dak}`);
  assert.ok(Math.abs(hy - 0.82) < 0.01, `HY ${hy}`);
  for (const [Tpr, Ppr] of [[1.3, 1], [1.5, 4], [2, 5], [2.5, 10], [1.2, 3]]) {
    const a = zFactorReduced({ Tpr, Ppr });
    const b = zFactorReduced({ Tpr, Ppr, method: 'hallYarborough' });
    assert.ok(rel(a, b) < 0.01, `Tpr ${Tpr} Ppr ${Ppr}: ${a} vs ${b}`);
  }
  assert.ok(Math.abs(zFactorReduced({ Tpr: 1.5, Ppr: 1e-4 }) - 1) < 1e-4);
});

test('DAK solution satisfies its own defining equation', () => {
  const Tr = 1.4;
  const Pr = 6;
  const Z = zFactorReduced({ Tpr: Tr, Ppr: Pr });
  const r = 0.27 * Pr / (Z * Tr);
  const A = [0.3265, -1.0700, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210];
  const rhs = 1 + (A[0] + A[1] / Tr + A[2] / Tr ** 3 + A[3] / Tr ** 4 + A[4] / Tr ** 5) * r
    + (A[5] + A[6] / Tr + A[7] / Tr ** 2) * r * r - A[8] * (A[6] / Tr + A[7] / Tr ** 2) * r ** 5
    + A[9] * (1 + A[10] * r * r) * (r * r / Tr ** 3) * Math.exp(-A[10] * r * r);
  assert.ok(Math.abs(Z - rhs) < 1e-10);
});

test('gas FVF and Lee–Gonzalez–Eakin viscosity', () => {
  assert.ok(rel(gasFvf({ P: P_SC, T: T_SC, Z: 1 }), 1) < 1e-12);
  // methane at 300 K and low density: K = (9.4 + 0.02·16.043)·540^1.5/(209 + 19·16.043 + 540) → 1.16e-5 Pa·s;
  // measured dilute-gas value 1.12e-5 Pa·s
  const mu = gasViscosityLee({ T: 300, rho: 0.65, M: 0.016043 });
  const K = (9.4 + 0.02 * 16.043) * 540 ** 1.5 / (209 + 19 * 16.043 + 540);
  assert.ok(rel(mu, 1e-7 * K) < 2e-3);
  assert.ok(rel(mu, 1.12e-5) < 0.05);
  assert.ok(gasViscosityLee({ T: 300, rho: 150, M: 0.016043 }) > mu);
});

test('Beggs–Robinson viscosities — hand calculation (30 °API, 200 °F, Rs = 350 scf/STB)', () => {
  // Z = 3.0324 − 0.02023·30 = 2.4255, Y = 10^Z = 266.4, X = Y·200^−1.163 = 0.5616, μod = 10^X − 1 = 2.64 cP
  const mud = deadOilViscosity({ T: fromF(200), api: 30 });
  assert.ok(rel(mud * 1e3, 10 ** (10 ** 2.4255 * 200 ** -1.163) - 1) < 1e-9);
  assert.ok(rel(mud * 1e3, 2.64) < 0.01);
  // A = 10.715·450^−0.515, B = 5.44·500^−0.338
  const live = liveOilViscosity({ muDead: mud, Rs: 350 / SCF_STB });
  assert.ok(rel(live * 1e3, 10.715 * 450 ** -0.515 * (mud * 1e3) ** (5.44 * 500 ** -0.338)) < 1e-9);
  assert.ok(live < mud);
  // undersaturated: m = 2.6 p^1.187 exp(−11.513 − 8.98e-5 p)
  const p = 4000;
  const m = 2.6 * p ** 1.187 * Math.exp(-11.513 - 8.98e-5 * p);
  assert.ok(rel(undersaturatedOilViscosity({ muOb: live, P: p * PSI, Pb: 2000 * PSI }), live * 2 ** m) < 1e-9);
});

test('oil compressibility (Vazquez–Beggs) — hand calculation', () => {
  // co = (−1433 + 5·500 + 17.2·200 − 1180·0.8 + 12.61·35)/(1e5·3000) 1/psi
  const co = (-1433 + 5 * 500 + 17.2 * 200 - 1180 * 0.8 + 12.61 * 35) / (1e5 * 3000);
  assert.ok(rel(oilCompressibility({ P: 3000 * PSI, T: fromF(200), Rsb: 500 / SCF_STB, api: 35, gammaG: 0.8 }) * PSI, co) < 1e-9);
});

test('interfacial tension correlations — hand calculation', () => {
  // dead oil at 68 °F, 30 °API: 39 − 0.2571·30 = 31.29 dyn/cm; factor at 1000 psia: 1 − 0.024·1000^0.45
  const s = gasOilIft({ P: 1000 * PSI, T: fromF(68), api: 30 });
  assert.ok(rel(s * 1e3, 31.287 * (1 - 0.024 * 1000 ** 0.45)) < 1e-9);
  // gas–water at 74 °F, 1000 psia: 75 − 1.108·1000^0.349
  assert.ok(rel(gasWaterIft({ P: 1000 * PSI, T: fromF(74) }) * 1e3, 75 - 1.108 * 1000 ** 0.349) < 1e-9);
  assert.ok(gasWaterIft({ P: 1000 * PSI, T: fromF(200) }) < gasWaterIft({ P: 1000 * PSI, T: fromF(74) }));
});

test('oil density follows the material balance ρo = (ρsto + Rs ρg,sc)/Bo and Bo is continuous at Pb', () => {
  const a = { T: 350, api: 32, gammaG: 0.75, Rsb: 120 };
  const o = oilProperties({ P: 8e6, ...a });
  assert.ok(o.saturated && o.Rs < 120);
  assert.ok(rel(o.rho, (apiToSg(32) * 999.016 + o.Rs * gasDensitySc(0.75)) / o.Bo) < 1e-12);
  const below = oilProperties({ P: o.Pb * (1 - 1e-9), ...a });
  const above = oilProperties({ P: o.Pb * (1 + 1e-9), ...a });
  assert.ok(rel(below.Bo, above.Bo) < 1e-6);
  assert.ok(rel(below.mu, above.mu) < 1e-6);
  // above the bubble point the oil shrinks and thickens with pressure
  const hi = oilProperties({ P: o.Pb * 1.5, ...a });
  assert.ok(hi.Bo < above.Bo && hi.mu > above.mu && hi.Rs === 120);
});

test('fluidAt conserves mass at every pressure and reports no free gas above the bubble point', () => {
  const stream = { qo: 0.02, GOR: 150, WC: 0.25, gammaG: 0.75, api: 32, T: 340, salinity: 35 };
  const qw = 0.02 * 0.25 / 0.75;
  const waterSc = fluidAt({ P: 1e7, ...stream });
  const mRef = 0.02 * apiToSg(32) * 999.016 + 0.02 * 150 * gasDensitySc(0.75) + waterSc.mWater;
  let lastGas = -1;
  for (const P of [6e7, 2e7, 1e7, 5e6, 2e6, 5e5]) {
    const f = fluidAt({ P, ...stream });
    assert.ok(rel(f.mTotal, mRef) < 1e-12, `mass at ${P}`);
    assert.ok(rel(f.qOil * f.rhoO + f.qGas * f.rhoG + f.qWater * f.rhoW, mRef) < 1e-10, `phase mass sum at ${P}`);
    assert.ok(f.qGas >= lastGas, 'free gas volume grows as pressure falls');
    lastGas = f.qGas;
    assert.ok(rel(f.lambdaL, f.qLiquid / f.qTotal) < 1e-12);
    assert.ok(rel(f.qWater, qw * f.Bw) < 1e-12);
  }
  const hi = fluidAt({ P: 6e7, ...stream });
  assert.equal(hi.qGas, 0);
  assert.equal(hi.lambdaL, 1);
  assert.ok(hi.Pb < 6e7);
});

test('fluidAt handles gas–water streams and superficial velocities', () => {
  const f = fluidAt({ P: 8e6, T: 300, qgSc: 20, qwSc: 0.002, gammaG: 0.65, D: 0.3 });
  const g = gasProperties({ P: 8e6, T: 300, gammaG: 0.65 });
  assert.ok(rel(f.qGas, 20 * g.Bg) < 1e-12);
  assert.ok(rel(f.vsg, f.qGas / (Math.PI * 0.09 / 4)) < 1e-12);
  assert.ok(rel(f.rhoL, f.rhoW) < 1e-12);
  assert.ok(rel(g.rho, 8e6 * 0.65 * 0.0289647 / (g.Z * 8.314462618 * 300)) < 1e-12);
  assert.throws(() => fluidAt({ P: 8e6, T: 300, qo: 0.01, WC: 1, gammaG: 0.7 }), RangeError);
  assert.throws(() => fluidAt({ P: -1, T: 300, qo: 0.01, gammaG: 0.7 }), RangeError);
});

test('zFactor from gas gravity equals the reduced-coordinate evaluation', () => {
  const pc = gasPseudoCriticals({ gammaG: 0.7 });
  assert.ok(rel(zFactor({ P: 1e7, T: 320, gammaG: 0.7 }), zFactorReduced({ Tpr: 320 / pc.Tpc, Ppr: 1e7 / pc.Ppc })) < 1e-12);
});
