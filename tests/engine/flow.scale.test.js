import test from 'node:test';
import assert from 'node:assert/strict';
import * as S from '../../src/engine/flow/scale.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const T25 = 298.15;

test('water analysis to molality: 1 mol/L NaCl at 1040 kg/m³ → 1/0.98156 = 1.0188 mol/kg, charge balanced; ionic strength', () => {
  const r = S.toMolal({ mgL: { Na: 22989.77, Cl: 35453 }, density: 1040 });
  rel(r.tds, 58442.77, 1e-9);
  rel(r.molal.Na, 1 / ((1040 - 58.44277) / 1000), 1e-9); rel(r.molal.Na, 1.0188, 1e-4); rel(r.molal.Cl, r.molal.Na, 1e-9);
  close(r.chargeBalance, 0, 1e-9); assert.equal(r.density, 1040);
  // an excess of cations shows up in the charge balance: Ca alone → +1
  close(S.toMolal({ mgL: { Ca: 400 }, density: 1000 }).chargeBalance, 1, 1e-12);
  // default density from the brine correlation is above that of fresh water
  assert.ok(S.toMolal({ mgL: { Na: 22989.77, Cl: 35453 } }).density > 1020);
  // I = ½ Σ m z²: 0.5 NaCl + 0.1 CaCl2 → ½ (0.5 + 0.4 + 0.7) = 0.8
  close(S.ionicStrength({ Na: 0.5, Ca: 0.1, Cl: 0.7 }), 0.8, 1e-12);
  // standard seawater per kg of solution: ionic strength about 0.7 mol/kg
  const sw = S.toMolal({ mgL: Object.fromEntries(Object.entries(S.SEAWATER_MG_PER_KG).map(([k, v]) => [k, v * 1.025])), density: 1025 });
  const I = S.ionicStrength(sw.molal);
  assert.ok(I > 0.65 && I < 0.75, `seawater ionic strength ${I}`);
  assert.ok(Math.abs(sw.chargeBalance) < 0.01, `seawater charge balance ${sw.chargeBalance}`);
  assert.equal(S.IONS.SO4.z, -2);
});

test('dielectric constant, Debye–Hückel parameters and activity coefficients at 25 °C', () => {
  // ε = 87.740 − 0.40008 × 25 + 9.398e-4 × 625 − 1.410e-6 × 15625 = 78.303
  close(S.waterDielectric(T25), 87.740 - 10.002 + 0.587375 - 0.02203125, 1e-9); close(S.waterDielectric(T25), 78.30, 0.01);
  const dh = S.debyeHuckelParameters(T25);
  close(dh.A, 0.5108, 1.5e-3); close(dh.B, 0.3286, 1e-3);
  // Davies, I = 0.1: log γ1 = −A (0.31623/1.31623 − 0.03) = −A × 0.21025 → γ1 = 0.781; γ2 = γ1⁴
  const g1 = S.activityCoefficient({ z: 1, I: 0.1 });
  rel(g1, 10 ** (-dh.A * (Math.sqrt(0.1) / (1 + Math.sqrt(0.1)) - 0.03)), 1e-12); close(g1, 0.781, 2e-3);
  rel(S.activityCoefficient({ z: 2, I: 0.1 }), g1 ** 4, 1e-12);
  close(S.activityCoefficient({ z: 2, I: 0 }), 1, 1e-15); close(S.activityCoefficient({ z: 0, I: 0.5 }), 1, 1e-15);
  // extended Debye–Hückel for Ca²⁺ (å = 6 Å) at I = 0.01: log γ = −4A × 0.1/(1 + 0.6 B)... = −0.1707 → 0.675
  const ge = S.activityCoefficient({ z: 2, I: 0.01, model: 'extendedDH', a: 6 });
  rel(ge, 10 ** (-dh.A * 4 * 0.1 / (1 + dh.B * 6 * 0.1)), 1e-12); close(ge, 0.675, 3e-3);
  // limiting law at very low ionic strength: both models agree
  rel(S.activityCoefficient({ z: 1, I: 1e-6 }), S.activityCoefficient({ z: 1, I: 1e-6, model: 'extendedDH', a: 4 }), 1e-5);
});

test('solubility products and carbonate constants at 25 °C match the tabulated values', () => {
  for (const [m, v] of [['calcite', -8.48], ['barite', -9.97], ['celestite', -6.63], ['gypsum', -4.58], ['anhydrite', -4.36]]) {
    close(S.solubilityProduct({ mineral: m, T: T25 }).logK, v, 0.01, m);
  }
  const b = S.solubilityProduct({ mineral: 'barite', T: T25 });
  rel(b.K, 10 ** b.logK, 1e-12); rel(b.molarMass, 0.23339, 1e-3);
  // retrograde calcite and anhydrite (less soluble hot), barite more soluble hot
  assert.ok(S.solubilityProduct({ mineral: 'calcite', T: 353.15 }).logK < -8.48);
  assert.ok(S.solubilityProduct({ mineral: 'anhydrite', T: 373.15 }).logK < -4.36);
  assert.ok(S.solubilityProduct({ mineral: 'barite', T: 353.15 }).logK > -9.97);
  // pressure correction: ΔVr = −50 cm³/mol at 300 bar: Δ ln K = 50e-6 × 299e5/(R T) = 0.6031 → Δ log K = 0.2619
  const p = S.solubilityProduct({ mineral: 'barite', T: T25, P: 300e5, dVr: -50e-6 });
  close(p.logK - b.logK, 50e-6 * (300e5 - 101325) / (8.314462 * T25) / Math.LN10, 1e-6);
  const c = S.carbonateConstants(T25);
  close(c.logKH, -1.468, 2e-3); close(c.logK1, -6.352, 2e-3); close(c.logK2, -10.329, 2e-3); close(c.logKw, -13.995, 5e-3);
  assert.throws(() => S.solubilityProduct({ mineral: 'halite', T: T25 }), RangeError);
});

test('pH of water under CO2: condensed water pH = ½(pK1 + pKH) − ½ log p; buffered water pH = pK1 + pKH + log(Alk/p)', () => {
  const c = S.carbonateConstants(T25);
  // 1 atm CO2, no alkalinity: pH = (6.352 + 1.468)/2 = 3.91; ten times the pressure: −0.5
  close(S.phFromCO2({ T: T25, pCO2: 101325 }), -0.5 * (c.logK1 + c.logKH), 1e-6); close(S.phFromCO2({ T: T25, pCO2: 101325 }), 3.91, 0.005);
  close(S.phFromCO2({ T: T25, pCO2: 1013250 }) - S.phFromCO2({ T: T25, pCO2: 101325 }), -0.5, 1e-9);
  // 0.01 mol/kg bicarbonate: a_H ≈ K1 KH p/Alk → pH = 7.82 − 2 = 5.82
  close(S.phFromCO2({ T: T25, pCO2: 101325, mHCO3: 0.01 }), -(c.logK1 + c.logKH) - 2, 1e-3);
  // exact quadratic
  const K = 10 ** (c.logK1 + c.logKH); const a = 0.5 * (-0.01 + Math.sqrt(1e-4 + 4 * K));
  close(S.phFromCO2({ T: T25, pCO2: 101325, mHCO3: 0.01 }), -Math.log10(a), 1e-9);
});

test('calcite saturation index from its definition; pH route and CO2 route agree', () => {
  const molal = { Ca: 0.002, HCO3: 0.004, Na: 0.01, Cl: 0.01 };
  const r = S.calciteSaturation({ molal, T: T25, pH: 8 });
  const I = S.ionicStrength(molal); close(r.I, I, 1e-12);
  const g1 = S.activityCoefficient({ z: 1, I }); const g2 = S.activityCoefficient({ z: 2, I });
  const c = S.carbonateConstants(T25); const logKsp = S.solubilityProduct({ mineral: 'calcite', T: T25 }).logK;
  // a_CO3 = K2 a_HCO3/a_H
  const SI = Math.log10(g2 * 0.002) + c.logK2 + Math.log10(g1 * 0.004) + 8 - logKsp;
  close(r.SI, SI, 1e-9); rel(r.SR, 10 ** SI, 1e-9); close(r.logKsp, logKsp, 1e-12);
  // one pH unit less → SI one lower; the CO2 route reproduces the pH route at the pH it computes
  close(S.calciteSaturation({ molal, T: T25, pH: 7 }).SI, SI - 1, 1e-9);
  const viaCO2 = S.calciteSaturation({ molal, T: T25, pCO2: 1000 });
  close(viaCO2.SI, S.calciteSaturation({ molal, T: T25, pH: viaCO2.pH }).SI, 1e-9);
  // more CO2 → lower pH → less scaling
  assert.ok(S.calciteSaturation({ molal, T: T25, pCO2: 1e5 }).SI < viaCO2.SI);
  assert.throws(() => S.calciteSaturation({ molal, T: T25 }), RangeError);
});

test('sulphate scale: saturation index, precipitation to equilibrium and solid mass', () => {
  // 1e-4 mol/kg Ba with 1e-3 mol/kg SO4 in 0.01 mol/kg NaCl
  const molal = { Ba: 1e-4, SO4: 1e-3, Na: 0.012, Cl: 0.0102 };
  const r = S.sulphateSaturation({ mineral: 'barite', molal, T: T25 });
  const I = S.ionicStrength(molal); const g = S.activityCoefficient({ z: 2, I });
  const K = S.solubilityProduct({ mineral: 'barite', T: T25 }).K;
  close(r.SI, Math.log10(g * g * 1e-7 / K), 1e-9); rel(r.SR, g * g * 1e-7 / K, 1e-9); assert.ok(r.SI > 2);
  // x = ½[(a + b) − √((a − b)² + 4K/γ²)]; after precipitation the product is back at K_sp
  const x = 0.5 * (1.1e-3 - Math.sqrt(0.9e-3 ** 2 + 4 * K / (g * g)));
  rel(r.precipitate, x, 1e-9); rel((1e-4 - x) * (1e-3 - x) * g * g, K, 1e-6);
  assert.ok(x > 0.99e-4 && x < 1e-4, 'nearly all barium precipitates');
  rel(r.mass, x * 0.23339, 1e-3);
  // undersaturated water precipitates nothing
  const u = S.sulphateSaturation({ mineral: 'gypsum', molal: { Ca: 1e-3, SO4: 1e-3 }, T: T25 });
  assert.ok(u.SI < 0); close(u.precipitate, 0, 1e-15); close(u.mass, 0, 1e-15);
});

test('mixing of waters and the seawater-injection scaling scan', () => {
  const mix = S.mixWaters({ Na: 0.5, SO4: 0.03 }, { Na: 1.0, Ba: 0.002 }, 0.25);
  close(mix.Na, 0.875, 1e-12); close(mix.SO4, 0.0075, 1e-12); close(mix.Ba, 0.0015, 1e-12);
  // sulphate-rich seawater into a barium-rich formation water
  const sea = { Na: 0.48, Cl: 0.56, Mg: 0.054, Ca: 0.0105, SO4: 0.029 };
  const fw = { Na: 1.2, Cl: 1.42, Ca: 0.05, Ba: 0.002, Sr: 0.004 };
  const r = S.seawaterMixingScale({ seawater: sea, formationWater: fw, T: 353.15, n: 21 });
  assert.equal(r.fraction.length, 21); close(r.fraction[0], 0, 1e-12); close(r.fraction[20], 1, 1e-12);
  const b = r.results.barite;
  // neither end member alone can precipitate barite (no sulphate in one, no barium in the other)
  close(b.mass[0], 0, 1e-12); close(b.mass[20], 0, 1e-12);
  assert.ok(b.worstFraction > 0 && b.worstFraction < 1 && b.worstMass > 0 && b.maxSI > 1);
  // mass balance: the largest deposit cannot exceed the barium available at that ratio (kg/m³)
  for (let i = 0; i < 21; i++) assert.ok(b.mass[i] <= (1 - r.fraction[i]) * 0.002 * 0.23339 * 1000 * (1 + 1e-9));
  // barium is the limiting ion at high seawater fraction: nearly complete precipitation there
  rel(b.mass[10], 0.5 * 0.002 * 233.39, 0.02);
  assert.equal(b.worstMass, Math.max(...b.mass));
  assert.ok(r.results.celestite && r.results.anhydrite);
});

test('Langelier / Ryznar and Stiff–Davis indices by hand; deposition-rate screen', () => {
  // pH 7.5, 25 °C, TDS 500, Ca 200 and alkalinity 150 mg/L as CaCO3:
  //   A = 0.1699, B = −13.12 log 298.15 + 34.55 = 2.0855, C = log 200 − 0.4 = 1.9010, D = log 150 = 2.1761 → pH_s = 7.478
  const l = S.langelierRyznar({ pH: 7.5, T: T25, tds: 500, calcium: 200, alkalinity: 150 });
  close(l.pHs, 9.3 + 0.16990 + 2.0855 - 1.9010 - 2.1761, 1e-3); close(l.LSI, 7.5 - l.pHs, 1e-12); close(l.RSI, 2 * l.pHs - 7.5, 1e-12);
  // warmer water is more scale-forming
  assert.ok(S.langelierRyznar({ pH: 7.5, T: 333.15, tds: 500, calcium: 200, alkalinity: 150 }).LSI > l.LSI);
  // SDI = 7 − 2 − 2.301 − 2.5 = 0.199
  const sd = S.stiffDavis({ pH: 7, Ca: 0.01, alkalinity: 0.005, K: 2.5 });
  close(sd.SDI, 7 - 2 - Math.log10(200) - 2.5, 1e-12); close(sd.pHs, 2 + Math.log10(200) + 2.5, 1e-12);
  // SR = 4, ν = 2 → (√4 − 1) = 1: rate = k; second order with SR = 9: (3 − 1)² = 4k
  const d = S.scaleDepositionRate({ SR: 4, k: 1e-7 });
  rel(d.rate, 1e-7, 1e-12); rel(d.thicknessRate, 1e-7 / 2710, 1e-12); assert.equal(d.inhibited, false);
  rel(S.scaleDepositionRate({ SR: 9, k: 1e-7, order: 2, rhoScale: 4480 }).rate, 4e-7, 1e-12);
  close(S.scaleDepositionRate({ SR: 0.8, k: 1e-7 }).rate, 0, 1e-20);
  const inh = S.scaleDepositionRate({ SR: 4, k: 1e-7, inhibitorDose: 10, mic: 5 });
  close(inh.rate, 0, 1e-20); assert.equal(inh.inhibited, true);
  assert.equal(S.scaleDepositionRate({ SR: 4, k: 1e-7, inhibitorDose: 2, mic: 5 }).inhibited, false);
});
