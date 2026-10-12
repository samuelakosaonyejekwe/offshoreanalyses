import test from 'node:test';
import assert from 'node:assert/strict';
import * as L from '../../src/engine/econ/lcoe.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const pvA = (i, n) => (1 - (1 + i) ** -n) / i; // present worth of an annuity (P/A, i, n)

test('LCOE with level costs and energy equals (CRF × capex + opex)/energy', () => {
  // capex 1000, opex 50/yr, 100 MWh/yr, 20 yr, 8 %: (P/A) = 9.8181, CRF = 0.101852 → LCOE = (101.852 + 50)/100 = 1.5185
  const r = L.lcoe({ capex: 1000, opex: 50, energy: 100, rate: 0.08, years: 20 });
  const crf = 1 / pvA(0.08, 20);
  close(r.lcoe, (1000 * crf + 50) / 100, 1e-10); close(r.lcoe, 1.5185, 1e-4);
  close(r.pvEnergy, 100 * pvA(0.08, 20), 1e-8); close(r.pvCost, 1000 + 50 * pvA(0.08, 20), 1e-8);
  close(r.components.capex, 10 * crf, 1e-10); close(r.components.opex, 0.5, 1e-10);
  close(Object.values(r.shares).reduce((s, v) => s + v, 0), 1, 1e-10);
  // the same with explicit series
  const n = 20; const ser = (v) => [0, ...new Array(n).fill(v)];
  close(L.lcoe({ capex: [1000, ...new Array(n).fill(0)], opex: ser(50), energy: ser(100), rate: 0.08 }).lcoe, r.lcoe, 1e-10);
  // zero discount rate: total cost over total energy
  close(L.lcoe({ capex: 1000, opex: 50, energy: 100, rate: 0, years: 20 }).lcoe, (1000 + 1000) / 2000, 1e-12);
  // fixed-charge-rate form with FCR = CRF is identical
  const f = L.lcoeFixedChargeRate({ capex: 1000, rate: 0.08, life: 20, fixedOm: 50, annualEnergy: 100 });
  close(f.lcoe, r.lcoe, 1e-10); close(f.fcr, crf, 1e-12);
  close(L.lcoeFixedChargeRate({ capex: 1000, fcr: 0.1, fixedOm: 50, annualEnergy: 100, variableOm: 0.2, fuel: 0.3 }).lcoe, 1.5 + 0.5, 1e-12);
});

test('LCOE with fuel, decommissioning, salvage and tax; real against nominal', () => {
  // two operating years, 10 %: capex 100 at t0; opex 10, fuel 5 each year; decom 20 and salvage 8 at t2; energy 50, 50
  const r = L.lcoe({ capex: [100, 0, 0], opex: [0, 10, 10], fuel: [0, 5, 5], decom: [0, 0, 20], salvage: [0, 0, 8], energy: [0, 50, 50], rate: 0.1 });
  const d1 = 1 / 1.1; const d2 = 1 / 1.21;
  const pvE = 50 * (d1 + d2);
  close(r.lcoe, (100 + 15 * (d1 + d2) + 12 * d2) / pvE, 1e-10);
  close(r.components.salvage, -8 * d2 / pvE, 1e-10); close(r.components.fuel, 5 * (d1 + d2) / pvE, 1e-10);
  // post-tax form: [capex − T PV(dep) + (1 − T) PV(opex)] / [(1 − T) PV(E)], T = 25 %, straight-line depreciation 50, 50
  const t = L.lcoe({ capex: [100, 0, 0], opex: [0, 10, 10], energy: [0, 50, 50], rate: 0.1, taxRate: 0.25, depreciation: [0, 50, 50] });
  close(t.lcoe, (100 - 0.25 * 50 * (d1 + d2) + 0.75 * 10 * (d1 + d2)) / (0.75 * pvE), 1e-10);
  assert.ok(t.components.taxShield < 0);
  // mid-year discounting moves every operating flow half a year earlier: energy and opex scale alike, capex does not
  const m = L.lcoe({ capex: 1000, opex: 50, energy: 100, rate: 0.08, years: 20, midYear: true });
  close(m.lcoe, (1000 / (100 * pvA(0.08, 20) * Math.sqrt(1.08))) + 0.5, 1e-9);
  // real 5 %, inflation 2 %: nominal rate 7.1 %; the flat nominal price exceeds the real LCOE
  const rn = L.lcoeRealNominal({ capex: 1000, opex: 50, energy: 100, rate: 0.05, years: 20, inflation: 0.02 });
  close(rn.nominalRate, 1.05 * 1.02 - 1, 1e-12);
  close(rn.real, (1000 / pvA(0.05, 20) + 50) / 100, 1e-10);
  // flat nominal price P with PV at 7.1 % of P × 100 = PV of real costs → P = real × (P/A, 5 %)/(P/A, 7.1 %)
  close(rn.nominal, rn.real * pvA(0.05, 20) / pvA(0.071, 20), 1e-9);
  assert.throws(() => L.lcoe({ capex: 1000, energy: 0, rate: 0.08, years: 5 }), RangeError);
});

test('levelised cost of storage and of hydrogen by hand', () => {
  // 10 MWh usable, 300 cycles/yr, 85 % round trip, charge at 40/MWh, capex 3 M, O&M 30 k/yr, 10 yr, 7 %, no fade
  const s = L.lcos({ capex: 3e6, omPerYear: 30e3, capacity_MWh: 10, cyclesPerYear: 300, roundTripEfficiency: 0.85, chargePrice: 40, life: 10, rate: 0.07 });
  const E = 3000; const a = pvA(0.07, 10);
  close(s.pvEnergy, E * a, 1e-6);
  close(s.lcos, (3e6 + (30e3 + 40 * E / 0.85) * a) / (E * a), 1e-8);
  close(s.lcos, 3e6 / (E * a) + 10 + 40 / 0.85, 1e-8); // capex + O&M + charging: 142.38 + 10 + 47.06
  // a replacement in year 5 and an end-of-life cost add their present values
  const s2 = L.lcos({ capex: 3e6, capacity_MWh: 10, cyclesPerYear: 300, roundTripEfficiency: 0.85, life: 10, rate: 0.07, replacements: [{ year: 5, cost: 5e5 }], endOfLife: 1e5 });
  close(s2.pvCost, 3e6 + 5e5 / 1.07 ** 5 + 1e5 / 1.07 ** 10, 1e-4);
  // degradation 2 %/yr reduces discharged energy geometrically
  const s3 = L.lcos({ capex: 3e6, capacity_MWh: 10, cyclesPerYear: 300, roundTripEfficiency: 0.85, life: 2, rate: 0, degradation: 0.02 });
  close(s3.pvEnergy, 3000 * (1 + 0.98), 1e-8);
  // 10 MW electrolyser, 50 % load factor, 50 kWh/kg: 10 × 8760 × 0.5 × 1000/50 = 876 000 kg/yr
  const h = L.lcoh({ capex: 10e6, capacityMW: 10, loadFactor: 0.5, specificConsumption_kWhPerKg: 50, electricityPrice: 40, fixedOmPerYear: 2e5, waterCostPerKg: 0.02, life: 20, rate: 0.06 });
  close(h.annualOutput_kg, 876000, 1e-6);
  // electricity 50 kWh/kg × 0.04 per kWh = 2.00 per kg; capex 10e6 CRF/876000; O&M 0.2283; water 0.02
  close(h.lcoh, 10e6 / pvA(0.06, 20) / 876000 + 2e5 / 876000 + 2 + 0.02, 1e-8);
  close(h.components.electricity ?? h.components.energy ?? 2, 2, 1e-8);
});

test('delivered energy, cost of service, LACE and blended revenue', () => {
  // generation 50/MWh, 1 TWh/yr, 4 % losses; transmission 300 M over 30 yr at 6 % plus 3 M/yr
  const d = L.lcoeDelivered({ lcoeGeneration: 50, annualEnergy: 1e6, lossesFrac: 0.04, transmissionCapex: 300e6, transmissionOmPerYear: 3e6, rate: 0.06, life: 30 });
  close(d.deliveredEnergy, 0.96e6, 1e-6); close(d.generation, 50 / 0.96, 1e-10);
  close(d.transmission, (300e6 / pvA(0.06, 30) + 3e6) / 0.96e6, 1e-8); close(d.lcode, d.generation + d.transmission, 1e-10);
  // costs 100 at t0 and 10, 10; quantities 0, 20, 20 at 10 %
  const s = L.levelisedCostOfService({ costs: [100, 10, 10], quantities: [0, 20, 20], rate: 0.1 });
  close(s.levelisedCost, (100 + 10 * (1 / 1.1 + 1 / 1.21)) / (20 * (1 / 1.1 + 1 / 1.21)), 1e-10);
  // LACE: flat marginal price 60 → energy value 60; capacity 100 MW × 20 % credit × 50 k/MW·yr = 1 M/yr over 200 000 MWh = 5
  const en = [0, 2e5, 2e5, 2e5];
  const la = L.lace({ energy: en, marginalPrice: 60, capacityMW: 100, capacityCredit: 0.2, capacityPaymentPerMWyr: 50e3, rate: 0.07, lcoe: 58 });
  close(la.energyValue, 60, 1e-9); close(la.capacityValue, 5, 1e-9); close(la.lace, 65, 1e-9); close(la.netValue, 7, 1e-9);
  // 100 MWh/yr for 3 years: 50 % PPA at 80, 30 % CfD strike 90 against a 60 market, remainder merchant at 60
  const b = L.blendedRevenue({ energy: [100, 100, 100], merchant: 60, ppa: { share: 0.5, price: 80 }, cfd: { share: 0.3, strike: 90 }, rate: 0.1 });
  // PPA 4000; CfD: market 30 × 60 = 1800 + top-up 30 × 30 = 900 → 2700; merchant 20 × 60 = 1200: total 7900 → 79/MWh
  close(b.revenue[0], 7900, 1e-9); close(b.averagePrice, 79, 1e-9); close(b.levelisedPrice, 79, 1e-9);
  close(b.byStream.ppa[1], 4000, 1e-9); close(b.byStream.cfd[1], 2700, 1e-9); close(b.byStream.merchant[2], 1200, 1e-9);
  close(b.pv, 7900 * pvA(0.1, 3), 1e-8); close(b.contractedShareOfPv, 6700 / 7900, 1e-10);
  // a one-year PPA reverts to merchant afterwards
  close(L.blendedRevenue({ energy: [100, 100], merchant: 60, ppa: { share: 1, price: 80, termYears: 1 }, rate: 0 }).revenue[1], 6000, 1e-9);
  assert.throws(() => L.blendedRevenue({ energy: [100], merchant: 60, ppa: { share: 0.8, price: 80 }, cfd: { share: 0.5, strike: 90 }, rate: 0 }), RangeError);
});

test('sensitivity grid, learning projection and capacity-factor energy', () => {
  const base = { capex: 1000, opex: 50, energy: 100, rate: 0.08, years: 20 };
  const g = L.lcoeSensitivityGrid({ base, capexFactors: [0.8, 1, 1.2], capacityFactorFactors: [0.9, 1], rates: [0.06, 0.08] });
  const b0 = L.lcoe(base).lcoe;
  close(g.base, b0, 1e-12); close(g.z[1][1], b0, 1e-10); assert.equal(g.z.length, 2); assert.equal(g.z[0].length, 3);
  // energy × 0.9 → LCOE / 0.9; capex × 1.2 raises only the capex component
  close(g.z[0][1], b0 / 0.9, 1e-10);
  close(g.z[1][2], (1200 / pvA(0.08, 20) + 50) / 100, 1e-10);
  close(g.lcoeVsRate[1], b0, 1e-10); assert.ok(g.lcoeVsRate[0] < g.lcoeVsRate[1]);
  // 15 % learning rate, capacity doubling twice (growth 100 %/yr for 2 years): 0.85² = 0.7225; half the cost learns
  const p = L.lcoeLearningProjection({ lcoe0: 100, year0: 2025, year: 2027, learningRate: 0.15, capacityGrowth: 1 });
  close(p.doublings, 2, 1e-12); close(p.capacityRatio, 4, 1e-12); close(p.lcoe, 72.25, 1e-9); assert.equal(p.series.length, 3);
  close(L.lcoeLearningProjection({ lcoe0: 100, year0: 2025, year: 2027, learningRate: 0.15, capacityGrowth: 1, learningShare: 0.5 }).lcoe, 50 * 0.7225 + 50, 1e-9);
  // 100 MW at 45 %: 394 200 MWh in year 1, falling 0.5 % per year
  const e = L.energyFromCapacityFactor({ capacityMW: 100, capacityFactor: 0.45, degradation: 0.005, years: 3 });
  close(e[0], 394200, 1e-6); close(e[2], 394200 * 0.995 ** 2, 1e-6); assert.equal(e.length, 3);
  assert.equal(typeof L.requiredStrikePrice, 'function');
});

test('required strike price: the goal-seek price delivers the target equity IRR; a higher target needs a higher price', () => {
  const inputs = { timeline: { constructionYears: 2, operatingYears: 12 }, inflation: 0.02, discountRate: 0.08,
    capex: [{ name: 'plant', amount: 1000, phasing: [0.4, 0.6] }], opex: { fixed: 30 }, production: 110, price: { type: 'indexed', price: 4 },
    tax: { rate: 0.25, depreciation: { method: 'sl', life: 10 } } };
  const debt = { gearing: 0.6, rate: 0.05, tenor: 8, type: 'annuity' };
  const a = L.requiredStrikePrice({ inputs, debt, targetEquityIrr: 0.10, lo: 0.5, hi: 20 });
  close(a.equityIrr, 0.10, 1e-5); close(a.gearing, 0.6, 1e-9); assert.ok(a.price > 0.5 && a.price < 20 && a.minDscr > 0);
  const b = L.requiredStrikePrice({ inputs, debt, targetEquityIrr: 0.14, lo: 0.5, hi: 20 });
  close(b.equityIrr, 0.14, 1e-5); assert.ok(b.price > a.price);
  // a flat nominal price must start higher than an inflation-indexed one for the same return
  const flat = L.requiredStrikePrice({ inputs, debt, targetEquityIrr: 0.10, indexed: false, lo: 0.5, hi: 20 });
  assert.ok(flat.price > a.price, `flat ${flat.price} vs indexed ${a.price}`);
});
