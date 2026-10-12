import test from 'node:test';
import assert from 'node:assert/strict';
import * as C from '../../src/engine/econ/carbon.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);

test('embodied carbon: Σ mass × factor with item factors, table overrides and ranges', () => {
  // 1000 t at 2.0 kg/kg + 500 t at 0.15 kg/kg = 2000 + 75 = 2075 t CO2e
  const r = C.embodiedCarbon({ items: [{ name: 'jacket', mass_t: 1000, factor: 2.0 }, { name: 'grout', mass_t: 500, factor: 0.15 }] });
  close(r.total_t, 2075, 1e-9); assert.equal(r.lines.length, 2); close(r.lines[1].emissions_t, 75, 1e-12);
  // a material key resolved through a supplied table (number or { value })
  const t = C.embodiedCarbon({ items: [{ material: 'mySteel', mass_t: 10 }, { material: 'myAlloy', mass_t: 2 }], factors: { mySteel: 1.5, myAlloy: { value: 8 } } });
  close(t.total_t, 15 + 16, 1e-9);
  // default table entries are indicative and carry a range that brackets the value
  const key = Object.keys(C.defaultFactors)[0];
  const d = C.embodiedCarbon({ items: [{ material: key, mass_t: 100 }] });
  close(d.total_t, 100 * C.defaultFactors[key].value, 1e-9); assert.equal(d.usesIndicativeFactors, true);
  assert.ok(d.low_t <= d.total_t && d.high_t >= d.total_t);
  assert.throws(() => C.embodiedCarbon({ items: [{ material: 'unobtainium', mass_t: 1 }] }), RangeError);
});

test('fabrication, transport and marine-operations emissions by hand', () => {
  // 2000 t × 500 kWh/t = 1000 MWh at 400 g/kWh = 400 t
  const f = C.fabricationEmissions({ mass_t: 2000, energy_kWhPerT: 500, gridIntensity_gPerKWh: 400 });
  close(f.energy_MWh, 1000, 1e-9); close(f.emissions_t, 400, 1e-9);
  // 5000 t × 2000 km × 0.01 kg/tkm = 100 t
  const tr = C.transportEmissions({ legs: [{ mass_t: 5000, distance_km: 2000, mode: 'ship', factor_kgPerTkm: 0.01 }, { mass_t: 100, distance_km: 300, mode: 'truck', factor_kgPerTkm: 0.1 }] });
  close(tr.tkm, 1e7 + 3e4, 1e-6); close(tr.total_t, 100 + 3, 1e-9);
  // IMO conversion factors are exact
  assert.equal(C.FUEL_CO2_FACTORS.MGO, 3.206); assert.equal(C.FUEL_CO2_FACTORS.HFO, 3.114); assert.equal(C.FUEL_CO2_FACTORS.LNG, 2.75);
  // 2 vessels × 30 days × 20 t/day MGO = 1200 t fuel → 3847.2 t CO2; 10 days × 15 t/day LNG = 150 t → 412.5 t
  const m = C.marineOpsEmissions({ vessels: [{ vessel: 'hlv', days: 30, fuel_tPerDay: 20, count: 2 }, { vessel: 'psv', days: 10, fuel_tPerDay: 15, fuel: 'LNG' }] });
  close(m.fuel_t, 1350, 1e-9); close(m.total_t, 1200 * 3.206 + 150 * 2.75, 1e-9); assert.equal(m.usesPlaceholderConsumption, false);
  close(C.marineOpsEmissions({ vessels: [{ vessel: 'hlv', days: 30, fuel_tPerDay: 20 }], co2eUplift: 0.02 }).total_t, 600 * 3.206 * 1.02, 1e-9);
});

test('operation emissions: fuel, grid import, flaring and methane with an explicit GWP', () => {
  // 100 t MGO → 320.6; 2000 MWh at 250 g/kWh → 500; 40 t flared at 2.75 → 110; 2 t CH4 × 29.8 → 59.6; other 10
  const o = C.operationEmissions({ logisticsFuel_t: 100, gridImport_MWh: 2000, gridIntensity_gPerKWh: 250, flaredGas_t: 40, fugitiveMethane_t: 2, other_t: 10 });
  close(o.total_t, 320.6 + 500 + 110 + 59.6 + 10, 1e-9);
  close(Object.values(o.breakdown).reduce((s, v) => s + v, 0), o.total_t, 1e-9);
  close(C.operationEmissions({ fugitiveMethane_t: 2, methaneGwp: 82.5 }).total_t, 165, 1e-9);
  close(C.operationEmissions({ flaredGas_t: 10, flareFactor_tCO2PerT: 3 }).total_t, 30, 1e-9);
});

test('end-of-life credit (net scrap), life-cycle roll-up and intensity', () => {
  // 1000 t, recovery 90 %, recycled content 20 %: net scrap 700 t × (2.3 − 0.5) = 1260 t credit
  const e = C.endOfLifeCredit({ mass_t: 1000, recoveryRate: 0.9, recycledContent: 0.2, primaryFactor: 2.3, secondaryFactor: 0.5 });
  close(e.netScrap_t, 700, 1e-9); close(e.credit_t, 1260, 1e-9);
  assert.ok(C.endOfLifeCredit({ mass_t: 1000, recoveryRate: 0.1, recycledContent: 0.6, primaryFactor: 2.3, secondaryFactor: 0.5 }).credit_t < 0);
  const l = C.lifeCycleGwp({ product_t: 5000, transport_t: 300, installation_t: 700, operationPerYear_t: 100, operatingYears: 25, endOfLife_t: 500, credit_t: 1260 });
  close(l.total_t, 5000 + 300 + 700 + 2500 + 500, 1e-9); close(l.totalWithCredit_t, 9000 - 1260, 1e-9);
  close(Object.values(l.shares).reduce((s, v) => s + v, 0), 1, 1e-12);
  close(C.lifeCycleGwp({ product_t: 10, operationTotal_t: 30 }).total_t, 40, 1e-12);
  // 9000 t over 1.5 TWh = 6 g/kWh; 9000 t over 3e6 boe = 3 kg/boe
  const i = C.carbonIntensity({ emissions_t: 9000, energy_MWh: 1.5e6 });
  close(i.gPerKWh, 6, 1e-12); assert.equal(i.kgPerBoe, null);
  close(C.carbonIntensity({ emissions_t: 9000, production_boe: 3e6 }).kgPerBoe, 3, 1e-12);
});

test('avoided emissions, carbon and energy payback, EROI', () => {
  // 60 000 MWh/yr against 400 g/kWh avoids 24 000 t/yr; life-cycle 9000 t → payback 0.375 yr; 25 yr: net 591 000 t
  const a = C.avoidedEmissions({ annualEnergy_MWh: 60000, referenceIntensity_gPerKWh: 400, lifeCycleEmissions_t: 9000, operatingYears: 25, embodiedEnergy_MWh: 30000, annualEnergyInput_MWh: 600 });
  close(a.annualAvoided_t, 24000, 1e-9); close(a.lifetimeAvoided_t, 600000, 1e-6); close(a.netAvoided_t, 591000, 1e-6);
  close(a.carbonPaybackYears, 0.375, 1e-12); close(a.projectIntensity_gPerKWh, 6, 1e-12);
  // energy payback = embodied / net annual output = 30000/59400; EROI = 1.5e6/(30000 + 15000) = 33.33
  close(a.energyPaybackYears, 30000 / 59400, 1e-9); close(a.eroi, 1.5e6 / 45000, 1e-9);
  assert.equal(C.avoidedEmissions({ annualEnergy_MWh: 1, referenceIntensity_gPerKWh: 400, lifeCycleEmissions_t: 1, operatingYears: 1 }).eroi, null);
});

test('carbon cost, marginal abatement cost and the abatement curve', () => {
  // emissions 100, 100, 100 t at 50 per t, 10 %: PV = 5000 (1 + 1/1.1 + 1/1.21) = 13677.69
  const c = C.carbonCost({ emissions_t: [100, 100, 100], price: 50, rate: 0.1, npv: 20000 });
  close(c.total, 15000, 1e-9); close(c.pv, 5000 * (1 + 1 / 1.1 + 1 / 1.21), 1e-9); close(c.npvWithCarbon, 20000 - c.pv, 1e-9);
  close(C.carbonCost({ emissions_t: [0, 10], price: [40, 60], rate: 0 }).pv, 600, 1e-12);
  // alternative gives up 2 M of NPV to avoid 40 000 t → 50 per tonne
  const m = C.marginalAbatementCost({ npvBase: 10e6, npvAlternative: 8e6, emissionsBase_t: 100000, emissionsAlternative_t: 60000 });
  close(m.mac, 50, 1e-9); close(m.abatement_t, 40000, 1e-9); close(m.cost, 2e6, 1e-6);
  const curve = C.abatementCurve({ measures: [{ name: 'electrify', abatement_t: 1000, cost: 80000 }, { name: 'led', abatement_t: 100, cost: -5000 }, { name: 'ccs', abatement_t: 5000, cost: 1e6 }], carbonPrice: 100 });
  assert.deepEqual(curve.curve.map((x) => x.name), ['led', 'electrify', 'ccs']);
  close(curve.curve[0].mac, -50, 1e-9); close(curve.curve[1].mac, 80, 1e-9); close(curve.curve[2].mac, 200, 1e-9);
  close(curve.curve[1].cumulativeFrom_t, 100, 1e-9); close(curve.curve[2].cumulativeTo_t, 6100, 1e-9);
  close(curve.totalAbatement_t, 6100, 1e-9); close(curve.worthwhileAbatement_t, 1100, 1e-9);
  assert.equal(curve.curve[2].worthwhile, false);
});

test('material circularity indicator: fully linear product 0.1, fully circular product 1', () => {
  // all virgin, nothing collected: V = M, W = M, LFI = 1, MCI = 1 − 0.9 = 0.1
  const lin = C.circularity({ mass_t: 100 });
  close(lin.lfi, 1, 1e-12); close(lin.mci, 0.1, 1e-12); close(lin.utility, 1, 1e-12);
  // all recycled feedstock, all collected, lossless recycling: V = 0, W = 0 → LFI = 0, MCI = 1
  const circ = C.circularity({ mass_t: 100, recycledFeedstock: 1, collectedForRecycling: 1, recyclingEfficiency: 1, feedstockRecyclingEfficiency: 1 });
  close(circ.lfi, 0, 1e-12); close(circ.mci, 1, 1e-12);
  // linear product lasting twice the average: MCI = 1 − 0.9/2 = 0.55
  close(C.circularity({ mass_t: 100, lifeRatio: 2 }).mci, 0.55, 1e-12);
  // recyclable mass fraction: (800 × 0.95 + 200 × 0.1)/1000 = 0.78
  close(C.circularity({ items: [{ mass_t: 800, recyclableFraction: 0.95 }, { mass_t: 200, recyclableFraction: 0.1 }] }).recyclableFraction, 0.78, 1e-12);
});

test('footprint, water, checklist flags, taxonomy threshold and weighted score', () => {
  const f = C.footprintAndWater({ footprint: [{ name: 'mudmat', area_m2: 100, count: 4 }, { name: 'anchor spread', area_m2: 5000, temporary: true }], water: [{ name: 'hydrotest', quantity: 2, m3PerUnit: 1500 }] });
  close(f.permanent_m2, 400, 1e-12); close(f.temporary_m2, 5000, 1e-12); close(f.total_m2, 5400, 1e-12); close(f.water_m3, 3000, 1e-12);
  const first = C.ENVIRONMENTAL_CHECKLIST[0];
  const flags = C.environmentalFlags({ [first.key]: true });
  assert.equal(flags.length, C.ENVIRONMENTAL_CHECKLIST.length);
  assert.equal(flags.find((x) => x.key === first.key).flagged, true);
  assert.equal(flags.filter((x) => x.flagged).length, 1);
  const tx = C.taxonomyThresholdCheck({ intensity_gPerKWh: 12 });
  assert.equal(tx.pass, true); close(tx.margin_gPerKWh, 88, 1e-12); close(tx.utilisation, 0.12, 1e-12);
  assert.equal(C.taxonomyThresholdCheck({ intensity_gPerKWh: 350, threshold_gPerKWh: 270 }).pass, false);
  assert.ok(Array.isArray(C.SDG_MAPPING) && C.SDG_MAPPING.length > 0);
  // carbon 30 between worst 100 and best 0 → 0.7 (weight 2); jobs 50 between 0 and 200 → 0.25 (weight 1): score 1.65/3
  const s = C.sustainabilityScore({ criteria: [{ name: 'carbon', value: 30, best: 0, worst: 100, weight: 2 }, { name: 'jobs', value: 50, best: 200, worst: 0, weight: 1 }] });
  close(s.items[0].normalised, 0.7, 1e-12); close(s.items[1].normalised, 0.25, 1e-12); close(s.score, 0.55, 1e-12);
  close(C.sustainabilityScore({ criteria: [{ name: 'x', value: 500, best: 0, worst: 100, weight: 1 }] }).score, 0, 1e-12);
});
