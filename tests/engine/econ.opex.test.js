import test from 'node:test';
import assert from 'node:assert/strict';
import * as opex from '../../src/engine/econ/opex.js';

const close = (a, b, tol = 1e-9, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);

test('OPEX build-up by hand against a controlled reference dataset', () => {
  const r = opex.opexBuildUp({ personnel: [{ count: 10, costPerYear: 100000 }], logistics: [{ days: 100, dayRate: 20000 }],
    maintenancePlanned: 500000, maintenanceUnplanned: 300000, insuranceRate: 0.01, capex: 1e8, leases: 200000, gridCharges: 150000,
    consumables: 50000, inspections: [{ cost: 600000, intervalYears: 3 }], variablePerUnit: 2, production: 1e6 });
  close(r.fixed, 1e6 + 2e6 + 5e5 + 3e5 + 1e6 + 2e5 + 1.5e5 + 5e4 + 2e5); close(r.variable, 2e6);
  close(r.total, r.fixed + r.variable); close(r.breakdown.inspection, 200000);
  close(Object.values(r.breakdown).reduce((s, v) => s + v, 0), r.total, 1e-12);
});

test('availability from MTBF/MTTR and downtime cost', () => {
  const a = opex.availability({ mtbf: 990, mttr: 10 });
  close(a.availability, 0.99); close(a.downtimeHoursPerYear, 87.6, 1e-12); close(a.failuresPerYear, 8.76, 1e-12);
  const d = opex.downtimeCost({ mtbf: 990, mttr: 10, capacityPerHour: 100, price: 50, variableCostSaved: 5 });
  close(d.lostProduction, 8760, 1e-12); close(d.lostRevenue, 438000, 1e-12); close(d.lostMargin, 8760 * 45, 1e-12);
  close(opex.downtimeCost({ availability: 0.95, capacityPerHour: 10, price: 1, recoveryFraction: 0.5 }).lostProduction, 8760 * 0.05 * 10 * 0.5, 1e-12);
});

test('cost of failure and insurance deductible / limit / coinsurance logic', () => {
  const f = opex.costOfFailure({ repairCost: 2e6, downtimeHours: 240, productionPerHour: 500, marginPerUnit: 40, consequenceCost: 1e6, annualFrequency: 0.05 });
  close(f.lostProduction, 4.8e6); close(f.costPerEvent, 7.8e6); close(f.expectedAnnualLoss, 390000);
  assert.deepEqual(opex.insuranceRecovery({ loss: 1000, deductible: 100, limit: 500 }), { insured: 500, retained: 500 });
  assert.deepEqual(opex.insuranceRecovery({ loss: 50, deductible: 100 }), { insured: 0, retained: 50 });
  close(opex.insuranceRecovery({ loss: 1100, deductible: 100, coinsurance: 0.2 }).insured, 800);
});

test('IMR campaign cost: interval x vessel days x day rate', () => {
  const r = opex.imrCost({ length_km: 100, surveyRate_kmPerDay: 10, structures: 5, daysPerStructure: 2, intervalYears: 5, dayRate: 100000,
    mobilisation: 500000, operability: 0.8, years: 20 });
  close(r.vesselDays, 25); close(r.campaignCost, 3e6);
  assert.deepEqual(r.series.map((v, k) => (v > 0 ? k + 1 : 0)).filter(Boolean), [5, 10, 15, 20]);
  close(r.annualAverage, 6e5);
});

test('oil and gas opex: fixed fraction of capex plus per-barrel variable cost', () => {
  const r = opex.oilGasOpex({ capex: 1e9, fixedFraction: 0.04, production_boe: [1e7, 5e6], variablePerBoe: 3, tariffPerBoe: 1 });
  close(r.fixed, 4e7); assert.deepEqual(r.variable, [4e7, 2e7]); close(r.unitCost[0], 8); close(r.unitCost[1], 12);
});

test('escalation series', () => {
  assert.deepEqual(opex.escalateSeries({ base: 100, rate: 0.1, years: 3 }).map((v) => Math.round(v * 100) / 100), [100, 110, 121]);
  close(opex.escalateSeries({ base: 100, rate: [0, 0.1, 0.2], years: 3 })[2], 132, 1e-12);
});

test('wind O&M model: transit erodes CTV working time, SOV is distance-insensitive, costs reconcile', () => {
  const base = { turbines: 80, ratedMW: 10 };
  const near = opex.windOmModel({ ...base, strategy: 'ctv', distance_km: 20 });
  const far = opex.windOmModel({ ...base, strategy: 'ctv', distance_km: 100 });
  // transit time = distance / (22 kn x 1.852)
  close(near.transitHours, 20 / (22 * 1.852), 1e-12);
  close(near.productiveHoursPerShift, 12 - 2 * near.transitHours, 1e-12);
  assert.ok(far.productiveHoursPerShift < near.productiveHoursPerShift && far.annualCost > near.annualCost);
  const sovNear = opex.windOmModel({ ...base, strategy: 'sov', distance_km: 20 }); const sovFar = opex.windOmModel({ ...base, strategy: 'sov', distance_km: 100 });
  close(sovNear.productiveHoursPerShift, sovFar.productiveHoursPerShift); close(sovNear.annualCost, sovFar.annualCost);
  close(Object.values(near.breakdown).reduce((s, v) => s + v, 0), near.annualCost, 1e-12);
  assert.ok(near.availability > 0.8 && near.availability < 1);
  close(near.lostEnergy_MWh, 80 * 10 * 0.45 * near.downtimeHoursPerTurbine, 1e-12);
  assert.equal(opex.windOmModel({ ...base, strategy: 'ctv', distance_km: 250 }).feasible, false);
  assert.equal(near.usesPlaceholderData, true);
  assert.match(opex.defaultWindFailureData.source, /indicative placeholder/);
  assert.match(opex.defaultOmRates.source, /indicative placeholder/);
});

test('wind O&M model responds to user failure data, and a strategy break-even distance exists', () => {
  const data = { subsystems: [{ name: 'only', rate: 2, repairCost: 1000, repairHours: 10, technicians: 2, heavyLift: false }] };
  const r = opex.windOmModel({ turbines: 10, ratedMW: 5, strategy: 'sov', distance_km: 50, failureData: data, meanWaitHours: 0,
    rates: { scheduledServiceHours: 0, scheduledServiceCostPerTurbine: 0 } });
  close(r.downtimeHoursPerTurbine, 20); close(r.breakdown.spares, 20000); close(r.availability, 1 - 20 / 8760, 1e-12);
  assert.equal(r.breakdown.heavyLift, 0);
  const be = opex.windOmStrategyBreakEven({ base: { turbines: 100, ratedMW: 12 }, energyPrice: 80, distances_km: [10, 30, 50, 70, 90, 110, 130, 150] });
  assert.ok(be.ctv[0] < be.sov[0]);                      // close to port the CTV strategy is cheaper with the placeholder data
  assert.ok(be.breakEven_km > 10 && be.breakEven_km <= 150);
});

test('life extension versus decommissioning now: hand NPV', () => {
  const r = opex.lifeExtensionCashflows({ upgradeCapex: 50, netOperatingCashflow: [40, 30, 20], decomCost: 100, rate: 0.10 });
  assert.deepEqual(r.extend, [-50, 40, 30, -80]); assert.deepEqual(r.decommissionNow, [-100, 0, 0, 0]);
  close(r.npvExtend, -50 + 40 / 1.1 + 30 / 1.21 - 80 / 1.331, 1e-13); close(r.npvDecommissionNow, -100);
  assert.equal(r.decision, 'extend'); close(r.npvIncremental, r.npvExtend + 100, 1e-13);
});
