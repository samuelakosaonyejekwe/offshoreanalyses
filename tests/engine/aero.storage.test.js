import test from 'node:test';
import assert from 'node:assert/strict';
import {
  batteryStep, runBattery, batteryDegradation, electrolyser, dispatchFirm, dispatchCurtailment,
  dispatchArbitrage, complementarity,
} from '../../src/engine/aero/storage.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);

test('battery step by hand: efficiency, power rating and state-of-charge limits', () => {
  const battery = { capacityWh: 10e6, cRate: 0.5, etaCharge: 0.9, etaDischarge: 0.95, socMin: 0.1, socMax: 0.9 };
  // charge request 8 MW is limited to 5 MW (0.5 C); stored +4.5 MWh, loss 0.5 MWh
  const c = batteryStep({ request: -8e6, soc: 0.3, battery });
  close(c.power, -5e6, 1e-12);
  close(c.storedDeltaWh, 4.5e6, 1e-12);
  close(c.soc, 0.75, 1e-12);
  close(c.lossWh, 0.5e6, 1e-12);
  // next hour: only 1.5 MWh of room → terminal charge 1.5/0.9 = 1.6667 MW
  const c2 = batteryStep({ request: -8e6, soc: c.soc, battery });
  close(c2.power, -1.5e6 / 0.9, 1e-12);
  close(c2.soc, 0.9, 1e-12);
  // discharge: 3 MW for 1 h draws 3/0.95 MWh from the store
  const d = batteryStep({ request: 3e6, soc: 0.9, battery });
  close(d.storedDeltaWh, -3e6 / 0.95, 1e-12);
  // limited by the lower SoC bound: (0.15 − 0.1)·10 MWh·0.95 = 0.475 MW
  close(batteryStep({ request: 3e6, soc: 0.15, battery }).power, 0.475e6, 1e-12);
  assert.equal(batteryStep({ request: 3e6, soc: 0.1, battery }).power, 0);
});

test('battery energy balance: energy in × η_c = energy out / η_d + Δ stored', () => {
  const n = 500;
  const request = Float64Array.from({ length: n }, (_, i) => 6e6 * Math.sin(i * 0.37) + 2e6 * Math.cos(i * 1.3));
  for (const battery of [
    { capacityWh: 20e6, cRate: 0.25, etaCharge: 0.92, etaDischarge: 0.94, soc0: 0.5 },
    { capacityWh: 20e6, cRate: 1, roundTripEfficiency: 0.85, soc0: 0.2 },
    { capacityWh: 20e6, cRate: 0.5, etaCharge: 0.88, etaDischarge: 1, soc0: 0.5 },
  ]) {
    const r = runBattery({ request, battery, dtHours: 0.5 });
    const etaC = battery.etaCharge ?? Math.sqrt(battery.roundTripEfficiency);
    const etaD = battery.etaDischarge ?? Math.sqrt(battery.roundTripEfficiency);
    close(r.chargedWh * etaC, r.dischargedWh / etaD + r.storedDeltaWh, 1e-10);
    // first law at the terminals: in = out + losses + Δ stored
    close(r.chargedWh, r.dischargedWh + r.lossWh + r.storedDeltaWh, 1e-10);
    for (let i = 0; i <= n; i++) assert.ok(r.soc[i] >= 0.1 - 1e-12 && r.soc[i] <= 0.9 + 1e-12);
    for (let i = 0; i < n; i++) assert.ok(Math.abs(r.power[i]) <= battery.cRate * 20e6 + 1e-6);
    close(r.equivalentFullCycles, (r.chargedWh * etaC) / 20e6, 1e-10);
  }
  // all loss on charging (η_d = 1): energy in × η = energy out + ΔSoC·capacity
  const r = runBattery({ request, battery: { capacityWh: 20e6, cRate: 0.5, etaCharge: 0.88, etaDischarge: 1, soc0: 0.5 } });
  close(r.chargedWh * 0.88, r.dischargedWh + (r.soc[n] - r.soc[0]) * 20e6, 1e-10);
});

test('degradation law and electrolyser by hand', () => {
  // 1.5 %/yr calendar + 0.004 %/cycle × 300 cycles/yr = 2.7 %/yr → 20 % fade after 7.41 years
  const d = batteryDegradation({ years: 5, cyclesPerYear: 300, calendarFadePerYear: 0.015, cycleFadePerCycle: 4e-5 });
  close(d.capacityFraction, 1 - 5 * 0.027, 1e-12);
  close(d.lifeYears, 0.2 / 0.027, 1e-12);
  close(d.lifeCycles, (0.2 / 0.027) * 300, 1e-12);
  // 10 MW electrolyser, 52 kWh/kg: full load 1 h → 192.3 kg; 0.5 MW is below the 10 % minimum load
  const e = electrolyser({ power: [10e6, 15e6, 5e6, 0.5e6], ratedPower: 10e6, specificConsumption: 52, waterMargin: 1.5 });
  close(e.hydrogenKg, (10000 + 10000 + 5000) / 52, 1e-12);
  close(e.energyUsedWh, 25e6, 1e-12);
  close(e.energyUnusedWh, 5e6 + 0.5e6, 1e-12);
  close(e.onFraction, 0.75, 1e-12);
  // water: stoichiometric 8.936 kg/kg × 1.5
  close(e.waterLitres / e.hydrogenKg, (18.01528 / 2.01588) * 1.5, 1e-12);
  close(e.waterLitres / e.hydrogenKg / 1.5, 8.937, 1e-3);
  // LHV efficiency 33.33/52 = 64.1 %
  close(e.efficiencyLhv, 33.33 / 52, 1e-12);
  // part-load penalty: +10 % consumption at half load
  const pl = electrolyser({ power: [5e6], ratedPower: 10e6, specificConsumption: 52, consumptionFactor: (l) => 1 + 0.2 * (1 - l) });
  close(pl.hydrogenKg, 5000 / (52 * 1.1), 1e-12);
  close(pl.meanSpecificConsumption, 57.2, 1e-12);
});

test('firm-power dispatch: storage raises firmness; delivered energy is conserved', () => {
  const n = 240;
  const generation = Float64Array.from({ length: n }, (_, i) => (i % 24 < 12 ? 9e6 : 1e6));
  const battery = { capacityWh: 60e6, cRate: 0.5, etaCharge: 0.95, etaDischarge: 0.95, socMin: 0, socMax: 1, soc0: 0 };
  const r = dispatchFirm({ generation, target: 4.5e6, battery });
  close(r.firmnessWithoutStorage, 0.5, 1e-12);
  // each day: 12 h × 4.5 MW surplus → 51.3 MWh stored; 12 h × 3.5 MW deficit needs 44.2 MWh from the store
  assert.ok(r.firmness > 0.95, `firmness ${r.firmness}`);
  assert.ok(r.delivered.every((p) => p <= 4.5e6 + 1e-6));
  // energy bookkeeping: generation = delivered + surplus spilled + battery net charge at the terminals
  let gen = 0;
  let net = 0;
  for (let i = 0; i < n; i++) {
    gen += generation[i];
    net += -r.battery[i];
  }
  close(gen, r.deliveredWh + r.surplusWh + net, 1e-9);
  close(r.energyShiftedWh, r.battery.reduce((s, p) => s + Math.max(0, p), 0), 1e-12);
  // a target above the mean generation cannot be firm
  const hard = dispatchFirm({ generation, target: 8e6, battery });
  assert.ok(hard.firmness < 0.6 && hard.shortfallWh > 0);
  const exp = dispatchFirm({ generation, target: 4.5e6, battery: { ...battery, capacityWh: 10e6 }, exportSurplus: true });
  assert.ok(Math.max(...exp.delivered) > 4.5e6);
});

test('curtailment absorption at an export limit', () => {
  const n = 96;
  const generation = Float64Array.from({ length: n }, (_, i) => (i % 24 < 6 ? 12e6 : 4e6));
  const battery = { capacityWh: 20e6, cRate: 0.5, etaCharge: 0.95, etaDischarge: 0.95, socMin: 0, socMax: 1, soc0: 0 };
  const r = dispatchCurtailment({ generation, exportLimit: 8e6, battery });
  // without storage: 4 MW × 6 h × 4 days = 96 MWh curtailed
  close(r.curtailedWithoutWh, 96e6, 1e-12);
  // the battery (20 MWh, 10 MW) absorbs 20/0.95 = 21.05 MWh per day at the terminals, then is full
  close(r.curtailmentAvoidedWh, 4 * (20e6 / 0.95), 1e-9);
  close(r.curtailedWithWh, 96e6 - 4 * (20e6 / 0.95), 1e-9);
  // and returns 20·0.95 = 19 MWh per day into the 4 MW headroom
  close(r.energyShiftedWh, 4 * 19e6, 1e-9);
  assert.ok(r.exported.every((p) => p <= 8e6 + 1e-6));
  close(r.equivalentFullCycles, 4, 1e-9);
  let gen = 0;
  for (let i = 0; i < n; i++) gen += generation[i];
  // generation = exported + curtailed + conversion losses (store returns to empty)
  close(gen, r.exportedWh + r.curtailedWithWh + 4 * (20e6 / 0.95 - 19e6), 1e-9);
  close(r.cableUtilisation, r.exportedWh / (8e6 * n), 1e-12);
});

test('price arbitrage: greedy threshold rule and linear programme', () => {
  const prices = Array.from({ length: 48 }, (_, i) => 60 + 40 * Math.sin((2 * Math.PI * (i - 9)) / 24) + (i % 5));
  const battery = { capacityWh: 8e6, cRate: 0.5, etaCharge: 0.95, etaDischarge: 0.95, socMin: 0.1, socMax: 0.9, soc0: 0.1 };
  const g = dispatchArbitrage({ prices, battery });
  assert.equal(g.method, 'threshold');
  assert.ok(g.revenue > 0);
  close(g.revenue, g.power.reduce((s, p, i) => s + (p / 1e6) * prices[i], 0), 1e-12);
  const lp = dispatchArbitrage({ prices, battery, method: 'lp' });
  assert.ok(lp.revenue >= g.revenue - 1e-9);
  for (const r of [g, lp]) {
    for (const s of r.soc) assert.ok(s >= 0.1 - 1e-9 && s <= 0.9 + 1e-9);
    for (const p of r.power) assert.ok(Math.abs(p) <= 4e6 + 1e-6);
  }
  // upper bound on revenue: usable energy 6.4 MWh cycled twice a day between the cheapest and dearest hours
  const pmax = Math.max(...prices);
  const pmin = Math.min(...prices);
  assert.ok(lp.revenue <= 2 * 2 * 6.4 * (pmax * 0.95 - pmin / 0.95) + 1e-6);
  // flat prices: nothing to gain; a lossy battery must not trade
  const flat = dispatchArbitrage({ prices: new Array(24).fill(50), battery, method: 'lp' });
  assert.ok(flat.revenue <= 1e-9);
  // the LP exploits a single spread that the quantile rule also finds; both obey the energy balance
  close(lp.chargedWh * 0.95, lp.dischargedWh / 0.95 + (lp.soc[48] - lp.soc[0]) * 8e6, 1e-9);
});

test('LP arbitrage on a two-price day attains the analytic optimum', () => {
  // 12 h at 20, then 12 h at 100; 8 MWh battery, 6.4 MWh usable, 4 MW rating:
  // buy 6.4/0.95 MWh at 20, sell 6.4·0.95 MWh at 100 → 608 − 134.74 = 473.26
  const prices = Array.from({ length: 24 }, (_, i) => (i < 12 ? 20 : 100));
  const battery = { capacityWh: 8e6, cRate: 0.5, etaCharge: 0.95, etaDischarge: 0.95, socMin: 0.1, socMax: 0.9, soc0: 0.1 };
  const best = 6.4 * 0.95 * 100 - (6.4 / 0.95) * 20;
  const lp = dispatchArbitrage({ prices, battery, method: 'lp' });
  close(lp.revenue, best, 1e-6);
  const greedy = dispatchArbitrage({ prices, battery });
  close(greedy.revenue, best, 1e-6);
});

test('complementarity statistics of hybrid generation', () => {
  const n = 1000;
  const wind = Float64Array.from({ length: n }, (_, i) => 50e6 * (0.5 + 0.5 * Math.sin((2 * Math.PI * i) / 100)));
  const solar = Float64Array.from({ length: n }, (_, i) => 40e6 * (0.5 - 0.5 * Math.sin((2 * Math.PI * i) / 100)));
  const r = complementarity({ series: { wind, solar }, rated: { wind: 50e6, solar: 40e6 }, exportLimit: 50e6 });
  close(r.correlation['wind|solar'], -1, 1e-9);
  close(r.capacityFactors.wind, 0.5, 1e-9);
  close(r.capacityFactors.solar, 0.5, 1e-9);
  close(r.combinedCapacityFactor, 0.5, 1e-9);
  // anti-correlated sources smooth the sum: combined = 45 MW + 5 MW sin → CoV = (5/√2)/45
  close(r.combinedCov, 5 / Math.SQRT2 / 45, 1e-6);
  assert.ok(r.combinedCov < r.cov.wind && r.combinedCov < r.cov.solar);
  // the sum never exceeds 50 MW: no curtailment, cable utilisation 45/50
  close(r.curtailedFraction, 0, 1e-12);
  close(r.cableUtilisation, 0.9, 1e-9);
  // fully correlated sources overload the same cable
  const same = complementarity({ series: { a: wind, b: wind }, rated: { a: 50e6, b: 50e6 }, exportLimit: 50e6 });
  close(same.correlation['a|b'], 1, 1e-9);
  assert.ok(same.curtailedFraction > 0.2);
});
