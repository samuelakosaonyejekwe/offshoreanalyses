// Energy storage and hybrid dispatch: battery model with efficiency, power and state-of-charge limits and a
// simple degradation law; electrolyser hydrogen production; hourly dispatch for firm power, price arbitrage and
// curtailment absorption at an export limit; complementarity statistics of hybrid generation.
//
// Units: power W, energy Wh, time step in hours (dtHours), prices in currency per MWh.

import { linprog } from '../math/optim.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function finitePos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

function normBattery(b) {
  need(b && typeof b === 'object', 'battery must be an object');
  finitePos(b.capacityWh, 'battery.capacityWh');
  const rte = b.roundTripEfficiency ?? 0.88;
  need(rte > 0 && rte <= 1, 'battery.roundTripEfficiency must be within (0, 1]');
  const etaC = b.etaCharge ?? Math.sqrt(rte);
  const etaD = b.etaDischarge ?? (b.etaCharge !== undefined ? rte / b.etaCharge : Math.sqrt(rte));
  need(etaC > 0 && etaC <= 1 && etaD > 0 && etaD <= 1, 'battery charge and discharge efficiencies must be within (0, 1]');
  const socMin = b.socMin ?? 0.1;
  const socMax = b.socMax ?? 0.9;
  need(socMin >= 0 && socMax <= 1 && socMax > socMin, 'battery SoC limits must satisfy 0 <= socMin < socMax <= 1');
  const cRate = b.cRate ?? 0.5;
  finitePos(cRate, 'battery.cRate');
  const pMax = b.powerW ?? cRate * b.capacityWh;
  const soc0 = b.soc0 ?? socMin;
  need(soc0 >= socMin - 1e-12 && soc0 <= socMax + 1e-12, 'battery.soc0 must be within the SoC limits');
  return { cap: b.capacityWh, etaC, etaD, socMin, socMax, pMax, soc0 };
}

/**
 * One time step of the battery model. Sign convention at the AC terminals: positive = discharge (export),
 * negative = charge. Stored energy changes by η_c P_charge Δt when charging and by P_discharge Δt / η_d when
 * discharging; the request is limited by the power rating (C-rate × capacity) and by the SoC window.
 * @param {object} o
 * @param {number} o.request requested terminal power (W), + discharge / − charge
 * @param {number} o.soc state of charge at the start of the step (0–1)
 * @param {number} [o.dtHours=1]
 * @param {object} o.battery { capacityWh, cRate=0.5 (1/h) or powerW, roundTripEfficiency=0.88 or etaCharge /
 *   etaDischarge, socMin=0.1, socMax=0.9 }
 * @returns {{power:number, soc:number, storedDeltaWh:number, lossWh:number}} delivered terminal power (W), new
 *   SoC, change of stored energy (Wh), conversion loss in the step (Wh)
 */
export function batteryStep({ request, soc, dtHours = 1, battery }) {
  finitePos(dtHours, 'dtHours');
  need(Number.isFinite(request), 'request must be finite');
  const b = normBattery({ ...battery, soc0: undefined });
  need(soc >= -1e-9 && soc <= 1 + 1e-9, 'soc must be within [0, 1]');
  return stepCore(b, request, soc, dtHours);
}

function stepCore(b, request, soc, dt) {
  let p = Math.max(-b.pMax, Math.min(b.pMax, request));
  let dE;
  if (p >= 0) {
    const avail = Math.max(0, (soc - b.socMin) * b.cap); // Wh that may leave the store
    p = Math.min(p, (avail * b.etaD) / dt);
    dE = (-p * dt) / b.etaD;
  } else {
    const room = Math.max(0, (b.socMax - soc) * b.cap);
    p = -Math.min(-p, room / (b.etaC * dt));
    dE = -p * dt * b.etaC;
  }
  const loss = p >= 0 ? (p * dt) / b.etaD - p * dt : -p * dt * (1 - b.etaC);
  return { power: p, soc: soc + dE / b.cap, storedDeltaWh: dE, lossWh: loss };
}

/**
 * Run the battery model over a series of requested terminal powers.
 * @param {object} o
 * @param {number[]|Float64Array} o.request requested power per step (W), + discharge / − charge
 * @param {object} o.battery see `batteryStep`; soc0 (default socMin) sets the initial state
 * @param {number} [o.dtHours=1]
 * @returns {{power:Float64Array, soc:Float64Array, chargedWh:number, dischargedWh:number, lossWh:number,
 *   storedDeltaWh:number, equivalentFullCycles:number}} soc has length n + 1; chargedWh / dischargedWh are
 *   terminal energies; cycles = stored-energy throughput (charge side) / capacity
 */
export function runBattery({ request, battery, dtHours = 1 }) {
  need(request && request.length > 0, 'request must be a non-empty series');
  finitePos(dtHours, 'dtHours');
  const b = normBattery(battery);
  const n = request.length;
  const power = new Float64Array(n);
  const soc = new Float64Array(n + 1);
  soc[0] = b.soc0;
  let cin = 0;
  let cout = 0;
  let loss = 0;
  let through = 0;
  for (let i = 0; i < n; i++) {
    const r = stepCore(b, request[i], soc[i], dtHours);
    power[i] = r.power;
    soc[i + 1] = r.soc;
    if (r.power < 0) {
      cin += -r.power * dtHours;
      through += r.storedDeltaWh;
    } else cout += r.power * dtHours;
    loss += r.lossWh;
  }
  return {
    power,
    soc,
    chargedWh: cin,
    dischargedWh: cout,
    lossWh: loss,
    storedDeltaWh: (soc[n] - soc[0]) * b.cap,
    equivalentFullCycles: through / b.cap,
  };
}

/**
 * Simple battery capacity fade: linear calendar ageing plus linear cycle ageing per equivalent full cycle,
 *   C(t)/C0 = 1 − r_cal t − r_cyc N_efc(t).
 * The two rates are technology data supplied by the user (from the cell warranty or test data); end of life is
 * reached at the `endOfLife` capacity fraction. A screening law: no temperature, depth-of-discharge or C-rate
 * dependence.
 * @param {object} o
 * @param {number} o.years service time (years)
 * @param {number} o.cyclesPerYear equivalent full cycles per year (-)
 * @param {number} o.calendarFadePerYear capacity fraction lost per year at rest (-)
 * @param {number} o.cycleFadePerCycle capacity fraction lost per equivalent full cycle (-)
 * @param {number} [o.endOfLife=0.8] end-of-life capacity fraction (-)
 * @returns {{capacityFraction:number, lifeYears:number, lifeCycles:number}}
 */
export function batteryDegradation({ years, cyclesPerYear, calendarFadePerYear, cycleFadePerCycle, endOfLife = 0.8 }) {
  need(years >= 0 && cyclesPerYear >= 0, 'years and cyclesPerYear must be >= 0');
  need(calendarFadePerYear >= 0 && cycleFadePerCycle >= 0, 'fade rates must be >= 0');
  need(endOfLife > 0 && endOfLife < 1, 'endOfLife must be within (0, 1)');
  const rate = calendarFadePerYear + cycleFadePerCycle * cyclesPerYear;
  const lifeYears = rate > 0 ? (1 - endOfLife) / rate : Infinity;
  return { capacityFraction: Math.max(0, 1 - rate * years), lifeYears, lifeCycles: lifeYears * cyclesPerYear };
}

/**
 * Green-hydrogen electrolyser fed by a power series.
 * Hydrogen rate = P / SEC(load), with the specific energy consumption SEC (kWh/kg) at rated load scaled by an
 * optional load-dependent factor; below the minimum load the stack is off; power above rating is not used.
 * Water demand: stoichiometric M_H2O/M_H2 = 18.015/2.016 = 8.94 kg per kg H2, times a margin for purification
 * reject and cooling (argument). Lower heating value of hydrogen 33.33 kWh/kg, higher 39.41 kWh/kg.
 * @param {object} o
 * @param {number[]|Float64Array} o.power available electrical power per step (W)
 * @param {number} o.ratedPower electrolyser rated input power (W)
 * @param {number} o.specificConsumption system specific energy consumption at rated load (kWh/kg), typically
 *   50–58 for alkaline/PEM systems
 * @param {(loadFraction:number)=>number} [o.consumptionFactor] multiplier on the specific consumption as a
 *   function of load fraction (1 at rated load); default constant 1
 * @param {number} [o.minLoad=0.1] minimum load fraction (-)
 * @param {number} [o.waterMargin=1.5] total-to-stoichiometric water ratio (-)
 * @param {number} [o.dtHours=1]
 * @returns {{hydrogenKg:number, waterLitres:number, energyUsedWh:number, energyUnusedWh:number,
 *   capacityFactor:number, meanSpecificConsumption:number, efficiencyLhv:number, onFraction:number,
 *   hydrogenRate:Float64Array}} hydrogenRate in kg per step
 */
export function electrolyser({ power, ratedPower, specificConsumption, consumptionFactor, minLoad = 0.1, waterMargin = 1.5, dtHours = 1 }) {
  need(power && power.length > 0, 'power must be a non-empty series');
  finitePos(ratedPower, 'ratedPower');
  finitePos(specificConsumption, 'specificConsumption');
  need(minLoad >= 0 && minLoad < 1, 'minLoad must be within [0, 1)');
  need(waterMargin >= 1, 'waterMargin must be >= 1');
  const n = power.length;
  const rate = new Float64Array(n);
  let kg = 0;
  let used = 0;
  let unused = 0;
  let on = 0;
  for (let i = 0; i < n; i++) {
    const avail = Math.max(0, power[i]);
    const p = Math.min(avail, ratedPower);
    const load = p / ratedPower;
    if (load >= minLoad && load > 0) {
      const sec = specificConsumption * (consumptionFactor ? consumptionFactor(load) : 1);
      need(sec > 0, 'consumptionFactor must return a positive value');
      rate[i] = ((p * dtHours) / 1000 / sec);
      kg += rate[i];
      used += p * dtHours;
      unused += (avail - p) * dtHours;
      on++;
    } else unused += avail * dtHours;
  }
  return {
    hydrogenKg: kg,
    waterLitres: kg * (18.01528 / 2.01588) * waterMargin,
    energyUsedWh: used,
    energyUnusedWh: unused,
    capacityFactor: used / (ratedPower * n * dtHours),
    meanSpecificConsumption: kg > 0 ? used / 1000 / kg : NaN,
    efficiencyLhv: used > 0 ? (kg * 33.33 * 1000) / used : NaN,
    onFraction: on / n,
    hydrogenRate: rate,
  };
}

/**
 * Dispatch against a firm-power target: surplus generation charges the battery, deficits are met by
 * discharging. Surplus that cannot be stored is exported above the target or curtailed.
 * @param {object} o
 * @param {number[]|Float64Array} o.generation generation per step (W)
 * @param {number} o.target firm power target (W)
 * @param {object} o.battery see `batteryStep`
 * @param {number} [o.dtHours=1]
 * @param {boolean} [o.exportSurplus=false] export unstorable surplus (true) or curtail it (false)
 * @param {number} [o.tolerance=0.01] relative shortfall still counted as firm (-)
 * @returns {{delivered:Float64Array, soc:Float64Array, battery:Float64Array, firmness:number,
 *   firmnessWithoutStorage:number, energyShiftedWh:number, shortfallWh:number, surplusWh:number,
 *   equivalentFullCycles:number, deliveredWh:number}} firmness = fraction of steps meeting the target
 */
export function dispatchFirm({ generation, target, battery, dtHours = 1, exportSurplus = false, tolerance = 0.01 }) {
  need(generation && generation.length > 0, 'generation must be a non-empty series');
  finitePos(target, 'target');
  const b = normBattery(battery);
  const n = generation.length;
  const delivered = new Float64Array(n);
  const bp = new Float64Array(n);
  const soc = new Float64Array(n + 1);
  soc[0] = b.soc0;
  let firm = 0;
  let firm0 = 0;
  let shifted = 0;
  let short = 0;
  let surplus = 0;
  let through = 0;
  let total = 0;
  for (let i = 0; i < n; i++) {
    const g = Math.max(0, generation[i]);
    const r = stepCore(b, target - g, soc[i], dtHours);
    soc[i + 1] = r.soc;
    bp[i] = r.power;
    let out = g + r.power;
    if (out > target) {
      surplus += (out - target) * dtHours;
      if (!exportSurplus) out = target;
    }
    delivered[i] = out;
    total += out * dtHours;
    if (out >= target * (1 - tolerance)) firm++;
    else short += (target - out) * dtHours;
    if (g >= target * (1 - tolerance)) firm0++;
    if (r.power > 0) shifted += r.power * dtHours;
    else through += r.storedDeltaWh;
  }
  return {
    delivered,
    soc,
    battery: bp,
    firmness: firm / n,
    firmnessWithoutStorage: firm0 / n,
    energyShiftedWh: shifted,
    shortfallWh: short,
    surplusWh: surplus,
    equivalentFullCycles: through / b.cap,
    deliveredWh: total,
  };
}

/**
 * Curtailment absorption at an export-cable limit: generation above the limit charges the battery; when there
 * is headroom the battery discharges into it.
 * @param {object} o
 * @param {number[]|Float64Array} o.generation generation per step (W)
 * @param {number} o.exportLimit export limit (W)
 * @param {object} o.battery see `batteryStep`
 * @param {number} [o.dtHours=1]
 * @returns {{exported:Float64Array, soc:Float64Array, curtailedWithoutWh:number, curtailedWithWh:number,
 *   curtailmentAvoidedWh:number, energyShiftedWh:number, exportedWh:number, equivalentFullCycles:number,
 *   cableUtilisation:number}} energyShifted = energy delivered from the battery (Wh)
 */
export function dispatchCurtailment({ generation, exportLimit, battery, dtHours = 1 }) {
  need(generation && generation.length > 0, 'generation must be a non-empty series');
  finitePos(exportLimit, 'exportLimit');
  const b = normBattery(battery);
  const n = generation.length;
  const exported = new Float64Array(n);
  const soc = new Float64Array(n + 1);
  soc[0] = b.soc0;
  let c0 = 0;
  let c1 = 0;
  let shifted = 0;
  let total = 0;
  let through = 0;
  for (let i = 0; i < n; i++) {
    const g = Math.max(0, generation[i]);
    const excess = g - exportLimit;
    if (excess > 0) c0 += excess * dtHours;
    const r = stepCore(b, -excess, soc[i], dtHours); // excess > 0 → charge request; headroom → discharge
    soc[i + 1] = r.soc;
    let out = g + r.power;
    if (out > exportLimit) {
      c1 += (out - exportLimit) * dtHours;
      out = exportLimit;
    }
    exported[i] = out;
    total += out * dtHours;
    if (r.power > 0) shifted += r.power * dtHours;
    else through += r.storedDeltaWh;
  }
  return {
    exported,
    soc,
    curtailedWithoutWh: c0,
    curtailedWithWh: c1,
    curtailmentAvoidedWh: c0 - c1,
    energyShiftedWh: shifted,
    exportedWh: total,
    equivalentFullCycles: through / b.cap,
    cableUtilisation: total / (exportLimit * n * dtHours),
  };
}

function quantile(arr, q) {
  const s = Float64Array.from(arr).sort();
  const pos = q * (s.length - 1);
  const lo = Math.floor(pos);
  const hi = Math.min(s.length - 1, lo + 1);
  return s[lo] + (pos - lo) * (s[hi] - s[lo]);
}

function revenueOf(power, prices, dt) {
  let r = 0;
  for (let i = 0; i < power.length; i++) r += ((power[i] * dt) / 1e6) * prices[i];
  return r;
}

/**
 * Price arbitrage with a grid-connected battery.
 *  - 'threshold': greedy rule — charge at full power when the price is at or below the `lowQuantile` of the
 *    series, discharge when at or above the `highQuantile`.
 *  - 'lp': linear programme maximising Σ price_t (d_t − c_t) Δt subject to the SoC recursion, power and SoC
 *    limits (solved with the shared simplex routine; horizons up to `lpMaxSteps`). The schedule is then replayed
 *    through the battery model. If the LP is unavailable, fails, or does not beat the greedy rule, the greedy
 *    result is returned and `method` reports 'threshold'.
 * No degradation cost and perfect price foresight are assumed.
 * @param {object} o
 * @param {number[]|Float64Array} o.prices price per step (currency/MWh)
 * @param {object} o.battery see `batteryStep`
 * @param {number} [o.dtHours=1]
 * @param {'threshold'|'lp'} [o.method='threshold']
 * @param {number} [o.lowQuantile=0.25]
 * @param {number} [o.highQuantile=0.75]
 * @param {number} [o.lpMaxSteps=96]
 * @returns {{power:Float64Array, soc:Float64Array, revenue:number, chargedWh:number, dischargedWh:number,
 *   equivalentFullCycles:number, method:string}} revenue in currency (positive = profit)
 */
export function dispatchArbitrage({ prices, battery, dtHours = 1, method = 'threshold', lowQuantile = 0.25, highQuantile = 0.75, lpMaxSteps = 96 }) {
  need(prices && prices.length > 0, 'prices must be a non-empty series');
  need(method === 'threshold' || method === 'lp', "method must be 'threshold' or 'lp'");
  need(lowQuantile >= 0 && highQuantile <= 1 && highQuantile > lowQuantile, 'quantiles must satisfy 0 <= low < high <= 1');
  const b = normBattery(battery);
  const n = prices.length;
  const lo = quantile(prices, lowQuantile);
  const hi = quantile(prices, highQuantile);
  const greedyReq = Float64Array.from(prices, (p) => (p <= lo ? -b.pMax : p >= hi ? b.pMax : 0));
  const pack = (run, used) => ({
    power: run.power,
    soc: run.soc,
    revenue: revenueOf(run.power, prices, dtHours),
    chargedWh: run.chargedWh,
    dischargedWh: run.dischargedWh,
    equivalentFullCycles: run.equivalentFullCycles,
    method: used,
  });
  const greedy = pack(runBattery({ request: greedyReq, battery, dtHours }), 'threshold');
  if (method !== 'lp' || n > lpMaxSteps) return greedy;
  try {
    // variables x = [c_0..c_(n−1), d_0..d_(n−1)] (terminal charge and discharge power, W), minimise −revenue
    const c = new Array(2 * n);
    for (let t = 0; t < n; t++) {
      c[t] = (prices[t] * dtHours) / 1e6;
      c[n + t] = (-prices[t] * dtHours) / 1e6;
    }
    const e0 = b.soc0 * b.cap;
    const Aub = [];
    const bub = [];
    for (let t = 0; t < n; t++) {
      const up = new Array(2 * n).fill(0);
      const dn = new Array(2 * n).fill(0);
      for (let s = 0; s <= t; s++) {
        up[s] = b.etaC * dtHours;
        up[n + s] = -dtHours / b.etaD;
        dn[s] = -b.etaC * dtHours;
        dn[n + s] = dtHours / b.etaD;
      }
      Aub.push(up);
      bub.push(b.socMax * b.cap - e0);
      Aub.push(dn);
      bub.push(e0 - b.socMin * b.cap);
    }
    for (let j = 0; j < 2 * n; j++) {
      const row = new Array(2 * n).fill(0);
      row[j] = 1;
      Aub.push(row);
      bub.push(b.pMax);
    }
    const sol = linprog({ c, Aub, bub, bounds: c.map(() => [0, b.pMax]) });
    if (!sol || sol.status !== 'optimal' || !sol.x) return greedy;
    const req = new Float64Array(n);
    for (let t = 0; t < n; t++) req[t] = sol.x[n + t] - sol.x[t];
    const lp = pack(runBattery({ request: req, battery, dtHours }), 'lp');
    return lp.revenue >= greedy.revenue ? lp : greedy;
  } catch {
    return greedy;
  }
}

function meanStd(x) {
  let m = 0;
  for (let i = 0; i < x.length; i++) m += x[i];
  m /= x.length;
  let v = 0;
  for (let i = 0; i < x.length; i++) v += (x[i] - m) ** 2;
  return { mean: m, std: Math.sqrt(v / x.length) };
}

/**
 * Complementarity statistics of hybrid generation (e.g. wind + solar + wave) sharing one export connection:
 * pairwise Pearson correlations, individual and combined capacity factors, variability (coefficient of
 * variation) of the sum against the capacity-weighted individual values, cable utilisation and curtailment at
 * the export limit.
 * @param {object} o
 * @param {Object<string, number[]|Float64Array>} o.series generation per step by source (W), equal lengths
 * @param {Object<string, number>} o.rated rated power by source (W)
 * @param {number} [o.exportLimit] export limit (W); default the sum of ratings
 * @returns {{capacityFactors:Object<string, number>, combinedCapacityFactor:number,
 *   correlation:Object<string, number>, cov:Object<string, number>, combinedCov:number,
 *   cableUtilisation:number, curtailedFraction:number, combined:Float64Array}} correlation keys are 'a|b';
 *   combinedCapacityFactor is relative to the sum of ratings; cableUtilisation = exported energy / (limit × time)
 */
export function complementarity({ series, rated, exportLimit }) {
  const names = Object.keys(series ?? {});
  need(names.length >= 1, 'series must contain at least one source');
  const n = series[names[0]].length;
  need(n > 1, 'series must contain at least two steps');
  let ratedSum = 0;
  const cf = {};
  const cov = {};
  const stats = {};
  const combined = new Float64Array(n);
  for (const k of names) {
    need(series[k].length === n, 'all series must have the same length');
    finitePos(rated?.[k], `rated.${k}`);
    ratedSum += rated[k];
    const s = meanStd(series[k]);
    stats[k] = s;
    cf[k] = s.mean / rated[k];
    cov[k] = s.mean > 0 ? s.std / s.mean : NaN;
    for (let i = 0; i < n; i++) combined[i] += series[k][i];
  }
  const correlation = {};
  for (let a = 0; a < names.length; a++) {
    for (let b = a + 1; b < names.length; b++) {
      const x = series[names[a]];
      const y = series[names[b]];
      let sxy = 0;
      for (let i = 0; i < n; i++) sxy += (x[i] - stats[names[a]].mean) * (y[i] - stats[names[b]].mean);
      const d = n * stats[names[a]].std * stats[names[b]].std;
      correlation[`${names[a]}|${names[b]}`] = d > 0 ? sxy / d : NaN;
    }
  }
  const limit = exportLimit ?? ratedSum;
  finitePos(limit, 'exportLimit');
  let exp = 0;
  let tot = 0;
  for (let i = 0; i < n; i++) {
    tot += combined[i];
    exp += Math.min(combined[i], limit);
  }
  const cs = meanStd(combined);
  return {
    capacityFactors: cf,
    combinedCapacityFactor: cs.mean / ratedSum,
    correlation,
    cov,
    combinedCov: cs.mean > 0 ? cs.std / cs.mean : NaN,
    cableUtilisation: exp / (limit * n),
    curtailedFraction: tot > 0 ? 1 - exp / tot : 0,
    combined,
  };
}
