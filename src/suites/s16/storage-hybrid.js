// Storage and hybrid generation behind one export connection: complementarity of wind, solar and wave output,
// battery dispatch for firm power, curtailment absorption at the export limit and price arbitrage, battery
// capacity fade, and hydrogen production from the energy that would otherwise be curtailed.
// Engineering units (MW, MWh, USD/MWh) are converted to SI (W, Wh) at the top of compute.
import { runBattery, batteryDegradation, electrolyser, dispatchFirm, dispatchCurtailment, dispatchArbitrage, complementarity } from '../../engine/aero/storage.js';
import { makeRng } from '../../engine/math/rng.js';
import { adv, check, worst, sev, fmt, round, run, need, thin, clamp, LIMIT_CODE, CONFIRM, CO2_GRID, energyCo2Text, costInput, money, COST_HELP } from './_common.js';

// standard normal distribution function (Abramowitz & Stegun 26.2.17, |error| < 7.5e-8)
function phi(z) {
  const t = 1 / (1 + 0.2316419 * Math.abs(z));
  const p = 1 - Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI) * t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
  return z >= 0 ? p : 1 - p;
}

/**
 * Synthetic hourly series (fractions of rated power and a price) for a representative period. Wind speed is a
 * first-order autoregressive Gaussian process mapped onto the Weibull distribution and a generic power curve; the
 * sea state lags the wind (waves keep running after the wind drops); solar follows the hour of day with a random
 * daily clearness; the price has a morning and an evening peak and falls when the wind output is high.
 */
function synthetic(inp, n) {
  const rng = makeRng(inp.seed);
  const phi1 = Math.exp(-1 / inp.windTau), wind = new Float64Array(n), solar = new Float64Array(n), wave = new Float64Array(n), price = new Float64Array(n);
  let z = rng.normal(), uLag = inp.weibA * 0.89, clear = 0.7, fMean = 0;
  for (let i = 0; i < n; i++) {
    z = phi1 * z + Math.sqrt(1 - phi1 * phi1) * rng.normal();
    const U = inp.weibA * (-Math.log(Math.max(1 - phi(z), 1e-12))) ** (1 / inp.weibK);
    wind[i] = U < 3 || U > 25 ? 0 : Math.min(1, ((U - 3) / (11.5 - 3)) ** 3);
    // sea state: exponential lag of the wind speed (10 h), fully developed Hs ≈ 0.022 U² and power ∝ Hs² up to the rating
    uLag += (U - uLag) / 10;
    const Hs = 0.022 * uLag * uLag;
    wave[i] = Hs > 7 ? 0 : Math.min(1, (Hs / inp.waveHsRated) ** 2);
    const h = i % 24;
    if (h === 0) clear = 0.25 + 0.75 * Math.sqrt(rng.next());
    solar[i] = h > 6 && h < 18 ? inp.solarPeak * clear * Math.sin(Math.PI * (h - 6) / 12) : 0;
    fMean += wind[i] / n;
  }
  for (let i = 0; i < n; i++) {
    const h = i % 24, shape = 0.6 * Math.cos(2 * Math.PI * (h - 19) / 24) + 0.4 * Math.cos(4 * Math.PI * (h - 8) / 24);
    price[i] = Math.max(0, inp.priceMean * (1 + inp.priceSwing * shape - inp.priceWind * (wind[i] - fMean) + 0.05 * rng.normal()));
  }
  return { wind, solar, wave, price };
}

export default {
  id: 's16-storage-hybrid', suite: 16, title: 'Storage and hybrid dispatch (battery, electrolyser, wind + solar + wave)',
  summary: 'Complementarity of wind, solar and wave generation on one export cable and what a battery and an electrolyser add: firm power, curtailment absorbed, arbitrage revenue, battery life and hydrogen output.',
  covers: [13], fidelity: 'reduced-order',
  standards: ['IEC 61400-26-1:2019 (availability and energy terms)', 'IEC 62933-2-1:2017 (electrical energy storage systems: unit parameters and testing)', 'IEC TS 62933-3-1:2018 (planning and performance assessment of storage systems)', 'ISO 22734:2019 (hydrogen generators using water electrolysis)', 'IEA Wind TCP Task 50 (hybrid power plants) terminology'],
  models: ['Battery energy balance with charge / discharge efficiency, power and state-of-charge limits', 'Linear calendar and cycle capacity fade', 'Electrolyser specific energy consumption with minimum load', 'Firm-power dispatch', 'Curtailment absorption at an export limit', 'Price arbitrage (quantile threshold rule or linear programme)', 'Pearson correlation and variability of hybrid generation', 'Autoregressive Weibull wind series with lagged sea state and diurnal solar (synthetic input)'],
  limits: `Hourly energy balance with perfect foresight of prices in the arbitrage case and no forecast error: no converter dynamics, grid-code response, state-of-charge-dependent efficiency or thermal behaviour, and battery fade is a linear screening law without temperature, depth-of-discharge or C-rate dependence. The built-in generation and price series are synthetic (statistically plausible, not a hindcast): import measured or modelled hourly series for decisions. Each service is dispatched on its own; stacking of services on one battery is not optimised. The time step is the hourly resolution of the data, so no discretisation study applies. ${LIMIT_CODE}`,
  publishes: ['hybrid.combinedCapacityFactor', 'hybrid.curtailmentAvoided', 'hybrid.hydrogenPerYear', 'hybrid.batteryLife'],
  searchTerms: 'hybrid power plant battery storage firming curtailment arbitrage electrolyser offshore wind solar wave complementarity',
  inputs: [
    { key: 'windMW', label: 'Wind capacity', unit: 'MW', type: 'number', default: 400, min: 0, max: 5000, group: 'Generation', role: 'geometry', bind: 'energy.ratedPower', design: true },
    { key: 'solarMW', label: 'Floating solar capacity (AC)', unit: 'MW', type: 'number', default: 100, min: 0, max: 5000, group: 'Generation', role: 'geometry', design: true },
    { key: 'waveMW', label: 'Wave capacity', unit: 'MW', type: 'number', default: 20, min: 0, max: 1000, group: 'Generation', role: 'geometry', design: true },
    { key: 'exportMW', label: 'Export connection limit', unit: 'MW', type: 'number', default: 420, min: 1, max: 10000, group: 'Generation', role: 'bc', design: true, help: 'Boundary condition at the point of connection: firm export capacity of the cable and grid agreement. Below the sum of the ratings the site is overplanted.' },
    { key: 'series', label: 'Hourly series (optional import)', type: 'table', group: 'Generation', role: 'bc', default: [],
      columns: [{ key: 'wind', label: 'Wind', unit: 'MW', type: 'number' }, { key: 'solar', label: 'Solar', unit: 'MW', type: 'number' }, { key: 'wave', label: 'Wave', unit: 'MW', type: 'number' }, { key: 'price', label: 'Price', unit: 'USD/MWh', type: 'number' }],
      help: 'One row per hour, at least 24 rows; results are scaled to a year. Leave empty to use the synthetic period below.' },
    { key: 'days', label: 'Synthetic period', unit: 'days', type: 'int', default: 56, min: 7, max: 366, group: 'Synthetic period', role: 'control', help: 'Length of the representative hourly period; results are scaled to a year.' },
    { key: 'weibA', label: 'Weibull scale A at hub height', unit: 'm/s', type: 'number', default: 10.6, min: 4, max: 16, group: 'Synthetic period', role: 'environment', bind: 'site.windWeibullA' },
    { key: 'weibK', label: 'Weibull shape k', unit: '-', type: 'number', default: 2.2, min: 1.2, max: 3.5, group: 'Synthetic period', role: 'environment', bind: 'site.windWeibullK' },
    { key: 'windTau', label: 'Wind persistence time', unit: 'h', type: 'number', default: 30, min: 3, max: 200, group: 'Synthetic period', role: 'environment', help: 'Decorrelation time of the hourly wind speed: 20–40 h for synoptic weather.' },
    { key: 'solarPeak', label: 'Clear-day solar peak / rated', unit: '-', type: 'number', default: 0.85, min: 0.1, max: 1, group: 'Synthetic period', role: 'environment', help: 'Noon output on a clear day of the season represented, as a fraction of the AC rating.' },
    { key: 'waveHsRated', label: 'Significant wave height at rated wave power', unit: 'm', type: 'number', default: 3.5, min: 1, max: 8, group: 'Synthetic period', role: 'environment' },
    { key: 'priceMean', label: 'Mean power price', unit: 'USD/MWh', type: 'number', default: 80, min: 1, max: 1000, group: 'Synthetic period', role: 'cost', help: `Level of the synthetic hourly price. ${COST_HELP}` },
    { key: 'priceSwing', label: 'Daily price swing (fraction of the mean)', unit: '-', type: 'number', default: 0.35, min: 0, max: 1, group: 'Synthetic period', role: 'bc' },
    { key: 'priceWind', label: 'Price sensitivity to wind output', unit: '-', type: 'number', default: 0.4, min: 0, max: 1.5, group: 'Synthetic period', role: 'bc', help: 'Fractional fall of the price for wind output one rated unit above its mean (merit-order effect).' },
    { key: 'seed', label: 'Random seed', type: 'int', default: 11, min: 1, max: 99999, group: 'Synthetic period', role: 'control' },
    { key: 'capMWh', label: 'Battery energy capacity', unit: 'MWh', type: 'number', default: 400, min: 0.1, max: 50000, group: 'Battery', role: 'geometry', design: true },
    { key: 'cRate', label: 'Power rating / energy capacity', unit: '1/h', type: 'number', default: 0.5, min: 0.05, max: 4, group: 'Battery', role: 'geometry', design: true, help: '0.5 = a two-hour battery, 0.25 = four hours.' },
    { key: 'rte', label: 'Round-trip efficiency (AC to AC)', unit: '-', type: 'number', default: 0.88, min: 0.5, max: 1, group: 'Battery', role: 'operation', uncertain: { dist: 'normal', cov: 0.02 } },
    { key: 'socMin', label: 'Lowest state of charge', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.5, group: 'Battery', role: 'bc' },
    { key: 'socMax', label: 'Highest state of charge', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1, group: 'Battery', role: 'bc' },
    { key: 'soc0', label: 'Initial state of charge', unit: '-', type: 'number', default: 0.1, min: 0, max: 1, group: 'Battery', role: 'ic', help: 'Initial condition of every dispatch; must lie within the limits.' },
    { key: 'service', label: 'Service the battery is committed to', type: 'select', default: 'curtail', group: 'Dispatch', role: 'control', options: [{ value: 'curtail', label: 'Curtailment absorption at the export limit' }, { value: 'firm', label: 'Firm power' }, { value: 'arbitrage', label: 'Price arbitrage' }],
      help: 'All three services are evaluated; the cycle count and the battery life are taken from this one.' },
    { key: 'firmMW', label: 'Firm power target', unit: 'MW', type: 'number', default: 60, min: 0.1, max: 10000, group: 'Dispatch', role: 'bc', design: true, help: 'Constant delivery the plant commits to; surplus above it is exported as non-firm energy.' },
    { key: 'arbMethod', label: 'Arbitrage schedule', type: 'select', default: 'threshold', group: 'Dispatch', role: 'control', options: [{ value: 'threshold', label: 'Quantile rule (charge in the cheapest quarter, discharge in the dearest)' }, { value: 'lp', label: 'Linear programme on the first four days (slower)' }] },
    { key: 'calFade', label: 'Calendar fade', unit: '%/yr', type: 'number', default: 1.0, min: 0, max: 10, group: 'Battery ageing', role: 'material', help: `Capacity lost per year at rest, from the cell warranty. ${CONFIRM}` },
    { key: 'cycFade', label: 'Cycle fade per 1000 equivalent full cycles', unit: '%', type: 'number', default: 2.5, min: 0, max: 40, group: 'Battery ageing', role: 'material', help: `From the cell warranty or test data. ${CONFIRM}` },
    { key: 'eol', label: 'End-of-life capacity', unit: '%', type: 'number', default: 80, min: 50, max: 95, group: 'Battery ageing', role: 'criteria' },
    { key: 'elMW', label: 'Electrolyser rated input', unit: 'MW', type: 'number', default: 20, min: 0, max: 5000, group: 'Electrolyser', role: 'geometry', design: true, help: 'Fed with the energy still curtailed after the battery. 0 = no electrolyser.' },
    { key: 'sec', label: 'Specific energy consumption', unit: 'kWh/kg', type: 'number', default: 54, min: 40, max: 80, group: 'Electrolyser', role: 'operation', uncertain: { dist: 'normal', cov: 0.04 }, help: 'System value at rated load: 50–58 for alkaline and PEM systems (the lower heating value of hydrogen is 33.3 kWh/kg).' },
    { key: 'elMin', label: 'Minimum load', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.6, group: 'Electrolyser', role: 'bc', help: 'Below this fraction of the rating the stack is off: about 0.05–0.1 PEM, 0.2–0.4 alkaline.' },
    { key: 'waterMargin', label: 'Water demand / stoichiometric', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Electrolyser', role: 'operation' },
    { key: 'firmReq', label: 'Required firmness', unit: '%', type: 'number', default: 80, min: 1, max: 100, group: 'Acceptance', role: 'criteria', help: 'Share of hours in which the firm target must be met (contract or grid requirement).' },
    { key: 'curtMax', label: 'Allowable curtailment after storage and electrolysis', unit: '%', type: 'number', default: 3, min: 0.1, max: 50, group: 'Acceptance', role: 'criteria' },
    { key: 'lifeReq', label: 'Required battery life', unit: 'yr', type: 'number', default: 10, min: 1, max: 30, group: 'Acceptance', role: 'criteria', help: 'Warranty period or planned replacement interval.' },
    costInput('battCost', 'Battery installed cost', 'USD/kWh', 300, 3000, 'Per kWh of energy capacity, installed offshore or at the landfall.'),
    costInput('h2Price', 'Value of hydrogen at the plant gate', 'USD/kg', 5, 50),
  ],
  compute(inp) {
    const MW = 1e6;
    need(inp.socMax > inp.socMin + 0.05, 'The highest state of charge must be above the lowest by at least 0.05.');
    need(inp.soc0 >= inp.socMin - 1e-9 && inp.soc0 <= inp.socMax + 1e-9, `The initial state of charge (${inp.soc0}) must lie between the limits ${inp.socMin} and ${inp.socMax}.`);
    need(inp.windMW + inp.solarMW + inp.waveMW > 0, 'Give at least one generation source a capacity above zero.');
    // series: imported hourly table, or the synthetic period
    const rows = (inp.series || []).filter((r) => r && [r.wind, r.solar, r.wave].some((v) => Number(v) > 0));
    const imported = rows.length > 0;
    let W, So, Wa, price, n;
    if (imported) {
      need(rows.length >= 24, `The imported series has ${rows.length} hourly rows: at least 24 are needed.`);
      need(rows.length <= 8784, 'The imported series is longer than one year of hours (8784 rows).');
      n = rows.length;
      const col = (k) => Float64Array.from(rows, (r) => { const v = Number(r[k]) || 0; need(v >= 0, `Imported series: ${k} must not be negative.`); return v * MW; });
      W = col('wind'); So = col('solar'); Wa = col('wave'); price = Float64Array.from(rows, (r) => Math.max(0, Number(r.price) || 0));
    } else {
      n = Math.round(inp.days) * 24;
      const s = synthetic(inp, n);
      W = s.wind.map((f) => f * inp.windMW * MW); So = s.solar.map((f) => f * inp.solarMW * MW); Wa = s.wave.map((f) => f * inp.waveMW * MW); price = s.price;
    }
    const ann = 8760 / n, limit = inp.exportMW * MW;
    const src = [['wind', W, inp.windMW, 'Wind'], ['solar', So, inp.solarMW, 'Solar'], ['wave', Wa, inp.waveMW, 'Wave']].filter(([, s, r]) => r > 0 || s.some((v) => v > 0));
    need(src.length > 0, 'No generation: every series is zero.');
    const series = {}, rated = {};
    for (const [k, s, r] of src) { series[k] = s; rated[k] = Math.max(r * MW, Math.max(...s), 1); }
    const cmp = run(() => complementarity({ series, rated, exportLimit: limit }), 'Complementarity');
    const gen = cmp.combined, genWh = gen.reduce((s, v) => s + v, 0);
    need(genWh > 0, 'The generation series sums to zero.');
    const corr = (a, b) => { const v = cmp.correlation[`${a}|${b}`] ?? cmp.correlation[`${b}|${a}`]; return Number.isFinite(v) ? v : 0; };
    const ratedSum = Object.values(rated).reduce((s, v) => s + v, 0);
    const covW = src.reduce((s, [k]) => s + (Number.isFinite(cmp.cov[k]) ? cmp.cov[k] : 0) * rated[k], 0) / ratedSum;

    const battery = { capacityWh: inp.capMWh * MW, cRate: inp.cRate, roundTripEfficiency: inp.rte, socMin: inp.socMin, socMax: inp.socMax, soc0: inp.soc0 };
    const etaC = Math.sqrt(inp.rte);
    const cur = run(() => dispatchCurtailment({ generation: gen, exportLimit: limit, battery }), 'Curtailment dispatch');
    const firm = run(() => dispatchFirm({ generation: gen, target: inp.firmMW * MW, battery, exportSurplus: true }), 'Firm-power dispatch');
    const meanPrice = price.reduce((s, v) => s + v, 0) / n;
    const hasPrice = meanPrice > 0;
    let arb = { revenue: 0, equivalentFullCycles: 0, method: 'no price series', dischargedWh: 0, chargedWh: 0, soc: new Float64Array(n + 1).fill(inp.soc0), power: new Float64Array(n) };
    if (hasPrice) {
      arb = run(() => dispatchArbitrage({ prices: price, battery, method: 'threshold' }), 'Arbitrage dispatch');
      if (inp.arbMethod === 'lp') {
        // the linear programme covers the first four days; its gain over the rule on the same hours scales the whole period
        const m = Math.min(96, n), p4 = price.slice(0, m);
        const lp = run(() => dispatchArbitrage({ prices: p4, battery, method: 'lp', lpMaxSteps: 96 }), 'Arbitrage linear programme'), th = dispatchArbitrage({ prices: p4, battery, method: 'threshold' });
        if (lp.method === 'lp' && th.revenue > 0) arb = { ...arb, revenue: arb.revenue * lp.revenue / th.revenue, equivalentFullCycles: arb.equivalentFullCycles * lp.equivalentFullCycles / Math.max(th.equivalentFullCycles, 1e-9), method: 'lp (first four days, scaled)' };
      }
    }
    // energy still curtailed after the battery, hour by hour: generation − export − battery charging
    const resid = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const charge = Math.max(0, cur.soc[i + 1] - cur.soc[i]) * battery.capacityWh / etaC;
      resid[i] = gen[i] > limit ? Math.max(0, gen[i] - cur.exported[i] - charge) : 0;
    }
    const el = inp.elMW > 0 ? run(() => electrolyser({ power: resid, ratedPower: inp.elMW * MW, specificConsumption: inp.sec, minLoad: inp.elMin, waterMargin: inp.waterMargin }), 'Electrolyser') : { hydrogenKg: 0, waterLitres: 0, energyUsedWh: 0, capacityFactor: 0, onFraction: 0, efficiencyLhv: 0 };
    const curtFinalWh = Math.max(0, cur.curtailedWithWh - el.energyUsedWh);
    const svc = inp.service === 'firm' ? firm : inp.service === 'arbitrage' ? arb : cur;
    const cyclesYr = svc.equivalentFullCycles * ann;
    const deg = batteryDegradation({ years: inp.lifeReq, cyclesPerYear: cyclesYr, calendarFadePerYear: inp.calFade / 100, cycleFadePerCycle: inp.cycFade / 100 / 1000, endOfLife: inp.eol / 100 });
    const life = Math.min(deg.lifeYears, 99);
    const h2t = el.hydrogenKg * ann / 1e3, avoidedMWh = cur.curtailmentAvoidedWh * ann / MW, arbYr = arb.revenue * ann;
    const curt0 = cur.curtailedWithoutWh / genWh * 100, curt1 = curtFinalWh / genWh * 100;

    const checks = [
      check('Curtailment after storage and electrolysis', curt1, inp.curtMax, '%', 'Project energy-yield criterion (IEC 61400-26-1 energy terms: curtailed energy as lost production)',
        'Share of the generated energy that still cannot be exported, stored or converted. Above the allowance the plant is overplanted for its connection: energy and revenue are thrown away every windy, sunny hour.'),
      { name: 'Firmness of the firm-power target', demand: round(inp.firmReq, 5), capacity: round(firm.firmness * 100, 5), unit: '%', uc: round(inp.firmReq / Math.max(firm.firmness * 100, 1e-9), 5), ref: 'IEC TS 62933-3-1:2018 §5 (firming application); contract requirement',
        explain: 'Required share of hours at the firm target against the share achieved with the battery. Below the requirement the plant cannot offer the firm product without more storage, a lower target or back-up capacity.' },
      { name: 'Battery life against the required life', demand: round(inp.lifeReq, 5), capacity: round(life, 5), unit: 'yr', uc: round(inp.lifeReq / Math.max(life, 1e-9), 5), ref: 'IEC 62933-2-1:2017 (capacity and end-of-life definition); cell warranty fade rates',
        explain: 'Years until the capacity fades to the end-of-life value at the cycling of the committed service. Shorter than required means an earlier battery replacement or augmentation than the business case assumes.' },
    ];
    const gov = worst(checks);
    const hrs = Array.from({ length: n }, (_, i) => i), nPlot = 600;
    const [th, tGen, tExp, tSoc, tPrice] = thin([hrs, gen, cur.exported, cur.soc.slice(0, n), price], nPlot);
    const [, tDel, tSocF] = thin([hrs, firm.delivered, firm.soc.slice(0, n)], nPlot);
    const sorted = (a) => Array.from(a).sort((x, y) => y - x), dur = thin([Array.from({ length: n }, (_, i) => 100 * (i + 0.5) / n), sorted(gen), sorted(cur.exported), sorted(W)], 300);
    const battCapex = inp.capMWh * 1e3 * inp.battCost, avoidedRev = avoidedMWh * (hasPrice ? meanPrice : inp.priceMean), h2Rev = h2t * 1e3 * inp.h2Price;
    const best = [['curtailment absorption', avoidedRev], ['price arbitrage', arbYr]].sort((p, q) => q[1] - p[1])[0];
    const labels = src.map(([, , , l]) => l);

    return {
      kpis: [
        { key: 'cfCombined', label: 'Combined capacity factor', value: round(cmp.combinedCapacityFactor * 100, 4), unit: '%', explain: `Mean output of all sources over the sum of their ratings (${labels.map((l, i) => `${l} ${fmt(cmp.capacityFactors[src[i][0]] * 100)} %`).join(', ')}) before the export limit.` },
        { key: 'corrWindSolar', label: 'Correlation wind – solar', value: round(corr('wind', 'solar'), 3), unit: '-', explain: 'Pearson correlation of the hourly outputs. Near zero or negative means the two fill each other\'s gaps and share the cable well; near +1 means they peak together and compete for it.' },
        { key: 'corrWindWave', label: 'Correlation wind – wave', value: round(corr('wind', 'wave'), 3), unit: '-', explain: 'Waves follow the wind with a lag of hours, so the correlation is positive but below 1: wave output carries on after the wind has dropped.' },
        { key: 'covCombined', label: 'Variability of the combined output', value: round(cmp.combinedCov, 4), unit: '-', explain: `Standard deviation over mean of the summed output. The capacity-weighted value of the sources taken alone is ${fmt(covW)}: the difference is the smoothing gained by combining them.` },
        { key: 'cableUtil', label: 'Export-cable utilisation with the battery', value: round(cur.cableUtilisation * 100, 4), unit: '%', explain: `Exported energy over the export limit running all the time (${fmt(cmp.cableUtilisation * 100)} % without storage). The cable is the fixed asset the hybrid is meant to use better.` },
        { key: 'curtNoStore', label: 'Curtailment without storage', value: round(curt0, 4), unit: '%', explain: 'Share of generated energy above the export limit with no battery or electrolyser.' },
        { key: 'curtFinal', label: 'Curtailment after storage and electrolysis', value: round(curt1, 4), unit: '%', limit: inp.curtMax, sense: 'max', explain: 'Share of generated energy that is still lost after the battery has absorbed what it can and the electrolyser has taken the rest up to its rating.' },
        { key: 'curtAvoided', label: 'Curtailment avoided by the battery', value: round(avoidedMWh, 5), unit: 'MWh/yr', explain: 'Energy that the battery takes above the export limit instead of it being curtailed, scaled to a year; it is re-exported, less the round-trip loss, when the cable has headroom.' },
        { key: 'firmness', label: 'Firmness with the battery', value: round(firm.firmness * 100, 4), unit: '%', limit: inp.firmReq, sense: 'min', explain: `Share of hours in which ${fmt(inp.firmMW)} MW is delivered with the battery covering the gaps.` },
        { key: 'firmNoStore', label: 'Firmness without storage', value: round(firm.firmnessWithoutStorage * 100, 4), unit: '%', explain: 'Share of hours in which generation alone meets the firm target: the baseline the battery improves on.' },
        { key: 'arbRevenue', label: 'Arbitrage margin (stand-alone battery)', value: round(arbYr / 1e3, 4), unit: 'k USD/yr', explain: `Revenue of discharging at high prices less the cost of charging at low prices (${arb.method}), with perfect foresight and no degradation cost: an upper bound.` },
        { key: 'cycles', label: 'Equivalent full cycles in the committed service', value: round(cyclesYr, 4), unit: '1/yr', explain: 'Stored-energy throughput over the capacity, per year, for the service selected. It drives the cycle part of the capacity fade.' },
        { key: 'battLife', label: 'Battery life to end-of-life capacity', value: round(life, 4), unit: 'yr', limit: inp.lifeReq, sense: 'min', explain: `Years until the capacity falls to ${inp.eol} % with the calendar and cycle fade rates entered (99 = no fade entered). After the required ${fmt(inp.lifeReq)} years the capacity is ${fmt(deg.capacityFraction * 100)} %.` },
        { key: 'h2', label: 'Hydrogen from curtailed energy', value: round(h2t, 5), unit: 't/yr', explain: `Hydrogen made from the energy still curtailed after the battery, at ${inp.sec} kWh/kg; needs about ${fmt(el.waterLitres * ann / 1e6)} thousand m³ of purified water a year.` },
        { key: 'elCF', label: 'Electrolyser capacity factor', value: round(el.capacityFactor * 100, 4), unit: '%', explain: 'Energy used over the rating running all the time. An electrolyser fed only by curtailed energy runs few hours: low values mean the stack cost is spread over little hydrogen.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Generation, export and state of charge (curtailment service)', x: { label: 'Time', unit: 'h' }, y: { label: 'Power (MW) / state of charge (%)', unit: '' },
          series: [{ name: 'Generation (MW)', x: th, y: tGen.map((v) => round(v / MW)) }, { name: 'Export with battery (MW)', x: th, y: tExp.map((v) => round(v / MW)) }, { name: 'Export limit (MW)', x: [0, n - 1], y: [inp.exportMW, inp.exportMW], dash: true }, { name: 'State of charge (%)', x: th, y: tSoc.map((v) => round(v * 100)) }],
          explain: 'Hours above the dashed limit are curtailed unless the battery (rising state of charge) takes the excess; it empties again into the next lull. A battery pinned at its upper limit during long windy spells is too small to matter.' },
        { type: 'line', title: 'Firm-power dispatch', x: { label: 'Time', unit: 'h' }, y: { label: 'Power (MW) / state of charge (%)', unit: '' },
          series: [{ name: 'Delivered (MW)', x: th, y: tDel.map((v) => round(Math.min(v, 3 * inp.firmMW * MW) / MW)) }, { name: 'Firm target (MW)', x: [0, n - 1], y: [inp.firmMW, inp.firmMW], dash: true }, { name: 'State of charge (%)', x: th, y: tSocF.map((v) => round(v * 100)) }],
          explain: 'Delivery against the firm target (surplus is clipped at three times the target for legibility). Dips below the dashed line are the hours the commitment is missed: they cluster in multi-day lulls that a few hours of storage cannot bridge.' },
        { type: 'line', title: 'Duration curves', x: { label: 'Share of time exceeded', unit: '%' }, y: { label: 'Power', unit: 'MW' },
          series: [{ name: 'Combined generation', x: dur[0].map((v) => round(v)), y: dur[1].map((v) => round(v / MW)) }, { name: 'Export with battery', x: dur[0].map((v) => round(v)), y: dur[2].map((v) => round(v / MW)) }, { name: 'Wind alone', x: dur[0].map((v) => round(v)), y: dur[3].map((v) => round(v / MW)) }],
          explain: 'Output sorted from highest to lowest. The area between the generation and export curves at the left is curtailment; the lift of the combined curve over wind alone at the right is what solar and wave add in calm hours.' },
        { type: 'line', title: 'Price series', x: { label: 'Time', unit: 'h' }, y: { label: 'Price', unit: 'USD/MWh' }, series: [{ name: 'Price', x: th, y: tPrice.map((v) => round(v)) }],
          explain: 'Hourly price used for the arbitrage case. The spread between the cheap and dear hours of a day, not the level, sets what a battery can earn.' },
      ],
      tables: [
        { title: 'Sources', columns: [{ key: 's', label: 'Source' }, { key: 'r', label: 'Rated', unit: 'MW' }, { key: 'cf', label: 'Capacity factor', unit: '%' }, { key: 'cov', label: 'Variability (σ/mean)' }, { key: 'e', label: 'Energy', unit: 'GWh/yr' }],
          rows: src.map(([k, s, , l]) => ({ s: l, r: round(rated[k] / MW), cf: round(cmp.capacityFactors[k] * 100), cov: round(Number.isFinite(cmp.cov[k]) ? cmp.cov[k] : 0), e: round(s.reduce((q, v) => q + v, 0) * ann / 1e9) })),
          explain: 'One row per source: rating, capacity factor, hour-to-hour variability and annual energy before the export limit.' },
        { title: 'Battery services compared', columns: [{ key: 's', label: 'Service' }, { key: 'e', label: 'Energy shifted', unit: 'MWh/yr' }, { key: 'c', label: 'Equivalent full cycles', unit: '1/yr' }, { key: 'v', label: 'Result' }],
          rows: [{ s: 'Curtailment absorption', e: round(cur.energyShiftedWh * ann / MW), c: round(cur.equivalentFullCycles * ann), v: `${fmt(avoidedMWh)} MWh/yr of curtailment avoided` }, { s: 'Firm power', e: round(firm.energyShiftedWh * ann / MW), c: round(firm.equivalentFullCycles * ann), v: `firmness ${fmt(firm.firmness * 100)} % (${fmt(firm.firmnessWithoutStorage * 100)} % without)` }, { s: 'Price arbitrage', e: round(arb.dischargedWh * ann / MW), c: round(arb.equivalentFullCycles * ann), v: `${money(arbYr)} a year` }],
          explain: 'The same battery dispatched for each service on its own. Compare the energy moved and the cycling each service costs with what it earns; the committed service sets the battery life.' },
      ],
      notes: [
        `${imported ? `Imported series of ${n} hours` : `Synthetic period of ${Math.round(inp.days)} days (seed ${inp.seed})`}, scaled to a year by ${fmt(ann)}. The plant is ${ratedSum > limit ? `overplanted by ${fmt((ratedSum / limit - 1) * 100)} %` : 'not overplanted'} against the ${fmt(inp.exportMW)} MW connection.`,
        `Battery ${fmt(inp.capMWh)} MWh / ${fmt(inp.capMWh * inp.cRate)} MW committed to ${inp.service === 'curtail' ? 'curtailment absorption' : inp.service === 'firm' ? 'firm power' : 'price arbitrage'}: ${fmt(cyclesYr)} equivalent full cycles a year, life ${fmt(life)} years.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit} (utilisation ${fmt(gov.uc)}). Curtailment falls from ${fmt(curt0)} % to ${fmt(curt1)} %; firmness of ${fmt(inp.firmMW)} MW rises from ${fmt(firm.firmnessWithoutStorage * 100)} % to ${fmt(firm.firmness * 100)} %.`,
          gov.uc > 1 ? (gov === checks[0] ? `Raise the electrolyser rating (the peak residual surplus is ${fmt(Math.max(...resid) / MW)} MW), enlarge the battery, or reduce the overplanting of the ${fmt(inp.exportMW)} MW connection.` : gov === checks[1] ? `A few hours of storage cannot bridge multi-day lulls: lower the firm target (generation alone meets ${fmt(inp.firmMW)} MW in ${fmt(firm.firmnessWithoutStorage * 100)} % of hours), add long-duration storage or contract back-up capacity.` : `Choose a cell with lower fade, limit the cycling (${fmt(cyclesYr)} a year now) or plan capacity augmentation in year ${fmt(Math.floor(life))}.`)
            : 'The storage and conversion plant meets the criteria entered. Repeat with at least one measured or hindcast year of hourly data and with forecast error before fixing the ratings, and check the grid-code response separately.',
          'IEC TS 62933-3-1:2018 (planning of electrical energy storage systems) and IEC 62933-2-1:2017 (unit parameters); IEA Wind TCP Task 50 hybrid-plant practice.',
          `Avoided curtailment and hydrogen recover ${fmt(avoidedMWh + el.energyUsedWh * ann / MW)} MWh a year of otherwise wasted renewable energy ≈ ${fmt((avoidedMWh * inp.rte + el.energyUsedWh * ann / MW * el.efficiencyLhv) * CO2_GRID)} t CO₂ a year displaced after conversion losses (screening factor ${CO2_GRID} t CO₂/MWh). Battery cells carry embodied emissions of the order of 60–100 kg CO₂e per kWh of capacity (${fmt(inp.capMWh * 0.08)} kt here at 80 kg/kWh): a second-life or recycling route should be fixed at procurement.`),
        adv('economic', 'info', `At the indicative rates the battery costs ${money(battCapex)}; curtailment absorption is worth ${money(avoidedRev)} a year, arbitrage up to ${money(arbYr)} a year and the hydrogen ${money(h2Rev)} a year (${fmt(h2t)} t).`,
          `${best[1] > 0 ? `The better single service is ${best[0]}: simple payback ${fmt(battCapex / best[1])} years against a battery life of ${fmt(life)} years${battCapex / best[1] > life ? ' — the battery does not pay back on one service alone' : ''}.` : 'Neither service earns on these inputs.'} Stack services (firming payments, grid services, avoided cable capacity) before deciding, and size the electrolyser on its capacity factor (${fmt(el.capacityFactor * 100)} % on curtailed energy alone).`,
          'Whole-life costing per ISO 15686-5 with battery replacement at end of life; LCOE / LCOS / LCOH conventions of IEA and IRENA; value stacking as in IEC TS 62933-3-1 application planning.',
          energyCo2Text(avoidedMWh)),
      ],
      publish: { 'hybrid.combinedCapacityFactor': cmp.combinedCapacityFactor, 'hybrid.curtailmentAvoided': cur.curtailmentAvoidedWh * ann, 'hybrid.hydrogenPerYear': h2t * 1e3, 'hybrid.batteryLife': life },
    };
  },
  verify: [
    { // 24 h of 60 MW against a 50 MW limit: 10 MW × 24 h = 240 MWh curtailed without storage = 16.667 % of 1440 MWh.
      // Loss-free 40 MWh battery from 10 % to 90 % state of charge stores 32 MWh and never finds headroom to discharge:
      // avoided = 32 MWh per day × 365 = 11 680 MWh/yr. The remaining 208 MWh feed the 20 MW electrolyser (always above its
      // 2 MW minimum load): 208 000 kWh / 52 kWh/kg = 4000 kg per day → 1460 t/yr, and nothing is left curtailed.
      name: 'Constant surplus above the export limit: battery fills once, electrolyser takes the rest', source: 'Hand calculation (energy balance)',
      inputs: { series: Array.from({ length: 24 }, () => ({ wind: 60, solar: 0, wave: 0, price: 50 })), windMW: 60, solarMW: 0, waveMW: 0, exportMW: 50, capMWh: 40, cRate: 0.5, rte: 1, socMin: 0.1, socMax: 0.9, soc0: 0.1, elMW: 20, sec: 52, elMin: 0.1, firmMW: 50 },
      expect: [{ kpi: 'curtNoStore', value: 16.6667, tol: 3e-4 }, { kpi: 'curtAvoided', value: 11680, tol: 1e-6 }, { kpi: 'h2', value: 1460, tol: 1e-6 }, { kpi: 'curtFinal', value: 0, abs: 1e-9 }] },
    { // 12 h at 80 MW then 12 h at 20 MW against a 50 MW firm target: generation alone is firm half the time. A loss-free
      // 450 MWh battery (usable 0.8 × 450 = 360 MWh, 225 MW) stores 30 MW × 12 h = 360 MWh and returns it: firm all 24 h.
      name: 'Firm power from a two-level generation day with a loss-free battery', source: 'Hand calculation (energy balance)',
      inputs: { series: Array.from({ length: 24 }, (_, i) => ({ wind: i < 12 ? 80 : 20, solar: 0, wave: 0, price: 50 })), windMW: 80, solarMW: 0, waveMW: 0, exportMW: 100, capMWh: 450, cRate: 0.5, rte: 1, socMin: 0.1, socMax: 0.9, soc0: 0.1, firmMW: 50, elMW: 0, service: 'firm' },
      expect: [{ kpi: 'firmNoStore', value: 50, tol: 1e-9 }, { kpi: 'firmness', value: 100, tol: 1e-9 }, { kpi: 'cycles', value: 0.8 * 365, tol: 1e-6 }] },
  ],
};
