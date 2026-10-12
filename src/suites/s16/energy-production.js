// Energy conversion efficiency and power production: gross annual energy from the power curve and the Weibull wind
// climate, site air-density correction, the gross-to-net loss cascade, P50 / P75 / P90 / P99 exceedance values,
// monthly and daily generation, capacity factor and resource utilisation.
import { aepWeibull, aepBins, airDensity, airDensityCorrection, weibull, weibullScaleFromMean, lossCascade, aepUncertainty, curvePower } from '../../engine/aero/aep.js';
import { curveInputs, getCurve } from './_curve.js';
import { site, adv, check, fmt, round, run, need, linspace, HOURS_YEAR, LIMIT_POLAR, LIMIT_CODE, PXX, CONFIRM, energyCo2Text, CO2_GRID } from './_common.js';

const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
const DAYS = [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
// Typical North-Sea seasonality of the monthly mean wind speed relative to the annual mean (indicative).
const SEASON = [1.18, 1.13, 1.08, 0.98, 0.9, 0.86, 0.83, 0.86, 0.96, 1.05, 1.11, 1.17];
const pct = (label, key, def, max, help, extra = {}) => ({ key, label, unit: '%', type: 'number', default: def, min: 0, max, group: 'Losses', role: 'operation', help, ...extra });
const unc = (label, key, def, help) => ({ key, label, unit: '%', type: 'number', default: def, min: 0, max: 20, group: 'Uncertainty', role: 'criteria', help });

export default {
  id: 's16-energy-production', suite: 16, title: 'Energy yield, losses and P50/P90',
  summary: 'Gross and net annual energy production of the plant from the power curve and the wind climate, with the loss cascade, exceedance values, monthly and daily generation and capacity factor.',
  covers: [13], fidelity: 'analytical',
  standards: ['IEC 61400-12-1:2017 (power curve, method of bins, air-density normalisation)', 'IEC 61400-15-1 (site energy yield assessment, loss and uncertainty categories)', 'MEASNET Evaluation of site-specific wind conditions (2016)', 'IEC TS 61400-26-1 (availability)'],
  models: ['Turbine power-coefficient curve and tip-speed-ratio schedule', 'Weibull wind-speed distribution', 'IEC method of bins and Simpson integration of P(u)·f(u)', 'Air-density correction of a pitch-regulated curve', 'Multiplicative loss cascade', 'Root-sum-square uncertainty with normal exceedance levels', 'Betz limit / resource utilisation'],
  limits: `${LIMIT_POLAR} A single Weibull distribution represents the long-term wind climate (no direction dependence here — wake losses come from the array module or are entered); monthly energies use a seasonal scaling of the Weibull scale at constant shape. ${PXX} The uncertainty model assumes independent, normally distributed components. ${LIMIT_CODE}`,
  publishes: ['energy.ratedPower', 'energy.aep', 'energy.capacityFactor', 'energy.availability', 'energy.aepGross', 'energy.aepP90'],
  searchTerms: 'wind farm energy yield assessment P50 P90 loss categories uncertainty',
  inputs: [
    ...curveInputs(),
    { key: 'nUnits', label: 'Number of turbines', type: 'int', default: 60, min: 1, max: 400, group: 'Plant', role: 'geometry', design: true },
    { key: 'windSpec', label: 'Wind climate given by', type: 'select', default: 'weibull', group: 'Wind resource', role: 'bc', options: [{ value: 'weibull', label: 'Weibull scale A and shape k' }, { value: 'mean', label: 'Mean wind speed and shape k' }], help: 'Time-varying resource boundary condition, expressed as the long-term distribution at hub height.' },
    site('A', { role: 'bc', showIf: { key: 'windSpec', in: ['weibull'] } }),
    { key: 'hubMean', label: 'Mean hub-height wind speed', unit: 'm/s', type: 'number', default: 9.4, min: 3, max: 15, group: 'Wind resource', role: 'bc', bind: 'site.hubWindMean', showIf: { key: 'windSpec', in: ['mean'] }, uncertain: { dist: 'normal', cov: 0.05 } },
    site('k', { role: 'bc' }),
    { key: 'Tair', label: 'Mean air temperature at hub height', unit: '°C', type: 'number', default: 10, min: -20, max: 40, group: 'Wind resource', role: 'environment', bind: 'site.Tair', help: 'Sets the site air density at standard pressure; the curve is corrected from 1.225 kg/m³ per IEC 61400-12-1.' },
    { key: 'monthly', label: 'Monthly mean wind speed / annual mean', type: 'series', default: { x: MONTHS.map((_, i) => i + 1), y: [...SEASON] }, xLabel: 'Month', xUnit: '-', yLabel: 'Ratio', yUnit: '-', group: 'Wind resource', role: 'bc', help: 'Twelve ratios (January to December). All ones = no seasonality.' },
    pct('Wake (internal array) loss', 'wake', 8, 40, 'From the wake and array module when it has been run; 5–12 % is usual for large offshore arrays.', { bind: 'energy.wakeLossFrac', uncertain: { dist: 'normal', cov: 0.15 } }),
    pct('Blockage loss', 'blockage', 1, 10, 'Upstream induction of the array as a whole (0.5–2 % in recent assessments).'),
    pct('Electrical loss to the metering point', 'elec', 2.5, 15, 'Array cables, transformers and export system.', { bind: 'cable.lossesFrac' }),
    { key: 'avail', label: 'Energy-weighted availability (turbines, balance of plant, grid)', unit: '%', type: 'number', default: 94, min: 50, max: 100, group: 'Losses', role: 'ic', bind: 'energy.availability', uncertain: { dist: 'normal', cov: 0.02 }, help: 'Operational availability state the production estimate starts from; from the installation and maintenance module when run.' },
    pct('Environmental (icing, degradation, high-wind hysteresis)', 'env', 1.5, 15, 'Blade degradation, high-wind hysteresis, site-specific environmental shutdowns.'),
    pct('Curtailment (grid, noise, wildlife)', 'curtail', 0.5, 30, 'Planned curtailments and grid constraints.'),
    unc('Wind-speed measurement uncertainty', 'uMeas', 2.5, 'Fraction of mean wind speed (cup / lidar calibration, mounting).'),
    unc('Long-term correction uncertainty', 'uLong', 2.5, 'Measure–correlate–predict and representativeness of the reference period.'),
    unc('Vertical and horizontal extrapolation uncertainty', 'uExtrap', 2.0, 'From the mast or lidar position and height to every hub.'),
    unc('Power-curve uncertainty (energy)', 'uCurve', 4.0, 'Fraction of energy; larger for a generic curve than for a warranted one.'),
    unc('Wake-model and loss uncertainty (energy)', 'uLoss', 3.0, 'Fraction of energy on wakes, electrical losses and availability.'),
    unc('Inter-annual variability of energy (one year)', 'uIav', 6.0, 'Standard deviation of annual energy from year to year.'),
    { key: 'years', label: 'Averaging period for the exceedance values', unit: 'yr', type: 'int', default: 10, min: 1, max: 40, group: 'Uncertainty', role: 'criteria', help: 'Lenders usually ask for 1-year and 10-year P90; inter-annual variability shrinks with √n.' },
    { key: 'cfTarget', label: 'Required net capacity factor', unit: '%', type: 'number', default: 40, min: 5, max: 70, group: 'Acceptance', role: 'criteria', help: `Project threshold for the P50 net capacity factor. ${CONFIRM}` },
  ],
  compute(inp) {
    const c0 = getCurve(inp);
    const A = inp.windSpec === 'mean' ? weibullScaleFromMean(inp.hubMean, inp.weibK) : inp.weibA, k = inp.weibK;
    const rho = airDensity({ pressure: 101325, temperatureC: inp.Tair });
    const curve = run(() => airDensityCorrection({ u: c0.u, p: c0.p, ct: c0.ct }, { rho, rho0: 1.225, regulation: 'pitch' }), 'Air-density correction');
    const wb = weibull(A, k);
    const one = run(() => aepWeibull({ curve, A, k }), 'Annual energy integration');
    const bins = run(() => aepBins({ curve, A, k }), 'Method of bins');
    const n = inp.nUnits, rated = c0.ratedPower * n, grossWh = one.aepWh * n;
    need(grossWh > 0, 'The gross energy is zero: the power curve does not overlap the wind-speed distribution.');
    const losses = [['Wake', inp.wake], ['Blockage', inp.blockage], ['Electrical', inp.elec], ['Availability', 100 - inp.avail], ['Environmental', inp.env], ['Curtailment', inp.curtail]].map(([name, l]) => ({ name, loss: l / 100 }));
    const cas = lossCascade({ grossWh, losses, ratedPower: rated });
    // sensitivity of energy to wind speed by central difference on the Weibull scale (±2 %)
    const eAt = (f) => aepWeibull({ curve, A: A * f, k, n: 1000 }).aepWh;
    // (floored: at very high mean wind speeds more wind adds time above cut-out and the sensitivity tends to zero or below)
    const sens = Math.max((eAt(1.02) - eAt(0.98)) / (0.04 * one.aepWh), 0.05);
    const mk = (yrs) => aepUncertainty({ p50Wh: cas.netWh, windSpeed: { 'Wind-speed measurement': inp.uMeas / 100, 'Long-term correction': inp.uLong / 100, 'Spatial extrapolation': inp.uExtrap / 100 }, sensitivity: sens,
      energy: { 'Power curve': inp.uCurve / 100, 'Wake model and losses': inp.uLoss / 100 }, interAnnual: inp.uIav / 100, years: yrs, levels: [75, 90, 99] });
    const un = mk(inp.years), un1 = mk(1);
    need(un.p90Wh > 0, 'The combined uncertainty is so large that the P90 energy is not positive: review the uncertainty inputs.');
    const cf = cas.capacityFactor;
    // monthly energy: Weibull scale scaled by the monthly ratio, then normalised so that the months sum to the annual net energy
    const ratios = (inp.monthly?.y || []).map(Number);
    need(ratios.length === 12 && ratios.every((r) => r > 0.2 && r < 3), 'The monthly series needs twelve ratios between 0.2 and 3.');
    const raw = ratios.map((r, i) => aepWeibull({ curve, A: A * r, k, n: 600 }).aepWh * DAYS[i] / 365.25);
    const sRaw = raw.reduce((a, b) => a + b, 0);
    const monthWh = raw.map((v) => v / sRaw * cas.netWh);
    // resource utilisation: net energy over the kinetic energy of the wind through the swept area
    const sweep = Math.PI * c0.D ** 2 / 4;
    const windWh = 0.5 * rho * sweep * wb.cubeMean * HOURS_YEAR * n;
    const us = linspace(0, 30, 121);
    const pdfMax = Math.max(...us.map((u) => wb.pdf(u)));
    const eDens = us.map((u) => wb.pdf(u) * curvePower(curve, u));
    const eMax = Math.max(...eDens);
    const gwh = (wh) => round(wh / 1e9, 5);
    const lostAvail = cas.steps.find((s) => s.name === 'Availability').lostWh, lostWake = cas.steps.find((s) => s.name === 'Wake').lostWh;
    const onePct = cas.netWh * 0.01 / 1e6;
    const checks = [check('Net capacity factor (P50) against the project requirement', inp.cfTarget, cf * 100, '%', 'Project criterion; IEC 61400-15-1 energy-yield reporting', 'Demand is the required capacity factor, capacity the predicted P50 value. Below the requirement the plant produces less energy per installed megawatt than the business case assumes.')];
    return {
      kpis: [
        { key: 'aepGross', label: 'Gross annual energy', value: gwh(grossWh), unit: 'GWh/yr', explain: 'Energy of all turbines in the undisturbed wind climate at site air density, before any loss. It is the starting point of the loss cascade.' },
        { key: 'aepNet', label: 'Net annual energy (P50)', value: gwh(cas.netWh), unit: 'GWh/yr', explain: `Gross energy after wake, blockage, electrical, availability, environmental and curtailment losses: the central estimate. ${PXX}` },
        { key: 'aepP90', label: `Net annual energy P90 (${inp.years}-year average)`, value: gwh(un.p90Wh), unit: 'GWh/yr', explain: `Energy exceeded with 90 % probability over the averaging period — the value lenders size debt on. ${PXX}` },
        { key: 'aepP75', label: `Net annual energy P75 (${inp.years}-year average)`, value: gwh(un.p75Wh), unit: 'GWh/yr', explain: `Energy exceeded with 75 % probability. ${PXX}` },
        { key: 'aepP90one', label: 'Net annual energy P90 (any single year)', value: gwh(un1.p90Wh), unit: 'GWh/yr', explain: `Energy exceeded with 90 % probability in any one year, carrying the full inter-annual variability; used for debt-service cover in the worst year. ${PXX}` },
        { key: 'cf', label: 'Net capacity factor (P50)', value: round(cf * 100, 4), unit: '%', limit: inp.cfTarget, sense: 'min', explain: 'Net energy over rated power times 8760 h. Offshore wind typically reaches 40–55 %; it is the single number that summarises site and turbine match.' },
        { key: 'cfGross', label: 'Gross capacity factor', value: round(grossWh / (rated * HOURS_YEAR) * 100, 4), unit: '%', explain: 'Capacity factor of an isolated, always-available turbine in this wind climate: the upper bound that losses eat into.' },
        { key: 'flh', label: 'Full-load hours', value: round(cas.fullLoadHours, 4), unit: 'h/yr', explain: 'Net energy divided by rated power: the hours at rated output that would give the same energy.' },
        { key: 'totalLoss', label: 'Total loss, gross to net', value: round(cas.totalLoss * 100, 4), unit: '%', explain: 'Combined effect of the multiplicative loss chain. Values of 15–22 % are typical for offshore wind including wakes and availability.' },
        { key: 'sigma', label: `Uncertainty of net energy (${inp.years}-year)`, value: round(un.sigma * 100, 4), unit: '%', explain: 'One standard deviation of the net energy as a fraction of P50, combining wind-speed, power-curve, loss and inter-annual components by root-sum-square.' },
        { key: 'sens', label: 'Energy sensitivity to wind speed', value: round(sens, 4), unit: '-', explain: 'Percent change of energy per percent change of mean wind speed. Values of 1.3–2 are normal; high-wind sites sit lower because more time is spent at rated power.' },
        { key: 'meanWind', label: 'Mean hub-height wind speed', value: round(wb.mean, 4), unit: 'm/s', explain: 'Mean of the Weibull distribution, A·Γ(1 + 1/k).' },
        { key: 'utilisation', label: 'Resource utilisation (net energy / wind energy through the rotors)', value: round(cas.netWh / windWh * 100, 4), unit: '%', explain: 'Overall conversion efficiency from the kinetic energy flowing through the swept area to energy at the meter. It is bounded by the Betz limit (59.3 %) and by the time spent above rated wind speed.' },
        { key: 'daily', label: 'Mean daily generation', value: round(cas.netWh / 365.25 / 1e6, 5), unit: 'MWh/day', explain: 'Net annual energy divided by the days in a year; actual days swing between zero and rated output for 24 h.' },
        { key: 'meanPower', label: 'Mean net power (instantaneous average)', value: round(cas.netWh / HOURS_YEAR / 1e6, 5), unit: 'MW', explain: 'Long-term average export power of the plant: the level a firm delivery contract can be compared with.' },
        { key: 'binsDiff', label: 'Method-of-bins minus integration', value: round((bins.aepWh / one.aepWh - 1) * 100, 3), unit: '%', explain: 'Difference between the IEC 61400-12-1 method of bins on the table points and the fine numerical integration. More than about 1 % means the curve table is too coarse.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Power curve, wind distribution and energy density', x: { label: 'Hub-height wind speed', unit: 'm/s' }, y: { label: 'Normalised', unit: '-' },
          series: [{ name: 'Power / rated', x: us.map((v) => round(v, 3)), y: us.map((u) => round(curvePower(curve, u) / c0.ratedPower, 4)) }, { name: 'Weibull density / peak', x: us.map((v) => round(v, 3)), y: us.map((u) => round(wb.pdf(u) / pdfMax, 4)) }, { name: 'Energy density / peak', x: us.map((v) => round(v, 3)), y: eDens.map((v) => round(v / eMax, 4)) }],
          explain: 'Site-density power curve, the wind-speed probability density and their product. The peak of the energy-density curve shows which wind speeds produce the energy: curve accuracy matters most there.' },
        { type: 'bar', title: 'Loss cascade from gross to net', categories: ['Gross', ...cas.steps.map((s) => `after ${s.name.toLowerCase()}`)], series: [{ name: 'Energy (GWh/yr)', y: [gwh(grossWh), ...cas.steps.map((s) => gwh(s.energyWh))] }], y: { label: 'Annual energy', unit: 'GWh/yr' },
          explain: 'Annual energy remaining after each loss in the chain. The largest steps (normally wakes and availability) are where design and O&M effort pays back.' },
        { type: 'bar', title: 'Monthly net generation', categories: MONTHS, series: [{ name: 'Energy (GWh)', y: monthWh.map(gwh) }], y: { label: 'Energy', unit: 'GWh' },
          explain: 'Net energy per calendar month from the seasonal wind-speed ratios. Winter months carry the yield; summer is the access window for planned maintenance.' },
        { type: 'bar', title: 'Exceedance values of net annual energy', categories: ['P50', 'P75', 'P90', 'P99', 'P90 single year'], series: [{ name: 'Energy (GWh/yr)', y: [gwh(cas.netWh), gwh(un.p75Wh), gwh(un.p90Wh), gwh(un.p99Wh), gwh(un1.p90Wh)] }], y: { label: 'Annual energy', unit: 'GWh/yr' },
          explain: `Net energy at rising confidence. ${PXX} The gap between P50 and P90 is the price of uncertainty in the financing.` },
        { type: 'tornado', title: 'Uncertainty contributions to net energy (± one standard deviation)', base: gwh(cas.netWh), items: un.components.map((q) => ({ name: q.name, low: gwh(cas.netWh * (1 - q.sigma)), high: gwh(cas.netWh * (1 + q.sigma)) })),
          explain: 'Each bar is the swing of net energy for one standard deviation of a single uncertainty component. The longest bars are the ones worth reducing with more measurement or a warranted curve.' },
      ],
      tables: [
        { title: 'Loss cascade', columns: [{ key: 'name', label: 'Loss' }, { key: 'loss', label: 'Loss', unit: '%' }, { key: 'lost', label: 'Energy lost', unit: 'GWh/yr' }, { key: 'left', label: 'Energy remaining', unit: 'GWh/yr' }],
          rows: cas.steps.map((s) => ({ name: s.name, loss: round(s.loss * 100, 4), lost: gwh(s.lostWh), left: gwh(s.energyWh) })),
          explain: 'One row per loss in the order applied (each acts on what the previous one left). Use it to reconcile the yield report with the financial model line by line.' },
        { title: 'Monthly and daily generation', columns: [{ key: 'm', label: 'Month' }, { key: 'ratio', label: 'Wind ratio', unit: '-' }, { key: 'e', label: 'Net energy', unit: 'GWh' }, { key: 'd', label: 'Mean daily energy', unit: 'MWh/day' }, { key: 'cf', label: 'Capacity factor', unit: '%' }],
          rows: MONTHS.map((m, i) => ({ m, ratio: ratios[i], e: gwh(monthWh[i]), d: round(monthWh[i] / DAYS[i] / 1e6, 5), cf: round(monthWh[i] / (rated * DAYS[i] * 24) * 100, 4) })),
          explain: 'Net energy, mean daily energy and capacity factor per month. Compare the low months with the maintenance plan and with seasonal price curves.' },
        { title: 'Overall system performance and operating limits', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Installed capacity', v: round(rated / 1e6, 5), u: 'MW' }, { q: 'Generating units', v: n, u: '-' }, { q: 'Cut-in / cut-out wind speed', v: `${fmt(c0.cutIn)} / ${fmt(c0.cutOut)}`, u: 'm/s' }, { q: 'Rated wind speed (first speed at rated power)', v: fmt(c0.u[c0.p.findIndex((p) => p >= 0.999 * c0.ratedPower)] ?? c0.cutOut), u: 'm/s' },
            { q: 'Time between cut-in and cut-out', v: round((wb.cdf(c0.cutOut) - wb.cdf(c0.cutIn)) * 100, 4), u: '%' }, { q: 'Time at rated power', v: round((wb.cdf(c0.cutOut) - wb.cdf(c0.u[c0.p.findIndex((p) => p >= 0.999 * c0.ratedPower)] ?? c0.cutOut)) * 100, 4), u: '%' },
            { q: 'Gross / net (P50) annual energy', v: `${fmt(grossWh / 1e9, 4)} / ${fmt(cas.netWh / 1e9, 4)}`, u: 'GWh/yr' }, { q: `P90 (${inp.years}-year) annual energy`, v: fmt(un.p90Wh / 1e9, 4), u: 'GWh/yr' }, { q: 'Net capacity factor against requirement', v: `${fmt(cf * 100)} / ${inp.cfTarget}`, u: '%' }, { q: 'Availability assumed', v: inp.avail, u: '%' }, { q: 'Design compliance (capacity-factor criterion)', v: cf * 100 >= inp.cfTarget ? 'met' : 'not met', u: '' }],
          explain: 'One-page summary of the plant: capacity, operating wind-speed limits, time in operation and at rated power, energy at P50 and P90 (exceedance values) and compliance with the yield requirement. Structural and electrical compliance is reported by the checks of the other modules of the suite.' },
        { title: 'Exceedance values', columns: [{ key: 'level', label: 'Level' }, { key: 'n', label: `${inp.years}-year average`, unit: 'GWh/yr' }, { key: 'one', label: 'Single year', unit: 'GWh/yr' }, { key: 'cf', label: `Capacity factor (${inp.years}-year)`, unit: '%' }],
          rows: [['P50', cas.netWh, cas.netWh], ['P75', un.p75Wh, un1.p75Wh], ['P90', un.p90Wh, un1.p90Wh], ['P99', un.p99Wh, un1.p99Wh]].map(([level, a, b]) => ({ level, n: gwh(a), one: gwh(b), cf: round(a / (rated * HOURS_YEAR) * 100, 4) })),
          explain: `Net energy at each exceedance level for the averaging period and for any single year. ${PXX}` },
      ],
      notes: [
        `Power curve: ${c0.label}, ${fmt(c0.ratedPower / 1e6)} MW × ${n} = ${fmt(rated / 1e6, 4)} MW installed; rotor diameter ${fmt(c0.D, 4)} m; site air density ${fmt(rho, 4)} kg/m³ (curve corrected from 1.225 kg/m³).`,
        `Wind climate: Weibull A = ${fmt(A, 4)} m/s, k = ${fmt(k)}, mean ${fmt(wb.mean, 4)} m/s.`,
        `Largest losses: wakes ${fmt(lostWake / 1e9)} GWh/yr and unavailability ${fmt(lostAvail / 1e9)} GWh/yr.`,
        `Exceedance convention: ${PXX} P90 (${inp.years}-year) = P50 × (1 − 1.2816 × ${fmt(un.sigma, 3)}).`,
      ],
      advice: [
        adv('engineering', cf * 100 >= inp.cfTarget ? 'ok' : 'fail', `Net P50 energy is ${fmt(cas.netWh / 1e9, 4)} GWh/yr (capacity factor ${fmt(cf * 100)} %), P90 over ${inp.years} years ${fmt(un.p90Wh / 1e9, 4)} GWh/yr; total losses ${fmt(cas.totalLoss * 100)} %.`,
          cf * 100 >= inp.cfTarget ? `Confirm the wake loss (${inp.wake} %) with the array module and replace the ${c0.label} by the warranted supplier curve; the method of bins and the integration differ by ${fmt((bins.aepWh / one.aepWh - 1) * 100, 2)} %.` : `The capacity factor is below the ${inp.cfTarget} % requirement: a larger rotor per megawatt (lower specific power), wider spacing or higher availability are the levers — one point of availability is worth ${fmt(grossWh * 0.01 * (cas.netWh / grossWh) / (inp.avail / 100) / 1e6)} MWh/yr.`,
          'IEC 61400-15-1 (energy-yield assessment and loss categories) and IEC 61400-12-1:2017 (power curve and air-density normalisation), latest editions known.',
          `The net yield displaces about ${fmt(cas.netWh / 1e6 * CO2_GRID / 1e3, 4)} kt CO₂ a year from gas-fired generation (screening factor ${CO2_GRID} t CO₂/MWh — use the grid factor of the project).`),
        adv('economic', 'info', `One percent of net energy is ${fmt(onePct, 4)} MWh a year; the financing haircut from P50 to P90 is ${fmt((cas.netWh - un.p90Wh) / 1e9, 3)} GWh/yr (${fmt((1 - un.p90Wh / cas.netWh) * 100)} %).`,
          `Reduce the largest uncertainty first (${un.components.reduce((a, b) => (b.sigma > a.sigma ? b : a)).name}, ${fmt(un.components.reduce((a, b) => (b.sigma > a.sigma ? b : a)).sigma * 100)} % of energy): a second year of on-site lidar or a warranted power curve typically pays back through a higher P90 debt capacity.`,
          'Project-finance practice sizes debt on P90 (or P99 one-year) energy with a debt-service cover ratio; LCOE per IEA / IRENA convention uses the P50 net energy as denominator; value-of-information reasoning for extra measurement.',
          energyCo2Text(onePct)),
      ],
      publish: { 'energy.ratedPower': rated, 'energy.aep': cas.netWh, 'energy.capacityFactor': cf, 'energy.availability': inp.avail / 100, 'energy.aepGross': grossWh, 'energy.aepP90': un.p90Wh },
    };
  },
  verify: [
    { name: 'Constant 1 MW between 3 and 25 m/s in a Rayleigh climate: AEP = P·8760·[exp(−(3/A)²) − exp(−(25/A)²)]',
      source: 'Closed-form integral of the Weibull density (k = 2, A = 10 m/s); Manwell, McGowan & Rogers, Wind Energy Explained §2.4',
      // exp(−0.09) − exp(−6.25) = 0.913931 − 0.001930 = 0.912001 → 1 MW × 8760 h × 0.912001 = 7989.1 MWh = 7.9891 GWh; capacity factor 91.20 %
      // 15 °C at 101325 Pa is the reference density 1.225 kg/m³, so no curve shift; all losses and unavailability set to zero
      inputs: { curveSource: 'table', curve: [{ u: 3, p: 1000, ct: 0.5 }, { u: 14, p: 1000, ct: 0.5 }, { u: 25, p: 1000, ct: 0.5 }], D: 60, nUnits: 1, windSpec: 'weibull', weibA: 10, weibK: 2, Tair: 15, wake: 0, blockage: 0, elec: 0, avail: 100, env: 0, curtail: 0 },
      expect: [{ kpi: 'aepGross', value: 7.9891, tol: 0.003 }, { kpi: 'aepNet', value: 7.9891, tol: 0.003 }, { kpi: 'cf', value: 91.2, tol: 0.003 }] },
    { name: 'Loss cascade: net = gross × Π(1 − lᵢ) and P90 = P50·(1 − 1.2816 σ)',
      source: 'Hand calculation: (1−0.10)(1−0.05)(1−0.06) = 0.8037; σ = 10 % energy only → P90/P50 = 0.87184',
      // gross 7.9891 GWh × 0.8037 = 6.4208 GWh; P90 = 6.4208 × (1 − 1.2816 × 0.10) = 5.5979 GWh
      inputs: { curveSource: 'table', curve: [{ u: 3, p: 1000, ct: 0.5 }, { u: 14, p: 1000, ct: 0.5 }, { u: 25, p: 1000, ct: 0.5 }], D: 60, nUnits: 1, weibA: 10, weibK: 2, Tair: 15, wake: 10, blockage: 0, elec: 5, avail: 94, env: 0, curtail: 0,
        uMeas: 0, uLong: 0, uExtrap: 0, uCurve: 10, uLoss: 0, uIav: 0, years: 1 },
      expect: [{ kpi: 'aepNet', value: 6.4208, tol: 0.003 }, { kpi: 'aepP90', value: 5.5979, tol: 0.003 }] },
  ],
};
