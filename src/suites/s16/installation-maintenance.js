// Installation and maintenance operability: weather windows and accessibility from an hourly wave and wind record
// (imported or synthesised from the site statistics), Monte Carlo duration of the installation campaign, and a
// time-stepping operations-and-maintenance simulation giving technical availability, downtime and vessel use.
import { operability, persistenceTable, campaignMonteCarlo } from '../../engine/metocean/windows.js';
import { simulateWindOM } from '../../engine/ram/maintenance.js';
import { weibullScaleFromMean } from '../../engine/aero/aep.js';
import { makeRng } from '../../engine/math/rng.js';
import { curveInputs, getCurve } from './_curve.js';
import { costInput, money, site, adv, check, worst, sev, fmt, round, run, need, interp, thin, LIMIT_CODE, CONFIRM, vesselCo2Text, energyCo2Text, CO2_MGO } from './_common.js';

const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
// Indicative installation cycle per position and failure data per turbine-year (order of magnitude of published
// offshore surveys, e.g. Carroll et al. 2016) — replace with the contractor's method statement and fleet data.
// Wave limits apply to positioning and jacking; once the vessel is elevated the lifts are wind-limited.
const OPS = [['Position and jack up', 6, 1.75, 20], ['Monopile upending and driving', 12, 2.0, 16], ['Transition piece and tower', 8, 3.0, 14], ['Nacelle lift', 5, 3.0, 12], ['Three blade lifts', 9, 3.0, 12], ['Mechanical completion and jack down', 5, 1.75, 20]];
const FAIL = [['Manual reset on site', 3.0, 3, 'ctv', 4], ['Minor repair', 2.5, 8, 'ctv', 24], ['Major repair', 0.3, 30, 'ctv', 96], ['Major replacement (heavy lift)', 0.06, 60, 'jackup', 240]];

/** Hourly synthetic year of hub wind speed and Hs: Gaussian AR(1) processes rank-mapped to Weibull marginals, with a seasonal factor. */
function synthWeather({ A, k, hsMean, hsShape, tau, season, corr, seed, n = 8760 }) {
  const rng = makeRng(seed), phi = Math.exp(-1 / tau), q = Math.sqrt(1 - phi * phi);
  const z1 = new Float64Array(n), z2 = new Float64Array(n);
  let a = rng.normal(), b = rng.normal();
  for (let i = 0; i < n; i++) { a = phi * a + q * rng.normal(); b = phi * b + q * rng.normal(); z1[i] = a; z2[i] = corr * a + Math.sqrt(1 - corr * corr) * b; }
  // rank mapping gives the target marginal exactly (plotting position (rank + 0.5)/n)
  const map = (z, scale, shape) => { const idx = Array.from(z, (_, i) => i).sort((p, r) => z[p] - z[r]), out = new Float64Array(n); idx.forEach((id, rank) => { out[id] = scale * (-Math.log(1 - (rank + 0.5) / n)) ** (1 / shape); }); return out; };
  const wind = map(z1, A, k), hs = map(z2, weibullScaleFromMean(hsMean, hsShape), hsShape);
  const t0 = Date.UTC(2023, 0, 1), timeMs = new Float64Array(n);
  for (let i = 0; i < n; i++) { const f = 1 + season * Math.cos(2 * Math.PI * (i / 24 - 15) / 365); wind[i] *= f; hs[i] *= f; timeMs[i] = t0 + i * 3600e3; }
  return { wind, hs, timeMs };
}

export default {
  id: 's16-installation-maintenance', suite: 16, title: 'Installation and maintenance operability',
  summary: 'Weather windows and accessibility at the site, duration of the installation campaign with weather downtime, and the technical availability, downtime and vessel demand of operations and maintenance.',
  covers: [16], fidelity: 'numerical',
  standards: ['DNV-ST-N001 (2021) Marine operations and marine warranty §2 (operation reference period, weather-restricted operations, alpha factors)', 'DNV-RP-N103 (2017) §9 (weather criteria and availability analysis)', 'IEC TS 61400-26-1:2019 (time-based and production-based availability)', 'IEC 61400-3-1:2019 §12–13 (assembly, installation, maintenance)', 'Dinwoodie et al. (2015) Wind Engineering 39(1) (O&M simulation reference cases)', 'Carroll, McDonald & McMillan (2016) Wind Energy 19 (failure rates and repair times)'],
  models: ['Weather-window (persistence) statistics of a time series', 'Operability with uninterrupted reference period', 'Monte Carlo replay of sequential weather-restricted operations', 'Gaussian autoregressive weather synthesis with Weibull marginals and seasonal modulation', 'Poisson failure processes per sub-assembly', 'Time-stepping O&M simulation with vessel weather limits, lead times and vessel counts', 'Time-based and energy-based availability'],
  limits: `The synthetic weather year reproduces the Weibull marginals, an exponential persistence time and a seasonal factor, not real storm sequences or swell: import a measured or hindcast series for decisions. Installation operations cannot be suspended once started; port logistics, technical downtime and learning effects are not modelled. The O&M simulation has no travel time, shifts or daylight restriction and treats repairs as pre-emptible in bad weather; failure data are indicative. ${LIMIT_CODE}`,
  publishes: ['install.vesselDays', 'install.operabilityFrac', 'energy.availability', 'om.vesselDaysPerYear', 'om.energyLostFrac'],
  searchTerms: 'offshore wind installation weather window operability O&M availability simulation',
  inputs: [
    ...curveInputs({ group: 'Wind farm' }),
    { key: 'nUnits', label: 'Number of turbines', type: 'int', default: 60, min: 1, max: 300, group: 'Wind farm', role: 'geometry', bind: 'array.count' },
    { key: 'hsSeries', label: 'Hourly significant wave height record (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Time', xUnit: 'h', yLabel: 'Hs', yUnit: 'm', group: 'Weather record', role: 'bc', help: 'Measured or hindcast hourly Hs over one or more whole years. With an imported record the wind speed is taken from the companion series below, or synthesised if that is empty.' },
    { key: 'windSeries', label: 'Hourly hub-height wind speed record (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Time', xUnit: 'h', yLabel: 'Wind speed', yUnit: 'm/s', group: 'Weather record', role: 'bc' },
    site('A', { group: 'Weather record', role: 'bc' }), site('k', { group: 'Weather record', role: 'bc' }),
    { key: 'hsMean', label: 'Annual mean significant wave height', unit: 'm', type: 'number', default: 1.4, min: 0.3, max: 4, group: 'Weather record', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Central North Sea about 1.5–2.0 m, Irish Sea 1.0–1.3 m, exposed Atlantic 2.5 m and more.' },
    { key: 'hsShape', label: 'Weibull shape of Hs', unit: '-', type: 'number', default: 1.5, min: 1, max: 2.5, group: 'Weather record', role: 'bc' },
    { key: 'tau', label: 'Weather persistence time', unit: 'h', type: 'number', default: 36, min: 6, max: 120, group: 'Weather record', role: 'bc', help: 'Decorrelation time of the weather systems; longer persistence gives longer but fewer windows.' },
    { key: 'season', label: 'Seasonal amplitude of wind and waves', unit: '-', type: 'number', default: 0.3, min: 0, max: 0.6, group: 'Weather record', role: 'bc', help: 'Winter maximum / annual mean − 1. 0 = no seasonality.' },
    { key: 'seed', label: 'Random seed', type: 'int', default: 5, min: 1, max: 99999, group: 'Weather record', role: 'ic' },
    { key: 'ops', label: 'Installation cycle per position', type: 'table', group: 'Installation', role: 'operation', columns: [{ key: 'name', label: 'Operation', type: 'text' }, { key: 'hours', label: 'Net duration', unit: 'h', type: 'number' }, { key: 'hs', label: 'Hs limit', unit: 'm', type: 'number' }, { key: 'wind', label: 'Hub wind limit', unit: 'm/s', type: 'number' }],
      default: OPS.map(([name, hours, hs, wind]) => ({ name, hours, hs, wind })), help: 'Operations in sequence, each needing an uninterrupted window of its net duration plus contingency within its limits (vessel capability and lifting criteria).' },
    { key: 'contingency', label: 'Contingency on operation duration', unit: '%', type: 'number', default: 50, min: 0, max: 100, group: 'Installation', role: 'criteria', help: `DNV-ST-N001: reference period TR = planned time + contingency, normally at least 50 % of the planned time. ${CONFIRM}` },
    { key: 'alpha', label: 'Alpha factor on the limits (forecast uncertainty)', unit: '-', type: 'number', default: 0.85, min: 0.5, max: 1, group: 'Installation', role: 'criteria', help: `DNV-ST-N001 §2.6: 0.65–0.9 depending on duration, forecast level and monitoring. ${CONFIRM}` },
    { key: 'transit', label: 'Transit and loading time per position', unit: 'h', type: 'number', default: 20, min: 0, max: 200, group: 'Installation', role: 'operation', help: 'Share of the port round trip and load-out per position (not weather-restricted here).' },
    { key: 'nVessels', label: 'Installation vessels working in parallel', type: 'int', default: 2, min: 1, max: 6, group: 'Installation', role: 'operation', design: true, help: 'For example a foundation vessel and a turbine installation vessel. Vessel days are unchanged; the elapsed campaign time divides by this number.' },
    { key: 'startMonth', label: 'Campaign start month', type: 'int', default: 4, min: 1, max: 12, group: 'Installation', role: 'ic', help: 'Initial condition of the campaign: 4 = April.' },
    { key: 'seasonDays', label: 'Available campaign duration', unit: 'days', type: 'number', default: 240, min: 20, max: 730, group: 'Installation', role: 'criteria', help: 'Vessel charter window or season: the P90 elapsed campaign duration should fit.' },
    { key: 'wowLimit', label: 'Acceptable waiting-on-weather share of the cycle', unit: '%', type: 'number', default: 60, min: 10, max: 90, group: 'Installation', role: 'criteria', help: 'Planning threshold: above it the vessel and its limits are mismatched with the site.' },
    { key: 'fail', label: 'Failure and repair data per turbine', type: 'table', group: 'Maintenance', role: 'operation', columns: [{ key: 'name', label: 'Failure class', type: 'text' }, { key: 'rate', label: 'Failures', unit: '1/yr', type: 'number' }, { key: 'hours', label: 'Repair time', unit: 'h', type: 'number' }, { key: 'vessel', label: 'Vessel (ctv / jackup)', type: 'text' }, { key: 'lead', label: 'Lead time', unit: 'h', type: 'number' }],
      default: FAIL.map(([name, rate, hours, vessel, lead]) => ({ name, rate, hours, vessel, lead })), help: 'Failure rate per turbine-year, active repair time, access vessel class and lead time for parts and planning (the jack-up mobilisation is added separately).' },
    { key: 'ctvHs', label: 'Access-vessel Hs limit', unit: 'm', type: 'number', default: 1.75, min: 0.5, max: 4, group: 'Maintenance', role: 'operation', design: true, help: 'CTV 1.25–1.75 m; service operation vessel with motion-compensated gangway 2.5–3.5 m.' },
    { key: 'ctvWind', label: 'Access wind limit (hub height)', unit: 'm/s', type: 'number', default: 18, min: 8, max: 30, group: 'Maintenance', role: 'operation' },
    { key: 'nCtv', label: 'Access vessels (teams working in parallel)', type: 'int', default: 3, min: 1, max: 20, group: 'Maintenance', role: 'operation', design: true },
    { key: 'juHs', label: 'Jack-up Hs limit', unit: 'm', type: 'number', default: 2.0, min: 0.5, max: 3.5, group: 'Maintenance', role: 'operation' },
    { key: 'juWind', label: 'Jack-up lifting wind limit (hub height)', unit: 'm/s', type: 'number', default: 11, min: 6, max: 20, group: 'Maintenance', role: 'operation' },
    { key: 'juMob', label: 'Jack-up mobilisation time', unit: 'days', type: 'number', default: 45, min: 0, max: 200, group: 'Maintenance', role: 'operation', help: 'On-demand charter: waiting for a vessel dominates the downtime of major replacements.' },
    { key: 'availTarget', label: 'Required energy-based availability', unit: '%', type: 'number', default: 93, min: 70, max: 99.5, group: 'Acceptance', role: 'criteria' },
    { key: 'omYears', label: 'Simulated O&M period', unit: 'yr', type: 'int', default: 3, min: 1, max: 25, group: 'Solver', role: 'control' },
    { key: 'omSims', label: 'O&M Monte Carlo histories', type: 'int', default: 2, min: 1, max: 50, group: 'Solver', role: 'control' },
    { key: 'nTrials', label: 'Installation Monte Carlo trials', type: 'int', default: 200, min: 20, max: 5000, group: 'Solver', role: 'control' },
    costInput('energyPrice', 'Value of delivered energy', 'USD/MWh', 80, 1000, 'Used to express downtime as lost revenue.'),
    costInput('instDayRate', 'Installation vessel spread day rate', 'USD/day', 250000, 2e6, 'Jack-up or heavy-lift vessel with crew and equipment.'),
    costInput('ctvDayRate', 'Crew-transfer vessel day rate', 'USD/day', 6000, 1e5),
    costInput('juDayRate', 'Maintenance jack-up day rate', 'USD/day', 150000, 2e6, 'For major-component exchange.'),
  ],
  compute(inp) {
    const c = getCurve(inp);
    // weather record: imported series or synthetic year
    const hsIn = (inp.hsSeries?.y || []).map(Number), wIn = (inp.windSeries?.y || []).map(Number);
    let wx;
    if (hsIn.length >= 720) {
      need(hsIn.every((v) => v >= 0 && v < 25), 'The imported Hs record contains values outside 0–25 m.');
      const n = hsIn.length, syn = synthWeather({ A: inp.weibA, k: inp.weibK, hsMean: inp.hsMean, hsShape: inp.hsShape, tau: inp.tau, season: inp.season, corr: 0.7, seed: inp.seed, n });
      if (wIn.length) need(wIn.length === n && wIn.every((v) => v >= 0 && v < 80), 'The wind record must have the same length as the Hs record and values between 0 and 80 m/s.');
      wx = { hs: Float64Array.from(hsIn), wind: wIn.length ? Float64Array.from(wIn) : syn.wind, timeMs: syn.timeMs, imported: true };
    } else {
      need(hsIn.length === 0, 'The imported Hs record is too short: at least 720 hourly values (ideally whole years) are needed.');
      wx = { ...synthWeather({ A: inp.weibA, k: inp.weibK, hsMean: inp.hsMean, hsShape: inp.hsShape, tau: inp.tau, season: inp.season, corr: 0.7, seed: inp.seed }), imported: false };
    }
    const series = { hs: wx.hs, wind: wx.wind }, dtHours = 1, n = wx.hs.length;
    // access statistics for the maintenance vessel
    const acc = run(() => operability({ series, limits: [{ key: 'hs', max: inp.ctvHs }, { key: 'wind', max: inp.ctvWind }], dtHours, requiredHours: 8, timeMs: wx.timeMs }), 'Accessibility analysis');
    const durs = [6, 12, 24, 48, 72];
    const pers = persistenceTable({ values: wx.hs, dtHours, thresholds: [1.0, 1.5, 2.0, 2.5], durations: durs, sense: 'below' });
    // installation campaign
    const ops = (inp.ops || []).filter((o) => o.hours > 0);
    need(ops.length >= 1, 'The installation table needs at least one operation with a positive duration.');
    need(ops.every((o) => o.hs > 0 && o.wind > 0), 'Every installation operation needs positive Hs and wind limits.');
    const operations = ops.map((o) => ({ name: o.name, netHours: o.hours, contingencyHours: o.hours * inp.contingency / 100, limits: [{ key: 'hs', max: o.hs, alpha: inp.alpha }, { key: 'wind', max: o.wind, alpha: inp.alpha }] }));
    const startMonths = [0, 1, 2, 3, 4, 5].map((d) => ((inp.startMonth - 1 + d) % 12) + 1);
    const mc = run(() => campaignMonteCarlo({ series, dtHours, operations, nTrials: inp.nTrials, seed: inp.seed, timeMs: wx.timeMs, startMonths }), 'Installation campaign simulation (no weather window exists for at least one operation — relax its limits or shorten it)');
    need(Number.isFinite(mc.completionHours.p50), 'No weather window was found for the installation cycle: an operation needs a longer window than the record contains within its limits.');
    const net = mc.netHours, cycle50 = mc.completionHours.p50 + inp.transit, cycle90 = mc.completionHours.p90 + inp.transit, cycleMean = mc.completionHours.mean + inp.transit;
    // campaign of N cycles: mean adds linearly; the spread of the sum of N independent cycles shrinks with √N
    const sd1 = (mc.completionHours.p90 - mc.completionHours.p50) / 1.2816;
    const camp50 = inp.nUnits * cycleMean / 24, camp90 = camp50 + 1.2816 * sd1 * Math.sqrt(inp.nUnits) / 24;
    const wow = (mc.completionHours.mean - net) / Math.max(mc.completionHours.mean, 1e-9);
    const instOper = net / Math.max(mc.completionHours.mean, 1e-9);
    // operations and maintenance
    const rows = (inp.fail || []).filter((f) => f.rate > 0 && f.hours > 0);
    need(rows.length >= 1, 'The failure table needs at least one class with a positive rate and repair time.');
    const components = rows.map((f, i) => { const v = String(f.vessel || 'ctv').trim().toLowerCase(); need(v === 'ctv' || v === 'jackup', `Failure class ${i + 1}: vessel must be 'ctv' or 'jackup'.`); return { name: f.name || `Class ${i + 1}`, failureRate: f.rate, repairHours: f.hours, vessel: v, leadTimeHours: f.lead || 0, stopsTurbine: true }; });
    const vessels = { ctv: { hsLimit: inp.ctvHs, windLimit: inp.ctvWind, count: inp.nCtv, charter: 'permanent' }, jackup: { hsLimit: inp.juHs, windLimit: inp.juWind, count: 1, charter: 'on-demand', mobilisationHours: inp.juMob * 24 } };
    const rated = c.ratedPower, pf = (u) => (u < c.cutIn || u > c.cutOut ? 0 : interp(c.u, c.p, u) / rated);
    const om = run(() => simulateWindOM({ turbines: inp.nUnits, ratedPower: rated, years: inp.omYears, weather: { hs: wx.hs, wind: wx.wind, dtHours }, components, vessels, powerFraction: pf, simulations: inp.omSims, seed: inp.seed }), 'O&M simulation');
    const availE = om.energyAvailability, availT = om.timeAvailability;
    const ctvDays = (om.vesselUtilisation.ctv ?? 0) * 365 * (vessels.ctv.count || 1), juDays = (om.vesselUtilisation.jackup ?? 0) * 365;
    const checks = [
      check('Energy-based availability against the requirement', inp.availTarget, availE * 100, '%', 'IEC TS 61400-26-1:2019 (production-based availability); project requirement', 'Demand is the required availability, capacity the simulated one. Below the requirement the energy yield and the revenue fall short of the business case.'),
      check('P90 elapsed installation campaign against the available season', camp90 / inp.nVessels, inp.seasonDays, 'days', 'DNV-ST-N001 (2021) §2 (weather-restricted operations, planning with P90 durations)', 'Campaign duration that is not exceeded in 90 % of the simulated weather sequences against the charter or season window. An overrun pushes work into winter or into the next season.'),
      check('Weather downtime share of the installation cycle', wow * 100, inp.wowLimit, '%', 'Marine-operations planning practice (project threshold on waiting on weather)', 'Fraction of the weather-restricted installation time spent waiting for a window. Above the threshold, the vessel and limits are mismatched with the site: a more capable vessel or relaxed criteria pay back.'),
    ];
    const gov = worst(checks);
    const lostMwh = om.energyLostMWhPerYear, mgoInst = camp50 * 35, vdYear = ctvDays + juDays;
    const [hx, hy, hw] = thin([Array.from({ length: n }, (_, i) => i / 24), wx.hs, wx.wind], 730);
    const samp = Array.from(mc.samples, (v) => round((v + inp.transit) / 24, 4));
    // accessibility sensitivity to the Hs limit
    const lims = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0], accL = lims.map((L) => operability({ series, limits: [{ key: 'hs', max: L }, { key: 'wind', max: inp.ctvWind }], dtHours, requiredHours: 8 }).startableFraction);
    return {
      kpis: [
        { key: 'accessWork', label: 'Accessibility: time within the access limits', value: round(acc.workableFraction * 100, 4), unit: '%', explain: 'Share of the year in which Hs and wind are within the limits of the access vessel. Central North Sea sites give 60–75 % for a 1.5 m crew-transfer limit.' },
        { key: 'accessStart', label: 'Accessibility for an 8-hour working window', value: round(acc.startableFraction * 100, 4), unit: '%', explain: 'Share of the time from which at least 8 uninterrupted hours within the limits follow: the accessibility that matters for planning a working day.' },
        { key: 'waitMean', label: 'Mean waiting time for an 8-hour window', value: round(acc.meanWaitingHours, 4), unit: 'h', explain: 'Average wait from a random moment to the next usable 8-hour access window. It adds directly to the downtime of every failure.' },
        { key: 'cycleP50', label: 'Installation cycle per position (P50)', value: round(cycle50 / 24, 4), unit: 'days', explain: `Median time to install one position including waiting on weather and transit; the net working time is ${fmt((net + inp.transit) / 24)} days. P50 here is the value not exceeded in half of the simulated weather sequences.` },
        { key: 'cycleP90', label: 'Installation cycle per position (P90)', value: round(cycle90 / 24, 4), unit: 'days', explain: 'Cycle time not exceeded in 90 % of the simulated weather sequences (non-exceedance convention: P90 is the conservative duration).' },
        { key: 'vesselDays', label: 'Installation campaign duration (mean)', value: round(camp50, 4), unit: 'days', explain: 'Expected installation-vessel days for all positions: number of positions times the mean cycle. This is the quantity passed to the cost estimate.' },
        { key: 'campP90', label: 'Elapsed installation campaign (P90)', value: round(camp90 / inp.nVessels, 4), unit: 'days', limit: inp.seasonDays, sense: 'max', explain: `Elapsed duration with ${inp.nVessels} vessel(s) in parallel that is not exceeded with 90 % probability (non-exceedance), assuming independent cycles. Plan charter and season against this value.` },
        { key: 'wow', label: 'Waiting on weather during installation', value: round(wow * 100, 4), unit: '%', limit: inp.wowLimit, sense: 'max', explain: 'Share of the weather-restricted cycle spent waiting for a window.' },
        { key: 'instOper', label: 'Installation operability (net / total time)', value: round(instOper, 4), unit: '-', explain: 'Net working time over total time on site: the weather operability factor used by schedule and cost models.' },
        { key: 'availE', label: 'Energy-based availability', value: round(availE * 100, 4), unit: '%', limit: inp.availTarget, sense: 'min', explain: 'Energy produced over the energy that would have been produced with no downtime. It is lower than time-based availability when failures cluster in windy, inaccessible periods.' },
        { key: 'availT', label: 'Time-based availability', value: round(availT * 100, 4), unit: '%', explain: 'Turbine up-time over total time (IEC TS 61400-26-1). Modern offshore farms achieve 94–97 %.' },
        { key: 'downtime', label: 'Mean downtime per failure', value: round(om.meanDowntimeHours, 4), unit: 'h', explain: 'Average time from failure to completed repair, including lead time, vessel mobilisation and waiting on weather.' },
        { key: 'failures', label: 'Failures per year in the farm', value: round(om.failuresPerYear, 4), unit: '1/yr', explain: 'Total number of failure events per year over all turbines: the workload of the maintenance organisation.' },
        { key: 'energyLost', label: 'Energy lost to downtime', value: round(lostMwh / 1e3, 5), unit: 'GWh/yr', explain: 'Production lost while turbines wait for repair. It is the yardstick for any investment in access capability.' },
        { key: 'omVesselDays', label: 'Maintenance vessel working days per year', value: round(vdYear, 4), unit: 'days/yr', explain: 'Access-vessel and jack-up days actually worked per year for corrective maintenance (scheduled service comes on top).' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Weather record: significant wave height and wind speed', x: { label: 'Day of the year', unit: 'days' }, y: { label: 'Hs (m) · wind/10 (m/s)', unit: '-' }, series: [{ name: 'Hs (m)', x: hx.map((v) => round(v, 3)), y: hy.map((v) => round(v, 3)) }, { name: 'Hub wind speed / 10', x: hx.map((v) => round(v, 3)), y: hw.map((v) => round(v / 10, 3)) }, { name: 'Access Hs limit', x: [0, n / 24], y: [inp.ctvHs, inp.ctvHs], dash: true }],
          explain: `${wx.imported ? 'Imported' : 'Synthetic'} hourly record (decimated). Access is possible while the wave trace stays under the dashed limit; the long winter exceedances are what drive downtime.` },
        { type: 'bar', title: 'Monthly accessibility for an 8-hour window', categories: acc.monthly.map((m) => MONTHS[m.month - 1]), series: [{ name: 'Accessibility (%)', y: acc.monthly.map((m) => round(m.startableFraction * 100, 3)) }], y: { label: 'Accessibility', unit: '%' },
          explain: 'Share of each month from which an 8-hour access window is available. Schedule campaigns and planned service in the high-accessibility months.' },
        { type: 'hist', title: 'Installation cycle time per position', values: samp, bins: 24, x: { label: 'Cycle time', unit: 'days' },
          explain: 'Distribution of the time to install one position over the simulated start dates. The long right tail is waiting on weather: it sets the P90 planning value.' },
        { type: 'line', title: 'Accessibility versus access-vessel wave limit', x: { label: 'Hs limit', unit: 'm' }, y: { label: 'Accessibility for 8 h', unit: '%' }, series: [{ name: 'Accessibility', x: lims, y: accL.map((v) => round(v * 100, 3)) }, { name: 'Selected limit', x: [inp.ctvHs, inp.ctvHs], y: [0, 100], dash: true }],
          explain: 'Gain in accessibility from a more capable access system. The slope at the selected limit shows what a service operation vessel with a gangway (2.5–3 m) would buy.' },
        { type: 'bar', title: 'Downtime per failure class', categories: om.byComponent.map((b) => b.name), series: [{ name: 'Mean downtime (h)', y: om.byComponent.map((b) => round(b.meanDowntimeHours, 3)) }, { name: 'Farm downtime (1000 h/yr)', y: om.byComponent.map((b) => round(b.meanDowntimeHours * b.failuresPerYear / 1e3, 4)) }], y: { label: 'Hours', unit: 'h' },
          explain: 'Mean downtime of one event and total downtime per year for each failure class. Frequent small faults and rare heavy-lift replacements usually share the loss about equally.' },
      ],
      tables: [
        { title: 'Weather windows: persistence of calm periods', columns: [{ key: 'thr', label: 'Hs threshold', unit: 'm' }, { key: 'occ', label: 'Time below', unit: '%' }, { key: 'mean', label: 'Mean window', unit: 'h' }, ...durs.map((d) => ({ key: `d${d}`, label: `Time in windows ≥ ${d} h`, unit: '%' }))],
          rows: pers.map((p) => ({ thr: p.threshold, occ: round(p.occurrence * 100, 3), mean: round(p.meanDurationHours, 3), ...Object.fromEntries(durs.map((d, i) => [`d${d}`, round(p.timeFraction[i] * 100, 3)])) })),
          explain: 'For each wave-height threshold: share of the time below it, mean length of the calm periods and share of the time lying in windows of at least the stated duration. Read off the usable time for an operation from its limit and reference period.' },
        { title: 'Installation operations', columns: [{ key: 'name', label: 'Operation' }, { key: 'tr', label: 'Reference period TR', unit: 'h' }, { key: 'start', label: 'Startable', unit: '%' }, { key: 'wait', label: 'Mean wait', unit: 'h' }, { key: 'w90', label: 'P90 wait', unit: 'h' }],
          rows: mc.operations.map((o) => ({ name: o.name, tr: round(o.requiredHours, 3), start: round(o.startableFraction * 100, 3), wait: round(o.meanWaitHours, 3), w90: round(o.p90WaitHours, 3) })),
          explain: 'Each operation with its reference period (net time plus contingency), the share of time it could start, and the waiting time before it. The operation with the longest wait is the bottleneck of the cycle.' },
        { title: 'Maintenance by failure class', columns: [{ key: 'name', label: 'Failure class' }, { key: 'n', label: 'Events', unit: '1/yr' }, { key: 'dt', label: 'Mean downtime', unit: 'h' }, { key: 'tot', label: 'Farm downtime', unit: 'h/yr' }],
          rows: om.byComponent.map((b) => ({ name: b.name, n: round(b.failuresPerYear, 4), dt: round(b.meanDowntimeHours, 4), tot: round(b.failuresPerYear * b.meanDowntimeHours, 5) })),
          explain: 'Events per year in the farm, mean downtime per event and total turbine-hours lost per year for each failure class.' },
      ],
      notes: [
        `${wx.imported ? 'Imported' : 'Synthetic'} weather record of ${fmt(n / 8760, 3)} year(s): mean Hs ${fmt(wx.hs.reduce((a, b) => a + b, 0) / n)} m, mean hub wind ${fmt(wx.wind.reduce((a, b) => a + b, 0) / n)} m/s.`,
        `Installation: ${ops.length} operations, ${fmt(net)} h net per position (TR with ${inp.contingency} % contingency, alpha ${inp.alpha}); cycle P50 ${fmt(cycle50 / 24)} days, P90 ${fmt(cycle90 / 24)} days; ${inp.nUnits} positions take ${fmt(camp50)} days on average (P90 ${fmt(camp90)}). P90 durations are non-exceedance values.`,
        `Maintenance: ${fmt(om.failuresPerYear)} failures a year, mean downtime ${fmt(om.meanDowntimeHours)} h; availability ${fmt(availT * 100, 4)} % by time, ${fmt(availE * 100, 4)} % by energy (${inp.omYears} years × ${inp.omSims} histories; standard deviation ${fmt((om.availabilityStd ?? 0) * 100, 2)} %).`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Accessibility for an 8-hour window is ${fmt(acc.startableFraction * 100)} % at an Hs limit of ${inp.ctvHs} m; energy-based availability ${fmt(availE * 100, 4)} % against ${inp.availTarget} % required; installation needs ${fmt(camp50)} vessel days (P90 ${fmt(camp90)}).`,
          gov === checks[0] && gov.uc > 1 ? `Raise availability: an access limit of ${fmt(Math.min(3, inp.ctvHs + 1))} m (service operation vessel) lifts accessibility from ${fmt(acc.startableFraction * 100)} % to ${fmt(interp(lims, accL, Math.min(3, inp.ctvHs + 1)) * 100)} %; a standing heavy-lift arrangement cuts the ${fmt(inp.juMob)}-day mobilisation that dominates major replacements.` : gov === checks[1] && gov.uc > 1 ? `The campaign overruns the season by ${fmt(camp90 / inp.nVessels - inp.seasonDays)} days at P90: start earlier, add a second vessel or relax the governing limit ('${mc.operations.reduce((a, o) => (o.meanWaitHours > a.meanWaitHours ? o : a)).name}').` : gov.uc > 1 ? 'Waiting on weather exceeds the planning threshold: choose a vessel with higher limits or split long operations into shorter weather-restricted steps.' : `Operability is adequate. Replace the ${wx.imported ? 'record' : 'synthetic weather year'} by a multi-year hindcast, and confirm limits with the vessel and lifting analyses (alpha factors per DNV-ST-N001).`,
          'DNV-ST-N001 (2021) §2 for weather-restricted operations; IEC TS 61400-26-1:2019 for availability; O&M simulation practice after Dinwoodie et al. (2015).',
          `Installation campaign fuel: ${fmt(camp50)} days ≈ ${fmt(mgoInst)} t marine gas oil ≈ ${fmt(mgoInst * CO2_MGO)} t CO₂ (screening: 35 t/day, ${CO2_MGO} t CO₂/t). Downtime costs ${fmt(lostMwh / 1e3)} GWh a year of clean energy.`),
        adv('economic', 'info', `Downtime loses ${fmt(lostMwh / 1e3, 4)} GWh a year (${fmt((1 - availE) * 100)} % of potential production); the installation campaign needs ${fmt(camp50)} vessel days, of which ${fmt(wow * 100)} % of the weather-restricted time is waiting. At the indicative rates entered: lost production ${money(lostMwh * inp.energyPrice)} a year, maintenance vessels ${money(ctvDays * inp.ctvDayRate + juDays * inp.juDayRate)} a year, installation campaign ${money(camp50 * inp.instDayRate)} at P50 and ${money(camp90 * inp.instDayRate)} at P90.`,
          `One percentage point of availability is worth ${fmt(lostMwh / Math.max((1 - availE) * 100, 1e-6) / 1e3, 3)} GWh a year (${money(lostMwh / Math.max((1 - availE) * 100, 1e-6) * inp.energyPrice)} at the indicative energy value): compare with the charter of a more capable access vessel (${fmt(vdYear)} corrective vessel days a year now). Waiting on weather in the installation campaign costs about ${money(Math.max(camp90 - camp50, 0) * inp.instDayRate)} between P50 and P90: carry it as contingency. For installation, each day saved per position saves ${inp.nUnits} vessel days.`,
          'LCOE convention (energy availability in the denominator) and marine-operations cost planning on P50 with P90 contingency (AACE RP 18R-97 Class 3–4); value of availability per IEC TS 61400-26.',
          `${vesselCo2Text(vdYear, 8)} ${energyCo2Text(lostMwh)}`),
      ],
      publish: { 'install.vesselDays': camp50, 'install.operabilityFrac': instOper, 'energy.availability': availE, 'om.vesselDaysPerYear': vdYear, 'om.energyLostFrac': 1 - availE },
      quantities: { vessel_days: { installation: round(camp50, 3), maintenance_per_year: round(vdYear, 3) } },
    };
  },
  verify: [
    { name: 'Accessibility of a Weibull-distributed sea state: P(Hs ≤ limit) = 1 − exp(−(limit/scale)^k)', source: 'Weibull cumulative distribution; scale = mean / Γ(1 + 1/k) (hand calculation)',
      // mean 1.6 m, k = 1.5: Γ(1.6667) = 0.90275 → scale 1.7724 m; P(Hs ≤ 1.5) = 1 − exp(−(0.84633)^1.5) = 1 − exp(−0.77860) = 54.10 % (no seasonality, wind limit not governing)
      inputs: { hsMean: 1.6, hsShape: 1.5, season: 0, ctvHs: 1.5, ctvWind: 30, omYears: 1, omSims: 1, nTrials: 20 }, expect: [{ kpi: 'accessWork', value: 54.10, tol: 0.003 }] },
    { name: 'Repairable turbines with unrestricted access: availability = 1 / (1 + λ·MTTR)', source: 'Steady-state availability of a two-state repairable unit (Rausand & Høyland, System Reliability Theory, §9.3)',
      // λ = 2 per operating year, repair 100 h, no lead time, limits never exceeded: A = 1/(1 + 2 × 100/8760) = 97.77 %
      inputs: { fail: [{ name: 'Test', rate: 2, hours: 100, vessel: 'ctv', lead: 0 }], ctvHs: 4, ctvWind: 30, nCtv: 20, hsMean: 0.5, season: 0, nUnits: 60, omYears: 5, omSims: 3, nTrials: 20 }, expect: [{ kpi: 'availT', value: 97.77, abs: 0.3 }] },
  ],
};
