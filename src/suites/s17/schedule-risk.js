// Installation and construction schedule: critical path, PERT, weather-adjusted Monte Carlo schedule risk,
// cost-loaded S-curve and the cost of delay (standby, deferred revenue, liquidated damages).
import { cpm, pertNetwork, scheduleRisk, costLoading, delayCost, liquidatedDamages, weatherAdjustedDuration } from '../../engine/econ/schedule.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, PCT_COST, CO2_MGO, MGO_T_PER_DAY, need, run, round, fmt, mil, fM, pct, adv, check, cost, tableRows, uniqueNames, seriesXY, percentiles, thinSample, range } from './_common.js';

const ACT = [
  { id: 'ENG', name: 'Detailed engineering', predecessors: '', a: 240, m: 300, b: 400, weather: 0, cost: 9.0e7 },
  { id: 'PROC', name: 'Procurement of long-lead items', predecessors: 'ENG', a: 300, m: 365, b: 520, weather: 0, cost: 7.2e8 },
  { id: 'FAB', name: 'Foundation fabrication', predecessors: 'ENG', a: 330, m: 400, b: 540, weather: 0, cost: 2.6e8 },
  { id: 'ONS', name: 'Onshore substation and grid works', predecessors: 'ENG', a: 380, m: 450, b: 600, weather: 0, cost: 8.0e7 },
  { id: 'TRN', name: 'Transport and marshalling', predecessors: 'FAB', a: 30, m: 45, b: 80, weather: 1, cost: 2.5e7 },
  { id: 'FND', name: 'Foundation installation', predecessors: 'TRN', a: 70, m: 95, b: 150, weather: 1, cost: 9.5e7 },
  { id: 'CBL', name: 'Cable lay and burial', predecessors: 'FND', a: 60, m: 80, b: 140, weather: 1, cost: 7.0e7 },
  { id: 'OSS', name: 'Offshore substation installation', predecessors: 'FND,PROC', a: 12, m: 18, b: 40, weather: 1, cost: 3.0e7 },
  { id: 'WTG', name: 'Turbine installation', predecessors: 'FND,PROC', a: 80, m: 105, b: 170, weather: 1, cost: 8.5e7 },
  { id: 'COM', name: 'Commissioning and energisation', predecessors: 'WTG,CBL,OSS,ONS', a: 45, m: 60, b: 110, weather: 0, cost: 2.5e7 },
];
const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

export default {
  id: 's17-schedule-risk', suite: 17, title: 'Installation schedule, schedule risk and cost of delay',
  summary: 'Analyses the fabrication-to-commissioning network for its critical path, simulates activity durations with seasonal weather operability to give the completion-date distribution, loads cost onto the schedule and prices delay.',
  covers: [2, 7], fidelity: 'numerical',
  standards: ['AACE International RP 57R-09 (integrated cost and schedule risk analysis using Monte Carlo simulation of a CPM model)', 'AACE International RP 49R-06 (identifying the critical path)', 'AACE International RP 55R-09 (S-curves)', 'DNV-ST-N001 (2021) §3 (weather-restricted operations and operability)', 'FIDIC Conditions of Contract (2017) Sub-Clause 8.8 (delay damages)', 'SCL Delay and Disruption Protocol, 2nd ed. (2017)'],
  models: ['Critical-path method', 'PERT three-point estimates and distributions', 'Monte Carlo schedule risk with merge bias', 'Seasonal weather operability stretch of marine activities', 'Criticality index', 'Cost loading and expenditure S-curve', 'Delay cost: standby, deferred revenue and liquidated damages', 'Probabilistic schedule distribution'],
  limits: `Activity-on-node network with finish-to-start links and beta-PERT durations sampled independently; weather acts through monthly mean operability, so persistence of individual weather windows and vessel availability clashes are not resolved (use a time-domain weather-window simulation for campaign planning). ${CUR_NOTE}`,
  publishes: ['econ.scheduleP50', 'econ.scheduleP80', 'econ.delayCost', 'econ.marineGrossDays'],
  searchTerms: 'schedule risk analysis critical path offshore installation weather downtime',
  mesh: { key: 'nSamples', label: 'Monte Carlo samples', levels: [250, 500, 1000, 2000], qoi: 'mean' },
  inputs: [
    { key: 'startDay', label: 'Project start (day of year)', unit: 'day of year', type: 'int', default: 0, min: 0, max: 365, group: 'Opening state', role: 'ic', help: 'Commencement date: 0 = 1 January. Sets which season each offshore activity falls in.' },
    { key: 'activities', label: 'Activity network', type: 'table', group: 'Installation schedule', role: 'operation',
      help: `Durations in days as optimistic / most likely / pessimistic net working time. Predecessors: comma-separated ids (finish-to-start). Weather = 1 for offshore work stretched by operability. Cost in ${CUR}. ${PLACEHOLDER}`,
      columns: [{ key: 'id', label: 'Id', type: 'text' }, { key: 'name', label: 'Activity', type: 'text' }, { key: 'predecessors', label: 'Predecessors', type: 'text' }, { key: 'a', label: 'Optimistic', unit: 'days', type: 'number' }, { key: 'm', label: 'Most likely', unit: 'days', type: 'number' }, { key: 'b', label: 'Pessimistic', unit: 'days', type: 'number' },
        { key: 'weather', label: 'Weather sensitive (0/1)', type: 'number' }, { key: 'cost', label: 'Cost', unit: CUR, type: 'number' }], default: ACT },
    { key: 'operabilityShape', label: 'Monthly operability', type: 'series', group: 'Weather', role: 'environment', xLabel: 'Month', yLabel: 'Operability', xUnit: '-', yUnit: '-',
      help: 'Fraction of time the limiting marine operation is workable in each month (1 = January), from the installation analyses.',
      default: { x: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], y: [0.42, 0.46, 0.55, 0.68, 0.80, 0.87, 0.90, 0.88, 0.76, 0.62, 0.50, 0.43] } },
    { key: 'meanOperability', label: 'Annual mean operability', unit: '-', type: 'number', default: 0.66, min: 0.05, max: 1, group: 'Weather', role: 'environment', bind: 'install.operabilityFrac', uncertain: { dist: 'normal', std: 0.06 },
      help: 'The monthly shape is scaled so that its mean equals this value; bound to the operability from the installation suites.' },
    { key: 'target', label: 'Required completion (from start)', unit: 'days (elapsed)', type: 'number', default: 1280, min: 1, max: 10000, group: 'Constraints', role: 'bc', help: 'Contractual or licence completion date.' },
    { key: 'confidence', label: 'Confidence level for the schedule', unit: '%', type: 'number', default: 80, min: 50, max: 99, group: 'Constraints', role: 'criteria' },
    cost({ key: 'standbyRate', label: 'Standby cost of the marine spread and site team', unit: `${CUR}/day`, default: 4.2e5, min: 0, max: 1e8, group: 'Cost of delay' }),
    cost({ key: 'revenuePerDay', label: 'Net revenue deferred per day of late start-up', unit: `${CUR}/day`, default: 4.6e5, min: 0, max: 1e9, group: 'Cost of delay', help: 'Net operating cash flow per day once in operation.' }),
    cost({ key: 'ldRate', label: 'Liquidated damages payable', unit: `${CUR}/day`, default: 2.5e5, min: 0, max: 1e8, group: 'Cost of delay' }),
    { key: 'ldCap', label: 'Cap on liquidated damages', unit: '% of contract value', type: 'number', default: 10, min: 0, max: 100, group: 'Cost of delay', role: 'bc' },
    cost({ key: 'contractValue', label: 'Contract value subject to damages', unit: CUR, default: 4.0e8, min: 0, max: 1e12, group: 'Cost of delay' }),
    { key: 'nSamples', label: 'Monte Carlo samples', unit: '-', type: 'int', default: 1000, min: 100, max: 20000, group: 'Analysis control', role: 'control' },
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 17, min: 1, max: 1e6, group: 'Analysis control', role: 'control' },
  ],
  compute(inp) {
    const rows = tableRows(inp.activities, 'Activity network', { a: { min: 0 }, m: { min: 0 }, b: { min: 0 }, weather: { min: 0, max: 1 }, cost: { min: 0 } });
    uniqueNames(rows, 'Activity network', 'id');
    const sh = seriesXY(inp.operabilityShape, 'Monthly operability', 12);
    need(sh.y.length === 12, 'Monthly operability needs exactly 12 values (January to December).');
    need(sh.y.every((v) => v > 0 && v <= 1), 'Monthly operability values must lie in (0, 1].');
    const avg = sh.y.reduce((a, c) => a + c, 0) / 12;
    const op = sh.y.map((v) => Math.min(1, Math.max(0.05, v * inp.meanOperability / avg)));
    const ids = new Set(rows.map((r) => String(r.id).trim()));
    const acts = rows.map((r) => {
      need(r.a <= r.m && r.m <= r.b, `Activity '${r.id}': durations must satisfy optimistic ≤ most likely ≤ pessimistic.`);
      const preds = String(r.predecessors ?? '').split(/[,;\s]+/).map((x) => x.trim()).filter(Boolean);
      for (const pr of preds) need(ids.has(pr), `Activity '${r.id}': predecessor '${pr}' is not in the table.`);
      return { id: String(r.id).trim(), name: r.name || r.id, predecessors: preds, a: r.a, m: r.m, b: r.b, weather: r.weather >= 0.5, cost: r.cost };
    });
    const n = Math.round(inp.nSamples);
    // deterministic network on the most-likely net durations, and the same with seasonal weather stretch
    const det = run(() => cpm(acts.map((x) => ({ id: x.id, duration: x.m, predecessors: x.predecessors }))), 'Critical path (check for dependency cycles)');
    const pert = run(() => pertNetwork(acts.map((x) => ({ id: x.id, a: x.a, m: x.m, b: x.b, predecessors: x.predecessors })), inp.target), 'PERT');
    const mc = run(() => scheduleRisk({ activities: acts.map((x) => (x.b > x.a ? { id: x.id, predecessors: x.predecessors, a: x.a, m: x.m, b: x.b, ...(x.weather ? { operabilityByMonth: op } : {}) } : { id: x.id, predecessors: x.predecessors, duration: x.m, ...(x.weather ? { operabilityByMonth: op } : {}) })), n, seed: inp.seed, projectStartDay: inp.startDay }), 'Schedule risk');
    const P = percentiles(mc.durations);
    const pConf = P(inp.confidence);
    let ok = 0; for (let k = 0; k < n; k++) if (mc.durations[k] <= inp.target) ok++;
    const pMeet = ok / n;
    // weather stretch of each marine activity at its mean start date (for the table and the vessel-day total)
    const sched = det.activities;
    const info = acts.map((x) => {
      const s = sched.find((q) => q.id === x.id);
      const start = mc.meanStart[x.id] ?? s.ES; const fin = mc.meanFinish[x.id] ?? s.EF;
      const wx = x.weather ? weatherAdjustedDuration({ netDays: x.m, operabilityByMonth: op, startDay: (inp.startDay + start) % 365 }) : { grossDays: x.m, weatherDays: 0 };
      return { ...x, ES: s.ES, EF: s.EF, float: s.totalFloat, critical: s.critical, ci: mc.criticalityIndex[x.id] ?? 0, meanStart: start, meanFinish: fin, gross: wx.grossDays, weatherDays: wx.weatherDays };
    });
    const marineNet = info.filter((x) => x.weather).reduce((s, x) => s + x.m, 0); const marineGross = info.filter((x) => x.weather).reduce((s, x) => s + x.gross, 0);
    // cost loading on the mean simulated dates
    const total = acts.reduce((s, x) => s + x.cost, 0);
    const simSched = { duration: mc.mean, activities: info.map((x) => ({ id: x.id, duration: Math.max(x.meanFinish - x.meanStart, 1e-6), ES: x.meanStart, EF: x.meanFinish, LS: x.meanStart, LF: x.meanFinish })) };
    const load = total > 0 ? run(() => costLoading({ schedule: simSched, costs: Object.fromEntries(acts.map((x) => [x.id, x.cost])), periodDays: 365.25 / 4 }), 'Cost loading') : null;
    const annual = total > 0 ? costLoading({ schedule: simSched, costs: Object.fromEntries(acts.map((x) => [x.id, x.cost])), periodDays: 365.25 }) : null;
    // delay beyond the required completion, at the chosen confidence and in expectation
    const delayConf = Math.max(pConf - inp.target, 0);
    let dsum = 0; for (let k = 0; k < n; k++) dsum += Math.max(mc.durations[k] - inp.target, 0);
    const delayExp = dsum / n;
    const dc = delayCost({ delayDays: delayConf, standbyDayRates: [inp.standbyRate], revenuePerDay: inp.revenuePerDay, mode: 'lost' });
    const ld = liquidatedDamages({ delayDays: delayConf, ratePerDay: inp.ldRate, contractValue: inp.contractValue, capFraction: inp.ldCap / 100 });
    const dcExp = (inp.standbyRate + inp.revenuePerDay) * delayExp;
    const critNames = det.criticalPath.join(' → ');

    const kpis = [
      { key: 'detDuration', label: 'Deterministic duration (most-likely, no weather)', value: round(det.duration), unit: 'days', explain: 'Critical-path length with every activity at its most-likely net duration and no weather downtime: the schedule usually shown on a plan, and an optimistic one.' },
      { key: 'pertMean', label: 'PERT expected duration', value: round(pert.expectedDuration), unit: 'days', explain: `Critical-path length on the PERT means (a + 4m + b)/6, standard deviation ${fmt(pert.std)} days. It ignores near-critical paths and weather.` },
      { key: 'mean', label: 'Simulated mean duration', value: round(mc.mean), unit: 'days', explain: 'Mean completion time from the Monte Carlo simulation with merge bias and seasonal weather. The excess over the deterministic figure is the schedule risk hidden in the plan.' },
      { key: 'p10', label: 'P10 duration (cumulative)', value: round(P(10)), unit: 'days', explain: 'Duration not exceeded with 10 % probability: the optimistic end.' },
      { key: 'p50', label: 'P50 duration (cumulative)', value: round(P(50)), unit: 'days', explain: 'Duration with an even chance of being exceeded.' },
      { key: 'p80', label: 'P80 duration (cumulative)', value: round(P(80)), unit: 'days', explain: 'Duration not exceeded with 80 % probability; a common basis for commitments to third parties.' },
      { key: 'p90', label: 'P90 duration (cumulative)', value: round(P(90)), unit: 'days', explain: 'Duration not exceeded with 90 % probability: the conservative, late case.' },
      { key: 'pConf', label: `Duration at P${inp.confidence}`, value: round(pConf), unit: 'days', limit: inp.target, sense: 'max', explain: 'Completion time at the chosen confidence level, compared with the required completion date.' },
      { key: 'pMeet', label: 'Probability of meeting the required date', value: round(100 * pMeet), unit: '%', limit: inp.confidence, sense: 'min', explain: 'Share of simulated outcomes finishing by the required completion date. It should be at least the chosen confidence level.' },
      { key: 'schedCont', label: 'Schedule contingency to the chosen confidence', value: round(pConf - det.duration), unit: 'days', explain: 'Days to add to the deterministic plan to reach the chosen confidence: the time allowance equivalent to cost contingency.' },
      { key: 'weatherDays', label: 'Expected waiting on weather (marine activities)', value: round(marineGross - marineNet), unit: 'days', explain: 'Gross less net duration of the weather-sensitive activities at their expected season. These are vessel days paid without progress.' },
      { key: 'delayCost', label: `Cost of delay at P${inp.confidence}`, value: mil(dc.total + ld.damages), unit: MU, explain: `Standby (${fM(dc.standby)}), revenue lost while not operating (${fM(dc.revenueImpact)}) and liquidated damages (${fM(ld.damages)}${ld.capped ? ', capped' : ''}) for the days beyond the required date at the chosen confidence.` },
      { key: 'delayCostExp', label: 'Expected cost of delay', value: mil(dcExp), unit: MU, explain: 'Probability-weighted standby and lost revenue over all simulated outcomes: the amount to carry in the cost risk analysis for schedule.' },
    ];
    const checks = [check(`P${inp.confidence} completion within the required date`, pConf, inp.target, 'days', 'AACE RP 57R-09 (schedule confidence level)', 'Guards against committing to a completion date the network can meet only with luck; above 1.0 the date is not achievable at the chosen confidence and float or resources must be added.')];
    const cdfN = 50;
    const plots = [
      { type: 'hist', title: 'Distribution of completion time', values: thinSample(mc.durations, 1500).map((v) => round(v)), bins: 36, x: { label: 'Duration from start', unit: 'days' }, y: { label: 'Frequency' }, explain: 'Simulated completion times. The right tail comes from offshore work slipping into the winter season; compare its extent with the required date.' },
      { type: 'line', title: 'Cumulative probability of completion', x: { label: 'Duration from start', unit: 'days' }, y: { label: 'Probability of completing by then', unit: '%' },
        series: [{ name: 'Completion', x: range(cdfN + 1, (k) => round(P(100 * k / cdfN))), y: range(cdfN + 1, (k) => round(100 * k / cdfN)) }, { name: 'Required date', x: [inp.target, inp.target], y: [0, 100], dash: true }, { name: 'Deterministic plan', x: [det.duration, det.duration], y: [0, 100], dash: true }],
        explain: 'Read up from a date to the chance of finishing by it. The height at which the deterministic plan cuts the curve is the real confidence of that plan — usually low.' },
      { type: 'bar', title: 'Criticality index by activity', categories: info.map((x) => x.id), x: { label: 'Activity' }, y: { label: 'On the critical path in', unit: '% of runs' }, series: [{ name: 'Criticality', y: info.map((x) => round(100 * x.ci, 4)) }],
        explain: 'How often each activity lies on the critical path across the simulation. Activities with a high index but float in the deterministic plan are the hidden schedule drivers.' },
      { type: 'bar', title: 'Monthly operability used', categories: MONTHS, x: { label: 'Month' }, y: { label: 'Operability', unit: '-' }, series: [{ name: 'Operability', y: op.map((v) => round(v, 4)) }], explain: 'Workable fraction of each month after scaling to the annual mean. Marine activities falling in the low months are stretched most.' },
    ];
    if (load) plots.push({ type: 'line', title: 'Expenditure S-curve', x: { label: 'Quarter from start', unit: '-' }, y: { label: 'Cumulative spend', unit: MU }, series: [{ name: 'Cumulative cost', x: load.periods, y: load.cumulative.map(mil) }],
      explain: 'Cumulative cost loaded on the mean simulated dates. The steep middle section is peak cash demand; the annual fractions feed the capex phasing of the cash-flow models.' });
    const tables = [
      { title: 'Activity schedule', columns: [{ key: 'id', label: 'Id' }, { key: 'name', label: 'Activity' }, { key: 'm', label: 'Net duration', unit: 'days' }, { key: 'gross', label: 'With weather', unit: 'days' }, { key: 'ES', label: 'Early start', unit: 'day' }, { key: 'float', label: 'Total float', unit: 'days' }, { key: 'crit', label: 'Critical (plan)' }, { key: 'ci', label: 'Criticality index', unit: '%' }, { key: 'ms', label: 'Mean simulated start', unit: 'day' }, { key: 'mf', label: 'Mean simulated finish', unit: 'day' }],
        rows: info.map((x) => ({ id: x.id, name: x.name, m: x.m, gross: round(x.gross, 4), ES: round(x.ES, 4), float: round(x.float, 4), crit: x.critical ? 'yes' : 'no', ci: round(100 * x.ci, 3), ms: round(x.meanStart, 4), mf: round(x.meanFinish, 4) })),
        explain: 'Each activity with its planned dates and float, its weather-stretched duration and how often it is critical in the simulation.' },
      { title: 'Schedule percentiles (cumulative)', columns: [{ key: 'p', label: 'Percentile' }, { key: 'd', label: 'Duration', unit: 'days' }, { key: 'slip', label: 'Beyond required date', unit: 'days' }, { key: 'c', label: 'Standby and lost revenue', unit: MU }],
        rows: [10, 20, 30, 40, 50, 60, 70, 80, 90, 95].map((q) => { const d = P(q); const sl = Math.max(d - inp.target, 0); return { p: `P${q}`, d: round(d, 5), slip: round(sl, 4), c: mil(sl * (inp.standbyRate + inp.revenuePerDay)) }; }),
        explain: 'Completion time and the cost of the resulting slip at each confidence level (Pxx = not exceeded with xx % probability).' },
    ];
    if (annual) tables.push({ title: 'Annual capex phasing', columns: [{ key: 'y', label: 'Year' }, { key: 'c', label: 'Spend', unit: MU }, { key: 'f', label: 'Fraction', unit: '%' }], rows: annual.periods.map((y, k) => ({ y, c: mil(annual.cost[k]), f: round(100 * annual.fractions[k], 4) })),
      explain: 'Share of the activity costs falling in each year of the programme: the phased capital spending used by the cash-flow models.' });
    const hidden = info.filter((x) => !x.critical && x.ci > 0.2).sort((a, b) => b.ci - a.ci)[0];
    const topCi = info.slice().sort((a, b) => b.ci - a.ci)[0];
    const advice = [
      adv('engineering', pMeet * 100 >= inp.confidence ? 'info' : 'warn', `Critical path of the plan: ${critNames} (${fmt(det.duration)} days). Simulated mean ${fmt(mc.mean)} days, P${inp.confidence} ${fmt(pConf)} days; "${topCi.name}" is critical in ${pct(topCi.ci)} of runs${hidden ? `, and "${hidden.name}" in ${pct(hidden.ci)} although the plan shows it with ${fmt(hidden.float)} days of float` : ''}. Marine work waits about ${fmt(marineGross - marineNet)} days on weather.`,
        `Move the weather-sensitive activities into the high-operability months (start date, earlier fabrication release or a second spread), and protect "${topCi.name}". Raising the limiting sea state of the marine operations lifts operability directly: each 0.05 on the annual mean shortens the marine work by roughly ${fmt(marineNet / inp.meanOperability - marineNet / Math.min(1, inp.meanOperability + 0.05))} days.`,
        'Schedule risk analysis on a CPM model, AACE RP 57R-09; weather-restricted operations and operability, DNV-ST-N001 (2021).',
        `Waiting on weather burns fuel without progress: ${fmt(marineGross - marineNet)} standby days at about ${MGO_T_PER_DAY * 0.4} t/day in standby mode is roughly ${fmt((marineGross - marineNet) * MGO_T_PER_DAY * 0.4 * CO2_MGO)} t CO₂ (screening: ${CO2_MGO} t CO₂ per t marine gas oil). Better seasonal placement cuts both cost and emissions.`),
      adv('economic', pConf <= inp.target ? 'ok' : 'warn', `Probability of meeting the required date of ${inp.target} days is ${pct(pMeet)}. At P${inp.confidence} the slip is ${fmt(delayConf)} days, costing ${fM(dc.total + ld.damages)} (standby ${fM(dc.standby)}, lost revenue ${fM(dc.revenueImpact)}, damages ${fM(ld.damages)}); expected delay cost ${fM(dcExp)}.`,
        pConf <= inp.target ? `The date is achievable at the chosen confidence with ${fmt(inp.target - pConf)} days in hand. Carry the expected delay cost of ${fM(dcExp)} in the cost risk register.` : `Add ${fmt(pConf - inp.target)} days of float, or accept and price the risk: one day of slip costs ${fM(inp.standbyRate + inp.revenuePerDay)}, so spending up to that per day saved on acceleration is justified.`,
        'Integrated cost–schedule risk, AACE RP 57R-09; delay costing per the SCL Delay and Disruption Protocol (2017); delay damages, FIDIC (2017) Sub-Clause 8.8.'),
    ];
    return { kpis, checks, plots, tables, advice,
      notes: [PCT_COST.replace('For a cost, P90 is the high-cost (conservative) case.', 'For a duration, P90 is the late (conservative) case.'), `${n} samples, seed ${inp.seed}; durations are beta-PERT between the optimistic and pessimistic values; weather-sensitive activities are stretched by the monthly operability from their sampled start date.`, `Deterministic plan ${fmt(det.duration)} days; PERT ${fmt(pert.expectedDuration)} ± ${fmt(pert.std)} days; simulation ${fmt(mc.mean)} ± ${fmt(mc.std)} days.`],
      publish: { 'econ.scheduleP50': P(50), 'econ.scheduleP80': P(80), 'econ.delayCost': dcExp, 'econ.marineGrossDays': marineGross },
      quantities: { vessel_days: { marine_gross: round(marineGross) } } };
  },
  verify: [
    // A (3) → B (5) and C (8) in parallel → D (4): critical path A–C–D = 15 days; B has 3 days of float. Fixed durations, no weather.
    { name: 'Critical path of a four-activity network', source: 'Hand forward pass: EF(A) = 3, EF(B) = 8, EF(C) = 11, ES(D) = 11, EF(D) = 15',
      inputs: { activities: [{ id: 'A', name: 'A', predecessors: '', a: 3, m: 3, b: 3, weather: 0, cost: 0 }, { id: 'B', name: 'B', predecessors: 'A', a: 5, m: 5, b: 5, weather: 0, cost: 0 }, { id: 'C', name: 'C', predecessors: 'A', a: 8, m: 8, b: 8, weather: 0, cost: 0 }, { id: 'D', name: 'D', predecessors: 'B,C', a: 4, m: 4, b: 4, weather: 0, cost: 0 }], target: 20, nSamples: 200 },
      expect: [{ kpi: 'detDuration', value: 15, tol: 1e-9 }, { kpi: 'mean', value: 15, tol: 1e-9 }, { kpi: 'p90', value: 15, tol: 1e-9 }, { kpi: 'pMeet', value: 100, tol: 1e-9 }] },
    // single activity a = 2, m = 5, b = 14: PERT mean (2 + 20 + 14)/6 = 6.0; the beta-PERT distribution has the same mean, so the simulation mean → 6.0
    { name: 'PERT mean of one activity (2, 5, 14)', source: 'Malcolm et al. (1959): (a + 4m + b)/6 = 6.0 days',
      inputs: { activities: [{ id: 'A', name: 'A', predecessors: '', a: 2, m: 5, b: 14, weather: 0, cost: 0 }], target: 20, nSamples: 4000 }, expect: [{ kpi: 'pertMean', value: 6, tol: 1e-9 }, { kpi: 'mean', value: 6, tol: 0.02 }] },
    // one weather-sensitive activity of 30 net days at a flat operability of 0.5 → 60 gross days
    { name: 'Weather stretch at constant operability 0.5', source: 'Hand calculation: gross = net / operability = 30 / 0.5 = 60 days',
      inputs: { activities: [{ id: 'A', name: 'A', predecessors: '', a: 30, m: 30, b: 30, weather: 1, cost: 0 }], operabilityShape: { x: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], y: new Array(12).fill(0.5) }, meanOperability: 0.5, target: 100, nSamples: 200 },
      expect: [{ kpi: 'mean', value: 60, tol: 1e-6 }, { kpi: 'weatherDays', value: 30, tol: 1e-6 }, { kpi: 'detDuration', value: 30, tol: 1e-9 }] },
  ],
};
