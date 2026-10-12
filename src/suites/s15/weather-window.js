// Installation weather windows and operability from a metocean time series: workable and startable fractions,
// waiting times, persistence, a Markov-chain sea-state model and a Monte Carlo campaign duration.
import { operationReferencePeriod, forecastLimit, operability, persistenceTable, markovFit, markovSimulate, campaignMonteCarlo } from '../../engine/metocean/windows.js';
import { fitWeibull } from '../../engine/reliability/extreme.js';
import { need, run, round, fmt, kpi, check, adv, line, bar, tbl, vesselCo2Text, syntheticWeather, CONFIRM, N001 } from './_common.js';

const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
const SEASONS = { all: null, summer: [4, 5, 6, 7, 8, 9], winter: [10, 11, 12, 1, 2, 3] };

export default {
  id: 's15-weather-window', suite: 15, title: 'Installation weather windows and operability',
  summary: 'Finds how often and for how long the operation limits on wave height, period and wind are met for the full reference period, the waiting time for a window by month, and the P10/P50/P90 duration of a campaign of repeated operations.',
  covers: [15], fidelity: 'numerical',
  standards: [`${N001} Sec. 2 (operation reference period, weather-restricted operations, alpha-factor, forecast limits)`, 'DNV-RP-N103 (2017) §9 (operational criteria and weather windows)', 'ISO 19901-1 (2015) §5–6 (metocean data for operations, persistence)', 'Graham (1982) persistence statistics, Coastal Eng. 6; Anastasiou & Tsekos (1996) Markov persistence, Applied Ocean Research 18'],
  models: ['Wave scatter / time-series statistics', 'Threshold and duration (persistence) logic', 'Markov-chain weather model (first order, fitted and simulated)', 'Monte Carlo simulation of sequential operations', 'Weibull distribution of significant wave height', 'Weather downtime model', 'Probabilistic installation schedule (P10/P50/P90)'],
  limits: 'Results are only as good as the series: without a pasted series or site hindcast statistics a seeded synthetic series with realistic persistence and seasonality is used, which is indicative only. Alpha-factors are recalled DNV-ST-N001 practice and must be read from the governing table for the planned duration and forecast level. The campaign simulation waits for a window covering each whole operation and cannot suspend an operation once started; technical downtime, mobilisation, crew change and learning effects are not modelled, and a first-order Markov chain under-predicts long calms. Extreme-value (generalised extreme value) analysis of design storms belongs to the metocean suite.',
  publishes: ['install.operabilityFrac', 'install.vesselDays', 'install.waitHoursMean'],
  searchTerms: 'weather window operability persistence waiting on weather Monte Carlo campaign duration alpha factor',
  inputs: [
    { key: 'series', label: 'Metocean time series (optional)', type: 'table', group: 'Metocean series', role: 'environment',
      columns: [{ key: 'hs', label: 'Hs', unit: 'm', type: 'number' }, { key: 'tp', label: 'Tp', unit: 's', type: 'number' }, { key: 'wind', label: 'Wind speed', unit: 'm/s', type: 'number' }],
      default: [], help: 'Paste a hindcast or measured series at a constant interval, starting on 1 January, ideally whole years. Leave empty to use the site hindcast statistics when the project has them, otherwise a seeded synthetic series.' },
    { key: 'dtHours', label: 'Sampling interval of the series', unit: 'h', type: 'number', default: 3, min: 0.25, max: 24, group: 'Metocean series', role: 'environment' },
    { key: 'hsMean', label: 'Synthetic series: annual mean Hs', unit: 'm', type: 'number', default: 1.8, min: 0.2, max: 6, group: 'Synthetic climate', role: 'environment', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Used only when no series is given and the site record has no monthly wave statistics. North Sea about 2 m, West Africa about 1.2 m, Gulf of Mexico about 1 m.' },
    { key: 'season', label: 'Synthetic series: seasonal variation of the mean', unit: '-', type: 'number', default: 0.35, min: 0, max: 0.8, group: 'Synthetic climate', role: 'environment', help: 'Winter mean = annual mean × (1 + this), summer mean = × (1 − this).' },
    { key: 'tau', label: 'Synthetic series: persistence time', unit: 'h', type: 'number', default: 36, min: 6, max: 240, group: 'Synthetic climate', role: 'environment', help: 'Decorrelation time of the sea state; 24–48 h for storm-dominated seas, longer for swell.' },
    { key: 'years', label: 'Synthetic series: length', unit: 'years', type: 'int', default: 3, min: 1, max: 20, group: 'Synthetic climate', role: 'control' },
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 7, min: 1, max: 99999, group: 'Synthetic climate', role: 'control' },
    { key: 'hsLim', label: 'Design limit on Hs (OPLIM)', unit: 'm', type: 'number', default: 2.0, min: 0.2, max: 10, group: 'Operational limits', role: 'bc', bind: 'install.splashHsLimit', design: true, help: 'From the operation analysis (lift, splash zone, mating, lay). The forecast must show less: OPWF = α × OPLIM.' },
    { key: 'alphaHs', label: 'Alpha-factor for waves', unit: '-', type: 'number', default: 0.76, min: 0.3, max: 1, group: 'Operational limits', role: 'criteria', help: `Depends on planned duration, Hs level and forecast/monitoring level: roughly 0.65–0.80 for 12–24 h operations with standard forecasting, up to 1.0 with on-site monitoring and a meteorologist. ${CONFIRM}` },
    { key: 'windLim', label: 'Design limit on wind speed', unit: 'm/s', type: 'number', default: 12, min: 2, max: 40, group: 'Operational limits', role: 'bc' },
    { key: 'alphaWind', label: 'Alpha-factor for wind', unit: '-', type: 'number', default: 0.8, min: 0.3, max: 1, group: 'Operational limits', role: 'criteria', help: CONFIRM },
    { key: 'tpMax', label: 'Upper limit on Tp (0 = none)', unit: 's', type: 'number', default: 0, min: 0, max: 25, group: 'Operational limits', role: 'bc', help: 'Swell-sensitive operations (float-over, heavy lift from a barge) often carry a period limit.' },
    { key: 'tpop', label: 'Planned operation period TPOP', unit: 'h', type: 'number', default: 18, min: 0.5, max: 400, group: 'Operation', role: 'operation', design: true },
    { key: 'tc', label: 'Contingency time TC', unit: 'h', type: 'number', default: 9, min: 0, max: 400, group: 'Operation', role: 'criteria', help: `Normally at least 50 % of TPOP (and not less than 6 h) unless justified. ${CONFIRM}` },
    { key: 'nOps', label: 'Number of operations in the campaign', unit: '-', type: 'int', default: 12, min: 1, max: 500, group: 'Campaign', role: 'operation', help: 'Repetitions of the same weather-restricted operation, e.g. foundations or lifts, executed one after another.' },
    { key: 'gapHours', label: 'Unrestricted work between operations', unit: 'h', type: 'number', default: 12, min: 0, max: 400, group: 'Campaign', role: 'operation', help: 'Transit, preparation and seafastening removal that can proceed in any weather; added to the campaign duration.' },
    { key: 'startSeason', label: 'Campaign start', type: 'select', default: 'all', group: 'Campaign', role: 'ic', options: [{ value: 'all', label: 'Any time of year' }, { value: 'summer', label: 'April–September' }, { value: 'winter', label: 'October–March' }], help: 'Initial condition of the schedule simulation: when the spread arrives on site.' },
    { key: 'nTrials', label: 'Monte Carlo trials', unit: '-', type: 'int', default: 500, min: 50, max: 20000, group: 'Campaign', role: 'control' },
    { key: 'operReq', label: 'Target operability', unit: '%', type: 'number', default: 50, min: 1, max: 100, group: 'Operation', role: 'criteria', help: 'Project planning threshold for the share of time at which the operation could start.' },
    { key: 'dayRate', label: 'Indicative spread day rate (0 = quantities only)', unit: 'USD', type: 'number', default: 0, min: 0, max: 5e6, group: 'Campaign', role: 'cost', help: 'Optional and indicative: used only to express waiting on weather as money in the advice.' },
  ],
  compute(inp, ctx) {
    const rows = (inp.series || []).filter((r) => Number.isFinite(Number(r.hs)) && r.hs !== '' && r.hs !== null);
    let ser, source;
    const dt = inp.dtHours;
    if (rows.length >= 8) {
      const hs = Float64Array.from(rows, (r) => Number(r.hs)), tp = Float64Array.from(rows, (r) => (Number.isFinite(Number(r.tp)) ? Number(r.tp) : 0)), wind = Float64Array.from(rows, (r) => (Number.isFinite(Number(r.wind)) ? Number(r.wind) : 0));
      const t0 = Date.UTC(2001, 0, 1);
      ser = { hs, tp, wind, timeMs: Float64Array.from(rows, (r, i) => t0 + i * dt * 3.6e6), dtHours: dt }; source = `series supplied (${rows.length} samples at ${dt} h)`;
    } else {
      need(rows.length === 0, 'The series has fewer than 8 valid rows: paste a longer record or clear the table to use the site statistics or the synthetic series.');
      const mm = ctx.site?.waves?.monthlyMean;
      const useSite = Array.isArray(mm) && mm.length === 12 && mm.every((v) => Number.isFinite(v) && v > 0);
      ser = syntheticWeather({ years: inp.years, dtHours: dt, hsMean: inp.hsMean, seasonality: inp.season, tauHours: inp.tau, monthlyMean: useSite ? mm : null, tpFactor: ctx.site?.waves?.tpFactor > 0 ? ctx.site.waves.tpFactor : 4.4, seed: inp.seed });
      source = useSite ? `synthetic series calibrated to the site hindcast monthly mean Hs (${inp.years} years at ${dt} h, seed ${inp.seed})` : `seeded synthetic series (${inp.years} years at ${dt} h, mean Hs ${inp.hsMean} m, seed ${inp.seed}) — indicative only`;
    }
    const n = ser.hs.length, spanDays = n * dt / 24;
    const ref = operationReferencePeriod({ plannedHours: inp.tpop, contingencyHours: inp.tc });
    const limits = [{ key: 'hs', max: inp.hsLim, alpha: inp.alphaHs }, { key: 'wind', max: inp.windLim, alpha: inp.alphaWind }];
    if (inp.tpMax > 0) limits.push({ key: 'tp', max: inp.tpMax });
    const series = { hs: ser.hs, tp: ser.tp, wind: ser.wind };
    const monthly = spanDays >= 300;
    const op = run(() => operability({ series, limits, dtHours: dt, requiredHours: ref.TR, timeMs: monthly ? ser.timeMs : undefined }), 'Operability');
    const opwf = forecastLimit(inp.hsLim, inp.alphaHs);
    const nWin = op.windows.length, winPerYear = nWin / (spanDays / 365.25);
    const meanWait = Number.isFinite(op.meanWaitingHours) ? op.meanWaitingHours : n * dt, p90Wait = Number.isFinite(op.p90WaitingHours) ? op.p90WaitingHours : n * dt;
    // Campaign: nOps identical operations, each waiting for a window of TR and then taking TPOP.
    let camp = null, campErr = '';
    const startMonths = monthly ? SEASONS[inp.startSeason] : null;
    try { camp = campaignMonteCarlo({ series, dtHours: dt, operations: Array.from({ length: inp.nOps }, (_, i) => ({ name: `Operation ${i + 1}`, limits, netHours: inp.tpop, contingencyHours: inp.tc })), nTrials: inp.nTrials, seed: inp.seed, timeMs: startMonths ? ser.timeMs : undefined, startMonths: startMonths || undefined }); } catch (e) { campErr = e.message; }
    const gap = inp.gapHours * inp.nOps;
    const netDays = (inp.nOps * inp.tpop + gap) / 24;
    const p = camp ? camp.completionHours : null;
    const d10 = p ? (p.p10 + gap) / 24 : netDays, d50 = p ? (p.p50 + gap) / 24 : netDays, d90 = p ? (p.p90 + gap) / 24 : netDays;
    const down50 = camp ? camp.downtimeHours.p50 / 24 : 0, down90 = camp ? camp.downtimeHours.p90 / 24 : 0;
    // Persistence of calms below the forecast limit, Markov chain on Hs and a Weibull fit, as supporting statistics.
    const durs = [6, 12, 24, 48, 72], ths = [0.5 * opwf, 0.75 * opwf, opwf, 1.25 * opwf, 1.5 * opwf].map((v) => round(v, 3));
    const pers = persistenceTable({ values: ser.hs, dtHours: dt, thresholds: ths, durations: durs, sense: 'below' });
    let hsMax = 0, hsSum = 0; for (let i = 0; i < n; i++) { if (ser.hs[i] > hsMax) hsMax = ser.hs[i]; hsSum += ser.hs[i]; }
    let operMarkov = op.workableFraction, markovOcc = [], markovEdges = [];
    try {
      const edges = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6].filter((e) => e < hsMax); edges.push(Math.max(hsMax * 1.001, edges[edges.length - 1] + 0.5));
      const mk = markovFit(Array.from(ser.hs), edges);
      const sim = markovSimulate(mk, { nSteps: Math.min(Math.round(10 * 8766 / dt), 40000), seed: inp.seed });
      let w = 0; for (let i = 0; i < sim.values.length; i++) if (sim.values[i] <= opwf) w++;
      operMarkov = w / sim.values.length; markovOcc = mk.stationary; markovEdges = edges;
    } catch { /* supporting statistic only */ }
    let wb = null; try { wb = fitWeibull(Array.from(ser.hs).filter((v) => v > 0)); } catch { wb = null; }
    const pWb = wb ? 1 - Math.exp(-((opwf / wb.scale) ** wb.shape)) : 0;
    let hsBelow = 0; for (let i = 0; i < n; i++) if (ser.hs[i] <= opwf) hsBelow++;
    const checks = [
      check('Operability against the planning target', inp.operReq, Math.max(op.startableFraction * 100, 1e-6), '% (required / achieved)', `${N001} Sec. 2 (weather-restricted operations: window of TR below the forecast limit); project planning target`, 'Share of time at which a window long enough for the whole operation plus contingency begins. Low operability means long standby and a schedule that depends on luck.'),
      check('Weather-restricted classification: TR ≤ 96 h and TPOP ≤ 72 h', Math.max(ref.TR / 96, inp.tpop / 72), 1, '-', `${N001} Sec. 2 (operation reference period)`, 'Longer operations cannot rely on a forecast and must be designed for seasonal extremes as weather-unrestricted operations.'),
    ];
    const ok = checks.every((c) => c.uc <= 1);
    const best = monthly ? op.monthly.reduce((a, m) => (m.startableFraction > a.startableFraction ? m : a), op.monthly[0]) : null;
    const worstM = monthly ? op.monthly.reduce((a, m) => (m.startableFraction < a.startableFraction ? m : a), op.monthly[0]) : null;
    // Sensitivity: operability against the design limit.
    const lims = [0.5, 0.75, 1, 1.25, 1.5, 2].map((f) => round(f * inp.hsLim, 3));
    const sens = lims.map((hl) => operability({ series, limits: [{ key: 'hs', max: hl, alpha: inp.alphaHs }, ...limits.slice(1)], dtHours: dt, requiredHours: ref.TR }).startableFraction * 100);
    const iUp = 3, gain = sens[iUp] - sens[2];
    const waitDays = d50 - netDays;
    return {
      kpis: [
        kpi('workable', 'Workable time (all limits met)', op.workableFraction * 100, '%', 'Share of the record with every limit below its forecast value, regardless of how long the calm lasts.'),
        kpi('operability', 'Operability (a window of TR can start)', op.startableFraction * 100, '%', 'Share of time from which the limits hold without interruption for the reference period TR = TPOP + TC. This is the figure for planning.', inp.operReq, 'min'),
        kpi('TR', 'Operation reference period TR', ref.TR, 'h', 'Planned operation period plus contingency: the length of weather window that must be forecast before starting.', 96, 'max'),
        kpi('opwf', 'Forecast limit on Hs (α × design limit)', opwf, 'm', 'The operation may start only when the forecast for the whole reference period is below this value.'),
        kpi('nWindows', 'Weather windows in the record', nWin, '-', 'Uninterrupted periods at least TR long.'),
        kpi('winPerYear', 'Weather windows per year', winPerYear, '1/yr', 'How many separate opportunities a year offers on average.'),
        kpi('meanWait', 'Mean waiting time for a window', meanWait, 'h', 'Average time from an arbitrary moment to the next start opportunity.'),
        kpi('p90Wait', 'P90 waiting time for a window', p90Wait, 'h', 'Waiting time exceeded one time in ten: use it for contingency in the schedule.'),
        kpi('d50', 'Campaign duration P50', d50, 'days', 'Median duration of the campaign including waiting on weather and the unrestricted work between operations.'),
        kpi('d90', 'Campaign duration P90', d90, 'days', 'Duration not exceeded in nine campaigns out of ten: the conservative schedule basis.'),
        kpi('d10', 'Campaign duration P10', d10, 'days', 'Optimistic outcome: one campaign in ten is faster.'),
        kpi('netDays', 'Campaign duration without weather', netDays, 'days', 'Net working time: operations plus unrestricted work.'),
        kpi('down50', 'Weather downtime P50', down50, 'days', 'Median waiting on weather over the campaign.'),
        kpi('down90', 'Weather downtime P90', down90, 'days', 'Weather downtime exceeded one time in ten.'),
        kpi('hsMeanRec', 'Mean Hs of the record', hsSum / n, 'm', 'Check that the series represents the site: compare with the metocean report.'),
        kpi('operMarkov', 'Workable time from the fitted Markov chain (Hs only)', operMarkov * 100, '%', 'Share of time below the Hs forecast limit in a ten-year series simulated from the fitted first-order Markov chain; it should agree with the record.'),
        kpi('pWeibull', 'Probability of Hs below the forecast limit (Weibull fit)', pWb * 100, '%', 'From a two-parameter Weibull distribution fitted to all Hs values; a cross-check of the empirical figure and a way to extrapolate to other limits.'),
        kpi('hsBelow', 'Share of time with Hs below the forecast limit', 100 * hsBelow / n, '%', 'Empirical non-exceedance of the Hs forecast limit alone, before wind, period and duration requirements.'),
      ],
      checks,
      plots: [
        ...(monthly ? [bar('Operability by month', MONTHS, 'Share of time', '%', [{ name: 'Workable', y: op.monthly.map((m) => round(m.workableFraction * 100)) }, { name: `Window of ${fmt(ref.TR)} h can start`, y: op.monthly.map((m) => round(m.startableFraction * 100)) }],
          'Seasonal pattern of the operation. The gap between the two bars is the cost of needing an uninterrupted window: short calms do not help. Plan the campaign in the months where the second bar is high.')] : []),
        line('Operability against the design limit', 'Design limit on Hs', 'm', 'Operability', '%', [{ name: 'Operability', x: lims, y: sens.map((v) => round(v)) }, { name: 'Target', x: [lims[0], lims[lims.length - 1]], y: [inp.operReq, inp.operReq], dash: true }],
          'What a higher (or lower) limiting sea state is worth. A steep curve at the present limit means design measures that raise the limit pay back quickly in vessel time.'),
        ...(camp ? [{ type: 'hist', title: 'Campaign duration distribution', values: Array.from(camp.samples).map((v) => round((v + gap) / 24)), bins: 24, x: { label: 'Campaign duration', unit: 'days' }, y: { label: 'Trials', unit: '-' },
          explain: 'Spread of the campaign duration over the random start times. A long right tail means a few bad-weather years dominate the risk: price the schedule on P90, not on the mean.' }] : []),
        line('Persistence of calms', 'Minimum duration', 'h', 'Share of time in calms at least this long', '%', pers.map((r) => ({ name: `Hs ≤ ${r.threshold} m`, x: durs, y: r.timeFraction.map((v) => round(v * 100)) })),
          'Share of the record occupied by calms that last at least the given duration, for several thresholds around the forecast limit. Read the curve at TR to see how the duration requirement eats into the workable time.'),
      ],
      tables: [
        tbl('Persistence of calms (Hs below threshold)', [['th', 'Threshold', 'm'], ['occ', 'Time below', '%'], ['ne', 'Events'], ['md', 'Mean duration', 'h'], ['mx', 'Longest', 'h'], ['f24', 'Time in calms ≥ 24 h', '%'], ['f48', 'Time in calms ≥ 48 h', '%']],
          pers.map((r) => ({ th: r.threshold, occ: round(r.occurrence * 100), ne: r.nEvents, md: round(r.meanDurationHours), mx: round(r.maxDurationHours), f24: round(r.timeFraction[2] * 100), f48: round(r.timeFraction[3] * 100) })),
          'One row per Hs threshold: how often the sea is below it, in how many separate calms, and how long they last. Use it to choose between a shorter operation and a higher limit.'),
        ...(monthly ? [tbl('Monthly statistics', [['m', 'Month'], ['w', 'Workable', '%'], ['s', 'Operability', '%'], ['wt', 'Mean wait', 'h']], op.monthly.map((m) => ({ m: MONTHS[m.month - 1], w: round(m.workableFraction * 100), s: round(m.startableFraction * 100), wt: round(Number.isFinite(m.meanWaitingHours) ? m.meanWaitingHours : 0) })),
          'Workable share, operability and mean waiting time by calendar month of the start.')] : []),
        ...(markovOcc.length ? [tbl('Fitted Markov chain: sea-state classes', [['c', 'Hs class', 'm'], ['p', 'Stationary probability', '%']], markovOcc.map((v, i) => ({ c: `${round(markovEdges[i], 3)}–${round(markovEdges[i + 1], 3)}`, p: round(v * 100) })),
          'Long-run share of time in each Hs class from the fitted transition matrix; it should reproduce the scatter of the record.')] : []),
      ],
      notes: [`Series: ${source}.`, `Forecast limit Hs ${fmt(opwf)} m (α = ${inp.alphaHs}), wind ${fmt(inp.windLim * inp.alphaWind)} m/s; TR = ${fmt(ref.TR)} h${ref.weatherRestricted ? '' : ' — outside the weather-restricted class'}.`, monthly && best ? `Best month ${MONTHS[best.month - 1]} (${fmt(best.startableFraction * 100)} %), worst ${MONTHS[worstM.month - 1]} (${fmt(worstM.startableFraction * 100)} %).` : 'The series is shorter than a year: no monthly breakdown.', campErr ? `Campaign simulation not available: ${campErr}.` : `Campaign of ${inp.nOps} operations: P50 ${fmt(d50)} days, P90 ${fmt(d90)} days against ${fmt(netDays)} days net.`],
      advice: [
        ok ? adv('engineering', op.startableFraction * 100 < 1.2 * inp.operReq ? 'warn' : 'ok', `Operability is ${fmt(op.startableFraction * 100)} % for a ${fmt(ref.TR)} h window below Hs ${fmt(opwf)} m (workable ${fmt(op.workableFraction * 100)} %); the mean wait for a window is ${fmt(meanWait)} h, P90 ${fmt(p90Wait)} h.`,
          `${monthly && best ? `Schedule the operation around ${MONTHS[best.month - 1]} where operability reaches ${fmt(best.startableFraction * 100)} %. ` : ''}Confirm the alpha-factor (${inp.alphaHs}) from the governing table for TPOP = ${inp.tpop} h and the forecast level, and repeat with the measured or hindcast series for the site.`, `${N001} Sec. 2; DNV-RP-N103 (2017) §9`)
          : adv('engineering', 'fail', `Operability is ${fmt(op.startableFraction * 100)} % against a target of ${inp.operReq} %${ref.weatherRestricted ? '' : `, and TR = ${fmt(ref.TR)} h is outside the weather-restricted class`}.`,
            `Raise the design limit (to ${fmt(lims[iUp])} m gives ${fmt(sens[iUp])} %), shorten the operation or split it at safe hold points so that TR falls, improve the forecast level to gain alpha, or move the campaign to the calm season.`, `${N001} Sec. 2; DNV-RP-N103 (2017) §9`),
        adv('economic', 'info', `The campaign of ${inp.nOps} operations takes ${fmt(netDays)} days net and ${fmt(d50)} days at P50 (${fmt(d90)} days at P90): about ${fmt(waitDays)} days of waiting on weather at P50.${inp.dayRate > 0 ? ` At the indicative spread rate entered, the P50 waiting is about ${fmt(waitDays * inp.dayRate / 1e6)} million USD and the P90–P50 exposure ${fmt((d90 - d50) * inp.dayRate / 1e6)} million USD (indicative).` : ''}`,
          `Raising the design limit from ${inp.hsLim} to ${fmt(lims[iUp])} m would add ${fmt(gain)} percentage points of operability: compare the vessel days saved with the cost of the design change. Carry the P90–P50 difference of ${fmt(d90 - d50)} days as schedule contingency.`,
          'Schedule risk analysis on P50/P90 (AACE RP 57R-09); waiting on weather priced at the spread day rate; DNV-ST-N001 Sec. 2 for the window definition',
          `Waiting on weather at P50: ${vesselCo2Text(Math.max(waitDays, 0), 25)} A higher limit or a calmer season cuts this directly.`),
      ],
      publish: { 'install.operabilityFrac': op.startableFraction, 'install.vesselDays': d50, 'install.waitHoursMean': meanWait },
      quantities: { vessel_days: { installationSpread: round(d50) } },
    };
  },
  verify: [
    { name: 'Threshold and duration logic on a hand-analysed series', source: 'Hand count: 24 samples at 3 h, limit Hs ≤ 1.5 m (α = 1), TR = 9 h = 3 steps. Calm runs of length 4, 2, 5 and 6 → workable 17/24 = 70.83 %; start possible at 2 + 0 + 3 + 4 = 9 samples → 37.5 %; windows of at least 9 h: 3',
      inputs: { series: [1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2].map((hs) => ({ hs, tp: 8, wind: 5 })), dtHours: 3, hsLim: 1.5, alphaHs: 1, tpop: 6, tc: 3, nOps: 2, nTrials: 50 },
      expect: [{ kpi: 'workable', value: 100 * 17 / 24, tol: 1e-9 }, { kpi: 'operability', value: 37.5, tol: 1e-9 }, { kpi: 'nWindows', value: 3, tol: 1e-9 }, { kpi: 'TR', value: 9, tol: 1e-9 }] },
    { name: 'Always-calm series: fully workable, no weather downtime, end-of-record effect on start opportunities', source: 'Limiting case: a constant sea state below every limit is 100 % workable and the campaign takes exactly its net time, 2 × 6 h + 2 × 12 h = 1.5 days. In a 40-sample record a 3-step window can start at samples 1–38 only (the last two cannot complete inside the record): 38/40 = 95 %',
      inputs: { series: Array.from({ length: 40 }, () => ({ hs: 0.5, tp: 6, wind: 4 })), dtHours: 3, tpop: 6, tc: 3, nOps: 2, nTrials: 50 },
      expect: [{ kpi: 'operability', value: 95, tol: 1e-9 }, { kpi: 'workable', value: 100, tol: 1e-9 }, { kpi: 'd50', value: 1.5, tol: 1e-9 }, { kpi: 'down90', value: 0, abs: 1e-9 }] },
  ],
};
