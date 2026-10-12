// Riser fatigue life and integrity: combination of fatigue damage sources with design fatigue factors,
// consumed and remaining life, fracture-mechanics crack growth from a weld flaw and inspection interval.
import { parisLawBs7910, surfaceCrackGrowth, parisLifeClosedForm } from '../../engine/fatigue/fracture.js';
import { fatigueLife } from '../../engine/fatigue/damage.js';
import { pod } from '../../engine/fatigue/inspection.js';
import {
  need, run, fmt, round, pipeInputs, adv, check, kpi, plot, finish, vesselCo2Text, steelCo2Text, CONFIRM,
} from './_common.js';

export default {
  id: 's09-fatigue-integrity', suite: 9, title: 'Fatigue life and integrity assessment',
  summary: 'Combined fatigue damage from waves, VIV, slugging and installation, factored and remaining life, crack growth from a girth-weld flaw and the in-service inspection interval.',
  covers: [14], fidelity: 'analytical',
  standards: ['DNV-ST-F201 (2018) Sec. 5 F fatigue limit state', 'DNV-RP-F204 (2019) riser fatigue', 'DNV-RP-C203 (2021) §2.2, §6 fracture mechanics', 'BS 7910:2019 clause 8 fatigue crack growth', 'DNV-RP-C210 (2019) inspection planning', 'API RP 2RD (1998, R2013) §5.6'],
  models: ['Palmgren–Miner damage summation with design fatigue factors', 'Fracture mechanics: Paris–Erdogan law with Newman–Raju surface-crack solution', 'Closed-form Paris life', 'Probability of detection curve', 'Damage-tolerant inspection interval'],
  limits: 'Damage rates are inputs (from the wave-fatigue, VIV and internal-flow modules or from measurements) and are assumed to occur at the same hot spot — conservative when the critical locations differ. Crack growth uses a constant-amplitude equivalent stress range, no weld-toe magnification unless folded into the stress, no crack-growth threshold (conservative with an equivalent range), no residual-stress or fracture (toughness) limit; the crack is taken to failure at through-wall. Probabilistic inspection updating is outside this module.',
  publishes: ['life.fatigueYears', 'riser.remainingLife', 'riser.inspectionInterval'],
  searchTerms: 'riser fatigue life integrity fracture mechanics inspection interval',
  inputs: [
    { key: 'dWave', label: 'Wave / vessel-motion damage rate', unit: '1/yr', type: 'number', default: 0.0013, min: 0, max: 10, group: 'Fatigue damage sources', role: 'load', bind: 'riser.waveDamagePerYear', uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Unfactored Miner damage per year at the critical weld.' },
    { key: 'dViv', label: 'VIV damage rate', unit: '1/yr', type: 'number', default: 0.0002, min: 0, max: 10, group: 'Fatigue damage sources', role: 'load', bind: 'riser.vivDamagePerYear', uncertain: { dist: 'lognormal', cov: 0.8 } },
    { key: 'dSlug', label: 'Slugging / operational damage rate', unit: '1/yr', type: 'number', default: 0.0001, min: 0, max: 10, group: 'Fatigue damage sources', role: 'load', bind: 'riser.slugDamagePerYear' },
    { key: 'dInstall', label: 'Installation damage (one-off)', unit: '-', type: 'number', default: 0.005, min: 0, max: 1, group: 'Fatigue damage sources', role: 'ic', help: 'Miner damage accumulated during lay, tow and hang-off: the initial condition of the fatigue budget.' },
    { key: 'yearsInService', label: 'Years in service', unit: 'yr', type: 'number', default: 8, min: 0, max: 80, group: 'Service history', role: 'ic' },
    { key: 'dMeasured', label: 'Additional accumulated damage from monitoring', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'Service history', role: 'ic', help: 'Damage already consumed that is not covered by the rates above (e.g. a recorded storm or VIV event).' },
    { key: 'designLife', label: 'Design life', unit: 'yr', type: 'number', default: 25, min: 1, max: 80, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
    { key: 'dff', label: 'Design fatigue factor (wave, slug, installation)', unit: '-', type: 'number', default: 10, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: `3 / 6 / 10 for low / normal / high safety class (DNV-ST-F201). ${CONFIRM}` },
    { key: 'dffViv', label: 'Design fatigue factor (VIV)', unit: '-', type: 'number', default: 10, min: 1, max: 30, group: 'Acceptance criteria', role: 'criteria', help: `Often taken higher than the wave factor because VIV prediction is more uncertain. ${CONFIRM}` },
    ...pipeInputs().filter((i) => ['D', 't'].includes(i.key)).map((i) => ({ ...i, group: 'Weld and flaw' })),
    { key: 'dSigma', label: 'Equivalent stress range at the weld', unit: 'MPa', type: 'number', default: 20, min: 1, max: 400, group: 'Weld and flaw', role: 'load', bind: 'riser.fatigueStressRange', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Constant-amplitude hot-spot range giving the same damage as the real spectrum.' },
    { key: 'cyclesYr', label: 'Stress cycles per year', unit: '1/yr', type: 'number', default: 4e6, min: 1e3, max: 5e7, group: 'Weld and flaw', role: 'load', bind: 'riser.fatigueCyclesPerYear' },
    { key: 'bendFrac', label: 'Bending share of the stress range', unit: '-', type: 'number', default: 0.3, min: 0, max: 1, group: 'Weld and flaw', role: 'load', help: 'Through-wall bending / total range at the weld (0 = pure membrane, the usual case for global riser bending on a thin wall).' },
    { key: 'a0', label: 'Initial flaw depth', unit: 'mm', type: 'number', default: 1.5, min: 0.1, max: 20, group: 'Weld and flaw', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Largest surface-breaking flaw the fabrication NDT may leave (acceptance criterion or ECA flaw size).' },
    { key: 'c0', label: 'Initial flaw half-length', unit: 'mm', type: 'number', default: 7.5, min: 0.2, max: 200, group: 'Weld and flaw', role: 'ic' },
    { key: 'aDetect', label: 'Reliably detectable flaw depth in service', unit: 'mm', type: 'number', default: 3, min: 0.3, max: 20, group: 'Inspection', role: 'operation', help: 'Depth found with high probability by the in-service method (e.g. internal UT tool).' },
    { key: 'podX0', label: 'Inspection POD median depth', unit: 'mm', type: 'number', default: 1.2, min: 0.1, max: 20, group: 'Inspection', role: 'operation', help: 'Depth detected with 50 % probability, POD(a) = 1 − 1/(1 + (a/x0)^b).' },
    { key: 'podB', label: 'Inspection POD exponent', unit: '-', type: 'number', default: 1.8, min: 0.5, max: 6, group: 'Inspection', role: 'operation' },
    { key: 'inspFactor', label: 'Safety factor on the inspection interval', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Inspection', role: 'criteria', help: 'Interval = time for a just-missed crack to grow through the wall divided by this factor.' },
    { key: 'environment', label: 'Crack-growth environment', type: 'select', default: 'marine-cp850', group: 'Weld and flaw', role: 'material', options: [{ value: 'air', label: 'Air / dry internal' }, { value: 'marine-cp850', label: 'Sea water with cathodic protection (−850 mV)' }, { value: 'marine-free', label: 'Sea water, free corrosion' }] },
    { key: 'curve', label: 'Crack-growth law', type: 'select', default: 'mean+2sd', group: 'Weld and flaw', role: 'material', options: [{ value: 'mean+2sd', label: 'BS 7910 two-stage, mean + 2 SD' }, { value: 'mean', label: 'BS 7910 two-stage, mean' }, { value: 'simplified', label: 'BS 7910 simplified (single stage)' }] },
    { key: 'inspDays', label: 'Vessel days per inspection campaign', unit: 'days', type: 'number', default: 4, min: 0, max: 60, group: 'Inspection', role: 'operation' },
    { key: 'dayRate', label: 'Indicative inspection vessel day rate', unit: 'USD', type: 'number', default: 120000, min: 0, max: 1e6, group: 'Inspection', role: 'cost', help: 'Indicative unit rate; replace with the contract rate.' },
  ],

  compute(inp) {
    // engineering units → SI: mm → m, MPa → Pa
    const rate = inp.dWave + inp.dViv + inp.dSlug;
    need(rate > 0, 'All damage rates are zero: there is nothing to assess. Enter at least one damage rate.');
    const rateF = inp.dff * (inp.dWave + inp.dSlug) + inp.dffViv * inp.dViv;
    const D0 = inp.dInstall + inp.dMeasured, D0F = inp.dff * D0;
    need(D0 < 1, 'Accumulated one-off damage already reaches 1.0: the weld is at the end of its characteristic fatigue life.');
    const fl = fatigueLife({ damageRate: rate, designLife: inp.designLife, dff: rateF / rate, accumulatedDamage: D0 });
    const lifeChar = (1 - D0) / rate, lifeF = Math.max(0, (1 - D0F) / rateF);
    const usedNow = D0F + rateF * inp.yearsInService, remaining = (1 - usedNow) / rateF;
    const ucEnd = D0F + rateF * inp.designLife;
    // fracture mechanics
    const t = inp.t / 1e3, a0 = inp.a0 / 1e3, c0 = inp.c0 / 1e3, ad = inp.aDetect / 1e3;
    need(a0 < 0.8 * t && ad < 0.9 * t, `Flaw depths must be well inside the wall (${inp.t} mm).`);
    need(2 * c0 * Math.max(1, inp.aDetect / inp.a0) < 0.5 * Math.PI * (inp.D - inp.t) / 1e3, 'The flaw is longer than half the pipe circumference (at the initial or detectable depth): treat it as a fully circumferential flaw in a dedicated assessment.');
    need(c0 >= a0, 'The flaw half-length must be at least its depth (semi-elliptical surface flaw, a/c ≤ 1).');
    // threshold set to zero: an equivalent constant range hides the larger cycles of the real spectrum that do exceed it
    const envS = inp.environment === 'air' ? 'air' : 'marine';
    const law = run(() => parisLawBs7910({ environment: inp.curve === 'simplified' ? envS : inp.environment, curve: inp.curve, threshold: 0 }), 'Crack-growth law');
    const lawS = parisLawBs7910({ environment: envS, curve: 'simplified', threshold: 0 }).stages[0];
    const ds = inp.dSigma * 1e6, W = Math.PI * (inp.D / 1e3 - t);
    const grow = (a, c) => run(() => surfaceCrackGrowth({ a0: a, c0: c, t, W, membraneRange: (1 - inp.bendFrac) * ds, bendingRange: inp.bendFrac * ds, law, maxCycles: 5e11 }), 'Crack-growth integration failed');
    const g0 = grow(a0, c0), gd = grow(ad, Math.max(c0, ad * c0 / a0));
    const CAP = 1e5, yrs = (g) => (Number.isFinite(g.cycles) && g.status !== 'arrested' && g.status !== 'max-cycles' ? Math.min(CAP, g.cycles / inp.cyclesYr) : CAP);
    const crackLife = yrs(g0), crackDet = yrs(gd), interval = Math.max(crackDet / inp.inspFactor, 1e-3);
    const closed = parisLifeClosedForm({ a0, af: t, stressRange: ds, Y: 1.12, A: lawS.A, m: lawS.m }) / inp.cyclesYr;
    const left = Math.max(0, inp.designLife - inp.yearsInService);
    const nIn = interval >= CAP ? 0 : Math.max(0, Math.floor((left - 1e-9) / interval));
    const podD = pod(ad, { x0: inp.podX0 / 1e3, b: inp.podB });
    const shares = [['Waves and vessel motion', inp.dWave * inp.dff], ['VIV', inp.dViv * inp.dffViv], ['Slugging / operations', inp.dSlug * inp.dff]];
    const govS = shares.reduce((a, s) => (s[1] > a[1] ? s : a), shares[0]);

    const kpis = [
      kpi('damageRate', 'Total annual damage (unfactored)', rate, '1/yr', 'Sum of the Miner damage rates of all cyclic sources at the critical weld, before design fatigue factors.'),
      kpi('lifeYears', 'Characteristic fatigue life', lifeChar, 'yr', 'Years to a Miner sum of 1.0 from new without design fatigue factors, after deducting installation damage. It is not the usable life.'),
      kpi('lifeFactored', 'Factored fatigue life', lifeF, 'yr', 'Years until the sum of damage × design fatigue factor reaches 1.0. This is the life that must cover the design life.', inp.designLife, 'min'),
      kpi('usedNow', 'Factored damage consumed to date', usedNow, '-', 'Share of the fatigue budget (Σ DFF·D ≤ 1) used after the years in service, including installation and any recorded events.', 1, 'max'),
      kpi('remainingLife', 'Remaining factored life', remaining, 'yr', 'Years of further service before the factored damage reaches 1.0 at the present damage rates. Negative means the budget is already exceeded.'),
      kpi('ucEnd', 'Factored damage at the end of the design life', ucEnd, '-', 'Fatigue utilisation projected to the end of the design life. Below 1.0 the fatigue limit state is satisfied; life extension is possible while it stays below 1.0.', 1, 'max'),
      kpi('crackLife', 'Crack growth life from the initial flaw', crackLife, 'yr', `Time for the assumed fabrication flaw to grow through the wall under the equivalent stress range (capped at ${CAP} when growth arrests below threshold). Compare with the S–N life as an independent estimate.`),
      kpi('crackLifeDetect', 'Growth time from a detectable crack to through-wall', crackDet, 'yr', 'Time for a crack of the reliably detectable depth to penetrate the wall: the window in which inspection can still find it.'),
      kpi('inspInterval', 'Maximum inspection interval', interval, 'yr', 'Growth time from the detectable size divided by the safety factor, so that a crack just missed at one inspection is found at the next before it is critical.'),
      kpi('nInspections', 'Inspections in the remaining design life', nIn, '-', 'Number of in-service inspections needed between now and the end of the design life at that interval.'),
      kpi('podDetect', 'Probability of detecting the stated detectable flaw', podD, '-', 'Probability of detection of the "detectable" depth from the POD curve. If it is well below 0.9 the detectable size or the method should be revised.'),
      kpi('parisYears', 'Closed-form Paris life (simplified law, Y = 1.12)', Math.min(CAP * 10, closed), 'yr', 'Hand-checkable estimate of growth from the initial flaw to through-wall with a constant geometry factor and the single-stage BS 7910 simplified law: a conservative order-of-magnitude check on the numerical integration.'),
    ];
    const checks = [
      check('Fatigue limit state at the end of the design life', ucEnd, 1, '-', 'DNV-ST-F201 (2018) Sec. 5 F: Σ D_i·DFF_i ≤ 1; DNV-RP-F204 (2019) §6', 'Requires the factored damage from all sources to stay below 1.0 over the design life. Otherwise a fatigue crack at the critical girth weld is credible in service.'),
      check('Fatigue budget consumed to date', usedNow, 1, '-', 'DNV-RP-F204 (2019) §6 (reassessment and life extension)', 'Compares the factored damage accumulated so far with the budget. Above 1.0 continued operation needs inspection-based justification or mitigation.'),
    ];
    if (interval < CAP) checks.push(check('Practical inspection interval (damage tolerance)', 1, interval, 'yr', 'DNV-RP-C210 (2019); DNV-RP-C203 (2021) §6', 'Compares a one-year minimum practical offshore inspection interval with the calculated maximum interval. Above 1.0 the weld is not damage-tolerant: it must be upgraded rather than managed by inspection.'));

    const tEnd = Math.max(inp.designLife * 1.2, Math.min(lifeF * 1.1, inp.designLife * 4), inp.yearsInService + 1), tx = [], dy = [];
    for (let i = 0; i <= 40; i++) { const y = tEnd * i / 40; tx.push(y); dy.push(D0F + rateF * y); }
    const hist = g0.history, hy = hist.cycles.map((n) => n / inp.cyclesYr), step = Math.max(1, Math.floor(hy.length / 200));
    const hx = hy.filter((_, i) => i % step === 0 || i === hy.length - 1), ha = hist.a.filter((_, i) => i % step === 0 || i === hy.length - 1).map((v) => v * 1e3), hc = hist.c.filter((_, i) => i % step === 0 || i === hy.length - 1).map((v) => v * 1e3);
    const as = [], ps = []; for (let i = 1; i <= 40; i++) { const a = t * i / 40; as.push(a * 1e3); ps.push(pod(a, { x0: inp.podX0 / 1e3, b: inp.podB })); }
    const plots = [
      plot('Factored damage accumulation', 'Years from start of service', 'yr', 'Σ DFF·D', '-', [{ name: 'Factored damage', x: tx, y: dy }, { name: 'Limit', x: [0, tEnd], y: [1, 1], dash: true }, { name: 'Today', x: [inp.yearsInService], y: [usedNow] }, { name: 'End of design life', x: [inp.designLife], y: [ucEnd] }],
        'Fatigue budget consumed against time, starting from installation damage. The line must stay below 1.0 at the design life; its crossing of the limit is the factored fatigue life.'),
      { type: 'bar', title: 'Factored damage rate by source', categories: shares.map((s) => s[0]), x: { label: 'Source', unit: '' }, y: { label: 'DFF × damage per year', unit: '1/yr' }, series: [{ name: 'Factored annual damage', y: shares.map((s) => s[1]) }],
        explain: 'Contribution of each loading source to the factored annual damage. The tallest bar is where analysis refinement or mitigation pays most.' },
      plot('Crack growth from the initial flaw', 'Years', 'yr', 'Crack size', 'mm', [{ name: 'Depth a', x: hx, y: ha }, { name: 'Half-length c', x: hx, y: hc }, { name: 'Wall thickness', x: [0, hx[hx.length - 1] || 1], y: [inp.t, inp.t], dash: true }, { name: 'Detectable depth', x: [0, hx[hx.length - 1] || 1], y: [inp.aDetect, inp.aDetect], dash: true }],
        'Depth and half-length of the surface crack against time. Growth is slow for most of the life and accelerates near the end; the time between the detectable depth and the wall thickness is the inspection window.'),
      plot('Probability of detection', 'Crack depth', 'mm', 'POD', '-', [{ name: 'POD curve', x: as, y: ps }, { name: 'Detectable depth', x: [inp.aDetect], y: [podD] }], 'Probability that the in-service inspection finds a crack of a given depth. The marker shows the depth assumed to be reliably detectable.'),
    ];
    const tables = [
      { title: 'Fatigue damage budget', columns: [{ key: 'src', label: 'Source' }, { key: 'd', label: 'Damage rate', unit: '1/yr' }, { key: 'dff', label: 'DFF', unit: '-' }, { key: 'df', label: 'Factored rate', unit: '1/yr' }, { key: 'life', label: 'Design-life damage', unit: '-' }, { key: 'share', label: 'Share', unit: '%' }],
        rows: [['Waves and vessel motion', inp.dWave, inp.dff], ['VIV', inp.dViv, inp.dffViv], ['Slugging / operations', inp.dSlug, inp.dff]].map(([src, d, f]) => ({ src, d: Number(d.toExponential(3)), dff: f, df: Number((d * f).toExponential(3)), life: round(d * f * inp.designLife, 4), share: round(100 * d * f / rateF, 3) }))
          .concat([{ src: 'Installation + recorded events (one-off)', d: '–', dff: inp.dff, df: '–', life: round(D0F, 4), share: '–' }]),
        explain: 'Each damage source with its design fatigue factor and its contribution over the design life. The total of the design-life column is the fatigue utilisation.' },
      { title: 'Integrity plan', columns: [{ key: 'item', label: 'Item' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
        rows: [['Inspection-critical location', `${govS[0]} hot spot (touchdown zone and hang-off welds)`, ''], ['Maximum inspection interval', round(interval, 4), 'yr'], ['First inspection no later than', round(Math.min(interval, left), 4), 'yr from now'], ['Inspections to end of design life', nIn, '-'], ['Vessel days for inspections', nIn * inp.inspDays, 'days'], ['Remaining factored life', round(remaining, 4), 'yr'], ['Crack-growth law', law.label ?? inp.curve, '']].map(([item, value, unit]) => ({ item, value, unit })),
        explain: 'Summary for the integrity management plan: where to look, how often, and the vessel time it implies over the remaining design life.' },
    ];
    const ucMax = Math.max(...checks.map((c) => c.uc));
    const redNeed = ucEnd > 1 ? (ucEnd - D0F) / (1 - D0F) : 1;
    const advice = [
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.8 ? 'warn' : 'ok',
        `Factored damage at ${inp.designLife} years is ${fmt(ucEnd)} (${fmt(usedNow)} consumed after ${inp.yearsInService} years); factored life ${fmt(lifeF)} yr, remaining ${fmt(remaining)} yr. ${govS[0]} contributes most. Crack growth from a ${inp.a0} mm flaw takes ${fmt(crackLife)} yr.`,
        ucEnd > 1 ? `The damage rate must fall by a factor of ${fmt(redNeed)} — a stress-range reduction of ${fmt((1 - redNeed ** (-1 / 3)) * 100)} % on an m = 3 curve. Target ${govS[0].toLowerCase()} first: better weld class, strakes, or a lower design fatigue factor justified by an inspection programme.`
          : `Inspect the critical welds at intervals of no more than ${fmt(interval)} years with a method whose POD at ${inp.aDetect} mm is at least 0.9 (now ${fmt(podD)}); replace the predicted damage rates with monitored ones when strain or motion data are available.`,
        `DNV-ST-F201 (2018) Sec. 5 F and DNV-RP-F204 (2019) §6; crack growth per BS 7910:2019 clause 8; inspection planning per DNV-RP-C210 (2019). ${CONFIRM}`),
      adv('economic', 'info',
        `${nIn} inspection campaign(s) of ${inp.inspDays} vessel days each are needed to the end of the design life: ${fmt(nIn * inp.inspDays)} vessel days, about ${fmt(nIn * inp.inspDays * inp.dayRate / 1e6)} million USD at the indicative day rate (unit rate × days). Remaining factored life beyond the design life is ${fmt(Math.max(0, remaining - left))} years of potential life extension.`,
        'Use risk-based inspection: concentrate on the critical welds and extend the interval where monitoring shows lower damage than predicted. A demonstrated margin supports life extension, which defers replacement of the whole riser.',
        'Risk-based inspection per DNV-RP-G101 / API RP 580; value of information of monitoring; whole-life costing ISO 15686-5.',
        `Inspection campaigns: ${vesselCo2Text(Math.max(1, nIn * inp.inspDays))} Extending the life of an existing riser instead of replacing it avoids the steel of a new one — for a 200 t riser, ${steelCo2Text(200)}`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`Factored damage rate ${fmt(rateF)} per year (${govS[0]} governs); one-off damage ${fmt(D0)} unfactored. Cross-check (1 − D0)/(DFF·ḋ) with the blended factor: ${fmt(fl.remainingLife)} yr.`,
        `Crack growth: ${law.label ?? inp.curve}, status "${g0.status}" from the initial flaw and "${gd.status}" from the detectable size; equivalent range ${inp.dSigma} MPa at ${fmt(inp.cyclesYr)} cycles per year.`],
      advice,
      publish: { 'life.fatigueYears': lifeF, 'riser.remainingLife': remaining, 'riser.inspectionInterval': Math.min(interval, CAP) },
      quantities: { vessel_days: { riser_inspection: nIn * inp.inspDays } },
    });
  },

  verify: [
    {
      // Miner: characteristic life = (1 − D0)/Σd = 1/(0.001 + 0.0005) = 666.7 yr; with DFF = 10 on both sources the
      // factored life is 66.67 yr; after 10 years in service 0.15 of the budget is used, 56.67 yr remain.
      name: 'Palmgren–Miner summation with design fatigue factors', source: 'DNV-RP-C203 (2021) Eq. 2.2.1; DNV-ST-F201 (2018) Sec. 5 F',
      inputs: { dWave: 0.001, dViv: 0.0005, dSlug: 0, dInstall: 0, dMeasured: 0, dff: 10, dffViv: 10, yearsInService: 10 },
      expect: [{ kpi: 'lifeYears', value: 1 / 0.0015, tol: 1e-9 }, { kpi: 'lifeFactored', value: 1 / 0.015, tol: 1e-9 }, { kpi: 'usedNow', value: 0.15, tol: 1e-9 }, { kpi: 'remainingLife', value: 1 / 0.015 - 10, tol: 1e-9 }],
    },
    {
      // Paris law, constant Y: N = [a_f^(1−m/2) − a_0^(1−m/2)] / [(1 − m/2)·A·(Y·Δσ·√π)^m].
      // BS 7910 simplified law in air: A = 5.21e-13, m = 3 (N, mm). a0 = 1.5 mm, af = 20.6 mm, Δσ = 30 MPa, Y = 1.12:
      // N = (20.6^−0.5 − 1.5^−0.5)/(−0.5·5.21e-13·(1.12·30·√π)³) ≈ 1.08e7 cycles → ≈ 2.7 yr at 4e6 cycles/yr.
      name: 'Closed-form Paris life with BS 7910 simplified law in air', source: 'Paris & Erdogan (1963); BS 7910:2019 clause 8 simplified law (A = 5.21e-13, m = 3, N–mm units)',
      inputs: { environment: 'air', curve: 'simplified', dSigma: 30 },
      expect: [{ kpi: 'parisYears', value: (20.6 ** -0.5 - 1.5 ** -0.5) / (-0.5 * 5.21e-13 * (1.12 * 30 * Math.sqrt(Math.PI)) ** 3) / 4e6, tol: 0.002 }],
    },
  ],
};
