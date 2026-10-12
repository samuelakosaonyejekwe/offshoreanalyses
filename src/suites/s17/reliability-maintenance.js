// Reliability and maintenance economics: Weibull life and replacement policies, repairable-system (NHPP)
// overhaul age, system availability by reliability block diagram and Markov repair model, fault and event trees,
// cost of downtime and failure, risk-based inspection interval and insurance split.
import { weibullLife, exponentialLife, lognormalLife, ageReplacement, blockReplacement, nhppPowerLaw, markovSteadyState, semiMarkovSteadyState, rbdReliability, faultTreeProbability, eventTree, rbiEconomics } from '../../engine/econ/lifecycle.js';
import { availability as unitAvailability, downtimeCost, costOfFailure, insuranceRecovery } from '../../engine/econ/opex.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, CO2_MGO, MGO_T_PER_DAY, need, run, round, fmt, mil, fM, pct, adv, check, checkMin, cost, tableRows, uniqueNames, linspace, range } from './_common.js';

const k$ = (v) => round(v / 1e3);
const KU = `k${CUR}`;
const SUBS = [
  { name: 'Turbine drivetrain and generator', mtbf: 9000, mttr: 160, n: 1, k: 1, repairCost: 3.2e5 },
  { name: 'Power converter and transformer', mtbf: 14000, mttr: 90, n: 1, k: 1, repairCost: 1.4e5 },
  { name: 'Array cable string', mtbf: 120000, mttr: 1100, n: 1, k: 1, repairCost: 3.5e6 },
  { name: 'Offshore substation transformers', mtbf: 70000, mttr: 1400, n: 2, k: 1, repairCost: 4.5e6 },
  { name: 'Export cable circuits', mtbf: 95000, mttr: 1900, n: 2, k: 1, repairCost: 9.0e6 },
  { name: 'Control, SCADA and communications', mtbf: 26000, mttr: 30, n: 2, k: 1, repairCost: 6.0e4 },
];
const BARRIERS = [
  { name: 'Condition monitoring fails to give warning', pFail: 25, cost: 4.0e5 }, { name: 'Protection fails to isolate the fault', pFail: 4, cost: 6.5e6 }, { name: 'Fire suppression fails', pFail: 8, cost: 2.4e7 },
];

/** Birth–death repair model: n units of failure rate λ, c repair crews of rate μ; states = number failed. */
function markovKofN(lambda, mu, n, k, crews) {
  const rates = range(n + 1, () => new Array(n + 1).fill(0));
  for (let j = 0; j <= n; j++) {
    if (j < n) rates[j][j + 1] = (n - j) * lambda;
    if (j > 0) rates[j][j - 1] = Math.min(j, crews) * mu;
  }
  return markovSteadyState({ rates, upStates: range(n - k + 1) }).availability;
}

export default {
  id: 's17-reliability-maintenance', suite: 17, title: 'Reliability, availability and maintenance economics',
  summary: 'Finds the cost-optimal preventive replacement age and overhaul age of a wearing component, the availability of the production system with redundancy and limited repair crews, the expected annual cost of failures and downtime, and the inspection interval with the lowest lifecycle cost.',
  covers: [8], fidelity: 'analytical',
  standards: ['IEC 60300-3-11:2009 (reliability-centred maintenance)', 'IEC 61078:2016 (reliability block diagrams) and IEC 61025:2006 (fault tree analysis)', 'IEC 61165:2006 (Markov techniques) and IEC 62502:2010 (event tree analysis)', 'ISO 14224:2016 (reliability and maintenance data) and ISO 20815:2018 (production assurance)', 'DNV-RP-G101 (2021) / API RP 580 (2016) (risk-based inspection)', 'Barlow & Hunter (1960), optimum preventive maintenance policies'],
  models: ['Weibull failure-time model', 'Exponential reliability model', 'Lognormal lifetime model', 'Renewal process theory (age and block replacement)', 'Non-homogeneous Poisson process (power law, minimal repair)', 'Markov chain availability (birth–death repair model)', 'Semi-Markov steady state', 'Reliability block diagram (series, parallel, k-out-of-n)', 'Fault-tree analysis', 'Event-tree analysis', 'Mean time between failures and mean time to repair', 'Expected downtime cost', 'Reliability-centred inspection planning', 'Insurance recovery (insured and retained loss)'],
  limits: `Steady-state availability with independent failures (no common cause), constant failure and repair rates in the system model, and perfect renewal at replacement. Weather-dependent access, spare-part logistics and queueing between subsystems are represented only through the repair time and the number of crews. Failure data are placeholders to be replaced by ISO 14224-style field data. ${CUR_NOTE}`,
  publishes: ['econ.availability', 'econ.correctiveAnnual', 'econ.downtimeCostAnnual', 'econ.optimalReplacementAge', 'econ.inspectionInterval'],
  searchTerms: 'age replacement Weibull optimal preventive maintenance availability Markov risk based inspection economics',
  inputs: [
    { key: 'shape', label: 'Weibull shape β of the wearing component', unit: '-', type: 'number', default: 2.6, min: 0.3, max: 10, group: 'Component life', role: 'material', uncertain: { dist: 'normal', cov: 0.15 }, help: 'β > 1 means wear-out (preventive replacement can pay); β ≤ 1 means random or early failures (it cannot).' },
    { key: 'scale', label: 'Weibull characteristic life η', unit: 'yr', type: 'number', default: 14, min: 0.1, max: 200, group: 'Component life', role: 'material', bind: 'life.fatigueYears', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Age by which 63 % have failed; may be taken from the fatigue or degradation analysis of the component.' },
    { key: 'units', label: 'Number of such components', unit: '-', type: 'int', default: 34, min: 1, max: 5000, group: 'Component life', role: 'operation' },
    cost({ key: 'cPrev', label: 'Planned replacement cost per component', unit: CUR, default: 4.2e5, min: 1, max: 1e10, group: 'Maintenance information', help: 'Parts, vessel and labour in a planned campaign.' }),
    cost({ key: 'cFail', label: 'Unplanned failure cost per component', unit: CUR, default: 1.6e6, min: 1, max: 1e10, group: 'Maintenance information', help: 'Emergency intervention, secondary damage and repair, excluding lost production.' }),
    { key: 'dtPrev', label: 'Downtime of a planned replacement', unit: 'hours', type: 'number', default: 36, min: 0, max: 20000, group: 'Maintenance information', role: 'operation' },
    { key: 'dtFail', label: 'Downtime of an unplanned failure', unit: 'hours', type: 'number', default: 520, min: 0, max: 50000, group: 'Maintenance information', role: 'operation', help: 'Includes waiting for parts, vessel and weather.' },
    { key: 'unitMW', label: 'Mean output lost per failed component', unit: 'MW (mean)', type: 'number', default: 7.2, min: 0, max: 5000, group: 'Production', role: 'operation', help: 'Average, not rated, output of the unit served by the component.' },
    { key: 'systemMW', label: 'Mean output of the whole system', unit: 'MW (mean)', type: 'number', default: 245, min: 0.001, max: 50000, group: 'Production', role: 'operation', help: 'Rated capacity × capacity factor.' },
    cost({ key: 'price', label: 'Net margin on lost production', unit: `${CUR}/MWh`, default: 95, min: 0, max: 1e5, group: 'Production', help: 'Price less variable cost saved.' }),
    { key: 'nhppLambda', label: 'Repairable system: power-law scale λ', unit: 'failures/yr^β', type: 'number', default: 0.12, min: 1e-6, max: 100, group: 'Repairable system (minimal repair)', role: 'material', help: 'Expected failures N(t) = λ t^β for a system restored to as-bad-as-old at each failure.' },
    { key: 'nhppBeta', label: 'Repairable system: power-law exponent β', unit: '-', type: 'number', default: 1.7, min: 0.3, max: 6, group: 'Repairable system (minimal repair)', role: 'material' },
    cost({ key: 'cMinimal', label: 'Cost of a minimal repair', unit: CUR, default: 1.8e5, min: 0, max: 1e10, group: 'Repairable system (minimal repair)' }),
    cost({ key: 'cOverhaul', label: 'Cost of overhaul or renewal', unit: CUR, default: 2.4e6, min: 1, max: 1e11, group: 'Repairable system (minimal repair)' }),
    { key: 'subsystems', label: 'Production system (series of subsystems)', type: 'table', group: 'System availability', role: 'ic',
      help: `Opening reliability state of the asset: the failure and repair data assumed at the start of the evaluation. Each row is a subsystem in series; n installed units of which k are needed. MTBF and MTTR in hours per unit; repair cost per failure. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Subsystem', type: 'text' }, { key: 'mtbf', label: 'MTBF', unit: 'h', type: 'number' }, { key: 'mttr', label: 'MTTR', unit: 'h', type: 'number' }, { key: 'n', label: 'Installed n', type: 'number' }, { key: 'k', label: 'Required k', type: 'number' }, { key: 'repairCost', label: 'Repair cost per failure', unit: CUR, type: 'number' }], default: SUBS },
    { key: 'crews', label: 'Repair crews available at once', unit: '-', type: 'int', default: 1, min: 1, max: 20, group: 'System availability', role: 'bc', help: 'Operational constraint: with fewer crews than failed units, repairs queue and redundancy is worth less than the block diagram suggests.' },
    { key: 'targetAvail', label: 'Required production availability', unit: '%', type: 'number', default: 96, min: 50, max: 99.99, group: 'Acceptance criteria', role: 'criteria', bind: 'energy.availability', help: 'Availability assumed in the energy-yield or production forecast.' },
    { key: 'initFreq', label: 'Initiating event frequency', unit: '1/yr', type: 'number', default: 0.35, min: 0, max: 100, group: 'Failure consequences (event tree)', role: 'load', help: 'For example electrical faults in the substation per year.' },
    cost({ key: 'baseEventCost', label: 'Cost when all barriers work', unit: CUR, default: 6.0e4, min: 0, max: 1e10, group: 'Failure consequences (event tree)' }),
    { key: 'barriers', label: 'Barriers', type: 'table', group: 'Failure consequences (event tree)', role: 'cost', help: `Probability that each barrier fails on demand and the extra consequence cost if it does. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Barrier failure', type: 'text' }, { key: 'pFail', label: 'Probability', unit: '%', type: 'number' }, { key: 'cost', label: 'Added consequence', unit: CUR, type: 'number' }], default: BARRIERS },
    cost({ key: 'deductible', label: 'Insurance deductible per event', unit: CUR, default: 1.0e6, min: 0, max: 1e10, group: 'Failure consequences (event tree)' }),
    cost({ key: 'insLimit', label: 'Insurance limit per event', unit: CUR, default: 5.0e7, min: 0, max: 1e12, group: 'Failure consequences (event tree)' }),
    cost({ key: 'inspCost', label: 'Cost of one inspection campaign', unit: CUR, default: 6.5e5, min: 0, max: 1e10, group: 'Inspection planning' }),
    cost({ key: 'inspFailCost', label: 'Cost of the failure the inspection guards against', unit: CUR, default: 6.0e7, min: 0, max: 1e12, group: 'Inspection planning', help: 'Repair, lost production and consequence of, for example, a through-thickness fatigue crack.' }),
    { key: 'inspEff', label: 'Inspection effectiveness', unit: '-', type: 'number', default: 0.8, min: 0.05, max: 1, group: 'Inspection planning', role: 'operation', help: 'Probability that an inspection finds the defect so that it is repaired (probability of detection).' },
    { key: 'inspShape', label: 'Degradation Weibull shape', unit: '-', type: 'number', default: 3, min: 0.5, max: 10, group: 'Inspection planning', role: 'material' },
    { key: 'inspScale', label: 'Degradation Weibull scale', unit: 'yr', type: 'number', default: 32, min: 1, max: 2000, group: 'Inspection planning', role: 'material', help: 'From the fatigue or corrosion assessment of the detail.' },
    { key: 'horizon', label: 'Planning horizon', unit: 'yr', type: 'int', default: 25, min: 2, max: 60, group: 'Inspection planning', role: 'bc', bind: 'life.designYears' },
    { key: 'rate', label: 'Real discount rate', unit: '%/yr', type: 'number', default: 6, min: 0, max: 40, group: 'Inspection planning', role: 'bc', help: PLACEHOLDER },
  ],
  compute(inp) {
    const subs = tableRows(inp.subsystems, 'Production system', { mtbf: { min: 1e-6 }, mttr: { min: 0 }, n: { min: 1, max: 12 }, k: { min: 1, max: 12 }, repairCost: { min: 0 } });
    uniqueNames(subs, 'Production system');
    for (const s of subs) { need(Number.isInteger(s.n) && Number.isInteger(s.k), `Subsystem '${s.name}': n and k must be whole numbers.`); need(s.k <= s.n, `Subsystem '${s.name}': ${s.k} units are required but only ${s.n} are installed — the system can never work.`); }
    const bars = tableRows(inp.barriers, 'Barriers', { pFail: { min: 0, max: 100 }, cost: { min: 0 } }, { min: 0 });
    need(inp.cFail > inp.cPrev, 'The unplanned failure cost must exceed the planned replacement cost; otherwise preventive replacement can never pay and the policy question does not arise.');
    const marginH = inp.unitMW * inp.price; // currency per hour of one unit down
    const cp = inp.cPrev + inp.dtPrev * marginH; const cf = inp.cFail + inp.dtFail * marginH;

    // (1) component life and replacement policy
    const w1 = weibullLife({ t: inp.scale, shape: inp.shape, scale: inp.scale });
    const age = run(() => ageReplacement({ shape: inp.shape, scale: inp.scale, costPreventive: cp, costFailure: cf, maxAge: 4 * inp.scale }), 'Age replacement');
    const Topt = age.optimalAge;
    const block = Topt ? run(() => blockReplacement({ interval: Topt, shape: inp.shape, scale: inp.scale, costPreventive: cp, costFailure: cf }), 'Block replacement') : null;
    const nh = run(() => nhppPowerLaw({ lambda: inp.nhppLambda, beta: inp.nhppBeta, t: inp.horizon, minimalRepairCost: inp.cMinimal, replacementCost: inp.cOverhaul }), 'Repairable system');
    // (2) system availability: block diagram (independent repair), Markov with limited crews, fault tree of the complement
    const blocks = subs.map((s) => {
      const u = unitAvailability({ mtbf: s.mtbf, mttr: s.mttr });
      const lambda = 1 / s.mtbf; const mu = s.mttr > 0 ? 1 / s.mttr : 1e9;
      const aRbd = rbdReliability(s.n === 1 ? u.availability : { type: 'koon', k: s.k, blocks: new Array(s.n).fill(u.availability) });
      const aMk = s.mttr > 0 ? markovKofN(lambda, mu, s.n, s.k, Math.round(inp.crews)) : 1;
      const failures = s.n * u.failuresPerYear;
      return { ...s, aUnit: u.availability, aRbd, aMk, failures, repair: failures * s.repairCost, down: (1 - aMk) * 8760 };
    });
    const aRbd = rbdReliability({ type: 'series', blocks: blocks.map((b) => b.aRbd) });
    const aSys = blocks.reduce((p, b) => p * b.aMk, 1);
    const qFt = faultTreeProbability({ gate: 'OR', inputs: blocks.map((b) => (b.n === 1 ? 1 - b.aUnit : b.k === 1 ? { gate: 'AND', inputs: new Array(b.n).fill(1 - b.aUnit) } : 1 - b.aRbd)) });
    // semi-Markov cross-check of the least available single unit (up/down alternation with general sojourn times)
    const worst = blocks.slice().sort((a, b) => a.aMk - b.aMk)[0];
    const sm = semiMarkovSteadyState({ P: [[0, 1], [1, 0]], meanSojourn: [worst.mtbf, Math.max(worst.mttr, 1e-9)], upStates: [0] });
    const dt = downtimeCost({ availability: aSys, capacityPerHour: inp.systemMW, price: inp.price });
    const repairAnnual = blocks.reduce((s, b) => s + b.repair, 0);
    // (3) event tree and insurance
    const et = inp.initFreq > 0 && bars.length ? eventTree({ frequency: inp.initFreq, barriers: bars.map((b) => ({ name: b.name, pFail: b.pFail / 100 })), cost: (failed) => inp.baseEventCost + failed.reduce((s, f, i) => s + (f ? bars[i].cost : 0), 0) }) : { sequences: [], expectedAnnualLoss: 0 };
    let insured = 0; let retained = 0;
    for (const s of et.sequences) { const r = insuranceRecovery({ loss: s.cost, deductible: inp.deductible, limit: inp.insLimit }); insured += s.frequency * r.insured; retained += s.frequency * r.retained; }
    const worstSeq = et.sequences.reduce((a, s) => (s.cost > (a?.cost ?? -1) ? s : a), null);
    // (4) inspection interval
    const cand = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15].filter((v) => v <= inp.horizon);
    const rbi = run(() => rbiEconomics({ intervals: cand, horizon: Math.round(inp.horizon), inspectionCost: inp.inspCost, failureCost: inp.inspFailCost, effectiveness: inp.inspEff, rate: inp.rate / 100, weibull: { shape: inp.inspShape, scale: inp.inspScale } }), 'Inspection planning');
    // fleet-level annual figures for the wearing component under the chosen policy
    const rateOpt = age.costRate; const rateRtf = age.runToFailureCostRate;
    const cof = costOfFailure({ repairCost: inp.cFail, downtimeHours: inp.dtFail, productionPerHour: inp.unitMW, marginPerUnit: inp.price, annualFrequency: 1 / w1.mttf });
    const correctiveAnnual = repairAnnual + inp.units * (Topt ? rateOpt : rateRtf);
    const targ = inp.targetAvail / 100;

    const kpis = [
      { key: 'mttf', label: 'Mean time to failure of the component', value: round(w1.mttf), unit: 'yr', explain: 'Mean of the Weibull life, η Γ(1 + 1/β). With β > 1 failures cluster around this age instead of being spread evenly, which is what makes planned replacement worthwhile.' },
      { key: 'optAge', label: 'Cost-optimal preventive replacement age', value: round(Topt ?? 0), unit: 'yr', explain: Topt ? 'Replacement age that minimises long-run cost per year, counting lost production. Replace earlier and good life is wasted; later and expensive failures dominate.' : 'Zero: no finite replacement age beats running to failure, because the hazard does not increase enough (β ≤ 1) or the failure premium is too small.' },
      { key: 'costRateOpt', label: 'Cost rate at the optimum, per component', value: k$(rateOpt), unit: `${KU}/yr`, explain: 'Long-run annual cost of the age-replacement policy for one component, including the margin lost during planned and unplanned downtime.' },
      { key: 'costRateRtf', label: 'Cost rate if run to failure, per component', value: k$(rateRtf), unit: `${KU}/yr`, explain: 'Annual cost with no preventive replacement: one failure cost per mean life.' },
      { key: 'policySaving', label: 'Fleet saving of the optimal policy', value: mil(inp.units * (rateRtf - rateOpt)), unit: `${MU}/yr`, explain: `Difference between run-to-failure and optimal age replacement over ${inp.units} components: the value of the maintenance strategy.` },
      { key: 'blockRate', label: 'Cost rate of block replacement at the same interval', value: k$(block ? block.costRate : rateRtf), unit: `${KU}/yr`, explain: 'Replacing every unit at fixed calendar intervals regardless of age. Simpler to plan offshore than age replacement, at a somewhat higher cost because recently replaced units are replaced again.' },
      { key: 'nhppFailures', label: `Expected failures of the repairable system in ${inp.horizon} years`, value: round(nh.expectedFailures), unit: '-', explain: 'Power-law (Crow–AMSAA) count of failures under minimal repair. An exponent above 1 means the failure intensity grows with age.' },
      { key: 'overhaulAge', label: 'Cost-optimal overhaul age (minimal repair)', value: round(nh.optimalReplacementAge ?? 0), unit: 'yr', explain: nh.optimalReplacementAge ? 'Age at which the growing cost of minimal repairs justifies renewing the system: the mid-life refurbishment date.' : 'Zero: the failure intensity does not grow (β ≤ 1), so overhaul on age never pays.' },
      { key: 'availability', label: 'System production availability', value: round(100 * aSys, 6), unit: '%', limit: inp.targetAvail, sense: 'min', explain: `Steady-state fraction of time the series system can produce, with ${inp.crews} repair crew(s). This is the availability to carry into the energy and revenue forecasts.` },
      { key: 'availRbd', label: 'Availability by block diagram (unlimited crews)', value: round(100 * aRbd, 6), unit: '%', explain: 'The same system assuming every failed unit is repaired independently. The gap to the Markov figure is the cost of repairs queuing for a crew.' },
      { key: 'unavailFt', label: 'Top-event probability by fault tree', value: round(100 * qFt, 6), unit: '%', explain: 'Probability of "no production" from the fault tree of the same system (OR of subsystem failures, AND across redundant units). It must equal one minus the block-diagram availability.' },
      { key: 'downHours', label: 'Expected downtime', value: round(dt.downtimeHours), unit: 'h/yr', explain: 'Equivalent full-system hours lost per year at the Markov availability.' },
      { key: 'downCost', label: 'Expected downtime cost', value: mil(dt.lostMargin), unit: `${MU}/yr`, explain: 'Lost production valued at the net margin: maintenance-related revenue loss. It is usually larger than the repair bill and is what redundancy and faster access buy back.' },
      { key: 'repairCost', label: 'Expected corrective repair cost', value: mil(repairAnnual), unit: `${MU}/yr`, explain: 'Failure frequency × repair cost summed over all installed units of the system table.' },
      { key: 'eventLoss', label: 'Expected annual loss from the event tree', value: mil(et.expectedAnnualLoss ?? 0), unit: `${MU}/yr`, explain: 'Initiating frequency × probability × consequence summed over all barrier outcomes: the cost-of-failure estimate for the escalating event.' },
      { key: 'retained', label: 'Retained (uninsured) part of that loss', value: mil(retained), unit: `${MU}/yr`, explain: `Expected loss falling on the owner after the deductible and above the limit; ${fM(insured)}/yr is expected to be recovered from insurers.` },
      { key: 'inspInterval', label: 'Lowest-cost inspection interval', value: rbi.best.interval, unit: 'yr', explain: 'Interval minimising the present value of inspection cost plus expected failure cost over the horizon. Shorter intervals buy little extra risk reduction; longer ones let the failure probability grow.' },
      { key: 'inspPV', label: 'PV of inspection and risk at that interval', value: mil(rbi.best.pvTotal), unit: MU, explain: `Discounted cost of the inspection programme plus residual failure risk, against ${fM(rbi.noInspection.pvRisk)} of risk with no inspection.` },
    ];
    const checks = [
      checkMin('Production availability ≥ requirement', 100 * aSys, inp.targetAvail, '%', 'ISO 20815:2018 (production assurance target)', 'Guards against a production forecast built on an availability the system cannot deliver; below the requirement every revenue figure downstream is overstated.'),
      check('Retained event loss ≤ 10 % of expected downtime and repair cost', retained, Math.max(0.1 * (dt.lostMargin + repairAnnual), 1e-9), `${CUR}/yr`, 'Risk retention review (insurance programme design)', 'Guards against an insurance structure that leaves a disproportionate share of rare, severe losses on the owner; above 1.0 review the deductible and limit.'),
    ];
    const ages = linspace(inp.scale * 0.02, inp.scale * 2.2, 80);
    const sg = Math.sqrt(Math.log(1 + 0.4 ** 2)); // lognormal σ for a coefficient of variation of 0.4
    const plots = [
      { type: 'line', title: 'Cost rate against preventive replacement age', x: { label: 'Replacement age', unit: 'yr' }, y: { label: 'Cost per component', unit: `${KU}/yr` },
        series: [{ name: 'Age replacement', x: age.curve.age.map((v) => round(v, 4)), y: age.curve.costRate.map(k$) }, { name: 'Run to failure', x: [age.curve.age[0], age.curve.age[age.curve.age.length - 1]], y: [k$(rateRtf), k$(rateRtf)], dash: true }, ...(Topt ? [{ name: 'Optimum', x: [round(Topt, 4)], y: [k$(rateOpt)] }] : [])],
        explain: 'Long-run annual cost for each replacement age. The minimum is the optimal policy; a flat bottom means the campaign can be timed to suit weather and vessel availability at little penalty.' },
      { type: 'line', title: 'Survival of the component', x: { label: 'Age', unit: 'yr' }, y: { label: 'Probability of surviving', unit: '-' },
        series: [{ name: `Weibull β = ${inp.shape}`, x: ages.map((v) => round(v, 4)), y: ages.map((t) => round(weibullLife({ t, shape: inp.shape, scale: inp.scale }).reliability, 5)) }, { name: 'Exponential, same MTTF', x: ages.map((v) => round(v, 4)), y: ages.map((t) => round(exponentialLife({ t, rate: 1 / w1.mttf }).reliability, 5)), dash: true },
          { name: 'Lognormal, same median, CoV 0.4', x: ages.map((v) => round(v, 4)), y: ages.map((t) => round(lognormalLife({ t, mu: Math.log(inp.scale * Math.LN2 ** (1 / inp.shape)), sigma: sg }).reliability, 5)), dash: true }],
        explain: 'Survival curve of the chosen Weibull life beside an exponential with the same mean and a lognormal with the same median. The exponential has no memory, so it can never justify preventive replacement.' },
      { type: 'bar', title: 'Unavailability by subsystem', categories: blocks.map((b) => b.name), x: { label: 'Subsystem' }, y: { label: 'Unavailability', unit: 'h/yr' },
        series: [{ name: `Markov, ${inp.crews} crew(s)`, y: blocks.map((b) => round(b.down, 5)) }, { name: 'Block diagram (independent repair)', y: blocks.map((b) => round((1 - b.aRbd) * 8760, 5)) }],
        explain: 'Hours per year each subsystem keeps the whole system down. The tallest bar is where redundancy, spares or faster access pays; a gap between the two series shows repair crews are a bottleneck.' },
      { type: 'line', title: 'Inspection interval: cost against risk', x: { label: 'Inspection interval', unit: 'yr' }, y: { label: 'Present value', unit: MU },
        series: [{ name: 'Inspection cost', x: rbi.results.map((r) => r.interval), y: rbi.results.map((r) => mil(r.pvInspection)) }, { name: 'Failure risk', x: rbi.results.map((r) => r.interval), y: rbi.results.map((r) => mil(r.pvRisk)) }, { name: 'Total', x: rbi.results.map((r) => r.interval), y: rbi.results.map((r) => mil(r.pvTotal)) }],
        explain: 'Inspection cost falls and failure risk rises as the interval lengthens; the minimum of the total is the economic interval. It must still satisfy any code-mandated maximum.' },
    ];
    const tables = [
      { title: 'Subsystem availability and cost', columns: [{ key: 'name', label: 'Subsystem' }, { key: 'cfg', label: 'Redundancy' }, { key: 'au', label: 'Unit availability', unit: '%' }, { key: 'am', label: 'Block availability (Markov)', unit: '%' }, { key: 'ar', label: 'Block availability (RBD)', unit: '%' }, { key: 'f', label: 'Failures', unit: '1/yr' }, { key: 'rc', label: 'Repair cost', unit: `${KU}/yr` }, { key: 'dc', label: 'Downtime cost', unit: `${KU}/yr` }],
        rows: blocks.map((b) => ({ name: b.name, cfg: `${b.k} of ${b.n}`, au: round(100 * b.aUnit, 6), am: round(100 * b.aMk, 6), ar: round(100 * b.aRbd, 6), f: round(b.failures, 4), rc: k$(b.repair), dc: k$(b.down * inp.systemMW * inp.price) })),
        explain: 'Each subsystem with its redundancy, availability by both methods, failure frequency and annual repair and downtime cost. Rank improvement options on the last two columns.' },
      { title: 'Event-tree sequences', columns: [{ key: 'label', label: 'Sequence' }, { key: 'f', label: 'Frequency', unit: '1/yr' }, { key: 'c', label: 'Consequence', unit: MU }, { key: 'r', label: 'Risk', unit: `${KU}/yr` }],
        rows: et.sequences.map((s) => ({ label: s.label, f: round(s.frequency, 5), c: mil(s.cost), r: k$(s.frequency * s.cost) })), explain: 'Every combination of barrier success and failure after the initiating event, with its frequency and consequence. Low-frequency, high-consequence rows drive the insurance need.' },
      { title: 'Inspection interval comparison', columns: [{ key: 'i', label: 'Interval', unit: 'yr' }, { key: 'pi', label: 'PV inspection', unit: MU }, { key: 'pr', label: 'PV failure risk', unit: MU }, { key: 'pt', label: 'PV total', unit: MU }, { key: 'ef', label: 'Expected failures' }],
        rows: rbi.results.map((r) => ({ i: r.interval, pi: mil(r.pvInspection), pr: mil(r.pvRisk), pt: mil(r.pvTotal), ef: round(r.expectedFailures, 4) })), explain: 'Present value of inspection cost and of expected failure cost for each candidate interval; choose the lowest total that also meets the mandatory maximum interval.' },
    ];
    const bn = blocks.slice().sort((a, b) => b.down - a.down)[0];
    const gap = aRbd - aSys;
    const advice = [
      adv('engineering', aSys >= targ ? 'info' : 'warn', `System availability is ${pct(aSys, 5)} against ${inp.targetAvail} % required; "${bn.name}" contributes ${fmt(bn.down)} of the ${fmt(dt.downtimeHours)} h/yr lost.${gap > 1e-5 ? ` Limited repair crews cost ${fmt(gap * 8760)} h/yr compared with independent repair.` : ''} The wearing component (β = ${inp.shape}, η = ${inp.scale} yr) is best replaced at ${Topt ? `${fmt(Topt)} years` : 'failure only'}.`,
        aSys >= targ ? `Availability has ${fmt((aSys - targ) * 8760)} h/yr of margin. Confirm MTBF and MTTR with ISO 14224 field data: repair time, driven by access and weather, usually matters more than failure rate offshore.` : `Recover ${fmt((targ - aSys) * 8760)} h/yr: shorten the repair time of "${bn.name}" (spares held, faster vessel access), add redundancy there, or add a repair crew.`,
        'Production assurance, ISO 20815:2018; reliability-centred maintenance, IEC 60300-3-11:2009; reliability data per ISO 14224:2016.',
        `Unplanned repairs need unscheduled vessel mobilisations: with about ${fmt(blocks.reduce((s, b) => s + b.failures, 0))} system failures a year at, say, 2 vessel days each, that is ${fmt(2 * blocks.reduce((s, b) => s + b.failures, 0) * MGO_T_PER_DAY * CO2_MGO)} t CO₂ a year (${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂/t, screening). Planned campaigns that bundle replacements cut mobilisations, and replacing parts before failure avoids secondary damage and scrapped components.`),
      adv('economic', Topt ? 'ok' : 'info', `Optimal age replacement costs ${fmt(rateOpt / 1e3)} ${KU}/yr per component against ${fmt(rateRtf / 1e3)} run to failure: ${fM(inp.units * (rateRtf - rateOpt))}/yr over the fleet. Downtime costs ${fM(dt.lostMargin)}/yr and corrective repairs ${fM(repairAnnual)}/yr; the event tree adds ${fM(et.expectedAnnualLoss ?? 0)}/yr of which ${fM(retained)} is retained.`,
        `Adopt ${Topt ? `replacement at about ${fmt(Topt)} years (block replacement at the same interval costs ${fmt((block.costRate / rateOpt - 1) * 100)} % more but is easier to plan)` : 'condition-based rather than age-based replacement'}, inspect every ${rbi.best.interval} years (saving ${fM(rbi.noInspection.pvRisk - rbi.best.pvTotal)} in present value against no inspection), and spend up to ${fM(dt.lostMargin / Math.max(dt.downtimeHours, 1e-9))} per hour of downtime avoided.`,
        'Optimum preventive maintenance, Barlow & Hunter (1960); risk-based inspection economics, DNV-RP-G101 (2021) and API RP 580 (2016); whole-life costing per ISO 15686-5:2017.'),
    ];
    return { kpis, checks, plots, tables, advice,
      notes: [`Replacement costs include lost margin: planned ${fM(cp)} and unplanned ${fM(cf)} per event (ratio ${fmt(cf / cp)}).`, `Semi-Markov cross-check of the least available unit ("${worst.name}"): ${pct(sm.availability, 6)}, equal to MTBF/(MTBF + MTTR) whatever the shape of the repair-time distribution.`,
        `Single-failure cost of the wearing component: ${fM(cof.costPerEvent)} (repair ${fM(cof.repair)}, lost production ${fM(cof.lostProduction)}).`, 'Failure rates, repair times and costs are indicative placeholders.'],
      publish: { 'econ.availability': aSys, 'econ.correctiveAnnual': correctiveAnnual, 'econ.downtimeCostAnnual': dt.lostMargin, 'econ.optimalReplacementAge': Topt ?? 0, 'econ.inspectionInterval': rbi.best.interval } };
  },
  verify: [
    // single unit, MTBF 1000 h, MTTR 10 h: A = 1000/1010 = 99.0099 %; downtime = 8760 × 10/1010 = 86.733 h/yr
    { name: 'Two-state Markov availability A = MTBF/(MTBF + MTTR)', source: 'Rausand & Høyland, System Reliability Theory, 2nd ed., eq. (8.14): A = μ/(λ + μ)',
      inputs: { subsystems: [{ name: 'U', mtbf: 1000, mttr: 10, n: 1, k: 1, repairCost: 0 }] }, expect: [{ kpi: 'availability', value: 99.0099, tol: 1e-6 }, { kpi: 'availRbd', value: 99.0099, tol: 1e-6 }, { kpi: 'unavailFt', value: 0.990099, tol: 1e-5 }, { kpi: 'downHours', value: 86.7327, tol: 1e-5 }] },
    // 1-of-2 redundant pair, two crews (independent repair): A = 1 − (10/1010)² = 99.990197 %;
    // with one crew (birth–death): P2/P0 = (2λ/μ)(λ/μ) = 2 × 0.01 × 0.01 = 2e-4, P1/P0 = 0.02 → P2 = 2e-4/1.0202 → A = 99.980396 %
    { name: 'Redundant pair with independent repair equals the block diagram', source: 'Parallel structure: A = 1 − (1 − A_unit)² = 1 − 0.00990099² = 0.99990197',
      inputs: { subsystems: [{ name: 'P', mtbf: 1000, mttr: 10, n: 2, k: 1, repairCost: 0 }], crews: 2 }, expect: [{ kpi: 'availability', value: 99.990197, tol: 1e-7 }, { kpi: 'availRbd', value: 99.990197, tol: 1e-7 }] },
    { name: 'Redundant pair with a single repair crew (birth–death chain)', source: 'Hand solution of the three-state chain: P2 = 2ρ²/(1 + 2ρ + 2ρ²), ρ = λ/μ = 0.01 → A = 1 − 1.96040e-4',
      inputs: { subsystems: [{ name: 'P', mtbf: 1000, mttr: 10, n: 2, k: 1, repairCost: 0 }], crews: 1 }, expect: [{ kpi: 'availability', value: 99.980396, tol: 1e-7 }] },
    // age replacement, β = 2.5, η = 10, cp = 1, cf = 10 (no downtime margin): numerical minimisation of [cp R(T) + cf F(T)]/∫R gives T* = 3.55 yr, C* = 0.4750;
    // run to failure: cf/MTTF = 10/(10 Γ(1.4)) = 10/8.8726 = 1.1271
    { name: 'Optimal age replacement, Weibull β = 2.5', source: 'Barlow & Hunter (1960) cost-rate function minimised by independent midpoint integration on a 0.01-yr grid: T* = 3.55, C(T*) = 0.47506; MTTF = η Γ(1.4) = 8.8726',
      inputs: { shape: 2.5, scale: 10, cPrev: 1000, cFail: 10000, dtPrev: 0, dtFail: 0, units: 1 }, expect: [{ kpi: 'optAge', value: 3.55, tol: 0.005 }, { kpi: 'costRateOpt', value: 0.47506, tol: 1e-3 }, { kpi: 'costRateRtf', value: 1.1271, tol: 1e-4 }, { kpi: 'mttf', value: 8.8726, tol: 1e-4 }] },
    // NHPP: N(10) = 0.1 × 10² = 10; T* = sqrt(cr/(cm λ (β − 1))) = sqrt(50/(1 × 0.1 × 1)) = 22.36 yr
    { name: 'Power-law process: expected failures and optimal overhaul age', source: 'Crow (1974) N(t) = λ t^β; T* = [c_r/(c_m λ (β − 1))]^(1/β) = √500 = 22.3607',
      inputs: { nhppLambda: 0.1, nhppBeta: 2, cMinimal: 1, cOverhaul: 50, horizon: 10 }, expect: [{ kpi: 'nhppFailures', value: 10, tol: 1e-9 }, { kpi: 'overhaulAge', value: 22.3607, tol: 1e-5 }] },
  ],
};
