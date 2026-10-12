// Asset replacement economics: economic service life of a new asset, defender–challenger analysis of an existing
// one, repair versus replace, and the keep/replace plan over a planning horizon by dynamic programming.
import { economicServiceLife, replacementAnalysis, repairVsReplace } from '../../engine/econ/lifecycle.js';
import { equipmentReplacementDp } from '../../engine/econ/optimise.js';
import { CUR, CUR_NOTE, PLACEHOLDER, CO2_STEEL, need, run, round, fmt, pct, adv, check, cost, tableRows, range } from './_common.js';

const k$ = (v) => round(v / 1e3, 6);
const KU = `k${CUR}`;
const CHALLENGER = [[1, 1.9e5, 2.9e6], [2, 2.3e5, 2.4e6], [3, 2.8e5, 2.0e6], [4, 3.5e5, 1.65e6], [5, 4.4e5, 1.35e6], [6, 5.6e5, 1.1e6], [7, 7.1e5, 9.0e5], [8, 9.0e5, 7.2e5], [9, 1.14e6, 5.8e5], [10, 1.45e6, 4.6e5], [11, 1.83e6, 3.6e5], [12, 2.3e6, 2.8e5]]
  .map(([year, om, salvage]) => ({ year, om, salvage }));
const DEFENDER = [[1, 8.2e5, 5.6e5], [2, 1.02e6, 4.2e5], [3, 1.28e6, 3.0e5], [4, 1.62e6, 2.0e5], [5, 2.05e6, 1.2e5], [6, 2.6e6, 6.0e4]].map(([year, om, salvage]) => ({ year, om, salvage }));

function profile(rows, what) {
  const r = tableRows(rows, what, { year: { min: 1 }, om: { min: 0 }, salvage: {} }, { min: 2 }).sort((a, b) => a.year - b.year);
  r.forEach((x, i) => need(x.year === i + 1, `${what}: years must run 1, 2, 3 … without gaps (row ${i + 1} is year ${x.year}).`));
  return r;
}

export default {
  id: 's17-replacement-service-life', suite: 17, title: 'Asset replacement and economic service life',
  summary: 'Finds the retention period that minimises the equivalent annual cost of a new asset, decides whether an existing asset should be kept, repaired or replaced now, and lays out the replacement plan over a planning horizon.',
  covers: [3], fidelity: 'analytical',
  standards: ['Blank & Tarquin, Engineering Economy, 8th ed., ch. 11 (replacement and retention decisions)', 'ISO 15686-5:2017 §5 (life-cycle costing of replacement options)', 'ISO 55000:2014 (asset management — value from assets)', 'Bellman (1955), equipment replacement by dynamic programming'],
  models: ['Asset replacement analysis (defender versus challenger)', 'Optimal economic service life', 'Equivalent annual cost and capital recovery', 'Marginal cost of retention', 'Repair-versus-replace comparison', 'Dynamic programming keep/replace plan', 'Residual value estimation'],
  limits: `Deterministic costs in constant (real) money, like-for-like replacement under the repeatability assumption, year-end cash flows and no tax. Technical obsolescence, production differences between old and new assets and failure consequences are not modelled (see the reliability module). ${CUR_NOTE}`,
  publishes: ['econ.economicLife', 'econ.replacementEac'],
  searchTerms: 'economic service life replacement analysis defender challenger equivalent annual cost',
  inputs: [
    { key: 'rate', label: 'Real discount rate (MARR)', unit: '%/yr', type: 'number', default: 8, min: 0, max: 50, group: 'Financial constraints', role: 'bc', help: `Minimum attractive rate of return in real terms. ${PLACEHOLDER}` },
    cost({ key: 'capex', label: 'Installed cost of a new asset (challenger)', unit: CUR, default: 4.2e6, min: 1, max: 1e11, group: 'Challenger (new asset)', design: true, uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Supply, installation and commissioning of the replacement, e.g. a pump package, crane or cable section.' }),
    { key: 'challenger', label: 'New asset: operating cost and salvage by age', type: 'table', group: 'Challenger (new asset)', role: 'cost', help: `Operating and maintenance cost during each year of age and resale or scrap value at the end of it (may be negative when removal costs more than the scrap). ${PLACEHOLDER}`,
      columns: [{ key: 'year', label: 'Year of age', type: 'number' }, { key: 'om', label: 'O&M cost in the year', unit: CUR, type: 'number' }, { key: 'salvage', label: 'Salvage at year end', unit: CUR, type: 'number' }], default: CHALLENGER },
    { key: 'defAge', label: 'Present age of the existing asset', unit: 'yr', type: 'number', default: 9, min: 0, max: 100, group: 'Opening state (existing asset)', role: 'ic', help: 'Initial condition of the defender. Its original cost is sunk and does not enter the decision.' },
    cost({ key: 'defValue', label: 'Market value of the existing asset now', unit: CUR, default: 7.2e5, min: -1e10, max: 1e11, group: 'Opening state (existing asset)', help: 'What it would realise today (resale or scrap less removal): the capital tied up by keeping it.' }),
    { key: 'defender', label: 'Existing asset: operating cost and salvage by further year', type: 'table', group: 'Opening state (existing asset)', role: 'cost', help: `Forecast for each further year of retention. ${PLACEHOLDER}`,
      columns: [{ key: 'year', label: 'Further year', type: 'number' }, { key: 'om', label: 'O&M cost in the year', unit: CUR, type: 'number' }, { key: 'salvage', label: 'Salvage at year end', unit: CUR, type: 'number' }], default: DEFENDER },
    { key: 'techLife', label: 'Remaining technical life of the existing asset', unit: 'yr', type: 'number', default: 5, min: 0, max: 100, group: 'Opening state (existing asset)', role: 'bc', bind: 'life.fatigueYears',
      help: 'Engineering limit from fatigue, corrosion or wear analysis. Retention beyond it is not permitted whatever the economics say.' },
    cost({ key: 'repairCost', label: 'Cost of a life-extending repair', unit: CUR, default: 1.1e6, min: 0, max: 1e11, group: 'Repair option', design: true }),
    { key: 'repairLife', label: 'Extra life bought by the repair', unit: 'yr', type: 'number', default: 5, min: 0.5, max: 60, group: 'Repair option', role: 'operation' },
    cost({ key: 'repairOm', label: 'O&M after repair', unit: `${CUR}/yr`, default: 6.5e5, min: 0, max: 1e10, group: 'Repair option' }),
    { key: 'horizon', label: 'Planning horizon', unit: 'yr', type: 'int', default: 20, min: 2, max: 60, group: 'Planning', role: 'bc', bind: 'life.designYears', help: 'Remaining field or facility life over which the service is needed.' },
  ],
  compute(inp) {
    const ch = profile(inp.challenger, 'New asset profile'); const df = profile(inp.defender, 'Existing asset profile');
    const i = inp.rate / 100;
    const esl = run(() => economicServiceLife({ capex: inp.capex, om: ch.map((r) => r.om), salvage: ch.map((r) => r.salvage), rate: i }), 'Economic service life');
    const rep = run(() => replacementAnalysis({ defender: { marketValue: inp.defValue, om: df.map((r) => r.om), salvage: df.map((r) => r.salvage) }, challenger: { capex: inp.capex, om: ch.map((r) => r.om), salvage: ch.map((r) => r.salvage) }, rate: i }), 'Replacement analysis');
    const life = esl.economicLife;
    const rr = run(() => repairVsReplace({ repairCost: inp.repairCost, repairLife: inp.repairLife, omOld: inp.repairOm, replaceCost: inp.capex, newLife: life, omNew: esl.omEquivalent[life - 1], salvageOld: Math.max(inp.defValue, 0), salvageNew: ch[life - 1].salvage, rate: i }), 'Repair versus replace');
    // dynamic programme over the horizon, starting from a new asset and replacing like for like
    const maxAge = ch.length;
    const omAt = (age) => ch[Math.min(Math.max(Math.round(age), 0), maxAge - 1)].om;
    const svAt = (age) => (age <= 0 ? inp.capex : ch[Math.min(Math.round(age), maxAge) - 1].salvage);
    const dp = run(() => equipmentReplacementDp({ horizon: Math.round(inp.horizon), maxAge, price: inp.capex, om: omAt, salvage: svAt, rate: i, initialAge: 0 }), 'Replacement plan');
    const cycle = dp.replacements.length ? dp.replacements.map((y, k) => (k === 0 ? y - 1 : y - dp.replacements[k - 1])) : [];
    const keep = Math.min(rep.keepYears, Math.floor(inp.techLife + 1e-9));
    const techGoverns = rep.keepYears > inp.techLife + 1e-9;
    const replaceNow = rep.replaceNow || inp.techLife < 1;
    const flat = (esl.eac[Math.min(life, maxAge - 1)] ?? esl.minEac) / esl.minEac - 1;

    const kpis = [
      { key: 'economicLife', label: 'Economic service life of a new asset', value: life, unit: 'yr', explain: 'Retention period that minimises equivalent annual cost: before it, capital recovery dominates; after it, rising operating cost does. Replace like-for-like at this age unless the technical life is shorter.' },
      { key: 'minEac', label: 'Minimum equivalent annual cost (new asset)', value: k$(esl.minEac), unit: `${KU}/yr`, explain: 'Annual cost of owning and operating a new asset when kept for its economic life. This is the benchmark the existing asset must beat each year.' },
      { key: 'defMarginal', label: 'Cost of keeping the existing asset one more year', value: k$(rep.defenderMarginalCost[0]), unit: `${KU}/yr`, limit: k$(esl.minEac), sense: 'max', explain: 'Loss of market value, interest on the capital tied up and next year’s operating cost. When this exceeds the new asset’s minimum annual cost, replacement is due.' },
      { key: 'keepYears', label: 'Further years to keep the existing asset', value: keep, unit: 'yr', explain: 'Economic retention capped by the remaining technical life. Zero means replace now.' },
      { key: 'replaceNow', label: 'Replace now (1 = yes)', value: replaceNow ? 1 : 0, unit: '-', explain: 'One when the existing asset’s next-year cost already exceeds the new asset’s minimum annual cost, or its technical life is exhausted.' },
      { key: 'eacRepair', label: 'Equivalent annual cost of repairing', value: k$(rr.eacRepair), unit: `${KU}/yr`, explain: 'Annual cost of the repair route: the repair cost recovered over the life it buys, the opportunity cost of not selling the old asset, and its operating cost.' },
      { key: 'eacReplace', label: 'Equivalent annual cost of replacing', value: k$(rr.eacReplace), unit: `${KU}/yr`, explain: 'Annual cost of the replace route at the new asset’s economic life; the lower of the two routes is preferred.' },
      { key: 'breakEvenRepair', label: 'Break-even repair cost', value: k$(rr.breakEvenRepairCost), unit: KU, explain: 'The most that is worth paying for the repair: above it, replacement has the lower annual cost.' },
      { key: 'planPV', label: `Present value of cost over ${Math.round(inp.horizon)} years (optimal plan)`, value: k$(dp.pvCost), unit: KU, explain: 'Discounted cost of the cheapest keep/replace sequence over the planning horizon starting with a new asset, final salvage credited. Useful as the budget for this service.' },
      { key: 'nReplacements', label: 'Replacements in the planning horizon', value: dp.replacements.length, unit: '-', explain: 'Number of renewals in the optimal plan. With a finite horizon the last cycles may be shorter or longer than the economic life.' },
    ];
    const checks = [
      check('Keeping cost ≤ new-asset annual cost', rep.defenderMarginalCost[0], esl.minEac, `${CUR}/yr`, 'Blank & Tarquin ch. 11 (marginal cost test)', 'Guards against running an asset whose next year costs more than the annual cost of a new one. Exceeding 1.0 means money is lost every year replacement is deferred.'),
      check('Planned retention within technical life', rep.keepYears, Math.max(inp.techLife, 1e-9), 'yr', 'Engineering life limit (fatigue, corrosion or wear analysis)', 'Guards against an economic retention period that outlasts the structural or mechanical life. Beyond 1.0 the economics ask for more years than the integrity assessment allows, so the technical limit governs.'),
    ];
    const ages = range(esl.eac.length, (k) => k + 1);
    const plots = [
      { type: 'line', title: 'Equivalent annual cost against retention period (new asset)', x: { label: 'Years retained', unit: 'yr' }, y: { label: 'Annual cost', unit: `${KU}/yr` },
        series: [{ name: 'Total EAC', x: ages, y: esl.eac.map(k$) }, { name: 'Capital recovery', x: ages, y: esl.capitalRecovery.map(k$) }, { name: 'Equivalent O&M', x: ages, y: esl.omEquivalent.map(k$) }, { name: 'Economic life', x: [life, life], y: [0, k$(esl.minEac)], dash: true }],
        explain: 'Falling capital recovery and rising equivalent operating cost add to a U-shaped total; its minimum is the economic service life. A flat bottom means the exact year matters little.' },
      { type: 'line', title: 'Existing asset: cost of each further year', x: { label: 'Further year', unit: 'yr' }, y: { label: 'Annual cost', unit: `${KU}/yr` },
        series: [{ name: 'Marginal cost of keeping', x: range(rep.defenderMarginalCost.length, (k) => k + 1), y: rep.defenderMarginalCost.map(k$) }, { name: 'New asset minimum EAC', x: [1, rep.defenderMarginalCost.length], y: [k$(esl.minEac), k$(esl.minEac)], dash: true },
          { name: 'Technical life limit', x: [inp.techLife, inp.techLife], y: [0, k$(Math.max(esl.minEac, ...rep.defenderMarginalCost))], dash: true }],
        explain: 'Replace in the first year where the marginal-cost line rises above the new-asset line, or at the technical-life limit if that comes first.' },
      { type: 'bar', title: 'Repair or replace', categories: ['Repair', 'Replace'], x: { label: 'Option' }, y: { label: 'Equivalent annual cost', unit: `${KU}/yr` }, series: [{ name: 'EAC', y: [k$(rr.eacRepair), k$(rr.eacReplace)] }],
        explain: 'Annual cost of the two routes on a like-for-like basis. The shorter bar wins; the break-even repair cost says how robust the choice is.' },
    ];
    const tables = [
      { title: 'New asset: retention analysis', columns: [{ key: 'n', label: 'Years retained' }, { key: 'om', label: 'O&M in year', unit: KU }, { key: 'sv', label: 'Salvage', unit: KU }, { key: 'cr', label: 'Capital recovery', unit: `${KU}/yr` }, { key: 'eo', label: 'Equivalent O&M', unit: `${KU}/yr` }, { key: 'eac', label: 'EAC', unit: `${KU}/yr` }],
        rows: ch.map((r, k) => ({ n: r.year, om: k$(r.om), sv: k$(r.salvage), cr: k$(esl.capitalRecovery[k]), eo: k$(esl.omEquivalent[k]), eac: k$(esl.eac[k]) })), explain: 'Equivalent annual cost for each possible retention period with its two components; the minimum row is the economic service life.' },
      { title: 'Optimal keep/replace plan', columns: [{ key: 'year', label: 'Year' }, { key: 'age', label: 'Age at start', unit: 'yr' }, { key: 'action', label: 'Action' }],
        rows: dp.plan.map((p) => ({ year: p.year, age: p.age, action: p.action })), explain: 'Year-by-year decision that minimises discounted cost over the planning horizon, starting from a new asset.' },
    ];
    const advice = [
      adv('economic', replaceNow ? 'warn' : 'ok', `New asset: economic life ${life} years at ${fmt(esl.minEac / 1e3)} ${KU}/yr. Existing asset: next year costs ${fmt(rep.defenderMarginalCost[0] / 1e3)} ${KU}. ${rr.decision === 'repair' ? 'Repair' : 'Replacement'} has the lower annual cost (${fmt(rr.eacRepair / 1e3)} vs ${fmt(rr.eacReplace / 1e3)} ${KU}/yr).`,
        replaceNow ? `Replace now: each year of deferral loses about ${fmt((rep.defenderMarginalCost[0] - esl.minEac) / 1e3)} ${KU}.` : `Keep the existing asset ${keep} more year(s), then replace${techGoverns ? ' (the technical life, not the economics, sets the date)' : ''}. Do not pay more than ${fmt(rr.breakEvenRepairCost / 1e3)} ${KU} for a repair giving ${inp.repairLife} years.`,
        'Defender–challenger analysis with sunk cost excluded, Blank & Tarquin ch. 11; whole-life costing per ISO 15686-5:2017.'),
      adv('engineering', techGoverns ? 'warn' : 'info', `The economic answer rests on how fast operating cost rises with age (from ${fmt(ch[0].om / 1e3)} to ${fmt(ch[ch.length - 1].om / 1e3)} ${KU}/yr over ${ch.length} years) and on a remaining technical life of ${fmt(inp.techLife)} years.`,
        techGoverns ? `Economics favour ${rep.keepYears} more years but only ${fmt(inp.techLife)} are technically available: either demonstrate more life (inspection, re-assessment of fatigue or corrosion margins) or plan the replacement at the technical limit.` : `Confirm the cost-growth curve with maintenance records; the EAC curve changes by ${pct(Math.abs(flat))} one year past the optimum, so the replacement date can be moved to suit a planned shutdown without material loss.`,
        'Asset management, ISO 55000:2014; life-extension assessment per NORSOK N-006 (2015) / ISO 19902 for structures.',
        `Extending life defers the manufacture of a replacement: for a 100 t steel package that is about ${fmt(100 * CO2_STEEL)} t CO₂e of embodied emissions postponed (screening factor ${CO2_STEEL} t CO₂e/t), set against the higher energy use and leak risk of ageing equipment. Recover the old unit for reuse or recycling rather than scrap-to-landfill.`),
    ];
    return { kpis, checks, plots, tables, advice,
      notes: [`Real terms at ${inp.rate} %/yr; the existing asset's original cost and book value are sunk and excluded.`, dp.replacements.length ? `Optimal plan over ${Math.round(inp.horizon)} years replaces at the start of year(s) ${dp.replacements.join(', ')} (cycle lengths ${cycle.join(', ')} yr).` : `No replacement falls inside the ${Math.round(inp.horizon)}-year horizon.`, 'The keep/replace plan starts from a new asset; the existing asset’s own decision is given by the marginal-cost test.'],
      publish: { 'econ.economicLife': life, 'econ.replacementEac': esl.minEac } };
  },
  verify: [
    // P = 20 000; O&M 2000, 3000, 4500, 6500, 9000; salvage 14 000, 10 000, 7000, 4500, 2500; i = 10 %.
    // EAC(3) = [20 000 − 7000 (P/F,10%,3) + 2000 (P/F,1) + 3000 (P/F,2) + 4500 (P/F,3)] (A/P,10%,3)
    //        = [20 000 − 5259.20 + 1818.18 + 2479.34 + 3380.92] × 0.402115 = 22 419.24 × 0.402115 = 9015.1;
    // EAC(1) = 20 000 × 1.1 − 14 000 + 2000 = 10 000; EAC(2) = 9238.1, EAC(4) = 9162.7 → minimum at n = 3
    { name: 'Economic service life by tabulated EAC', source: 'Hand calculation with (P/F, 10 %, n) and (A/P, 10 %, n) factors, method of Blank & Tarquin ch. 11.2',
      inputs: { rate: 10, capex: 20000, challenger: [[1, 2000, 14000], [2, 3000, 10000], [3, 4500, 7000], [4, 6500, 4500], [5, 9000, 2500]].map(([year, om, salvage]) => ({ year, om, salvage })),
        defValue: 8000, defender: [[1, 3000, 6000], [2, 4000, 4500], [3, 5000, 3000]].map(([year, om, salvage]) => ({ year, om, salvage })), techLife: 10, horizon: 10 },
      // defender marginal cost year 1 = 8000 × 1.1 − 6000 + 3000 = 5800
      expect: [{ kpi: 'economicLife', value: 3, abs: 0 }, { kpi: 'minEac', value: 9.0151, tol: 1e-4 }, { kpi: 'defMarginal', value: 5.8, tol: 1e-9 }, { kpi: 'replaceNow', value: 0, abs: 0 }] },
  ],
};
