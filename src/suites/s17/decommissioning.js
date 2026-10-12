// Decommissioning cost and liability: bottom-up removal estimate with recycling credit, the asset retirement
// obligation and its accretion, funding of a decommissioning security, and the timing of removal.
import { decomEstimate, assetRetirementObligation, decomFunding, decomTiming, residualValue } from '../../engine/econ/lifecycle.js';
import { endOfLifeCredit, circularity } from '../../engine/econ/carbon.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, CO2_STEEL, CO2_MGO, MGO_T_PER_DAY, SCREEN, need, run, round, fmt, mil, fM, pct, adv, check, cost, tableRows, uniqueNames, range } from './_common.js';

const ITEMS = [
  { name: 'Structure and foundation removal', cost: 5.6e7 }, { name: 'Topsides / turbine removal', cost: 2.4e7 }, { name: 'Subsea equipment and cable recovery', cost: 1.5e7 },
  { name: 'Pipeline or cable abandonment in place', cost: 4.0e6 }, { name: 'Transport to shore', cost: 9.0e6 }, { name: 'Onshore dismantling and waste disposal', cost: 1.1e7 }, { name: 'Site clearance, restoration and verification survey', cost: 5.0e6 },
];

export default {
  id: 's17-decommissioning', suite: 17, title: 'Decommissioning cost, liability and funding',
  summary: 'Estimates the cost of removal and restoration, nets off the recycling credit, converts it to a discounted abandonment liability with its accretion, sizes the annual funding of a security and compares removal dates.',
  covers: [3], fidelity: 'analytical',
  standards: ['IAS 37 §45–47 and IFRIC 1 (provisions for decommissioning, present value and unwinding)', 'US GAAP ASC 410-20 (asset retirement obligations)', 'OSPAR Decision 98/3 and IMO Resolution A.672(16) (removal of offshore installations)', 'ISO 15686-5:2017 §4 (end-of-life costs in whole-life costing)', 'EN 15804:2012+A2:2019 module D (recycling credit)'],
  models: ['Decommissioning liability modelling', 'Bottom-up removal cost estimate with contingency', 'Asset retirement obligation (inflate, discount, accrete)', 'Sinking-fund and unit-of-production funding', 'Deferral (timing) valuation', 'Residual and scrap value estimation', 'End-of-life recycling credit and material circularity indicator'],
  limits: `Single removal campaign at a known date; cost today is inflated at one escalation rate and discounted at one rate. Regulatory derogations, residual liability for items left in place and tax relief on decommissioning spend are outside the model. Emission factors are ${SCREEN}. ${CUR_NOTE}`,
  publishes: ['econ.decomCost', 'econ.decomLiability', 'econ.decomFunding'],
  searchTerms: 'offshore decommissioning cost estimate asset retirement obligation provision',
  inputs: [
    { key: 'baseYear', label: 'Price base year of the estimate', unit: 'yr', type: 'int', default: 2026, min: 1990, max: 2100, group: 'Opening financial state', role: 'ic' },
    { key: 'fundOpening', label: 'Decommissioning fund already held', unit: CUR, type: 'number', default: 0, min: 0, max: 1e12, group: 'Opening financial state', role: 'ic', help: 'Opening balance of the escrow or security at the base year.' },
    { key: 'yearsToDecom', label: 'Years until decommissioning', unit: 'yr', type: 'int', default: 25, min: 1, max: 80, group: 'Terminal conditions', role: 'bc', bind: 'life.designYears', design: true, help: 'Terminal condition: remaining operating life from the design-life or fatigue-life assessment.' },
    { key: 'items', label: 'Decommissioning scope', type: 'table', group: 'Decommissioning costs', role: 'cost', help: `Cost of each activity in base-year money. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Activity', type: 'text' }, { key: 'cost', label: 'Cost today', unit: CUR, type: 'number' }], default: ITEMS },
    { key: 'wells', label: 'Wells to plug and abandon', unit: '-', type: 'int', default: 0, min: 0, max: 500, group: 'Decommissioning costs', role: 'operation' },
    cost({ key: 'costPerWell', label: 'Plug-and-abandon cost per well', unit: CUR, default: 8.0e6, min: 0, max: 1e9, group: 'Decommissioning costs' }),
    { key: 'monitorYears', label: 'Post-decommissioning monitoring period', unit: 'yr', type: 'int', default: 5, min: 0, max: 50, group: 'Decommissioning costs', role: 'operation' },
    cost({ key: 'monitorCost', label: 'Environmental monitoring cost', unit: `${CUR}/yr`, default: 4.0e5, min: 0, max: 1e9, group: 'Decommissioning costs' }),
    cost({ key: 'contingency', label: 'Contingency on gross cost', unit: '%', default: 20, min: 0, max: 200, group: 'Decommissioning costs', help: 'Decommissioning estimates made decades ahead carry wide uncertainty.' }),
    { key: 'vesselDays', label: 'Removal vessel days', unit: 'vessel-days', type: 'number', default: 150, min: 0, max: 5000, group: 'Decommissioning costs', role: 'operation', help: 'For the emissions of the removal campaign only; the vessel cost is inside the scope table.' },
    { key: 'steel', label: 'Steel to be removed', unit: 'tonnes', type: 'number', default: 62000, min: 0, max: 5e6, group: 'Recycling and circularity', role: 'geometry', bind: 'quantities.steel_t' },
    { key: 'recoverable', label: 'Fraction of steel recovered for recycling', unit: '-', type: 'number', default: 0.9, min: 0, max: 1, group: 'Recycling and circularity', role: 'operation', design: true },
    cost({ key: 'scrapPrice', label: 'Scrap steel price (delivered to recycler)', unit: `${CUR}/t`, default: 260, min: 0, max: 5000, group: 'Recycling and circularity', uncertain: { dist: 'lognormal', cov: 0.3 } }),
    { key: 'recycledContent', label: 'Recycled content of the original steel', unit: '-', type: 'number', default: 0.25, min: 0, max: 1, group: 'Recycling and circularity', role: 'material' },
    cost({ key: 'escalation', label: 'Escalation of decommissioning cost', unit: '%/yr', default: 2.5, min: -5, max: 20, group: 'Liability', help: 'Nominal growth of removal cost to the settlement date.' }),
    { key: 'discountRate', label: 'Discount rate for the liability (nominal)', unit: '%/yr', type: 'number', default: 4.5, min: 0, max: 30, group: 'Liability', role: 'bc', help: `Rate prescribed for the provision (risk-free or credit-adjusted), in nominal terms to match the escalated cost. ${PLACEHOLDER}` },
    { key: 'fundRate', label: 'Return earned by the fund (nominal)', unit: '%/yr', type: 'number', default: 3.0, min: 0, max: 20, group: 'Funding', role: 'bc', help: `Net return on the escrow account. ${PLACEHOLDER}` },
    { key: 'securityFrac', label: 'Required security as share of the future cost', unit: '%', type: 'number', default: 100, min: 0, max: 200, group: 'Funding', role: 'bc', help: 'Regulatory or lease requirement for financial security at the removal date.' },
    { key: 'maxDeferral', label: 'Longest deferral of removal considered', unit: 'yr', type: 'int', default: 8, min: 0, max: 30, group: 'Timing', role: 'bc', help: 'Regulatory limit on leaving the asset idle after cessation.' },
    cost({ key: 'holdingCost', label: 'Cost of holding the idle asset', unit: `${CUR}/yr`, default: 3.5e6, min: 0, max: 1e10, group: 'Timing', help: 'Monitoring, navigation aids, insurance and integrity maintenance after cessation (base-year money).' }),
  ],
  compute(inp) {
    const items = tableRows(inp.items, 'Decommissioning scope', { cost: { min: 0 } });
    uniqueNames(items, 'Decommissioning scope');
    const n = Math.round(inp.yearsToDecom); const f = inp.escalation / 100; const r = inp.discountRate / 100; const fr = inp.fundRate / 100;
    const lump = items.reduce((s, i) => s + i.cost, 0);
    const est = run(() => decomEstimate({ removal: lump, wells: { count: inp.wells, costPerWell: inp.costPerWell }, monitoring: { years: inp.monitorYears, costPerYear: inp.monitorCost },
      scrap: { mass_t: inp.steel, recoverableFraction: inp.recoverable, pricePerT: inp.scrapPrice }, contingency: inp.contingency / 100 }), 'Decommissioning estimate');
    need(est.gross > 0, 'The decommissioning scope has no cost: enter at least one activity.');
    const netToday = Math.max(est.net, 0);
    const aro = run(() => assetRetirementObligation({ costToday: netToday, inflation: f, years: n, discountRate: r }), 'Asset retirement obligation');
    const target = Math.max(aro.futureCost * inp.securityFrac / 100 - inp.fundOpening * (1 + fr) ** n, 0);
    const fund = target > 0 ? run(() => decomFunding({ target, years: n, fundRate: fr, method: 'level' }), 'Funding schedule') : { contributions: new Array(n).fill(0), balance: new Array(n).fill(0), totalContributed: 0, interestEarned: 0 };
    // timing is evaluated from the cessation date (cost and holding cost in money of that year), then discounted
    // to today, so that holding cost is counted only for the years the asset stands idle
    const g = (1 + f) ** n; const d0 = (1 + r) ** -n;
    const timC = run(() => decomTiming({ costToday: netToday * g, escalation: f, rate: r, holdingCostPerYear: inp.holdingCost * g, candidateYears: range(Math.round(inp.maxDeferral) + 1) }), 'Timing');
    const tim = { options: timC.options.map((o) => ({ year: n + o.year, pvCost: o.pvCost * d0 })), best: { year: n + timC.best.year, pvCost: timC.best.pvCost * d0 } };
    const scrap = residualValue({ method: 'scrap', mass_t: inp.steel, scrapPricePerT: inp.scrapPrice, recoverableFraction: inp.recoverable, removalCost: 0 });
    // end-of-life carbon: recycling credit (net-scrap approach) against the removal campaign's fuel
    const credit = inp.steel > 0 ? endOfLifeCredit({ mass_t: inp.steel, recoveryRate: inp.recoverable, recycledContent: inp.recycledContent, primaryFactor: 2.3, secondaryFactor: 0.5 }) : { credit_t: 0, netScrap_t: 0 };
    const fuelCO2 = inp.vesselDays * MGO_T_PER_DAY * CO2_MGO;
    const circ = inp.steel > 0 ? circularity({ mass_t: inp.steel, recycledFeedstock: inp.recycledContent, collectedForRecycling: inp.recoverable, recyclingEfficiency: 0.95 }) : { mci: 0 };
    const balOpen = range(n, (k) => inp.fundOpening * (1 + fr) ** (k + 1) + fund.balance[k]);

    const kpis = [
      { key: 'gross', label: 'Gross decommissioning cost (today’s money)', value: mil(est.gross), unit: MU, explain: `Removal, recovery, disposal, restoration, well abandonment and monitoring priced in ${inp.baseYear} money before contingency and recycling credit.` },
      { key: 'net', label: 'Net cost with contingency and scrap credit', value: mil(netToday), unit: MU, explain: 'Gross cost plus contingency less the recycling credit for recovered steel. This is the figure carried as the terminal cost in the project cash flow.' },
      { key: 'credit', label: 'Recycling credit', value: mil(est.recyclingCredit), unit: MU, explain: 'Recovered steel tonnage at the scrap price. It offsets only a small part of removal cost and moves with a volatile scrap market, so it should not be relied on to fund removal.' },
      { key: 'future', label: `Cost at settlement (${inp.baseYear + n}, money of the day)`, value: mil(aro.futureCost), unit: MU, explain: `Net cost escalated ${n} years at ${inp.escalation} %/yr: the cash that will actually be needed.` },
      { key: 'liability', label: 'Discounted abandonment liability today', value: mil(aro.presentValue), unit: MU, explain: `Future cost discounted at the nominal liability rate of ${inp.discountRate} %: the provision recognised today. A lower discount rate raises it.` },
      { key: 'accretion1', label: 'First-year unwinding of the discount', value: mil(aro.accretion[0].expense), unit: `${MU}/yr`, explain: 'Annual finance charge that grows the provision towards the settlement amount; it increases every year.' },
      { key: 'funding', label: 'Level annual contribution to the security', value: mil(fund.contributions[0]), unit: `${MU}/yr`, explain: `Sinking-fund payment that, with the fund's return of ${inp.fundRate} %, reaches the required security at the removal date. An operating cost in all but name.` },
      { key: 'fundingPerSteel', label: 'Net cost per tonne of steel removed', value: inp.steel > 0 ? round(netToday / inp.steel) : 0, unit: `${CUR}/t`, explain: 'Benchmark denominator: removal cost per tonne is set by weight per lift, water depth and cutting method rather than by total tonnage.' },
      { key: 'bestYear', label: 'Lowest-cost removal year (years from now)', value: tim.best.year, unit: 'yr', explain: 'Removal date with the lowest present value of escalated cost plus the cost of holding the idle asset, within the permitted deferral.' },
      { key: 'co2Credit', label: 'Net end-of-life carbon balance', value: round(credit.credit_t - fuelCO2), unit: 't CO₂e', explain: 'Recycling credit of the recovered steel beyond its recycled content, less the fuel emissions of the removal campaign. Positive means removal and recycling avoid more than they emit.' },
      { key: 'mci', label: 'Material circularity indicator (steel)', value: round(circ.mci, 4), unit: '-', explain: 'Ellen MacArthur Foundation indicator from 0.1 (fully linear) to 1 (fully circular) for the steel, from its recycled content in and its recovery out.' },
    ];
    const checks = [check('Security held at removal ≥ required', aro.futureCost * inp.securityFrac / 100, balOpen[n - 1] || 1e-9, CUR, 'Financial security for decommissioning (lease or regulatory requirement)',
      'Guards against reaching the removal date without the funds to pay for it; a shortfall falls on the owner’s balance sheet or, on default, on the state.')];
    const yrs = range(n + 1, (k) => inp.baseYear + k);
    const plots = [
      { type: 'line', title: 'Provision and fund build-up', x: { label: 'Year', unit: '' }, y: { label: 'Amount', unit: MU },
        series: [{ name: 'Provision (accreting liability)', x: yrs, y: [aro.presentValue, ...aro.accretion.map((a) => a.closing)].map(mil) }, { name: 'Fund balance', x: yrs, y: [inp.fundOpening, ...balOpen].map(mil) }, { name: 'Cost at settlement', x: [yrs[0], yrs[n]], y: [mil(aro.futureCost), mil(aro.futureCost)], dash: true }],
        explain: 'The liability unwinds upward to the settlement cost while contributions build the fund. Both should meet at the removal date; a gap at the right-hand end is unfunded liability.' },
      { type: 'bar', title: 'Decommissioning cost build-up', categories: [...items.map((i) => i.name), ...(inp.wells > 0 ? ['Well plug and abandonment'] : []), 'Monitoring', 'Contingency', 'Recycling credit'], x: { label: 'Item' }, y: { label: 'Cost', unit: MU },
        series: [{ name: 'Cost', y: [...items.map((i) => mil(i.cost)), ...(inp.wells > 0 ? [mil(inp.wells * inp.costPerWell)] : []), mil(inp.monitorYears * inp.monitorCost), mil(est.contingency), -mil(est.recyclingCredit)] }],
        explain: 'Each scope item with contingency and the (negative) recycling credit. Offshore removal dominates; the credit is small in comparison.' },
      { type: 'line', title: 'Present value of cost against removal date', x: { label: 'Years from now', unit: 'yr' }, y: { label: 'PV of removal and holding cost', unit: MU }, series: [{ name: 'PV of cost', x: tim.options.map((o) => o.year), y: tim.options.map((o) => mil(o.pvCost)) }],
        explain: 'Deferral lowers the present value of removal only if the discount rate exceeds escalation and the saving outweighs the cost of holding the idle asset. The minimum is the economic removal date.' },
    ];
    const tables = [
      { title: 'Liability accretion schedule', columns: [{ key: 'year', label: 'Year' }, { key: 'opening', label: 'Opening provision', unit: MU }, { key: 'expense', label: 'Unwinding of discount', unit: MU }, { key: 'closing', label: 'Closing provision', unit: MU }, { key: 'contrib', label: 'Fund contribution', unit: MU }, { key: 'fund', label: 'Fund balance', unit: MU }],
        rows: aro.accretion.map((a, k) => ({ year: inp.baseYear + a.year, opening: mil(a.opening), expense: mil(a.expense), closing: mil(a.closing), contrib: mil(fund.contributions[k]), fund: mil(balOpen[k]) })),
        explain: 'Year-by-year provision with its finance charge, beside the contributions and balance of the security. Use it for the forecast of decommissioning expenditure and liabilities.' },
      { title: 'Estimate summary', columns: [{ key: 'item', label: 'Item' }, { key: 'v', label: 'Amount', unit: MU }],
        rows: [{ item: 'Scope items', v: mil(lump) }, { item: 'Well plug and abandonment', v: mil(est.breakdown.wellPA) }, { item: 'Monitoring', v: mil(est.breakdown.monitoring) }, { item: 'Gross', v: mil(est.gross) }, { item: `Contingency (${inp.contingency} %)`, v: mil(est.contingency) }, { item: 'Recycling credit', v: -mil(est.recyclingCredit) }, { item: 'Net, today’s money', v: mil(netToday) }, { item: 'At settlement, money of the day', v: mil(aro.futureCost) }, { item: 'Present value (provision)', v: mil(aro.presentValue) }],
        explain: 'From scope to provision: gross cost, contingency, credit, escalation to the removal date and discounting back to today.' },
    ];
    const perDay = inp.vesselDays > 0 ? lump / inp.vesselDays : 0;
    const funded = balOpen[n - 1] >= aro.futureCost * inp.securityFrac / 100 * 0.999;
    const advice = [
      adv('economic', funded ? 'ok' : 'warn', `Net cost today ${fM(netToday)}; ${fM(aro.futureCost)} at settlement in ${inp.baseYear + n}; provision today ${fM(aro.presentValue)}. Funding the security takes ${fM(fund.contributions[0])} a year.`,
        `Carry the net cost as the terminal cash flow of the project and the contribution as an operating cost. A one-point fall in the liability discount rate raises the provision by ${fM(netToday * (1 + f) ** n / (1 + r - 0.01) ** n - aro.presentValue)}; ${tim.best.year > n ? `deferring removal to year ${tim.best.year} lowers the present value by ${fM(tim.options[0].pvCost - tim.best.pvCost)}` : 'deferral beyond the end of operations does not pay once holding costs are counted'}.`,
        'Provisions at present value with unwinding of the discount, IAS 37 §45–47 and IFRIC 1; end-of-life cost in whole-life costing, ISO 15686-5:2017 §4.'),
      adv('engineering', 'info', `Removal is driven by ${fmt(inp.steel)} t of steel and about ${fmt(inp.vesselDays)} vessel days${perDay > 0 ? ` (scope cost per vessel day ${fM(perDay)})` : ''}; net cost is ${inp.steel > 0 ? fmt(netToday / inp.steel) : '–'} ${CUR} per tonne removed.`,
        'Design for removal now: lifting points and cut lines sized for the removal vessel, piles cut below mudline with access for internal cutting, and connections that can be released subsea reduce vessel days decades later at almost no fabrication cost.',
        'OSPAR Decision 98/3 (complete removal as the base case); IMO Resolution A.672(16); DNV-RP-H102 (2004, reissued) for marine operations during removal.',
        `Recovering ${pct(inp.recoverable)} of the steel returns ${fmt(credit.netScrap_t)} t of net scrap and a credit of ${fmt(credit.credit_t)} t CO₂e (EN 15804 module D, primary 2.3 and secondary 0.5 t CO₂e/t as screening values), against ${fmt(fuelCO2)} t CO₂ from the removal vessels (${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂/t). Material circularity indicator ${fmt(circ.mci)}; reuse of sound members scores higher than remelting.`),
    ];
    return { kpis, checks, plots, tables, advice,
      notes: [`Nominal terms for the liability: cost escalated at ${inp.escalation} %/yr and discounted at ${inp.discountRate} %/yr nominal; the net cost today is in real ${inp.baseYear} money.`, `Scrap value of the recovered steel ${fM(scrap)}; it is credited in the net cost.`, 'All removal costs and the scrap price are indicative placeholders.'],
      publish: { 'econ.decomCost': netToday, 'econ.decomLiability': aro.presentValue, 'econ.decomFunding': fund.contributions[0] },
      quantities: { vessel_days: { decommissioning: inp.vesselDays } } };
  },
  verify: [
    // 100 M today, 2 %/yr for 10 years → 121.899 M; discounted at 5 % → 121.899 / 1.628895 = 74.836 M; first-year accretion 5 % × 74.836 = 3.7418 M;
    // sinking fund at 5 %: 121.899 × 0.05 / (1.05^10 − 1) = 121.899 × 0.0795046 = 9.6916 M/yr
    { name: 'Asset retirement obligation: 100 M, 2 % escalation, 5 % discount, 10 years', source: 'Hand calculation: 100 × 1.02^10 / 1.05^10 = 74.836; sinking-fund factor (A/F, 5 %, 10) = 0.0795046',
      inputs: { items: [{ name: 'Removal', cost: 1e8 }], wells: 0, monitorYears: 0, contingency: 0, steel: 0, escalation: 2, discountRate: 5, fundRate: 5, yearsToDecom: 10, securityFrac: 100, fundOpening: 0 },
      expect: [{ kpi: 'net', value: 100, tol: 1e-9 }, { kpi: 'future', value: 121.8994, tol: 1e-6 }, { kpi: 'liability', value: 74.8363, tol: 1e-5 }, { kpi: 'accretion1', value: 3.74181, tol: 1e-5 }, { kpi: 'funding', value: 9.69156, tol: 1e-5 }] },
    // gross = 100 + 2 × 20 + 3 × 1 = 143; contingency 10 % = 14.3; credit = 1000 t × 0.9 × 200 = 0.18 M → net 157.12 M
    { name: 'Estimate build-up with wells, monitoring, contingency and scrap credit', source: 'Hand calculation: (100 + 40 + 3) × 1.10 − 1000 × 0.9 × 200/1e6 = 157.12 M',
      inputs: { items: [{ name: 'Removal', cost: 1e8 }], wells: 2, costPerWell: 2e7, monitorYears: 3, monitorCost: 1e6, contingency: 10, steel: 1000, recoverable: 0.9, scrapPrice: 200 },
      expect: [{ kpi: 'gross', value: 143, tol: 1e-9 }, { kpi: 'net', value: 157.12, tol: 1e-9 }, { kpi: 'credit', value: 0.18, tol: 1e-9 }] },
  ],
};
