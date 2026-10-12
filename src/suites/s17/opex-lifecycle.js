// Operating expenditure build-up and whole-life cost: personnel, logistics, maintenance, insurance, energy demand,
// inspections and major replacements, rolled into annual OPEX, discounted life-cycle cost and equivalent annual cost.
import { opexBuildUp } from '../../engine/econ/opex.js';
import { lifeCycleCost, totalCostOfOwnership } from '../../engine/econ/lifecycle.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, CO2_MGO, SCREEN, need, run, round, fmt, mil, fM, pct, adv, check, cost, tableRows, uniqueNames, range, crf } from './_common.js';

const PERSONNEL = [
  { name: 'Site management and engineering', count: 14, costPerYear: 165000 }, { name: 'Turbine and BoP technicians', count: 36, costPerYear: 112000 },
  { name: 'Control room and marine coordination', count: 12, costPerYear: 105000 }, { name: 'HSE, warehouse and administration', count: 9, costPerYear: 92000 },
];
const LOGISTICS = [
  { name: 'Service operation vessel', days: 330, dayRate: 36000 }, { name: 'Crew transfer vessels (2)', days: 520, dayRate: 4200 },
  { name: 'Jack-up for major component exchange', days: 18, dayRate: 190000 }, { name: 'Helicopter support', days: 30, dayRate: 16000 },
];
const INSPECTIONS = [
  { name: 'Subsea structural and scour survey (ROV)', cost: 2.4e6, intervalYears: 2 }, { name: 'Cable route and burial survey', cost: 1.8e6, intervalYears: 3 },
  { name: 'Statutory lifting, electrical and safety inspection', cost: 9.0e5, intervalYears: 1 }, { name: 'Blade and tower inspection campaign', cost: 1.5e6, intervalYears: 2 },
];
const RENEWALS = [
  { name: 'Major component exchange programme (gearboxes / main bearings)', year: 12, cost: 4.8e7 }, { name: 'Export cable repair allowance', year: 9, cost: 1.4e7 },
  { name: 'Coating and cathodic protection retrofit', year: 15, cost: 1.1e7 }, { name: 'Control system and SCADA upgrade', year: 14, cost: 6.0e6 },
];

export default {
  id: 's17-opex-lifecycle', suite: 17, title: 'Operating cost build-up and whole-life cost',
  summary: 'Builds annual OPEX from personnel, logistics, maintenance, insurance, energy demand and inspections, adds major renewals and decommissioning, and reports undiscounted and discounted life-cycle cost and its annual equivalent.',
  covers: [3], fidelity: 'analytical',
  standards: ['ISO 15686-5:2017 (life-cycle costing)', 'ISO 15663:2021 (life-cycle costing, petroleum and natural gas industries)', 'IEC 60300-3-3:2017 (life cycle costing)', 'NORSOK O-CR-001 / Z-008 (2017) (maintenance cost and criticality)'],
  models: ['Lifecycle cost analysis', 'Whole-life costing', 'Total cost of ownership', 'Discounted maintenance cost', 'Fixed and variable operating cost build-up', 'Energy-demand cost (fuel and imported power)', 'Equivalent annual cost and capital recovery factor'],
  limits: `Annual cost model in constant (real) base-year money with one real escalation rate on operating cost; costs fall at year ends, capital at year 0 and decommissioning in the year after the last operating year. Failure-driven cost is entered as an annual expectation (see the reliability module); tax and financing are excluded. ${CUR_NOTE}`,
  publishes: ['econ.opexAnnual', 'econ.opexPerMW', 'econ.lifecycleCost', 'econ.lifecycleCostPV', 'econ.eac'],
  searchTerms: 'offshore operating expenditure life cycle cost whole life costing ISO 15686-5',
  inputs: [
    { key: 'baseYear', label: 'Price base year', unit: 'yr', type: 'int', default: 2026, min: 1990, max: 2100, group: 'Opening financial state', role: 'ic', help: 'All costs below are in constant money of this year (real terms).' },
    cost({ key: 'capex', label: 'Capital expenditure (acquisition)', unit: CUR, default: 1.5e9, min: 0, max: 1e12, group: 'Opening financial state', bind: 'econ.capex', help: 'Initial investment committed at time zero; from the capital cost estimate when available.' }),
    { key: 'ratedMW', label: 'Installed capacity', unit: 'MW', type: 'number', default: 500, min: 0, max: 20000, group: 'Project information', role: 'operation', bind: 'energy.ratedPower' },
    { key: 'production', label: 'Annual production', unit: 'GWh/yr', type: 'number', default: 2150, min: 0, max: 2e5, group: 'Project information', role: 'operation', bind: 'energy.aep', help: 'For variable cost and cost per unit of output; 0 if not applicable.' },
    { key: 'opLife', label: 'Operating life', unit: 'yr', type: 'int', default: 25, min: 1, max: 60, group: 'Project information', role: 'bc', bind: 'life.designYears', design: true, help: 'Operational constraint from the design or fatigue life supplied by the engineering suites.' },
    { key: 'discountRate', label: 'Real discount rate', unit: '%/yr', type: 'number', default: 6, min: -5, max: 40, group: 'Financial constraints', role: 'bc', help: `Real rate matching the constant-money costs. ${PLACEHOLDER}` },
    { key: 'opexCap', label: 'Annual operating budget', unit: `${CUR}/yr`, type: 'number', default: 6.5e7, min: 1, max: 1e11, group: 'Financial constraints', role: 'bc', help: 'Affordability limit for first-year OPEX.' },

    { key: 'personnel', label: 'Personnel', type: 'table', group: 'Operational expenditure', role: 'cost', help: `Head count and fully-loaded annual cost per head. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Role', type: 'text' }, { key: 'count', label: 'Head count', type: 'number' }, { key: 'costPerYear', label: 'Cost per head', unit: `${CUR}/yr`, type: 'number' }], default: PERSONNEL },
    { key: 'logistics', label: 'Logistics and offshore support', type: 'table', group: 'Operational expenditure', role: 'cost', help: `Vessel and aircraft days per year at their day rates. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Service', type: 'text' }, { key: 'days', label: 'Days per year', type: 'number' }, { key: 'dayRate', label: 'Day rate', unit: `${CUR}/day`, type: 'number' }], default: LOGISTICS },
    cost({ key: 'maintPlanned', label: 'Planned maintenance (parts and contracts)', unit: `${CUR}/yr`, default: 6.5e6, min: 0, max: 1e11, group: 'Maintenance', uncertain: { dist: 'normal', cov: 0.15 } }),
    cost({ key: 'maintCorrective', label: 'Corrective maintenance (expected)', unit: `${CUR}/yr`, default: 5.5e6, min: 0, max: 1e11, group: 'Maintenance', bind: 'econ.correctiveAnnual', uncertain: { dist: 'lognormal', cov: 0.35 }, help: 'Expected annual repair cost; from the reliability and maintenance module when available.' }),
    cost({ key: 'insuranceRate', label: 'Insurance premium', unit: '% of capex/yr', default: 0.45, min: 0, max: 10, group: 'Operational expenditure' }),
    cost({ key: 'leases', label: 'Seabed lease, licences and levies', unit: `${CUR}/yr`, default: 3.0e6, min: 0, max: 1e11, group: 'Operational expenditure' }),
    cost({ key: 'gridCharges', label: 'Grid, transmission or tariff charges', unit: `${CUR}/yr`, default: 2.5e6, min: 0, max: 1e11, group: 'Operational expenditure' }),
    cost({ key: 'consumables', label: 'Chemicals, utilities and consumables', unit: `${CUR}/yr`, default: 1.2e6, min: 0, max: 1e11, group: 'Operational expenditure' }),
    cost({ key: 'variablePerMWh', label: 'Variable cost per unit of output', unit: `${CUR}/MWh`, default: 1.0, min: 0, max: 1e4, group: 'Operational expenditure' }),
    { key: 'inspections', label: 'Inspection campaigns', type: 'table', group: 'Maintenance', role: 'cost', help: `Campaign cost and interval; annualised as cost / interval. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Campaign', type: 'text' }, { key: 'cost', label: 'Cost per campaign', unit: CUR, type: 'number' }, { key: 'intervalYears', label: 'Interval', unit: 'yr', type: 'number' }], default: INSPECTIONS },

    { key: 'auxPower', label: 'Average imported electrical power', unit: 'kW', type: 'number', default: 900, min: 0, max: 1e6, group: 'Energy demand', role: 'load', help: 'Auxiliary systems, pumps, compressors and heating drawn from the grid or generators.' },
    cost({ key: 'powerPrice', label: 'Price of imported power', unit: `${CUR}/MWh`, default: 110, min: 0, max: 5000, group: 'Energy demand' }),
    { key: 'fuel', label: 'Fuel consumption (vessels and generators)', unit: 't/yr', type: 'number', default: 3600, min: 0, max: 1e6, group: 'Energy demand', role: 'load', design: true },
    cost({ key: 'fuelPrice', label: 'Marine gas oil price', unit: `${CUR}/t`, default: 820, min: 0, max: 1e4, group: 'Energy demand', uncertain: { dist: 'lognormal', cov: 0.3 } }),
    cost({ key: 'carbonPrice', label: 'Carbon price on operational emissions', unit: `${CUR}/t CO₂`, default: 90, min: 0, max: 2000, group: 'Energy demand', help: 'Emission allowance or internal carbon price applied to the fuel burnt.' }),

    cost({ key: 'opexEscalation', label: 'Real escalation of operating cost', unit: '%/yr', default: 1.0, min: -5, max: 15, group: 'Lifecycle', help: 'Growth above general inflation, e.g. ageing-related maintenance and wage drift.' }),
    { key: 'renewals', label: 'Replacement and refurbishment', type: 'table', group: 'Lifecycle', role: 'cost', help: `Major renewals by operating year. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Renewal', type: 'text' }, { key: 'year', label: 'Operating year', type: 'number' }, { key: 'cost', label: 'Cost', unit: CUR, type: 'number' }], default: RENEWALS },
    cost({ key: 'decom', label: 'Decommissioning cost', unit: CUR, default: 1.2e8, min: 0, max: 1e12, group: 'Lifecycle', bind: 'econ.decomCost', help: 'Terminal condition, spent in the year after the last operating year.' }),
    cost({ key: 'residual', label: 'Residual value at end of life', unit: CUR, default: 1.0e7, min: 0, max: 1e12, group: 'Lifecycle', help: 'Scrap or resale value credited with the decommissioning.' }),
  ],
  compute(inp) {
    const pers = tableRows(inp.personnel, 'Personnel', { count: { min: 0 }, costPerYear: { min: 0 } }, { min: 0 });
    const logi = tableRows(inp.logistics, 'Logistics', { days: { min: 0, max: 366 * 50 }, dayRate: { min: 0 } }, { min: 0 });
    const insp = tableRows(inp.inspections, 'Inspection campaigns', { cost: { min: 0 }, intervalYears: { min: 0.1 } }, { min: 0 });
    const ren = tableRows(inp.renewals, 'Replacement and refurbishment', { year: { min: 1 }, cost: { min: 0 } }, { min: 0 });
    uniqueNames(pers, 'Personnel'); uniqueNames(logi, 'Logistics'); uniqueNames(insp, 'Inspection campaigns'); uniqueNames(ren, 'Replacement and refurbishment');
    const n = Math.round(inp.opLife); const r = inp.discountRate / 100; const g = 1 + inp.opexEscalation / 100;
    need(r > -1, 'The discount rate must be above −100 %.');
    for (const x of ren) need(x.year <= n, `Renewal '${x.name}' is scheduled in operating year ${x.year}, after the ${n}-year operating life.`);
    const mwh = inp.production * 1e3;

    const energyCost = inp.auxPower / 1e3 * 8760 * inp.powerPrice + inp.fuel * inp.fuelPrice;
    const co2 = inp.fuel * CO2_MGO;
    const carbon = co2 * inp.carbonPrice;
    const bu = run(() => opexBuildUp({ personnel: pers.map((p) => ({ role: p.name, count: p.count, costPerYear: p.costPerYear })), logistics: logi.map((l) => ({ name: l.name, days: l.days, dayRate: l.dayRate })),
      maintenancePlanned: inp.maintPlanned, maintenanceUnplanned: inp.maintCorrective, insuranceRate: inp.insuranceRate / 100, capex: inp.capex, leases: inp.leases, gridCharges: inp.gridCharges,
      consumables: inp.consumables + energyCost + carbon, inspections: insp.map((i) => ({ name: i.name, cost: i.cost, intervalYears: i.intervalYears })), variablePerUnit: inp.variablePerMWh, production: mwh }), 'OPEX build-up');
    const opex1 = bu.total;
    need(opex1 > 0 || inp.capex > 0, 'Nothing to cost: capital and operating costs are all zero.');
    const persCost = pers.reduce((s, p) => s + p.count * p.costPerYear, 0); const logiCost = logi.reduce((s, l) => s + l.days * l.dayRate, 0); const inspCost = insp.reduce((s, i) => s + i.cost / i.intervalYears, 0);
    const lines = [['Personnel', persCost], ['Logistics and offshore support', logiCost], ['Planned maintenance', inp.maintPlanned], ['Corrective maintenance', inp.maintCorrective], ['Inspection (annualised)', inspCost],
      ['Insurance', inp.insuranceRate / 100 * inp.capex], ['Leases and levies', inp.leases], ['Grid and tariff charges', inp.gridCharges], ['Chemicals and consumables', inp.consumables], ['Energy (power and fuel)', energyCost],
      ['Carbon cost of fuel', carbon], ['Variable cost on output', inp.variablePerMWh * mwh]];
    const chk = lines.reduce((s, l) => s + l[1], 0);
    need(Math.abs(chk - opex1) <= 1e-6 * Math.max(1, opex1), 'Internal reconciliation failed: the OPEX lines do not add up to the build-up total.');

    // life-cycle phases, index 0 = time zero; operating years 1..n; decommissioning in year n + 1
    const N = n + 2;
    const cap = new Array(N).fill(0); cap[0] = inp.capex;
    const op = range(N, (t) => (t >= 1 && t <= n ? (opex1 - inp.maintPlanned - inp.maintCorrective - inspCost) * g ** (t - 1) : 0));
    const mp = range(N, (t) => (t >= 1 && t <= n ? (inp.maintPlanned + inspCost) * g ** (t - 1) : 0));
    const mc = range(N, (t) => (t >= 1 && t <= n ? inp.maintCorrective * g ** (t - 1) : 0));
    const rn = new Array(N).fill(0); for (const x of ren) rn[Math.round(x.year)] += x.cost;
    const dc = new Array(N).fill(0); dc[n + 1] = inp.decom - inp.residual;
    const lcc = run(() => lifeCycleCost({ phases: [{ name: 'Acquisition (CAPEX)', costs: cap }, { name: 'Operations', costs: op }, { name: 'Planned maintenance and inspection', costs: mp }, { name: 'Corrective maintenance', costs: mc }, { name: 'Replacement and refurbishment', costs: rn }, { name: 'Decommissioning less residual value', costs: dc }], rate: r }), 'Life-cycle cost');
    const tco = run(() => totalCostOfOwnership({ acquisition: inp.capex, operating: opex1 - inp.maintPlanned - inp.maintCorrective - inspCost, maintenance: inp.maintPlanned + inspCost + inp.maintCorrective, disposal: inp.decom, residual: inp.residual, life: n, rate: r, escalation: inp.opexEscalation / 100 }), 'Total cost of ownership');
    const eac = lcc.discounted * crf(r, n);
    const pvE = mwh > 0 ? range(n, (k) => mwh / (1 + r) ** (k + 1)).reduce((a, b) => a + b, 0) : 0;
    const phase = (name) => lcc.byPhase.find((p) => p.name === name);
    const maintLife = phase('Planned maintenance and inspection').undiscounted; const corrLife = phase('Corrective maintenance').undiscounted;
    const yearly = range(N, (t) => cap[t] + op[t] + mp[t] + mc[t] + rn[t] + dc[t]);
    let cum = 0; const cumPV = yearly.map((v, t) => { cum += v / (1 + r) ** t; return cum; });

    const kpis = [
      { key: 'opex', label: 'Annual OPEX, first year', value: mil(opex1), unit: `${MU}/yr`, limit: mil(inp.opexCap), sense: 'max', explain: `All operating cost lines in ${inp.baseYear} money for the first operating year. It rises afterwards at the real escalation rate.` },
      { key: 'fixed', label: 'Fixed OPEX', value: mil(bu.fixed), unit: `${MU}/yr`, explain: 'Costs incurred whether or not the asset produces: people, vessels, insurance, leases, maintenance contracts and inspection.' },
      { key: 'variable', label: 'Variable OPEX', value: mil(bu.variable), unit: `${MU}/yr`, explain: 'Costs proportional to output. A low share means lost production is not compensated by saved cost.' },
      { key: 'opexPerMW', label: 'OPEX per MW installed', value: inp.ratedMW > 0 ? round(opex1 / inp.ratedMW / 1e3) : 0, unit: `k${CUR}/MW/yr`, explain: 'Benchmark denominator for comparison with operating assets of the same type and distance to port.' },
      { key: 'opexPerMWh', label: 'OPEX per unit of output', value: mwh > 0 ? round(opex1 / mwh) : 0, unit: `${CUR}/MWh`, explain: 'First-year operating cost divided by annual production: the operating part of the levelised cost.' },
      { key: 'opexRatio', label: 'OPEX as share of CAPEX', value: inp.capex > 0 ? round(100 * opex1 / inp.capex) : 0, unit: '%/yr', explain: 'Annual operating cost relative to the investment; a quick plausibility check against comparable assets.' },
      { key: 'lifeOpex', label: 'Lifecycle OPEX (undiscounted)', value: mil(phase('Operations').undiscounted + maintLife + corrLife), unit: MU, explain: `Sum of operating and maintenance cost over ${n} years with real escalation, before renewals and decommissioning.` },
      { key: 'maintPlannedLife', label: 'Planned maintenance over the life', value: mil(maintLife), unit: MU, explain: 'Planned maintenance and annualised inspection summed over the operating life (undiscounted).' },
      { key: 'maintCorrLife', label: 'Corrective maintenance over the life', value: mil(corrLife), unit: MU, explain: 'Expected failure-driven repair cost summed over the operating life (undiscounted). It is an expectation: individual years will differ widely.' },
      { key: 'lcc', label: 'Total lifecycle cost (undiscounted)', value: mil(lcc.total), unit: MU, explain: 'Every cost from acquisition to decommissioning, less residual value, added without discounting: the cash the asset consumes over its life.' },
      { key: 'lccPV', label: 'Discounted lifecycle cost', value: mil(lcc.discounted), unit: MU, explain: `Present value at ${inp.baseYear} of the same cost stream at the real discount rate: the figure for comparing design options with different cost timing.` },
      { key: 'eac', label: 'Equivalent annual cost', value: mil(eac), unit: `${MU}/yr`, explain: 'Discounted lifecycle cost spread as a level annual amount over the operating life (× capital recovery factor); compare options with different lives on this basis.' },
      { key: 'tco', label: 'Total cost of ownership (formula check)', value: mil(tco.tco), unit: MU, explain: 'Acquisition plus present value of operating, maintenance and disposal cost less residual value from the closed-form annuity expressions, excluding scheduled renewals. It reconciles with the year-by-year figure.' },
      { key: 'levelised', label: 'Levelised cost of output (pre-tax, real)', value: pvE > 0 ? round(lcc.discounted / pvE) : 0, unit: `${CUR}/MWh`, explain: 'Discounted lifecycle cost over discounted production: the constant real price that just recovers all costs at the discount rate.' },
      { key: 'co2', label: 'Operational CO₂ from fuel', value: round(co2), unit: 't CO₂/yr', explain: `Fuel burnt by vessels and generators at ${CO2_MGO} t CO₂ per tonne of marine gas oil (screening factor).` },
    ];
    const checks = [check('First-year OPEX within operating budget', opex1, inp.opexCap, `${CUR}/yr`, 'Financial constraint (operating budget)', 'Guards against an operating concept whose running cost exceeds what the owner has budgeted; above 1.0 the strategy or the budget must change.')];
    const yrs = range(N, (t) => inp.baseYear + t);
    const plots = [
      { type: 'bar', title: 'Annual OPEX breakdown (first year)', categories: lines.map((l) => l[0]), x: { label: 'Cost line' }, y: { label: 'Cost', unit: `${MU}/yr` }, series: [{ name: 'OPEX', y: lines.map((l) => mil(l[1])) }],
        explain: 'Each operating cost line. The tallest bars are where the maintenance and logistics strategy has most leverage.' },
      { type: 'bar', title: 'Lifecycle cost profile by year', categories: yrs.map(String), x: { label: 'Year' }, y: { label: 'Cost', unit: MU },
        series: [{ name: 'Acquisition', y: cap.map(mil) }, { name: 'Operations', y: op.map(mil) }, { name: 'Maintenance and inspection', y: mp.map((v, t) => mil(v + mc[t])) }, { name: 'Renewals', y: rn.map(mil) }, { name: 'Decommissioning net of residual', y: dc.map(mil) }],
        explain: 'When each cost falls. Look for the renewal spikes in mid-life and the decommissioning bar at the end: both are easy to under-provide for.' },
      { type: 'line', title: 'Cumulative discounted lifecycle cost', x: { label: 'Year', unit: '' }, y: { label: 'Present value', unit: MU }, series: [{ name: 'Cumulative PV', x: yrs, y: cumPV.map(mil) }, { name: 'Acquisition only', x: [yrs[0], yrs[N - 1]], y: [mil(inp.capex), mil(inp.capex)], dash: true }],
        explain: 'Growth of the present value of cost over the life. The gap above the acquisition line is the discounted cost of owning and operating the asset.' },
    ];
    const tables = [
      { title: 'OPEX build-up', columns: [{ key: 'line', label: 'Line' }, { key: 'cost', label: 'Annual cost', unit: `${MU}/yr` }, { key: 'share', label: 'Share', unit: '%' }],
        rows: [...lines.map(([line, c]) => ({ line, cost: mil(c), share: round(100 * c / (opex1 || 1), 3) })), { line: 'Total', cost: mil(opex1), share: 100 }], explain: 'Every operating cost line and its share of the total; reconcile against recorded expenditure of comparable operating assets.' },
      { title: 'Lifecycle cost by phase', columns: [{ key: 'name', label: 'Phase' }, { key: 'u', label: 'Undiscounted', unit: MU }, { key: 'd', label: 'Discounted', unit: MU }, { key: 's', label: 'Share of discounted', unit: '%' }],
        rows: [...lcc.byPhase.map((p) => ({ name: p.name, u: mil(p.undiscounted), d: mil(p.discounted), s: round(100 * p.share, 3) })), { name: 'Total', u: mil(lcc.total), d: mil(lcc.discounted), s: 100 }],
        explain: 'Undiscounted and discounted cost of each life-cycle phase. Discounting shrinks late costs: the discounted share of decommissioning is far smaller than its cash amount.' },
    ];
    const top = lines.slice().sort((a, b) => b[1] - a[1])[0];
    const opShare = 1 - phase('Acquisition (CAPEX)').share;
    const vesselShare = logiCost / (opex1 || 1);
    const advice = [
      adv('engineering', 'info', `"${top[0]}" is the largest operating line at ${fM(top[1])}/yr (${pct(top[1] / (opex1 || 1))}); logistics take ${pct(vesselShare)} and maintenance plus inspection ${pct((inp.maintPlanned + inp.maintCorrective + inspCost) / (opex1 || 1))}.`,
        'Reduce the driver behind it: access strategy and distance to port set vessel days, component reliability sets corrective cost, and inspection intervals follow from fatigue and corrosion margins — a design with more margin can justify longer intervals.',
        'Risk-based maintenance and inspection planning, NORSOK Z-008 (2017) and DNV-RP-G101 (2021).',
        `Operations burn ${fmt(inp.fuel)} t of fuel a year, ${fmt(co2)} t CO₂/yr and ${fmt(co2 * n)} t over the life (${SCREEN}); at ${inp.carbonPrice} ${CUR}/t this is ${fM(carbon)}/yr. Cutting vessel days by 10 % saves about ${fmt(0.1 * co2)} t CO₂ and ${fM(0.1 * (logiCost + inp.fuel * inp.fuelPrice + carbon))} a year.`),
      adv('economic', opex1 <= inp.opexCap ? 'ok' : 'fail', `First-year OPEX ${fM(opex1)} (${inp.capex > 0 ? pct(opex1 / inp.capex) : '–'} of CAPEX). Discounted lifecycle cost ${fM(lcc.discounted)}, of which ${pct(opShare)} arises after acquisition; equivalent annual cost ${fM(eac)}/yr.`,
        opex1 <= inp.opexCap ? `Compare design options on discounted lifecycle cost, not CAPEX: at ${inp.discountRate} % real, one unit of annual OPEX saved for ${n} years is worth ${fmt(1 / crf(r, n))} units of CAPEX today.` : `OPEX exceeds the ${fM(inp.opexCap)}/yr budget by ${fM(opex1 - inp.opexCap)}/yr: revisit the logistics concept or the budget.`,
        'Whole-life costing per ISO 15686-5:2017 §5 and ISO 15663:2021; costs in real terms discounted at a real rate.'),
    ];
    return { kpis, checks, plots, tables, advice,
      notes: [`Real terms: constant ${inp.baseYear} money discounted at a real rate of ${inp.discountRate} %/yr; operating cost escalates ${inp.opexEscalation} %/yr above inflation.`, 'Timing: acquisition in year 0, operating cost at the end of operating years 1 to n, decommissioning net of residual value in year n + 1.', 'All rates, salaries and day rates are indicative placeholders.'],
      publish: { 'econ.opexAnnual': opex1, 'econ.opexPerMW': inp.ratedMW > 0 ? opex1 / inp.ratedMW : 0, 'econ.lifecycleCost': lcc.total, 'econ.lifecycleCostPV': lcc.discounted, 'econ.eac': eac },
      quantities: { vessel_days: { operations_per_year: round(logi.reduce((s, l) => s + l.days, 0)) } } };
  },
  verify: [
    // OPEX = 10 × 100 000 = 1.0 M/yr; LCC_PV = 10 + 1.0 × (P/A, 10 %, 10) = 10 + 6.144567 = 16.144567 M; EAC = 16.144567 × (A/P, 10 %, 10) = 16.144567 × 0.162745 = 2.627454 M/yr
    { name: 'Capital plus level annual cost: 10 M + 1 M/yr for 10 years at 10 %', source: 'Hand calculation with (P/A, 10 %, 10) = 6.144567 and (A/P, 10 %, 10) = 0.162745 (Blank & Tarquin tables)',
      inputs: { capex: 1e7, opLife: 10, discountRate: 10, personnel: [{ name: 'crew', count: 10, costPerYear: 1e5 }], logistics: [], inspections: [], renewals: [], maintPlanned: 0, maintCorrective: 0, insuranceRate: 0, leases: 0, gridCharges: 0, consumables: 0, variablePerMWh: 0, auxPower: 0, fuel: 0, opexEscalation: 0, decom: 0, residual: 0 },
      expect: [{ kpi: 'opex', value: 1.0, tol: 1e-9 }, { kpi: 'lccPV', value: 16.144567, tol: 1e-6 }, { kpi: 'eac', value: 2.627454, tol: 1e-5 }, { kpi: 'tco', value: 16.144567, tol: 1e-6 }, { kpi: 'lcc', value: 20, tol: 1e-9 }] },
    // 1000 t/yr fuel → 3206 t CO2; energy cost = 1000 × 800 + 1 MW × 8760 h × 100 = 0.8 M + 0.876 M; carbon 3206 × 50 = 0.1603 M → OPEX 1.8363 M
    { name: 'Energy demand cost and emissions', source: 'Hand calculation: 1000 t × 800 + 8760 MWh × 100 + 3206 t CO₂ × 50 = 1.8363 M; IMO factor 3.206 t CO₂/t MGO',
      inputs: { personnel: [], logistics: [], inspections: [], renewals: [], maintPlanned: 0, maintCorrective: 0, insuranceRate: 0, leases: 0, gridCharges: 0, consumables: 0, variablePerMWh: 0, auxPower: 1000, powerPrice: 100, fuel: 1000, fuelPrice: 800, carbonPrice: 50 },
      expect: [{ kpi: 'opex', value: 1.8363, tol: 1e-6 }, { kpi: 'co2', value: 3206, tol: 1e-9 }] },
  ],
};
