// Bottom-up capital cost estimate: engineering quantities × unit rates, equipment quotations, fabrication labour,
// marine installation spread × weather-adjusted days, indirects, contingency and escalation to the spend date.
import { bottomUpEstimate, installationCampaignCost, estimateRange, benchmarkMetrics, AACE_CLASSES } from '../../engine/econ/cost.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, CO2_STEEL, CO2_CONC, CO2_MGO, MGO_T_PER_DAY, SCREEN, need, run, round, fmt, mil, fM, pct, adv, check, cost, tableRows, uniqueNames } from './_common.js';

const EQUIP = [
  { name: 'Wind turbine generators (supply)', discipline: 'Generation', qty: 34, unitCost: 1.7e7, fx: 1 },
  { name: 'Offshore substation topside', discipline: 'Electrical', qty: 1, unitCost: 9.5e7, fx: 1 },
  { name: 'Onshore substation and grid connection', discipline: 'Electrical', qty: 1, unitCost: 6.0e7, fx: 1 },
  { name: 'SCADA, control and communications', discipline: 'Control', qty: 1, unitCost: 9.0e6, fx: 1 },
  { name: 'Platform crane and auxiliary systems', discipline: 'Mechanical', qty: 1, unitCost: 6.5e6, fx: 1 },
];
const FAB = [
  { name: 'Secondary steel and outfitting labour', hours: 240000, rate: 68, lump: 0 },
  { name: 'Welding and NDT inspection', hours: 60000, rate: 85, lump: 0 },
  { name: 'Coating and corrosion protection', hours: 0, rate: 0, lump: 1.4e7 },
  { name: 'Quality assurance and yard overheads', hours: 0, rate: 0, lump: 9.0e6 },
  { name: 'Load-out and sea transport', hours: 0, rate: 0, lump: 1.8e7 },
];
const SPREAD = [
  { vessel: 'Heavy-lift / installation jack-up', count: 1, dayRate: 280000, share: 1 },
  { vessel: 'Cable-lay vessel', count: 1, dayRate: 150000, share: 0.45 },
  { vessel: 'Offshore construction / ROV support vessel', count: 1, dayRate: 110000, share: 0.6 },
  { vessel: 'Tugs and feeder barges', count: 3, dayRate: 22000, share: 1 },
  { vessel: 'Crew transfer and guard vessels', count: 2, dayRate: 6000, share: 1 },
];

export default {
  id: 's17-capex-estimate', suite: 17, title: 'Capital cost estimate (bottom-up, quantity based)',
  summary: 'Converts the material quantities, equipment list, fabrication effort and marine installation days delivered by the engineering suites into a classified CAPEX estimate with its breakdown and accuracy range.',
  covers: [2], fidelity: 'analytical',
  standards: ['AACE International RP 18R-97 (cost estimate classification, process industries)', 'AACE International RP 17R-97 (cost estimate classification system)', 'AACE International RP 34R-05 (basis of estimate)', 'AACE International RP 58R-10 (escalation)', 'DNV-ST-N001 (2021) §3 (weather-restricted marine operations, for the operability allowance)'],
  models: ['Bottom-up engineering cost estimation', 'Quantity-based unit-rate cost model', 'Vessel day-rate installation cost model', 'Weather-operability allowance on marine durations', 'Indirect-cost and contingency allowances', 'Procurement cost escalation', 'Currency conversion of quotations', 'AACE estimate classes and accuracy ranges'],
  limits: `Deterministic point estimate: quantities come from the bus or the form, every rate is a user input, and indirects are percentages of direct cost. Market conditions, contract strategy and schedule-driven cost are not modelled here (see the cost-risk and schedule-risk modules). ${CUR_NOTE}`,
  publishes: ['econ.capex', 'econ.capexDirect', 'econ.installCost', 'econ.materialCost', 'econ.capexPerMW', 'econ.embodiedCO2'],
  searchTerms: 'offshore capital cost estimate unit rates vessel day rates AACE estimate class',
  inputs: [
    { key: 'projectType', label: 'Project type', type: 'select', default: 'wind-fixed', group: 'Project information', role: 'operation',
      options: [{ value: 'wind-fixed', label: 'Offshore wind, fixed foundations' }, { value: 'wind-floating', label: 'Offshore wind, floating' }, { value: 'platform', label: 'Oil and gas platform' }, { value: 'subsea', label: 'Subsea tie-back / pipeline' }, { value: 'terminal', label: 'Port or marine terminal' }, { value: 'other', label: 'Other marine asset' }],
      help: 'Label for the report and the benchmark denominators; it does not change any rate.' },
    { key: 'ratedMW', label: 'Design capacity', unit: 'MW', type: 'number', default: 510, min: 0, max: 20000, group: 'Project information', role: 'operation', bind: 'energy.ratedPower', help: 'For the cost-per-MW benchmark. Enter 0 for assets without a power rating.' },
    { key: 'baseYear', label: 'Price base year of the rates', unit: 'yr', type: 'int', default: 2026, min: 1990, max: 2100, group: 'Opening financial state', role: 'ic', help: 'All rates and quotations below are in money of this year.' },
    { key: 'spendYear', label: 'Mid-point of expenditure', unit: 'yr', type: 'int', default: 2028, min: 1990, max: 2110, group: 'Opening financial state', role: 'ic', help: 'The estimate is escalated from the base year to this year.' },
    { key: 'aaceClass', label: 'Estimate class (AACE 18R-97)', type: 'select', default: '3', group: 'Opening financial state', role: 'ic',
      options: [5, 4, 3, 2, 1].map((c) => ({ value: String(c), label: `Class ${c} — ${AACE_CLASSES[c].usage}` })), help: 'Maturity of project definition; sets the expected accuracy range around the point estimate.' },
    { key: 'budgetCap', label: 'Available capital budget', unit: CUR, type: 'number', default: 2.75e9, min: 1, max: 1e12, group: 'Financial constraints', role: 'bc', help: 'Funding available including headroom (committed facilities plus contingent equity). The escalated estimate and the high end of its accuracy range are both checked against it.' },

    { key: 'steel', label: 'Structural steel', unit: 'tonnes', type: 'number', default: 62000, min: 0, max: 5e6, group: 'Material quantities', role: 'geometry', bind: 'quantities.steel_t', design: true, uncertain: { dist: 'normal', cov: 0.07 }, help: 'Primary and secondary steel from the structural suites.' },
    { key: 'concrete', label: 'Concrete', unit: 'm³', type: 'number', default: 4500, min: 0, max: 5e6, group: 'Material quantities', role: 'geometry', bind: 'quantities.concrete_m3' },
    { key: 'rock', label: 'Rock and scour protection', unit: 'tonnes of rock', type: 'number', default: 150000, min: 0, max: 5e7, group: 'Material quantities', role: 'geometry', bind: 'quantities.rock_t' },
    { key: 'cable', label: 'Subsea cable route length', unit: 'km of cable', type: 'number', default: 165, min: 0, max: 5000, group: 'Material quantities', role: 'geometry', bind: 'quantities.cable_km' },
    { key: 'pipe', label: 'Pipeline length', unit: 'km of pipe', type: 'number', default: 0, min: 0, max: 5000, group: 'Material quantities', role: 'geometry', bind: 'quantities.pipe_km' },
    { key: 'piles', label: 'Pile length driven or drilled', unit: 'm of pile', type: 'number', default: 1300, min: 0, max: 1e6, group: 'Material quantities', role: 'geometry' },
    { key: 'coating', label: 'Coated surface area', unit: 'm² coated', type: 'number', default: 210000, min: 0, max: 1e8, group: 'Material quantities', role: 'geometry' },

    cost({ key: 'rSteel', label: 'Fabricated steel', unit: `${CUR}/t`, default: 3200, min: 0, max: 1e5, group: 'Material unit costs', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Supplied and fabricated at yard.' }),
    cost({ key: 'rConcrete', label: 'Concrete in place', unit: `${CUR}/m³`, default: 650, min: 0, max: 1e5, group: 'Material unit costs' }),
    cost({ key: 'rRock', label: 'Rock supplied and placed', unit: `${CUR}/t`, default: 55, min: 0, max: 1e4, group: 'Material unit costs' }),
    cost({ key: 'rCable', label: 'Cable supply', unit: `${CUR}/km`, default: 7.0e5, min: 0, max: 1e8, group: 'Material unit costs', uncertain: { dist: 'lognormal', cov: 0.2 } }),
    cost({ key: 'rPipe', label: 'Linepipe supplied, coated and welded', unit: `${CUR}/km`, default: 1.1e6, min: 0, max: 1e8, group: 'Material unit costs' }),
    cost({ key: 'rPile', label: 'Pile installation', unit: `${CUR}/m`, default: 9000, min: 0, max: 1e6, group: 'Material unit costs' }),
    cost({ key: 'rCoating', label: 'Coating', unit: `${CUR}/m²`, default: 45, min: 0, max: 1e4, group: 'Material unit costs' }),

    { key: 'equipment', label: 'Equipment and packages', type: 'table', group: 'Equipment costs', role: 'cost',
      help: `Supplier quotations. FX = project currency per unit of the quotation currency (1 when quoted in ${CUR}). ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Item', type: 'text' }, { key: 'discipline', label: 'Discipline', type: 'text' }, { key: 'qty', label: 'Quantity', type: 'number' }, { key: 'unitCost', label: 'Unit cost (quote currency)', type: 'number' }, { key: 'fx', label: 'FX to project currency', type: 'number' }],
      default: EQUIP },
    { key: 'fabrication', label: 'Fabrication, QA and transport', type: 'table', group: 'Fabrication costs', role: 'cost',
      help: `Labour as hours × all-in rate, plus lump sums for subcontracts and overheads. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Activity', type: 'text' }, { key: 'hours', label: 'Labour hours', unit: 'h', type: 'number' }, { key: 'rate', label: 'All-in rate', unit: `${CUR}/h`, type: 'number' }, { key: 'lump', label: 'Lump sum', unit: CUR, type: 'number' }],
      default: FAB },

    { key: 'vesselDays', label: 'Net offshore working days of the campaign', unit: 'vessel-days', type: 'number', default: 210, min: 0, max: 5000, group: 'Installation schedule', role: 'operation', bind: 'install.vesselDays', design: true, uncertain: { dist: 'lognormal', cov: 0.2 },
      help: 'Productive days of the main spread, before weather. From the installation analyses when available.' },
    { key: 'operability', label: 'Weather operability', unit: '-', type: 'number', default: 0.68, min: 0.05, max: 1, group: 'Installation schedule', role: 'environment', bind: 'install.operabilityFrac', uncertain: { dist: 'normal', std: 0.07 },
      help: 'Fraction of time the limiting operation is workable; gross days = net days / operability.' },
    { key: 'transitDays', label: 'Transit and standby days', unit: 'vessel-days', type: 'number', default: 14, min: 0, max: 1000, group: 'Installation schedule', role: 'operation' },
    { key: 'spread', label: 'Marine spread', type: 'table', group: 'Marine installation costs', role: 'cost',
      help: `Share = fraction of the campaign each vessel is on hire. ${PLACEHOLDER}`,
      columns: [{ key: 'vessel', label: 'Vessel', type: 'text' }, { key: 'count', label: 'Number', type: 'number' }, { key: 'dayRate', label: 'Day rate', unit: `${CUR}/day`, type: 'number' }, { key: 'share', label: 'Share of campaign', unit: '-', type: 'number' }],
      default: SPREAD },
    cost({ key: 'mobDemob', label: 'Mobilisation and demobilisation', unit: CUR, default: 1.2e7, min: 0, max: 1e10, group: 'Marine installation costs' }),
    cost({ key: 'rCableLay', label: 'Cable or pipe lay and burial (non-vessel)', unit: `${CUR}/km`, default: 1.2e5, min: 0, max: 1e8, group: 'Marine installation costs', help: 'Trenching tools, survey, pull-ins and terminations per kilometre.' }),

    cost({ key: 'pEngineering', label: 'Engineering and design', unit: '% of direct', default: 6, min: 0, max: 50, group: 'Indirects and allowances' }),
    cost({ key: 'pManagement', label: 'Project management and supervision', unit: '% of direct', default: 4, min: 0, max: 50, group: 'Indirects and allowances' }),
    cost({ key: 'pInsurance', label: 'Construction insurance', unit: '% of direct', default: 1.5, min: 0, max: 20, group: 'Indirects and allowances' }),
    cost({ key: 'pOwner', label: "Owner's costs, permits and surveys", unit: '% of direct', default: 4, min: 0, max: 50, group: 'Indirects and allowances' }),
    cost({ key: 'pCommissioning', label: 'Commissioning and start-up', unit: '% of direct', default: 2, min: 0, max: 30, group: 'Indirects and allowances' }),
    cost({ key: 'pContingency', label: 'Contingency', unit: '% of base estimate', default: 12, min: 0, max: 100, group: 'Indirects and allowances', help: 'Allowance for items expected to be spent but not yet defined; not a reserve for scope change.' }),
    cost({ key: 'escalation', label: 'Cost escalation', unit: '%/yr', default: 3, min: -10, max: 40, group: 'Indirects and allowances', help: 'Composite escalation index growth between the base year and the spend year.' }),
    cost({ key: 'carbonPrice', label: 'Internal carbon price', unit: `${CUR}/t CO₂e`, default: 90, min: 0, max: 2000, group: 'Indirects and allowances', help: 'Used only to express embodied emissions as a shadow cost.' }),
  ],
  compute(inp) {
    const eq = tableRows(inp.equipment, 'Equipment', { qty: { min: 0 }, unitCost: { min: 0 }, fx: { min: 1e-9 } }, { min: 0 });
    const fab = tableRows(inp.fabrication, 'Fabrication', { hours: { min: 0 }, rate: { min: 0 }, lump: { min: 0 } }, { min: 0 });
    const sp = tableRows(inp.spread, 'Marine spread', { count: { min: 0 }, dayRate: { min: 0 }, share: { min: 0, max: 1 } }, { min: 0 });
    uniqueNames(eq, 'Equipment'); uniqueNames(fab, 'Fabrication'); uniqueNames(sp, 'Marine spread', 'vessel');
    need(inp.spendYear >= inp.baseYear, 'The mid-point of expenditure cannot be earlier than the price base year of the rates.');
    need(inp.vesselDays === 0 || sp.length > 0, 'Offshore working days are given but the marine spread table is empty: add the vessels or set the days to zero.');

    const items = [];
    const add = (name, quantity, unitRate, unit, discipline, category) => { if (quantity > 0 && unitRate > 0) items.push({ name, quantity, unitRate, unit, discipline, category }); };
    add('Structural steel, fabricated', inp.steel, inp.rSteel, 't', 'Structural', 'Materials');
    add('Concrete', inp.concrete, inp.rConcrete, 'm³', 'Civil', 'Materials');
    add('Rock and scour protection', inp.rock, inp.rRock, 't', 'Civil', 'Materials');
    add('Subsea cable supply', inp.cable, inp.rCable, 'km', 'Electrical', 'Materials');
    add('Linepipe', inp.pipe, inp.rPipe, 'km', 'Pipeline', 'Materials');
    add('Coating', inp.coating, inp.rCoating, 'm²', 'Structural', 'Materials');
    for (const e of eq) add(e.name, e.qty, e.unitCost * e.fx, 'each', e.discipline || 'Equipment', 'Equipment');
    for (const f of fab) { add(`${f.name} (labour)`, f.hours, f.rate, 'h', 'Fabrication', 'Fabrication'); add(f.lump > 0 && f.hours > 0 ? `${f.name} (lump sum)` : f.name, f.lump > 0 ? 1 : 0, f.lump, 'lump', 'Fabrication', 'Fabrication'); }
    // marine campaign: weather-adjusted gross days per vessel from the engine's day-rate model
    let grossMain = 0; let vesselDayTotal = 0; let fuel = 0;
    for (const s of sp) {
      if (!(s.count > 0 && s.share > 0 && inp.vesselDays > 0)) continue;
      const c = run(() => installationCampaignCost({ spread: [{ vessel: s.vessel, count: s.count, dayRate: s.dayRate }], netDays: inp.vesselDays * s.share, operability: inp.operability, extraDays: inp.transitDays * s.share, mobDemob: 0 }), `Marine spread '${s.vessel}'`);
      grossMain = Math.max(grossMain, c.grossDays); vesselDayTotal += c.grossDays * s.count;
      // screening fuel burn: construction-class vessels at the main-vessel rate, small craft at 3 t/day
      fuel += c.grossDays * s.count * (s.dayRate >= 50000 ? MGO_T_PER_DAY : 3);
      add(s.vessel, c.grossDays * s.count, s.dayRate, 'vessel-day', 'Marine installation', 'Installation');
    }
    add('Mobilisation and demobilisation', inp.mobDemob > 0 && inp.vesselDays > 0 ? 1 : 0, inp.mobDemob, 'lump', 'Marine installation', 'Installation');
    add('Cable / pipe lay, burial and terminations', inp.cable + inp.pipe, inp.rCableLay, 'km', 'Marine installation', 'Installation');
    add('Pile installation', inp.piles, inp.rPile, 'm', 'Marine installation', 'Installation');
    need(items.length > 0, 'Nothing to estimate: every quantity or rate is zero.');

    const years = inp.spendYear - inp.baseYear;
    const est = run(() => bottomUpEstimate({ items,
      indirects: { 'Engineering and design': inp.pEngineering / 100, 'Project management': inp.pManagement / 100, 'Construction insurance': inp.pInsurance / 100, "Owner's costs": inp.pOwner / 100, Commissioning: inp.pCommissioning / 100 },
      contingency: inp.pContingency / 100, escalation: { rate: inp.escalation / 100, years } }), 'Estimate');
    const cat = (k) => est.byCategory[k] || 0;
    const cls = Number(inp.aaceClass);
    const rng = estimateRange(est.total, cls);
    const bm = benchmarkMetrics({ cost: est.total, MW: inp.ratedMW > 0 ? inp.ratedMW : undefined, tonnes: inp.steel > 0 ? inp.steel : undefined, km: inp.cable + inp.pipe > 0 ? inp.cable + inp.pipe : undefined });

    // screening embodied carbon of the estimate's physical content
    const co2 = inp.steel * CO2_STEEL + inp.concrete * CO2_CONC + fuel * CO2_MGO;
    const carbonCost = co2 * inp.carbonPrice;

    const kpis = [
      { key: 'capex', label: `Total CAPEX (money of ${inp.spendYear})`, value: mil(est.total), unit: MU, limit: mil(inp.budgetCap), sense: 'max',
        explain: `Direct cost plus indirects, contingency and escalation from ${inp.baseYear} to ${inp.spendYear}. This is the Class ${cls} point estimate; the true outcome is expected inside the accuracy range below.` },
      { key: 'direct', label: 'Direct cost', value: mil(est.direct), unit: MU, explain: 'Sum of quantity × unit rate for materials, equipment, fabrication and installation in base-year money, before any percentage allowance.' },
      { key: 'materials', label: 'Material procurement', value: mil(cat('Materials')), unit: MU, explain: 'Steel, concrete, rock, cable, linepipe and coating at their unit rates. Scales directly with the quantities from the engineering suites.' },
      { key: 'equipment', label: 'Equipment and packages', value: mil(cat('Equipment')), unit: MU, explain: 'Quoted equipment converted to the project currency. Usually the least design-sensitive and the most market-sensitive block.' },
      { key: 'fabrication', label: 'Fabrication, QA and transport', value: mil(cat('Fabrication')), unit: MU, explain: 'Yard labour hours at the all-in rate plus subcontract lump sums; it grows with weld volume and outfitting complexity rather than with tonnage alone.' },
      { key: 'installation', label: 'Offshore installation and marine vessels', value: mil(cat('Installation')), unit: MU, explain: 'Vessel day rates over weather-adjusted days plus mobilisation, lay and piling. The block most exposed to weather and schedule.' },
      { key: 'indirects', label: 'Indirects and commissioning', value: mil(est.indirectTotal), unit: MU, explain: 'Engineering, management, insurance, owner’s costs and commissioning as percentages of direct cost.' },
      { key: 'contingency', label: 'Contingency', value: mil(est.contingency), unit: MU, explain: 'Allowance on the base estimate for undefined items within scope. Compare it with the cost-risk module’s P50–P80 result (cumulative percentiles: P80 is the cost not exceeded with 80 % probability) rather than a rule of thumb.' },
      { key: 'escalation', label: 'Escalation', value: mil(est.escalation), unit: MU, explain: `Growth of the whole estimate over ${years} year(s) at the composite escalation rate; it is the difference between base-year and spend-year money.` },
      { key: 'rangeLow', label: `Class ${cls} accuracy range, low`, value: mil(rng.wide.low), unit: MU, explain: 'Lower end of the wide AACE 18R-97 accuracy band (about 80 % confidence after contingency). Outcomes below this are unlikely.' },
      { key: 'rangeHigh', label: `Class ${cls} accuracy range, high`, value: mil(rng.wide.high), unit: MU, explain: 'Upper end of the wide AACE 18R-97 accuracy band. Financing and budget headroom should cover this, not just the point estimate.' },
      { key: 'grossDays', label: 'Main-spread gross days offshore', value: round(grossMain), unit: 'days', explain: 'Net working days divided by operability plus transit and standby. The gap to the net days is paid waiting on weather.' },
      { key: 'co2', label: 'Embodied and installation emissions (screening)', value: round(co2), unit: 't CO₂e', explain: `Steel, concrete and marine fuel only, with ${SCREEN}. An order-of-magnitude figure for comparing options, not a product declaration.` },
    ];
    if (bm.perMW != null) kpis.push({ key: 'perMW', label: 'CAPEX per MW', value: round(bm.perMW / 1e6), unit: `${MU}/MW`, explain: 'Benchmark denominator for comparison with completed projects of the same type, vintage and water depth; differences in scope make raw comparisons unreliable.' });

    const checks = [check('CAPEX within available budget', est.total, inp.budgetCap, CUR, 'Financial constraint (affordability limit)',
      'Guards against sanctioning a project whose point estimate already exceeds the funds available; above 1.0 scope, phasing or funding must change.'),
    check('High end of accuracy range within budget', rng.wide.high, inp.budgetCap, CUR, 'AACE RP 18R-97 Table 1 (expected accuracy range)',
      'Guards against a budget that covers the point estimate but not the credible overrun for this class of estimate; exceeding it means funding headroom or a more mature estimate is needed before commitment.')];

    const disc = Object.entries(est.byDiscipline).sort((a, b) => b[1] - a[1]);
    const waterfall = [['Materials', cat('Materials')], ['Equipment', cat('Equipment')], ['Fabrication', cat('Fabrication')], ['Installation', cat('Installation')], ...Object.entries(est.indirects), ['Contingency', est.contingency], ['Escalation', est.escalation]];
    const plots = [
      { type: 'bar', title: 'CAPEX build-up', categories: waterfall.map((w) => w[0]), x: { label: 'Cost block' }, y: { label: 'Cost', unit: MU }, series: [{ name: 'Cost', y: waterfall.map((w) => mil(w[1])) }],
        explain: 'Each block of the estimate from materials to escalation. The tallest bars are where a better quotation or a design change pays most.' },
      { type: 'bar', title: 'Direct cost by engineering discipline', categories: disc.map((d) => d[0]), x: { label: 'Discipline' }, y: { label: 'Direct cost', unit: MU }, series: [{ name: 'Direct cost', y: disc.map((d) => mil(d[1])) }],
        explain: 'Direct cost grouped by the discipline that owns the quantity. Use it to see which engineering suite drives the money.' },
      { type: 'bar', title: `Point estimate and AACE Class ${cls} range`, categories: ['Wide low', 'Narrow low', 'Estimate', 'Narrow high', 'Wide high', 'Budget'], x: { label: '' }, y: { label: 'CAPEX', unit: MU },
        series: [{ name: 'CAPEX', y: [rng.wide.low, rng.narrow.low, est.total, rng.narrow.high, rng.wide.high, inp.budgetCap].map(mil) }],
        explain: 'The estimate with the narrow and wide ends of its class accuracy band, against the budget. The band is asymmetric: overruns are larger than underruns.' },
    ];
    const tables = [
      { title: 'Estimate lines', columns: [{ key: 'name', label: 'Line' }, { key: 'category', label: 'Block' }, { key: 'discipline', label: 'Discipline' }, { key: 'quantity', label: 'Quantity' }, { key: 'unit', label: 'Unit' }, { key: 'unitRate', label: 'Unit rate', unit: CUR }, { key: 'cost', label: 'Cost', unit: MU }, { key: 'share', label: 'Share of direct', unit: '%' }],
        rows: items.map((i) => ({ name: i.name, category: i.category, discipline: i.discipline, quantity: round(i.quantity), unit: i.unit, unitRate: round(i.unitRate), cost: mil(i.quantity * i.unitRate), share: round(100 * i.quantity * i.unitRate / est.direct, 3) })).sort((a, b) => b.cost - a.cost),
        explain: 'Every priced line, largest first, with its quantity, rate and share of direct cost: the basis of estimate to be reconciled against quotations.' },
      { title: 'Estimate summary', columns: [{ key: 'item', label: 'Item' }, { key: 'cost', label: 'Cost', unit: MU }, { key: 'share', label: 'Share of total', unit: '%' }],
        rows: [...waterfall, ['Total CAPEX', est.total]].map(([item, c]) => ({ item, cost: mil(c), share: round(100 * c / est.total, 3) })),
        explain: 'Roll-up of the estimate by block with each block’s share of the total; the percentages show whether the allowances are in a credible proportion to the direct cost.' },
    ];

    const top = items.map((i) => ({ n: i.name, c: i.quantity * i.unitRate })).sort((a, b) => b.c - a.c)[0];
    const loading = est.total / est.direct;
    const perT = inp.rSteel * loading; const perDay = sp.reduce((s, v) => s + v.count * v.dayRate * v.share, 0) * loading;
    const weatherDays = inp.vesselDays > 0 ? inp.vesselDays * (1 / inp.operability - 1) : 0;
    const contShare = est.contingency / est.base;
    const advice = [
      adv('engineering', 'info', `The largest line is "${top.n}" at ${fM(top.c)} (${pct(top.c / est.direct)} of direct cost). With allowances, each tonne of steel carries ${fmt(perT)} ${CUR} and each net working day offshore ${fM(perDay / inp.operability)}.`,
        `Target the quantity behind the largest line first. ${weatherDays > 0 ? `About ${fmt(weatherDays)} of the ${fmt(grossMain)} gross days are waiting on weather: raising operability from ${fmt(inp.operability)} to ${fmt(Math.min(1, inp.operability + 0.1))} (higher limiting sea state or a more capable vessel) would release ${fmt(inp.vesselDays / inp.operability - inp.vesselDays / Math.min(1, inp.operability + 0.1))} days of the main spread.` : 'No offshore days are priced.'}`,
        'Quantity-based estimating, AACE RP 34R-05; weather-restricted operations, DNV-ST-N001 (2021).',
        `Screening embodied emissions are ${fmt(co2)} t CO₂e (${fmt(inp.steel * CO2_STEEL)} t from steel, ${fmt(inp.concrete * CO2_CONC)} t from concrete, ${fmt(fuel * CO2_MGO)} t from ${fmt(fuel)} t of marine fuel; ${SCREEN}). At ${inp.carbonPrice} ${CUR}/t that is a shadow cost of ${fM(carbonCost)}, ${pct(carbonCost / est.total)} of CAPEX — every 1000 t of steel saved avoids about ${fmt(1000 * CO2_STEEL)} t CO₂e.`),
      adv('economic', est.total <= inp.budgetCap ? (rng.wide.high <= inp.budgetCap ? 'ok' : 'warn') : 'fail',
        `Class ${cls} estimate of ${fM(est.total)} with an expected range of ${fM(rng.wide.low)} to ${fM(rng.wide.high)}; contingency is ${pct(contShare)} of the base estimate. Budget ${fM(inp.budgetCap)}.`,
        est.total > inp.budgetCap ? `The point estimate exceeds the budget by ${fM(est.total - inp.budgetCap)}: reduce scope or raise funding before further commitment.`
          : rng.wide.high > inp.budgetCap ? `The budget covers the estimate but not the upper range (short by ${fM(rng.wide.high - inp.budgetCap)}): mature the estimate to Class ${Math.max(1, cls - 1)} or secure that headroom before sanction.`
            : 'Budget covers the upper accuracy range. Replace the placeholder rates with quotations and re-run before using the figure for commitment.',
        'AACE International RP 18R-97 Table 1 (estimate classes and accuracy ranges) and RP 17R-97; contingency per AACE RP 41R-08.'),
    ];
    if (contShare < Math.abs(AACE_CLASSES[cls].high[0]) * 0.5) advice.push(adv('economic', 'warn', `Contingency of ${pct(contShare)} is small against the +${pct(AACE_CLASSES[cls].high[0])} to +${pct(AACE_CLASSES[cls].high[1])} high-side accuracy of a Class ${cls} estimate.`,
      'Derive contingency from the probabilistic cost-risk module (cumulative P50 to P80 of total cost, P80 being the cost not exceeded with 80 % probability) rather than a flat percentage.', 'AACE RP 41R-08 (range estimating) and RP 18R-97.'));

    return { kpis, checks, plots, tables, advice,
      notes: [`Rates are in ${inp.baseYear} money; the total is escalated ${years} year(s) to ${inp.spendYear} at ${inp.escalation} %/yr.`,
        `Marine durations: ${fmt(inp.vesselDays)} net days at ${fmt(inp.operability)} operability plus ${fmt(inp.transitDays)} transit days give ${fmt(grossMain)} gross days for the longest-hired vessel.`,
        'Every unit rate, day rate and percentage is an indicative placeholder entered by the user, not a market fact.'],
      publish: { 'econ.capex': est.total, 'econ.capexDirect': est.direct, 'econ.installCost': cat('Installation'), 'econ.materialCost': cat('Materials'), 'econ.capexPerMW': bm.perMW ?? 0, 'econ.embodiedCO2': co2 },
      quantities: { steel_t: inp.steel, concrete_m3: inp.concrete, rock_t: inp.rock, cable_km: inp.cable, pipe_km: inp.pipe, vessel_days: { installation: round(vesselDayTotal) } } };
  },
  verify: [
    // 1000 t × 4000 = 4.000 M direct; engineering 10 % → base 4.400 M; contingency 10 % → 4.840 M; escalation 3 %/yr × 2 yr → × 1.0609 = 5.134756 M
    { name: 'Single line: 1000 t × 4000/t, 10 % indirect, 10 % contingency, 3 %/yr over 2 years', source: 'Hand calculation: 4.0 × 1.10 × 1.10 × 1.03² = 5.134756 M',
      inputs: { steel: 1000, rSteel: 4000, concrete: 0, rock: 0, cable: 0, pipe: 0, piles: 0, coating: 0, equipment: [], fabrication: [], spread: [], vesselDays: 0, mobDemob: 0, pEngineering: 10, pManagement: 0, pInsurance: 0, pOwner: 0, pCommissioning: 0, pContingency: 10, escalation: 3, baseYear: 2026, spendYear: 2028 },
      expect: [{ kpi: 'capex', value: 5.134756, tol: 1e-5 }, { kpi: 'direct', value: 4.0, tol: 1e-6 }, { kpi: 'co2', value: 1900, tol: 1e-6 }] },
    // one vessel, 100 net days at 0.8 operability + 5 transit days = 130 gross days × 200 000 = 26.0 M, + 2.0 M mobilisation = 28.0 M direct
    { name: 'Vessel day-rate model: 100 net days at 80 % operability', source: 'Hand calculation: (100/0.8 + 5) × 0.2 M + 2 M = 28 M',
      inputs: { steel: 0, concrete: 0, rock: 0, cable: 0, pipe: 0, piles: 0, coating: 0, equipment: [], fabrication: [], spread: [{ vessel: 'HLV', count: 1, dayRate: 200000, share: 1 }], vesselDays: 100, operability: 0.8, transitDays: 5, mobDemob: 2e6, pEngineering: 0, pManagement: 0, pInsurance: 0, pOwner: 0, pCommissioning: 0, pContingency: 0, escalation: 0 },
      expect: [{ kpi: 'installation', value: 28, tol: 1e-6 }, { kpi: 'grossDays', value: 130, tol: 1e-6 }, { kpi: 'capex', value: 28, tol: 1e-6 }] },
  ],
};
