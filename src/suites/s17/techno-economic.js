// Integrated techno-economic model: the engineering results on the project bus (material quantities, vessel days
// and operability, rated power, capacity factor, availability, design and fatigue life) are carried through cost,
// yield, reliability-informed maintenance, after-tax discounted cash flow, debt cover, a seeded stochastic
// simulation, a physics-based life-cycle carbon balance, a Bayesian update of the cost estimate and a coupled
// design–cost optimisation of the structural section against fatigue life.
import { bottomUpEstimate, AACE_CLASSES } from '../../engine/econ/cost.js';
import { energyFromCapacityFactor } from '../../engine/econ/lcoe.js';
import { embodiedCarbon, marineOpsEmissions, endOfLifeCredit, lifeCycleGwp, carbonIntensity, avoidedEmissions, taxonomyThresholdCheck } from '../../engine/econ/carbon.js';
import { bayesNormalUpdate, priorFromEstimateClass } from '../../engine/econ/forecast.js';
import { surrogateOptimise, multiObjective } from '../../engine/econ/optimise.js';
import { optimalAgeReplacement } from '../../engine/ram/maintenance.js';
import { goldenSection } from '../../engine/math/optim.js';
import { normCdf } from '../../engine/math/stats.js';
import { dcf, fastFlow, lcoeOf } from './_project.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, PCT_VALUE, SCREEN, need, run, round, fmt, mil, fM, pct, adv, check, checkMin, cost, tableRows, uniqueNames, threePoint, sampleLHS, percentiles, thinSample, irrReport, linspace, range, crf, sum } from './_common.js';

const QUOTES = [{ name: 'Fabrication tender A', value: 106 }, { name: 'Fabrication tender B', value: 113 }];
const S_LO = 0.75; const S_HI = 1.35; // range of the section scale factor explored by the design optimisation

export default {
  id: 's17-techno-economic', suite: 17, title: 'Integrated techno-economic model (design to value)',
  summary: 'Chains the engineering results of the other suites — quantities, vessel days, yield, availability, design and fatigue life — into capital cost, energy, operating cost, after-tax value, debt cover, risk, life-cycle carbon and the cost-optimal structural section, on one consistent basis.',
  covers: [12], fidelity: 'reduced-order',
  standards: ['ISO 15686-5:2017 (whole-life costing)', 'AACE International RP 18R-97 (estimate classes and accuracy ranges)', 'IEA/NEA (2020), Projected Costs of Generating Electricity (LCOE convention)', 'Yescombe, Principles of Project Finance, 2nd ed. (2014), ch. 12 (cover ratios)', 'EN 15804:2012+A2:2019 and ISO 14040/14044 (life-cycle assessment modules A–D)', 'Commission Delegated Regulation (EU) 2021/2139, Annex I (life-cycle threshold for electricity generation)', 'DNV-RP-C203 (2021) §2.4 and DNV-ST-0126 (2021) §4.11 (S-N fatigue and design fatigue factors)', 'Jones, Schonlau & Welch (1998) efficient global optimisation; Deb et al. (2002) NSGA-II'],
  models: ['Quantity-based capital cost from engineering quantities', 'Vessel day-rate installation cost with weather operability', 'Capacity-factor and availability-adjusted energy yield with degradation', 'Reliability-informed maintenance planning (Weibull age replacement)', 'After-tax discounted cash flow, real terms', 'Levelised cost of energy', 'Debt-service cover ratio (screening)', 'Stochastic project-finance simulation (Latin hypercube Monte Carlo)', 'Physics-based life-cycle carbon assessment', 'Carbon-adjusted NPV and taxonomy threshold', 'Bayesian updating of the cost estimate', 'Coupled design–cost optimisation against an S-N fatigue constraint', 'Probabilistic fatigue life (lognormal) and expected failure cost', 'Surrogate-assisted (Gaussian-process) optimisation', 'Multi-objective optimisation and Pareto front (return, carbon, safety)', 'One-at-a-time sensitivity to engineering inputs'],
  limits: `A reduced-order chain: each link is the simplest model that keeps the coupling (cost linear in quantities, yield linear in capacity factor and availability, fatigue life proportional to the section factor to the power of the S-N slope, lognormal fatigue life). It shows how an engineering change propagates to value and ranks the drivers; the detailed modules of this suite and the engineering suites supply the numbers. Debt cover is a screening ratio on an annuity loan, not a lender model. Uncertain drivers are sampled independently. ${CUR_NOTE}`,
  publishes: ['econ.technoCapex', 'econ.technoOpexAnnual', 'econ.technoNpv', 'econ.technoLcoe', 'econ.technoMinDscr', 'econ.carbonIntensity', 'econ.carbonPayback', 'econ.lifecycleCO2', 'econ.optimalSectionFactor', 'econ.steelSaving'],
  searchTerms: 'techno-economic analysis offshore wind coupled design cost optimisation life cycle assessment LCOE',
  mesh: { key: 'nSamples', label: 'Monte Carlo samples', levels: [250, 500, 1000, 2000], qoi: 'npvMean' },
  inputs: [
    { key: 'baseYear', label: 'Price base year (time zero)', unit: 'yr', type: 'int', default: 2026, min: 1990, max: 2100, group: 'Opening financial state', role: 'ic', help: 'All money is constant money of this year; the discount rate is real.' },
    { key: 'inflation', label: 'General inflation', unit: '%/yr', type: 'number', default: 2, min: -2, max: 30, group: 'Opening financial state', role: 'ic', help: 'Used for tax depreciation and debt service, which are fixed in money of the day.' },
    { key: 'aaceClass', label: 'Class of the quantity-based estimate (AACE 18R-97)', type: 'select', default: '4', group: 'Opening financial state', role: 'ic', options: [5, 4, 3, 2, 1].map((c) => ({ value: String(c), label: `Class ${c} — ${AACE_CLASSES[c].usage}` })),
      help: 'Sets the uncertainty of the capital cost in the simulation and the prior of the Bayesian update.' },

    { key: 'steel', label: 'Structural steel', unit: 'tonnes', type: 'number', default: 62000, min: 0, max: 5e6, group: 'Engineering quantities', role: 'geometry', bind: 'quantities.steel_t', design: true, uncertain: { dist: 'normal', cov: 0.07 }, help: 'From the structural suites.' },
    { key: 'concrete', label: 'Concrete', unit: 'm³', type: 'number', default: 4500, min: 0, max: 5e6, group: 'Engineering quantities', role: 'geometry', bind: 'quantities.concrete_m3' },
    { key: 'rock', label: 'Rock and scour protection', unit: 'tonnes of rock', type: 'number', default: 150000, min: 0, max: 5e7, group: 'Engineering quantities', role: 'geometry', bind: 'quantities.rock_t' },
    { key: 'cable', label: 'Subsea cable', unit: 'km of cable', type: 'number', default: 165, min: 0, max: 5000, group: 'Engineering quantities', role: 'geometry', bind: 'quantities.cable_km' },
    { key: 'pipe', label: 'Pipeline', unit: 'km of pipe', type: 'number', default: 0, min: 0, max: 5000, group: 'Engineering quantities', role: 'geometry', bind: 'quantities.pipe_km' },
    { key: 'vesselDays', label: 'Net offshore working days', unit: 'vessel-days', type: 'number', default: 210, min: 0, max: 5000, group: 'Engineering quantities', role: 'operation', bind: 'install.vesselDays', design: true, uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Productive days of the installation spread, from the marine-operations suite.' },
    { key: 'operability', label: 'Weather operability', unit: '-', type: 'number', default: 0.68, min: 0.05, max: 1, group: 'Engineering quantities', role: 'environment', bind: 'install.operabilityFrac', help: 'Workable fraction of time; gross days = net days / operability.' },

    { key: 'ratedMW', label: 'Installed capacity', unit: 'MW', type: 'number', default: 500, min: 0.1, max: 20000, group: 'Energy yield', role: 'operation', bind: 'energy.ratedPower' },
    { key: 'yieldSource', label: 'Energy yield from', type: 'select', default: 'cf', group: 'Energy yield', role: 'control', options: [{ value: 'cf', label: 'Capacity factor × availability' }, { value: 'aep', label: 'Net annual production' }] },
    { key: 'capacityFactor', label: 'Capacity factor before unavailability', unit: '-', type: 'number', default: 0.51, min: 0.01, max: 1, group: 'Energy yield', role: 'operation', bind: 'energy.capacityFactor', design: true, uncertain: { dist: 'normal', cov: 0.07 }, showIf: { key: 'yieldSource', in: ['cf'] },
      help: 'Net of wake and electrical losses, at 100 % availability. From the energy-yield suite.' },
    { key: 'availability', label: 'Production availability', unit: '%', type: 'number', default: 96, min: 10, max: 100, group: 'Energy yield', role: 'operation', bind: 'energy.availability', design: true, showIf: { key: 'yieldSource', in: ['cf'] }, help: 'From the reliability and maintenance analysis.' },
    { key: 'aep', label: 'Net annual production, first year', unit: 'GWh/yr', type: 'number', default: 2145, min: 0.001, max: 2e5, group: 'Energy yield', role: 'operation', bind: 'energy.aep', showIf: { key: 'yieldSource', in: ['aep'] } },
    { key: 'degradation', label: 'Annual output degradation', unit: '%/yr', type: 'number', default: 0.3, min: 0, max: 10, group: 'Energy yield', role: 'operation' },

    { key: 'designLife', label: 'Design (operating) life', unit: 'yr', type: 'int', default: 25, min: 1, max: 60, group: 'Design life and fatigue', role: 'bc', bind: 'life.designYears', design: true, help: 'Operational constraint from the engineering suites: the years of production the value rests on.' },
    { key: 'fatigueLife', label: 'Calculated fatigue life of the governing detail', unit: 'yr', type: 'number', default: 95, min: 0.1, max: 1e5, group: 'Design life and fatigue', role: 'bc', bind: 'life.fatigueYears', help: 'Engineering performance limit at the present design, with the design S-N curve.' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Design life and fatigue', role: 'criteria', help: 'Required ratio of calculated fatigue life to design life for the inspection and consequence class.' },
    { key: 'snSlope', label: 'S-N curve slope m', unit: '-', type: 'number', default: 4, min: 2, max: 6, group: 'Design life and fatigue', role: 'material', help: 'Fatigue life varies with stress range to the power −m; 3 for the upper branch of welded-joint curves, 5 for the lower.' },
    { key: 'lifeBias', label: 'Median fatigue life over calculated life', unit: '-', type: 'number', default: 2.5, min: 1, max: 20, group: 'Design life and fatigue', role: 'material', help: 'The design S-N curve is the mean minus two standard deviations, so the median life is longer than the calculated one.' },
    { key: 'lifeCov', label: 'Coefficient of variation of fatigue life', unit: '-', type: 'number', default: 0.9, min: 0.05, max: 3, group: 'Design life and fatigue', role: 'material', help: 'Scatter of life from S-N data and load uncertainty together.' },
    { key: 'steelScaled', label: 'Share of steel that scales with the section', unit: '%', type: 'number', default: 70, min: 1, max: 100, group: 'Design life and fatigue', role: 'geometry', help: 'Primary steel whose mass follows the section scale factor; secondary steel does not.' },
    cost({ key: 'failCost', label: 'Cost of a fatigue failure of the governing detail', unit: CUR, default: 6.0e7, min: 0, max: 1e12, group: 'Design life and fatigue', help: 'Offshore repair, lost production and consequence.' }),

    cost({ key: 'rSteel', label: 'Fabricated steel', unit: `${CUR}/t`, default: 3200, min: 0, max: 1e5, group: 'Unit rates', uncertain: { dist: 'lognormal', cov: 0.2 } }),
    cost({ key: 'rConcrete', label: 'Concrete in place', unit: `${CUR}/m³`, default: 650, min: 0, max: 1e5, group: 'Unit rates' }),
    cost({ key: 'rRock', label: 'Rock supplied and placed', unit: `${CUR}/t`, default: 55, min: 0, max: 1e4, group: 'Unit rates' }),
    cost({ key: 'rCable', label: 'Cable supplied and laid', unit: `${CUR}/km`, default: 8.2e5, min: 0, max: 1e8, group: 'Unit rates' }),
    cost({ key: 'rPipe', label: 'Pipeline supplied and laid', unit: `${CUR}/km`, default: 1.25e6, min: 0, max: 1e8, group: 'Unit rates' }),
    cost({ key: 'dayRate', label: 'Installation spread day rate', unit: `${CUR}/day`, default: 5.4e5, min: 0, max: 1e8, group: 'Unit rates', help: 'All vessels of the spread together, per gross day.' }),
    cost({ key: 'equipPerMW', label: 'Generation and electrical equipment', unit: `${CUR}/MW`, default: 1.5e6, min: 0, max: 1e8, group: 'Unit rates', help: 'Turbines or process equipment and electrical plant per megawatt installed.' }),
    cost({ key: 'pIndirect', label: 'Indirect costs', unit: '% of direct', default: 17.5, min: 0, max: 100, group: 'Unit rates', help: 'Engineering, management, insurance, owner’s costs and commissioning.' }),
    cost({ key: 'pContingency', label: 'Contingency', unit: '% of base estimate', default: 12, min: 0, max: 100, group: 'Unit rates' }),
    cost({ key: 'opexPerMW', label: 'Fixed operating cost', unit: `${CUR}/MW/yr`, default: 9.0e4, min: 0, max: 1e7, group: 'Unit rates', help: 'Excluding the component replacement policy below.' }),
    cost({ key: 'opexVar', label: 'Variable operating cost', unit: `${CUR}/MWh`, default: 1.5, min: 0, max: 1e4, group: 'Unit rates' }),
    cost({ key: 'pDecom', label: 'Decommissioning cost', unit: '% of CAPEX', default: 8, min: 0, max: 100, group: 'Unit rates', help: 'Terminal condition: spent in the year after the last operating year.' }),
    cost({ key: 'price', label: 'Unit price (inflation-indexed)', unit: `${CUR}/MWh`, default: 108, min: 0.01, max: 1e5, group: 'Unit rates', uncertain: { dist: 'normal', cov: 0.12 } }),

    { key: 'units', label: 'Wearing components of one kind', unit: '-', type: 'int', default: 34, min: 0, max: 5000, group: 'Reliability and maintenance', role: 'operation', help: 'For example main bearings or gearboxes. 0 leaves the replacement policy out.' },
    { key: 'mShape', label: 'Weibull shape of component life', unit: '-', type: 'number', default: 2.6, min: 0.3, max: 10, group: 'Reliability and maintenance', role: 'ic', help: 'Opening reliability assumption; above 1 the component wears out and planned replacement can pay.' },
    { key: 'mScale', label: 'Weibull characteristic life', unit: 'yr', type: 'number', default: 14, min: 0.1, max: 200, group: 'Reliability and maintenance', role: 'ic' },
    cost({ key: 'cPrev', label: 'Planned replacement cost per component', unit: CUR, default: 4.2e5, min: 1, max: 1e10, group: 'Reliability and maintenance' }),
    cost({ key: 'cFail', label: 'Unplanned failure cost per component', unit: CUR, default: 1.6e6, min: 1, max: 1e10, group: 'Reliability and maintenance', help: 'Including emergency intervention and lost production.' }),

    { key: 'constructionYears', label: 'Construction period', unit: 'yr', type: 'int', default: 3, min: 1, max: 10, group: 'Fiscal and financial constraints', role: 'operation' },
    { key: 'discountRate', label: 'Real discount rate (hurdle rate)', unit: '%/yr', type: 'number', default: 6, min: -5, max: 50, group: 'Fiscal and financial constraints', role: 'bc', help: `Minimum acceptable real return. ${PLACEHOLDER}` },
    { key: 'taxRate', label: 'Corporate tax rate', unit: '%', type: 'number', default: 25, min: 0, max: 90, group: 'Fiscal and financial constraints', role: 'bc' },
    { key: 'depLife', label: 'Tax depreciation life (straight line)', unit: 'yr', type: 'int', default: 15, min: 1, max: 60, group: 'Fiscal and financial constraints', role: 'bc' },
    { key: 'gearing', label: 'Debt share of funding', unit: '%', type: 'number', default: 65, min: 0, max: 95, group: 'Fiscal and financial constraints', role: 'bc' },
    cost({ key: 'debtRate', label: 'Interest rate on debt (nominal)', unit: '%/yr', default: 5.5, min: 0, max: 30, group: 'Fiscal and financial constraints' }),
    { key: 'tenor', label: 'Debt tenor', unit: 'yr', type: 'int', default: 18, min: 1, max: 40, group: 'Fiscal and financial constraints', role: 'bc' },
    { key: 'minDscr', label: 'Minimum DSCR covenant', unit: '-', type: 'number', default: 1.2, min: 1, max: 3, group: 'Fiscal and financial constraints', role: 'bc' },
    { key: 'budgetCap', label: 'Available capital budget', unit: CUR, type: 'number', default: 2.0e9, min: 1, max: 1e13, group: 'Fiscal and financial constraints', role: 'bc' },

    cost({ key: 'efSteel', label: 'Emission factor, fabricated steel', unit: 't CO₂e/t', default: 1.9, min: 0, max: 10, group: 'Life-cycle carbon', help: 'Screening cradle-to-gate factor.' }),
    cost({ key: 'efConcrete', label: 'Emission factor, concrete', unit: 't CO₂e/m³', default: 0.3, min: 0, max: 2, group: 'Life-cycle carbon', help: 'Screening cradle-to-gate factor.' }),
    cost({ key: 'efCable', label: 'Emission factor, subsea cable', unit: 't CO₂e/km', default: 95, min: 0, max: 5000, group: 'Life-cycle carbon', help: 'Conductor, insulation and armour per kilometre.' }),
    cost({ key: 'efEquip', label: 'Emission factor, generation equipment', unit: 't CO₂e/MW', default: 420, min: 0, max: 1e4, group: 'Life-cycle carbon' }),
    cost({ key: 'fuelPerDay', label: 'Fuel burnt by the spread', unit: 't MGO/day', default: 55, min: 0, max: 1000, group: 'Life-cycle carbon', help: 'Marine gas oil per gross day; 3.206 t CO₂ per tonne burnt (IMO MEPC.308(73)).' }),
    cost({ key: 'opEmissions', label: 'Operational emissions', unit: 't CO₂e/yr', default: 9000, min: 0, max: 1e8, group: 'Life-cycle carbon', help: 'Service vessels, helicopters and imported power.' }),
    cost({ key: 'recovery', label: 'Steel recovered for recycling at end of life', unit: '%', default: 90, min: 0, max: 100, group: 'Life-cycle carbon' }),
    cost({ key: 'recycled', label: 'Recycled content of the steel supplied', unit: '%', default: 20, min: 0, max: 100, group: 'Life-cycle carbon' }),
    cost({ key: 'efPrimary', label: 'Emission factor, primary steel', unit: 't CO₂e/t', default: 2.3, min: 0, max: 10, group: 'Life-cycle carbon' }),
    cost({ key: 'efSecondary', label: 'Emission factor, secondary steel', unit: 't CO₂e/t', default: 0.5, min: 0, max: 10, group: 'Life-cycle carbon' }),
    cost({ key: 'gridIntensity', label: 'Carbon intensity of displaced generation', unit: 'g CO₂e/kWh', default: 350, min: 0, max: 1500, group: 'Life-cycle carbon', help: 'Benchmark for avoided emissions; location and year dependent.' }),
    cost({ key: 'taxonomy', label: 'Life-cycle intensity threshold', unit: 'g CO₂e/kWh', default: 100, min: 1, max: 1500, group: 'Life-cycle carbon', help: 'Benchmark of the substantial-contribution test for electricity generation.' }),
    cost({ key: 'carbonPrice', label: 'Carbon price', unit: `${CUR}/t CO₂e`, default: 90, min: 0, max: 5000, group: 'Life-cycle carbon' }),

    { key: 'quotes', label: 'Quotations received (percent of the estimate)', type: 'table', group: 'Uncertainty', role: 'cost', help: `Each quotation for the whole scope as a percentage of the estimate above; they update the estimate by Bayes’ rule. Leave empty for none. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Source', type: 'text' }, { key: 'value', label: 'Quotation', unit: '% of estimate', type: 'number' }], default: QUOTES },
    { key: 'quoteStd', label: 'Scatter of a quotation about the true cost', unit: '%', type: 'number', default: 8, min: 0.5, max: 100, group: 'Uncertainty', role: 'ic' },
    { key: 'covPrice', label: 'Uncertainty of the price (CoV)', unit: '%', type: 'number', default: 12, min: 0, max: 80, group: 'Uncertainty', role: 'ic', help: 'Opening uncertainty state of the whole-life average price.' },
    { key: 'covYield', label: 'Uncertainty of the energy yield (CoV)', unit: '%', type: 'number', default: 7, min: 0, max: 60, group: 'Uncertainty', role: 'ic' },
    { key: 'covOpex', label: 'Uncertainty of the operating cost (CoV)', unit: '%', type: 'number', default: 15, min: 0, max: 80, group: 'Uncertainty', role: 'ic' },
    { key: 'maxProbLoss', label: 'Tolerated probability of negative NPV', unit: '%', type: 'number', default: 30, min: 0.1, max: 100, group: 'Uncertainty', role: 'criteria' },
    { key: 'nSamples', label: 'Monte Carlo samples', unit: '-', type: 'int', default: 1000, min: 100, max: 20000, group: 'Analysis control', role: 'control' },
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 17, min: 1, max: 1e6, group: 'Analysis control', role: 'control', help: 'Fixed seed: identical inputs reproduce identical results.' },
  ],
  compute(inp) {
    const N = Math.round(inp.designLife); const nc = Math.round(inp.constructionYears); const T = inp.taxRate / 100;
    const infl = inp.inflation / 100; const rReal = inp.discountRate / 100; const rNom = (1 + rReal) * (1 + infl) - 1;
    need(rReal > -1, 'The discount rate must be above −100 %.');
    need(inp.tenor <= N, `The debt tenor of ${inp.tenor} years exceeds the ${N}-year operating life.`);
    const quotes = tableRows(inp.quotes, 'Quotations', { value: { min: 1, max: 1000 } }, { min: 0 });
    uniqueNames(quotes, 'Quotations');
    const loading = (1 + inp.pIndirect / 100) * (1 + inp.pContingency / 100);

    /* ---- reliability-informed maintenance policy for the wearing components ---- */
    let policy = null;
    if (inp.units > 0) {
      need(inp.cFail >= inp.cPrev, 'Reliability and maintenance: an unplanned failure cannot cost less than a planned replacement.');
      policy = run(() => optimalAgeReplacement({ shape: inp.mShape, scale: inp.mScale, costPreventive: inp.cPrev, costFailure: inp.cFail }), 'Replacement policy');
    }
    const maintAnnual = policy ? inp.units * policy.costRate : 0;
    const maintRtf = policy ? inp.units * policy.runToFailureCostRate : 0;

    /* ---- the chain: engineering values → capex, energy, opex → after-tax cash flow ---- */
    // `o` overrides engineering inputs (for the sensitivity and the design optimisation); returns the whole case
    const build = (o = {}) => {
      const v = { ...inp, ...o };
      const gross = v.vesselDays > 0 ? v.vesselDays / v.operability : 0;
      const items = [];
      const add = (name, quantity, unitRate, unit, category) => { if (quantity > 0 && unitRate > 0) items.push({ name, quantity, unitRate, unit, discipline: category, category }); };
      add('Structural steel', v.steel, v.rSteel, 't', 'Structure'); add('Concrete', v.concrete, v.rConcrete, 'm³', 'Structure'); add('Rock and scour protection', v.rock, v.rRock, 't', 'Structure');
      add('Subsea cable', v.cable, v.rCable, 'km', 'Cables and pipelines'); add('Pipeline', v.pipe, v.rPipe, 'km', 'Cables and pipelines');
      add('Marine installation', gross, v.dayRate, 'vessel-day', 'Installation'); add('Generation and electrical equipment', v.ratedMW, v.equipPerMW, 'MW', 'Equipment');
      need(items.length > 0, 'Nothing to cost: every quantity or unit rate is zero.');
      const est = bottomUpEstimate({ items, indirects: { Indirects: v.pIndirect / 100 }, contingency: v.pContingency / 100 });
      const e0 = v.yieldSource === 'cf' ? v.ratedMW * 8760 * v.capacityFactor * Math.min(v.availability, 100) / 100 : v.aep * 1e3; // MWh in the first year
      const nOp = Math.round(v.designLife);
      const prod = Float64Array.from(v.yieldSource === 'cf'
        ? energyFromCapacityFactor({ capacityMW: v.ratedMW, capacityFactor: Math.min(1, v.capacityFactor * Math.min(v.availability, 100) / 100), degradation: v.degradation / 100, years: nOp })
        : range(nOp, (k) => e0 * (1 - v.degradation / 100) ** k));
      const capex = est.total;
      const opexFixed = v.opexPerMW * v.ratedMW + maintAnnual;
      const mk = (real) => ({ nc, nOp, capex, phasing: new Float64Array(nc).fill(1 / nc), sunk: 0, prod, price: v.price, indexed: true, opexFixed, opexVar: v.opexVar, opexEsc: 0,
        decom: capex * v.pDecom / 100, taxRate: T, depLife: Math.round(v.depLife), royalty: 0, infl, rate: real ? rReal : rNom, real, delay: 0, baseYear: v.baseYear });
      return { est, gross, e0, prod, capex, opexFixed, nOp, mk };
    };
    const base = run(() => build(), 'Techno-economic chain');
    const p = base.mk(true); const pNom = base.mk(false);
    const cf = dcf(p); const NPV = cf.npv;
    const L = lcoeOf(p, cf); const ir = irrReport(cf.fcff);
    const opex1 = base.opexFixed + inp.opexVar * base.e0;
    const netCF = base.e0 / (inp.ratedMW * 8760);
    // discounted payback from the first row
    let cum = 0; let payback = null;
    for (let t = 0; t < cf.rows; t++) { const d = cf.fcff[t] * (1 + rReal) ** -t; if (payback == null && t > 0 && cum < 0 && cum + d >= 0) payback = t - 1 + (-cum) / d; cum += d; }
    let pvIn = 0; let pvOut = 0;
    for (let t = 0; t < cf.rows; t++) { const d = (1 + rReal) ** -t; pvIn += cf.revenue[t] * d; pvOut += (cf.capex[t] + cf.opex[t] + cf.tax[t] + cf.decom[t] + cf.royalty[t]) * d; }

    /* ---- debt cover (screening): annuity loan on the nominal construction cost ---- */
    const cfN = dcf(pNom);
    const capexNom = sum(cfN.capex); const debt = inp.gearing / 100 * capexNom; const tenor = Math.round(inp.tenor);
    const ds = debt > 0 ? debt * crf(inp.debtRate / 100, tenor) : 0;
    const dscrOf = (flows, service) => { let mn = Infinity; for (let k = 0; k < tenor; k++) mn = Math.min(mn, flows[nc + k] / service); return mn; };
    const minDscr = ds > 0 ? dscrOf(cfN.fcff, ds) : null;

    /* ---- stochastic simulation: price, yield, capital and operating cost ---- */
    const cls = Number(inp.aaceClass); const band = AACE_CLASSES[cls];
    const n = Math.round(inp.nSamples); const vars = {};
    if (inp.covPrice > 0) vars.price = { dist: 'normal', mean: 1, std: inp.covPrice / 100 };
    if (inp.covYield > 0) vars.prod = { dist: 'normal', mean: 1, std: inp.covYield / 100 };
    // the wide accuracy band of the estimate class read as the cumulative P10–P90 of a lognormal
    vars.capex = threePoint('lognormal', 1 + band.low[1], 1, 1 + band.high[1], 'Capital cost range');
    if (inp.covOpex > 0) vars.opex = { dist: 'normal', mean: 1, std: inp.covOpex / 100 };
    const S = run(() => sampleLHS(vars, n, inp.seed, 0), 'Sampling');
    const npvS = new Float64Array(n); const lcoeS = new Float64Array(n); const dscrS = new Float64Array(n);
    const buf = new Float64Array(nc + N + 16); const res = {}; const f = { price: 1, prod: 1, capex: 1, opex: 1 };
    let breach = 0;
    for (let k = 0; k < n; k++) {
      f.price = S.price ? Math.max(S.price[k], 0) : 1; f.prod = S.prod ? Math.max(S.prod[k], 0) : 1; f.capex = S.capex[k]; f.opex = S.opex ? Math.max(S.opex[k], 0) : 1;
      // nominal flows discounted at the nominal rate give the same NPV as real flows at the real rate (Fisher)
      fastFlow(pNom, f, 0, buf, res);
      npvS[k] = res.npv; lcoeS[k] = res.lcoe;
      dscrS[k] = ds > 0 ? dscrOf(buf, ds * f.capex) : Infinity;
      if (ds > 0 && dscrS[k] < inp.minDscr) breach++;
    }
    const PN = percentiles(npvS); const PL = percentiles(lcoeS);
    let mean = 0; for (let k = 0; k < n; k++) mean += npvS[k]; mean /= n;
    let sd = 0; let loss = 0; for (let k = 0; k < n; k++) { sd += (npvS[k] - mean) ** 2; if (npvS[k] < 0) loss++; } sd = Math.sqrt(sd / (n - 1));
    const probLoss = loss / n; const probBreach = ds > 0 ? breach / n : 0;

    /* ---- sensitivity of NPV to the engineering inputs (±10 %, one at a time) ---- */
    const npvOf = (o) => dcf(build(o).mk(true)).npv;
    const drivers = [['Structural steel', 'steel'], ['Offshore vessel days', 'vesselDays'], ['Weather operability', 'operability'], ['Design life', 'designLife'],
      ...(inp.yieldSource === 'cf' ? [['Capacity factor', 'capacityFactor'], ['Availability', 'availability']] : [['Annual production', 'aep']]), ['Cable length', 'cable']]
      .filter(([, k]) => inp[k] > 0)
      .map(([name, key]) => {
        const cap = key === 'operability' || key === 'capacityFactor' ? 1 : key === 'availability' ? 100 : Infinity;
        const lo = key === 'designLife' ? Math.max(inp.tenor, Math.round(inp[key] * 0.9)) : inp[key] * 0.9; const hi = key === 'designLife' ? Math.round(inp[key] * 1.1) : Math.min(inp[key] * 1.1, cap);
        const a = npvOf({ [key]: lo }); const b = npvOf({ [key]: hi });
        return { name, key, lo, hi, low: a, high: b, swing: Math.abs(b - a) };
      }).sort((a, b) => b.swing - a.swing);

    /* ---- coupled design–cost optimisation of the structural section ---- */
    // Section scale factor s on the primary steel: mass ∝ s, stress range ∝ 1/s, so fatigue life L(s) = L0 s^m.
    // Fatigue life is lognormal about the median bias·L(s); a failure in year t costs failCost.
    const m = inp.snSlope; const L0 = inp.fatigueLife; const sig = Math.sqrt(Math.log(1 + inp.lifeCov ** 2));
    const steelS = inp.steel * inp.steelScaled / 100;
    const lifeOf = (s) => L0 * s ** m;
    const pfAt = (s, t) => normCdf((Math.log(t) - Math.log(inp.lifeBias * lifeOf(s))) / sig);
    const efc = (s) => { let pv = 0; let prev = 0; for (let k = 1; k <= N; k++) { const F = pfAt(s, k); pv += (F - prev) * inp.failCost * (1 + rReal) ** -(nc - 1 + k); prev = F; } return pv; };
    const capexOf = (s) => base.capex + steelS * (s - 1) * inp.rSteel * loading;
    const sBuf = new Float64Array(nc + N + 16);
    const pS = base.mk(true);
    const npvS1 = (s) => { const c = capexOf(s); pS.capex = c; pS.decom = c * inp.pDecom / 100; return fastFlow(pS, null, 0, sBuf) - (1 - T) * efc(s); };
    const sReq = (inp.dff * N / L0) ** (1 / m);
    need(sReq <= S_HI, `The calculated fatigue life of ${fmt(L0)} years cannot reach ${fmt(inp.dff)} × ${N} years within a section increase of ${pct(S_HI - 1)}: the detail needs redesign (lower stress concentration, better S-N class or grinding), not more thickness.`);
    const sA = Math.max(S_LO, sReq);
    const gs = goldenSection((s) => -npvS1(s) / 1e6, sA, S_HI);
    const sOpt = -npvS1(sA) / 1e6 <= gs.f ? sA : gs.x; // the optimum may sit on the fatigue limit itself
    const npvOpt = npvS1(sOpt); const npvAsIs = npvS1(1);
    const sur = run(() => surrogateOptimise({ f: (x) => -npvS1(x[0]) / 1e6, bounds: [[sA, S_HI]], nInit: 6, nIter: 10, seed: inp.seed }), 'Surrogate optimisation');
    const steelSaved = steelS * (1 - sOpt);
    const co2Of = (s) => (inp.steel + steelS * (s - 1)) * inp.efSteel;
    const mo = run(() => multiObjective({ objectives: (x) => [-npvS1(x[0]) / 1e6, co2Of(x[0]) / 1e3, Math.log10(Math.max(pfAt(x[0], N), 1e-300))], bounds: [[S_LO, S_HI]], names: ['Value', 'Carbon', 'Failure probability'], pop: 24, generations: 16, seed: inp.seed }), 'Multi-objective optimisation');
    const front = mo.front.slice().sort((a, b) => a.x[0] - b.x[0]);

    /* ---- physics-based life-cycle carbon ---- */
    const emb = embodiedCarbon({ items: [{ name: 'Structural steel', mass_t: inp.steel, factor: inp.efSteel }, { name: 'Concrete', mass_t: inp.concrete, factor: inp.efConcrete }, { name: 'Subsea cable', mass_t: inp.cable, factor: inp.efCable }, { name: 'Generation equipment', mass_t: inp.ratedMW, factor: inp.efEquip }] });
    const marine = marineOpsEmissions({ vessels: [{ vessel: 'Installation spread', days: base.gross, fuel_tPerDay: inp.fuelPerDay }] });
    const eol = endOfLifeCredit({ mass_t: inp.steel, recoveryRate: inp.recovery / 100, recycledContent: inp.recycled / 100, primaryFactor: inp.efPrimary, secondaryFactor: inp.efSecondary });
    // removal is taken as a marine campaign of half the installation effort
    const gwp = lifeCycleGwp({ product_t: emb.total_t, installation_t: marine.total_t, operationPerYear_t: inp.opEmissions, operatingYears: N, endOfLife_t: 0.5 * marine.total_t, credit_t: Math.max(eol.credit_t, 0) });
    const lifeEnergy = sum(base.prod);
    const ci = carbonIntensity({ emissions_t: gwp.total_t, energy_MWh: lifeEnergy });
    const av = avoidedEmissions({ annualEnergy_MWh: lifeEnergy / N, referenceIntensity_gPerKWh: inp.gridIntensity, lifeCycleEmissions_t: gwp.total_t, operatingYears: N });
    const tx = taxonomyThresholdCheck({ intensity_gPerKWh: ci.gPerKWh, threshold_gPerKWh: inp.taxonomy });
    // carbon priced when emitted: construction emissions over the construction years, operation each year, removal at the end
    let pvCarbon = 0;
    for (let t = 0; t < nc; t++) pvCarbon += (emb.total_t + marine.total_t) / nc * inp.carbonPrice * (1 + rReal) ** -t;
    for (let k = 0; k < N; k++) pvCarbon += inp.opEmissions * inp.carbonPrice * (1 + rReal) ** -(nc + k);
    pvCarbon += 0.5 * marine.total_t * inp.carbonPrice * (1 + rReal) ** -(nc + N);
    const npvCarbon = NPV - pvCarbon;

    /* ---- Bayesian update of the cost estimate with quotations ---- */
    const prior = priorFromEstimateClass({ estimate: base.capex, cls });
    const post = bayesNormalUpdate({ priorMean: prior.priorMean, priorStd: prior.priorStd, observations: quotes.map((q) => q.value / 100 * base.capex), obsStd: inp.quoteStd / 100 * base.capex });

    const cat = (k) => base.est.byCategory[k] || 0;
    const kpis = [
      { key: 'capex', label: 'CAPEX from the engineering quantities', value: mil(base.capex), unit: MU, limit: mil(inp.budgetCap), sense: 'max', explain: `Quantities and vessel days at their unit rates, with ${inp.pIndirect} % indirects and ${inp.pContingency} % contingency, in ${inp.baseYear} money. A Class ${cls} figure: it moves one-for-one with the engineering quantities.` },
      { key: 'capexPerMW', label: 'CAPEX per MW', value: round(base.capex / inp.ratedMW / 1e6), unit: `${MU}/MW`, explain: 'Benchmark denominator for comparison with completed projects of the same type and water depth.' },
      { key: 'aep', label: 'Net annual energy, first year', value: round(base.e0 / 1e3), unit: 'GWh/yr', explain: inp.yieldSource === 'cf' ? `Capacity × 8760 h × capacity factor × availability. Each point of availability is worth ${fmt(base.e0 / 1e3 / Math.min(inp.availability, 100))} GWh a year.` : 'Net annual production as supplied by the energy-yield analysis.' },
      { key: 'netCF', label: 'Net capacity factor', value: round(100 * netCF), unit: '%', explain: 'First-year energy over capacity × 8760 h, after availability. The single number that summarises site, layout and reliability.' },
      { key: 'opexAnnual', label: 'Annual OPEX, first year', value: mil(opex1), unit: `${MU}/yr`, explain: `Fixed cost per MW, variable cost per MWh and ${fM(maintAnnual)}/yr for the planned replacement policy of the wearing components.` },
      { key: 'npv', label: 'NPV (real, after tax)', value: mil(NPV), unit: MU, limit: 0, sense: 'min', explain: `Present value at ${inp.baseYear} of the after-tax cash flow at ${inp.discountRate} % real. Positive means the engineering design earns more than the hurdle rate.` },
    ];
    if (ir.value != null) kpis.push({ key: 'irr', label: 'IRR (real, after tax)', value: round(100 * ir.value), unit: '%/yr', limit: inp.discountRate, sense: 'min', explain: `Real return on the whole investment. ${ir.count > 1 ? `The cash flow has ${ir.count} mathematical rates because of the decommissioning outflow; the single positive one is quoted.` : 'Compare with the real hurdle rate.'}` });
    kpis.push(
      { key: 'lcoe', label: 'LCOE (real, pre-tax)', value: round(L.lcoe), unit: `${CUR}/MWh`, limit: inp.price, sense: 'max', explain: 'Discounted capital, operating and decommissioning cost over discounted energy, constant money and real rate (IEA convention). Below the price the project covers its cost of capital before tax.' },
      { key: 'payback', label: 'Discounted payback', value: round(payback ?? cf.rows - 1), unit: 'yr', explain: payback == null ? 'The discounted cumulative cash flow does not turn positive within the life; the horizon is shown.' : 'Years from time zero until the discounted after-tax cash flow has recovered the investment.' });
    if (minDscr != null) kpis.push({ key: 'minDscr', label: 'Minimum DSCR (screening)', value: round(minDscr), unit: '-', limit: inp.minDscr, sense: 'min', explain: `Lowest ratio of after-tax operating cash flow to the service of a ${inp.gearing} % annuity loan over ${tenor} years, in money of the day. Below the covenant the design does not carry this debt.` });
    kpis.push(
      { key: 'npvMean', label: 'Expected NPV (simulation)', value: mil(mean), unit: MU, explain: `Mean of ${n} sampled outcomes with price, yield, capital and operating cost uncertain. Usually below the deterministic NPV because the Class ${cls} cost range is skewed to overrun.` },
      { key: 'npvP10', label: 'NPV P10 (cumulative, downside)', value: mil(PN(10)), unit: MU, explain: 'NPV not exceeded with 10 % probability: the downside case. In the exceedance convention this same number is called P90.' },
      { key: 'npvP90', label: 'NPV P90 (cumulative, upside)', value: mil(PN(90)), unit: MU, explain: 'NPV not exceeded with 90 % probability: the upside case (exceedance P10).' },
      { key: 'probLoss', label: 'Probability of negative NPV', value: round(100 * probLoss), unit: '%', limit: inp.maxProbLoss, sense: 'max', explain: 'Share of sampled outcomes that destroy value at the hurdle rate, against the tolerated probability.' },
      { key: 'probBreach', label: 'Probability of breaching the DSCR covenant', value: round(100 * probBreach), unit: '%', explain: `Share of outcomes in which the minimum cover falls below ${inp.minDscr}. Lenders look at this before the NPV.` },
      { key: 'pfDesign', label: 'Probability of fatigue failure in the design life', value: round(100 * pfAt(1, N), 4), unit: '%', explain: `Chance that the governing detail fails within ${N} years at the present design, from a lognormal life with median ${fmt(inp.lifeBias * L0)} years. It prices the fatigue reserve rather than only checking it.` },
      { key: 'sReq', label: 'Section factor required by the fatigue criterion', value: round(sReq), unit: '-', explain: `Smallest scale of the primary section at which calculated life reaches ${fmt(inp.dff)} × ${N} years, from life ∝ s^${fmt(m)}. Below 1.0 the present design has fatigue reserve that could be given back.` },
      { key: 'sOpt', label: 'Cost-optimal section factor', value: round(sOpt), unit: '-', explain: 'Scale of the primary section that maximises NPV, counting steel cost and the expected cost of fatigue failure, subject to the fatigue criterion. Equal to the required factor when the criterion governs.' },
      { key: 'steelSaved', label: 'Steel change at the optimum', value: round(-steelSaved), unit: 't', explain: 'Change in primary steel from the present design to the cost-optimal section (negative = steel saved). The quantity to take back to the structural design.' },
      { key: 'npvGain', label: 'NPV gain of the optimal section', value: mil(npvOpt - npvAsIs), unit: MU, explain: 'After-tax value added by moving from the present section to the optimal one, net of the change in expected failure cost.' },
      { key: 'surGap', label: 'Surrogate optimum, difference in section factor', value: round(sur.x[0] - sOpt, 4), unit: '-', explain: `The same optimum found by a Gaussian-process surrogate in ${sur.evaluations} model runs, against the exact line search. A small difference shows the surrogate method can be trusted when each run is a full engineering analysis.` },
      { key: 'maintSaving', label: 'Saving of planned replacement over run-to-failure', value: mil(maintRtf - maintAnnual), unit: `${MU}/yr`, explain: policy && Number.isFinite(policy.interval) ? `Fleet saving from replacing the wearing component at ${fmt(policy.interval)} years instead of on failure. It is already inside the operating cost above.` : 'No planned replacement age beats run-to-failure for this component (no wear-out or no cost advantage), so none is assumed.' },
      { key: 'co2', label: 'Life-cycle emissions', value: round(gwp.total_t / 1e3), unit: 'kt CO₂e', explain: `Materials, installation, ${N} years of operation and removal, before the recycling credit of ${fmt(Math.max(eol.credit_t, 0) / 1e3)} kt. Computed from the same quantities as the cost.` },
      { key: 'intensity', label: 'Life-cycle carbon intensity', value: round(ci.gPerKWh), unit: 'g CO₂e/kWh', limit: inp.taxonomy, sense: 'max', explain: 'Life-cycle emissions over life-time energy. Compared with the threshold benchmark entered; passing it is necessary, not sufficient, for taxonomy alignment.' },
      { key: 'npvCarbon', label: 'NPV with carbon priced', value: mil(npvCarbon), unit: MU, explain: `NPV less the present value of life-cycle emissions at ${inp.carbonPrice} ${CUR}/t (${fM(pvCarbon)}). The value if embodied and operational carbon were charged in full.` },
      { key: 'capexPosterior', label: 'CAPEX after the quotations (posterior mean)', value: mil(post.mean), unit: MU, explain: `The Class ${cls} estimate combined with ${quotes.length} quotation(s) by Bayes’ rule; the prior carries ${pct(post.priorWeight)} of the weight. It replaces the estimate as the expected cost.` },
      { key: 'capexPostP90', label: 'CAPEX after the quotations, P90 (cumulative)', value: mil(post.p90), unit: MU, explain: 'Cost not exceeded with 90 % probability after the update: the high case. Quotations narrow the range far more than they move the mean.' },
    );

    // with nothing displaced (reference intensity zero) there is no carbon payback to report
    const paysBack = Number.isFinite(av.carbonPaybackYears);
    if (paysBack) kpis.splice(kpis.findIndex((k) => k.key === 'npvCarbon'), 0, { key: 'carbonPayback', label: 'Carbon payback', value: round(av.carbonPaybackYears), unit: 'yr', explain: `Years of generation needed to avoid as much as the asset emits over its life, against displaced generation at ${inp.gridIntensity} g/kWh.` });
    const checks = [
      check('Design life × design fatigue factor ≤ calculated fatigue life', N * inp.dff, L0, 'yr', 'DNV-RP-C203 (2021) §2.4; DNV-ST-0126 (2021) §4.11 (design fatigue factors)',
        'Guards against fatigue cracking of the governing detail before the end of the life the revenue is counted over. Above 1.0 the value is built on years the structure is not shown to survive.'),
      check('Present value of cost ≤ present value of revenue', pvOut, pvIn, CUR, 'ISO 15686-5:2017 §5.4 (net present value)', 'Guards against a design whose whole-life cost, including tax and decommissioning, is not recovered at the hurdle rate. Above 1.0 the NPV is negative.'),
      check('CAPEX within the available budget', base.capex, inp.budgetCap, CUR, 'Financial constraint (affordability limit)', 'Guards against a design the available capital cannot fund; above 1.0 the quantities or the scope must come down.'),
      check('Probability of loss ≤ tolerance', 100 * probLoss, inp.maxProbLoss, '%', 'Investment risk tolerance (corporate criterion); AACE RP 57R-09', 'Guards against sanctioning a design with a higher chance of destroying value than the investor accepts.'),
      check('Life-cycle carbon intensity ≤ threshold', ci.gPerKWh, inp.taxonomy, 'g CO₂e/kWh', 'Commission Delegated Regulation (EU) 2021/2139, Annex I', 'Guards against losing access to taxonomy-aligned finance. Above 1.0 the asset emits more per unit of energy over its life than the threshold allows.'),
    ];
    if (minDscr != null) checks.push(checkMin('Minimum DSCR ≥ covenant', minDscr, inp.minDscr, '-', 'Yescombe (2014) §12.3 (annual debt-service cover ratio)', 'Guards against cash flow that cannot service the debt in the weakest year; below the covenant the gearing assumed is not bankable.'));

    const years = range(cf.rows, (t) => inp.baseYear + t);
    const ss = linspace(S_LO, S_HI, 41);
    let c2 = 0; const cumDisc = Array.from(cf.fcff, (v, t) => { c2 += v * (1 + rReal) ** -t; return c2; });
    const plots = [
      { type: 'bar', title: 'CAPEX from the engineering quantities', categories: ['Structure', 'Cables and pipelines', 'Installation', 'Equipment', 'Indirects', 'Contingency'], x: { label: 'Block' }, y: { label: 'Cost', unit: MU },
        series: [{ name: 'Cost', y: [cat('Structure'), cat('Cables and pipelines'), cat('Installation'), cat('Equipment'), base.est.indirectTotal, base.est.contingency].map(mil) }],
        explain: 'Where the capital goes, by the engineering quantity that causes it. The tallest block is the suite whose results matter most to cost.' },
      { type: 'tornado', title: 'NPV sensitivity to the engineering inputs (± 10 %)', base: mil(NPV), items: drivers.map((d) => ({ name: d.name, low: mil(Math.min(d.low, d.high)), high: mil(Math.max(d.low, d.high)) })), x: { label: 'NPV', unit: MU },
        explain: 'NPV with each engineering input alone moved by 10 % (design life in whole years, fractions capped at their physical limit). The longest bar is the technical result worth the most further analysis.' },
      { type: 'line', title: 'Cumulative discounted cash flow', x: { label: 'Year', unit: '' }, y: { label: 'Cumulative', unit: MU }, series: [{ name: 'Discounted, after tax', x: years, y: cumDisc.map(mil) }, { name: 'Break-even', x: [years[0], years[years.length - 1]], y: [0, 0], dash: true }],
        explain: 'Running present value of the after-tax cash flow. The lowest point is the capital at risk, the zero crossing the discounted payback and the end value the NPV.' },
      { type: 'hist', title: 'Distribution of NPV', values: thinSample(npvS, 1500).map(mil), bins: 36, x: { label: 'NPV', unit: MU }, y: { label: 'Frequency' }, explain: 'Simulated NPV outcomes. The area to the left of zero is the probability of loss; the long left tail comes from the overrun side of the estimate class.' },
      { type: 'line', title: 'Value and fatigue safety against the section factor', x: { label: 'Section scale factor s', unit: '-' }, y: { label: 'NPV', unit: MU },
        series: [{ name: 'NPV net of expected failure cost', x: ss.map((s) => round(s, 4)), y: ss.map((s) => mil(npvS1(s))) }, { name: 'Fatigue criterion (minimum s)', x: [round(sReq, 4), round(sReq, 4)], y: [mil(npvS1(S_HI)), mil(Math.max(npvS1(S_LO), npvOpt))], dash: true },
          { name: 'Optimum', x: [round(sOpt, 4), round(sOpt, 4)], y: [mil(npvS1(S_HI)), mil(npvOpt)], dash: true }, { name: 'Present design', x: [1, 1], y: [mil(npvS1(S_HI)), mil(npvAsIs)], dash: true }],
        explain: 'Value falls with added steel and, to the left, with the rising expected cost of fatigue failure. The feasible region is to the right of the fatigue-criterion line; the optimum is the highest feasible point.' },
      { type: 'scatter', title: 'Pareto front: value against embodied carbon of the steel', x: { label: 'Embodied carbon of structural steel', unit: 'kt CO₂e' }, y: { label: 'NPV', unit: MU },
        series: [{ name: 'Non-dominated designs', x: front.map((q) => round(q.f[1])), y: front.map((q) => round(-q.f[0])) }, { name: 'Knee', x: [round(mo.knee.f[1])], y: [round(-mo.knee.f[0])] }],
        explain: 'Designs for which no other is better in value, carbon and failure probability together. Along the front more steel buys safety at the price of value and carbon; the knee is the balanced choice.' },
      { type: 'bar', title: 'Life-cycle emissions by phase', categories: Object.keys(gwp.byPhase), x: { label: 'Life-cycle module (EN 15804)' }, y: { label: 'Emissions', unit: 'kt CO₂e' }, series: [{ name: 'Emissions', y: Object.values(gwp.byPhase).map((v) => round(v / 1e3)) }],
        explain: 'Emissions by life-cycle module, with the recycling credit shown negative. Materials usually dominate, so the steel tonnage drives carbon as it drives cost.' },
    ];
    const tables = [
      { title: 'Project economic dashboard', columns: [{ key: 'area', label: 'Area' }, { key: 'item', label: 'Measure' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }, { key: 'src', label: 'Driven by' }],
        rows: [['Cost', 'CAPEX', mil(base.capex), MU, 'quantities, vessel days, operability'], ['Cost', 'OPEX, first year', mil(opex1), `${MU}/yr`, 'capacity, maintenance policy'], ['Cost', 'Decommissioning', mil(p.decom), MU, 'share of CAPEX'],
          ['Yield', 'Net energy, first year', round(base.e0 / 1e3), 'GWh/yr', 'capacity factor, availability'], ['Yield', 'Net capacity factor', round(100 * netCF, 4), '%', 'site, layout, reliability'],
          ['Value', 'NPV, real after tax', mil(NPV), MU, 'all of the above'], ['Value', 'IRR, real after tax', ir.value != null ? round(100 * ir.value, 4) : '–', '%/yr', ''], ['Value', 'LCOE, real pre-tax', round(L.lcoe, 5), `${CUR}/MWh`, ''], ['Value', 'Discounted payback', payback != null ? round(payback, 4) : '–', 'yr', ''],
          ['Finance', 'Minimum DSCR', minDscr != null ? round(minDscr, 4) : '–', '-', `${inp.gearing} % debt, ${tenor} yr`], ['Risk', 'Expected NPV', mil(mean), MU, `Class ${cls} cost range, price, yield`], ['Risk', 'NPV P10 / P90 (cumulative)', `${fmt(PN(10) / 1e6)} / ${fmt(PN(90) / 1e6)}`, MU, ''],
          ['Risk', 'Probability of loss', round(100 * probLoss, 4), '%', ''], ['Risk', 'LCOE P10 / P90 (cumulative)', `${fmt(PL(10))} / ${fmt(PL(90))}`, `${CUR}/MWh`, 'P90 is the high-cost case'], ['Integrity', 'Fatigue utilisation of life', round(N * inp.dff / L0, 4), '-', 'design life × DFF / fatigue life'],
          ['Integrity', 'Optimal section factor', round(sOpt, 5), '-', 'steel rate, failure cost, S-N slope'], ['Carbon', 'Life-cycle emissions', round(gwp.total_t / 1e3, 5), 'kt CO₂e', 'steel, concrete, fuel'], ['Carbon', 'Intensity', round(ci.gPerKWh, 4), 'g CO₂e/kWh', ''], ['Carbon', 'Carbon payback', paysBack ? round(av.carbonPaybackYears, 4) : '–', 'yr', 'displaced generation']]
          .map(([area, item, value, unit, src]) => ({ area, item, value, unit, src })),
        explain: 'The headline results of the whole chain on one page, each with the engineering input that drives it: the summary for an investment decision report.' },
      { title: 'Lifecycle financial forecast (real, after tax)', columns: [{ key: 'y', label: 'Year' }, { key: 'e', label: 'Energy', unit: 'GWh' }, { key: 'rev', label: 'Revenue', unit: MU }, { key: 'opex', label: 'Opex', unit: MU }, { key: 'tax', label: 'Tax', unit: MU }, { key: 'capex', label: 'Capex and decommissioning', unit: MU }, { key: 'fcf', label: 'Free cash flow', unit: MU }, { key: 'cum', label: 'Cumulative discounted', unit: MU }],
        rows: years.map((yv, t) => ({ y: yv, e: round(cf.energy[t] / 1e3, 5), rev: mil(cf.revenue[t]), opex: mil(cf.opex[t]), tax: mil(cf.tax[t]), capex: mil(cf.capex[t] + cf.decom[t]), fcf: mil(cf.fcff[t]), cum: mil(cumDisc[t]) })),
        explain: 'Year-by-year forecast in constant money: use it to check that capital, production, tax and decommissioning fall in their intended years.' },
      { title: 'Engineering drivers of value', columns: [{ key: 'name', label: 'Engineering input' }, { key: 'base', label: 'Value' }, { key: 'lo', label: 'NPV at −10 %', unit: MU }, { key: 'hi', label: 'NPV at +10 %', unit: MU }, { key: 'el', label: 'NPV change per 1 % of input', unit: MU }],
        rows: drivers.map((d) => ({ name: d.name, base: round(inp[d.key], 5), lo: mil(d.low), hi: mil(d.high), el: mil((d.high - d.low) / (100 * (d.hi - d.lo) / inp[d.key])) })),
        explain: 'Each engineering input with the NPV at its low and high value and the value of one per cent of it. The last column is the budget for the analysis or design change that would secure that one per cent.' },
      { title: 'Life-cycle carbon inventory', columns: [{ key: 'name', label: 'Item' }, { key: 'q', label: 'Quantity' }, { key: 'f', label: 'Factor' }, { key: 'e', label: 'Emissions', unit: 't CO₂e' }],
        rows: [...emb.lines.map((l) => ({ name: l.name, q: round(l.mass_t), f: l.factor, e: round(l.emissions_t) })), { name: 'Installation fuel', q: round(marine.fuel_t), f: 3.206, e: round(marine.total_t) }, { name: `Operation, ${N} years`, q: N, f: inp.opEmissions, e: round(inp.opEmissions * N) },
          { name: 'Removal (half the installation effort)', q: round(0.5 * marine.fuel_t), f: 3.206, e: round(0.5 * marine.total_t) }, { name: 'Recycling credit (module D)', q: round(eol.netScrap_t), f: round(inp.efPrimary - inp.efSecondary, 4), e: -round(Math.max(eol.credit_t, 0)) }],
        explain: 'Every emission source as quantity times factor, from the same quantities that are costed. Replace the screening factors with supplier declarations as they become available.' },
    ];

    const top = drivers[0];
    const co2Saved = steelSaved * inp.efSteel;
    const fatigueGoverns = Math.abs(sOpt - sReq) < 1e-3 * sReq;
    const advice = [
      adv('economic', NPV >= 0 && probLoss * 100 <= inp.maxProbLoss && (minDscr == null || minDscr >= inp.minDscr) ? 'ok' : NPV < 0 ? 'fail' : 'warn',
        `From the engineering results: CAPEX ${fM(base.capex)} (Class ${cls}), OPEX ${fM(opex1)}/yr, ${fmt(base.e0 / 1e3)} GWh/yr; NPV ${fM(NPV)} real after tax, LCOE ${fmt(L.lcoe)} ${CUR}/MWh against a price of ${inp.price}${minDscr != null ? `, minimum DSCR ${fmt(minDscr)}` : ''}. Simulation: expected NPV ${fM(mean)}, cumulative P10 / P90 ${fM(PN(10))} / ${fM(PN(90))}, probability of loss ${pct(probLoss)}.`,
        NPV < 0 ? `The design does not recover its cost: LCOE must fall by ${fmt(L.lcoe - inp.price)} ${CUR}/MWh or the price rise to ${fmt(L.lcoe)}. Start with "${top.name}", worth ${fM(Math.abs(top.high - top.low) / 20)} per percent.`
          : probLoss * 100 > inp.maxProbLoss ? `Value is positive but the chance of loss exceeds the ${inp.maxProbLoss} % tolerance: mature the estimate to Class ${Math.max(1, cls - 1)} or fix the price before sanction.`
            : `Carry these figures into the finance and risk modules. With the quotations the expected CAPEX is ${fM(post.mean)} (P90 ${fM(post.p90)}); use that, not the estimate, as the budget basis.`,
        `Whole-life costing, ISO 15686-5:2017; estimate classes, AACE RP 18R-97; LCOE convention of IEA/NEA (2020); cover-ratio covenants, Yescombe (2014) ch. 12. ${PCT_VALUE}`),
      adv('engineering', N * inp.dff > L0 ? 'fail' : 'info',
        `"${top.name}" is the technical driver: ± 10 % moves NPV between ${fM(Math.min(top.low, top.high))} and ${fM(Math.max(top.low, top.high))}. Fatigue: calculated life ${fmt(L0)} years against ${fmt(N * inp.dff)} required; the section could be scaled to ${fmt(sReq, 4)} on the criterion and the cost-optimal factor is ${fmt(sOpt, 4)} (${fatigueGoverns ? 'the fatigue criterion governs' : 'the expected failure cost, not the criterion, governs'}).`,
        N * inp.dff > L0 ? `Raise the fatigue life before relying on ${N} years of revenue: the section factor needed is ${fmt(sReq, 4)}, about ${fmt(steelS * (sReq - 1))} t more primary steel, or improve the detail class.`
          : `${steelSaved > 0 ? `Test a section factor of ${fmt(sOpt, 4)} in the structural and fatigue analyses: about ${fmt(steelSaved)} t less primary steel and ${fM(npvOpt - npvAsIs)} more NPV, if ultimate strength, buckling and stiffness (natural frequency) still pass — those limit states are not in this model.` : `Keep or increase the section: the optimum is ${fmt(-steelSaved)} t above the present design because the expected cost of fatigue failure outweighs the steel.`} Then secure "${top.name}" at its source suite.`,
        'S-N fatigue and design fatigue factors, DNV-RP-C203 (2021) and DNV-ST-0126 (2021); risk-based structural optimisation (minimum expected life-cycle cost), ISO 2394:2015 §4.4 and Annex G.',
        `Life-cycle emissions are ${fmt(gwp.total_t / 1e3)} kt CO₂e, ${fmt(ci.gPerKWh)} g/kWh, ${tx.pass ? 'inside' : 'above'} the ${inp.taxonomy} g/kWh threshold (${pct(tx.utilisation)} of it)${paysBack ? `, repaid in ${fmt(av.carbonPaybackYears)} years of generation` : ''}. ${steelSaved > 0 ? `The optimal section saves about ${fmt(co2Saved)} t CO₂e of steel` : `The optimal section adds about ${fmt(-co2Saved)} t CO₂e of steel`} at ${inp.efSteel} t CO₂e/t; priced at ${inp.carbonPrice} ${CUR}/t the whole life cycle costs ${fM(pvCarbon)} in present value, leaving an NPV of ${fM(npvCarbon)}. Recycling ${inp.recovery} % of the steel returns ${fmt(Math.max(eol.credit_t, 0))} t CO₂e (${SCREEN} are the defaults; all factors here are inputs).`),
    ];
    if (policy && Number.isFinite(policy.interval)) advice.push(adv('economic', 'info', `Planned replacement of the wearing component at ${fmt(policy.interval)} years costs ${fM(maintAnnual)}/yr for ${inp.units} units against ${fM(maintRtf)}/yr run to failure.`,
      `Plan the replacement campaign at that age and budget ${fM(maintAnnual)}/yr; the saving of ${fM(maintRtf - maintAnnual)}/yr is worth ${fM((maintRtf - maintAnnual) * (1 - (1 + rReal) ** -N) / (rReal || 1e-12))} over the life before tax.`, 'Age-replacement policy, Barlow & Proschan (1965); maintenance within whole-life cost, ISO 15686-5:2017 §5.'));
    if (Math.abs(sur.x[0] - sOpt) > 0.02) advice.push(adv('engineering', 'warn', `The surrogate optimiser stopped at a section factor of ${fmt(sur.x[0], 4)}, ${fmt(Math.abs(sur.x[0] - sOpt), 3)} from the exact optimum.`, 'Allow more adaptive evaluations before using the surrogate method on a model that is expensive to run.', 'Efficient global optimisation, Jones, Schonlau & Welch (1998).'));

    return { kpis, checks, plots, tables, advice,
      notes: [PCT_VALUE, `Real terms: constant ${inp.baseYear} money at ${inp.discountRate} % real (${pct(rNom)} nominal with ${inp.inflation} % inflation); tax depreciation and debt service are computed in money of the day and deflated.`,
        `Chain: ${fmt(inp.steel)} t steel, ${fmt(base.gross)} gross vessel days (${fmt(inp.vesselDays)} net at ${fmt(inp.operability)} operability), ${fmt(inp.ratedMW)} MW at ${pct(netCF)} net capacity factor over ${N} years.`,
        `${n} Latin-hypercube samples, seed ${inp.seed}; price, yield, capital and operating cost sampled independently, capital cost lognormal over the Class ${cls} band read as cumulative P10–P90.`, ir.text],
      publish: { 'econ.technoCapex': base.capex, 'econ.technoOpexAnnual': opex1, 'econ.technoNpv': NPV, 'econ.technoLcoe': L.lcoe, ...(minDscr != null ? { 'econ.technoMinDscr': minDscr } : {}), 'econ.carbonIntensity': ci.gPerKWh, ...(paysBack ? { 'econ.carbonPayback': av.carbonPaybackYears } : {}),
        'econ.lifecycleCO2': gwp.total_t, 'econ.optimalSectionFactor': sOpt, 'econ.steelSaving': steelSaved },
      quantities: { steel_t: inp.steel, concrete_m3: inp.concrete, rock_t: inp.rock, cable_km: inp.cable, pipe_km: inp.pipe, vessel_days: { installation: round(base.gross) } } };
  },
  verify: [
    // 1000 t × 3000 = 3.0 M capex in year 0; 10 MW × 8760 h × 0.5 = 43 800 MWh/yr × 100 = 4.38 M/yr for 20 years, no opex, tax or decommissioning:
    // NPV = −3 + 4.38 × (P/A, 8 %, 20) = −3 + 4.38 × 9.818147 = 40.00348 M; LCOE = 3e6 / (43 800 × 9.818147) = 6.97618 per MWh
    { name: 'Chain from quantities to NPV and LCOE', source: 'Hand calculation with (P/A, 8 %, 20) = 9.818147: capex 3.0 M, energy 43.8 GWh/yr, NPV 40.00348 M, LCOE 6.97618',
      inputs: { steel: 1000, rSteel: 3000, concrete: 0, rock: 0, cable: 0, pipe: 0, vesselDays: 0, equipPerMW: 0, pIndirect: 0, pContingency: 0, ratedMW: 10, capacityFactor: 0.5, availability: 100, degradation: 0, designLife: 20, tenor: 10, fatigueLife: 1e4, failCost: 0,
        constructionYears: 1, opexPerMW: 0, opexVar: 0, pDecom: 0, price: 100, taxRate: 0, discountRate: 8, inflation: 0, units: 0, quotes: [], budgetCap: 1e9 },
      expect: [{ kpi: 'capex', value: 3.0, tol: 1e-9 }, { kpi: 'aep', value: 43.8, tol: 1e-9 }, { kpi: 'npv', value: 40.00348, tol: 1e-6 }, { kpi: 'lcoe', value: 6.97618, tol: 1e-6 }] },
    // same asset: 1000 t × 1.9 = 1900 t CO2e and nothing else; 43 800 MWh × 20 yr = 876 GWh → 1900e6 g / 876e6 kWh = 2.16895 g/kWh;
    // avoided 43 800 MWh × 0.4 t/MWh = 17 520 t/yr → payback 1900 / 17 520 = 0.108447 yr
    { name: 'Life-cycle carbon intensity and payback', source: 'Hand calculation: 1900 t CO₂e over 876 GWh = 2.16895 g/kWh; 1900 / (43 800 × 0.4) = 0.108447 yr',
      inputs: { steel: 1000, rSteel: 3000, concrete: 0, rock: 0, cable: 0, pipe: 0, vesselDays: 0, equipPerMW: 0, pIndirect: 0, pContingency: 0, ratedMW: 10, capacityFactor: 0.5, availability: 100, degradation: 0, designLife: 20, tenor: 10, fatigueLife: 1e4, failCost: 0,
        constructionYears: 1, opexPerMW: 0, opexVar: 0, pDecom: 0, price: 100, taxRate: 0, discountRate: 8, inflation: 0, units: 0, quotes: [], budgetCap: 1e9, efSteel: 1.9, efEquip: 0, opEmissions: 0, recovery: 0, recycled: 0, gridIntensity: 400 },
      expect: [{ kpi: 'co2', value: 1.9, tol: 1e-9 }, { kpi: 'intensity', value: 2.16895, tol: 1e-5 }, { kpi: 'carbonPayback', value: 0.108447, tol: 1e-5 }] },
    // life ∝ s^m: s = (DFF × N / L0)^(1/m) = (3 × 25 / 100)^(1/3) = 0.9085603; with no failure cost less steel is always cheaper, so the optimum is on the limit
    { name: 'Fatigue-limited optimum of the section', source: 'Closed form from the S-N relation N ∝ S^(−m), S ∝ 1/s: s = (3 × 25 / 100)^(1/3) = 0.9085603',
      inputs: { fatigueLife: 100, dff: 3, designLife: 25, snSlope: 3, failCost: 0 }, expect: [{ kpi: 'sReq', value: 0.9085603, tol: 1e-6 }, { kpi: 'sOpt', value: 0.9085603, tol: 1e-4 }] },
    // prior from Class 4 (−30 % / +50 % read as P10–P90) on 3.0 M: mean 3.3, std 2.4 / (2 × 1.2815516) = 0.936365; one quotation at 3.0 M with std 0.3 M:
    // precision 1/0.876779 + 1/0.09 = 1.140539 + 11.111111 = 12.251650; mean = (3.3 × 1.140539 + 3.0 × 11.111111) / 12.251650 = 3.027928
    { name: 'Bayesian update of the estimate with one quotation', source: 'Conjugate normal model by hand (Gelman et al., BDA3 §2.5): posterior mean 3.027928 M from a Class 4 prior N(3.3, 0.936365²) and one quotation of 3.0 M ± 0.3 M',
      inputs: { steel: 1000, rSteel: 3000, concrete: 0, rock: 0, cable: 0, pipe: 0, vesselDays: 0, equipPerMW: 0, pIndirect: 0, pContingency: 0, aaceClass: '4', quotes: [{ name: 'q', value: 100 }], quoteStd: 10, units: 0, tenor: 10, budgetCap: 1e9 },
      expect: [{ kpi: 'capexPosterior', value: 3.027928, tol: 1e-5 }] },
  ],
};
