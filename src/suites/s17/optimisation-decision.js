// Optimisation and decision models: economic comparison and multi-criteria ranking of candidate designs,
// decision under scenario uncertainty, project selection under a budget, and a coupled design–cost optimisation
// (economic conductor sizing) with constrained, multi-objective, chance-constrained and stochastic variants.
import { constrainedMinimise, multiObjective, portfolioMilp, capitalRationing, knapsack, minMaxRegret, chanceConstrained, twoStageStochastic, surrogateOptimise } from '../../engine/econ/optimise.js';
import { topsis, weightedSum, paretoFilter } from '../../engine/econ/decision.js';
import { incrementalAnalysis, compareUnequalLives } from '../../engine/econ/cashflow.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, need, run, round, fmt, mil, fM, pct, adv, check, checkMin, cost, tableRows, uniqueNames, linspace, range, crf } from './_common.js';

const ALTS = [
  { name: 'Monopile', capex: 1.46e9, opex: 4.6e7, energy: 2150, life: 25, co2: 131000, availability: 96.2 },
  { name: 'Jacket', capex: 1.58e9, opex: 4.4e7, energy: 2165, life: 30, co2: 118000, availability: 96.6 },
  { name: 'Gravity base', capex: 1.52e9, opex: 4.3e7, energy: 2150, life: 35, co2: 142000, availability: 96.4 },
  { name: 'Suction bucket jacket', capex: 1.55e9, opex: 4.4e7, energy: 2160, life: 30, co2: 109000, availability: 96.5 },
  { name: 'Floating semi-submersible', capex: 2.05e9, opex: 5.6e7, energy: 2290, life: 25, co2: 176000, availability: 94.8 },
];
const PROJ = [
  { name: 'Scour protection upgrade', capex: 1.8e7, npv: 9.0e6, group: '' }, { name: 'Second export circuit', capex: 9.5e7, npv: 4.1e7, group: 'export' }, { name: 'Dynamic rating of export cable', capex: 1.2e7, npv: 1.9e7, group: 'export' },
  { name: 'Condition monitoring retrofit', capex: 8.0e6, npv: 1.3e7, group: '' }, { name: 'Service operation vessel purchase', capex: 6.0e7, npv: 2.2e7, group: '' }, { name: 'Battery storage 100 MWh', capex: 5.5e7, npv: 1.5e7, group: '' }, { name: 'Blade leading-edge protection', capex: 6.0e6, npv: 7.0e6, group: '' },
];
const METHODS = [{ value: 'de', label: 'Differential evolution' }, { value: 'pso', label: 'Particle swarm' }, { value: 'sa', label: 'Simulated annealing' }, { value: 'ga', label: 'Genetic algorithm' }, { value: 'nm', label: 'Nelder–Mead (local)' }];

export default {
  id: 's17-optimisation-decision', suite: 17, title: 'Design selection, portfolio and cost optimisation',
  summary: 'Ranks candidate engineering designs on value, cost, carbon and performance, tests the choice against price scenarios, selects the best set of projects within a budget, and optimises a design variable against lifecycle cost under an engineering constraint.',
  covers: [9], fidelity: 'numerical',
  standards: ['ISO 15686-5:2017 §5 (comparison of options on whole-life cost)', 'Hwang & Yoon (1981) TOPSIS; Deb et al. (2002) NSGA-II', 'Blank & Tarquin, Engineering Economy, ch. 8 (incremental analysis) and ch. 12 (capital budgeting)', 'IEC 60287-3-2:2012 (economic optimisation of power cable size)', 'Savage (1951) minimax regret; Birge & Louveaux (2011), Introduction to Stochastic Programming'],
  models: ['Linear programming (two-stage recourse in extensive form)', 'Mixed-integer linear programming (portfolio selection)', 'Nonlinear programming with constraints (augmented Lagrangian)', 'Dynamic programming (0/1 knapsack)', 'Stochastic programming (two-stage with recourse, EVPI)', 'Robust decision (minimax regret, maximin, Hurwicz)', 'Chance-constrained optimisation', 'Multi-objective optimisation and Pareto-front analysis (NSGA-II)', 'Genetic algorithm, differential evolution, particle swarm, simulated annealing', 'Bayesian (surrogate-assisted) optimisation', 'Multi-criteria ranking (TOPSIS, weighted sum)', 'Incremental and equivalent-annual-worth comparison'],
  limits: `Candidate designs are compared on the figures supplied for them: the module does not re-engineer the candidates. The conductor-sizing problem is single-variable with a current-density limit standing in for the thermal rating from the cable suites. Mixed-integer nonlinear problems are handled by enumeration of the supplied candidates, not by a general MINLP solver. ${CUR_NOTE}`,
  publishes: ['econ.bestAlternativeNpv', 'econ.portfolioNpv', 'econ.optimalConductor'],
  searchTerms: 'multi-criteria decision foundation concept selection economic cable sizing optimisation',
  inputs: [
    { key: 'alternatives', label: 'Candidate designs', type: 'table', group: 'Alternative design options', role: 'cost', help: `One row per candidate with its cost, yield, life, embodied carbon and availability from the engineering suites. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Design', type: 'text' }, { key: 'capex', label: 'CAPEX', unit: CUR, type: 'number' }, { key: 'opex', label: 'OPEX', unit: `${CUR}/yr`, type: 'number' }, { key: 'energy', label: 'Annual output', unit: 'GWh/yr', type: 'number' }, { key: 'life', label: 'Design life', unit: 'yr', type: 'number' }, { key: 'co2', label: 'Embodied carbon', unit: 't CO₂e', type: 'number' }, { key: 'availability', label: 'Availability', unit: '%', type: 'number' }], default: ALTS },
    cost({ key: 'price', label: 'Base price', unit: `${CUR}/MWh`, default: 96, min: 0.01, max: 1e5, group: 'Economic scenarios' }),
    cost({ key: 'priceLow', label: 'Low price scenario', unit: `${CUR}/MWh`, default: 78, min: 0.01, max: 1e5, group: 'Economic scenarios' }),
    cost({ key: 'priceHigh', label: 'High price scenario', unit: `${CUR}/MWh`, default: 115, min: 0.01, max: 1e5, group: 'Economic scenarios' }),
    cost({ key: 'carbonPrice', label: 'Carbon price on embodied emissions', unit: `${CUR}/t CO₂e`, default: 90, min: 0, max: 2000, group: 'Economic scenarios' }),
    { key: 'rate', label: 'Real discount rate (MARR)', unit: '%/yr', type: 'number', default: 6, min: 0, max: 40, group: 'Financial constraints', role: 'bc', help: PLACEHOLDER },
    { key: 'budget', label: 'Capital budget for the design', unit: CUR, type: 'number', default: 1.8e9, min: 1, max: 1e12, group: 'Financial constraints', role: 'bc', help: 'Candidates above this are infeasible.' },
    { key: 'minAvail', label: 'Minimum availability', unit: '%', type: 'number', default: 95, min: 0, max: 100, group: 'Financial constraints', role: 'bc', help: 'Minimum performance requirement; candidates below it are infeasible.' },
    { key: 'wValue', label: 'Weight: value (NPV)', unit: '-', type: 'number', default: 5, min: 0, max: 10, group: 'Decision weights', role: 'criteria' },
    { key: 'wCapex', label: 'Weight: capital exposure', unit: '-', type: 'number', default: 2, min: 0, max: 10, group: 'Decision weights', role: 'criteria' },
    { key: 'wCarbon', label: 'Weight: embodied carbon', unit: '-', type: 'number', default: 2, min: 0, max: 10, group: 'Decision weights', role: 'criteria' },
    { key: 'wPerf', label: 'Weight: technical performance (availability)', unit: '-', type: 'number', default: 1, min: 0, max: 10, group: 'Decision weights', role: 'criteria' },

    { key: 'projects', label: 'Candidate projects', type: 'table', group: 'Project portfolio', role: 'cost', help: `Independent improvement projects with their capital and NPV. Rows sharing a group name are mutually exclusive. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Project', type: 'text' }, { key: 'capex', label: 'Capital', unit: CUR, type: 'number' }, { key: 'npv', label: 'NPV', unit: CUR, type: 'number' }, { key: 'group', label: 'Exclusive group', type: 'text' }], default: PROJ },
    { key: 'portfolioBudget', label: 'Portfolio capital budget', unit: CUR, type: 'number', default: 1.2e8, min: 0, max: 1e12, group: 'Project portfolio', role: 'bc' },
    { key: 'committed', label: 'Capital already committed from that budget', unit: CUR, type: 'number', default: 8.0e6, min: 0, max: 1e12, group: 'Project portfolio', role: 'ic',
      help: 'Opening financial state: commitments already made against the portfolio budget. They are sunk for the selection and only reduce the capital still available.' },

    { key: 'current', label: 'Maximum load current per phase', unit: 'A', type: 'number', default: 800, min: 1, max: 10000, group: 'Design optimisation (conductor size)', role: 'load', uncertain: { dist: 'normal', cov: 0.08 } },
    { key: 'length', label: 'Cable route length', unit: 'km of cable', type: 'number', default: 60, min: 0.01, max: 2000, group: 'Design optimisation (conductor size)', role: 'geometry', bind: 'quantities.cable_km' },
    { key: 'resistivity', label: 'Conductor AC resistivity at operating temperature', unit: 'Ω·mm²/m', type: 'number', default: 0.0225, min: 0.01, max: 0.08, group: 'Design optimisation (conductor size)', role: 'material', help: 'Copper ≈ 0.0225, aluminium ≈ 0.037 at 90 °C including skin and proximity effects.' },
    { key: 'lossLoadFactor', label: 'Loss load factor', unit: '-', type: 'number', default: 0.35, min: 0.01, max: 1, group: 'Design optimisation (conductor size)', role: 'operation', help: 'Mean of (I/Imax)² over the year; about 0.3–0.4 for wind export.' },
    { key: 'jMax', label: 'Permissible current density (thermal limit)', unit: 'A/mm²', type: 'number', default: 1.7, min: 0.2, max: 6, group: 'Design optimisation (conductor size)', role: 'bc', help: 'Engineering performance limit standing in for the ampacity from the cable thermal analysis.' },
    cost({ key: 'cPerArea', label: 'Cable cost per conductor area', unit: `${CUR}/m/mm²`, default: 0.55, min: 0.001, max: 100, group: 'Design optimisation (conductor size)', help: 'Incremental supply cost of three cores per metre for each extra mm² of conductor.' }),
    cost({ key: 'cFixed', label: 'Area-independent cable and lay cost', unit: `${CUR}/m`, default: 520, min: 0, max: 1e5, group: 'Design optimisation (conductor size)' }),
    cost({ key: 'lossPrice', label: 'Value of lost energy', unit: `${CUR}/MWh`, default: 90, min: 0, max: 1e5, group: 'Design optimisation (conductor size)' }),
    { key: 'cableLife', label: 'Cable service life', unit: 'yr', type: 'int', default: 25, min: 1, max: 60, group: 'Design optimisation (conductor size)', role: 'bc', bind: 'life.designYears' },
    cost({ key: 'efConductor', label: 'Embodied carbon of conductor metal', unit: 't CO₂e/t', default: 4.0, min: 0, max: 50, group: 'Design optimisation (conductor size)', help: 'Screening emission factor for copper (aluminium is higher per tonne, lower per ampere).' }),
    cost({ key: 'gridIntensity', label: 'Carbon intensity of lost energy', unit: 'g CO₂/kWh', default: 350, min: 0, max: 1200, group: 'Design optimisation (conductor size)' }),
    { key: 'method', label: 'Global optimiser', type: 'select', default: 'de', options: METHODS, group: 'Analysis control', role: 'control' },
    { key: 'loadConfidence', label: 'Required probability of staying within the thermal limit', unit: '%', type: 'number', default: 95, min: 50, max: 99.9, group: 'Design optimisation (conductor size)', role: 'criteria' },
    { key: 'loadCov', label: 'Uncertainty of the maximum load (CoV)', unit: '%', type: 'number', default: 8, min: 0.1, max: 50, group: 'Design optimisation (conductor size)', role: 'load' },

    { key: 'daysExpected', label: 'Expected offshore vessel days needed', unit: 'vessel-days', type: 'number', default: 210, min: 1, max: 5000, group: 'Vessel booking under uncertainty', role: 'operation', bind: 'install.vesselDays' },
    { key: 'daysSpread', label: 'Low/high spread of vessel days', unit: '± %', type: 'number', default: 25, min: 0, max: 90, group: 'Vessel booking under uncertainty', role: 'operation' },
    cost({ key: 'rateReserved', label: 'Day rate if reserved in advance', unit: `${CUR}/day`, default: 2.6e5, min: 1, max: 1e8, group: 'Vessel booking under uncertainty' }),
    cost({ key: 'rateSpot', label: 'Day rate for extra days at short notice', unit: `${CUR}/day`, default: 4.2e5, min: 1, max: 1e8, group: 'Vessel booking under uncertainty' }),
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 17, min: 1, max: 1e6, group: 'Analysis control', role: 'control' },
  ],
  compute(inp) {
    const alts = tableRows(inp.alternatives, 'Candidate designs', { capex: { min: 1 }, opex: { min: 0 }, energy: { min: 1e-9 }, life: { min: 1, max: 80 }, co2: { min: 0 }, availability: { min: 0, max: 100 } }, { min: 2 });
    uniqueNames(alts, 'Candidate designs');
    const projs = tableRows(inp.projects, 'Candidate projects', { capex: { min: 0 }, npv: {} }, { min: 1 });
    uniqueNames(projs, 'Candidate projects');
    need(inp.priceLow <= inp.price && inp.price <= inp.priceHigh, 'Price scenarios must satisfy low ≤ base ≤ high.');
    need(inp.wValue + inp.wCapex + inp.wCarbon + inp.wPerf > 0, 'At least one decision weight must be positive.');
    need(inp.rateSpot >= inp.rateReserved, 'The short-notice day rate must be at least the reserved rate; otherwise there is nothing to decide.');
    const r = inp.rate / 100;
    const pa = (n) => 1 / crf(r, n);

    /* ---- (A) candidate designs ---- */
    const flows = (a, price) => [-a.capex - a.co2 * inp.carbonPrice, ...range(Math.round(a.life), () => a.energy * 1e3 * price - a.opex)];
    const ev = alts.map((a) => {
      const n = Math.round(a.life);
      const npvAt = (p) => -a.capex - a.co2 * inp.carbonPrice + (a.energy * 1e3 * p - a.opex) * pa(n);
      const npv = npvAt(inp.price);
      return { ...a, n, npv, npvLow: npvAt(inp.priceLow), npvHigh: npvAt(inp.priceHigh), eaw: npv * crf(r, n), lcoe: (a.capex * crf(r, n) + a.opex) / (a.energy * 1e3), intensity: a.co2 * 1e6 / (a.energy * 1e6 * n),
        feasible: a.capex <= inp.budget && a.availability >= inp.minAvail };
    });
    const feas = ev.filter((a) => a.feasible);
    need(feas.length >= 1, `No candidate design satisfies both the capital budget (${fM(inp.budget)}) and the minimum availability (${inp.minAvail} %). Relax a constraint or add candidates.`);
    const crit = [{ name: 'Annual worth', weight: inp.wValue, sense: 'max' }, { name: 'CAPEX', weight: inp.wCapex, sense: 'min' }, { name: 'Embodied carbon', weight: inp.wCarbon, sense: 'min' }, { name: 'Availability', weight: inp.wPerf, sense: 'max' }].filter((c) => c.weight > 0);
    const col = { 'Annual worth': 'eaw', CAPEX: 'capex', 'Embodied carbon': 'co2', Availability: 'availability' };
    const matrix = feas.map((a) => crit.map((c) => a[col[c.name]]));
    const tp = feas.length > 1 ? run(() => topsis({ alternatives: feas.map((a) => a.name), criteria: crit, matrix }), 'TOPSIS') : { best: feas[0].name, closeness: [1], ranking: [{ name: feas[0].name, closeness: 1, rank: 1 }] };
    const ws = feas.length > 1 ? run(() => weightedSum({ alternatives: feas.map((a) => a.name), criteria: crit, matrix }), 'Weighted sum') : { best: feas[0].name, scores: [1] };
    const par = paretoFilter(feas.map((a) => [a.eaw, a.capex, a.co2]), ['max', 'min', 'min']);
    // unequal lives: equivalent annual worth under like-for-like renewal; equal lives: incremental NPV analysis
    const sameLife = feas.every((a) => a.n === feas[0].n);
    const ul = run(() => compareUnequalLives({ alternatives: feas.map((a) => ({ name: a.name, cashflows: flows(a, inp.price) })), rate: r }), 'Annual-worth comparison');
    const inc = sameLife && feas.length > 1 ? run(() => incrementalAnalysis({ alternatives: feas.map((a) => ({ name: a.name, cashflows: flows(a, inp.price) })), marr: r }), 'Incremental analysis') : null;
    const byEaw = feas.slice().sort((a, b) => b.eaw - a.eaw);
    const best = byEaw[0];
    const reg = feas.length > 1 ? minMaxRegret({ alternatives: feas.map((a) => a.name), scenarios: ['Low price', 'Base price', 'High price'], payoffs: feas.map((a) => [a.npvLow * crf(r, a.n), a.eaw, a.npvHigh * crf(r, a.n)]), sense: 'max', probabilities: [0.25, 0.5, 0.25] }) : null;
    const chosen = feas.find((a) => a.name === tp.best);

    /* ---- (B) portfolio under one budget ---- */
    const groups = {}; projs.forEach((p, i) => { const g = String(p.group ?? '').trim(); if (g) (groups[g] = groups[g] || []).push(i); });
    const exclusive = Object.values(groups).filter((g) => g.length > 1);
    const pos = projs.map((p) => Math.max(p.npv, 0));
    // opening commitments are sunk: they reduce the capital still available but do not enter any project's value
    need(inp.committed <= inp.portfolioBudget, `Capital already committed (${fM(inp.committed)}) exceeds the portfolio budget (${fM(inp.portfolioBudget)}): no capital is left to allocate.`);
    const avail = inp.portfolioBudget - inp.committed;
    const milp = run(() => portfolioMilp({ values: projs.map((p) => p.npv), budgets: [{ use: projs.map((p) => p.capex), limit: avail }], exclusive }), 'Portfolio selection');
    const ration = run(() => capitalRationing({ projects: projs.map((p) => ({ name: p.name, cost: Math.max(p.capex, 1e-9), npv: p.npv })), budget: avail }), 'Capital rationing');
    const unit = Math.max(1, 10 ** Math.floor(Math.log10(Math.max(avail, 1) / 2000)));
    const kn = run(() => knapsack({ values: pos, weights: projs.map((p) => p.capex), capacity: avail, resolution: unit }), 'Knapsack');
    const sel = new Set(milp.selected);
    const used = milp.selected.reduce((s, i) => s + projs[i].capex, 0);

    /* ---- (C) economic conductor size ---- */
    const L = inp.length * 1e3; const n = Math.round(inp.cableLife);
    // PV of joule losses per unit of 1/A: 3 I² ρ L [W·mm²] × h × LLF × price/1e6 × (P/A)
    const lossMWhPerInvA = 3 * inp.current ** 2 * inp.resistivity * L * 8760 * inp.lossLoadFactor / 1e6;
    const K = lossMWhPerInvA * inp.lossPrice * pa(n);
    const lcc = (A) => (inp.cFixed + inp.cPerArea * A) * L + K / A;
    const aMin = inp.current / inp.jMax;
    const aKelvin = Math.sqrt(K / (inp.cPerArea * L)); // unconstrained closed-form optimum (Kelvin's law)
    const hi = Math.max(4 * aKelvin, 2 * aMin, 100);
    const opt = run(() => constrainedMinimise({ objective: (x) => lcc(x[0]) / 1e6, bounds: [[Math.min(aMin, aKelvin) * 0.25, hi]], ineq: [(x) => (aMin - x[0]) / aMin], method: inp.method, seed: inp.seed, maxGen: 60, outerIter: 25 }), 'Conductor optimisation');
    // The multiplier method approaches an active limit from the infeasible side and the engine flags any residual above
    // 1e-6 as infeasible. A residual of 1e-4 of the thermal minimum is far below a conductor size step, so it is
    // accepted and the size is taken on the limit itself.
    need(Number.isFinite(opt.x[0]) && opt.maxViolation <= 1e-4, 'The conductor optimisation did not converge to a feasible size; try another optimiser.');
    const aOpt = Math.max(opt.x[0], aMin);
    const sur = run(() => surrogateOptimise({ f: (x) => lcc(Math.max(x[0], aMin)) / 1e6, bounds: [[aMin, hi]], nInit: 6, nIter: 12, seed: inp.seed }), 'Surrogate optimisation');
    // carbon of a size: conductor metal (8.9 t/m³ copper-equivalent, three cores) + emissions of the lifetime losses
    const co2Of = (A) => 3 * A * 1e-6 * L * 8.9 * inp.efConductor + lossMWhPerInvA / A * n * inp.gridIntensity / 1e3;
    const mo = run(() => multiObjective({ objectives: (x) => [lcc(x[0]) / 1e6, co2Of(x[0]) / 1e3], bounds: [[aMin, hi]], names: ['Lifecycle cost', 'Carbon'], pop: 32, generations: 30, seed: inp.seed }), 'Multi-objective optimisation');
    const front = mo.front.slice().sort((a, b) => a.x[0] - b.x[0]);
    const aCarbon = Math.sqrt(lossMWhPerInvA * n * inp.gridIntensity / 1e3 / (3 * 1e-6 * L * 8.9 * inp.efConductor || 1e-12));
    // chance constraint: P[load ≤ jMax × A] ≥ p with a normally distributed peak load
    const cc = run(() => chanceConstrained({ objective: (x) => lcc(x[0]) / 1e6, constraint: (x, xi) => (xi.load - inp.jMax * x[0]) / inp.current, vars: { load: { dist: 'normal', mean: inp.current, std: inp.current * inp.loadCov / 100 } }, bounds: [[aMin * 0.5, hi]], p: inp.loadConfidence / 100, n: 200, seed: inp.seed, maxGen: 40 }), 'Chance-constrained optimisation');

    /* ---- (D) two-stage stochastic booking of vessel days ---- */
    const dLo = inp.daysExpected * (1 - inp.daysSpread / 100); const dHi = inp.daysExpected * (1 + inp.daysSpread / 100);
    const sc = [[0.25, dLo], [0.5, inp.daysExpected], [0.25, dHi]];
    // minimise c x + Σ p q y_s subject to x + y_s ≥ d_s  (written −x − y_s ≤ −d_s)
    const ts = run(() => twoStageStochastic({ c: [inp.rateReserved], scenarios: sc.map(([p, d]) => ({ p, q: [inp.rateSpot], T: [[-1]], W: [[-1]], h: [-d] })) }), 'Stochastic programme');
    const detCost = inp.rateReserved * inp.daysExpected + sc.reduce((s, [p, d]) => s + p * inp.rateSpot * Math.max(d - inp.daysExpected, 0), 0);

    const kpis = [
      { key: 'bestNpv', label: `NPV of the preferred design (${chosen.name})`, value: mil(chosen.npv), unit: MU, explain: 'Net present value at the base price, including the carbon cost of embodied emissions, of the design ranked first by the weighted multi-criteria method.' },
      { key: 'bestEaw', label: 'Highest equivalent annual worth among feasible designs', value: mil(best.eaw), unit: `${MU}/yr`, explain: `Annual worth of "${best.name}". With different design lives, annual worth (not NPV) is the like-for-like economic measure.` },
      { key: 'topsis', label: 'TOPSIS closeness of the preferred design', value: round(Math.max(...tp.closeness)), unit: '-', explain: 'Relative closeness to the ideal solution on all weighted criteria (1 = best on everything, 0 = worst on everything). A clear gap to the runner-up means the choice is robust to the weights.' },
      { key: 'nPareto', label: 'Designs on the Pareto front', value: par.indices.length, unit: '-', explain: 'Feasible designs not beaten on value, capital and carbon together. Only these deserve further study; the rest are dominated.' },
      { key: 'nFeasible', label: 'Feasible designs', value: feas.length, unit: '-', explain: 'Candidates meeting both the capital budget and the minimum availability.' },
      { key: 'lcoeBest', label: 'Lowest LCOE among feasible designs', value: round(Math.min(...feas.map((a) => a.lcoe))), unit: `${CUR}/MWh`, explain: 'Capital recovery plus OPEX per MWh at the real discount rate; a value-independent cost measure for the comparison.' },
      { key: 'portfolioNpv', label: 'NPV of the optimal project portfolio', value: mil(milp.value), unit: MU, explain: 'Largest total NPV obtainable within the portfolio budget, respecting mutually exclusive groups (mixed-integer programme).' },
      { key: 'portfolioUsed', label: 'Budget used by the portfolio', value: mil(used), unit: MU, limit: mil(avail), sense: 'max', explain: 'Capital committed by the selected projects, against the budget left after opening commitments; the unused remainder is too small for any further positive-NPV project.' },
      { key: 'greedyNpv', label: 'NPV by profitability-index ranking', value: mil(ration.greedy.npv), unit: MU, explain: 'Total NPV if projects are simply funded in order of NPV per unit of capital and exclusivity is ignored. It can fall short of the optimum when projects are lumpy.' },
      { key: 'optA', label: 'Cost-optimal conductor cross-section', value: round(aOpt), unit: 'mm²', explain: 'Cross-section minimising cable cost plus the present value of joule losses over the service life, subject to the thermal limit. Choose the next standard size above it.' },
      { key: 'kelvinA', label: 'Unconstrained economic cross-section (Kelvin)', value: round(aKelvin), unit: 'mm²', explain: 'Closed-form optimum where the marginal cost of conductor equals the marginal saving in losses. When it is below the thermal minimum, the thermal limit governs.' },
      { key: 'minA', label: 'Thermal minimum cross-section', value: round(aMin), unit: 'mm²', explain: 'Smallest cross-section that carries the maximum load at the permissible current density.' },
      { key: 'lccOpt', label: 'Lifecycle cost at the optimum', value: mil(lcc(aOpt)), unit: MU, explain: `Cable supply and lay cost plus discounted losses. Sizing at the thermal minimum instead would cost ${fM(lcc(aMin))}.` },
      { key: 'ccA', label: `Cross-section for ${inp.loadConfidence} % thermal confidence`, value: round(cc.x[0]), unit: 'mm²', explain: 'Cost-optimal size that keeps the uncertain peak load within the thermal limit with the required probability (chance constraint).' },
      { key: 'kneeA', label: 'Knee of the cost–carbon front', value: round(mo.knee.x[0]), unit: 'mm²', explain: 'Compromise size on the Pareto front between lifecycle cost and lifecycle carbon: beyond it, further carbon savings cost disproportionately more.' },
      { key: 'reserveDays', label: 'Vessel days to reserve in advance', value: round(ts.x[0]), unit: 'days', explain: 'First-stage decision of the two-stage stochastic programme: days booked now at the reserved rate, with shortfalls bought at the short-notice rate in each scenario.' },
      { key: 'evpi', label: 'Expected value of perfect information on vessel days', value: mil(ts.evpi ?? 0), unit: MU, explain: 'How much less the campaign would cost on average if the required days were known before booking: the ceiling on what better duration forecasts are worth.' },
    ];
    const checks = [
      check('Preferred design within capital budget', chosen.capex, inp.budget, CUR, 'Financial constraint (budget limit)', 'Guards against selecting a concept that cannot be funded.'),
      checkMin('Preferred design meets minimum availability', chosen.availability, inp.minAvail, '%', 'Minimum performance requirement (ISO 20815:2018 production assurance)', 'Guards against trading away technical performance for cost; below the minimum the revenue forecast behind the NPV is not deliverable.'),
      check('Portfolio within budget', used, Math.max(avail, 1e-9), CUR, 'Capital rationing constraint', 'Guards against over-committing capital across projects.'),
      check('Conductor current density ≤ thermal limit', inp.current / aOpt, inp.jMax, 'A/mm²', 'IEC 60287 thermal rating (current-density proxy)', 'Guards against a conductor that overheats at maximum load, which ages the insulation and shortens cable life.'),
    ];
    const as = linspace(Math.max(aMin * 0.6, 20), hi * 0.8, 70);
    const plots = [
      { type: 'scatter', title: 'Candidate designs: annual worth against embodied carbon', x: { label: 'Embodied carbon', unit: 'kt CO₂e' }, y: { label: 'Equivalent annual worth', unit: `${MU}/yr` },
        series: [{ name: 'Pareto-optimal', x: par.indices.map((i) => round(feas[i].co2 / 1e3)), y: par.indices.map((i) => mil(feas[i].eaw)) }, { name: 'Dominated', x: par.dominated.map((i) => round(feas[i].co2 / 1e3)), y: par.dominated.map((i) => mil(feas[i].eaw)) },
          { name: 'Infeasible', x: ev.filter((a) => !a.feasible).map((a) => round(a.co2 / 1e3)), y: ev.filter((a) => !a.feasible).map((a) => mil(a.eaw)) }],
        explain: 'Each candidate by value and carbon. The best designs are top-left; points labelled dominated are beaten by another candidate on value, capital and carbon at once.' },
      { type: 'bar', title: 'Multi-criteria ranking', categories: feas.map((a) => a.name), x: { label: 'Design' }, y: { label: 'Score', unit: '-' }, series: [{ name: 'TOPSIS closeness', y: tp.closeness.map((v) => round(v, 4)) }, { name: 'Weighted sum', y: ws.scores.map((v) => round(v, 4)) }],
        explain: 'Scores of the feasible designs by two ranking methods. Agreement between them, and a clear leader, indicate a choice that does not hinge on the method.' },
      { type: 'line', title: 'Lifecycle cost against conductor size', x: { label: 'Conductor cross-section', unit: 'mm²' }, y: { label: 'Cost', unit: MU },
        series: [{ name: 'Total', x: as.map((v) => round(v, 4)), y: as.map((A) => mil(lcc(A))) }, { name: 'Cable supply and lay', x: as.map((v) => round(v, 4)), y: as.map((A) => mil((inp.cFixed + inp.cPerArea * A) * L)) }, { name: 'PV of losses', x: as.map((v) => round(v, 4)), y: as.map((A) => mil(K / A)) },
          { name: 'Thermal minimum', x: [aMin, aMin], y: [0, mil(lcc(aMin))], dash: true }, { name: 'Optimum', x: [round(aOpt, 4)], y: [mil(lcc(aOpt))] }],
        explain: 'Cable cost rises and the present value of losses falls with cross-section; the total has a shallow minimum. Sizes left of the dashed line violate the thermal limit.' },
      { type: 'scatter', title: 'Cost–carbon trade-off of the conductor size', x: { label: 'Lifecycle cost', unit: MU }, y: { label: 'Lifecycle carbon', unit: 'kt CO₂e' }, series: [{ name: 'Pareto front', x: front.map((f) => round(f.f[0], 5)), y: front.map((f) => round(f.f[1], 5)) }, { name: 'Knee', x: [round(mo.knee.f[0], 5)], y: [round(mo.knee.f[1], 5)] }],
        explain: 'Non-dominated sizes: moving along the front trades money for carbon. The knee is the last point where carbon is cheap to remove.' },
      { type: 'bar', title: 'Project portfolio', categories: projs.map((p) => p.name), x: { label: 'Project' }, y: { label: 'NPV', unit: MU }, series: [{ name: 'Selected', y: projs.map((p, i) => (sel.has(i) ? mil(p.npv) : 0)) }, { name: 'Not selected', y: projs.map((p, i) => (sel.has(i) ? 0 : mil(p.npv))) }],
        explain: 'NPV of each candidate project, split by whether the optimal portfolio funds it. A high-NPV project left out is either too large for the remaining budget or excluded by a better alternative in its group.' },
    ];
    const tables = [
      { title: 'Economic comparison of candidate designs', columns: [{ key: 'name', label: 'Design' }, { key: 'feas', label: 'Feasible' }, { key: 'capex', label: 'CAPEX', unit: MU }, { key: 'npv', label: 'NPV', unit: MU }, { key: 'eaw', label: 'Annual worth', unit: `${MU}/yr` }, { key: 'lcoe', label: 'LCOE', unit: `${CUR}/MWh` }, { key: 'ci', label: 'Embodied intensity', unit: 'g CO₂e/kWh' }, { key: 'av', label: 'Availability', unit: '%' }, { key: 'par', label: 'Pareto' }, { key: 'tp', label: 'TOPSIS rank' }],
        rows: ev.map((a) => { const fi = feas.indexOf(a); return { name: a.name, feas: a.feasible ? 'yes' : (a.capex > inp.budget ? 'over budget' : 'availability'), capex: mil(a.capex), npv: mil(a.npv), eaw: mil(a.eaw), lcoe: round(a.lcoe, 5), ci: round(a.intensity, 4), av: a.availability, par: fi >= 0 ? (par.indices.includes(fi) ? 'yes' : 'no') : '–', tp: fi >= 0 ? tp.ranking.find((q) => q.name === a.name).rank : '–' }; }),
        explain: 'All candidates with feasibility, economic measures, carbon intensity and rank. Annual worth is the comparison measure when lives differ.' },
      { title: 'Decision under price scenarios', columns: [{ key: 'rule', label: 'Decision rule' }, { key: 'choice', label: 'Chosen design' }],
        rows: reg ? [{ rule: 'Expected value (25/50/25 %)', choice: reg.expectedValue }, { rule: 'Minimax regret (Savage)', choice: reg.minimaxRegret }, { rule: 'Maximin (Wald, pessimist)', choice: reg.maximin }, { rule: 'Maximax (optimist)', choice: reg.maximax }, { rule: 'Hurwicz α = 0.5', choice: reg.hurwicz }, { rule: 'Multi-criteria (TOPSIS)', choice: tp.best }, { rule: 'Highest annual worth', choice: ul.bestByEaw }, ...(inc ? [{ rule: 'Incremental NPV analysis', choice: inc.selected ?? 'do nothing' }] : [])] : [{ rule: 'Only one feasible design', choice: feas[0].name }],
        explain: 'The design each decision rule selects across the low, base and high price scenarios. If all rules agree, the choice is robust; disagreement shows which attitude to risk drives it.' },
      { title: 'Portfolio selection', columns: [{ key: 'name', label: 'Project' }, { key: 'capex', label: 'Capital', unit: MU }, { key: 'npv', label: 'NPV', unit: MU }, { key: 'pi', label: 'NPV per unit capital' }, { key: 'grp', label: 'Exclusive group' }, { key: 'sel', label: 'Optimal portfolio' }],
        rows: projs.map((p, i) => ({ name: p.name, capex: mil(p.capex), npv: mil(p.npv), pi: p.capex > 0 ? round(p.npv / p.capex, 4) : '–', grp: String(p.group ?? '') || '–', sel: sel.has(i) ? 'yes' : 'no' })), explain: 'Each project with its capital, NPV and profitability index, and whether the optimal portfolio includes it.' },
    ];
    const runner = tp.ranking.length > 1 ? tp.ranking[1] : null;
    const kg = inp.gridIntensity; const lossMWh = lossMWhPerInvA / aOpt;
    const advice = [
      adv('engineering', aKelvin < aMin ? 'warn' : 'info', `Economic conductor size ${fmt(aOpt)} mm² (current density ${fmt(inp.current / aOpt)} A/mm²) against a thermal minimum of ${fmt(aMin)} mm²: ${aKelvin < aMin ? 'the thermal limit governs' : `losses, not ampacity, govern — the economic size is ${pct(aOpt / aMin - 1)} larger than the thermal minimum`}. With load uncertainty of ${inp.loadCov} %, ${fmt(cc.x[0])} mm² keeps the cable within its rating with ${inp.loadConfidence} % probability.`,
        `Select the standard size at or above ${fmt(Math.max(aOpt, cc.x[0]))} mm² and confirm its ampacity with the thermal analysis of the cable suite. The cost curve is flat near the optimum: ±10 % on area changes lifecycle cost by only ${pct(lcc(aOpt * 1.1) / lcc(aOpt) - 1, 2)}, so round up.`,
        'Economic optimisation of cable size, IEC 60287-3-2:2012; thermal rating per IEC 60287-1-1:2023.',
        `At the optimum the cable loses about ${fmt(lossMWh)} MWh a year (${fmt(lossMWh * kg / 1e3)} t CO₂ at ${kg} g/kWh) and embodies ${fmt(3 * aOpt * 1e-6 * L * 8.9 * inp.efConductor)} t CO₂e of conductor at ${inp.efConductor} t CO₂e/t (screening factor). The carbon-optimal size would be ${fmt(aCarbon)} mm²; the knee of the cost–carbon front is ${fmt(mo.knee.x[0])} mm².`),
      adv('economic', 'ok', `"${tp.best}" ranks first on the weighted criteria (closeness ${fmt(Math.max(...tp.closeness))}${runner ? `, runner-up "${runner.name}" ${fmt(runner.closeness)}` : ''}); "${best.name}" has the highest annual worth at ${fM(best.eaw)}/yr. ${par.indices.length} of ${feas.length} feasible designs are Pareto-optimal. The optimal portfolio funds ${milp.selected.length} projects for ${fM(milp.value)} of NPV using ${fM(used)} of the ${fM(avail)} still available (${fM(inp.committed)} of the ${fM(inp.portfolioBudget)} budget is already committed).`,
        `${reg && reg.minimaxRegret !== tp.best ? `The minimax-regret choice ("${reg.minimaxRegret}") differs from the multi-criteria choice: decide which price view governs before freezing the concept. ` : ''}${ration.rankingIsOptimal ? 'Ranking by profitability index happens to give the optimal portfolio here.' : `Ranking by profitability index would give ${fM(ration.greedy.npv)}, less than the optimum: use the integer optimum.`} Reserve ${fmt(ts.x[0])} vessel days now (expected cost ${fM(ts.objective)} against ${fM(detCost)} if the expected ${fmt(inp.daysExpected)} days were booked).`,
        'Comparison of alternatives on equivalent annual worth and incremental analysis, Blank & Tarquin ch. 6 and 8; capital rationing, Brealey, Myers & Allen ch. 5; whole-life option appraisal, ISO 15686-5:2017 §5.'),
    ];
    return { kpis, checks, plots, tables, advice,
      notes: [`Real terms at ${inp.rate} %/yr; design NPVs include the carbon cost of embodied emissions at ${inp.carbonPrice} ${CUR}/t.`, `Conductor optimisation by ${METHODS.find((m) => m.value === inp.method).label.toLowerCase()} with an augmented-Lagrangian constraint: ${fmt(aOpt)} mm² (closed-form check max(Kelvin, thermal) = ${fmt(Math.max(aKelvin, aMin))} mm²; surrogate-assisted search ${fmt(Math.max(sur.x[0], aMin))} mm² in ${sur.evaluations} evaluations).`,
        `Knapsack cross-check of the portfolio without exclusivity: ${fM(kn.value)}.`, 'Seeded optimisers: identical inputs reproduce identical results.'],
      publish: { 'econ.bestAlternativeNpv': chosen.npv, 'econ.portfolioNpv': milp.value, 'econ.optimalConductor': aOpt } };
  },
  verify: [
    // K/(c1 L) = 3 ρ h LLF p (P/A)/1e6 / c1 × I² : 3 × 0.0225 × 8760 × 0.35 × 90/1e6 × 12.78336/0.55 = 0.432914 → A* = 800 × 0.657962 = 526.37 mm² (> thermal minimum 800/1.7 = 470.6)
    { name: 'Economic conductor size (Kelvin’s law)', source: 'IEC 60287-3-2 principle, closed form A* = I √(3 ρ h LLF p (P/A, 6 %, 25) / (10⁶ c₁)) = 800 × √0.432914 = 526.37 mm²',
      inputs: { current: 800, resistivity: 0.0225, lossLoadFactor: 0.35, lossPrice: 90, rate: 6, cableLife: 25, cPerArea: 0.55, jMax: 1.7 }, expect: [{ kpi: 'kelvinA', value: 526.37, tol: 1e-4 }, { kpi: 'optA', value: 526.37, tol: 0.01 }, { kpi: 'minA', value: 470.588, tol: 1e-5 }] },
    // thermal limit governs when jMax is low: A = 800/1.2 = 666.67 mm²
    { name: 'Thermal constraint active', source: 'Constrained optimum of a convex cost: A = max(A_Kelvin, I/J_max) = 800/1.2 = 666.67 mm²', inputs: { current: 800, jMax: 1.2 }, expect: [{ kpi: 'optA', value: 666.667, tol: 0.005 }] },
    // 0/1 knapsack: values 60, 100, 120; costs 10, 20, 30; budget 50 → take the last two, value 220 (ranking by value/cost would take the first two: 160)
    { name: 'Portfolio optimum of the textbook knapsack', source: 'Classic 0/1 knapsack instance (values 60/100/120, weights 10/20/30, capacity 50): optimum 220, greedy 160',
      inputs: { projects: [{ name: 'a', capex: 10e6, npv: 60e6, group: '' }, { name: 'b', capex: 20e6, npv: 100e6, group: '' }, { name: 'c', capex: 30e6, npv: 120e6, group: '' }], portfolioBudget: 50e6, committed: 0 }, expect: [{ kpi: 'portfolioNpv', value: 220, tol: 1e-9 }, { kpi: 'greedyNpv', value: 160, tol: 1e-9 }] },
    // newsvendor: reserve up to the critical fractile (q − c)/q = (300 − 200)/300 = 1/3 → the smallest scenario with cumulative probability ≥ 1/3 is the middle one: 100 days
    { name: 'Two-stage booking equals the newsvendor fractile', source: 'Newsvendor solution: reserve the (q − c)/q = 1/3 quantile of demand = 100 days for scenarios 75/100/125 with probabilities 0.25/0.5/0.25',
      inputs: { daysExpected: 100, daysSpread: 25, rateReserved: 2e5, rateSpot: 3e5 }, expect: [{ kpi: 'reserveDays', value: 100, tol: 1e-6 }] },
  ],
};
