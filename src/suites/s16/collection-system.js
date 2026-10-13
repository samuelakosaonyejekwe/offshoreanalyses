// Array collection system of an offshore wind farm: string layout (capacitated minimum spanning tree or simple
// rows), cable section per branch, radial AC load flow at rated power and at any production level, annual energy
// loss from the wind-speed distribution, energy not delivered after a single cable fault (N-1) and the comparison
// of 33 kV with 66 kV collection. Engineering units (MW, kV, km, Ω/km) are converted to SI at the top of compute.
import { genericCatalogue, buildNetwork, esauWilliams, radialStrings, sizeCables, radialLoadFlow, annualEnergyLoss, energyNotDelivered, cableQuantities, compareCollectionVoltages } from '../../engine/electrical/collection.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, HOURS_YEAR, LIMIT_CODE, CONFIRM, CO2_COPPER, energyCo2Text, vesselCo2Text, costInput, money } from './_common.js';

const RHO_CU = 8960;   // kg/m³
// Sections commonly stocked for array cables; rating each one with the IEC 60287 engine is the slow step, so the
// ratings are kept per voltage and frequency (a pure function of both) and only the costs are applied per run.
const AREAS = [150, 240, 400, 630, 800].map((a) => a * 1e-6);
const catCache = new Map();
function catalogueFor(voltage, f, costFn) {
  const key = `${voltage}:${f}`;
  if (!catCache.has(key)) catCache.set(key, genericCatalogue({ voltage, f, areas: AREAS }));
  return catCache.get(key).map((c) => ({ ...c, cost: costFn(c.area) }));
}

/** Production levels and their probabilities from a Weibull wind-speed distribution and a generic power curve (cubic between cut-in and rated). */
function productionLevels({ A, k, nBins, cutIn = 3, ratedU = 11.5, cutOut = 25 }) {
  const F = (u) => 1 - Math.exp(-((u / A) ** k)), frac = (u) => (u < cutIn || u > cutOut ? 0 : Math.min(1, ((u - cutIn) / (ratedU - cutIn)) ** 3));
  const out = [], du = (ratedU - cutIn) / nBins;
  out.push({ fraction: 0, probability: F(cutIn) + 1 - F(cutOut) });
  for (let i = 0; i < nBins; i++) { const u0 = cutIn + i * du; out.push({ fraction: frac(u0 + du / 2), probability: F(u0 + du) - F(u0) }); }
  out.push({ fraction: 1, probability: F(cutOut) - F(ratedU) });
  return out;
}

export default {
  id: 's16-collection-system', suite: 16, title: 'Array collection system (load flow, losses, string sizing, N-1)',
  summary: 'Lays out the array cables in strings, sizes each section, solves the radial load flow at rated and part load, and reports voltage rise, annual losses, cable quantities and the energy lost to a single cable fault.',
  covers: [15], fidelity: 'numerical',
  standards: ['IEC 60287-1-1:2023 (current rating)', 'IEC 60502-2:2014 (cables 6–30 kV: construction)', 'IEC 63026:2019 (submarine cables 6–60 kV)', 'IEC 61400-3-1:2019 §9 (electrical system)', 'DNV-ST-0359 (2021) §4 (subsea power cables for wind power plants)', 'CIGRE TB 379 (2009) and TB 815 (2020) (submarine cable service experience)'],
  models: ['Capacitated minimum spanning tree (Esau–Williams heuristic)', 'Cable selection by continuous current rating (IEC 60287)', 'Backward/forward-sweep AC load flow of a radial network (nominal-π branches)', 'Annual energy loss by production-level weighting (Weibull × power curve)', 'Expected energy not delivered for single cable faults (N-1)', 'Collection-voltage comparison'],
  limits: `Balanced steady-state load flow of the array cables only, with the substation busbar as the slack: no transformer tap control, harmonics, switching transients, short-circuit or protection study. The generic cable catalogue (three-core copper XLPE, wire armoured, buried 1 m in 0.7 K·m/W soil at 15 °C) is a planning aid: replace ratings, impedances and all costs with manufacturer data. The layout heuristic gives a good, not a proven optimal, cable length and ignores seabed constraints and exclusion zones. Reliability counts one cable fault at a time on radial strings; cable dynamics and the thermal rating at J-tubes are in the array-cable module. ${LIMIT_CODE}`,
  publishes: ['array.cableLength', 'cable.lossesFrac', 'array.cableCost', 'array.energyNotDeliveredFrac', 'quantities.cable_km'],
  searchTerms: 'offshore wind farm collection system array cable layout load flow losses 66 kV reliability',
  inputs: [
    { key: 'layout', label: 'Imported turbine coordinates (optional)', type: 'geometry', kinds: ['points'], default: null, group: 'Layout', role: 'geometry', help: 'Point set with x, y in metres (up to 150 positions). When empty the regular grid below is used.' },
    { key: 'nRows', label: 'Rows', type: 'int', default: 5, min: 1, max: 15, group: 'Layout', role: 'geometry', design: true },
    { key: 'nCols', label: 'Turbines per row', type: 'int', default: 8, min: 1, max: 15, group: 'Layout', role: 'geometry', design: true },
    { key: 'sx', label: 'Spacing along the rows', unit: 'm', type: 'number', default: 1000, min: 200, max: 5000, group: 'Layout', role: 'geometry', design: true },
    { key: 'sy', label: 'Spacing between rows', unit: 'm', type: 'number', default: 1400, min: 200, max: 5000, group: 'Layout', role: 'geometry', design: true },
    { key: 'subPos', label: 'Offshore substation position', type: 'select', default: 'centre', group: 'Layout', role: 'bc', options: [{ value: 'centre', label: 'Centre of the array' }, { value: 'custom', label: 'Coordinates below' }], help: 'The substation busbar is the slack node of the load flow.' },
    { key: 'subX', label: 'Substation x', unit: 'm', type: 'number', default: 0, min: -1e5, max: 1e5, group: 'Layout', role: 'bc', showIf: { key: 'subPos', in: ['custom'] }, help: 'In the frame of the layout; the grid starts with its first turbine at (0, 0).' },
    { key: 'subY', label: 'Substation y', unit: 'm', type: 'number', default: 0, min: -1e5, max: 1e5, group: 'Layout', role: 'bc', showIf: { key: 'subPos', in: ['custom'] } },
    { key: 'stringing', label: 'String layout', type: 'select', default: 'ew', group: 'Layout', role: 'control', options: [{ value: 'ew', label: 'Shortest cable (branching allowed)' }, { value: 'chain', label: 'Shortest cable, daisy chains only' }], help: 'Capacitated minimum spanning tree without cable crossings. Daisy chains suit two-way switchgear at every turbine.' },
    { key: 'maxPerString', label: 'Largest string (0 = what the biggest cable allows)', type: 'int', default: 8, min: 0, max: 20, group: 'Layout', role: 'criteria', design: true },
    { key: 'slack', label: 'Laid length / plan distance', unit: '-', type: 'number', default: 1.03, min: 1, max: 1.3, group: 'Layout', role: 'geometry', help: 'Seabed slope and micro-routing allowance.' },
    { key: 'endLength', label: 'Cable at each end (J-tube, deck, burial transition)', unit: 'm', type: 'number', default: 60, min: 0, max: 300, group: 'Layout', role: 'geometry' },
    { key: 'turbineMW', label: 'Rated power per turbine', unit: 'MW', type: 'number', default: 8, min: 1, max: 25, group: 'Turbines', role: 'bc', design: true, help: 'Constant-power injection at each turbine node.' },
    { key: 'pf', label: 'Power factor for cable sizing', unit: '-', type: 'number', default: 0.95, min: 0.8, max: 1, group: 'Turbines', role: 'criteria', help: `Sizing current I = n·P/(√3·U·cos φ). ${CONFIRM}` },
    { key: 'tanPhi', label: 'Reactive / active power of the turbines in operation', unit: '-', type: 'number', default: 0, min: -0.4, max: 0.4, group: 'Turbines', role: 'bc', help: 'Q/P injected (positive = exporting reactive power). Zero is unity power factor at the turbine terminals.' },
    { key: 'level', label: 'Production level for the part-load case', unit: '%', type: 'number', default: 50, min: 0, max: 100, group: 'Turbines', role: 'bc', help: 'The load flow is solved at rated power and at this share of it.' },
    { key: 'weibA', label: 'Weibull scale A at hub height', unit: 'm/s', type: 'number', default: 10.6, min: 4, max: 16, group: 'Wind resource', role: 'environment', bind: 'site.windWeibullA', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'weibK', label: 'Weibull shape k', unit: '-', type: 'number', default: 2.2, min: 1.2, max: 3.5, group: 'Wind resource', role: 'environment', bind: 'site.windWeibullK' },
    { key: 'voltage', label: 'Collection voltage', type: 'select', default: '66', group: 'Cables', role: 'bc', options: [{ value: '33', label: '33 kV' }, { value: '66', label: '66 kV' }], help: 'Slack voltage at the substation busbar.' },
    { key: 'freq', label: 'System frequency', unit: 'Hz', type: 'number', default: 50, min: 45, max: 65, group: 'Cables', role: 'bc' },
    { key: 'cableSource', label: 'Cable data', type: 'select', default: 'generic', group: 'Cables', role: 'material', options: [{ value: 'generic', label: 'Generic catalogue (150–800 mm² copper XLPE)' }, { value: 'table', label: 'Own cable table' }] },
    { key: 'cables', label: 'Own cables', type: 'table', group: 'Cables', role: 'material', showIf: { key: 'cableSource', in: ['table'] }, default: [{ area: 240, I: 480, r: 0.1, x: 0.12, c: 0.2, cost: 400 }, { area: 630, I: 780, r: 0.045, x: 0.1, c: 0.3, cost: 700 }],
      columns: [{ key: 'area', label: 'Section', unit: 'mm²', type: 'number' }, { key: 'I', label: 'Rating', unit: 'A', type: 'number' }, { key: 'r', label: 'AC resistance', unit: 'Ω/km', type: 'number' }, { key: 'x', label: 'Reactance', unit: 'Ω/km', type: 'number' }, { key: 'c', label: 'Capacitance', unit: 'µF/km', type: 'number' }, { key: 'cost', label: 'Supplied and installed', unit: 'USD/m', type: 'number' }],
      help: 'Manufacturer data per phase at maximum temperature, including screen and armour losses in the resistance.' },
    { key: 'compare', label: 'Compare 33 kV with 66 kV (slower)', type: 'bool', default: false, group: 'Cables', role: 'control', showIf: { key: 'cableSource', in: ['generic'] }, help: 'Repeats the layout, sizing, load flow and annual loss at both voltages with the generic cables.' },
    { key: 'maxTypes', label: 'Different cable sections allowed', type: 'int', default: 2, min: 1, max: 6, group: 'Cables', role: 'criteria', help: 'Two is usual for procurement and spares.' },
    { key: 'failRate', label: 'Cable failure rate', unit: '1/(km·yr)', type: 'number', default: 0.003, min: 0, max: 0.05, group: 'Reliability', role: 'criteria', uncertain: { dist: 'lognormal', cov: 0.5 }, help: `Internal and external faults per kilometre per year: of the order 0.001–0.01 in service statistics (CIGRE TB 379 / TB 815). ${CONFIRM}` },
    { key: 'mttr', label: 'Mean time to repair a cable', unit: 'days', type: 'number', default: 45, min: 1, max: 365, group: 'Reliability', role: 'operation', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Fault location, vessel mobilisation, weather window and jointing: one to three months is typical.' },
    { key: 'dUmax', label: 'Allowable voltage rise at the turbines', unit: '%', type: 'number', default: 3, min: 0.5, max: 10, group: 'Acceptance', role: 'criteria', help: `Above the busbar voltage at rated power. ${CONFIRM}` },
    { key: 'lossMax', label: 'Allowable annual array-cable loss', unit: '%', type: 'number', default: 1.0, min: 0.1, max: 5, group: 'Acceptance', role: 'criteria' },
    { key: 'endMax', label: 'Allowable energy not delivered (cable faults)', unit: '%', type: 'number', default: 0.3, min: 0.01, max: 5, group: 'Acceptance', role: 'criteria' },
    costInput('cableBase', 'Array cable, base cost supplied and installed', 'USD/m', 280, 5000, 'Part of the cost per metre that does not depend on the section (lay, burial, armour, terminations).'),
    costInput('cablePerMm2', 'Array cable, cost per mm² of copper section', 'USD/m per mm²', 0.6, 20),
    costInput('energyPrice', 'Value of delivered energy', 'USD/MWh', 80, 1000, 'Used to express losses and energy not delivered as money.'),
    { key: 'nBins', label: 'Production levels for the annual loss', type: 'int', default: 16, min: 2, max: 256, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nBins', label: 'Production levels', levels: [2, 4, 8, 16, 32], qoi: 'lossAnnual' },
  compute(inp) {
    const MW = 1e6, U = Number(inp.voltage) * 1e3, P = inp.turbineMW * MW;
    // turbine positions: imported points or the regular grid
    let turbines = [];
    const g = inp.layout;
    if (g && g.positions && g.positions.length >= 3) {
      const n = Math.floor(g.positions.length / 3);
      need(n <= 150, `The imported layout has ${n} positions; the layout heuristic in this module is limited to 150 turbines per run.`);
      for (let i = 0; i < n; i++) { const x = g.positions[3 * i], y = g.positions[3 * i + 1]; need(Number.isFinite(x) && Number.isFinite(y), 'The imported layout contains non-numeric coordinates.'); turbines.push({ x, y }); }
    } else {
      need(inp.nRows * inp.nCols <= 150, 'The grid is limited to 150 turbines per run.');
      for (let j = 0; j < inp.nRows; j++) for (let i = 0; i < inp.nCols; i++) turbines.push({ x: i * inp.sx, y: j * inp.sy });
    }
    const n = turbines.length;
    let sub;
    if (inp.subPos === 'custom') sub = { x: inp.subX, y: inp.subY };
    else {
      sub = { x: turbines.reduce((s, t) => s + t.x, 0) / n, y: turbines.reduce((s, t) => s + t.y, 0) / n };
      // a platform cannot stand on a turbine position: step half a spacing aside
      const dMin = Math.min(...turbines.map((t) => Math.hypot(t.x - sub.x, t.y - sub.y)));
      if (dMin < 100) sub = { x: sub.x + 0.5 * (n > 1 ? Math.max(200, Math.min(...turbines.slice(1).map((t) => Math.hypot(t.x - turbines[0].x, t.y - turbines[0].y)))) : 500), y: sub.y };
    }
    for (let i = 0; i < n; i++) {
      need(Math.hypot(turbines[i].x - sub.x, turbines[i].y - sub.y) > 1, `Turbine ${i + 1} coincides with the substation: move the substation.`);
      for (let j = i + 1; j < n; j++) need(Math.hypot(turbines[i].x - turbines[j].x, turbines[i].y - turbines[j].y) > 1, `Turbines ${i + 1} and ${j + 1} share one position.`);
    }
    // cable catalogue
    const costFn = (area) => inp.cableBase + inp.cablePerMm2 * area * 1e6;
    let cat;
    if (inp.cableSource === 'table') {
      const rows = (inp.cables || []).filter((r) => r && r.I > 0);
      need(rows.length >= 1, 'The cable table needs at least one row with a positive rating.');
      for (const r of rows) need(r.area > 0 && r.r > 0 && r.x >= 0 && r.c >= 0 && r.cost >= 0, 'Cable table: section and resistance must be positive; reactance, capacitance and cost must not be negative.');
      cat = rows.map((r) => ({ area: r.area * 1e-6, voltage: U, I: Number(r.I), r: r.r / 1e3, x: r.x / 1e3, c: r.c * 1e-9, cost: Number(r.cost) }));   // Ω/km → Ω/m, µF/km → F/m
    } else cat = run(() => catalogueFor(U, inp.freq, costFn), 'Cable catalogue');
    const Imax = Math.max(...cat.map((c) => c.I)), I1 = P / (Math.sqrt(3) * U * inp.pf), capMax = Math.floor(Imax / I1 + 1e-9);
    need(capMax >= 1, `The largest cable (${fmt(Imax)} A) cannot carry one ${inp.turbineMW} MW turbine (${fmt(I1)} A) at ${inp.voltage} kV: raise the voltage or add a larger cable.`);
    const mps = inp.maxPerString > 0 ? Math.min(inp.maxPerString, capMax) : capMax;
    const ew = run(() => esauWilliams({ substation: sub, turbines, maxPerString: mps, noCrossing: true, noBranching: inp.stringing === 'chain' }), 'String layout');
    const net = run(() => buildNetwork({ substation: sub, turbines, parent: ew.parent, slackFactor: inp.slack, endLength: inp.endLength }), 'Network');
    const sizing = run(() => sizeCables({ network: net, catalogue: cat, turbinePower: P, voltage: U, powerFactor: inp.pf, maxTypes: inp.maxTypes }), 'Cable sizing');
    const lf = { network: net, catalogue: cat, type: sizing.type, voltage: U, f: inp.freq };
    const flow = (frac) => run(() => radialLoadFlow({ ...lf, P: frac * P, Q: frac * P * inp.tanPhi }), 'Load flow');
    const rated = flow(1), part = flow(inp.level / 100);
    const dist = productionLevels({ A: inp.weibA, k: inp.weibK, nBins: Math.round(inp.nBins) });
    const annual = run(() => annualEnergyLoss({ ...lf, turbinePower: P, distribution: dist, powerFactorTan: inp.tanPhi }), 'Annual losses');
    const cf = annual.energyGenerated / (n * P * HOURS_YEAR);
    const rel = run(() => energyNotDelivered({ network: net, turbinePower: P, capacityFactor: cf, failureRate: inp.failRate / 1e3, mttr: inp.mttr * 24 }), 'Reliability');   // per km → per m; days → h
    const nStrings = ew.strings.length, longest = Math.max(...ew.strings.map((s) => s.length));
    const dU = (rated.Umax / U - 1) * 100, utilMax = Math.max(...Array.from(sizing.utilisation).slice(1));
    const km = net.totalLength / 1e3, qty = cableQuantities(sizing.quantities);
    const lossAnnualPct = annual.lossFraction * 100, endPct = rel.unavailability * 100;
    const lossMWh = annual.energyLoss / MW, endMWh = rel.energyNotDelivered / MW;
    // conductor copper: three cores per cable
    const copperT = sizing.quantities.reduce((s, q) => s + 3 * q.area * q.length * RHO_CU, 0) / 1e3;

    const checks = [
      check('Cable current against the continuous rating (most loaded section)', utilMax, 1, '-', 'IEC 60287-1-1:2023 §1.4; DNV-ST-0359 §4',
        'Sizing current of each section over the rating of the cable chosen for it. Above 1 the conductor exceeds its maximum temperature at full output and the insulation ages faster.'),
      check('Voltage rise at the most remote turbine (rated power)', Math.max(dU, 0), inp.dUmax, '%', 'IEC 61400-3-1:2019 §9; IEC 60038 voltage bands',
        'Rise of the turbine terminal voltage above the substation busbar at full output. Too much rise pushes the turbine transformers and converters out of their voltage band and trips them.'),
      check('Annual array-cable loss', lossAnnualPct, inp.lossMax, '%', 'IEC 60287-3-2:2012 (economic optimisation of cable size); project loss budget',
        'Energy dissipated in the array cables over a year as a share of the energy generated. Exceeding the budget lowers the net yield on which the project is financed.'),
      check('Energy not delivered after single cable faults', endPct, inp.endMax, '%', 'CIGRE TB 379 / TB 815 failure statistics; N-1 assessment of radial strings',
        'Expected production stranded behind cable faults while they are repaired, as a share of the annual energy. Above the allowance, shorter strings, ring connections or a faster repair contract are needed.'),
    ];
    const gov = worst(checks);

    // 33 kV against 66 kV with the generic catalogue (same layout rule and sizing)
    let cmpRows = [];
    if (inp.cableSource === 'generic' && inp.compare) {
      const cmp = run(() => compareCollectionVoltages({ substation: sub, turbines, turbinePower: P, f: inp.freq, powerFactor: inp.pf, maxTypes: inp.maxTypes, maxPerString: inp.maxPerString > 0 ? inp.maxPerString : undefined, distribution: dist, slackFactor: inp.slack, endLength: inp.endLength,
        catalogueFor: (v) => catalogueFor(v, inp.freq, costFn) }), 'Voltage comparison');
      cmpRows = cmp.map((c) => ({ v: c.voltage / 1e3, strings: c.nStrings, perString: c.maxPerString, km: round(c.cableLength / 1e3), cost: round(c.cableCost / 1e6), lossRated: round(c.lossFractionAtRated * 100), lossAnnual: c.annualLossFraction != null ? round(c.annualLossFraction * 100) : 0, dU: round((c.Umax / c.voltage - 1) * 100) }));
    }
    // worst string: the path from the node with the highest voltage back to the substation
    let far = 1; for (let i = 1; i <= n; i++) if (rated.U[i] > rated.U[far]) far = i;
    const path = []; for (let i = far; i > 0; i = net.parent[i]) path.push(i); path.reverse();
    let dist0 = 0; const prof = [{ d: 0, u: rated.U[0] / U }];
    for (const i of path) { dist0 += net.length[i]; prof.push({ d: dist0 / 1e3, u: rated.U[i] / U }); }
    const fr = linspace(0, 1, 11), lossCurve = fr.map((f) => flow(f).losses / 1e3);
    const lossValue = lossMWh * inp.energyPrice, endValue = endMWh * inp.energyPrice, layDays = km / 4;

    return {
      kpis: [
        { key: 'nTurbines', label: 'Turbines connected', value: n, unit: '-', explain: `Installed capacity ${fmt(n * inp.turbineMW)} MW collected at ${inp.voltage} kV.` },
        { key: 'nStrings', label: 'Strings', value: nStrings, unit: '-', explain: `Feeders leaving the substation, each needing a switchgear bay. The longest carries ${longest} turbine(s); the largest cable allows ${capMax}.` },
        { key: 'cableKm', label: 'Array cable length', value: round(km, 5), unit: 'km', explain: 'Laid length of all sections including the slack allowance and the cable in the J-tubes at each end.' },
        { key: 'utilMax', label: 'Highest cable utilisation', value: round(utilMax, 4), unit: '-', limit: 1, sense: 'max', explain: 'Sizing current over the continuous rating for the most loaded section. Well below 1 on many sections means the cable sizes could be stepped down along the string.' },
        { key: 'lossRated', label: 'Array-cable loss at rated power', value: round(rated.losses / 1e3, 5), unit: 'kW', explain: `Three-phase loss in all sections with every turbine at full output: ${fmt(rated.lossFraction * 100)} % of the generated power.` },
        { key: 'lossLevel', label: 'Array-cable loss at the part-load level', value: round(part.losses / 1e3, 5), unit: 'kW', explain: `Loss with every turbine at ${fmt(inp.level)} % of rated power. Resistive loss falls with the square of the output; the dielectric charging current keeps it from reaching zero.` },
        { key: 'lossAnnual', label: 'Annual array-cable loss', value: round(lossAnnualPct, 5), unit: '%', limit: inp.lossMax, sense: 'max', explain: `Loss energy over generated energy across the wind-speed distribution: ${fmt(lossMWh)} MWh a year. Lower than the rated-power figure because the farm spends most hours at part load.` },
        { key: 'dU', label: 'Voltage rise at the most remote turbine', value: round(dU, 4), unit: '%', limit: inp.dUmax, sense: 'max', explain: 'Highest turbine terminal voltage above the substation busbar at rated power.' },
        { key: 'qSub', label: 'Reactive power into the substation at rated power', value: round(rated.Qsubstation / 1e6, 4), unit: 'Mvar', explain: 'Positive values flow into the busbar. Cable charging dominates at low output, the series reactance at high output; the difference sizes the reactive compensation.' },
        { key: 'end', label: 'Energy not delivered (single cable faults)', value: round(endPct, 4), unit: '%', limit: inp.endMax, sense: 'max', explain: `Expected ${fmt(endMWh)} MWh a year stranded behind ${fmt(rel.expectedFailures, 2)} cable fault(s) a year, each isolating the turbines downstream for ${fmt(inp.mttr)} days.` },
        { key: 'endMWh', label: 'Energy not delivered', value: round(endMWh, 5), unit: 'MWh/yr', explain: 'Σ over the sections of failure rate × length × repair time × turbines downstream × mean output. Sections next to the substation count most: they strand the whole string.' },
        { key: 'cableCost', label: 'Array cable cost (indicative)', value: round(sizing.totalCost / 1e6, 4), unit: 'M USD', explain: 'Length of each section times the indicative supplied-and-installed rate for its conductor size. A comparison figure, not an estimate.' },
        { key: 'cf', label: 'Capacity factor of the production distribution', value: round(cf * 100, 4), unit: '%', explain: 'Mean output over rated from the Weibull distribution and the generic power curve (before wake and availability losses); it weights both the losses and the energy not delivered.' },
      ],
      checks,
      plots: [
        { type: 'model3d', title: 'String layout and voltage at rated power', nodes: net.nodes.map((p) => [round(p.x, 6), round(p.y, 6), 0]), lines: Array.from({ length: n }, (_, i) => [i + 1, net.parent[i + 1]]), tris: [], nodeValues: Array.from(rated.U, (v) => round(v / U, 5)), valueLabel: 'Voltage', valueUnit: 'p.u.',
          explain: 'Plan view of the turbines, the substation (node at the hub of the strings) and the cables, coloured by voltage. The highest voltage is at the far end of the longest, most heavily loaded string.' },
        { type: 'line', title: 'Array-cable loss against production level', x: { label: 'Production level', unit: '% of rated' }, y: { label: 'Loss', unit: 'kW' }, series: [{ name: 'Total loss', x: fr.map((f) => round(f * 100)), y: lossCurve.map((v) => round(v)) }],
          explain: 'Loss in all array cables as the farm output rises. The curve is close to a parabola (I²R); its value at zero output is the loss from the charging current alone.' },
        { type: 'line', title: 'Voltage along the string with the highest voltage (rated power)', x: { label: 'Cable distance from the substation', unit: 'km' }, y: { label: 'Voltage', unit: 'p.u.' },
          series: [{ name: 'Node voltage', x: prof.map((p) => round(p.d)), y: prof.map((p) => round(p.u, 6)) }, { name: 'Limit', x: [0, round(prof[prof.length - 1].d)], y: [1 + inp.dUmax / 100, 1 + inp.dUmax / 100], dash: true }],
          explain: 'Voltage at each turbine from the substation outwards. The steep part next to the substation carries the whole string: a larger section there buys the most voltage and loss reduction.' },
        { type: 'bar', title: 'Probability of the production levels', categories: dist.map((d) => `${fmt(d.fraction * 100, 2)} %`), series: [{ name: 'Share of the year (%)', y: dist.map((d) => round(d.probability * 100, 3)) }], x: { label: 'Production level', unit: '' }, y: { label: 'Share of the year', unit: '%' },
          explain: 'Weights of the load-flow cases behind the annual loss. The bars at zero and at rated power are the hours below cut-in or above cut-out, and between rated wind speed and cut-out.' },
      ],
      tables: [
        { title: 'Cable schedule', columns: [{ key: 'size', label: 'Section', unit: 'mm²' }, { key: 'km', label: 'Length', unit: 'km' }, { key: 'I', label: 'Rating', unit: 'A' }, { key: 'cost', label: 'Cost', unit: 'M USD' }],
          rows: sizing.quantities.map((q, i) => ({ size: qty[i].sizeMm2, km: round(qty[i].km), I: round(cat.find((c) => Math.abs(c.area - q.area) < 1e-9)?.I ?? 0), cost: round(q.cost / 1e6) })),
          explain: 'One row per conductor size used: the bill of cable for procurement and the cost model.' },
        { title: 'Strings', columns: [{ key: 's', label: 'String' }, { key: 'n', label: 'Turbines' }, { key: 'km', label: 'Cable', unit: 'km' }, { key: 'I', label: 'Current at the substation end', unit: 'A' }, { key: 'u', label: 'Highest voltage', unit: 'p.u.' }, { key: 'loss', label: 'Loss at rated', unit: 'kW' }],
          rows: ew.strings.map((s, k) => { const ids = s.map((t) => t + 1), head = ids.find((i) => net.parent[i] === 0); return { s: k + 1, n: s.length, km: round(ids.reduce((q, i) => q + net.length[i], 0) / 1e3), I: round(rated.I[head]), u: round(Math.max(...ids.map((i) => rated.U[i])) / U, 5), loss: round(ids.reduce((q, i) => q + rated.branchLoss[i], 0) / 1e3) }; }),
          explain: 'Each feeder with its turbine count, cable length, current in the first section, highest node voltage and loss at full output. Uneven strings show where rebalancing would help.' },
        { title: '33 kV against 66 kV (generic cables, same layout rule)', columns: [{ key: 'v', label: 'Voltage', unit: 'kV' }, { key: 'strings', label: 'Strings' }, { key: 'perString', label: 'Largest string' }, { key: 'km', label: 'Cable', unit: 'km' }, { key: 'cost', label: 'Cable cost', unit: 'M USD' }, { key: 'lossRated', label: 'Loss at rated', unit: '%' }, { key: 'lossAnnual', label: 'Annual loss', unit: '%' }, { key: 'dU', label: 'Voltage rise', unit: '%' }],
          rows: cmpRows, explain: 'The same farm collected at each voltage. 66 kV carries about twice the power per cable, so fewer strings and less cable; the table is filled when the comparison is switched on with the generic cables.' },
      ],
      notes: [
        `${n} turbines of ${inp.turbineMW} MW in ${nStrings} string(s) of up to ${longest}; ${fmt(km)} km of cable in ${sizing.quantities.length} section size(s); the load flow converged in ${rated.iterations} sweep(s).`,
        `${gov.name} governs at ${fmt(gov.uc)}. ${ew.crossings > 0 ? `${ew.crossings} cable crossing(s) remain in the layout.` : 'The layout has no cable crossings.'}`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)}${gov.unit === '-' ? '' : ` ${gov.unit}`} (utilisation ${fmt(gov.uc)}). Annual loss ${fmt(lossAnnualPct)} %, voltage rise ${fmt(dU)} %, energy not delivered ${fmt(endPct)} %.`,
          gov.uc > 1 ? (gov === checks[1] || gov === checks[2] ? `Shorten the strings (largest now ${longest}; try ${Math.max(1, longest - 1)}), allow a larger section next to the substation, or collect at ${inp.voltage === '33' ? '66' : 'a higher'} kV.` : gov === checks[3] ? `Limit the strings to ${Math.max(1, longest - 2)} turbines, close pairs of strings into rings (the outage then lasts only the switching time), or contract a repair time below ${fmt(inp.mttr * inp.endMax / Math.max(endPct, 1e-9))} days.` : 'Reduce the string size or add a larger cable to the catalogue.')
            : `The collection system meets the criteria entered. Confirm the cable ratings at the J-tubes and landfall hot spots with the array-cable module, and check reactive compensation for ${fmt(rated.Qsubstation / 1e6)} Mvar at rated power.`,
          'IEC 60287 (ratings) and IEC 60287-3-2 (economic sizing); IEC 61400-3-1:2019 §9; DNV-ST-0359 (2021) §4; CIGRE TB 815 (2020) for failure statistics.',
          `The conductors hold about ${fmt(copperT)} t of copper ≈ ${fmt(copperT * CO2_COPPER)} t CO₂e (screening factor ${CO2_COPPER} t/t); the array-cable losses are ${fmt(lossMWh)} MWh a year of renewable energy not delivered. ${vesselCo2Text(layDays, 25)} Shorter total cable length also means less seabed disturbed by burial.`),
        adv('economic', 'info', `At the indicative rates the array cable costs ${money(sizing.totalCost)}; losses are worth ${money(lossValue)} a year and energy not delivered after cable faults ${money(endValue)} a year.`,
          `${cmpRows.length === 2 ? `Collecting at ${cmpRows[0].cost < cmpRows[1].cost ? cmpRows[0].v : cmpRows[1].v} kV has the lower cable cost (${fmt(cmpRows[0].cost)} against ${fmt(cmpRows[1].cost)} M USD for 33 and 66 kV) before switchgear and transformer differences. ` : ''}Compare a larger section on the first section of each string with the capitalised value of its losses, and a ring connection with the ${money(endValue)} a year it would largely remove.`,
          'IEC 60287-3-2 economic optimisation of cable size; whole-life costing per ISO 15686-5 (capitalised losses and unavailability against capital cost); AACE RP 18R-97 for the class of the cost figures.',
          energyCo2Text(lossMWh + endMWh)),
      ],
      publish: { 'array.cableLength': net.totalLength, 'cable.lossesFrac': annual.lossFraction, 'array.cableCost': sizing.totalCost, 'array.energyNotDeliveredFrac': rel.unavailability, 'quantities.cable_km': km },
      quantities: { cable_km: round(km, 5), copper_t: round(copperT, 4), vessel_days: { cable_lay: round(layDays, 3) } },
    };
  },
  verify: [
    { // One 8 MW turbine on 1 km of cable, r = 0.1 Ω/km, unity power factor at 33 kV: I = 8e6/(√3·33e3) = 139.96 A;
      // loss = 3·I²·R = 3·139.96²·0.1 = 5877 W (the 24 V of voltage rise lowers the current by under 0.1 %).
      name: 'Single turbine on one cable: I²R loss by hand', source: 'Hand calculation: I = P/(√3 U), P_loss = 3 I² R',
      inputs: { nRows: 1, nCols: 1, subPos: 'custom', subX: -1000, subY: 0, slack: 1, endLength: 0, turbineMW: 8, pf: 1, tanPhi: 0, voltage: '33', cableSource: 'table', cables: [{ area: 240, I: 480, r: 0.1, x: 0, c: 0, cost: 400 }], level: 50 },
      expect: [{ kpi: 'lossRated', value: 5.877, tol: 0.003 }, { kpi: 'lossLevel', value: 5.877 / 4, tol: 0.003 }, { kpi: 'cableKm', value: 1, tol: 1e-9 }, { kpi: 'cableCost', value: 0.4, tol: 1e-9 }] },
    { // Three turbines 1 km apart in one chain on the same cable: sections carry 3, 2 and 1 turbine currents, so the loss
      // is (9 + 4 + 1) × 5877 W = 82.28 kW (voltage rise under 0.3 % lowers it slightly). Voltage rise at the far end:
      // √3·R·I·(3 + 2 + 1) = 1.732 × 0.1 × 139.96 × 6 = 145.4 V = 0.441 % of 33 kV.
      name: 'Chain of three turbines: loss ∝ Σ n² and voltage rise ∝ Σ n', source: 'Hand calculation (Kersting, Distribution System Modeling and Analysis, ch. 10: series voltage drop of a radial feeder)',
      inputs: { nRows: 1, nCols: 3, sx: 1000, subPos: 'custom', subX: -1000, subY: 0, slack: 1, endLength: 0, turbineMW: 8, pf: 1, tanPhi: 0, voltage: '33', cableSource: 'table', cables: [{ area: 240, I: 480, r: 0.1, x: 0, c: 0, cost: 400 }] },
      expect: [{ kpi: 'lossRated', value: 82.28, tol: 0.01 }, { kpi: 'dU', value: 0.441, tol: 0.01 }, { kpi: 'nStrings', value: 1, abs: 0 }, { kpi: 'cableKm', value: 3, tol: 1e-9 }] },
    { // N-1 for the same chain. The capacity factor multiplies both the stranded and the annual energy and cancels:
      // Σ n_b = 3 + 2 + 1 = 6; stranded = 0.01 /km/yr × 1 km × 240 h × 6 × 8 MW = 115.2 MWh/yr per unit capacity factor;
      // annual = 3 × 8 MW × 8760 h = 210 240 MWh per unit capacity factor; ratio 0.05479 %.
      name: 'Energy not delivered after single faults on a chain: λ·L·MTTR·Σn / (N·8760)', source: 'Hand calculation (expected-value N-1 assessment of a radial feeder)',
      inputs: { nRows: 1, nCols: 3, sx: 1000, subPos: 'custom', subX: -1000, subY: 0, slack: 1, endLength: 0, turbineMW: 8, pf: 1, voltage: '33', cableSource: 'table', cables: [{ area: 240, I: 480, r: 0.1, x: 0, c: 0, cost: 400 }], failRate: 0.01, mttr: 10 },
      expect: [{ kpi: 'end', value: 0.05479, tol: 2e-3 }] },
  ],
};
