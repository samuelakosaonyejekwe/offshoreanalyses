// Renewable-energy economics: levelised cost of energy (real and nominal), delivered energy, storage, avoided
// cost, contracted and merchant revenue, break-even tariff and the LCOE sensitivity map.
import { lcoe, lcoeRealNominal, lcoeFixedChargeRate, lcos, lcoeDelivered, lace, blendedRevenue } from '../../engine/econ/lcoe.js';
import { renewableGeneration } from '../../engine/econ/production.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, CO2_STEEL, need, run, round, fmt, mil, fM, pct, adv, check, cost, range, crf } from './_common.js';

const PMWH = `${CUR}/MWh`;

export default {
  id: 's17-renewable-lcoe', suite: 17, title: 'Levelised cost of energy and renewable revenue',
  summary: 'Computes the real and nominal LCOE of a generating project with degradation, losses, replacement and decommissioning, the cost of delivered and stored energy, the value of the energy displaced, and revenue under a PPA, a contract for difference and merchant sales.',
  covers: [5], fidelity: 'analytical',
  standards: ['IEA/NEA, Projected Costs of Generating Electricity (2020 ed.), methodology annex (discounted LCOE)', 'NREL/TP-462-5173, Short et al. (1995), Manual for the Economic Evaluation of Energy Efficiency and Renewable Energy Technologies', 'IRENA, Renewable Power Generation Costs (annual; LCOE conventions, real pre-tax)', 'Schmidt et al. (2019), Joule 3(1) (levelised cost of storage)', 'US EIA, Levelized Costs of New Generation Resources (LACE definition)'],
  models: ['Levelised cost of energy (discounting method)', 'Fixed-charge-rate LCOE cross-check', 'Real versus nominal LCOE (Fisher)', 'Levelised cost of delivered energy', 'Levelised cost of storage', 'Levelised avoided cost of energy and net value', 'Capacity-factor and degradation-adjusted energy yield', 'Replacement cost model', 'Power purchase agreement, contract-for-difference and merchant revenue', 'Break-even tariff'],
  limits: `Pre-tax LCOE in real terms (the IEA/IRENA/NREL convention) with annual energy and year-end cash flows; taxes, financing structure and intra-year price shape are excluded (see the project-finance module). Market price paths are flat real inputs. ${CUR_NOTE}`,
  publishes: ['econ.lcoe', 'econ.lcoeNominal', 'econ.lace', 'econ.breakEvenTariff'],
  searchTerms: 'levelized cost of energy offshore wind LCOE convention real nominal degradation',
  inputs: [
    { key: 'baseYear', label: 'Price base year', unit: 'yr', type: 'int', default: 2026, min: 1990, max: 2100, group: 'Opening financial state', role: 'ic' },
    { key: 'inflation', label: 'General inflation', unit: '%/yr', type: 'number', default: 2, min: -2, max: 30, group: 'Opening financial state', role: 'ic', help: 'Used only to express the nominal LCOE and nominal contract prices.' },
    { key: 'ratedMW', label: 'Installed capacity', unit: 'MW', type: 'number', default: 500, min: 0.1, max: 20000, group: 'Energy generation', role: 'operation', bind: 'energy.ratedPower' },
    { key: 'aep', label: 'Annual energy production, year 1', unit: 'GWh/yr', type: 'number', default: 2210, min: 0.001, max: 2e5, group: 'Energy generation', role: 'operation', bind: 'energy.aep', uncertain: { dist: 'normal', cov: 0.08 },
      help: 'As delivered by the energy-yield suite (net of wake and availability losses). Export losses and curtailment are applied below.' },
    { key: 'lossesFrac', label: 'Export transmission losses', unit: '%', type: 'number', default: 2.5, min: 0, max: 40, group: 'Energy generation', role: 'operation', bind: 'cable.lossesFrac', help: 'From the cable electrical analysis when available.' },
    { key: 'curtailment', label: 'Curtailment', unit: '%', type: 'number', default: 1, min: 0, max: 60, group: 'Energy generation', role: 'bc', help: 'Grid or market constraint on output.' },
    { key: 'degradation', label: 'Annual output degradation', unit: '%/yr', type: 'number', default: 0.3, min: 0, max: 10, group: 'Energy generation', role: 'operation' },
    { key: 'sigmaYield', label: 'Uncertainty of long-term net yield (1σ)', unit: '%', type: 'number', default: 6, min: 0, max: 40, group: 'Energy generation', role: 'operation', help: 'Standard deviation of the yield estimate, for the P50/P90 energy cases.' },
    { key: 'constructionYears', label: 'Construction period', unit: 'yr', type: 'int', default: 3, min: 1, max: 10, group: 'Timeline', role: 'operation' },
    { key: 'opLife', label: 'Operating life', unit: 'yr', type: 'int', default: 25, min: 1, max: 60, group: 'Timeline', role: 'bc', bind: 'life.designYears', design: true, help: 'Operational constraint from the design or fatigue life.' },
    cost({ key: 'capex', label: 'Capital expenditure', unit: CUR, default: 1.5e9, min: 1, max: 1e12, group: 'Costs (base-year money)', bind: 'econ.capex', design: true, uncertain: { dist: 'lognormal', cov: 0.15 } }),
    cost({ key: 'opexFixed', label: 'Fixed operating cost', unit: `${CUR}/yr`, default: 4.5e7, min: 0, max: 1e11, group: 'Costs (base-year money)', bind: 'econ.opexAnnual', uncertain: { dist: 'normal', cov: 0.15 } }),
    cost({ key: 'opexVar', label: 'Variable operating cost', unit: PMWH, default: 1.5, min: 0, max: 1e4, group: 'Costs (base-year money)' }),
    cost({ key: 'replCost', label: 'Mid-life replacement cost', unit: CUR, default: 6.0e7, min: 0, max: 1e12, group: 'Costs (base-year money)', help: 'Major component exchange or repowering allowance.' }),
    { key: 'replYear', label: 'Operating year of the replacement', unit: 'yr', type: 'int', default: 13, min: 1, max: 60, group: 'Costs (base-year money)', role: 'operation' },
    cost({ key: 'decom', label: 'Decommissioning cost', unit: CUR, default: 1.2e8, min: 0, max: 1e12, group: 'Costs (base-year money)', bind: 'econ.decomCost', help: 'Terminal cost in the year after the last operating year.' }),
    { key: 'discountRate', label: 'Real discount rate (pre-tax WACC)', unit: '%/yr', type: 'number', default: 6, min: 0, max: 40, group: 'Financial constraints', role: 'bc', uncertain: { dist: 'normal', std: 0.75 }, help: `Real rate matching the constant-money costs. ${PLACEHOLDER}` },
    cost({ key: 'merchant', label: 'Expected market (merchant) price', unit: PMWH, default: 88, min: 0, max: 5000, group: 'Revenue', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Externally varying market price, flat in real terms; captured price after shape effects.' }),
    { key: 'ppaShare', label: 'Share of output under the PPA', unit: '%', type: 'number', default: 40, min: 0, max: 100, group: 'Revenue', role: 'bc', help: 'Fixed contractual price for this share.' },
    cost({ key: 'ppaPrice', label: 'PPA price, year 1', unit: PMWH, default: 106, min: 0, max: 5000, group: 'Revenue' }),
    { key: 'ppaEscalation', label: 'PPA price escalation (real)', unit: '%/yr', type: 'number', default: 0, min: -10, max: 10, group: 'Revenue', role: 'bc', help: 'Zero means fully inflation-indexed; a fixed nominal PPA has a real escalation of about minus the inflation rate.' },
    { key: 'ppaTerm', label: 'PPA term', unit: 'yr', type: 'int', default: 15, min: 0, max: 60, group: 'Revenue', role: 'bc' },
    { key: 'cfdShare', label: 'Share of output under the contract for difference', unit: '%', type: 'number', default: 40, min: 0, max: 100, group: 'Revenue', role: 'bc' },
    cost({ key: 'cfdStrike', label: 'CfD strike price', unit: PMWH, default: 114, min: 0, max: 5000, group: 'Revenue' }),
    { key: 'cfdTerm', label: 'CfD term', unit: 'yr', type: 'int', default: 15, min: 0, max: 60, group: 'Revenue', role: 'bc' },
    { key: 'capacityCredit', label: 'Capacity credit', unit: '%', type: 'number', default: 12, min: 0, max: 100, group: 'Avoided cost', role: 'operation', help: 'Share of installed capacity counted as firm by the system operator.' },
    cost({ key: 'capacityPayment', label: 'Capacity payment', unit: `${CUR}/MW/yr`, default: 45000, min: 0, max: 1e7, group: 'Avoided cost' }),
    cost({ key: 'txCapex', label: 'Transmission asset capital cost', unit: CUR, default: 2.6e8, min: 0, max: 1e12, group: 'Delivered energy', help: 'Export system not already in the project CAPEX (enter 0 if it is).' }),
    cost({ key: 'txOpex', label: 'Transmission operating cost', unit: `${CUR}/yr`, default: 5.0e6, min: 0, max: 1e11, group: 'Delivered energy' }),
    { key: 'storage', label: 'Include a storage system', type: 'bool', default: true, group: 'Storage', role: 'control' },
    cost({ key: 'stCapex', label: 'Storage capital cost', unit: CUR, default: 1.1e8, min: 1, max: 1e12, group: 'Storage', showIf: { key: 'storage', in: [true] } }),
    { key: 'stEnergy', label: 'Usable storage capacity', unit: 'MWh (storage)', type: 'number', default: 400, min: 0.1, max: 1e6, group: 'Storage', role: 'operation', showIf: { key: 'storage', in: [true] } },
    { key: 'stCycles', label: 'Full cycles per year', unit: '1/yr', type: 'number', default: 300, min: 1, max: 3000, group: 'Storage', role: 'operation', showIf: { key: 'storage', in: [true] } },
    { key: 'stEff', label: 'Round-trip efficiency', unit: '%', type: 'number', default: 86, min: 20, max: 100, group: 'Storage', role: 'operation', showIf: { key: 'storage', in: [true] } },
    { key: 'stLife', label: 'Storage life', unit: 'yr', type: 'int', default: 15, min: 1, max: 50, group: 'Storage', role: 'operation', showIf: { key: 'storage', in: [true] } },
    cost({ key: 'stOpex', label: 'Storage operating cost', unit: `${CUR}/yr`, default: 2.2e6, min: 0, max: 1e10, group: 'Storage', showIf: { key: 'storage', in: [true] } }),
    cost({ key: 'gridIntensity', label: 'Carbon intensity of displaced generation', unit: 'g CO₂/kWh', default: 350, min: 0, max: 1200, group: 'Avoided cost', help: 'Marginal or average grid emission factor of the receiving system.' }),
  ],
  compute(inp) {
    const n = Math.round(inp.opLife); const nc = Math.round(inp.constructionYears); const r = inp.discountRate / 100; const f = inp.inflation / 100;
    need(inp.ppaShare + inp.cfdShare <= 100 + 1e-9, `The PPA and CfD shares add up to ${inp.ppaShare + inp.cfdShare} % of output: contracted shares cannot exceed 100 %.`);
    need(inp.replYear <= n || inp.replCost === 0, `The replacement is scheduled in operating year ${inp.replYear}, after the ${n}-year operating life.`);
    const gen = run(() => renewableGeneration({ grossAep_MWh: inp.aep * 1e3, availability: 1, lossesFrac: inp.lossesFrac / 100, curtailment: inp.curtailment / 100, degradation: inp.degradation / 100, years: n, sigma: inp.sigmaYield / 100, exceedance: [50, 75, 90, 99] }), 'Energy yield');
    const e = gen.p50; // MWh per operating year
    const N = nc + n + 1; // rows; operations in rows nc … nc+n−1; decommissioning in row nc+n
    const capex = range(N, (t) => (t < nc ? inp.capex / nc : 0));
    capex[nc - 1 + Math.round(inp.replYear)] += inp.replCost;
    const energy = range(N, (t) => (t >= nc && t < nc + n ? e[t - nc] : 0));
    const opex = range(N, (t) => (t >= nc && t < nc + n ? inp.opexFixed + inp.opexVar * e[t - nc] : 0));
    const decom = range(N, (t) => (t === nc + n ? inp.decom : 0));
    const L = run(() => lcoe({ capex, opex, decom, energy, rate: r }), 'LCOE');
    const rn = run(() => lcoeRealNominal({ capex, opex, decom, energy, rate: r, inflation: f }), 'Real/nominal LCOE');
    const eMean = e.reduce((a, c) => a + c, 0) / n;
    const fcr = lcoeFixedChargeRate({ capex: inp.capex, rate: r, life: n, fixedOm: inp.opexFixed, annualEnergy: eMean, variableOm: inp.opexVar });
    const p90e = gen.cases.P90 || gen.cases.p90 || gen.cases['90'];
    need(Array.isArray(p90e), 'Energy exceedance cases could not be evaluated.');
    const L90 = lcoe({ capex, opex: range(N, (t) => (t >= nc && t < nc + n ? inp.opexFixed + inp.opexVar * p90e[t - nc] : 0)), decom, energy: range(N, (t) => (t >= nc && t < nc + n ? p90e[t - nc] : 0)), rate: r }).lcoe;
    const del = inp.txCapex > 0 || inp.txOpex > 0 ? lcoeDelivered({ lcoeGeneration: L.lcoe, annualEnergy: eMean / (1 - inp.lossesFrac / 100), lossesFrac: inp.lossesFrac / 100, transmissionCapex: inp.txCapex, transmissionOmPerYear: inp.txOpex, rate: r, life: n }) : null;
    const av = lace({ energy, marginalPrice: inp.merchant, capacityMW: inp.ratedMW, capacityCredit: inp.capacityCredit / 100, capacityPaymentPerMWyr: inp.capacityPayment, rate: r, lcoe: L.lcoe });
    const rev = run(() => blendedRevenue({ energy: e, merchant: inp.merchant, ppa: inp.ppaShare > 0 && inp.ppaTerm > 0 ? { share: inp.ppaShare / 100, price: inp.ppaPrice, escalation: inp.ppaEscalation / 100, termYears: Math.round(inp.ppaTerm) } : undefined,
      cfd: inp.cfdShare > 0 && inp.cfdTerm > 0 ? { share: inp.cfdShare / 100, strike: inp.cfdStrike, termYears: Math.round(inp.cfdTerm), reference: inp.merchant } : undefined, rate: r }), 'Revenue');
    const st = inp.storage ? run(() => lcos({ capex: inp.stCapex, omPerYear: inp.stOpex, capacity_MWh: inp.stEnergy, cyclesPerYear: inp.stCycles, roundTripEfficiency: inp.stEff / 100, chargePrice: L.lcoe, degradation: 0.01, life: Math.round(inp.stLife), rate: r }), 'LCOS') : null;
    const margin = rev.levelisedPrice - L.lcoe;
    const ncf = eMean / (inp.ratedMW * 8760);
    const avoided = eMean * inp.gridIntensity / 1e3; // t CO2 per year
    // sensitivity map: LCOE against capex and energy multipliers
    const fx = [0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3]; const fy = [0.8, 0.9, 1, 1.1, 1.2];
    const grid = fy.map((ey) => fx.map((cx) => lcoe({ capex: capex.map((v) => v * cx), opex: range(N, (t) => (t >= nc && t < nc + n ? inp.opexFixed + inp.opexVar * e[t - nc] * ey : 0)), decom, energy: energy.map((v) => v * ey), rate: r }).lcoe));
    const rates = [3, 4, 5, 6, 7, 8, 9, 10].filter((x) => x > 0);
    const byRate = rates.map((x) => lcoe({ capex, opex, decom, energy, rate: x / 100 }).lcoe);

    const kpis = [
      { key: 'lcoe', label: 'LCOE (real, pre-tax)', value: round(L.lcoe), unit: PMWH, limit: round(rev.levelisedPrice), sense: 'max', explain: `Discounted lifetime cost over discounted lifetime energy in constant ${inp.baseYear} money at a ${inp.discountRate} % real rate: the constant real price that just recovers all costs. This is the IEA/IRENA/NREL convention.` },
      { key: 'lcoeNominal', label: 'LCOE (nominal)', value: round(rn.nominal), unit: PMWH, explain: `Flat money-of-the-day price with the same present value, at the nominal rate of ${pct(rn.nominalRate)}. Always higher than the real LCOE under positive inflation; never compare a nominal LCOE with a real one.` },
      { key: 'lcoeP90', label: 'LCOE at P90 energy (exceedance)', value: round(L90), unit: PMWH, explain: 'LCOE if yield is the P90 exceedance case — the energy exceeded with 90 % probability, the low-yield case lenders use. Note the convention: exceedance P90 equals the cumulative P10 of energy.' },
      { key: 'lcoeFcr', label: 'LCOE by fixed-charge rate (cross-check)', value: round(fcr.lcoe), unit: PMWH, explain: 'Capital recovery factor × overnight capex plus O&M over mean annual energy. It ignores construction time, replacement and decommissioning, so it sits below the discounted figure.' },
      { key: 'capexShare', label: 'Capital share of LCOE', value: round(100 * L.shares.capex), unit: '%', explain: 'Fraction of the levelised cost that is capital. A high share makes LCOE sensitive to the discount rate and to yield rather than to operating cost.' },
      { key: 'ncf', label: 'Net capacity factor', value: round(100 * ncf), unit: '%', explain: 'Mean delivered energy over rated capacity × 8760 h, after export losses, curtailment and degradation.' },
      { key: 'price', label: 'Levelised realised price', value: round(rev.levelisedPrice), unit: PMWH, explain: 'Discounted revenue over discounted energy from the PPA, CfD and merchant shares: directly comparable with the real LCOE.' },
      { key: 'margin', label: 'Margin over LCOE', value: round(margin), unit: PMWH, limit: 0, sense: 'min', explain: 'Levelised realised price less LCOE. Positive means the project earns more than the discount rate before tax.' },
      { key: 'contracted', label: 'Contracted share of revenue (PV)', value: round(100 * rev.contractedShareOfPv), unit: '%', explain: 'Share of discounted revenue coming from fixed-price contracts rather than the market; lenders size debt mainly on this part.' },
      { key: 'breakEven', label: 'Break-even tariff', value: round(L.lcoe), unit: PMWH, explain: 'Constant real, inflation-indexed tariff on all output at which pre-tax NPV is zero at the discount rate; by definition equal to the real LCOE.' },
      { key: 'lace', label: 'Levelised avoided cost of energy', value: round(av.lace), unit: PMWH, explain: 'Value to the system of the energy and firm capacity displaced, per MWh. A project is economically attractive to the system when this exceeds its LCOE.' },
      { key: 'netValue', label: 'Net system value (LACE − LCOE)', value: round(av.netValue), unit: PMWH, explain: 'Avoided cost less levelised cost. Positive: the project is worth more to the system than it costs.' },
      { key: 'avoidedCO2', label: 'Emissions displaced', value: round(avoided / 1e3), unit: 'kt CO₂/yr', explain: 'Mean annual energy times the emission intensity of the generation displaced; an indicator for carbon revenue and taxonomy reporting, not a life-cycle footprint.' },
    ];
    if (del) kpis.push({ key: 'lcode', label: 'Levelised cost of delivered energy', value: round(del.lcode), unit: PMWH, explain: 'Generation LCOE grossed up for transmission losses plus the levelised cost of the transmission asset, per MWh delivered at the grid connection.' });
    if (st) kpis.push({ key: 'lcos', label: 'Levelised cost of storage', value: round(st.lcos), unit: PMWH, explain: 'Discounted storage capital, O&M and charging cost per MWh discharged, charging at the project LCOE. Storage pays only if the price spread it captures exceeds this less the charging cost.' });
    const checks = [check('LCOE ≤ levelised realised price', L.lcoe, rev.levelisedPrice, PMWH, 'IEA/NEA Projected Costs of Generating Electricity (LCOE as break-even price)', 'Guards against a project whose lifetime cost per MWh exceeds what it will be paid; above 1.0 it does not earn the discount rate.'),
      check('LCOE at P90 energy ≤ realised price', L90, rev.levelisedPrice, PMWH, 'Lender downside case on P90 (exceedance) energy yield', 'Guards against a project that only works at the median yield; above 1.0 a low-wind outcome would not cover costs at the required return.')];
    const comp = [['Capital', L.components.capex], ['Operation and maintenance', L.components.opex], ['Decommissioning', L.components.decom]];
    const yrs = range(n, (k) => inp.baseYear + nc + k);
    const plots = [
      { type: 'bar', title: 'LCOE components', categories: comp.map((c) => c[0]), x: { label: 'Component' }, y: { label: 'Contribution', unit: PMWH }, series: [{ name: 'LCOE', y: comp.map((c) => round(c[1])) }],
        explain: 'Contribution of each cost block to the LCOE. Capital dominates for wind, wave and solar; decommissioning is small because it is discounted over the whole life.' },
      { type: 'heatmap', title: 'LCOE against capital cost and energy yield', xs: fx.map((v) => Math.round(v * 100)), ys: fy.map((v) => Math.round(v * 100)), z: grid.map((row) => row.map((v) => round(v))), x: { label: 'CAPEX', unit: '% of base' }, y: { label: 'Energy', unit: '% of base' }, zLabel: 'LCOE', zUnit: PMWH,
        explain: 'LCOE for combinations of capital cost and energy yield. Read along a row for the cost of a capital overrun and down a column for the cost of a yield shortfall; the two are nearly interchangeable per percent.' },
      { type: 'line', title: 'LCOE against discount rate', x: { label: 'Real discount rate', unit: '%/yr' }, y: { label: 'LCOE', unit: PMWH }, series: [{ name: 'LCOE', x: rates, y: byRate.map((v) => round(v)) }, { name: 'Realised price', x: [rates[0], rates[rates.length - 1]], y: [round(rev.levelisedPrice), round(rev.levelisedPrice)], dash: true }],
        explain: 'How the cost of capital moves the LCOE. The rate at which the curve crosses the realised price is the pre-tax real return the project earns.' },
      { type: 'line', title: 'Annual energy and revenue', x: { label: 'Year', unit: '' }, y: { label: 'Energy (GWh) / revenue (M)', unit: '' },
        series: [{ name: 'P50 energy (GWh)', x: yrs, y: e.map((v) => round(v / 1e3)) }, { name: 'P90 exceedance energy (GWh)', x: yrs, y: p90e.map((v) => round(v / 1e3)), dash: true }, { name: `Revenue (${MU})`, x: yrs, y: rev.revenue.map(mil) }],
        explain: 'Degrading energy output with its P50 and P90 (exceedance) cases, and annual revenue. The step in revenue marks the end of a contract term and the return to merchant prices.' },
    ];
    const tables = [
      { title: 'Annual energy and revenue by stream', columns: [{ key: 'year', label: 'Year' }, { key: 'e', label: 'Energy P50', unit: 'GWh' }, { key: 'e90', label: 'Energy P90 (exceedance)', unit: 'GWh' }, { key: 'ppa', label: 'PPA', unit: MU }, { key: 'cfd', label: 'CfD', unit: MU }, { key: 'mer', label: 'Merchant', unit: MU }, { key: 'rev', label: 'Revenue', unit: MU }],
        rows: yrs.map((y, k) => ({ year: y, e: round(e[k] / 1e3), e90: round(p90e[k] / 1e3), ppa: mil(rev.byStream.ppa[k]), cfd: mil(rev.byStream.cfd[k]), mer: mil(rev.byStream.merchant[k]), rev: mil(rev.revenue[k]) })),
        explain: 'Yearly energy and the revenue from each contract stream; fixed contractual prices are kept separate from the externally varying market price.' },
      { title: 'Levelised cost summary', columns: [{ key: 'item', label: 'Measure' }, { key: 'v', label: 'Value', unit: PMWH }, { key: 'note', label: 'Convention' }],
        rows: [{ item: 'LCOE real', v: round(L.lcoe), note: `constant ${inp.baseYear} money, real rate ${inp.discountRate} %, pre-tax` }, { item: 'LCOE nominal', v: round(rn.nominal), note: `money of the day, nominal rate ${pct(rn.nominalRate)}` }, { item: 'LCOE at P90 energy', v: round(L90), note: 'exceedance P90 yield' },
          ...(del ? [{ item: 'Delivered energy (LCODE)', v: round(del.lcode), note: `generation ${fmt(del.generation)} + transmission ${fmt(del.transmission)}` }] : []), ...(st ? [{ item: 'Storage (LCOS)', v: round(st.lcos), note: 'per MWh discharged, charging at LCOE' }] : []),
          { item: 'Avoided cost (LACE)', v: round(av.lace), note: `energy ${fmt(av.energyValue)} + capacity ${fmt(av.capacityValue)}` }, { item: 'Levelised realised price', v: round(rev.levelisedPrice), note: 'PPA + CfD + merchant, real' }],
        explain: 'All levelised measures side by side with the convention each uses, so that like is compared with like.' },
    ];
    const dCapex = L.lcoe * L.shares.capex * 0.1; const dE = L.lcoe * (1 / 0.99 - 1);
    const advice = [
      adv('engineering', 'info', `Capital is ${pct(L.shares.capex)} of the LCOE and the net capacity factor is ${pct(ncf)}. A 10 % CAPEX change moves LCOE by ${fmt(dCapex)} ${PMWH}; each 1 % of energy lost (wake, availability, cable losses, curtailment) adds ${fmt(dE)} ${PMWH}.`,
        `Weigh design choices on LCOE, not CAPEX: a measure costing ${fM(0.01 * inp.capex)} (1 % of CAPEX) pays if it raises lifetime energy by more than ${pct(0.01 * L.shares.capex)}. Export losses of ${inp.lossesFrac} % alone cost ${fmt(L.lcoe * (1 / (1 - inp.lossesFrac / 100) - 1))} ${PMWH}. Extending the life from ${n} to ${n + 5} years lowers LCOE by ${fmt(L.lcoe - lcoe({ capex: [...capex.slice(0, nc + n), 0, 0, 0, 0, 0, 0], opex: range(N + 5, (t) => (t >= nc && t < nc + n + 5 ? inp.opexFixed + inp.opexVar * eMean : 0)), decom: range(N + 5, (t) => (t === nc + n + 5 ? inp.decom : 0)), energy: range(N + 5, (t) => (t >= nc && t < nc + n + 5 ? e[Math.min(t - nc, n - 1)] : 0)), rate: r }).lcoe)} ${PMWH} if the fatigue life allows it.`,
        'LCOE as the design objective, IEA Wind TCP Task 26 and NREL/TP-462-5173.',
        `The project displaces about ${fmt(avoided / 1e3)} kt CO₂ a year at ${inp.gridIntensity} g/kWh (${fmt(avoided * n / 1e6)} Mt over the life). Its own embodied emissions are typically 6–15 g CO₂e/kWh for offshore wind, well below the 100 g CO₂e/kWh life-cycle threshold used by the EU taxonomy for electricity generation; confirm with the techno-economic module, where steel at ${CO2_STEEL} t CO₂e/t dominates.`),
      adv('economic', margin >= 0 ? (L90 <= rev.levelisedPrice ? 'ok' : 'warn') : 'fail', `Real LCOE ${fmt(L.lcoe)} ${PMWH} (nominal ${fmt(rn.nominal)}); levelised realised price ${fmt(rev.levelisedPrice)} ${PMWH}, margin ${fmt(margin)}; ${pct(rev.contractedShareOfPv)} of revenue is contracted. LACE ${fmt(av.lace)} ${PMWH}.`,
        margin < 0 ? `The break-even tariff is ${fmt(L.lcoe)} ${PMWH}: the realised price is ${fmt(-margin)} short. Raise the contracted price, cut cost or accept a lower return.` : L90 > rev.levelisedPrice ? `The margin disappears at P90 (exceedance) energy (LCOE ${fmt(L90)} ${PMWH}): reduce yield uncertainty with a longer measurement campaign or hold a debt-service reserve.` : `The margin survives the P90 energy case. Extending contract cover beyond year ${Math.max(inp.ppaTerm, inp.cfdTerm)} would cut the merchant tail exposure.`,
        'LCOE conventions of IEA/NEA Projected Costs of Generating Electricity (2020) and IRENA (real, pre-tax, discounted); P50/P90 energy as exceedance probabilities in lender practice.'),
    ];
    return { kpis, checks, plots, tables, advice,
      notes: [`Real terms: constant ${inp.baseYear} money and a real discount rate of ${inp.discountRate} %; the nominal LCOE uses ${pct(rn.nominalRate)} with ${inp.inflation} % inflation.`, `Percentile convention for energy: exceedance — ${gen.convention} The cost percentiles elsewhere in the suite are cumulative.`, `Timing: capex over rows 0–${nc - 1}, operation for ${n} years from row ${nc}, decommissioning in row ${nc + n}.`],
      publish: { 'econ.lcoe': L.lcoe, 'econ.lcoeNominal': rn.nominal, 'econ.lace': av.lace, 'econ.breakEvenTariff': L.lcoe } };
  },
  verify: [
    // one construction year (capex at t = 0), 20 operating years, no degradation or losses: LCOE = (1000 × (A/P, 8 %, 20) + 50)/100 GWh
    //   = (1000 × 0.1018522 + 50) M / 100 000 MWh = 1518.52 per MWh; fixed-charge-rate form identical
    { name: 'Level costs and energy: LCOE = (capex × CRF + O&M) / E', source: 'NREL/TP-462-5173 §4 (fixed-charge-rate identity) with (A/P, 8 %, 20) = 0.1018522',
      inputs: { capex: 1e9, constructionYears: 1, opLife: 20, discountRate: 8, aep: 100, lossesFrac: 0, curtailment: 0, degradation: 0, opexFixed: 5e7, opexVar: 0, replCost: 0, decom: 0, storage: false, txCapex: 0, txOpex: 0 },
      expect: [{ kpi: 'lcoe', value: 1518.522, tol: 1e-5 }, { kpi: 'lcoeFcr', value: 1518.522, tol: 1e-5 }, { kpi: 'capexShare', value: 100 * 101.8522 / 151.8522, tol: 1e-4 }] },
    // delivered energy: LCODE = LCOE/(1 − l) + (tx capex × CRF + tx O&M)/(E_gen (1 − l)); with l = 0: 1518.522 + (200 × 0.1018522 + 2) M / 100 000 MWh = 1518.522 + 223.704
    { name: 'Delivered-energy adder of a transmission asset', source: 'Hand calculation: (200 M × 0.1018522 + 2 M)/100 000 MWh = 223.704 per MWh added to the generation LCOE',
      inputs: { capex: 1e9, constructionYears: 1, opLife: 20, discountRate: 8, aep: 100, lossesFrac: 0, curtailment: 0, degradation: 0, opexFixed: 5e7, opexVar: 0, replCost: 0, decom: 0, storage: false, txCapex: 2e8, txOpex: 2e6 },
      expect: [{ kpi: 'lcode', value: 1742.226, tol: 1e-5 }] },
    // revenue: 50 % PPA at 50 and 50 % merchant at 40, flat → levelised price 45
    { name: 'Blended price of a half-contracted project', source: 'Hand calculation: 0.5 × 50 + 0.5 × 40 = 45 per MWh', inputs: { ppaShare: 50, ppaPrice: 50, ppaEscalation: 0, ppaTerm: 60, cfdShare: 0, merchant: 40 }, expect: [{ kpi: 'price', value: 45, tol: 1e-6 }] },
  ],
};
