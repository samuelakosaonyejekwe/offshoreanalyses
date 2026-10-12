// Hydrocarbon production and revenue: plateau and Arps decline profile, associated gas and water, price
// scenarios with probabilities, royalties and tariffs, economic limit, unit technical cost and break-even price.
import { plateauDeclineProfile, associatedStreams, gasSales, arpsFit, arpsRate, arpsEur, declineConversion } from '../../engine/econ/production.js';
import { levelisedCostOfService } from '../../engine/econ/lcoe.js';
import { expectedMonetaryValue } from '../../engine/econ/risk.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, need, run, round, fmt, mil, fM, pct, adv, check, checkMin, cost, tableRows, uniqueNames, seriesXY, range, pv, linspace } from './_common.js';

const DAYS = 365.25;
const SCEN = [{ name: 'Low', oil: 52, gas: 2.6, probability: 25 }, { name: 'Base', oil: 72, gas: 3.6, probability: 50 }, { name: 'High', oil: 95, gas: 5.0, probability: 25 }];
const HIST = { x: [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4], y: [24800, 22300, 20400, 18600, 17300, 15900, 14900, 13800, 13000] };

export default {
  id: 's17-production-revenue', suite: 17, title: 'Production forecast and revenue (decline-curve based)',
  summary: 'Builds the oil, gas and water production profile from a plateau and an Arps decline, prices it under probability-weighted scenarios with royalties and tariffs, and reports revenue, cash flow, economic limit, unit technical cost and break-even price.',
  covers: [4], fidelity: 'analytical',
  standards: ['Arps, J.J. (1945), Analysis of decline curves, Trans. AIME 160', 'SPE Petroleum Resources Management System (PRMS, 2018) §3.1 (economic limit and forecasts)', 'SPEE Recommended Evaluation Practice (discounting and price cases)', 'ISO 15663:2021 (life-cycle costing in the petroleum industry)'],
  models: ['Production decline curve analysis', 'Arps exponential decline', 'Arps hyperbolic decline', 'Arps harmonic decline', 'Plateau and ramp-up production profile', 'Availability-adjusted production and production-loss cost', 'Hydrocarbon price scenarios', 'Expected monetary value over price scenarios', 'Unit technical cost (levelised cost per barrel of oil equivalent)', 'Economic limit and break-even price'],
  limits: `Single-segment field profile with constant uptime; water cut rises linearly with cumulative recovery; pre-tax cash flow with royalty only (corporate tax and financing are in the project-finance module). Reservoir simulation, well scheduling and facility constraints are outside the module. Prices are real (constant base-year money) and discounted at a real rate. ${CUR_NOTE}`,
  publishes: ['econ.eurMMboe', 'econ.revenueAnnual', 'econ.breakEvenPrice', 'econ.unitTechnicalCost'],
  searchTerms: 'Arps decline curve production forecast economic limit unit technical cost',
  inputs: [
    { key: 'baseYear', label: 'Price base year (first production the year after)', unit: 'yr', type: 'int', default: 2026, min: 1990, max: 2100, group: 'Opening state', role: 'ic' },
    { key: 'plateau', label: 'Initial (plateau) oil rate', unit: 'bbl/d', type: 'number', default: 25000, min: 1, max: 5e6, group: 'Opening state', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Initial production capability from the reservoir and flow-assurance studies.' },
    { key: 'rampYears', label: 'Ramp-up period', unit: 'yr', type: 'number', default: 1, min: 0, max: 10, group: 'Production information', role: 'operation' },
    { key: 'plateauYears', label: 'Plateau duration', unit: 'yr', type: 'number', default: 3, min: 0, max: 40, group: 'Production information', role: 'operation', design: true },
    { key: 'decline', label: 'Nominal initial decline rate', unit: '%/yr', type: 'number', default: 18, min: 0.1, max: 200, group: 'Production information', role: 'operation', uncertain: { dist: 'normal', cov: 0.2 } },
    { key: 'bExp', label: 'Arps exponent b', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Production information', role: 'operation', help: '0 = exponential, 1 = harmonic, between = hyperbolic.' },
    { key: 'years', label: 'Field life considered', unit: 'yr', type: 'int', default: 20, min: 2, max: 60, group: 'Production information', role: 'bc', bind: 'life.designYears', help: 'Operational constraint: facility design life from the engineering suites.' },
    { key: 'uptime', label: 'Production availability (uptime)', unit: '%', type: 'number', default: 94, min: 10, max: 100, group: 'Production information', role: 'bc', bind: 'energy.availability', design: true, help: 'From the reliability and maintenance analysis; applied to all volumes.' },
    { key: 'wc0', label: 'Initial water cut', unit: '%', type: 'number', default: 5, min: 0, max: 98, group: 'Production information', role: 'ic' },
    { key: 'wcEnd', label: 'Water cut at end of life', unit: '%', type: 'number', default: 85, min: 0, max: 98, group: 'Production information', role: 'operation' },
    { key: 'gor', label: 'Gas–oil ratio', unit: 'scf/bbl', type: 'number', default: 650, min: 0, max: 1e5, group: 'Production information', role: 'operation' },
    { key: 'fuelFlare', label: 'Fuel and flare share of gas', unit: '%', type: 'number', default: 12, min: 0, max: 100, group: 'Production information', role: 'operation' },
    { key: 'history', label: 'Measured rate history (optional fit)', type: 'series', group: 'Decline-curve fit', role: 'ic', xLabel: 'Time on decline', yLabel: 'Oil rate', xUnit: 'yr', yUnit: 'bbl/d', help: 'Historical production of this or an analogue field. The fitted decline is reported beside the entered one; it replaces it only when selected.', default: HIST },
    { key: 'useFit', label: 'Use the fitted decline parameters', type: 'bool', default: false, group: 'Decline-curve fit', role: 'control' },
    { key: 'scenarios', label: 'Price scenarios', type: 'table', group: 'Product prices', role: 'cost', help: `Real prices and scenario probabilities (sum 100 %). ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Scenario', type: 'text' }, { key: 'oil', label: 'Oil price', unit: `${CUR}/bbl`, type: 'number' }, { key: 'gas', label: 'Gas price', unit: `${CUR}/Mscf`, type: 'number' }, { key: 'probability', label: 'Probability', unit: '%', type: 'number' }], default: SCEN },
    cost({ key: 'tariff', label: 'Transport and processing tariff', unit: `${CUR}/boe`, default: 4.5, min: 0, max: 200, group: 'Product prices' }),
    { key: 'royalty', label: 'Royalty on gross revenue', unit: '%', type: 'number', default: 10, min: 0, max: 90, group: 'Fiscal regime', role: 'bc' },
    cost({ key: 'capex', label: 'Development capital (spent in the base year)', unit: CUR, default: 1.05e9, min: 0, max: 1e12, group: 'Costs', bind: 'econ.capex' }),
    cost({ key: 'opexFixed', label: 'Fixed operating cost', unit: `${CUR}/yr`, default: 6.5e7, min: 0, max: 1e11, group: 'Costs', bind: 'econ.opexAnnual' }),
    cost({ key: 'opexVar', label: 'Variable operating cost', unit: `${CUR}/boe`, default: 6, min: 0, max: 500, group: 'Costs' }),
    cost({ key: 'decom', label: 'Abandonment cost', unit: CUR, default: 1.8e8, min: 0, max: 1e12, group: 'Costs', bind: 'econ.decomCost', help: 'Terminal condition: spent in the year after cessation.' }),
    { key: 'discountRate', label: 'Real discount rate', unit: '%/yr', type: 'number', default: 10, min: 0, max: 50, group: 'Financial constraints', role: 'bc', help: `Hurdle rate in real terms, matching the real prices. ${PLACEHOLDER}` },
    { key: 'maxBreakEven', label: 'Maximum acceptable break-even oil price', unit: `${CUR}/bbl`, type: 'number', default: 55, min: 1, max: 500, group: 'Acceptance criteria', role: 'criteria' },
  ],
  compute(inp) {
    const sc = tableRows(inp.scenarios, 'Price scenarios', { oil: { min: 0 }, gas: { min: 0 }, probability: { min: 0, max: 100 } });
    uniqueNames(sc, 'Price scenarios');
    const pSum = sc.reduce((s, x) => s + x.probability, 0);
    need(Math.abs(pSum - 100) < 0.01, `Scenario probabilities must sum to 100 % (they sum to ${fmt(pSum, 5)} %): the scenarios are meant to be exclusive and exhaustive, not simultaneous extremes.`);
    need(inp.wcEnd >= inp.wc0, 'The end-of-life water cut cannot be below the initial water cut.');
    const n = Math.round(inp.years); const r = inp.discountRate / 100; const up = inp.uptime / 100;

    // optional decline fit to the measured history
    let fit = null;
    if (inp.history && (inp.history.x || []).length >= 3) {
      const h = seriesXY(inp.history, 'Rate history');
      need(h.y.every((v) => v > 0), 'Rate history: rates must be positive.');
      fit = run(() => arpsFit({ t: h.x, q: h.y, model: 'hyperbolic', bMax: 1 }), 'Decline-curve fit');
      fit.t = h.x; fit.q = h.y;
    }
    need(!inp.useFit || fit, 'No rate history is available to fit: enter at least three points or switch off the fitted decline.');
    const Di = inp.useFit ? fit.Di : inp.decline / 100; const b = inp.useFit ? fit.b : inp.bExp;
    const prof = run(() => plateauDeclineProfile({ plateauRate: inp.plateau * DAYS, rampYears: inp.rampYears, plateauYears: inp.plateauYears, Di, b, years: n, uptime: up }), 'Production profile');
    const ideal = run(() => plateauDeclineProfile({ plateauRate: inp.plateau * DAYS, rampYears: inp.rampYears, plateauYears: inp.plateauYears, Di, b, years: n, uptime: 1 }), 'Production profile');
    const oil = prof.annual; // bbl per year
    need(prof.total > 0, 'The production profile is empty.');
    // water cut rising linearly with recovery
    let cum = 0; const wc = oil.map((q) => { const f = (cum + 0.5 * q) / prof.total; cum += q; return (inp.wc0 + (inp.wcEnd - inp.wc0) * f) / 100; });
    const st = associatedStreams({ oil, waterCut: wc, gor: inp.gor / 1000 }); // gas in Mscf
    const gs = gasSales({ gross: st.gas, fuelFraction: inp.fuelFlare / 100 });
    const boe = oil.map((q, k) => q + gs.sales[k] / 6); // 6 Mscf per boe (energy-equivalence convention)
    const eurBoe = boe.reduce((a, c) => a + c, 0);
    const conv = declineConversion({ Di, b });
    const tail = arpsEur({ qi: inp.plateau * DAYS, Di, b, qLimit: Math.max(inp.plateau * DAYS * 0.02, 1) });

    // cash flow for one price case; rows 0..n+1 (capex at 0, production years 1..n, abandonment the year after cessation)
    const flow = (po, pg) => {
      const net = new Float64Array(n + 2); const rev = new Float64Array(n + 2); const op = new Float64Array(n + 2);
      net[0] = -inp.capex;
      let last = 0;
      for (let k = 0; k < n; k++) {
        rev[k + 1] = oil[k] * po + gs.sales[k] * pg;
        op[k + 1] = inp.opexFixed + (inp.opexVar + inp.tariff) * boe[k];
        net[k + 1] = rev[k + 1] * (1 - inp.royalty / 100) - op[k + 1];
        if (net[k + 1] > 0) last = k + 1;
      }
      // economic limit: stop after the last year with positive operating cash flow, then abandon
      const trunc = new Float64Array(n + 2); trunc[0] = net[0];
      for (let t = 1; t <= last; t++) trunc[t] = net[t];
      trunc[Math.min(last + 1, n + 1)] -= inp.decom;
      return { net: trunc, rev, op, last, npv: pv(trunc, r), opcf: net };
    };
    const cases = sc.map((s) => ({ ...s, ...flow(s.oil, s.gas) }));
    const emv = expectedMonetaryValue(cases.map((c) => ({ p: c.probability / 100, value: c.npv })));
    const base = cases.reduce((a, c) => (c.probability > a.probability ? c : a), cases[0]);
    // break-even oil price at the base gas price (NPV is piecewise linear in price because of the economic limit: bisection)
    let lo = 0; let hi = 1000;
    need(flow(hi, base.gas).npv > 0, 'No oil price up to 1000 per barrel makes the development break even: costs exceed any plausible revenue.');
    for (let it = 0; it < 60; it++) { const mid = 0.5 * (lo + hi); if (flow(mid, base.gas).npv > 0) hi = mid; else lo = mid; }
    const breakEven = 0.5 * (lo + hi);
    // unit technical cost: PV of capex + opex + tariff + abandonment over PV of boe, to the base-case economic limit
    const costs = range(n + 2, (t) => (t === 0 ? inp.capex : t <= base.last ? base.op[t] : 0)); costs[Math.min(base.last + 1, n + 1)] += inp.decom;
    const qty = range(n + 2, (t) => (t >= 1 && t <= base.last ? boe[t - 1] : 0));
    const utc = qty.some((v) => v > 0) ? levelisedCostOfService({ costs, quantities: qty, rate: r }).levelisedCost : 0;
    const lostOil = ideal.annual.map((q, k) => q - oil[k]);
    const downtimeCost1 = lostOil[Math.min(Math.ceil(inp.rampYears), n - 1)] * base.oil * (1 - inp.royalty / 100);
    const downtimePV = pv([0, ...lostOil.map((q) => q * (base.oil * (1 - inp.royalty / 100) - inp.opexVar - inp.tariff))], r);
    const peakRev = Math.max(...base.rev);

    const kpis = [
      { key: 'eurOil', label: `Oil recovered in ${n} years`, value: round(prof.total / 1e6), unit: 'MMbbl', explain: 'Cumulative oil over the field life considered, after uptime. Volumes beyond the economic limit of a price case are not sold in that case.' },
      { key: 'eurBoe', label: 'Oil-equivalent recovery', value: round(eurBoe / 1e6), unit: 'MMboe', explain: 'Oil plus sales gas at 6 Mscf per barrel of oil equivalent (energy-equivalence convention, not value equivalence).' },
      { key: 'peakRate', label: 'Peak annual average oil rate', value: round(Math.max(...oil) / DAYS), unit: 'bbl/d', explain: 'Highest yearly average rate after uptime; it sizes the processing and export capacity that must be paid for.' },
      { key: 'effDecline', label: 'Effective first-year decline', value: round(100 * conv.De), unit: '%/yr', explain: 'Fraction of rate lost over the first year of decline, the figure usually quoted, derived from the nominal decline and the exponent b.' },
      { key: 'revenue', label: `Peak annual revenue (${base.name} case)`, value: mil(peakRev), unit: `${MU}/yr`, explain: 'Gross oil and gas sales revenue in the best production year at the most probable price case, before royalty.' },
      { key: 'npv', label: `NPV, ${base.name} price case (real, pre-tax)`, value: mil(base.npv), unit: MU, limit: 0, sense: 'min', explain: `Net present value at ${inp.baseYear} of revenue after royalty less operating cost, tariff, capital and abandonment, at the real discount rate, produced to the economic limit.` },
      { key: 'emv', label: 'Expected NPV over the price scenarios', value: mil(emv.emv), unit: MU, explain: 'Probability-weighted NPV of the scenarios (expected monetary value). The scenarios are alternatives with stated probabilities, not simultaneous extremes.' },
      { key: 'npvStd', label: 'Standard deviation of NPV over scenarios', value: mil(emv.std), unit: MU, explain: 'Spread of NPV caused by price alone; compare it with the expected NPV to judge how price-dependent the decision is.' },
      { key: 'breakEven', label: 'Break-even oil price', value: round(breakEven), unit: `${CUR}/bbl`, limit: inp.maxBreakEven, sense: 'max', explain: 'Constant real oil price at which NPV is zero at the hurdle rate, with gas at its base-case price. The lower it is, the more price falls the development survives.' },
      { key: 'utc', label: 'Unit technical cost', value: round(utc), unit: `${CUR}/boe`, explain: 'Discounted capital, operating, tariff and abandonment cost per discounted barrel of oil equivalent: the levelised cost of the produced service, before royalty.' },
      { key: 'econLimit', label: `Economic limit (${base.name} case)`, value: base.last, unit: 'yr', explain: 'Last production year with positive operating cash flow. Producing beyond it loses money, so cessation and abandonment follow.' },
      { key: 'downtime', label: 'Revenue lost to downtime, plateau year', value: mil(downtimeCost1), unit: `${MU}/yr`, explain: `Net revenue of the production not achieved because uptime is ${inp.uptime} % rather than 100 %. Part is deferred rather than lost in a real reservoir, so this is an upper bound.` },
      { key: 'downtimePV', label: 'Present value of downtime losses', value: mil(downtimePV), unit: MU, explain: 'Discounted net margin of all production lost to downtime over the life: the most that reliability improvements could be worth.' },
    ];
    const checks = [
      check('Break-even price ≤ acceptance limit', breakEven, inp.maxBreakEven, `${CUR}/bbl`, 'Investment criterion (price resilience)', 'Guards against sanctioning a development that only works at high prices; above 1.0 the project fails the owner’s low-price test.'),
      checkMin('Economic limit covers the field life', base.last, n, 'yr', 'SPE PRMS (2018) §3.1.2 (economic limit)', 'Guards against counting production that would be sold at a loss. Above 1.0 the field stops paying before the design life is reached and the late volumes are not reserves.'),
    ];
    const yrs = range(n, (k) => inp.baseYear + 1 + k);
    const plots = [
      { type: 'area', title: 'Production profile', x: { label: 'Year', unit: '' }, y: { label: 'Rate', unit: 'bbl/d or boe/d' },
        series: [{ name: 'Oil', x: yrs, y: oil.map((q) => round(q / DAYS)) }, { name: 'Sales gas (boe)', x: yrs, y: gs.sales.map((q) => round(q / 6 / DAYS)) }, { name: 'Water', x: yrs, y: st.water.map((q) => round(q / DAYS)) }],
        explain: 'Yearly average oil, sales gas and produced water. Watch the water curve: it loads the facilities and operating cost just as the oil that pays for them declines.' },
      { type: 'line', title: 'Net cash flow by price scenario', x: { label: 'Year', unit: '' }, y: { label: 'Net cash flow', unit: MU },
        series: cases.map((c) => ({ name: `${c.name} (${c.oil} ${CUR}/bbl)`, x: [inp.baseYear, ...yrs, inp.baseYear + n + 1], y: Array.from(c.net, mil) })),
        explain: 'Annual pre-tax net cash flow for each price case, truncated at its economic limit with abandonment the year after. The lines separate most in the plateau years.' },
      { type: 'bar', title: 'NPV by price scenario', categories: [...cases.map((c) => c.name), 'Expected'], x: { label: 'Scenario' }, y: { label: 'NPV', unit: MU }, series: [{ name: 'NPV', y: [...cases.map((c) => mil(c.npv)), mil(emv.emv)] }],
        explain: 'NPV in each price case and the probability-weighted expectation. A negative low case with a positive expectation is a price bet, not a robust project.' },
    ];
    if (fit) {
      const tt = linspace(0, Math.max(...fit.t) * 1.6, 50);
      plots.push({ type: 'line', title: 'Decline-curve fit to the rate history', x: { label: 'Time on decline', unit: 'yr' }, y: { label: 'Oil rate', unit: 'bbl/d' },
        series: [{ name: 'Measured', x: fit.t, y: fit.q }, { name: `Fitted (Di ${fmt(fit.Di * 100)} %/yr, b ${fmt(fit.b)})`, x: tt.map((v) => round(v)), y: tt.map((t) => round(arpsRate({ qi: fit.qi, Di: fit.Di, b: fit.b, t }))) },
          { name: `Entered (Di ${inp.decline} %/yr, b ${inp.bExp})`, x: tt.map((v) => round(v)), y: tt.map((t) => round(arpsRate({ qi: fit.qi, Di: inp.decline / 100, b: inp.bExp, t }))), dash: true }],
        explain: 'Measured rates with the least-squares Arps fit and the entered decline from the same starting rate. Divergence beyond the data shows how much the forecast depends on the exponent b.' });
    }
    const tables = [
      { title: `Annual forecast (${base.name} price case)`, columns: [{ key: 'year', label: 'Year' }, { key: 'oil', label: 'Oil', unit: 'bbl/d' }, { key: 'gas', label: 'Sales gas', unit: 'MMscf/d' }, { key: 'water', label: 'Water', unit: 'bbl/d' }, { key: 'wc', label: 'Water cut', unit: '%' }, { key: 'rev', label: 'Revenue', unit: MU }, { key: 'opex', label: 'Opex and tariff', unit: MU }, { key: 'net', label: 'Operating cash flow', unit: MU }],
        rows: yrs.map((y, k) => ({ year: y, oil: round(oil[k] / DAYS), gas: round(gs.sales[k] / 1000 / DAYS), water: round(st.water[k] / DAYS), wc: round(100 * wc[k], 3), rev: mil(base.rev[k + 1]), opex: mil(base.op[k + 1]), net: mil(base.opcf[k + 1]) })),
        explain: 'Yearly volumes, revenue and operating cash flow before truncation. The first year with a negative operating cash flow marks the economic limit.' },
      { title: 'Price scenarios', columns: [{ key: 'name', label: 'Scenario' }, { key: 'p', label: 'Probability', unit: '%' }, { key: 'oil', label: 'Oil', unit: `${CUR}/bbl` }, { key: 'gas', label: 'Gas', unit: `${CUR}/Mscf` }, { key: 'limit', label: 'Economic limit', unit: 'yr' }, { key: 'npv', label: 'NPV', unit: MU }],
        rows: cases.map((c) => ({ name: c.name, p: c.probability, oil: c.oil, gas: c.gas, limit: c.last, npv: mil(c.npv) })), explain: 'Each price case with its probability, the year production stops paying and its NPV; the expected value weights the last column by the probabilities.' },
    ];
    const low = cases.reduce((a, c) => (c.npv < a.npv ? c : a), cases[0]);
    const advice = [
      adv('economic', breakEven <= inp.maxBreakEven && emv.emv > 0 ? 'ok' : 'warn', `Expected NPV ${fM(emv.emv)} (${cases.map((c) => `${c.name} ${fM(c.npv)}`).join(', ')}); break-even ${fmt(breakEven)} ${CUR}/bbl against a limit of ${inp.maxBreakEven}; unit technical cost ${fmt(utc)} ${CUR}/boe.`,
        low.npv < 0 ? `The ${low.name} case loses ${fM(-low.npv)}: phase capital, secure a lower tariff or hedge so the project survives ${low.oil} ${CUR}/bbl, or accept the exposure explicitly.` : 'All price cases are positive: the development is robust to the stated price range. Test it next against schedule delay and capital overrun.',
        'Price-scenario evaluation and economic limit, SPE PRMS (2018) §3.1 and SPEE recommended evaluation practices; expected monetary value, Clemen & Reilly ch. 4.'),
      adv('engineering', 'info', `Uptime of ${inp.uptime} % costs ${fM(downtimeCost1)} of net revenue in a plateau year and ${fM(downtimePV)} in present value; the decline is ${fmt(100 * conv.De)} %/yr effective${fit ? ` (history fits Di ${fmt(100 * fit.Di)} %/yr, b ${fmt(fit.b)}, R² ${fmt(fit.r2)})` : ''}.`,
        `Each percentage point of availability is worth about ${fM(downtimePV / Math.max(100 - inp.uptime, 1e-9))} in present value: compare that with the cost of sparing, redundancy or faster intervention from the reliability module. Water handling must reach ${fmt(Math.max(...st.water) / DAYS)} bbl/d.`,
        'Decline-curve analysis, Arps (1945); production availability targets per ISO 20815:2018.',
        `Fuel and flare take ${inp.fuelFlare} % of produced gas — ${fmt(gs.fuel.reduce((a, c) => a + c, 0) / 1e6 + gs.flare.reduce((a, c) => a + c, 0) / 1e6)} Bscf over the life, roughly ${fmt((gs.fuel.reduce((a, c) => a + c, 0) + gs.flare.reduce((a, c) => a + c, 0)) * 0.0537 / 1e3)} kt CO₂ at about 53.7 kg CO₂ per Mscf burnt (screening value for pipeline-quality gas). Flare recovery and electrification reduce both the emission charge and the lost sales volume; late-life high-water-cut barrels carry the highest emissions per barrel.`),
    ];
    return { kpis, checks, plots, tables, advice,
      notes: [`Real terms: constant ${inp.baseYear} prices discounted at ${inp.discountRate} % real; pre-tax, after ${inp.royalty} % royalty.`, `Decline used: Di = ${fmt(100 * Di)} %/yr nominal, b = ${fmt(b)} (${inp.useFit ? 'fitted to the history' : 'as entered'}); unconstrained decline to 2 % of the initial rate would take ${fmt(tail.timeToLimit)} years.`, 'Scenario probabilities describe alternative price paths; prices are indicative placeholders.'],
      publish: { 'econ.eurMMboe': eurBoe / 1e6, 'econ.revenueAnnual': peakRev, 'econ.breakEvenPrice': breakEven, 'econ.unitTechnicalCost': utc } };
  },
  verify: [
    // exponential decline from time zero: Np(10 yr) = qi/Di (1 − e^(−Di t)) = 10 000 × 365.25 / 0.2 × (1 − e^−2) = 18.2625e6 × 0.864665 = 15.7909 MMbbl
    { name: 'Arps exponential cumulative, Di = 20 %/yr over 10 years', source: 'Arps (1945): Np = (qi − q)/Di = 10 000 × 365.25/0.2 × (1 − e⁻²) = 15.7909 MMbbl',
      inputs: { plateau: 10000, rampYears: 0, plateauYears: 0, decline: 20, bExp: 0, years: 10, uptime: 100, gor: 0, history: { x: [], y: [] } }, expect: [{ kpi: 'eurOil', value: 15.7909, tol: 1e-4 }, { kpi: 'effDecline', value: 18.1269, tol: 1e-4 }] },
    // hyperbolic b = 0.5: Np = qi/((1 − b) Di) [1 − (1 + b Di t)^−(1−b)/b] = 3.6525e6/0.1 × (1 − 1/2) = 18.2625 MMbbl
    { name: 'Arps hyperbolic cumulative, b = 0.5', source: 'Arps (1945): Np = qi^b (qi^(1−b) − q^(1−b))/((1 − b) Di) with q(10) = qi/4 → 18.2625 MMbbl',
      inputs: { plateau: 10000, rampYears: 0, plateauYears: 0, decline: 20, bExp: 0.5, years: 10, uptime: 100, gor: 0, history: { x: [], y: [] } }, expect: [{ kpi: 'eurOil', value: 18.2625, tol: 1e-4 }] },
    // harmonic b = 1: Np = (qi/Di) ln(1 + Di t) = 18.2625e6 × ln 3 = 20.0633 MMbbl
    { name: 'Arps harmonic cumulative, b = 1', source: 'Arps (1945): Np = (qi/Di) ln(qi/q) = 18.2625 × ln 3 = 20.0633 MMbbl',
      inputs: { plateau: 10000, rampYears: 0, plateauYears: 0, decline: 20, bExp: 1, years: 10, uptime: 100, gor: 0, history: { x: [], y: [] } }, expect: [{ kpi: 'eurOil', value: 20.0633, tol: 1e-4 }] },
  ],
};
