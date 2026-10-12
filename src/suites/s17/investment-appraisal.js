// Classical engineering-economy appraisal of a user-supplied cash-flow table: time-value factors, NPV, IRR with
// multiple/absent-root diagnostics, MIRR, profitability index, benefit–cost ratio, ROI, payback, annual worth.
import { npv, npvProfile, mirr, profitabilityIndex, benefitCostRatio, roi, payback, discountSeries } from '../../engine/econ/cashflow.js';
import { factorTable, effectiveRate, equivalentAnnualWorth, fvContinuous, realToNominal, nominalToReal } from '../../engine/econ/tvm.js';
import { CUR, MU, CUR_NOTE, CO2_STEEL, need, run, round, fmt, mil, fM, pct, adv, check, checkMin, tableRows, irrReport, linspace } from './_common.js';

const flow = (year, investment, benefit, cost, salvage = 0) => ({ year, investment, benefit, cost, salvage });
const DEFAULT_FLOWS = [flow(0, 6.0e7, 0, 0), flow(1, 6.0e7, 0, 0), ...Array.from({ length: 14 }, (_, i) => flow(i + 2, 0, 2.4e7, 6.5e6)), flow(15, 0, 2.4e7, 6.5e6, 8e6), flow(16, 0, 0, 1.4e7)];

export default {
  id: 's17-investment-appraisal', suite: 17, title: 'Investment appraisal of a cash-flow series',
  summary: 'Turns a table of yearly investments, benefits and costs into NPV, IRR (with multiple-root diagnostics), MIRR, profitability index, benefit–cost ratio, ROI, payback and equivalent annual worth.',
  covers: [1], fidelity: 'analytical',
  standards: ['Blank & Tarquin, Engineering Economy, 8th ed., ch. 2–9 and 13', 'ISO 15686-5:2017 §5 (discounting and NPV)', 'Norstrom (1972) and Descartes sign tests for rate-of-return uniqueness'],
  models: ['Time value of money and compound-interest factors', 'Continuous and periodic compounding', 'Discounted cash flow and net present value', 'Internal rate of return (all real roots)', 'Modified internal rate of return', 'Profitability index and benefit–cost ratio', 'Return on investment', 'Simple and discounted payback', 'Equivalent annual worth and capital recovery factor', 'Fisher relation (real versus nominal)'],
  limits: `End-of-year discrete cash flows on whole years (optional mid-year discounting); taxes, financing and inflation must already be inside the table and the table must be wholly real or wholly nominal, matching the rate. ${CUR_NOTE}`,
  publishes: ['econ.appraisalNpv', 'econ.appraisalIrr'],
  searchTerms: 'engineering economy net present value internal rate of return multiple roots',
  inputs: [
    { key: 'baseYear', label: 'Price base year (year 0)', unit: 'yr', type: 'int', default: 2026, min: 1990, max: 2100, group: 'Opening financial state', role: 'ic', help: 'Calendar year of time zero: flows in year 0 are not discounted.' },
    { key: 'basis', label: 'Terms of the table and the rate', type: 'select', default: 'real', group: 'Opening financial state', role: 'ic',
      options: [{ value: 'real', label: 'Real (constant money)' }, { value: 'nominal', label: 'Nominal (money of the day)' }], help: 'Table and discount rate must be in the same terms; the other rate is reported through the Fisher relation.' },
    { key: 'inflation', label: 'General inflation', unit: '%/yr', type: 'number', default: 2, min: -2, max: 30, group: 'Opening financial state', role: 'ic' },
    { key: 'flows', label: 'Cash-flow lines', type: 'table', group: 'Cash-flow timing', role: 'cost',
      help: `One row per year and line; rows of the same year are added. All amounts positive, in ${CUR}. Indicative placeholder — replace with current quotations or indices.`,
      columns: [{ key: 'year', label: 'Year', type: 'number' }, { key: 'investment', label: 'Investment', unit: CUR, type: 'number' }, { key: 'benefit', label: 'Revenue / benefit', unit: CUR, type: 'number' },
        { key: 'cost', label: 'Operating cost', unit: CUR, type: 'number' }, { key: 'salvage', label: 'Salvage / residual', unit: CUR, type: 'number' }],
      default: DEFAULT_FLOWS },
    { key: 'rate', label: 'Discount rate (hurdle rate, MARR)', unit: '%/yr', type: 'number', default: 8, min: -20, max: 100, group: 'Financial constraints', role: 'bc', uncertain: { dist: 'normal', std: 1 },
      help: 'Minimum attractive rate of return. Indicative placeholder — replace with current quotations or indices.' },
    { key: 'financeRate', label: 'Finance rate for MIRR', unit: '%/yr', type: 'number', default: 6, min: -20, max: 100, group: 'Financial constraints', role: 'bc', help: 'Rate at which outflows are funded.' },
    { key: 'reinvestRate', label: 'Reinvestment rate for MIRR', unit: '%/yr', type: 'number', default: 8, min: -20, max: 100, group: 'Financial constraints', role: 'bc', help: 'Rate earned on released cash; MIRR removes the IRR assumption that inflows are reinvested at the IRR itself.' },
    { key: 'compounding', label: 'Compounding of the quoted rate', type: 'select', default: '1', group: 'Financial constraints', role: 'control',
      options: [{ value: '1', label: 'Annual' }, { value: '2', label: 'Semi-annual' }, { value: '4', label: 'Quarterly' }, { value: '12', label: 'Monthly' }, { value: 'inf', label: 'Continuous' }],
      help: 'The quoted rate is converted to the effective annual rate used for discounting.' },
    { key: 'midYear', label: 'Mid-year discounting', type: 'bool', default: false, group: 'Financial constraints', role: 'control', help: 'Flows of year k ≥ 1 taken at k − 0.5 (continuous receipts).' },
    { key: 'maxPayback', label: 'Maximum discounted payback', unit: 'yr', type: 'number', default: 12, min: 0.5, max: 100, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'minPI', label: 'Minimum profitability index', unit: '-', type: 'number', default: 1, min: 0, max: 10, group: 'Acceptance criteria', role: 'criteria' },
  ],
  compute(inp) {
    const rows = tableRows(inp.flows, 'Cash-flow lines', { year: { min: 0, max: 100 }, investment: { min: 0 }, benefit: { min: 0 }, cost: { min: 0 }, salvage: { min: 0 } });
    for (const r of rows) need(Number.isInteger(r.year), `Cash-flow lines: year ${r.year} is not a whole year.`);
    const n = Math.max(...rows.map((r) => r.year));
    need(n >= 1, 'The cash-flow table must span at least years 0 and 1.');
    const inv = new Array(n + 1).fill(0); const ben = new Array(n + 1).fill(0); const om = new Array(n + 1).fill(0); const sal = new Array(n + 1).fill(0);
    for (const r of rows) { inv[r.year] += r.investment; ben[r.year] += r.benefit; om[r.year] += r.cost; sal[r.year] += r.salvage; }
    const net = inv.map((_, t) => ben[t] + sal[t] - om[t] - inv[t]);
    need(inv.some((v) => v > 0), 'The table contains no investment: at least one row needs an investment amount for the return measures to be defined.');
    need(net.some((v) => v !== 0), 'All net cash flows are zero.');

    // quoted rate → effective annual rate
    const m = inp.compounding === 'inf' ? Infinity : Number(inp.compounding);
    const i = run(() => effectiveRate(inp.rate / 100, m), 'Discount rate');
    need(i > -0.99, 'The effective discount rate must be above −99 %.');
    const f = inp.inflation / 100;
    const other = inp.basis === 'real' ? realToNominal(i, f) : nominalToReal(i, f);
    const o = { midYear: !!inp.midYear };

    const NPV = npv(net, i, o);
    const ir = irrReport(net);
    const MIRR = mirr(net, inp.financeRate / 100, inp.reinvestRate / 100);
    const PI = profitabilityIndex(net, i, { investment: inv, ...o });
    const BC = benefitCostRatio({ benefits: ben, investment: inv, om, salvage: sal, rate: i });
    const R = roi(net, { investment: inv });
    const pb = payback(net, i, o);
    const eaw = equivalentAnnualWorth(net, i);
    const ft = factorTable(i, n);
    const ds = discountSeries(net, i, o);
    const cont = fvContinuous(1, inp.rate / 100, 1) - 1;

    const kpis = [
      { key: 'npv', label: `Net present value (${inp.basis})`, value: mil(NPV), unit: MU, limit: 0, sense: 'min',
        explain: `Sum of all yearly net flows discounted to ${inp.baseYear} at ${pct(i)}. Positive means the project earns more than the hurdle rate; it is the primary decision measure.` },
    ];
    if (ir.value != null) kpis.push({ key: 'irr', label: 'Internal rate of return', value: round(ir.value * 100), unit: '%/yr', limit: round(i * 100), sense: 'min',
      explain: 'Discount rate at which NPV is zero. Comparable with the hurdle rate only when it is unique or the single positive root; see the IRR-count KPI.' });
    kpis.push(
      { key: 'irrCount', label: 'Number of real IRRs', value: ir.count, unit: '-',
        explain: 'Count of real rates above −100 % that zero the NPV. One is the normal case; zero or more than one means IRR must not be used alone and NPV/MIRR govern.' });
    if (Number.isFinite(MIRR)) kpis.push({ key: 'mirr', label: 'Modified internal rate of return', value: round(MIRR * 100), unit: '%/yr',
      explain: 'Single-valued return with outflows financed at the finance rate and inflows reinvested at the reinvestment rate; always unique, so it replaces IRR when IRR is ambiguous.' });
    kpis.push(
      { key: 'pi', label: 'Profitability index', value: round(PI.pi), unit: '-', limit: inp.minPI, sense: 'min',
        explain: 'Present value of the returns per unit of present value invested (1 + NPV / PV of investment). Above 1 the project adds value; used to rank projects when capital is rationed.' },
      { key: 'bcr', label: 'Benefit–cost ratio (conventional)', value: round(BC.conventional), unit: '-',
        explain: 'PV of benefits over PV of investment plus operating cost less salvage. Above 1 the benefits outweigh all costs; the modified ratio nets operating cost off the benefits instead.' },
      { key: 'bcrMod', label: 'Benefit–cost ratio (modified)', value: round(BC.modified), unit: '-',
        explain: 'PV of benefits net of operating cost over PV of investment less salvage. Gives the same accept/reject answer as the conventional ratio but a different magnitude.' },
      { key: 'roi', label: 'Return on investment (undiscounted)', value: round(R.roi * 100), unit: '%',
        explain: 'Total net return less investment, over investment, with no discounting. A quick accounting measure that overstates long-dated projects; never use it alone.' },
    );
    if (pb.simple != null) kpis.push({ key: 'paybackSimple', label: 'Simple payback', value: round(pb.simple), unit: 'yr',
      explain: 'Years from time zero until cumulative undiscounted net flow reaches zero. A liquidity measure; it ignores everything after recovery.' });
    if (pb.discounted != null) kpis.push({ key: 'paybackDisc', label: 'Discounted payback', value: round(pb.discounted), unit: 'yr', limit: inp.maxPayback, sense: 'max',
      explain: 'Years until cumulative discounted net flow reaches zero. Longer than simple payback; if it does not exist the NPV is negative.' });
    kpis.push(
      { key: 'eaw', label: 'Equivalent annual worth', value: mil(eaw), unit: `${MU}/yr`,
        explain: 'The NPV spread as a level annual amount over the study period (NPV × capital recovery factor). Use it to compare alternatives of different lives.' },
      { key: 'crf', label: 'Capital recovery factor (A/P)', value: round(ft.AP, 6), unit: '1/yr',
        explain: `Level annual payment that repays one unit of capital over ${n} years at the discount rate; multiplies capex to give its annual equivalent.` },
      { key: 'effRate', label: 'Effective annual discount rate', value: round(i * 100, 6), unit: '%/yr',
        explain: `Quoted rate converted for its compounding frequency; with continuous compounding the same quoted rate gives ${pct(cont, 5)}. All discounting here uses this effective rate.` },
      { key: 'otherRate', label: inp.basis === 'real' ? 'Equivalent nominal rate' : 'Equivalent real rate', value: round(other * 100), unit: '%/yr',
        explain: 'The same hurdle expressed in the other terms through the Fisher relation (1 + nominal) = (1 + real)(1 + inflation). Discounting a nominal table with a real rate, or the reverse, is a common and material error.' },
    );

    const checks = [
      check('Value test: PV of costs ≤ PV of benefits', BC.pvInvestment + BC.pvOm - BC.pvSalvage, BC.pvBenefits, CUR, 'Blank & Tarquin ch. 9 (B/C ≥ 1 ⇔ NPV ≥ 0)',
        'Guards against committing capital that does not recover its cost at the hurdle rate. A ratio above 1 means every unit of benefit costs more than one unit in present value.'),
      checkMin('Profitability index', PI.pi, inp.minPI, '-', 'Brealey, Myers & Allen ch. 5 (capital rationing)',
        'Guards against low value per unit of scarce capital. Below the minimum the project crowds out better uses of the budget.'),
    ];
    if (ir.value != null) checks.push(checkMin('IRR ≥ hurdle rate', ir.value * 100, Math.max(i * 100, 1e-9), '%/yr', 'Blank & Tarquin ch. 7',
      'Guards against a return below the cost of capital. Valid only for the unique or single positive root; with several roots rely on the NPV check.'));
    checks.push(pb.discounted != null
      ? check('Discounted payback ≤ limit', pb.discounted, inp.maxPayback, 'yr', 'Blank & Tarquin ch. 13', 'Guards against capital being exposed for longer than the investor accepts; a long payback raises exposure to late-life price and integrity risk.')
      : { name: 'Discounted payback ≤ limit', demand: round(n), capacity: round(inp.maxPayback), unit: 'yr', uc: 99, ref: 'Blank & Tarquin ch. 13', explain: 'The discounted cumulative flow never reaches zero within the table, so the investment is not recovered at the hurdle rate.' });

    const years = net.map((_, t) => inp.baseYear + t);
    const hiRate = Math.max(0.3, (ir.roots.length ? Math.max(...ir.roots) : 0) * 1.5, i * 2);
    const prof = npvProfile(net, linspace(Math.max(-0.5, Math.min(0, ...ir.roots) - 0.05), Math.min(hiRate, 5), 121));
    const plots = [
      { type: 'bar', title: 'Net cash flow by year', categories: years.map(String), x: { label: 'Year' }, y: { label: 'Cash flow', unit: MU },
        series: [{ name: 'Benefit + salvage', y: ben.map((b, t) => mil(b + sal[t])) }, { name: 'Operating cost', y: om.map((v) => -mil(v)) }, { name: 'Investment', y: inv.map((v) => -mil(v)) }],
        explain: 'Yearly inflows above the axis and outflows below. Look for late outflows (overhauls, decommissioning): they are what create a second IRR.' },
      { type: 'line', title: 'Cumulative cash position', x: { label: 'Year', unit: '' }, y: { label: 'Cumulative', unit: MU },
        series: [{ name: 'Undiscounted', x: years, y: ds.cumulative.map(mil) }, { name: `Discounted at ${pct(i)}`, x: years, y: ds.cumulativeDiscounted.map(mil) }, { name: 'Break-even', x: [years[0], years[n]], y: [0, 0], dash: true }],
        explain: 'Where each curve crosses zero is the simple and the discounted payback; the end value of the discounted curve is the NPV. The deepest point is the maximum capital exposure.' },
      { type: 'line', title: 'NPV profile', x: { label: 'Discount rate', unit: '%/yr' }, y: { label: 'NPV', unit: MU },
        series: [{ name: 'NPV', x: prof.rates.map((r) => round(r * 100)), y: prof.npv.map(mil) }, { name: 'Zero', x: [prof.rates[0] * 100, prof.rates[prof.rates.length - 1] * 100], y: [0, 0], dash: true },
          { name: 'Hurdle rate', x: [i * 100, i * 100], y: [Math.min(0, ...prof.npv.map(mil)), Math.max(0, ...prof.npv.map(mil))], dash: true }],
        explain: 'NPV against discount rate. Each crossing of zero is an IRR: one crossing is the conventional case, none or several mean IRR is not a usable criterion. Read the NPV at the hurdle-rate line.' },
    ];
    const tables = [
      { title: 'Discounted cash flow', columns: [{ key: 'year', label: 'Year' }, { key: 'inv', label: 'Investment', unit: MU }, { key: 'ben', label: 'Benefit', unit: MU }, { key: 'om', label: 'Operating cost', unit: MU },
        { key: 'sal', label: 'Salvage', unit: MU }, { key: 'net', label: 'Net', unit: MU }, { key: 'df', label: 'Discount factor' }, { key: 'disc', label: 'Discounted', unit: MU }, { key: 'cum', label: 'Cumulative discounted', unit: MU }],
        rows: net.map((v, t) => ({ year: years[t], inv: mil(inv[t]), ben: mil(ben[t]), om: mil(om[t]), sal: mil(sal[t]), net: mil(v), df: round(ds.factors[t], 5), disc: mil(ds.discounted[t]), cum: mil(ds.cumulativeDiscounted[t]) })),
        explain: 'One row per year with the timing of every flow, its discount factor and the running discounted total; use it to audit that each cost and revenue falls in its intended year.' },
      { title: `Compound-interest factors at ${pct(i)} and n = ${n}`, columns: [{ key: 'f', label: 'Factor' }, { key: 'v', label: 'Value' }, { key: 'use', label: 'Converts' }],
        rows: [['F/P', ft.FP, 'present sum → future sum'], ['P/F', ft.PF, 'future sum → present sum'], ['P/A', ft.PA, 'level annual amount → present sum'], ['A/P', ft.AP, 'present sum → level annual amount (capital recovery)'],
          ['F/A', ft.FA, 'level annual amount → future sum'], ['A/F', ft.AF, 'future sum → level annual amount (sinking fund)'], ['P/G', ft.PG, 'arithmetic gradient → present sum'], ['A/G', ft.AG, 'arithmetic gradient → level annual amount']]
          .map(([fn, v, use]) => ({ f: fn, v: round(v, 6), use })),
        explain: 'The standard engineering-economy factors for the effective rate and the study period, for hand checks of any figure above.' },
      { title: 'IRR diagnostics', columns: [{ key: 'item', label: 'Item' }, { key: 'value', label: 'Value' }],
        rows: [{ item: 'Sign changes in the net flow (Descartes bound)', value: ir.signChanges }, { item: 'Real roots found', value: ir.count }, ...ir.roots.map((r, k) => ({ item: `Root ${k + 1}`, value: pct(r, 5) })), { item: 'Verdict', value: ir.text }],
        explain: 'Every real internal rate of return of the series and the test that says whether one of them can be used. More sign changes than one is the warning sign.' },
    ];

    const pvInv = PI.pvInvestment;
    const advice = [
      adv('economic', NPV >= 0 ? 'ok' : 'fail', `NPV is ${fM(NPV)} at ${pct(i)} ${inp.basis} (${pct(other)} ${inp.basis === 'real' ? 'nominal' : 'real'}); profitability index ${fmt(PI.pi)}, equivalent annual worth ${fM(eaw)}/yr.`,
        NPV >= 0 ? `Proceed on NPV; the present value of investment could rise by ${fM(NPV)} (${pct(NPV / pvInv)}) before value is destroyed — hold that as the ceiling on scope growth.`
          : `Reject or re-scope: the present value of investment must fall by ${fM(-NPV)} (${pct(-NPV / pvInv)}) or annual net benefit rise by ${fM(-eaw)}/yr to break even.`,
        'Net-present-value rule, Blank & Tarquin ch. 5; ISO 15686-5:2017 §5.4 (NPV as the primary whole-life measure).'),
      adv('economic', ir.count === 1 ? 'info' : 'warn', ir.text,
        ir.count === 1 ? 'IRR may be quoted alongside NPV; rank competing options on NPV or incremental IRR, not on IRR alone.'
          : ir.value != null ? `Quote the positive root together with this diagnostic and the MIRR${Number.isFinite(MIRR) ? ` (${pct(MIRR)})` : ''}; the negative root has no economic meaning.`
            : `Do not quote a single IRR. Report NPV and MIRR${Number.isFinite(MIRR) ? ` (${pct(MIRR)})` : ''}, and show the NPV profile.`,
        'Descartes rule of signs and Norstrom criterion, Blank & Tarquin ch. 7 (multiple rates of return).'),
      adv('engineering', 'info', `The value rests on ${fM(pvInv)} of discounted investment against ${fM(BC.pvBenefits)} of discounted benefit and ${fM(BC.pvOm)} of discounted operating cost.`,
        'Trace each line back to its technical driver: investment to installed quantities and vessel days, benefit to yield and availability, operating cost to the maintenance strategy. A one-year slip of first benefit costs about one year of discounting on the whole benefit stream.',
        'AACE RP 34R-05 (basis of estimate traceability).',
        `Late-life rows such as removal and disposal are discounted heavily (factor ${fmt(ds.factors[n])} in year ${n}) but are physical obligations: the steel removed can displace primary steel at about ${CO2_STEEL} t CO₂e per tonne (screening factor), so keep the liability and the recycling credit visible rather than letting discounting hide them.`),
    ];
    if (pb.discounted == null) advice.push(adv('economic', 'warn', 'The discounted cumulative cash flow never turns positive.', 'Extend the benefit period, reduce early capital or lower the cost of capital; otherwise the capital is not recovered.', 'Discounted payback, Blank & Tarquin ch. 13.'));
    else if (!pb.staysPositive) advice.push(adv('economic', 'warn', 'The cumulative position turns negative again after payback because of late outflows.', 'Provide for the late outflow (reserve or fund) and quote payback together with the end-of-life position.', 'Whole-life costing, ISO 15686-5:2017 §4.'));

    const publish = { 'econ.appraisalNpv': NPV };
    if (ir.value != null) publish['econ.appraisalIrr'] = ir.value;
    return { kpis, checks, plots, tables, advice, publish,
      notes: [`${inp.basis === 'real' ? 'Real' : 'Nominal'} terms throughout: the table and the ${pct(i)} rate are both ${inp.basis}.`, ir.text,
        `Discounting is ${inp.midYear ? 'mid-year for years ≥ 1' : 'end-of-year'}; year 0 is ${inp.baseYear}.`] };
  },
  verify: [
    // NPV = −1000 + 300 (P/A, 10 %, 5) = −1000 + 300 × 3.790787 = 137.236; IRR solves 300 (P/A, i, 5) = 1000 → 15.238 %;
    // PI = 1137.236/1000; simple payback 1000/300 = 3.3333; discounted: 950.960 after 4 yr, (1000 − 950.960)/186.276 = 0.2633 → 4.2633
    { name: 'Uniform series: −1000 then 300 × 5 at 10 %', source: 'Hand calculation with (P/A, 10 %, 5) = 3.790787 (Blank & Tarquin compound-interest tables)',
      inputs: { rate: 10, flows: [flow(0, 1000e6, 0, 0), ...[1, 2, 3, 4, 5].map((y) => flow(y, 0, 300e6, 0))] },
      expect: [{ kpi: 'npv', value: 137.236, tol: 1e-4 }, { kpi: 'irr', value: 15.238, tol: 1e-4 }, { kpi: 'pi', value: 1.137236, tol: 1e-4 }, { kpi: 'paybackSimple', value: 3.3333, tol: 1e-4 }, { kpi: 'paybackDisc', value: 4.2633, tol: 1e-3 }, { kpi: 'irrCount', value: 1, abs: 0 }] },
    // −1000 + 2300/(1+i) − 1320/(1+i)² = 0 → (1+i) = 1.1 or 1.2: two IRRs, 10 % and 20 %
    { name: 'Two IRRs: −1000, +2300, −1320 has roots 10 % and 20 %', source: 'Quadratic 1000 x² − 2300 x + 1320 = 0, x = 1 + i → x = 1.1, 1.2',
      inputs: { rate: 15, flows: [flow(0, 1000e6, 0, 0), flow(1, 0, 2300e6, 0), flow(2, 1320e6, 0, 0)] },
      expect: [{ kpi: 'irrCount', value: 2, abs: 0 }, { kpi: 'npv', value: -1000 + 2300 / 1.15 - 1320 / 1.15 ** 2, tol: 1e-4 }] },
    // 1000 − 2000 x + 2000 x² (x = 1/(1+i)) has discriminant 4 − 8 < 0: no real root
    { name: 'No IRR: +1000, −2000, +2000', source: 'Quadratic with negative discriminant: no real rate zeroes the NPV',
      inputs: { rate: 10, flows: [flow(0, 0, 1000e6, 0), flow(1, 2000e6, 0, 0), flow(2, 0, 2000e6, 0)] },
      expect: [{ kpi: 'irrCount', value: 0, abs: 0 }] },
    // Monthly compounding of 12 %: (1.01)^12 − 1 = 12.6825 %
    { name: 'Effective rate of 12 % compounded monthly', source: '(1 + 0.12/12)^12 − 1 = 0.126825', inputs: { rate: 12, compounding: '12' }, expect: [{ kpi: 'effRate', value: 12.6825, tol: 1e-5 }] },
  ],
};
