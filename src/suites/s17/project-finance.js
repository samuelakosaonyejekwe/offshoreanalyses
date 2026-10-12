// Project financing and after-tax investment appraisal: cost of capital (CAPM, WACC), unlevered after-tax cash
// flow with tax depreciation, debt sizing and repayment, coverage ratios, equity and project returns, and the
// one-at-a-time sensitivity of NPV and IRR.
import { buildProjectCashflow, cashflowMetrics, npv as npvOf } from '../../engine/econ/cashflow.js';
import { wacc, capm, financeProject } from '../../engine/econ/finance.js';
import { taxShieldPV } from '../../engine/econ/depreciation.js';
import { realToNominal, nominalToReal } from '../../engine/econ/tvm.js';
import { projectInputs } from './_project.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, CO2_STEEL, need, run, round, fmt, mil, fM, pct, adv, check, checkMin, cost, irrReport, irrFast, range } from './_common.js';

const DEP = [{ value: 'sl', label: 'Straight line' }, { value: 'ddb', label: 'Double declining balance (switch to straight line)' }, { value: 'syd', label: 'Sum of years’ digits' }];

export default {
  id: 's17-project-finance', suite: 17, title: 'Project finance and after-tax investment appraisal',
  summary: 'Derives the cost of capital, builds the after-tax project cash flow with tax depreciation and inflation, sizes and schedules the debt, and reports project and equity returns, coverage ratios, funding requirement and the sensitivity of NPV and IRR.',
  covers: [6], fidelity: 'analytical',
  standards: ['Yescombe, Principles of Project Finance, 2nd ed. (2014), ch. 12–13 (debt sizing, cover ratios, covenants)', 'Gatti, Project Finance in Theory and Practice, 3rd ed. (2018)', 'Brealey, Myers & Allen, Principles of Corporate Finance, ch. 6 and 19 (after-tax cash flow, WACC)', 'Sharpe (1964) capital asset pricing model', 'IAS 16 / IAS 12 (depreciation and tax treatment, for the timing of allowances)'],
  models: ['Weighted average cost of capital', 'Capital asset pricing model', 'Debt–equity financing with interest during construction', 'Amortising (annuity), equal-principal and sculpted debt repayment', 'Debt-service coverage ratio', 'Loan life coverage ratio', 'Project life coverage ratio', 'Equity internal rate of return', 'Project internal rate of return', 'Tax depreciation and loss carry-forward', 'After-tax cash flow', 'Inflation-adjusted valuation, real versus nominal discounting', 'Break-even price, volume and capital cost', 'One-at-a-time sensitivity (tornado)'],
  limits: `Annual model in nominal money (money of the day), which is required for tax depreciation and debt service to be right; real equivalents are quoted through the Fisher relation. One senior debt tranche, tax shield on interest taken when it arises, no refinancing, hedging or intra-year cash sweep. Lender models need semi-annual periods and the actual term sheet. ${CUR_NOTE}`,
  publishes: ['econ.npv', 'econ.wacc', 'econ.projectIrr', 'econ.equityIrr', 'econ.minDscr', 'econ.debt', 'econ.equity', 'econ.npvAfterTax'],
  searchTerms: 'project finance debt sizing DSCR LLCR equity IRR offshore wind',
  inputs: [
    ...projectInputs({}, ['basis', 'discountRate']),
    { key: 'depMethod', label: 'Tax depreciation method', type: 'select', default: 'sl', options: DEP, group: 'Fiscal regime', role: 'bc' },
    { key: 'workingCapital', label: 'Working capital', unit: '% of revenue', type: 'number', default: 2, min: 0, max: 50, group: 'Fiscal regime', role: 'ic', help: 'Opening working-capital requirement, recovered in the last operating year.' },
    { key: 'riskFree', label: 'Risk-free rate (nominal)', unit: '%/yr', type: 'number', default: 3.5, min: -2, max: 30, group: 'Cost of capital', role: 'bc', help: `Long government bond yield in the project currency. ${PLACEHOLDER}` },
    { key: 'beta', label: 'Equity beta (relevered)', unit: '-', type: 'number', default: 1.1, min: 0, max: 5, group: 'Cost of capital', role: 'bc' },
    { key: 'marketPremium', label: 'Equity market risk premium', unit: '%', type: 'number', default: 5.5, min: 0, max: 20, group: 'Cost of capital', role: 'bc', help: PLACEHOLDER },
    { key: 'countryPremium', label: 'Country and project risk premium', unit: '%', type: 'number', default: 1.0, min: 0, max: 30, group: 'Cost of capital', role: 'bc', help: PLACEHOLDER },
    { key: 'rateSource', label: 'Discount rate for project NPV', type: 'select', default: 'wacc', group: 'Cost of capital', role: 'bc', options: [{ value: 'wacc', label: 'Computed post-tax WACC' }, { value: 'input', label: 'Entered hurdle rate' }] },
    { key: 'hurdle', label: 'Hurdle rate (nominal, post-tax)', unit: '%/yr', type: 'number', default: 7.5, min: 0, max: 50, group: 'Cost of capital', role: 'bc', showIf: { key: 'rateSource', in: ['input'] }, help: PLACEHOLDER },
    { key: 'equityAvailable', label: 'Equity available from sponsors', unit: CUR, type: 'number', default: 6.5e8, min: 0, max: 1e12, group: 'Opening financial state', role: 'ic', help: 'Committed equity at financial close; the equity requirement is checked against it.' },
    { key: 'sizing', label: 'Debt sizing', type: 'select', default: 'dscr', group: 'Debt', role: 'bc', options: [{ value: 'dscr', label: 'To a target DSCR, capped by maximum gearing' }, { value: 'gearing', label: 'Fixed gearing' }] },
    { key: 'gearing', label: 'Gearing (debt share of funding)', unit: '%', type: 'number', default: 70, min: 0, max: 95, group: 'Debt', role: 'bc', help: 'Fixed gearing, or the cap when sizing to DSCR (borrowing limit).' },
    { key: 'targetDscr', label: 'Target DSCR for sizing', unit: '-', type: 'number', default: 1.30, min: 1.01, max: 5, group: 'Debt', role: 'bc', showIf: { key: 'sizing', in: ['dscr'] } },
    { key: 'debtType', label: 'Repayment profile', type: 'select', default: 'sculpted', group: 'Debt', role: 'bc', options: [{ value: 'sculpted', label: 'Sculpted to cash flow (constant DSCR)' }, { value: 'annuity', label: 'Annuity (level payment)' }, { value: 'equalPrincipal', label: 'Equal principal' }] },
    cost({ key: 'debtRate', label: 'All-in interest rate on debt', unit: '%/yr', default: 5.5, min: 0, max: 30, group: 'Debt', uncertain: { dist: 'normal', std: 0.75 }, help: 'Base rate plus margin, nominal.' }),
    { key: 'tenor', label: 'Debt tenor from first operating year', unit: 'yr', type: 'int', default: 18, min: 1, max: 40, group: 'Debt', role: 'bc' },
    cost({ key: 'arrangementFee', label: 'Arrangement fee', unit: '% of facility', default: 2, min: 0, max: 10, group: 'Debt' }),
    cost({ key: 'commitmentFee', label: 'Commitment fee on undrawn facility', unit: '%/yr', default: 0.8, min: 0, max: 5, group: 'Debt' }),
    { key: 'dsraMonths', label: 'Debt-service reserve account', unit: 'months', type: 'number', default: 6, min: 0, max: 24, group: 'Debt', role: 'bc', help: 'Liquidity requirement: months of debt service held in reserve.' },
    { key: 'minDscr', label: 'Minimum DSCR covenant', unit: '-', type: 'number', default: 1.15, min: 1, max: 3, group: 'Acceptance criteria', role: 'bc', help: 'Default covenant in the loan agreement.' },
    { key: 'lockupDscr', label: 'Distribution lock-up DSCR', unit: '-', type: 'number', default: 1.20, min: 1, max: 3, group: 'Acceptance criteria', role: 'bc' },
    { key: 'minLlcr', label: 'Minimum loan life coverage ratio', unit: '-', type: 'number', default: 1.25, min: 1, max: 3, group: 'Acceptance criteria', role: 'bc' },
    { key: 'minEquityIrr', label: 'Required equity IRR (nominal)', unit: '%/yr', type: 'number', default: 10, min: 0, max: 50, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'sensRange', label: 'Sensitivity range', unit: '± %', type: 'number', default: 10, min: 1, max: 50, group: 'Analysis control', role: 'control' },
  ],
  compute(inp) {
    const n = Math.round(inp.opLife); const nc = Math.round(inp.constructionYears); const f = inp.inflation / 100; const T = inp.taxRate / 100;
    need(inp.tenor <= n, `The debt tenor of ${inp.tenor} years exceeds the ${n}-year operating life: debt must be repaid from operating cash flow.`);
    // cost of capital (nominal)
    const kE = capm({ riskFree: inp.riskFree / 100, beta: inp.beta, marketPremium: inp.marketPremium / 100, countryRiskPremium: inp.countryPremium / 100 });
    const kD = inp.debtRate / 100;
    const deg = inp.degradation / 100; const e0 = inp.aep * 1e3;
    const production = range(n, (k) => e0 * (1 - deg) ** k);
    const model = (rate, factors) => buildProjectCashflow({ timeline: { constructionYears: nc, operatingYears: n, baseYear: inp.baseYear }, mode: 'nominal', inflation: f, discountRate: rate, factors,
      capex: [{ name: 'Capital expenditure', amount: inp.capex, phasing: range(nc, () => 1 / nc), startYear: 0 }, ...(inp.sunk > 0 ? [{ name: 'Sunk cost', amount: inp.sunk, year: 0, sunk: true }] : [])],
      opex: { fixed: inp.opexFixed, variablePerUnit: inp.opexVar }, production, price: { type: inp.priceType, price: inp.price }, royaltyRate: inp.royalty / 100,
      tax: { rate: T, lossCarryForward: true, depreciation: { method: inp.depMethod, life: Math.round(inp.depLife) } }, workingCapital: { fractionOfRevenue: inp.workingCapital / 100 }, decommissioning: { cost: inp.decom } });
    // first pass at a provisional rate to obtain the cash flow, then finance it and derive the WACC at the achieved gearing
    const m0 = run(() => model(0.07), 'Project cash flow');
    const debtSpec = { rate: kD, tenor: Math.round(inp.tenor), type: inp.debtType, arrangementFeeRate: inp.arrangementFee / 100, commitmentFeeRate: inp.commitmentFee / 100, dsraMonths: inp.dsraMonths, lockupDscr: inp.lockupDscr };
    if (inp.sizing === 'dscr') { debtSpec.targetDscr = inp.targetDscr; debtSpec.maxGearing = inp.gearing / 100; } else debtSpec.gearing = inp.gearing / 100;
    need(m0.fcff.slice(m0.timeline.firstOperatingRow).some((v) => v > 0), 'The project never generates positive operating cash flow, so it cannot carry debt: check price, production and operating cost.');
    const fin = run(() => financeProject({ fcff: m0.fcff, firstOperatingRow: m0.timeline.firstOperatingRow, debt: debtSpec, taxRate: T, equityDiscountRate: kE }), 'Financing');
    const W = wacc({ equity: Math.max(fin.equity, 1e-9), debt: Math.max(fin.debt, 0), costOfEquity: kE, costOfDebt: kD, taxRate: T });
    const rNom = inp.rateSource === 'wacc' ? W.postTax : inp.hurdle / 100;
    const rReal = nominalToReal(rNom, f);
    const m = run(() => model(rNom), 'Project cash flow');
    const NPV = m.metrics.npv;
    const pIrr = irrReport(m.fcff); const eIrr = irrReport(fin.equityCashflows);
    const met = cashflowMetrics(m.fcff, rNom);
    const shield = taxShieldPV({ charges: m.table.depreciation.slice(m.timeline.firstOperatingRow), taxRate: T, rate: rNom, startPeriod: m.timeline.firstOperatingRow });
    const dsraInit = Math.max(...fin.dsraBalance, 0);
    const funding = fin.debt + fin.equity;
    // one-at-a-time sensitivity of NPV and IRR
    const s = inp.sensRange / 100;
    const sens = [['Price', 'price'], ['Production', 'production'], ['Capital cost', 'capex'], ['Operating cost', 'opex']].map(([name, key]) => {
      const lo = model(rNom, { [key]: 1 - s }); const hi = model(rNom, { [key]: 1 + s });
      const il = irrFast(lo.fcff); const ih = irrFast(hi.fcff);
      return { name, key, low: lo.metrics.npv, high: hi.metrics.npv, irrLow: Number.isFinite(il) ? il : null, irrHigh: Number.isFinite(ih) ? ih : null };
    });
    // break-even multipliers on price and production by secant iteration started from the sensitivity runs
    const be = {};
    for (const key of ['price', 'production']) {
      const sv = sens.find((x) => x.key === key);
      let f0 = 1 - s; let g0 = sv.low; let f1 = 1 + s; let g1 = sv.high; let ok = false;
      for (let it = 0; it < 8 && g1 !== g0; it++) {
        const f2 = Math.min(Math.max(f1 - g1 * (f1 - f0) / (g1 - g0), 0.02), 20);
        const g2 = model(rNom, { [key]: f2 }).metrics.npv;
        f0 = f1; g0 = g1; f1 = f2; g1 = g2;
        if (Math.abs(g2) < 1e-7 * Math.max(Math.abs(inp.capex), 1)) { ok = true; break; }
      }
      be[key === 'price' ? 'price' : 'volume'] = ok ? { value: key === 'price' ? f1 * inp.price : f1 } : null;
    }

    const dr = Math.max(rNom * s, 0.0025);
    sens.push({ name: 'Discount rate', low: npvOf(m.fcff, rNom - dr), high: npvOf(m.fcff, rNom + dr), irrLow: pIrr.value, irrHigh: pIrr.value });
    sens.sort((a, b) => Math.abs(b.high - b.low) - Math.abs(a.high - a.low));

    const first = m.timeline.firstOperatingRow;
    // an all-equity case (no debt raised) has no debt service: cover ratios and their covenants do not apply
    const hasDebt = fin.debt > 0 && Number.isFinite(fin.minDscr);
    const dscrVals = fin.dscr.map((v, t) => ({ t, v })).filter((x) => x.v != null && Number.isFinite(x.v));
    const kpis = [
      { key: 'wacc', label: 'WACC (post-tax, nominal)', value: round(100 * W.postTax), unit: '%/yr', explain: `Weighted cost of ${pct(W.weightEquity)} equity at ${pct(kE)} and ${pct(W.weightDebt)} debt at ${pct(kD)} after tax relief on interest. The real equivalent is ${pct(nominalToReal(W.postTax, f))}.` },
      { key: 'costEquity', label: 'Cost of equity (CAPM)', value: round(100 * kE), unit: '%/yr', explain: 'Risk-free rate plus beta times the market premium plus the country premium: the return shareholders require for this risk.' },
      { key: 'npv', label: 'Project NPV (after tax, unlevered)', value: mil(NPV), unit: MU, limit: 0, sense: 'min', explain: `Present value at ${inp.baseYear} of the after-tax free cash flow at ${pct(rNom)} nominal (${pct(rReal)} real). The same value results from real flows at the real rate.` },
    ];
    if (pIrr.value != null) kpis.push({ key: 'projectIrr', label: 'Project IRR (nominal, after tax)', value: round(100 * pIrr.value), unit: '%/yr', limit: round(100 * rNom), sense: 'min', explain: `Return on the whole investment regardless of financing; ${pct(nominalToReal(pIrr.value, f))} in real terms. Compare with the WACC.` });
    kpis.push({ key: 'projectIrrCount', label: 'Number of project IRRs', value: pIrr.count, unit: '-', explain: 'Real roots of the project cash flow. Decommissioning at the end adds a second, negative root; the single positive root is the economic return.' });
    if (eIrr.value != null) kpis.push({ key: 'equityIrr', label: 'Equity IRR (nominal)', value: round(100 * eIrr.value), unit: '%/yr', limit: inp.minEquityIrr, sense: 'min', explain: 'Return on the sponsors’ cash after debt service and tax. Leverage raises it above the project IRR when the project earns more than the cost of debt.' });
    kpis.push(
      { key: 'equityNpv', label: 'Equity NPV at the cost of equity', value: mil(fin.equityNpv), unit: MU, explain: 'Present value of the equity cash flows discounted at the CAPM cost of equity: value created for the sponsors.' },
      { key: 'debt', label: 'Senior debt', value: mil(fin.debt), unit: MU, explain: `Debt raised including capitalised interest and fees; sized by ${fin.binding === 'dscr' ? 'the target DSCR' : 'the gearing limit'}.` },
      { key: 'gearing', label: 'Gearing achieved', value: round(100 * fin.gearing), unit: '%', explain: 'Debt as a share of total funding. Lower than the cap when cash flow, not the cap, limits borrowing.' },
      { key: 'equity', label: 'Equity required', value: mil(fin.equity), unit: MU, limit: mil(inp.equityAvailable), sense: 'max', explain: 'Sponsor cash needed during construction after debt drawdown.' },
      { key: 'funding', label: 'Total funding requirement', value: mil(funding), unit: MU, explain: `Capital cost in money of the day plus interest during construction of ${fM(fin.idc)} and financing fees of ${fM(fin.fees)}.` },
      { key: 'mirr', label: 'Project MIRR', value: round(100 * met.mirr), unit: '%/yr', explain: 'Return with outflows financed and inflows reinvested at the discount rate; unique where IRR is not.' },
      { key: 'payback', label: 'Discounted payback', value: round(met.paybackDiscounted ?? m.n), unit: 'yr', explain: met.paybackDiscounted == null ? 'The discounted cumulative cash flow does not turn positive within the model; the horizon is shown.' : 'Years from time zero until the discounted cumulative after-tax cash flow reaches zero.' },
      { key: 'taxShield', label: 'PV of depreciation tax shield', value: mil(shield.pv), unit: MU, explain: 'Present value of the tax saved by depreciation allowances. Faster depreciation methods raise it without changing the total allowance.' },
      { key: 'taxPaid', label: 'Total corporate tax paid', value: mil(m.totals.tax), unit: MU, explain: 'Undiscounted tax over the life after depreciation and loss carry-forward.' },
    );
    if (hasDebt) kpis.push(
      { key: 'minDscr', label: 'Minimum DSCR', value: round(fin.minDscr), unit: '-', limit: inp.minDscr, sense: 'min', explain: 'Lowest ratio of cash available for debt service to debt service in any year. Below the covenant the loan is in default; below the lock-up level dividends stop.' },
      { key: 'avgDscr', label: 'Average DSCR', value: round(fin.averageDscr), unit: '-', explain: 'Mean coverage over the loan life; with sculpted repayment it equals the minimum by construction.' },
      { key: 'llcr', label: 'Loan life coverage ratio', value: round(fin.llcr ?? 0), unit: '-', limit: inp.minLlcr, sense: 'min', explain: 'Present value of cash flow over the loan life divided by the debt outstanding at the start of operations: the lender’s measure of overall repayment capacity.' },
      { key: 'plcr', label: 'Project life coverage ratio', value: round(fin.plcr ?? 0), unit: '-', explain: 'As LLCR but over the whole project life; the excess over LLCR is the cash-flow “tail” after the loan matures.' },
    );
    if (be.price) kpis.push({ key: 'breakEvenPrice', label: 'Break-even price (after tax)', value: round(be.price.value), unit: `${CUR}/MWh`, explain: 'Unit price (base-year money, same contract type) at which the after-tax NPV is zero at the discount rate.' });
    if (be.volume) kpis.push({ key: 'breakEvenVolume', label: 'Break-even production', value: round(100 * be.volume.value), unit: '% of forecast', explain: 'Fraction of forecast production at which the after-tax NPV is zero: the yield shortfall the project can absorb.' });

    const checks = [
      ...(hasDebt ? [checkMin('Minimum DSCR ≥ covenant', fin.minDscr, inp.minDscr, '-', 'Yescombe (2014) §12.3 (annual debt-service cover ratio)', 'Guards against cash flow that cannot service the debt in a weak year; below the covenant the lenders may accelerate the loan.'),
      checkMin('LLCR ≥ minimum', fin.llcr ?? 0, inp.minLlcr, '-', 'Yescombe (2014) §12.3 (loan-life cover ratio)', 'Guards against insufficient total cash flow over the loan life to repay the debt, even if single years pass.')] : []),
      check('Equity required ≤ equity available', fin.equity, Math.max(inp.equityAvailable, 1e-9), CUR, 'Financial constraint (committed equity)', 'Guards against a funding gap at financial close; above 1.0 the sponsors must commit more equity or find mezzanine finance.'),
      eIrr.value != null ? checkMin('Equity IRR ≥ required return', 100 * eIrr.value, inp.minEquityIrr, '%/yr', 'Investment criterion (sponsor hurdle rate)', 'Guards against an equity return below what sponsors need for the risk; below it the project will not attract equity.')
        : { name: 'Equity IRR ≥ required return', demand: inp.minEquityIrr, capacity: 0, unit: '%/yr', uc: 99, ref: 'Investment criterion (sponsor hurdle rate)', explain: 'No usable equity IRR exists for this cash flow, so the return cannot be shown to meet the hurdle; decide on equity NPV.' },
    ];
    const yrs = m.table.calendarYear;
    const plots = [
      { type: 'bar', title: 'After-tax project cash flow', categories: yrs.map(String), x: { label: 'Year' }, y: { label: 'Cash flow (nominal)', unit: MU },
        series: [{ name: 'Revenue', y: m.table.revenue.map(mil) }, { name: 'Operating cost and royalty', y: m.table.opex.map((v, t) => -mil(v + m.table.royalty[t])) }, { name: 'Tax', y: m.table.tax.map((v) => -mil(v)) }, { name: 'Capital and decommissioning', y: m.table.capex.map((v, t) => -mil(v + m.table.decomSpend[t])) }],
        explain: 'Yearly nominal inflows and outflows of the unlevered project. Tax starts once depreciation and carried-forward losses are used up.' },
      { type: 'line', title: 'Debt service and coverage', x: { label: 'Year', unit: '' }, y: { label: 'Amount', unit: MU },
        series: [{ name: 'Cash available for debt service', x: yrs.slice(first), y: fin.cfads.slice(first).map(mil) }, { name: 'Debt service', x: yrs.slice(first), y: fin.debtService.slice(first).map(mil) }, { name: 'Debt outstanding', x: yrs, y: fin.balance.map(mil), dash: true }],
        explain: 'The gap between cash available and debt service is the lenders’ cushion; where the two lines approach, the DSCR is at its minimum. The dashed line is the loan balance.' },
      { type: 'line', title: 'DSCR by year', x: { label: 'Year', unit: '' }, y: { label: 'DSCR', unit: '-' },
        series: [{ name: 'DSCR', x: dscrVals.map((d) => yrs[d.t]), y: dscrVals.map((d) => round(d.v)) }, { name: 'Covenant', x: [yrs[first], yrs[yrs.length - 1]], y: [inp.minDscr, inp.minDscr], dash: true }, { name: 'Lock-up', x: [yrs[first], yrs[yrs.length - 1]], y: [inp.lockupDscr, inp.lockupDscr], dash: true }],
        explain: 'Coverage in each year against the default covenant and the dividend lock-up level. Any point below a dashed line is a breach.' },
      { type: 'tornado', title: `NPV sensitivity (± ${inp.sensRange} %)`, base: mil(NPV), items: sens.map((x) => ({ name: x.name, low: mil(Math.min(x.low, x.high)), high: mil(Math.max(x.low, x.high)) })), x: { label: 'NPV', unit: MU },
        explain: 'NPV when each driver alone moves by the stated percentage (the discount rate by the same relative amount). The longest bar is the variable to pin down first.' },
      { type: 'bar', title: 'Equity cash flow', categories: yrs.map(String), x: { label: 'Year' }, y: { label: 'Cash flow to equity', unit: MU }, series: [{ name: 'Equity', y: fin.equityCashflows.map(mil) }],
        explain: 'Sponsor contributions during construction and distributions afterwards. Distributions rise when the debt is repaid; a year at zero is a lock-up.' },
    ];
    const tables = [
      { title: 'Cash-flow statement (nominal)', columns: [{ key: 'y', label: 'Year' }, { key: 'rev', label: 'Revenue', unit: MU }, { key: 'opex', label: 'Opex', unit: MU }, { key: 'dep', label: 'Depreciation', unit: MU }, { key: 'tax', label: 'Tax', unit: MU }, { key: 'capex', label: 'Capex', unit: MU }, { key: 'fcff', label: 'Free cash flow', unit: MU }, { key: 'ds', label: 'Debt service', unit: MU }, { key: 'dscr', label: 'DSCR' }, { key: 'eq', label: 'Equity cash flow', unit: MU }],
        rows: yrs.map((y, t) => ({ y, rev: mil(m.table.revenue[t]), opex: mil(m.table.opex[t]), dep: mil(m.table.depreciation[t]), tax: mil(m.table.tax[t]), capex: mil(m.table.capex[t] + m.table.decomSpend[t]), fcff: mil(m.fcff[t]), ds: mil(fin.debtService[t]), dscr: fin.dscr[t] != null && Number.isFinite(fin.dscr[t]) ? round(fin.dscr[t], 4) : '–', eq: mil(fin.equityCashflows[t]) })),
        explain: 'Year-by-year statement from revenue to equity cash flow; use it to check that capital, operating cost, tax, financing and decommissioning fall in their intended periods.' },
      { title: 'Sources and uses of funds', columns: [{ key: 'item', label: 'Item' }, { key: 'v', label: 'Amount', unit: MU }, { key: 's', label: 'Share', unit: '%' }],
        rows: [{ item: 'Uses: construction cost (money of the day)', v: mil(m.totals.capex), s: round(100 * m.totals.capex / funding, 4) }, { item: 'Uses: interest during construction', v: mil(fin.idc), s: round(100 * fin.idc / funding, 4) }, { item: 'Uses: financing fees', v: mil(fin.fees), s: round(100 * fin.fees / funding, 4) },
          { item: 'Sources: senior debt', v: mil(fin.debt), s: round(100 * fin.debt / funding, 4) }, { item: 'Sources: equity', v: mil(fin.equity), s: round(100 * fin.equity / funding, 4) }, { item: 'Reserve: peak DSRA balance', v: mil(dsraInit), s: round(100 * dsraInit / funding, 4) }],
        explain: 'What the funding pays for and where it comes from at financial close; uses and sources must balance.' },
      { title: `Sensitivity of NPV and IRR (± ${inp.sensRange} %)`, columns: [{ key: 'name', label: 'Driver' }, { key: 'lo', label: 'NPV low', unit: MU }, { key: 'hi', label: 'NPV high', unit: MU }, { key: 'il', label: 'IRR low', unit: '%' }, { key: 'ih', label: 'IRR high', unit: '%' }],
        rows: sens.map((x) => ({ name: x.name, lo: mil(x.low), hi: mil(x.high), il: x.irrLow != null ? round(100 * x.irrLow, 4) : '–', ih: x.irrHigh != null ? round(100 * x.irrHigh, 4) : '–' })),
        explain: 'NPV and project IRR with each driver at the low and high end of the range, one at a time and reproducible. Low and high refer to the driver value, not to the result.' },
    ];
    const top = sens[0];
    const lever = eIrr.value != null && pIrr.value != null ? eIrr.value - pIrr.value : null;
    const advice = [
      adv('economic', NPV >= 0 && (!hasDebt || fin.minDscr >= inp.minDscr) ? 'ok' : 'warn', `After-tax NPV ${fM(NPV)} at ${pct(rNom)} nominal (${pct(rReal)} real); project IRR ${pIrr.value != null ? pct(pIrr.value) : 'not defined'}, equity IRR ${eIrr.value != null ? pct(eIrr.value) : 'not defined'}. ${hasDebt ? `Debt ${fM(fin.debt)} at ${pct(fin.gearing)} gearing, minimum DSCR ${fmt(fin.minDscr)}, LLCR ${fmt(fin.llcr ?? 0)}.` : 'No debt is raised: the project is funded entirely by equity and no cover ratio applies.'}`,
        hasDebt && fin.minDscr < inp.minDscr ? `Coverage is below the ${inp.minDscr} covenant: reduce gearing, lengthen the tenor or sculpt repayments to cash flow.` : NPV < 0 ? `Value is negative: the break-even price is ${be.price ? fmt(be.price.value) : '–'} ${CUR}/MWh against ${inp.price}.` : `${fin.binding === 'dscr' ? `Cash flow limits debt to ${pct(fin.gearing)} gearing, below the ${inp.gearing} % cap; a higher contracted price share or longer tenor would release more debt.` : 'The gearing cap, not cash flow, limits debt; coverage has headroom.'}${lever != null ? ` Leverage adds ${fmt(100 * lever)} points to the equity return.` : ''}`,
        'Debt sizing and cover-ratio covenants, Yescombe, Principles of Project Finance (2014) ch. 12; WACC and after-tax cash-flow rules, Brealey, Myers & Allen ch. 6 and 19.'),
      adv('engineering', 'info', `"${top.name}" is the strongest driver: ± ${inp.sensRange} % moves NPV between ${fM(Math.min(top.low, top.high))} and ${fM(Math.max(top.low, top.high))}.${be.volume ? ` The project breaks even at ${pct(be.volume.value)} of forecast production.` : ''}`,
        `Tie that driver to its engineering source: production to the energy-yield and availability analyses, capital cost to quantities and vessel days, operating cost to the maintenance strategy. Lenders will run the P90 yield (exceedance convention: the yield exceeded with 90 % probability, equal to the cumulative P10) through this model — check that the minimum DSCR stays above 1.0 at ${be.volume ? pct(be.volume.value) : 'the downside'} production.`,
        'Lender sensitivity cases on energy yield and cost overrun, Gatti (2018) ch. 5; AACE RP 18R-97 for the capital estimate class behind the base case.',
        `Taxonomy-aligned or green-loan financing can lower the margin: 25 basis points on ${fM(fin.debt)} of debt is about ${fM(0.0025 * fin.debt)} a year. Eligibility needs evidence of life-cycle emissions (steel at about ${CO2_STEEL} t CO₂e/t dominates the embodied part) and a funded decommissioning plan, both of which this suite quantifies.`),
    ];
    if (pIrr.count !== 1 || eIrr.count !== 1) advice.push(adv('economic', 'info', `IRR diagnostics — project: ${pIrr.text} Equity: ${eIrr.text}`, 'Quote NPV and MIRR alongside any IRR; do not rank options on IRR when more than one root exists.', 'Blank & Tarquin ch. 7 (multiple rates of return), Norstrom criterion.'));
    const publish = { 'econ.npv': NPV, 'econ.wacc': W.postTax, 'econ.debt': fin.debt, 'econ.equity': fin.equity, 'econ.npvAfterTax': NPV };
    if (hasDebt) publish['econ.minDscr'] = fin.minDscr;
    if (pIrr.value != null) publish['econ.projectIrr'] = pIrr.value;
    if (eIrr.value != null) publish['econ.equityIrr'] = eIrr.value;
    return { kpis, checks, plots, tables, advice, publish,
      notes: [`Nominal terms throughout (money of the day, inflation ${inp.inflation} %/yr): the ${pct(rNom)} nominal discount rate equals ${pct(rReal)} real by the Fisher relation, and real flows at the real rate give the same NPV.`,
        `Sunk cost of ${fM(inp.sunk)} is excluded from the decision cash flow and kept only in the depreciation basis.`, `Debt: ${inp.debtType} repayment over ${inp.tenor} years at ${inp.debtRate} %; ${fin.lockedUp.filter(Boolean).length} year(s) in distribution lock-up.`, pIrr.text] };
  },
  verify: [
    // CAPM: 3 + 1.2 × 5 + 0 = 9 %
    { name: 'CAPM cost of equity', source: 'Sharpe (1964): kE = rf + β × MRP = 3 % + 1.2 × 5 % = 9 %', inputs: { riskFree: 3, beta: 1.2, marketPremium: 5, countryPremium: 0 }, expect: [{ kpi: 'costEquity', value: 9, tol: 1e-9 }] },
    // no tax, no inflation, capex 1000 M in year 0, EBITDA = 2000 GWh × 100 − 50 M = 150 M for 20 years, no decommissioning:
    // NPV at 8 % = −1000 + 150 × (P/A, 8 %, 20) = −1000 + 150 × 9.818147 = 472.722 M; IRR: (P/A, i, 20) = 6.6667 → i = 13.8866 % (bisection on the annuity factor)
    // debt 60 % annuity over 10 years at 5 %, no fees, draws at end of year 0 so no IDC: 600 M × (A/P, 5 %, 10) = 600 × 0.1295046 = 77.703 M/yr → DSCR = 150/77.703 = 1.9304
    // WACC at 60/40 with kE = 9 %: 0.4 × 9 + 0.6 × 5 = 6.6 %
    { name: 'Level cash flow, 60 % annuity debt: NPV, IRR, DSCR and WACC', source: 'Hand calculation with (P/A, 8 %, 20) = 9.818147 and (A/P, 5 %, 10) = 0.1295046',
      inputs: { inflation: 0, sunk: 0, aep: 2000, degradation: 0, constructionYears: 1, opLife: 20, capex: 1e9, opexFixed: 5e7, opexVar: 0, decom: 0, price: 100, taxRate: 0, royalty: 0, workingCapital: 0, rateSource: 'input', hurdle: 8,
        riskFree: 3, beta: 1.2, marketPremium: 5, countryPremium: 0, sizing: 'gearing', gearing: 60, debtType: 'annuity', debtRate: 5, tenor: 10, arrangementFee: 0, commitmentFee: 0, dsraMonths: 0, equityAvailable: 1e9 },
      expect: [{ kpi: 'npv', value: 472.722, tol: 1e-5 }, { kpi: 'projectIrr', value: 13.8866, tol: 1e-4 }, { kpi: 'minDscr', value: 1.9304, tol: 1e-4 }, { kpi: 'wacc', value: 6.6, tol: 1e-6 }, { kpi: 'debt', value: 600, tol: 1e-6 }] },
  ],
};
