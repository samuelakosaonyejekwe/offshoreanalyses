// Real options and strategic investment: value of waiting (deferral), expansion and abandonment options on one
// lattice, staged (compound) development, Black–Scholes and least-squares Monte Carlo cross-checks, and a decision
// tree for buying information before commitment.
// Unit trap: the engines work with continuously compounded rates. The risk-free rate is entered as an effective
// annual rate and converted with r = ln(1 + rf); the value leakage is entered directly as a continuous yield.
import { blackScholes, valueLattice, sequentialCompoundOption, projectVolatility, expandedNpv, decisionTreeRollback, longstaffSchwartz } from '../../engine/econ/realoptions.js';
import { valueOfFlexibility } from '../../engine/econ/decision.js';
import { evpi } from '../../engine/econ/risk.js';
import { assetOptions } from './_lattice.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, CO2_STEEL, need, run, round, fmt, mil, fM, pct, adv, check, cost, tableRows, uniqueNames, linspace, sampleLHS } from './_common.js';

const PHASES = [
  { name: 'Consent, surveys and concept', time: 0, cost: 3.0e7 },
  { name: 'Site investigation and FEED', time: 1, cost: 6.0e7 },
  { name: 'Construction commitment (FID)', time: 3, cost: 1.41e9 },
];

export default {
  id: 's17-real-options', suite: 17, title: 'Real options and staged investment',
  summary: 'Values the flexibility a static NPV ignores: the option to wait, to expand, to abandon and to develop in stages, with the investment trigger, an expanded NPV and the value of buying information before commitment.',
  covers: [10], fidelity: 'numerical',
  standards: ['Dixit & Pindyck (1994), Investment under Uncertainty, ch. 5 (value of waiting, investment trigger)', 'Trigeorgis (1996), Real Options, ch. 4–5 and 7 (expanded NPV, option interactions)', 'Copeland & Antikarov (2001), Real Options: A Practitioner’s Guide, ch. 5–6 and 9 (lattices, compound options, consolidated volatility)', 'Cox, Ross & Rubinstein (1979) binomial lattice; Black & Scholes (1973) and Merton (1973)', 'Longstaff & Schwartz (2001) least-squares Monte Carlo', 'AACE International RP 18R-97 (class of the estimate used as the exercise price)'],
  models: ['Real options valuation', 'Binomial lattice (Cox–Ross–Rubinstein)', 'Black–Scholes–Merton option value with Greeks', 'Investment deferral option and investment trigger', 'Expansion option', 'Abandonment option', 'Interacting options on one lattice', 'Staged development (sequential compound option)', 'Decision-tree analysis with expected-value and exponential-utility rollback', 'Expected value of perfect information', 'Project volatility by the consolidated (market-asset-disclaimer) approach', 'Least-squares Monte Carlo (Longstaff–Schwartz)', 'Expanded (strategic) NPV'],
  limits: `Project value follows a lognormal diffusion with constant volatility, risk-free rate and leakage; jumps, mean reversion of prices and competitive pre-emption are not modelled. Black–Scholes figures are a European reference only and are valid where the lognormal, constant-volatility assumption is defensible. The deferral option is valued on the bare project; expansion and abandonment are valued for the owner of the operating project, who receives its cash flow while holding it, and exercising one ends the other. The appraisal in the decision tree is assumed to reveal the state perfectly. ${CUR_NOTE}`,
  publishes: ['econ.optionValueDefer', 'econ.expandedNpv', 'econ.investmentTrigger', 'econ.stagedValue', 'econ.valueOfInformation'],
  searchTerms: 'real options offshore investment timing deferral option staged development binomial lattice',
  mesh: { key: 'steps', label: 'Lattice time steps', levels: [50, 100, 200, 400], qoi: 'deferValue' },
  inputs: [
    { key: 'baseYear', label: 'Valuation date (year 0)', unit: 'yr', type: 'int', default: 2026, min: 1990, max: 2100, group: 'Opening financial state', role: 'ic', help: 'All present values are at this date, in real money of this year.' },
    { key: 'V', label: 'Present value of operating cash flows if built today', unit: CUR, type: 'number', default: 1.65e9, min: 1, max: 1e13, group: 'Opening financial state', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.15 },
      help: 'Starting project valuation: PV of revenue less operating cost, tax and end-of-life cost, excluding the investment below.' },
    { key: 'sunk', label: 'Development cost already incurred (sunk)', unit: CUR, type: 'number', default: 4.0e7, min: 0, max: 1e12, group: 'Opening financial state', role: 'ic', help: 'Reported for the record and excluded from every value: only future spend enters the decision.' },
    cost({ key: 'invest', label: 'Investment to build (exercise price)', unit: CUR, default: 1.5e9, min: 1, max: 1e13, group: 'Opening financial state', bind: 'econ.capex', design: true, uncertain: { dist: 'lognormal', cov: 0.12 }, help: 'From the capital cost estimate when available; use a Class 3 or better estimate for a sanction decision.' }),

    { key: 'volSource', label: 'Volatility of project value', type: 'select', default: 'input', group: 'Uncertainty', role: 'ic', options: [{ value: 'input', label: 'Entered below' }, { value: 'consolidated', label: 'Derived from price volatility and operating leverage' }],
      help: 'The consolidated approach simulates the first-year return on project value from the price volatility and the cost share.' },
    { key: 'sigma', label: 'Volatility of project value', unit: '%/√yr', type: 'number', default: 22, min: 1, max: 150, group: 'Uncertainty', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.25 }, showIf: { key: 'volSource', in: ['input'] }, help: 'Annual standard deviation of the logarithmic return on project value.' },
    { key: 'priceVol', label: 'Volatility of the realised price', unit: '%/√yr', type: 'number', default: 15, min: 0.5, max: 150, group: 'Uncertainty', role: 'ic', help: 'Annual volatility of the uncontracted price; a fixed-price contract lowers it.' },
    { key: 'costShare', label: 'Operating cost as a share of revenue', unit: '%', type: 'number', default: 35, min: 0, max: 90, group: 'Uncertainty', role: 'operation', help: 'Operating leverage: fixed costs amplify price swings into larger swings of value.' },

    { key: 'riskFree', label: 'Risk-free rate (effective annual)', unit: '%/yr', type: 'number', default: 3.5, min: -2, max: 30, group: 'Financial constraints', role: 'bc', help: `Government bond yield over the option life. ${PLACEHOLDER}` },
    { key: 'hurdle', label: 'Risk-adjusted discount rate (decision tree)', unit: '%/yr', type: 'number', default: 8, min: 0, max: 50, group: 'Financial constraints', role: 'bc', help: `Hurdle rate applied to values deferred by the appraisal. ${PLACEHOLDER}` },
    { key: 'leakage', label: 'Value lost per year of waiting (cash yield)', unit: '%/yr', type: 'number', default: 5, min: 0, max: 40, group: 'Financial constraints', role: 'bc', design: true,
      help: 'Cash flow forgone while the project is not built, as a continuous yield on project value (about annual free cash flow / value). Zero means waiting costs nothing and early investment is never optimal.' },
    { key: 'window', label: 'Years the opportunity stays open', unit: 'yr', type: 'number', default: 5, min: 0.25, max: 30, group: 'Financial constraints', role: 'bc', help: 'Lease, licence or grid-connection expiry: the latest date the investment can be made.' },
    { key: 'exercise', label: 'Exercise of the options', type: 'select', default: 'american', group: 'Financial constraints', role: 'control', options: [{ value: 'american', label: 'Any time within the window (American)' }, { value: 'european', label: 'Only at the end of the window (European)' }] },

    { key: 'expFactor', label: 'Expansion: scale of the enlarged project', unit: '-', type: 'number', default: 1.3, min: 1.01, max: 5, group: 'Options on the operating asset', role: 'operation', design: true, help: '1.3 = 30 % more value-generating capacity after expansion.' },
    cost({ key: 'expCost', label: 'Expansion: cost', unit: CUR, default: 5.6e8, min: 0, max: 1e13, group: 'Options on the operating asset' }),
    cost({ key: 'salvage', label: 'Abandonment: net salvage value', unit: CUR, default: 4.5e8, min: 0, max: 1e13, group: 'Options on the operating asset', help: 'Resale, redeployment and scrap value less the cost of bringing removal forward.' }),
    { key: 'opWindow', label: 'Years in which the asset options can be used', unit: 'yr', type: 'number', default: 10, min: 0.5, max: 40, group: 'Options on the operating asset', role: 'bc', help: 'Terminal condition: after this date the expansion consent and the resale market are assumed gone.' },
    { key: 'steel', label: 'Recoverable structural steel', unit: 'tonnes', type: 'number', default: 62000, min: 0, max: 5e6, group: 'Options on the operating asset', role: 'geometry', bind: 'quantities.steel_t', help: 'For the scrap and carbon content of the salvage value.' },
    cost({ key: 'scrapPrice', label: 'Scrap steel price', unit: `${CUR}/t`, default: 260, min: 0, max: 5000, group: 'Options on the operating asset' }),

    { key: 'phases', label: 'Development stages', type: 'table', group: 'Staged development', role: 'cost', help: `Each stage is a gate: its cost is committed only if the project is still worth continuing. Time in years from the valuation date; the last stage buys the operating project. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Stage', type: 'text' }, { key: 'time', label: 'Decision time', unit: 'yr', type: 'number' }, { key: 'cost', label: 'Commitment', unit: CUR, type: 'number' }], default: PHASES },

    { key: 'pGood', label: 'Probability of the favourable outcome', unit: '%', type: 'number', default: 60, min: 1, max: 99, group: 'Information before commitment (decision tree)', role: 'ic', help: 'Opening uncertainty state, for example ground conditions or yield as assumed.' },
    { key: 'upFactor', label: 'Project value if favourable', unit: '× V', type: 'number', default: 1.25, min: 0.1, max: 10, group: 'Information before commitment (decision tree)', role: 'ic' },
    { key: 'downFactor', label: 'Project value if unfavourable', unit: '× V', type: 'number', default: 0.62, min: 0, max: 10, group: 'Information before commitment (decision tree)', role: 'ic' },
    cost({ key: 'appraisalCost', label: 'Cost of the appraisal', unit: CUR, default: 2.5e7, min: 0, max: 1e12, group: 'Information before commitment (decision tree)', help: 'Site investigation, measurement campaign or pilot that resolves the uncertainty.' }),
    { key: 'appraisalYears', label: 'Delay caused by the appraisal', unit: 'yr', type: 'number', default: 1, min: 0, max: 10, group: 'Information before commitment (decision tree)', role: 'operation' },
    { key: 'riskTolerance', label: 'Corporate risk tolerance', unit: CUR, type: 'number', default: 1.5e9, min: 1e3, max: 1e14, group: 'Information before commitment (decision tree)', role: 'bc', help: `Exponential-utility risk tolerance used for the certainty equivalent of the tree. ${PLACEHOLDER}` },

    { key: 'steps', label: 'Lattice time steps', unit: '-', type: 'int', default: 200, min: 20, max: 1000, group: 'Analysis control', role: 'control' },
    { key: 'lsmPaths', label: 'Least-squares Monte Carlo paths', unit: '-', type: 'int', default: 2000, min: 0, max: 50000, group: 'Analysis control', role: 'control', help: 'Independent simulation cross-check of the deferral option. 0 skips it.' },
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 17, min: 1, max: 1e6, group: 'Analysis control', role: 'control', help: 'Fixed seed: identical inputs reproduce identical results.' },
  ],
  compute(inp) {
    const V = inp.V; const I = inp.invest;
    const r = Math.log(1 + inp.riskFree / 100); // effective annual → continuous
    const q = inp.leakage / 100;
    const american = inp.exercise === 'american';
    need(inp.downFactor <= inp.upFactor, 'Decision tree: the unfavourable project value cannot exceed the favourable one.');

    // volatility: entered, or simulated first-year return on project value (consolidated approach).
    // Value is taken proportional to the operating margin (f − c)/(1 − c) for a price factor f with E[f] = 1.
    const c = inp.costShare / 100; const sp = inp.priceVol / 100;
    const f = sampleLHS({ f: { dist: 'lognormal', mu: -0.5 * sp * sp, sigma: sp } }, 2000, inp.seed).f;
    let sigCons = null; let consNote = '';
    try { sigCons = projectVolatility({ pv0: 1, pv1: Float64Array.from(f, (x) => (x - c) / (1 - c)) }).sigma; } catch (e) {
      consNote = 'The simulated operating margin turns negative at this price volatility and cost share, so the consolidated volatility is not defined.';
    }
    need(inp.volSource === 'input' || sigCons != null, `${consNote} Enter the volatility of project value directly.`);
    const sigma = inp.volSource === 'input' ? inp.sigma / 100 : sigCons;

    // a recombining lattice needs a risk-neutral probability inside (0, 1): σ√Δt > |r − q|Δt
    const minSteps = (T) => Math.ceil(1.2 * (r - q) ** 2 * T / (sigma * sigma)) + 1;
    const stepsFor = (T, n) => {
      const s = Math.max(Math.round(n), minSteps(T));
      need(s <= 20000, `Volatility of ${pct(sigma)} is too low for a lattice with a ${fmt(100 * (r - q))} %/yr drift over ${fmt(T)} years: the outcome is effectively certain, so use the static NPV.`);
      return s;
    };
    const N = stepsFor(inp.window, inp.steps); const Nop = stepsFor(inp.opWindow, inp.steps);
    const callLattice = (v, sg, n, am = american) => valueLattice({ V0: v, sigma: sg, r, T: inp.window, q, steps: n,
      terminal: (x) => Math.max(x - I, 0), exercise: am ? (x, cont) => Math.max(cont, x - I) : undefined });

    /* ---- (A) option to defer ---- */
    const staticNpv = V - I;
    const lat = run(() => callLattice(V, sigma, N), 'Deferral lattice');
    const defer = lat.value;
    const waiting = defer - Math.max(staticNpv, 0);
    const bs = run(() => blackScholes({ S: V, K: I, r, sigma, T: inp.window, q }), 'Black–Scholes');
    // immediate investment is optimal when the right to wait adds nothing to a positive static NPV (McDonald & Siegel 1986)
    const investNow = american && staticNpv > 0 && defer - staticNpv < 1e-9 * Math.max(1, V);
    // investment trigger V*: lowest project value at which immediate investment is optimal (bisection on the lattice)
    let trigger = null;
    if (american && q > 0) {
      const nT = Math.min(N, 100); const tol = 1e-7 * I;
      const g = (v) => callLattice(v, sigma, nT, true).value - (v - I);
      let lo = I; let hi = 2 * I; let k = 0;
      while (g(hi) > tol && k < 8) { lo = hi; hi *= 2; k++; }
      if (g(hi) <= tol) { for (let it = 0; it < 26; it++) { const m = 0.5 * (lo + hi); if (g(m) > tol) lo = m; else hi = m; } trigger = hi; }
    }
    let lsm = null;
    if (inp.lsmPaths >= 100 && american) lsm = run(() => longstaffSchwartz({ S: V, K: I, r, sigma, T: inp.window, q, steps: 50, paths: Math.round(inp.lsmPaths), seed: inp.seed, type: 'call' }), 'Least-squares Monte Carlo');

    /* ---- (B) options on the operating asset ---- */
    // Valued with the owner's lattice of this suite (cash flow received while the asset is held); the engine's
    // real-option lattice carries the terminal value only and overstates early exercise when the leakage is not zero.
    const oExp = { type: 'expand', factor: inp.expFactor, cost: inp.expCost }; const oAbn = { type: 'abandon', salvage: inp.salvage };
    const owner = (sg, n, options) => assetOptions({ V0: V, sigma: sg, r, q, T: inp.opWindow, steps: n, american, options });
    const exp = owner(sigma, Nop, [oExp]);
    const abn = owner(sigma, Nop, [oAbn]);
    const both = owner(sigma, Nop, [oExp, oAbn]);
    const interaction = both.optionValue - exp.optionValue - abn.optionValue;
    const en = expandedNpv({ staticNpv, optionValue: Math.max(both.optionValue, 0) });
    const buildNowValue = en.expandedNpv;

    /* ---- (C) staged development ---- */
    const ph = tableRows(inp.phases, 'Development stages', { time: { min: 0, max: 40 }, cost: { min: 0 } }).sort((a, b) => a.time - b.time);
    uniqueNames(ph, 'Development stages');
    const Tst = ph[ph.length - 1].time;
    const nSt = Tst > 0 ? stepsFor(Tst, Math.max(inp.steps, 120)) : 1;
    const staged = run(() => sequentialCompoundOption({ V, phases: ph.map((p) => ({ time: p.time, cost: p.cost })), r, sigma, q, steps: nSt }), 'Staged development');
    // continuation value before the first gate, for the "worth starting" test
    const c0 = ph[0].time === 0 ? ph[0].cost : 0;
    const cont0 = c0 > 0 ? run(() => sequentialCompoundOption({ V, phases: ph.map((p, k) => ({ time: p.time, cost: k === 0 ? 0 : p.cost })), r, sigma, q, steps: nSt }), 'Staged development').value : staged.value;
    const pvPhases = ph.map((p) => p.cost * Math.exp(-r * p.time));

    /* ---- (D) decision tree: commit, appraise first, or walk away ---- */
    const p = inp.pGood / 100; const df = (1 + inp.hurdle / 100) ** -inp.appraisalYears;
    const vG = V * inp.upFactor; const vB = V * inp.downFactor;
    const T0 = (v) => ({ type: 'terminal', value: v });
    const commit = { name: 'Commit now', payoff: -I, node: { type: 'chance', name: 'Outcome', branches: [{ name: 'favourable', p, node: T0(vG) }, { name: 'unfavourable', p: 1 - p, node: T0(vB) }] } };
    const walk = { name: 'Walk away', node: T0(0) };
    const after = (v, nm) => ({ type: 'decision', name: `After appraisal (${nm})`, options: [{ name: 'Build', node: T0((v - I) * df) }, { name: 'Do not build', node: T0(0) }] });
    const appraise = { name: 'Appraise first', payoff: -inp.appraisalCost, node: { type: 'chance', name: 'Appraisal result', branches: [{ name: 'favourable', p, node: after(vG, 'favourable') }, { name: 'unfavourable', p: 1 - p, node: after(vB, 'unfavourable') }] } };
    const flexible = { type: 'decision', name: 'Today', options: [commit, appraise, walk] };
    const rigid = { type: 'decision', name: 'Today', options: [commit, walk] };
    const tree = run(() => decisionTreeRollback(flexible), 'Decision tree');
    let treeCE = run(() => decisionTreeRollback(flexible, { riskTolerance: inp.riskTolerance }), 'Decision tree');
    // outcomes many times larger than the risk tolerance overflow the exponential utility: the decision maker is then
    // effectively risk-neutral on this scale or the tolerance is mis-stated, and the expected value is reported
    const ceOk = Number.isFinite(treeCE.certaintyEquivalent);
    if (!ceOk) treeCE = { ...tree, certaintyEquivalent: tree.expectedValue };
    const vof = run(() => valueOfFlexibility({ flexible, rigid }), 'Decision tree');
    const info = evpi([{ p, payoffs: [vG - I, 0] }, { p: 1 - p, payoffs: [vB - I, 0] }]);
    const evCommit = p * vG + (1 - p) * vB - I;
    const evAppraise = -inp.appraisalCost + (p * Math.max(vG - I, 0) + (1 - p) * Math.max(vB - I, 0)) * df;
    // largest appraisal cost at which appraising first still beats the better of committing and walking away
    const appraisalBreakEven = Math.max(evAppraise + inp.appraisalCost - Math.max(evCommit, 0), 0);

    const waitBetter = defer > buildNowValue + 1e-9 * I && !investNow;
    const kpis = [
      { key: 'staticNpv', label: 'Static NPV (build now, no flexibility)', value: mil(staticNpv), unit: MU, explain: `Project value less investment at ${inp.baseYear}, ignoring every later choice. The conventional rule accepts it when positive; real options ask what the choices are worth on top.` },
      { key: 'deferValue', label: 'Value of the opportunity with the right to wait', value: mil(defer), unit: MU, explain: `Value today of holding the right, not the obligation, to invest within ${fmt(inp.window)} years (${american ? 'American' : 'European'} call on project value). Above the static NPV, waiting has value; equal to it, invest now.` },
      { key: 'waitingValue', label: 'Value of waiting', value: mil(waiting), unit: MU, explain: 'Excess of the opportunity value over the larger of static NPV and zero. It is what is destroyed by committing today; zero means nothing is gained by delay.' },
      { key: 'bsValue', label: 'Black–Scholes reference (European)', value: mil(bs.call), unit: MU, explain: 'Closed-form value if investment were possible only at the end of the window. The lattice value exceeds it by the early-exercise premium when waiting costs cash flow.' },
      { key: 'delta', label: 'Sensitivity of the opportunity to project value (delta)', value: round(lat.delta), unit: '-', explain: 'Change of the opportunity value per unit change of project value. Near 1 the option is deep in the money and behaves like the project itself; near 0 it is an out-of-the-money bet on volatility.' },
      { key: 'sigmaUsed', label: 'Volatility of project value used', value: round(100 * sigma), unit: '%/√yr', explain: `Annual volatility in the lattice.${sigCons != null ? ` The consolidated estimate from ${inp.priceVol} % price volatility and ${inp.costShare} % cost share is ${pct(sigCons)} (operating leverage alone gives ${pct(sp / (1 - c))}).` : ` ${consNote}`}` },
    ];
    if (trigger != null) kpis.push({ key: 'trigger', label: 'Investment trigger V*', value: mil(trigger), unit: MU, explain: `Project value at which immediate investment becomes optimal: ${fmt(trigger / I)} times the investment instead of the 1.0 of the NPV rule. Below it the option to wait is worth more than the cash flow forgone.` });
    if (lsm) kpis.push({ key: 'lsmValue', label: 'Least-squares Monte Carlo cross-check', value: mil(lsm.value), unit: MU, explain: `Independent simulation estimate of the same deferral option (standard error ${fM(lsm.standardError)}). Agreement with the lattice within a few standard errors confirms the numerics, not the assumptions.` });
    kpis.push(
      { key: 'expandValue', label: 'Value of the option to expand', value: mil(exp.optionValue), unit: MU, explain: `Worth today of being able to scale the asset to ${fmt(inp.expFactor)} times for ${fM(inp.expCost)} within ${fmt(inp.opWindow)} years. It is the most that enabling pre-investment (spare capacity, consent) is worth.` },
      { key: 'abandonValue', label: 'Value of the option to abandon', value: mil(abn.optionValue), unit: MU, explain: `Worth today of being able to exit for ${fM(inp.salvage)} of net salvage: an insurance against the low-value outcomes. Larger for redeployable, removable assets.` },
      { key: 'combinedValue', label: 'Combined value of expansion and abandonment', value: mil(both.optionValue), unit: MU, explain: `Both options valued together on one lattice. It differs from the sum by ${fM(interaction)} because the two options interact; adding separately valued options double counts.` },
      { key: 'expandedNpv', label: 'Expanded NPV (build now with asset options)', value: mil(buildNowValue), unit: MU, limit: 0, sense: 'min', explain: 'Static NPV plus the combined option value: the worth of committing today with the later choices kept open. Compare with the value of the opportunity with the right to wait.' },
      { key: 'stagedValue', label: 'Value of the staged programme', value: mil(staged.value), unit: MU, explain: `Value today of developing in ${ph.length} gated stages, each committed only if the project still justifies it. Zero means the first stage is not worth starting.` },
      { key: 'stagedStatic', label: 'Static NPV of the staged programme', value: mil(staged.staticNpv), unit: MU, explain: 'Project value at completion less all stage costs discounted at the risk-free rate, as if every stage were committed in advance.' },
      { key: 'stagedFlex', label: 'Value of the stage gates', value: mil(staged.flexibilityValue), unit: MU, explain: 'What the right to stop at each gate adds over the committed programme. It grows with volatility and with the share of spend that sits behind later gates.' },
      { key: 'treeValue', label: 'Expected value of the best first decision', value: mil(tree.expectedValue), unit: MU, explain: `Rollback of the tree "commit now / appraise first / walk away": the optimal first move is "${tree.policy[0].choice}". Risk-neutral expected value.` },
      { key: 'treeCE', label: 'Certainty equivalent of the tree', value: mil(treeCE.certaintyEquivalent), unit: MU, explain: `Value of the tree to a decision maker with a risk tolerance of ${fM(inp.riskTolerance)}; the optimal first move is then "${treeCE.policy[0].choice}". ${ceOk ? 'The gap to the expected value is the risk premium.' : 'The outcomes are too large relative to that tolerance for the utility to be evaluated, so the expected value is shown: check the tolerance.'}` },
      { key: 'infoValue', label: 'Expected value of perfect information', value: mil(info.evpi), unit: MU, explain: 'Most that any appraisal can be worth before delay and cost: the expected gain from being able to decide after the uncertainty is resolved. An appraisal costing more than this cannot pay.' },
      { key: 'appraisalGain', label: 'Net gain from appraising first', value: mil(vof.valueOfFlexibility), unit: MU, explain: 'Value of the tree with the appraisal available less the value without it, after the appraisal cost and the delay. Zero means the appraisal does not change the decision enough to pay for itself.' },
    );

    const checks = [
      check('Investment ≤ project value with flexibility (expanded NPV ≥ 0)', I, V + Math.max(both.optionValue, 0), CUR, 'Trigeorgis (1996) ch. 4 (expanded NPV rule)',
        'Guards against committing capital that neither the cash flows nor the later choices recover. Above 1.0 the project destroys value even with its options counted.'),
    ];
    if (trigger != null) checks.push(check('Investment trigger reached (V ≥ V*)', trigger, V, CUR, 'Dixit & Pindyck (1994) ch. 5 (optimal investment rule)',
      'Guards against killing a valuable option by investing too early. Above 1.0 the project value is still below the trigger, so the right to wait is worth more than the cash flow lost by waiting.'));
    if (c0 > 0) checks.push(check('First stage worth starting (stage cost ≤ continuation value)', c0, cont0, CUR, 'Copeland & Antikarov (2001) ch. 6 (sequential compound options)',
      'Guards against starting a development whose first commitment exceeds the value of the right to continue. Above 1.0 the first stage should not be funded.'));

    // sweeps for the plots use a lighter lattice: 60 steps resolve the curve shape
    const nP = Math.min(N, 60); const nPo = Math.min(Nop, 60);
    const vs = linspace(0.4 * I, Math.max(2.0 * I, trigger ? 1.15 * trigger : 0, 1.2 * V), 25);
    const sgs = linspace(Math.max(0.05, 0.4 * sigma), Math.max(0.6, 1.8 * sigma), 9).filter((s) => 1.2 * (r - q) ** 2 * Math.max(inp.window, inp.opWindow) / (s * s) < 400);
    const atSig = (s) => {
      const n1 = Math.max(nP, Math.ceil(1.2 * (r - q) ** 2 * inp.window / (s * s)) + 1); const n2 = Math.max(nPo, Math.ceil(1.2 * (r - q) ** 2 * inp.opWindow / (s * s)) + 1);
      return { d: callLattice(V, s, n1).value - Math.max(staticNpv, 0),
        e: owner(s, n2, [oExp]).optionValue, a: owner(s, n2, [oAbn]).optionValue };
    };
    const sw = sgs.map(atSig);
    const plots = [
      { type: 'line', title: 'Opportunity value against project value', x: { label: 'Project value V', unit: MU }, y: { label: 'Value', unit: MU },
        series: [{ name: 'With the right to wait', x: vs.map(mil), y: vs.map((v) => mil(callLattice(v, sigma, nP).value)) }, { name: 'Invest now: max(V − I, 0)', x: vs.map(mil), y: vs.map((v) => mil(Math.max(v - I, 0))), dash: true },
          { name: 'This project', x: [mil(V), mil(V)], y: [0, mil(Math.max(defer, staticNpv, 0))], dash: true }, ...(trigger != null ? [{ name: 'Trigger V*', x: [mil(trigger), mil(trigger)], y: [0, mil(trigger - I)], dash: true }] : [])],
        explain: 'The gap between the curve and the invest-now line is the value of waiting. Where the two meet is the investment trigger; to its left, wait, to its right, invest.' },
      { type: 'line', title: 'Option values against volatility', x: { label: 'Volatility of project value', unit: '%/√yr' }, y: { label: 'Option value', unit: MU },
        series: [{ name: 'Value of waiting', x: sgs.map((s) => round(100 * s)), y: sw.map((o) => mil(o.d)) }, { name: 'Expand', x: sgs.map((s) => round(100 * s)), y: sw.map((o) => mil(o.e)) }, { name: 'Abandon', x: sgs.map((s) => round(100 * s)), y: sw.map((o) => mil(o.a)) }],
        explain: 'Flexibility is worth more the less certain the future: every curve rises with volatility. Read how far the conclusion depends on the volatility assumed, the least observable input.' },
      { type: 'bar', title: 'From static to expanded NPV', categories: ['Static NPV', 'Expand', 'Abandon', 'Interaction', 'Expanded NPV', 'Right to wait'], x: { label: '' }, y: { label: 'Value', unit: MU },
        series: [{ name: 'Value', y: [staticNpv, exp.optionValue, abn.optionValue, interaction, buildNowValue, defer].map(mil) }],
        explain: 'Static NPV, the contribution of each asset option and their interaction, the resulting value of building now, and beside it the value of keeping the investment decision open. The taller of the last two is the better course.' },
      { type: 'bar', title: 'Staged development: commitments and value', categories: [...ph.map((s) => `${s.name} (yr ${fmt(s.time)})`), 'Programme value'], x: { label: 'Stage' }, y: { label: 'Amount', unit: MU },
        series: [{ name: 'PV of commitment', y: [...pvPhases.map(mil), 0] }, { name: 'Value', y: [...ph.map(() => 0), mil(staged.value)] }],
        explain: 'Present value of each gated commitment and the value of the programme today. Small early stages that buy the right to a large later decision are where staging pays.' },
      { type: 'bar', title: 'First decision: expected value by strategy', categories: ['Commit now', 'Appraise first', 'Walk away'], x: { label: 'Strategy' }, y: { label: 'Expected value', unit: MU }, series: [{ name: 'Expected value', y: [mil(evCommit), mil(evAppraise), 0] }],
        explain: 'Expected value of each first move. The appraisal wins when the chance of avoiding a bad investment outweighs its cost and the delay.' },
    ];
    const tables = [
      { title: 'Option summary', columns: [{ key: 'opt', label: 'Option' }, { key: 'value', label: 'Value', unit: MU }, { key: 'share', label: 'Share of investment', unit: '%' }, { key: 'basis', label: 'Basis' }],
        rows: [['Right to wait (opportunity value)', defer, `${american ? 'American' : 'European'} call, ${N} steps, ${fmt(inp.window)} yr`], ['Value of waiting', waiting, 'Opportunity value − max(static NPV, 0)'], ['Black–Scholes reference', bs.call, 'European call with yield q'],
          ...(lsm ? [['Least-squares Monte Carlo', lsm.value, `${lsm.paths} paths, 50 steps`]] : []), ['Expand', exp.optionValue, `× ${fmt(inp.expFactor)} for ${fM(inp.expCost)}, ${fmt(inp.opWindow)} yr`], ['Abandon', abn.optionValue, `Salvage ${fM(inp.salvage)}, ${fmt(inp.opWindow)} yr`],
          ['Expand + abandon together', both.optionValue, 'One lattice (interaction included)'], ['Stage gates', staged.flexibilityValue, `${ph.length} stages to year ${fmt(Tst)}`], ['Perfect information', info.evpi, 'Upper bound on any appraisal']]
          .map(([opt, v, basis]) => ({ opt, value: mil(v), share: round(100 * v / I, 4), basis })),
        explain: 'Every flexibility valued here with its basis and its size relative to the investment. Values on different rows are not additive unless stated.' },
      { title: 'Sensitivities of the European reference (Black–Scholes Greeks)', columns: [{ key: 'g', label: 'Measure' }, { key: 'v', label: 'Value' }, { key: 'unit', label: 'Unit' }, { key: 'read', label: 'Reading' }],
        rows: [{ g: 'Delta', v: round(bs.greeks.call.delta, 5), unit: '-', read: 'Change in option value per unit of project value' }, { g: 'Vega', v: mil(bs.greeks.call.vega / 100), unit: `${MU} per vol point`, read: 'Change per one percentage point of volatility' },
          { g: 'Theta', v: mil(bs.greeks.call.theta), unit: `${MU}/yr`, read: 'Change per year as the window shortens (negative = decays)' }, { g: 'Rho', v: mil(bs.greeks.call.rho / 100), unit: `${MU} per rate point`, read: 'Change per percentage point of the risk-free rate' }],
        explain: 'Local sensitivities of the option to its inputs, for hand checks and for judging which assumption deserves scrutiny. They are exact for the European case only.' },
      { title: 'Staged development', columns: [{ key: 'name', label: 'Stage' }, { key: 'time', label: 'Decision time', unit: 'yr' }, { key: 'cost', label: 'Commitment', unit: MU }, { key: 'pv', label: 'PV at risk-free rate', unit: MU }, { key: 'cum', label: 'Cumulative PV', unit: MU }],
        rows: ph.map((s, k) => ({ name: s.name, time: s.time, cost: mil(s.cost), pv: mil(pvPhases[k]), cum: mil(pvPhases.slice(0, k + 1).reduce((a, b) => a + b, 0)) })),
        explain: 'Each gate with the capital it commits. The cumulative column is the exposure reached if every gate is passed; the option value comes from being able to stop before the large rows.' },
      { title: 'Decision tree policy', columns: [{ key: 'node', label: 'Decision point' }, { key: 'choice', label: 'Best choice' }, { key: 'ev', label: 'Expected value', unit: MU }],
        rows: tree.policy.map((s) => ({ node: s.node, choice: s.choice, ev: mil(s.expectedValue) })),
        explain: 'The optimal choice at every decision point reached under the best strategy, with the expected value from that point on.' },
    ];

    const vega1 = bs.greeks.call.vega / 100;
    const scrap = inp.steel * 0.9 * inp.scrapPrice;
    const advice = [
      adv('economic', buildNowValue < 0 && defer <= 1e-9 * I ? 'fail' : waitBetter ? 'warn' : 'ok',
        `Static NPV ${fM(staticNpv)}; opportunity value with the right to wait ${fM(defer)} (value of waiting ${fM(waiting)}); expanded NPV of building now ${fM(buildNowValue)} with ${fM(both.optionValue)} of expansion and abandonment value.${trigger != null ? ` Investment trigger ${fM(trigger)} against a project value of ${fM(V)}.` : ''}`,
        waitBetter ? `Keep the option open: waiting is worth ${fM(defer - buildNowValue)} more than committing today. Invest when project value reaches ${trigger != null ? fM(trigger) : 'the trigger'}, or when the window is about to close with V above ${fM(I)}.`
          : buildNowValue >= 0 ? `Commit: building now is worth at least as much as waiting (${fM(buildNowValue)} against ${fM(defer)}). Secure the expansion and exit rights in the contracts and consents, since ${fM(both.optionValue)} of the value rests on them.`
            : 'Do not commit on these figures: neither the cash flows nor the flexibility recover the investment. Hold the opportunity only if holding it is free.',
        'Expanded NPV rule, Trigeorgis (1996) ch. 4; optimal investment timing, Dixit & Pindyck (1994) ch. 5; the exercise price should be a Class 3 or better estimate, AACE RP 18R-97.'),
      adv('engineering', 'info', `Volatility of ${pct(sigma)} and a cash yield of ${inp.leakage} %/yr drive the result: one more point of volatility adds about ${fM(vega1)} to the opportunity value, and the option to expand is worth ${fM(exp.optionValue)}.`,
        `Spend up to ${fM(exp.optionValue)} on enabling works that make expansion possible (spare J-tubes and switchgear bays, export-cable and foundation capacity, consent envelope) and no more. Narrow the volatility with site data: yield measurement and ground investigation are worth up to ${fM(info.evpi)} here. The leakage equals annual free cash flow over value, so it rises with yield and availability from the energy analyses.`,
        'Design for flexibility, de Neufville & Scholtes (2011), Flexibility in Engineering Design; value of information, Clemen & Reilly, Making Hard Decisions, ch. 12.',
        `The abandonment value presumes the asset can be removed and its materials recovered: ${fmt(inp.steel)} t of steel at 90 % recovery is about ${fM(scrap)} of scrap and displaces up to ${fmt(inp.steel * 0.9 * CO2_STEEL)} t CO₂e of primary steel (screening factor ${CO2_STEEL} t CO₂e/t). A funded decommissioning liability is what makes the ${fM(inp.salvage)} salvage figure real; staging also avoids fabricating steel for a project that is later stopped.`),
      adv('economic', staged.value > 0 ? 'info' : 'warn', `Developed in ${ph.length} stages the programme is worth ${fM(staged.value)} against a committed-in-advance NPV of ${fM(staged.staticNpv)}; the gates add ${fM(staged.flexibilityValue)}. Best first decision in the tree: "${tree.policy[0].choice}" (${fM(tree.expectedValue)}); the appraisal nets ${fM(vof.valueOfFlexibility)} and perfect information would be worth ${fM(info.evpi)}.`,
        staged.value > 0 ? `Fund the first stage (${fM(c0)}) and re-value at each gate with the updated estimate and price outlook. ${vof.valueOfFlexibility > 0 ? `Run the appraisal: it pays for itself by ${fM(vof.valueOfFlexibility)}.` : `Skip the appraisal at ${fM(inp.appraisalCost)}: it does not change the decision enough to pay; it would have to cost less than ${fM(appraisalBreakEven)} to be worth running.`}`
          : 'The first stage is not worth funding at these values: reduce the early commitment, or wait for a higher project value.',
        'Stage-gate valuation as a sequential compound option, Copeland & Antikarov (2001) ch. 6; AACE RP 17R-97 (estimate maturity by gate); decision-tree rollback, Raiffa (1968).'),
    ];
    if (Math.abs(p * inp.upFactor + (1 - p) * inp.downFactor - 1) > 0.05) advice.push(adv('economic', 'warn', `The decision-tree outcomes average ${fmt(p * inp.upFactor + (1 - p) * inp.downFactor)} × V, not 1.0, so the tree and the lattice start from different expected project values.`,
      'Make the two outcome factors and the probability consistent with the project value entered, or accept that the tree answers a different question.', 'Consistency of scenario and distribution inputs, IEC 31010:2019 §B.2.'));

    const notes = [`Real terms at ${inp.baseYear}. Risk-free rate ${inp.riskFree} % effective (${fmt(100 * r)} % continuous); leakage ${inp.leakage} %/yr continuous; volatility ${pct(sigma)} (${inp.volSource === 'input' ? 'entered' : 'consolidated approach'}).`,
      `Sunk development cost of ${fM(inp.sunk)} is excluded: only future spend enters the option values.`,
      `Lattice: ${N} steps for the deferral option, ${Nop} for the asset options, ${nSt} for the staged programme${N > inp.steps || Nop > inp.steps ? ' (raised above the entered value to keep the risk-neutral probability inside 0–1)' : ''}.`,
      waitBetter ? 'Waiting governs: the right to defer is worth more than building now with the asset options.' : 'Building now is worth at least as much as holding the option to wait.'];
    if (consNote) notes.push(consNote);
    const publish = { 'econ.optionValueDefer': defer, 'econ.expandedNpv': buildNowValue, 'econ.stagedValue': staged.value, 'econ.valueOfInformation': info.evpi };
    if (trigger != null) publish['econ.investmentTrigger'] = trigger;
    return { kpis, checks, plots, tables, advice, notes, publish };
  },
  verify: [
    // Hull's standard example: S = K = 100, r = 5 % continuous (effective e^0.05 − 1 = 5.12711 %), σ = 20 %, T = 1, no yield
    // d1 = (0 + (0.05 + 0.02))/0.2 = 0.35, d2 = 0.15; c = 100 N(0.35) − 100 e^−0.05 N(0.15) = 63.6831 − 53.2325 = 10.4506.
    // With no leakage an American call is never exercised early, so the lattice converges to the same value.
    { name: 'Black–Scholes call, S = K = 100, r = 5 %, σ = 20 %, T = 1', source: 'Hull, Options, Futures and Other Derivatives: c = S N(d1) − K e^(−rT) N(d2) = 10.4506 with d1 = 0.35, d2 = 0.15',
      inputs: { V: 100e6, invest: 100e6, volSource: 'input', sigma: 20, riskFree: 5.127109637602412, leakage: 0, window: 1, steps: 400, lsmPaths: 0 },
      expect: [{ kpi: 'bsValue', value: 10.4506, tol: 1e-4 }, { kpi: 'deferValue', value: 10.4506, tol: 2e-3 }, { kpi: 'staticNpv', value: 0, abs: 1e-9 }] },
    // European abandonment = European put with strike = salvage on project value: S = 100, K = 90, r = 5 %, σ = 30 %, T = 2
    // d1 = [ln(100/90) + (0.05 + 0.045)·2]/(0.3·√2) = 0.696172, d2 = 0.271908; p = 90 e^−0.1 N(−d2) − 100 N(−d1) = 31.9917 − 24.3161 = 7.676
    { name: 'European abandonment option equals a Black–Scholes put', source: 'Put on project value with the salvage as strike (Myers & Majd 1990): p = K e^(−rT) N(−d2) − S N(−d1) = 7.676 for S = 100, K = 90, r = 5 %, σ = 30 %, T = 2',
      inputs: { V: 100e6, invest: 100e6, volSource: 'input', sigma: 30, riskFree: 5.127109637602412, leakage: 0, salvage: 90e6, opWindow: 2, exercise: 'european', steps: 400, lsmPaths: 0 },
      expect: [{ kpi: 'abandonValue', value: 7.676, tol: 3e-3 }] },
    // single-stage programme with the whole cost at year 1 is the European call of the first case
    { name: 'One-stage programme reduces to a call', source: 'Limiting case of the compound option: one gate at T = 1 with cost 100 on V = 100 is the call above, 10.4506',
      inputs: { V: 100e6, invest: 100e6, volSource: 'input', sigma: 20, riskFree: 5.127109637602412, leakage: 0, window: 1, steps: 400, lsmPaths: 0, phases: [{ name: 'Build', time: 1, cost: 100e6 }] },
      expect: [{ kpi: 'stagedValue', value: 10.4506, tol: 2e-3 }] },
    // tree: commit = 0.6·1.3·100 + 0.4·0.5·100 − 100 = −2; appraise (free, no delay) = 0.6·(130 − 100) + 0.4·0 = 18; walk = 0 → 18.
    // EVPI = 0.6·30 + 0.4·0 − max(−2, 0) = 18
    { name: 'Decision tree: appraise before committing', source: 'Hand rollback: commit 0.6 × 130 + 0.4 × 50 − 100 = −2; appraise 0.6 × 30 = 18; walk 0 → 18; EVPI = 18 − 0',
      inputs: { V: 100e6, invest: 100e6, pGood: 60, upFactor: 1.3, downFactor: 0.5, appraisalCost: 0, appraisalYears: 0, lsmPaths: 0 },
      expect: [{ kpi: 'treeValue', value: 18, tol: 1e-9 }, { kpi: 'infoValue', value: 18, tol: 1e-9 }, { kpi: 'appraisalGain', value: 18, tol: 1e-9 }] },
    // no operating cost: project value moves one-for-one with price, so the consolidated volatility equals the price volatility
    { name: 'Consolidated volatility without operating leverage', source: 'Copeland & Antikarov (2001) ch. 9: with zero cost share the log return on value is the log return on price, σ_V = σ_P = 15 %',
      inputs: { volSource: 'consolidated', priceVol: 15, costShare: 0, lsmPaths: 0 }, expect: [{ kpi: 'sigmaUsed', value: 15, tol: 0.01 }] },
  ],
};
