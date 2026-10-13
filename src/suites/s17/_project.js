// Generic project definition shared by the appraisal modules, and a fast after-tax cash-flow kernel.
// The kernel reproduces the unlevered free-cash-flow logic of engine/econ/cashflow.js buildProjectCashflow for
// the feature set used here (phased capex, indexed or fixed price, fixed + variable opex with real escalation,
// royalty, straight-line tax depreciation with loss carry-forward, decommissioning in the year after the last
// operating year, whole-year start-up delay). It exists because sampling studies need thousands of evaluations
// per run; the engine model is used for the single deterministic cases.
import { cost, need, PLACEHOLDER, CUR, pv } from './_common.js';

/**
 * Form fields of the generic project. `d` overrides defaults; `skip` removes keys a module does not need.
 * Role test: opening financial state (base year, terms basis, inflation expectation, sunk cost) is `ic`;
 * constraints (hurdle rate, price contract, tax regime, operating life) are `bc`; prices and costs are `cost`.
 */
export function projectInputs(d = {}, skip = []) {
  const v = (k, def) => (d[k] !== undefined ? d[k] : def);
  const list = [
    { key: 'baseYear', label: 'Price base year (time zero)', unit: 'yr', type: 'int', default: v('baseYear', 2026), min: 1990, max: 2100, group: 'Opening financial state', role: 'ic',
      help: 'Calendar year of row 0. All constant-money inputs are in money of this year.' },
    { key: 'basis', label: 'Terms of the discount rate and results', type: 'select', default: v('basis', 'real'), group: 'Opening financial state', role: 'ic',
      options: [{ value: 'real', label: 'Real (constant base-year money)' }, { value: 'nominal', label: 'Nominal (money of the day)' }],
      help: 'The discount rate must be stated in the same terms. Real and nominal give the same NPV when the rates are linked by the Fisher relation.' },
    { key: 'inflation', label: 'General inflation', unit: '%/yr', type: 'number', default: v('inflation', 2), min: -2, max: 30, group: 'Opening financial state', role: 'ic',
      help: 'Expectation used to move between real and nominal money.' },
    { key: 'sunk', label: 'Cost already incurred (sunk)', unit: CUR, type: 'number', default: v('sunk', 2.5e7), min: 0, max: 1e12, group: 'Opening financial state', role: 'ic',
      help: 'Historical spend. Excluded from the decision cash flow; kept only in the tax depreciation basis.' },
    { key: 'ratedMW', label: 'Installed capacity', unit: 'MW', type: 'number', default: v('ratedMW', 500), min: 0.1, max: 20000, group: 'Production', role: 'operation', bind: 'energy.ratedPower',
      help: 'From the energy-yield suite when available.' },
    { key: 'aep', label: 'Net annual production, first year', unit: 'GWh/yr', type: 'number', default: v('aep', 2150), min: 0.001, max: 2e5, group: 'Production', role: 'operation', bind: 'energy.aep',
      uncertain: { dist: 'normal', cov: 0.08 }, help: 'Net of wake, electrical and availability losses. From the energy-yield suite when available.' },
    { key: 'degradation', label: 'Annual output degradation', unit: '%/yr', type: 'number', default: v('degradation', 0.3), min: 0, max: 10, group: 'Production', role: 'operation' },
    { key: 'constructionYears', label: 'Construction period', unit: 'yr', type: 'int', default: v('constructionYears', 3), min: 1, max: 10, group: 'Timeline', role: 'operation', design: true,
      help: 'Whole years before first production; capex is spread evenly unless the module has a phasing table.' },
    { key: 'opLife', label: 'Operating life', unit: 'yr', type: 'int', default: v('opLife', 25), min: 1, max: 60, group: 'Timeline', role: 'bc', bind: 'life.designYears', design: true,
      help: 'Operational constraint from the design life or the fatigue life supplied by the engineering suites.' },
    cost({ key: 'capex', label: 'Capital expenditure (base-year money)', unit: CUR, default: v('capex', 1.5e9), min: 1, max: 1e12, group: 'Costs', bind: 'econ.capex', design: true,
      uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'From the capital cost estimate when available.' }),
    cost({ key: 'opexFixed', label: 'Fixed operating cost (base-year money)', unit: `${CUR}/yr`, default: v('opexFixed', 4.5e7), min: 0, max: 1e11, group: 'Costs', bind: 'econ.opexAnnual',
      uncertain: { dist: 'normal', cov: 0.15 }, help: 'From the operating cost build-up when available.' }),
    cost({ key: 'opexVar', label: 'Variable operating cost', unit: `${CUR}/MWh`, default: v('opexVar', 1.5), min: 0, max: 1e4, group: 'Costs' }),
    cost({ key: 'decom', label: 'Decommissioning cost (base-year money)', unit: CUR, default: v('decom', 1.2e8), min: 0, max: 1e12, group: 'Costs', bind: 'econ.decomCost',
      help: 'Terminal condition: spent in the year after the last operating year.' }),
    { key: 'priceType', label: 'Price contract', type: 'select', default: v('priceType', 'indexed'), group: 'Revenue', role: 'bc',
      options: [{ value: 'indexed', label: 'Inflation-indexed contract price' }, { value: 'fixed', label: 'Fixed nominal contract price' }],
      help: 'Revenue boundary: an indexed price keeps its real value; a fixed nominal price loses real value with inflation.' },
    cost({ key: 'price', label: 'Unit price (base-year money)', unit: `${CUR}/MWh`, default: v('price', 116), min: 0.01, max: 1e5, group: 'Revenue',
      uncertain: { dist: 'normal', cov: 0.15 }, help: 'Tariff, strike or expected market price.' }),
    { key: 'taxRate', label: 'Corporate tax rate', unit: '%', type: 'number', default: v('taxRate', 25), min: 0, max: 90, group: 'Fiscal regime', role: 'bc' },
    { key: 'royalty', label: 'Royalty or levy on revenue', unit: '%', type: 'number', default: v('royalty', 0), min: 0, max: 90, group: 'Fiscal regime', role: 'bc' },
    { key: 'depLife', label: 'Tax depreciation life (straight line)', unit: 'yr', type: 'int', default: v('depLife', 15), min: 1, max: 60, group: 'Fiscal regime', role: 'bc' },
    { key: 'discountRate', label: 'Discount rate (hurdle rate)', unit: '%/yr', type: 'number', default: v('discountRate', 6), min: -5, max: 50, group: 'Financial constraints', role: 'bc',
      help: `Minimum acceptable return, in the terms selected above. ${PLACEHOLDER}` },
  ];
  return list.filter((i) => !skip.includes(i.key));
}

/** Form values (engineering units) → kernel parameters (currency, MWh, fractions). */
export function projectParams(inp) {
  const nOp = Math.round(inp.opLife); const nc = Math.round(inp.constructionYears);
  need(nOp >= 1 && nc >= 1, 'The project needs at least one construction year and one operating year.');
  const infl = (inp.inflation ?? 0) / 100;
  const rate = inp.discountRate / 100;
  need(rate > -1 && infl > -1, 'Rates must be above −100 %.');
  const deg = (inp.degradation ?? 0) / 100;
  const e0 = inp.aep * 1e3; // GWh → MWh
  const prod = new Float64Array(nOp);
  for (let k = 0; k < nOp; k++) prod[k] = e0 * (1 - deg) ** k;
  const phasing = new Float64Array(nc).fill(1 / nc);
  return {
    nc, nOp, capex: inp.capex, phasing, sunk: inp.sunk ?? 0, prod, price: inp.price, indexed: (inp.priceType ?? 'indexed') === 'indexed',
    opexFixed: inp.opexFixed, opexVar: inp.opexVar ?? 0, opexEsc: 0, decom: inp.decom ?? 0, taxRate: (inp.taxRate ?? 0) / 100,
    depLife: Math.round(inp.depLife ?? nOp), royalty: (inp.royalty ?? 0) / 100, infl, rate, real: (inp.basis ?? 'real') === 'real', delay: 0,
    baseYear: inp.baseYear ?? 0,
  };
}

/** Index and discount tables of a parameter set (inflation index, discount factors at the stated and the real rate). */
function tables(p) {
  if (p.cache) return p.cache;
  const T = p.nc + p.nOp + 16;
  const I = new Float64Array(T); const D = new Float64Array(T); const DR = new Float64Array(T);
  const g = 1 + p.infl; const rr = p.real ? p.rate : (1 + p.rate) / g - 1;
  for (let t = 0; t < T; t++) { I[t] = g ** t; D[t] = (1 + p.rate) ** -t; DR[t] = (1 + rr) ** -t; }
  Object.defineProperty(p, 'cache', { value: { I, D, DR, T }, enumerable: true, writable: true });
  return p.cache;
}

/**
 * Allocation-free evaluation for sampling: writes the reporting-terms free cash flow into `buf` and returns the
 * NPV; `res` receives { npv, rows, lcoe } (LCOE pre-tax, real). Same logic as dcf() for a whole-year delay.
 */
export function fastFlow(p, f, delay, buf, res) {
  const { I, D, DR, T } = tables(p);
  const fp = f?.price ?? 1; const fq = f?.prod ?? 1; const fc = f?.capex ?? 1; const fo = f?.opex ?? 1;
  const first = p.nc + delay; const N = first + p.nOp;
  if (N >= T) throw new Error('Start-up delay is beyond the range of the cash-flow model (about 14 years).');
  let basis = p.sunk || 0; let npv = 0; let pvC = 0; let pvE = 0;
  for (let t = 0; t < p.nc; t++) {
    const c = p.capex * fc * p.phasing[t] * I[t]; basis += c;
    const x = p.real ? -c / I[t] : -c;
    buf[t] = x; npv += x * D[t]; pvC += c / I[t] * DR[t];
  }
  for (let t = p.nc; t < first; t++) buf[t] = 0;
  const depYr = basis / p.depLife; let loss = 0;
  for (let k = 0; k < p.nOp; k++) {
    const t = first + k; const It = I[t];
    const q = p.prod[k] * fq;
    const rev = q * p.price * fp * (p.indexed ? It : 1);
    const ox = (p.opexFixed + p.opexVar * q) * fo * It;
    const ti = rev * (1 - p.royalty) - ox - (k < p.depLife ? depYr : 0);
    let tax = 0;
    if (ti < 0) loss -= ti; else { const use = loss < ti ? loss : ti; loss -= use; tax = p.taxRate * (ti - use); }
    let x = rev * (1 - p.royalty) - ox - tax; if (p.real) x /= It;
    buf[t] = x; npv += x * D[t]; pvC += ox / It * DR[t]; pvE += q * DR[t];
  }
  buf[N] = p.real ? -p.decom : -p.decom * I[N];
  npv += buf[N] * D[N]; pvC += p.decom * DR[N];
  if (res) { res.npv = npv; res.rows = N + 1; res.lcoe = pvE > 0 ? pvC / pvE : NaN; }
  return npv;
}

/**
 * After-tax unlevered cash flow. Evaluated in nominal money (so depreciation and fixed prices are right) and
 * deflated to base-year money when p.real. Returns per-row arrays (row 0 = time zero) in the reporting terms.
 * f = { price, prod, capex, opex } multipliers for sensitivity and sampling.
 */
export function dcf(p, f) {
  const fp = f?.price ?? 1; const fq = f?.prod ?? 1; const fc = f?.capex ?? 1; const fo = f?.opex ?? 1;
  const delay = Math.max(0, Math.round(p.delay || 0));
  const first = p.nc + delay; const N = first + p.nOp; const rows = N + 1;
  const g = 1 + p.infl; const ge = 1 + (p.opexEsc || 0);
  const revenue = new Float64Array(rows); const opex = new Float64Array(rows); const capex = new Float64Array(rows);
  const tax = new Float64Array(rows); const dep = new Float64Array(rows); const decom = new Float64Array(rows);
  const royalty = new Float64Array(rows); const fcff = new Float64Array(rows); const energy = new Float64Array(rows);
  let basis = p.sunk || 0;
  for (let t = 0; t < p.nc; t++) { capex[t] = p.capex * fc * p.phasing[t] * g ** t; basis += capex[t]; }
  const depYr = basis / p.depLife;
  let loss = 0;
  for (let k = 0; k < p.nOp; k++) {
    const t = first + k; const I = g ** t;
    const q = p.prod[k] * fq;
    energy[t] = q;
    revenue[t] = q * p.price * fp * (p.indexed ? I : 1);
    royalty[t] = revenue[t] * p.royalty;
    opex[t] = (p.opexFixed + p.opexVar * q) * fo * I * ge ** t;
    dep[t] = k < p.depLife ? depYr : 0;
    const ti = revenue[t] - royalty[t] - opex[t] - dep[t];
    if (ti < 0) loss -= ti;
    else { const use = Math.min(loss, ti); loss -= use; tax[t] = p.taxRate * (ti - use); }
  }
  decom[N] = p.decom * g ** N;
  for (let t = 0; t < rows; t++) fcff[t] = revenue[t] - royalty[t] - opex[t] - tax[t] - capex[t] - decom[t];
  if (p.real) {
    for (let t = 0; t < rows; t++) {
      const I = g ** t;
      revenue[t] /= I; royalty[t] /= I; opex[t] /= I; tax[t] /= I; dep[t] /= I; capex[t] /= I; decom[t] /= I; fcff[t] /= I;
    }
  }
  return { rows, first, last: N - 1, decomRow: N, revenue, royalty, opex, tax, dep, capex, decom, fcff, energy, npv: pv(fcff, p.rate), unusedLoss: loss };
}

const SCRATCH = new Float64Array(160);
/** NPV with a fractional start-up delay: linear interpolation between the neighbouring whole-year cases. */
export function npvDelayed(p, f, delayYears) {
  const d = Math.max(0, delayYears); const d0 = Math.floor(d); const w = d - d0;
  const buf = p.nc + p.nOp + d0 + 3 <= SCRATCH.length ? SCRATCH : new Float64Array(p.nc + p.nOp + d0 + 3);
  const a = fastFlow(p, f, d0, buf);
  return w < 1e-9 ? a : a + w * (fastFlow(p, f, d0 + 1, buf) - a);
}

/** Real and nominal discount rates of the project (Fisher). */
export function rates(p) {
  const real = p.real ? p.rate : (1 + p.rate) / (1 + p.infl) - 1;
  return { real, nominal: (1 + real) * (1 + p.infl) - 1 };
}

/**
 * Pre-tax real levelised cost: PV(capex + opex + decommissioning) / PV(energy), constant money, real rate
 * (IEA/NREL convention). `cf` is a dcf() result in either terms.
 */
export function lcoeOf(p, cf) {
  const r = rates(p).real; const g = 1 + p.infl;
  let c = 0; let e = 0;
  for (let t = 0; t < cf.rows; t++) {
    const defl = p.real ? 1 : g ** t;
    const d = (1 + r) ** -t;
    c += (cf.capex[t] + cf.opex[t] + cf.decom[t]) / defl * d;
    e += cf.energy[t] * d;
  }
  return { lcoe: e > 0 ? c / e : NaN, pvCost: c, pvEnergy: e };
}
