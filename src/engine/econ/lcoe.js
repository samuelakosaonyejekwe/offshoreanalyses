// Levelised cost measures: LCOE (discounted, real/nominal, pre/post-tax, fixed-charge-rate), storage, hydrogen,
// delivered energy, generic levelised cost of service, avoided cost (LACE-style), component breakdown,
// sensitivity grids, required strike price for a target equity IRR, blended contract revenue, learning projection.
// Costs in currency, energy in MWh (or any consistent unit), so LCOE is currency per MWh.

import { requireNumber, requireSeries, discountFactors, factorAP, realToNominal } from './tvm.js';
import { npv, solveFor, buildProjectCashflow } from './cashflow.js';
import { financeProject } from './finance.js';
import { learningExponent } from './cost.js';

const pvOf = (series, rate, o) => (series ? npv(series, rate, o) : 0);
const asSeries = (v, n, name) => {
  if (v === undefined) return new Array(n).fill(0);
  if (typeof v === 'number') { const a = new Array(n).fill(v); a[0] = 0; return a; }   // a scalar means "every year after time zero"
  requireSeries(name, v);
  const a = new Array(n).fill(0);
  for (let k = 0; k < Math.min(n, v.length); k++) a[k] = v[k];
  return a;
};

/**
 * Levelised cost of energy by the discounting method: LCOE = Σ (capex_t + opex_t + fuel_t + decom_t − salvage_t)
 * DF_t / Σ E_t DF_t. Index 0 is time zero. Use a real discount rate with constant-money costs for the real LCOE;
 * a nominal rate with nominal costs gives the nominal LCOE (the flat nominal price). Pre-tax unless tax terms are
 * given, in which case the post-tax form is used: [PV(capex) − T PV(depreciation) + (1 − T) PV(opex + fuel +
 * decom)] / [(1 − T) PV(E)] with a post-tax discount rate.
 * A scalar opex/fuel/energy is applied to every period after time zero.
 * Ref: IEA/NEA (2020), Projected Costs of Generating Electricity, ch. 2 and Annex; Short, Packey & Holt (1995),
 * A Manual for the Economic Evaluation of Energy Efficiency and Renewable Energy Technologies, NREL/TP-462-5173.
 * @param {{capex: number[] | number, opex?: number[] | number, fuel?: number[] | number, decom?: number[],
 *   salvage?: number[], energy: number[] | number, rate: number, years?: number, midYear?: boolean,
 *   taxRate?: number, depreciation?: number[]}} o per-period series (currency; energy in MWh); a scalar capex is
 *   spent at time zero; `years` sets the horizon when only scalars are given
 * @returns {{lcoe: number, pvCost: number, pvEnergy: number, components: {capex: number, opex: number, fuel: number,
 *   decom: number, salvage: number, taxShield: number}, shares: Record<string, number>, basis: string}} lcoe and
 *   components in currency per MWh; shares are fractions of LCOE
 */
export function lcoe({ capex, opex, fuel, decom, salvage, energy, rate, years, midYear = false, taxRate, depreciation }) {
  requireNumber('rate', rate, { gt: -1 });
  const lens = [capex, opex, fuel, decom, salvage, energy].filter((s) => s && typeof s !== 'number').map((s) => s.length);
  const n = Math.max(years !== undefined ? years + 1 : 0, ...lens);
  if (!(n >= 2)) throw new RangeError('give at least one series or `years` to set the horizon');
  const cx = typeof capex === 'number' ? [requireNumber('capex', capex, { min: 0 }), ...new Array(n - 1).fill(0)] : asSeries(capex, n, 'capex');
  const ox = asSeries(opex, n, 'opex'); const fx = asSeries(fuel, n, 'fuel'); const dx = asSeries(decom, n, 'decom');
  const sx = asSeries(salvage, n, 'salvage'); const en = asSeries(energy, n, 'energy');
  const o = { midYear };
  const pvE = npv(en, rate, o);
  if (!(pvE > 0)) throw new RangeError('discounted energy must be > 0');
  const pv = { capex: npv(cx, rate, o), opex: npv(ox, rate, o), fuel: npv(fx, rate, o), decom: npv(dx, rate, o),
    salvage: -npv(sx, rate, o), taxShield: 0 };
  let denom = pvE; let basis = 'pre-tax';
  if (taxRate !== undefined && taxRate > 0) {
    requireNumber('taxRate', taxRate, { gt: 0, lt: 1 });
    basis = 'post-tax';
    pv.taxShield = -taxRate * pvOf(depreciation, rate, o);
    for (const k of ['opex', 'fuel', 'decom', 'salvage']) pv[k] *= 1 - taxRate;
    denom = (1 - taxRate) * pvE;
  }
  const pvCost = Object.values(pv).reduce((s, v) => s + v, 0);
  const components = {}; const shares = {};
  for (const k of Object.keys(pv)) { components[k] = pv[k] / denom; shares[k] = pvCost !== 0 ? pv[k] / pvCost : 0; }
  return { lcoe: pvCost / denom, pvCost, pvEnergy: pvE, components, shares, basis };
}

/**
 * Real and nominal LCOE from constant-money costs: real LCOE uses the real discount rate; the nominal LCOE is the
 * flat nominal price with the same present value, obtained by inflating the costs and discounting costs and energy
 * at the nominal rate (Fisher).
 * Ref: Short, Packey & Holt (1995), NREL/TP-462-5173 (real vs nominal levelised cost).
 * @param {object} o arguments of lcoe() in constant money with `rate` = real rate, plus `inflation` (fraction/yr)
 * @returns {{real: number, nominal: number, nominalRate: number}} currency per MWh
 */
export function lcoeRealNominal(o) {
  const { inflation, ...rest } = o;
  requireNumber('inflation', inflation, { gt: -1 });
  const real = lcoe(rest);
  const nominalRate = realToNominal(rest.rate, inflation);
  const n = Math.max(...[rest.capex, rest.opex, rest.fuel, rest.decom, rest.salvage, rest.energy]
    .filter((s) => s && typeof s !== 'number').map((s) => s.length), rest.years !== undefined ? rest.years + 1 : 0);
  const infl = (v, name) => {
    if (v === undefined) return undefined;
    const a = name === 'capex' && typeof v === 'number' ? [v, ...new Array(n - 1).fill(0)] : asSeries(v, n, name);
    return a.map((x, k) => x * Math.pow(1 + inflation, k));
  };
  const nominal = lcoe({ ...rest, capex: infl(rest.capex, 'capex'), opex: infl(rest.opex, 'opex'), fuel: infl(rest.fuel, 'fuel'),
    decom: infl(rest.decom, 'decom'), salvage: infl(rest.salvage, 'salvage'), depreciation: rest.depreciation, rate: nominalRate });
  return { real: real.lcoe, nominal: nominal.lcoe, nominalRate };
}

/**
 * Simple fixed-charge-rate LCOE: (FCR × overnight capex + fixed O&M) / annual energy + variable O&M + fuel.
 * With FCR = capital-recovery factor (A/P, i, n) this equals the discounted LCOE for level annual costs and energy.
 * Ref: NREL, Simple Levelized Cost of Energy calculator documentation; Short, Packey & Holt (1995).
 * @param {{capex: number, fcr?: number, rate?: number, life?: number, fixedOm: number, annualEnergy: number,
 *   variableOm?: number, fuel?: number}} o capex (currency), FCR (1/yr) or rate + life to use the CRF, fixed O&M
 *   (currency/yr), energy (MWh/yr), variable O&M and fuel (currency/MWh)
 * @returns {{lcoe: number, fcr: number, components: {capex: number, fixedOm: number, variableOm: number,
 *   fuel: number}}} currency per MWh
 */
export function lcoeFixedChargeRate({ capex, fcr, rate, life, fixedOm, annualEnergy, variableOm = 0, fuel = 0 }) {
  requireNumber('capex', capex, { min: 0 }); requireNumber('fixedOm', fixedOm, { min: 0 });
  requireNumber('annualEnergy', annualEnergy, { gt: 0 });
  const f = fcr ?? factorAP(rate, life);
  requireNumber('fcr', f, { gt: 0 });
  const c = { capex: f * capex / annualEnergy, fixedOm: fixedOm / annualEnergy, variableOm, fuel };
  return { lcoe: c.capex + c.fixedOm + c.variableOm + c.fuel, fcr: f, components: c };
}

/**
 * Levelised cost of storage: [capex + Σ (O&M_t + charging cost_t + replacement_t) DF_t] / Σ discharged energy_t
 * DF_t, with discharged energy = usable capacity × cycles × (1 − degradation)^(t−1) and charging energy =
 * discharged / round-trip efficiency.
 * Ref: Schmidt, O., Melchior, S., Hawkes, A. & Staffell, I. (2019), Projecting the future levelized cost of
 * electricity storage technologies, Joule 3(1), 81–100; Lazard LCOS methodology.
 * @param {{capex: number, omPerYear?: number, capacity_MWh: number, cyclesPerYear: number, roundTripEfficiency: number,
 *   chargePrice?: number, degradation?: number, life: number, rate: number, replacements?: {year: number,
 *   cost: number}[], endOfLife?: number}} o capacity = usable energy per cycle (MWh), efficiency 0–1, charging price
 *   (currency/MWh), degradation per year, life (years), end-of-life cost (+) or value (−) at the end
 * @returns {{lcos: number, pvCost: number, pvEnergy: number, components: Record<string, number>}} currency per MWh
 */
export function lcos({ capex, omPerYear = 0, capacity_MWh, cyclesPerYear, roundTripEfficiency, chargePrice = 0,
  degradation = 0, life, rate, replacements = [], endOfLife = 0 }) {
  requireNumber('capex', capex, { min: 0 }); requireNumber('capacity_MWh', capacity_MWh, { gt: 0 });
  requireNumber('cyclesPerYear', cyclesPerYear, { gt: 0 }); requireNumber('roundTripEfficiency', roundTripEfficiency, { gt: 0, max: 1 });
  requireNumber('life', life, { min: 1, integer: true }); requireNumber('degradation', degradation, { min: 0, lt: 1 });
  const df = discountFactors(rate, life + 1);
  let pvE = 0; let pvOm = 0; let pvCh = 0; let pvRep = 0;
  for (let t = 1; t <= life; t++) {
    const e = capacity_MWh * cyclesPerYear * Math.pow(1 - degradation, t - 1);
    pvE += e * df[t]; pvOm += omPerYear * df[t]; pvCh += e / roundTripEfficiency * chargePrice * df[t];
  }
  for (const r of replacements) { requireNumber('replacement year', r.year, { min: 1, max: life, integer: true }); pvRep += r.cost * df[r.year]; }
  const pvEol = endOfLife * df[life];
  const pvCost = capex + pvOm + pvCh + pvRep + pvEol;
  return { lcos: pvCost / pvE, pvCost, pvEnergy: pvE,
    components: { capex: capex / pvE, om: pvOm / pvE, charging: pvCh / pvE, replacement: pvRep / pvE, endOfLife: pvEol / pvE } };
}

/**
 * Levelised cost of hydrogen: [capex + Σ (fixed O&M + electricity + water + stack replacement) DF] / Σ kg H2 DF,
 * with annual output = electrolyser MW × 8760 × load factor × 1000 / specific consumption (kWh/kg).
 * Ref: IEA (2019), The Future of Hydrogen, Annex (LCOH definition); IRENA (2020), Green Hydrogen Cost Reduction.
 * @param {{capex: number, capacityMW: number, loadFactor: number, specificConsumption_kWhPerKg: number,
 *   electricityPrice: number, fixedOmPerYear?: number, waterCostPerKg?: number, stackReplacements?: {year: number,
 *   cost: number}[], life: number, rate: number, degradation?: number}} o electricity price in currency/MWh;
 *   degradation = annual increase of specific consumption (fraction)
 * @returns {{lcoh: number, annualOutput_kg: number, components: Record<string, number>}} currency per kg H2
 */
export function lcoh({ capex, capacityMW, loadFactor, specificConsumption_kWhPerKg, electricityPrice, fixedOmPerYear = 0,
  waterCostPerKg = 0, stackReplacements = [], life, rate, degradation = 0 }) {
  requireNumber('capex', capex, { min: 0 }); requireNumber('capacityMW', capacityMW, { gt: 0 });
  requireNumber('loadFactor', loadFactor, { gt: 0, max: 1 }); requireNumber('specificConsumption_kWhPerKg', specificConsumption_kWhPerKg, { gt: 0 });
  requireNumber('electricityPrice', electricityPrice, { min: 0 }); requireNumber('life', life, { min: 1, integer: true });
  const df = discountFactors(rate, life + 1);
  const mwh = capacityMW * 8760 * loadFactor;
  let pvKg = 0; let pvEl = 0; let pvOm = 0; let pvW = 0; let pvStack = 0;
  for (let t = 1; t <= life; t++) {
    const kg = mwh * 1000 / (specificConsumption_kWhPerKg * Math.pow(1 + degradation, t - 1));
    pvKg += kg * df[t]; pvEl += mwh * electricityPrice * df[t]; pvOm += fixedOmPerYear * df[t]; pvW += kg * waterCostPerKg * df[t];
  }
  for (const r of stackReplacements) { requireNumber('stack replacement year', r.year, { min: 1, max: life, integer: true }); pvStack += r.cost * df[r.year]; }
  const total = capex + pvEl + pvOm + pvW + pvStack;
  return { lcoh: total / pvKg, annualOutput_kg: mwh * 1000 / specificConsumption_kWhPerKg,
    components: { capex: capex / pvKg, electricity: pvEl / pvKg, om: pvOm / pvKg, water: pvW / pvKg, stack: pvStack / pvKg } };
}

/**
 * Levelised cost of delivered energy: generation LCOE grossed up for transmission losses plus the levelised cost
 * of the transmission asset per delivered MWh: LCODE = LCOE_gen / (1 − losses) + [capex_tx CRF + O&M_tx] /
 * [E_gen (1 − losses)].
 * Ref: delivered-cost extension of the LCOE definition, IEA/NEA (2020) ch. 2 (system and grid costs kept explicit).
 * @param {{lcoeGeneration: number, annualEnergy: number, lossesFrac: number, transmissionCapex: number,
 *   transmissionOmPerYear?: number, rate: number, life: number}} o currency/MWh, MWh/yr generated, loss fraction
 *   (0–1), currency, currency/yr, discount rate, transmission asset life (years)
 * @returns {{lcode: number, generation: number, transmission: number, deliveredEnergy: number}} currency per
 *   delivered MWh
 */
export function lcoeDelivered({ lcoeGeneration, annualEnergy, lossesFrac, transmissionCapex, transmissionOmPerYear = 0, rate, life }) {
  requireNumber('lcoeGeneration', lcoeGeneration, { min: 0 }); requireNumber('annualEnergy', annualEnergy, { gt: 0 });
  requireNumber('lossesFrac', lossesFrac, { min: 0, lt: 1 }); requireNumber('transmissionCapex', transmissionCapex, { min: 0 });
  const delivered = annualEnergy * (1 - lossesFrac);
  const generation = lcoeGeneration / (1 - lossesFrac);
  const transmission = (transmissionCapex * factorAP(rate, life) + transmissionOmPerYear) / delivered;
  return { lcode: generation + transmission, generation, transmission, deliveredEnergy: delivered };
}

/**
 * Levelised cost of a service or product (per barrel produced, per tonne handled, per m³ transported):
 * Σ cost_t DF_t / Σ quantity_t DF_t.
 * Ref: unit technical cost definition used in upstream economics; same discounting identity as LCOE.
 * @param {{costs: number[], quantities: number[], rate: number, midYear?: boolean}} o per-period cost (currency) and
 *   quantity (any unit); index 0 = time zero
 * @returns {{levelisedCost: number, pvCost: number, pvQuantity: number}} currency per unit
 */
export function levelisedCostOfService({ costs, quantities, rate, midYear = false }) {
  const pvCost = npv(costs, rate, { midYear }); const pvQuantity = npv(quantities, rate, { midYear });
  if (!(pvQuantity > 0)) throw new RangeError('discounted quantity must be > 0');
  return { levelisedCost: pvCost / pvQuantity, pvCost, pvQuantity };
}

/**
 * Levelised avoided cost of energy (LACE-style) and net value: LACE = Σ (marginal energy price_t × E_t + capacity
 * payment × capacity credit × capacity) DF_t / Σ E_t DF_t; net value = LACE − LCOE (positive = economically
 * attractive against the displaced system cost).
 * Ref: US EIA (2013–2023), Levelized Costs of New Generation Resources — LCOE and LACE methodology.
 * @param {{energy: number[], marginalPrice: number[] | number, capacityMW?: number, capacityCredit?: number,
 *   capacityPaymentPerMWyr?: number, rate: number, lcoe?: number}} o energy per period (MWh, index 0 = time zero),
 *   marginal price of displaced energy (currency/MWh), capacity value terms, discount rate
 * @returns {{lace: number, energyValue: number, capacityValue: number, netValue: number | null}} currency per MWh
 */
export function lace({ energy, marginalPrice, capacityMW = 0, capacityCredit = 0, capacityPaymentPerMWyr = 0, rate, lcoe: lc }) {
  requireSeries('energy', energy, { minLength: 2 });
  const df = discountFactors(rate, energy.length);
  let pvE = 0; let pvV = 0; let pvC = 0;
  for (let t = 0; t < energy.length; t++) {
    const p = typeof marginalPrice === 'number' ? marginalPrice : marginalPrice[Math.min(t, marginalPrice.length - 1)];
    pvE += energy[t] * df[t]; pvV += energy[t] * p * df[t];
    if (energy[t] > 0) pvC += capacityMW * capacityCredit * capacityPaymentPerMWyr * df[t];
  }
  if (!(pvE > 0)) throw new RangeError('discounted energy must be > 0');
  const l = (pvV + pvC) / pvE;
  return { lace: l, energyValue: pvV / pvE, capacityValue: pvC / pvE, netValue: lc === undefined ? null : l - lc };
}

/**
 * Sensitivity grid of LCOE against capex multiplier × capacity-factor multiplier (energy scales with capacity
 * factor), and LCOE against discount rate.
 * Ref: sensitivity presentation of LCOE, IEA/NEA (2020) ch. 3.
 * @param {{base: object, capexFactors: number[], capacityFactorFactors: number[], rates?: number[]}} o base =
 *   arguments of lcoe(); factors are multipliers on capex and on energy
 * @returns {{xs: number[], ys: number[], z: number[][], rates: number[], lcoeVsRate: number[], base: number}}
 *   z[iy][ix] with x = capex factor and y = capacity-factor factor (currency per MWh)
 */
export function lcoeSensitivityGrid({ base, capexFactors, capacityFactorFactors, rates = [] }) {
  requireSeries('capexFactors', capexFactors); requireSeries('capacityFactorFactors', capacityFactorFactors);
  const scale = (v, f) => (v === undefined ? undefined : typeof v === 'number' ? v * f : Array.from(v, (x) => x * f));
  const z = Array.from(capacityFactorFactors, (fe) => Array.from(capexFactors, (fc) => lcoe({ ...base,
    capex: scale(base.capex, fc), energy: scale(base.energy, fe) }).lcoe));
  return { xs: Array.from(capexFactors), ys: Array.from(capacityFactorFactors), z, rates: Array.from(rates),
    lcoeVsRate: Array.from(rates, (r) => lcoe({ ...base, rate: r }).lcoe), base: lcoe(base).lcoe };
}

/**
 * Required strike (constant real, inflation-indexed) or flat nominal price that delivers a target equity IRR,
 * found by goal seek on the cash-flow model plus the project-finance model.
 * Ref: goal seek on the flow-to-equity model; Yescombe, Principles of Project Finance.
 * @param {{inputs: object, debt: object, targetEquityIrr: number, indexed?: boolean, lo?: number, hi?: number}} o
 *   inputs as for buildProjectCashflow (its price is replaced), debt as for financeProject, target IRR (fraction),
 *   indexed true = constant real price (default), false = flat nominal
 * @returns {{price: number, equityIrr: number, projectIrr: number | null, minDscr: number, gearing: number}}
 *   price in currency per unit of production
 */
export function requiredStrikePrice({ inputs, debt, targetEquityIrr, indexed = true, lo = 1, hi = 200 }) {
  requireNumber('targetEquityIrr', targetEquityIrr, { gt: -1 });
  const run = (p) => {
    const m = buildProjectCashflow({ ...inputs, price: { type: indexed ? 'indexed' : 'fixed', price: p } });
    return financeProject({ fcff: m.fcff, firstOperatingRow: m.timeline.firstOperatingRow, debt, taxRate: inputs.tax?.rate ?? 0 });
  };
  // equity NPV at the target rate is monotonic in price and defined even where the IRR is not unique
  const f = (p) => npv(run(p).equityCashflows, targetEquityIrr);
  const price = solveFor(f, { lo, hi, min: 0, max: 1e9 });
  const r = run(price);
  return { price, equityIrr: r.equityIrr ?? targetEquityIrr, projectIrr: r.projectIrr, minDscr: r.minDscr, gearing: r.gearing };
}

/**
 * Blended revenue of a generator selling shares of output under a PPA (fixed, escalating), a two-way CfD and the
 * merchant market; shares must sum to at most 1 and the remainder is merchant. Contract terms shorter than the
 * horizon revert to merchant.
 * Ref: revenue-stack valuation; CfD two-way settlement as in the UK CfD standard terms (difference payment =
 * strike − reference price).
 * @param {{energy: number[], merchant: number[] | number, ppa?: {share: number, price: number, escalation?: number,
 *   termYears?: number}, cfd?: {share: number, strike: number, indexation?: number, termYears?: number,
 *   reference?: number[] | number}, rate: number}} o energy per year (MWh, index 0 = first operating year,
 *   discounted from the end of that year), prices in currency/MWh
 * @returns {{revenue: number[], byStream: {ppa: number[], cfd: number[], merchant: number[]}, averagePrice: number,
 *   pv: number, levelisedPrice: number, contractedShareOfPv: number}}
 */
export function blendedRevenue({ energy, merchant, ppa, cfd, rate }) {
  requireSeries('energy', energy); requireNumber('rate', rate, { gt: -1 });
  const sp = ppa?.share ?? 0; const sc = cfd?.share ?? 0;
  if (sp < 0 || sc < 0 || sp + sc > 1 + 1e-12) throw new RangeError('ppa.share + cfd.share must be between 0 and 1');
  const at = (s, k) => (typeof s === 'number' ? s : s[Math.min(k, s.length - 1)]);
  const rev = []; const by = { ppa: [], cfd: [], merchant: [] };
  let pv = 0; let pvE = 0; let pvContract = 0; let sumRev = 0; let sumE = 0;
  for (let k = 0; k < energy.length; k++) {
    const m = at(merchant, k); const e = energy[k];
    const inPpa = ppa && (ppa.termYears === undefined || k < ppa.termYears);
    const inCfd = cfd && (cfd.termYears === undefined || k < cfd.termYears);
    const rp = inPpa ? sp * e * ppa.price * Math.pow(1 + (ppa.escalation ?? 0), k) : 0;
    // two-way CfD: market capture plus the difference payment against the reference price
    const rc = inCfd ? sc * e * (m + cfd.strike * Math.pow(1 + (cfd.indexation ?? 0), k) - at(cfd.reference ?? merchant, k)) : 0;
    const rm = e * m * (1 - (inPpa ? sp : 0) - (inCfd ? sc : 0));
    const d = Math.pow(1 + rate, -(k + 1));
    by.ppa.push(rp); by.cfd.push(rc); by.merchant.push(rm); rev.push(rp + rc + rm);
    pv += (rp + rc + rm) * d; pvE += e * d; pvContract += (rp + rc) * d; sumRev += rp + rc + rm; sumE += e;
  }
  return { revenue: rev, byStream: by, averagePrice: sumE > 0 ? sumRev / sumE : NaN, pv, levelisedPrice: pvE > 0 ? pv / pvE : NaN,
    contractedShareOfPv: pv !== 0 ? pvContract / pv : 0 };
}

/**
 * Learning-rate projection of LCOE: the capex-driven share falls with cumulative installed capacity,
 * LCOE(year) = LCOE0 [share (Q/Q0)^b + (1 − share)], Q = Q0 (1 + growth)^(year − year0), b = log2(1 − LR).
 * Ref: experience-curve projection, Rubin, Azevedo, Jaramillo & Yeh (2015), A review of learning rates for
 * electricity supply technologies, Energy Policy 86; IRENA renewable cost reports for the method.
 * @param {{lcoe0: number, year0: number, year: number, learningRate: number, capacityGrowth: number,
 *   learningShare?: number}} o LCOE at year0 (currency/MWh), learning rate per doubling (fraction), annual growth
 *   of cumulative capacity (fraction), share of LCOE subject to learning (0–1, default 1)
 * @returns {{lcoe: number, doublings: number, capacityRatio: number, series: {year: number, lcoe: number}[]}}
 */
export function lcoeLearningProjection({ lcoe0, year0, year, learningRate, capacityGrowth, learningShare = 1 }) {
  requireNumber('lcoe0', lcoe0, { min: 0 }); requireNumber('year0', year0); requireNumber('year', year, { min: year0 });
  requireNumber('capacityGrowth', capacityGrowth, { gt: -1 }); requireNumber('learningShare', learningShare, { min: 0, max: 1 });
  const { b } = learningExponent({ learningRate });
  const at = (y) => { const ratio = Math.pow(1 + capacityGrowth, y - year0); return lcoe0 * (learningShare * Math.pow(ratio, b) + 1 - learningShare); };
  const series = [];
  for (let y = Math.ceil(year0); y <= year; y++) series.push({ year: y, lcoe: at(y) });
  const capacityRatio = Math.pow(1 + capacityGrowth, year - year0);
  return { lcoe: at(year), doublings: Math.log2(capacityRatio), capacityRatio, series };
}

/**
 * Capacity-factor and degradation adjustment of annual energy: E_t = capacity × 8760 × CF × (1 − degradation)^(t−1).
 * Ref: definition of capacity factor; IEA/NEA (2020) ch. 2.
 * @param {{capacityMW: number, capacityFactor: number, degradation?: number, years: number}} o MW, fraction, fraction
 *   per year, count
 * @returns {number[]} MWh per operating year (length `years`)
 */
export function energyFromCapacityFactor({ capacityMW, capacityFactor, degradation = 0, years }) {
  requireNumber('capacityMW', capacityMW, { min: 0 }); requireNumber('capacityFactor', capacityFactor, { min: 0, max: 1 });
  requireNumber('degradation', degradation, { min: 0, lt: 1 }); requireNumber('years', years, { min: 1, integer: true });
  return Array.from({ length: years }, (_, k) => capacityMW * 8760 * capacityFactor * Math.pow(1 - degradation, k));
}
