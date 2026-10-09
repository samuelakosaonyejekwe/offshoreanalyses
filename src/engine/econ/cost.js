// Capital-cost estimation: bottom-up, capacity scaling, factorial, transparent parametric offshore models,
// learning curves, AACE estimate classes with probabilistic contingency, escalation, roll-ups and benchmarks.
//
// IMPORTANT — every unit rate and day rate in `defaultRates` is an indicative placeholder so that the models run
// out of the box. They are NOT market data and must be replaced with current quotations or indices. All cost
// functions accept a `rates` argument ({ key: number | { value } }) that overrides the defaults.
// Quantities use the customary cost-engineering units stated per argument (t = tonne, km, MW, days).

import { requireNumber, requireSeries } from './tvm.js';
import { sampleInputs, quantile } from './risk.js';

const PLACEHOLDER = 'indicative placeholder — replace with current quotations/indices';
const rate = (value, unit, note = '') => ({ value, unit, currency: 'USD', basisYear: 2023, source: PLACEHOLDER, note });

/**
 * Editable table of default unit rates, day rates and percentage allowances. Each entry:
 * { value, unit, currency, basisYear, source, note }. Every value is an indicative placeholder (see `source`);
 * the application is expected to overwrite them from live indices, quotations or user data.
 * Ref: none — placeholders only, deliberately not attributed to any publication.
 */
export const defaultRates = {
  // fabricated steel and concrete, supplied at yard (currency per tonne or per m³)
  'steel.jacket': rate(4500, 'USD/t', 'fabricated jacket steel incl. coating and appurtenances'),
  'steel.piles': rate(1800, 'USD/t', 'rolled and welded tubular piles'),
  'steel.topsidesStructural': rate(9000, 'USD/t', 'topsides structural steel fabricated'),
  'steel.topsidesIntegrated': rate(35000, 'USD/t', 'topsides all-in (structure, equipment, bulks) per dry tonne'),
  'steel.monopile': rate(2200, 'USD/t', 'monopile, fabricated and coated'),
  'steel.transitionPiece': rate(3500, 'USD/t', 'transition piece incl. secondary steel'),
  'steel.suctionBucket': rate(4000, 'USD/t', 'suction bucket/caisson fabricated'),
  'steel.subseaStructure': rate(12000, 'USD/t', 'manifolds, PLETs, templates incl. piping'),
  'steel.floaterSemi': rate(4500, 'USD/t', 'semi-submersible hull steel'),
  'steel.floaterSpar': rate(3500, 'USD/t', 'spar hull steel'),
  'steel.floaterTlp': rate(5500, 'USD/t', 'TLP hull steel'),
  'concrete.gbs': rate(700, 'USD/m3', 'slip-formed structural concrete in place incl. formwork'),
  'concrete.floater': rate(900, 'USD/m3', 'floating substructure concrete in place'),
  'concrete.marine': rate(450, 'USD/m3', 'marine civil concrete in place'),
  'steel.rebar': rate(1500, 'USD/t', 'reinforcement cut, bent and fixed'),
  'ballast.solid': rate(60, 'USD/t', 'solid ballast supplied and placed'),
  // mooring
  'mooring.chain': rate(3000, 'USD/t', 'studless chain'),
  'mooring.wire': rate(6000, 'USD/t', 'spiral strand / six-strand wire'),
  'mooring.polyester': rate(9000, 'USD/t', 'polyester rope'),
  'anchor.drag': rate(250000, 'USD/each', 'drag embedment anchor'),
  'anchor.suction': rate(900000, 'USD/each', 'suction anchor'),
  'anchor.pile': rate(700000, 'USD/each', 'driven anchor pile'),
  // pipelines, umbilicals, flexibles
  'pipe.linepipe': rate(1500, 'USD/t', 'carbon-steel linepipe'),
  'pipe.coatingAnticorrosion': rate(60, 'USD/m', 'anti-corrosion coating and anodes'),
  'pipe.coatingConcrete': rate(150, 'USD/m', 'concrete weight coating'),
  'pipe.trenching': rate(120, 'USD/m', 'post-lay trenching'),
  'pipe.rockDump': rate(60, 'USD/t', 'rock supplied and placed'),
  'umbilical.static': rate(900, 'USD/m', 'static umbilical supply'),
  'umbilical.dynamic': rate(1500, 'USD/m', 'dynamic umbilical supply'),
  'flexible.small': rate(2500, 'USD/m', 'flexible pipe up to about 8 in'),
  'flexible.large': rate(5000, 'USD/m', 'flexible pipe above about 8 in'),
  // cables and electrical
  'cable.array.33kV': rate(350000, 'USD/km', 'array cable supply'),
  'cable.array.66kV': rate(450000, 'USD/km', 'array cable supply'),
  'cable.array.132kV': rate(600000, 'USD/km', 'array cable supply'),
  'cable.export.220kV_HVAC': rate(1200000, 'USD/km', 'export cable supply, three-core'),
  'cable.export.275kV_HVAC': rate(1400000, 'USD/km', 'export cable supply, three-core'),
  'cable.export.320kV_HVDC': rate(1000000, 'USD/km', 'export cable supply, bipole pair'),
  'cable.export.525kV_HVDC': rate(1600000, 'USD/km', 'export cable supply, bipole pair'),
  'cable.layBurial.array': rate(250000, 'USD/km', 'lay and burial'),
  'cable.layBurial.export': rate(500000, 'USD/km', 'lay and burial'),
  'cable.termination.array': rate(120000, 'USD/each', 'pull-in, termination and testing'),
  'cable.termination.export': rate(400000, 'USD/each', 'pull-in, jointing/termination and testing'),
  'substation.hvac': rate(180000, 'USD/MW', 'offshore HVAC substation incl. topside and foundation'),
  'substation.hvdc': rate(450000, 'USD/MW', 'offshore HVDC converter platform'),
  'turbine.supply': rate(1400000, 'USD/MW', 'wind turbine supply incl. tower'),
  // marine civil works
  'civil.deck': rate(3500, 'USD/m2', 'suspended deck on piles'),
  'civil.pileSteel': rate(2500, 'USD/t', 'tubular pile supply'),
  'civil.pileDriving': rate(800, 'USD/m', 'pile installation per metre driven'),
  'civil.rock': rate(45, 'USD/t', 'armour/core rock supplied and placed'),
  'civil.dredging': rate(12, 'USD/m3', 'dredging and disposal'),
  // vessel and spread day rates
  'vessel.heavyLift': rate(400000, 'USD/day', 'heavy-lift crane vessel spread'),
  'vessel.wtiv': rate(250000, 'USD/day', 'wind turbine installation jack-up'),
  'vessel.pipelay': rate(450000, 'USD/day', 'pipelay vessel spread'),
  'vessel.cableLay': rate(150000, 'USD/day', 'cable-lay vessel spread'),
  'vessel.ocv': rate(120000, 'USD/day', 'offshore construction / ROV support vessel'),
  'vessel.dsv': rate(150000, 'USD/day', 'diving support vessel'),
  'vessel.ahts': rate(60000, 'USD/day', 'anchor-handling tug supply'),
  'vessel.tug': rate(15000, 'USD/day', 'ocean-going tug'),
  'vessel.psv': rate(25000, 'USD/day', 'platform supply vessel'),
  'vessel.barge': rate(8000, 'USD/day', 'cargo/launch barge'),
  'vessel.rockDump': rate(90000, 'USD/day', 'fall-pipe rock installation vessel'),
  'vessel.sov': rate(45000, 'USD/day', 'service operation vessel'),
  'vessel.ctv': rate(4000, 'USD/day', 'crew transfer vessel'),
  'vessel.rig': rate(300000, 'USD/day', 'mobile drilling unit'),
  'vessel.mobDemob': rate(3000000, 'USD/campaign', 'mobilisation and demobilisation of a main spread'),
  // percentage allowances (fractions)
  'indirect.engineering': rate(0.08, 'fraction of direct cost', 'engineering and design'),
  'indirect.projectManagement': rate(0.05, 'fraction of direct cost', 'project management and supervision'),
  'indirect.insurance': rate(0.015, 'fraction of direct cost', 'construction all-risk insurance'),
  'indirect.ownersCost': rate(0.05, 'fraction of direct cost', "owner's team, permits, studies"),
  'allowance.contingency': rate(0.15, 'fraction of base estimate', 'contingency'),
  'allowance.commissioning': rate(0.03, 'fraction of installed cost', 'pre-commissioning and commissioning'),
  'allowance.decommissioning': rate(0.6, 'fraction of installation cost', 'removal cost relative to installation cost'),
};

/** Density of structural steel used for tonnage from dimensions, kg/m³ (EN 1991-1-1 Table A.4 nominal value). */
export const STEEL_DENSITY = 7850;

/**
 * Resolve a rate: the caller's override ({ key: number } or { key: { value } }) or the placeholder default.
 * Ref: none (lookup helper).
 * @param {string} key key of defaultRates
 * @param {Record<string, number | {value: number}>} [rates] overrides
 * @returns {{value: number, unit: string, source: string, basisYear: number | null, overridden: boolean}}
 */
export function resolveRate(key, rates) {
  const o = rates ? rates[key] : undefined;
  const d = defaultRates[key];
  if (o !== undefined) {
    const value = typeof o === 'number' ? o : o.value;
    requireNumber(`rate '${key}'`, value, { min: 0 });
    return { value, unit: (typeof o === 'object' && o.unit) || (d ? d.unit : ''), basisYear: (typeof o === 'object' && o.basisYear) || null,
      source: (typeof o === 'object' && o.source) || 'user supplied', overridden: true };
  }
  if (!d) throw new RangeError(`rate '${key}' is not defined; supply it in rates`);
  return { value: d.value, unit: d.unit, basisYear: d.basisYear, source: d.source, overridden: false };
}

const mkLine = (name, quantity, unit, key, rates, category) => {
  requireNumber(`${name} quantity`, quantity, { min: 0 });
  const r = resolveRate(key, rates);
  return { name, quantity, unit, rateKey: key, rate: r.value, rateSource: r.source, cost: quantity * r.value, category };
};
const pack = (lines, extra = {}) => {
  let total = 0; let placeholders = 0;
  for (const l of lines) { total += l.cost; if (l.rateSource === PLACEHOLDER) placeholders++; }
  return { total, lines, usesPlaceholderRates: placeholders > 0, ...extra };
};
const tubeTonnes = (D, t, L) => {
  requireNumber('D', D, { gt: 0 }); requireNumber('t', t, { gt: 0, lt: D / 2 }); requireNumber('L', L, { min: 0 });
  return Math.PI * (D - t) * t * L * STEEL_DENSITY / 1000;
};

// ---------------------------------------------------------------------------------------------------------------
// (1) Bottom-up
// ---------------------------------------------------------------------------------------------------------------

/**
 * Bottom-up estimate: Σ quantity × unit rate, plus indirects as fractions of direct cost, contingency on the base
 * estimate (direct + indirect) and escalation of the whole to the spend date.
 * Ref: AACE International RP 34R-05 (basis of estimate); Humphreys, Project and Cost Engineers' Handbook, ch. 1–2.
 * @param {{items: {name: string, quantity: number, unitRate: number, unit?: string, discipline?: string,
 *   category?: string, code?: string}[], indirects?: Record<string, number>, contingency?: number,
 *   escalation?: {rate: number, years: number}, rates?: object}} o unitRate in currency per unit; indirects are
 *   fractions of direct cost (defaults from defaultRates 'indirect.*'); contingency fraction of base estimate
 *   (default from 'allowance.contingency'); escalation rate per year (fraction) over `years`
 * @returns {{direct: number, indirects: Record<string, number>, indirectTotal: number, base: number,
 *   contingency: number, escalation: number, total: number, byDiscipline: Record<string, number>,
 *   byCategory: Record<string, number>, lines: object[]}} currency
 */
export function bottomUpEstimate({ items, indirects, contingency, escalation, rates }) {
  if (!Array.isArray(items) || items.length === 0) throw new RangeError('items must be a non-empty array');
  const byDiscipline = {}; const byCategory = {}; const lines = [];
  let direct = 0;
  for (const it of items) {
    requireNumber(`quantity of '${it.name}'`, it.quantity, { min: 0 });
    requireNumber(`unitRate of '${it.name}'`, it.unitRate, { min: 0 });
    const cost = it.quantity * it.unitRate;
    direct += cost;
    const d = it.discipline ?? 'unassigned'; const c = it.category ?? 'unassigned';
    byDiscipline[d] = (byDiscipline[d] ?? 0) + cost; byCategory[c] = (byCategory[c] ?? 0) + cost;
    lines.push({ ...it, cost });
  }
  const frac = indirects ?? { engineering: resolveRate('indirect.engineering', rates).value,
    projectManagement: resolveRate('indirect.projectManagement', rates).value,
    insurance: resolveRate('indirect.insurance', rates).value, ownersCost: resolveRate('indirect.ownersCost', rates).value };
  const ind = {}; let indirectTotal = 0;
  for (const k of Object.keys(frac)) { requireNumber(`indirects.${k}`, frac[k], { min: 0 }); ind[k] = direct * frac[k]; indirectTotal += ind[k]; }
  const base = direct + indirectTotal;
  const cf = contingency ?? resolveRate('allowance.contingency', rates).value;
  requireNumber('contingency', cf, { min: 0 });
  const cont = base * cf;
  const esc = escalation ? (base + cont) * (Math.pow(1 + requireNumber('escalation.rate', escalation.rate, { gt: -1 }),
    requireNumber('escalation.years', escalation.years)) - 1) : 0;
  return { direct, indirects: ind, indirectTotal, base, contingency: cont, escalation: esc, total: base + cont + esc,
    byDiscipline, byCategory, lines };
}

// ---------------------------------------------------------------------------------------------------------------
// (2) Analogous / capacity scaling, (3) factorial
// ---------------------------------------------------------------------------------------------------------------

/**
 * Cost–capacity scaling C2 = C1 (S2/S1)^n × (index2/index1) × location factor × complexity factor.
 * n = 0.6 is the classic "six-tenths rule"; typical exponents lie between about 0.4 and 0.9 by equipment type.
 * Ref: Williams, R. (1947), Six-tenths factor aids in approximating costs, Chem. Eng. 54(12); Chilton (1950).
 * @param {{C1: number, S1: number, S2: number, n?: number, index1?: number, index2?: number,
 *   locationFactor?: number, complexityFactor?: number}} o C1 known cost (currency), S1/S2 capacities (same unit),
 *   n exponent (-), cost-index values at the dates of C1 and of the estimate, dimensionless factors
 * @returns {number} C2 (currency)
 */
export function capacityScaling({ C1, S1, S2, n = 0.6, index1 = 1, index2 = 1, locationFactor = 1, complexityFactor = 1 }) {
  requireNumber('C1', C1, { min: 0 }); requireNumber('S1', S1, { gt: 0 }); requireNumber('S2', S2, { gt: 0 });
  requireNumber('n', n, { min: 0, max: 1.5 }); requireNumber('index1', index1, { gt: 0 }); requireNumber('index2', index2, { gt: 0 });
  requireNumber('locationFactor', locationFactor, { gt: 0 }); requireNumber('complexityFactor', complexityFactor, { gt: 0 });
  return C1 * Math.pow(S2 / S1, n) * (index2 / index1) * locationFactor * complexityFactor;
}

/**
 * Scaling exponent implied by two known (capacity, cost) points: n = ln(C2/C1) / ln(S2/S1).
 * Ref: Williams (1947) cost–capacity relation, inverted.
 * @param {{C1: number, S1: number, C2: number, S2: number}} o costs (currency) and capacities (same unit)
 * @returns {number} exponent (-)
 */
export function scalingExponent({ C1, S1, C2, S2 }) {
  requireNumber('C1', C1, { gt: 0 }); requireNumber('C2', C2, { gt: 0 }); requireNumber('S1', S1, { gt: 0 }); requireNumber('S2', S2, { gt: 0 });
  if (S1 === S2) throw new RangeError('S1 and S2 must differ');
  return Math.log(C2 / C1) / Math.log(S2 / S1);
}

/**
 * Value of a cost index at a date by linear interpolation in a supplied index series, and the escalation factor
 * between two dates.
 * Ref: AACE International RP 58R-10 (escalation estimating using indices).
 * @param {{years: number[], values: number[]}} series ascending years (decimal) and index values (> 0)
 * @param {number} fromYear date of the known cost
 * @param {number} toYear date of the estimate
 * @returns {{from: number, to: number, factor: number}} index values and to/from ratio
 */
export function indexFactor(series, fromYear, toYear) {
  if (!series) throw new RangeError('series is required');
  requireSeries('series.years', series.years, { minLength: 2 }); requireSeries('series.values', series.values, { minLength: 2 });
  const at = (y) => {
    const xs = series.years; const ys = series.values;
    if (y < xs[0] || y > xs[xs.length - 1]) throw new RangeError(`year ${y} is outside the index series`);
    let k = 1; while (k < xs.length - 1 && xs[k] < y) k++;
    return ys[k - 1] + (ys[k] - ys[k - 1]) * (y - xs[k - 1]) / (xs[k] - xs[k - 1]);
  };
  const from = at(fromYear); const to = at(toYear);
  return { from, to, factor: to / from };
}

/**
 * Lang factors (ratio of total plant capital cost to delivered equipment cost) and Hand factors (installation
 * multipliers by equipment type). Published values; editable copies may be passed to factorialEstimate.
 * Ref: Lang, H.J. (1947/1948), Chem. Eng. 54(10), 55(6): 3.10 solid, 3.63 solid–fluid, 4.74 fluid processing;
 * Hand, W.E. (1958), Petroleum Refiner 37(9).
 */
export const LANG_FACTORS = Object.freeze({ solid: 3.10, solidFluid: 3.63, fluid: 4.74 });
export const HAND_FACTORS = Object.freeze({ fractionatingColumn: 4, pressureVessel: 4, heatExchanger: 3.5,
  firedHeater: 2, pump: 4, compressor: 2.5, instrument: 4, miscellaneous: 2.5 });

/**
 * Factorial estimate. 'lang': total = Lang factor × Σ delivered equipment cost. 'hand': total = Σ equipment cost ×
 * Hand factor of its type. Factors are for onshore process plant; offshore installation needs project-specific
 * factors, which is why `factors` is an argument.
 * Ref: Lang (1947/1948); Hand (1958); Towler & Sinnott, Chemical Engineering Design, ch. 7.
 * @param {{method?: 'lang' | 'hand', equipment: {name: string, cost: number, type?: string}[], plantType?: string,
 *   factors?: Record<string, number>}} o equipment delivered costs (currency), plantType key for Lang, type keys
 *   for Hand, optional factor overrides
 * @returns {{total: number, equipment: number, factor: number, lines: object[]}} currency; factor = total/equipment
 */
export function factorialEstimate({ method = 'lang', equipment, plantType = 'fluid', factors }) {
  if (!Array.isArray(equipment) || equipment.length === 0) throw new RangeError('equipment must be a non-empty array');
  let eq = 0; let total = 0; const lines = [];
  const tab = { ...(method === 'lang' ? LANG_FACTORS : HAND_FACTORS), ...(factors || {}) };
  for (const e of equipment) {
    requireNumber(`cost of '${e.name}'`, e.cost, { min: 0 });
    const key = method === 'lang' ? plantType : (e.type ?? 'miscellaneous');
    const f = tab[key];
    if (!(f > 0)) throw new RangeError(`no factor for '${key}'`);
    eq += e.cost; total += e.cost * f;
    lines.push({ name: e.name, cost: e.cost, factor: f, installed: e.cost * f });
  }
  return { total, equipment: eq, factor: eq > 0 ? total / eq : NaN, lines };
}

// ---------------------------------------------------------------------------------------------------------------
// (4) Parametric offshore models — every equation is quantity × rate and is listed line by line in the result
// ---------------------------------------------------------------------------------------------------------------

/**
 * Jacket substructure supply cost: jacket tonnage × fabricated rate + pile tonnage × pile rate.
 * Ref: transparent tonnage × rate model (placeholder rates in defaultRates).
 * @param {{jacket_t: number, piles_t?: number, rates?: object}} o tonnes
 * @returns {{total: number, lines: object[], usesPlaceholderRates: boolean}} currency
 */
export function jacketCost({ jacket_t, piles_t = 0, rates }) {
  return pack([mkLine('Jacket steel', jacket_t, 't', 'steel.jacket', rates, 'substructure'),
    mkLine('Foundation piles', piles_t, 't', 'steel.piles', rates, 'substructure')]);
}

/**
 * Topsides supply cost: either structural tonnage × structural rate, or all-in dry weight × integrated rate.
 * Ref: transparent tonnage × rate model (placeholder rates).
 * @param {{structural_t?: number, integrated_t?: number, rates?: object}} o tonnes; give one or both
 * @returns {{total: number, lines: object[], usesPlaceholderRates: boolean}}
 */
export function topsidesCost({ structural_t = 0, integrated_t = 0, rates }) {
  return pack([mkLine('Topsides structural steel', structural_t, 't', 'steel.topsidesStructural', rates, 'topsides'),
    mkLine('Topsides integrated (all-in)', integrated_t, 't', 'steel.topsidesIntegrated', rates, 'topsides')]);
}

/**
 * Monopile and transition-piece supply cost. Tonnage from geometry: m = π (D − t) t L ρ_steel for each can
 * (thin-walled tube, 7850 kg/m³); a transition-piece mass may be given directly instead of its dimensions.
 * Ref: tube mass from geometry; placeholder rates.
 * @param {{D: number, t: number, L: number, tp?: {D?: number, t?: number, L?: number, mass_t?: number},
 *   rates?: object}} o monopile outer diameter, wall thickness and length (m); transition piece likewise
 * @returns {{total: number, lines: object[], monopile_t: number, tp_t: number, usesPlaceholderRates: boolean}}
 */
export function monopileCost({ D, t, L, tp, rates }) {
  const mp = tubeTonnes(D, t, L);
  const tpT = tp ? (tp.mass_t ?? tubeTonnes(tp.D, tp.t, tp.L)) : 0;
  return pack([mkLine('Monopile', mp, 't', 'steel.monopile', rates, 'substructure'),
    mkLine('Transition piece', tpT, 't', 'steel.transitionPiece', rates, 'substructure')], { monopile_t: mp, tp_t: tpT });
}

/**
 * Suction bucket (caisson) supply cost for n buckets: skirt π (D − t) t L ρ plus lid plate π D²/4 × tLid × ρ ×
 * stiffening factor.
 * Ref: steel mass from geometry; placeholder rate.
 * @param {{D: number, L: number, tSkirt: number, tLid: number, n?: number, lidStiffeningFactor?: number,
 *   rates?: object}} o diameter, skirt length, thicknesses (m), number of buckets, lid stiffener mass multiplier (>= 1)
 * @returns {{total: number, lines: object[], steel_t: number, usesPlaceholderRates: boolean}}
 */
export function suctionBucketCost({ D, L, tSkirt, tLid, n = 1, lidStiffeningFactor = 1.5, rates }) {
  requireNumber('n', n, { min: 1, integer: true }); requireNumber('tLid', tLid, { gt: 0 });
  requireNumber('lidStiffeningFactor', lidStiffeningFactor, { min: 1 });
  const one = tubeTonnes(D, tSkirt, L) + Math.PI * D * D / 4 * tLid * STEEL_DENSITY / 1000 * lidStiffeningFactor;
  return pack([mkLine('Suction buckets', one * n, 't', 'steel.suctionBucket', rates, 'substructure')], { steel_t: one * n });
}

/**
 * Concrete gravity-base structure cost: concrete volume × in-place rate + reinforcement (given tonnage, or volume ×
 * reinforcement ratio in kg/m³) × rate + solid ballast × rate.
 * Ref: transparent quantity × rate model (placeholder rates).
 * @param {{concrete_m3: number, rebar_t?: number, rebarRatio_kgm3?: number, ballast_t?: number, rates?: object}} o
 * @returns {{total: number, lines: object[], rebar_t: number, usesPlaceholderRates: boolean}}
 */
export function gbsCost({ concrete_m3, rebar_t, rebarRatio_kgm3 = 200, ballast_t = 0, rates }) {
  requireNumber('concrete_m3', concrete_m3, { min: 0 });
  const rb = rebar_t ?? concrete_m3 * requireNumber('rebarRatio_kgm3', rebarRatio_kgm3, { min: 0 }) / 1000;
  return pack([mkLine('Structural concrete', concrete_m3, 'm3', 'concrete.gbs', rates, 'substructure'),
    mkLine('Reinforcement', rb, 't', 'steel.rebar', rates, 'substructure'),
    mkLine('Solid ballast', ballast_t, 't', 'ballast.solid', rates, 'substructure')], { rebar_t: rb });
}

/**
 * Floating substructure cost: steel tonnage × rate for the hull type (semi, spar, TLP) and/or concrete volume ×
 * rate, plus solid ballast.
 * Ref: transparent quantity × rate model (placeholder rates).
 * @param {{type: 'semi' | 'spar' | 'tlp', steel_t?: number, concrete_m3?: number, ballast_t?: number, rates?: object}} o
 * @returns {{total: number, lines: object[], usesPlaceholderRates: boolean}}
 */
export function floaterCost({ type, steel_t = 0, concrete_m3 = 0, ballast_t = 0, rates }) {
  const key = { semi: 'steel.floaterSemi', spar: 'steel.floaterSpar', tlp: 'steel.floaterTlp' }[type];
  if (!key) throw new RangeError("type must be 'semi', 'spar' or 'tlp'");
  return pack([mkLine(`Hull steel (${type})`, steel_t, 't', key, rates, 'substructure'),
    mkLine('Hull concrete', concrete_m3, 'm3', 'concrete.floater', rates, 'substructure'),
    mkLine('Solid ballast', ballast_t, 't', 'ballast.solid', rates, 'substructure')]);
}

/**
 * Mooring system supply cost: Σ segment length × mass per metre × rate for its material + anchors × unit cost.
 * Ref: transparent quantity × rate model (placeholder rates).
 * @param {{segments: {type: 'chain' | 'wire' | 'polyester', length_m: number, massPerM_kg: number, count?: number}[],
 *   anchors?: {type: 'drag' | 'suction' | 'pile', count: number}[], rates?: object}} o count = number of identical lines
 * @returns {{total: number, lines: object[], mass_t: number, usesPlaceholderRates: boolean}}
 */
export function mooringCost({ segments, anchors = [], rates }) {
  if (!Array.isArray(segments)) throw new RangeError('segments must be an array');
  const lines = []; let mass = 0;
  for (const s of segments) {
    requireNumber('segment length_m', s.length_m, { min: 0 }); requireNumber('segment massPerM_kg', s.massPerM_kg, { min: 0 });
    const t = s.length_m * s.massPerM_kg / 1000 * (s.count ?? 1);
    mass += t;
    lines.push(mkLine(`Mooring ${s.type}`, t, 't', `mooring.${s.type}`, rates, 'mooring'));
  }
  for (const a of anchors) lines.push(mkLine(`Anchor (${a.type})`, a.count, 'each', `anchor.${a.type}`, rates, 'mooring'));
  return pack(lines, { mass_t: mass });
}

/**
 * Marine installation campaign cost: Σ vessel day rate × count × gross days + mobilisation/demobilisation, with
 * gross days = net working days / weather operability + fixed transit/standby days.
 * Ref: vessel day-rate model; operability from workability analysis (e.g. DNV-ST-N001 weather-restricted
 * operations) supplied by the caller.
 * @param {{spread: {vessel: string, count?: number, dayRate?: number}[], netDays: number, operability?: number,
 *   extraDays?: number, mobDemob?: number, rates?: object}} o vessel = key suffix of 'vessel.*' in defaultRates (or
 *   give dayRate, currency/day), netDays working days, operability fraction (0–1], extra days, mob/demob lump sum
 *   (default 'vessel.mobDemob')
 * @returns {{total: number, lines: object[], grossDays: number, weatherDays: number, spreadDayRate: number,
 *   usesPlaceholderRates: boolean}}
 */
export function installationCampaignCost({ spread, netDays, operability = 1, extraDays = 0, mobDemob, rates }) {
  if (!Array.isArray(spread) || spread.length === 0) throw new RangeError('spread must be a non-empty array');
  requireNumber('netDays', netDays, { min: 0 }); requireNumber('operability', operability, { gt: 0, max: 1 });
  requireNumber('extraDays', extraDays, { min: 0 });
  const grossDays = netDays / operability + extraDays;
  const lines = []; let dayRate = 0;
  for (const v of spread) {
    const key = `vessel.${v.vessel}`;
    const ov = v.dayRate !== undefined ? { ...(rates || {}), [key]: v.dayRate } : rates;
    const l = mkLine(`${v.vessel} × ${v.count ?? 1}`, grossDays * (v.count ?? 1), 'days', key, ov, 'installation');
    dayRate += l.rate * (v.count ?? 1);
    lines.push(l);
  }
  const mobKey = 'vessel.mobDemob';
  lines.push(mkLine('Mobilisation / demobilisation', 1, 'campaign', mobKey,
    mobDemob !== undefined ? { ...(rates || {}), [mobKey]: mobDemob } : rates, 'installation'));
  return pack(lines, { grossDays, weatherDays: grossDays - netDays - extraDays, spreadDayRate: dayRate });
}

/**
 * Rigid pipeline cost: linepipe tonnage (π (D − t) t ρ L) × rate + coatings per metre + installation spread
 * (lay days = length / lay rate / operability) + trenching per metre + rock tonnage × rate.
 * Ref: transparent quantity × rate model (placeholder rates).
 * @param {{length_m: number, D: number, t: number, concreteCoated?: boolean, layRate_mPerDay: number,
 *   operability?: number, trenchLength_m?: number, rock_t?: number, mobDemob?: number, rates?: object}} o SI lengths
 *   (m), lay rate in m/day, operability fraction
 * @returns {{total: number, lines: object[], linepipe_t: number, layDays: number, supply: number,
 *   installation: number, usesPlaceholderRates: boolean}}
 */
export function pipelineCost({ length_m, D, t, concreteCoated = false, layRate_mPerDay, operability = 1,
  trenchLength_m = 0, rock_t = 0, mobDemob, rates }) {
  requireNumber('length_m', length_m, { min: 0 }); requireNumber('layRate_mPerDay', layRate_mPerDay, { gt: 0 });
  const tonnes = tubeTonnes(D, t, length_m);
  const supply = [mkLine('Linepipe', tonnes, 't', 'pipe.linepipe', rates, 'supply'),
    mkLine('Anti-corrosion coating and anodes', length_m, 'm', 'pipe.coatingAnticorrosion', rates, 'supply')];
  if (concreteCoated) supply.push(mkLine('Concrete weight coating', length_m, 'm', 'pipe.coatingConcrete', rates, 'supply'));
  const inst = installationCampaignCost({ spread: [{ vessel: 'pipelay' }], netDays: length_m / layRate_mPerDay,
    operability, mobDemob, rates });
  const post = [mkLine('Trenching', trenchLength_m, 'm', 'pipe.trenching', rates, 'installation'),
    mkLine('Rock placement', rock_t, 't', 'pipe.rockDump', rates, 'installation')];
  const lines = [...supply, ...inst.lines, ...post];
  const s = supply.reduce((a, l) => a + l.cost, 0);
  const all = pack(lines);
  return { ...all, linepipe_t: tonnes, layDays: inst.grossDays, supply: s, installation: all.total - s };
}

/**
 * Umbilical or flexible pipe supply cost per metre by class.
 * Ref: transparent length × rate model (placeholder rates).
 * @param {{kind: 'umbilical' | 'flexible', cls: string, length_m: number, rates?: object}} o cls: 'static' |
 *   'dynamic' for umbilicals, 'small' | 'large' for flexibles
 * @returns {{total: number, lines: object[], usesPlaceholderRates: boolean}}
 */
export function umbilicalFlexibleCost({ kind, cls, length_m, rates }) {
  if (kind !== 'umbilical' && kind !== 'flexible') throw new RangeError("kind must be 'umbilical' or 'flexible'");
  return pack([mkLine(`${kind} (${cls})`, length_m, 'm', `${kind}.${cls}`, rates, 'supply')]);
}

/**
 * Subsea structure supply cost by tonnage.
 * Ref: transparent tonnage × rate model (placeholder rate).
 * @param {{mass_t: number, rates?: object}} o tonnes
 * @returns {{total: number, lines: object[], usesPlaceholderRates: boolean}}
 */
export function subseaStructureCost({ mass_t, rates }) {
  return pack([mkLine('Subsea structures', mass_t, 't', 'steel.subseaStructure', rates, 'supply')]);
}

/**
 * Array or export cable cost: supply per km by voltage/cross-section class + lay and burial per km + terminations.
 * Ref: transparent length × rate model (placeholder rates).
 * @param {{kind: 'array' | 'export', cls: string, length_km: number, terminations?: number, rates?: object}} o cls
 *   e.g. '66kV' (array) or '220kV_HVAC' (export); any class may be added through `rates`
 * @returns {{total: number, lines: object[], perKm: number, usesPlaceholderRates: boolean}}
 */
export function cableCost({ kind, cls, length_km, terminations = 0, rates }) {
  if (kind !== 'array' && kind !== 'export') throw new RangeError("kind must be 'array' or 'export'");
  const r = pack([mkLine(`${kind} cable supply (${cls})`, length_km, 'km', `cable.${kind}.${cls}`, rates, 'supply'),
    mkLine('Lay and burial', length_km, 'km', `cable.layBurial.${kind}`, rates, 'installation'),
    mkLine('Terminations', terminations, 'each', `cable.termination.${kind}`, rates, 'installation')]);
  return { ...r, perKm: length_km > 0 ? r.total / length_km : NaN };
}

/**
 * Offshore substation cost per MW of export capacity (HVAC or HVDC), and turbine supply per MW.
 * Ref: transparent capacity × rate model (placeholder rates).
 * @param {{substationMW?: number, substationType?: 'hvac' | 'hvdc', turbineMW?: number, rates?: object}} o MW
 * @returns {{total: number, lines: object[], usesPlaceholderRates: boolean}}
 */
export function generationPlantCost({ substationMW = 0, substationType = 'hvac', turbineMW = 0, rates }) {
  return pack([mkLine(`Offshore substation (${substationType})`, substationMW, 'MW', `substation.${substationType}`, rates, 'electrical'),
    mkLine('Wind turbines', turbineMW, 'MW', 'turbine.supply', rates, 'turbines')]);
}

/**
 * Marine civil structure cost (wharf, jetty, dolphin, breakwater): deck area, pile tonnage and driven length, rock
 * tonnage, concrete volume and dredged volume, each × its rate.
 * Ref: transparent quantity × rate model (placeholder rates).
 * @param {{deckArea_m2?: number, pile_t?: number, pileLength_m?: number, rock_t?: number, concrete_m3?: number,
 *   dredge_m3?: number, rates?: object}} o
 * @returns {{total: number, lines: object[], usesPlaceholderRates: boolean}}
 */
export function marineCivilCost({ deckArea_m2 = 0, pile_t = 0, pileLength_m = 0, rock_t = 0, concrete_m3 = 0, dredge_m3 = 0, rates }) {
  return pack([mkLine('Deck', deckArea_m2, 'm2', 'civil.deck', rates, 'civil'),
    mkLine('Pile supply', pile_t, 't', 'civil.pileSteel', rates, 'civil'),
    mkLine('Pile installation', pileLength_m, 'm', 'civil.pileDriving', rates, 'civil'),
    mkLine('Rock', rock_t, 't', 'civil.rock', rates, 'civil'),
    mkLine('Concrete', concrete_m3, 'm3', 'concrete.marine', rates, 'civil'),
    mkLine('Dredging', dredge_m3, 'm3', 'civil.dredging', rates, 'civil')]);
}

/**
 * Commissioning allowance as a fraction of installed cost, and decommissioning as a fraction of installation
 * cost or as a bottom-up list of items.
 * Ref: percentage-allowance models (placeholder fractions in defaultRates).
 * @param {{installedCost?: number, installationCost?: number, decomItems?: {name: string, cost: number}[],
 *   rates?: object}} o currency
 * @returns {{commissioning: number, decommissioning: number, decomBasis: 'fraction' | 'bottom-up'}}
 */
export function commissioningDecomCost({ installedCost = 0, installationCost = 0, decomItems, rates }) {
  requireNumber('installedCost', installedCost, { min: 0 }); requireNumber('installationCost', installationCost, { min: 0 });
  const commissioning = installedCost * resolveRate('allowance.commissioning', rates).value;
  if (decomItems) {
    let d = 0; for (const it of decomItems) d += requireNumber(`cost of '${it.name}'`, it.cost, { min: 0 });
    return { commissioning, decommissioning: d, decomBasis: 'bottom-up' };
  }
  return { commissioning, decommissioning: installationCost * resolveRate('allowance.decommissioning', rates).value,
    decomBasis: 'fraction' };
}

// ---------------------------------------------------------------------------------------------------------------
// (5) Learning / experience curves
// ---------------------------------------------------------------------------------------------------------------

/**
 * Learning exponent b = ln(progress ratio) / ln 2 (negative for learning); progress ratio = 1 − learning rate.
 * Ref: Wright, T.P. (1936), Factors affecting the cost of airplanes, J. Aeronautical Sciences 3(4).
 * @param {{progressRatio?: number, learningRate?: number}} o give one (fractions; PR in (0, 1.5])
 * @returns {{b: number, progressRatio: number, learningRate: number}}
 */
export function learningExponent({ progressRatio, learningRate }) {
  const pr = progressRatio ?? 1 - requireNumber('learningRate', learningRate, { lt: 1 });
  requireNumber('progressRatio', pr, { gt: 0, max: 1.5 });
  return { b: Math.log(pr) / Math.LN2, progressRatio: pr, learningRate: 1 - pr };
}

/**
 * Wright cumulative-average learning curve: average cost of the first n units Ȳ(n) = T1 n^b; total = T1 n^(b+1);
 * cost of unit n = total(n) − total(n − 1). Doubling n multiplies the cumulative average by the progress ratio.
 * Ref: Wright (1936).
 * @param {{T1: number, n: number, progressRatio?: number, learningRate?: number}} o first-unit cost (currency),
 *   number of units (>= 1)
 * @returns {{average: number, total: number, unit: number}} currency
 */
export function wrightCurve({ T1, n, progressRatio, learningRate }) {
  requireNumber('T1', T1, { min: 0 }); requireNumber('n', n, { min: 1 });
  const { b } = learningExponent({ progressRatio, learningRate });
  const total = T1 * Math.pow(n, b + 1);
  return { average: T1 * Math.pow(n, b), total, unit: total - (n > 1 ? T1 * Math.pow(n - 1, b + 1) : 0) };
}

/**
 * Crawford unit learning curve: cost of unit n is T1 n^b; total is the exact sum over units 1..n. Doubling n
 * multiplies the unit cost by the progress ratio.
 * Ref: Crawford, J.R. (1944), Learning curve, ship curve, ratios, related data, Lockheed Aircraft Corporation.
 * @param {{T1: number, n: number, progressRatio?: number, learningRate?: number}} o first-unit cost (currency),
 *   integer number of units (>= 1)
 * @returns {{unit: number, total: number, average: number}} currency
 */
export function crawfordCurve({ T1, n, progressRatio, learningRate }) {
  requireNumber('T1', T1, { min: 0 }); requireNumber('n', n, { min: 1, integer: true });
  const { b } = learningExponent({ progressRatio, learningRate });
  let total = 0;
  for (let k = 1; k <= n; k++) total += Math.pow(k, b);
  return { unit: T1 * Math.pow(n, b), total: T1 * total, average: T1 * total / n };
}

/**
 * Experience curve on cumulative volume: C(Q) = C0 (Q/Q0)^b — cost falls by the learning rate for each doubling
 * of cumulative production or installed capacity.
 * Ref: Boston Consulting Group (1968), Perspectives on Experience; Wright (1936).
 * @param {{c0: number, q0: number, q: number, learningRate: number}} o reference cost (currency per unit) at
 *   cumulative volume q0, target cumulative volume q (same unit), learning rate (fraction per doubling)
 * @returns {number} cost at cumulative volume q (currency per unit)
 */
export function costAtCumulativeVolume({ c0, q0, q, learningRate }) {
  requireNumber('c0', c0, { min: 0 }); requireNumber('q0', q0, { gt: 0 }); requireNumber('q', q, { gt: 0 });
  return c0 * Math.pow(q / q0, learningExponent({ learningRate }).b);
}

/**
 * Serial-production saving for n identical units (e.g. foundations) relative to n × first-unit cost.
 * Ref: Wright (1936) / Crawford (1944) as selected.
 * @param {{firstUnitCost: number, n: number, progressRatio?: number, learningRate?: number,
 *   model?: 'wright' | 'crawford'}} o currency, integer count
 * @returns {{total: number, withoutLearning: number, saving: number, savingFraction: number, average: number}}
 */
export function serialProductionSavings({ firstUnitCost, n, progressRatio, learningRate, model = 'crawford' }) {
  const c = model === 'wright' ? wrightCurve({ T1: firstUnitCost, n, progressRatio, learningRate })
    : crawfordCurve({ T1: firstUnitCost, n, progressRatio, learningRate });
  const flat = firstUnitCost * n;
  return { total: c.total, withoutLearning: flat, saving: flat - c.total, savingFraction: flat > 0 ? 1 - c.total / flat : 0,
    average: c.average };
}

// ---------------------------------------------------------------------------------------------------------------
// (6) Estimate classes and probabilistic contingency
// ---------------------------------------------------------------------------------------------------------------

/**
 * AACE International cost-estimate classification for the process industries: maturity of project definition
 * and typical expected accuracy range (low and high side, as fractions of the estimate, each given as the
 * narrow–wide band at about 80 % confidence after contingency).
 * Ref: AACE International Recommended Practice 18R-97, Table 1.
 */
export const AACE_CLASSES = Object.freeze({
  5: Object.freeze({ definition: [0, 0.02], usage: 'Concept screening', low: [-0.20, -0.50], high: [0.30, 1.00] }),
  4: Object.freeze({ definition: [0.01, 0.15], usage: 'Study or feasibility', low: [-0.15, -0.30], high: [0.20, 0.50] }),
  3: Object.freeze({ definition: [0.10, 0.40], usage: 'Budget authorisation or control', low: [-0.10, -0.20], high: [0.10, 0.30] }),
  2: Object.freeze({ definition: [0.30, 0.75], usage: 'Control or bid/tender', low: [-0.05, -0.15], high: [0.05, 0.20] }),
  1: Object.freeze({ definition: [0.65, 1.00], usage: 'Check estimate or bid/tender', low: [-0.03, -0.10], high: [0.03, 0.15] }),
});

/**
 * Accuracy range of an estimate for an AACE class: low/high bounds in currency for the narrow and the wide end of
 * the typical band.
 * Ref: AACE International RP 18R-97, Table 1.
 * @param {number} estimate point estimate (currency)
 * @param {1 | 2 | 3 | 4 | 5} cls estimate class
 * @returns {{cls: number, narrow: {low: number, high: number}, wide: {low: number, high: number}, usage: string}}
 */
export function estimateRange(estimate, cls) {
  requireNumber('estimate', estimate, { min: 0 });
  const c = AACE_CLASSES[cls];
  if (!c) throw new RangeError('cls must be 1, 2, 3, 4 or 5');
  return { cls, usage: c.usage, narrow: { low: estimate * (1 + c.low[0]), high: estimate * (1 + c.high[0]) },
    wide: { low: estimate * (1 + c.low[1]), high: estimate * (1 + c.high[1]) } };
}

/**
 * Probabilistic contingency: Monte Carlo (Latin hypercube) of line-item costs with optional rank correlations.
 * Contingency at a confidence level = percentile of total cost − base estimate. Percentiles are cumulative
 * (non-exceedance): P80 is the cost not exceeded with 80 % probability.
 * Ref: AACE International RP 41R-08 (range estimating) and 57R-09 (integrated cost and schedule risk analysis).
 * @param {{items: {name: string, dist: object}[], base?: number, correlations?: {a: string, b: string, rho: number}[],
 *   n?: number, seed?: number}} o dist = distribution spec of math/stats dist() in currency; base = deterministic
 *   base estimate (default: sum of distribution means); correlations = Spearman rank correlations between items
 * @returns {{base: number, mean: number, std: number, p10: number, p50: number, p80: number, p90: number,
 *   contingencyP50: number, contingencyP80: number, contingencyP90: number, samples: Float64Array,
 *   convention: string}} currency
 */
export function probabilisticContingency({ items, base, correlations, n = 5000, seed = 1 }) {
  if (!Array.isArray(items) || items.length === 0) throw new RangeError('items must be a non-empty array');
  const vars = {};
  for (const it of items) {
    if (vars[it.name]) throw new RangeError(`duplicate item name '${it.name}'`);
    vars[it.name] = it.dist;
  }
  const s = sampleInputs({ vars, n, seed, rankCorr: correlations });
  const total = new Float64Array(n);
  for (const name of Object.keys(vars)) { const col = s.samples[name]; for (let k = 0; k < n; k++) total[k] += col[k]; }
  let mean = 0; for (let k = 0; k < n; k++) mean += total[k]; mean /= n;
  let v = 0; for (let k = 0; k < n; k++) v += (total[k] - mean) ** 2;
  const sorted = Float64Array.from(total).sort();
  const b = base ?? s.means.reduce((a, m) => a + m, 0);
  const p = (q) => quantile(sorted, q);
  return { base: b, mean, std: Math.sqrt(v / (n - 1)), p10: p(0.1), p50: p(0.5), p80: p(0.8), p90: p(0.9),
    contingencyP50: p(0.5) - b, contingencyP80: p(0.8) - b, contingencyP90: p(0.9) - b, samples: total,
    convention: 'Cumulative (non-exceedance): Pxx is the cost not exceeded with xx % probability.' };
}

// ---------------------------------------------------------------------------------------------------------------
// (7) Escalation, currency exposure, roll-ups, benchmarks
// ---------------------------------------------------------------------------------------------------------------

/**
 * Composite escalation factor Σ w_i (I_i,to / I_i,from) from weighted cost indices (weights are normalised).
 * Ref: AACE International RP 58R-10 (escalation estimating principles and methods using indices).
 * @param {{weights: Record<string, number>, indices: Record<string, {from: number, to: number}>}} o weights by
 *   cost component (any positive scale), index values at the base date and at the target date
 * @returns {{factor: number, contributions: Record<string, number>}} factor to multiply base-date cost by
 */
export function compositeEscalation({ weights, indices }) {
  let sw = 0; for (const k of Object.keys(weights)) sw += requireNumber(`weights.${k}`, weights[k], { min: 0 });
  if (!(sw > 0)) throw new RangeError('weights must sum to a positive number');
  let factor = 0; const contributions = {};
  for (const k of Object.keys(weights)) {
    const i = indices[k];
    if (!i || !(i.from > 0) || !(i.to > 0)) throw new RangeError(`indices.${k} needs positive from and to`);
    contributions[k] = weights[k] / sw * (i.to / i.from);
    factor += contributions[k];
  }
  return { factor, contributions };
}

/**
 * Currency exposure split of an estimate: total by currency and converted to the project currency.
 * Ref: AACE International RP 58R-10; exposure reporting practice.
 * @param {{items: {name: string, amount: number, currency: string}[], rates: Record<string, number>,
 *   projectCurrency: string}} o rates[code] = units of `code` per unit of a common reference currency
 * @returns {{byCurrency: Record<string, {native: number, converted: number, share: number}>, total: number}}
 */
export function currencyExposure({ items, rates, projectCurrency }) {
  if (!(rates?.[projectCurrency] > 0)) throw new RangeError('rates must include the project currency');
  const by = {}; let total = 0;
  for (const it of items) {
    if (!(rates[it.currency] > 0)) throw new RangeError(`rates has no entry for '${it.currency}'`);
    const conv = requireNumber(`amount of '${it.name}'`, it.amount) * rates[projectCurrency] / rates[it.currency];
    by[it.currency] = by[it.currency] ?? { native: 0, converted: 0, share: 0 };
    by[it.currency].native += it.amount; by[it.currency].converted += conv; total += conv;
  }
  for (const k of Object.keys(by)) by[k].share = total !== 0 ? by[k].converted / total : 0;
  return { byCurrency: by, total };
}

/**
 * Cost-breakdown-structure roll-up: sums leaf costs up a dotted code hierarchy ('1.2.3' rolls into '1.2' and '1')
 * and summarises by discipline. Duplicate codes and missing costs raise an error so aggregation cannot silently
 * double count.
 * Ref: AACE International RP 20R-98 (project code of accounts).
 * @param {{code: string, name?: string, cost: number, discipline?: string}[]} items leaf items (currency)
 * @returns {{total: number, nodes: Record<string, number>, byDiscipline: Record<string, number>,
 *   levels: {code: string, level: number, cost: number}[]}}
 */
export function cbsRollup(items) {
  if (!Array.isArray(items) || items.length === 0) throw new RangeError('items must be a non-empty array');
  const nodes = {}; const byDiscipline = {}; const seen = new Set(); let total = 0;
  for (const it of items) {
    if (typeof it.code !== 'string' || !it.code) throw new RangeError('every item needs a code');
    if (seen.has(it.code)) throw new RangeError(`duplicate cost code '${it.code}'`);
    seen.add(it.code);
    requireNumber(`cost of '${it.code}'`, it.cost);
    const parts = it.code.split('.');
    for (let k = 1; k <= parts.length; k++) { const c = parts.slice(0, k).join('.'); nodes[c] = (nodes[c] ?? 0) + it.cost; }
    const d = it.discipline ?? 'unassigned';
    byDiscipline[d] = (byDiscipline[d] ?? 0) + it.cost; total += it.cost;
  }
  const levels = Object.keys(nodes).sort().map((code) => ({ code, level: code.split('.').length, cost: nodes[code] }));
  return { total, nodes, byDiscipline, levels };
}

/**
 * Benchmark metrics of an estimate: cost per MW, per tonne and per km (null where the divisor is not given).
 * Ref: unit-cost benchmarking practice (cost per unit of capacity, mass or length).
 * @param {{cost: number, MW?: number, tonnes?: number, km?: number}} o currency and normalising quantities
 * @returns {{perMW: number | null, perTonne: number | null, perKm: number | null}}
 */
export function benchmarkMetrics({ cost, MW, tonnes, km }) {
  requireNumber('cost', cost);
  const d = (q) => (q > 0 ? cost / q : null);
  return { perMW: d(MW), perTonne: d(tonnes), perKm: d(km) };
}
