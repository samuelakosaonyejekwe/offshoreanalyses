// Sustainability and life-cycle assessment: embodied carbon from quantities × emission factors, fabrication
// energy, transport, marine operations, operation, decommissioning and recycling credit, life-cycle GWP by phase,
// carbon intensity, avoided emissions, carbon and energy payback, EROI, carbon cost and abatement economics,
// circularity indicators, footprint data, taxonomy threshold check, SDG mapping and a weighted score.
//
// IMPORTANT — every emission factor, energy factor and vessel fuel consumption in the default tables is an
// indicative mid-range value of the order reported in open literature. They are NOT product data: replace them
// with EPDs, supplier declarations or project measurements. Only the fuel carbon factors (t CO2 per t fuel) are
// exact regulatory values.
// Units: masses in tonnes (t), emissions in t CO2e unless stated.

import { requireNumber, requireSeries } from './tvm.js';
import { npv } from './cashflow.js';

const NOTE = 'indicative — replace with EPD/project data';
const SRC = 'order of magnitude of values widely reported in open literature (e.g. ICE database, worldsteel LCI); not a quotation of any single entry';
const ef = (value, range, note = '') => ({ value, unit: 'kg CO2e/kg', range, basisYear: 2019, source: SRC, status: NOTE, note });

/**
 * Editable table of default cradle-to-gate emission factors (EN 15804 modules A1–A3) by material:
 * { value, unit, range: [low, high], basisYear, source, status, note }. Every entry is indicative.
 * Ref: none authoritative — see `source` and `status` on each entry.
 */
export const defaultFactors = {
  steelPlate: ef(2.5, [1.5, 3.1], 'structural plate, mainly primary (BF-BOF) route'),
  steelSections: ef(1.6, [0.7, 2.9], 'sections and tubulars; strongly dependent on recycled content'),
  rebar: ef(1.4, [0.5, 2.3], 'reinforcing bar; EAF route at the low end'),
  concreteC20: ef(0.10, [0.07, 0.13], 'ready-mix C20/25, unreinforced'),
  concreteC30: ef(0.13, [0.09, 0.16], 'ready-mix C30/37, unreinforced'),
  concreteC40: ef(0.16, [0.11, 0.20], 'ready-mix C40/50, unreinforced'),
  concreteC50: ef(0.18, [0.13, 0.23], 'ready-mix C50/60, unreinforced'),
  cement: ef(0.83, [0.60, 0.95], 'Portland cement; lower with clinker substitution'),
  grout: ef(0.60, [0.35, 0.85], 'high-strength cementitious grout'),
  aggregate: ef(0.005, [0.003, 0.010], 'crushed rock / sand and gravel at quarry gate'),
  armourRock: ef(0.010, [0.003, 0.080], 'quarried armour stone; excludes transport'),
  copper: ef(3.0, [2.5, 4.5], 'conductor-grade copper, typical recycled share'),
  aluminium: ef(9.0, [6.0, 17.0], 'primary/recycled mix; smelter electricity dominates'),
  lead: ef(1.7, [1.3, 3.4], 'cable sheath lead'),
  zincAnode: ef(3.1, [2.6, 4.2], 'zinc sacrificial anodes'),
  aluminiumAnode: ef(12.0, [8.0, 17.0], 'aluminium alloy sacrificial anodes (mainly primary metal)'),
  xlpe: ef(2.5, [1.9, 3.5], 'cross-linked polyethylene insulation'),
  polyethylene: ef(2.0, [1.6, 2.6], 'HDPE/MDPE sheathing and pipe coating'),
  polyurethane: ef(4.5, [3.5, 5.5], 'polyurethane coating / elastomer'),
  gfrp: ef(8.0, [1.5, 9.0], 'glass-fibre reinforced polymer'),
  cfrp: ef(30.0, [20.0, 45.0], 'carbon-fibre reinforced polymer'),
  syntacticFoam: ef(6.0, [3.0, 10.0], 'syntactic buoyancy foam; low confidence'),
};

/**
 * CO2 emitted per tonne of marine fuel burned (t CO2 / t fuel), the non-dimensional conversion factors C_F used
 * for EEDI/EEXI/CII. Exact regulatory values.
 * Ref: IMO Resolution MEPC.308(73), 2018 Guidelines on the method of calculation of the attained EEDI, §2.2.1.
 */
export const FUEL_CO2_FACTORS = Object.freeze({ MGO: 3.206, LFO: 3.151, HFO: 3.114, LNG: 2.750, LPG_PROPANE: 3.000,
  LPG_BUTANE: 3.030, METHANOL: 1.375, ETHANOL: 1.913 });

/**
 * Editable placeholder fuel consumption of vessel classes while working (t fuel per day) and placeholder
 * transport emission factors (kg CO2e per tonne-kilometre). Indicative only.
 * Ref: none authoritative — replace with vessel-specific consumption and operator/carrier data.
 */
export const defaultOperationalFactors = {
  source: 'indicative placeholder — replace with vessel specifications, charter-party data or carrier factors',
  basisYear: 2023,
  vesselFuel_tPerDay: { heavyLift: 40, wtiv: 25, pipelay: 45, cableLay: 20, ocv: 15, dsv: 15, ahts: 15, tug: 8, psv: 8,
    rockDump: 20, sov: 6, ctv: 1.5, rig: 20, barge: 0 },
  transport_kgPerTkm: { ship: 0.012, barge: 0.030, rail: 0.030, road: 0.100 },
  fabricationEnergy_kWhPerT: { steel: 400, concrete: 15 },
};

/**
 * Embodied carbon inventory (cradle to gate): Σ mass × emission factor. Items may reference a key of the factor
 * table or carry their own factor (e.g. from an EPD).
 * Ref: EN 15804:2012+A2:2019 modules A1–A3; ISO 14040/14044 (LCA framework).
 * @param {{items: {name?: string, material?: string, mass_t: number, factor?: number}[], factors?: object}} o mass in
 *   tonnes; factor in kg CO2e/kg (= t CO2e/t); `factors` overrides or extends the default table (number or { value })
 * @returns {{total_t: number, lines: {name: string, material: string, mass_t: number, factor: number,
 *   emissions_t: number, low_t: number | null, high_t: number | null, indicative: boolean}[], low_t: number,
 *   high_t: number, usesIndicativeFactors: boolean}} t CO2e; low/high use the factor ranges where available
 */
export function embodiedCarbon({ items, factors }) {
  if (!Array.isArray(items) || items.length === 0) throw new RangeError('items must be a non-empty array');
  const lines = []; let total = 0; let low = 0; let high = 0; let ind = false;
  for (const it of items) {
    requireNumber(`mass_t of '${it.name ?? it.material}'`, it.mass_t, { min: 0 });
    let f; let range = null; let indicative = false;
    if (it.factor !== undefined) f = requireNumber('factor', it.factor, { min: 0 });
    else {
      const o = factors ? factors[it.material] : undefined;
      if (o !== undefined) { f = typeof o === 'number' ? o : o.value; range = typeof o === 'object' ? o.range ?? null : null; }
      else {
        const d = defaultFactors[it.material];
        if (!d) throw new RangeError(`no emission factor for material '${it.material}'`);
        f = d.value; range = d.range; indicative = true;
      }
      requireNumber(`factor of '${it.material}'`, f, { min: 0 });
    }
    const e = it.mass_t * f;
    total += e; low += it.mass_t * (range ? range[0] : f); high += it.mass_t * (range ? range[1] : f); ind = ind || indicative;
    lines.push({ name: it.name ?? it.material, material: it.material ?? 'user', mass_t: it.mass_t, factor: f, emissions_t: e,
      low_t: range ? it.mass_t * range[0] : null, high_t: range ? it.mass_t * range[1] : null, indicative });
  }
  return { total_t: total, lines, low_t: low, high_t: high, usesIndicativeFactors: ind };
}

/**
 * Fabrication emissions from process energy: mass × specific energy × grid (or fuel) carbon intensity.
 * Ref: EN 15804 module A3; energy × emission-factor method of the GHG Protocol (Scope 2 guidance).
 * @param {{mass_t: number, energy_kWhPerT: number, gridIntensity_gPerKWh: number}} o tonnes, kWh per tonne, g CO2e/kWh
 * @returns {{energy_MWh: number, emissions_t: number}}
 */
export function fabricationEmissions({ mass_t, energy_kWhPerT, gridIntensity_gPerKWh }) {
  requireNumber('mass_t', mass_t, { min: 0 }); requireNumber('energy_kWhPerT', energy_kWhPerT, { min: 0 });
  requireNumber('gridIntensity_gPerKWh', gridIntensity_gPerKWh, { min: 0 });
  const kWh = mass_t * energy_kWhPerT;
  return { energy_MWh: kWh / 1000, emissions_t: kWh * gridIntensity_gPerKWh / 1e6 };
}

/**
 * Transport emissions by tonne-kilometre: Σ mass × distance × mode factor.
 * Ref: EN 15804 module A4; ISO 14083:2023 (quantification of transport-chain GHG emissions) for the t·km method.
 * @param {{legs: {name?: string, mass_t: number, distance_km: number, mode: string, factor_kgPerTkm?: number}[],
 *   factors?: Record<string, number>}} o mode = 'ship' | 'barge' | 'rail' | 'road' (or any key supplied in factors)
 * @returns {{total_t: number, tkm: number, lines: object[], usesIndicativeFactors: boolean}} t CO2e
 */
export function transportEmissions({ legs, factors }) {
  if (!Array.isArray(legs)) throw new RangeError('legs must be an array');
  let total = 0; let tkm = 0; let ind = false; const lines = [];
  for (const l of legs) {
    requireNumber('mass_t', l.mass_t, { min: 0 }); requireNumber('distance_km', l.distance_km, { min: 0 });
    let f = l.factor_kgPerTkm ?? (factors ? factors[l.mode] : undefined);
    if (f === undefined) { f = defaultOperationalFactors.transport_kgPerTkm[l.mode]; ind = true; }
    if (f === undefined) throw new RangeError(`no transport factor for mode '${l.mode}'`);
    requireNumber('transport factor', f, { min: 0 });
    const e = l.mass_t * l.distance_km * f / 1000;
    total += e; tkm += l.mass_t * l.distance_km;
    lines.push({ name: l.name ?? l.mode, mode: l.mode, tkm: l.mass_t * l.distance_km, factor_kgPerTkm: f, emissions_t: e });
  }
  return { total_t: total, tkm, lines, usesIndicativeFactors: ind };
}

/**
 * Marine-operations emissions: Σ vessel days × fuel consumption (t/day) × fuel CO2 factor (t CO2 / t fuel).
 * The CO2 factors are the exact IMO values (MGO 3.206, HFO 3.114, LNG 2.750); consumption defaults are placeholders.
 * Methane slip and N2O are not included (CO2 only) unless added through `co2eUplift`.
 * Ref: IMO Resolution MEPC.308(73) for C_F; EN 15804 modules A5 / B / C for the life-cycle allocation.
 * @param {{vessels: {vessel: string, days: number, fuel?: string, fuel_tPerDay?: number, count?: number}[],
 *   co2eUplift?: number}} o vessel = class key of the placeholder table (or give fuel_tPerDay); fuel = key of
 *   FUEL_CO2_FACTORS (default 'MGO'); co2eUplift = fractional addition for non-CO2 gases (default 0)
 * @returns {{total_t: number, fuel_t: number, lines: object[], usesPlaceholderConsumption: boolean}} t CO2(e), t fuel
 */
export function marineOpsEmissions({ vessels, co2eUplift = 0 }) {
  if (!Array.isArray(vessels)) throw new RangeError('vessels must be an array');
  requireNumber('co2eUplift', co2eUplift, { min: 0 });
  let total = 0; let fuelT = 0; let ph = false; const lines = [];
  for (const v of vessels) {
    requireNumber(`days of '${v.vessel}'`, v.days, { min: 0 });
    const cf = FUEL_CO2_FACTORS[v.fuel ?? 'MGO'];
    if (!cf) throw new RangeError(`fuel '${v.fuel}' is not in FUEL_CO2_FACTORS`);
    let rate = v.fuel_tPerDay;
    if (rate === undefined) { rate = defaultOperationalFactors.vesselFuel_tPerDay[v.vessel]; ph = true; }
    if (rate === undefined) throw new RangeError(`no fuel consumption for vessel class '${v.vessel}'`);
    requireNumber('fuel_tPerDay', rate, { min: 0 });
    const fuel = v.days * rate * (v.count ?? 1); const e = fuel * cf * (1 + co2eUplift);
    total += e; fuelT += fuel;
    lines.push({ vessel: v.vessel, days: v.days, fuel: v.fuel ?? 'MGO', fuel_t: fuel, co2Factor: cf, emissions_t: e });
  }
  return { total_t: total, fuel_t: fuelT, lines, usesPlaceholderConsumption: ph };
}

/**
 * Operational emissions per year: logistics fuel, imported grid electricity, flaring and fugitive methane (for
 * oil and gas), each as an explicit argument.
 * Ref: GHG Protocol Corporate Standard (Scopes 1 and 2); IPCC AR6 WG1 Table 7.15 for the methane GWP100 default
 * (29.8 for fossil methane) — pass another value to use a different assessment report.
 * @param {{logisticsFuel_t?: number, fuel?: string, gridImport_MWh?: number, gridIntensity_gPerKWh?: number,
 *   flaredGas_t?: number, flareFactor_tCO2PerT?: number, fugitiveMethane_t?: number, methaneGwp?: number,
 *   other_t?: number}} o annual quantities; flareFactor = t CO2 per t gas flared (composition dependent, typical
 *   2.5–3.0, no default assumed beyond 2.75)
 * @returns {{total_t: number, breakdown: Record<string, number>}} t CO2e per year
 */
export function operationEmissions({ logisticsFuel_t = 0, fuel = 'MGO', gridImport_MWh = 0, gridIntensity_gPerKWh = 0,
  flaredGas_t = 0, flareFactor_tCO2PerT = 2.75, fugitiveMethane_t = 0, methaneGwp = 29.8, other_t = 0 }) {
  const cf = FUEL_CO2_FACTORS[fuel];
  if (!cf) throw new RangeError(`fuel '${fuel}' is not in FUEL_CO2_FACTORS`);
  for (const [k, v] of Object.entries({ logisticsFuel_t, gridImport_MWh, gridIntensity_gPerKWh, flaredGas_t, fugitiveMethane_t, other_t })) requireNumber(k, v, { min: 0 });
  const b = { logistics: logisticsFuel_t * cf, gridImport: gridImport_MWh * gridIntensity_gPerKWh / 1000,
    flaring: flaredGas_t * flareFactor_tCO2PerT, fugitive: fugitiveMethane_t * methaneGwp, other: other_t };
  return { total_t: Object.values(b).reduce((s, v) => s + v, 0), breakdown: b };
}

/**
 * End-of-life recycling credit (EN 15804 module D, net-scrap approach): only the scrap recovered beyond the
 * recycled content already used earns a credit, credit = mass × (recovery rate − recycled content) × (primary
 * factor − secondary factor). A negative result is a burden (more scrap used than returned).
 * Ref: EN 15804:2012+A2:2019 Annex D (net output flow); worldsteel LCA methodology report (2017) (scrap
 * burden / credit).
 * @param {{mass_t: number, recoveryRate: number, recycledContent: number, primaryFactor: number,
 *   secondaryFactor: number}} o tonnes, fractions 0–1, factors in kg CO2e/kg of primary and of secondary production
 * @returns {{credit_t: number, netScrap_t: number}} t CO2e (positive = emissions avoided) and net scrap (t)
 */
export function endOfLifeCredit({ mass_t, recoveryRate, recycledContent, primaryFactor, secondaryFactor }) {
  requireNumber('mass_t', mass_t, { min: 0 }); requireNumber('recoveryRate', recoveryRate, { min: 0, max: 1 });
  requireNumber('recycledContent', recycledContent, { min: 0, max: 1 }); requireNumber('primaryFactor', primaryFactor, { min: 0 });
  requireNumber('secondaryFactor', secondaryFactor, { min: 0 });
  const net = mass_t * (recoveryRate - recycledContent);
  return { credit_t: net * (primaryFactor - secondaryFactor), netScrap_t: net };
}

/**
 * Life-cycle global-warming potential by phase and in total. Phases follow the EN 15804 modules: product (A1–A3),
 * transport (A4), installation (A5), operation (B, annual × years or a total), end of life (C) and the recycling
 * credit beyond the system boundary (D), which is reported separately and included in `totalWithCredit` only.
 * Ref: EN 15804:2012+A2:2019 (module structure); ISO 14040:2006 / ISO 14044:2006.
 * @param {{product_t?: number, transport_t?: number, installation_t?: number, operationPerYear_t?: number,
 *   operatingYears?: number, operationTotal_t?: number, endOfLife_t?: number, credit_t?: number}} o t CO2e
 * @returns {{total_t: number, totalWithCredit_t: number, byPhase: Record<string, number>,
 *   shares: Record<string, number>}} shares of the total without credit
 */
export function lifeCycleGwp({ product_t = 0, transport_t = 0, installation_t = 0, operationPerYear_t = 0, operatingYears = 0,
  operationTotal_t, endOfLife_t = 0, credit_t = 0 }) {
  const op = operationTotal_t ?? operationPerYear_t * operatingYears;
  const byPhase = { 'A1-A3 product': product_t, 'A4 transport': transport_t, 'A5 installation': installation_t,
    'B operation': op, 'C end of life': endOfLife_t };
  for (const [k, v] of Object.entries(byPhase)) requireNumber(k, v, { min: 0 });
  requireNumber('credit_t', credit_t);
  const total = Object.values(byPhase).reduce((s, v) => s + v, 0);
  const shares = {}; for (const k of Object.keys(byPhase)) shares[k] = total > 0 ? byPhase[k] / total : 0;
  return { total_t: total, totalWithCredit_t: total - credit_t, byPhase: { ...byPhase, 'D credit': -credit_t }, shares };
}

/**
 * Carbon intensity of the product: g CO2e per kWh for energy, kg CO2e per boe for hydrocarbons.
 * Ref: functional-unit normalisation, ISO 14044 §4.2.3.2.
 * @param {{emissions_t: number, energy_MWh?: number, production_boe?: number}} o life-cycle emissions (t CO2e) and
 *   life-cycle output
 * @returns {{gPerKWh: number | null, kgPerBoe: number | null}}
 */
export function carbonIntensity({ emissions_t, energy_MWh, production_boe }) {
  requireNumber('emissions_t', emissions_t);
  return { gPerKWh: energy_MWh > 0 ? emissions_t / energy_MWh * 1000 : null,
    kgPerBoe: production_boe > 0 ? emissions_t * 1000 / production_boe : null };
}

/**
 * Avoided emissions against a reference intensity (grid average, marginal plant or displaced fuel), carbon payback
 * time, energy payback time and energy return on investment.
 * The reference intensity is an argument because it is location and time dependent (it may come from live data).
 * Ref: avoided-emissions accounting, GHG Protocol Guidelines for Grid-Connected Electricity Projects (2007);
 * EROI after Hall, Lambert & Balogh (2014), Energy Policy 64.
 * @param {{annualEnergy_MWh: number, referenceIntensity_gPerKWh: number, lifeCycleEmissions_t: number,
 *   operatingYears: number, embodiedEnergy_MWh?: number, annualEnergyInput_MWh?: number}} o embodiedEnergy = primary
 *   energy invested in manufacture, installation and decommissioning; annualEnergyInput = operating energy use
 * @returns {{annualAvoided_t: number, lifetimeAvoided_t: number, netAvoided_t: number, carbonPaybackYears: number,
 *   projectIntensity_gPerKWh: number, energyPaybackYears: number | null, eroi: number | null}}
 */
export function avoidedEmissions({ annualEnergy_MWh, referenceIntensity_gPerKWh, lifeCycleEmissions_t, operatingYears,
  embodiedEnergy_MWh, annualEnergyInput_MWh = 0 }) {
  requireNumber('annualEnergy_MWh', annualEnergy_MWh, { gt: 0 }); requireNumber('referenceIntensity_gPerKWh', referenceIntensity_gPerKWh, { min: 0 });
  requireNumber('lifeCycleEmissions_t', lifeCycleEmissions_t, { min: 0 }); requireNumber('operatingYears', operatingYears, { gt: 0 });
  const annual = annualEnergy_MWh * referenceIntensity_gPerKWh / 1000;
  const life = annual * operatingYears;
  const net = annualEnergy_MWh - annualEnergyInput_MWh;
  return { annualAvoided_t: annual, lifetimeAvoided_t: life, netAvoided_t: life - lifeCycleEmissions_t,
    carbonPaybackYears: annual > 0 ? lifeCycleEmissions_t / annual : Infinity,
    projectIntensity_gPerKWh: lifeCycleEmissions_t / (annualEnergy_MWh * operatingYears) * 1000,
    energyPaybackYears: embodiedEnergy_MWh !== undefined && net > 0 ? embodiedEnergy_MWh / net : null,
    eroi: embodiedEnergy_MWh !== undefined ? annualEnergy_MWh * operatingYears / (embodiedEnergy_MWh + annualEnergyInput_MWh * operatingYears) : null };
}

/**
 * Carbon cost along a carbon-price path and its present value (internal carbon pricing): cost_t = emissions_t ×
 * price_t; NPV including carbon = NPV − PV(carbon cost).
 * Ref: internal carbon pricing practice, e.g. CDP (2021), Putting a Price on Carbon; discounting as elsewhere.
 * @param {{emissions_t: number[], price: number[] | number, rate: number, npv?: number}} o emissions per period
 *   (index 0 = time zero), price in currency per t CO2e, discount rate; optional project NPV to adjust
 * @returns {{cost: number[], pv: number, total: number, npvWithCarbon: number | null}} currency
 */
export function carbonCost({ emissions_t, price, rate, npv: base }) {
  requireSeries('emissions_t', emissions_t);
  const cost = Array.from(emissions_t, (e, k) => e * requireNumber('price', typeof price === 'number' ? price : price[Math.min(k, price.length - 1)], { min: 0 }));
  const pv = npv(cost, rate);
  return { cost, pv, total: cost.reduce((s, v) => s + v, 0), npvWithCarbon: base === undefined ? null : base - pv };
}

/**
 * Marginal abatement cost of a design alternative: MAC = (NPV_base − NPV_alternative) / (emissions_base −
 * emissions_alternative), the cost per tonne of CO2e avoided (negative = the alternative saves money and carbon).
 * Ref: marginal abatement cost definition, e.g. Kesicki & Strachan (2011), Environmental Science & Policy 14(8).
 * @param {{npvBase: number, npvAlternative: number, emissionsBase_t: number, emissionsAlternative_t: number}} o
 * @returns {{mac: number, abatement_t: number, cost: number}} currency per t CO2e, t CO2e, currency (NPV given up)
 */
export function marginalAbatementCost({ npvBase, npvAlternative, emissionsBase_t, emissionsAlternative_t }) {
  requireNumber('npvBase', npvBase); requireNumber('npvAlternative', npvAlternative);
  requireNumber('emissionsBase_t', emissionsBase_t); requireNumber('emissionsAlternative_t', emissionsAlternative_t);
  const abatement = emissionsBase_t - emissionsAlternative_t; const cost = npvBase - npvAlternative;
  return { mac: abatement !== 0 ? cost / abatement : NaN, abatement_t: abatement, cost };
}

/**
 * Marginal abatement cost curve data: measures sorted by cost per tonne with cumulative abatement, and the set
 * that is worthwhile at a given carbon price.
 * Ref: MAC curve construction, Kesicki & Strachan (2011); measures are assumed independent and additive.
 * @param {{measures: {name: string, abatement_t: number, cost: number}[], carbonPrice?: number}} o abatement over
 *   the life (t CO2e, > 0) and NPV cost (currency; negative = saving) of each measure
 * @returns {{curve: {name: string, mac: number, abatement_t: number, cumulativeFrom_t: number, cumulativeTo_t: number,
 *   worthwhile: boolean | null}[], totalAbatement_t: number, worthwhileAbatement_t: number | null}}
 */
export function abatementCurve({ measures, carbonPrice }) {
  if (!Array.isArray(measures) || measures.length === 0) throw new RangeError('measures must be a non-empty array');
  const rows = measures.map((m) => ({ name: m.name, mac: requireNumber('cost', m.cost) / requireNumber('abatement_t', m.abatement_t, { gt: 0 }),
    abatement_t: m.abatement_t })).sort((a, b) => a.mac - b.mac);
  let cum = 0; let ok = carbonPrice === undefined ? null : 0;
  const curve = rows.map((r) => {
    const w = carbonPrice === undefined ? null : r.mac <= carbonPrice;
    const row = { ...r, cumulativeFrom_t: cum, cumulativeTo_t: cum + r.abatement_t, worthwhile: w };
    cum += r.abatement_t; if (w) ok += r.abatement_t;
    return row;
  });
  return { curve, totalAbatement_t: cum, worthwhileAbatement_t: ok };
}

/**
 * Circularity indicators: recyclable mass fraction of the bill of materials and the material circularity
 * indicator in its basic form: virgin feedstock V = M(1 − F_R − F_U); unrecoverable waste W = W_0 + (W_F + W_C)/2
 * with W_0 = M(1 − C_R − C_U), W_C = M(1 − E_C)C_R, W_F = M(1 − E_F)F_R/E_F; linear flow index
 * LFI = (V + W) / (2M + (W_F − W_C)/2); utility X = (L/L_av)(U/U_av); MCI = max(0, 1 − LFI × 0.9/X).
 * Ref: Ellen MacArthur Foundation & Granta Design (2015/2019), Circularity Indicators — Methodology.
 * @param {{items?: {name?: string, mass_t: number, recyclableFraction: number}[], mass_t?: number,
 *   recycledFeedstock?: number, reusedFeedstock?: number, collectedForRecycling?: number, collectedForReuse?: number,
 *   recyclingEfficiency?: number, feedstockRecyclingEfficiency?: number, lifeRatio?: number,
 *   intensityRatio?: number}} o fractions 0–1; lifeRatio = L/L_av, intensityRatio = U/U_av (default 1)
 * @returns {{recyclableFraction: number | null, mci: number, lfi: number, utility: number}}
 */
export function circularity({ items, mass_t, recycledFeedstock = 0, reusedFeedstock = 0, collectedForRecycling = 0,
  collectedForReuse = 0, recyclingEfficiency = 1, feedstockRecyclingEfficiency = 1, lifeRatio = 1, intensityRatio = 1 }) {
  let rec = null; let M = mass_t;
  if (items) {
    let m = 0; let r = 0;
    for (const it of items) { m += requireNumber('mass_t', it.mass_t, { min: 0 }); r += it.mass_t * requireNumber('recyclableFraction', it.recyclableFraction, { min: 0, max: 1 }); }
    rec = m > 0 ? r / m : null; M = M ?? m;
  }
  requireNumber('mass_t', M, { gt: 0 });
  for (const [k, v] of Object.entries({ recycledFeedstock, reusedFeedstock, collectedForRecycling, collectedForReuse })) requireNumber(k, v, { min: 0, max: 1 });
  if (recycledFeedstock + reusedFeedstock > 1 || collectedForRecycling + collectedForReuse > 1) throw new RangeError('feedstock and collection fractions must each sum to <= 1');
  requireNumber('recyclingEfficiency', recyclingEfficiency, { gt: 0, max: 1 }); requireNumber('feedstockRecyclingEfficiency', feedstockRecyclingEfficiency, { gt: 0, max: 1 });
  const V = M * (1 - recycledFeedstock - reusedFeedstock);
  const W0 = M * (1 - collectedForRecycling - collectedForReuse);
  const Wc = M * (1 - recyclingEfficiency) * collectedForRecycling;
  const Wf = M * (1 - feedstockRecyclingEfficiency) * recycledFeedstock / feedstockRecyclingEfficiency;
  const W = W0 + (Wf + Wc) / 2;
  const lfi = (V + W) / (2 * M + (Wf - Wc) / 2);
  const X = requireNumber('lifeRatio', lifeRatio, { gt: 0 }) * requireNumber('intensityRatio', intensityRatio, { gt: 0 });
  return { recyclableFraction: rec, mci: Math.max(0, 1 - lfi * 0.9 / X), lfi, utility: X };
}

/**
 * Seabed footprint and water use roll-ups: Σ count × area (permanent and temporary separately) and Σ quantity ×
 * specific water use.
 * Ref: environmental-statement quantity schedules; simple summations.
 * @param {{footprint?: {name: string, area_m2: number, count?: number, temporary?: boolean}[],
 *   water?: {name: string, quantity: number, m3PerUnit: number}[]}} o areas in m², water in m³ per unit of quantity
 * @returns {{permanent_m2: number, temporary_m2: number, total_m2: number, water_m3: number}}
 */
export function footprintAndWater({ footprint = [], water = [] }) {
  let perm = 0; let temp = 0; let w = 0;
  for (const f of footprint) { const a = requireNumber(`area of '${f.name}'`, f.area_m2, { min: 0 }) * (f.count ?? 1); if (f.temporary) temp += a; else perm += a; }
  for (const x of water) w += requireNumber(`quantity of '${x.name}'`, x.quantity, { min: 0 }) * requireNumber(`m3PerUnit of '${x.name}'`, x.m3PerUnit, { min: 0 });
  return { permanent_m2: perm, temporary_m2: temp, total_m2: perm + temp, water_m3: w };
}

/**
 * Qualitative environmental checklist (data only): activity triggers and the topic each raises. The application
 * flags an item when the project includes the trigger; no quantitative assessment is implied.
 * Ref: topic list typical of offshore environmental impact assessments (e.g. OSPAR and national EIA scoping
 * guidance); qualitative.
 */
export const ENVIRONMENTAL_CHECKLIST = Object.freeze([
  { key: 'impactPiling', topic: 'Underwater noise', note: 'Impact piling: marine-mammal and fish noise exposure; consider soft start, bubble curtains, seasonal windows.' },
  { key: 'uxoClearance', topic: 'Underwater noise', note: 'UXO clearance: high-order detonation noise; consider low-order deflagration.' },
  { key: 'seabedDisturbance', topic: 'Benthic habitat', note: 'Trenching, dredging, rock placement: habitat loss and sediment plumes.' },
  { key: 'scourProtection', topic: 'Benthic habitat', note: 'Rock or mattresses change the substrate (reef effect and habitat conversion).' },
  { key: 'emf', topic: 'Electromagnetic fields', note: 'Power cables: EMF for electro-sensitive species; burial depth mitigates.' },
  { key: 'dischargeProducedWater', topic: 'Discharges', note: 'Produced water and drill cuttings: oil-in-water and chemical discharge limits.' },
  { key: 'flaring', topic: 'Air emissions', note: 'Flaring and venting: CO2, methane and black carbon.' },
  { key: 'vesselTraffic', topic: 'Collision and disturbance', note: 'Construction and service vessels: collision risk, disturbance, invasive species via ballast.' },
  { key: 'rotatingBlades', topic: 'Birds and bats', note: 'Wind turbines: collision risk and displacement.' },
  { key: 'decommissioningInSitu', topic: 'Legacy', note: 'Leaving infrastructure in place: long-term degradation, snagging risk, monitoring liability.' },
]);

/**
 * Flag checklist items for a project description (object of booleans keyed as in ENVIRONMENTAL_CHECKLIST).
 * Ref: ENVIRONMENTAL_CHECKLIST (qualitative).
 * @param {Record<string, boolean>} project trigger flags
 * @returns {{key: string, topic: string, note: string, flagged: boolean}[]}
 */
export function environmentalFlags(project) {
  return ENVIRONMENTAL_CHECKLIST.map((c) => ({ ...c, flagged: !!(project && project[c.key]) }));
}

/**
 * Substantial-contribution threshold check on life-cycle carbon intensity (EU-taxonomy style): passes when the
 * intensity is below the threshold. The threshold is an argument (default 100 g CO2e/kWh) because it is set by
 * regulation and may change; passing the check is not a statement of taxonomy alignment, which also needs the
 * do-no-significant-harm and minimum-safeguard tests.
 * Ref: Commission Delegated Regulation (EU) 2021/2139, Annex I (life-cycle GHG threshold of 100 g CO2e/kWh for
 * several electricity-generation activities).
 * @param {{intensity_gPerKWh: number, threshold_gPerKWh?: number}} o
 * @returns {{pass: boolean, margin_gPerKWh: number, utilisation: number, threshold_gPerKWh: number}}
 */
export function taxonomyThresholdCheck({ intensity_gPerKWh, threshold_gPerKWh = 100 }) {
  requireNumber('intensity_gPerKWh', intensity_gPerKWh, { min: 0 }); requireNumber('threshold_gPerKWh', threshold_gPerKWh, { gt: 0 });
  return { pass: intensity_gPerKWh < threshold_gPerKWh, margin_gPerKWh: threshold_gPerKWh - intensity_gPerKWh,
    utilisation: intensity_gPerKWh / threshold_gPerKWh, threshold_gPerKWh };
}

/**
 * Qualitative mapping of project attributes to UN Sustainable Development Goals (data only; a prompt for
 * reporting, not an assessment).
 * Ref: United Nations (2015), Transforming our world: the 2030 Agenda for Sustainable Development (goal titles).
 */
export const SDG_MAPPING = Object.freeze([
  { sdg: 7, title: 'Affordable and clean energy', relevance: 'Energy generated, levelised cost, capacity added.' },
  { sdg: 8, title: 'Decent work and economic growth', relevance: 'Local content, jobs during construction and operation, safety performance.' },
  { sdg: 9, title: 'Industry, innovation and infrastructure', relevance: 'Port, grid and fabrication infrastructure; technology readiness.' },
  { sdg: 12, title: 'Responsible consumption and production', relevance: 'Material efficiency, recycled content, end-of-life recovery, circularity indicator.' },
  { sdg: 13, title: 'Climate action', relevance: 'Life-cycle GWP, carbon intensity, avoided emissions, carbon payback.' },
  { sdg: 14, title: 'Life below water', relevance: 'Seabed footprint, underwater noise, discharges, decommissioning strategy.' },
]);

/**
 * Multi-criteria sustainability score: each criterion is scaled linearly between its worst (0) and best (1)
 * values (clamped) and combined with editable weights.
 * Ref: weighted-sum aggregation (see decision.weightedSum); the scales and weights are the user's value
 * judgements and are reported back with the result.
 * @param {{criteria: {name: string, value: number, best: number, worst: number, weight: number}[]}} o best and
 *   worst in the unit of the value (best may be smaller than worst for "lower is better")
 * @returns {{score: number, items: {name: string, normalised: number, weight: number, contribution: number}[]}}
 *   score 0–1
 */
export function sustainabilityScore({ criteria }) {
  if (!Array.isArray(criteria) || criteria.length === 0) throw new RangeError('criteria must be a non-empty array');
  let sw = 0; for (const c of criteria) sw += requireNumber(`weight of '${c.name}'`, c.weight, { min: 0 });
  if (!(sw > 0)) throw new RangeError('weights must sum to a positive number');
  let score = 0;
  const items = criteria.map((c) => {
    requireNumber(`value of '${c.name}'`, c.value);
    if (c.best === c.worst) throw new RangeError(`best and worst of '${c.name}' must differ`);
    const nrm = Math.min(1, Math.max(0, (c.value - c.worst) / (c.best - c.worst)));
    const contribution = nrm * c.weight / sw; score += contribution;
    return { name: c.name, normalised: nrm, weight: c.weight / sw, contribution };
  });
  return { score, items };
}
