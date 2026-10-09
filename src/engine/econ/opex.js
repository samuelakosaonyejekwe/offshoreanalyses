// Operating-expenditure models: fixed/variable build-up, offshore-wind O&M (port-based vs SOV), oil and gas
// facilities, pipeline/subsea IMR, escalation, availability and downtime cost, cost of failure, insurance recovery
// and life-extension versus decommissioning cash flows.
// Generic failure data and day rates are indicative placeholders to be replaced with operator or OEM records.

import { requireNumber, requireSeries } from './tvm.js';
import { npv } from './cashflow.js';

const PLACEHOLDER = 'indicative placeholder — replace with operator/OEM reliability records and current rates';

/**
 * Editable generic failure data for an offshore wind turbine, by subsystem: failures per turbine-year, mean
 * repair cost (materials, USD), mean active repair time (hours), technicians required and whether a heavy-lift
 * (jack-up) vessel is needed. Values only have the right order of magnitude; they are placeholders, not a dataset.
 * Ref: none — placeholders; replace with site/OEM data (the structure follows the usual subsystem breakdown of
 * published offshore-wind reliability surveys).
 */
export const defaultWindFailureData = {
  basisYear: 2023, currency: 'USD', source: PLACEHOLDER,
  subsystems: [
    { name: 'Pitch and hydraulics', rate: 1.0, repairCost: 3000, repairHours: 9, technicians: 2, heavyLift: false },
    { name: 'Generator', rate: 0.6, repairCost: 8000, repairHours: 12, technicians: 3, heavyLift: false },
    { name: 'Generator (major replacement)', rate: 0.03, repairCost: 300000, repairHours: 80, technicians: 6, heavyLift: true },
    { name: 'Gearbox / drivetrain', rate: 0.4, repairCost: 6000, repairHours: 12, technicians: 3, heavyLift: false },
    { name: 'Gearbox (major replacement)', rate: 0.04, repairCost: 500000, repairHours: 100, technicians: 6, heavyLift: true },
    { name: 'Blades', rate: 0.5, repairCost: 5000, repairHours: 15, technicians: 3, heavyLift: false },
    { name: 'Blade (major replacement)', rate: 0.01, repairCost: 350000, repairHours: 90, technicians: 6, heavyLift: true },
    { name: 'Power converter and electrical', rate: 0.9, repairCost: 4000, repairHours: 8, technicians: 2, heavyLift: false },
    { name: 'Controls and sensors', rate: 0.9, repairCost: 1500, repairHours: 6, technicians: 2, heavyLift: false },
    { name: 'Yaw system', rate: 0.2, repairCost: 3000, repairHours: 8, technicians: 2, heavyLift: false },
    { name: 'Other (tower, safety, auxiliaries)', rate: 1.2, repairCost: 1000, repairHours: 5, technicians: 2, heavyLift: false },
  ],
};

/**
 * Editable placeholder rates for the offshore-wind O&M model (USD, basis 2023): vessel day rates, technician
 * cost, vessel speeds and heavy-lift campaign parameters.
 * Ref: none — placeholders only.
 */
export const defaultOmRates = {
  basisYear: 2023, currency: 'USD', source: PLACEHOLDER,
  ctvDayRate: 4000, sovDayRate: 45000, heavyLiftDayRate: 200000, heavyLiftMobilisation: 800000,
  heavyLiftDaysPerRepair: 4, technicianCostPerYear: 110000, ctvSpeed_kn: 22, sovSpeed_kn: 12,
  onshoreBaseCostPerYear: 1500000, scheduledServiceCostPerTurbine: 40000, scheduledServiceHours: 50,
};

/**
 * Annual OPEX build-up from first principles: personnel, logistics (vessels, helicopters), planned and unplanned
 * maintenance, insurance as a fraction of capex, leases, grid/transmission charges, consumables/chemicals/fuel,
 * inspection campaigns (annualised by interval) and a per-unit variable cost.
 * Ref: cost build-up practice, AACE International RP 34R-05 (basis of estimate) applied to operating cost.
 * @param {{personnel?: {role?: string, count: number, costPerYear: number}[], logistics?: {name?: string,
 *   days: number, dayRate: number}[], maintenancePlanned?: number, maintenanceUnplanned?: number,
 *   insuranceRate?: number, capex?: number, leases?: number, gridCharges?: number, consumables?: number,
 *   inspections?: {name?: string, cost: number, intervalYears: number}[], variablePerUnit?: number,
 *   production?: number}} o currency per year unless stated; insuranceRate is a fraction of capex per year;
 *   production in units per year
 * @returns {{fixed: number, variable: number, total: number, breakdown: Record<string, number>}} currency per year
 */
export function opexBuildUp(o) {
  const nn = (name, v) => requireNumber(name, v ?? 0, { min: 0 });
  const b = {};
  b.personnel = (o.personnel || []).reduce((s, p) => s + nn('personnel.count', p.count) * nn('personnel.costPerYear', p.costPerYear), 0);
  b.logistics = (o.logistics || []).reduce((s, l) => s + nn('logistics.days', l.days) * nn('logistics.dayRate', l.dayRate), 0);
  b.maintenancePlanned = nn('maintenancePlanned', o.maintenancePlanned);
  b.maintenanceUnplanned = nn('maintenanceUnplanned', o.maintenanceUnplanned);
  b.insurance = nn('insuranceRate', o.insuranceRate) * nn('capex', o.capex);
  b.leases = nn('leases', o.leases);
  b.gridCharges = nn('gridCharges', o.gridCharges);
  b.consumables = nn('consumables', o.consumables);
  b.inspection = (o.inspections || []).reduce((s, i) => s + nn('inspections.cost', i.cost) / requireNumber('inspections.intervalYears', i.intervalYears, { gt: 0 }), 0);
  const fixed = Object.values(b).reduce((s, v) => s + v, 0);
  const variable = nn('variablePerUnit', o.variablePerUnit) * nn('production', o.production);
  b.variable = variable;
  return { fixed, variable, total: fixed + variable, breakdown: b };
}

/**
 * Parametric offshore-wind O&M model comparing a port-based crew-transfer-vessel (CTV) strategy with a service
 * operation vessel (SOV) stationed in the field.
 * Logic (every step visible in the result): transit time = distance / speed; productive hours per technician-shift
 * = shift length − 2 × transit (CTV) or shift length − in-field transfer allowance (SOV); corrective workload =
 * Σ failure rate × (repair hours × technicians) per turbine; technician-shifts needed = workload / productive
 * hours / weather accessibility; vessel days follow from technicians per vessel; downtime per failure = repair
 * hours + mean waiting time for access and logistics; availability = 1 − downtime hours / 8760.
 * Ref: structure follows published strategy models, e.g. Dalgic, Lazakis, Dinwoodie, McMillan & Revie (2015),
 * Advanced logistics planning for offshore wind farm operation and maintenance activities, Ocean Engineering 101;
 * all numerical data are placeholders supplied through `failureData` and `rates`.
 * @param {{turbines: number, ratedMW: number, distance_km: number, strategy: 'ctv' | 'sov', capacityFactor?: number,
 *   accessibility?: number, shiftHours?: number, techniciansPerVessel?: number, meanWaitHours?: number,
 *   heavyLiftWaitHours?: number, failureData?: object, rates?: object}} o accessibility = fraction of days
 *   workable for the chosen access system (0–1]; meanWaitHours = logistic + weather waiting per minor failure
 * @returns {{annualCost: number, costPerMW: number, breakdown: Record<string, number>, availability: number,
 *   downtimeHoursPerTurbine: number, transitHours: number, productiveHoursPerShift: number, vesselDays: number,
 *   technicians: number, vessels: number, lostEnergy_MWh: number, feasible: boolean, usesPlaceholderData: boolean}}
 */
export function windOmModel({ turbines, ratedMW, distance_km, strategy, capacityFactor = 0.45, accessibility,
  shiftHours = 12, techniciansPerVessel, meanWaitHours, heavyLiftWaitHours = 720, failureData = defaultWindFailureData,
  rates }) {
  requireNumber('turbines', turbines, { min: 1, integer: true }); requireNumber('ratedMW', ratedMW, { gt: 0 });
  requireNumber('distance_km', distance_km, { min: 0 }); requireNumber('capacityFactor', capacityFactor, { gt: 0, max: 1 });
  if (strategy !== 'ctv' && strategy !== 'sov') throw new RangeError("strategy must be 'ctv' or 'sov'");
  const R = { ...defaultOmRates, ...(rates || {}) };
  const sov = strategy === 'sov';
  const acc = requireNumber('accessibility', accessibility ?? (sov ? 0.85 : 0.65), { gt: 0, max: 1 });
  const perVessel = requireNumber('techniciansPerVessel', techniciansPerVessel ?? (sov ? 40 : 12), { min: 1 });
  const wait = requireNumber('meanWaitHours', meanWaitHours ?? (sov ? 24 : 48), { min: 0 });
  const speed = (sov ? R.sovSpeed_kn : R.ctvSpeed_kn) * 1.852;       // knots to km/h
  const transitHours = distance_km / speed;
  // CTV crews commute every shift; SOV crews live in the field and lose only an in-field transfer allowance.
  const productive = sov ? shiftHours - 1 : shiftHours - 2 * transitHours;
  const feasible = productive > 1;
  let techHours = 0; let materials = 0; let downtime = 0; let heavyEvents = 0;
  for (const s of failureData.subsystems) {
    techHours += s.rate * s.repairHours * s.technicians;
    materials += s.rate * s.repairCost;
    downtime += s.rate * (s.repairHours + (s.heavyLift ? heavyLiftWaitHours : wait));
    if (s.heavyLift) heavyEvents += s.rate;
  }
  techHours += R.scheduledServiceHours * 2;                             // annual service by a two-person team
  downtime += R.scheduledServiceHours;
  const totalTechHours = techHours * turbines;
  const shifts = feasible ? totalTechHours / productive : Infinity;      // technician-shifts per year
  const workableDays = 365 * acc;
  const technicians = Math.ceil(shifts / workableDays / 0.6);           // 0.6 = rota utilisation (leave, training)
  const vesselDays = shifts / perVessel;
  const vessels = Math.max(1, Math.ceil(vesselDays / workableDays));
  const breakdown = {
    vessels: sov ? vessels * 365 * R.sovDayRate : vessels * 365 * R.ctvDayRate,
    technicians: technicians * R.technicianCostPerYear,
    spares: materials * turbines,
    scheduledService: R.scheduledServiceCostPerTurbine * turbines,
    heavyLift: heavyEvents * turbines * (R.heavyLiftDaysPerRepair * R.heavyLiftDayRate) + (heavyEvents * turbines > 0 ? R.heavyLiftMobilisation * Math.max(1, Math.round(heavyEvents * turbines / 4)) : 0),
    onshoreBase: R.onshoreBaseCostPerYear,
  };
  const annualCost = Object.values(breakdown).reduce((s, v) => s + v, 0);
  const availability = Math.max(0, 1 - downtime / 8760);
  return { annualCost, costPerMW: annualCost / (turbines * ratedMW), breakdown, availability,
    downtimeHoursPerTurbine: downtime, transitHours, productiveHoursPerShift: productive, vesselDays, technicians,
    vessels, lostEnergy_MWh: turbines * ratedMW * capacityFactor * downtime, feasible,
    usesPlaceholderData: failureData === defaultWindFailureData || !rates };
}

/**
 * Distance at which the SOV strategy becomes cheaper than CTV, including the value of lost energy, by scanning
 * the distance range (the comparison is not guaranteed monotonic, so the first crossing is reported).
 * Ref: strategy break-even on the windOmModel equations.
 * @param {{base: object, energyPrice: number, distances_km: number[]}} o base = windOmModel arguments without
 *   strategy/distance; energyPrice in currency per MWh
 * @returns {{distances_km: number[], ctv: number[], sov: number[], breakEven_km: number | null}} total annual cost
 *   (O&M + lost revenue), currency per year
 */
export function windOmStrategyBreakEven({ base, energyPrice, distances_km }) {
  requireNumber('energyPrice', energyPrice, { min: 0 }); requireSeries('distances_km', distances_km, { minLength: 2 });
  const tot = (strategy, d) => {
    const r = windOmModel({ ...base, strategy, distance_km: d });
    return r.feasible ? r.annualCost + r.lostEnergy_MWh * energyPrice : Infinity;
  };
  const ctv = Array.from(distances_km, (d) => tot('ctv', d)); const sov = Array.from(distances_km, (d) => tot('sov', d));
  let be = null;
  for (let k = 0; k < ctv.length && be === null; k++) {
    if (sov[k] <= ctv[k]) {
      if (k === 0 || !Number.isFinite(ctv[k]) || !Number.isFinite(ctv[k - 1])) be = distances_km[k];
      else {
        const d0 = ctv[k - 1] - sov[k - 1]; const d1 = ctv[k] - sov[k];
        be = distances_km[k - 1] + (distances_km[k] - distances_km[k - 1]) * d0 / (d0 - d1);
      }
    }
  }
  return { distances_km: Array.from(distances_km), ctv, sov, breakEven_km: be };
}

/**
 * Oil and gas facilities OPEX: fixed cost as a fraction of capex per year plus a variable cost per barrel of oil
 * equivalent, plus optional tariffs.
 * Ref: screening-level opex model widely used in field-development economics, e.g. Jahn, Cook & Graham,
 * Hydrocarbon Exploration and Production, 2nd ed.; the fraction and unit cost are inputs.
 * @param {{capex: number, fixedFraction: number, production_boe: number | number[], variablePerBoe: number,
 *   tariffPerBoe?: number}} o capex (currency), fixed fraction of capex per year (typically a few per cent),
 *   production per year (boe), unit costs (currency/boe)
 * @returns {{fixed: number, variable: number[], total: number[], unitCost: number[]}} currency per year and
 *   currency per boe
 */
export function oilGasOpex({ capex, fixedFraction, production_boe, variablePerBoe, tariffPerBoe = 0 }) {
  requireNumber('capex', capex, { min: 0 }); requireNumber('fixedFraction', fixedFraction, { min: 0, max: 1 });
  requireNumber('variablePerBoe', variablePerBoe, { min: 0 });
  const prod = typeof production_boe === 'number' ? [production_boe] : Array.from(requireSeries('production_boe', production_boe));
  const fixed = capex * fixedFraction;
  const variable = prod.map((q) => q * (variablePerBoe + tariffPerBoe));
  const total = variable.map((v) => v + fixed);
  return { fixed, variable, total, unitCost: total.map((t, k) => (prod[k] > 0 ? t / prod[k] : NaN)) };
}

/**
 * Pipeline / subsea inspection, maintenance and repair (IMR) cost: campaigns every `intervalYears`, each needing
 * vessel days = length / survey rate + days per structure × structures, at a vessel day rate plus mobilisation;
 * returned as a year-by-year series and its annual average.
 * Ref: transparent interval × vessel-days model; inspection intervals from the integrity-management plan
 * (e.g. DNV-RP-F116).
 * @param {{length_km?: number, surveyRate_kmPerDay?: number, structures?: number, daysPerStructure?: number,
 *   intervalYears: number, dayRate: number, mobilisation?: number, operability?: number, years: number,
 *   firstYear?: number}} o operability fraction (0–1], years of operation, year of first campaign (1-based)
 * @returns {{series: number[], campaignCost: number, vesselDays: number, annualAverage: number}} currency
 */
export function imrCost({ length_km = 0, surveyRate_kmPerDay = 10, structures = 0, daysPerStructure = 1, intervalYears,
  dayRate, mobilisation = 0, operability = 1, years, firstYear }) {
  requireNumber('intervalYears', intervalYears, { min: 1, integer: true }); requireNumber('dayRate', dayRate, { min: 0 });
  requireNumber('years', years, { min: 1, integer: true }); requireNumber('operability', operability, { gt: 0, max: 1 });
  requireNumber('surveyRate_kmPerDay', surveyRate_kmPerDay, { gt: 0 });
  const vesselDays = (length_km / surveyRate_kmPerDay + structures * daysPerStructure) / operability;
  const campaignCost = vesselDays * dayRate + mobilisation;
  const series = new Array(years).fill(0);
  for (let y = firstYear ?? intervalYears; y <= years; y += intervalYears) series[y - 1] = campaignCost;
  return { series, campaignCost, vesselDays, annualAverage: series.reduce((s, v) => s + v, 0) / years };
}

/**
 * Escalate a base-year annual cost over a number of years: cost_k = base (1 + rate)^k, k = 0 … years − 1; a
 * rate series may be given for year-specific escalation.
 * Ref: AACE International RP 58R-10 (escalation).
 * @param {{base: number, rate: number | number[], years: number}} o currency per year, fraction per year, count
 * @returns {number[]} currency per year
 */
export function escalateSeries({ base, rate, years }) {
  requireNumber('base', base); requireNumber('years', years, { min: 1, integer: true });
  const out = []; let f = 1;
  for (let k = 0; k < years; k++) {
    if (k > 0) f *= 1 + requireNumber('rate', typeof rate === 'number' ? rate : rate[Math.min(k, rate.length - 1)], { gt: -1 });
    out.push(base * f);
  }
  return out;
}

/**
 * Steady-state (inherent) availability A = MTBF / (MTBF + MTTR), with the expected failures and downtime per year.
 * Ref: IEC 60050-192 (dependability vocabulary); Rausand & Høyland, System Reliability Theory, 2nd ed., ch. 9.
 * @param {{mtbf: number, mttr: number}} o mean time between failures and mean time to restore (hours)
 * @returns {{availability: number, failuresPerYear: number, downtimeHoursPerYear: number}}
 */
export function availability({ mtbf, mttr }) {
  requireNumber('mtbf', mtbf, { gt: 0 }); requireNumber('mttr', mttr, { min: 0 });
  const a = mtbf / (mtbf + mttr);
  return { availability: a, failuresPerYear: 8760 / (mtbf + mttr), downtimeHoursPerYear: 8760 * (1 - a) };
}

/**
 * Downtime and lost-production cost per year. Availability is taken from the caller or from MTBF/MTTR. Lost
 * margin = lost production × (price − variable cost saved); deferred production may be credited back with
 * `recoveryFraction` (e.g. reservoir volumes produced later).
 * Ref: production-loss accounting, ISO 20815:2018 (production assurance and reliability management).
 * @param {{availability?: number, mtbf?: number, mttr?: number, capacityPerHour: number, price: number,
 *   variableCostSaved?: number, recoveryFraction?: number}} o capacity in units per hour, price and saved cost in
 *   currency per unit, recoveryFraction 0–1
 * @returns {{availability: number, downtimeHours: number, lostProduction: number, lostRevenue: number,
 *   lostMargin: number}} per year
 */
export function downtimeCost({ availability: a, mtbf, mttr, capacityPerHour, price, variableCostSaved = 0, recoveryFraction = 0 }) {
  const av = a ?? availability({ mtbf, mttr }).availability;
  requireNumber('availability', av, { min: 0, max: 1 }); requireNumber('capacityPerHour', capacityPerHour, { min: 0 });
  requireNumber('price', price, { min: 0 }); requireNumber('recoveryFraction', recoveryFraction, { min: 0, max: 1 });
  const hours = 8760 * (1 - av);
  const lost = hours * capacityPerHour * (1 - recoveryFraction);
  return { availability: av, downtimeHours: hours, lostProduction: lost, lostRevenue: lost * price,
    lostMargin: lost * (price - variableCostSaved) };
}

/**
 * Cost of a failure event and its annual expectation: repair + lost production margin over the outage +
 * consequence costs (safety, environment, penalties), times the annual probability or frequency.
 * Ref: risk = probability × consequence, ISO 31000:2018; DNV-RP-G101 (risk-based inspection) for the consequence
 * categories.
 * @param {{repairCost: number, downtimeHours?: number, productionPerHour?: number, marginPerUnit?: number,
 *   consequenceCost?: number, annualFrequency?: number}} o currency, hours, units/hour, currency/unit, currency,
 *   expected events per year
 * @returns {{repair: number, lostProduction: number, consequence: number, costPerEvent: number,
 *   expectedAnnualLoss: number}} currency
 */
export function costOfFailure({ repairCost, downtimeHours = 0, productionPerHour = 0, marginPerUnit = 0, consequenceCost = 0, annualFrequency = 1 }) {
  requireNumber('repairCost', repairCost, { min: 0 }); requireNumber('downtimeHours', downtimeHours, { min: 0 });
  requireNumber('annualFrequency', annualFrequency, { min: 0 }); requireNumber('consequenceCost', consequenceCost, { min: 0 });
  const lost = downtimeHours * productionPerHour * marginPerUnit;
  const per = repairCost + lost + consequenceCost;
  return { repair: repairCost, lostProduction: lost, consequence: consequenceCost, costPerEvent: per,
    expectedAnnualLoss: per * annualFrequency };
}

/**
 * Insurance recovery for one loss: the insurer pays (loss − deductible) × (1 − coinsurance) up to the limit; the
 * remainder is retained (uninsured).
 * Ref: standard deductible / limit / coinsurance policy mechanics, e.g. Rejda & McNamara, Principles of Risk
 * Management and Insurance.
 * @param {{loss: number, deductible?: number, limit?: number, coinsurance?: number}} o currency; coinsurance =
 *   insured's share of the loss above the deductible (fraction 0–1)
 * @returns {{insured: number, retained: number}} currency
 */
export function insuranceRecovery({ loss, deductible = 0, limit = Infinity, coinsurance = 0 }) {
  requireNumber('loss', loss, { min: 0 }); requireNumber('deductible', deductible, { min: 0 });
  requireNumber('coinsurance', coinsurance, { min: 0, max: 1 });
  if (!(limit >= 0)) throw new RangeError('limit must be >= 0');
  const insured = Math.min(limit, Math.max(0, loss - deductible) * (1 - coinsurance));
  return { insured, retained: loss - insured };
}

/**
 * Life extension versus decommissioning now: incremental cash flows of extending (upgrade capex now, net operating
 * cash flows during the extension, decommissioning deferred to the end) against decommissioning today.
 * Costs are in constant money; `decomEscalation` is the real escalation of the decommissioning cost.
 * Ref: incremental cash-flow principle, Brealey, Myers & Allen ch. 6; NORSOK N-006 for the engineering side of
 * life extension.
 * @param {{upgradeCapex: number, netOperatingCashflow: number[], decomCost: number, decomEscalation?: number,
 *   rate: number}} o currency; net operating cash flow per extension year (revenue − opex, may be negative); real
 *   discount rate (fraction)
 * @returns {{extend: number[], decommissionNow: number[], incremental: number[], npvExtend: number,
 *   npvDecommissionNow: number, npvIncremental: number, decision: 'extend' | 'decommission'}}
 */
export function lifeExtensionCashflows({ upgradeCapex, netOperatingCashflow, decomCost, decomEscalation = 0, rate }) {
  requireNumber('upgradeCapex', upgradeCapex, { min: 0 }); requireSeries('netOperatingCashflow', netOperatingCashflow);
  requireNumber('decomCost', decomCost, { min: 0 }); requireNumber('rate', rate, { gt: -1 });
  const n = netOperatingCashflow.length;
  const extend = [-upgradeCapex, ...netOperatingCashflow];
  extend[n] -= decomCost * Math.pow(1 + decomEscalation, n);
  const now = new Array(n + 1).fill(0); now[0] = -decomCost;
  const incremental = extend.map((v, k) => v - now[k]);
  const ne = npv(extend, rate); const nd = npv(now, rate);
  return { extend, decommissionNow: now, incremental, npvExtend: ne, npvDecommissionNow: nd, npvIncremental: ne - nd,
    decision: ne >= nd ? 'extend' : 'decommission' };
}
