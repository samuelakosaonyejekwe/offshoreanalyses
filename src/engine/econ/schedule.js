// Project schedule: critical-path method, PERT, Monte Carlo schedule risk, weather-downtime adjustment,
// cost loading to a capex S-curve, delay cost and liquidated damages. Durations and times are in days.

import { makeRng } from '../math/rng.js';
import { dist } from '../math/stats.js';
import { requireNumber } from './tvm.js';
import { normCdf, quantile } from './risk.js';

const DAYS_PER_MONTH = 365.25 / 12;

// Validates the network once and returns a topological order with resolved links.
function compile(activities) {
  if (!Array.isArray(activities) || activities.length === 0) throw new RangeError('activities must be a non-empty array');
  const index = new Map();
  activities.forEach((a, i) => {
    if (a.id === undefined || index.has(a.id)) throw new RangeError(`activity id '${a.id}' is missing or duplicated`);
    index.set(a.id, i);
  });
  const n = activities.length;
  const preds = activities.map(() => []); const succs = activities.map(() => []);
  activities.forEach((a, i) => {
    for (const p of a.predecessors || []) {
      const link = typeof p === 'object' ? p : { id: p };
      const j = index.get(link.id);
      if (j === undefined) throw new RangeError(`activity '${a.id}' has unknown predecessor '${link.id}'`);
      const type = link.type ?? 'FS';
      if (!['FS', 'SS', 'FF', 'SF'].includes(type)) throw new RangeError(`link type '${type}' is not FS, SS, FF or SF`);
      const l = { from: j, to: i, type, lag: link.lag ?? 0 };
      preds[i].push(l); succs[j].push(l);
    }
  });
  const indeg = preds.map((p) => p.length); const order = [];
  const queue = []; for (let i = 0; i < n; i++) if (indeg[i] === 0) queue.push(i);
  while (queue.length) {
    const i = queue.shift(); order.push(i);
    for (const l of succs[i]) if (--indeg[l.to] === 0) queue.push(l.to);
  }
  if (order.length !== n) throw new RangeError('activities contain a dependency cycle');
  return { n, preds, succs, order };
}

// Forward and backward pass for given durations. adjust(i, ES, d) may stretch a duration once its start is known.
function pass(net, dur, adjust) {
  const { n, preds, succs, order } = net;
  const ES = new Float64Array(n); const EF = new Float64Array(n); const d = Float64Array.from(dur);
  for (const i of order) {
    let es = 0;
    for (const l of preds[i]) {
      const j = l.from;
      const c = l.type === 'FS' ? EF[j] + l.lag : l.type === 'SS' ? ES[j] + l.lag
        : l.type === 'FF' ? EF[j] + l.lag - d[i] : ES[j] + l.lag - d[i];
      if (c > es) es = c;
    }
    ES[i] = es;
    if (adjust) d[i] = adjust(i, es, d[i]);
    EF[i] = es + d[i];
  }
  let T = 0; for (let i = 0; i < n; i++) if (EF[i] > T) T = EF[i];
  const LS = new Float64Array(n); const LF = new Float64Array(n);
  for (let q = n - 1; q >= 0; q--) {
    const i = order[q];
    let lf = T;
    for (const l of succs[i]) {
      const s = l.to;
      const c = l.type === 'FS' ? LS[s] - l.lag : l.type === 'SS' ? LS[s] - l.lag + d[i]
        : l.type === 'FF' ? LF[s] - l.lag : LF[s] - l.lag + d[i];
      if (c < lf) lf = c;
    }
    LF[i] = lf; LS[i] = lf - d[i];
  }
  return { ES, EF, LS, LF, d, T };
}

/**
 * Critical-path method: forward and backward pass over an activity-on-node network with FS/SS/FF/SF links and
 * lags, giving early/late dates, total and free float and the critical path.
 * Ref: Kelley, J.E. & Walker, M.R. (1959), Critical-path planning and scheduling, Proc. Eastern Joint Computer
 * Conf.; Kerzner, Project Management; AACE International RP 49R-06 (identifying the critical path).
 * @param {{id: string, duration: number, predecessors?: (string | {id: string, type?: string, lag?: number})[]}[]} activities
 *   durations and lags in days (duration >= 0)
 * @returns {{duration: number, activities: {id: string, duration: number, ES: number, EF: number, LS: number,
 *   LF: number, totalFloat: number, freeFloat: number, critical: boolean}[], criticalPath: string[]}} days
 */
export function cpm(activities) {
  const net = compile(activities);
  const dur = activities.map((a) => requireNumber(`duration of '${a.id}'`, a.duration, { min: 0 }));
  return report(activities, net, pass(net, dur));
}

function report(activities, net, r) {
  const eps = 1e-9 * Math.max(1, r.T);
  const out = activities.map((a, i) => {
    let ff = r.T - r.EF[i];
    for (const l of net.succs[i]) {
      const s = l.to;
      const slack = l.type === 'FS' ? r.ES[s] - l.lag - r.EF[i] : l.type === 'SS' ? r.ES[s] - l.lag - r.ES[i]
        : l.type === 'FF' ? r.EF[s] - l.lag - r.EF[i] : r.EF[s] - l.lag - r.ES[i];
      if (slack < ff) ff = slack;
    }
    const tf = r.LS[i] - r.ES[i];
    return { id: a.id, duration: r.d[i], ES: r.ES[i], EF: r.EF[i], LS: r.LS[i], LF: r.LF[i], totalFloat: Math.abs(tf) < eps ? 0 : tf,
      freeFloat: Math.abs(ff) < eps ? 0 : Math.max(0, ff), critical: tf <= eps };
  });
  // trace one critical chain: start at the earliest critical activity and follow driving critical links
  const path = [];
  let cur = -1;
  for (const i of net.order) if (out[i].critical && (cur < 0 || r.ES[i] < r.ES[cur])) cur = i;
  while (cur >= 0) {
    path.push(activities[cur].id);
    let next = -1;
    for (const l of net.succs[cur]) {
      const s = l.to;
      if (!out[s].critical) continue;
      const drive = l.type === 'FS' ? r.EF[cur] + l.lag : l.type === 'SS' ? r.ES[cur] + l.lag
        : l.type === 'FF' ? r.EF[cur] + l.lag - r.d[s] : r.ES[cur] + l.lag - r.d[s];
      if (Math.abs(drive - r.ES[s]) <= eps && (next < 0 || r.EF[s] > r.EF[next])) next = s;
    }
    cur = next;
  }
  return { duration: r.T, activities: out, criticalPath: path };
}

/**
 * PERT three-point estimate: mean (a + 4m + b)/6, standard deviation (b − a)/6.
 * Ref: Malcolm, Roseboom, Clark & Fazar (1959), Application of a technique for research and development program
 * evaluation, Operations Research 7(5).
 * @param {{a: number, m: number, b: number}} o optimistic, most likely and pessimistic durations (days), a <= m <= b
 * @returns {{mean: number, std: number, variance: number}} days, days, days²
 */
export function pertEstimate({ a, m, b }) {
  requireNumber('a', a, { min: 0 }); requireNumber('m', m, { min: a }); requireNumber('b', b, { min: m });
  const std = (b - a) / 6;
  return { mean: (a + 4 * m + b) / 6, std, variance: std * std };
}

/**
 * PERT network analysis: CPM on the expected durations, variance of the critical path as the sum of its
 * activities' variances, and the normal-approximation probability of finishing by a target date.
 * The approximation ignores near-critical paths (merge bias); use scheduleRisk for that.
 * Ref: Malcolm et al. (1959); Kerzner, Project Management.
 * @param {{id: string, a: number, m: number, b: number, predecessors?: object[]}[]} activities three-point durations (days)
 * @param {number} [target] target completion (days)
 * @returns {{expectedDuration: number, std: number, criticalPath: string[], probabilityByTarget: number | null,
 *   schedule: object}}
 */
export function pertNetwork(activities, target) {
  const est = activities.map((a) => pertEstimate(a));
  const sched = cpm(activities.map((a, i) => ({ ...a, duration: est[i].mean })));
  let v = 0;
  for (const id of sched.criticalPath) v += est[activities.findIndex((a) => a.id === id)].variance;
  const std = Math.sqrt(v);
  const p = target === undefined ? null : std > 0 ? normCdf((target - sched.duration) / std) : (target >= sched.duration ? 1 : 0);
  return { expectedDuration: sched.duration, std, criticalPath: sched.criticalPath, probabilityByTarget: p, schedule: sched };
}

/**
 * Gross duration of weather-sensitive work: the net working time is consumed month by month at that month's
 * operability (fraction of time workable), starting from a given day of the year.
 * Ref: weather-restricted operations planning, DNV-ST-N001 (weather-restricted operations); operability from the caller's
 * workability analysis.
 * @param {{netDays: number, operabilityByMonth: number[], startDay?: number}} o net working days, 12 fractions
 *   (January … December, each in (0, 1]), start as day of year (0 = 1 January)
 * @returns {{grossDays: number, weatherDays: number}} days
 */
export function weatherAdjustedDuration({ netDays, operabilityByMonth, startDay = 0 }) {
  requireNumber('netDays', netDays, { min: 0 }); requireNumber('startDay', startDay, { min: 0 });
  if (!operabilityByMonth || operabilityByMonth.length !== 12) throw new RangeError('operabilityByMonth must have 12 values');
  for (const f of operabilityByMonth) requireNumber('operabilityByMonth', f, { gt: 0, max: 1 });
  let t = startDay; let left = netDays;
  for (let guard = 0; left > 1e-12 && guard < 100000; guard++) {
    const mIdx = Math.floor(t / DAYS_PER_MONTH);
    const f = operabilityByMonth[((mIdx % 12) + 12) % 12];
    const monthEnd = (mIdx + 1) * DAYS_PER_MONTH;
    const possible = (monthEnd - t) * f;
    if (possible >= left) { t += left / f; left = 0; } else { left -= possible; t = monthEnd; }
  }
  return { grossDays: t - startDay, weatherDays: t - startDay - netDays };
}

/**
 * Monte Carlo schedule risk analysis (seeded): activity durations are sampled from their distributions, the
 * network is recalculated for every sample, and the completion distribution, per-activity criticality index and
 * mean dates are returned. Weather hook: an activity with `operabilityByMonth` (12 fractions) is stretched from
 * its sampled start date; or pass `operability(activity, monthIndex)` for a general rule. `projectStartDay` is
 * the day of year of time zero.
 * Percentiles are cumulative: P80 is the duration not exceeded with 80 % probability.
 * Ref: AACE International RP 57R-09 (integrated cost and schedule risk analysis using Monte Carlo simulation of
 * a CPM model); Van Slyke, R.M. (1963), Monte Carlo methods and the PERT problem, Operations Research 11(5).
 * @param {{activities: object[], n?: number, seed?: number, projectStartDay?: number,
 *   operability?: (activity: object, monthIndex: number) => number}} o each activity: { id, predecessors, and one
 *   of: dist (spec of math/stats dist(), days) | a, m, b (PERT) | duration (fixed) , operabilityByMonth? }
 * @returns {{durations: Float64Array, mean: number, std: number, p10: number, p50: number, p80: number, p90: number,
 *   deterministic: number, criticalityIndex: Record<string, number>, meanStart: Record<string, number>,
 *   meanFinish: Record<string, number>, convention: string, n: number, seed: number}} days
 */
export function scheduleRisk({ activities, n = 2000, seed = 1, projectStartDay = 0, operability }) {
  requireNumber('n', n, { min: 2, integer: true });
  const net = compile(activities);
  const rng = makeRng(seed);
  const samplers = activities.map((a) => {
    if (a.dist) { const d = dist(a.dist); return { mean: d.mean, draw: () => Math.max(0, d.sample(rng)) }; }
    if (a.m !== undefined) {
      if (a.b === a.a) return { mean: a.m, draw: () => a.m };
      const d = dist({ dist: 'pert', min: a.a, mode: a.m, max: a.b });
      return { mean: d.mean, draw: () => d.sample(rng) };
    }
    requireNumber(`duration of '${a.id}'`, a.duration, { min: 0 });
    return { mean: a.duration, draw: () => a.duration };
  });
  const adjust = (i, es, d) => {
    const a = activities[i];
    if (a.operabilityByMonth) return weatherAdjustedDuration({ netDays: d, operabilityByMonth: a.operabilityByMonth, startDay: projectStartDay + es }).grossDays;
    if (operability) {
      const tab = []; for (let m = 0; m < 12; m++) tab.push(operability(a, m));
      return weatherAdjustedDuration({ netDays: d, operabilityByMonth: tab, startDay: projectStartDay + es }).grossDays;
    }
    return d;
  };
  const needAdjust = !!operability || activities.some((a) => a.operabilityByMonth);
  const N = net.n;
  const durations = new Float64Array(n); const crit = new Float64Array(N); const sES = new Float64Array(N); const sEF = new Float64Array(N);
  const d = new Float64Array(N);
  for (let k = 0; k < n; k++) {
    for (let i = 0; i < N; i++) d[i] = samplers[i].draw();
    const r = pass(net, d, needAdjust ? adjust : undefined);
    durations[k] = r.T;
    const eps = 1e-9 * Math.max(1, r.T);
    for (let i = 0; i < N; i++) { if (r.LS[i] - r.ES[i] <= eps) crit[i]++; sES[i] += r.ES[i]; sEF[i] += r.EF[i]; }
  }
  const sorted = Float64Array.from(durations).sort();
  let mean = 0; for (let k = 0; k < n; k++) mean += durations[k]; mean /= n;
  let v = 0; for (let k = 0; k < n; k++) v += (durations[k] - mean) ** 2;
  const criticalityIndex = {}; const meanStart = {}; const meanFinish = {};
  activities.forEach((a, i) => { criticalityIndex[a.id] = crit[i] / n; meanStart[a.id] = sES[i] / n; meanFinish[a.id] = sEF[i] / n; });
  const det = pass(net, samplers.map((s) => s.mean), needAdjust ? adjust : undefined).T;
  return { durations, mean, std: Math.sqrt(v / (n - 1)), p10: quantile(sorted, 0.1), p50: quantile(sorted, 0.5),
    p80: quantile(sorted, 0.8), p90: quantile(sorted, 0.9), deterministic: det, criticalityIndex, meanStart, meanFinish,
    convention: 'Cumulative: Pxx is the duration not exceeded with xx % probability.', n, seed };
}

/**
 * Cost loading of a schedule: each activity's cost is spread uniformly over its early (or late) dates and summed
 * per period, giving the expenditure S-curve and the phasing fractions used by the cash-flow model.
 * Ref: AACE International RP 55R-09 (analysing S-curves); Kerzner, Project Management.
 * @param {{schedule: {duration: number, activities: object[]}, costs: Record<string, number>, periodDays?: number,
 *   timing?: 'early' | 'late'}} o schedule from cpm(), cost per activity id (currency), period length in days
 *   (default one month; use 365.25 for annual capex phasing)
 * @returns {{periods: number[], cost: number[], cumulative: number[], fractions: number[], total: number}}
 *   periods are 1-based period numbers; cumulative is the S-curve (currency)
 */
export function costLoading({ schedule, costs, periodDays = DAYS_PER_MONTH, timing = 'early' }) {
  if (!schedule || !Array.isArray(schedule.activities)) throw new RangeError('schedule must be a cpm() result');
  requireNumber('periodDays', periodDays, { gt: 0 });
  const nP = Math.max(1, Math.ceil(schedule.duration / periodDays - 1e-9));
  const cost = new Array(nP).fill(0); let total = 0;
  for (const a of schedule.activities) {
    const c = costs[a.id] ?? 0;
    requireNumber(`cost of '${a.id}'`, c, { min: 0 });
    if (c === 0) continue;
    total += c;
    const s = timing === 'late' ? a.LS : a.ES; const f = timing === 'late' ? a.LF : a.EF;
    if (f - s < 1e-12) { cost[Math.min(nP - 1, Math.floor(s / periodDays))] += c; continue; }
    for (let p = Math.floor(s / periodDays); p < nP && p * periodDays < f; p++) {
      const overlap = Math.min(f, (p + 1) * periodDays) - Math.max(s, p * periodDays);
      if (overlap > 0) cost[p] += c * overlap / (f - s);
    }
  }
  const cumulative = []; let run = 0;
  for (const v of cost) { run += v; cumulative.push(run); }
  return { periods: cost.map((_, i) => i + 1), cost, cumulative, fractions: cost.map((v) => (total > 0 ? v / total : 0)), total };
}

/**
 * Cost of a delay: standby of the marine spread and site team, plus revenue either lost for the delay period or
 * deferred (in which case the loss is the fall in present value of the operating cash flows shifted by the delay).
 * Ref: delay and disruption costing practice, Society of Construction Law Delay and Disruption Protocol (2017).
 * @param {{delayDays: number, standbyDayRates?: number[], revenuePerDay?: number, mode?: 'lost' | 'deferred',
 *   pvOfOperations?: number, discountRate?: number}} o standby rates (currency/day), revenue (currency/day) for
 *   mode 'lost'; for 'deferred' the PV of operating cash flows (currency) and the annual discount rate
 * @returns {{standby: number, revenueImpact: number, total: number}} currency
 */
export function delayCost({ delayDays, standbyDayRates = [], revenuePerDay = 0, mode = 'lost', pvOfOperations, discountRate }) {
  requireNumber('delayDays', delayDays, { min: 0 });
  let rate = 0; for (const r of standbyDayRates) rate += requireNumber('standbyDayRates', r, { min: 0 });
  const standby = rate * delayDays;
  let revenueImpact;
  if (mode === 'lost') revenueImpact = requireNumber('revenuePerDay', revenuePerDay, { min: 0 }) * delayDays;
  else if (mode === 'deferred') {
    requireNumber('pvOfOperations', pvOfOperations); requireNumber('discountRate', discountRate, { gt: -1 });
    revenueImpact = pvOfOperations * (1 - Math.pow(1 + discountRate, -delayDays / 365.25));
  } else throw new RangeError("mode must be 'lost' or 'deferred'");
  return { standby, revenueImpact, total: standby + revenueImpact };
}

/**
 * Liquidated damages for late completion: daily rate × days beyond the grace period, capped at a fraction of the
 * contract value.
 * Ref: typical delay-damages clause, e.g. FIDIC Conditions of Contract (2017) Sub-Clause 8.8.
 * @param {{delayDays: number, ratePerDay: number, graceDays?: number, contractValue?: number, capFraction?: number}} o
 *   days, currency/day, days, currency, cap as a fraction of contract value
 * @returns {{damages: number, cap: number, capped: boolean, chargeableDays: number}} currency
 */
export function liquidatedDamages({ delayDays, ratePerDay, graceDays = 0, contractValue, capFraction }) {
  requireNumber('delayDays', delayDays, { min: 0 }); requireNumber('ratePerDay', ratePerDay, { min: 0 });
  requireNumber('graceDays', graceDays, { min: 0 });
  const days = Math.max(0, delayDays - graceDays);
  const cap = contractValue !== undefined && capFraction !== undefined
    ? requireNumber('contractValue', contractValue, { min: 0 }) * requireNumber('capFraction', capFraction, { min: 0 }) : Infinity;
  const raw = days * ratePerDay;
  return { damages: Math.min(raw, cap), cap, capped: raw > cap, chargeableDays: days };
}
