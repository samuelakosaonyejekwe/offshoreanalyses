import test from 'node:test';
import assert from 'node:assert/strict';
import * as L from '../../src/engine/econ/lifecycle.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const pvA = (i, n) => (1 - (1 + i) ** -n) / i;

test('life-cycle cost by phase, total cost of ownership and equivalent annual cost', () => {
  // capex 1000 at t0; opex 100 in years 1–3; decommissioning 200 less residual 50 in year 3; 10 %
  const r = L.lifeCycleCost({ phases: [{ name: 'capex', costs: [1000, 0, 0, 0] }, { name: 'opex', costs: [0, 100, 100, 100] }, { name: 'decom', costs: [0, 0, 0, 150] }], rate: 0.1 });
  close(r.total, 1450, 1e-9); close(r.discounted, 1000 + 100 * pvA(0.1, 3) + 150 / 1.331, 1e-9);
  close(r.eac, r.discounted / pvA(0.1, 3), 1e-9); assert.equal(r.years, 3);
  close(r.byPhase[1].discounted, 100 * pvA(0.1, 3), 1e-9); close(r.byPhase.reduce((s, p) => s + p.share, 0), 1, 1e-12);
  // TCO: 500 + (60 + 30 + 10) (P/A, 8 %, 5) + (20 − 80)/1.08⁵
  const t = L.totalCostOfOwnership({ acquisition: 500, operating: 60, maintenance: 30, downtime: 10, disposal: 20, residual: 80, life: 5, rate: 0.08 });
  close(t.tco, 500 + 100 * pvA(0.08, 5) - 60 / 1.08 ** 5, 1e-9); close(t.eac, t.tco / pvA(0.08, 5), 1e-9);
  // escalation equal to the discount rate: every year's PV is cost/(1 + i)
  close(L.totalCostOfOwnership({ acquisition: 0, operating: 100, life: 4, rate: 0.05, escalation: 0.05 }).tco, 400 / 1.05, 1e-9);
});

test('EAC of unequal lives and economic service life (textbook cases)', () => {
  // A: 10 000, 5 yr, 2000/yr, salvage 1000; B: 15 000, 10 yr, 1500/yr, salvage 0; 10 %
  //   EAC_A = 9000 × 0.263797 + 100 + 2000 = 4474.18; EAC_B = 15000 × 0.162745 + 1500 = 3941.18
  const e = L.eacComparison([{ name: 'A', capex: 10000, annualCost: 2000, salvage: 1000, life: 5 }, { name: 'B', capex: 15000, annualCost: 1500, life: 10 }], 0.1);
  close(e.results[0].eac, 4474.18, 0.01); close(e.results[1].eac, 3941.18, 0.01); assert.equal(e.best, 'B');
  close(e.results[0].capitalRecovery, 9000 / pvA(0.1, 5) + 100, 1e-8);
  // economic life: P = 20 000; O&M 2000, 3000, 4620, 8000, 12000; salvage 15000, 11250, 8500, 6500, 4750; 10 %
  // (Park, Contemporary Engineering Economics, lift-truck example pattern) — computed from the definition
  const om = [2000, 3000, 4620, 8000, 12000]; const sv = [15000, 11250, 8500, 6500, 4750];
  const s = L.economicServiceLife({ capex: 20000, om, salvage: sv, rate: 0.1 });
  const eac = om.map((_, k) => { const n = k + 1; let pv = 20000 - sv[k] / 1.1 ** n; for (let j = 0; j < n; j++) pv += om[j] / 1.1 ** (j + 1); return pv / pvA(0.1, n); });
  for (let k = 0; k < 5; k++) close(s.eac[k], eac[k], 1e-6);
  close(s.eac[0], 20000 * 1.1 - 15000 + 2000, 1e-9); // one year: 9000
  assert.equal(s.economicLife, eac.indexOf(Math.min(...eac)) + 1); close(s.minEac, Math.min(...eac), 1e-6);
  close(s.capitalRecovery[2] + s.omEquivalent[2], s.eac[2], 1e-8);
});

test('replacement analysis: defender marginal cost against the challenger minimum EAC', () => {
  // defender worth 5000 now; next years: market values 4000, 3000; O&M 2000, 3500; 10 %
  //   MC1 = 5000 × 1.1 − 4000 + 2000 = 3500; MC2 = 4000 × 1.1 − 3000 + 3500 = 4900
  const ch = { capex: 12000, om: [1000, 1000, 1000, 1000], salvage: [9000, 7000, 5500, 4000] };
  const r = L.replacementAnalysis({ defender: { marketValue: 5000, om: [2000, 3500], salvage: [4000, 3000] }, challenger: ch, rate: 0.1 });
  close(r.defenderMarginalCost[0], 3500, 1e-9); close(r.defenderMarginalCost[1], 4900, 1e-9);
  const cEac = L.economicServiceLife({ ...ch, rate: 0.1 });
  close(r.challengerMinEac, cEac.minEac, 1e-9); assert.equal(r.challengerLife, cEac.economicLife);
  // challenger minimum EAC lies between 3500 and 4900 → keep one more year, then replace
  assert.ok(r.challengerMinEac > 3500 && r.challengerMinEac < 4900, `challenger EAC ${r.challengerMinEac}`);
  assert.equal(r.replaceNow, false); assert.equal(r.keepYears, 1);
  const now = L.replacementAnalysis({ defender: { marketValue: 5000, om: [6000], salvage: [1000] }, challenger: ch, rate: 0.1 });
  assert.equal(now.replaceNow, true); assert.equal(now.keepYears, 0);
});

test('lifetime distributions: Weibull, exponential and lognormal closed forms', () => {
  // β = 2, η = 1000 h at t = 500: R = e^(−0.25) = 0.7788, h = 2 × 500/10⁶ = 0.001, MTTF = 1000 Γ(1.5) = 886.23
  const w = L.weibullLife({ t: 500, shape: 2, scale: 1000 });
  close(w.reliability, Math.exp(-0.25), 1e-12); close(w.hazard, 0.001, 1e-12); close(w.mttf, 886.2269, 1e-3);
  close(w.reliability + w.unreliability, 1, 1e-15);
  // β = 1 is the exponential law
  const e = L.exponentialLife({ t: 500, rate: 0.001 });
  close(L.weibullLife({ t: 500, shape: 1, scale: 1000 }).reliability, e.reliability, 1e-12);
  close(e.reliability, Math.exp(-0.5), 1e-15); close(e.mttf, 1000, 1e-12); close(e.hazard, 0.001, 1e-15);
  // lognormal: median e^μ has R = 0.5; MTTF = e^(μ + σ²/2)
  const l = L.lognormalLife({ t: Math.exp(6), mu: 6, sigma: 0.5 });
  close(l.reliability, 0.5, 1e-9); close(l.mttf, Math.exp(6.125), 1e-9);
  // hazard = pdf/R at the median: φ(0)/(σ t)/0.5
  close(l.hazard, 0.3989423 / (0.5 * Math.exp(6)) / 0.5, 1e-9);
  assert.throws(() => L.weibullLife({ t: 1, shape: 0, scale: 1 }), RangeError);
});

test('maintenance policies: age replacement optimum, renewal function limits, block replacement, power-law NHPP', () => {
  // exponential life (β = 1): no ageing, preventive replacement never pays — run to failure, cost rate cf/MTTF
  const ex = L.ageReplacement({ shape: 1, scale: 100, costPreventive: 1, costFailure: 5 });
  assert.equal(ex.optimalAge, null); close(ex.runToFailureCostRate, 0.05, 1e-9); close(ex.costRate, 0.05, 1e-6);
  // β = 2, η = 100, cp = 1, cf = 10: optimum from dC/dT = 0, h(T) ∫R − F = cp/(cf − cp) → checked on the returned curve
  const a = L.ageReplacement({ shape: 2, scale: 100, costPreventive: 1, costFailure: 10 });
  assert.ok(a.optimalAge > 25 && a.optimalAge < 45, `optimal age ${a.optimalAge}`);
  close(a.runToFailureCostRate, 10 / (100 * 0.8862269), 1e-6);
  assert.ok(a.costRate < a.runToFailureCostRate && Math.abs(a.saving - (a.runToFailureCostRate - a.costRate)) < 1e-12);
  // independent evaluation of C(T) at the optimum by Simpson integration of R(t) = exp(−(t/100)²)
  const T = a.optimalAge; const n = 2000; let integ = 0;
  for (let i = 0; i <= n; i++) { const t = T * i / n; integ += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * Math.exp(-((t / 100) ** 2)); }
  integ *= T / n / 3;
  const RT = Math.exp(-((T / 100) ** 2));
  rel(a.costRate, (1 * RT + 10 * (1 - RT)) / integ, 1e-4);
  // renewal function of an exponential life is exactly t/MTTF; for any life M(t) → t/MTTF − ½ + σ²/(2 MTTF²) eventually
  rel(L.renewalFunction({ t: 300, shape: 1, scale: 100 }).expectedFailures, 3, 2e-3);
  const m = L.renewalFunction({ t: 1000, shape: 2, scale: 100, steps: 2000 });
  const mu = 88.62269; const var2 = 100 ** 2 * (1 - Math.PI / 4);
  rel(m.expectedFailures, 1000 / mu + var2 / (2 * mu * mu) - 0.5, 5e-3);
  // short horizon: M(t) ≈ F(t)
  rel(L.renewalFunction({ t: 10, shape: 2, scale: 100 }).expectedFailures, 1 - Math.exp(-0.01), 1e-2);
  const b = L.blockReplacement({ interval: 50, shape: 2, scale: 100, costPreventive: 1, costFailure: 10 });
  close(b.costRate, (1 + 10 * b.expectedFailures) / 50, 1e-12);
  rel(b.expectedFailures, L.renewalFunction({ t: 50, shape: 2, scale: 100 }).expectedFailures, 1e-9);
  // NHPP: N(t) = 0.01 t^1.5 at t = 100 → 10 failures; intensity 0.015 √100 = 0.15; T* = [cr/(cm λ (β − 1))]^(1/β)
  const p = L.nhppPowerLaw({ lambda: 0.01, beta: 1.5, t: 100, minimalRepairCost: 2, replacementCost: 50 });
  close(p.expectedFailures, 10, 1e-9); close(p.intensity, 0.15, 1e-12); close(p.expectedRepairCost, 20, 1e-9);
  const Ts = (50 / (2 * 0.01 * 0.5)) ** (1 / 1.5);
  close(p.optimalReplacementAge, Ts, 1e-6); close(p.minCostRate, (50 + 2 * 0.01 * Ts ** 1.5) / Ts, 1e-9);
  assert.equal(L.nhppPowerLaw({ lambda: 0.01, beta: 0.8, t: 100, minimalRepairCost: 2, replacementCost: 50 }).optimalReplacementAge, null);
});

test('Markov models: two-state availability μ/(λ + μ), chain evolution, semi-Markov time fractions', () => {
  // failure rate 0.01/h, repair rate 0.5/h → availability 0.5/0.51 = 0.98039
  const s = L.markovSteadyState({ rates: [[0, 0.01], [0.5, 0]], upStates: [0] });
  close(s.availability, 0.5 / 0.51, 1e-12); close(s.probabilities[0] + s.probabilities[1], 1, 1e-12);
  // three-state deterioration chain
  const P = [[0.9, 0.1, 0], [0, 0.8, 0.2], [0, 0, 1]];
  const c = L.markovChainEvolve({ P, p0: [1, 0, 0], steps: 2 });
  assert.equal(c.history.length, 3);
  close(c.final[0], 0.81, 1e-12); close(c.final[1], 0.09 + 0.08, 1e-12); close(c.final[2], 0.02, 1e-12);
  // alternating up/down with mean sojourns 98 h and 2 h: availability 0.98 whatever the distributions
  const sm = L.semiMarkovSteadyState({ P: [[0, 1], [1, 0]], meanSojourn: [98, 2], upStates: [0] });
  close(sm.availability, 0.98, 1e-12); close(sm.embedded[0], 0.5, 1e-9);
});

test('system reliability: block diagram, fault tree and event tree', () => {
  // series 0.9 × 0.95; parallel 1 − 0.1 × 0.2; 2-out-of-3 of 0.9: 3 × 0.81 × 0.1 + 0.729 = 0.972
  close(L.rbdReliability({ type: 'series', blocks: [0.9, 0.95] }), 0.855, 1e-12);
  close(L.rbdReliability({ type: 'parallel', blocks: [0.9, 0.8] }), 0.98, 1e-12);
  close(L.rbdReliability({ type: 'koon', k: 2, blocks: [0.9, 0.9, 0.9] }), 0.972, 1e-12);
  close(L.rbdReliability({ type: 'series', blocks: [0.99, { type: 'parallel', blocks: [0.9, { type: 'koon', k: 1, blocks: [0.5, 0.5] }] }] }), 0.99 * (1 - 0.1 * 0.25), 1e-12);
  // non-identical 2-out-of-3: 0.9, 0.8, 0.7 → abc + ab(1−c) + a(1−b)c + (1−a)bc = 0.504 + 0.216 + 0.126 + 0.056 = 0.902
  close(L.rbdReliability({ type: 'koon', k: 2, blocks: [0.9, 0.8, 0.7] }), 0.902, 1e-12);
  // fault tree: OR(0.01, AND(0.1, 0.2)) = 1 − 0.99 × 0.98 = 0.0298; 2-of-3 vote of 0.1 = 0.028
  close(L.faultTreeProbability({ gate: 'OR', inputs: [0.01, { gate: 'AND', inputs: [0.1, 0.2] }] }), 0.0298, 1e-12);
  close(L.faultTreeProbability({ gate: 'VOTE', k: 2, inputs: [0.1, 0.1, 0.1] }), 0.028, 1e-12);
  // event tree: 0.1 /yr leak; detection fails 0.1, isolation fails 0.2
  const et = L.eventTree({ frequency: 0.1, barriers: [{ name: 'detection', pFail: 0.1 }, { name: 'isolation', pFail: 0.2 }],
    cost: (f) => (f[0] && f[1] ? 1e7 : f[0] || f[1] ? 1e6 : 0) });
  assert.equal(et.sequences.length, 4);
  close(et.sequences.reduce((s, q) => s + q.frequency, 0), 0.1, 1e-15);
  close(et.sequences.find((q) => q.failed[0] && q.failed[1]).frequency, 0.002, 1e-15);
  // EAL = 0.002 × 1e7 + (0.1 × 0.1 × 0.8 + 0.1 × 0.9 × 0.2) × 1e6 = 20 000 + 26 000
  close(et.expectedAnnualLoss, 46000, 1e-6);
});

test('inspection economics, repair or replace, and the life-extension case', () => {
  const base = { horizon: 20, inspectionCost: 50, failureCost: 5000, effectiveness: 0.8, rate: 0.05, weibull: { shape: 3, scale: 30 } };
  const r = L.rbiEconomics({ ...base, intervals: [2, 5, 10, 20] });
  assert.equal(r.results.length, 4);
  // more frequent inspection costs more and leaves less risk
  for (let i = 1; i < 4; i++) { assert.ok(r.results[i].pvInspection < r.results[i - 1].pvInspection); assert.ok(r.results[i].pvRisk > r.results[i - 1].pvRisk - 1e-9); }
  for (const x of r.results) close(x.pvTotal, x.pvInspection + x.pvRisk, 1e-9);
  assert.equal(r.best.pvTotal, Math.min(...r.results.map((x) => x.pvTotal)));
  // without inspection the yearly failure probability is the conditional Weibull increment: Σ ΔF/R × cost × DF
  let pv = 0; const R = (t) => Math.exp(-((t / 30) ** 3));
  for (let y = 1; y <= 20; y++) pv += (1 - R(y) / R(y - 1)) * 5000 / 1.05 ** y;
  rel(r.noInspection.pvRisk, pv, 1e-6);
  // useless inspection (effectiveness 0) only adds its cost
  const z = L.rbiEconomics({ ...base, effectiveness: 0, intervals: [5] });
  rel(z.results[0].pvRisk, pv, 1e-6);
  // repair 300 for 3 more years at O&M 100 against a new unit 1000 for 10 years at O&M 40; 8 %
  const rr = L.repairVsReplace({ repairCost: 300, repairLife: 3, omOld: 100, replaceCost: 1000, newLife: 10, omNew: 40, rate: 0.08 });
  close(rr.eacRepair, 300 / pvA(0.08, 3) + 100, 1e-9); close(rr.eacReplace, 1000 / pvA(0.08, 10) + 40, 1e-9);
  assert.equal(rr.decision, rr.eacRepair <= rr.eacReplace ? 'repair' : 'replace');
  close(rr.breakEvenRepairCost, (rr.eacReplace - 100) * pvA(0.08, 3), 1e-8);
  // life extension without ageing risk: −capex + revenue (P/A) − decom/(1 + i)ⁿ against −decom now
  const le = L.lifeExtensionCase({ currentAge: 25, extensionYears: 5, upgradeCapex: 100, netRevenue: 60, decomCost: 200, rate: 0.07 });
  close(le.npvDecommissionNow, -200, 1e-9); close(le.npvExtend, -100 + 60 * pvA(0.07, 5) - 200 / 1.07 ** 5, 1e-8);
  close(le.npvIncremental, le.npvExtend - le.npvDecommissionNow, 1e-9); assert.equal(le.decision, 'extend'); assert.equal(le.optimalExtension, 5);
  // with ageing risk the annual failure probability is the conditional Weibull probability at the current age
  const ag = L.lifeExtensionCase({ currentAge: 25, extensionYears: 3, upgradeCapex: 100, netRevenue: 60, failureCost: 1000, weibull: { shape: 3, scale: 30 }, decomCost: 200, rate: 0.07 });
  close(ag.annualFailureProbability[0], 1 - R(26) / R(25), 1e-9); close(ag.expectedFailureCost[1], 1000 * (1 - R(27) / R(26)), 1e-6);
  assert.ok(ag.npvExtend < L.lifeExtensionCase({ currentAge: 25, extensionYears: 3, upgradeCapex: 100, netRevenue: 60, decomCost: 200, rate: 0.07 }).npvExtend);
});

test('decommissioning: estimate build-up, retirement obligation, funding schedules, timing and residual value', () => {
  const d = L.decomEstimate({ removal: 50, wells: { count: 4, costPerWell: 10 }, subseaRecovery: 15, disposal: 8, restoration: 2, monitoring: { years: 5, costPerYear: 1 },
    scrap: { mass_t: 10000, recoverableFraction: 0.9, pricePerT: 0.0002 }, contingency: 0.2 });
  close(d.gross, 50 + 40 + 15 + 8 + 2 + 5, 1e-9); close(d.contingency, 24, 1e-9); close(d.recyclingCredit, 1.8, 1e-9);
  close(d.net, 120 + 24 - 1.8, 1e-9);
  // ARO: 100 today, 2 % inflation, 10 years, 5 %: future 121.899, PV 74.836; accretion opens at PV and closes at the future cost
  const a = L.assetRetirementObligation({ costToday: 100, inflation: 0.02, years: 10, discountRate: 0.05 });
  close(a.futureCost, 100 * 1.02 ** 10, 1e-9); close(a.presentValue, a.futureCost / 1.05 ** 10, 1e-9);
  close(a.accretion[0].opening, a.presentValue, 1e-9); close(a.accretion[0].expense, 0.05 * a.presentValue, 1e-9);
  close(a.accretion[9].closing, a.futureCost, 1e-8); assert.equal(a.accretion.length, 10);
  // sinking fund: target 1000 in 5 years at 4 % → A = 1000 × 0.04/(1.04⁵ − 1) = 184.627
  const f = L.decomFunding({ target: 1000, years: 5, fundRate: 0.04 });
  close(f.contributions[0], 1000 * 0.04 / (1.04 ** 5 - 1), 1e-9); close(f.balance[4], 1000, 1e-8);
  close(f.totalContributed + f.interestEarned, 1000, 1e-8);
  // unit of production: contributions proportional to production and compounding to the target
  const u = L.decomFunding({ target: 1000, years: 3, fundRate: 0.04, method: 'unitOfProduction', production: [30, 20, 10] });
  close(u.contributions[0] / u.contributions[2], 3, 1e-9); close(u.balance[2], 1000, 1e-8);
  // timing: escalation 3 % below the discount rate 8 % favours waiting unless holding costs are high
  const t = L.decomTiming({ costToday: 100, escalation: 0.03, rate: 0.08, candidateYears: [0, 5, 10] });
  close(t.options[1].pvCost, 100 * (1.03 / 1.08) ** 5, 1e-9); assert.equal(t.best.year, 10);
  const th = L.decomTiming({ costToday: 100, escalation: 0.03, rate: 0.08, holdingCostPerYear: 10, candidateYears: [0, 5, 10] });
  close(th.options[1].pvCost, 100 * (1.03 / 1.08) ** 5 + 10 * pvA(0.08, 5), 1e-9); assert.equal(th.best.year, 0);
  // residual values
  close(L.residualValue({ method: 'linear', cost: 1000, life: 20, age: 5, endValue: 100 }), 1000 - 900 * 0.25, 1e-9);
  close(L.residualValue({ method: 'declining', cost: 1000, age: 3, declineRate: 0.2 }), 512, 1e-9);
  close(L.residualValue({ method: 'scrap', mass_t: 5000, scrapPricePerT: 250, recoverableFraction: 0.9, removalCost: 1.5e6 }), 5000 * 250 * 0.9 - 1.5e6, 1e-6);
});
