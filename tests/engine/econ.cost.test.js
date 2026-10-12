import test from 'node:test';
import assert from 'node:assert/strict';
import * as cost from '../../src/engine/econ/cost.js';

const close = (a, b, tol = 1e-9, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);

test('every default rate is a labelled placeholder with basis year, unit and source text', () => {
  const keys = Object.keys(cost.defaultRates);
  assert.ok(keys.length > 50);
  for (const k of keys) {
    const r = cost.defaultRates[k];
    assert.ok(Number.isFinite(r.value) && r.value >= 0, k);
    assert.match(r.source, /indicative placeholder — replace with current quotations\/indices/, k);
    assert.equal(typeof r.basisYear, 'number', k); assert.ok(r.unit.length > 0, k);
  }
});

test('rates can be overridden by number or object, and unknown keys are rejected', () => {
  assert.equal(cost.resolveRate('steel.jacket').overridden, false);
  assert.equal(cost.resolveRate('steel.jacket', { 'steel.jacket': 5200 }).value, 5200);
  const o = cost.resolveRate('steel.jacket', { 'steel.jacket': { value: 6000, source: 'quote 2025-03', basisYear: 2025 } });
  assert.equal(o.value, 6000); assert.equal(o.source, 'quote 2025-03'); assert.equal(o.overridden, true);
  assert.throws(() => cost.resolveRate('steel.unobtainium'), RangeError);
});

test('bottom-up estimate by hand: direct 250 000, indirects 10 %, contingency 10 %, 3 % escalation for 2 years', () => {
  const r = cost.bottomUpEstimate({ items: [{ name: 'steel', quantity: 100, unitRate: 2000, discipline: 'structural', category: 'materials' },
    { name: 'welding', quantity: 50, unitRate: 1000, discipline: 'structural', category: 'fabrication' }],
  indirects: { engineering: 0.06, projectManagement: 0.04 }, contingency: 0.10, escalation: { rate: 0.03, years: 2 } });
  close(r.direct, 250000); close(r.indirectTotal, 25000); close(r.base, 275000); close(r.contingency, 27500);
  close(r.escalation, 302500 * (1.03 ** 2 - 1), 1e-12); close(r.total, 302500 * 1.0609, 1e-12);
  close(r.byDiscipline.structural, 250000); close(r.byCategory.fabrication, 50000);
  // consistency: total = direct + indirect + contingency + escalation
  close(r.total, r.direct + r.indirectTotal + r.contingency + r.escalation, 1e-12);
});

test('capacity scaling: six-tenths rule, doubling capacity raises cost by 2^0.6 = 1.5157', () => {
  close(cost.capacityScaling({ C1: 1e6, S1: 100, S2: 200 }), 1e6 * Math.pow(2, 0.6), 1e-14);
  close(cost.capacityScaling({ C1: 1e6, S1: 100, S2: 200 }), 1515716.57, 1e-8);
  close(cost.capacityScaling({ C1: 1e6, S1: 100, S2: 200, index1: 500, index2: 600, locationFactor: 1.1, complexityFactor: 1.2 }),
    1e6 * Math.pow(2, 0.6) * 1.2 * 1.1 * 1.2, 1e-13);
  close(cost.scalingExponent({ C1: 1e6, S1: 100, C2: 1515716.57, S2: 200 }), 0.6, 1e-8);
});

test('cost index interpolation and factor', () => {
  const f = cost.indexFactor({ years: [2018, 2020, 2022], values: [100, 110, 132] }, 2019, 2021);
  close(f.from, 105); close(f.to, 121); close(f.factor, 121 / 105, 1e-14);
  assert.throws(() => cost.indexFactor({ years: [2018, 2020], values: [100, 110] }, 2017, 2019), /outside/);
});

test('factorial estimates: Lang 4.74 for fluid processing, Hand factors by equipment type', () => {
  assert.deepEqual({ ...cost.LANG_FACTORS }, { solid: 3.10, solidFluid: 3.63, fluid: 4.74 });
  close(cost.factorialEstimate({ equipment: [{ name: 'a', cost: 100 }, { name: 'b', cost: 50 }] }).total, 150 * 4.74, 1e-12);
  const h = cost.factorialEstimate({ method: 'hand', equipment: [{ name: 'pump', cost: 100, type: 'pump' }, { name: 'comp', cost: 200, type: 'compressor' }] });
  close(h.total, 100 * 4 + 200 * 2.5); close(h.factor, 3);
  close(cost.factorialEstimate({ equipment: [{ name: 'a', cost: 100 }], factors: { fluid: 6 } }).total, 600);
});

test('monopile tonnage from geometry: pi (D - t) t L rho, D = 8 m, t = 80 mm, L = 70 m -> 1093.8 t', () => {
  const r = cost.monopileCost({ D: 8, t: 0.08, L: 70, tp: { mass_t: 300 }, rates: { 'steel.monopile': 2000, 'steel.transitionPiece': 3000 } });
  close(r.monopile_t, Math.PI * 7.92 * 0.08 * 70 * 7.85, 1e-12); close(r.monopile_t, 1093.79, 1e-5);
  close(r.total, r.monopile_t * 2000 + 900000, 1e-12);
  assert.equal(r.usesPlaceholderRates, false);
  assert.equal(cost.monopileCost({ D: 8, t: 0.08, L: 70 }).usesPlaceholderRates, true);
});

test('parametric models are quantity x rate, line by line', () => {
  const R = { 'steel.jacket': 4000, 'steel.piles': 2000 };
  close(cost.jacketCost({ jacket_t: 1000, piles_t: 500, rates: R }).total, 5e6);
  close(cost.gbsCost({ concrete_m3: 10000, rebarRatio_kgm3: 250, rates: { 'concrete.gbs': 600, 'steel.rebar': 1200, 'ballast.solid': 0 } }).total, 6e6 + 2500 * 1200);
  close(cost.mooringCost({ segments: [{ type: 'chain', length_m: 1000, massPerM_kg: 200, count: 3 }], anchors: [{ type: 'drag', count: 3 }],
    rates: { 'mooring.chain': 3000, 'anchor.drag': 100000 } }).total, 600 * 3000 + 300000);
  close(cost.cableCost({ kind: 'array', cls: '66kV', length_km: 10, terminations: 4,
    rates: { 'cable.array.66kV': 400000, 'cable.layBurial.array': 200000, 'cable.termination.array': 100000 } }).total, 4e6 + 2e6 + 4e5);
  close(cost.generationPlantCost({ substationMW: 500, turbineMW: 500, rates: { 'substation.hvac': 2e5, 'turbine.supply': 1.5e6 } }).total, 1e8 + 7.5e8);
  close(cost.marineCivilCost({ deckArea_m2: 1000, rock_t: 5000, rates: { 'civil.deck': 3000, 'civil.rock': 40 } }).total, 3e6 + 2e5);
  close(cost.subseaStructureCost({ mass_t: 100, rates: { 'steel.subseaStructure': 10000 } }).total, 1e6);
  close(cost.umbilicalFlexibleCost({ kind: 'umbilical', cls: 'static', length_m: 2000, rates: { 'umbilical.static': 800 } }).total, 1.6e6);
  close(cost.floaterCost({ type: 'semi', steel_t: 3000, rates: { 'steel.floaterSemi': 4000 } }).total, 1.2e7);
  close(cost.topsidesCost({ structural_t: 100, rates: { 'steel.topsidesStructural': 9000 } }).total, 9e5);
  const sb = cost.suctionBucketCost({ D: 10, L: 10, tSkirt: 0.03, tLid: 0.04, n: 3, lidStiffeningFactor: 1, rates: { 'steel.suctionBucket': 1 } });
  close(sb.steel_t, 3 * (Math.PI * 9.97 * 0.03 * 10 * 7.85 + Math.PI * 25 * 0.04 * 7.85), 1e-12);
});

test('installation campaign: gross days = net days / operability, plus mobilisation', () => {
  const r = cost.installationCampaignCost({ spread: [{ vessel: 'heavyLift' }, { vessel: 'tug', count: 2 }], netDays: 20, operability: 0.8,
    mobDemob: 1e6, rates: { 'vessel.heavyLift': 300000, 'vessel.tug': 10000 } });
  close(r.grossDays, 25); close(r.weatherDays, 5); close(r.spreadDayRate, 320000);
  close(r.total, 25 * 320000 + 1e6);
});

test('pipeline cost: linepipe tonnage, lay days from lay rate and operability, supply/installation split', () => {
  const r = cost.pipelineCost({ length_m: 30000, D: 0.4064, t: 0.0159, layRate_mPerDay: 3000, operability: 0.8, trenchLength_m: 5000, mobDemob: 0,
    rates: { 'pipe.linepipe': 1500, 'pipe.coatingAnticorrosion': 50, 'vessel.pipelay': 400000, 'pipe.trenching': 100, 'pipe.rockDump': 50 } });
  const tonnes = Math.PI * (0.4064 - 0.0159) * 0.0159 * 30000 * 7.85;
  close(r.linepipe_t, tonnes, 1e-12); close(r.layDays, 12.5);
  close(r.supply, tonnes * 1500 + 30000 * 50, 1e-12);
  close(r.installation, 12.5 * 400000 + 5000 * 100, 1e-12);
});

test('commissioning and decommissioning allowances', () => {
  const r = cost.commissioningDecomCost({ installedCost: 1e8, installationCost: 2e7, rates: { 'allowance.commissioning': 0.02, 'allowance.decommissioning': 0.5 } });
  close(r.commissioning, 2e6); close(r.decommissioning, 1e7);
  assert.equal(cost.commissioningDecomCost({ decomItems: [{ name: 'removal', cost: 5 }, { name: 'disposal', cost: 2 }] }).decommissioning, 7);
});

test('Wright learning curve: doubling quantity multiplies the cumulative average by the progress ratio', () => {
  for (const pr of [0.8, 0.9, 0.95]) for (const n of [1, 4, 25]) {
    const a = cost.wrightCurve({ T1: 100, n, progressRatio: pr }); const b = cost.wrightCurve({ T1: 100, n: 2 * n, progressRatio: pr });
    close(b.average / a.average, pr, 1e-13);
  }
  const w = cost.wrightCurve({ T1: 100, n: 4, progressRatio: 0.8 });
  close(w.average, 64, 1e-12); close(w.total, 256, 1e-12);
});

test('Crawford unit curve: unit 4 on an 80 % curve costs 0.64 T1; total by exact summation', () => {
  const c = cost.crawfordCurve({ T1: 100, n: 4, progressRatio: 0.8 });
  close(c.unit, 64, 1e-12);
  close(c.total, 100 * (1 + 0.8 + Math.pow(3, Math.log2(0.8)) + 0.64), 1e-12);
  close(cost.crawfordCurve({ T1: 100, n: 8, learningRate: 0.2 }).unit / c.unit, 0.8, 1e-13);
  const e = cost.learningExponent({ learningRate: 0.2 });
  close(e.b, Math.log2(0.8), 1e-15); close(e.progressRatio, 0.8);
});

test('experience curve on cumulative volume and serial-production saving', () => {
  close(cost.costAtCumulativeVolume({ c0: 100, q0: 10, q: 40, learningRate: 0.15 }), 100 * 0.85 ** 2, 1e-12);
  const s = cost.serialProductionSavings({ firstUnitCost: 100, n: 4, progressRatio: 0.8 });
  close(s.withoutLearning, 400); close(s.saving, 400 - s.total, 1e-12); assert.ok(s.savingFraction > 0.2 && s.savingFraction < 0.22);
  close(cost.serialProductionSavings({ firstUnitCost: 100, n: 4, progressRatio: 0.8, model: 'wright' }).total, 256, 1e-12);
});

test('AACE 18R-97 classes: accuracy bands widen from class 1 to class 5', () => {
  assert.deepEqual(cost.AACE_CLASSES[5].low, [-0.20, -0.50]); assert.deepEqual(cost.AACE_CLASSES[5].high, [0.30, 1.00]);
  assert.deepEqual(cost.AACE_CLASSES[1].low, [-0.03, -0.10]); assert.deepEqual(cost.AACE_CLASSES[1].high, [0.03, 0.15]);
  assert.deepEqual(cost.AACE_CLASSES[3].definition, [0.10, 0.40]);
  for (let c = 1; c < 5; c++) assert.ok(cost.AACE_CLASSES[c + 1].high[1] > cost.AACE_CLASSES[c].high[1]);
  const r = cost.estimateRange(100, 4);
  close(r.wide.low, 70); close(r.wide.high, 150); close(r.narrow.low, 85);
});

test('probabilistic contingency: independent normals add in variance; P80 = mean + 0.8416 sigma; correlation widens it', () => {
  const items = [{ name: 'a', dist: { dist: 'normal', mean: 100, std: 10 } }, { name: 'b', dist: { dist: 'normal', mean: 200, std: 20 } }];
  const r = cost.probabilisticContingency({ items, n: 4000, seed: 3 });
  const sd = Math.sqrt(500);
  close(r.mean, 300, 1e-3); close(r.std, sd, 0.05); close(r.p80, 300 + 0.8416 * sd, 5e-3); close(r.p50, 300, 3e-3);
  close(r.contingencyP80, r.p80 - 300, 1e-12);
  const c = cost.probabilisticContingency({ items, correlations: [{ a: 'a', b: 'b', rho: 0.9 }], n: 4000, seed: 3 });
  assert.ok(c.std > 1.2 * r.std && c.std < 30.5);          // perfectly correlated limit is 30
  close(c.mean, 300, 1e-3);
  const again = cost.probabilisticContingency({ items, n: 4000, seed: 3 });
  assert.equal(again.p80, r.p80);
  assert.match(r.convention, /not exceeded/);
});

test('composite escalation, currency exposure, CBS roll-up and benchmarks', () => {
  const e = cost.compositeEscalation({ weights: { steel: 60, labour: 40 }, indices: { steel: { from: 100, to: 110 }, labour: { from: 200, to: 240 } } });
  close(e.factor, 0.6 * 1.1 + 0.4 * 1.2, 1e-14);
  const x = cost.currencyExposure({ items: [{ name: 'a', amount: 110, currency: 'USD' }, { name: 'b', amount: 100, currency: 'EUR' }], rates: { EUR: 1, USD: 1.1 }, projectCurrency: 'EUR' });
  close(x.total, 200); close(x.byCurrency.USD.share, 0.5); close(x.byCurrency.USD.native, 110);
  const cbs = cost.cbsRollup([{ code: '1.1', cost: 10, discipline: 'structural' }, { code: '1.2', cost: 20, discipline: 'marine' }, { code: '2.1.1', cost: 5, discipline: 'marine' }]);
  assert.equal(cbs.nodes['1'], 30); assert.equal(cbs.nodes['2'], 5); assert.equal(cbs.nodes['2.1'], 5); assert.equal(cbs.total, 35);
  assert.equal(cbs.byDiscipline.marine, 25);
  assert.throws(() => cost.cbsRollup([{ code: '1', cost: 1 }, { code: '1', cost: 2 }]), /duplicate/);
  assert.throws(() => cost.cbsRollup([{ code: '1', cost: NaN }]), RangeError);
  assert.deepEqual(cost.benchmarkMetrics({ cost: 3e9, MW: 1000, km: 100 }), { perMW: 3e6, perTonne: null, perKm: 3e7 });
});
