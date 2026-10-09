import test from 'node:test';
import assert from 'node:assert/strict';
import {
  GENERIC_HAMMERS, impactVelocity, waveEquation, rodImpactForce, clayPileCapacityFactor, srdProfile, bearingGraph,
  drivabilityAnalysis, selfWeightPenetration, pileSubmergedWeight,
} from '../../src/engine/geotech/driving.js';
import { axialCapacity } from '../../src/engine/geotech/axialcapacity.js';
import { createProfile, exampleProfile } from '../../src/engine/geotech/soilprofile.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-12), `${a} vs ${b}`);

test('wave equation benchmark: rod struck by a rigid mass, sigma = rho c v0 exp(-rho c A t / M)', () => {
  const area = 0.3, M = 30e3, v0 = 5, E = 210e9, rho = 7850, c = Math.sqrt(E / rho);
  // 200 m rod: the reflection returns to the head after 2L/c = 77 ms; compare over the first 50 ms
  const run = (nSegments, dtFactor) => waveEquation({ hammer: { ramMass: M, impactVelocity: v0 }, pile: { length: 200, area, nSegments }, dtFactor, duration: 0.06, history: true });
  const r = run(400, 1);
  close(r.waveSpeed, c, 1e-12); close(r.impedance, rho * c * area, 1e-12);
  const Z = r.impedance;
  let maxErr = 0;
  r.history.t.forEach((t, i) => {
    if (t < 5e-4 || t > 0.05) return;
    const ref = rodImpactForce({ ramMass: M, impactVelocity: v0, area, t });
    close(ref.force, Z * v0 * Math.exp((-Z * t) / M), 1e-12);
    maxErr = Math.max(maxErr, Math.abs(r.history.headForce[i] - ref.force) / (Z * v0));
  });
  assert.ok(maxErr < 2e-3, `max head-force error ${maxErr}`);
  // peak stress rho c v0 (the lumped head mass removes m_head/(M + m_head) = 1.9 % of it)
  close(rodImpactForce({ ramMass: M, impactVelocity: v0, area, t: 0 }).peakStress, rho * c * v0, 1e-12);
  close(r.maxCompressiveStress, rho * c * v0, 0.025);
  // before the front reaches the toe (L/c = 39 ms) there is no tension anywhere in the rod
  const early = waveEquation({ hammer: { ramMass: M, impactVelocity: v0 }, pile: { length: 200, area, nSegments: 400 }, dtFactor: 1, duration: 0.03 });
  assert.ok(early.maxTensileStress < 1e-6 * early.maxCompressiveStress);
  assert.ok(r.maxTensileStress > 0); // free-end reflection after 39 ms
  // decay time constant: force falls to 1/e after t = M/Z
  const tau = M / Z, iTau = r.history.t.findIndex((t) => t >= tau), h = r.history;
  const w = (tau - h.t[iTau - 1]) / (h.t[iTau] - h.t[iTau - 1]);
  close(h.headForce[iTau - 1] + w * (h.headForce[iTau] - h.headForce[iTau - 1]), (Z * v0) / Math.E, 5e-3);
  // a reduced time step keeps the head force accurate (dispersion only affects the wave front)
  const r9 = run(400, 0.9);
  let err9 = 0;
  r9.history.t.forEach((t, i) => { if (t > 1e-3 && t < 0.05) err9 = Math.max(err9, Math.abs(r9.history.headForce[i] - rodImpactForce({ ramMass: M, impactVelocity: v0, area, t }).force) / (Z * v0)); });
  assert.ok(err9 < 0.01, `max head-force error at dtFactor 0.9: ${err9}`);
});

test('energy conservation without soil: kinetic + strain energy constant; only the impact merge loss m_head/(M + m_head)', () => {
  const area = 0.3, M = 30e3, v0 = 5, rho = 7850;
  const loss = (n) => { const mHead = 0.5 * rho * area * (200 / n); return mHead / (M + mHead); };
  for (const [n, dtFactor] of [[200, 1], [400, 1], [400, 0.7]]) {
    // long enough for several reflections and ram rebound
    const r = waveEquation({ hammer: { ramMass: M, impactVelocity: v0 }, pile: { length: 200, area, nSegments: n }, dtFactor, duration: 0.07 });
    assert.ok(Math.abs(r.energyDrift - loss(n)) < 2e-4, `n=${n}: drift ${r.energyDrift} vs merge loss ${loss(n)}`);
    close(r.impactEnergy, 0.5 * M * v0 * v0, 1e-12);
    assert.ok(r.enthru <= r.impactEnergy && r.enthru > 0.9 * r.impactEnergy);
  }
  // with a cushion spring there is no merge: energy is conserved to integration accuracy
  const c = waveEquation({ hammer: { ramMass: M, impactVelocity: v0 }, cushionStiffness: 3e9, helmetMass: 2e3, pile: { length: 200, area, nSegments: 400 }, duration: 0.07 });
  assert.ok(c.energyDrift < 2e-3, `drift ${c.energyDrift}`);
  // the cushion lengthens the pulse and lowers the peak stress below rho c v0
  assert.ok(c.maxCompressiveStress < 0.9 * rho * Math.sqrt(210e9 / rho) * v0);
});

test('hammer library is labelled generic; impact velocity from delivered energy', () => {
  for (const h of Object.values(GENERIC_HAMMERS)) {
    assert.match(h.label, /^Generic/);
    close(h.ratedEnergy / (h.ramMass * 9.80665), 2.0394, 1e-3); // equivalent stroke
  }
  close(impactVelocity({ ratedEnergy: 500e3, ramMass: 25e3, efficiency: 0.9 }), 6, 1e-12); // sqrt(2·0.9·500e3/25e3)
});

test('bearing graph: blow count rises with SRD, refusal when the hammer cannot move the toe; static sanity', () => {
  const o = { hammer: GENERIC_HAMMERS['generic-1200'], pile: { length: 60, D: 2.5, t: 0.05 }, penetration: 40 };
  const g = bearingGraph({ ...o, srd: [5e6, 10e6, 20e6, 40e6, 80e6, 2000e6] });
  for (let i = 1; i < 5; i++) assert.ok(g.blowsPerMetre[i] > g.blowsPerMetre[i - 1]);
  assert.equal(g.blowsPerMetre[5], Infinity);
  assert.equal(g.set[5], 0);
  for (let i = 0; i < 5; i++) close(g.blowsPer025m[i], g.blowsPerMetre[i] / 4);
  // energy bound: work done against the static resistance over the set cannot exceed the transferred energy
  for (let i = 0; i < 5; i++) assert.ok(g.srd[i] * g.set[i] < g.enthru[i], `SRD·set ${g.srd[i] * g.set[i]} vs ENTHRU ${g.enthru[i]}`);
  // mesh and time-step independence of the blow count (within 5 %)
  const a = bearingGraph({ ...o, srd: [40e6] }).blowsPerMetre[0];
  const b = bearingGraph({ ...o, pile: { ...o.pile, nSegments: 200 }, dtFactor: 0.5, srd: [40e6] }).blowsPerMetre[0];
  close(a, b, 0.05);
  // no soil: stress at the head is rho c v0 for a rigid ram on a prismatic pile
  const v0 = impactVelocity(o.hammer);
  const free = waveEquation({ hammer: o.hammer, pile: o.pile, dtFactor: 1 });
  close(free.maxCompressiveStress, 7850 * Math.sqrt(210e9 / 7850) * v0, 0.04);
  // free toe: the compression wave reflects in tension of almost the same magnitude
  assert.ok(free.maxTensileStress > 0.7 * free.maxCompressiveStress);
});

test('SRD (Stevens et al. 1982 cases): pile capacity factor, hand-checked components, ordering of the bounds', () => {
  close(clayPileCapacityFactor(1), 0.5); close(clayPileCapacityFactor(4), 0.5 * 4 ** 0.3);
  const D = 2, t = 0.05, Di = 1.9, L = 30;
  const clay = createProfile([{ type: 'clay', top: 0, thickness: 40, gammaSub: 8000, suTop: 60e3, suGrad: 0, OCR: 4 }]);
  const s = srdProfile({ profile: clay, D, t, Lmax: L, dz: 0.25 });
  const i = s.z.length - 1;
  // shaft: Fp × static API outer friction; trapezoidal versus Gauss integration
  const Qfo = axialCapacity({ profile: clay, L, D, t }).Qfo, Fp = 0.5 * 4 ** 0.3;
  close(s.cases.pluggedLower.shaft[i], Fp * Qfo, 5e-3);
  close(s.cases.coringUpper.shaft[i], s.cases.pluggedLower.shaft[i] * (1 + Di / D), 1e-12);
  close(s.cases.coringLower.shaft[i], s.cases.pluggedLower.shaft[i] * (1 + (0.5 * Di) / D), 1e-12);
  // toe: 9 su on the annulus (coring) or the gross area (plugged)
  close(s.cases.coringLower.toe[i], (9 * 60e3 * Math.PI * (D * D - Di * Di)) / 4, 1e-12);
  close(s.cases.pluggedLower.toe[i], (9 * 60e3 * Math.PI * D * D) / 4, 1e-12);
  close(s.cases.pluggedUpper.total[i], s.cases.pluggedLower.total[i], 1e-12); // clay upper-bound factors default to 1
  assert.ok(s.cases.coringLower.total[i] < s.cases.coringUpper.total[i]);
  // sand: upper-bound plugged multipliers 1.3 (friction) and 1.5 (end bearing) by default, both overridable
  const sand = exampleProfile('denseSand');
  const ss = srdProfile({ profile: sand, D, t, Lmax: L });
  const j = ss.z.length - 1;
  close(ss.cases.pluggedUpper.shaft[j], 1.3 * ss.cases.pluggedLower.shaft[j], 1e-12);
  close(ss.cases.pluggedUpper.toe[j], 1.5 * ss.cases.pluggedLower.toe[j], 1e-12);
  const s2 = srdProfile({ profile: sand, D, t, Lmax: L, upperPlugged: { sandFriction: 1, sandBearing: 2 } });
  close(s2.cases.pluggedUpper.toe[j], 2 * ss.cases.pluggedLower.toe[j], 1e-12);
});

test('drivability with depth: SRD and blow counts grow, totals integrate, refusal limit is applied', () => {
  const profile = exampleProfile('layeredSandClay');
  const o = { profile, hammer: GENERIC_HAMMERS['generic-1200'], pile: { length: 60, D: 2.5, t: 0.05, nSegments: 60 }, depths: [5, 10, 20, 30, 40, 50] };
  const d = drivabilityAnalysis(o);
  for (let i = 1; i < 6; i++) { assert.ok(d.srd[i] > d.srd[i - 1]); assert.ok(d.blowsPerMetre[i] > d.blowsPerMetre[i - 1]); }
  assert.equal(d.refusalDepth, null);
  // trapezoidal integration of blows per metre, written out
  const bpm = d.blowsPerMetre, z = o.depths;
  let total = bpm[0] * z[0];
  for (let i = 1; i < 6; i++) total += 0.5 * (bpm[i] + bpm[i - 1]) * (z[i] - z[i - 1]);
  close(d.totalBlows, total, 1e-12);
  close(d.fatigueBlocks.reduce((s, b) => s + b.blows, 0), total, 1e-12);
  assert.equal(d.fatigueBlocks.length, 6);
  close(d.fatigueBlocks[3].stressRange, d.maxCompressiveStress[3] + d.maxTensileStress[3], 1e-12);
  // a small hammer meets refusal: the same pile with a strict limit
  const small = drivabilityAnalysis({ ...o, hammer: GENERIC_HAMMERS['generic-500'], refusalBlowsPer025m: 20 });
  assert.ok(small.refusalDepth !== null);
  assert.ok(small.refusal[5]);
  assert.ok(small.blowsPer025m[small.depth.indexOf(small.refusalDepth)] > 20);
  assert.ok(small.totalBlows > 0);
  // plugged upper bound is harder to drive than coring lower bound
  const hard = drivabilityAnalysis({ ...o, srdCase: 'pluggedUpper', depths: [40] }), easy = drivabilityAnalysis({ ...o, srdCase: 'coringLower', depths: [40] });
  assert.ok(hard.srd[0] > easy.srd[0] && hard.blowsPerMetre[0] > easy.blowsPerMetre[0], `${hard.srd[0]} ${easy.srd[0]} ${hard.blowsPerMetre[0]} ${easy.blowsPerMetre[0]}`);
  for (const v of [...hard.maxCompressiveStress, ...d.maxCompressiveStress]) assert.ok(v < 600e6); // bounded (stable)
});

test('self-weight penetration and pile-run screening by hand', () => {
  // SRD rising 1 MN/m to 10 m, a crust, then a weak layer: 0→10 MN at 10 m, drops to 4 MN at 12 m, recovers by 20 m
  const z = [0, 10, 12, 16, 20, 30], srd = [0, 10e6, 4e6, 5e6, 12e6, 30e6];
  const r = selfWeightPenetration({ z, srd, pileWeight: 3e6, hammerWeight: 3e6 });
  close(r.swpPile, 3, 1e-12); close(r.swpPileHammer, 6, 1e-12);
  assert.equal(r.runRisk, true);
  assert.deepEqual(r.runIntervals, [{ from: 12, to: 16, minSrd: 4e6 }]);
  assert.equal(selfWeightPenetration({ z, srd, pileWeight: 1e6, hammerWeight: 1e6 }).runRisk, false);
  assert.equal(selfWeightPenetration({ z, srd, pileWeight: 50e6 }).swpPile, null);
  // submerged weight: (7850 − 1025) g A L
  const A = (Math.PI * (2.5 ** 2 - 2.4 ** 2)) / 4;
  close(pileSubmergedWeight({ D: 2.5, t: 0.05, length: 60 }), 6825 * 9.80665 * A * 60, 1e-12);
});
