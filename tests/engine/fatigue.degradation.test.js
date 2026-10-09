import test from 'node:test';
import assert from 'node:assert/strict';
import {
  corrosionWastage, remainingWall, coatingBreakdown, coatingBreakdownConstants, anodeResistance, anodeDesign,
  marineGrowth, archardWear, chainInterlinkWear, weibullDeterioration, markovDeterioration,
} from '../../src/engine/fatigue/degradation.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('corrosion wastage models: hand values', () => {
  // Linear: 0.1 mm/yr after a 10-year coating life → 1.5 mm at 25 years.
  close(corrosionWastage({ time: 25, coatingLife: 10, rate: 0.1e-3 }).depth, 1.5e-3, 1e-12);
  assert.equal(corrosionWastage({ time: 5, coatingLife: 10, rate: 0.1e-3 }).depth, 0);
  // Guedes Soares–Garbatov: d∞ = 2 mm, τc = 5, τt = 10: at t = 15, d = 2·(1 − e^-1) = 1.26424 mm
  const g = corrosionWastage({ time: 15, model: 'guedes-soares-garbatov', coatingLife: 5, dInf: 2e-3, transitionTime: 10 });
  close(g.depth, 2e-3 * (1 - Math.exp(-1)), 1e-12);
  close(g.rate, (2e-3 / 10) * Math.exp(-1), 1e-12);
  assert.ok(corrosionWastage({ time: 1e4, model: 'guedes-soares-garbatov', coatingLife: 5, dInf: 2e-3, transitionTime: 10 }).depth <= 2e-3);
  // Paik: C1 = 0.17 mm, C2 = 0.67, 16 years of exposure: 0.17·16^0.67 = 1.0893 mm
  close(corrosionWastage({ time: 21, model: 'paik', coatingLife: 5, c1: 0.17e-3, c2: 0.67 }).depth, 0.17e-3 * 16 ** 0.67, 1e-12);
});

test('remaining wall and remaining life against a minimum thickness', () => {
  // 20 mm plate, 15 mm minimum, both sides at 0.1 mm/yr after 5 years of coating: 5 mm lost after 25 years of
  // exposure → end of life at age 30; at age 12 the loss is 1.4 mm.
  const r = remainingWall({ thickness: 0.02, minThickness: 0.015, age: 12, surfaces: 2, corrosion: { coatingLife: 5, rate: 0.1e-3 } });
  close(r.wastage, 1.4e-3, 1e-12);
  close(r.remainingThickness, 18.6e-3, 1e-12);
  close(r.endOfLifeAge, 30, 1e-9);
  close(r.remainingLife, 18, 1e-9);
  // The saturating model never reaches the limit when d∞ is below the allowance.
  const s = remainingWall({ thickness: 0.02, minThickness: 0.015, age: 12, corrosion: { model: 'guedes-soares-garbatov', coatingLife: 5, dInf: 3e-3, transitionTime: 8 } });
  assert.equal(s.endOfLifeAge, Infinity);
});

test('DNV-RP-B401 coating breakdown factors', () => {
  // Category III below 30 m, 25 years: mean 0.02 + 0.008·12.5 = 0.12, final 0.02 + 0.008·25 = 0.22
  const f = coatingBreakdown({ ...coatingBreakdownConstants('III', 100), designLife: 25 });
  close(f.initial, 0.02, 1e-12);
  close(f.mean, 0.12, 1e-12);
  close(f.final, 0.22, 1e-12);
  // Category I in shallow water for 25 years exceeds 1 at the end: final capped, mean = 1 − (1 − a)²/(2·b·tf)
  const c = coatingBreakdown({ ...coatingBreakdownConstants('I', 10), designLife: 25 });
  assert.equal(c.final, 1);
  close(c.mean, 1 - 0.81 / 5, 1e-12);
  // The capped mean equals the time average of min(1, a + b·t).
  let avg = 0;
  for (let i = 0; i < 100000; i++) avg += Math.min(1, 0.1 + 0.1 * ((i + 0.5) * 25) / 100000) / 100000;
  close(c.mean, avg, 1e-6);
});

test('DNV-RP-B401 anode resistance formulas: hand values', () => {
  // Long slender stand-off, ρ = 0.3 Ω·m, L = 2.5 m, r = 0.1 m: 0.3/(2π·2.5)·(ln 100 − 1) = 0.068854 Ω
  close(anodeResistance({ type: 'long-slender-standoff', resistivity: 0.3, length: 2.5, radius: 0.1 }), 0.068854, 1e-5);
  // Short slender: L = 0.3, r = 0.1: x = r/2L = 1/6
  const x = 1 / 6;
  const root = Math.sqrt(1 + x * x);
  close(anodeResistance({ type: 'short-slender-standoff', resistivity: 0.3, length: 0.3, radius: 0.1 }), (0.3 / (2 * Math.PI * 0.3)) * (Math.log(6 * (1 + root)) + x - root), 1e-12);
  // Long flush: ρ/(2S), S = (1.0 + 0.2)/2 → 0.25 Ω; bracelet: 0.315·ρ/√A, A = 0.5 m² → 0.13364 Ω
  close(anodeResistance({ type: 'long-flush', resistivity: 0.3, length: 1.0, width: 0.2 }), 0.25, 1e-12);
  close(anodeResistance({ type: 'bracelet', resistivity: 0.3, area: 0.5 }), 0.133643, 1e-5);
  // At the changeover L = 4r the two slender formulas differ by 7 % (ln 16 − 1 = 1.773 against 1.894);
  // for a long anode the short-anode expression tends to the long one.
  const a = anodeResistance({ type: 'long-slender-standoff', resistivity: 0.3, length: 0.4, radius: 0.1 });
  const b = anodeResistance({ type: 'short-slender-standoff', resistivity: 0.3, length: 0.4, radius: 0.1 });
  close(a, b, 0.08);
  close(a, (0.3 / (2 * Math.PI * 0.4)) * (Math.log(16) - 1), 1e-12);
  close(anodeResistance({ type: 'short-slender-standoff', resistivity: 0.3, length: 50, radius: 0.1 }),
    (0.3 / (2 * Math.PI * 50)) * (Math.log(2000) - 1), 1e-3);
});

test('DNV-RP-B401 anode design: current demand, mass M = Icm·tf·8760/(u·ε) and current checks, by hand', () => {
  // Bare steel 500 m², densities 0.15 / 0.07 / 0.10 A/m²: Ici = 75 A, Icm = 35 A, Icf = 50 A.
  // M = 35·25·8760/(0.9·2000) = 4258.3 kg → with 200 kg anodes: 22 anodes for mass.
  const anode = { type: 'long-slender-standoff', netMass: 200, length: 2.5, radius: 0.1, coreRadius: 0.03 };
  const r = anodeDesign({ zones: [{ area: 500, initialDensity: 0.15, meanDensity: 0.07, finalDensity: 0.10 }], designLife: 25, anode, resistivity: 0.3 });
  close(r.currentInitial, 75, 1e-12);
  close(r.currentMean, 35, 1e-12);
  close(r.currentFinal, 50, 1e-12);
  close(r.requiredMass, (35 * 25 * 8760) / 1800, 1e-12);
  assert.equal(r.anodesForMass, 22);
  // Initial output per anode: ΔE/Ra = 0.25/0.068854 = 3.6309 A → 75/3.6309 = 20.66 → 21 anodes.
  close(r.anodeCurrentInitial, 0.25 / 0.068854, 1e-4);
  assert.equal(r.anodesForInitialCurrent, 21);
  // Final: length 2.25 m, remaining volume 20 kg/2750 = 7.2727e-3 m³ → r = √(V/(πL) + 0.03²) = 0.043936 m,
  // Ra = 0.3/(2π·2.25)·(ln(4·2.25/0.043936) − 1) = 0.091722 Ω → 2.7256 A → 50/2.7256 = 18.3 → 19 anodes.
  const rf = Math.sqrt(20 / 2750 / (Math.PI * 2.25) + 0.0009);
  close(r.resistanceFinal, (0.3 / (2 * Math.PI * 2.25)) * (Math.log(9 / rf) - 1), 1e-12);
  assert.equal(r.anodesForFinalCurrent, 19);
  assert.equal(r.anodes, 22);
  assert.equal(r.governing, 'mass');
  close(r.totalNetMass, 4400, 1e-12);
  // Coated zone: Category III below 30 m scales the mean demand by 0.12.
  const c = anodeDesign({ zones: [{ area: 500, initialDensity: 0.15, meanDensity: 0.07, finalDensity: 0.10, ...coatingBreakdownConstants('III', 100) }], designLife: 25, anode, resistivity: 0.3 });
  close(c.currentMean, 35 * 0.12, 1e-12);
  close(c.currentFinal, 50 * 0.22, 1e-12);
});

test('marine growth default profile (NORSOK N-003)', () => {
  assert.deepEqual(marineGrowth({ elevation: -10, latitudeDeg: 57 }), { thickness: 0.1, density: 1325 });
  assert.equal(marineGrowth({ elevation: -60, latitudeDeg: 57 }).thickness, 0.05);
  assert.equal(marineGrowth({ elevation: -10, latitudeDeg: 65 }).thickness, 0.06);
  assert.equal(marineGrowth({ elevation: -60, latitudeDeg: 65 }).thickness, 0.03);
  assert.equal(marineGrowth({ elevation: 5, latitudeDeg: 65 }).thickness, 0);
  assert.equal(marineGrowth({ elevation: -20, profile: [{ zTop: 0, zBottom: -30, thickness: 0.15 }] }).thickness, 0.15);
  assert.throws(() => marineGrowth({ elevation: -10, latitudeDeg: 30 }), /latitudeDeg/);
});

test('Archard wear: V = k·F·s/H and chain interlink application', () => {
  // k = 1e-4, F = 1 kN, s = 1 km, H = 2 GPa → V = 1e-4·1e3·1e3/2e9 = 5e-8 m³
  const w = archardWear({ force: 1e3, slidingDistance: 1e3, k: 1e-4, hardness: 2e9, contactArea: 1e-4 });
  close(w.volume, 5e-8, 1e-12);
  close(w.depth, 5e-4, 1e-12);
  // Chain: d = 0.1 m, Δθ = 0.01 rad, 1e6 cycles → s = 2·0.01·0.05·1e6 = 1000 m; K = 1e-15 1/Pa, T = 1 MN →
  // V = 1e-15·1e6·1e3 = 1e-6 m³ over 0.25·d² = 2.5e-3 m² → 0.4 mm
  const c = chainInterlinkWear({ tension: 1e6, diameter: 0.1, rotationRange: 0.01, cycles: 1e6, wearCoefficient: 1e-15 });
  close(c.slidingDistance, 1000, 1e-12);
  close(c.depth, 4e-4, 1e-12);
  close(c.remainingDiameter, 0.0996, 1e-12);
});

test('Weibull deterioration and condition-state Markov chain', () => {
  const d = weibullDeterioration({ time: 10, eta: 20, beta: 2, threshold: 0.5 });
  close(d.condition, Math.exp(-0.25), 1e-12);
  close(d.timeToThreshold, 20 * Math.sqrt(Math.log(2)), 1e-12);
  close(d.remainingLife, 20 * Math.sqrt(Math.log(2)) - 10, 1e-12);
  // Two-state chain with annual failure probability p: P[worst] after n years = 1 − (1 − p)^n, mean time 1/p.
  const m = markovDeterioration({ P: [[0.9, 0.1], [0, 1]], initial: [1, 0], years: 10, failureProbabilityLimit: 0.5 });
  close(m.pWorst[10], 1 - 0.9 ** 10, 1e-12);
  close(m.meanTimeToWorst[0], 10, 1e-9);
  assert.equal(m.yearLimitExceeded, 7); // 1 − 0.9^7 = 0.522 is the first value above 0.5
  // Three-state progressive chain: mean time = 1/0.2 + 1/0.1 = 15 years.
  const c = markovDeterioration({ P: [[0.8, 0.2, 0], [0, 0.9, 0.1], [0, 0, 1]], initial: [1, 0, 0], years: 200 });
  close(c.meanTimeToWorst[0], 15, 1e-9);
  close(c.meanTimeToWorst[1], 10, 1e-9);
  close(c.distribution[200].reduce((a, b) => a + b, 0), 1, 1e-12);
  close(c.expectedState[0], 1, 1e-12);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => corrosionWastage({ time: -1, rate: 1 }), /time/);
  assert.throws(() => anodeResistance({ type: 'long-slender-standoff', resistivity: 0.3, length: 0.2, radius: 0.1 }), /length/);
  assert.throws(() => markovDeterioration({ P: [[0.5, 0.4], [0, 1]], initial: [1, 0], years: 1 }), /row 0/);
});
