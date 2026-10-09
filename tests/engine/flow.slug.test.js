import test from 'node:test';
import assert from 'node:assert/strict';
import {
  R_GAS, slugTranslationalVelocity, slugBodyHoldupGregory, slugBodyHoldupBarneaBrauner, slugFrequencyGregoryScott,
  slugFrequencyZabaras, meanSlugLength, slugLengthDistribution, slugUnitCell, slugTrain, kelvinHelmholtzStability,
  slugGrowthCriterion, boeCriterion, potsNumber, taitelStability, severeSluggingCycle, severeSluggingMitigation,
  surgeVolume, pigSlugVolume,
} from '../../src/engine/flow/slug.js';
import { G } from '../../src/engine/flow/driftflux.js';
import { stratifiedGeometry } from '../../src/engine/flow/flowregime.js';
import { stratifiedTwoFluid } from '../../src/engine/flow/multiphase.js';
import { makeRng } from '../../src/engine/math/rng.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg ?? ''} got ${a}, expected ${b}`);
const oilGas = { D: 0.1, rhoL: 850, rhoG: 20, muL: 2e-3, muG: 1.5e-5 };

test('slug translational velocity: Nicklin vertical 1.2 vm + 0.35 √(gD), Bendiksen horizontal 1.05 vm + 0.54 √(gD)', () => {
  const v = slugTranslationalVelocity({ vm: 2, D: 0.1, theta: Math.PI / 2 });
  close(v.vt, 2.4 + 0.35 * Math.sqrt(G * 0.1), 1e-13);
  const h = slugTranslationalVelocity({ vm: 2, D: 0.1 });
  close(h.vt, 2.1 + 0.54 * Math.sqrt(G * 0.1), 1e-13);
  close(slugTranslationalVelocity({ vm: 2, D: 0.1, C0: 1.3, vd: 0.1 }).vt, 2.7, 1e-13);
});

test('Gregory, Nicholson & Aziz (1978) slug-body holdup by hand', () => {
  close(slugBodyHoldupGregory(8.66), 0.5, 1e-14);
  // vm = 3 m/s: (3/8.66)^1.39 = exp(1.39 × ln 0.34642) = exp(−1.47355) = 0.22911 → 1/1.22911 = 0.81360
  close(slugBodyHoldupGregory(3), 0.81360, 1e-4);
  assert.equal(slugBodyHoldupGregory(0), 1);
});

test('Barnea & Brauner (1985) slug-body holdup: unaerated at low velocity, bounded by 0.48 at high velocity', () => {
  const p = { ...oilGas, sigma: 0.025 };
  assert.equal(slugBodyHoldupBarneaBrauner({ ...p, vm: 0.2 }), 1);
  const mid = slugBodyHoldupBarneaBrauner({ ...p, vm: 3 });
  assert.ok(mid < 1 && mid > 0.48);
  assert.ok(slugBodyHoldupBarneaBrauner({ ...p, vm: 6 }) <= mid);
  close(slugBodyHoldupBarneaBrauner({ ...p, vm: 50 }), 0.48, 1e-12);
});

test('Gregory & Scott (1969) slug frequency by hand: vsl = 1, vm = 3 m/s, D = 0.1 m', () => {
  // vsl/(gD) = 1/0.980665 = 1.019716; 19.75/3 + 3 = 9.58333; product 9.77228; ^1.2 = exp(1.2 × 2.27955) = 15.417; × 0.0226 = 0.3484 Hz
  close(slugFrequencyGregoryScott({ vsl: 1, vm: 3, D: 0.1 }), 0.3484, 5e-4);
  // minimum of the frequency with respect to vm at vm = √19.75 = 4.44 m/s for fixed vsl
  const f = (vm) => slugFrequencyGregoryScott({ vsl: 1, vm, D: 0.1 });
  assert.ok(f(Math.sqrt(19.75)) < f(4.0) && f(Math.sqrt(19.75)) < f(5.0));
});

test('Zabaras (2000) frequency: SI form equals the published US-unit form; horizontal factor 0.836', () => {
  const vsl = 1.2; const vm = 3.5; const D = 0.15; const theta = 5 * Math.PI / 180;
  const ft = 1 / 0.3048;
  const us = 0.0226 * (vsl * ft / (32.174 * D * ft)) ** 1.2 * (212.6 / (vm * ft) + vm * ft) ** 1.2 * (0.836 + 2.75 * Math.sin(theta) ** 0.25);
  close(slugFrequencyZabaras({ vsl, vm, D, theta }), us, 1e-3);
  close(slugFrequencyZabaras({ vsl, vm, D }) / slugFrequencyGregoryScott({ vsl, vm, D }), 0.836, 1e-12);
  assert.ok(slugFrequencyZabaras({ vsl, vm, D, theta }) > slugFrequencyZabaras({ vsl, vm, D }));
  assert.throws(() => slugFrequencyZabaras({ vsl, vm, D, theta: -0.1 }), /theta/);
});

test('mean slug length: 32 D, and Scott, Shoham & Brill (1989) for a 24 in line by hand', () => {
  close(meanSlugLength({ D: 0.1 }), 3.2, 1e-14);
  close(meanSlugLength({ D: 0.1, factor: 18 }), 1.8, 1e-14);
  // ln 24 = 3.17805; ^0.1 = 1.122577; ln Ls = −25.4 + 28.5 × 1.122577 = 6.5934 → Ls = 730.3 ft = 222.6 m
  close(meanSlugLength({ D: 24 * 0.0254, method: 'scott' }), 222.6, 2e-3);
  assert.throws(() => meanSlugLength({ D: 0.02, method: 'scott' }), RangeError);
});

test('lognormal slug-length distribution: median, 1/1000 maximum and sample mean', () => {
  const d = slugLengthDistribution({ meanLength: 10, sigmaLn: 0.5 });
  close(d.median, 10 * Math.exp(-0.125), 1e-13);
  // Φ⁻¹(0.999) = 3.0902
  close(d.maxLength, d.median * Math.exp(0.5 * 3.0902), 1e-4);
  close(d.quantile(0.5), d.median, 1e-9);
  const rng = makeRng(7);
  let s = 0;
  const n = 40000;
  for (let i = 0; i < n; i++) s += d.sample(rng);
  close(s / n, 10, 0.02);
});

test('slug unit cell: phase mass balances, unit-holdup identity and film momentum balance', () => {
  for (const [vsl, vsg, theta] of [[1, 2, 0], [0.5, 3, 0.03], [1.5, 1, -0.01]]) {
    const u = slugUnitCell({ ...oilGas, vsl, vsg, theta });
    assert.equal(u.feasible, true);
    const b = u.slugFraction;
    // liquid and gas volumetric balances over the unit
    close(u.vm * u.slugHoldup * b + u.vLF * u.filmHoldup * (1 - b), vsl, 1e-10);
    close(u.vm * (1 - u.slugHoldup) * b + u.vGF * (1 - u.filmHoldup) * (1 - b), vsg, 1e-10);
    // Taitel & Barnea (1990): HLU = [vt HLS + vGS (1 − HLS) − vsg]/vt, independent of the zone lengths
    close(u.holdup, (u.vt * u.slugHoldup + u.vm * (1 - u.slugHoldup) - vsg) / u.vt, 1e-10);
    assert.ok(Math.abs(u.residual) < 1e-6 * Math.max(1, u.gradFrictionSlug));
    assert.ok(u.filmHoldup < u.slugHoldup && u.vLF < u.vm && u.vGF < u.vt);
    close(u.frequency, u.vt / u.unitLength, 1e-13);
    close(u.unitLength, u.slugLength / b, 1e-13);
    close(u.gradFriction, b * u.gradFrictionSlug + (1 - b) * u.gradFrictionFilm, 1e-13);
    close(u.unitPressureDrop, u.gradTotal * u.unitLength, 1e-13);
    assert.equal(stratifiedGeometry(u.filmLevel).holdup, u.filmHoldup);
  }
  close(slugUnitCell({ ...oilGas, vsl: 1, vsg: 2 }).slugLength, 3.2, 1e-13);
});

test('synthetic slug train: deterministic for a seed and consistent with the liquid and gas throughput', () => {
  const vsl = 1; const vsg = 2;
  const u = slugUnitCell({ ...oilGas, vsl, vsg });
  const args = { D: 0.1, vm: u.vm, vt: u.vt, slugHoldup: u.slugHoldup, filmHoldup: u.filmHoldup, vLF: u.vLF, vGF: u.vGF,
    slugFraction: u.slugFraction, meanSlugLength: u.slugLength, sigmaLn: 0.5, duration: 4000, nSamples: 4000, seed: 11 };
  const a = slugTrain(args);
  const b = slugTrain(args);
  assert.deepEqual(a.slugs.slice(0, 5), b.slugs.slice(0, 5));
  assert.notDeepEqual(a.slugs.slice(0, 5), slugTrain({ ...args, seed: 12 }).slugs.slice(0, 5));
  const A = Math.PI * 0.01 / 4;
  close(a.meanQL, vsl * A, 0.05);
  close(a.meanQG, vsg * A, 0.05);
  close(a.meanFrequency, u.frequency, 0.08);
  let s = 0;
  for (const q of a.qL) s += q;
  close(s / a.qL.length, a.meanQL, 1e-9);
  close(a.peakQL, A * u.vm * u.slugHoldup, 1e-13);
  assert.ok(a.maxSlugVolume > A * u.slugLength * u.slugHoldup);
});

test('Kelvin–Helmholtz stability of stratified flow: inviscid limit by hand, viscous criterion more restrictive', () => {
  const pts = [[0.02, 1], [0.05, 2], [0.2, 3], [0.3, 6]];
  const res = pts.map(([vsl, vsg]) => {
    const s = stratifiedTwoFluid({ ...oilGas, vsl, vsg });
    const k = kelvinHelmholtzStability({ ...oilGas, h: s.hL, vsl, vsg });
    const g = stratifiedGeometry(s.hL);
    // IKH: (vG − vL)² ≤ (αL/ρL + αG/ρG) (ρL − ρG) g A/(dAL/dhL)
    const lim = Math.sqrt((g.holdup / 850 + (1 - g.holdup) / 20) * 830 * G * (Math.PI * 0.1 / 4) / g.dAL);
    close(k.ikhLimit, lim, 1e-12);
    close(k.relativeVelocity, s.vG - s.vL, 1e-12);
    // the kinematic-wave term is bounded below by the inviscid one
    assert.ok(k.vkhRatio >= k.ikhRatio * (1 - 1e-9));
    return k;
  });
  assert.equal(res[0].vkhStable, true); assert.equal(res[0].ikhStable, true);
  assert.equal(res[1].vkhStable, true);
  assert.equal(res[2].vkhStable, false); assert.equal(res[2].ikhStable, true);
  assert.equal(res[3].ikhStable, false); assert.equal(res[3].vkhStable, false);
  // stable flow: the kinematic wave speed lies between the dynamic wave speeds
  assert.ok(res[0].cMinus <= res[0].kinematicVelocity && res[0].kinematicVelocity <= res[0].cPlus);
});

test('minimum-slip slug growth criterion by hand', () => {
  // HLS = 0.9, HLF = 0.3, vm = 3, vLF = 0.5: vF = (2.7 − 0.15)/0.6 = 4.25 m/s; vt = 1.2 × 3 + 0.2 = 3.8 m/s
  const g = slugGrowthCriterion({ vm: 3, D: 0.1, filmHoldup: 0.3, vLF: 0.5, slugHoldup: 0.9, C0: 1.2, vd: 0.2 });
  close(g.frontVelocity, 4.25, 1e-13); close(g.tailVelocity, 3.8, 1e-13); close(g.growthRate, 0.45, 1e-12);
  assert.equal(g.slugSurvives, true);
  // a thin film cannot feed the slug: vF = (2.7 − 0.025)/0.85 = 3.147 < 3.8
  assert.equal(slugGrowthCriterion({ vm: 3, D: 0.1, filmHoldup: 0.05, vLF: 0.5, slugHoldup: 0.9, C0: 1.2, vd: 0.2 }).slugSurvives, false);
});

test('Bøe (1981) criterion by hand and equivalence with the Pots ΠSS number', () => {
  // vsl,crit = P vsg / (ρL g α L) = 2e5 × 0.5 / (1000 × 9.80665 × 0.9 × 100) = 0.113302 m/s
  const b = boeCriterion({ vsl: 0.2, vsg: 0.5, P: 2e5, rhoL: 1000, alphaPipe: 0.9, L: 100 });
  close(b.vslCritical, 0.113302, 1e-5);
  assert.equal(b.severeSlugging, true);
  close(b.piSS, 0.113302 / 0.2, 1e-5);
  assert.equal(boeCriterion({ vsl: 0.05, vsg: 0.5, P: 2e5, rhoL: 1000, alphaPipe: 0.9, L: 100 }).severeSlugging, false);
  // same state expressed as mass rates of an ideal gas: ΠSS = (Z R T/M) WG/(g α L WL)
  const T = 293; const M = 0.029; const A = 0.01;
  const rhoG = 2e5 * M / (R_GAS * T);
  const p = potsNumber({ massRateG: rhoG * 0.5 * A, massRateL: 1000 * 0.2 * A, T, M, alphaPipe: 0.9, L: 100 });
  close(p.piSS, b.piSS, 1e-12);
  assert.equal(p.severeSlugging, true);
  // inclined riser: head builds more slowly
  close(boeCriterion({ vsl: 0.2, vsg: 0.5, P: 2e5, rhoL: 1000, alphaPipe: 0.9, L: 100, riserAngle: Math.PI / 6 }).vslCritical, 2 * 0.113302, 1e-5);
});

test('Taitel (1986) stability criterion by hand', () => {
  // (α/α') L − h = 0.9/0.89 × 100 − 10 = 91.1236 m of liquid → Ps,min = 1000 × 9.80665 × 91.1236 = 893 617 Pa
  const t = taitelStability({ Psep: 5e5, rhoL: 1000, alphaPipe: 0.9, L: 100, riserHeight: 10 });
  close(t.PsepMin, 893617, 1e-5);
  assert.equal(t.stable, false);
  close(t.rhs, 893617 / 101325, 1e-5);
  assert.equal(taitelStability({ Psep: 9e5, rhoL: 1000, alphaPipe: 0.9, L: 100, riserHeight: 10 }).stable, true);
  // a short flowline (α L/α' < h) is always stable
  assert.equal(taitelStability({ Psep: 1e5, rhoL: 1000, alphaPipe: 0.5, L: 10, riserHeight: 10 }).stable, true);
});

test('severe-slugging cycle estimate satisfies its liquid-volume and gas-mass balances', () => {
  const D = 0.05; const A = Math.PI * D * D / 4;
  const o = { vsl: 0.1, massRateG: 1.2 * 0.1 * A, gasRT: 287 * 293, D, L: 60, riserHeight: 10, Psep: 1e5, rhoL: 1000, alphaPipe: 0.83 };
  const c = severeSluggingCycle(o);
  assert.equal(c.severeSlugging, true);
  const P1 = 1e5 + 1000 * G * 10;
  close(c.peakPressure, P1, 1e-13);
  // liquid: inflow during build-up fills the riser and the penetration length
  close(o.vsl * c.buildUpTime, 10 + o.alphaPipe * c.penetrationLength, 1e-12);
  close(c.slugVolume, A * o.vsl * c.buildUpTime, 1e-12);
  // gas: P1 α (L − x) A = P0 α L A + ṁ RT t
  close(P1 * o.alphaPipe * (60 - c.penetrationLength) * A, 1e5 * o.alphaPipe * 60 * A + o.massRateG * o.gasRT * c.buildUpTime, 1e-12);
  // production stage: the tail travels x at vsg1/α
  const vsg1 = o.massRateG * o.gasRT / (P1 * A);
  close(c.productionTime, o.alphaPipe * c.penetrationLength / vsg1, 1e-12);
  close(c.productionLiquidRate, A * (o.vsl + vsg1), 1e-12);
  close(c.period, c.buildUpTime + c.productionTime + c.blowoutTime, 1e-13);
  close(c.blowoutVelocity, Math.sqrt(2 * G * D / 0.02), 1e-12);
  assert.ok(c.peakLiquidRate > c.meanLiquidRate);
  // plenty of gas: the flowline pressure reaches the full-riser head before the riser is full
  const stable = severeSluggingCycle({ ...o, massRateG: 1.2 * 5 * A });
  assert.equal(stable.severeSlugging, false);
  assert.equal(stable.penetrationLength, 0);
});

test('mitigation: choking and gas lift shift both criteria towards stability', () => {
  const base = { vsl: 0.2, vsg: 0.5, Psep: 2e5, rhoL: 1000, alphaPipe: 0.9, L: 100, riserHeight: 10, D: 0.1 };
  const none = severeSluggingMitigation(base);
  assert.equal(none.base.taitel.stable, false);
  assert.deepEqual(none.mitigated.taitel, none.base.taitel);
  close(none.chokeDpForStability, none.base.taitel.PsepMin - 2e5, 1e-12);
  const choked = severeSluggingMitigation({ ...base, chokeDp: none.chokeDpForStability * 1.01 });
  assert.equal(choked.mitigated.taitel.stable, true);
  const lifted = severeSluggingMitigation({ ...base, vsgLift: 3 });
  close(lifted.riserVoidLift, 3 / (1.2 * 3.2 + 0.35 * Math.sqrt(G * 0.1)), 1e-13);
  assert.ok(lifted.mitigated.boe.piSS > lifted.base.boe.piSS);
  assert.ok(lifted.chokeDpForStability < none.chokeDpForStability);
});

test('surge volume and pigging slug volume by hand', () => {
  // 20 m³ slug arriving at 0.5 m³/s (40 s) drained at 0.2 m³/s: 20 − 8 = 12 m³
  const s = surgeVolume({ slugVolume: 20, inflowRate: 0.5, drainRate: 0.2 });
  close(s.slugDuration, 40, 1e-13); close(s.surgeVolume, 12, 1e-13); close(s.drainedDuringSlug, 8, 1e-13);
  assert.equal(surgeVolume({ slugVolume: 20, inflowRate: 0.5, drainRate: 0.7 }).surgeVolume, 0);
  close(surgeVolume({ slugVolume: 20, inflowRate: 0.5, drainRate: 0.2, margin: 1.2 }).surgeVolume, 14.4, 1e-13);
  // 10 km of 0.3 m line, holdup 0.25, no-slip holdup 0.05: V = π/4 × 0.09 × 10000 × 0.20 = 141.37 m³
  const p = pigSlugVolume({ D: 0.3, holdup: 0.25, lambdaL: 0.05, length: 10000, vm: 4 });
  close(p.slugVolume, Math.PI / 4 * 0.09 * 10000 * 0.2, 1e-13);
  close(p.inventory, Math.PI / 4 * 0.09 * 10000 * 0.25, 1e-13);
  close(p.slugLength, 2000, 1e-13); close(p.slugDuration, 500, 1e-13);
  const seg = pigSlugVolume({ D: 0.3, holdup: [0.25, 0.5], lambdaL: 0.05, length: [5000, 5000] });
  close(seg.slugVolume, Math.PI / 4 * 0.09 * 5000 * (0.2 + 0.45), 1e-13);
});

test('invalid arguments throw RangeError', () => {
  assert.throws(() => slugBodyHoldupGregory(-1), /vm/);
  assert.throws(() => boeCriterion({ vsl: 0.1, vsg: 1, P: 1e5, rhoL: 1000, alphaPipe: 0, L: 10 }), /alphaPipe/);
  assert.throws(() => slugGrowthCriterion({ vm: 1, D: 0.1, filmHoldup: 0.95, vLF: 0, slugHoldup: 0.9 }), RangeError);
  assert.throws(() => surgeVolume({ slugVolume: 1, inflowRate: 0, drainRate: 0 }), /inflowRate/);
});
