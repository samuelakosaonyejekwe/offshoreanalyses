import test from 'node:test';
import assert from 'node:assert/strict';
import * as X from '../../src/engine/pipeline/expansion.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${a} vs ${b}`);

const D = 0.3239; const t = 0.0127; const E = 207e9; const alpha = 11.7e-6;
const As = Math.PI / 4 * (D * D - (D - 2 * t) ** 2);       // 0.0124163 m2
const Ai = Math.PI / 4 * (D - 2 * t) ** 2;                 // 0.0699807 m2
const EA = E * As;

test('restrained axial force S0 = H - dp Ai (1 - 2 nu) - As E alpha dT: hand calculation', () => {
  // pressure term: 10e6*0.0699807*0.4 = 279 923 N
  // thermal term: 0.0124163*207e9*11.7e-6*60 = 0.0124163*145.314e6 = 1 804 266 N
  // S0 = 100 000 - 279 923 - 1 804 266 = -1 984 189 N
  const r = X.restrainedAxialForce({ H: 100e3, deltaPi: 10e6, Ai, As, deltaT: 60 });
  close(r.pressureTerm, -279923, 2e-5);
  close(r.thermalTerm, -1804266, 2e-5);
  close(r.S0, -1984189, 2e-5);
  close(r.freeStrain, alpha * 60 + 10e6 * Ai * 0.4 / EA, 1e-12);
});

test('anchor length = force / friction and end expansion = P^2/(2 f EA) for uniform heating (long line)', () => {
  const dT = 80; const f = 500;
  const P = EA * alpha * dT;
  const r = X.endExpansion({ L: 20000, As, Ai, friction: f, deltaT: dT });
  assert.equal(r.regime, 'long');
  close(r.anchorHot, P / f, 1e-6);
  close(20000 - r.anchorCold, P / f, 1e-6);
  close(r.expansionHot, P * P / (2 * f * EA), 1e-6);
  close(r.expansionCold, r.expansionHot, 1e-9);
  close(r.maxCompression, P, 1e-12);
  const u = X.endExpansionUniform({ P, f, EA, L: 20000 });
  close(u.anchorLength, P / f, 1e-12);
  close(u.expansion, P * P / (2 * f * EA), 1e-12);
  // force profile: linear build-up at f per metre from the free end
  close(-r.S[100], f * r.x[100], 1e-9);
});

test('short pipeline: single virtual anchor at mid-length, expansion = (P L/2 - f L^2/8)/EA', () => {
  const dT = 80; const f = 500; const L = 4000;
  const P = EA * alpha * dT;
  assert.ok(L < 2 * P / f);
  const r = X.endExpansion({ L, As, Ai, friction: f, deltaT: dT });
  assert.equal(r.regime, 'short');
  close(r.anchorHot, L / 2, 1e-9);
  close(r.expansionHot, (P * L / 2 - f * L * L / 8) / EA, 1e-6);
  close(r.maxCompression, f * L / 2, 1e-9);
  close(X.endExpansionUniform({ P, f, EA, L }).expansion, r.expansionHot, 1e-6);
});

test('end expansion with decaying temperature profile equals the direct integral of the released strain', () => {
  const f = 400; const L = 30000; const lam = 8000;
  const dT = (x) => 90 * Math.exp(-x / lam);
  const r = X.endExpansion({ L, As, Ai, friction: f, deltaT: dT, deltaPi: 12e6, n: 6000 });
  // hot end: solve f*xa = P(xa) by bisection, then integrate (P(x) - f x)/EA with Simpson
  const P = (x) => 12e6 * Ai * 0.4 + EA * alpha * dT(x);
  let lo = 0; let hi = L;
  for (let i = 0; i < 100; i++) { const m = 0.5 * (lo + hi); if (f * m < P(m)) lo = m; else hi = m; }
  const xa = lo; const n = 2000; const h = xa / n; let s = 0;
  for (let i = 0; i <= n; i++) s += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * (P(i * h) - f * i * h) / EA;
  close(r.anchorHot, xa, 1e-5);
  close(r.expansionHot, s * h / 3, 1e-5);
  assert.ok(r.expansionHot > r.expansionCold && r.expansionCold > 0);
  const T = X.temperatureProfile({ x: [0, lam, 1e9], TinC: 95, TambC: 5, lambda: lam });
  close(T[0], 95); close(T[1], 5 + 90 / Math.E, 1e-12); close(T[2], 5, 1e-12);
  close(X.temperatureProfile({ x: 1000, TinC: 95, TambC: 5, massFlow: 50, cp: 2000, U: 5, Dref: 0.3 }), 5 + 90 * Math.exp(-1000 * 5 * Math.PI * 0.3 / 1e5), 1e-12);
});

test('tie-in spring: reaction R = k delta and delta = (P - R)^2/(2 f EA)', () => {
  const dT = 80; const f = 500; const k = 2e6;
  const P = EA * alpha * dT;
  const r = X.endExpansion({ L: 20000, As, Ai, friction: f, deltaT: dT, kHot: k });
  close(r.Rhot, k * r.expansionHot, 1e-6);
  close(r.expansionHot, (P - r.Rhot) ** 2 / (2 * f * EA), 1e-5);
  assert.ok(r.expansionHot < X.endExpansionUniform({ P, f, EA }).expansion);
  close(X.endExpansionUniform({ P, f, EA, R: r.Rhot }).expansion, r.expansionHot, 1e-5);
  const p = X.pressureProfile({ x: [0, 1000], pIn: 20e6, frictionGradient: 100, rho: 800, z: [0, -50] });
  close(p[1], 20e6 - 1e5 + 800 * 9.80665 * 50, 1e-12);
});

test('guided cantilever spool: sigma = 3 E D delta / L^2', () => {
  const I = Math.PI / 64 * (D ** 4 - (D - 2 * t) ** 4);
  const r = X.spoolGuidedCantilever({ delta: 0.5, D, I, Lleg: 30, sigmaAllow: 300e6 });
  close(r.stress, 3 * E * D * 0.5 / 900, 1e-12);
  close(r.force, 12 * E * I * 0.5 / 27000, 1e-12);
  close(r.moment, r.force * 30 / 2, 1e-12);
  close(X.spoolGuidedCantilever({ delta: 0.5, D, I, Lleg: r.requiredLeg }).stress, 300e6, 1e-12);
});

const W = 1000; const mu = 0.5; const L = 2000; const nElem = 200;
const freeStrain = alpha * 80;
const dP = EA * freeStrain;

test('walking model: no walking for symmetric loading without slope, tension or thermal gradient', () => {
  const r = X.axialWalkingModel({ L, EA, W, mu, steps: X.uniformCycleSteps({ nElem, freeStrain }), cycles: 4, uMob: 1e-4 });
  for (const w of r.walkPerCycle) assert.ok(Math.abs(w) < 1e-9, `walk ${w}`);
  // ends breathe symmetrically; cold-state residual tension from friction is symmetric about mid-length
  close(r.endDisplacementStart[3], -r.endDisplacementEnd[3], 1e-6);
  close(r.coldForce[10], r.coldForce[nElem - 11], 1e-6);
  // hot state: fully mobilised short line, compression at mid-length = mu W L/2
  close(-Math.min(...r.hotForce), mu * W * L / 2, 2e-2);
});

test('walking model: downhill walking on a slope matches the Carr-Sinclair-Bruton expression', () => {
  const phi = 2 * Math.PI / 180;
  const r = X.axialWalkingModel({ L, EA, W, mu, steps: X.uniformCycleSteps({ nElem, freeStrain }), cycles: 4, slopeRad: phi, uMob: 1e-4 });
  // [dP + W L sin(phi) - mu W L cos(phi)] L tan(phi)/(EA mu)
  const carr = (dP + W * L * Math.sin(phi) - mu * W * L * Math.cos(phi)) * L * Math.tan(phi) / (EA * mu);
  assert.ok(r.walkRate > 0);                                  // towards +x, i.e. downhill
  close(r.walkRate, carr, 5e-3);
  close(r.walkPerCycle[2], r.walkPerCycle[3], 1e-9);          // steady rate after the first cycle
  const a = X.walkingSlope({ deltaP: dP, W, L, mu, slopeRad: phi, EA });
  close(a.walkPerCycle, carr, 1e-12);
  close(a.anchorSeparation, L * Math.tan(phi) / mu, 1e-12);
  // reversing the slope reverses the walk
  const up = X.axialWalkingModel({ L, EA, W, mu, steps: X.uniformCycleSteps({ nElem, freeStrain }), cycles: 3, slopeRad: -phi, uMob: 1e-4 });
  close(up.walkRate, -r.walkRate, 1e-6);
});

test('walking model: end tension (SCR) matches (dP + SR - mu W L) SR/(EA mu W)', () => {
  const SR = 2e5;
  const r = X.axialWalkingModel({ L, EA, W, mu, steps: X.uniformCycleSteps({ nElem, freeStrain }), cycles: 4, tensionAtEnd: SR, uMob: 1e-4 });
  const carr = (dP + SR - mu * W * L) * SR / (EA * mu * W);
  close(r.walkRate, carr, 5e-3);
  close(X.walkingEndTension({ deltaP: dP, SR, W, L, mu, EA }).walkPerCycle, carr, 1e-12);
  assert.equal(X.walkingEndTension({ deltaP: 1e4, SR, W, L, mu, EA }).walkPerCycle, 0);   // not fully mobilised
});

test('walking model: thermal transient walks towards the cold end; a stiff anchor arrests slope walking', () => {
  const tr = X.axialWalkingModel({ L, EA, W, mu, steps: X.thermalTransientSteps({ L, nElem, alphaDeltaT: freeStrain, lambda: 3000, nHeat: 40 }), cycles: 4, uMob: 1e-4 });
  assert.ok(tr.walkRate > 1e-3, `transient walk ${tr.walkRate}`);
  close(tr.walkPerCycle[2], tr.walkPerCycle[3], 1e-6);
  const phi = 2 * Math.PI / 180;
  const an = X.axialWalkingModel({ L, EA, W, mu, steps: X.uniformCycleSteps({ nElem, freeStrain }), cycles: 4, slopeRad: phi, uMob: 1e-4,
    anchor: { x: 0, stiffness: 1e10 } });
  assert.ok(Math.abs(an.walkRate) < 1e-9);
  const lim = X.walkingArrestForce({ W, L, mu, slopeRad: phi });
  close(an.anchorForceMax, lim.rigidAnchorTension, 1e-2);       // cool-down: (mu cos + sin) W L
  close(-an.anchorForceMin, lim.rigidAnchorCompression, 1e-2);  // heat-up: (mu cos - sin) W L
  close(lim.holdBackTension, W * L * Math.sin(phi), 1e-12);
  // a constant hold-back tension W L sin(phi) at the upper end also stops the walk
  const hb = X.axialWalkingModel({ L, EA, W, mu, steps: X.uniformCycleSteps({ nElem, freeStrain }), cycles: 4, slopeRad: phi, uMob: 1e-4,
    tensionAtStart: lim.holdBackTension });
  assert.ok(Math.abs(hb.walkRate) < 2e-4, `held ${hb.walkRate}`);
});
