import test from 'node:test';
import assert from 'node:assert/strict';
import {
  beggsBrill, stratifiedTwoFluid, pressureGradient, gradientAt, marchPipeline, turndownCurve, hydraulicCapacity,
  erosionalVelocity,
} from '../../src/engine/flow/multiphase.js';
import { G, darcyFrictionFactor } from '../../src/engine/flow/driftflux.js';
import { stratifiedForces } from '../../src/engine/flow/flowregime.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg ?? ''} got ${a}, expected ${b}`);
const oilGas = { D: 0.1, rhoL: 850, rhoG: 20, muL: 2e-3, muG: 1.5e-5, sigma: 0.025 };
const rough = 4.5e-5;

test('Beggs–Brill reduces exactly to the single-phase Darcy gradient (liquid and gas)', () => {
  for (const theta of [0, 0.3, -0.3, Math.PI / 2]) {
    const l = beggsBrill({ ...oilGas, vsl: 2, vsg: 0, theta, roughness: rough });
    const f = darcyFrictionFactor(850 * 2 * 0.1 / 2e-3, rough / 0.1);
    assert.equal(l.holdup, 1);
    close(l.gradFriction, f * 850 * 4 / (2 * 0.1), 1e-13);
    close(l.gradGravity, 850 * G * Math.sin(theta), 1e-13);
    close(l.gradTotal, l.gradFriction + l.gradGravity, 1e-13);
    const g = beggsBrill({ ...oilGas, vsl: 0, vsg: 8, theta, roughness: rough });
    assert.equal(g.holdup, 0);
    close(g.gradFriction, darcyFrictionFactor(20 * 8 * 0.1 / 1.5e-5, rough / 0.1) * 20 * 64 / (2 * 0.1), 1e-13);
  }
});

test('Beggs–Brill tends continuously to the single-phase liquid gradient as the gas rate vanishes', () => {
  for (const theta of [0, 0.2, Math.PI / 2]) {
    const one = beggsBrill({ ...oilGas, vsl: 2, vsg: 0, theta, roughness: rough });
    const two = beggsBrill({ ...oilGas, vsl: 2, vsg: 2e-6, theta, roughness: rough });
    close(two.holdup, 1, 1e-5, `holdup θ=${theta}`);
    close(two.gradTotal, one.gradTotal, 1e-4, `gradient θ=${theta}`);
  }
});

test('Beggs–Brill hand calculation: horizontal intermittent flow, D = 0.1 m, vsl = 1, vsg = 2 m/s', () => {
  // λ = 1/3; Fr = 9/(9.80665 × 0.1) = 9.1774; L3 = 0.1 × 3^1.4516 = 0.4927 < Fr ≤ L1 = 316 × 3^-0.302 = 226.8 → intermittent
  // HL(0) = 0.845 λ^0.5351 / Fr^0.0173 = 0.845 × 0.55551 / 1.03909 = 0.45175
  // y = λ/HL² = 1.6334, ln y = 0.49066, s = 0.49066/(−0.0523 + 1.56128 − 0.21005 + 0.00107) = 0.37743, e^s = 1.4585
  const r = beggsBrill({ ...oilGas, vsl: 1, vsg: 2, roughness: rough });
  assert.equal(r.regime, 'intermittent');
  close(r.detail.Fr, 9.1774, 1e-4);
  close(r.detail.L3, 0.4927, 1e-3);
  close(r.detail.L1, 226.8, 1e-3);
  close(r.holdup, 0.45175, 2e-4);
  close(r.detail.y, 1.6334, 5e-4);
  close(r.detail.ftp / r.detail.fn, 1.4585, 1e-3);
  const rhoN = 850 / 3 + 20 * 2 / 3;
  const muN = 2e-3 / 3 + 1.5e-5 * 2 / 3;
  close(r.detail.fn, darcyFrictionFactor(rhoN * 3 * 0.1 / muN, rough / 0.1), 1e-13);
  close(r.gradFriction, r.detail.ftp * rhoN * 9 / 0.2, 1e-13);
  assert.equal(r.gradGravity, 0);
});

test('Beggs–Brill hand calculation: vertical upflow inclination correction and Payne factor 0.924', () => {
  // Nlv = vsl (ρL/(g σ))^¼ = (850/0.245166)^¼ = 7.6734
  // C = (1 − λ) ln(2.96 λ^0.305 Nlv^−0.4473 Fr^0.0978) = 0.66667 × ln(2.96 × 0.71528 × 0.40192 × 1.24209) = 0.03694
  // ψ = 1 + C [sin 162° − sin³162°/3] = 1 + 0.03694 × 0.29918 = 1.01105; HL = 0.45175 × 1.01105 × 0.924 = 0.4220
  const r = beggsBrill({ ...oilGas, vsl: 1, vsg: 2, theta: Math.PI / 2, roughness: rough });
  close(r.detail.Nlv, 7.6734, 1e-4);
  close(r.detail.C, 0.03694, 2e-3);
  close(r.detail.psi, 1.01105, 1e-5);
  close(r.holdup, 0.4220, 3e-4);
  close(r.gradGravity, (r.holdup * 850 + (1 - r.holdup) * 20) * G, 1e-13);
  const noPayne = beggsBrill({ ...oilGas, vsl: 1, vsg: 2, theta: Math.PI / 2, roughness: rough, payne: false });
  close(r.holdup / noPayne.holdup, 0.924, 1e-12);
  // downhill: all-pattern coefficients and factor 0.685; holdup below the horizontal value
  const down = beggsBrill({ ...oilGas, vsl: 1, vsg: 2, theta: -0.3, roughness: rough });
  const downNo = beggsBrill({ ...oilGas, vsl: 1, vsg: 2, theta: -0.3, roughness: rough, payne: false });
  close(down.holdup / downNo.holdup, 0.685, 1e-12);
  const Cd = (2 / 3) * Math.log(4.70 * (1 / 3) ** -0.3692 * 7.6734 ** 0.1244 * 9.1774 ** -0.5056);
  close(downNo.detail.C, Cd, 1e-3);
  assert.ok(down.holdup < 0.45175);
});

test('Beggs–Brill: transition pattern interpolates continuously between segregated and intermittent', () => {
  // λ = 0.5: L2 = 0.00512, L3 = 0.2735 (Fr = vm²/(gD))
  const at = (Fr) => { const vm = Math.sqrt(Fr * G * 0.1); return beggsBrill({ ...oilGas, vsl: vm / 2, vsg: vm / 2, roughness: rough }); };
  assert.equal(at(0.2734).regime, 'transition');
  assert.equal(at(0.2737).regime, 'intermittent');
  close(at(0.2734).holdup, at(0.2737).holdup, 2e-3);
  assert.equal(at(0.00511).regime, 'segregated');
  assert.equal(at(0.00513).regime, 'transition');
  close(at(0.00511).holdup, at(0.00513).holdup, 2e-3);
});

test('Beggs–Brill acceleration term: Ek = ρs vm vsg / P and total = (friction + gravity)/(1 − Ek)', () => {
  const r = beggsBrill({ ...oilGas, vsl: 0.5, vsg: 15, theta: 0.1, roughness: rough, P: 5e5 });
  close(r.Ek, r.rhoM * 15.5 * 15 / 5e5, 1e-13);
  close(r.gradTotal, (r.gradFriction + r.gradGravity) / (1 - r.Ek), 1e-13);
  close(r.gradAcceleration, r.gradTotal * r.Ek, 1e-12);
});

test('no-slip limit: drift flux with C0 = 1, vd = 0 equals the homogeneous gradient', () => {
  const p = { ...oilGas, vsl: 1.2, vsg: 3, theta: 0.4, roughness: rough, P: 2e6 };
  const h = pressureGradient({ ...p, method: 'homogeneous' });
  const d = pressureGradient({ ...p, method: 'drift-flux', closure: 'slug', C0: 1, vd: 0 });
  close(d.holdup, 1.2 / 4.2, 1e-13);
  close(d.gradTotal, h.gradTotal, 1e-13);
  close(h.gradGravity, h.rhoNs * G * Math.sin(0.4), 1e-13);
  // every method returns the single-phase gradient without gas
  for (const method of ['homogeneous', 'drift-flux', 'mechanistic']) {
    const s = pressureGradient({ ...oilGas, method, vsl: 2, vsg: 0, theta: 0.1, roughness: rough });
    assert.equal(s.holdup, 1);
    close(s.gradTotal, darcyFrictionFactor(850 * 2 * 0.1 / 2e-3, rough / 0.1) * 850 * 4 / 0.2 + 850 * G * Math.sin(0.1), 1e-13);
  }
});

test('stratified two-fluid model: combined balance is satisfied and the gradient equals each phase balance', () => {
  for (const [vsl, vsg, theta] of [[0.02, 1, 0], [0.05, 2, -0.02], [0.03, 4, 0.005]]) {
    const s = stratifiedTwoFluid({ ...oilGas, vsl, vsg, theta });
    const f = stratifiedForces({ ...oilGas, h: s.hL, vL: s.vL, vG: s.vG, theta });
    assert.ok(Math.abs(s.residual) < 1e-6 * Math.max(1, Math.abs(f.tauWL * f.SL / f.AL)), `residual ${s.residual}`);
    // gas-phase momentum balance: −dP/dx = (τwG SG + τi Si)/AG + ρG g sinθ
    const gas = (f.tauWG * f.SG + f.tauI * f.Si) / f.AG + 20 * G * Math.sin(theta);
    // liquid-phase momentum balance: −dP/dx = (τwL SL − τi Si)/AL + ρL g sinθ
    const liq = (f.tauWL * f.SL - f.tauI * f.Si) / f.AL + 850 * G * Math.sin(theta);
    assert.ok(Math.abs(gas - liq) < 1e-6 * Math.max(1, Math.abs(gas)));
    assert.ok(Math.abs(s.gradTotal - gas) < 1e-6 * Math.max(1, Math.abs(gas)), `${s.gradTotal} vs ${gas}`);
    close(s.vL * s.holdup, vsl, 1e-12);
  }
  // a downward slope thins the liquid layer
  assert.ok(stratifiedTwoFluid({ ...oilGas, vsl: 0.05, vsg: 2, theta: -0.05 }).holdup < stratifiedTwoFluid({ ...oilGas, vsl: 0.05, vsg: 2 }).holdup);
});

test('mechanistic method selects models by flow pattern', () => {
  const s = pressureGradient({ ...oilGas, method: 'mechanistic', vsl: 0.05, vsg: 2, roughness: rough });
  assert.match(s.regime, /^stratified/); assert.equal(s.detail.model, 'two-fluid stratified');
  const i = pressureGradient({ ...oilGas, method: 'mechanistic', vsl: 1, vsg: 2, roughness: rough });
  assert.equal(i.regime, 'intermittent'); assert.equal(i.detail.model, 'slug unit cell');
  const a = pressureGradient({ ...oilGas, method: 'mechanistic', vsl: 0.1, vsg: 15, roughness: rough });
  assert.equal(a.regime, 'annular'); assert.match(a.detail.fallback, /annular/);
  assert.throws(() => pressureGradient({ ...oilGas, method: 'nope', vsl: 1, vsg: 1 }), RangeError);
});

test('API RP 14E erosional velocity: C = 100 and ρ = 1 lb/ft³ gives 100 ft/s', () => {
  close(erosionalVelocity(16.018463), 30.48, 1e-9);
  close(erosionalVelocity(100, 125), 1.2199 * 125 / 10, 1e-4);
  const r = pressureGradient({ ...oilGas, method: 'homogeneous', vsl: 1, vsg: 3 });
  close(r.erosionalRatio, 4 / erosionalVelocity(r.rhoNs), 1e-13);
});

const liquidOnly = () => ({ rhoL: 850, rhoG: 20, muL: 2e-3, muG: 1.5e-5, sigma: 0.025 });
const gasRT = 8.314462618 * 300 / 0.018;
const fluidAt = (P) => ({ rhoL: 850, rhoG: P / gasRT, muL: 2e-3, muG: 1.5e-5, sigma: 0.025 });
const route = [{ x: 0, z: 0 }, { x: 2000, z: -20 }, { x: 4000, z: 10 }, { x: 4000, z: 110 }];

test('marching solver, incompressible single-phase liquid on an inclined line: ΔP = (f ρ v²/2D + ρ g sinθ) L', () => {
  const D = 0.1;
  const A = Math.PI * D * D / 4;
  const v = 1.5;
  const profile = [{ x: 0, z: 0 }, { x: 800, z: 40 }, { x: 1500, z: 10 }];
  const L = Math.hypot(800, 40) + Math.hypot(700, 30);
  const f = darcyFrictionFactor(850 * v * D / 2e-3, rough / D);
  const expected = f * 850 * v * v / (2 * D) * L + 850 * G * 10;
  for (const method of ['beggs-brill', 'drift-flux', 'homogeneous', 'mechanistic']) {
    const r = marchPipeline({ profile, D, roughness: rough, fluidAt: liquidOnly, temperature: 300, pIn: 50e5, massRateL: 850 * v * A, massRateG: 0, method });
    close(r.dP, expected, 1e-10, method);
    close(r.liquidInventory, A * L, 1e-12);
    close(r.dPGravity, 850 * G * 10, 1e-10);
    assert.equal(r.solvedFor, 'pOut');
    assert.equal(r.stations[0].regime, 'single-phase-liquid');
  }
});

test('marching solver, isothermal ideal gas in a horizontal line: P1² − P2² = G² RT [f L/D + 2 ln(P1/P2)]', () => {
  const D = 0.2;
  const A = Math.PI * D * D / 4;
  const L = 20000;
  const m = 6;
  const Gm = m / A;
  const r = marchPipeline({ profile: [{ x: 0, z: 0 }, { x: L, z: 0 }], D, roughness: rough, fluidAt, temperature: 300, pIn: 60e5, massRateL: 0, massRateG: m, maxStep: 100 });
  // Reynolds number G D/μ is constant along the line, so f is constant.
  const f = darcyFrictionFactor(Gm * D / 1.5e-5, rough / D);
  const lhs = 60e5 ** 2 - r.pOut ** 2;
  const rhs = Gm * Gm * gasRT * (f * L / D + 2 * Math.log(60e5 / r.pOut));
  close(lhs, rhs, 1e-7);
  assert.ok(r.dP > 1e5);
  assert.equal(r.liquidInventory, 0);
});

test('marching solver: inlet- and outlet-pressure specifications are mutually consistent and the rate solve recovers the rate', () => {
  const base = { profile: route, D: 0.1, roughness: rough, fluidAt, temperature: 300 };
  const fwd = marchPipeline({ ...base, pIn: 30e5, massRateL: 5, massRateG: 0.3 });
  assert.equal(fwd.solvedFor, 'pOut');
  assert.ok(fwd.pOut < 30e5 && fwd.pOut > 0);
  const back = marchPipeline({ ...base, pOut: fwd.pOut, massRateL: 5, massRateG: 0.3 });
  assert.equal(back.solvedFor, 'pIn');
  close(back.pIn, 30e5, 1e-6);
  const rate = marchPipeline({ ...base, pIn: 30e5, pOut: fwd.pOut, gasMassFraction: 0.3 / 5.3 });
  assert.equal(rate.solvedFor, 'massRate');
  close(rate.massRate, 5.3, 1e-5);
  close(rate.massRateG, 0.3, 1e-5);
  // components add up to the total drop
  close(fwd.dPFriction + fwd.dPGravity + fwd.dPAcceleration, fwd.dP, 2e-3);
  // station output
  const st = fwd.stations[10];
  for (const k of ['s', 'x', 'z', 'P', 'T', 'holdup', 'vsl', 'vsg', 'vm', 'vL', 'vG', 'gradFriction', 'gradGravity', 'gradAcceleration', 'erosionalRatio']) assert.ok(Number.isFinite(st[k]), k);
  assert.equal(typeof st.regime, 'string');
  close(st.vsl, 5 / (850 * Math.PI * 0.01 / 4), 1e-12);
  close(st.vsg, 0.3 / (st.P / gasRT * Math.PI * 0.01 / 4), 1e-12);
  assert.ok(fwd.liquidInventory > 0 && fwd.liquidInventory < Math.PI * 0.01 / 4 * fwd.length);
  close(fwd.meanHoldup, fwd.liquidInventory / (Math.PI * 0.01 / 4 * fwd.length), 1e-13);
  assert.ok(fwd.maxErosionalRatio > 0);
});

test('marching solver: step refinement converges (4th-order integration of a smooth gradient)', () => {
  const base = { profile: [{ x: 0, z: 0 }, { x: 5000, z: 50 }], D: 0.1, roughness: rough, fluidAt, temperature: 300, pIn: 40e5, massRateL: 3, massRateG: 0.5, method: 'drift-flux', closure: 'slug', C0: 1.2, vd: 0.2 };
  const p = [1000, 500, 250].map((maxStep) => marchPipeline({ ...base, maxStep }).pOut);
  const ref = marchPipeline({ ...base, maxStep: 20 }).pOut;
  const e = p.map((v) => Math.abs(v - ref));
  assert.ok(e[2] < 1, `fine-grid error ${e[2]} Pa`);
  assert.ok(e[0] > e[1] && e[1] > e[2]);
  assert.ok(Math.log2(e[0] / e[1]) > 3, `observed order ${Math.log2(e[0] / e[1])}`);
});

test('marching solver guards against over- and under-constrained boundary specification', () => {
  const base = { profile: route, D: 0.1, fluidAt, temperature: 300 };
  assert.throws(() => marchPipeline({ ...base, pIn: 30e5, pOut: 10e5, massRateL: 5, massRateG: 0.3 }), /Over-constrained/);
  assert.throws(() => marchPipeline({ ...base, pIn: 30e5 }), /Under-constrained/);
  assert.throws(() => marchPipeline({ ...base, massRateL: 5, massRateG: 0.3 }), /Under-constrained/);
  assert.throws(() => marchPipeline({ ...base, pIn: 30e5, massRateL: 5, massRateG: 0.3, massRate: 5 }), /Over-constrained/);
  assert.throws(() => marchPipeline({ ...base, pIn: 30e5, massRateL: 5 }), /Under-constrained/);
  assert.throws(() => marchPipeline({ ...base, pIn: 2e5, massRateL: 50, massRateG: 3 }), /cannot be delivered|Pressure falls|choking/);
  assert.throws(() => marchPipeline({ ...base, pIn: 10e5, pOut: 30e5, gasMassFraction: 0.05 }), /No flow solution/);
  assert.throws(() => turndownCurve({ ...base, pOut: 10e5, pIn: 30e5, massRates: [1, 2] }), /Over-constrained/);
});

test('turndown curve of an uphill line has a gravity/friction minimum; capacity is the friction-dominated root', () => {
  const base = { profile: route, D: 0.1, roughness: rough, fluidAt, temperature: 300, gasMassFraction: 0.3 / 5.3 };
  const td = turndownCurve({ ...base, pOut: 10e5, massRates: [0.05, 0.1, 0.2, 0.5, 1, 2, 4, 8] });
  assert.equal(td.hasInteriorMinimum, true);
  assert.ok(td.dP[0] > td.minimum.dP && td.dP[7] > td.minimum.dP);
  // low rates: gravity dominates and liquid accumulates; high rates: friction dominates
  assert.ok(td.dPGravity[0] > td.dPFriction[0]);
  assert.ok(td.dPFriction[7] > td.dPGravity[7]);
  assert.ok(td.liquidInventory[0] > td.liquidInventory[7]);
  // capacity for a pressure budget above the minimum: two rates satisfy it, the higher one is returned
  const pIn = 10e5 + 1.5 * td.minimum.dP;
  const cap = hydraulicCapacity({ ...base, pIn, pOut: 10e5 });
  assert.equal(cap.multipleSolutions, true);
  assert.ok(cap.massRate > td.minimum.massRate);
  close(cap.result.pIn, pIn, 1e-6);
  const check = marchPipeline({ ...base, pOut: 10e5, massRate: cap.massRate });
  close(check.pIn, pIn, 1e-6);
});

test('gradientAt wrapper: mass rates and fluidAt(P, T) give the same result as explicit superficial velocities', () => {
  const A = Math.PI * 0.01 / 4;
  const w = gradientAt({ fluidAt, P: 20e5, T: 300, massRateL: 4, massRateG: 0.2, D: 0.1, theta: 0.05, roughness: rough });
  const fl = fluidAt(20e5);
  const d = pressureGradient({ ...fl, D: 0.1, theta: 0.05, roughness: rough, P: 20e5, vsl: 4 / (850 * A), vsg: 0.2 / (fl.rhoG * A) });
  close(w.gradTotal, d.gradTotal, 1e-13);
  const q = gradientAt({ fluidAt: (P) => ({ ...fluidAt(P), quality: 0.05 }), P: 20e5, T: 300, massRate: 4, D: 0.1 });
  close(q.massRateG, 0.2, 1e-13);
});
