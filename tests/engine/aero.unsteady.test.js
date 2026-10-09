import test from 'node:test';
import assert from 'node:assert/strict';
import {
  oyeTimeConstants, oyeDynamicInflowStep, pittPetersRate, createSeparationModel, oyeDynamicStallStep,
  beddoesLeishmanInit, beddoesLeishmanStep, towerShadowPotential, shearedWind, harmonicAmplitudes,
  rotationalSampling, unsteadyBem,
} from '../../src/engine/aero/unsteady.js';
import { bemSolve } from '../../src/engine/aero/bem.js';
import { nrel5mw } from '../../src/engine/aero/rotors.js';
import { genericPolar } from '../../src/engine/aero/polar.js';

const D2R = Math.PI / 180;
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} +- ${tol}`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${tol})`);

const rotor = nrel5mw();
const polar = genericPolar({ thickness: 0.18 });
const sep = createSeparationModel(polar);

test('oyeTimeConstants: hand values (Hansen Eqs. 9.11-9.12), induction capped at 0.5', () => {
  // a = 0.3, R = 63, V0 = 10, r = 0.7 R: tau1 = 1.1/(1 - 0.39) * 6.3 = 11.3607 s, tau2 = (0.39 - 0.26*0.49) tau1 = 2.9833 s
  const t = oyeTimeConstants({ a: 0.3, R: 63, V0: 10, r: 44.1 });
  close(t.tau1, 11.3607, 1e-4); close(t.tau2, 2.9833, 1e-4);
  close(oyeTimeConstants({ a: 0.9, R: 63, V0: 10, r: 0 }).tau1, 1.1 / 0.35 * 6.3, 1e-9);
  assert.throws(() => oyeTimeConstants({ a: 0.3, R: 63, V0: 0, r: 1 }), RangeError);
});

test('Oye dynamic-inflow filter: unit-step response equals the analytic two-time-constant solution', () => {
  const tau1 = 10; const tau2 = 3; const k = 0.6; const dt = 0.005;
  const st = { xPrev: 0, y: 0, W: 0 };
  const C = -1 + (1 - k) * tau1 / (tau1 - tau2);
  const exact = (t) => 1 - (1 - k) * tau1 / (tau1 - tau2) * Math.exp(-t / tau1) + C * Math.exp(-t / tau2);
  let t = 0;
  for (let i = 1; i <= 8000; i++) {
    t = i * dt;
    oyeDynamicInflowStep(st, { x: 1, dt, tau1, tau2 });
    if (i === 1) close(st.y, k, 2e-3, 'instant 60 % response of the first filter');
    if (i % 1000 === 0) close(st.W, exact(t), 2e-3, `t = ${t}`);
  }
  close(st.W, 1, 0.02);
});

test('Pitt-Peters uniform inflow: momentum theory is the steady state; apparent mass 8/3 rho R^3', () => {
  const rho = 1.225; const R = 40; const V = 10; const a = 0.25;
  const T = 2 * rho * Math.PI * R * R * a * V * V * (1 - a);
  close(pittPetersRate({ v: a * V, thrust: T, V, R, rho }), 0, 1e-12);
  close(pittPetersRate({ v: 0, thrust: T, V, R, rho }), T / (8 / 3 * rho * R ** 3), 1e-12);
});

test('Kirchhoff separation model reproduces the static lift: Cl = m (a - a0) ((1 + sqrt f)/2)^2', () => {
  close(sep.alpha0, -4 * D2R, 1e-6);
  close(sep.fAt(2 * D2R), 1, 1e-6, 'attached');
  assert.equal(sep.fAt(80 * D2R), 0);
  let partly = 0;
  for (let deg = -3; deg <= 22; deg += 0.5) {
    const al = deg * D2R;
    const f = sep.fAt(al);
    if (f > 0.02 && f < 0.98) partly++;
    const lin = sep.slope * (al - sep.alpha0);
    // both representations of the static curve
    close(lin * ((1 + Math.sqrt(f)) / 2) ** 2, polar.clAt(al), f > 0 ? 6e-3 : 10, `Kirchhoff at ${deg}`);
    close(lin * f + sep.clFsAt(al) * (1 - f), polar.clAt(al), 6e-3, `decomposition at ${deg}`);
  }
  assert.ok(partly >= 5);
  close(sep.clMax, 1.45, 1e-3);
});

test('Oye dynamic stall: steady state is the static polar; a sudden incidence step delays the lift loss', () => {
  const st = { fs: sep.fAt(8 * D2R) };
  const args = { sep, alpha: 18 * D2R, W: 60, chord: 3, dt: 0.01 };
  oyeDynamicStallStep(st, args);
  assert.ok(st.cl > polar.clAt(18 * D2R) + 0.2, 'lift overshoot right after the step');
  for (let i = 0; i < 2000; i++) oyeDynamicStallStep(st, args);
  close(st.cl, polar.clAt(18 * D2R), 6e-3);
  // time constant tau = 4 c / W: after one tau the separation state has covered 1 - 1/e of the gap
  const s2 = { fs: 1 };
  const tau = 4 * 3 / 60;
  for (let i = 0; i < 20; i++) oyeDynamicStallStep(s2, { ...args, dt: tau / 20 });
  close((1 - s2.fs) / (1 - sep.fAt(18 * D2R)), 1 - Math.exp(-1), 1e-9);
});

test('Beddoes-Leishman attached flow: step response follows the Jones approximation of the Wagner function', () => {
  const U = 50; const c = 2;
  const a0 = 0; const a1 = 2 * D2R;
  const st = beddoesLeishmanInit(sep, a0);
  const dt = 0.0002;                                  // ds = 2 U dt / c = 0.01 semi-chords
  const lin = (al) => sep.slope * (al - sep.alpha0);
  const wag = (s) => 1 - 0.165 * Math.exp(-0.0455 * s) - 0.335 * Math.exp(-0.3 * s);
  for (let i = 1; i <= 8000; i++) {
    beddoesLeishmanStep(st, { sep, polar, alpha: a1, W: U, chord: c, dt });
    const s = i * 0.01;
    if (i === 200 || i === 1000 || i === 4000 || i === 8000) {
      const frac = (st.cl - lin(a0)) / (lin(a1) - lin(a0));
      close(frac, wag(s), 6e-3, `s = ${s}`);
    }
  }
  // half of the circulatory lift appears instantly: Wagner(0) = 0.5
  close(wag(0), 0.5, 1e-12);
});

test('Beddoes-Leishman: steady state equals the static polar before and beyond stall; drag returns to static', () => {
  for (const deg of [6, 12, 17, 24]) {
    const st = beddoesLeishmanInit(sep, 2 * D2R);
    for (let i = 0; i < 6000; i++) beddoesLeishmanStep(st, { sep, polar, alpha: deg * D2R, W: 60, chord: 3, dt: 0.01 });
    close(st.cl, polar.clAt(deg * D2R), 8e-3, `Cl at ${deg} deg`);
    close(st.cd, polar.cdAt(deg * D2R), 2e-3, `Cd at ${deg} deg`);
  }
});

test('Beddoes-Leishman: pitch oscillation through stall gives a hysteresis loop (upstroke lift above downstroke and above static)', () => {
  const U = 60; const c = 3; const kRed = 0.1;        // reduced frequency k = omega c / (2U)
  const om = kRed * 2 * U / c;
  const dt = 2 * Math.PI / om / 400;
  const run = (opts) => {
    const st = beddoesLeishmanInit(sep, 12 * D2R);
    let up = null; let down = null; let clMax = -Infinity;
    for (let i = 1; i <= 1200; i++) {
      const al = (12 + 6 * Math.sin(om * i * dt)) * D2R;
      beddoesLeishmanStep(st, { sep, polar, alpha: al, W: U, chord: c, dt, ...opts });
      if (i > 800) {
        clMax = Math.max(clMax, st.cl);
        const ph = (om * i * dt) % (2 * Math.PI);
        if (Math.abs(ph - Math.asin(0.5)) < 0.009) up = st.cl;               // 15 deg, rising
        if (Math.abs(ph - (Math.PI - Math.asin(0.5))) < 0.009) down = st.cl; // 15 deg, falling
      }
    }
    return { up, down, clMax };
  };
  const r = run({});
  const stat = polar.clAt(15 * D2R);
  assert.ok(r.up > stat + 0.1, `upstroke ${r.up} vs static ${stat}`);
  assert.ok(r.down < r.up - 0.2, `downstroke ${r.down}`);
  assert.ok(r.clMax > 1.45, 'dynamic maximum lift exceeds the static maximum');
  // the optional vortex-lift state adds lift during the stall event
  const v = run({ vortexLift: true });
  assert.ok(v.clMax > r.clMax);
});

test('towerShadowPotential and shearedWind: hand values', () => {
  close(towerShadowPotential({ x: -8, y: 0, towerDiameter: 4 }), 1 - (2 / 8) ** 2, 1e-14);
  close(shearedWind({ Vhub: 10, z: 45, hubHeight: 90, shearExp: 0.2 }).speed, 10 * 0.5 ** 0.2, 1e-12);
  close(shearedWind({ Vhub: 10, z: 45, hubHeight: 90, z0: 0.001 }).speed, 10 * Math.log(45000) / Math.log(90000), 1e-12);
  close(shearedWind({ Vhub: 10, z: 120, hubHeight: 90, veerDegPerM: 0.1 }).directionDeg, 3, 1e-12);
});

test('harmonicAmplitudes: recovers mean, amplitude and phase of a known signal', () => {
  const N = 48; const y = [];
  for (let i = 0; i < N; i++) { const p = 2 * Math.PI * i / N; y.push(5 + 3 * Math.cos(p - 0.4) + 1.5 * Math.cos(3 * p + 1.0)); }
  const h = harmonicAmplitudes(y, 4);
  close(h.mean, 5, 1e-12);
  close(h.harmonics[0].amp, 3, 1e-12); close(h.harmonics[0].phase, 0.4, 1e-12);
  close(h.harmonics[1].amp, 0, 1e-12);
  close(h.harmonics[2].amp, 1.5, 1e-12); close(h.harmonics[2].phase, -1.0, 1e-12);
});

test('rotational sampling: one blade sees 1P, the three-bladed rotor only 3P, 6P, 9P', () => {
  const rs = rotationalSampling({ rotor, U: 10, omega: 1.1, shearExp: 0.2, towerShadow: true, tiltDeg: 0, preconeDeg: 0 });
  const bh = rs.bladeFlapHarmonics.harmonics;
  assert.ok(bh[0].amp > 0.02 * rs.bladeFlapHarmonics.mean, '1P on the blade');
  assert.ok(bh[0].amp > bh[1].amp && bh[0].amp > bh[2].amp);
  const th = rs.rotorThrustHarmonics.harmonics;
  for (const n of [1, 2, 4, 5, 7, 8]) close(th[n - 1].amp / rs.rotorThrustHarmonics.mean, 0, 1e-10, `${n}P cancels`);
  assert.ok(th[2].amp > 1e-4 * rs.rotorThrustHarmonics.mean, '3P present');
  close(rs.frequency1P, 1.1 / (2 * Math.PI), 1e-12);
  // mean rotor thrust equals the sector-averaged BEM thrust
  rel(rs.rotorThrustHarmonics.mean, bemSolve({ rotor, U: 10, omega: 1.1, shearExp: 0.2, towerShadow: true, tiltDeg: 0, preconeDeg: 0, nSectors: 72 }).thrust, 1e-10);
});

test('unsteadyBem: steady uniform inflow reproduces the steady BEM exactly', () => {
  const u = unsteadyBem({ rotor, dt: 0.05, nSteps: 60, wind: 8, omega: 0.9, pitchDeg: 0, tiltDeg: 0, preconeDeg: 0 });
  const s = bemSolve({ rotor, U: 8, omega: 0.9, tiltDeg: 0, preconeDeg: 0 });
  rel(u.thrust[60], s.thrust, 1e-8); rel(u.torque[60], s.torque, 1e-8);
  rel(u.flapMoment[1][60], s.flapMoment, 1e-8);
  rel(u.power[60], s.power, 1e-8);
  close(u.azimuth[60], 0.9 * 3, 1e-9);
});

test('unsteadyBem with dynamic inflow: a pitch step overshoots the load change and relaxes to the new steady state', () => {
  const dt = 0.1; const nSteps = 900;
  const pitch = (t) => (t < 5 ? 0 : 3);
  const u = unsteadyBem({ rotor, dt, nSteps, wind: 9, omega: 1.0, pitchDeg: pitch, tiltDeg: 0, preconeDeg: 0 });
  const s0 = bemSolve({ rotor, U: 9, omega: 1.0, pitchDeg: 0, tiltDeg: 0, preconeDeg: 0 });
  const s1 = bemSolve({ rotor, U: 9, omega: 1.0, pitchDeg: 3, tiltDeg: 0, preconeDeg: 0 });
  rel(u.thrust[40], s0.thrust, 1e-6);
  const just = u.thrust[51];
  assert.ok(just < s1.thrust - 0.15 * (s0.thrust - s1.thrust), `undershoot ${just} below the new steady thrust ${s1.thrust}`);
  rel(u.thrust[nSteps], s1.thrust, 5e-3);
  // without dynamic inflow the response is quasi-steady
  const q = unsteadyBem({ rotor, dt, nSteps: 60, wind: 9, omega: 1.0, pitchDeg: pitch, tiltDeg: 0, preconeDeg: 0, dynamicInflow: 'none' });
  rel(q.thrust[51], s1.thrust, 1e-8);
});

test('unsteadyBem with shear: rotor thrust contains the blade-passing (3P) frequency; dynamic stall options run', () => {
  const om = 1.2; const dt = 2 * Math.PI / om / 120;
  const u = unsteadyBem({ rotor, dt, nSteps: 720, wind: 11, omega: om, shearExp: 0.2, towerShadow: true, tiltDeg: 0, preconeDeg: 0, dynamicInflow: 'none' });
  const last = Array.from(u.thrust.subarray(600, 720));
  const h = harmonicAmplitudes(last, 6);
  assert.ok(h.harmonics[2].amp > 20 * h.harmonics[0].amp && h.harmonics[2].amp > 20 * h.harmonics[1].amp, '3P dominates 1P and 2P');
  for (const mode of ['oye', 'bl']) {
    const d = unsteadyBem({ rotor, dt, nSteps: 240, wind: 11, omega: om, shearExp: 0.2, tiltDeg: 0, preconeDeg: 0, dynamicStall: mode });
    rel(d.thrust[240], u.thrust[240], 0.05, mode);
  }
  assert.throws(() => unsteadyBem({ rotor, dt: 0, nSteps: 10, wind: 8, omega: 1 }), RangeError);
});
