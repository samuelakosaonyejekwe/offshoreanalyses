import test from 'node:test';
import assert from 'node:assert/strict';
import {
  reynolds, frictionFactor, minorLoss, liquidLinePressureDrop, gasFlowGeneral, gasTransmission, gasJouleThomson,
  lineProfile, noiseVelocityLimit, velocityChecks, sizeLine, pumpPower, compressorPolytropic, waveSpeed, joukowsky,
  waterHammerMOC,
} from '../../src/engine/flow/singlephase.js';
import { zFactor, M_AIR, R_GAS } from '../../src/engine/flow/blackoil.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const PSI = 6894.757293168;
const G = 9.80665;

test('Colebrook–White by hand iteration at Re = 1e5, ε/D = 0.002 (Moody chart ≈ 0.025)', () => {
  // fixed-point iteration 1/√f = −2 log10(ε/(3.7D) + 2.51/(Re√f)), written out
  let f = 0.02;
  for (let i = 0; i < 100; i++) f = (-2 * Math.log10(0.002 / 3.7 + 2.51 / (1e5 * Math.sqrt(f)))) ** -2;
  const got = frictionFactor({ Re: 1e5, relRough: 0.002 });
  assert.ok(rel(got, f) < 1e-12);
  assert.ok(Math.abs(got - 0.0251) < 1e-4, `f = ${got}`);
  // fully rough limit (von Kármán): 1/√f = 2 log10(3.7 D/ε) → ε/D = 0.01: f = 0.0379
  assert.ok(rel(frictionFactor({ Re: 1e9, relRough: 0.01 }), (2 * Math.log10(370)) ** -2) < 1e-4);
  // smooth pipe, Prandtl law 1/√f = 2 log10(Re√f) − 0.8 at Re = 1e5: f = 0.0180
  const fs = frictionFactor({ Re: 1e5, relRough: 0 });
  assert.ok(Math.abs(1 / Math.sqrt(fs) - (2 * Math.log10(1e5 * Math.sqrt(fs)) - 0.8)) < 2e-3);
  assert.ok(Math.abs(fs - 0.0180) < 1e-4);
});

test('explicit friction-factor formulas track Colebrook: Haaland within 2 %, Churchill within 3 %', () => {
  for (const Re of [5e3, 1e4, 1e5, 1e6, 1e7, 1e8]) {
    for (const rr of [0, 1e-5, 1e-4, 1e-3, 1e-2]) {
      const c = frictionFactor({ Re, relRough: rr });
      assert.ok(rel(frictionFactor({ Re, relRough: rr, method: 'haaland' }), c) < 0.02, `Haaland Re ${Re} rr ${rr}`);
      assert.ok(rel(frictionFactor({ Re, relRough: rr, method: 'churchill' }), c) < 0.03, `Churchill Re ${Re} rr ${rr}`);
    }
  }
  assert.ok(rel(frictionFactor({ Re: 500, method: 'churchill' }), 64 / 500) < 1e-3);
  assert.equal(frictionFactor({ Re: 1000 }), 0.064);
});

test('Hagen–Poiseuille laminar pressure drop ΔP = 32 μ L v/D²', () => {
  const D = 0.05;
  const mu = 0.1;
  const rho = 900;
  const v = 0.5; // Re = 225
  const r = liquidLinePressureDrop({ Q: v * Math.PI * D * D / 4, D, rho, mu, L: 100 });
  assert.ok(rel(r.Re, 225) < 1e-12);
  assert.ok(rel(r.dP, 32 * mu * 100 * v / (D * D)) < 1e-12);
});

test('liquid line with elevation profile and fittings — hand calculation', () => {
  // water, D = 0.2 m, v = 2 m/s: Re = 4e5, ε/D = 2.25e-4
  const D = 0.2;
  const rho = 1000;
  const mu = 1e-3;
  const Q = 2 * Math.PI * D * D / 4;
  const f = frictionFactor({ Re: 4e5, relRough: 4.5e-5 / D });
  const r = liquidLinePressureDrop({ Q, D, rho, mu, K: 3, Pin: 5e6,
    profile: [{ x: 0, z: 0 }, { x: 400, z: -30 }, { x: 1000, z: 25 }] });
  const fric = f * (1000 / D) * rho * 4 / 2;
  assert.ok(rel(r.dPfriction, fric) < 1e-12);
  assert.ok(rel(r.dPminor, 3 * rho * 4 / 2) < 1e-12);
  assert.ok(rel(r.dPelevation, rho * G * 25) < 1e-12);
  assert.ok(rel(r.P[2], 5e6 - fric - 6000 - rho * G * 25) < 1e-12);
  assert.ok(r.P[1] > r.P[0] - fric * 0.4 - 6000 * 0.4); // downhill leg gains static head
  assert.ok(rel(minorLoss({ K: 3, rho, v: 2, D, f }).Leq, 3 * D / f) < 1e-12);
  assert.ok(rel(reynolds({ rho, v: 2, D, mu }), 4e5) < 1e-12);
});

test('Weymouth equation — hand calculation in oilfield units and equivalence with the general equation', () => {
  // 15.5 in ID, 50 miles, P1 = 1000 psia, P2 = 800 psia, γ = 0.6, T = 520 °R, Z = 0.9, Tb = 520 °R, Pb = 14.696 psia
  // Q = 433.5 (Tb/Pb) [(P1² − P2²)/(γ T L Z)]^0.5 D^2.667
  const Tb = 519.67;
  const Pb = 101325 / PSI;
  const Q = 433.5 * (Tb / Pb) * Math.sqrt((1000 ** 2 - 800 ** 2) / (0.6 * 520 * 50 * 0.9)) * 15.5 ** 2.667; // scf/d
  const a = { P1: 1000 * PSI, P2: 800 * PSI, D: 15.5 * 0.0254, L: 50 * 1609.344, T: 520 / 1.8, Z: 0.9, gammaG: 0.6 };
  const w = gasTransmission({ equation: 'weymouth', ...a });
  assert.ok(rel(w.qStd, Q * 0.028316846592 / 86400) < 1e-9);
  assert.ok(Q > 110e6 && Q < 120e6); // ≈ 116 MMscfd
  // the general flow equation with the Weymouth friction factor f = 0.032/D^(1/3) gives the same flow (the
  // published constant 433.5 and exponent 2.667 are rounded: agreement to 0.2 %)
  const g = gasFlowGeneral({ ...a, equation: 'weymouth' });
  assert.ok(rel(g.qStd, w.qStd) < 2e-3);
  assert.ok(rel(g.f, 0.032 / 15.5 ** (1 / 3)) < 1e-12);
  // inverse: outlet pressure for the computed flow
  assert.ok(rel(gasTransmission({ equation: 'weymouth', ...a, P2: undefined, qStd: w.qStd }).P2, 800 * PSI) < 1e-9);
});

test('Panhandle A and B — hand calculation in oilfield units', () => {
  const Tb = 519.67;
  const Pb = 101325 / PSI;
  const X = (1000 ** 2 - 800 ** 2);
  const a = { P1: 1000 * PSI, P2: 800 * PSI, D: 15.5 * 0.0254, L: 50 * 1609.344, T: 520 / 1.8, Z: 0.9, gammaG: 0.6, E: 0.92 };
  const QA = 435.87 * 0.92 * (Tb / Pb) ** 1.0788 * (X / (0.6 ** 0.8539 * 520 * 50 * 0.9)) ** 0.5394 * 15.5 ** 2.6182;
  const QB = 737 * 0.92 * (Tb / Pb) ** 1.02 * (X / (0.6 ** 0.961 * 520 * 50 * 0.9)) ** 0.51 * 15.5 ** 2.53;
  assert.ok(rel(gasTransmission({ equation: 'panhandleA', ...a }).qStd, QA * 0.028316846592 / 86400) < 1e-9);
  assert.ok(rel(gasTransmission({ equation: 'panhandleB', ...a }).qStd, QB * 0.028316846592 / 86400) < 1e-9);
  const inv = gasTransmission({ equation: 'panhandleB', ...a, P2: undefined, qStd: QB * 0.028316846592 / 86400 });
  assert.ok(rel(inv.P2, 800 * PSI) < 1e-9);
  assert.throws(() => gasTransmission({ equation: 'panhandleA', ...a, P2: undefined, qStd: 1e6 }), RangeError);
});

test('general gas flow equation matches numerical integration of the momentum balance, including elevation', () => {
  // dP/dx = −f ρ v²/(2D) − ρ g dz/dx, ρ = P M/(Z R T), ρ v = G constant → integrate with RK4
  const D = 0.3;
  const L = 60000;
  const dz = 400;
  const T = 300;
  const Z = 0.88;
  const M = 0.65 * M_AIR;
  const f = 0.012;
  const P1 = 9e6;
  const r = gasFlowGeneral({ P1, P2: 6.5e6, D, L, T, Z, M, dz, f });
  const Gm = r.mDot / (Math.PI * D * D / 4);
  const rhs = (P) => {
    const rho = P * M / (Z * R_GAS * T);
    return -f * Gm * Gm / (2 * D * rho) - rho * G * dz / L;
  };
  let P = P1;
  const n = 4000;
  const h = L / n;
  for (let i = 0; i < n; i++) {
    const k1 = rhs(P);
    const k2 = rhs(P + 0.5 * h * k1);
    const k3 = rhs(P + 0.5 * h * k2);
    const k4 = rhs(P + h * k3);
    P += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  }
  assert.ok(rel(P, 6.5e6) < 1e-8, `integrated outlet ${P}`);
  // flow → pressure inverse, and a rising line delivers less than a level one
  assert.ok(rel(gasFlowGeneral({ P1, mDot: r.mDot, D, L, T, Z, M, dz, f }).P2, 6.5e6) < 1e-10);
  assert.ok(gasFlowGeneral({ P1, P2: 6.5e6, D, L, T, Z, M, f }).mDot > r.mDot);
  // AGA fully turbulent: F = 4 log10(3.7 D/ε) = 2/√f_Darcy
  const aga = gasFlowGeneral({ P1, P2: 6.5e6, D, L, T, Z, M, equation: 'aga', rough: 2e-5 });
  assert.ok(rel(2 / Math.sqrt(aga.f), 4 * Math.log10(3.7 * D / 2e-5)) < 1e-12);
  // friction factor from Colebrook when viscosity is given: converged Reynolds number is self-consistent
  const col = gasFlowGeneral({ P1, P2: 6.5e6, D, L, T, Z, M, mu: 1.3e-5, rough: 2e-5 });
  assert.ok(rel(col.f, frictionFactor({ Re: 4 * col.mDot / (Math.PI * D * 1.3e-5), relRough: 2e-5 / D })) < 1e-10);
});

test('marching solver: adiabatic liquid line reproduces Darcy–Weisbach and frictional heating ΔT = ΔP/(ρ cp)', () => {
  const D = 0.25;
  const rho = 850;
  const mu = 4e-3;
  const cp = 2000;
  const mDot = 80;
  const prof = [{ x: 0, z: 0 }, { x: 5000, z: 0 }, { x: 12000, z: 0 }];
  const r = lineProfile({ fluid: { type: 'liquid', rho, mu, cp }, mDot, Pin: 8e6, Tin: 320, D, profile: prof, U: 0, Tamb: 277 });
  const ref = liquidLinePressureDrop({ Q: mDot / rho, D, rho, mu, L: 12000, method: 'churchill' });
  assert.ok(rel(r.dP, ref.dP) < 1e-10);
  assert.ok(rel(r.Tout - 320, ref.dP / (rho * cp)) < 1e-9);
  // uphill: static head added, potential-energy cooling g Δz/cp removed from the frictional heating
  const up = lineProfile({ fluid: { type: 'liquid', rho, mu, cp }, mDot, Pin: 8e6, Tin: 320, D,
    profile: [{ x: 0, z: 0 }, { x: 12000, z: 300 }], U: 0, Tamb: 277 });
  const len = Math.hypot(12000, 300);
  const dPu = ref.gradient * len + rho * G * 300;
  assert.ok(rel(up.dP, dPu) < 1e-10);
  assert.ok(rel(up.Tout - 320, dPu / (rho * cp) - G * 300 / cp) < 1e-9);
});

test('marching solver: gas line agrees with the isothermal general equation and cools by heat loss + Joule–Thomson', () => {
  const D = 0.4;
  const L = 40000;
  const gammaG = 0.65;
  const mDot = 60;
  const Pin = 1.0e7;
  // (a) near-isothermal: no heat exchange and μJT forced to zero → compare with the average-Z closed form
  const iso = lineProfile({ fluid: { type: 'gas', gammaG, cp: 2300, muJT: 0 }, mDot, Pin, Tin: 300, D, rough: 2e-5,
    profile: [{ x: 0, z: 0 }, { x: L, z: 0 }], U: 0, Tamb: 300, nSub: 200 });
  assert.ok(Math.abs(iso.Tout - 300) < 1e-9);
  const Pavg = (2 / 3) * (Pin + iso.Pout - Pin * iso.Pout / (Pin + iso.Pout));
  const Zavg = zFactor({ P: Pavg, T: 300, gammaG });
  // back out the mean friction factor from the marching result and check the closed form with it
  const Mg = gammaG * M_AIR;
  const A = Math.PI * D * D / 4;
  const fEff = (Pin ** 2 - iso.Pout ** 2) * Mg * D * A * A / (mDot ** 2 * Zavg * R_GAS * 300 * L);
  assert.ok(fEff > 0.0105 && fEff < 0.0125, `effective f ${fEff}`); // Colebrook at Re ≈ 1.5e7, ε/D = 5e-5: 0.0109
  const closed = gasFlowGeneral({ P1: Pin, mDot, D, L, T: 300, Z: Zavg, M: Mg, mu: 1.35e-5, rough: 2e-5 });
  assert.ok(rel(closed.P2, iso.Pout) < 0.01, `closed ${closed.P2} vs marched ${iso.Pout}`);
  // (b) with heat loss the gas approaches and (through Joule–Thomson cooling) can fall below ambient
  const cold = lineProfile({ fluid: { type: 'gas', gammaG, cp: 2300 }, mDot, Pin, Tin: 330, D, rough: 2e-5,
    profile: [{ x: 0, z: 0 }, { x: L / 2, z: 0 }, { x: L, z: 0 }], U: 3, Tamb: 278, nSub: 50 });
  assert.ok(cold.T[1] < 330 && cold.T[2] < cold.T[1]);
  // without the Joule–Thomson term the arrival temperature is the pure exponential; JT cooling lowers it by less
  // than μJT·ΔP (part of the cooling is recovered from the surroundings)
  const noJT = 278 + (330 - 278) * Math.exp(-L * 3 * Math.PI * D / (mDot * 2300));
  assert.ok(cold.Tout < noJT - 0.3, `${cold.Tout} vs ${noJT}`);
  assert.ok(noJT - cold.Tout < 6e-6 * cold.dP);
  const jt = gasJouleThomson({ P: 8e6, T: 300, gammaG, cp: 2300 });
  assert.ok(jt > 2e-6 && jt < 6e-6, `μJT = ${jt} K/Pa`);
  // continuity: ρ v is constant along the line
  assert.ok(rel(cold.rho[0] * cold.v[0], cold.rho[2] * cold.v[2]) < 1e-12);
  assert.throws(() => lineProfile({ fluid: { type: 'gas', gammaG, cp: 2300 }, mDot: 5000, Pin: 2e6, Tin: 300, D: 0.1,
    profile: [{ x: 0, z: 0 }, { x: 5e4, z: 0 }], U: 0, Tamb: 300 }), /Error/);
});

test('line sizing, velocity checks, pump and compressor power — hand calculations', () => {
  const s = sizeLine({ Q: 0.1, rho: 850, mu: 5e-3, vMax: 3, gradMax: 100, candidates: [0.15, 0.2, 0.25, 0.3] });
  assert.ok(rel(s.Dvelocity, Math.sqrt(4 * 0.1 / (Math.PI * 3))) < 1e-12);
  // at the gradient-governed diameter the Darcy–Weisbach gradient equals the limit
  const v = 0.1 / (Math.PI * s.Dgradient ** 2 / 4);
  const f = frictionFactor({ Re: 850 * v * s.Dgradient / 5e-3, relRough: 4.5e-5 / s.Dgradient });
  assert.ok(rel(f * 850 * v * v / (2 * s.Dgradient), 100) < 1e-8);
  assert.equal(s.governing, 'gradient');
  assert.equal(s.selected, 0.3);
  // noise limit 175 (1/ρ)^0.43 capped at 60 m/s; erosional velocity from API RP 14E
  assert.ok(rel(noiseVelocityLimit({ rho: 50 }), 175 * 50 ** -0.43) < 1e-12);
  assert.equal(noiseVelocityLimit({ rho: 1 }), 60);
  const vc = velocityChecks({ v: 10, rho: 100, gas: true });
  assert.ok(rel(vc.erosional, 12.199) < 1e-3 && vc.ok);
  // pump: 0.05 m³/s against 2 MPa at 75 % → 133.3 kW
  assert.ok(rel(pumpPower({ Q: 0.05, dP: 2e6, efficiency: 0.75 }).shaft, 0.05 * 2e6 / 0.75) < 1e-12);
  assert.ok(rel(pumpPower({ Q: 0.05, head: 100, rho: 1000 }).hydraulic, 0.05 * 1000 * G * 100) < 1e-12);
  // compressor: k = 1.3, ηp = 0.8 → (n−1)/n = 0.3/(1.3·0.8) = 0.28846; r = 3; T1 = 300 K
  const c = compressorPolytropic({ mDot: 10, P1: 2e6, P2: 6e6, T1: 300, M: 0.018, k: 1.3, Z: 0.95, polytropicEfficiency: 0.8, mechanicalEfficiency: 1 });
  const e = 0.3 / (1.3 * 0.8);
  assert.ok(rel(c.T2, 300 * 3 ** e) < 1e-12);
  assert.ok(rel(c.head, 0.95 * R_GAS * 300 / 0.018 * (3 ** e - 1) / e) < 1e-12);
  assert.ok(rel(c.shaftPower, 10 * c.head / 0.8) < 1e-12);
  assert.ok(rel(c.n, 1 / (1 - e)) < 1e-12);
});

test('Joukowsky surge and elastic wave speed — hand calculations', () => {
  // ΔP = ρ a Δv = 1000·1200·2 = 2.4 MPa
  const j = joukowsky({ rho: 1000, a: 1200, dv: 2, L: 3000, tClose: 2 });
  assert.equal(j.dP, 2.4e6);
  assert.ok(rel(j.dHead, 1200 * 2 / G) < 1e-12);
  assert.equal(j.period, 5);
  assert.equal(j.rapid, true);
  assert.ok(rel(joukowsky({ rho: 1000, a: 1200, dv: 2, L: 3000, tClose: 20 }).dPslow, 2.4e6 * 5 / 20) < 1e-12);
  // water in steel pipe, D/e = 50, anchored throughout: a = √(K/ρ)/√(1 + (K/E)(D/e)(1 − ν²))
  const a = waveSpeed({ K: 2.2e9, rho: 1000, D: 0.5, e: 0.01, E: 2.07e11, nu: 0.3 });
  assert.ok(rel(a, Math.sqrt(2.2e6) / Math.sqrt(1 + (2.2e9 / 2.07e11) * 50 * 0.91)) < 1e-12);
  assert.ok(a > 1200 && a < 1240);
  // rigid-pipe limit → acoustic speed of the liquid
  assert.ok(rel(waveSpeed({ K: 2.2e9, rho: 1000, D: 0.5, e: 0.01, E: 1e15 }), Math.sqrt(2.2e6)) < 1e-4);
});

test('MOC: frictionless instantaneous closure gives the Joukowsky head and a 4L/a period', () => {
  const L = 1000;
  const a = 1000;
  const D = 0.5;
  const v0 = 1;
  const H0 = 100;
  const r = waterHammerMOC({ L, D, a, H0, Q0: v0 * Math.PI * D * D / 4, tEnd: 4 * 4 * L / a, n: 20 });
  const dH = a * v0 / G;
  assert.ok(Math.abs(r.HmaxAll - (H0 + dH)) < 1e-9);
  assert.ok(Math.abs(r.HminAll - (H0 - dH)) < 1e-9);
  // square wave at the valve: +ΔH for 0 < t < 2L/a, −ΔH for 2L/a < t < 4L/a, repeating with period 4L/a
  const at = (t) => r.Hend[Math.round(t / r.dt)];
  const T = 4 * L / a;
  for (const k of [0, 1, 2]) {
    assert.ok(Math.abs(at(k * T + 0.25 * T) - (H0 + dH)) < 1e-9, `first half-period ${k}`);
    assert.ok(Math.abs(at(k * T + 0.75 * T) - (H0 - dH)) < 1e-9, `second half-period ${k}`);
  }
  // flow at the reservoir reverses after L/a < t < 3L/a
  assert.ok(r.Qin[Math.round(2 * L / a / r.dt)] < 0);
  assert.equal(r.Qend[r.Qend.length - 1], 0);
});

test('MOC: steady state is preserved with friction, slow closure reduces the surge, dead end doubles a pressure step', () => {
  const L = 2000;
  const a = 1200;
  const D = 0.3;
  const Q0 = 0.1;
  const base = { L, D, a, f: 0.02, H0: 150, Q0, n: 40 };
  const steady = waterHammerMOC({ ...base, tau: () => 1, tEnd: 20 });
  const A = Math.PI * D * D / 4;
  const hf = 0.02 * (L / D) * (Q0 / A) ** 2 / (2 * G);
  assert.ok(Math.abs(steady.Hend[steady.Hend.length - 1] - (150 - hf)) < 1e-6);
  assert.ok(rel(steady.Qend[steady.Qend.length - 1], Q0) < 1e-9);
  const fast = waterHammerMOC({ ...base, tClose: 0, tEnd: 30 });
  const slow = waterHammerMOC({ ...base, tClose: 10 * 2 * L / a, tEnd: 80 });
  const jouk = a * (Q0 / A) / G;
  assert.ok(fast.HmaxAll - (150 - hf) >= jouk - 1e-6, 'instantaneous closure reaches at least the Joukowsky rise');
  assert.ok(fast.HmaxAll - (150 - hf) < jouk + hf + 1e-6, 'line packing adds at most the friction head');
  assert.ok(slow.HmaxAll - 150 < 0.35 * jouk, 'slow closure');
  // long after closure the line settles at the reservoir head (friction damps the oscillation)
  const late = waterHammerMOC({ ...base, tClose: 0, tEnd: 600 });
  const tail = late.Hend.slice(Math.floor(0.95 * late.Hend.length)).map((h) => Math.abs(h - 150));
  const head = late.Hend.slice(0, Math.floor(0.05 * late.Hend.length)).map((h) => Math.abs(h - 150));
  assert.ok(Math.max(...tail) < 0.5 * Math.max(...head), `${Math.max(...tail)} vs ${Math.max(...head)}`);
  // dead end: a head step ΔH at the reservoir reflects as 2ΔH at the closed end (frictionless)
  const de = waterHammerMOC({ L, D, a, H0: (t) => (t > 0 ? 120 : 100), Hinit: 100, downstream: 'deadEnd', n: 40, tEnd: 1.5 * L / a });
  assert.ok(Math.abs(de.Hend[de.Hend.length - 1] - 140) < 1e-9);
  assert.throws(() => waterHammerMOC({ ...base, downstream: 'pump', tEnd: 1 }), RangeError);
});
