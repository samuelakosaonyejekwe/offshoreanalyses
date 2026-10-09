import test from 'node:test';
import assert from 'node:assert/strict';
import * as V from '../../src/engine/lines/viv.js';
import { solveRodStatic, rodModes } from '../../src/engine/lines/rod.js';
import { besselJ0, besselY0 } from '../../src/engine/math/special.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('Strouhal curve, shedding frequency and reduced velocity: hand values', () => {
  assert.ok(rel(V.strouhalNumber(1e4), 0.198 * (1 - 19.7 / 1e4)) < 1e-12);
  assert.ok(rel(V.strouhalNumber(100), 0.212 * (1 - 21.2 / 100)) < 1e-12);
  assert.equal(V.strouhalNumber(30), 0);
  assert.equal(V.strouhalNumber(1e6, { stHighRe: 0.25 }), 0.25);
  // continuity of the blend at Re = 150 and 250
  assert.ok(Math.abs(V.strouhalNumber(150) - V.strouhalNumber(149.999)) < 1e-6);
  assert.ok(Math.abs(V.strouhalNumber(250) - V.strouhalNumber(249.999)) < 1e-6);
  assert.ok(rel(V.sheddingFrequency({ U: 1.5, D: 0.3, St: 0.2 }), 1.0) < 1e-12);
  assert.ok(rel(V.reducedVelocity({ U: 1.5, fn: 0.5, D: 0.3 }), 10) < 1e-12);
  const li = V.lockIn({ Vr: 4, Ks: 0.5 });
  assert.ok(li.inline && li.crossflow && !li.crossflowPeak);
  assert.ok(!V.lockIn({ Vr: 4, Ks: 2.5 }).inline);
  // Ks = 4π me ζ/(ρD²); SG = 2π St² Ks
  const sp = V.stabilityParameter({ me: 200, zeta: 0.01, rho: 1025, D: 0.3 });
  assert.ok(rel(sp.Ks, 4 * Math.PI * 200 * 0.01 / (1025 * 0.09)) < 1e-12 && rel(sp.SG, 2 * Math.PI * 0.04 * sp.Ks) < 1e-12);
  assert.ok(rel(V.inlineOnsetF105({ Ks: 0.2 }), 1 / 1.1) < 1e-12);
  assert.ok(rel(V.inlineOnsetF105({ Ks: 1.0 }), 1.6 / 1.1) < 1e-12);
  assert.ok(rel(V.inlineOnsetF105({ Ks: 3 }), 2.2 / 1.1) < 1e-12);
  assert.equal(V.keuleganCarpenter({ Um: 2, T: 10, D: 0.5 }), 40);
});

test('response amplitude curves and drag amplification: limiting and hand values', () => {
  assert.ok(rel(V.responseAmplitude({ SG: 0 }), 1.29) < 1e-12);
  assert.ok(rel(V.responseAmplitude({ SG: 1, gamma: 1.155 }), 1.29 * 1.155 / 1.43 ** 3.35) < 1e-12);
  assert.ok(rel(V.responseAmplitude({ SG: 0, model: 'sarpkaya' }), 0.32 / Math.sqrt(0.06)) < 1e-12);
  assert.equal(V.responseAmplitude({ SG: 0.3, override: 0.9 }), 0.9);
  // Vandiver: 2 Arms/D = 1 → factor 2.043
  assert.ok(rel(V.dragAmplification({ Cd0: 1.2, aRmsOverD: 0.5 }).Cd, 1.2 * 2.043) < 1e-12);
  assert.ok(rel(V.dragAmplification({ Cd0: 1, aOverD: 1, model: 'blevins' }).factor, 3.1) < 1e-12);
});

test('tensioned beam modes vs closed form fn = (n/2L) sqrt(T/m) sqrt(1 + n²π²EI/(T L²)) (analytic path and Rayleigh–Ritz path)', () => {
  const L = 300, EI = 4e7, m = 250, T = 8e5;
  const fn = (n) => (n / (2 * L)) * Math.sqrt(T / m) * Math.sqrt(1 + n * n * Math.PI ** 2 * EI / (T * L * L));
  const a = V.tensionedBeamModes({ L, EI, m, T, nModes: 6 });
  const r = V.tensionedBeamModes({ L, EI, m, T, nModes: 6, forceRitz: true });
  for (let n = 1; n <= 6; n++) {
    assert.ok(rel(a.modes[n - 1].fHz, fn(n)) < 1e-12);
    assert.ok(rel(r.modes[n - 1].fHz, fn(n)) < 1e-9);
  }
  // curvature of the sine mode at the antinode: (nπ/L)²
  assert.ok(rel(Math.abs(a.modes[0].curvature[100]), (Math.PI / L) ** 2) < 1e-12);
});

test('string with linearly varying tension: Rayleigh–Ritz vs the exact Bessel solution J0(a)Y0(b) − J0(b)Y0(a) = 0 and vs WKB', () => {
  const L = 1000, m = 300, T0 = 4e5, T1 = 2.4e6;
  const md = V.tensionedBeamModes({ L, EI: 0, m, T0, T1, nModes: 5, nBasis: 60 });
  for (let n = 1; n <= 5; n++) {
    const exact = V.hangingStringFrequency({ n, L, m, T0, T1, besselJ0, besselY0 });
    assert.ok(rel(md.modes[n - 1].fHz, exact) < 2e-4, `mode ${n}: ${md.modes[n - 1].fHz} vs ${exact}`);
    assert.ok(rel(md.modes[n - 1].wkbHz, exact) < 0.02);
  }
});

test('modes supplied by the rod solver are interchangeable with the analytic ones in the screening', () => {
  const L = 300, EI = 4e7, m = 250, T = 8e5, D = 0.3;
  const rod = solveRodStatic({ sections: [{ length: L, EA: 5e9, EI, w: 0, m }], nElem: 120, endA: { x: 0, z: 0 }, endB: { z: 0, Fx: T } });
  const rm = rodModes(rod, { nModes: 4 });
  const an = V.tensionedBeamModes({ L, EI, m, T, nModes: 4, nPoints: 121 });
  const args = { U: 0.6, D, mass: m, zeta: 0.003 };
  const s1 = V.vivModalScreening({ s: rm.s, modes: rm.modes, ...args });
  const s2 = V.vivModalScreening({ s: an.s, modes: an.modes, ...args });
  assert.equal(s1.dominant, s2.dominant);
  assert.ok(rel(s1.modes[s1.dominant].aOverD, s2.modes[s2.dominant].aOverD) < 0.01);
  assert.ok(rel(rm.modes[2].fHz, an.modes[2].fHz) < 2e-3);
});

test('uniform current: the locked-in mode has power-in length = L, γ = 2/√3, me = m; sheared current shortens it', () => {
  const L = 300, EI = 4e7, m = 250, T = 8e5, D = 0.3;
  const md = V.tensionedBeamModes({ L, EI, m, T, nModes: 6, nPoints: 601 });
  const U = 6 * md.modes[2].fHz * D;                 // Vr = 6 for mode 3
  const sc = V.vivModalScreening({ s: md.s, modes: md.modes, U, D, mass: m, zeta: 0.003 });
  const m3 = sc.modes[2];
  assert.ok(m3.excited && rel(m3.powerInLength, L) < 1e-9);
  assert.ok(rel(m3.gamma, 2 / Math.sqrt(3)) < 1e-4);
  assert.ok(rel(m3.me, m) < 1e-6);
  assert.ok(rel(m3.zetaEff, 0.003) < 1e-9);          // no hydrodynamic damping when the whole length is power-in
  assert.ok(rel(m3.Ks, 4 * Math.PI * m * 0.003 / (1025 * D * D)) < 1e-6);
  assert.ok(rel(m3.aOverD, 1.29 * m3.gamma / (1 + 0.43 * m3.SG) ** 3.35) < 1e-9);
  // energy balance with full-length power-in: A = ∫½ρDU²Cl|φ| / (ω·2ζω·∫mφ²) = (½ρDU²Cl·2L/π)/(2ζω²·mL/2) (n odd)
  const om = 2 * Math.PI * m3.fHz;
  assert.ok(rel(m3.aOverDEnergy * D, (0.5 * 1025 * D * U * U * 0.7 * 2 * L / Math.PI) / (2 * 0.003 * om * om * m * L / 2)) < 1e-3);
  // linear shear 0 → 2U: only part of the length is inside the window and damping rises
  const sh = V.vivModalScreening({ s: md.s, modes: md.modes, U: (s) => 2 * U * s / L, D, mass: m, zeta: 0.003 });
  const h3 = sh.modes[2];
  // window Vr 4…8 ↔ s/L from 4/12 to 8/12
  assert.ok(rel(h3.powerInLength, L / 3) < 0.01);
  assert.ok(h3.zetaEff > 0.003 && h3.aOverD < m3.aOverD);
});

test('VIV fatigue: hand calculation σ = E·(D/2)·(nπ/L)²·A, damage/yr = fn·31 557 600/N(2σ); scatter weighting and suppression', () => {
  const L = 300, EI = 4e7, m = 250, T = 8e5, D = 0.3, E = 207e9;
  const md = V.tensionedBeamModes({ L, EI, m, T, nModes: 6, nPoints: 601 });
  const U = 6 * md.modes[2].fHz * D;
  const sn = (S) => 1e12 * (S / 1e6) ** -3;
  const base = { s: md.s, modes: md.modes, D, mass: m, zeta: 0.003, stressDiameter: D, snCyclesToFailure: sn, E };
  const sc = V.vivModalScreening({ ...base, U });
  const dm = sc.modes[sc.dominant];
  const fat = V.vivFatigue({ ...base, screening: sc });
  const n = dm.index + 1, sig = E * (D / 2) * (n * Math.PI / L) ** 2 * dm.aOverD * D;
  let smax = 0; for (const v of fat.stressAmplitude) smax = Math.max(smax, v);
  assert.ok(rel(smax, sig) < 1e-4);
  assert.ok(rel(fat.maxDamagePerYear, dm.fHz * 365.25 * 86400 / sn(2 * sig)) < 1e-3);
  assert.ok(rel(fat.lifeYears, 1 / fat.maxDamagePerYear) < 1e-12);
  // current occurring 10 % of the time → ten times the life
  const life = V.vivFatigueLife({ ...base, cases: [{ probability: 0.1, U }, { probability: 0.9, U: 0 }] });
  assert.ok(rel(life.lifeYears, 10 * fat.lifeYears) < 1e-9);
  assert.equal(life.cases[1].dominantMode, -1);
  // strakes with 80 % efficiency: coverage to halve the amplitude
  const cov = V.suppressionCoverage({ ...base, U, efficiency: 0.8, targetAOverD: 0.5 * dm.aOverD });
  assert.ok(cov.achievable && cov.coverageFraction > 0 && cov.coverageFraction <= 1 && cov.aOverD <= 0.5 * dm.aOverD);
  const none = V.suppressionCoverage({ ...base, U, efficiency: 0.8, targetAOverD: 5 });
  assert.equal(none.coveredLength, 0);
});

test('hydrodynamic damping forms and span screening hook: hand values', () => {
  const om = 2, D = 0.3, rho = 1025, nu = 1.19e-6;
  const rsw = (om * Math.PI * rho * D * D / 2) * (2 * Math.SQRT2 / Math.sqrt(om * D * D / nu) + 0.25 * 0.25);
  assert.ok(rel(V.hydrodynamicDamping({ zone: 'still', omega: om, D, aOverD: 0.5 }), rsw) < 1e-12);
  assert.ok(rel(V.hydrodynamicDamping({ zone: 'low', omega: om, D, V: 0.4, aOverD: 0.5 }), rsw + 0.18 * rho * D * 0.4) < 1e-12);
  assert.ok(rel(V.hydrodynamicDamping({ zone: 'high', omega: om, D, V: 2 }), 0.2 * rho * 4 / om) < 1e-12);
  const sp = V.spanVivScreening({ frequenciesHz: [0.5, 2], U: 1.0, D: 0.4 });
  assert.ok(rel(sp[0].Vr, 5) < 1e-12 && sp[0].crossflow && sp[0].crossflowPeak && !sp[0].inline);
  assert.ok(rel(sp[1].Vr, 1.25) < 1e-12 && sp[1].inline && !sp[1].crossflow);
});

test('wake oscillator (Facchinetti et al. 2004, A = 12, ε = 0.3): lock-in peak just above Ur = 1/St = 5 and amplitude vs the harmonic-balance closed form', () => {
  const Ur = []; for (let v = 3; v <= 9.01; v += 0.25) Ur.push(v);
  const w = V.wakeOscillator({ massRatio: 2.0, zeta: 0.001, Ur, cycles: 200, stepsPerCycle: 30 });
  // published qualitative behaviour: resonance of wake and structure near Ur·St = 1 with a lock-in plateau above it
  assert.ok(w.peakUr >= 4.5 && w.peakUr <= 7.5, `peak at Ur = ${w.peakUr}`);
  // parameters: μ = (m* + CM)π/4, γ = CD/(4π St) = 0.796, M = CL0/(16π² St² μ)
  assert.ok(rel(w.mu, 3 * Math.PI / 4) < 1e-12 && rel(w.gamma, 2 / (4 * Math.PI * 0.2)) < 1e-12);
  assert.ok(rel(w.M, 0.3 / (16 * Math.PI ** 2 * 0.04 * w.mu)) < 1e-12);
  // simulated amplitude at Ur = 5 (δ = 1) against the model's own lock-in solution (sinusoidal approximation)
  const at5 = w.aOverD[Ur.findIndex((u) => Math.abs(u - 5) < 1e-9)];
  assert.ok(rel(at5, w.lockInAmplitude) < 0.15, `${at5} vs ${w.lockInAmplitude}`);
  assert.ok(w.peakAOverD > 0.1 && w.peakAOverD < 1.5, `A/D = ${w.peakAOverD}`);
  // away from lock-in the response collapses and the lift returns to its fixed-cylinder level (q ≈ 2)
  assert.ok(w.aOverD[0] < 0.3 * w.peakAOverD);
  assert.ok(Math.abs(w.liftMagnification[0] - 1) < 0.15);
  const k = Ur.findIndex((u) => Math.abs(u - 5.5) < 1e-9);
  assert.ok(w.liftMagnification[k] > 1.3);            // lift magnification inside lock-in
  assert.ok(Math.abs(w.frequencyRatio[k] - 1) < 0.25);  // response close to the structural frequency
});
