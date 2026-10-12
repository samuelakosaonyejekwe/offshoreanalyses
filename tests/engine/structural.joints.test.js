import test from 'node:test';
import assert from 'node:assert/strict';
import { jointGeometry, qBeta, qGap, quApi, qfApi, classifyJoint, apiJointCheck } from '../../src/engine/structural/joints.js';

const near = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('joint geometry parameters and validity range', () => {
  const g = jointGeometry({ D: 1.0, T: 0.025, d: 0.5, t: 0.0125, thetaDeg: 45, gap: 0.1 });
  near(g.beta, 0.5, 1e-12); near(g.gamma, 20, 1e-12); near(g.tau, 0.5, 1e-12); near(g.theta, Math.PI / 4, 1e-12); near(g.gapRatio, 0.1, 1e-12);
  assert.equal(g.warnings.length, 0);
  assert.ok(jointGeometry({ D: 1.0, T: 0.008, d: 0.15, t: 0.006, thetaDeg: 25 }).warnings.length === 3, 'β < 0.2, γ > 50, θ < 30°');
  assert.throws(() => jointGeometry({ D: 1.0, T: 0.02, d: 1.2, t: 0.01, thetaDeg: 45 }), RangeError);
});

test('Qu strength factors, API RP 2A-WSD Table 4.3-1 (hand values for β = 0.5, γ = 20)', () => {
  const beta = 0.5; const gamma = 20;
  // 0.5^1.2 = 0.435275; 0.5^1.6 = 0.329877; 0.5^2.6 = 0.164938
  near(quApi({ type: 'K', load: 'compression', beta, gamma, Qg: 1 }), 40 * 0.435275, 1e-5, 'K: (16 + 1.2γ)β^1.2 = 40 β^1.2 at γ = 20');
  near(quApi({ type: 'K', load: 'tension', beta, gamma: 30, Qg: 1.1 }), 40 * 0.435275 * 1.1, 1e-5, 'K: capped at 40 β^1.2 Qg');
  near(quApi({ type: 'K', load: 'tension', beta, gamma: 12, Qg: 1 }), (16 + 14.4) * 0.435275, 1e-5);
  near(quApi({ type: 'Y', load: 'tension', beta, gamma }), 15, 1e-12);
  near(quApi({ type: 'T', load: 'compression', beta, gamma }), 2.8 + 36 * 0.329877, 1e-5);
  near(quApi({ type: 'Y', load: 'compression', beta, gamma: 10 }), 2.8 + 28 * 0.329877, 1e-5);
  near(quApi({ type: 'Y', load: 'compression', beta, gamma: 40 }), 2.8 + 36 * 0.329877, 1e-5, 'cap 2.8 + 36 β^1.6');
  near(quApi({ type: 'X', load: 'tension', beta, gamma }), 11.5, 1e-12);
  near(quApi({ type: 'X', load: 'tension', beta: 1.0, gamma }), 20.7 + 0.1 * (340 - 220), 1e-12);
  near(quApi({ type: 'X', load: 'compression', beta, gamma }), 2.8 + 14 * 0.5, 1e-12);
  near(quApi({ type: 'X', load: 'compression', beta: 0.8, gamma }), (2.8 + 14 * 0.8) * 0.3 / (0.8 * (1 - 0.6664)), 1e-12);
  near(quApi({ type: 'K', load: 'ipb', beta, gamma }), 19 * 0.435275, 1e-5);
  near(quApi({ type: 'X', load: 'opb', beta, gamma }), 2.5 + 8.5 * 0.164938, 1e-5);
  near(qBeta(0.8), 1.12410, 1e-5);
  assert.equal(qBeta(0.6), 1);
});

test('gap factor Qg', () => {
  // g/D = 0.1: 1 + 0.2 (1 − 0.28)³ = 1.074650
  near(qGap({ gapRatio: 0.1, gamma: 20 }), 1.07465, 1e-5);
  assert.equal(qGap({ gapRatio: 0.5, gamma: 20 }), 1, 'not less than 1.0');
  // overlap g/D = −0.2: 0.13 + 0.65 φ γ^0.5 with φ = τ Fyb/Fy = 0.5
  near(qGap({ gapRatio: -0.2, gamma: 16, tau: 0.5, fyBrace: 355e6, fyChord: 355e6 }), 0.13 + 0.65 * 0.5 * 4, 1e-12);
  // linear interpolation between g/D = −0.05 and +0.05
  const lo = 0.13 + 0.65 * 0.5 * 4; const hi = 1 + 0.2 * (1 - 0.14) ** 3;
  near(qGap({ gapRatio: 0, gamma: 16, tau: 0.5, fyBrace: 355e6, fyChord: 355e6 }), 0.5 * (lo + hi), 1e-12);
});

test('chord load factor Qf, Table 4.3-2 coefficients', () => {
  const Py = 1e7; const Mp = 2e6; const FS = 1.6;
  // T/Y axial, chord compression FS·Pc/Py = −0.4: Qf = 1 + 0.3(−0.4) − 0.8(0.16) = 0.752
  near(qfApi({ type: 'Y', load: 'axial', beta: 0.5, Pc: -0.4 * Py / FS, Py, Mp }).Qf, 0.752, 1e-12);
  // K axial with chord in-plane moment: FS Pc/Py = −0.3, FS M/Mp = 0.2: 1 − 0.06 − 0.04 − 0.3 (0.09 + 0.04) = 0.861
  near(qfApi({ type: 'K', load: 'axial', beta: 0.5, Pc: -0.3 * Py / FS, Mipb: 0.2 * Mp / FS, Py, Mp }).Qf, 0.861, 1e-12);
  // X axial: β ≤ 0.9 → (0.2, 0, 0.5); β = 1.0 → (−0.2, 0, 0.2); β = 0.95 → (0, 0, 0.35)
  near(qfApi({ type: 'X', load: 'axial', beta: 0.7, Pc: -0.5 * Py / FS, Py, Mp }).Qf, 1 - 0.1 - 0.5 * 0.25, 1e-12);
  near(qfApi({ type: 'X', load: 'axial', beta: 1.0, Pc: -0.5 * Py / FS, Py, Mp }).Qf, 1 + 0.1 - 0.2 * 0.25, 1e-12);
  const mid = qfApi({ type: 'X', load: 'axial', beta: 0.95, Pc: -0.5 * Py / FS, Py, Mp });
  near(mid.C1 + 1, 1, 1e-12); near(mid.C3, 0.35, 1e-12);
  // moment loading: (0.2, 0, 0.4)
  near(qfApi({ type: 'K', load: 'moment', beta: 0.5, Pc: -0.5 * Py / FS, Py, Mp }).Qf, 1 - 0.1 - 0.4 * 0.25, 1e-12);
  assert.equal(qfApi({ type: 'Y', load: 'axial', beta: 0.5, Pc: 0, Py, Mp }).Qf, 1);
});

test('joint classification from the load path', () => {
  // balanced K: one brace in tension, one in compression on the same side
  const k = classifyJoint([{ P: 1000, thetaDeg: 45, side: 1 }, { P: -1000, thetaDeg: 45, side: 1 }]);
  near(k[0].K, 1, 1e-12); near(k[1].K, 1, 1e-12); assert.equal(k[0].Y, 0);
  // X: equal compression braces on opposite sides
  const x = classifyJoint([{ P: -800, thetaDeg: 90, side: 1 }, { P: -800, thetaDeg: 90, side: -1 }]);
  near(x[0].X, 1, 1e-12); near(x[1].X, 1, 1e-12);
  // single brace: Y
  assert.deepEqual(classifyJoint([{ P: -500, thetaDeg: 60, side: 1 }])[0], { K: 0, X: 0, Y: 1 });
  // unbalanced K: 1000 against 600 (same angle) → the larger brace is 60 % K and 40 % Y
  const u = classifyJoint([{ P: -1000, thetaDeg: 45, side: 1 }, { P: 600, thetaDeg: 45, side: 1 }]);
  near(u[0].K, 0.6, 1e-12); near(u[0].Y, 0.4, 1e-12); near(u[1].K, 1, 1e-12);
});

test('allowable capacities, interaction and the 50 % rule (hand calculation, T joint)', () => {
  // D = 1.0, T = 0.04, d = 0.5, t = 0.02, θ = 90°, Fy = 345 MPa: β = 0.5, γ = 12.5
  const Fy = 345e6; const T = 0.04; const d = 0.5;
  const QuC = 2.8 + (20 + 0.8 * 12.5) * 0.329877; // = 12.6963
  const Pa = QuC * Fy * T * T / 1.6;
  const r = apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 90, fyChord: Fy, P: -1e6, Mipb: 1e5, Mopb: 5e4 });
  near(r.Qu.Y, QuC, 1e-5);
  near(r.Pa, Pa, 1e-5);
  const MaI = (5 + 0.7 * 12.5) * 0.435275 * Fy * T * T * d / 1.6;
  const MaO = (2.5 + (4.5 + 2.5) * 0.164938) * Fy * T * T * d / 1.6;
  near(r.MaIpb, MaI, 1e-5); near(r.MaOpb, MaO, 1e-5);
  near(r.uc, 1e6 / Pa + (1e5 / MaI) ** 2 + 5e4 / MaO, 1e-5);
  // 50 % rule: brace yield strength π/4 (0.5² − 0.46²) Fy
  const Ab = Math.PI / 4 * (0.25 - 0.46 ** 2);
  near(r.minimumCapacity.required, 0.5 * Ab * Fy, 1e-12);
  near(r.minimumCapacity.provided, Pa * 1.6, 1e-5);
  // chord compression lowers the capacity through Qf; tension load uses Qu = 30β
  const loaded = apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 90, fyChord: Fy, P: -1e6, chord: { P: -1.5e7 } });
  assert.ok(loaded.Pa < r.Pa);
  near(apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 90, fyChord: Fy, P: 1e6 }).Pa, 15 * Fy * T * T / 1.6, 1e-12);
  // inclined brace: capacity ∝ 1/sinθ; mixed 50 % K / 50 % Y classification averages the capacities
  const y = apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 45, fyChord: Fy, P: -1e6 });
  near(y.Pa, Pa / Math.sin(Math.PI / 4), 1e-5);
  const kk = apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 45, fyChord: Fy, P: -1e6, classification: { K: 1 }, gap: 0.1 });
  const mix = apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 45, fyChord: Fy, P: -1e6, classification: { K: 0.5, Y: 0.5 }, gap: 0.1 });
  near(mix.Pa, 0.5 * (kk.Pa + y.Pa), 1e-12);
  near(kk.Qu.K, (16 + 1.2 * 12.5) * 0.435275 * 1.07465, 1e-5);
  // Fyc limited to 0.8 × tensile strength
  near(apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 90, fyChord: Fy, fuChord: 400e6, P: 1e6 }).Pa, 15 * 320e6 * T * T / 1.6, 1e-12);
});

test('overstressed chord (Qf ≤ 0): zero capacity, large finite utilisation and a flag instead of a negative utilisation', () => {
  // T joint as above. Chord yield capacity Py = π/4 (1.0² − 0.92²) × 345e6 = 41.62 MN.
  // Y axial: Qf = 1 + 0.3 (FS Pc/Py) − 0.8 (FS Pc/Py)². With Pc = −38 MN: FS Pc/Py = −1.4608 → Qf = 1 − 0.438 − 1.707 = −1.145
  const Fy = 345e6; const T = 0.04; const d = 0.5;
  const Py = Math.PI / 4 * (1 - 0.92 ** 2) * Fy;
  const x = 1.6 * -38e6 / Py;
  assert.ok(1 + 0.3 * x - 0.8 * x * x < -1);
  const r = apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 90, fyChord: Fy, P: -1e6, Mipb: 1e5, Mopb: 5e4, chord: { P: -38e6 } });
  assert.equal(r.chordOverstressed, true);
  assert.equal(r.Pa, 0); assert.equal(r.MaIpb, 0); assert.equal(r.MaOpb, 0);
  assert.equal(r.Qf.Y, 0); assert.equal(r.Qf.moment, 0);
  assert.equal(r.uc, 1000); assert.equal(r.ucAxial, 1000);
  for (const c of r.checks) assert.ok(Number.isFinite(c.uc) && c.uc >= 1000, `${c.name}: ${c.uc}`);
  assert.equal(r.minimumCapacity.ok, false);
  assert.ok(r.warnings.some((w) => /Qf is not positive/.test(w)));
  // even an unloaded brace is flagged: the chord itself has failed
  const idle = apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 90, fyChord: Fy, chord: { P: -38e6 } });
  assert.equal(idle.chordOverstressed, true); assert.equal(idle.uc, 1000);
  // moderately loaded chord: reduced but positive capacity, ordinary utilisation, no flag
  const ok = apiJointCheck({ D: 1.0, T, d, t: 0.02, thetaDeg: 90, fyChord: Fy, P: -1e6, chord: { P: -1.5e7 } });
  const xo = 1.6 * -1.5e7 / Py;
  assert.equal(ok.chordOverstressed, false);
  near(ok.Qf.Y, 1 + 0.3 * xo - 0.8 * xo * xo, 1e-12);
  assert.ok(ok.Pa > 0 && ok.uc > 0 && ok.uc < 1000);
  near(ok.uc, 1e6 / ok.Pa, 1e-12);
});
