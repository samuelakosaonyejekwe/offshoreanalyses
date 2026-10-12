import test from 'node:test';
import assert from 'node:assert/strict';
import * as B from '../../src/engine/aero/blade.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const uniform = (L = 40, EI = 4e8, m = 200, n = 21) => ({ x: Array.from({ length: n }, (_, i) => L * i / (n - 1)), EI: new Array(n).fill(EI), massPerLength: new Array(n).fill(m) });

test('uniform cantilever at standstill: f_n = (β_n L)²/(2π) √(EI/(m L⁴)) with β L = 1.8751, 4.6941, 7.8548', () => {
  const u = uniform();
  const r = B.bladeFrequencies({ ...u, nModes: 3 });
  const base = Math.sqrt(4e8 / (200 * 40 ** 4)) / (2 * Math.PI);
  rel(r.frequencies[0], 1.875104 ** 2 * base, 1e-5); rel(r.frequencies[1], 4.694091 ** 2 * base, 1e-4); rel(r.frequencies[2], 7.854757 ** 2 * base, 5e-4);
  rel(r.omegaN[0], 2 * Math.PI * r.frequencies[0], 1e-12);
  // mode shape normalised to unit tip deflection, clamped at the root
  close(r.shapes[0][r.shapes[0].length - 1], 1, 1e-9); close(r.shapes[0][0], 0, 1e-12);
  // a tip mass equal to the blade mass: first frequency from the Rayleigh estimate with m_eff = 0.2427 mL → f ≈ f0 √(0.2427/1.2427)
  const tm = B.bladeFrequencies({ ...u, tipMass: 200 * 40, nModes: 1 }).frequencies[0];
  rel(tm, r.frequencies[0] * Math.sqrt(0.2427 / 1.2427), 0.01);
  // very large shear stiffness reproduces Euler–Bernoulli; a finite one lowers the frequency
  rel(B.bladeFrequencies({ ...u, GAs: new Array(21).fill(1e14), nModes: 1 }).frequencies[0], r.frequencies[0], 1e-4);
  assert.ok(B.bladeFrequencies({ ...u, GAs: new Array(21).fill(2e7), nModes: 1 }).frequencies[0] < r.frequencies[0]);
});

test('rotating uniform blade: Southwell coefficient 1.19 for flap, one less for edge (ω_edge² = ω_flap² − Ω²)', () => {
  const u = uniform();
  const w0 = B.bladeFrequencies({ ...u, nModes: 1 }).omegaN[0];
  const Om = 2.0;
  const flap = B.bladeFrequencies({ ...u, omega: Om, nModes: 1 }).omegaN[0];
  const edge = B.bladeFrequencies({ ...u, omega: Om, direction: 'edge', nModes: 1 }).omegaN[0];
  // Southwell: ω² = ω0² + K Ω², K ≈ 1.19 for the first flap mode of a uniform cantilever without hub offset
  close((flap * flap - w0 * w0) / (Om * Om), 1.19, 0.02);
  rel(edge * edge, flap * flap - Om * Om, 1e-7);
  // a hub offset increases the centrifugal tension and the stiffening
  assert.ok(B.bladeFrequencies({ ...u, omega: Om, hubRadius: 3, nModes: 1 }).omegaN[0] > flap);
  const c = B.campbellData({ ...u, EIflap: u.EI, EIedge: u.EI.map((v) => 4 * v), omegas: [0, 0.5, 1, 1.5, 2], nModes: 2, harmonics: [1, 3] });
  assert.equal(c.omega.length, 5); assert.equal(c.flap.length, 2); assert.equal(c.flap[0].length, 5);
  close(c.southwellFlap[0], 1.19, 0.02); close(c.southwellEdge[0], 0.19, 0.02);
  rel(c.edge[0][0], 2 * c.flap[0][0], 1e-6); // four times the stiffness, twice the frequency
  rel(c.excitation[1].f[4], 3 * 2 / (2 * Math.PI), 1e-12); assert.equal(c.excitation[0].n, 1);
  // 3P at 2 rad/s is 0.955 Hz: the first flap mode (about 0.49 Hz at standstill) is crossed by 3P near Ω = 1.05 rad/s
  const x3 = c.crossings.find((k) => k.harmonic === 3 && /flap/.test(k.mode));
  assert.ok(x3 && x3.omega > 0.9 && x3.omega < 1.2, `3P crossing ${JSON.stringify(c.crossings)}`);
  rel(x3.f, 3 * x3.omega / (2 * Math.PI), 1e-6);
});

test('static deflection of the cantilever: uniform load tip = qL⁴/(8EI), root moment qL²/2; triangular load', () => {
  const u = uniform();
  const r = B.bladeStaticDeflection({ ...u, load: new Array(21).fill(1000) });
  rel(r.tipDeflection, 1000 * 40 ** 4 / (8 * 4e8), 1e-6); rel(r.rootMoment, 1000 * 1600 / 2, 1e-9); rel(r.rootShear, 4e4, 1e-9);
  rel(r.slope[r.slope.length - 1], 1000 * 40 ** 3 / (6 * 4e8), 1e-6);
  // moment at mid-span: q (L/2)²/2
  rel(r.moment[(r.moment.length - 1) / 2], 1000 * 400 / 2, 1e-6);
  // load rising linearly from 0 at the root to q at the tip: tip = 11 q L⁴/(120 EI), root moment q L²/3
  const t = B.bladeStaticDeflection({ ...u, load: u.x.map((x) => 1000 * x / 40) });
  rel(t.tipDeflection, 11 * 1000 * 40 ** 4 / (120 * 4e8), 1e-4); rel(t.rootMoment, 1000 * 1600 / 3, 1e-6);
  // centrifugal stiffening reduces the deflection but not the statically determinate root moment
  const s = B.bladeStaticDeflection({ ...u, load: new Array(21).fill(1000), omega: 2 });
  assert.ok(s.tipDeflection < r.tipDeflection); rel(s.rootMoment, r.rootMoment, 1e-9);
});

test('tip clearance, bending strain, mass scaling, tip speed and gravity cycles', () => {
  // 5 cos 5° + 63 sin 7.5° = 4.981 + 8.223 = 13.204 m to the tower axis; tower radius 2 m; deflection 6 m
  const c = B.tipClearance({ tipDeflection: 6, overhang: 5, tiltDeg: 5, preconeDeg: 2.5, tipRadius: 63, towerRadius: 2 });
  const und = 5 * Math.cos(5 * Math.PI / 180) + 63 * Math.sin(7.5 * Math.PI / 180);
  rel(c.undeflected, und - 2, 1e-9); rel(c.clearance, und - 2 - 6, 1e-9); rel(c.usedFraction, 6 / (und - 2), 1e-9); assert.equal(c.ok, true);
  assert.equal(B.tipClearance({ tipDeflection: 12, overhang: 5, tiltDeg: 5, preconeDeg: 2.5, tipRadius: 63, towerRadius: 2 }).ok, false);
  rel(B.tipClearance({ tipDeflection: 6, overhang: 5, tiltDeg: 5, preconeDeg: 2.5, tipRadius: 63, towerRadius: 2, prebend: 3 }).clearance, und + 3 - 8, 1e-9);
  close(B.bendingStrain(2e6, 4e8, 0.5), 0.0025, 1e-15);
  rel(B.bladeMassScaling({ massRef: 17740, radiusRef: 63, radius: 126 }), 8 * 17740, 1e-12);
  rel(B.bladeMassScaling({ massRef: 17740, radiusRef: 63, radius: 126, exponent: 2.2 }), 17740 * 2 ** 2.2, 1e-12);
  const ts = B.tipSpeedCheck({ omega: 1.267, radius: 63, limit: 90 });
  rel(ts.tipSpeed, 79.82, 1e-3); assert.equal(ts.ok, true); assert.equal(B.tipSpeedCheck({ omega: 1.267, radius: 63 }).ok, null);
  assert.equal(B.tipSpeedCheck({ omega: 1.6, radius: 63, limit: 90 }).ok, false);
  // 17.74 t blade, cg 20.5 m from the root: range 2 m g r = 7.13 MN·m; 1 rad/s for 25 years at 95 %: 1.193e8 revolutions
  const g = B.gravityEdgewiseCycles({ bladeMass: 17740, cgRadius: 20.5, meanRotorSpeed: 1 });
  rel(g.momentRange, 2 * 17740 * 9.80665 * 20.5, 1e-9); rel(g.cycles, 25 * 365.25 * 86400 * 0.95 / (2 * Math.PI), 1e-3);
});

test('generic blade distributions hit the requested mass and first frequencies', () => {
  const d = B.genericBladeDistributions({ length: 61.5, mass: 17740, flapFrequency: 0.70, edgeFrequency: 1.08 });
  assert.equal(d.x.length, 25); close(d.x[24], 61.5, 1e-9); assert.match(d.label, /generic/i);
  let m = 0; for (let i = 1; i < 25; i++) m += 0.5 * (d.massPerLength[i] + d.massPerLength[i - 1]) * (d.x[i] - d.x[i - 1]);
  rel(m, 17740, 1e-6); rel(d.massPerLength[24] / d.massPerLength[0], 0.1, 1e-9); rel(d.EIflap[24] / d.EIflap[0], 0.002, 1e-9);
  rel(B.bladeFrequencies({ x: d.x, EI: d.EIflap, massPerLength: d.massPerLength, nModes: 1 }).frequencies[0], 0.70, 5e-3);
  rel(B.bladeFrequencies({ x: d.x, EI: d.EIedge, massPerLength: d.massPerLength, nModes: 1 }).frequencies[0], 1.08, 5e-3);
});

test('classical laminate theory: unidirectional ply, symmetric cross-ply and ±45 lay-up', () => {
  const ply = { E1: 40e9, E2: 10e9, G12: 4e9, nu12: 0.3, thickness: 1e-3, density: 1900 };
  // Q11 = E1/(1 − ν12 ν21), ν21 = ν12 E2/E1 = 0.075 → Q11 = 40.92 GPa, Q22 = 10.23 GPa, Q12 = 3.069 GPa
  const den = 1 - 0.3 * 0.075; const Q11 = 40e9 / den; const Q22 = 10e9 / den; const Q12 = 0.3 * 10e9 / den;
  const ud = B.laminateABD([{ ply, angleDeg: 0 }]);
  rel(ud.A[0][0], Q11 * 1e-3, 1e-12); rel(ud.A[1][1], Q22 * 1e-3, 1e-12); rel(ud.A[0][1], Q12 * 1e-3, 1e-12); rel(ud.A[2][2], 4e9 * 1e-3, 1e-12);
  rel(ud.D[0][0], Q11 * 1e-9 / 12, 1e-12); close(ud.B[0][0], 0, 1e-6);
  rel(ud.Ex, 40e9, 1e-9); rel(ud.Ey, 10e9, 1e-9); rel(ud.Gxy, 4e9, 1e-9); rel(ud.nuxy, 0.3, 1e-9); rel(ud.arealMass, 1.9, 1e-12);
  // 90° ply swaps the axes
  rel(B.laminateABD([{ ply, angleDeg: 90 }]).Ex, 10e9, 1e-9);
  // [0/90]s: A11 = A22 = (Q11 + Q22) h/2, no coupling; outer 0° plies make D11 > D22
  const cp = B.laminateABD([{ ply, angleDeg: 0 }, { ply, angleDeg: 90 }, { ply, angleDeg: 90 }, { ply, angleDeg: 0 }]);
  rel(cp.A[0][0], (Q11 + Q22) * 2e-3, 1e-12); rel(cp.A[1][1], cp.A[0][0], 1e-12); close(cp.B[0][0] / cp.A[0][0], 0, 1e-12);
  // D11 = (1/3)[Q22 × 2 × 1e-9 + Q11 × 2 × (8e-9 − 1e-9)]
  rel(cp.D[0][0], (Q22 * 2e-9 + Q11 * 14e-9) / 3, 1e-12); assert.ok(cp.D[0][0] > cp.D[1][1]);
  close(cp.thickness, 4e-3, 1e-15); assert.deepEqual(cp.z.map((v) => Math.round(v * 1e6)), [-2000, -1000, 0, 1000, 2000]);
  // an unsymmetric [0/90] laminate couples extension and bending: B11 = ½ (Q22 − Q11) t²... with 0° at the bottom
  const us = B.laminateABD([{ ply, angleDeg: 0 }, { ply, angleDeg: 90 }]);
  rel(us.B[0][0], 0.5 * (Q22 - Q11) * 1e-6, 1e-9);
  // ±45: Gxy = (Q11 + Q22 − 2 Q12)/4, far above G12
  const pm = B.laminateABD([{ ply, angleDeg: 45 }, { ply, angleDeg: -45 }, { ply, angleDeg: -45 }, { ply, angleDeg: 45 }]);
  rel(pm.Gxy, (Q11 + Q22 - 2 * Q12) / 4, 1e-9); close(pm.A[0][2] / pm.A[0][0], 0, 1e-12);
  // generic lay-ups build symmetric stacks
  const tri = B.genericLayup('triax');
  assert.equal(tri.length % 2, 0); assert.equal(tri[0].angleDeg, tri[tri.length - 1].angleDeg);
  close(B.laminateABD(tri).B[0][0] / B.laminateABD(tri).A[0][0], 0, 1e-9);
  assert.ok(B.GENERIC_PLIES.glassEpoxyUD.E1 > B.GENERIC_PLIES.glassEpoxyUD.E2);
});

test('first-ply failure: a unidirectional ply fails at Xt, Xc and Yt under uniaxial load by both criteria', () => {
  const ply = { E1: 40e9, E2: 10e9, G12: 4e9, nu12: 0.3, thickness: 1e-3, Xt: 1000e6, Xc: 600e6, Yt: 40e6, Yc: 140e6, S: 60e6 };
  const lay = [{ ply, angleDeg: 0 }];
  for (const criterion of ['max-stress', 'tsai-wu']) {
    // Nx = 1e5 N/m → σ1 = 100 MPa → R = 10 in tension; −100 MPa → R = 6 in compression; Ny = 1e4 → σ2 = 10 MPa → R = 4
    rel(B.firstPlyFailure(lay, { N: [1e5, 0, 0] }, { criterion }).strengthRatio, 10, 1e-9, criterion);
    rel(B.firstPlyFailure(lay, { N: [-1e5, 0, 0] }, { criterion }).strengthRatio, 6, 1e-9, criterion);
    rel(B.firstPlyFailure(lay, { N: [0, 1e4, 0] }, { criterion }).strengthRatio, 4, 1e-9, criterion);
    rel(B.firstPlyFailure(lay, { N: [0, 0, 1e4] }, { criterion }).strengthRatio, 6, 1e-9, criterion);
  }
  const r = B.firstPlyFailure(lay, { N: [1e5, 0, 0] }, { criterion: 'max-stress' });
  rel(r.failureIndex, 0.1, 1e-9); rel(r.plies[0].sigma1, 100e6, 1e-9); rel(r.strain0[0], 100e6 / 40e9, 1e-9); assert.match(r.mode, /fibre|1|tension/i);
  // cross-ply [0/90]s under Nx: the 90° plies fail first in transverse tension
  const cp = [{ ply, angleDeg: 0 }, { ply, angleDeg: 90 }, { ply, angleDeg: 90 }, { ply, angleDeg: 0 }];
  const f = B.firstPlyFailure(cp, { N: [1e5, 0, 0] }, { criterion: 'max-stress' });
  assert.equal(f.criticalAngleDeg, 90); assert.ok(f.strengthRatio < 10);
  // pure bending of the UD ply: outer-fibre stress 6 M/t²: M = 100 N·m/m → 600 MPa → R = min(1000/600, 600/600) = 1 (compression side)
  rel(B.firstPlyFailure(lay, { N: [0, 0, 0], M: [100, 0, 0] }, { criterion: 'max-stress' }).strengthRatio, 1, 1e-9);
});

test('box spar stiffness, stress ratio, Goodman life and Miner damage', () => {
  // b = 0.6 m, h = 0.4 m; caps E 40 GPa, G 4 GPa, 20 mm; webs E 12 GPa, G 10 GPa, 8 mm
  const s = B.boxSparSection({ width: 0.6, height: 0.4, cap: { E: 40e9, G: 4e9, thickness: 0.02, density: 1900 }, web: { E: 12e9, G: 10e9, thickness: 0.008, density: 1800 } });
  rel(s.EIflap, 2 * 40e9 * 0.02 * 0.6 * 0.04 + 2 * 12e9 * 0.008 * 0.4 ** 3 / 12, 1e-12);
  rel(s.EIedge, 2 * 12e9 * 0.008 * 0.4 * 0.09 + 2 * 40e9 * 0.02 * 0.6 ** 3 / 12, 1e-12);
  rel(s.GJ, 4 * (0.24) ** 2 / (2 * 0.6 / (4e9 * 0.02) + 2 * 0.4 / (10e9 * 0.008)), 1e-12);
  rel(s.EA, 2 * 40e9 * 0.02 * 0.6 + 2 * 12e9 * 0.008 * 0.4, 1e-12); rel(s.massPerLength, 2 * 1900 * 0.02 * 0.6 + 2 * 1800 * 0.008 * 0.4, 1e-12);
  close(s.cFlap, 0.2, 0.011);
  // amplitude 45, mean 55: min 10, max 100 → R = 0.1; fully reversed → −1; compression–compression 10
  close(B.stressRatio(45, 55), 0.1, 1e-12); close(B.stressRatio(30, 0), -1, 1e-12); close(B.stressRatio(45, -55), 10, 1e-12);
  // symmetric strengths 500, fully reversed amplitude 100, m = 10: N = 5^10 = 9.766e6
  const g = B.goodmanCycles({ amplitude: 100, mean: 0, Rt: 500, Rc: 500 });
  rel(g.N, 5 ** 10, 1e-9); close(g.allowableAmplitude, 500, 1e-9); close(g.R, -1, 1e-12);
  // tensile mean 200 with Rt 500, Rc 300: envelope min(300, 500) = 300 → N = 3^10
  rel(B.goodmanCycles({ amplitude: 100, mean: 200, Rt: 500, Rc: 300 }).N, 3 ** 10, 1e-9);
  // compressive mean −100: min(600, 200) = 200 → N = 2^10
  rel(B.goodmanCycles({ amplitude: 100, mean: -100, Rt: 500, Rc: 300 }).N, 1024, 1e-9);
  assert.equal(B.goodmanCycles({ amplitude: 400, mean: 200, Rt: 500, Rc: 300 }).N, 0);
  // Miner: 1e5 cycles of range 200 (amplitude 100) at zero mean, life 9.766e6 → damage 0.01024; scaled ×10
  const d = B.goodmanDamage({ cycles: { ranges: [200], means: [0], counts: [1e5] }, Rt: 500, Rc: 500 });
  rel(d.damage, 1e5 / 5 ** 10, 1e-9); rel(d.life, 5 ** 10 / 1e5, 1e-9);
  rel(B.goodmanDamage({ cycles: { ranges: [200], means: [0], counts: [1e5] }, Rt: 500, Rc: 500, scale: 10 }).damage, 1e6 / 5 ** 10, 1e-9);
  // a partial factor 1.25 on the amplitude raises the damage by 1.25^10
  rel(B.goodmanDamage({ cycles: { ranges: [200], means: [0], counts: [1e5] }, Rt: 500, Rc: 500, gamma: 1.25 }).damage, 1e5 / 4 ** 10, 1e-9);
});

test('tubular section stresses and the classical shell-buckling stress 0.605 E t/R', () => {
  // D = 6 m, t = 60 mm: A = π/4 (36 − 34.5744) = 1.1196 m²; W = π (D⁴ − d⁴)/(32 D) = 1.6468 m³
  const o = { D: 6, t: 0.06, N: 5e6, M: 150e6, V: 2e6, T: 10e6 };
  const r = B.tubularSectionStress(o);
  const d = 5.88; const A = Math.PI / 4 * (36 - d * d); const W = Math.PI * (6 ** 4 - d ** 4) / (32 * 6); const J = Math.PI * (6 ** 4 - d ** 4) / 32;
  rel(r.area, A, 1e-12); rel(r.W, W, 1e-12); rel(r.sigmaAxial, 5e6 / A, 1e-12); rel(r.sigmaBending, 150e6 / W, 1e-12);
  rel(r.sigmaCompMax, 5e6 / A + 150e6 / W, 1e-12); rel(r.tauShear, 2 * 2e6 / A, 1e-12); rel(r.tauTorsion, 10e6 * 3 / J, 1e-12);
  const vmFibre = Math.hypot(r.sigmaCompMax, Math.sqrt(3) * r.tauTorsion); const vmNa = Math.hypot(r.sigmaAxial, Math.sqrt(3) * (r.tauShear + r.tauTorsion));
  rel(r.vonMises, Math.max(vmFibre, vmNa), 1e-9);
  // σ_cr = E t/(R √(3 (1 − ν²))) with the mid-surface radius R = (D − t)/2 = 2.97 m:
  //   210e9 × 0.06/(2.97 × 1.6523) = 2568 MPa; knock-down 0.3 → 770 MPa, capped at fy = 355 MPa
  const b = B.shellBucklingCheck({ D: 6, t: 0.06, sigmaComp: 100e6 });
  rel(b.sigmaCr, 210e9 * 0.06 / (2.97 * Math.sqrt(3 * 0.91)), 1e-12); rel(b.sigmaCr, 0.605 * 210e9 * 0.06 / 2.97, 1e-3);
  rel(b.utilisation, 100e6 / b.sigmaCr, 1e-12);
  const k = B.shellBucklingCheck({ D: 6, t: 0.06, sigmaComp: 100e6, knockdown: 0.3, fy: 355e6 });
  rel(k.capacity, 355e6, 1e-12); rel(k.utilisation, 100 / 355, 1e-12);
  rel(B.shellBucklingCheck({ D: 6, t: 0.06, sigmaComp: 100e6, knockdown: 0.1 }).capacity, 0.1 * b.sigmaCr, 1e-12);
  rel(B.shellBucklingCheck({ D: 6, t: 0.06, sigmaComp: 100e6, capacity: () => 200e6 }).utilisation, 0.5, 1e-12);
});
