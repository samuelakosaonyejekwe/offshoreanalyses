import test from 'node:test';
import assert from 'node:assert/strict';
import {
  scfParameters, shortChordFactors, scfTY, scfX, scfK, scfKT, hotSpotStresses, scfButtWeldPlate, scfGirthWeld, tyJointHotSpots,
} from '../../src/engine/structural/scf.js';

const near = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

// T joint: D = 1.0 m, T = 40 mm, d = 0.5 m, t = 20 mm, L = 10 m, θ = 90° → β = 0.5, γ = 12.5, τ = 0.5, α = 20
const P = { beta: 0.5, gamma: 12.5, tau: 0.5, alpha: 20, theta: Math.PI / 2 };

test('Efthymiou T-joint SCFs, chord ends fixed — hand calculation from the equations', () => {
  const p = scfParameters({ D: 1.0, T: 0.04, d: 0.5, t: 0.02, L: 10, thetaDeg: 90 });
  near(p.beta, 0.5, 1e-12); near(p.gamma, 12.5, 1e-12); near(p.tau, 0.5, 1e-12); near(p.alpha, 20, 1e-12);
  assert.equal(p.warnings.length, 0);
  const s = scfTY({ ...P, chordEnds: 'fixed' });
  // (1) chord saddle: γ τ^1.1 [1.11 − 3(β − 0.52)²] = 12.5 × 0.46652 × 1.1088 = 6.4660
  near(s.axial.chordSaddle, 6.4660, 1e-4);
  // (2) chord crown: γ^0.2 τ [2.65 + 5(β − 0.65)²] + τβ(0.25α − 3) = 1.65723 × 0.5 × 2.7625 + 0.25 × 2 = 2.78905
  near(s.axial.chordCrown, 2.78905, 1e-4);
  // (3) brace saddle: 1.3 + γ τ^0.52 α^0.1 [0.187 − 1.25 β^1.1 (β − 0.96)] = 1.3 + 12.5 × 0.69737 × 1.34928 × 0.455249 = 6.6546
  near(s.axial.braceSaddle, 6.6546, 1e-4);
  // (4) brace crown: 3 + γ^1.2 [0.12 e^(−4β) + 0.011β² − 0.045] + βτ(0.1α − 1.2) = 3 + 20.7153 × (−0.0260098) + 0.2 = 2.6612
  near(s.axial.braceCrown, 2.6612, 1e-4);
  // (8) IPB chord crown: 1.45 β τ^0.85 γ^(1 − 0.68β) = 1.45 × 0.5 × 0.55479 × 5.29615 = 2.1302
  near(s.ipb.chordCrown, 2.1302, 1e-4);
  // (9) IPB brace crown: 1 + 0.65 β τ^0.4 γ^(1.09 − 0.77β) = 1 + 0.65 × 0.5 × 0.75786 × 5.93365 = 2.4615
  near(s.ipb.braceCrown, 2.4615, 1e-4);
  // (10) OPB chord saddle: γ τ β (1.7 − 1.05β³) = 3.125 × 1.56875 = 4.90234
  near(s.opb.chordSaddle, 4.90234, 1e-5);
  // (11) OPB brace saddle: τ^−0.54 γ^−0.05 (0.99 − 0.47β + 0.08β⁴) × (10) = 1.45397 × 0.88136 × 0.76 × 4.90234 = 4.7745
  near(s.opb.braceSaddle, 4.7745, 1e-4);
});

test('general chord fixity (C = 0.7), inclination and short-chord factors', () => {
  const s = scfTY({ ...P, chordEnds: 'general', C: 0.7 });
  // θ = 90°: sin 2θ = 0 → the C1 term of Eq. (5) vanishes
  near(s.axial.chordSaddle, 6.4660, 1e-4);
  // (6): crown with C2 = 0.35: 2.28905 + 0.25 (0.35 × 20 − 3) = 3.28905; (7): C3 = 0.14: 2.46120 + 0.25 (2.8 − 1.2) = 2.8612
  near(s.axial.chordCrown, 3.28905, 1e-4);
  near(s.axial.braceCrown, 2.8612, 1e-4);
  // C = 1 (fully fixed) gives C2 = 0.5, larger than the fixed-end coefficient 0.25 of Eq. (2) by design of the general formula
  // Y joint at 45°: Eq. (1) scales with (sin θ)^1.6 and Eq. (5) adds C1 (0.8α − 6) τ β² √(1 − β²) sin²2θ
  const y = scfTY({ ...P, theta: Math.PI / 4, chordEnds: 'general', C: 0.7 });
  near(y.axial.chordSaddle, 6.4660 * Math.sin(Math.PI / 4) ** 1.6 + 0.4 * 10 * 0.5 * 0.25 * Math.sqrt(0.75), 1e-4);
  // short chord α = 8: F1 = 1 − (0.83β − 0.56β² − 0.02) γ^0.23 exp(−0.21 γ^−1.16 α^2.5)
  const F = shortChordFactors({ beta: 0.5, gamma: 12.5, alpha: 8 });
  const f1 = 1 - (0.415 - 0.14 - 0.02) * 12.5 ** 0.23 * Math.exp(-0.21 * 12.5 ** -1.16 * 8 ** 2.5);
  near(F.F1, f1, 1e-12);
  assert.ok(F.F1 < 1 && F.F2 < 1 && F.F3 < 1 && F.F1 > 0.5);
  near(F.F3, 1 - 0.55 * 0.5 ** 1.8 * 12.5 ** 0.16 * Math.exp(-0.49 * 12.5 ** -0.89 * 8 ** 1.8), 1e-12);
  const short = scfTY({ ...P, alpha: 8, chordEnds: 'fixed' });
  near(short.axial.chordSaddle, 12.5 * 0.5 ** 1.1 * 1.1088 * f1, 1e-10);
  assert.deepEqual(shortChordFactors({ beta: 0.5, gamma: 12.5, alpha: 12 }), { F1: 1, F2: 1, F3: 1, F4: 1 });
});

test('validity range is reported, values still returned', () => {
  const s = scfTY({ beta: 0.15, gamma: 40, tau: 1.2, alpha: 3, theta: 15 * Math.PI / 180 });
  assert.equal(s.warnings.length, 5);
  assert.ok(Number.isFinite(s.axial.chordSaddle));
  assert.throws(() => scfTY({ ...P, C: 0.3 }), RangeError);
});

test('X-joint SCFs (hand calculation)', () => {
  const s = scfX(P);
  // (12) 3.87 γ τ β (1.10 − β^1.8) = 3.87 × 3.125 × (1.10 − 0.287175) = 9.8301
  near(s.axial.chordSaddle, 3.87 * 3.125 * (1.10 - 0.5 ** 1.8), 1e-12);
  near(s.axial.chordSaddle, 9.8301, 1e-4);
  // (13) 2.28905 − 3 τ β = 1.53905
  near(s.axial.chordCrown, 1.53905, 1e-4);
  // (14) 1 + 1.9 γ τ^0.5 β^0.9 (1.09 − β^1.7) = 1 + 1.9 × 12.5 × 0.70711 × 0.53589 × (1.09 − 0.30779) = 8.0397
  near(s.axial.braceSaddle, 1 + 1.9 * 12.5 * Math.SQRT1_2 * 0.5 ** 0.9 * (1.09 - 0.5 ** 1.7), 1e-12);
  near(s.axial.braceSaddle, 8.0397, 2e-4);
  // (15) brace crown: 3 + 20.7153 × (−0.0260098) = 2.4612
  near(s.axial.braceCrown, 2.4612, 1e-4);
  // (16) OPB: γ τ β (1.56 − 1.34 β⁴) = 3.125 × 1.47625 = 4.61328
  near(s.opb.chordSaddle, 4.61328, 1e-5);
  near(s.ipb.chordCrown, 2.1302, 1e-4);
});

test('K-joint SCFs: balanced axial Eq. (20)/(21) and OPB Eq. (22) (hand calculation, identical braces)', () => {
  const br = { beta: 0.5, tau: 0.5, theta: Math.PI / 4 };
  const zeta = 0.1;
  const k = scfK({ gamma: 12.5, alpha: 20, zeta, A: br, B: br });
  // (20): τ^0.9 γ^0.5 (0.67 − β² + 1.16β) sinθ [1.64 + 0.29 β^−0.38 atan(8ζ)]
  const e20 = 0.5 ** 0.9 * Math.sqrt(12.5) * (0.67 - 0.25 + 0.58) * Math.SQRT1_2 * (1.64 + 0.29 * 0.5 ** -0.38 * Math.atan(0.8));
  near(k.axial.chord, e20, 1e-12);
  // numerically: 0.53589 × 3.53553 × 1.0 × 0.70711 × (1.64 + 0.29 × 1.30134 × 0.67474) = 1.33972 × 1.89464 = 2.5383
  near(k.axial.chord, 2.5383, 2e-4);
  // (21) gap joint (C = 0): 1 + (1.97 − 1.57 β^0.25) τ^−0.14 (sinθ)^0.7 × (20)
  near(k.axial.brace, 1 + (1.97 - 1.57 * 0.5 ** 0.25) * 0.5 ** -0.14 * Math.SQRT1_2 ** 0.7 * e20, 1e-12);
  // (22): x = 1 + ζ sinθ/β; both braces identical
  const x = 1 + zeta * Math.SQRT1_2 / 0.5;
  const e10 = 12.5 * 0.5 * 0.5 * (1.7 - 1.05 * 0.125) * Math.SQRT1_2 ** 1.6;
  const red = 1 - 0.08 * Math.sqrt(0.5 * 12.5) * Math.exp(-0.8 * x);
  near(k.opb.chordSaddle, e10 * red + e10 * red * 2.05 * Math.sqrt(0.5) * Math.exp(-1.3 * x), 1e-12);
  near(k.opb.braceSaddle, 0.5 ** -0.54 * 12.5 ** -0.05 * 0.76 * k.opb.chordSaddle, 1e-12);
  // overlapping through brace adds the C term
  assert.ok(scfK({ gamma: 12.5, alpha: 20, zeta: -0.1, A: br, B: br, C: 1 }).axial.brace !== scfK({ gamma: 12.5, alpha: 20, zeta: -0.1, A: br, B: br, C: 0 }).axial.brace);
  // KT joint with a negligible central brace load path: diagonal SCFs use ζ = ζAB + ζBC + βB
  const kt = scfKT({ gamma: 12.5, alpha: 20, zetaAB: 0.1, zetaBC: 0.1, A: br, B: { beta: 0.4, tau: 0.4, theta: Math.PI / 2 }, C: br });
  near(kt.A.axial.chord, scfK({ gamma: 12.5, alpha: 20, zeta: 0.6, A: br, B: br }).axial.chord, 1e-12);
  near(kt.A.opb.chordSaddle, kt.C.opb.chordSaddle, 1e-12, 'symmetric KT');
  assert.ok(kt.B.opb.chordSaddle > 0);
});

test('hot-spot stress superposition at eight points (DNV-RP-C203 Eq. 3.3.1)', () => {
  const r = hotSpotStresses({ sigmaX: 10e6, sigmaMy: 20e6, sigmaMz: 5e6, scfAC: 3, scfAS: 6, scfMIP: 2, scfMOP: 4 });
  const h = Math.SQRT1_2;
  const exp = [
    30 + 40, 45 + h * 40 - h * 20, 60 - 20, 45 - h * 40 - h * 20, 30 - 40, 45 - h * 40 + h * 20, 60 + 20, 45 + h * 40 + h * 20,
  ].map((v) => v * 1e6);
  r.sigma.forEach((v, i) => near(v, exp[i], 1e-12, `point ${i + 1}`));
  // point 8 combines the mean axial term with both bending terms at 45°: 45 + (40 + 20)/√2 = 87.4 MPa
  assert.equal(r.governingPoint, 8);
  near(r.max, exp[7], 1e-12);
  // pure axial: crown and saddle values, mean in between
  const a = hotSpotStresses({ sigmaX: 1, scfAC: 2, scfAS: 8, scfMIP: 3, scfMOP: 5 });
  assert.deepEqual(a.sigma, [2, 5, 8, 5, 2, 5, 8, 5]);
  const t = tyJointHotSpots({ D: 1.0, T: 0.04, d: 0.5, t: 0.02, L: 10, thetaDeg: 90, sigmaX: 50e6, chordEnds: 'fixed' });
  near(t.chord.sigma[2], 6.4660 * 50e6, 1e-4, 'chord saddle hot spot');
  near(t.brace.sigma[0], 2.6612 * 50e6, 1e-4, 'brace crown hot spot');
});

test('butt-weld SCFs: plate misalignment and tubular girth welds (DNV-RP-C203 §3.3.7)', () => {
  // plate: 1 + 3 (δm − 0.1 t)/t with t = 20 mm, δm = 4 mm → 1 + 3 × 0.1 = 1.3
  near(scfButtWeldPlate({ t: 0.02, deltaM: 0.004 }), 1.3, 1e-12);
  assert.equal(scfButtWeldPlate({ t: 0.02, deltaM: 0.001 }), 1);
  // tube, equal thickness: 1 + 3 δm/t e^(−√(t/D)) = 1 + 0.45 × e^(−0.141421) = 1.39065
  near(scfGirthWeld({ D: 1.0, t: 0.02, deltaM: 0.003 }).scf, 1.39065, 1e-5);
  // thickness transition T = 30, t = 20, D = 1.0, L = 30 mm, δm = 2 mm, δ0 = 0.1 t:
  // β = 1.5 − 1/log10(50) + 3/log10(50)² = 1.5 − 0.588592 + 1.039322 = 1.950730
  // (T/t)^β = 1.5^1.95073 = 2.20562; 1/(1 + 2.20562) = 0.311952; α = 1.82 × 0.03/√0.02 × 0.311952 = 0.120437
  // SCF = 1 + 6 (0.005 + 0.002 − 0.002)/0.02 × 0.311952 × e^(−0.120437) = 1 + 1.5 × 0.311952 × 0.886533 = 1.41483
  const g = scfGirthWeld({ D: 1.0, t: 0.02, T: 0.03, deltaM: 0.002, delta0: 0.002, L: 0.03 });
  near(g.beta, 1.95073, 1e-5);
  near(g.alpha, 0.120437, 1e-4);
  near(g.scf, 1.41483, 1e-4);
  assert.throws(() => scfGirthWeld({ D: 1.0, t: 0.02, T: 0.03, deltaM: 0.002 }), RangeError);
});
