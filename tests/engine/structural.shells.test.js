import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  lameCylinder, thinCylinder, lameSphere, classicalShellBuckling, c202ElasticBuckling, c202ShellCheck,
  coneEquivalentCylinder, pipeCollapsePressure, ovalTubeYieldPressure, brazierMoment,
} from '../../src/engine/structural/shells.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);

test('Lamé cylinder: boundary conditions, hand values and the thin-wall limit', () => {
  // ri = 100 mm, ro = 200 mm, pi = 60 MPa: σθ(ri) = pi (ro² + ri²)/(ro² − ri²) = 100 MPa, σθ(ro) = 40 MPa,
  // closed-end axial stress = pi ri²/(ro² − ri²) = 20 MPa.
  const a = lameCylinder({ ri: 0.1, ro: 0.2, pi: 60e6 });
  close(a.sigmaTheta, 100e6, 1e-12); close(a.sigmaR, -60e6, 1e-12); close(a.sigmaZ, 20e6, 1e-12);
  close(lameCylinder({ ri: 0.1, ro: 0.2, pi: 60e6, r: 0.2 }).sigmaTheta, 40e6, 1e-12);
  assert.ok(Math.abs(lameCylinder({ ri: 0.1, ro: 0.2, pi: 60e6, r: 0.2 }).sigmaR) < 1e-3);
  close(a.tresca, 160e6, 1e-12);
  // External pressure only: σr(ro) = −po, σθ(ri) = −2 po ro²/(ro² − ri²)
  const e = lameCylinder({ ri: 0.1, ro: 0.2, po: 30e6 });
  close(e.sigmaTheta, (-2 * 30e6 * 0.04) / 0.03, 1e-12);
  close(lameCylinder({ ri: 0.1, ro: 0.2, po: 30e6, r: 0.2 }).sigmaR, -30e6, 1e-12);
  // Equal pressures inside and out: hydrostatic state.
  const h = lameCylinder({ ri: 0.1, ro: 0.2, pi: 5e6, po: 5e6, r: 0.15 });
  close(h.sigmaR, -5e6, 1e-12); close(h.sigmaTheta, -5e6, 1e-12);
  // Radial equilibrium dσr/dr = (σθ − σr)/r checked numerically.
  const f = (r) => lameCylinder({ ri: 0.1, ro: 0.2, pi: 60e6, po: 10e6, r });
  const d = (f(0.1501).sigmaR - f(0.1499).sigmaR) / 2e-4;
  close(d, (f(0.15).sigmaTheta - f(0.15).sigmaR) / 0.15, 1e-6);
  // Thin wall r/t = 100: mean hoop stress → p r/t, axial → p r/(2t)
  const thin = lameCylinder({ ri: 0.995, ro: 1.005, pi: 1e6, r: 1.0 });
  const m = thinCylinder({ r: 1, t: 0.01, p: 1e6 });
  close(thin.sigmaTheta, m.sigmaTheta, 6e-3);
  // Exact end-cap stress is p ri²/(2 r t): one wall-thickness ratio (1 %) below the mean-radius formula.
  close(thin.sigmaZ, m.sigmaZ, 1.1e-2);
  close(m.sigmaTheta, 100e6, 1e-12); close(m.vonMises, 100e6 * Math.sqrt(0.75), 1e-12);
  // Displacement, open ends: u(ri) = (ri/E)(σθ − ν σr)
  close(lameCylinder({ ri: 0.1, ro: 0.2, pi: 60e6, ends: 'open', E: 200e9, nu: 0.3 }).u, (0.1 / 200e9) * (100e6 + 0.3 * 60e6), 1e-12);
  assert.throws(() => lameCylinder({ ri: 0.2, ro: 0.1 }), RangeError);
});

test('Lamé sphere: boundary conditions and thin-wall limit p r / 2t', () => {
  const s = lameSphere({ ri: 0.5, ro: 0.6, pi: 10e6 });
  close(s.sigmaR, -10e6, 1e-12);
  // σθ(ri) = pi (ro³ + 2 ri³) / (2 (ro³ − ri³))
  close(s.sigmaTheta, (10e6 * (0.216 + 0.25)) / (2 * (0.216 - 0.125)), 1e-12);
  assert.ok(Math.abs(lameSphere({ ri: 0.5, ro: 0.6, pi: 10e6, r: 0.6 }).sigmaR) < 1e-3);
  const thin = lameSphere({ ri: 0.999, ro: 1.001, pi: 1e6, r: 1 });
  close(thin.sigmaTheta, (1e6 * 1) / (2 * 0.002), 3e-3);
});

test('classical buckling: σ = 0.605 E t/r (ν = 0.3), ring/tube and sphere pressures', () => {
  const c = classicalShellBuckling({ r: 2, t: 0.02, E: 210e9, nu: 0.3 });
  close(c.sigmaAxial, 0.605 * 210e9 * 0.01, 1e-3);
  close(c.pCylinder, (210e9 * 1e-6) / (4 * 0.91), 1e-12);
  close(c.pSphere, 1.21 * 210e9 * 1e-4, 1.5e-3, 'sphere 1.21 E (t/r)²');
  close(classicalShellBuckling({ r: 2, t: 0.02, E: 210e9, knockDownAxial: 0.3 }).sigmaAxial, 0.3 * c.sigmaAxial, 1e-12);
});

test('DNV-RP-C202 elastic buckling coefficients: limits and consistency with classical theory', () => {
  const E = 210e9; const nu = 0.3; const r = 5; const t = 0.03;
  // Medium-length shell: axial C ≈ ρ ξ, i.e. f_Ea = ρ · 0.605 E t/r · (0.702 π² /(12 · 0.605 √(1−ν²))) → the
  // classical stress times the knock-down ρ = 0.5/√(1 + r/150t).
  const l = 6;
  const b = c202ElasticBuckling({ r, t, l, E, nu });
  close(b.Zl, ((l * l) / (r * t)) * Math.sqrt(0.91), 1e-12);
  const rho = 0.5 / Math.sqrt(1 + r / (150 * t));
  const classical = (E * t) / (r * Math.sqrt(3 * 0.91));
  close(b.fEa, rho * classical, 2e-3, 'axial: knock-down × classical');
  assert.ok(b.fEm > b.fEa, 'bending knock-down is milder than axial');
  // Short shell (Z → 0): plate limits C → ψ: 1 (axial, wide plate strip), 5.34 (shear), 4 (lateral pressure)
  const s = c202ElasticBuckling({ r: 50, t: 0.03, l: 0.05, E, nu });
  close(s.C.axial, 1, 1e-3); close(s.C.torsion, 5.34, 1e-3); close(s.C.lateral, 4, 1e-3); close(s.C.hydrostatic, 2, 2e-3);
  // Hand value, lateral pressure: Z = 228.9, ξ = 1.04 √Z = 15.74, C = 4 √(1 + (0.6·15.74/4)²) = 10.26
  close(b.C.lateral, 4 * Math.sqrt(1 + ((0.6 * 1.04 * Math.sqrt(b.Zl)) / 4) ** 2), 1e-12);
  close(b.Zl, 228.95, 1e-3); close(b.C.lateral, 10.26, 2e-3);
  close(b.fEhLateral, (10.26 * Math.PI ** 2 * E * (t / l) ** 2) / (12 * 0.91), 2e-3);
  // Long shell: f_Eh = 0.25 E (t/r)²
  const long = c202ElasticBuckling({ r: 0.5, t: 0.02, l: 40, E, nu });
  assert.ok(long.long);
  close(long.fEhLateral, 0.25 * E * (0.02 / 0.5) ** 2, 1e-12);
  close(long.fEtau, 0.25 * E * (0.02 / 0.5) ** 1.5, 1e-12);
});

test('DNV-RP-C202 shell check: slenderness, strength and material factor (hand calculation)', () => {
  const E = 210e9; const fy = 355e6; const r = 5; const t = 0.03; const l = 6;
  const el = c202ElasticBuckling({ r, t, l, E });
  // Pure axial compression of 100 MPa: λs² = fy/σj · σa0/f_Ea = fy/f_Ea
  const c = c202ShellCheck({ r, t, l, E, fy, sigmaA: -100e6 });
  close(c.sigmaJ, 100e6, 1e-12);
  close(c.lambdaS, Math.sqrt(fy / el.fEa), 1e-12);
  close(c.fks, fy / Math.sqrt(1 + (fy / el.fEa) ** 2), 1e-12);
  const g = c.lambdaS < 0.5 ? 1.15 : c.lambdaS <= 1 ? 0.85 + 0.6 * c.lambdaS : 1.45;
  close(c.gammaM, g, 1e-12);
  close(c.utilisation, 100e6 / (c.fks / g), 1e-12);
  // Axial tension does not contribute to the slenderness; external pressure does.
  const p = 0.4e6;
  const h = c202ShellCheck({ r, t, l, E, fy, sigmaA: 50e6, sigmaH: (-p * r) / t });
  const sj = Math.sqrt(50e6 ** 2 + 50e6 * ((p * r) / t) + ((p * r) / t) ** 2);
  close(h.sigmaJ, sj, 1e-12);
  close(h.lambdaS ** 2, (fy / sj) * (((p * r) / t) / el.fEhLateral), 1e-12);
  // Hydrostatic pressure uses the lower ψ = 2 coefficient: more slender.
  const hy = c202ShellCheck({ r, t, l, E, fy, sigmaH: (-p * r) / t, hydrostatic: true });
  const la = c202ShellCheck({ r, t, l, E, fy, sigmaH: (-p * r) / t });
  assert.ok(hy.lambdaS > la.lambdaS);
  // Very stocky shell: f_ks → fy, γM = 1.15
  const st = c202ShellCheck({ r: 0.2, t: 0.04, l: 0.3, E, fy, sigmaA: -10e6 });
  close(st.fks, fy, 5e-3); assert.equal(st.gammaM, 1.15);
  assert.throws(() => c202ShellCheck({ r, t, l, E, fy }), RangeError);
});

test('cone: equivalent cylinder', () => {
  const c = coneEquivalentCylinder({ r1: 2, r2: 3, l: 4 });
  close(c.alpha, Math.atan(0.25), 1e-12);
  close(c.rEquivalent, 2.5 / Math.cos(c.alpha), 1e-12);
  close(c.lEquivalent, Math.hypot(4, 1), 1e-12, 'slant length');
  assert.ok(c.withinCode);
  assert.ok(!coneEquivalentCylinder({ r1: 1, r2: 4, l: 2 }).withinCode);
});

test('pipe collapse: cubic limits (elastic, plastic, zero ovality) and residual', () => {
  const E = 207e9; const fy = 450e6;
  // Thin pipe, D/t = 60: elastic collapse governs; with f0 → 0 the root is min(pel, pp)
  const thin = pipeCollapsePressure({ D: 0.6, t: 0.01, E, fy, f0: 0 });
  close(thin.pel, (2 * E * (1 / 60) ** 3) / 0.91, 1e-12);
  close(thin.pc, Math.min(thin.pel, thin.pp), 1e-12);
  const thick = pipeCollapsePressure({ D: 0.3, t: 0.03, E, fy, f0: 0 });
  close(thick.pc, thick.pp, 1e-12); close(thick.pp, fy * 0.2, 1e-12);
  // General case satisfies the cubic and lies below both limits; ovality reduces the collapse pressure.
  const D = 0.508; const t = 0.0206;
  const a = pipeCollapsePressure({ D, t, E, fy, f0: 0.005 });
  const res = (a.pc - a.pel) * (a.pc ** 2 - a.pp ** 2) - (a.pc * a.pel * a.pp * 0.005 * D) / t;
  assert.ok(Math.abs(res) < 1e-9 * a.pel * a.pp * a.pp);
  assert.ok(a.pc < Math.min(a.pel, a.pp));
  const b = pipeCollapsePressure({ D, t, E, fy, f0: 0.02 });
  assert.ok(b.pc < a.pc);
  close(pipeCollapsePressure({ D, t, E, fy, alphaFab: 0.85 }).pp, 0.85 * a.pp, 1e-12);
});

test('Timoshenko oval tube: first-yield pressure limits and stress identity', () => {
  const E = 200e9; const fy = 300e6; const r = 0.5; const t = 0.01;
  // Perfect tube: the smaller of buckling and squash pressures
  const p0 = ovalTubeYieldPressure({ r, t, E, fy, w0: 0 });
  close(p0.pYield, Math.min(p0.pCr, p0.pSquash), 1e-9);
  close(p0.pCr, (E * (t / r) ** 3) / (4 * 0.91), 1e-12);
  // With ovality the extreme-fibre stress equals fy: p m + 6 p m² n / (1 − p/pcr) = fy
  const w0 = 0.002;
  const q = ovalTubeYieldPressure({ r, t, E, fy, w0 });
  const m = r / t; const n = w0 / r;
  close(q.pYield * m + (6 * q.pYield * m * m * n) / (1 - q.pYield / q.pCr), fy, 1e-9);
  assert.ok(q.pYield < p0.pYield);
});

test('Brazier moment: maximum of M(κ) = π E r³ t κ (1 − 3ζ/2) at ζ = 2/9; 0.987 E r t²/√(1 − ν²)', () => {
  const o = { r: 0.5, t: 0.005, E: 210e9, nu: 0.3 };
  const b = brazierMoment(o);
  close(b.momentLimit, (0.987 * 210e9 * 0.5 * 0.005 ** 2) / Math.sqrt(0.91), 1e-3);
  const at = brazierMoment({ ...o, kappa: b.curvatureLimit });
  close(at.moment, b.momentLimit, 1e-12); close(at.flattening, 2 / 9, 1e-12);
  for (const f of [0.9, 1.1]) assert.ok(brazierMoment({ ...o, kappa: f * b.curvatureLimit }).moment < b.momentLimit);
  // Small curvature: beam theory M = E I κ with I = π r³ t
  const k = 1e-4 * b.curvatureLimit;
  close(brazierMoment({ ...o, kappa: k }).moment, 210e9 * Math.PI * 0.5 ** 3 * 0.005 * k, 1e-6);
});
