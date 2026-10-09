import test from 'node:test';
import assert from 'node:assert/strict';
import {
  columnResonance, chamberCompliance, columnBody, owcRegular, optimalTurbineDamping, owcIrregular,
  owcPowerMatrix,
} from '../../src/engine/aero/owc.js';
import { pointAbsorber, waveSpectrum, regularWavePower } from '../../src/engine/aero/wec.js';

const G = 9.80665;
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);
const grid = (a, b, n) => Float64Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));

test('water-column resonance T = 2π√((d + l_a)/g) and chamber compliance by hand', () => {
  // d = 4 m, l_a = 1 m: 2π √(5/9.80665) = 4.4865 s
  close(columnResonance({ draft: 4, addedLength: 1 }).period, 2 * Math.PI * Math.sqrt(5 / G), 1e-12);
  close(columnResonance({ draft: 4, addedLength: 1 }).period, 4.4865, 1e-4);
  // V0 = 500 m³: C = 500/(1.4·101325) = 3.5247e-3 m³/Pa
  close(chamberCompliance({ volume: 500 }), 500 / (1.4 * 101325), 1e-15);
  close(chamberCompliance({ volume: 500 }), 3.5247e-3, 1e-4);
});

const body = columnBody({ radius: 3, draft: 4, addedLength: 1 });

test('incompressible chamber with a Wells turbine is a piston with a linear damper S²/G_t', () => {
  const w = 1.1;
  const col = body.column(w);
  const Gt = 4e-3;
  const o = owcRegular({ column: body.column, area: body.area, omega: w, amplitude: 0.8, compliance: 0, conductance: Gt });
  // equivalent mechanical damper
  const Bpto = body.area ** 2 / Gt;
  const pa = pointAbsorber({ mass: col.M, A: 0, B: col.B, K: col.K, Fe: col.Fe, omega: w, amplitude: 0.8, control: 'fixed', Bpto });
  close(o.pneumaticPower, pa.power, 1e-10);
  close(o.stroke, pa.stroke, 1e-10);
  close(o.flow, body.area * w * o.stroke, 1e-10);
  close(o.pressure, o.flow / Gt, 1e-10);
  close(o.pneumaticPower, 0.5 * o.pressure * o.flow, 1e-10);
});

test('optimum turbine damping at resonance reaches the axisymmetric limit J/k', () => {
  const w0 = columnResonance({ draft: 4, addedLength: 1 }).omega;
  const col = body.column(w0);
  const opt = optimalTurbineDamping({ column: body.column, area: body.area, omega: w0, amplitude: 0.5, compliance: 0 });
  close(opt.pneumaticPower, (col.Fe * 0.5) ** 2 / (8 * col.B), 1e-6);
  close(opt.damping, col.B / body.area ** 2, 1e-3);
  const J = regularWavePower({ amplitude: 0.5, omega: w0 });
  close(opt.captureWidth, G / (w0 * w0), 1e-6);
  close(opt.pneumaticPower / J, opt.captureWidth, 1e-12);
  // air compressibility acts as a spring in series with the turbine: it cannot raise the optimum above the limit
  const C = chamberCompliance({ volume: 400 });
  const comp = optimalTurbineDamping({ column: body.column, area: body.area, omega: w0, amplitude: 0.5, compliance: C });
  assert.ok(comp.pneumaticPower <= opt.pneumaticPower * (1 + 1e-9));
  // with the same turbine, compressibility changes the power (spring effect)
  const a = owcRegular({ column: body.column, area: body.area, omega: 0.9, amplitude: 0.5, compliance: 0, conductance: 3e-3 });
  const b = owcRegular({ column: body.column, area: body.area, omega: 0.9, amplitude: 0.5, compliance: C, conductance: 3e-3 });
  assert.ok(Math.abs(a.pneumaticPower - b.pneumaticPower) / a.pneumaticPower > 1e-3);
  // chamber mass balance amplitude: |S ẋ| = |q + iωC p| = p √(G_t² + ω²C²)
  close(body.area * 0.9 * b.stroke, b.pressure * Math.hypot(3e-3, 0.9 * C), 1e-10);
});

test('impulse turbine: harmonic linearisation is self-consistent and energy-equivalent', () => {
  const ki = 40;
  const r = owcRegular({ column: body.column, area: body.area, omega: 1.0, amplitude: 0.6, compliance: chamberCompliance({ volume: 300 }), turbine: 'impulse', ki });
  // G_t = 1/((8/3π) k_i q̂) and q̂ = G_t p̂
  close(r.conductance, 1 / ((8 / (3 * Math.PI)) * ki * r.flow), 1e-9);
  close(r.flow, r.conductance * r.pressure, 1e-12);
  // mean power of p = k q|q| with q = q̂ sin: k q̂³ · (4/3π)
  close(r.pneumaticPower, ki * r.flow ** 3 * (4 / (3 * Math.PI)), 1e-9);
  // stiffer turbine → smaller flow
  assert.ok(owcRegular({ column: body.column, area: body.area, omega: 1.0, amplitude: 0.6, compliance: 0, turbine: 'impulse', ki: 400 }).flow < r.flow);
});

test('turbine efficiency: constant value and flow-dependent curve', () => {
  const o = { column: body.column, area: body.area, omega: 1.0, amplitude: 0.6, compliance: 0, conductance: 4e-3 };
  const c = owcRegular({ ...o, efficiency: 0.6 });
  close(c.electricalPower, 0.6 * c.pneumaticPower, 1e-12);
  const f = owcRegular({ ...o, efficiency: () => 0.6 });
  close(f.electricalPower, c.electricalPower, 1e-12);
  // efficiency proportional to |q|/q̂ (capped): power-weighted mean over a sinusoid = ∫sin³/∫sin² = 8/(3π)
  const lin = owcRegular({ ...o, efficiency: (q) => q / c.flow });
  close(lin.efficiency, 8 / (3 * Math.PI), 1e-4);
});

test('irregular waves: P = G_t σp², quadratic in Hs for the Wells turbine; stochastic linearisation for the impulse turbine', () => {
  const om = grid(0.2, 3, 300);
  const C = chamberCompliance({ volume: 300 });
  const run = (Hs) => owcIrregular({ column: body.column, area: body.area, omega: om, S: waveSpectrum({ Hs, Tp: 8, omega: om }), compliance: C, conductance: 4e-3 });
  const r1 = run(1);
  const r2 = run(2);
  close(r1.pneumaticPower, 4e-3 * r1.sigmaP ** 2, 1e-12);
  close(r2.pneumaticPower, 4 * r1.pneumaticPower, 1e-10);
  close(r1.sigmaQ, 4e-3 * r1.sigmaP, 1e-12);
  const imp = owcIrregular({ column: body.column, area: body.area, omega: om, S: waveSpectrum({ Hs: 2, Tp: 8, omega: om }), compliance: C, turbine: 'impulse', ki: 40 });
  close(imp.conductance, 1 / (Math.sqrt(8 / Math.PI) * 40 * imp.sigmaQ), 1e-8);
  // E[k |q|³] for Gaussian q: k · 2√(2/π) σ³
  close(imp.pneumaticPower, 40 * 2 * Math.sqrt(2 / Math.PI) * imp.sigmaQ ** 3, 1e-8);
});

test('power matrix and annual energy', () => {
  const C = chamberCompliance({ volume: 300 });
  const o = { hs: [1, 2, 4], tp: [6, 9], column: body.column, area: body.area, compliance: C, conductance: 4e-3, efficiency: 0.55, ratedPower: 60e3 };
  const m = owcPowerMatrix(o);
  close(m.power[1][0] / m.power[0][0], 4, 1e-9);
  assert.ok(m.power.flat().every((p) => p <= 60e3));
  assert.ok(Number.isNaN(m.aepWh));
  const prob = [[0, 0], [0, 1], [0, 0]];
  const e = owcPowerMatrix({ ...o, prob, availability: 0.9 });
  close(e.aepWh, m.power[1][1] * 0.9 * 8760, 1e-12);
  close(e.capacityFactor, (m.power[1][1] * 0.9) / 60e3, 1e-12);
  assert.throws(() => owcRegular({ column: body.column, area: body.area, omega: 1, amplitude: 1, compliance: 0 }), RangeError);
});
