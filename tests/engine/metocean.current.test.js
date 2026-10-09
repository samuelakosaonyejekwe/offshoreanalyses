import test from 'node:test';
import assert from 'node:assert/strict';
import * as C from '../../src/engine/metocean/current.js';
import { wavenumber } from '../../src/engine/metocean/waves.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('1/7 power-law profile: surface and bed values, mid-depth hand value, depth average = 7/8 v0', () => {
  const o = { v0: 1.2, d: 40 };
  assert.equal(C.currentPowerLaw(0, o), 1.2);
  assert.equal(C.currentPowerLaw(-40, o), 0);
  assert.ok(rel(C.currentPowerLaw(-20, o), 1.2 * 0.5 ** (1 / 7)) < 1e-15); // 1.087 m/s
  assert.ok(rel(C.depthAveragedCurrent((z) => C.currentPowerLaw(z, o), { d: 40 }), (1.2 * 7) / 8) < 1e-7);
  assert.ok(rel(C.depthAveragedCurrent((z) => C.currentPowerLaw(z, { ...o, alpha: 0.25 }), { d: 40 }), 1.2 / 1.25) < 1e-6);
  assert.throws(() => C.currentPowerLaw(-41, o), RangeError);
});

test('linear wind-generated profile vanishes at the reference depth; depth average of a triangle', () => {
  const o = { v0: 0.6, d0: 50, d: 120 };
  assert.equal(C.currentLinear(0, o), 0.6);
  assert.ok(rel(C.currentLinear(-25, o), 0.3) < 1e-15);
  assert.equal(C.currentLinear(-50, o), 0);
  assert.equal(C.currentLinear(-100, o), 0);
  // in water shallower than d0 the profile reaches zero at the bed
  assert.ok(rel(C.currentLinear(-10, { v0: 0.6, d0: 50, d: 20 }), 0.3) < 1e-15);
  assert.ok(rel(C.depthAveragedCurrent((z) => C.currentLinear(z, o), { d: 120, n: 4000 }), (0.5 * 0.6 * 50) / 120) < 1e-4);
  // tide + wind superposition at the surface
  assert.ok(rel(C.currentTideAndWind(0, { vTide: 1.0, vWind: 0.6, d: 120 }), 1.6) < 1e-15);
  assert.ok(rel(C.currentTideAndWind(-60, { vTide: 1.0, vWind: 0.6, d: 120 }), 0.5 ** (1 / 7)) < 1e-15);
});

test('bilinear and tabulated profiles pass through their defining points', () => {
  const o = { vSurface: 1.5, vBreak: 0.9, vBed: 0.3, zBreak: -30, d: 100 };
  assert.equal(C.currentBilinear(0, o), 1.5);
  assert.ok(rel(C.currentBilinear(-30, o), 0.9) < 1e-15);
  assert.ok(rel(C.currentBilinear(-100, o), 0.3) < 1e-15);
  assert.ok(rel(C.currentBilinear(-15, o), 1.2) < 1e-15);
  assert.ok(rel(C.currentBilinear(-65, o), 0.6) < 1e-15);
  const t = { zs: [-100, -30, 0], vs: [0.3, 0.9, 1.5] };
  for (const z of [-100, -80, -65, -30, -15, 0]) assert.ok(rel(C.currentTabulated(z, t), C.currentBilinear(z, o)) < 1e-14);
  assert.equal(C.currentTabulated(5, t), 1.5);
  assert.equal(C.currentTabulated(-500, t), 0.3);
  // area under the bilinear profile: 30·(1.5+0.9)/2 + 70·(0.9+0.3)/2 = 78 m²/s
  assert.ok(rel(C.depthAveragedCurrent((z) => C.currentBilinear(z, o), { d: 100, n: 20000 }), 0.78) < 1e-5);
});

test('logarithmic near-bed profile (DNV-RP-F109 Eq. 3.2) and mean over a pipe diameter (Eq. 3.3)', () => {
  const o = { vRef: 0.8, zRef: 3, z0: C.SEABED_ROUGHNESS.mediumSand };
  assert.equal(C.SEABED_ROUGHNESS.mediumSand, 4e-5);
  assert.ok(rel(C.currentLogNearBed(3, o), 0.8) < 1e-15);
  assert.equal(C.currentLogNearBed(0, o), 0);
  assert.ok(rel(C.currentLogNearBed(1, o), (0.8 * Math.log(1 / 4e-5 + 1)) / Math.log(3 / 4e-5 + 1)) < 1e-15); // 0.722 m/s
  // Eq. 3.3 is the exact average of the profile over 0..D: compare with numerical quadrature
  const D = 0.5; let s = 0; const n = 200000;
  for (let i = 0; i < n; i++) s += C.currentLogNearBed((D * (i + 0.5)) / n, o) / n;
  assert.ok(rel(C.meanCurrentOverHeight({ ...o, D }), s) < 1e-6);
  assert.ok(rel(C.meanCurrentOverHeight({ ...o, D, theta: Math.PI / 6 }), 0.5 * s) < 1e-6);
});

test('wind-driven surface current by the percentage rule', () => {
  assert.ok(rel(C.windDrivenSurfaceCurrent(30), 0.9) < 1e-15);
  assert.ok(rel(C.windDrivenSurfaceCurrent(30, 0.015), 0.45) < 1e-15);
  assert.throws(() => C.windDrivenSurfaceCurrent(30, 0.5), RangeError);
});

test('current stretching: end points, defining equations, and reduction to no stretching for η = 0', () => {
  const d = 50; const k = 0.05;
  for (const method of ['linear', 'nonlinear']) {
    for (const eta of [4, -3]) {
      assert.ok(Math.abs(C.stretchedElevation({ z: eta, eta, d, method, k })) < 1e-10);
      assert.ok(Math.abs(C.stretchedElevation({ z: -d, eta, d, method, k }) + d) < 1e-10);
    }
    assert.ok(Math.abs(C.stretchedElevation({ z: -12, eta: 0, d, method, k }) + 12) < 1e-10);
  }
  // linear: z' = (z − η) d / (d + η)
  assert.ok(rel(C.stretchedElevation({ z: -10, eta: 4, d }), (-14 * 50) / 54) < 1e-15);
  // nonlinear: z = z' + η sinh(k(z' + d)) / sinh(kd)
  for (const z of [-40, -20, -5, 2]) {
    const zp = C.stretchedElevation({ z, eta: 4, d, method: 'nonlinear', k });
    assert.ok(Math.abs(zp + (4 * Math.sinh(k * (zp + d))) / Math.sinh(k * d) - z) < 1e-10);
  }
  // nonlinear stretching displaces deep water less than linear stretching does
  const zl = C.stretchedElevation({ z: -40, eta: 4, d, method: 'linear' });
  const zn = C.stretchedElevation({ z: -40, eta: 4, d, method: 'nonlinear', k });
  assert.ok(Math.abs(zn + 40) < Math.abs(zl + 40));
  // a slab profile is unchanged; a sheared profile keeps its surface value at the crest
  const prof = (z) => C.currentPowerLaw(z, { v0: 1, d });
  assert.ok(rel(C.stretchedCurrent(prof, { z: 4, eta: 4, d }), 1) < 1e-12);
  assert.equal(C.stretchedElevation({ z: 3, eta: 4, d, method: 'none' }), 0);
  assert.throws(() => C.stretchedElevation({ z: 5, eta: 4, d }), RangeError);
  assert.throws(() => C.stretchedElevation({ z: 0, eta: 4, d, method: 'nonlinear' }), RangeError);
});

test('linear stretching conserves the depth-integrated flux per unit still-water depth scaling', () => {
  // ∫_{−d}^{η} v(z'(z)) dz = (d + η)/d · ∫_{−d}^{0} v dz' for the linear mapping
  const d = 30; const eta = 3;
  const prof = (z) => C.currentBilinear(z, { vSurface: 1.2, vBreak: 0.8, vBed: 0.2, zBreak: -10, d });
  let s = 0; const n = 20000;
  for (let i = 0; i < n; i++) {
    const z = -d + ((d + eta) * (i + 0.5)) / n;
    s += (C.stretchedCurrent(prof, { z, eta, d }) * (d + eta)) / n;
  }
  const base = C.depthAveragedCurrent(prof, { d, n: 20000 }) * d;
  assert.ok(rel(s, (base * (d + eta)) / d) < 1e-5);
});

test('combined wave and current near the bed: linear orbital velocity π H / (T sinh kd) and vector combination', () => {
  const H = 4; const T = 10; const d = 25;
  const k = wavenumber({ T, d });
  const w = C.nearBedWaveCurrent({ H, T, d });
  assert.ok(rel(w.Uw, (Math.PI * H) / (T * Math.sinh(k * d))) < 1e-13);
  assert.ok(rel(w.Aw, (w.Uw * T) / (2 * Math.PI)) < 1e-13);
  assert.ok(rel(w.Umax, w.Uw) < 1e-15 && w.Umin === 0);
  // collinear current adds directly; perpendicular current adds in quadrature (no Doppler shift)
  const co = C.nearBedWaveCurrent({ H, T, d, Uc: 0.5, angle: 0 });
  assert.ok(rel(co.Umax, 0.5 + co.Uw) < 1e-14);
  assert.ok(co.k < k && co.Tr > T); // following current: longer wave, longer intrinsic period
  const pe = C.nearBedWaveCurrent({ H, T, d, Uc: 0.5, angle: Math.PI / 2 });
  assert.ok(rel(pe.Uw, w.Uw) < 1e-12);
  assert.ok(rel(pe.Umax, Math.hypot(0.5, w.Uw)) < 1e-12);
  assert.ok(rel(pe.Umin, 0.5) < 1e-12);
  // elevation above the bed: cosh(k z) growth
  assert.ok(rel(C.nearBedWaveCurrent({ H, T, d, zBed: 2 }).Uw, w.Uw * Math.cosh(2 * k)) < 1e-13);
});
