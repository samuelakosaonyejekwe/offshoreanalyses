import test from 'node:test';
import assert from 'node:assert/strict';
import {
  seawaterDensity, seawaterSecantBulkModulus, brineDensity, waterDensity, waterViscosity, seawaterViscosity,
  brineViscosity, seawaterConductivity, seawaterHeatCapacity, seawaterSoundSpeed, waterVapourPressure,
  hydrostaticPressure,
} from '../../src/engine/flow/water.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const ATM = 101325;

test('EOS-80 density reproduces the UNESCO Technical Paper 44 check values', () => {
  // UNESCO 44 (Fofonoff & Millard 1983) check table: (S, t °C, p bar) → ρ kg/m³
  const cases = [
    [0, 5, 0, 999.96675], [35, 5, 0, 1027.67547], [35, 5, 1000, 1069.48914],
    [0, 25, 1000, 1037.90204], [35, 25, 1000, 1062.53817],
  ];
  for (const [S, t, p, rho] of cases) {
    assert.ok(Math.abs(seawaterDensity({ T: t + 273.15, S, P: ATM + p * 1e5 }) - rho) < 2e-5, `S=${S} t=${t} p=${p}`);
  }
});

test('EOS-80 secant bulk modulus is consistent with the density ratio', () => {
  const T = 283.15;
  const P = ATM + 3e7;
  const K = seawaterSecantBulkModulus({ T, S: 35, P });
  const r0 = seawaterDensity({ T, S: 35 });
  assert.ok(rel(seawaterDensity({ T, S: 35, P }), r0 / (1 - 3e7 / K)) < 1e-12);
});

test('IAPWS-IF97 saturation pressure check value at 300 K and Antoine at the normal boiling point', () => {
  assert.ok(rel(waterVapourPressure({ T: 300 }), 3536.58941) < 1e-8); // IF97 Table 35
  assert.ok(rel(waterVapourPressure({ T: 500 }), 2.63889776e6) < 1e-8);
  assert.ok(rel(waterVapourPressure({ T: 373.15, method: 'antoine' }), 101325) < 5e-3);
  assert.ok(rel(waterVapourPressure({ T: 323.15, method: 'antoine' }), waterVapourPressure({ T: 323.15 })) < 5e-3);
  // salt lowers the vapour pressure (Raoult): factor 1/(1 + 0.57357·S/(1000 − S))
  assert.ok(rel(waterVapourPressure({ T: 300, S: 35 }), 3536.58941 / (1 + 0.57357 * 35 / 965)) < 1e-8);
});

test('Mackenzie sound speed check value (25 °C, S = 35, 1000 m → 1550.744 m/s)', () => {
  assert.ok(Math.abs(seawaterSoundSpeed({ T: 298.15, S: 35, depth: 1000 }) - 1550.744) < 1e-3);
});

test('viscosity: pure water 20 °C = 1.002 mPa·s (IAPWS) and seawater 35 g/kg ≈ 1.08 mPa·s', () => {
  assert.ok(rel(waterViscosity(293.15), 1.0016e-3) < 3e-3);
  assert.ok(rel(seawaterViscosity({ T: 293.15, S: 35 }), 1.077e-3) < 0.02);
  // hand evaluation of the salinity factor at 20 °C, s = 0.035: A = 1.541 + 0.3996 − 0.03808, B = 7.974 − 1.5122 + 0.18896
  const A = 1.541 + 1.998e-2 * 20 - 9.52e-5 * 400;
  const B = 7.974 - 7.561e-2 * 20 + 4.724e-4 * 400;
  assert.ok(rel(seawaterViscosity({ T: 293.15, S: 35 }) / waterViscosity(293.15), 1 + A * 0.035 + B * 0.035 ** 2) < 1e-12);
});

test('heat capacity and conductivity of water and seawater at 20 °C', () => {
  assert.ok(rel(seawaterHeatCapacity({ T: 293.15, S: 0 }), 4184) < 3e-3);
  assert.ok(rel(seawaterHeatCapacity({ T: 293.15, S: 35 }), 3993) < 5e-3);
  assert.ok(rel(seawaterConductivity({ T: 293.15, S: 0 }), 0.598) < 0.015);
  assert.ok(seawaterConductivity({ T: 293.15, S: 35 }) < seawaterConductivity({ T: 293.15, S: 0 }));
});

test('McCain brine density: fresh-water standard density and agreement with EOS-80 for seawater', () => {
  const d0 = brineDensity({ T: 288.7056, S: 0 });
  assert.ok(rel(d0.rhoSc, 62.368 * 16.01846337) < 1e-12);
  assert.ok(rel(brineDensity({ T: 293.15, S: 35 }).rho, seawaterDensity({ T: 293.15, S: 35 })) < 3e-3);
  // dispatcher: EOS-80 inside its range, McCain outside
  assert.equal(waterDensity({ T: 283.15, S: 35 }), seawaterDensity({ T: 283.15, S: 35 }));
  assert.equal(waterDensity({ T: 363.15, S: 100 }), brineDensity({ T: 363.15, S: 100 }).rho);
  assert.ok(brineViscosity({ T: 350, S: 100 }) > brineViscosity({ T: 350, S: 0 }));
});

test('hydrostatic pressure lies between the incompressible bounds and converges at second order', () => {
  const depth = 3000;
  const T = 277.15;
  const h = hydrostaticPressure({ depth, T, S: 35, n: 200 });
  const rhoTop = seawaterDensity({ T, S: 35 });
  assert.ok(h.gauge > rhoTop * 9.80665 * depth);
  assert.ok(h.gauge < h.rhoBottom * 9.80665 * depth);
  const e1 = Math.abs(hydrostaticPressure({ depth, T, S: 35, n: 2 }).P - h.P);
  const e2 = Math.abs(hydrostaticPressure({ depth, T, S: 35, n: 4 }).P - h.P);
  assert.ok(e1 / e2 > 3 && e1 / e2 < 5, `observed ratio ${e1 / e2}`);
});

test('invalid input throws RangeError naming the argument', () => {
  assert.throws(() => seawaterDensity({ T: 400, S: 35 }), /RangeError|T =/);
  assert.throws(() => seawaterDensity({ T: 280, S: -1 }), RangeError);
  assert.throws(() => waterVapourPressure({ T: 300, method: 'nope' }), RangeError);
});
