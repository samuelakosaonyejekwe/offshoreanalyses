import test from 'node:test';
import assert from 'node:assert/strict';
import {
  METALS, INSULATION, SKIN_PROXIMITY, dcResistance20, resistanceAtTemperature, conductorDiameter, sheathResistance20,
  armourResistance20, armourWireCount, skinProximityCoefficients, standardSizes, metalProperties, insulationProperties,
} from '../../src/engine/electrical/conductors.js';

const rel = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${a} vs ${b} (tol ${tol})`);

test('IEC 60228 Table 2 class-2 maximum DC resistance at 20 °C, copper (Ω/km)', () => {
  const expected = { 240: 0.0754, 500: 0.0366, 630: 0.0283, 800: 0.0221, 1000: 0.0176, 1200: 0.0151, 1600: 0.0113, 2000: 0.0090 };
  for (const [mm2, r] of Object.entries(expected)) {
    const out = dcResistance20({ material: 'Cu', area: Number(mm2) * 1e-6 });
    assert.equal(out.tabulated, true);
    rel(out.R20 * 1e3, r, 1e-12);
  }
});

test('IEC 60228 Table 2 class-2 maximum DC resistance at 20 °C, aluminium (Ω/km)', () => {
  const expected = { 240: 0.125, 500: 0.0605, 630: 0.0469, 800: 0.0367, 1000: 0.0291, 1200: 0.0247, 1600: 0.0186, 2000: 0.0149 };
  for (const [mm2, r] of Object.entries(expected)) rel(dcResistance20({ material: 'Al', area: Number(mm2) * 1e-6 }).R20 * 1e3, r, 1e-12);
});

test('tabulated resistances are consistent with ρ20/A_nominal: ratio between 1.00 and 1.14 for every size (no table typos)', () => {
  for (const material of ['Cu', 'Al']) {
    for (const area of standardSizes()) {
      const ratio = dcResistance20({ material, area }).R20 / (METALS[material].rho20 / area);
      assert.ok(ratio > 1.0 && ratio < 1.14, `${material} ${area * 1e6} mm²: ratio ${ratio}`);
    }
  }
});

test('non-tabulated size (3000 mm²) is computed from resistivity with the stated lay factor', () => {
  const out = dcResistance20({ material: 'Cu', area: 3000e-6, layFactor: 1.03 });
  assert.equal(out.tabulated, false);
  rel(out.R20, (1.7241e-8 * 1.03) / 3000e-6, 1e-12);
  assert.throws(() => dcResistance20({ material: 'Cu', area: 5000e-6 }), RangeError);
  assert.throws(() => dcResistance20({ material: 'lead', area: 1e-4 }), RangeError);
});

test('temperature correction by hand: 1000 mm² Cu at 90 °C = 0.0176·(1 + 0.00393·70) = 0.022442 Ω/km', () => {
  const R = resistanceAtTemperature({ R20: 0.0176e-3, thetaC: 90, material: 'Cu' });
  rel(R * 1e3, 0.0224418, 1e-5);
  rel(resistanceAtTemperature({ R20: 1, thetaC: 20, alpha20: 0.004 }), 1, 1e-15);
});

test('IEC 60287-1-1 Table 1 resistivities and temperature coefficients', () => {
  assert.equal(METALS.Cu.rho20, 1.7241e-8); assert.equal(METALS.Cu.alpha20, 3.93e-3);
  assert.equal(METALS.Al.rho20, 2.8264e-8); assert.equal(METALS.Al.alpha20, 4.03e-3);
  assert.equal(METALS.lead.rho20, 21.4e-8); assert.equal(METALS.lead.alpha20, 4.0e-3);
  assert.equal(METALS.steel.rho20, 13.8e-8); assert.equal(METALS.steel.alpha20, 4.5e-3);
});

test('IEC 60949 constant K follows from first principles: K = sqrt(σc·(β + 20)/ρ20)·1e-6 within 1 %', () => {
  for (const m of ['Cu', 'Al', 'lead', 'steel']) {
    const p = metalProperties(m);
    const K = Math.sqrt((p.heatCapacity * (p.beta + 20)) / p.rho20) * 1e-6;
    rel(K, p.K, 0.01);
  }
});

test('thermal resistivities and dielectric constants of insulations (IEC 60287-2-1 Table 1, 60287-1-1 Table 3)', () => {
  assert.equal(INSULATION.xlpe.rhoT, 3.5); assert.equal(INSULATION.pe.rhoT, 3.5);
  assert.equal(INSULATION.epr.rhoT, 5.0); assert.equal(INSULATION.epr.rhoTLowVoltage, 3.5);
  assert.equal(INSULATION.pvc.rhoT, 6.0); assert.equal(INSULATION.pvc.rhoTLowVoltage, 5.0);
  assert.equal(INSULATION.paper.rhoT, 6.0);
  assert.equal(INSULATION.xlpe.epsR, 2.5); assert.equal(INSULATION.xlpe.tanDelta, 0.001); assert.equal(INSULATION.xlpe.tanDeltaMv, 0.004);
  assert.equal(insulationProperties('paper').epsR, 4.0);
  assert.throws(() => insulationProperties('unobtainium'), RangeError);
});

test('skin and proximity coefficients (IEC 60287-1-1 Table 2)', () => {
  assert.deepEqual(skinProximityCoefficients('roundStranded'), { ks: 1, kp: 1 });
  assert.deepEqual(skinProximityCoefficients('roundStrandedImpregnated'), { ks: 1, kp: 0.8 });
  assert.deepEqual(skinProximityCoefficients('milliken'), { ks: 0.435, kp: 0.37 });
  assert.ok(Object.isFrozen(SKIN_PROXIMITY));
  assert.throws(() => skinProximityCoefficients('square'), RangeError);
});

test('conductor diameter by hand: 1000 mm², fill 0.9 → sqrt(4e-3/(0.9π)) = 37.61 mm; solid → 35.68 mm', () => {
  rel(conductorDiameter({ area: 1000e-6, fillFactor: 0.9 }), 0.037613, 1e-4);
  rel(conductorDiameter({ area: 1000e-6, fillFactor: 1 }), 0.0356825, 1e-5);
  assert.throws(() => conductorDiameter({ area: 1e-3, fillFactor: 1.2 }), RangeError);
});

test('lead sheath resistance by hand: ρ/(π·d·t) = 21.4e-8/(π·0.08·0.003) = 2.838e-4 Ω/m', () => {
  rel(sheathResistance20({ material: 'lead', meanDiameter: 0.08, thickness: 0.003 }), 2.8382e-4, 1e-4);
});

test('armour resistance by hand: 60 steel wires of 5 mm → 13.8e-8/(60·π·0.005²/4) = 1.1714e-4 Ω/m; lay factor scales', () => {
  rel(armourResistance20({ material: 'steel', wireDiameter: 0.005, nWires: 60 }), 1.17137e-4, 1e-4);
  rel(armourResistance20({ material: 'steel', wireDiameter: 0.005, nWires: 60, layFactor: 1.05 }), 1.17137e-4 * 1.05, 1e-4);
  assert.equal(armourWireCount({ meanDiameter: 0.1, wireDiameter: 0.005, coverage: 0.95 }), 59);
});
