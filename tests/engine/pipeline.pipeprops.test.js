import test from 'node:test';
import assert from 'node:assert/strict';
import {
  pipeSection, pipeWeights, steelGrade, listSteelGrades, temperatureDerating, fabricationFactor, materialStrengthFactor,
  characteristicStrength, npsOuterDiameter, standardPipeTable, scheduleWallThickness,
} from '../../src/engine/pipeline/pipeprops.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${a} vs ${b}`);

test('pipeSection: hand calculation for 323.9 x 12.7 mm', () => {
  // Di = 0.3239 - 2*0.0127 = 0.2985
  // As = pi/4 (0.3239^2 - 0.2985^2) = pi/4 (0.10491121 - 0.08910225) = pi/4 * 0.01580896 = 0.0124163 m2
  // Ai = pi/4 * 0.08910225 = 0.0699807 m2
  // I  = pi/64 (0.3239^4 - 0.2985^4) = pi/64 (0.011006362 - 0.007939211) = 1.50558e-4 m4
  const s = pipeSection({ D: 0.3239, t: 0.0127 });
  close(s.Di, 0.2985, 1e-12);
  close(s.As, 0.0124163, 2e-5);
  close(s.Ai, 0.0699807, 2e-5);
  close(s.I, 1.50558e-4, 2e-4);
  close(s.As + s.Ai, s.Ae, 1e-12);
  close(s.Z, s.I / (0.3239 / 2), 1e-12);
  // thin-wall check: As ~ pi (D - t) t exactly for an annulus
  close(s.As, Math.PI * (0.3239 - 0.0127) * 0.0127, 1e-12);
  assert.throws(() => pipeSection({ D: 0.3, t: 0.15 }), RangeError);
  assert.throws(() => pipeSection({ D: -1, t: 0.01 }), RangeError);
});

test('pipeWeights: submerged weight of a coated pipe, hand calculation', () => {
  // 323.9 x 12.7 mm steel, 3 mm coating at 1300 kg/m3, 40 mm concrete at 3040 kg/m3 with 3 % absorption,
  // contents 800 kg/m3, seawater 1025 kg/m3.
  const g = 9.80665;
  const steel = 7850 * (Math.PI / 4) * (0.3239 ** 2 - 0.2985 ** 2);            // 97.468 kg/m
  const coat = 1300 * (Math.PI / 4) * (0.3299 ** 2 - 0.3239 ** 2);             // 4.005 kg/m
  const conc = 3040 * (Math.PI / 4) * (0.4099 ** 2 - 0.3299 ** 2);             // 141.308 kg/m
  const absorbed = 0.03 * conc;                                                 // 4.239 kg/m
  const contents = 800 * (Math.PI / 4) * 0.2985 ** 2;                           // 55.985 kg/m
  const displaced = 1025 * (Math.PI / 4) * 0.4099 ** 2;                         // 135.260 kg/m
  const w = pipeWeights({ D: 0.3239, t: 0.0127, coatings: [{ thickness: 0.003, density: 1300 }, { thickness: 0.04, density: 3040, absorption: 0.03 }], rhoContents: 800 });
  close(w.Dtotal, 0.4099, 1e-12);
  close(w.massSteel, 97.468, 1e-4);
  close(w.massLayers[1], 141.308, 1e-4);
  close(w.buoyancy, 135.260 * g, 1e-4);
  close(w.weightDry, (steel + coat + conc) * g, 1e-12);
  close(w.submergedEmpty, (steel + coat + conc + absorbed - displaced) * g, 1e-12);
  close(w.submergedOperating, (steel + coat + conc + absorbed + contents - displaced) * g, 1e-12);
  close(w.submergedFlooded - w.submergedEmpty, 1025 * (Math.PI / 4) * 0.2985 ** 2 * g, 1e-10);
  // literal from the arithmetic above: (97.468 + 4.005 + 141.308 + 4.239 - 135.260) * 9.80665 = 1096.0 N/m
  close(w.submergedEmpty, 1096.0, 2e-4);
  close(w.sgEmpty, (w.submergedEmpty + w.buoyancy) / w.buoyancy, 1e-12);
});

test('pipeWeights: bare pipe floats when lighter than the displaced water; marine growth adds diameter', () => {
  const w = pipeWeights({ D: 0.5, t: 0.005 });
  assert.ok(w.submergedEmpty < 0 && w.sgEmpty < 1);
  const mg = pipeWeights({ D: 0.5, t: 0.02, marineGrowth: { thickness: 0.05, density: 1325 } });
  close(mg.Dhydro, 0.6, 1e-12);
  close(mg.massMarineGrowth, 1325 * (Math.PI / 4) * (0.36 - 0.25), 1e-12);
  const no = pipeWeights({ D: 0.5, t: 0.02 });
  close(mg.submergedEmpty - no.submergedEmpty, (1325 - 1025) * (Math.PI / 4) * (0.36 - 0.25) * 9.80665, 1e-9);
});

test('steel grades: API 5L PSL2 and DNV-ST-F101 SMYS / SMTS', () => {
  assert.deepEqual([steelGrade('X52').smys, steelGrade('X52').smts], [360e6, 460e6]);
  assert.deepEqual([steelGrade('X65').smys, steelGrade('X65').smts], [450e6, 535e6]);
  assert.deepEqual([steelGrade('x70').smys, steelGrade('X70').smts], [485e6, 570e6]);
  assert.deepEqual([steelGrade('X80').smys, steelGrade('X80').smts], [555e6, 625e6]);
  assert.deepEqual([steelGrade('DNV 450').smys, steelGrade('L450').smts], [450e6, 535e6]);
  assert.equal(steelGrade('DNV 245').smts, 370e6);
  assert.ok(listSteelGrades().includes('X60'));
  assert.throws(() => steelGrade('X999'), RangeError);
});

test('temperature de-rating: 0 to 50 C, 30 MPa at 100 C, 70 MPa at 200 C, linear between', () => {
  assert.equal(temperatureDerating(20), 0);
  assert.equal(temperatureDerating(50), 0);
  close(temperatureDerating(75), 15e6);
  close(temperatureDerating(100), 30e6);
  close(temperatureDerating(150), 50e6);
  assert.throws(() => temperatureDerating(250), RangeError);
  close(temperatureDerating(60, { curve: [[20, 0], [100, 80e6]] }), 40e6);
});

test('characteristic strength fy = (SMYS - fy,temp) alphaU', () => {
  const c = characteristicStrength({ smys: 450e6, smts: 535e6, tempC: 100 });
  close(c.fy, (450e6 - 30e6) * 0.96);
  close(c.fu, (535e6 - 30e6) * 0.96);
  assert.equal(materialStrengthFactor(true), 1);
  assert.equal(fabricationFactor('UOE'), 0.85);
  assert.equal(fabricationFactor('seamless'), 1);
  assert.equal(fabricationFactor('erw'), 0.93);
});

test('NPS table: outside diameters are exact (ASME B36.10M)', () => {
  close(npsOuterDiameter(12), 12.75 * 0.0254, 1e-14);
  close(npsOuterDiameter('6'), 0.168275, 1e-12);
  close(npsOuterDiameter(24), 0.6096, 1e-12);
  close(npsOuterDiameter('1 1/2'), 1.9 * 0.0254, 1e-14);
  close(npsOuterDiameter(0.5), 0.021336, 1e-12);
  close(scheduleWallThickness(8, 40), 0.322 * 0.0254, 1e-14);
  close(scheduleWallThickness(12, 'XS'), 0.0127, 1e-12);
  const tab = standardPipeTable({ sizes: [10] });
  assert.equal(tab.length, 1);
  close(tab[0].D, 0.27305, 1e-12);
  assert.ok(tab[0].walls.some((x) => x.schedule === '80' && Math.abs(x.t - 0.594 * 0.0254) < 1e-12));
  assert.ok(standardPipeTable().length >= 30);
  assert.throws(() => npsOuterDiameter(13), RangeError);
});
