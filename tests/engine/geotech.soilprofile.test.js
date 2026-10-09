import test from 'node:test';
import assert from 'node:assert/strict';
import {
  createProfile, validateLayers, verticalEffectiveStress, propertiesAt, layerIndexAt, applyScour, scourBaseDepth,
  sampleProfile, EXAMPLE_PROFILES, exampleProfile,
} from '../../src/engine/geotech/soilprofile.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-12), `${a} vs ${b}`);

const three = () => createProfile([
  { type: 'sand', top: 0, thickness: 5, gammaSub: 9000, phiDeg: 32 },
  { type: 'clay', top: 5, thickness: 10, gammaSub: 7000, suTop: 30e3, suGrad: 2000, eps50: 0.01 },
  { type: 'sand', top: 15, thickness: 25, gammaSub: 10000, phiDeg: 36, Vs: 250 },
]);

test('vertical effective stress: hand integration through three layers', () => {
  const p = three();
  // 0–5 m: 9 kPa/m; 5–15 m: 7 kPa/m; below: 10 kPa/m
  close(verticalEffectiveStress(p, 0), 0);
  close(verticalEffectiveStress(p, 3), 27e3);
  close(verticalEffectiveStress(p, 5), 45e3);
  close(verticalEffectiveStress(p, 12), 45e3 + 7 * 7e3);
  close(verticalEffectiveStress(p, 15), 115e3);
  close(verticalEffectiveStress(p, 40), 115e3 + 250e3);
});

test('property lookup: linear su, interface belongs to the lower layer, G0 = rho Vs^2', () => {
  const p = three();
  assert.equal(layerIndexAt(p, 5), 1);
  assert.equal(layerIndexAt(p, 4.999), 0);
  assert.equal(layerIndexAt(p, 40), 2);
  close(propertiesAt(p, 9).su, 30e3 + 4 * 2000);
  assert.equal(propertiesAt(p, 2).su, undefined);
  // saturated density (10000 + 10052)/9.80665 = 2044.7 kg/m3 → G0 = 2044.7 × 250² = 127.8 MPa
  close(propertiesAt(p, 20).G0, ((10000 + 10052) / 9.80665) * 250 * 250);
  assert.throws(() => propertiesAt(p, 41), RangeError);
});

test('scour: zero stress in the scour hole, linear overburden reduction to 6 D, general-scour stress below', () => {
  const p0 = createProfile([{ type: 'sand', top: 0, thickness: 50, gammaSub: 10e3, phiDeg: 35 }]);
  const p = applyScour(p0, { general: 1, D: 2 }); // local 1.5 D = 3 m, reduction depth 6 D = 12 m
  assert.equal(scourBaseDepth(p), 4);
  close(verticalEffectiveStress(p, 4), 0);
  close(verticalEffectiveStress(p, 2), 0);
  close(verticalEffectiveStress(p, 13), 120e3); // 12 m below the general scour level
  close(verticalEffectiveStress(p, 8.5), 60e3); // halfway between the hole base (4 m) and 13 m
  close(verticalEffectiveStress(p, 20), 190e3);
  close(verticalEffectiveStress(p0, 20), 200e3); // original profile untouched
  assert.equal(propertiesAt(p, 10).zBelowSurface, 6);
  assert.throws(() => applyScour(p0, { general: 60, local: 0, reductionDepth: 0 }), RangeError);
});

test('validation reports gaps, missing strength and bad values', () => {
  const v = validateLayers([
    { type: 'clay', top: 0, thickness: 5, gammaSub: 8000 },
    { type: 'sand', top: 6, thickness: 5, gammaSub: -1, phiDeg: 33 },
    { type: 'peat', top: 11, thickness: 1, gammaSub: 5000 },
  ]);
  assert.equal(v.valid, false);
  assert.ok(v.errors.some((e) => e.includes('suTop')));
  assert.ok(v.errors.some((e) => e.includes('contiguous')));
  assert.ok(v.errors.some((e) => e.includes('gammaSub')));
  assert.ok(v.errors.some((e) => e.includes('type')));
  assert.throws(() => createProfile([]), RangeError);
});

test('example profiles are labelled generic and build without error', () => {
  for (const key of Object.keys(EXAMPLE_PROFILES)) {
    assert.match(EXAMPLE_PROFILES[key].label, /^Generic example/);
    const p = exampleProfile(key);
    assert.ok(p.depth >= 50);
    const s = sampleProfile(p, { dz: 5 });
    for (let i = 1; i < s.z.length; i++) assert.ok(s.sigmaV[i] >= s.sigmaV[i - 1]);
  }
  assert.throws(() => exampleProfile('nope'), RangeError);
});
