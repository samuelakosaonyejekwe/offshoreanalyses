import { test } from 'node:test';
import assert from 'node:assert/strict';
import { classifySeabed, reduceSeabedSamples, SEABED_CLASSES } from '../src/data/site.js';

test('seabed descriptions: Folk class names', () => {
  for (const [t, c] of [['Sand', 'sand'], ['Mud', 'mud'], ['sandy Mud', 'sandyMud'], ['muddy Sand', 'muddySand'], ['sandy Gravel', 'gravel'], ['Mud to muddy Sand', 'sandyMud'],
    ['Coarse-grained sediment', 'gravel'], ['Mixed sediment', 'mixed'], ['Rock & boulders', 'rock'], ['gravelly muddy Sand', 'muddySand'], ['gravelly Sand', 'sand']]) assert.equal(classifySeabed(t), c, t);
});
test('seabed descriptions: chart abbreviations take the first constituent', () => {
  for (const [t, c] of [['SD,BRN,CRS,BRKN SHELLS', 'sand'], ['M', 'mud'], ['Cy,S', 'mud'], ['RK', 'rock'], ['SH,SD', 'biogenic'], ['G,S', 'gravel'], ['BRN,FNE,SD', 'sand'], ['hard', 'rock']]) assert.equal(classifySeabed(t), c, t);
  for (const t of ['', null, 'No data at this level', 'BRN,FNE']) assert.equal(classifySeabed(t), null);
});
test('seabed samples: nearest described samples decide, far ones are ignored', () => {
  const at = (dLat, text, year = 1990) => ({ lat: 30 + dLat, lon: -80, text, year });
  const r = reduceSeabedSamples([at(0.01, 'SD,GY'), at(0.02, 'SD'), at(0.03, 'M'), at(0.5, 'RK'), at(0.04, 'BRN')], 30, -80);
  assert.equal(r.cls, 'sand'); assert.equal(r.n, 3); assert.equal(r.agree, 2);
  assert.ok(Math.abs(r.nearestKm - 1.112) < 0.01);
  assert.equal(reduceSeabedSamples([at(0.5, 'RK')], 30, -80), null);
});
test('every class has a label, a meaning and a grain size in the range of its name', () => {
  for (const [k, c] of Object.entries(SEABED_CLASSES)) { assert.ok(c.label && c.means.length > 40, k); if (c.d50 != null) assert.ok(c.d50 > 1e-6 && c.d50 < 0.1, k); }
  assert.ok(SEABED_CLASSES.mud.d50 < 6.3e-5 && SEABED_CLASSES.sand.d50 > 6.3e-5 && SEABED_CLASSES.sand.d50 < 2e-3 && SEABED_CLASSES.gravel.d50 > 2e-3);
});
