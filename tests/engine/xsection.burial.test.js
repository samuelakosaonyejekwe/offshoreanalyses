import test from 'node:test';
import assert from 'node:assert/strict';
import {
  BPI_CLASSES, burialProtectionIndex, anchorPenetration, cbraStrikeProbability, depthOfLowering,
  GENERIC_FISHING_PENETRATION, fishingGearPenetration, trenchingRate, ploughTowForce, naturalInfillTime, rockBerm,
  rockSizeShields, mattressCount, burialThermalRecord, sandwaveExposure,
} from '../../src/engine/xsection/burial.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);

test('burial protection index is depth over the BPI = 1 reference depth', () => {
  const r = burialProtectionIndex({ referenceDepth: 0.6, burialDepth: 1.5, targetBpi: 3 });
  close(r.bpi, 2.5); close(r.requiredDepth, 1.8);
  assert.equal(r.protectedAgainst, BPI_CLASSES[1].description);
  assert.equal(burialProtectionIndex({ referenceDepth: 0.6, burialDepth: 0.3 }).protectedAgainst, 'Below BPI 1');
});

test('anchor penetration: fluke geometry and cube-root scaling', () => {
  close(anchorPenetration({ flukeLength: 1.6, flukeAngleDeg: 40, soilFactor: 1.5 }).penetration, 1.6 * Math.sin(40 * Math.PI / 180) * 1.5);
  const a = anchorPenetration({ anchorMass: 8000, flukeCoefficient: 0.08 });
  close(a.flukeLength, 0.08 * 20); close(a.penetration, 1.6 * Math.SQRT1_2);
  // geometric similarity: 8x the mass doubles the fluke length
  close(anchorPenetration({ anchorMass: 64000, flukeCoefficient: 0.08 }).flukeLength, 2 * a.flukeLength, 1e-12);
});

test('CBRA strike probability and return period by hand', () => {
  const classes = [
    // 2000 crossings/yr, 500 m hazard zone at 5 m/s -> 100 s each -> 55.56 h/yr; x 1e-2 /8760
    { name: 'small', crossings: 2000, speed: 5, hazardDistance: 500, pIncident: 1e-2, penetration: 0.8 },
    { name: 'large', crossings: 300, speed: 6, hazardDistance: 1200, pIncident: 2e-2, penetration: 2.0 },
  ];
  const r = cbraStrikeProbability({ classes, burialDepth: 0, pTraffic: 1, pWd: 0.5 });
  const h1 = 2000 * 100 / 3600, h2 = 300 * 200 / 3600;
  close(r.classes[0].exposureHours, h1, 1e-12);
  const p = 0.5 * (h1 / 8760 * 1e-2 + h2 / 8760 * 2e-2);
  close(r.pStrike, p, 1e-12);
  close(r.returnPeriod, 1 / p, 1e-12);
  // burial at 1 m removes the small class only
  const b = cbraStrikeProbability({ classes, burialDepth: 1.0, pWd: 0.5 });
  close(b.pStrike, 0.5 * h2 / 8760 * 2e-2, 1e-12);
  assert.equal(b.classes[0].threatens, false);
  assert.equal(cbraStrikeProbability({ classes, burialDepth: 2.0 }).returnPeriod, Infinity);
  // depth of lowering: target between the two levels -> bury below the small anchors
  const d = depthOfLowering({ classes, targetProbability: 2.5e-5, pWd: 0.5, margin: 0.2, fishingPenetration: 0.3 });
  close(d.recommendedDepth, 1.0); assert.equal(d.governing, 'anchor');
  close(d.pStrike, b.pStrike, 1e-12); close(d.pStrikeUnburied, p, 1e-12);
  // very strict target: below the deepest anchor
  close(depthOfLowering({ classes, targetProbability: 1e-9, pWd: 0.5 }).recommendedDepth, 2.0);
  // lenient target: fishing gear governs
  const f = depthOfLowering({ classes, targetProbability: 1, fishingPenetration: 0.3, margin: 0.2 });
  close(f.recommendedDepth, 0.5); assert.equal(f.governing, 'fishing');
});

test('fishing gear table is generic and editable', () => {
  assert.equal(GENERIC_FISHING_PENETRATION.generic, true);
  const r = fishingGearPenetration({ gear: 'beamTrawl', soil: 'mud', safetyFactor: 2 });
  close(r.design, 2 * GENERIC_FISHING_PENETRATION.beamTrawl.mud); assert.equal(r.generic, true);
  close(fishingGearPenetration({ gear: 'x', soil: 'y', table: { x: { y: 0.42 } } }).penetration, 0.42);
  assert.throws(() => fishingGearPenetration({ gear: 'nope', soil: 'sand' }), /gear/);
});

test('trenching rate from the energy balance; plough tow force forms by hand', () => {
  // 600 kW at 50 % into a 0.4 m x 2 m trench with e = 10 su, su = 40 kPa: v = 3e5/(4e5*0.8) = 0.9375 m/s
  const t = trenchingRate({ power: 600e3, efficiency: 0.5, strength: 40e3, strengthFactor: 10, trenchWidth: 0.4, trenchDepth: 2 });
  close(t.rate, 0.9375); close(t.ratePerHour, 3375);
  close(trenchingRate({ power: 600e3, specificEnergy: 4e5, trenchArea: 0.8 }).rate, 1.875);
  const s = ploughTowForce({ soil: 'sand', depth: 1.5, ploughWeight: 300e3, gammaSub: 10e3, speed: 0.1, Cw: 0.4, Cs: 10, Cd: 2e5, phiDeg: 35, shareWidth: 1.0 });
  close(s.staticPart, 0.4 * 300e3 + 10 * 10e3 * 1.5 ** 3); close(s.ratePart, 2e5 * 0.1 * 2.25); close(s.empirical, s.staticPart + s.ratePart);
  const Kp = Math.tan(62.5 * Math.PI / 180) ** 2;
  close(Kp, (1 + Math.sin(35 * Math.PI / 180)) / (1 - Math.sin(35 * Math.PI / 180)), 1e-12);
  close(s.rankine, 0.5 * Kp * 10e3 * 2.25 + 0.4 * 300e3, 1e-12);
  const c = ploughTowForce({ soil: 'clay', depth: 1.5, ploughWeight: 300e3, gammaSub: 8e3, su: 30e3, Cw: 0.4, Cc: 15, shareWidth: 1.0 });
  close(c.empirical, 120e3 + 15 * 30e3 * 2.25);
  close(c.rankine, (0.5 * 8e3 * 2.25 + 2 * 30e3 * 1.5) + 120e3);
  assert.throws(() => ploughTowForce({ soil: 'sand', depth: 1, ploughWeight: 1, gammaSub: 1e4 }), RangeError);
  close(naturalInfillTime({ trenchArea: 1.2, transportRate: 1e-6, trapEfficiency: 0.5 }).timeDays, 2.4e6 / 86400);
});

test('rock berm geometry and tonnage take-off by hand', () => {
  // D = 0.3 m, cover 0.7 m -> h = 1.0 m; crest 1 m, slopes 1:3 -> A = 1 x (1 + 3) = 4 m2, base 7 m
  const r = rockBerm({ productDiameter: 0.3, cover: 0.7, crestWidth: 1, slope: 3, length: 500, rockDensity: 2650, porosity: 0.4, overdump: 0.1, requiredAbsorptionCover: 0.5 });
  close(r.height, 1.0); close(r.baseWidth, 7.0);
  close(r.area, 4 - Math.PI * 0.09 / 4);
  close(r.tonnes, (4 - Math.PI * 0.0225) * 500 * 1.1 * 1590 / 1000);
  close(r.coverUtilisation, 0.5 / 0.7);
});

test('Shields rock size, mattress count and thermal record', () => {
  // tau = 1025 * 0.1^2 = 10.25 Pa; D50 = 10.25/(0.035 * 1625 * 9.80665) = 18.38 mm
  const s = rockSizeShields({ shearVelocity: 0.1, shieldsCritical: 0.035, d50: 0.05 });
  close(s.bedShear, 10.25); close(s.requiredD50, 10.25 / (0.035 * 1625 * 9.80665)); close(s.utilisation, s.requiredD50 / 0.05);
  const m = mattressCount({ length: 100, mattressLength: 6, overlap: -0.2, endOverrun: 2, coverWidth: 5, mattressWidth: 3, mattressMass: 9000 });
  assert.equal(m.perRow, Math.ceil(102 / 6.2)); assert.equal(m.rows, 2); assert.equal(m.count, 34); close(m.totalMass, 34 * 9000);
  assert.equal(mattressCount({ length: 60, mattressLength: 6 }).count, 10);
  const t = burialThermalRecord({ depthToTop: 1.5, outerDiameter: 0.2, soilResistivity: 0.7 });
  close(t.burialDepthAxis, 1.6); close(t.u, 16);
  close(t.T4, 0.7 / (2 * Math.PI) * Math.acosh(16), 1e-12);
});

test('sandwave exposure: crest case in closed form, cover history reproduced, never-exposed case', () => {
  const H = 4, lam = 200, c = 10, d0 = 1.5;
  const r = sandwaveExposure({ height: H, wavelength: lam, migrationRate: c, initialCover: d0, targetCover: 0.5 });
  close(r.time, lam / (2 * Math.PI * c) * Math.acos(1 - 2 * 1.0 / H), 1e-12);
  close(r.minimumCover, d0 - H);
  close(r.period, 20);
  // general phase: evaluate the cover history at the returned time and check no earlier crossing
  const phase = 2.0, target = 0.6;
  const g = sandwaveExposure({ height: H, wavelength: lam, migrationRate: c, initialCover: d0, targetCover: target, phase });
  const cover = (t) => d0 + H / 2 * (Math.cos(phase - 2 * Math.PI * c * t / lam) - Math.cos(phase));
  close(cover(g.time), target, 1e-9);
  for (let t = 0; t < g.time * 0.999; t += g.time / 500) assert.ok(cover(t) > target);
  // phase 2 rad is on the lee side: the cover first increases as the crest approaches
  assert.ok(cover(g.time / 10) > d0);
  // cable under the trough with cover larger than nothing to lose
  const n = sandwaveExposure({ height: H, wavelength: lam, migrationRate: c, initialCover: d0, targetCover: 0.5, phase: Math.PI });
  assert.equal(n.time, Infinity); assert.equal(n.exposed, false); close(n.minimumCover, d0);
  assert.equal(sandwaveExposure({ height: H, wavelength: lam, migrationRate: c, initialCover: 0.4, targetCover: 0.5 }).time, 0);
});
