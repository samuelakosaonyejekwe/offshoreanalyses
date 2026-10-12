import test from 'node:test';
import assert from 'node:assert/strict';
import {
  API_SAND_TABLE, densityDescription, apiSandParameters, clayAlpha, unitResistanceAt, axialCapacity,
  capacityVsPenetration, requiredPenetration, cptShaftFriction, CPT_SHAFT_COEFFICIENTS, cptEndBearingUWA,
  groupBlockCapacityClay,
} from '../../src/engine/geotech/axialcapacity.js';
import { createProfile } from '../../src/engine/geotech/soilprofile.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-12), `${a} vs ${b}`);

test('API alpha: 0.5 psi^-0.5 (psi <= 1), 0.5 psi^-0.25 (psi > 1), capped at 1', () => {
  close(clayAlpha(1), 0.5);
  close(clayAlpha(0.49), 0.5 / 0.7);
  close(clayAlpha(0.25), 1);
  close(clayAlpha(0.1), 1);
  close(clayAlpha(16), 0.25);
});

test('API RP 2GEO Table 1 design parameters for cohesionless soil', () => {
  assert.deepEqual({ ...API_SAND_TABLE['medium-dense sand-silt'] }, { beta: 0.29, fLim: 67e3, Nq: 12, qLim: 3e6 });
  assert.deepEqual({ ...API_SAND_TABLE['medium-dense sand'] }, { beta: 0.37, fLim: 81e3, Nq: 20, qLim: 5e6 });
  assert.deepEqual({ ...API_SAND_TABLE['dense sand'] }, { beta: 0.46, fLim: 96e3, Nq: 40, qLim: 10e6 });
  assert.deepEqual({ ...API_SAND_TABLE['very dense sand'] }, { beta: 0.56, fLim: 115e3, Nq: 50, qLim: 12e6 });
  assert.deepEqual(API_SAND_TABLE['dense sand-silt'], API_SAND_TABLE['medium-dense sand']);
  assert.equal(densityDescription(0.5), 'medium-dense');
  assert.equal(densityDescription(0.7), 'dense');
  assert.equal(apiSandParameters({ type: 'sand', Dr: 0.9 }).beta, 0.56);
  assert.equal(apiSandParameters({ type: 'silt', Dr: 0.7 }).beta, 0.37);
  assert.throws(() => apiSandParameters({ type: 'sand', Dr: 0.2 }), RangeError); // loose sand: no API values
});

const twoLayer = () => createProfile([
  { type: 'clay', top: 0, thickness: 20, gammaSub: 8000, suTop: 50e3, suGrad: 0 },
  { type: 'sand', top: 20, thickness: 20, gammaSub: 10000, phiDeg: 33, apiSand: 'medium-dense sand' },
]);

test('hand calculation: 1.0 m pile, 30 m penetration, clay (su 50 kPa) over medium-dense sand', () => {
  const D = 1.0, t = 0.025, L = 30;
  // Clay, sigma'v = 8 z kPa, psi = 6.25 / z:
  //   z < 6.25 m (psi > 1): f = 25 (z/6.25)^0.25 kPa  → ∫ = 25 · 6.25 / 1.25 = 125 kPa·m
  //   6.25–20 m (psi <= 1, alpha < 1 since psi > 0.25): f = 0.5 sqrt(50 · 8 z) = 10 sqrt(z)
  //     → ∫ = (20/3)(20^1.5 − 6.25^1.5) = 492.118 kPa·m
  const clay = 125 + (20 / 3) * (20 ** 1.5 - 6.25 ** 1.5);
  // Sand, sigma'v = 160 + 10 (z − 20) kPa, f = 0.37 sigma'v <= 81 kPa: limit reached at sigma'v = 218.92 kPa
  const zl = (81 / 0.37 - 160) / 10; // 5.892 m into the sand
  const sand = 0.5 * (0.37 * 160 + 81) * zl + 81 * (10 - zl);
  const Qfo = Math.PI * D * (clay + sand) * 1e3;
  // End bearing at 30 m: Nq sigma'v = 20 × 260 kPa = 5.2 MPa > 5 MPa limit
  const QpPlug = (5e6 * Math.PI * D * D) / 4;
  const r = axialCapacity({ profile: twoLayer(), L, D, t });
  close(r.Qfo, Qfo, 2e-3);
  close(r.qTip, 5e6);
  close(r.QpPlugged, QpPlug);
  close(r.Qfi, r.Qfo * (D - 2 * t) / D, 1e-12);
  close(r.QpUnplugged, (5e6 * Math.PI * (D * D - 0.95 * 0.95)) / 4 + r.Qfi, 1e-12);
  assert.equal(r.plugged, true); // plug resistance exceeds gross end bearing
  close(r.Qc, r.Qfo + QpPlug, 1e-12);
  close(r.Qt, r.Qfo);
  // unit values at spot depths
  close(unitResistanceAt(twoLayer(), 10).f, 10e3 * Math.sqrt(10), 1e-12);
  close(unitResistanceAt(twoLayer(), 10).q, 450e3);
  close(unitResistanceAt(twoLayer(), 22).f, 0.37 * 180e3, 1e-12);
});

test('coring governs for a short pile in clay; options scale the components', () => {
  const p = twoLayer();
  const r = axialCapacity({ profile: p, L: 5, D: 2, t: 0.04 });
  assert.equal(r.plugged, false);
  close(r.Qp, (450e3 * Math.PI * (4 - 1.92 ** 2)) / 4 + r.Qfi, 1e-12);
  const half = axialCapacity({ profile: p, L: 5, D: 2, t: 0.04, internalFrictionFactor: 0.5, setupFactor: 1.2, tensionFactor: 0.7 });
  close(half.Qfo, 1.2 * r.Qfo, 1e-12);
  close(half.Qfi, 0.6 * r.Qfi, 1e-12);
  close(half.Qt, 0.7 * half.Qfo, 1e-12);
  assert.equal(axialCapacity({ profile: p, L: 5, D: 2, t: 0.04, plugMode: 'plugged' }).plugged, true);
  close(axialCapacity({ profile: p, L: 5, D: 2, t: 0.04, closedEnded: true }).Qfi, 0);
});

test('capacity curve is consistent and required penetration meets the factored load (WSD 1.5 / 2.0, LRFD factor)', () => {
  const p = twoLayer(), o = { profile: p, D: 1, t: 0.025 };
  const cv = capacityVsPenetration({ ...o, Lmin: 5, Lmax: 40, step: 5 });
  assert.deepEqual(cv.L, [5, 10, 15, 20, 25, 30, 35, 40]);
  close(cv.Qc[5], axialCapacity({ ...o, L: 30 }).Qc);
  for (let i = 1; i < cv.L.length; i++) assert.ok(cv.Qfo[i] > cv.Qfo[i - 1]);
  const load = 3e6;
  for (const sf of [1.5, 2.0]) {
    const r = requiredPenetration({ ...o, load, safetyFactor: sf, mode: 'tension', tol: 0.001 });
    assert.equal(r.found, true);
    close(axialCapacity({ ...o, L: r.L }).Qt / sf, load, 2e-3);
  }
  const a = requiredPenetration({ ...o, load, safetyFactor: 1.5, mode: 'tension' }).L;
  const b = requiredPenetration({ ...o, load, safetyFactor: 2.0, mode: 'tension' }).L;
  assert.ok(b > a);
  const lrfd = requiredPenetration({ ...o, load, resistanceFactor: 0.8, mode: 'tension', tol: 0.001 });
  close(axialCapacity({ ...o, L: lrfd.L }).Qt * 0.8, load, 2e-3);
  assert.equal(requiredPenetration({ ...o, load: 1e9 }).found, false);
});

test('unified CPT shaft friction: hand evaluation for the UWA-05 and simplified ICP-05 coefficient sets', () => {
  const D = 2, t = 0.05, L = 40, z = 30, qc = 20e6, sigmaV = 300e3, deltaCvDeg = 28;
  const Ar = 1 - (1.9 / 2) ** 2; // 0.0975
  const h = (L - z) / D; // 5
  // UWA-05 compression: f = 0.030 qc Ar^0.3 [max(h, 2)]^-0.5 tan(delta)
  const uwa = 0.03 * qc * Ar ** 0.3 * h ** -0.5 * Math.tan((28 * Math.PI) / 180);
  close(cptShaftFriction({ qc, sigmaV, z, L, D, t, deltaCvDeg, coeff: CPT_SHAFT_COEFFICIENTS['UWA-05 compression'] }), uwa, 1e-12);
  // ICP-05 compression: f = 0.023 qc (sigma/pa)^0.1 Ar^0.2 [max(h, 4 sqrt(Ar))]^-0.4 tan(delta)
  const icp = 0.023 * qc * 3 ** 0.1 * Ar ** 0.2 * Math.max(h, 4 * Math.sqrt(Ar)) ** -0.4 * Math.tan((28 * Math.PI) / 180);
  close(cptShaftFriction({ qc, sigmaV, z, L, D, t, deltaCvDeg, coeff: CPT_SHAFT_COEFFICIENTS['ICP-05 compression'] }), icp, 1e-12);
  // tension / compression ratios 0.022/0.030 and 0.016/0.023; friction fatigue: higher friction near the tip
  const ft = cptShaftFriction({ qc, sigmaV, z, L, D, t, deltaCvDeg, coeff: CPT_SHAFT_COEFFICIENTS['UWA-05 tension'] });
  close(ft / uwa, 0.022 / 0.03, 1e-12);
  const tip = cptShaftFriction({ qc, sigmaV, z: 39.9, L, D, t, deltaCvDeg, coeff: CPT_SHAFT_COEFFICIENTS['UWA-05 compression'] });
  close(tip, 0.03 * qc * Ar ** 0.3 * 2 ** -0.5 * Math.tan((28 * Math.PI) / 180), 1e-12);
  const eb = cptEndBearingUWA({ qcAvg: 30e6, D, t });
  close(eb.q, 30e6 * (0.15 + 0.45 * Ar));
});

test('group block failure in clay: perimeter shear plus base bearing versus sum of piles', () => {
  // 3 × 3 group, block 5 m × 5 m × 20 m, su 50 kPa: 2·10·20·50 + 9·50·25 = 20 000 + 11 250 kN
  const r = groupBlockCapacityClay({ B: 5, W: 5, L: 20, suShaft: 50e3, suBase: 50e3, nPiles: 9, Qsingle: 4e6 });
  close(r.Qblock, 31.25e6);
  close(r.Qsum, 36e6);
  assert.equal(r.governs, 'block');
  close(r.efficiency, 31.25 / 36);
  assert.equal(groupBlockCapacityClay({ B: 5, W: 5, L: 20, suShaft: 50e3, suBase: 50e3, nPiles: 9, Qsingle: 2e6 }).governs, 'individual');
});

test('loose sand: outside API Table 1 by default; allowLooseSand uses the lowest class and flags the result as screening', () => {
  const profile = createProfile([{ type: 'sand', top: 0, thickness: 40, gammaSub: 9000, phiDeg: 28, Dr: 0.25 }]);
  const pile = { profile, L: 20, D: 1.0, t: 0.025 };
  assert.throws(() => axialCapacity(pile), (e) => e instanceof RangeError && /loose sand/.test(e.message) && /allowLooseSand/.test(e.message));
  assert.throws(() => apiSandParameters({ type: 'sand', Dr: 0.1 }), RangeError);
  const p = apiSandParameters({ type: 'sand', Dr: 0.25 }, { allowLooseSand: true });
  assert.equal(p.beta, 0.29); assert.equal(p.fLim, 67e3); assert.equal(p.Nq, 12); assert.equal(p.qLim, 3e6);
  assert.equal(p.screening, true); assert.match(p.warning, /loose sand/); assert.match(p.category, /medium-dense sand-silt/);
  // hand calculation, lowest class: f = 0.29 σ′v = 0.29 × 9 z kPa = 2.61 z kPa (67 kPa limit not reached: 52.2 kPa at 20 m)
  //   Qfo = π D ∫ f dz = π × 1.0 × 2.61 × 20²/2 = 1639.9 kN; q = 12 × 180 = 2160 kPa (< 3 MPa)
  const c = axialCapacity({ ...pile, allowLooseSand: true });
  assert.equal(c.screening, true); assert.equal(c.warnings.length, 1); assert.match(c.warnings[0], /Screening value only/);
  close(c.Qfo, Math.PI * 1.0 * 2.61e3 * 200, 1e-9);
  close(c.qTip, 12 * 180e3, 1e-12);
  // medium-dense sand is unaffected by the option and carries no flag
  const md = createProfile([{ type: 'sand', top: 0, thickness: 40, gammaSub: 9000, phiDeg: 33, Dr: 0.5 }]);
  const a = axialCapacity({ ...pile, profile: md }); const b = axialCapacity({ ...pile, profile: md, allowLooseSand: true });
  assert.equal(a.screening, false); assert.deepEqual(a.warnings, []); close(b.Qc, a.Qc, 1e-15); assert.equal(b.screening, false);
  assert.ok(c.Qfo < a.Qfo, 'the substituted class is weaker than medium-dense sand');
  // the flag follows through the wrappers
  const curve = capacityVsPenetration({ profile, D: 1.0, t: 0.025, Lmin: 5, Lmax: 20, step: 5, allowLooseSand: true });
  assert.equal(curve.screening, true); assert.equal(curve.L.length, 4);
  const req = requiredPenetration({ profile, D: 1.0, t: 0.025, load: 1e6, allowLooseSand: true });
  assert.ok(req.found && req.screening && req.warnings.length === 1);
  assert.equal(unitResistanceAt(profile, 10, { allowLooseSand: true }).screening, true);
  assert.equal(apiSandParameters({ type: 'sand', Dr: 0.5 }, { allowLooseSand: true }).screening, undefined);
});
