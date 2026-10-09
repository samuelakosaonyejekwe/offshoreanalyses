import test from 'node:test';
import assert from 'node:assert/strict';
import {
  UMBILICAL_UTILISATION, umbilicalComponent, assembleUmbilical, umbilicalLoadSharing, umbilicalUtilisation,
  umbilicalCapacityCurve, umbilicalMbr, ringRadialLoadMoment, umbilicalCrushCheck, umbilicalCoreTemperature,
  umbilicalTubeFatigue, tubePressureDrop, hydraulicResponseTime,
} from '../../src/engine/xsection/umbilical.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);
const DEG = Math.PI / 180;
const tube = { kind: 'tube', name: 'SD tube', od: 0.0191, id: 0.0161, E: 2.0e11, rho: 7800, yieldStress: 550e6, fluidDensity: 1000 };
const wire = (name) => ({ kind: 'wire', name, od: 0.004, E: 2.0e11, rho: 7850, yieldStress: 1200e6 });
const tubeA = Math.PI * (0.0191 ** 2 - 0.0161 ** 2) / 4;

const demo = () => assembleUmbilical({
  centre: { ...tube, name: 'centre tube' },
  layers: [
    { R: 0.02, alphaDeg: 8, lay: 1, elements: [{ component: tube, count: 6 }] },
    { R: 0.035, alphaDeg: 15, lay: -1, elements: [{ component: wire('armour 1'), count: 50 }] },
    { R: 0.039, alphaDeg: 15, lay: 1, elements: [{ component: wire('armour 2'), count: 56 }] },
  ],
  sheaths: [{ ri: 0.041, ro: 0.046, E: 8e8, nu: 0.42, rho: 950 }],
});

test('assembly: mass, stiffness and weight by hand sums', () => {
  const s = demo();
  const wA = Math.PI * 0.004 ** 2 / 4, shA = Math.PI * (0.046 ** 2 - 0.041 ** 2);
  const c8 = Math.cos(8 * DEG), c15 = Math.cos(15 * DEG), s8 = Math.sin(8 * DEG), s15 = Math.sin(15 * DEG);
  close(s.massEmpty, 7800 * tubeA * (1 + 6 / c8) + 7850 * wA * 106 / c15 + 950 * shA, 1e-12);
  const bore = Math.PI * 0.0161 ** 2 / 4;
  close(s.massFilled - s.massEmpty, 1000 * bore * (1 + 6 / c8), 1e-12);
  close(s.outerDiameter, 0.092, 1e-12);
  close(s.submergedFilled, (s.massFilled - 1025 * Math.PI * 0.046 ** 2) * 9.80665, 1e-12);
  const EAt = 2e11 * tubeA, EAw = 2e11 * wA;
  close(s.K11, EAt * (1 + 6 * c8 ** 3) + EAw * 106 * c15 ** 3 + 8e8 * shA, 1e-12);
  close(s.K12, 6 * EAt * 0.02 * s8 * c8 ** 2 + EAw * s15 * c15 ** 2 * (-50 * 0.035 + 56 * 0.039), 1e-10);
  assert.ok(s.EIstick > s.EIslip && s.EIslip > 0);
  assert.ok(s.EA < s.EAfixed);
  assert.equal(UMBILICAL_UTILISATION.normal, 0.67);
  assert.throws(() => umbilicalComponent({ od: 0.01, id: 0.02 }), /id/);
});

test('cross-wound symmetric armour pair: zero residual torque per unit tension', () => {
  const mk = (lay2) => assembleUmbilical({ layers: [
    { R: 0.03, alphaDeg: 12, lay: 1, elements: [{ component: wire('a'), count: 40 }] },
    { R: 0.03, alphaDeg: 12, lay: lay2, elements: [{ component: wire('b'), count: 40 }] },
  ] });
  const bal = mk(-1), unbal = mk(1);
  assert.equal(bal.torquePerTension, 0);
  assert.equal(bal.torqueBalanced, true);
  assert.equal(unbal.torqueBalanced, false);
  close(unbal.torquePerTension, 0.03 * Math.tan(12 * DEG), 1e-12); // K12/K11 = R tan(alpha)
  // torque-free unbalanced bundle is softer axially; balanced one is not
  close(bal.EA, bal.EAfixed, 1e-12);
  assert.ok(unbal.EA < 0.2 * unbal.EAfixed);
});

test('load sharing: shares sum to one, element forces and contact line loads by hand', () => {
  const s = demo();
  const T = 2e5;
  const r = umbilicalLoadSharing(s, { tension: T });
  close(r.axialStrain, T / s.K11, 1e-12);
  // the polymer sheath takes the (small) remainder
  const sum = r.elements.reduce((a, e) => a + e.share, 0);
  const sheathShare = 8e8 * Math.PI * (0.046 ** 2 - 0.041 ** 2) / s.K11;
  close(sum + sheathShare, 1, 1e-12);
  const el = r.elements.find((e) => e.name === 'SD tube');
  const c8 = Math.cos(8 * DEG), s8 = Math.sin(8 * DEG), c15 = Math.cos(15 * DEG), s15 = Math.sin(15 * DEG);
  close(el.force, 2e11 * tubeA * c8 ** 2 * r.axialStrain, 1e-12);
  close(el.lineLoad, el.force * s8 ** 2 / 0.02, 1e-12);
  // radial load of the two armour layers per metre of umbilical, shared by the six tubes
  const a1 = r.elements.find((e) => e.name === 'armour 1'), a2 = r.elements.find((e) => e.name === 'armour 2');
  const outerRadial = 50 * a1.force * s15 ** 2 / 0.035 / c15 + 56 * a2.force * s15 ** 2 / 0.039 / c15;
  close(el.contactOuter, outerRadial * c8 / 6, 1e-12);
  close(el.contactInner, el.lineLoad + el.contactOuter, 1e-12);
  assert.equal(a2.contactOuter, 0);
  close(r.torque, s.K12 * r.axialStrain, 1e-12);
  // torque-free: twist such that the net torque vanishes
  const f = umbilicalLoadSharing(s, { tension: T, twistRestrained: false });
  close(s.K12 * f.axialStrain + s.K22 * f.twist, 0, 0, 1e-6);
  close(s.K11 * f.axialStrain + s.K12 * f.twist, T, 1e-10);
});

test('capacity curve of a single straight tube: T = (eta sy - E r kappa) A, kappa limit = eta sy/(E r)', () => {
  const s = assembleUmbilical({ centre: tube, layers: [] });
  const eta = 0.67, A = tubeA, r = 0.0191 / 2;
  const cv = umbilicalCapacityCurve(s, { eta, nPoints: 5 });
  close(cv.tensionLimit, eta * 550e6 * A, 1e-8);
  close(cv.kappaLimit, eta * 550e6 / (2e11 * r), 1e-8);
  for (let i = 0; i < 4; i++) close(cv.allowableTension[i], (eta * 550e6 - 2e11 * r * cv.curvature[i]) * A, 1e-7);
  close(umbilicalMbr(s, { eta }).mbr, 2e11 * r / (eta * 550e6), 1e-8);
  // with tension the MBR grows
  const m = umbilicalMbr(s, { eta, tension: 0.5 * cv.tensionLimit });
  close(m.mbr, 2 * 2e11 * r / (eta * 550e6), 1e-7);
});

test('capacity of the bundle: pressure, friction and tension reduce the allowable curvature', () => {
  const s = demo();
  const base = umbilicalUtilisation(s, { tension: 1e5, curvature: 0.05, mu: 0.2 });
  assert.ok(base.elements.every((e) => e.frictionStress >= 0 && e.vonMises > 0));
  const noFric = umbilicalUtilisation(s, { tension: 1e5, curvature: 0.05, mu: 0 });
  assert.ok(base.utilisation >= noFric.utilisation);
  const cv = umbilicalCapacityCurve(s, { eta: 0.8, nPoints: 6 });
  for (let i = 1; i < cv.curvature.length; i++) assert.ok(cv.allowableTension[i] <= cv.allowableTension[i - 1] + 1e-6);
  close(umbilicalUtilisation(s, { tension: cv.allowableTension[2], curvature: cv.curvature[2], eta: 0.8 }).utilisation, 1, 1e-6);
  // a pressurised tube alone: zero-tension utilisation is the Lame von Mises stress at the bore
  const pt = assembleUmbilical({ centre: { ...tube, pressure: 69e6 }, layers: [] });
  const ri = 0.00805, ro = 0.00955, p = 69e6;
  const hoop = p * (ro ** 2 + ri ** 2) / (ro ** 2 - ri ** 2), ax = p * ri ** 2 / (ro ** 2 - ri ** 2);
  const vm = Math.sqrt(0.5 * ((hoop + p) ** 2 + (hoop - ax) ** 2 + (ax + p) ** 2));
  close(umbilicalUtilisation(pt, { tension: 0, curvature: 0, eta: 1 }).utilisation, vm / 550e6, 1e-10);
});

test('ring under equally spaced radial loads: M = 0.3183 P R for two opposing loads (statically indeterminate ring solved by quadrature)', () => {
  const P = 1000, R = 0.05;
  for (const n of [2, 3, 4]) {
    const th0 = Math.PI / n;
    // cut at a load point: shear P/2, hoop force N0 from zero shear at mid-span, M0 from zero rotation there:
    // M(th) = M0 + N0 R (1 - cos th) - (P/2) R sin th, integral over [0, th0] must vanish
    const N0 = P / 2 / Math.tan(th0);
    let int = 0; const m = 200000;
    for (let i = 0; i < m; i++) { const th = (i + 0.5) * th0 / m; int += (N0 * R * (1 - Math.cos(th)) - P / 2 * R * Math.sin(th)) * th0 / m; }
    const M0 = -int / th0;
    close(ringRadialLoadMoment({ P, R, nLoads: n }).momentAtLoad, M0, 1e-8);
  }
  close(ringRadialLoadMoment({ P, R }).coefficient, 1 / Math.PI, 1e-14);
  close(ringRadialLoadMoment({ P, R }).coefficient, 0.3183, 1e-4);
  close(ringRadialLoadMoment({ P, R }).momentBetween, -P * R * (0.5 - 1 / Math.PI), 1e-12); // -0.1817 P R
});

test('crush: tube line-load capacity, tensioner grip length and chute radius by hand', () => {
  const D = 0.0191, t = 0.0015, Rm = (D - t) / 2, sa = 0.8 * 550e6;
  const r = umbilicalCrushCheck({ tubeOd: D, tubeWall: t, allowableStress: sa, padLoad: 60e3, nTracks: 4, elementsUnderPad: 2, tension: 300e3, mu: 0.1, safetyFactor: 1.5, chuteRadius: 3 });
  const cap = sa * t * t / (6 * Rm / Math.PI);
  close(r.tubeLineLoadCapacity, cap, 1e-12);
  close(r.utilisationTensioner, 30e3 / cap, 1e-12);
  close(r.requiredTensionerLength, 1.5 * 300e3 / (0.1 * 4 * 60e3), 1e-12); // 18.75 m
  close(r.chuteLineLoad, 1e5, 1e-12);
  close(r.utilisationChute, 0.5e5 / cap, 1e-12);
  close(r.minimumChuteRadius, 300e3 / 2 / cap, 1e-12);
  assert.throws(() => umbilicalCrushCheck({ tubeOd: D, tubeWall: t, allowableStress: sa, nTracks: 5 }), /nTracks/);
});

test('core temperature: series thermal resistances by hand', () => {
  const r = umbilicalCoreTemperature({
    cores: [{ loss: 10, conductorDiameter: 0.006, insulationDiameter: 0.012, kInsulation: 0.3 }, { loss: 5, conductorDiameter: 0.006, insulationDiameter: 0.012, kInsulation: 0.3 }],
    ambientC: 10, bundleResistance: 0.5, sheath: { ri: 0.04, ro: 0.045, k: 0.4 }, hExternal: 200, outerDiameter: 0.09,
  });
  const rExt = 1 / (Math.PI * 0.09 * 200), rSh = Math.log(0.045 / 0.04) / (2 * Math.PI * 0.4), rIns = Math.log(2) / (2 * Math.PI * 0.3);
  close(r.surfaceTempC, 10 + 15 * rExt);
  close(r.coreTempsC[0], 10 + 15 * (rExt + rSh + 0.5) + 10 * rIns);
  close(r.coreTempsC[1], 10 + 15 * (rExt + rSh + 0.5) + 5 * rIns);
});

test('tube fatigue: stress range components and Miner damage by hand', () => {
  const s = demo();
  const unit = umbilicalLoadSharing(s, { tension: 1 }).elements.find((e) => e.name === 'SD tube');
  const curve = { a: 1e36, m: 3 };
  const bin = { tensionRange: 5e4, curvatureRange: 0.01, cycles: 1e6, meanTension: 2e5 };
  const r = umbilicalTubeFatigue(s, { element: 'SD tube', bins: [bin], curve, mu: 0.2, scf: 1.3 });
  const tPart = unit.force * 5e4 / tubeA, bPart = 2e11 * 0.00955 * 0.01;
  const f = 0.2 * (unit.contactInner + unit.contactOuter) * 2e5;
  const fPart = Math.min(2e11 * Math.cos(8 * DEG) ** 2 * 0.02 * 0.01, 2 * Math.PI * 0.02 * f / (2 * Math.sin(8 * DEG) * tubeA));
  close(r.ranges[0].tensionPart, tPart, 1e-12); close(r.ranges[0].bendingPart, bPart, 1e-12); close(r.ranges[0].frictionPart, fPart, 1e-12);
  const S = 1.3 * (tPart + bPart + fPart);
  close(r.damage, 1e6 / (1e36 * S ** -3), 1e-9);
  close(r.lifeWithDff, 1 / (10 * r.damage), 1e-9);
  assert.throws(() => umbilicalTubeFatigue(s, { element: 'nope', bins: [], curve }), /not found/);
});

test('Darcy pressure drop: laminar equals Hagen-Poiseuille; Haaland satisfies Colebrook within 2 %', () => {
  const lam = tubePressureDrop({ flow: 2e-6, d: 0.0127, length: 10000, rho: 1050, viscosity: 3e-3 });
  assert.equal(lam.regime, 'laminar');
  close(lam.dp, 128 * 3e-3 * 10000 * 2e-6 / (Math.PI * 0.0127 ** 4), 1e-12);
  const tur = tubePressureDrop({ flow: 2e-4, d: 0.0127, length: 100, rho: 1000, viscosity: 1e-3, roughness: 1.5e-6 });
  assert.equal(tur.regime, 'turbulent');
  const f = tur.frictionFactor;
  const colebrook = -2 * Math.log10(1.5e-6 / 0.0127 / 3.7 + 2.51 / (tur.reynolds * Math.sqrt(f)));
  close(1 / Math.sqrt(f), colebrook, 0.02);
  close(tur.dp, f * 100 / 0.0127 * 1000 * tur.velocity ** 2 / 2, 1e-12);
});

test('hydraulic response time: first-mode formula against a finite-difference solution of the line diffusion equation', () => {
  const o = { length: 20000, d: 0.0127, viscosity: 3e-3, bulkModulus: 2e9, volumetricExpansion: 1e-14 };
  const r = hydraulicResponseTime({ ...o, fraction: 0.9 });
  const rc = 128 * 3e-3 / (Math.PI * 0.0127 ** 4) * (Math.PI * 0.0127 ** 2 / 4 / 2e9 + 1e-14);
  close(r.tau1, 4 / Math.PI ** 2 * rc * 20000 ** 2, 1e-12);
  // explicit FD in normalised form p_t = p_xx on [0, 1], p(0) = 1, p_x(1) = 0
  const N = 60, dx = 1 / N, dt = 0.25 * dx * dx;
  let p = new Float64Array(N + 1); p[0] = 1;
  const tEnd = r.time / (rc * 20000 ** 2);
  const steps = Math.round(tEnd / dt);
  for (let k = 0; k < steps; k++) {
    const q = p.slice();
    for (let i = 1; i < N; i++) q[i] = p[i] + dt / (dx * dx) * (p[i + 1] - 2 * p[i] + p[i - 1]);
    q[N] = p[N] + dt / (dx * dx) * 2 * (p[N - 1] - p[N]);
    p = q;
  }
  close(p[N], 0.9, 5e-3);
});
