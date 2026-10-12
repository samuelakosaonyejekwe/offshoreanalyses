import test from 'node:test';
import assert from 'node:assert/strict';
import {
  skinEffect, proximityEffect, acResistance, capacitance, dielectricLoss, sheathLossSingleCore, sheathEddyLossFactor,
  sheathLossThreeCoreSL, armourLossThreeCoreSteelWire, layerThermalResistance, geometricFactorSL, t4Buried, t4Group,
  t4TrefoilTouching, t4FlatTouching, mutualHeatingRise, backfillCorrection, ductAirResistance,
  freeAirDissipationCoefficient, freeAirSurfaceRise, jTubeSurfaceTemperature, ratingEquation, cableGeometry,
  externalThermalModel, temperaturesAtCurrent, rateCable, sweepRating, routeHotspot,
} from '../../src/engine/electrical/ampacity.js';

const rel = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${a} vs ${b} (tol ${tol})`);

// 132 kV single-core cable used in several tests.
const single = {
  cores: 1,
  conductor: { material: 'Cu', area: 1000e-6, diameter: 0.0376 },
  conductorScreenThickness: 0.0012,
  insulation: { material: 'xlpe', thickness: 0.016 },
  insulationScreenThickness: 0.0012,
  sheath: { material: 'lead', thickness: 0.003 },
  serving: { thickness: 0.004, rhoT: 3.5 },
};
// 132 kV three-core SL-type export cable.
const threeCore = {
  cores: 3,
  conductor: { material: 'Cu', area: 800e-6, diameter: 0.0337 },
  conductorScreenThickness: 0.0015,
  insulation: { material: 'xlpe', thickness: 0.015 },
  insulationScreenThickness: 0.0015,
  sheath: { material: 'lead', thickness: 0.0025 },
  coreJacket: { thickness: 0.0025, rhoT: 3.5 },
  bedding: { thickness: 0.004, rhoT: 6 },
  armour: { material: 'steel', wireDiameter: 0.006 },
  serving: { thickness: 0.004, rhoT: 6 },
};
const ac = { f: 50, U: 132e3, thetaMaxC: 90 };

test('skin effect by hand, 1000 mm² Cu, 50 Hz, 90 °C: xs² = 8π·50·1e-7/2.24418e-5 = 5.5996, ys = 0.14444', () => {
  const Rdc = 0.0176e-3 * (1 + 0.00393 * 70); // 2.244176e-5 Ω/m
  const { xs, ys } = skinEffect({ Rdc, f: 50, ks: 1 });
  rel(xs * xs, 5.59955, 1e-5);
  // xs⁴ = 31.35496; ys = 31.35496/(192 + 0.8·31.35496) = 31.35496/217.08397
  rel(ys, 0.144437, 1e-5);
});

test('skin effect: Milliken coefficient ks = 0.435 scales xs², and the three ranges join continuously', () => {
  const Rdc = 2.244176e-5;
  rel(skinEffect({ Rdc, f: 50, ks: 0.435 }).xs ** 2, 0.435 * 5.59955, 1e-5);
  // continuity of the piecewise definition at xs = 2.8 and 3.8 (IEC 60287-1-1 A1:2014) within 0.5 %
  const at = (xs) => skinEffect({ Rdc: (8 * Math.PI * 50 * 1e-7) / (xs * xs), f: 50 }).ys;
  rel(at(2.8 + 1e-9), at(2.8 - 1e-9), 5e-3);
  rel(at(3.8 + 1e-9), at(3.8 - 1e-9), 5e-3);
});

test('proximity effect by hand, three cables dc = 37.6 mm at s = 90 mm: yp = F·(dc/s)²·[0.312·(dc/s)² + 1.18/(F + 0.27)] = 0.07315', () => {
  const Rdc = 2.244176e-5;
  // F = 0.144437 (same x as the skin effect with kp = 1); (dc/s)² = 0.174538
  // bracket = 0.312·0.174538 + 1.18/0.414437 = 0.054456 + 2.847236 = 2.901692
  const { yp } = proximityEffect({ Rdc, f: 50, kp: 1, dc: 0.0376, s: 0.09, conductors: 3 });
  rel(yp, 0.144437 * 0.174538 * 2.901692, 1e-5);
  rel(yp, 0.07315, 1e-3);
  // two conductors: yp = 2.9·F·(dc/s)²
  rel(proximityEffect({ Rdc, f: 50, dc: 0.0376, s: 0.09, conductors: 2 }).yp, 2.9 * 0.144437 * 0.174538, 1e-5);
});

test('AC resistance R = R′(1 + ys + yp); DC returns R′', () => {
  const r = acResistance({ R20: 0.0176e-3, thetaC: 90, alpha20: 0.00393, f: 50, dc: 0.0376, s: 0.09 });
  rel(r.Rac, 2.244176e-5 * (1 + 0.144437 + 0.07315), 1e-4);
  const d = acResistance({ R20: 0.0176e-3, thetaC: 90, alpha20: 0.00393, f: 0 });
  assert.equal(d.Rac, d.Rdc); assert.equal(d.ys, 0);
});

test('capacitance and dielectric loss by hand, 132 kV XLPE: C = 2.5/(18·ln(72/40))·1e-9 = 0.23629 nF/m, Wd = 0.4311 W/m', () => {
  const C = capacitance({ epsR: 2.5, Di: 0.072, dc: 0.040 });
  rel(C, 0.236292e-9, 1e-5);
  // Wd = 2π·50·C·(132e3/√3)²·0.001 = 314.159·0.236292e-9·5.808e9·0.001
  rel(dielectricLoss({ f: 50, C, U0: 132e3 / Math.sqrt(3), tanDelta: 0.001 }), 0.43115, 1e-4);
});

test('sheath circulating loss, trefoil, by hand: X = 2ω·1e-7·ln(2.5) = 5.7572e-5, λ1′ = 12/(1 + 27.153) = 0.42624', () => {
  const out = sheathLossSingleCore({ R: 2.5e-5, Rs: 3e-4, f: 50, d: 0.08, s: 0.1, formation: 'trefoil', bonding: 'solid' });
  rel(out.X, 5.75722e-5, 1e-5);
  rel(out.lambda1, 0.42624, 1e-4);
  assert.equal(out.lambda1Eddy, 0);
});

test('sheath loss, flat formation: transposed uses ln(2·∛2·s/d); untransposed three-cable factors reduce to it on average', () => {
  const base = { R: 2.5e-5, Rs: 3e-4, f: 50, d: 0.08, s: 0.2, formation: 'flat', bonding: 'solid' };
  const tr = sheathLossSingleCore({ ...base, transposed: true });
  const X1 = 2 * 2 * Math.PI * 50 * 1e-7 * Math.log((2 * Math.cbrt(2) * 0.2) / 0.08);
  rel(tr.lambda1, (3e-4 / 2.5e-5) / (1 + (3e-4 / X1) ** 2), 1e-12);
  const un = sheathLossSingleCore(base);
  // Leading + lagging outer cables: the ± cross term cancels, mean = Rs/R·(0.75P²/(Rs²+P²) + 0.25Q²/(Rs²+Q²))
  const X = 2 * 2 * Math.PI * 50 * 1e-7 * Math.log(0.4 / 0.08); const Xm = 2 * 2 * Math.PI * 50 * 1e-7 * Math.LN2;
  const P = X + Xm; const Q = X - Xm / 3; const Rs = 3e-4;
  rel(0.5 * (un.perCable.leading + un.perCable.lagging), 12 * ((0.75 * P * P) / (Rs * Rs + P * P) + (0.25 * Q * Q) / (Rs * Rs + Q * Q)), 1e-12);
  rel(un.perCable.centre, (12 * Q * Q) / (Rs * Rs + Q * Q), 1e-12);
  assert.equal(un.lambda1, Math.max(un.perCable.leading, un.perCable.centre, un.perCable.lagging));
  // The three untransposed cables lose, in total, about what three transposed cables lose (within 5 %).
  rel(un.perCable.leading + un.perCable.centre + un.perCable.lagging, 3 * tr.lambda1, 0.05);
});

test('single-point bonding: no circulating loss; eddy factor by hand for a lead sheath in trefoil', () => {
  // Rs = 3e-4 Ω/m → m = ω/Rs·1e-7 = 0.10472; d/2s = 0.4; λ0 = 3·m²/(1+m²)·0.16 = 0.0052066
  // Δ1 = (1.14·m^2.45 + 0.33)·0.4^(0.92m + 1.66) ; β1 = sqrt(4πω/(1e7·ρs)), ρs = 21.4e-8·1.2 = 2.568e-7
  const Rs = 3e-4; const R = 2.5e-5; const w = 2 * Math.PI * 50;
  const m = (w / Rs) * 1e-7;
  rel(m, 0.104720, 1e-5);
  const lambda0 = 3 * (m * m / (1 + m * m)) * 0.16;
  const delta1 = (1.14 * m ** 2.45 + 0.33) * 0.4 ** (0.92 * m + 1.66);
  const beta1 = Math.sqrt((4 * Math.PI * w) / (1e7 * 2.568e-7));
  const gs = 1 + (0.003 / 0.083) ** 1.74 * (beta1 * 0.083 - 1.6);
  const expected = (Rs / R) * (gs * lambda0 * (1 + delta1) + (beta1 * 0.003) ** 4 / 12);
  const e = sheathEddyLossFactor({ R, Rs, f: 50, d: 0.08, s: 0.1, thickness: 0.003, outerDiameter: 0.083, rhoS: 2.568e-7, position: 'trefoil' });
  rel(e.lambda1Eddy, expected, 1e-12);
  assert.ok(e.lambda1Eddy > 0.05 && e.lambda1Eddy < 0.09, `eddy factor ${e.lambda1Eddy}`);
  assert.ok(Math.abs(e.gs - 1) < 0.01);
  const sp = sheathLossSingleCore({ R, Rs, f: 50, d: 0.08, s: 0.1, bonding: 'singlePoint',
    sheath: { thickness: 0.003, outerDiameter: 0.083, rhoS: 2.568e-7 } });
  assert.equal(sp.lambda1Circulating, 0);
  rel(sp.lambda1, e.lambda1Eddy, 1e-12);
});

test('three-core SL sheath loss carries the factor 1.5 when armoured', () => {
  const a = sheathLossThreeCoreSL({ R: 3e-5, Rs: 4e-4, f: 50, d: 0.072, s: 0.08 });
  const b = sheathLossThreeCoreSL({ R: 3e-5, Rs: 4e-4, f: 50, d: 0.072, s: 0.08, armoured: false });
  rel(a.lambda1 / b.lambda1, 1.5, 1e-12);
  rel(b.lambda1, sheathLossSingleCore({ R: 3e-5, Rs: 4e-4, f: 50, d: 0.072, s: 0.08 }).lambda1, 1e-12);
});

test('steel-wire armour loss by hand: λ2 = 1.23·(6e-5/3e-5)·(0.09/0.18)²/((2.77·6e-5·1e6/314.159)² + 1) = 0.48052', () => {
  const out = armourLossThreeCoreSteelWire({ R: 3e-5, RA: 6e-5, f: 50, c: 0.045, dA: 0.18 });
  rel(out.lambda2, 0.480516, 1e-5);
  // SL screening factor (1 − R/Rs·λ1′) and a user reduction factor
  const scr = armourLossThreeCoreSteelWire({ R: 3e-5, RA: 6e-5, f: 50, c: 0.045, dA: 0.18, Rs: 4e-4, lambda1Trefoil: 0.2, reductionFactor: 0.5 });
  rel(scr.lambda2, 0.480516 * (1 - (3e-5 / 4e-4) * 0.2) * 0.5, 1e-5);
  assert.throws(() => armourLossThreeCoreSteelWire({ R: 3e-5, RA: 6e-5, f: 50, c: 0.045, dA: 0.18, reductionFactor: 1.2 }), RangeError);
});

test('T1 by hand: 3.5/(2π)·ln(1 + 2·18/40) = 0.55704·0.64185 = 0.35754 K·m/W', () => {
  rel(layerThermalResistance({ rhoT: 3.5, innerDiameter: 0.040, thickness: 0.018 }), 0.357541, 1e-5);
});

test('SL-type geometric factor: Ḡ(0.1) = 2π·0.0974053 = 0.61202; the two fitted branches agree within 1.5 % at X = 0.03', () => {
  rel(geometricFactorSL(0.1), 0.612016, 1e-5);
  rel(geometricFactorSL(0.03), geometricFactorSL(0.03 + 1e-12), 0.015);
  assert.equal(geometricFactorSL(0), 0);
  assert.throws(() => geometricFactorSL(0.2), RangeError);
});

test('T4 by hand (Kennelly): ρ = 1, L = 1 m, De = 0.1 m → u = 20, ln(20 + √399)/(2π) = 0.58700 K·m/W', () => {
  rel(t4Buried({ rhoT: 1, depth: 1, De: 0.1 }), 0.587004, 1e-5);
  assert.throws(() => t4Buried({ rhoT: 1, depth: 0.04, De: 0.1 }), RangeError);
});

test('T4 touching trefoil, metallic sheath: 1.5/π·(ln 40 − 0.630) = 1.46051 K·m/W; touching flat: 0.475·ln 40 − 0.346', () => {
  rel(t4TrefoilTouching({ rhoT: 1, depth: 1, De: 0.1 }), 1.460507, 1e-5);
  rel(t4FlatTouching({ rhoT: 1, depth: 1, De: 0.1 }), 0.475 * Math.log(40) - 0.346, 1e-12);
});

test('group T4 by images, three cables flat at s: centre cable = ρ/(2π)·[ln(u + √(u²−1)) + ln(1 + (2L/s)²)] (IEC 60287-2-1 §4.2.3.3.1)', () => {
  const positions = [{ x: -0.3, depth: 1 }, { x: 0, depth: 1 }, { x: 0.3, depth: 1 }];
  const t = t4Group({ rhoT: 1.2, positions, index: 1, De: 0.1 });
  rel(t, (1.2 / (2 * Math.PI)) * (Math.log(20 + Math.sqrt(399)) + Math.log(1 + (2 / 0.3) ** 2)), 1e-12);
  // outer cable is cooler than the centre cable
  assert.ok(t4Group({ rhoT: 1.2, positions, index: 0, De: 0.1 }) < t);
  // the touching-trefoil formula and the image method differ by a few per cent only
  const s = 0.1; const h = s / Math.sqrt(3);
  const tre = [{ x: 0, depth: 1 - h }, { x: -s / 2, depth: 1 + h / 2 }, { x: s / 2, depth: 1 + h / 2 }];
  const worst = Math.max(...tre.map((_, i) => t4Group({ rhoT: 1, positions: tre, index: i, De: 0.1 })));
  rel(worst, t4TrefoilTouching({ rhoT: 1, depth: 1, De: 0.1 }), 0.1);
});

test('mutual heating by images: one source 50 W/m at 0.5 m spacing, both 1 m deep → ρ/(2π)·50·ln(√(0.25+4)/0.5)', () => {
  const r = mutualHeatingRise({ rhoT: 1, target: { x: 0, depth: 1 }, sources: [{ x: 0.5, depth: 1, W: 50 }] });
  rel(r, (50 / (2 * Math.PI)) * Math.log(Math.sqrt(4.25) / 0.5), 1e-12);
});

test('backfill equivalent radius by hand (Neher–McGrath): 0.5 m × 1.0 m → ln rb = 0.25·(4/π − 0.5)·ln 5 + ln 0.25, rb = 0.34124 m', () => {
  const b = backfillCorrection({ width: 1.0, height: 0.5, depthCentre: 1.2, rhoSoil: 1.5, rhoBackfill: 0.8, nCables: 3 });
  rel(b.rb, 0.341239, 1e-5);
  const u = 1.2 / 0.341239;
  rel(b.deltaT4, (3 / (2 * Math.PI)) * 0.7 * Math.log(u + Math.sqrt(u * u - 1)), 1e-5);
  assert.throws(() => backfillCorrection({ width: 3.5, height: 1, depthCentre: 2 }), RangeError);
});

test('duct air-space resistance: plastic duct, De = 80 mm, θm = 50 °C → 1.87/(1 + 0.1·(0.312 + 0.0037·50)·80) = 1.87/4.976 = 0.37580', () => {
  rel(ductAirResistance({ De: 0.08, thetaMeanC: 50, kind: 'plastic' }), 1.87 / (1 + 0.1 * 0.497 * 80), 1e-12);
  rel(ductAirResistance({ De: 0.08, thetaMeanC: 50, kind: 'plastic' }), 0.375804, 1e-5);
});

test('free air: h = Z/De^g + E for a single cable (0.21, 3.94, 0.60) and the surface energy balance closes', () => {
  const h = freeAirDissipationCoefficient({ De: 0.1 });
  rel(h, 0.21 / 0.1 ** 0.6 + 3.94, 1e-12);
  const rise = freeAirSurfaceRise({ W: 60, De: 0.1, h, absorptivity: 0.4, solarIntensity: 1000 });
  rel(Math.PI * 0.1 * h * rise ** 1.25, 60 + 0.4 * 0.1 * 1000, 1e-12);
});

test('J-tube energy balance: heat leaving the cable equals heat crossing the gap and (with solar gain) leaving the tube', () => {
  const p = { W: 80, De: 0.15, tubeInnerDiameter: 0.3, tubeOuterDiameter: 0.32, ambientC: 25, hGap: 3, hOut: 10,
    emissivityCable: 0.9, emissivityTubeInner: 0.8, emissivityTubeOuter: 0.8, absorptivity: 0.6, solarIntensity: 800 };
  const r = jTubeSurfaceTemperature(p);
  const sb = 5.670374419e-8; const Ts = r.surfaceC + 273.15; const Tt = r.tubeC + 273.15; const Ta = 298.15;
  const fr = 1 / (1 / 0.9 + (0.15 / 0.3) * (1 / 0.8 - 1));
  rel(Math.PI * 0.15 * (3 * (Ts - Tt) + fr * sb * (Ts ** 4 - Tt ** 4)), 80, 1e-8);
  rel(Math.PI * 0.32 * (10 * (Tt - Ta) + 0.8 * sb * (Tt ** 4 - Ta ** 4)), 80 + 0.6 * 0.32 * 800, 1e-8);
  // solar radiation makes the cable hotter
  assert.ok(r.surfaceC > jTubeSurfaceTemperature({ ...p, solarIntensity: 0 }).surfaceC);
});

test('rating equation, DC closed form: I = sqrt(Δθ/(R·ΣT))', () => {
  const I = ratingEquation({ dTheta: 55, R: 1.2e-5, T1: 0.4, T3: 0.05, T4: 0.7 });
  rel(I, Math.sqrt(55 / (1.2e-5 * 1.15)), 1e-12);
});

test('DC cable rating reproduces I = sqrt(Δθ/(R(θmax)·(T1 + T3 + T4))) with R = R20·(1 + α·(θmax − 20))', () => {
  const g = cableGeometry(single);
  const inst = { type: 'buried', ambientC: 15, depth: 1.2, rhoSoil: 0.9 };
  const r = rateCable({ cable: single, installation: inst, system: { dc: true, thetaMaxC: 70 } });
  const R = 0.0176e-3 * (1 + 0.00393 * 50);
  const T4 = t4Buried({ rhoT: 0.9, depth: 1.2, De: g.De });
  rel(r.I, Math.sqrt(55 / (R * (g.T1 + g.T3 + T4))), 1e-7);
  assert.equal(r.lambda1, 0); assert.equal(r.Wd, 0); assert.equal(r.ys, 0);
  rel(r.conductorC, 70, 1e-7);
});

test('DC cable: insulation temperature-drop limit governs when I²·R·T1 would exceed it', () => {
  const inst = { type: 'buried', ambientC: 15, depth: 1.2, rhoSoil: 0.9 };
  const g = cableGeometry(single);
  const r = rateCable({ cable: single, installation: inst, system: { dc: true, thetaMaxC: 70, dTInsulationMaxK: 8 } });
  assert.equal(r.governedBy, 'insulationTemperatureDrop');
  rel(r.I * r.I * r.R * g.T1, 8, 1e-6);
  assert.ok(r.I < r.IThermal && r.conductorC < 70);
});

test('cable geometry: three touching cores → laid-up diameter (1 + 2/√3)·Dcore, c = Dcore/√3; T1, T2, T3 from the layers', () => {
  const g = cableGeometry(threeCore);
  const core = 0.0337 + 2 * (0.0015 + 0.015 + 0.0015 + 0.0025 + 0.0025);
  rel(g.coreDiameter, core, 1e-12);
  rel(g.laidUpDiameter, core * (1 + 2 / Math.sqrt(3)), 1e-12);
  rel(g.c, core / Math.sqrt(3), 1e-12);
  rel(g.T1, (3.5 / (2 * Math.PI)) * Math.log((0.0337 + 0.036) / 0.0337), 1e-12);
  const Ds = core - 0.005;
  const Tj = (3.5 / (2 * Math.PI)) * Math.log(core / Ds);
  rel(g.T2, (6 / (6 * Math.PI)) * geometricFactorSL(0.004 / core) + Tj / 3, 1e-12);
  rel(g.T3, (6 / (2 * Math.PI)) * Math.log(g.De / g.armourOuterDiameter), 1e-12);
  rel(g.De, g.laidUpDiameter + 2 * (0.004 + 0.006 + 0.004), 1e-12);
  rel(g.Rs20, 21.4e-8 / (Math.PI * g.sheathMeanDiameter * 0.0025), 1e-12);
});

test('energy balance at the converged AC rating: Σ losses × thermal resistances = Δθ, and the IEC rating equation returns the same current', () => {
  const inst = { type: 'buried', ambientC: 12, depth: 1.5, rhoSoil: 0.7 };
  const r = rateCable({ cable: threeCore, installation: inst, system: ac });
  const g = cableGeometry(threeCore);
  const T4 = t4Buried({ rhoT: 0.7, depth: 1.5, De: g.De });
  const n = 3; const Wc = r.I * r.I * r.R;
  const dTheta = (Wc + 0.5 * r.Wd) * g.T1 + (Wc * (1 + r.lambda1) + r.Wd) * n * g.T2
    + (Wc * (1 + r.lambda1 + r.lambda2) + r.Wd) * n * (g.T3 + T4);
  rel(dTheta, 90 - 12, 1e-7);
  const I = ratingEquation({ dTheta: 78, R: r.R, T1: g.T1, T2: g.T2, T3: g.T3, T4, n, Wd: r.Wd, lambda1: r.lambda1, lambda2: r.lambda2 });
  rel(I, r.I, 1e-7);
  rel(r.losses.total, n * (Wc * (1 + r.lambda1 + r.lambda2) + r.Wd), 1e-12);
  assert.ok(r.lambda1 > 0 && r.lambda2 > 0 && r.sheathC < 90 && r.armourC < r.sheathC && r.surfaceC < r.armourC);
  // plausibility: an 800 mm² Cu 132 kV three-core cable buried at 1.5 m carries several hundred amperes
  assert.ok(r.I > 600 && r.I < 1100, `rating ${r.I}`);
});

test('inverse problem: temperatures at the rated current return θmax; a lower current gives a cooler conductor', () => {
  const inst = { type: 'buried', ambientC: 12, depth: 1.5, rhoSoil: 0.7 };
  const r = rateCable({ cable: threeCore, installation: inst, system: ac });
  rel(temperaturesAtCurrent({ cable: threeCore, installation: inst, system: ac, I: r.I }).conductorC, 90, 1e-7);
  const half = temperaturesAtCurrent({ cable: threeCore, installation: inst, system: ac, I: 0.5 * r.I });
  assert.ok(half.conductorC < 40 && half.conductorC > 12);
  // no load: only dielectric loss heats the cable, θc − θa = Wd·(0.5·T1 + n·(T2 + T3 + T4))
  const g = cableGeometry(threeCore);
  const zero = temperaturesAtCurrent({ cable: threeCore, installation: inst, system: ac, I: 0 });
  rel(zero.conductorC - 12, zero.Wd * (0.5 * g.T1 + 3 * (g.T2 + g.T3 + t4Buried({ rhoT: 0.7, depth: 1.5, De: g.De }))), 1e-9);
});

test('armour-loss reduction factor (CIGRE TB 880 guidance) raises the rating; bonding and formation change it as expected', () => {
  const inst = { type: 'buried', ambientC: 12, depth: 1.5, rhoSoil: 0.7 };
  const base = rateCable({ cable: threeCore, installation: inst, system: ac });
  const red = rateCable({ cable: threeCore, installation: inst, system: { ...ac, armourLossReduction: 0.5 } });
  rel(red.lambda2 / base.lambda2, 0.5, 0.02);
  assert.ok(red.I > base.I);
  const b = { type: 'buried', ambientC: 15, depth: 1.0, rhoSoil: 1.0 };
  const solid = rateCable({ cable: single, installation: { ...b, formation: 'trefoil', bonding: 'solid' }, system: ac });
  const cross = rateCable({ cable: single, installation: { ...b, formation: 'trefoil', bonding: 'cross' }, system: ac });
  const flat = rateCable({ cable: single, installation: { ...b, formation: 'flat', spacing: 0.3, bonding: 'cross' }, system: ac });
  assert.ok(cross.I > solid.I && flat.I > cross.I);
  assert.equal(cross.lambda1Circulating, 0);
});

test('partial drying-out of the soil (IEC 60287-1-1 §1.4.2): rating equals the two-zone equation with v = ρd/ρw', () => {
  const g = cableGeometry(single);
  const inst = { type: 'buried', ambientC: 20, depth: 1.0, rhoSoil: 1.0, dryZone: { rhoDry: 2.5, criticalRiseK: 30 } };
  const r = rateCable({ cable: single, installation: inst, system: { dc: true, thetaMaxC: 90 } });
  assert.equal(r.dried, true);
  const R = 0.0176e-3 * (1 + 0.00393 * 70);
  const T4 = t4Buried({ rhoT: 1.0, depth: 1.0, De: g.De });
  rel(r.I, ratingEquation({ dTheta: 70, R, T1: g.T1, T3: g.T3, T4, v: 2.5, dThetaX: 30 }), 1e-7);
  const wet = rateCable({ cable: single, installation: { ...inst, dryZone: undefined }, system: { dc: true, thetaMaxC: 90 } });
  assert.ok(r.I < wet.I);
  // a critical rise that is never reached leaves the rating unchanged
  const never = rateCable({ cable: single, installation: { ...inst, dryZone: { rhoDry: 2.5, criticalRiseK: 80 } }, system: { dc: true, thetaMaxC: 90 } });
  rel(never.I, wet.I, 1e-9); assert.equal(never.dried, false);
});

test('free-air rating with solar radiation: surface energy balance closes and the IEC 60287-1-1 §1.4.4 equation returns the same current', () => {
  const g = cableGeometry(single);
  const inst = { type: 'air', ambientC: 30, absorptivity: 0.4, solarIntensity: 1000 };
  const r = rateCable({ cable: single, installation: inst, system: { dc: true, thetaMaxC: 90 } });
  const h = freeAirDissipationCoefficient({ De: g.De });
  const W = r.losses.total;
  const dThetaS = r.surfaceC - 30;
  rel(Math.PI * g.De * h * dThetaS ** 1.25, W + 0.4 * g.De * 1000, 1e-7);
  // IEC 60287-1-1 §1.4.4: I² = (Δθ − σ·De·H·T4*)/(R·(T1 + T3 + T4*)), T4* = 1/(π·De·h·Δθs^0.25)
  const T4s = 1 / (Math.PI * g.De * h * dThetaS ** 0.25);
  rel(r.I, Math.sqrt((60 - 0.4 * g.De * 1000 * T4s) / (r.R * (g.T1 + g.T3 + T4s))), 1e-7);
  const shade = rateCable({ cable: single, installation: { ...inst, solarIntensity: 0 }, system: { dc: true, thetaMaxC: 90 } });
  assert.ok(shade.I > r.I);
});

test('seabed-laid, backfilled, ducted and externally heated installations', () => {
  const g = cableGeometry(single);
  const dc = { dc: true, thetaMaxC: 90 };
  const sea = externalThermalModel({ type: 'seabed', ambientC: 10, h: 200 }, g);
  rel(sea.T4, 1 / (Math.PI * g.De * 200), 1e-12);
  const buried = rateCable({ cable: single, installation: { type: 'buried', ambientC: 15, depth: 1.2, rhoSoil: 1.5 }, system: dc });
  const backfilled = rateCable({ cable: single, installation: { type: 'buried', ambientC: 15, depth: 1.2, rhoSoil: 1.5,
    backfill: { width: 0.8, height: 0.6, depthCentre: 1.2, rhoBackfill: 0.8 } }, system: dc });
  const allGood = rateCable({ cable: single, installation: { type: 'buried', ambientC: 15, depth: 1.2, rhoSoil: 0.8 }, system: dc });
  assert.ok(backfilled.I > buried.I && backfilled.I < allGood.I);
  // HDD duct filled with bentonite: T4 = ρfill/(2π)·ln(Di/De) + ρduct/(2π)·ln(Do/Di) + soil part on the duct diameter
  const duct = { innerDiameter: 0.2, outerDiameter: 0.225, rhoDuct: 3.5, fill: 'solid', rhoFill: 1.2 };
  const d = rateCable({ cable: single, installation: { type: 'duct', ambientC: 15, depth: 8, rhoSoil: 1.0, duct }, system: dc });
  const T4 = (1.2 / (2 * Math.PI)) * Math.log(0.2 / g.De) + (3.5 / (2 * Math.PI)) * Math.log(0.225 / 0.2) + t4Buried({ rhoT: 1, depth: 8, De: 0.225 });
  rel(d.T4eff, T4, 1e-9);
  rel(d.I, Math.sqrt(75 / (d.R * (g.T1 + g.T3 + T4))), 1e-7);
  const airDuct = rateCable({ cable: single, installation: { type: 'duct', ambientC: 15, depth: 8, rhoSoil: 1.0, duct: { ...duct, fill: 'air' } }, system: dc });
  assert.ok(airDuct.I > 0 && Number.isFinite(airDuct.T4eff));
  // an external heat source (crossing cable) reduces the permissible rise by its image-method contribution
  const ext = rateCable({ cable: single, installation: { type: 'buried', ambientC: 15, depth: 1.2, rhoSoil: 1.5,
    externalSources: [{ x: 0.4, depth: 0.8, W: 40 }] }, system: dc });
  rel(ext.externalRise, mutualHeatingRise({ rhoT: 1.5, target: { x: 0, depth: 1.2 }, sources: [{ x: 0.4, depth: 0.8, W: 40 }] }), 1e-12);
  rel(ext.I / buried.I, Math.sqrt((75 - ext.externalRise) / 75), 1e-6);
});

test('sweeps: rating falls monotonically with burial depth and with soil thermal resistivity', () => {
  const inst = { type: 'buried', ambientC: 12, depth: 1.0, rhoSoil: 0.7 };
  const d = sweepRating({ cable: threeCore, installation: inst, system: ac, parameter: 'depth', values: [0.5, 1, 2, 4, 8] });
  const r = sweepRating({ cable: threeCore, installation: inst, system: ac, parameter: 'rhoSoil', values: [0.5, 0.7, 1.0, 1.5] });
  for (let i = 1; i < d.I.length; i++) assert.ok(d.I[i] < d.I[i - 1]);
  for (let i = 1; i < r.I.length; i++) assert.ok(r.I[i] < r.I[i - 1]);
});

test('route hotspot: the section with the worst thermal environment limits the route', () => {
  const sections = [
    { name: 'Offshore buried', length: 40e3, installation: { type: 'buried', ambientC: 10, depth: 1.5, rhoSoil: 0.7 } },
    { name: 'Landfall HDD', length: 800, installation: { type: 'duct', ambientC: 15, depth: 12, rhoSoil: 1.2,
      duct: { innerDiameter: 0.35, outerDiameter: 0.4, rhoDuct: 3.5, fill: 'solid', rhoFill: 1.0 } } },
    { name: 'J-tube', length: 30, installation: { type: 'jtube', ambientC: 25, tubeInnerDiameter: 0.45, tubeOuterDiameter: 0.48, solarIntensity: 600 } },
  ];
  const out = routeHotspot({ cable: threeCore, system: ac, sections, I: 600 });
  assert.equal(out.sections.length, 3);
  assert.equal(out.I, Math.min(...out.sections.map((s) => s.I)));
  assert.equal(out.limitingName, 'Landfall HDD');
  assert.equal(out.sections[2].kpStart, 40800);
  assert.ok(out.sections.every((s) => Number.isFinite(s.conductorAtLoadC)));
});

test('invalid input is rejected with RangeError naming the argument', () => {
  assert.throws(() => skinEffect({ Rdc: -1, f: 50 }), /Rdc/);
  assert.throws(() => capacitance({ epsR: 2.5, Di: 0.03, dc: 0.04 }), /Di/);
  assert.throws(() => rateCable({ cable: single, installation: { type: 'buried', ambientC: 95, depth: 1, rhoSoil: 1 }, system: ac }), /thetaMaxC/);
  assert.throws(() => rateCable({ cable: single, installation: { type: 'buried', depth: 1, rhoSoil: 1 }, system: ac }), /ambientC/);
  assert.throws(() => cableGeometry({ cores: 2, conductor: { material: 'Cu', area: 1e-4 }, insulation: { material: 'xlpe', thickness: 0.01 } }), /cores/);
});

test('three-core bedding thicker than 0.15 × core diameter: the error names the bedding thickness and the limit; T2 override still works', () => {
  // core diameter = 33.7 + 2·(1.5 + 15 + 1.5 + 2.5 + 2.5) = 79.7 mm → limit 0.15 × 79.7 = 11.96 mm
  const thick = { ...threeCore, bedding: { thickness: 0.015, rhoT: 6 } };
  assert.throws(() => cableGeometry(thick), (e) => e instanceof RangeError && /cable\.bedding\.thickness/.test(e.message)
    && /0\.15/.test(e.message) && /12\.0 mm/.test(e.message) && /T2/.test(e.message));
  assert.throws(() => geometricFactorSL(0.2), /0–0\.15.*IEC 60287-2-1/);
  // just inside the limit, and with T2 supplied directly, the geometry is built
  assert.ok(cableGeometry({ ...threeCore, bedding: { thickness: 0.0119, rhoT: 6 } }).T2 > 0);
  rel(cableGeometry({ ...thick, T2: 0.05 }).T2, 0.05, 1e-12);
});
