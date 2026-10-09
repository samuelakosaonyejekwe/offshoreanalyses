import test from 'node:test';
import assert from 'node:assert/strict';
import * as C from '../../src/engine/pipeline/crossing.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${a} vs ${b}`);

const D = 0.3239; const t = 0.0127; const E = 207e9;
const EI = E * Math.PI / 64 * (D ** 4 - (D - 2 * t) ** 4);

test('beam solver: simply supported span, M = w L^2/8, deflection 5 w L^4/(384 EI)', () => {
  const L = 30; const w = 1500;
  const r = C.beamOnSupports({ L, EI, w, ends: 'pinned', nElem: 60, D });
  close(r.maxMoment, w * L * L / 8, 1e-6);
  close(r.xMaxMoment, 15, 1e-12);
  close(-r.deflection[30], 5 * w * L ** 4 / (384 * EI), 1e-5);
  close(r.totalReaction, w * L, 1e-6);
  close(r.maxStress, (w * L * L / 8) * D / (2 * EI / E), 1e-6);
  close(Math.abs(r.shear[0]), w * L / 2, 1e-5);
  assert.equal(r.spans.length, 1); close(r.spans[0].length, L, 1e-12);
});

test('beam solver: two-span continuous beam, central reaction 1.25 w l and hogging moment w l^2/8', () => {
  const l = 20; const w = 1000;
  const r = C.beamOnSupports({ L: 2 * l, EI, w, ends: 'pinned', supports: [{ x: l, unilateral: false }], nElem: 80 });
  close(r.reactions[0].R, 1.25 * w * l, 1e-5);
  close(r.moment[40], -w * l * l / 8, 1e-5);
  close(r.totalReaction, 2 * w * l, 1e-6);
  // point load on a fixed-fixed beam: M = F L/8 at the ends and at mid-span, deflection F L^3/(192 EI)
  const p = C.beamOnSupports({ L: 20, EI, ends: 'fixed', pointLoads: [{ x: 10, F: 5e4 }], nElem: 40 });
  close(p.moment[20], 5e4 * 20 / 8, 1e-5);
  close(-p.deflection[20], 5e4 * 20 ** 3 / (192 * EI), 1e-5);
});

test('beam solver: pipe over a prop on a stiff seabed — touchdown length (72 EI d/w)^(1/4), reaction 4 w l/3, moment w l^2/6', () => {
  const w = 1000; const del = 0.5;
  const l = (72 * EI * del / w) ** 0.25;
  const r = C.beamOnSupports({ L: 8 * l, EI, w, supports: [{ x: 4 * l, height: del }], seabed: { stiffness: 1e9 }, nElem: 400 });
  close(r.reactions[0].R, 4 * w * l / 3, 5e-3);
  close(r.maxMoment, w * l * l / 6, 5e-3);
  close(r.xMaxMoment, 4 * l, 1e-12);
  assert.equal(r.spans.length, 2);
  close(r.spans[0].length, l, 0.03);                 // resolved to one element (l/50)
  close(r.totalReaction, r.totalLoad, 1e-4);
  close(r.deflection[200], del, 1e-4);
  // settlement of the support lowers the reaction
  const s = C.beamOnSupports({ L: 8 * l, EI, w, supports: [{ x: 4 * l, height: del - 0.1 }], seabed: { stiffness: 1e9 }, nElem: 400 });
  assert.ok(s.reactions[0].R < r.reactions[0].R);
  close(s.reactions[0].R, 4 * w * (72 * EI * 0.4 / w) ** 0.25 / 3, 5e-3);
});

test('beam solver: unilateral support below the deflected shape carries nothing; mechanism is reported', () => {
  const r = C.beamOnSupports({ L: 30, EI, w: 1500, ends: 'pinned', supports: [{ x: 15, height: -1 }], nElem: 60 });
  assert.equal(r.reactions[0].active, false);
  assert.equal(r.reactions[0].R, 0);
  // a support just under the unsupported deflection picks up load
  const sag = 5 * 1500 * 30 ** 4 / (384 * EI);
  const q = C.beamOnSupports({ L: 30, EI, w: 1500, ends: 'pinned', supports: [{ x: 15, height: -0.5 * sag }], nElem: 60 });
  // superposition: R L^3/(48 EI) = sag/2  ->  R = 24 EI sag/L^3 = (5/16) w L
  close(q.reactions[0].R, (5 / 16) * 1500 * 30, 1e-5);
  assert.throws(() => C.beamOnSupports({ L: 30, EI, w: 1500, nElem: 20 }), Error);
});

test('crossing separation and quantity take-off', () => {
  const s = C.crossingSeparation({ axisElevation: 1.2, Dcross: 0.41, topOfExisting: 0.6, settlement: 0.05 });
  close(s.separation, 1.2 - 0.205 - 0.05 - 0.6, 1e-12);
  assert.equal(s.ok, true);
  assert.equal(C.crossingSeparation({ axisElevation: 1.0, Dcross: 0.41, topOfExisting: 0.6 }).ok, false);
  const q = C.crossingQuantities({ mattresses: [{ count: 6, length: 6, width: 3, thickness: 0.3 }],
    rockBerms: [{ D: 0.4, H: 0.6, crestWidth: 1, sideSlope: 3, length: 50, rhoSolid: 2650 }] });
  close(q.mattressVolume, 6 * 5.4, 1e-12);
  close(q.mattressMass, 6 * 5.4 * 2400, 1e-12);
  close(q.mattressSubmergedWeight, 6 * 5.4 * (2400 - 1025) * 9.80665, 1e-12);
  close(q.rockVolume, 50 * (4 - Math.PI * 0.04), 1e-12);
});

test('dropped object: terminal velocity and impact energy by hand', () => {
  // m = 2000 kg, V = 0.5 m3, CD = 1.0, A = 1.0 m2: vT = sqrt(2*9.80665*(2000 - 512.5)/(1025*1*1)) = sqrt(28.4632) = 5.3351 m/s
  const r = C.droppedObjectImpact({ mass: 2000, volume: 0.5, CD: 1.0, area: 1.0, Ca: 1.0 });
  close(r.terminalVelocity, 5.3351, 1e-4);
  close(r.kineticEnergy, 0.5 * 2000 * 28.4632, 1e-4);
  close(r.impactEnergy, 0.5 * (2000 + 512.5) * 28.4632, 1e-4);
  // at terminal velocity drag balances the submerged weight
  close(0.5 * 1025 * 1.0 * 1.0 * r.terminalVelocity ** 2, (2000 - 512.5) * 9.80665, 1e-12);
  assert.throws(() => C.droppedObjectImpact({ mass: 100, volume: 0.5, CD: 1, area: 1 }), RangeError);
  close(C.impactEnergy({ mass: 3000, addedMass: 1500, velocity: 2.5 }), 0.5 * 4500 * 6.25, 1e-12);
});

test('pipe dent energy and protection: hand calculation and inversion', () => {
  // mp = 0.25*450e6*0.0127^2 = 18 145 N ; E = 16*sqrt(2 pi/9)*mp*sqrt(D/t)*D*(0.05)^1.5
  //   = 16*0.835548*18145*5.05016*0.3239*0.0111803 = 4 436 J
  const r = C.steelPipeDentEnergy({ D, t, fy: 450e6, dentRatio: 0.05 });
  close(r.mp, 18145, 1e-4);
  close(r.energy, 16 * Math.sqrt(2 * Math.PI / 9) * 18145.125 * Math.sqrt(D / t) * D * 0.05 ** 1.5, 1e-12);
  close(r.energy, 4436, 2e-3);
  close(C.steelPipeDentEnergy({ D, t, fy: 450e6, energy: r.energy }).dentRatio, 0.05, 1e-12);
  close(C.concreteCoatingEnergy({ Y: 120e6, b: 0.03, h: 0.2, x0: 0.04 }), 120e6 * 0.03 * 0.2 * 0.04, 1e-12);
  close(C.concreteCoatingEnergy({ mode: 'flat', Y: 120e6, b: 0.03, D: 0.41, x0: 0.04 }), 120e6 * 0.03 * (4 / 3) * Math.sqrt(0.41 * 0.04 ** 3), 1e-12);
  const p = C.impactProtection({ impactEnergy: 40e3, absorbed: [28.8e3], D, t, fy: 450e6 });
  close(p.energyToSteel, 11.2e3, 1e-12);
  assert.equal(p.ok, false);
  close(p.additionalProtection, 11.2e3 - r.energy, 1e-9);
  assert.equal(C.impactProtection({ impactEnergy: 30e3, absorbed: [28.8e3], D, t, fy: 450e6 }).ok, true);
});

test('trawl pull-over: force format and quasi-static beam response', () => {
  // Fp = CF V sqrt(m k) = 4*2.5*sqrt(3000*1e5) = 10*17320.5 = 173 205 N
  const r = C.trawlPullover({ CF: 4, V: 2.5, mass: 3000, kWarp: 1e5, L: 40, EI, D });
  close(r.Fp, 173205, 1e-5);
  close(r.momentMax, r.Fp * 40 / 8, 1e-12);
  close(r.deflection, r.Fp * 40 ** 3 / (192 * EI), 1e-12);
  // agrees with the beam solver for a fixed-fixed span with a central point load
  const b = C.beamOnSupports({ L: 40, EI, ends: 'fixed', pointLoads: [{ x: 20, F: r.Fp }], nElem: 40 });
  close(b.maxMoment, r.momentMax, 1e-5);
  const s = C.trawlPullover({ Fp: 1e5, L: 40, EI, q: 1000, ends: 'pinned' });
  close(s.momentMax, 1e5 * 10 - 1000 * 200, 1e-12);
  assert.equal(C.trawlPullover({ Fp: 1e4, L: 40, EI, q: 1000 }).soilHoldsLoad, true);
});

test('rock armour size from the Shields criterion', () => {
  // tau = 0.5*1025*0.02*2^2 = 41 Pa ; D50 = 41/(0.035*(2650 - 1025)*9.80665) = 0.07351 m
  const r = C.rockSizeShields({ thetaCr: 0.035, U: 2, f: 0.02 });
  close(r.tau, 41, 1e-12);
  close(r.D50, 0.07351, 1e-3);
  close(C.rockSizeShields({ thetaCr: 0.035, tau: 41 }).D50, r.D50, 1e-12);
});
