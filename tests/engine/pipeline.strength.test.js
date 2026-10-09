import test from 'node:test';
import assert from 'node:assert/strict';
import * as S from '../../src/engine/pipeline/strength.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${a} vs ${b}`);

test('burst pressure: hand calculation, DNV-ST-F101 Eq. 5.8', () => {
  // D = 323.9 mm, t1 = 10.0 mm, X65 at ambient: fy = 450*0.96 = 432 MPa, fu = 535*0.96 = 513.6 MPa
  // fcb = min(432, 513.6/1.15 = 446.6) = 432 MPa
  // pb = 2*0.010/(0.3239 - 0.010) * 432 * 2/sqrt(3) = 0.0637145 * 498.831 = 31.783 MPa
  const r = S.burstPressure({ D: 0.3239, t: 0.010, fy: 432e6, fu: 513.6e6 });
  close(r.fcb, 432e6);
  close(r.pb, 31.783e6, 2e-5);
  // tensile-governed case: fy = 500, fu = 550 -> fcb = 478.26 MPa
  close(S.burstPressure({ D: 0.3239, t: 0.010, fy: 500e6, fu: 550e6 }).fcb, 550e6 / 1.15);
});

test('pressure containment: t1 definition, resistance factors and utilisation', () => {
  // t = 12.7 mm, 12.5 % fabrication tolerance, 3 mm corrosion allowance:
  // t1 = 12.7 - 1.5875 - 3 = 8.1125 mm (operation), 11.1125 mm (system test)
  // pb(t1) = 2*0.0081125/(0.3239-0.0081125)*432e6*1.154701 = 25.6297 MPa
  // medium safety class: gamma_m*gamma_SC = 1.15*1.138 = 1.3087 -> 19.584 MPa
  const r = S.pressureContainment({ D: 0.3239, t: 0.0127, tFabFrac: 0.125, tCorr: 0.003, fy: 432e6, fu: 513.6e6, pLocalInc: 16.5e6, pe: 1e6 });
  close(r.t1, 0.0081125, 1e-12);
  close(r.t1Test, 0.0111125, 1e-12);
  close(r.pb, 25.6297e6, 1e-5);
  close(r.pbDesign, r.pb / (1.15 * 1.138), 1e-12);
  close(r.uc, 15.5e6 / (r.pb / 1.3087), 1e-12);
  assert.equal(r.governing, 'burst');
  const f = S.dnvResistanceFactors('high');
  assert.deepEqual([f.gammaM, f.gammaScPc, f.gammaSc, f.gammaEps], [1.15, 1.308, 1.26, 3.3]);
  assert.equal(S.dnvResistanceFactors('low', 'FLS').gammaM, 1.0);
  const w = S.characteristicWallThickness({ t: 0.02, tFab: 0.001, tCorr: 0.003 });
  close(w.t1, 0.016); close(w.t2, 0.017);
});

test('pressure containment: system-test and mill-test terms', () => {
  // mill test: 2*tmin/(D - tmin)*min(0.96*450, 0.84*535 = 449.4) = 2*0.0111125/0.3127875*432e6 = 30.696 MPa
  const pm = S.millTestPressure({ D: 0.3239, tMin: 0.0111125, smys: 450e6, smts: 535e6 });
  close(pm, (2 * 0.0111125 / 0.3127875) * 432e6, 1e-12);
  close(pm, 30.696e6, 1e-4);
  const r = S.pressureContainment({ D: 0.3239, t: 0.0127, tFabFrac: 0.125, tCorr: 0.003, fy: 432e6, fu: 513.6e6, pLocalInc: 16.5e6, pe: 1e6,
    pLocalTest: 17.4e6, pMill: pm });
  close(r.ucSystemTest, 15.5e6 / (17.4e6 / 1.05 - 1e6), 1e-12);
  close(r.ucMillTest, 15.5e6 / (pm * 0.96 / 1.088), 1e-12);
  assert.equal(r.governing, 'system test');
  assert.ok(r.testCheck.uc > 0 && r.testCheck.capacity <= pm);
});

test('required wall thickness for burst inverts the check exactly', () => {
  const a = { D: 0.6096, fy: 432e6, fu: 513.6e6, pLocalInc: 20e6, pe: 2e6, safetyClass: 'high', tCorr: 0.003, tFabFrac: 0.1 };
  const r = S.requiredWallThicknessBurst(a);
  const chk = S.pressureContainment({ D: a.D, t: r.tNominal, tFabFrac: 0.1, tCorr: 0.003, fy: a.fy, fu: a.fu, pLocalInc: 20e6, pe: 2e6, safetyClass: 'high' });
  close(chk.uc, 1, 1e-10);
  close(chk.t1, r.t1, 1e-10);
});

test('collapse pressure: root satisfies the cubic; thin wall -> pel, thick wall -> pp', () => {
  for (const [D, t] of [[0.3239, 0.0127], [0.6096, 0.0191], [1.0, 0.01], [0.2, 0.03]]) {
    const r = S.collapsePressure({ D, t, fy: 450e6, alphaFab: 0.93, f0: 0.015 });
    const lhs = (r.pc - r.pel) * (r.pc ** 2 - r.pp ** 2);
    const rhs = r.pc * r.pel * r.pp * 0.015 * D / t;
    close(lhs, rhs, 1e-9);
    assert.ok(r.pc > 0 && r.pc < Math.min(r.pel, r.pp));
    // hand formulae for the two limit pressures
    close(r.pel, 2 * 207e9 * (t / D) ** 3 / (1 - 0.09), 1e-12);
    close(r.pp, 450e6 * 0.93 * 2 * t / D, 1e-12);
  }
  const thin = S.collapsePressure({ D: 1, t: 0.005, fy: 450e6 });       // D/t = 200: elastic
  assert.ok(thin.pc / thin.pel > 0.97 && thin.pc < thin.pel, `thin ${thin.pc / thin.pel}`);
  const thick = S.collapsePressure({ D: 1, t: 0.2, fy: 450e6 });        // D/t = 5: plastic
  assert.ok(thick.pc / thick.pp > 0.98 && thick.pc < thick.pp, `thick ${thick.pc / thick.pp}`);
  // ovality is floored at 0.5 % and more ovality lowers the capacity
  close(S.collapsePressure({ D: 0.3239, t: 0.0127, fy: 450e6, f0: 0.001 }).pc, S.collapsePressure({ D: 0.3239, t: 0.0127, fy: 450e6, f0: 0.005 }).pc);
  assert.ok(S.collapsePressure({ D: 0.3239, t: 0.0127, fy: 450e6, f0: 0.03 }).pc < S.collapsePressure({ D: 0.3239, t: 0.0127, fy: 450e6 }).pc);
});

test('collapse check and required thickness are consistent', () => {
  const a = { D: 0.6096, fy: 432e6, pe: 15e6, safetyClass: 'medium', alphaFab: 0.85, f0: 0.01 };
  const t1 = S.requiredWallThicknessCollapse(a);
  close(S.collapseCheck({ ...a, t1 }).uc, 1, 1e-8);
  const c = S.collapseCheck({ ...a, t1: 0.02 });
  close(c.capacity, c.pc / (1.15 * 1.14), 1e-12);
});

test('propagation pressure: hand calculation ppr = 35 fy alpha_fab (t/D)^2.5', () => {
  // t/D = 0.0127/0.3239 = 0.0392096; (t/D)^2.5 = 0.0392096^2 * sqrt(0.0392096) = 1.537393e-3 * 0.198014 = 3.04425e-4
  // ppr = 35 * 450e6 * 3.04425e-4 = 4.7947 MPa
  const r = S.propagationPressure({ D: 0.3239, t2: 0.0127, fy: 450e6, pe: 5e6 });
  close(r.ppr, 4.7947e6, 1e-4);
  close(r.capacity, r.ppr / (1.15 * 1.14), 1e-12);
  assert.equal(r.arrestorsRequired, true);
  assert.equal(S.propagationPressure({ D: 0.3239, t2: 0.0127, fy: 450e6, pe: 1e6 }).arrestorsRequired, false);
  // arrestor: crossover pressure tends to the arrestor propagation pressure for a long arrestor
  const ba = S.buckleArrestorCapacity({ D: 0.3239, t2: 0.0127, tBA: 0.03, LBA: 50, fy: 450e6 });
  close(ba.pX, ba.pprBA, 1e-9);
  assert.ok(S.buckleArrestorCapacity({ D: 0.3239, t2: 0.0127, tBA: 0.03, LBA: 0.3, fy: 450e6 }).pX < ba.pX);
});

test('plastic capacities and load-controlled local buckling terms (hand calculation)', () => {
  // D = 0.3239, t2 = 0.0127, fy = 432 MPa, fu = 513.6 MPa
  // Sp = 432e6*pi*0.3112*0.0127 = 5.3639e6 N ; Mp = 432e6*0.3112^2*0.0127 = 531.33e3 N.m
  const { Sp, Mp } = S.plasticCapacities({ D: 0.3239, t: 0.0127, fy: 432e6 });
  close(Sp, 5.3639e6, 1e-4);
  close(Mp, 531.33e3, 1e-4);
  // beta = (60 - 25.504)/90 = 0.38329 ; alpha_c = 1 - beta + beta*513.6/432 = 1.07240
  const r = S.localBucklingLoadControlled({ D: 0.3239, t2: 0.0127, fy: 432e6, fu: 513.6e6, MSd: 150e3, SSd: -600e3, pi: 12e6, pe: 2e6 });
  close(r.beta, (60 - 0.3239 / 0.0127) / 90, 1e-12);
  close(r.alphaC, 1.0724, 1e-4);
  const g = 1.15 * 1.14;
  const m = g * 150e3 / (r.alphaC * Mp);
  const s = (g * 600e3 / (r.alphaC * Sp)) ** 2;
  const q = 10e6 / r.pb;                                   // < 2/3 -> alpha_p = 1 - beta
  assert.ok(q < 2 / 3);
  close(r.alphaP, 1 - r.beta, 1e-12);
  close(r.uc, (m + s) ** 2 + ((1 - r.beta) * q / r.alphaC) ** 2, 1e-12);
  assert.equal(r.mode, 'internal');
  // the allowable moment brings the interaction exactly to unity
  const lim = S.localBucklingLoadControlled({ D: 0.3239, t2: 0.0127, fy: 432e6, fu: 513.6e6, MSd: r.allowableMoment, SSd: -600e3, pi: 12e6, pe: 2e6 });
  close(lim.uc, 1, 1e-10);
  // high pressure branch of alpha_p
  const hp = S.localBucklingLoadControlled({ D: 0.3239, t2: 0.0127, fy: 432e6, fu: 513.6e6, MSd: 0, SSd: 0, pi: 0.8 * r.pb, pe: 0 });
  close(hp.alphaP, 1 - 3 * r.beta * (1 - 0.8), 1e-12);
});

test('load-controlled local buckling, external overpressure: reduces to the collapse check without M and S', () => {
  const a = { D: 0.3239, t2: 0.0127, fy: 432e6, fu: 513.6e6, pi: 0, pe: 8e6, alphaFab: 0.93 };
  const r = S.localBucklingLoadControlled({ ...a, MSd: 0, SSd: 0 });
  const pc = S.collapsePressure({ D: a.D, t: a.t2, fy: a.fy, alphaFab: 0.93 }).pc;
  assert.equal(r.mode, 'external');
  close(r.uc, (1.15 * 1.14 * 8e6 / pc) ** 2, 1e-12);
  assert.ok(S.localBucklingLoadControlled({ ...a, MSd: 200e3, SSd: 0 }).uc > r.uc);
  assert.ok(S.localBucklingLoadControlled({ D: 1, t2: 0.01, fy: 432e6, fu: 513.6e6, MSd: 0, SSd: 0 }).warnings.length > 0);
});

test('displacement-controlled local buckling: strain capacity hand calculation', () => {
  // t/D = 0.0392096; eps_c(0) = 0.78*(0.0392096 - 0.01)*0.93^-1.5 = 0.78*0.0292096*1.11500 = 0.025403
  const r = S.localBucklingDisplacementControlled({ D: 0.3239, t2: 0.0127, fy: 432e6, fu: 513.6e6, epsSd: 0.004 });
  close(r.epsC, 0.025403, 1e-4);
  close(r.epsRd, r.epsC / 2.5, 1e-12);
  close(r.uc, 0.004 / r.epsRd, 1e-12);
  const p = S.localBucklingDisplacementControlled({ D: 0.3239, t2: 0.0127, fy: 432e6, fu: 513.6e6, epsSd: 0.004, pMin: 10e6 });
  const pb = S.burstPressure({ D: 0.3239, t: 0.0127, fy: 432e6, fu: 513.6e6 }).pb;
  close(p.epsC, r.epsC * (1 + 5.75 * 10e6 / pb), 1e-12);
  const e = S.localBucklingDisplacementControlled({ D: 0.3239, t2: 0.0127, fy: 432e6, fu: 513.6e6, epsSd: 0.004, pe: 5e6 });
  assert.equal(e.mode, 'external');
  close(e.uc, (0.004 / r.epsRd) ** 0.8 + 5e6 / (e.pc / (1.15 * 1.14)), 1e-12);
});

test('ASME B31.8 / B31.4 stresses: hand calculation, restrained line', () => {
  // SH = 10e6*0.3239/(2*0.0127) = 127.52 MPa ; allowable hoop 0.72*450 = 324 MPa -> 0.3936
  // SL = 0.3*127.52 - 207e9*11.7e-6*50 = 38.256 - 121.095 = -82.839 MPa
  // Tresca: 127.52 + 82.839 = 210.36 MPa ; allowable 0.9*450 = 405 MPa
  const r = S.asmeStressCheck({ D: 0.3239, t: 0.0127, smys: 450e6, pi: 10e6, deltaT: 50 });
  close(r.SH, 127.52e6, 1e-4);
  close(r.ucHoop, 0.3936, 1e-3);
  close(r.SLmax, -82.839e6, 1e-4);
  close(r.Scombined, 210.36e6, 1e-4);
  close(r.ucCombined, 210.36 / 405, 1e-4);
  const vm = S.asmeStressCheck({ D: 0.3239, t: 0.0127, smys: 450e6, pi: 10e6, deltaT: 50, criterion: 'vonMises' });
  close(vm.Scombined, Math.sqrt(r.SH ** 2 - r.SH * r.SLmax + r.SLmax ** 2), 1e-12);
  assert.ok(vm.Scombined < r.Scombined);
  const free = S.asmeStressCheck({ D: 0.3239, t: 0.0127, smys: 450e6, pi: 10e6, restrained: false });
  close(free.SLmax, r.SH / 2, 1e-12);
  close(r.tRequiredHoop, 10e6 * 0.3239 / (2 * 324e6), 1e-12);
  assert.equal(S.b318DesignFactor('3'), 0.5);
  assert.equal(S.b318DesignFactor('1-1'), 0.8);
  close(S.b318TemperatureFactor(100), 1);
  close(S.b318TemperatureFactor((300 - 32) / 1.8), 0.967, 1e-9);
  close(S.b318TemperatureFactor((425 - 32) / 1.8), 0.8835, 1e-9);
});

test('Lame thick-wall stresses: boundary conditions and closed forms', () => {
  const D = 0.4; const t = 0.05; const pi = 30e6; const pe = 5e6;
  const ro = 0.2; const ri = 0.15;
  const inner = S.lameStresses({ D, t, pi, pe });
  const outer = S.lameStresses({ D, t, pi, pe, r: ro });
  close(inner.sigmaR, -pi, 1e-12);
  close(outer.sigmaR, -pe, 1e-12);
  // internal pressure only: hoop at bore = p (ro^2 + ri^2)/(ro^2 - ri^2) = 30*(0.04+0.0225)/(0.0175) = 107.14 MPa
  close(S.lameStresses({ D, t, pi }).sigmaH, 30e6 * 0.0625 / 0.0175, 1e-12);
  // sigma_h + sigma_r is constant through the wall (= 2A)
  close(inner.sigmaH + inner.sigmaR, outer.sigmaH + outer.sigmaR, 1e-12);
  close(inner.sigmaAxialCapped, (pi * ri * ri - pe * ro * ro) / (ro * ro - ri * ri), 1e-12);
  // thin-wall limit: mean hoop stress ~ p Dmean / 2t
  const thin = S.lameStresses({ D: 1, t: 0.005, pi: 1e6, r: 0.4975 });
  close(thin.sigmaH, 1e6 * 0.995 / (2 * 0.005), 6e-3);
});

test('equivalent stress, hydrotest, bend radius', () => {
  const r = S.equivalentStressCheck({ D: 0.3239, t: 0.0127, fy: 450e6, eta: 0.9, pi: 10e6 });
  close(r.sigmaH, 10e6 * 0.3112 / 0.0254, 1e-12);
  close(r.sigmaE, r.sigmaH, 1e-12);                       // no axial stress: von Mises = hoop
  const h = S.hydrotestPressure({ pDesign: 10e6, D: 0.3239, t: 0.0127, smys: 450e6 });
  close(h.pTest, 12.5e6);
  close(h.hoopFractionSmys, 12.5e6 * 0.3239 / 0.0254 / 450e6, 1e-12);
  const b = S.bendingStrainRadius({ D: 0.3239, R: 100, epsAllow: 0.002, sigmaAllow: 0.72 * 450e6 });
  close(b.strain, 0.0016195, 1e-12);
  close(b.minRadiusStrain, 80.975, 1e-12);
  close(b.minRadiusElastic, 207e9 * 0.3239 / (2 * 324e6), 1e-12);
});

test('Brazier ovalisation: limit moment 0.987 E r t^2/sqrt(1-nu^2) at zeta = 2/9', () => {
  const D = 0.5; const t = 0.01; const r = 0.245;
  const b0 = S.brazierOvalisation({ D, t, curvature: 0 });
  const b = S.brazierOvalisation({ D, t, curvature: b0.criticalCurvature });
  close(b.zeta, 2 / 9, 1e-12);
  close(b.moment, b.brazierMoment, 1e-12);
  close(b.brazierMoment, 0.98738 * 207e9 * r * t * t / Math.sqrt(0.91), 1e-4);
  // the moment-curvature curve peaks there
  assert.ok(S.brazierOvalisation({ D, t, curvature: 0.9 * b0.criticalCurvature }).moment < b.moment);
  assert.ok(S.brazierOvalisation({ D, t, curvature: 1.1 * b0.criticalCurvature }).moment < b.moment);
  const o = S.ovalisationBendingPressure({ D, t, strain: 0.002, f0: 0.005, pe: 1e6, pc: 4e6 });
  close(o.ovality, (0.005 + 0.03 * (1 + 50 / 120) * (2 * 0.002 * 50) ** 2) / 0.75, 1e-12);
});

test('B31 elbow factors h = tR/r^2, k = 1.65/h, i = 0.9/h^(2/3); bend thinning', () => {
  // t = 12.7 mm, R = 1.5 m, r = 0.1556 m: h = 0.01905/0.0242114 = 0.78682 ; k = 2.0970 ; i = 1.0560
  const e = S.elbowFactors({ t: 0.0127, R: 1.5, r: 0.1556 });
  close(e.h, 0.78682, 1e-4);
  close(e.k, 2.0970, 1e-4);
  close(e.i, 1.0560, 1e-4);
  assert.equal(S.elbowFactors({ t: 0.03, R: 3, r: 0.1 }).i, 1);       // floor at 1
  const b = S.bendWallThickness({ D: 0.3239, t: 0.0127, R: 5 * 0.3239 });
  close(b.tExtrados, 0.0127 * 5 / 5.5, 1e-12);
  close(b.lorenzIntrados, 19 / 18, 1e-12);
  close(b.lorenzExtrados, 21 / 22, 1e-12);
});
