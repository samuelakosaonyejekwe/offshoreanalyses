import test from 'node:test';
import assert from 'node:assert/strict';
import {
  vonMises, lameStresses, thinWallStresses, thickWallLimitPressures, b313PressureDesign, dnvBurstPressure,
  tubeCollapsePressure, tubeCombinedUtilisation, proofTestPressure, tubeBendThinning, hoseRating,
  insulationResistance, looseTubeStrainWindow,
} from '../../src/engine/xsection/tubes.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);

test('Lame: boundary conditions, textbook thick cylinder ro/ri = 2', () => {
  const ri = 0.01, ro = 0.02, pInt = 60e6, pExt = 5e6;
  close(lameStresses({ pInt, pExt, ri, ro, r: ri }).radial, -pInt);
  close(lameStresses({ pInt, pExt, ri, ro, r: ro }).radial, -pExt);
  // ro/ri = 2, internal pressure only: hoop at bore = p (ro^2 + ri^2)/(ro^2 - ri^2) = 5p/3, at outer wall 2p/3
  close(lameStresses({ pInt: 30e6, ri, ro, r: ri }).hoop, 50e6);
  close(lameStresses({ pInt: 30e6, ri, ro, r: ro }).hoop, 20e6);
  close(lameStresses({ pInt: 30e6, ri, ro }).axial, 10e6); // p ri^2/(ro^2 - ri^2)
  // sigma_r + sigma_theta is constant through the wall (= 2A)
  const a = lameStresses({ pInt, pExt, ri, ro, r: 0.012 }), b = lameStresses({ pInt, pExt, ri, ro, r: 0.018 });
  close(a.radial + a.hoop, b.radial + b.hoop);
  // equilibrium: integral of hoop stress over the wall = pi ri - po ro
  let s = 0; const n = 20000;
  for (let i = 0; i < n; i++) s += lameStresses({ pInt, pExt, ri, ro, r: ri + (i + 0.5) * (ro - ri) / n }).hoop * (ro - ri) / n;
  close(s, pInt * ri - pExt * ro, 1e-8);
});

test('Lame converges to the thin-wall formula as D/t grows', () => {
  const p = 10e6, D = 0.1;
  let prev = Infinity;
  for (const t of [0.01, 0.005, 0.001, 0.0002]) {
    const ro = D / 2, ri = ro - t;
    // the mean-diameter formula matches the Lame bore stress to second order: ratio = 1 + (t/(D - t))^2
    const bore = lameStresses({ pInt: p, ri, ro, r: ri }).hoop;
    const thin = thinWallStresses({ pInt: p, D, t }).hoop;
    const err = Math.abs(bore / thin - 1);
    close(bore / thin, 1 + (t / (D - t)) ** 2, 1e-12);
    assert.ok(err < prev); prev = err;
  }
  assert.ok(prev < 1e-4);
  close(thinWallStresses({ pInt: p, D, t: 0.002, barlow: true }).hoop, p * D / 0.004);
  close(thinWallStresses({ pInt: p, D, t: 0.002 }).axial, p * (D - 0.002) / 0.008);
});

test('thick-wall limits: von Mises at the bore equals yield at pYield; plastic pressure (2/sqrt3) sy ln(ro/ri)', () => {
  const ri = 0.00635, ro = 0.00935, sy = 550e6;
  const l = thickWallLimitPressures({ ri, ro, yieldStress: sy, uts: 800e6 });
  close(lameStresses({ pInt: l.pYield, ri, ro }).vonMises, sy, 1e-12);
  close(l.pPlastic, 2 / Math.sqrt(3) * sy * Math.log(ro / ri));
  assert.ok(l.pPlastic > l.pYield && l.pBurstEstimate > l.pPlastic);
  close(vonMises({ s1: 100, s2: 0, s3: 0, tau: 0 }).vonMises, 100);
  close(vonMises({ s1: 0, s2: 0, s3: 0, tau: 100 }).vonMises, 100 * Math.sqrt(3));
});

test('ASME B31.3 Eq. (3a) round trip and hand value', () => {
  // D = 25.4 mm, S = 138 MPa, P = 20 MPa, Y = 0.4: t = 20*25.4/(2*(138 + 8)) = 1.7397 mm
  const r = b313PressureDesign({ D: 0.0254, S: 138e6, P: 20e6 });
  close(r.t, 20 * 25.4 / (2 * 146) / 1000, 1e-12);
  close(b313PressureDesign({ D: 0.0254, S: 138e6, t: r.t }).P, 20e6, 1e-12);
  assert.equal(r.Y, 0.4);
  assert.throws(() => b313PressureDesign({ D: 0.0254, S: 138e6 }), RangeError);
});

test('DNV burst form by hand', () => {
  const r = dnvBurstPressure({ D: 0.0191, t: 0.0015, fy: 550e6, fu: 800e6 });
  close(r.pb, 2 * 0.0015 / (0.0191 - 0.0015) * 550e6 * 2 / Math.sqrt(3));
  assert.equal(r.fcb, 550e6);
  assert.equal(dnvBurstPressure({ D: 0.0191, t: 0.0015, fy: 550e6, fu: 600e6 }).fcb, 600e6 / 1.15);
});

test('collapse cubic: limits (elastic, plastic) and the residual of the root', () => {
  const E = 2e11, fy = 550e6;
  // thin tube without ovality -> elastic collapse; thick -> plastic
  const thin = tubeCollapsePressure({ D: 0.1, t: 0.001, E, fy, ovality: 0, enforceMinimum: false });
  close(thin.pc, 2 * E * 1e-6 / 0.91, 1e-12);
  const thick = tubeCollapsePressure({ D: 0.1, t: 0.02, E, fy, ovality: 0, enforceMinimum: false });
  close(thick.pc, fy * 0.4, 1e-12);
  const r = tubeCollapsePressure({ D: 0.1, t: 0.004, E, fy, ovality: 0.01 });
  close((r.pc - r.pel) * (r.pc ** 2 - r.pp ** 2), r.pc * r.pel * r.pp * 0.01 * 25, 1e-7);
  assert.ok(r.pc < Math.min(r.pel, r.pp));
  // code minimum ovality enforced
  assert.equal(tubeCollapsePressure({ D: 0.1, t: 0.004, E, fy, ovality: 0.001 }).f0, 0.005);
});

test('combined utilisation: pure tension, pure bending and pressure cases by hand', () => {
  const D = 0.02, t = 0.002, A = Math.PI * (0.01 ** 2 - 0.008 ** 2);
  const ten = tubeCombinedUtilisation({ D, t, tension: 20e3, yieldStress: 550e6 });
  close(ten.vonMisesMax, 20e3 / A);
  const bend = tubeCombinedUtilisation({ D, t, curvature: 0.5, E: 2e11, yieldStress: 550e6, eta: 0.67 });
  close(bend.vonMisesMax, 2e11 * 0.01 * 0.5); close(bend.utilisation, 1e9 / (0.67 * 550e6));
  const I = Math.PI * (0.01 ** 4 - 0.008 ** 4) / 4;
  close(tubeCombinedUtilisation({ D, t, moment: 2e11 * I * 0.5, yieldStress: 550e6 }).vonMisesMax, 1e9, 1e-12);
  const pr = tubeCombinedUtilisation({ D, t, pInt: 50e6, yieldStress: 550e6 });
  close(pr.vonMisesMax, lameStresses({ pInt: 50e6, ri: 0.008, ro: 0.01 }).vonMises);
  assert.match(pr.location, /bore/);
});

test('proof pressure, bend thinning, hose rating, insulation resistance', () => {
  const p = proofTestPressure({ designPressure: 69e6, ri: 0.00635, ro: 0.00935, yieldStress: 550e6 });
  close(p.testPressure, 103.5e6);
  close(p.utilisation, lameStresses({ pInt: 103.5e6, ri: 0.00635, ro: 0.00935 }).vonMises / 550e6);
  const b = tubeBendThinning({ t: 0.002, D: 0.02, bendRadius: 0.06 });
  close(b.tExtrados, 0.002 * 6 / 7); close(b.tIntrados, 0.002 * 6 / 5); close(b.bendingStrain, 1 / 6);
  const h = hoseRating({ burstPressure: 138e6, volumetricExpansion: 2e-15, length: 10000, pressureChange: 20e6 });
  close(h.designPressure, 34.5e6); close(h.expansionVolume, 4e-4);
  const ir = insulationResistance({ resistivity: 1e14, conductorDiameter: 0.004, insulationDiameter: 0.008, length: 1000 });
  close(ir.resistance, 1e14 * Math.log(2) / (2 * Math.PI * 1000));
  close(ir.resistanceKm, ir.resistance);
});

test('loose-tube fibre strain window: central tube window equals the excess fibre length; stranded tube by exact geometry', () => {
  const c = looseTubeStrainWindow({ layLength: 0.1, clearance: 0.001, efl: 0.003 });
  close(c.tensileWindow, 0.003, 1e-12);
  // stranded: p = 100 mm, Rp = 5 mm, a = 0.5 mm, EFL = 0
  const p = 0.1, Rp = 0.005, a = 0.0005;
  const s = looseTubeStrainWindow({ layLength: p, pitchRadius: Rp, clearance: a, efl: 0, cableStrain: 0.012 });
  const Lf = Math.hypot(p, 2 * Math.PI * Rp);
  // after stretching by the window the fibre at radius Rp - a has exactly its free length
  close(Math.hypot(p * (1 + s.tensileWindow), 2 * Math.PI * (Rp - a)), Lf, 1e-12);
  close(Math.hypot(p * (1 - s.contractionWindow), 2 * Math.PI * (Rp + a)), Lf, 1e-12);
  // first-order check: window ~ (2 pi)^2 Rp a / p^2
  close(s.tensileWindow, 4 * Math.PI ** 2 * Rp * a / p ** 2, 0.06);
  close(s.fibreStrain, 0.012 - s.tensileWindow, 1e-12);
  close(s.utilisation, s.fibreStrain / 0.002, 1e-12);
  // more excess length widens the tensile window and narrows the contraction window
  const e = looseTubeStrainWindow({ layLength: p, pitchRadius: Rp, clearance: a, efl: 0.001 });
  assert.ok(e.tensileWindow > s.tensileWindow && e.contractionWindow < s.contractionWindow);
});

test('invalid arguments are named', () => {
  assert.throws(() => lameStresses({ pInt: 1, ri: 0.02, ro: 0.01 }), /ro/);
  assert.throws(() => tubeBendThinning({ t: 0.002, D: 0.02, bendRadius: 0.005 }), /bendRadius/);
  assert.throws(() => looseTubeStrainWindow({ layLength: -1, clearance: 0 }), /layLength/);
});
