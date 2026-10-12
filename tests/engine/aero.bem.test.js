import test from 'node:test';
import assert from 'node:assert/strict';
import {
  BETZ_LIMIT, actuatorDisc, thrustCoefficientBuhl, inductionFromCt, glauertIdealCp, prandtlLoss, towerPotentialFactor,
  bemElement, bemSolve, cpSurface, surfaceLookup, maxCp, powerCurve,
} from '../../src/engine/aero/bem.js';
import { nrel5mw, designOptimalRotor } from '../../src/engine/aero/rotors.js';
import { createPolar } from '../../src/engine/aero/polar.js';

const D2R = Math.PI / 180;
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} +- ${tol}`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${tol})`);

const rotor = nrel5mw();
const flat = { tiltDeg: 0, preconeDeg: 0 };

test('actuator disc: Betz limit 16/27 at a = 1/3 with CT = 8/9', () => {
  close(BETZ_LIMIT, 0.592593, 1e-6);
  const d = actuatorDisc(1 / 3);
  close(d.cp, 16 / 27, 1e-14); close(d.ct, 8 / 9, 1e-14); close(d.wakeSpeedRatio, 1 / 3, 1e-14);
  // maximum of 4a(1-a)^2 is at a = 1/3
  assert.ok(actuatorDisc(0.32).cp < d.cp && actuatorDisc(0.35).cp < d.cp);
  assert.throws(() => actuatorDisc(0.6), RangeError);
});

test('Buhl (2005) high-thrust closure: value and slope continuous at a = 0.4, CT = 2 at a = 1, inverse round trip', () => {
  for (const F of [1, 0.8, 0.5]) {
    close(thrustCoefficientBuhl(0.4, F), 0.96 * F, 1e-12);
    const h = 1e-6;
    const sl = (thrustCoefficientBuhl(0.4 + h, F) - thrustCoefficientBuhl(0.4, F)) / h;
    close(sl, 4 * F * (1 - 0.8), 1e-4, 'slope');
    close(thrustCoefficientBuhl(1, F), 2, 1e-12);
    for (const a of [0.1, 0.3, 0.45, 0.7, 0.95]) close(inductionFromCt(thrustCoefficientBuhl(a, F), F), a, 1e-9);
  }
});

test('glauertIdealCp: Manwell et al., Wind Energy Explained, Table 3.2 and the Betz limit at high tip-speed ratio', () => {
  const tab = [[0.5, 0.289], [1.0, 0.416], [1.5, 0.477], [2.0, 0.511], [2.5, 0.533], [5.0, 0.570], [7.5, 0.581], [10.0, 0.585]];
  for (const [l, cp] of tab) close(glauertIdealCp(l), cp, 1.5e-3, `lambda ${l}`);   // table quoted to three decimals
  close(glauertIdealCp(200), 16 / 27, 1e-4);
});

test('prandtlLoss: hand value, F -> 1 for many blades / mid-span, F -> 0 at the tip and hub', () => {
  // B = 3, r/R = 0.9, phi = 10 deg: f = 1.5*0.1/(0.9*sin 10) = 0.9598, F = (2/pi) acos(exp(-0.9598)) = 0.7498
  close(prandtlLoss({ blades: 3, r: 9, tipRadius: 10, phi: 10 * D2R }).F, 0.7498, 2e-4);
  close(prandtlLoss({ blades: 3, r: 5, tipRadius: 10, phi: 5 * D2R }).F, 1, 1e-6);
  close(prandtlLoss({ blades: 300, r: 9.5, tipRadius: 10, phi: 10 * D2R }).F, 1, 1e-6);
  close(prandtlLoss({ blades: 3, r: 10, tipRadius: 10, phi: 10 * D2R }).F, 0, 1e-12);
  close(prandtlLoss({ blades: 3, r: 1, tipRadius: 10, hubRadius: 1, phi: 30 * D2R }).Fhub, 0, 1e-12);
  const both = prandtlLoss({ blades: 3, r: 2, tipRadius: 10, hubRadius: 1, phi: 30 * D2R });
  close(both.F, both.Ftip * both.Fhub, 1e-15);
});

test('towerPotentialFactor: 1 - (a/x)^2 on the stagnation line, speed-up abeam, zero inside', () => {
  close(towerPotentialFactor(-5, 0, 2), 1 - 0.16, 1e-14);
  close(towerPotentialFactor(0, 4, 2), 1.25, 1e-14);
  close(towerPotentialFactor(-500, 3, 2), 1, 1e-4);
  assert.equal(towerPotentialFactor(1, 0, 2), 0);
});

test('BEM on a Glauert-optimum rotor without drag or tip loss approaches the ideal Cp with wake rotation', () => {
  // ideal section: Cl = 2 pi alpha everywhere, no drag
  const polar = createPolar({ alphaDeg: [-180, 180], cl: [-2 * Math.PI * Math.PI, 2 * Math.PI * Math.PI], cd: [0, 0] });
  for (const tsr of [4, 7, 10]) {
    const ro = designOptimalRotor({ tipRadius: 50, hubRadius: 0.5, blades: 200, tsr, clDesign: 2 * Math.PI * 5 * D2R, alphaDesignDeg: 5, polar, nStations: 80 });
    const s = bemSolve({ rotor: ro, U: 10, omega: tsr * 10 / 50, tipLoss: false, hubLoss: false });
    rel(s.Cp, glauertIdealCp(tsr), 0.015, `Cp at lambda ${tsr}`);
    assert.ok(s.Cp < BETZ_LIMIT);
    // outboard annuli run at the Betz induction a = 1/3 and CT -> 8/9
    close(s.spanwise[60].a, 1 / 3, 5e-3);
    close(s.spanwise[60].alpha, 5 * D2R, 2e-3);
  }
  const ro = designOptimalRotor({ tipRadius: 50, hubRadius: 0.5, blades: 200, tsr: 12, clDesign: 2 * Math.PI * 5 * D2R, alphaDesignDeg: 5, polar, nStations: 120 });
  const s = bemSolve({ rotor: ro, U: 10, omega: 12 * 10 / 50, tipLoss: false, hubLoss: false });
  close(s.Cp, 16 / 27, 0.012, 'Cp tends to the Betz limit at high tip-speed ratio');
  close(s.Ct, 8 / 9, 0.01);
});

test('momentum / blade-element consistency: annulus thrust and torque balance at every station', () => {
  const U = 9; const om = 1.0; const rho = 1.225;
  const s = bemSolve({ rotor, U, omega: om, pitchDeg: 1, rho, ...flat });
  assert.equal(s.converged, true);
  let checked = 0;
  for (const e of s.spanwise) {
    if (e.a > 0.4) continue;
    checked++;
    // dT/dr = 4 pi r rho U^2 a (1 - a) F = B Np ;  dQ/dr = 4 pi r^3 rho U Omega a' (1 - a) F = B Tp r
    rel(3 * e.Np, 4 * Math.PI * e.r * rho * U * U * e.a * (1 - e.a) * e.F, 1e-7, `thrust r=${e.r}`);
    rel(3 * e.Tp * e.r, 4 * Math.PI * e.r ** 3 * rho * U * om * e.ap * (1 - e.a) * e.F, 1e-7, `torque r=${e.r}`);
    // velocity triangle
    close(Math.tan(e.phi), U * (1 - e.a) / (om * e.r * (1 + e.ap)), 1e-8);
  }
  assert.ok(checked >= 14);
});

test('heavily loaded annuli follow the Buhl thrust relation (a > 0.4)', () => {
  const U = 6; const om = 1.25; const rho = 1.225;
  const s = bemSolve({ rotor, U, omega: om, pitchDeg: 0, rho, ...flat });
  const heavy = s.spanwise.filter((e) => e.a > 0.42);
  assert.ok(heavy.length >= 3, 'case must contain heavily loaded annuli');
  for (const e of heavy) {
    const ctLocal = 3 * e.Np / (0.5 * rho * U * U * 2 * Math.PI * e.r);
    rel(ctLocal, thrustCoefficientBuhl(e.a, e.F), 1e-6);
  }
});

test('bemElement: residual of the inflow-angle equation is zero and the solution is reproducible', () => {
  const st = rotor.stations[10];
  const o = { Vx: 10, Vy: 1.2 * st.r, r: st.r, chord: st.chord, theta: st.twistDeg * D2R, polar: st.polar, blades: 3, tipRadius: 63, hubRadius: 1.5 };
  const e = bemElement(o);
  assert.equal(e.converged, true);
  close(Math.sin(e.phi) / (1 - e.a) - (o.Vx / o.Vy) * Math.cos(e.phi) / (1 + e.ap), 0, 1e-8);
  const out = {};
  assert.equal(bemElement(o, out), out);
  assert.equal(out.phi, e.phi);
});

test('NREL 5-MW planform with generic polars (plausibility band, not a validation): peak Cp 0.44-0.50 at tip-speed ratio 7-8.5', () => {
  const m = maxCp({ rotor, bem: flat });
  assert.ok(m.cpMax > 0.44 && m.cpMax < 0.50, `cpMax ${m.cpMax}`);
  assert.ok(m.tsrOpt > 7 && m.tsrOpt < 8.5, `tsrOpt ${m.tsrOpt}`);
  assert.ok(Math.abs(m.pitchOptDeg) < 2, `pitch ${m.pitchOptDeg}`);
  assert.ok(m.ct > 0.6 && m.ct < 0.95);
  // it is a maximum
  const R = rotor.tipRadius;
  for (const [dl, dp] of [[0.3, 0], [-0.3, 0], [0, 0.5], [0, -0.5]]) {
    assert.ok(bemSolve({ rotor, U: 10, omega: (m.tsrOpt + dl) * 10 / R, pitchDeg: m.pitchOptDeg + dp, ...flat }).Cp <= m.cpMax + 1e-6);
  }
});

test('NREL 5-MW power curve with generic polars (plausibility): rated power reached between 11 and 12 m/s and held to cut-out', () => {
  const pc = powerCurve({ rotor, speeds: [2, 3, 6, 9, 11, 12, 13, 18, 25, 26] });
  assert.ok(pc.ratedWind > 11 && pc.ratedWind < 12, `rated wind ${pc.ratedWind}`);
  assert.equal(pc.power[0], 0); assert.equal(pc.region[0], 'off');
  assert.equal(pc.power[9], 0); assert.equal(pc.region[9], 'off');
  assert.ok(pc.power[4] < 5e6);
  for (const i of [5, 6, 7, 8]) { rel(pc.power[i], 5e6, 1e-4); assert.equal(pc.region[i], '3'); }
  // monotonic below rated, pitch increasing above rated, thrust falling above rated
  assert.ok(pc.power[1] < pc.power[2] && pc.power[2] < pc.power[3] && pc.power[3] < pc.power[4]);
  assert.ok(pc.pitchDeg[5] < pc.pitchDeg[6] && pc.pitchDeg[6] < pc.pitchDeg[7] && pc.pitchDeg[7] < pc.pitchDeg[8]);
  assert.ok(pc.thrust[8] < pc.thrust[5]);
  assert.ok(pc.pitchDeg[8] > 18 && pc.pitchDeg[8] < 28, `pitch at 25 m/s ${pc.pitchDeg[8]}`);
  // electrical = efficiency x aerodynamic
  rel(pc.power[3], 0.944 * pc.aeroPower[3], 1e-12);
  // rotor speed limits
  rel(pc.omega[1], rotor.minRotorSpeed, 1e-12); rel(pc.omega[8], rotor.ratedRotorSpeed, 1e-12);
});

test('same rotor in water: loads scale with density, coefficients unchanged (tidal use)', () => {
  const a = bemSolve({ rotor, U: 2.5, omega: 0.3, rho: 1.225, ...flat });
  const w = bemSolve({ rotor, U: 2.5, omega: 0.3, rho: 1025, mu: 1.4e-3, ...flat });
  rel(w.thrust / a.thrust, 1025 / 1.225, 1e-12);
  rel(w.Cp, a.Cp, 1e-12);
  assert.ok(w.spanwise[8].Re > a.spanwise[8].Re);
  rel(w.spanwise[8].Re, 1025 * w.spanwise[8].W * rotor.stations[8].chord / 1.4e-3, 1e-12);
});

test('coefficients are consistent: Cp = Cq * tsr, power = torque * omega', () => {
  const s = bemSolve({ rotor, U: 8, omega: 0.95, ...flat });
  rel(s.Cp, s.Cq * s.tsr, 1e-12);
  rel(s.power, s.torque * 0.95, 1e-12);
  rel(s.thrust, s.Ct * 0.5 * 1.225 * Math.PI * 63 * 63 * 64, 1e-12);
  assert.ok(s.flapMoment > 0 && s.edgeMoment > 0);
  // root flap moment = integral Np (r - Rhub) dr: bounded by thrust per blade times the lever of the tip
  assert.ok(s.flapMoment < s.thrust / 3 * (63 - 1.5));
});

test('tilt, yaw, shear and tower shadow: sector averaging behaves physically', () => {
  const base = bemSolve({ rotor, U: 9, omega: 1.0, ...flat });
  const tilted = bemSolve({ rotor, U: 9, omega: 1.0 });
  assert.ok(tilted.power < base.power && tilted.power > 0.96 * base.power, 'tilt 5 deg + cone 2.5 deg cost a few per cent');
  const yawed = bemSolve({ rotor, U: 9, omega: 1.0, yawDeg: 20, ...flat });
  const c = Math.cos(20 * D2R);
  assert.ok(yawed.power / base.power < c && yawed.power / base.power > c ** 4, `yaw power ratio ${yawed.power / base.power}`);
  // yaw is symmetric for the rotor-average loads
  const yawedNeg = bemSolve({ rotor, U: 9, omega: 1.0, yawDeg: -20, ...flat });
  rel(yawedNeg.power, yawed.power, 2e-3);
  const sheared = bemSolve({ rotor, U: 9, omega: 1.0, shearExp: 0.2, ...flat, nSectors: 24 });
  rel(sheared.thrust, base.thrust, 0.02);
  // blade pointing up sees more wind than pointing down
  assert.ok(sheared.azimuthal.flapMoment[0] > sheared.azimuthal.flapMoment[12]);
  const shadow = bemSolve({ rotor, U: 9, omega: 1.0, towerShadow: true, ...flat, nSectors: 72 });
  assert.ok(shadow.azimuthal.bladeThrust[36] < shadow.azimuthal.bladeThrust[0], 'thrust dip in front of the tower');
  rel(shadow.azimuthal.bladeThrust[18], base.azimuthal.bladeThrust[0], 0.02);
});

test('parked rotor (omega = 0) is evaluated without induction and gives finite loads', () => {
  const s = bemSolve({ rotor, U: 50, omega: 0, pitchDeg: 90, ...flat });
  assert.ok(Number.isFinite(s.thrust) && Number.isFinite(s.torque));
  close(s.power, 0, 1e-12);
  // feathered blades: thrust far below the operating thrust at rated
  assert.ok(Math.abs(s.thrust) < 3e5);
  const z = bemSolve({ rotor, U: 0, omega: 1 });
  assert.equal(z.Cp, 0);
});

test('cpSurface and surfaceLookup: nodes reproduced, bilinear inside, clamped outside', () => {
  const sf = cpSurface({ rotor, tsr: [6, 7, 8], pitchDeg: [0, 2], bem: flat });
  const direct = bemSolve({ rotor, U: 10, omega: 7 * 10 / 63, pitchDeg: 2, ...flat });
  rel(sf.cp[1][1], direct.Cp, 1e-12); rel(sf.ct[1][1], direct.Ct, 1e-12);
  const l = surfaceLookup(sf, 7, 2);
  rel(l.cp, direct.Cp, 1e-12);
  const mid = surfaceLookup(sf, 6.5, 1);
  close(mid.cp, 0.25 * (sf.cp[0][0] + sf.cp[0][1] + sf.cp[1][0] + sf.cp[1][1]), 1e-12);
  close(surfaceLookup(sf, 20, 1).cp, 0.5 * (sf.cp[0][2] + sf.cp[1][2]), 1e-12);
  rel(mid.cq * 6.5, 0.25 * (sf.cq[0][0] + sf.cq[0][1] + sf.cq[1][0] + sf.cq[1][1]) * 6.5, 1e-12);
  assert.throws(() => cpSurface({ rotor, tsr: [3, 2], pitchDeg: [0, 1] }), RangeError);
});

test('input validation', () => {
  assert.throws(() => bemSolve({ rotor, U: -1, omega: 1 }), RangeError);
  assert.throws(() => bemSolve({ rotor, U: 10, omega: 1, rho: 0 }), RangeError);
  assert.throws(() => bemSolve({ rotor: {}, U: 10, omega: 1 }), RangeError);
  assert.throws(() => glauertIdealCp(0), RangeError);
});

test('parked, feathered rotor: no induction, drag-only thrust at zero yaw by strip integration, smooth variation with yaw', () => {
  const rotor = nrel5mw();
  const U = 70; const rho = 1.225;
  const r0 = bemSolve({ rotor, U, omega: 0, pitchDeg: 90 });
  assert.equal(r0.parked, true); assert.equal(r0.tsr, 0);
  // zero yaw, tilt and cone: the wind is normal to the rotor plane (phi = 90 deg), each element sits at
  // alpha = 90 deg − (twist + 90 deg) = −twist, and the force along the shaft is its drag plus nothing from lift:
  //   T = B ∫ ½ ρ U² c cd(−twist) dr  (trapezoidal through the stations, zero at hub and tip, as the solver integrates)
  const flat = bemSolve({ rotor, U, omega: 0, pitchDeg: 90, tiltDeg: 0, preconeDeg: 0 });
  const st = rotor.stations; let T = 0; let rp = rotor.hubRadius; let fp = 0;
  for (let i = 0; i <= st.length; i++) {
    const r = i < st.length ? st[i].r : rotor.tipRadius;
    const f = i < st.length ? 0.5 * rho * U * U * st[i].chord * st[i].polar.cdAt(-st[i].twistDeg * Math.PI / 180) : 0;
    T += 0.5 * (f + fp) * (r - rp); rp = r; fp = f;
  }
  rel(flat.thrust, rotor.blades * T, 1e-9, 'drag-only thrust');
  for (const s of flat.spanwise) { assert.equal(s.a, 0); assert.equal(s.ap, 0); assert.equal(s.F, 1); }
  // a feathered rotor is two orders of magnitude below the operating thrust coefficient
  assert.ok(flat.Ct > 0 && flat.Ct < 0.02, `Ct ${flat.Ct}`);
  // smooth in yaw: 2 deg steps never change Ct by more than a small fraction of the operating value (the momentum
  // iteration used to jump from 123 kN at 0 deg to 1990 kN at 8 deg and −96 kN at 30 deg)
  let prev = r0.Ct; let big = 0; const ct = [];
  for (let yaw = 2; yaw <= 90; yaw += 2) {
    const c = bemSolve({ rotor, U, omega: 0, pitchDeg: 90, yawDeg: yaw }).Ct;
    big = Math.max(big, Math.abs(c - prev)); prev = c; ct.push(c);
    assert.ok(Math.abs(c) < 0.02, `yaw ${yaw}: Ct ${c}`);
  }
  assert.ok(big < 0.003, `largest 2 deg step in Ct ${big}`);
  // thrust scales with U² and the answer is symmetric in the sign of the yaw error for an untilted, unconed rotor
  rel(bemSolve({ rotor, U: 35, omega: 0, pitchDeg: 90, yawDeg: 20 }).thrust, bemSolve({ rotor, U, omega: 0, pitchDeg: 90, yawDeg: 20 }).thrust / 4, 1e-9);
  rel(bemSolve({ rotor, U, omega: 0, pitchDeg: 90, yawDeg: -20, tiltDeg: 0, preconeDeg: 0 }).thrust, bemSolve({ rotor, U, omega: 0, pitchDeg: 90, yawDeg: 20, tiltDeg: 0, preconeDeg: 0 }).thrust, 1e-9);
  // slow idling stays on the same branch and close to the standstill value; an operating point is unaffected
  const idle = bemSolve({ rotor, U, omega: 0.02, pitchDeg: 90, yawDeg: 8 });
  assert.equal(idle.parked, true);
  assert.ok(Math.abs(idle.Ct - bemSolve({ rotor, U, omega: 0, pitchDeg: 90, yawDeg: 8 }).Ct) < 1e-3);
  const op = bemSolve({ rotor, U: 11.4, omega: 12.1 * Math.PI / 30, pitchDeg: 0 });
  assert.equal(op.parked, false); assert.ok(op.Ct > 0.6 && op.Ct < 0.95);
  // the branch can be switched off
  assert.equal(bemSolve({ rotor, U, omega: 0, pitchDeg: 90, parkedTsr: 0 }).parked, false);
});
