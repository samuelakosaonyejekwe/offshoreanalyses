import test from 'node:test';
import assert from 'node:assert/strict';
import {
  fluxDensity, arrangement, twistAttenuation, fluxDensityProfile, distanceToThreshold, dcFieldWithGeomagnetic,
  motionalElectricField, faradayElectricField, mutualImpedance, inducedVoltage, sheathStandingVoltage, svlSelection,
  acCorrosionCurrentDensity, circulatingCurrent,
} from '../../src/engine/electrical/emf.js';

const rel = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${a} vs ${b} (tol ${tol})`);
const MU0 = 4e-7 * Math.PI;

test('Biot–Savart, single long conductor: B = μ0·I/(2πr): 1000 A at 1 m → 200 µT, direction by the right-hand rule', () => {
  const b = fluxDensity({ conductors: [{ x: 0, z: 0, I: 1000 }], point: { x: 0, z: 1 } });
  rel(b.rms, 200e-6, 1e-12);
  // current along +y, point above (+z): field along +x
  rel(b.Bx.re, 200e-6, 1e-12); assert.ok(Math.abs(b.Bz.re) < 1e-20);
  const side = fluxDensity({ conductors: [{ x: 0, z: 0, I: 1000 }], point: { x: 2, z: 0 } });
  rel(side.Bz.re, -100e-6, 1e-12);
  rel(fluxDensity({ conductors: [{ x: 0, z: 0, I: 1000 }], point: { x: 3, z: 4 } }).rms, (MU0 * 1000) / (2 * Math.PI * 5), 1e-12);
});

test('go-and-return pair (DC bipole): far field μ0·I·s/(2πr²); exact on the perpendicular bisector μ0·I·s/(2π(r² + s²/4))', () => {
  const c = arrangement({ arrangement: 'bipole', I: 1500, s: 0.2, z: 0 });
  const b = fluxDensity({ conductors: c, point: { x: 0, z: 5 } });
  rel(b.rms, (MU0 * 1500 * 0.2) / (2 * Math.PI * (25 + 0.01)), 1e-12);
  rel(fluxDensity({ conductors: c, point: { x: 0, z: 50 } }).rms, (MU0 * 1500 * 0.2) / (2 * Math.PI * 2500), 1e-4);
  // separated poles produce a larger field than bundled ones
  const wide = arrangement({ arrangement: 'bipole', I: 1500, s: 20, z: 0 });
  assert.ok(fluxDensity({ conductors: wide, point: { x: 0, z: 5 } }).rms > 10 * b.rms);
});

test('balanced three-phase trefoil: far field decays as 1/r² and equals √6·μ0·I·s/(4πr²) (rms) in every direction', () => {
  const c = arrangement({ arrangement: 'trefoil', I: 1000, s: 0.1, z: 0 });
  const at = (r, ang) => fluxDensity({ conductors: c, point: { x: r * Math.cos(ang), z: r * Math.sin(ang) } }).rms;
  rel(at(20, 0.3) / at(10, 0.3), 0.25, 1e-3);
  for (const ang of [0, 0.7, 1.9, 3.3]) rel(at(20, ang), (Math.sqrt(6) * MU0 * 1000 * 0.1) / (4 * Math.PI * 400), 1e-3);
  // the phasors sum to zero: no monopole term, so the field falls faster than that of one conductor
  const sum = c.reduce((s, k) => ({ re: s.re + k.I.re, im: s.im + k.I.im }), { re: 0, im: 0 });
  assert.ok(Math.hypot(sum.re, sum.im) < 1e-9);
  assert.ok(at(10, 1) < 0.02 * (MU0 * 1000) / (2 * Math.PI * 10));
});

test('balanced three-phase flat formation: far field √3·μ0·I·s/(2πr²) (rms), √2 larger than trefoil at the same spacing', () => {
  const c = arrangement({ arrangement: 'flat', I: 1000, s: 0.1, z: 0 });
  const b = fluxDensity({ conductors: c, point: { x: 0, z: 30 } });
  rel(b.rms, (Math.sqrt(3) * MU0 * 1000 * 0.1) / (2 * Math.PI * 900), 1e-3);
  assert.ok(b.major >= b.minor && b.major <= b.rms + 1e-20);
  rel(Math.hypot(b.major, b.minor), b.rms, 1e-12);
  const t = fluxDensity({ conductors: arrangement({ arrangement: 'trefoil', I: 1000, s: 0.1 }), point: { x: 0, z: 30 } });
  rel(b.rms / t.rms, Math.SQRT2, 2e-3);
});

test('profile at the seabed for a buried cable peaks above the cable; shielding factor scales the field', () => {
  const c = arrangement({ arrangement: 'trefoil', I: 800, s: 0.08, z: -1.5 });
  const p = fluxDensityProfile({ conductors: c, height: 0, xMin: -10, xMax: 10, n: 201 });
  assert.ok(Math.abs(p.xMax) < 0.3);
  rel(p.Bmax, (Math.sqrt(6) * MU0 * 800 * 0.08) / (4 * Math.PI * 1.5 ** 2), 0.02);
  const sh = fluxDensityProfile({ conductors: c, height: 0, xMin: -10, xMax: 10, n: 201, shieldingFactor: 0.8 });
  rel(sh.Bmax / p.Bmax, 0.8, 1e-12);
  // 1 m above the seabed the field is lower
  assert.ok(fluxDensityProfile({ conductors: c, height: 1 }).Bmax < p.Bmax);
});

test('distance to a threshold: single conductor r = μ0·I/(2πB): 1000 A, 100 µT → 2 m; trefoil from the far-field law', () => {
  rel(distanceToThreshold({ conductors: [{ x: 0, z: 0, I: 1000 }], threshold: 100e-6 }), 2, 1e-6);
  const c = arrangement({ arrangement: 'trefoil', I: 1000, s: 0.1 });
  const d = distanceToThreshold({ conductors: c, threshold: 1e-6 });
  rel(d, Math.sqrt((Math.sqrt(6) * MU0 * 1000 * 0.1) / (4 * Math.PI * 1e-6)), 2e-3);
  // a threshold never reached outside the group envelope (ray pointing down, between the two lower cables)
  assert.equal(distanceToThreshold({ conductors: c, threshold: 1, directionDeg: -90 }), 0);
});

test('twist attenuation factor: 1 at the core circle, exp(−2π) one pitch further out', () => {
  rel(twistAttenuation({ r: 0.05, pitch: 2, coreRadius: 0.05 }), 1, 1e-12);
  rel(twistAttenuation({ r: 2.05, pitch: 2, coreRadius: 0.05 }), Math.exp(-2 * Math.PI), 1e-12);
});

test('DC cable and geomagnetic field: 1000 A monopole heading north, 5 m above → 40 µT eastwards; compass deviation atan(40/21.13) = 62.2°', () => {
  const r = dcFieldWithGeomagnetic({ conductors: [{ x: 0, z: 0, I: 1000 }], point: { x: 0, z: 5 }, headingDeg: 0, geomagneticTotal: 50e-6, inclinationDeg: 65 });
  rel(r.cable, 40e-6, 1e-12);
  const H = 50e-6 * Math.cos((65 * Math.PI) / 180);
  rel(r.compassDeviationDeg, (Math.atan2(40e-6, H) * 180) / Math.PI, 1e-9);
  rel(r.compassDeviationDeg, 62.15, 1e-3);
  rel(r.total, Math.hypot(H, 40e-6, 50e-6 * Math.sin((65 * Math.PI) / 180)), 1e-12);
  // a cable running east–west adds to or subtracts from the horizontal field without deviating the compass
  const ew = dcFieldWithGeomagnetic({ conductors: [{ x: 0, z: 0, I: 100 }], point: { x: 0, z: 5 }, headingDeg: 90 });
  assert.ok(Math.abs(ew.compassDeviationDeg) < 1e-9 || Math.abs(Math.abs(ew.compassDeviationDeg) - 180) < 1e-9);
  // bundled bipole: negligible deviation at the same distance
  const bi = dcFieldWithGeomagnetic({ conductors: arrangement({ arrangement: 'bipole', I: 1000, s: 0.15 }).map((c) => ({ x: c.x, z: c.z, I: c.I.re })), point: { x: 0, z: 5 }, headingDeg: 0 });
  assert.ok(Math.abs(bi.compassDeviationDeg) < 1 || Math.abs(bi.compassDeviationDeg) < r.compassDeviationDeg / 10);
});

test('induced electric fields: E = v·B (1 m/s in 50 µT → 50 µV/m); Faraday loop E = π·f·a·B', () => {
  rel(motionalElectricField({ velocity: 1, B: 50e-6 }), 50e-6, 1e-12);
  rel(motionalElectricField({ velocity: 2, B: 50e-6, angleDeg: 30 }), 50e-6, 1e-12);
  // loop of radius a: EMF = ω·B·π·a², divided by the circumference 2πa
  rel(faradayElectricField({ B: 10e-6, f: 50, radius: 0.1 }), (2 * Math.PI * 50 * 10e-6 * Math.PI * 0.01) / (2 * Math.PI * 0.1), 1e-12);
});

test('Carson mutual impedance by hand: f = 50 Hz, ρ = 100 Ω·m, d = 10 m → 4.935e-5 + j2.849e-4 Ω/m', () => {
  const { Zm, De } = mutualImpedance({ d: 10, f: 50, rhoEarth: 100 });
  rel(De, 931.78, 1e-4);
  rel(Zm.re, 4.9348e-5, 1e-5);
  rel(Zm.im, 2 * Math.PI * 50 * 2e-7 * Math.log(93.178), 1e-5);
  rel(Zm.im, 2.84912e-4, 1e-4);
  assert.throws(() => mutualImpedance({ d: 2000, f: 50, rhoEarth: 100 }), RangeError);
});

test('induced voltage on a parallel pipeline: single current uses the full Zm; a balanced three-phase group largely cancels', () => {
  const one = inducedVoltage({ conductors: [{ d: 10, I: 1000 }], f: 50, rhoEarth: 100, length: 5000 });
  rel(one.E, 1000 * Math.hypot(4.9348e-5, 2.84912e-4), 1e-4);
  rel(one.voltageInsulated, one.E * 5000, 1e-12);
  const a = 2 * Math.PI / 3;
  const three = inducedVoltage({ conductors: [
    { d: 9.9, I: { re: 1000, im: 0 } }, { d: 10, I: { re: 1000 * Math.cos(-a), im: 1000 * Math.sin(-a) } },
    { d: 10.1, I: { re: 1000 * Math.cos(a), im: 1000 * Math.sin(a) } }], f: 50, rhoEarth: 100, length: 5000 });
  assert.ok(three.E < 0.01 * one.E);
  // long-line limit and screening factor
  const long = inducedVoltage({ conductors: [{ d: 10, I: 1000 }], f: 50, rhoEarth: 100, length: 50e3, gamma: 2e-4, screeningFactor: 0.5 });
  rel(long.voltageLongLine, (0.5 * one.E) / 2e-4, 1e-12);
  assert.equal(long.voltage, Math.min(long.voltageInsulated, long.voltageLongLine));
  rel(circulatingCurrent({ voltage: 100, Zloop: { re: 3, im: 4 } }), 20, 1e-12);
});

test('sheath standing voltage, trefoil, equals X·I of IEC 60287: ω·2e-7·I·ln(2s/d) = 57.6 V/km at 1000 A, s = 100 mm, d = 80 mm', () => {
  const c = arrangement({ arrangement: 'trefoil', I: 1000, s: 0.1 });
  const v = sheathStandingVoltage({ conductors: c, sheathMeanDiameter: 0.08, f: 50, length: 1000 });
  const expected = 2 * Math.PI * 50 * 2e-7 * 1000 * Math.log(0.2 / 0.08);
  for (const e of v.perMetre) rel(e, expected, 1e-10);
  rel(v.voltage, 57.57, 1e-3);
  // flat formation: the outer cables see a higher voltage than the centre one, which sees X·I
  const f = sheathStandingVoltage({ conductors: arrangement({ arrangement: 'flat', I: 1000, s: 0.1 }), sheathMeanDiameter: 0.08, f: 50 });
  rel(f.perMetre[1], expected, 1e-10);
  assert.ok(f.perMetre[0] > f.perMetre[1]);
  const svl = svlSelection({ faultGradient: 20 * v.max, length: 800, margin: 1.2, ratings: [3000, 6000, 9000] });
  rel(svl.required, 20 * v.max * 800 * 1.2, 1e-12); assert.equal(svl.selected, 3000);
});

test('AC corrosion criterion: J = 8·V/(ρ·π·d): 10 V, 10 Ω·m, 1 cm² holiday → 225.7 A/m² > 30 A/m²', () => {
  const c = acCorrosionCurrentDensity({ Vac: 10, rhoSoil: 10 });
  rel(c.J, 80 / (10 * Math.PI * 0.01128), 1e-12);
  rel(c.J, 225.75, 1e-3);
  assert.equal(c.ok, false);
  rel(c.VacAllowable, 10 * 30 / c.J, 1e-12);
  // spreading resistance of the disc holiday R = ρ/(2d): J = V/(R·A)
  rel(c.J, 10 / ((10 / (2 * 0.01128)) * (Math.PI * 0.01128 ** 2 / 4)), 1e-12);
});
