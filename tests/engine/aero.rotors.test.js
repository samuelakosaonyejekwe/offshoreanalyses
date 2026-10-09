import test from 'node:test';
import assert from 'node:assert/strict';
import {
  createRotor, nrel5mw, nrel5mwGenericPolars, scaleRotor, PUBLIC_REFERENCE_TURBINES, optimalSection,
  designOptimalRotor, tubularDistribution, towerDistribution,
} from '../../src/engine/aero/rotors.js';
import { cylinderPolar } from '../../src/engine/aero/polar.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} +- ${tol}`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${tol})`);

test('NREL 5-MW (Jonkman et al. 2009): gross properties and 17-node planform', () => {
  const r = nrel5mw();
  assert.equal(r.blades, 3);
  assert.equal(r.tipRadius, 63); assert.equal(r.hubRadius, 1.5); assert.equal(r.hubHeight, 90);
  assert.equal(r.ratedPower, 5e6); assert.equal(r.ratedWind, 11.4); assert.equal(r.cutIn, 3); assert.equal(r.cutOut, 25);
  assert.equal(r.gearboxRatio, 97); assert.equal(r.drivetrainEfficiency, 0.944);
  assert.equal(r.tiltDeg, 5); assert.equal(r.preconeDeg, 2.5);
  assert.equal(r.rotorMass + r.nacelleMass, 350000);
  assert.equal(r.stations.length, 17);
  // node widths 2.7333 x3, 4.1 x11, 2.7333 x3 tile the blade from 1.5 m to 63 m
  const dr = [2.7333, 2.7333, 2.7333, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 2.7333, 2.7333, 2.7333];
  let edge = 1.5;
  r.stations.forEach((s, i) => { close(s.r, edge + dr[i] / 2, 2e-3, `node ${i}`); edge += dr[i]; });
  close(edge, 63, 1e-3);
  assert.equal(r.stations[4].chord, 4.652);           // maximum chord
  assert.equal(r.stations[16].chord, 1.419);
  assert.equal(r.stations[0].twistDeg, 13.308); assert.equal(r.stations[16].twistDeg, 0.106);
  assert.deepEqual(r.stations.map((s) => s.thickness), [1, 1, 1, 0.4, 0.35, 0.35, 0.3, 0.25, 0.25, 0.21, 0.21, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18]);
  // rated tip speed 12.1 rpm x 63 m = 79.8 m/s ("80 m/s" in the report)
  close(r.ratedRotorSpeed * 63, 12.1 * Math.PI / 30 * 63, 1e-9);
  close(r.ratedRotorSpeed * 63, 80, 0.2);
  // drivetrain inertia = rotor + generator referred through the gearbox
  close(r.rotorInertia + 534.116 * 97 * 97, 43784700, 1e-3);
  assert.match(r.label, /GENERIC/);
});

test('NREL 5-MW tower: integrated mass 347,460 kg and base EI 614.3e9 N m2 with the effective wall thickness', () => {
  const r = nrel5mw();
  const d = towerDistribution(r, 41);
  rel(d.mass, 347460, 1e-3);
  rel(d.massPerLength[0], 5590.87, 1e-3);
  rel(d.EI[0], 614.34e9, 2e-3);
  close(d.cogZ, 38.234, 0.1);        // report: centre of mass at 38.234 m
  assert.equal(r.tower.nominalBaseThickness, 0.027);
});

test('generic polars for the reference rotor are generic by thickness, cylinders use Cd 0.50 / 0.35', () => {
  const p = nrel5mwGenericPolars();
  assert.equal(p.Cylinder1.cdAt(0.3), 0.5); assert.equal(p.Cylinder2.cdAt(0.3), 0.35);
  for (const id of ['DU40', 'DU35', 'DU30', 'DU25', 'DU21', 'NACA64']) assert.equal(p[id].generic, true);
  assert.equal(p.DU40.thickness, 0.40); assert.equal(p.NACA64.thickness, 0.18);
});

test('scaleRotor: sqrt(P) length scaling at constant tip speed and specific power, labelled as scaled', () => {
  const ref = nrel5mw();
  const r = scaleRotor(ref, { ratedPower: 10e6 });
  const s = Math.SQRT2;
  rel(r.tipRadius, 63 * s, 1e-12); rel(r.stations[5].chord, ref.stations[5].chord * s, 1e-12);
  rel(r.ratedRotorSpeed * r.tipRadius, ref.ratedRotorSpeed * ref.tipRadius, 1e-12);
  rel(r.ratedPower / r.sweptArea, ref.ratedPower / ref.sweptArea, 1e-12);
  rel(r.rnaMass, 350000 * s ** 3, 1e-12);
  rel(r.drivetrainInertia, 43784700 * s ** 5, 1e-12);
  rel(r.solidity, ref.solidity, 1e-12);
  assert.equal(r.stations[5].twistDeg, ref.stations[5].twistDeg);
  assert.match(r.label, /generic scaled rotor/);
  assert.match(r.label, /not a published reference turbine/);
  // blade-tip clearance to the datum is preserved by default
  close(r.hubHeight - r.tipRadius, 90 - 63, 1e-9);
  const r2 = scaleRotor(ref, { ratedPower: 15e6, massExponent: 2.2 });
  rel(r2.bladeMass, 17740 * 3 ** 1.1, 1e-12);
  assert.throws(() => scaleRotor(ref, { ratedPower: -1 }), RangeError);
});

test('public headline parameters of the IEA 15-MW and DTU 10-MW reference turbines are recorded', () => {
  const iea = PUBLIC_REFERENCE_TURBINES['IEA-15MW'];
  assert.equal(iea.rotorDiameter, 240); assert.equal(iea.hubHeight, 150); assert.equal(iea.ratedWind, 10.59);
  assert.equal(iea.ratedRotorSpeedRpm, 7.56); assert.equal(iea.rnaMass, 1017000);
  const dtu = PUBLIC_REFERENCE_TURBINES['DTU-10MW'];
  assert.equal(dtu.rotorDiameter, 178.3); assert.equal(dtu.hubHeight, 119); assert.equal(dtu.ratedWind, 11.4); assert.equal(dtu.ratedRotorSpeedRpm, 9.6);
  // a 15-MW rotor scaled from the 5-MW machine is a different machine from the IEA 15-MW (218 m vs 240 m)
  const s15 = scaleRotor(nrel5mw(), { ratedPower: 15e6 });
  close(2 * s15.tipRadius, 126 * Math.sqrt(3), 1e-9);
  assert.ok(Math.abs(2 * s15.tipRadius - iea.rotorDiameter) > 15);
});

test('optimalSection: Betz and Glauert/Schmitz optimum by hand (lambda_r = 4, B = 3, Cl = 1, r = 10 m)', () => {
  const b = optimalSection({ r: 10, lambdaR: 4, blades: 3, clDesign: 1, wakeRotation: false });
  // phi = atan(2/12) = 9.4623 deg; c = 8 pi 10 sin(phi) / (3*3*1*4) = 1.1477 m
  close(b.phi * 180 / Math.PI, 9.4623, 1e-4);
  close(b.chord, 1.1477, 1e-4);
  const g = optimalSection({ r: 10, lambdaR: 4, blades: 3, clDesign: 1 });
  // phi = (2/3) atan(1/4) = 9.3575 deg; c = 8 pi 10 (1 - cos phi)/3 = 1.1148 m
  close(g.phi * 180 / Math.PI, 9.3575, 1e-4);
  close(g.chord, 1.1148, 2e-4);
  assert.throws(() => optimalSection({ r: 1, lambdaR: 0, blades: 3, clDesign: 1 }), RangeError);
});

test('designOptimalRotor: twist = phi - design alpha, stations at annulus mid-radius', () => {
  const polar = cylinderPolar({ cd: 0.5 });
  const r = designOptimalRotor({ tipRadius: 10, hubRadius: 1, blades: 3, tsr: 5, clDesign: 1, alphaDesignDeg: 6, polar, nStations: 9 });
  assert.equal(r.stations.length, 9);
  close(r.stations[0].r, 1.5, 1e-12); close(r.stations[8].r, 9.5, 1e-12);
  const phi = (2 / 3) * Math.atan(1 / (5 * 9.5 / 10)) * 180 / Math.PI;
  close(r.stations[8].twistDeg, phi - 6, 1e-9);
  assert.equal(r.designTsr, 5);
});

test('tubularDistribution: uniform tube mass rho pi (D - t) t L and I = pi/64 (D^4 - d^4)', () => {
  const d = tubularDistribution({ zBase: -30, zTop: 10, dBase: 6, dTop: 6, tBase: 0.06, tTop: 0.06, n: 5 });
  rel(d.mass, 7850 * Math.PI * 5.94 * 0.06 * 40, 1e-12);
  rel(d.I[2], Math.PI / 64 * (6 ** 4 - 5.88 ** 4), 1e-12);
  close(d.cogZ, -10, 1e-9);
  assert.throws(() => tubularDistribution({ zBase: 0, zTop: 10, dBase: 1, dTop: 1, tBase: 0.6, tTop: 0.01 }), RangeError);
});

test('createRotor: validation of stations and polars', () => {
  const p = cylinderPolar();
  assert.throws(() => createRotor({ blades: 3, tipRadius: 10, hubRadius: 1, stations: [{ r: 5, chord: 1, polar: p }] }), RangeError);
  assert.throws(() => createRotor({ blades: 3, tipRadius: 10, hubRadius: 1, stations: [{ r: 5, chord: 1, polar: p }, { r: 4, chord: 1, polar: p }] }), RangeError);
  assert.throws(() => createRotor({ blades: 3, tipRadius: 10, hubRadius: 1, stations: [{ r: 5, chord: 1, airfoil: 'x' }, { r: 6, chord: 1, airfoil: 'x' }] }), RangeError);
  const r = createRotor({ blades: 2, tipRadius: 10, hubRadius: 0, stations: [{ r: 3, chord: 1, polar: p }, { r: 7, chord: 1, polar: p }] });
  // constant chord 1 m: solidity = 2 * 10 / (pi * 100)
  close(r.solidity, 20 / (Math.PI * 100), 1e-12);
});
