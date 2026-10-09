import test from 'node:test';
import assert from 'node:assert/strict';
import {
  GENERIC_CABLE_MATERIALS, GENERIC_TELECOM_CABLES, buildCable, cablePullingTension, cableHandlingLimits,
  cableImpactClass, cableCoiling, cableLayTension, cableOnBottomStability, conductorThermalThrust,
  telecomTensionUtilisation, telecomSinkAngle, telecomSlack,
} from '../../src/engine/xsection/cablemech.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);
const DEG = Math.PI / 180;
const M = GENERIC_CABLE_MATERIALS;

const core = [
  { type: 'conductor', area: 630e-6, fillFactor: 0.9, ...M.copper },
  { type: 'insulation', thickness: 0.018, ...M.xlpe },
  { type: 'sheath', thickness: 0.003, ...M.lead },
  { type: 'oversheath', thickness: 0.003, ...M.pe },
];
const armour = (lay, d = 0.005) => ({ type: 'armour', wireDiameter: d, alphaDeg: 12, lay, fill: 0.95, ...M.steel, uts: 1200e6 });

test('single-core cable: diameter, mass, EA and submerged weight by hand sums', () => {
  const c = buildCable({ cores: 1, core, outer: [armour(1)] });
  const rc = Math.sqrt(630e-6 / (Math.PI * 0.9));
  const r1 = rc + 0.018, r2 = r1 + 0.003, r3 = r2 + 0.003;
  const R = r3 + 0.0025, ca = Math.cos(12 * DEG);
  const n = Math.floor(0.95 * 2 * Math.PI * R * ca / 0.005 + 1e-9), Aw = Math.PI * 0.005 ** 2 / 4;
  close(c.outerDiameter, 2 * (r3 + 0.005), 1e-12);
  const mass = 8900 * 630e-6 + 920 * Math.PI * (r1 ** 2 - rc ** 2) + 11340 * Math.PI * (r2 ** 2 - r1 ** 2) + 950 * Math.PI * (r3 ** 2 - r2 ** 2) + 7850 * n * Aw / ca;
  close(c.mass, mass, 1e-12);
  const EA = 1.15e11 * 630e-6 + 2e8 * Math.PI * (r1 ** 2 - rc ** 2) + 1.6e10 * Math.PI * (r2 ** 2 - r1 ** 2) + 8e8 * Math.PI * (r3 ** 2 - r2 ** 2) + n * 2e11 * Aw * ca ** 3;
  close(c.EAfixed, EA, 1e-12);
  close(c.submergedWeight, (mass - 1025 * Math.PI * c.outerDiameter ** 2 / 4) * 9.80665, 1e-12);
  close(c.specificGravity, mass / (1025 * Math.PI * c.outerDiameter ** 2 / 4), 1e-12);
  close(c.coupling, n * 2e11 * Aw * R * Math.sin(12 * DEG) * ca ** 2, 1e-12);
  assert.equal(c.torqueBalanced, false);
  assert.ok(c.EImax > c.EImin && c.EImin > 0 && c.EA < c.EAfixed);
  assert.equal(M.copper.generic, true);
});

test('cross-wound double armour is close to torque balance; identical opposite layers at one radius balance exactly', () => {
  const single = buildCable({ cores: 1, core, outer: [armour(1)] });
  const dbl = buildCable({ cores: 1, core, outer: [armour(1), armour(-1)] });
  assert.ok(Math.abs(dbl.torquePerTension) < 0.25 * Math.abs(single.torquePerTension));
  // by hand: coupling of the pair = k sin a cos^2 a (n1 R1 - n2 R2)
  const a1 = dbl.layers.find((L) => L.type === 'armour' && L.lay === 1), a2 = dbl.layers.find((L) => L.lay === -1);
  const kk = 2e11 * a1.wireArea * Math.sin(12 * DEG) * Math.cos(12 * DEG) ** 2;
  close(dbl.coupling, kk * (a1.nWires * a1.R - a2.nWires * a2.R), 1e-10);
  // sum of the standard coupling terms with opposite hand cancels when n R is equal
  const n = 40, R = 0.03;
  const term = (lay) => lay * n * 2e11 * 2e-5 * R * Math.sin(0.2) * Math.cos(0.2) ** 2;
  assert.equal(term(1) + term(-1), 0);
});

test('three-core lay-up geometry: pitch radius d/sqrt(3), envelope d(1 + 2/sqrt(3))', () => {
  const c = buildCable({ cores: 3, core, coreLayLength: 2.0, fillerDensity: 600, fibreUnit: { mass: 0.1 }, outer: [{ type: 'bedding', thickness: 0.002, ...M.pe }, armour(-1, 0.006)] });
  const d = c.coreDiameter;
  // three mutually touching circles: centre distance from the axis = d/sqrt(3); check the centres are d apart
  const Rc = d / Math.sqrt(3);
  close(2 * Rc * Math.sin(Math.PI / 3), d, 1e-12);
  const env = d * (1 + 2 / Math.sqrt(3));
  close(c.outerDiameter, env + 2 * (0.002 + 0.006), 1e-12);
  close(c.conductorArea, 3 * 630e-6, 1e-12);
  const a = Math.atan(2 * Math.PI * Rc / 2.0);
  const one = buildCable({ cores: 1, core });
  const filler = 600 * (Math.PI * env ** 2 / 4 - 3 * Math.PI * d ** 2 / 4);
  const outer = c.layers.filter((L) => L.group === 'outer').reduce((s, L) => s + L.mass, 0);
  close(c.mass, 3 * one.mass / Math.cos(a) + filler + 0.1 + outer, 1e-12);
});

test('pulling tension: 50 N/mm2 on copper, 30 N/mm2 on aluminium, armour-based by hand', () => {
  const c = buildCable({ cores: 3, core, coreLayLength: 2.0, outer: [armour(1)] });
  close(cablePullingTension(c).conductorLimit, 50e6 * 3 * 630e-6, 1e-12); // 94.5 kN
  close(cablePullingTension(c, { conductor: 'aluminium' }).conductorLimit, 30e6 * 3 * 630e-6, 1e-12);
  const arm = c.layers.find((L) => L.type === 'armour');
  const p = cablePullingTension(c, { armourUtilisation: 0.3 });
  close(p.armourLimit, 0.3 * arm.nWires * arm.wireArea * 1200e6 * Math.cos(12 * DEG), 1e-12);
});

test('handling limits: MBR multiples, sidewall pressure T/R, grip length', () => {
  const h = cableHandlingLimits({ outerDiameter: 0.15, tension: 80e3, bendRadius: 4, sidewallLimit: 25e3, squeeze: 20e3, squeezeLimit: 30e3, nTracks: 4, mu: 0.15 });
  close(h.mbrInstall, 3.0); close(h.mbrStatic, 1.8); close(h.bendUtilisation, 0.75);
  close(h.sidewallPressure, 20e3); close(h.sidewallUtilisation, 0.8); close(h.minRadiusForSidewall, 3.2);
  close(h.gripLength, 1.5 * 80e3 / (0.15 * 4 * 20e3)); // 10 m
  close(cableHandlingLimits({ outerDiameter: 0.15, installFactor: 15, staticFactor: 10 }).mbrInstall, 2.25);
  const bands = [{ maxEnergy: 500, label: 'no damage' }, { maxEnergy: 2000, label: 'sheath damage' }];
  assert.equal(cableImpactClass({ energy: 400, bands }).label, 'no damage');
  assert.equal(cableImpactClass({ energy: 1500, bands }).band, 1);
  assert.equal(cableImpactClass({ energy: 5000, bands }).exceedsAll, true);
});

test('coiling: unidirectional armour is coilable, cross-wound is not; twist rate 1/R', () => {
  const s = buildCable({ cores: 1, core, outer: [armour(1)] });
  const d = buildCable({ cores: 1, core, outer: [armour(1), armour(-1)] });
  const cs = cableCoiling(s, { coilRadius: 4 });
  assert.equal(cs.coilable, true); close(cs.twistRate, 0.25);
  const L = s.layers.find((x) => x.type === 'armour');
  close(Math.abs(cs.armour[0].strain), L.R * Math.sin(12 * DEG) * Math.cos(12 * DEG) / 4, 1e-12);
  assert.equal(cableCoiling(d, { coilRadius: 4 }).coilable, false);
});

test('laying tension: catenary identities T = H + w d, s^2 = d^2 + 2 d H/w, and the 1.3 w d + H test tension', () => {
  const w = 300, d = 120, H = 15e3;
  const r = cableLayTension({ submergedWeight: w, depth: d, bottomTension: H });
  close(r.staticTopTension, H + w * d);
  // independent catenary check: a = H/w, y = a (cosh(x/a) - 1) reaches d at the layback x, with s = a sinh(x/a)
  const a = H / w;
  close(a * (Math.cosh(r.layback / a) - 1), d, 1e-12);
  close(a * Math.sinh(r.layback / a), r.suspendedLength, 1e-12);
  // top tension equals sqrt(H^2 + (w s)^2) and the departure angle satisfies tan = w s/H
  close(Math.hypot(H, w * r.suspendedLength), r.staticTopTension, 1e-12);
  close(Math.tan(r.departureAngle), w * r.suspendedLength / H, 1e-12);
  close(r.designTension, 1.3 * w * d + H);
  close(r.touchdownRadius, 50);
});

test('on-bottom stability by hand and specific gravity', () => {
  const r = cableOnBottomStability({ submergedWeight: 400, outerDiameter: 0.15, velocity: 1.0, mu: 0.5, Cd: 1.0, Cl: 0.9, safetyFactor: 1.1, mass: 60 });
  const q = 0.5 * 1025 * 0.15; // 76.875 N/m
  close(r.drag, q); close(r.lift, 0.9 * q); close(r.resistance, 0.5 * (400 - 0.9 * q));
  close(r.requiredWeight, 0.9 * q + 1.1 * q / 0.5);
  assert.equal(r.stable, true);
  // at the critical velocity the utilisation is exactly one
  const u = cableOnBottomStability({ submergedWeight: 400, outerDiameter: 0.15, velocity: r.criticalVelocity, mu: 0.5, Cd: 1.0, Cl: 0.9 });
  close(u.utilisation, 1, 1e-12);
  close(r.specificGravity, 60 / (1025 * Math.PI * 0.15 ** 2 / 4));
});

test('thermal thrust F = E A alpha dT and snake amplitude from the sine arc length', () => {
  const r = conductorThermalThrust({ area: 1000e-6, Eeff: 4e10, alphaT: 1.7e-5, deltaT: 70, snakeWavelength: 6 });
  close(r.thrust, 4e10 * 1e-3 * 1.7e-5 * 70); // 47.6 kN
  // numerical arc length of y = a sin(2 pi x / L) over one wavelength exceeds L by the free strain
  const a = r.snakeAmplitude, L = 6, m = 200000;
  let s = 0;
  for (let i = 0; i < m; i++) { const x = (i + 0.5) * L / m; s += Math.hypot(1, a * 2 * Math.PI / L * Math.cos(2 * Math.PI * x / L)) * L / m; }
  close(s / L - 1, 1.7e-5 * 70, 2e-3);
  close(r.snakeCurvature, a * (2 * Math.PI / 6) ** 2);
});

test('telecom: tension utilisation, hydrodynamic constant and sink angle (normal force balance)', () => {
  const c = GENERIC_TELECOM_CABLES.LW;
  assert.equal(c.generic, true); assert.match(c.label, /GENERIC/);
  const u = telecomTensionUtilisation({ transientTension: 35e3, operatingTension: 30e3, permanentTension: 6e3, ntts: c.ntts, nots: c.nots, npts: c.npts });
  close(u.transient, 0.5); close(u.operating, 0.6); close(u.permanent, 0.2); close(u.governing, 0.6);
  const V = 3.0, Cd = 1.8;
  const s = telecomSinkAngle({ weightWater: c.weightWater, outerDiameter: c.outerDiameter, shipSpeed: V, Cd, depth: 3000 });
  close(s.hydrodynamicConstant, Math.sqrt(2 * 3.0 / (1025 * 1.8 * 0.017)), 1e-12);
  // exact angle satisfies w cos(a) = 0.5 rho Cd D (V sin a)^2
  close(c.weightWater * Math.cos(s.sinkAngle), 0.5 * 1025 * Cd * c.outerDiameter * (V * Math.sin(s.sinkAngle)) ** 2, 1e-10);
  close(s.sinkAngle, s.sinkAngleSmall, 0.02); // small-angle form alpha = H/V
  close(s.layback, 3000 / Math.tan(s.sinkAngle));
  close(s.hydrodynamicConstantDegKnots, s.hydrodynamicConstant * 180 / Math.PI / 0.514444, 1e-12);
});

test('telecom slack: flat seabed Vc/Vs = 1 + e; sloping seabed against explicit triangle geometry', () => {
  close(telecomSlack({ sinkAngle: 0.15, shipSpeed: 3, slack: 0.03 }).payoutSpeed, 3.09, 1e-12);
  close(telecomSlack({ sinkAngle: 0.15, shipSpeed: 3, payoutSpeed: 3.09 }).slack, 0.03, 1e-10);
  // explicit geometry: ship at x, cable straight at angle a, seabed y = -(h0 + x tan b)
  const a = 0.2, b = 0.08, h0 = 1000;
  const touch = (xs) => { const l = (h0 + xs * Math.tan(b)) / (Math.sin(a) + Math.cos(a) * Math.tan(b)); return { l, xt: xs - l * Math.cos(a) }; };
  const p0 = touch(0), p1 = touch(100);
  const seabedLength = (p1.xt - p0.xt) / Math.cos(b), dl = p1.l - p0.l;
  const e = 0.02;
  const need = ((1 + e) * seabedLength + dl) / 100;
  const r = telecomSlack({ sinkAngle: a, seabedSlope: b, shipSpeed: 2, slack: e });
  close(r.speedRatio, need, 1e-10);
  close(r.zeroSlackRatio, (seabedLength + dl) / 100, 1e-10);
  // rising seabed needs less cable than ship distance
  assert.ok(telecomSlack({ sinkAngle: a, seabedSlope: -0.08, shipSpeed: 2, slack: 0 }).speedRatio < 1);
  assert.throws(() => telecomSlack({ sinkAngle: a, shipSpeed: 2 }), RangeError);
});
