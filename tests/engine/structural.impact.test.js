import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  collisionEnergy, collisionStrainEnergy, tubeDentResistance, tubeDentDepth, beamMechanismResistance, beamMechanismDeflection,
  sharedEnergy, classifyCollisionDesign, dropInAir, terminalVelocityWater, dropThroughWater, plateImpactCapacity, pipeDentEnergy,
  pipeDentFromEnergy, pileImpact,
} from '../../src/engine/structural/impact.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);

test('collision energy ½(m + a)v²: 5000 t supply vessel at 2 m/s — 14 MJ broadside (a = 0.4 m), 11 MJ bow/stern (a = 0.1 m)', () => {
  // the classic design values of NORSOK N-003 / DNV-RP-C204
  const sway = collisionEnergy({ mass: 5e6, velocity: 2, addedMassFactor: 1.4 });
  const surge = collisionEnergy({ mass: 5e6, velocity: 2, addedMassFactor: 1.1 });
  close(sway.energy, 14e6, 1e-12); close(surge.energy, 11e6, 1e-12);
  close(sway.addedMass, 2e6, 1e-12); close(surge.momentum, 11e6, 1e-12);
  assert.throws(() => collisionEnergy({ mass: 5e6, velocity: 2, addedMassFactor: 0.9 }), RangeError);
  assert.throws(() => collisionEnergy({ mass: 5e6, velocity: 2 }), RangeError);
});

test('collision strain energy: fixed and compliant installation from conservation of momentum and energy', () => {
  const fixed = collisionStrainEnergy({ shipMass: 7e6, shipVelocity: 2 });
  close(fixed.strainEnergy, 14e6, 1e-12);
  // compliant: ms = 7000 t, mi = 21 000 t at rest → common velocity 0.5 m/s; Es = KE before − KE after = 14 − ½·28e6·0.25 = 10.5 MJ
  const c = collisionStrainEnergy({ shipMass: 7e6, shipVelocity: 2, installationMass: 21e6 });
  close(c.commonVelocity, 0.5, 1e-12); close(c.strainEnergy, 10.5e6, 1e-12);
  close(c.strainEnergy, 14e6 / (1 + 1 / 3), 1e-12);
  // moving installation: energy balance closes for any velocities
  const m = collisionStrainEnergy({ shipMass: 7e6, shipVelocity: 2, installationMass: 4e6, installationVelocity: 0.5 });
  close(m.strainEnergy, m.kineticEnergyBefore - 0.5 * 11e6 * m.commonVelocity ** 2, 1e-12);
  close(m.strainEnergy, 0.5 * 7e6 * 4 * (1 - 0.25) ** 2 / (1 + 7 / 4), 1e-12);
});

test('tubular denting, Furnes–Amdahl form R/Rc = c (wd/D)^n: force and energy by hand; API RP 2A constants reproduce E = 14.14 Mp X^1.5/t^0.5', () => {
  // D = 1.0 m, t = 25 mm, fy = 355 MPa: Rc = 355e6·0.025²/4·√40 = 350 812 N; wd = 0.1 m → R = 22·Rc·√0.1 = 2.4406 MN
  const r = tubeDentResistance({ D: 1, t: 0.025, fy: 355e6, dent: 0.1 });
  close(r.Rc, 350812, 1e-5); close(r.R, 2.4406e6, 1e-4);
  // E = ∫R dw = Rc·D·22·(2/3)·0.1^1.5 = 162.7 kJ
  close(r.energy, 350812 * 22 * 2 / 3 * 0.1 ** 1.5, 1e-5); close(r.energy, 162.7e3, 1e-3);
  // numerical integration of the force
  let e = 0; const n = 20000;
  for (let i = 0; i < n; i++) e += tubeDentResistance({ D: 1, t: 0.025, fy: 355e6, dent: 0.1 * (i + 0.5) / n }).R * 0.1 / n;
  close(e, r.energy, 1e-5);
  // API RP 2A commentary: Pd = 15 Mp (D/t)^½ (X/R)^½ ≡ c = 15√2; energy 14.14 Mp X^1.5/t^0.5 with Mp = fy t²/4
  const X = 0.08; const Mp = 355e6 * 0.025 ** 2 / 4;
  const api = tubeDentResistance({ D: 1, t: 0.025, fy: 355e6, dent: X, coefficient: 15 * Math.SQRT2 });
  close(api.R, 15 * Mp * Math.sqrt(40) * Math.sqrt(X / 0.5), 1e-12);
  close(api.energy, 14.14 * Mp * X ** 1.5 / Math.sqrt(0.025), 2e-4);
  // inverses
  close(tubeDentDepth({ D: 1, t: 0.025, fy: 355e6, energy: r.energy }).dent, 0.1, 1e-12);
  close(tubeDentDepth({ D: 1, t: 0.025, fy: 355e6, force: r.R }).dent, 0.1, 1e-12);
  // contact-width form of DNV-RP-C204 as arguments (B/D = 1): c = 23.2, n = 1.925/4.5
  const b = tubeDentResistance({ D: 1, t: 0.025, fy: 355e6, dent: 0.1, coefficient: 23.2, exponent: 1.925 / 4.5 });
  close(b.R, 350812 * 23.2 * 0.1 ** (1.925 / 4.5), 1e-5);
  assert.throws(() => tubeDentResistance({ D: 1, t: 0.025, fy: 355e6, dent: 1.2 }), RangeError);
});

test('three-hinge mechanism with membrane action: clamped rectangular beam P/Pc = 1 + (W/H)², 2W/H (Haythornthwaite; Jones, Structural Impact §7)', () => {
  const o = { section: 'rect', H: 0.1, B: 0.05, fy: 300e6, L: 2 };
  const Mp = 300e6 * 0.05 * 0.01 / 4;
  const r0 = beamMechanismResistance({ ...o, w: 0 });
  close(r0.R0, 8 * Mp / 2, 1e-12); close(r0.R, r0.R0, 1e-12);
  close(beamMechanismResistance({ ...o, w: 0.05 }).ratio, 1.25, 1e-12);
  close(beamMechanismResistance({ ...o, w: 0.1 }).ratio, 2, 1e-12);
  close(beamMechanismResistance({ ...o, w: 0.25 }).ratio, 5, 1e-12);
  // pinned ends: Pc = 4Mp/L, P/Pc = 1 + 4W²/H² (W ≤ H/2), 4W/H beyond
  const p = beamMechanismResistance({ ...o, ends: 'pinned', w: 0.025 });
  close(p.R0, 4 * Mp / 2, 1e-12); close(p.ratio, 1.25, 1e-12);
  close(beamMechanismResistance({ ...o, ends: 'pinned', w: 0.1 }).ratio, 4, 1e-12);
  // fully membrane state: R = 2 Np w/(L/2) for both end conditions (string solution)
  const Np = 300e6 * 0.05 * 0.1;
  close(beamMechanismResistance({ ...o, w: 0.3 }).R, 2 * Np * 0.3 / 1, 1e-12);
  close(beamMechanismResistance({ ...o, ends: 'pinned', w: 0.3 }).R, 2 * Np * 0.3 / 1, 1e-12);
  // energy: ∫(1 + w̄²) = w̄ + w̄³/3 up to 1, then w̄² − 1 beyond; clamped scale = H
  close(beamMechanismResistance({ ...o, w: 0.1 }).energy, r0.R0 * 0.1 * (4 / 3), 1e-9);
  close(beamMechanismResistance({ ...o, w: 0.2 }).energy, r0.R0 * 0.1 * (4 / 3 + 3), 1e-9);
});

test('three-hinge mechanism, tube: R/R0 = √(1 − w̄²) + w̄ asin w̄, (π/2) w̄ beyond w̄ = 1, w̄ = w/(D − t) for clamped ends', () => {
  const o = { section: 'tube', D: 0.6, t: 0.02, fy: 355e6, L: 12 };
  const Dm = 0.58; const Mp = 355e6 * Dm * Dm * 0.02; const Np = 355e6 * Math.PI * Dm * 0.02;
  const r = beamMechanismResistance({ ...o, w: 0.29 });
  close(r.R0, 8 * Mp / 12, 1e-12); close(r.wBar, 0.5, 1e-12);
  close(r.ratio, Math.sqrt(0.75) + 0.5 * Math.PI / 6, 1e-12); // 1.1278
  close(r.axialRatio, 1 / 3, 1e-12);
  // yield condition M/Mp = cos(π n/2) is satisfied by the solution: R = R0 [cos(πn/2) + n (π/2) w̄]
  close(r.ratio, Math.cos(Math.PI / 6) + (1 / 3) * (Math.PI / 2) * 0.5, 1e-12);
  // continuity and pure membrane beyond w̄ = 1: R = 2 Np w/(L/2)
  close(beamMechanismResistance({ ...o, w: Dm }).ratio, Math.PI / 2, 1e-12);
  close(beamMechanismResistance({ ...o, w: 1.0 }).R, 2 * Np * 1.0 / 6, 1e-12);
  // energy inverse
  const E = beamMechanismResistance({ ...o, w: 0.45 }).energy;
  close(beamMechanismDeflection({ ...o, energy: E }).w, 0.45, 1e-8);
  // with membrane action the deflection for a given energy is smaller than the bending-only estimate E/R0
  assert.ok(beamMechanismDeflection({ ...o, energy: E }).w < E / r.R0);
  assert.throws(() => beamMechanismResistance({ ...o, section: 'I', w: 0 }), RangeError);
});

test('energy sharing between ship and structure: linear springs split as the flexibilities; plateau takes the remainder; capacity', () => {
  // two linear springs in series: E_i = F²/(2k_i) → shares ∝ 1/k
  const ship = { d: [0, 1], F: [0, 10e6] }; const str = { d: [0, 0.5], F: [0, 20e6] };
  const r = sharedEnergy({ shipCurve: ship, structureCurve: str, energy: 1e6 });
  // F²(1/(2·10e6) + 1/(2·40e6)) = 1e6 → F = 4 MN; ship 0.8 MJ, structure 0.2 MJ
  close(r.force, 4e6, 1e-9); close(r.shipEnergy, 0.8e6, 1e-9); close(r.structureEnergy, 0.2e6, 1e-9);
  close(r.structureFraction, 0.2, 1e-9); close(r.shipDeformation, 0.4, 1e-9); close(r.structureDeformation, 0.1, 1e-9);
  assert.equal(r.exceeded, false);
  // rigid-plastic structure (plateau 3 MN over 1 m) against a ship stiffening linearly to 8 MN at 2 m:
  // the ship deforms until 3 MN (d = 0.75 m, 1.125 MJ), the structure then takes the rest at 3 MN
  const ship2 = { d: [0, 2], F: [0, 8e6] }; const str2 = { d: [0, 1], F: [3e6, 3e6] };
  const p = sharedEnergy({ shipCurve: ship2, structureCurve: str2, energy: 2e6 });
  close(p.force, 3e6, 1e-9); close(p.shipEnergy, 1.125e6, 1e-9); close(p.structureEnergy, 0.875e6, 1e-9);
  close(p.structureDeformation, 0.875 / 3, 1e-9); close(p.shipDeformation, 0.75, 1e-9);
  close(p.capacity, 1.125e6 + 3e6, 1e-9);
  // small energy: only the ship deforms (structure rigid below its plateau) — strength design
  const s = sharedEnergy({ shipCurve: ship2, structureCurve: str2, energy: 0.5e6 });
  close(s.shipEnergy, 0.5e6, 1e-9); assert.ok(s.structureEnergy < 1); close(s.force, Math.sqrt(2 * 0.5e6 * 4e6), 1e-9);
  assert.equal(classifyCollisionDesign({ structureEnergy: s.structureEnergy, shipEnergy: s.shipEnergy }).category, 'strength');
  assert.equal(classifyCollisionDesign({ structureEnergy: p.structureEnergy, shipEnergy: p.shipEnergy }).category, 'shared');
  assert.equal(classifyCollisionDesign({ structureEnergy: 9.5, shipEnergy: 0.5 }).category, 'ductile');
  // demand beyond the end of the structure curve
  const x = sharedEnergy({ shipCurve: ship2, structureCurve: str2, energy: 6e6 });
  assert.ok(x.exceeded); close(x.structureDeformation, 1, 1e-9);
  assert.throws(() => sharedEnergy({ shipCurve: { d: [0, 1], F: [5, 1] }, structureCurve: str, energy: 1 }), RangeError);
});

test('dropped object in air and in water: v = √(2gh); terminal velocity from weight–drag balance; ODE against v² = vT² + (v0² − vT²)e^(−ρCdA z/(m+ma))', () => {
  const a = dropInAir({ mass: 2000, height: 20 });
  close(a.velocity, Math.sqrt(2 * 9.80665 * 20), 1e-12); close(a.energy, 2000 * 9.80665 * 20, 1e-12);
  // steel container-like object: m = 2000 kg, V = 0.5 m³, Cd = 1.0, A = 1.5 m², Ca = 0.6
  const o = { mass: 2000, volume: 0.5, Cd: 1.0, area: 1.5, Ca: 0.6 };
  const tv = terminalVelocityWater(o);
  // W' = (2000 − 512.5)·9.80665 = 14 587 N; vT = √(2·14 587/(1025·1.0·1.5)) = 4.356 m/s
  close(tv.submergedWeight, 14587.4, 1e-5); close(tv.velocity, 4.356, 1e-3);
  close(tv.addedMass, 307.5, 1e-12); close(tv.energy, 0.5 * 2307.5 * tv.velocity ** 2, 1e-12);
  for (const [v0, depth] of [[0, 3], [a.velocity, 5], [a.velocity, 100], [0, 100]]) {
    const r = dropThroughWater({ ...o, depth, vEntry: v0 });
    const v2 = tv.velocity ** 2 + (v0 * v0 - tv.velocity ** 2) * Math.exp(-1025 * 1.0 * 1.5 * depth / 2307.5);
    assert.ok(r.reached);
    close(r.velocity, Math.sqrt(v2), 1e-6, `v0 = ${v0}, depth = ${depth}`);
    close(r.energy, 0.5 * 2307.5 * v2, 2e-6);
  }
  // deep water: terminal velocity regardless of the entry velocity
  close(dropThroughWater({ ...o, depth: 100, vEntry: a.velocity }).velocity, tv.velocity, 1e-6);
  // a buoyant object entering at speed stops: depth from the same closed form with vT² < 0
  const b = dropThroughWater({ mass: 300, volume: 0.5, Cd: 1.0, area: 1.5, Ca: 0.6, depth: 50, vEntry: 10 });
  assert.equal(b.reached, false); assert.equal(b.terminalVelocity, 0);
  const me = 300 + 307.5; const W = (512.5 - 300) * 9.80665; const c = 0.5 * 1025 * 1.5;
  // 0 = −W/c + (v0² + W/c) e^(−2cz/me) → z = me/(2c) ln(1 + c v0²/W)
  close(b.maxDepth, me / (2 * c) * Math.log(1 + c * 100 / W), 2e-3);
});

test('plate impact: yield-line collapse loads 2π mp R/(R − a), 4π mp R/(R − a); membrane 2π fy t δ/ln(R/a); energy of the envelope', () => {
  const t = 0.012; const fy = 355e6; const mp = fy * t * t / 4;
  const c = plateImpactCapacity({ t, fy, radius: 0.4, punchRadius: 0.05, deflection: 0 });
  close(c.Pc, 4 * Math.PI * mp * 0.4 / 0.35, 1e-12); close(c.mp, 12780, 1e-12);
  close(plateImpactCapacity({ t, fy, radius: 0.4, punchRadius: 0.05, edge: 'simple', deflection: 0 }).Pc, c.Pc / 2, 1e-12);
  // a → 0: classical point-load values 4π mp (clamped) and 2π mp (simply supported)
  close(plateImpactCapacity({ t, fy, radius: 0.4, punchRadius: 1e-9, membrane: false, deflection: 0 }).Pc, 4 * Math.PI * mp, 1e-8);
  // membrane stiffness and transition deflection
  const km = 2 * Math.PI * fy * t / Math.log(8);
  close(c.transitionDeflection, c.Pc / km, 1e-12);
  const d = plateImpactCapacity({ t, fy, radius: 0.4, punchRadius: 0.05, deflection: 0.1 });
  close(d.membraneLoad, km * 0.1, 1e-12); close(d.resistance, km * 0.1, 1e-12);
  const wt = c.transitionDeflection;
  close(d.energy, c.Pc * wt + 0.5 * km * (0.01 - wt * wt), 1e-12);
  // bending only: E = Pc δ
  close(plateImpactCapacity({ t, fy, radius: 0.4, punchRadius: 0.05, membrane: false, deflection: 0.1 }).energy, c.Pc * 0.1, 1e-12);
  assert.throws(() => plateImpactCapacity({ t, fy, radius: 0.4, punchRadius: 0.5, deflection: 0 }), RangeError);
});

test('DNV-RP-F107 pipe dent energy E = 16 (2π/9)^½ mp (D/t)^½ D (δ/D)^1.5: 12" pipe, 5 % dent, hand calculation', () => {
  // D = 323.9 mm, t = 12.7 mm, fy = 450 MPa: mp = 450e6·0.0127²/4 = 18 145 N; (D/t)^½ = 5.050; (0.05)^1.5 = 0.011180
  // E = 16·0.83554·18 145·5.050·0.3239·0.011180 = 4.436 kJ
  const r = pipeDentEnergy({ D: 0.3239, t: 0.0127, fy: 450e6, dent: 0.05 * 0.3239 });
  close(r.mp, 18145, 1e-4); close(r.energy, 4436, 2e-3);
  close(r.energy, 16 * Math.sqrt(2 * Math.PI / 9) * 18145.125 * Math.sqrt(0.3239 / 0.0127) * 0.3239 * 0.05 ** 1.5, 1e-12);
  close(pipeDentFromEnergy({ D: 0.3239, t: 0.0127, fy: 450e6, energy: r.energy }).dentRatio, 0.05, 1e-12);
  // energy scales with δ^1.5
  close(pipeDentEnergy({ D: 0.3239, t: 0.0127, fy: 450e6, dent: 0.1 * 0.3239 }).energy / r.energy, 2 ** 1.5, 1e-12);
});

test('vessel impact on a pile: elastic energy Mp² h/(6EI) and plastic hinge rotation, hand calculation', () => {
  // pile 1.2 m × 30 mm: I = 9.0e-3 m⁴ (approx), take EI = 1.9e9 N·m², Mp = 14.6 MN·m, h = 20 m
  const o = { Mp: 14.6e6, EI: 1.9e9, height: 20, rotationCapacity: 0.05 };
  const el = pileImpact({ ...o, energy: 100e3 });
  // Fc = 730 kN; δe = Fc h³/(3EI) = 1.0246 m; Ee = ½ Fc δe = 373.97 kJ
  close(el.collapseLoad, 730e3, 1e-12); close(el.elasticDeflection, 1.02456, 1e-5); close(el.elasticEnergy, 373.97e3, 1e-4);
  close(el.elasticEnergy, 14.6e6 ** 2 * 20 / (6 * 1.9e9), 1e-12);
  assert.ok(el.elastic); close(el.deflection, Math.sqrt(2 * 100e3 / (3 * 1.9e9 / 8000)), 1e-12);
  close(0.5 * el.force * el.deflection, 100e3, 1e-12);
  // 800 kJ: θp = (800 − 373.97)/14 600 = 0.02918 rad; capacity = 373.97 + 730 = 1103.97 kJ
  const pl = pileImpact({ ...o, energy: 800e3 });
  close(pl.plasticRotation, 0.029180, 1e-4); close(pl.capacity, 1103.97e3, 1e-5); close(pl.uc, 800 / 1103.97, 1e-5);
  close(pl.deflection, 1.02456 + 0.02918 * 20, 1e-4); assert.equal(pl.force, 730e3);
});
