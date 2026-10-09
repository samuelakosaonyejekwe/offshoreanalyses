import test from 'node:test';
import assert from 'node:assert/strict';
import * as mo from '../../src/engine/hydro/morison.js';
import { G, dispersion } from '../../src/engine/hydro/_wavehelp.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);

test('Morison per unit length: flow along the axis gives no force; normal flow gives ½ rho Cd D u|u| + rho Cm A a', () => {
  const along = mo.morisonForcePerLength({ D: 1, Cd: 1, Cm: 2, axis: [0, 0, 2], u: [0, 0, 3], a: [0, 0, 1] });
  assert.deepEqual(along.f, [0, 0, 0]);
  const r = mo.morisonForcePerLength({ rho: 1000, D: 2, Cd: 0.7, Cm: 2, axis: [0, 0, 1], u: [1.5, 0, 9], a: [0.4, 0, 5] });
  close(r.f[0], 0.5 * 1000 * 0.7 * 2 * 1.5 * 1.5 + 1000 * 2 * Math.PI * 0.4, 1e-12);
  // 45° inclined member in the x–z plane: only the normal component (u/sqrt2) loads the member
  const inc = mo.morisonForcePerLength({ rho: 1000, D: 1, Cd: 1, Cm: 0, axis: [1, 0, 1], u: [2, 0, 0], a: [0, 0, 0] });
  close(Math.hypot(...inc.f), 0.5 * 1000 * 1 * (2 / Math.SQRT2) ** 2, 1e-12);
});

test('Relative-velocity form: a member moving with the fluid feels only the Froude–Krylov part rho A a', () => {
  const r = mo.morisonForcePerLength({ rho: 1000, D: 1, Cd: 1, Cm: 2, axis: [0, 0, 1], u: [1, 0, 0], a: [0.5, 0, 0], vBody: [1, 0, 0], aBody: [0.5, 0, 0] });
  close(r.f[0], 1000 * (Math.PI / 4) * 0.5, 1e-12);
});

test('Vertical surface-piercing pile vs closed-form Morison integral (hand calculation of F_I and F_D)', () => {
  const H = 6; const T = 10; const d = 30; const D = 1.2; const Cd = 1.0; const Cm = 2.0; const rho = 1025;
  const w = (2 * Math.PI) / T; const k = dispersion(w, d);
  // hand formulas, written out independently of the library
  const A = (Math.PI * D * D) / 4;
  const FI = rho * Cm * A * (H / 2) * w * w * (1 / k);                       // ∫ cosh k(z+d)/sinh kd dz = 1/k
  const FD = 0.5 * rho * Cd * D * ((w * H) / 2) ** 2 * (Math.sinh(2 * k * d) / (4 * k) + d / 2) / Math.sinh(k * d) ** 2;
  const cf = mo.pileForceClosedForm({ H, T, depth: d, D, Cd, Cm, rho });
  close(cf.FI, FI, 1e-10); close(cf.FD, FD, 1e-10);
  const base = { p1: [0, 0, -d], p2: [0, 0, 10], D, rho, wave: { H, T, depth: d }, about: [0, 0, -d], nGauss: 12, nSeg: 4 };
  // drag only at the crest (t = 0); inertia only a quarter period later (|a| maximum)
  const drag = mo.memberWaveForce({ ...base, Cd, Cm: 0, t: 0 });
  close(drag.F[0], FD, 1e-9, 'drag force'); close(drag.M[1], cf.MD, 1e-9, 'drag moment');
  const iner = mo.memberWaveForce({ ...base, Cd: 0, Cm, t: T / 4 });
  close(Math.abs(iner.F[0]), FI, 1e-9, 'inertia force'); close(Math.abs(iner.M[1]), cf.MI, 1e-9, 'inertia moment');
  // maximum of the combined force over the cycle
  const hist = mo.memberForceHistory({ ...base, Cd, Cm, nSteps: 720 });
  close(hist.FhMax, cf.Fmax, 2e-4, 'max base shear'); close(hist.MhMax, cf.Mmax, 2e-4, 'max overturning moment');
});

test('Drag-dominated pile: F_max = F_D + F_I²/(4 F_D)', () => {
  const cf = mo.pileForceClosedForm({ H: 12, T: 12, depth: 20, D: 0.3, Cd: 1.05, Cm: 1.2 });
  assert.ok(cf.FI < 2 * cf.FD);
  close(cf.Fmax, cf.FD + (cf.FI * cf.FI) / (4 * cf.FD), 1e-12);
  const h = mo.memberForceHistory({ p1: [0, 0, -20], p2: [0, 0, 5], D: 0.3, Cd: 1.05, Cm: 1.2, wave: { H: 12, T: 12, depth: 20 }, nSteps: 1440 });
  close(h.FhMax, cf.Fmax, 1e-4);
});

test('Stretching: Wheeler and vertical extrapolation add load above still water under the crest', () => {
  const base = { p1: [0, 0, -30], p2: [0, 0, 15], D: 1, Cd: 1, Cm: 0, wave: { H: 8, T: 10, depth: 30 }, t: 0 };
  const none = mo.memberWaveForce(base).F[0];
  const vert = mo.memberWaveForce({ ...base, stretching: 'vertical' }).F[0];
  const whe = mo.memberWaveForce({ ...base, stretching: 'wheeler' }).F[0];
  // vertical extrapolation adds exactly eta × (force per unit length at z = 0)
  const w = (2 * Math.PI) / 10; const k = dispersion(w, 30);
  const f0 = 0.5 * 1025 * (4 * w / Math.tanh(k * 30)) ** 2;
  close(vert - none, 4 * f0, 1e-6);
  assert.ok(whe > none && whe < vert);
});

test('Horizontal member normal to deep-water waves: |F| = rho Cm A L a omega² exp(kz), constant in time', () => {
  const T = 8; const w = (2 * Math.PI) / T; const k = w * w / G; const z = -12;
  const exp = 1025 * 2 * (Math.PI / 4) * 20 * 1.5 * w * w * Math.exp(k * z);
  for (const t of [0, 1.3, 2.9, 5]) {
    const r = mo.memberWaveForce({ p1: [0, -10, z], p2: [0, 10, z], D: 1, Cd: 0, Cm: 2, wave: { H: 3, T }, t });
    close(Math.hypot(r.F[0], r.F[2]), exp, 1e-9);
    assert.ok(Math.abs(r.F[1]) < 1e-9 * exp);
  }
});

test('DNV-RP-C205 §6.7: C_DS(roughness) = 0.65 / (29 + 4 log10 Δ)/20 / 1.05, continuous at the limits', () => {
  assert.equal(mo.dragCoefficientSteady(0), 0.65);
  assert.equal(mo.dragCoefficientSteady(5e-2), 1.05);
  close(mo.dragCoefficientSteady(1e-3), (29 - 12) / 20, 1e-12);
  close(mo.dragCoefficientSteady(1e-4), 0.65, 1e-12);
  close(mo.dragCoefficientSteady(1e-2), 1.05, 1e-12);
});

test('DNV-RP-C205 §6.7.2 wake amplification and §6.9.1 inertia coefficient (hand values)', () => {
  const cpi = 1.5 - 0.024 * (12 / 1.05 - 10);
  close(cpi, 1.4657142857, 1e-9);
  close(mo.wakeAmplification(6, 1.05), cpi - 0.6, 1e-12);
  close(mo.wakeAmplification(1, 1.05), cpi - 1, 1e-12);
  close(mo.wakeAmplification(0.5, 1.05), cpi - 1 + 0.5, 1e-12);
  assert.throws(() => mo.wakeAmplification(15, 1.05), RangeError);
  assert.equal(mo.inertiaCoefficient(2, 0.65).Cm, 2);
  close(mo.inertiaCoefficient(10, 0.65).Ca, 1 - 0.044 * 7, 1e-12);
  close(mo.inertiaCoefficient(60, 0.65).Ca, 0.6, 1e-12);
  close(mo.inertiaCoefficient(60, 1.05).Ca, 0.2, 1e-12);
});

test('MacCamy–Fuchs tends to the Morison inertia force with C_M = 2 as ka -> 0', () => {
  const H = 1; const T = 15; const d = 40; const D = 0.5;
  const mf = mo.maccamyFuchs({ H, T, depth: d, D });
  assert.ok(mf.ka < 0.01);
  close(mf.CmEff, 2, 1e-3);
  const cf = mo.pileForceClosedForm({ H, T, depth: d, D, Cd: 0, Cm: 2 });
  close(mf.F, cf.FI, 1e-3); close(mf.M, cf.MI, 1e-3);
  // diffraction reduces the effective inertia coefficient for large ka
  const big = mo.maccamyFuchs({ H, T: 6, depth: d, D: 30 });
  assert.ok(big.CmEff < 1.2 && big.ka > 1.5);
  // distributed force integrates to the total
  let s = 0; const n = 2000;
  for (let i = 0; i < n; i++) s += big.forcePerLength(-d + ((i + 0.5) * d) / n) * (d / n);
  close(s, big.F, 1e-5);
});

test('Added mass: circle rho pi a², square 1.51, plate 1.0, sphere ½ displaced mass, disc 8/3 rho a³', () => {
  close(mo.addedMass2D({ shape: 'circle', a: 2, rho: 1000 }).addedMassPerLength, 1000 * Math.PI * 4, 1e-12);
  assert.equal(mo.addedMass2D({ shape: 'square', a: 1 }).Ca, 1.51);
  assert.equal(mo.addedMass2D({ shape: 'plate', a: 1 }).Ca, 1.0);
  close(mo.addedMass2D({ shape: 'rectangle', a: 1, b: 1 }).Ca, 1.51, 1e-12);
  close(mo.addedMass2D({ shape: 'rectangle', a: 1, b: 0.5 }).Ca, 1.36, 1e-12);
  close(mo.addedMass2D({ shape: 'rectangle', a: 1, b: 10 }).Ca, 2.23, 1e-12);
  close(mo.addedMass3D({ shape: 'sphere', a: 3, rho: 1000 }).addedMass, (2 / 3) * Math.PI * 27 * 1000, 1e-12);
  close(mo.addedMass3D({ shape: 'disc', a: 3, rho: 1000 }).addedMass, (8 / 3) * 1000 * 27, 1e-12);
});

test('Member buoyancy, marine growth, slam force and shedding frequency (hand calculations)', () => {
  const b = mo.memberBuoyancy({ p1: [0, 0, -10], p2: [0, 0, 10], D: 2, rho: 1000, g: 10 });
  close(b.volume, Math.PI * 10, 1e-12); close(b.force, 1e5 * Math.PI, 1e-12); assert.deepEqual(b.centre, [0, 0, -5]);
  const fl = mo.memberBuoyancy({ p1: [0, 0, -10], p2: [5, 0, -10], D: 1, t: 0.02, flooded: true });
  close(fl.volume, (Math.PI / 4) * (1 - 0.96 ** 2) * 5, 1e-12);
  const mg = mo.marineGrowth({ D: 1, thickness: 0.1 });
  close(mg.De, 1.2, 1e-12); close(mg.massPerLength, 1325 * (Math.PI / 4) * (1.44 - 1), 1e-12);
  close(mo.slamForcePerLength({ D: 0.8, v: 5, Cs: Math.PI, rho: 1025 }), 0.5 * 1025 * Math.PI * 0.8 * 25, 1e-12);
  close(mo.vortexSheddingFrequency({ U: 2, D: 0.5 }), 0.8, 1e-12);
  close(mo.keuleganCarpenter(2, 10, 1), 20, 1e-12);
});
