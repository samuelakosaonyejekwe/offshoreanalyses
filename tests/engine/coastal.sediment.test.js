import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  dimensionlessGrainSize, settlingVelocitySoulsby, settlingVelocityVanRijn, criticalShields,
  criticalShearStress, shieldsParameter, mobilityCheck, grainRoughness, bedformRoughness, currentRipples,
  bedShearCurrent, waveOrbital, waveFrictionFactor, bedShearWave, combinedShearSoulsby, meyerPeterMuller,
  engelundHansen, vanRijnBedLoad, vanRijnSuspendedLoad, vanRijnSimplified, soulsbyVanRijn, rouseProfile,
  isbashStoneSize, shieldsRockSize,
} from '../../src/engine/coastal/sediment.js';

const g = 9.80665;
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const abs = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b}`);

// Sea water at 10 °C, 35 ppt, quartz sand — the reference conditions of the worked examples in Soulsby (1997).
const sea = { rhoS: 2650, rho: 1027, nu: 1.36e-6, g: 9.81 };

test('Soulsby (1997) worked values for 0.2 mm sand in sea water: D* = 4.06, ws = 0.0202 m/s, θcr = 0.0553, τcr = 0.176 Pa', () => {
  const Ds = dimensionlessGrainSize({ d: 0.2e-3, ...sea });
  rel(Ds, 0.2e-3 * ((9.81 * (2650 / 1027 - 1)) / 1.36e-6 ** 2) ** (1 / 3), 1e-12);
  abs(Ds, 4.06, 0.005);
  abs(settlingVelocitySoulsby({ d: 0.2e-3, ...sea }), 0.0202, 1e-4);
  const c = criticalShearStress({ d: 0.2e-3, ...sea });
  abs(c.thetaCr, 0.0553, 1e-4);
  abs(c.tauCr, 0.176, 1e-3);
  rel(c.uStarCr, Math.sqrt(c.tauCr / 1027), 1e-12);
});

test('Soulsby–Whitehouse threshold: formula, coarse-grain asymptote 0.055 and fine-grain limit 0.30', () => {
  rel(criticalShields(10), 0.3 / 13 + 0.055 * (1 - Math.exp(-0.2)), 1e-12);
  abs(criticalShields(1e4), 0.055, 1e-4);
  abs(criticalShields(1e-6), 0.3, 1e-5);
  assert.throws(() => criticalShields(0), RangeError);
});

test('settling velocity: Stokes limit, van Rijn ranges and coarse-grain scaling', () => {
  const o = { rhoS: 2650, rho: 1000, nu: 1e-6 };
  const stokes = (d) => (1.65 * g * d * d) / (18 * 1e-6);
  // Fine silt: both formulae tend to Stokes (Soulsby's constant gives 1.049/(2·10.36) = 1/19.75 instead of 1/18).
  rel(settlingVelocityVanRijn({ d: 50e-6, ...o }), stokes(50e-6), 1e-12);
  rel(settlingVelocitySoulsby({ d: 5e-6, ...o }), stokes(5e-6) * (18 / 19.75), 2e-3);
  // van Rijn intermediate range, d = 0.3 mm, by hand
  const d = 0.3e-3;
  rel(settlingVelocityVanRijn({ d, ...o }), ((10 * 1e-6) / d) * (Math.sqrt(1 + (0.01 * 1.65 * g * d ** 3) / 1e-12) - 1), 1e-12);
  rel(settlingVelocityVanRijn({ d: 2e-3, ...o }), 1.1 * Math.sqrt(1.65 * g * 2e-3), 1e-12);
  // The two formulae agree within 15 % for sand.
  for (const dd of [0.15e-3, 0.3e-3, 0.6e-3, 0.9e-3]) {
    rel(settlingVelocityVanRijn({ d: dd, ...o }), settlingVelocitySoulsby({ d: dd, ...o }), 0.15, `d=${dd}`);
  }
});

test('current bed shear stress: depth-averaged log law verified by integrating the velocity profile', () => {
  const h = 5;
  const z0 = 0.001;
  const U = 1.2;
  const r = bedShearCurrent({ U, h, z0, rho: 1027 });
  rel(r.Cd, (0.4 / (Math.log(h / z0) - 1)) ** 2, 1e-12);
  // u(z) = (u*/κ) ln(z/z0); its depth average must return U (to O(z0/h)).
  const n = 200000;
  let s = 0;
  const lz0 = Math.log(z0);
  const dl = (Math.log(h) - lz0) / n;
  for (let i = 0; i < n; i++) {
    const z = Math.exp(lz0 + (i + 0.5) * dl);
    s += (r.uStar / 0.4) * Math.log(z / z0) * z * dl;
  }
  rel(s / h, U, 5e-4);
  rel(r.tau, 1027 * r.uStar ** 2, 1e-12);
  // Manning: τ = ρ g n² U²/h^(1/3); Colebrook–White: C = 18 log10(12h/ks)
  rel(bedShearCurrent({ U, h, method: 'manning', n: 0.025, rho: 1000 }).tau, (1000 * g * 0.025 ** 2 * U * U) / h ** (1 / 3), 1e-12);
  const cw = bedShearCurrent({ U, h, method: 'colebrook', ks: 0.05, rho: 1000 });
  rel(cw.chezy, 18 * Math.log10(1200), 1e-12);
  rel(cw.tau, (1000 * g * U * U) / cw.chezy ** 2, 1e-12);
  // The log law and Colebrook–White describe the same rough-wall profile: with ks = 30 z0 they agree within 5 %.
  rel(bedShearCurrent({ U, h, method: 'colebrook', ks: 30 * z0, rho: 1027 }).tau, r.tau, 0.05);
  assert.throws(() => bedShearCurrent({ U, h: 0.002, z0 }), RangeError);
});

test('wave orbital velocity and friction factors: linear theory, Soulsby and Swart', () => {
  const o = waveOrbital({ H: 2, T: 8, h: 10 });
  rel((2 * Math.PI / 8) ** 2, g * o.k * Math.tanh(o.k * 10), 1e-12);
  rel(o.Uw, (Math.PI * 2) / (8 * Math.sinh(o.k * 10)), 1e-12);
  rel(o.A, (o.Uw * 8) / (2 * Math.PI), 1e-12);
  // Shallow-water limit: Uw = (H/2) √(g/h)
  rel(waveOrbital({ H: 0.2, T: 60, h: 2 }).Uw, 0.1 * Math.sqrt(g / 2), 5e-3);
  rel(waveFrictionFactor({ A: 1, z0: 0.0005 }), 1.39 * 2000 ** -0.52, 1e-12);
  rel(waveFrictionFactor({ A: 1, ks: 0.015, method: 'swart' }), Math.exp(5.213 * (1 / 0.015) ** -0.194 - 5.977), 1e-12);
  assert.equal(waveFrictionFactor({ A: 0.01, ks: 0.01, method: 'swart' }), 0.3);
  // Swart's curve is continuous at the cap A/ks = 1.57 (fw = 0.30).
  abs(Math.exp(5.213 * 1.57 ** -0.194 - 5.977), 0.3, 2e-3);
  const w = bedShearWave({ Uw: o.Uw, T: 8, z0: 0.0005, rho: 1027 });
  rel(w.tau, 0.5 * 1027 * w.fw * o.Uw ** 2, 1e-12);
  assert.equal(bedShearWave({ Uw: 0, T: 8, z0: 0.0005 }).tau, 0);
});

test('Soulsby DATA2 wave–current stresses: limits and hand calculation', () => {
  const c = combinedShearSoulsby({ tauC: 0.5, tauW: 0 });
  rel(c.tauMean, 0.5, 1e-12);
  rel(c.tauMax, 0.5, 1e-12);
  const w = combinedShearSoulsby({ tauC: 0, tauW: 2 });
  assert.equal(w.tauMean, 0);
  rel(w.tauMax, 2, 1e-12);
  // τc = 0.5, τw = 1.5, φ = 30°: τm = 0.5 (1 + 1.2·0.75^3.2) = 0.739
  const r = combinedShearSoulsby({ tauC: 0.5, tauW: 1.5, phiDeg: 30 });
  const tm = 0.5 * (1 + 1.2 * 0.75 ** 3.2);
  rel(r.tauMean, tm, 1e-12);
  abs(r.tauMean, 0.739, 1e-3);
  rel(r.tauMax, Math.hypot(tm + 1.5 * Math.cos(Math.PI / 6), 1.5 * Math.sin(Math.PI / 6)), 1e-12);
  // Collinear: τmax = τm + τw; perpendicular: √(τm² + τw²)
  rel(combinedShearSoulsby({ tauC: 0.5, tauW: 1.5 }).tauMax, tm + 1.5, 1e-12);
  rel(combinedShearSoulsby({ tauC: 0.5, tauW: 1.5, phiDeg: 90 }).tauMax, Math.hypot(tm, 1.5), 1e-12);
  assert.ok(r.tauMean > 0.5); // waves enhance the mean stress
});

test('mobility check and roughness helpers', () => {
  const c = criticalShearStress({ d: 0.2e-3, ...sea });
  assert.equal(mobilityCheck({ tau: 0.9 * c.tauCr, d: 0.2e-3, ...sea }).mobile, false);
  const m = mobilityCheck({ tau: 2 * c.tauCr, d: 0.2e-3, ...sea });
  assert.equal(m.mobile, true);
  rel(m.ratio, 2, 1e-12);
  rel(shieldsParameter({ tau: 1, d: 1e-3, rhoS: 2650, rho: 1000 }), 1 / (g * 1650 * 1e-3), 1e-12);
  rel(grainRoughness(0.3e-3).ks, 0.75e-3, 1e-12);
  rel(grainRoughness(0.3e-3).z0, 0.3e-3 / 12, 1e-12);
  const rp = currentRipples(0.2e-3);
  rel(rp.length, 0.2, 1e-12);
  rel(rp.height, 0.2 / 7, 1e-12);
  rel(bedformRoughness({ height: rp.height, length: rp.length }).z0, rp.height ** 2 / rp.length, 1e-12);
});

test('Meyer-Peter & Müller: Φ = 8 (θ − 0.047)^1.5 by hand; zero below threshold', () => {
  assert.equal(meyerPeterMuller({ theta: 0.04, d: 2e-3 }).qb, 0);
  // θ = 0.147: Φ = 8·0.1^1.5 = 0.2530; d = 2 mm in fresh water: √(g·1.65·8e-9) = 3.598e-4 → qb = 9.10e-5 m²/s
  const r = meyerPeterMuller({ theta: 0.147, d: 2e-3, rhoS: 2650, rho: 1000 });
  rel(r.phi, 8 * 0.1 ** 1.5, 1e-10);
  rel(r.qb, r.phi * Math.sqrt(g * 1.65 * 8e-9), 1e-12);
  abs(r.qb, 9.1e-5, 1e-6);
});

test('Engelund–Hansen: the f·Φ = 0.1 θ^2.5 form equals 0.05 U⁵ /((s−1)² √g d50 C³)', () => {
  const U = 1.4;
  const h = 4;
  const d50 = 0.4e-3;
  const cw = bedShearCurrent({ U, h, method: 'colebrook', ks: 0.1, rho: 1000 });
  const r = engelundHansen({ U, tau: cw.tau, d50, rhoS: 2650, rho: 1000 });
  rel(r.qt, (0.05 * U ** 5) / (1.65 ** 2 * Math.sqrt(g) * d50 * cw.chezy ** 3), 1e-10);
  rel(r.theta, cw.tau / (1000 * g * 1.65 * d50), 1e-12);
  // qt ∝ U⁵ at fixed Chézy coefficient
  const r2 = engelundHansen({ U: 2 * U, tau: 4 * cw.tau, d50, rhoS: 2650, rho: 1000 });
  rel(r2.qt / r.qt, 32, 1e-10);
  assert.equal(engelundHansen({ U: 0, tau: 0, d50 }).qt, 0);
});

test('van Rijn (1984) bed load: transport stage and rate by hand', () => {
  const o = { U: 1.0, h: 3, d50: 0.5e-3, d90: 1e-3, rhoS: 2650, rho: 1000, nu: 1e-6 };
  const r = vanRijnBedLoad(o);
  const Dstar = 0.5e-3 * ((g * 1.65) / 1e-12) ** (1 / 3);
  rel(r.Dstar, Dstar, 1e-12);
  const C = 18 * Math.log10((12 * 3) / (3 * 1e-3));
  const us = (Math.sqrt(g) / C) * 1.0;
  rel(r.uStarGrain, us, 1e-12);
  const ucr2 = criticalShields(Dstar) * g * 1.65 * 0.5e-3;
  const T = (us * us - ucr2) / ucr2;
  rel(r.T, T, 1e-12);
  rel(r.qb, (0.053 * Math.sqrt(1.65 * g) * (0.5e-3) ** 1.5 * T ** 2.1) / Dstar ** 0.3, 1e-12);
  assert.ok(r.T > 1 && r.T < 15);
  assert.equal(vanRijnBedLoad({ ...o, U: 0.2 }).qb, 0);
});

test('van Rijn (1984) suspended load: reference concentration, suspension number and F factor by hand', () => {
  const o = { U: 1.5, h: 5, d50: 0.25e-3, d90: 0.5e-3, ks: 0.05, sigmaS: 2, rhoS: 2650, rho: 1000, nu: 1e-6 };
  const r = vanRijnSuspendedLoad(o);
  const Dstar = 0.25e-3 * ((g * 1.65) / 1e-12) ** (1 / 3);
  const a = Math.max(0.05, 0.01 * 5);
  rel(r.a, a, 1e-12);
  rel(r.ca, (0.015 * (0.25e-3 / a) * r.T ** 1.5) / Dstar ** 0.3, 1e-12);
  rel(r.uStar, (Math.sqrt(g) / (18 * Math.log10(1200))) * 1.5, 1e-12);
  const ds = r.T < 25 ? 0.25e-3 * (1 + 0.011 * (r.T - 25)) : 0.25e-3;
  rel(r.ds, ds, 1e-12);
  rel(r.ws, settlingVelocityVanRijn({ d: ds, rhoS: 2650, rho: 1000, nu: 1e-6 }), 1e-12);
  const beta = Math.min(2, 1 + 2 * (r.ws / r.uStar) ** 2);
  rel(r.Z, r.ws / (beta * 0.4 * r.uStar), 1e-12);
  const Zp = r.Z + 2.5 * (r.ws / r.uStar) ** 0.8 * (r.ca / 0.65) ** 0.4;
  rel(r.Zprime, Zp, 1e-12);
  const ah = a / 5;
  rel(r.F, (ah ** Zp - ah ** 1.2) / ((1 - ah) ** Zp * (1.2 - Zp)), 1e-12);
  rel(r.qs, r.F * 1.5 * 5 * r.ca, 1e-12);
  assert.ok(r.F > 0 && r.F < 1 && r.qs > 0);
  // Stronger flow suspends more; no transport below threshold.
  assert.ok(vanRijnSuspendedLoad({ ...o, U: 2 }).qs > r.qs);
  assert.equal(vanRijnSuspendedLoad({ ...o, U: 0.15 }).qs, 0);
});

test('van Rijn simplified formulae and Soulsby–van Rijn: identical for a current alone; waves enhance; slope term', () => {
  const o = { U: 1.2, h: 6, d50: 0.3e-3, d90: 0.6e-3, rhoS: 2650, rho: 1027, nu: 1.36e-6 };
  const v = vanRijnSimplified(o);
  rel(v.Ucr, 0.19 * (0.3e-3) ** 0.1 * Math.log10((4 * 6) / 0.6e-3), 1e-12);
  const s1 = 2650 / 1027 - 1;
  const m = ((1.2 - v.Ucr) / Math.sqrt(s1 * g * 0.3e-3)) ** 2.4;
  rel(v.qb, 0.005 * 1.2 * 6 * m * (0.3e-3 / 6) ** 1.2, 1e-12);
  const Dstar = dimensionlessGrainSize({ d: 0.3e-3, rhoS: 2650, rho: 1027, nu: 1.36e-6 });
  rel(v.qs, 0.012 * 1.2 * 6 * m * (0.3e-3 / 6) * Dstar ** -0.6, 1e-12);
  rel(v.qt, v.qb + v.qs, 1e-12);
  // Coarse branch of the threshold velocity
  rel(vanRijnSimplified({ ...o, d50: 1e-3, d90: 2e-3 }).Ucr, 8.5 * (1e-3) ** 0.6 * Math.log10((4 * 6) / 2e-3), 1e-12);
  const c = soulsbyVanRijn({ ...o, Urms: 0 });
  rel(c.qb, v.qb, 1e-10);
  rel(c.qs, v.qs, 1e-10);
  rel(c.Cd, (0.4 / (Math.log(6 / 0.006) - 1)) ** 2, 1e-12);
  const w = soulsbyVanRijn({ ...o, Urms: 0.4 });
  const Ue = Math.sqrt(1.2 ** 2 + (0.018 / c.Cd) * 0.16);
  rel(w.qt / c.qt, ((Ue - v.Ucr) / (1.2 - v.Ucr)) ** 2.4, 1e-10);
  rel(soulsbyVanRijn({ ...o, Urms: 0.4, bedSlope: 0.05 }).qt, w.qt * (1 - 0.08), 1e-12);
  // Waves can mobilise sediment under a current that is itself below threshold.
  assert.equal(vanRijnSimplified({ ...o, U: 0.3 }).qt, 0);
  assert.ok(soulsbyVanRijn({ ...o, U: 0.3, Urms: 0.5 }).qt > 0);
  assert.throws(() => vanRijnSimplified({ ...o, d50: 5e-3 }), RangeError);
});

test('Rouse profile: reference value, surface limit, power-law form and the governing diffusion balance', () => {
  const o = { h: 4, a: 0.05, ca: 1.0, ws: 0.02, uStar: 0.06 };
  rel(rouseProfile({ ...o, z: 0.05 }).c, 1, 1e-12);
  assert.equal(rouseProfile({ ...o, z: 4 }).c, 0);
  rel(rouseProfile({ ...o, z: 1 }).rouse, 0.02 / (0.4 * 0.06), 1e-12);
  // ws c + εs dc/dz = 0 with the parabolic diffusivity εs = κ u* z (1 − z/h)
  for (const z of [0.2, 1, 2.5, 3.6]) {
    const dz = 1e-5;
    const dc = (rouseProfile({ ...o, z: z + dz }).c - rouseProfile({ ...o, z: z - dz }).c) / (2 * dz);
    const eps = 0.4 * 0.06 * z * (1 - z / 4);
    abs(0.02 * rouseProfile({ ...o, z }).c + eps * dc, 0, 1e-8);
  }
  assert.throws(() => rouseProfile({ ...o, z: 0.01 }), RangeError);
});

test('Isbash and Shields-based rock sizing against currents', () => {
  // U = 3 m/s, Δ = 1.65, C = 0.86: D = 9/(2 g·1.65·0.7396) = 0.376 m
  const D = isbashStoneSize({ U: 3, delta: 1.65 });
  rel(D, 9 / (2 * g * 1.65 * 0.86 ** 2), 1e-12);
  abs(D, 0.376, 1e-3);
  rel(isbashStoneSize({ U: 3, delta: 1.65, C: 1.2 }) / D, (0.86 / 1.2) ** 2, 1e-12);
  // The Shields-based size satisfies θ = U²/(C² Δ d) = ψcr with C = 18 log10(12h/(2d)).
  const r = shieldsRockSize({ U: 3, h: 6, delta: 1.65, psiCr: 0.035 });
  rel(r.chezy, 18 * Math.log10((12 * 6) / (2 * r.d)), 1e-10);
  rel((3 * 3) / (r.chezy ** 2 * 1.65 * r.d), 0.035, 1e-8);
  // Deeper water at the same mean velocity means less shear and smaller rock; side slopes need larger rock.
  assert.ok(shieldsRockSize({ U: 3, h: 12, delta: 1.65 }).d < r.d);
  const sl = shieldsRockSize({ U: 3, h: 6, delta: 1.65, slopeDeg: (Math.atan(0.5) * 180) / Math.PI, reposeDeg: 40 });
  rel(sl.ksl, Math.sqrt(1 - (Math.sin(Math.atan(0.5)) / Math.sin((40 * Math.PI) / 180)) ** 2), 1e-12);
  assert.ok(sl.d > r.d);
  assert.throws(() => shieldsRockSize({ U: 8, h: 0.3, delta: 1.65 }), RangeError);
});
