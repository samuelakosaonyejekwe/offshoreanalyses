import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  beamTransformationFactors, equivalentSdof, loadHistory, sdofBlastResponse, strainEnergyEpp, pressureImpulseCurve, supportRotation,
  beamBlastResponse, dynamicIncreaseFactor, meanStrainRate, blastDragLoad,
} from '../../src/engine/structural/blast.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b} (abs ${abs})`);

test('transformation factors derived from the deflected shapes: exact integrals and Biggs Tables 5.1–5.3 (two decimals)', () => {
  const f = (support, load, range) => beamTransformationFactors({ support, load, range });
  // simple, uniform, elastic: KL = 16/25, KM = 3968/7875 (∫[16/5 (ξ − 2ξ³ + ξ⁴)]² dξ), V = 0.39R + 0.11F
  let r = f('simple', 'uniform', 'elastic');
  close(r.KL, 0.64, 1e-12); close(r.KM, 3968 / 7875, 1e-12); close(r.KLM, 0.7873, 1e-4);
  close(r.inertiaCentroid, 61 / 192, 1e-12);
  near(r.reaction.a, 0.39, 0.005); near(r.reaction.b, 0.11, 0.005);
  // plastic: KL = 1/2, KM = 1/3, V = 0.38R + 0.12F
  r = f('simple', 'uniform', 'plastic');
  close(r.KL, 0.5, 1e-12); close(r.KM, 1 / 3, 1e-12); close(r.KLM, 2 / 3, 1e-12);
  close(r.reaction.a, 0.375, 1e-12); close(r.reaction.b, 0.125, 1e-12);
  // simple, point load at midspan: KL = 1, KM = 17/35 = 0.49, V = 0.78R − 0.28F; plastic 0.33, 0.75R − 0.25F
  r = f('simple', 'point', 'elastic');
  assert.equal(r.KL, 1); close(r.KM, 17 / 35, 1e-12); close(r.reaction.a, 0.78125, 1e-12); close(r.reaction.b, -0.28125, 1e-12);
  r = f('simple', 'point', 'plastic');
  close(r.KM, 1 / 3, 1e-12); close(r.reaction.a, 0.75, 1e-12); close(r.reaction.b, -0.25, 1e-12);
  // fixed, uniform, elastic: KL = 8/15 = 0.53, KM = 256·B(5,5) = 128/315 = 0.41, KLM = 0.77, V = 0.36R + 0.14F
  r = f('fixed', 'uniform', 'elastic');
  close(r.KL, 8 / 15, 1e-12); close(r.KM, 128 / 315, 1e-12); near(r.KLM, 0.77, 0.01);
  near(r.reaction.a, 0.36, 0.005); near(r.reaction.b, 0.14, 0.005);
  // fixed, point: KM = 13/35 = 0.37, V = 0.71R − 0.21F
  r = f('fixed', 'point', 'elastic');
  close(r.KM, 13 / 35, 1e-12); near(r.reaction.a, 0.71, 0.005); near(r.reaction.b, -0.21, 0.005);
  // cantilever, uniform: KL = ∫(ξ⁴ − 4ξ³ + 6ξ²)/3 dξ = 0.4, KM = 104/405 = 0.2568
  r = f('cantilever', 'uniform', 'elastic');
  close(r.KL, 0.4, 1e-12); close(r.KM, 104 / 405, 1e-12);
  // statics is recovered when F = R: a + b = ½ per support (1 for the single cantilever support)
  for (const s of ['simple', 'fixed']) for (const l of ['uniform', 'point']) for (const g of ['elastic', 'plastic']) close(f(s, l, g).reaction.a + f(s, l, g).reaction.b, 0.5, 1e-12);
  for (const l of ['uniform', 'point']) for (const g of ['elastic', 'plastic']) close(f('cantilever', l, g).reaction.a + f('cantilever', l, g).reaction.b, 1, 1e-12);
  // shapes are normalised at the reference point
  close(f('fixed', 'point', 'elastic').shape(0.5), 1, 1e-12); close(f('cantilever', 'point', 'elastic').shape(1), 1, 1e-12);
  assert.throws(() => f('propped', 'uniform', 'elastic'), RangeError);
});

test('equivalent SDOF: stiffness, resistance, period vs the exact beam frequency, fixed-beam bilinear k_E = 307 EI/L³ (Biggs Table 5.2)', () => {
  const L = 6; const EI = 4e7; const Mp = 300e3; const mass = 1800;
  const ss = equivalentSdof({ support: 'simple', L, EI, Mp, mass });
  close(ss.k, 384 * EI / (5 * L ** 3), 1e-12); close(ss.Rm, 8 * Mp / L, 1e-12); close(ss.ue, ss.Rm / ss.k, 1e-12);
  // exact fundamental period of a simply supported beam: T = (2/π) L² √(m̄/EI); the equivalent system is within 0.1 %
  close(ss.period, 2 / Math.PI * L * L * Math.sqrt(mass / L / EI), 1e-3);
  const fx = equivalentSdof({ support: 'fixed', L, EI, Mp, mass });
  close(fx.kElastic, 384 * EI / L ** 3, 1e-12); close(fx.Re, 12 * Mp / L, 1e-12); close(fx.Rm, 16 * Mp / L, 1e-12);
  close(fx.k, 307.2 * EI / L ** 3, 1e-12); // equal-area bilinear: 16/(2(1/12 − 0.916667/16)) = 307.2
  assert.equal(fx.stages.length, 2); close(fx.stages[1].k, 384 * EI / (5 * L ** 3), 1e-12);
  close(fx.uMechanism, Mp * L * L / (12 * EI), 1e-12);
  // exact clamped–clamped period: T = 2π L²/(4.73004²) √(m̄/EI); elastic-stage equivalent within 1 %
  close(2 * Math.PI * Math.sqrt(fx.KLM.elastic * mass / fx.kElastic), 2 * Math.PI * L * L / 4.73004074 ** 2 * Math.sqrt(mass / L / EI), 0.01);
  // weak midspan (Mpm < Mps/2): midspan hinge first, then two cantilevers 128 EI/L³
  const wk = equivalentSdof({ support: 'fixed', L, EI, Mp: 100e3, MpSupport: 300e3, mass });
  close(wk.Re, 24 * 100e3 / L, 1e-12); close(wk.stages[1].k, 128 * EI / L ** 3, 1e-12); close(wk.Rm, 8 * 400e3 / L, 1e-12);
  // point loads and cantilever
  close(equivalentSdof({ support: 'simple', load: 'point', L, EI, Mp, mass }).k, 48 * EI / L ** 3, 1e-12);
  const fp = equivalentSdof({ support: 'fixed', load: 'point', L, EI, Mp, mass });
  close(fp.k, 192 * EI / L ** 3, 1e-12); close(fp.Rm, 8 * Mp / L, 1e-12); assert.equal(fp.stages.length, 1);
  const ca = equivalentSdof({ support: 'cantilever', L, EI, Mp, mass });
  close(ca.k, 8 * EI / L ** 3, 1e-12); close(ca.Rm, 2 * Mp / L, 1e-12);
  // exact cantilever period 2π L²/1.8751² √(m̄/EI): static-shape Rayleigh estimate within 1 %
  close(ca.period, 2 * Math.PI * L * L / 1.87510407 ** 2 * Math.sqrt(mass / L / EI), 0.01);
});

test('elastic SDOF: suddenly applied constant load gives DLF = 2; rectangular pulse DLF = 2 sin(π td/T) for td < T/2 (Biggs §2.3)', () => {
  const mass = 500; const k = 2e6; const T = 2 * Math.PI * Math.sqrt(mass / k);
  const step = sdofBlastResponse({ mass, k, load: { type: 'rectangular', peak: 1e4, duration: 5 * T } });
  close(step.dlf, 2, 1e-4); close(step.uMax, 2 * 1e4 / k, 1e-4); close(step.tMax, T / 2, 5e-3);
  close(step.period, T, 1e-12);
  for (const ratio of [0.1, 0.25, 0.4]) {
    const r = sdofBlastResponse({ mass, k, load: { type: 'rectangular', peak: 1e4, duration: ratio * T } });
    close(r.dlf, 2 * Math.sin(Math.PI * ratio), 2e-4, `td/T = ${ratio}`);
  }
});

test('elastic SDOF: impulsive limit x_max = I/(m ω); triangular pulse against the closed-form Duhamel solution', () => {
  const mass = 800; const k = 5e6; const w = Math.sqrt(k / mass); const T = 2 * Math.PI / w;
  // very short pulse: td = T/200
  const P = 4e5; const td = T / 200;
  const r = sdofBlastResponse({ mass, k, load: { type: 'triangular', peak: P, duration: td } });
  close(r.uMax, 0.5 * P * td / (mass * w), 2e-4, 'I/(mω)');
  // ideal impulse as an initial velocity
  const I = 300;
  close(sdofBlastResponse({ mass, k, load: () => 0, v0: I / mass, tEnd: T }).uMax, I / (mass * w), 1e-4);
  // triangular load, td = T: for t ≤ td  u = (P/k)[1 − cos ωt + sin ωt/(ω td) − t/td]; maximum of the closed form by dense sampling
  const tri = sdofBlastResponse({ mass, k, load: { type: 'triangular', peak: P, duration: T } });
  let best = 0;
  for (let i = 0; i <= 20000; i++) { const t = T * i / 20000; best = Math.max(best, 1 - Math.cos(w * t) + Math.sin(w * t) / (w * T) - t / T); }
  close(tri.dlf, best, 1e-4);
  near(tri.dlf, 1.55, 0.01, 'Biggs Fig. 2.7: DLF ≈ 1.55 at td/T = 1');
  // long triangular pulse tends to 2
  assert.ok(sdofBlastResponse({ mass, k, load: { type: 'triangular', peak: P, duration: 60 * T } }).dlf > 1.97);
  // energy is conserved in free vibration (velocity Verlet): ½mv² + ½ku² at the end equals the initial kinetic energy
  const fv = sdofBlastResponse({ mass, k, load: () => 0, v0: 1, tEnd: 10 * T });
  const n = fv.t.length - 1;
  close(0.5 * mass * fv.v[n] ** 2 + 0.5 * k * fv.u[n] ** 2, 0.5 * mass, 1e-4);
  // damping reduces the peak: logarithmic decrement over one cycle exp(−2πζ/√(1 − ζ²))
  const d = sdofBlastResponse({ mass, k, zeta: 0.05, load: () => 0, v0: 1, tEnd: 3 * T });
  close(d.uMax, Math.exp(-0.05 / Math.sqrt(1 - 0.0025) * Math.atan(Math.sqrt(1 - 0.0025) / 0.05)) / (w * 1), 1e-3);
});

test('elastic–plastic SDOF: energy solutions — step load μ = 1/(2(1 − F/Rm)); impulse μ = ½[1 + (Iω/Rm)²] (Biggs §2.7)', () => {
  const mass = 1000; const k = 4e6; const Rm = 1e5; const w = Math.sqrt(k / mass); const T = 2 * Math.PI / w;
  for (const ratio of [0.6, 0.8, 0.9]) {
    const r = sdofBlastResponse({ mass, k, Rm, load: { type: 'rectangular', peak: ratio * Rm, duration: 400 * T } });
    close(r.ductility, 1 / (2 * (1 - ratio)), 2e-3, `F/Rm = ${ratio}`);
  }
  // below half the resistance the response stays elastic: μ = 2F/Rm
  close(sdofBlastResponse({ mass, k, Rm, load: { type: 'rectangular', peak: 0.4 * Rm, duration: 20 * T } }).ductility, 0.8, 1e-3);
  const I = 3 * Rm / w;
  const imp = sdofBlastResponse({ mass, k, Rm, load: () => 0, v0: I / mass, tEnd: 4 * T });
  close(imp.ductility, 0.5 * (1 + 9), 1e-3);
  // permanent set = u_max − u_e, rebound oscillates elastically about it
  close(imp.plasticDisplacement, imp.uMax - Rm / k, 1e-3);
  assert.ok(imp.uMin >= 0);
  close(strainEnergyEpp({ k, Rm, u: imp.uMax }), I * I / (2 * mass), 2e-3);
  close(strainEnergyEpp({ k, Rm, u: 0.01 }), 0.5 * k * 1e-4, 1e-12);
  // load–mass factor: plastic K_LM applies to the kinetic energy of the mechanism
  const kl = sdofBlastResponse({ mass, k, Rm, klm: 0.66, load: () => 0, v0: 2, tEnd: 4 * T });
  close(strainEnergyEpp({ k, Rm, u: kl.uMax }), 0.5 * 0.66 * mass * 4, 3e-3);
});

test('pressure–impulse curve: points reproduce the target ductility and approach the energy asymptotes', () => {
  const mass = 1000; const k = 4e6; const Rm = 1e5; const mu = 5;
  const pi = pressureImpulseCurve({ mass, k, Rm, ductility: mu, nPoints: 12, durationRatio: [0.01, 200] });
  // impulsive asymptote I = (Rm/ω)√(2μ − 1); quasi-static F = Rm(1 − 1/(2μ))
  const w = Math.sqrt(k / mass);
  close(pi.impulsiveAsymptote, Rm / w * Math.sqrt(2 * mu - 1), 1e-12);
  close(pi.quasiStaticAsymptote, Rm * (1 - 1 / (2 * mu)), 1e-12);
  close(pi.points[0].I, pi.impulsiveAsymptote, 5e-3, 'short duration → impulsive asymptote');
  close(pi.points.at(-1).F, pi.quasiStaticAsymptote, 0.02, 'long duration → quasi-static asymptote');
  // every point really produces μ = 5
  for (const p of [pi.points[3], pi.points[6], pi.points[9]]) {
    const r = sdofBlastResponse({ mass, k, Rm, load: { type: 'triangular', peak: p.F, duration: p.td } });
    close(r.ductility, mu, 5e-3, `td = ${p.td}`);
  }
  // the curve is monotonic: load falls and impulse rises with duration; all points lie outside both asymptotes
  for (let i = 1; i < pi.points.length; i++) {
    assert.ok(pi.points[i].F < pi.points[i - 1].F && pi.points[i].I > pi.points[i - 1].I);
    assert.ok(pi.points[i].F > pi.quasiStaticAsymptote * 0.999 && pi.points[i].I > pi.impulsiveAsymptote * 0.999);
  }
  // elastic target with a rectangular pulse: asymptotes I = u√(km), F = k u/2
  const el = pressureImpulseCurve({ mass, k, uTarget: 0.01, shape: 'rectangular', nPoints: 6 });
  close(el.impulsiveAsymptote, 0.01 * Math.sqrt(k * mass), 1e-12); close(el.quasiStaticAsymptote, 0.5 * k * 0.01, 1e-12);
  close(el.points.at(-1).F, el.quasiStaticAsymptote, 1e-3);
  assert.throws(() => pressureImpulseCurve({ mass, k, ductility: 3 }), RangeError);
});

test('beam blast response: simply supported steel beam, ductility, support rotation, reactions bounded by statics', () => {
  // span 5 m, EI = 2.1e11·8.36e-5, Mp = 1.1·355e6·6.28e-4 (dynamic), 3 m tributary width carrying 0.5 bar over 50 ms, mass 2500 kg
  const L = 5; const EI = 2.1e11 * 8.356e-5; const Mp = 1.1 * 355e6 * 6.284e-4; const mass = 2500;
  const F0 = 0.5e5 * 3 * L;
  const r = beamBlastResponse({ support: 'simple', L, EI, Mp, mass, history: { type: 'triangular', peak: F0, duration: 0.05 } });
  assert.ok(r.ductility > 1 && r.ductility < 50, `μ = ${r.ductility}`);
  close(r.rotationDeg, Math.atan(r.uMax / 2.5) * 180 / Math.PI, 1e-12);
  // V = aR + bF ≤ 0.39 Rm + 0.11 F0 (elastic coefficients) and at least the plastic value 0.375 Rm at some instant after the load has decayed
  assert.ok(r.reactionMax <= 0.3935 * r.sdof.Rm + 0.1066 * F0 + 1);
  assert.ok(r.reactionMax >= 0.375 * r.sdof.Rm * 0.999);
  // a stronger beam of the same stiffness deflects less
  const strong = beamBlastResponse({ support: 'simple', L, EI, Mp: 3 * Mp, mass, history: { type: 'triangular', peak: F0, duration: 0.05 } });
  assert.ok(strong.uMax < r.uMax);
  close(supportRotation({ deflection: 0.1, L: 4 }).rotationDeg, 2.8624, 1e-4);
  close(supportRotation({ deflection: 0.1, L: 4, support: 'cantilever' }).rotation, Math.atan(0.025), 1e-12);
});

test('load histories: triangular with rise time, table interpolation, impulse', () => {
  const tri = loadHistory({ type: 'triangular', peak: 100, duration: 0.1, riseTime: 0.02 });
  close(tri.f(0.01), 50, 1e-12); close(tri.f(0.06), 50, 1e-12); assert.equal(tri.f(0.2), 0); close(tri.impulse, 5, 1e-12);
  const tab = loadHistory({ type: 'table', t: [0, 0.01, 0.03], F: [0, 200, 0] });
  close(tab.f(0.005), 100, 1e-12); close(tab.f(0.02), 100, 1e-12); close(tab.impulse, 3, 1e-12); assert.equal(tab.peak, 200);
  // table and triangular descriptions of the same pulse give the same response
  const a = sdofBlastResponse({ mass: 100, k: 1e6, load: { type: 'triangular', peak: 200, duration: 0.03, riseTime: 0.01 }, dt: 1e-5 });
  const b = sdofBlastResponse({ mass: 100, k: 1e6, load: { type: 'table', t: [0, 0.01, 0.03], F: [0, 200, 0] }, dt: 1e-5 });
  close(a.uMax, b.uMax, 1e-9);
  assert.throws(() => loadHistory({ type: 'table', t: [0, 0], F: [1, 1] }), RangeError);
});

test('Cowper–Symonds dynamic increase factor: mild steel constants D = 40.4 1/s, q = 5 give DIF = 2 at 40.4 1/s; drag load by hand', () => {
  close(dynamicIncreaseFactor({ strainRate: 40.4, D: 40.4, q: 5 }), 2, 1e-12);
  // ε̇ = 1 1/s: 1 + (1/40.4)^0.2 = 1.477
  close(dynamicIncreaseFactor({ strainRate: 1, D: 40.4, q: 5 }), 1.477, 1e-3);
  assert.equal(dynamicIncreaseFactor({ strainRate: 0, D: 40.4, q: 5 }), 1);
  close(meanStrainRate({ fy: 355e6, E: 2.1e11, timeToYield: 0.01 }), 0.169, 1e-3);
  // 219 mm pipe, gas 1.5 kg/m³ at 150 m/s, Cd = 1.2: q = 16 875 Pa, F = 1.2·16 875·0.219 = 4434.8 N/m
  const d = blastDragLoad({ rho: 1.5, velocity: 150, Cd: 1.2, D: 0.219, length: 4 });
  close(d.dynamicPressure, 16875, 1e-12); close(d.forcePerLength, 4434.75, 1e-9); close(d.force, 17739, 1e-9);
  close(blastDragLoad({ dynamicPressure: 1e4, Cd: 2, area: 3 }).force, 6e4, 1e-12);
});
