import test from 'node:test';
import assert from 'node:assert/strict';
import {
  dispersion, wavePower, regularWavePower, waveSpectrum, spectralPeriods, scatterResource, monthlyResource,
  cylinderHeave, pointAbsorber, captureWidth, heaveLimits, cylinderBody, spectralResponse, powerMatrix,
  matrixEnergy, radiationKernel, infiniteAddedMass, simulateCummins, twoBodyAbsorber, flapHydrodynamics2D,
  flapStiffness, flapConverter, captureWidthDevice, overtoppingPower, arrayPower, survivalResponse,
  slammingForce,
} from '../../src/engine/aero/wec.js';

const G = 9.80665;
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);
const grid = (a, b, n) => Float64Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));

test('wave power: Hs = 2 m, Te = 8 s in deep water is 15.7 kW/m (ρ = 1025)', () => {
  // ρ g² Hs² Te/(64π) = 1025 · 96.17038 · 4 · 8 / 201.0619 = 15 689 W/m
  const J = wavePower({ Hm0: 2, Te: 8 });
  close(J, (1025 * G * G * 32) / (64 * Math.PI), 1e-12);
  close(J / 1000, 15.69, 1e-3);
  // finite depth converges to deep water and exceeds it at intermediate depth (c_g > deep value for kh ≈ 1)
  close(wavePower({ Hm0: 2, Te: 8, depth: 2000 }), J, 1e-9);
  assert.ok(wavePower({ Hm0: 2, Te: 8, depth: 30 }) > J);
  // regular wave H = 2 m, T = 8 s: ρ g² H² T/(32π) = 2× the irregular formula with Hs = H
  close(regularWavePower({ amplitude: 1, omega: (2 * Math.PI) / 8 }), 2 * J, 1e-12);
});

test('dispersion relation and group velocity limits', () => {
  const w = (2 * Math.PI) / 10;
  const d = dispersion(w, 20);
  close(w * w, G * d.k * Math.tanh(d.k * 20), 1e-12);
  close(dispersion(w).wavelength, (G * 100) / (2 * Math.PI), 1e-12);
  // shallow water: c = c_g = √(g h)
  const sh = dispersion((2 * Math.PI) / 200, 5);
  close(sh.cg, Math.sqrt(G * 5), 2e-3);
  close(sh.c, sh.cg, 5e-3);
});

test('wave spectra: m0 = Hs²/16 and Te/Tp = 0.857 for Bretschneider', () => {
  const om = grid(0.05, 12, 24000);
  const S = waveSpectrum({ Hs: 3, Tp: 10, gamma: 1, omega: om });
  const p = spectralPeriods({ omega: om, S });
  close(p.Hm0, 3, 1e-3);
  // Te = Tp (4/5)^(1/4) Γ(5/4) = 0.85726 Tp;  Tz = Tp (5/4)^(-1/4) ... = 0.7104 Tp (m2 by integration)
  close(p.Te / 10, 0.8573, 1e-3);
  close(p.Tz / 10, 0.7104, 5e-3);
  // JONSWAP γ = 3.3 with the 1 − 0.287 ln γ normalisation keeps Hs within a few per cent and raises Te/Tp
  const Sj = waveSpectrum({ Hs: 3, Tp: 10, gamma: 3.3, omega: om });
  const pj = spectralPeriods({ omega: om, S: Sj });
  close(pj.Hm0, 3, 0.03);
  assert.ok(pj.Te / 10 > 0.88 && pj.Te / 10 < 0.93);
});

test('scatter-diagram resource', () => {
  const hs = [1, 2, 3];
  const te = [6, 8];
  const prob = [[0, 0], [0, 5], [0, 0]];
  close(scatterResource({ hs, te, prob }).meanJ, wavePower({ Hm0: 2, Te: 8 }), 1e-12);
  const mixed = scatterResource({ hs, te, prob: [[1, 1], [1, 1], [0, 0]] });
  // mean of J ∝ Hs² Te: (1·6 + 1·8 + 4·6 + 4·8)/4 = 17.5 against 32 for the (2 m, 8 s) cell
  close(mixed.meanJ / wavePower({ Hm0: 2, Te: 8 }), 17.5 / 32, 1e-12);
  close(mixed.contribution.flat().reduce((s, x) => s + x, 0), 1, 1e-12);
  const m = monthlyResource({ hs, te, monthly: [prob, [[1, 0], [0, 0], [0, 0]]] });
  close(m[1], wavePower({ Hm0: 1, Te: 6 }), 1e-12);
});

test('truncated-cylinder coefficients: hydrostatics, long-wave excitation, Haskind damping', () => {
  const c = cylinderHeave({ radius: 5, draft: 4, omega: 0.9 });
  close(c.K, 1025 * G * Math.PI * 25, 1e-12);
  close(c.A, (1025 * 1000) / 6, 1e-12);
  close(c.displacedMass, 1025 * Math.PI * 25 * 4, 1e-12);
  // deep water Haskind: B = ω k |Fe|²/(2 ρ g²)
  close(c.B, (0.9 * c.k * c.Fe ** 2) / (2 * 1025 * G * G), 1e-10);
  // long waves: Fe → ρ g Aw
  const lw = cylinderHeave({ radius: 5, draft: 4, omega: 0.05 });
  close(lw.Fe, lw.K, 2e-3);
  // short waves see less of the bottom pressure
  assert.ok(cylinderHeave({ radius: 5, draft: 4, omega: 1.5 }).Fe < c.Fe);
  // finite depth reduces to deep water
  close(cylinderHeave({ radius: 5, draft: 4, omega: 0.9, depth: 500 }).Fe, c.Fe, 1e-6);
});

test('point absorber: reactive control gives |Fe|²/(8B); at resonance the passive optimum does too; capture width = λ/2π', () => {
  const radius = 5;
  const draft = 4;
  const body = cylinderBody({ radius, draft });
  const probe = cylinderHeave({ radius, draft, omega: 1 });
  const mass = probe.displacedMass;
  // resonance (added mass is frequency independent here)
  const w0 = Math.sqrt(probe.K / (mass + probe.A));
  const c = cylinderHeave({ radius, draft, omega: w0 });
  const zeta = 0.5;
  const pMax = (c.Fe * zeta) ** 2 / (8 * c.B);
  const passive = pointAbsorber({ mass, A: c.A, B: c.B, K: c.K, Fe: c.Fe, omega: w0, amplitude: zeta });
  close(passive.power, pMax, 1e-10);
  close(passive.Bpto, c.B, 1e-9);
  // theoretical maximum for a heaving axisymmetric body: J/k, i.e. capture width λ/2π
  const J = regularWavePower({ amplitude: zeta, omega: w0 });
  close(captureWidth({ power: passive.power, J }).captureWidth, 1 / c.k, 1e-9);
  close(captureWidth({ power: passive.power, J }).captureWidth, dispersion(w0).wavelength / (2 * Math.PI), 1e-9);
  const lim = heaveLimits({ omega: w0, H: 2 * zeta, volume: Math.PI * 25 * 4 });
  close(lim.radiationLimit, pMax, 1e-9);
  // off resonance: reactive control still reaches the maximum, the passive damper does not
  for (const w of [0.6, 0.9, 1.6]) {
    const cc = cylinderHeave({ radius, draft, omega: w });
    const o = { mass, A: cc.A, B: cc.B, K: cc.K, Fe: cc.Fe, omega: w, amplitude: zeta };
    const re = pointAbsorber({ ...o, control: 'reactive' });
    const pa = pointAbsorber({ ...o, control: 'passive' });
    close(re.power, (cc.Fe * zeta) ** 2 / (8 * cc.B), 1e-10);
    close(re.power, re.powerMax, 1e-12);
    const k = regularWavePower({ amplitude: zeta, omega: w });
    assert.ok(re.power / k <= 1 / cc.k + 1e-9);
    // passive optimum closed form: P = ¼ |F|² / (B + √(B² + X²)), X = ω M − K/ω
    const X = w * (mass + cc.A) - cc.K / w;
    close(pa.power, (0.25 * (cc.Fe * zeta) ** 2) / (cc.B + Math.hypot(cc.B, X)), 1e-10);
    assert.ok(pa.power < re.power);
    // any other damping absorbs less
    for (const f of [0.5, 2]) {
      assert.ok(pointAbsorber({ ...o, control: 'fixed', Bpto: pa.Bpto * f }).power < pa.power);
    }
  }
  assert.ok(body(w0).M > mass);
});

test('point absorber with stroke and force limits', () => {
  const c = cylinderHeave({ radius: 5, draft: 4, omega: 0.8 });
  const o = { mass: c.displacedMass, A: c.A, B: c.B, K: c.K, Fe: c.Fe, omega: 0.8, amplitude: 1 };
  const free = pointAbsorber({ ...o, control: 'reactive' });
  assert.ok(free.stroke > 3);
  const lim = pointAbsorber({ ...o, control: 'reactive', limits: { stroke: 2 } });
  assert.equal(lim.constrained, true);
  assert.equal(lim.feasible, true);
  close(lim.stroke, 2, 1e-6);
  assert.ok(lim.power < free.power && lim.power > 0);
  // at resonance of the controlled system |X| = F/(ω (B + Bpto)) → Bpto = F/(ω x_max) − B
  close(lim.Bpto, c.Fe / (0.8 * 2) - c.B, 1e-5);
  close(lim.power, 0.5 * lim.Bpto * (0.8 * 2) ** 2, 1e-6);
  const impossible = pointAbsorber({ ...o, control: 'passive', limits: { stroke: 0.5, force: 1 } });
  assert.equal(impossible.feasible, false);
  assert.equal(impossible.power, 0);
});

test('Budal upper bound by hand', () => {
  // V = 100 m³, H = 2 m, T = 8 s: (π/4)·1025·9.80665·100·2/8 = 197 364 W
  const l = heaveLimits({ omega: (2 * Math.PI) / 8, H: 2, volume: 100 });
  close(l.budalLimit, ((Math.PI / 4) * 1025 * G * 100 * 2) / 8, 1e-12);
  close(l.budalLimit, 197364, 1e-4);
  close(l.maxCaptureWidth, (G * 64) / (4 * Math.PI * Math.PI), 1e-12);
  assert.equal(l.limit, Math.min(l.budalLimit, l.radiationLimit));
});

test('irregular waves: P = Bpto σv², quadratic in Hs; power matrix and annual energy', () => {
  const body = cylinderBody({ radius: 5, draft: 4 });
  const om = grid(0.2, 3, 400);
  const S1 = waveSpectrum({ Hs: 1, Tp: 8, omega: om });
  const S2 = waveSpectrum({ Hs: 2, Tp: 8, omega: om });
  const r1 = spectralResponse({ omega: om, S: S1, body, Bpto: 2e5 });
  const r2 = spectralResponse({ omega: om, S: S2, body, Bpto: 2e5 });
  close(r1.power, 2e5 * r1.sigmaV ** 2, 1e-12);
  close(r2.power, 4 * r1.power, 1e-10);
  close(r1.sigmaF, 2e5 * r1.sigmaV, 1e-10);
  // mean power cannot exceed the spectral sum of the axisymmetric limit J(ω)/k(ω) = ρ g³ S(ω)/(2 ω³) dω
  let bound = 0;
  for (let i = 1; i < om.length; i++) {
    const f = (j) => (1025 * G ** 3 * S1[j]) / (2 * om[j] ** 3);
    bound += 0.5 * (f(i) + f(i - 1)) * (om[i] - om[i - 1]);
  }
  assert.ok(r1.power < bound);
  const pm = powerMatrix({ hs: [1, 2, 4], tp: [6, 8, 10], body, ratedPower: 190e3 });
  close(pm.power[1][1] / pm.power[0][1], 4, 1e-9);
  assert.ok(pm.power[2].every((p) => p <= 190e3));
  assert.equal(Math.max(...pm.power[2]), 190e3);
  const prob = [[2, 2, 0], [0, 4, 0], [0, 1, 1]];
  const e = matrixEnergy({ power: pm.power, prob, ratedPower: 190e3, availability: 0.9, devices: 10, arrayFactor: 0.95 });
  const mean = (2 * pm.power[0][0] + 2 * pm.power[0][1] + 4 * pm.power[1][1] + pm.power[2][1] + pm.power[2][2]) / 10;
  close(e.meanPower, mean * 0.9 * 10 * 0.95, 1e-12);
  close(e.aepWh, e.meanPower * 8760, 1e-12);
  close(e.capacityFactor, (mean * 0.9 * 0.95) / 190e3, 1e-12);
  const sv = survivalResponse({ omega: om, S: waveSpectrum({ Hs: 8, Tp: 12, omega: om }), body, Bpto: 1e6, strokeLimit: 3 });
  close(sv.maxStroke, sv.sigmaX * Math.sqrt(2 * Math.log(sv.cycles)), 1e-12);
  assert.ok(sv.maxStroke > 2 * sv.sigmaX && sv.strokeUtilisation > 0);
});

test('radiation kernel and Ogilvie relation against the analytic pair K(t) = c e^(−βt)', () => {
  // B(ω) = c β/(β² + ω²), A(ω) = A∞ − c/(β² + ω²)
  const c = 4e4;
  const beta = 1;
  const om = grid(0, 400, 80001);
  const B = om.map((w) => (c * beta) / (beta * beta + w * w));
  const ker = radiationKernel({ omega: om, B, dt: 0.05, tMax: 12 });
  for (const t of [0.5, 1, 2, 4]) close(ker.K[Math.round(t / 0.05)], c * Math.exp(-beta * t), 5e-3 * c / Math.max(1, c * Math.exp(-beta * t)) + 1e-3, `t = ${t}`);
  const w = 1.3;
  const Ainf = 5e4;
  const Aw = Ainf - c / (beta * beta + w * w);
  close(infiniteAddedMass({ A: Aw, omega: w, kernel: ker }), Ainf, 2e-3);
});

test('Cummins equation: free-decay period equals the damped eigenfrequency (exponential kernel) and 2π√(M/K) without radiation', () => {
  const mass = 1e5;
  const Ainf = 5e4;
  const M = mass + Ainf;
  const Kh = 3e5;
  const c = 4e4;
  const beta = 1;
  const period = (r) => {
    const zc = [];
    for (let i = 1; i < r.x.length; i++) {
      if (r.x[i - 1] > 0 && r.x[i] <= 0) zc.push(r.t[i - 1] + ((r.t[i] - r.t[i - 1]) * r.x[i - 1]) / (r.x[i - 1] - r.x[i]));
    }
    return (zc[zc.length - 1] - zc[0]) / (zc.length - 1);
  };
  const und = simulateCummins({ mass, Ainf, stiffness: Kh, dt: 0.01, duration: 60, x0: 1 });
  close(period(und), 2 * Math.PI * Math.sqrt(M / Kh), 1e-6);
  close(und.x[und.x.length - 1], Math.cos(Math.sqrt(Kh / M) * 60), 1e-6);
  // characteristic equation (M s² + K)(s + β) + c s = 0: real root by bisection, then the complex pair
  const poly = (s) => M * s ** 3 + M * beta * s * s + (Kh + c) * s + Kh * beta;
  let lo = -beta;
  let hi = 0;
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (poly(lo) * poly(mid) <= 0) hi = mid;
    else lo = mid;
  }
  const s1 = 0.5 * (lo + hi);
  // M s³ + Mβ s² + (K + c) s + Kβ = M (s − s1)(s² + p s + q): p = β + s1, q = −Kβ/(M s1)
  const p = beta + s1;
  const q = (-Kh * beta) / (M * s1);
  const wd = Math.sqrt(q - (p * p) / 4);
  const dec = simulateCummins({ mass, Ainf, stiffness: Kh, kernel: (t) => c * Math.exp(-beta * t), dt: 0.005, duration: 40, x0: 1, kernelLength: 20 });
  close(period(dec), (2 * Math.PI) / wd, 2e-3);
  // decay envelope: successive positive peaks ratio exp(−p/2 · T)
  assert.ok(Math.abs(dec.x[dec.x.length - 1]) < 0.6);
  assert.ok((2 * Math.PI) / wd !== 2 * Math.PI * Math.sqrt(M / Kh));
});

test('Cummins time domain with a linear PTO reproduces the frequency-domain absorbed power', () => {
  const mass = 1e5;
  const Ainf = 5e4;
  const Kh = 3e5;
  const c = 4e4;
  const beta = 1;
  const w = 1.2;
  const F0 = 5e4;
  const Bp = 6e4;
  const Bw = (c * beta) / (beta * beta + w * w);
  const Aw = Ainf - c / (beta * beta + w * w);
  const fd = pointAbsorber({ mass, A: Aw, B: Bw, K: Kh, Fe: F0, omega: w, amplitude: 1, control: 'fixed', Bpto: Bp });
  const dt = 0.005;
  const run = (duration) => simulateCummins({ mass, Ainf, stiffness: Kh, kernel: (t) => c * Math.exp(-beta * t), excitation: (t) => F0 * Math.cos(w * t), dt, duration, pto: { type: 'linear', B: Bp }, kernelLength: 15 });
  const T = (2 * Math.PI) / w;
  const a = run(30 * T);
  const b = run(40 * T);
  // mean power over the last ten periods (transient decayed)
  const pSteady = (b.absorbedEnergy - a.absorbedEnergy) / (10 * T);
  close(pSteady, fd.power, 1e-2);
  let xmax = 0;
  for (let i = Math.round((30 * T) / dt); i < b.x.length; i++) xmax = Math.max(xmax, Math.abs(b.x[i]));
  close(xmax, fd.stroke, 1e-2);
});

test('nonlinear PTO, end stops, latching and declutching in the time domain', () => {
  const base = { mass: 1e5, Ainf: 5e4, stiffness: 3e5, dt: 0.005, duration: 120, excitation: (t) => 8e4 * Math.cos(0.9 * t) };
  // Coulomb PTO: absorbed energy = Fc ∫|v| dt (up to the tanh smoothing)
  const cou = simulateCummins({ ...base, pto: { type: 'coulomb', force: 2e4, vSmooth: 1e-3 } });
  let path = 0;
  for (let i = 1; i < cou.v.length; i++) path += 0.5 * (Math.abs(cou.v[i]) + Math.abs(cou.v[i - 1])) * 0.005;
  close(cou.absorbedEnergy, 2e4 * path, 2e-2);
  assert.ok(cou.maxPtoForce <= 2e4 * (1 + 1e-9));
  // end stops bound the stroke of the lightly damped body
  const free = simulateCummins({ ...base, pto: { type: 'linear', B: 2e3 } });
  const stop = simulateCummins({ ...base, pto: { type: 'linear', B: 2e3 }, endStop: { stroke: 0.4, stiffness: 5e7 } });
  assert.ok(free.maxStroke > 0.6);
  assert.ok(stop.maxStroke < 0.45);
  // hydraulic PTO: pressure stays positive and generated energy is consistent with absorbed energy
  const hyd = simulateCummins({ ...base, duration: 300, pto: { type: 'hydraulic', pistonArea: 0.01, p0: 2e6, V0: 0.5, motorCoeff: 2e-9 } });
  assert.ok(hyd.pressure.every((p) => p > 0));
  assert.ok(hyd.meanGeneratedPower > 0);
  // energy absorbed at the piston = energy through the motor + energy stored in the gas (isentropic, γ = 1.4)
  const pEnd = hyd.pressure[hyd.pressure.length - 1];
  const vEnd = 0.5 * (2e6 / pEnd) ** (1 / 1.4);
  const stored = (pEnd * vEnd - 2e6 * 0.5) / 0.4;
  close(hyd.absorbedEnergy, hyd.meanGeneratedPower * 300 + stored, 3e-2);
  // latching holds the body at the motion extremes; declutching removes the PTO force at low speed
  const lat = simulateCummins({ ...base, pto: { type: 'linear', B: 5e4 }, latching: { duration: 1.0 } });
  let held = 0;
  for (let i = 2000; i < lat.v.length; i++) if (lat.v[i] === 0) held++;
  assert.ok(held * 0.005 > 20, `latched time ${held * 0.005}`);
  const dec = simulateCummins({ ...base, pto: { type: 'linear', B: 5e4 }, declutch: { threshold: 0.05 } });
  for (let i = 0; i < dec.v.length; i++) if (Math.abs(dec.v[i]) < 0.05) assert.equal(dec.ptoForce[i], 0);
  assert.ok(dec.meanAbsorbedPower > 0);
  const moor = simulateCummins({ mass: 1e5, Ainf: 5e4, stiffness: 3e5, mooringStiffness: 3e5, dt: 0.01, duration: 20, x0: 1 });
  close(moor.x[moor.x.length - 1], Math.cos(2 * 20), 1e-5);
});

test('two-body absorber reduces to the single body when the reaction body is immovable', () => {
  const c = cylinderHeave({ radius: 5, draft: 4, omega: 1.0 });
  const b1 = { mass: c.displacedMass, A: c.A, B: c.B, K: c.K, Fe: c.Fe };
  const single = pointAbsorber({ ...b1, omega: 1.0, amplitude: 1, control: 'fixed', Bpto: 3e5 });
  const two = twoBodyAbsorber({ omega: 1.0, amplitude: 1, body1: b1, body2: { mass: 1e15, A: 0, B: 0, K: 1e18, Fe: 0 }, Bpto: 3e5 });
  close(two.power, single.power, 1e-6);
  close(two.x1, single.stroke, 1e-6);
  assert.ok(two.x2 < 1e-9);
  // a free reaction body moves and changes the absorbed power; the optimiser beats an arbitrary damping
  const b2 = { mass: 4e5, A: 6e5, B: 2e3, K: 1e5, Fe: { re: -2e4, im: 5e3 } };
  const fixed = twoBodyAbsorber({ omega: 1.0, amplitude: 1, body1: b1, body2: b2, Bpto: 3e5 });
  const opt = twoBodyAbsorber({ omega: 1.0, amplitude: 1, body1: b1, body2: b2, optimise: true });
  assert.ok(fixed.x2 > 0.01);
  assert.ok(opt.power >= fixed.power);
  for (const f of [0.7, 1.4]) assert.ok(twoBodyAbsorber({ omega: 1.0, amplitude: 1, body1: b1, body2: b2, Bpto: opt.Bpto * f }).power <= opt.power * (1 + 1e-9));
  // power balance: PTO power = excitation power − radiated power
  close(fixed.ptoForce, fixed.relative * 1.0 * 3e5, 1e-9);
});

test('bottom-hinged flap: wavemaker ratio limit, Haskind consistency and the 50 % two-dimensional limit', () => {
  // deep-water asymptote of the flap wavemaker: H/S → 2 (1 − 1/kh)
  const deep = flapHydrodynamics2D({ omega: 3, depth: 30 });
  close(deep.waveRatio, 2 * (1 - 1 / (deep.k * 30)), 1e-9);
  // shallow-water limit: H/S → kh/2 (Dean & Dalrymple, simplified wavemaker theory for a flap)
  const shallow = flapHydrodynamics2D({ omega: 0.05, depth: 5 });
  close(shallow.waveRatio / (shallow.k * 5), 0.5, 1e-3);
  const w = 0.8;
  const depth = 12;
  const h = flapHydrodynamics2D({ omega: w, depth });
  // 2-D Haskind: b = 2 |X|²/(4 ρ g c_g)
  close(h.damping, (2 * h.excitation ** 2) / (4 * 1025 * G * h.cg), 1e-12);
  // tuned flap (resonant, optimal damping): P = |Me ζ|²/(8B) = J w/2 → capture-width ratio 0.5
  const I = 4e6;
  const Ia = 6e6;
  const r = flapConverter({ width: 18, depth, omega: w, amplitude: 0.6, inertia: I, addedInertia: Ia, stiffness: w * w * (I + Ia) });
  close(r.captureWidthRatio, 0.5, 1e-9);
  close(r.Bpto, r.radiationDamping, 1e-9);
  const J = regularWavePower({ amplitude: 0.6, omega: w, depth });
  close(r.power, (J * 18) / 2, 1e-9);
  // detuned flap absorbs less; an angle limit constrains the motion
  const off = flapConverter({ width: 18, depth, omega: w, amplitude: 0.6, inertia: I, addedInertia: Ia, stiffness: 2e7 });
  assert.ok(off.captureWidthRatio < 0.5);
  const lim = flapConverter({ width: 18, depth, omega: w, amplitude: 0.6, inertia: I, addedInertia: Ia, stiffness: w * w * (I + Ia), maxAngle: r.angle / 2 });
  close(lim.angle, r.angle / 2, 1e-6);
  assert.ok(lim.power < r.power);
  // hydrostatic stiffness by hand: g (ρ V z_b − m z_g) = 9.80665 (1025·100·5 − 40000·4) = 3.4568e6
  close(flapStiffness({ volume: 100, buoyancyArm: 5, mass: 4e4, gravityArm: 4 }), G * (1025 * 500 - 160000), 1e-12);
});

test('parameterised devices, array factor and slamming hook by hand', () => {
  close(captureWidthDevice({ J: 30e3, captureWidthRatio: 0.2, length: 120, efficiency: 0.8 }), 576e3, 1e-12);
  assert.equal(captureWidthDevice({ J: 30e3, captureWidthRatio: 0.2, length: 120, ratedPower: 500e3 }), 500e3);
  // ρ g q Rc w η = 1025·9.80665·0.05·2·100·0.9
  close(overtoppingPower({ q: 0.05, crestFreeboard: 2, width: 100, efficiency: 0.9 }), 1025 * G * 0.05 * 2 * 100 * 0.9, 1e-12);
  assert.equal(arrayPower({ devices: 20, power: 100e3, q: 0.9 }), 1.8e6);
  close(slammingForce({ Cs: 3, area: 2, velocity: 4 }), 0.5 * 1025 * 3 * 2 * 16, 1e-12);
});
