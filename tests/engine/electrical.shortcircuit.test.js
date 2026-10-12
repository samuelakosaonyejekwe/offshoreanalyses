import test from 'node:test';
import assert from 'node:assert/strict';
import {
  adiabaticShortCircuit, adiabaticMinimumArea, nonAdiabaticFactorConductor, nonAdiabaticFactorSheath, faultWithstand,
  gridImpedance, transformerImpedance, peakFactor, heatFactorM, threePhaseFault, breakingFactorMu, busFault,
  inductionMachineImpedance, synchronousGeneratorImpedance, earthFaultSplit, shortCircuitForce, seriesImpedance,
} from '../../src/engine/electrical/shortcircuit.js';

const rel = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${a} vs ${b} (tol ${tol})`);

test('IEC 60949 adiabatic rating by hand: 630 mm² Cu, 90 → 250 °C, 1 s: 226·630·sqrt(ln(484.5/324.5)) = 90.14 kA', () => {
  // ln(484.5/324.5) = 0.400830, sqrt = 0.633111, 226·630 = 142380 → 90142 A
  const r = adiabaticShortCircuit({ material: 'Cu', area: 630e-6, thetaInitialC: 90, thetaFinalC: 250, duration: 1 });
  rel(r.I, 90142, 1e-4);
  rel(r.I2t, 90142 ** 2, 2e-4);
  rel(r.currentDensity, 143.08e6, 1e-3);
  // 0.5 s allows √2 more current; aluminium of the same size carries 148/226 · sqrt(ln(478/318)/ln(484.5/324.5))
  rel(adiabaticShortCircuit({ material: 'Cu', area: 630e-6, thetaInitialC: 90, thetaFinalC: 250, duration: 0.5 }).I, 90142 * Math.SQRT2, 1e-4);
  const al = adiabaticShortCircuit({ material: 'Al', area: 630e-6, thetaInitialC: 90, thetaFinalC: 250, duration: 1 });
  rel(al.I, 148 * 630 * Math.sqrt(Math.log(478 / 318)), 1e-12);
});

test('lead sheath and steel armour constants: K, β = 41, 230 and 78, 202', () => {
  // lead sheath 500 mm², 80 → 200 °C, 1 s: 41·500·sqrt(ln(430/310)) = 11.73 kA
  rel(adiabaticShortCircuit({ material: 'lead', area: 500e-6, thetaInitialC: 80, thetaFinalC: 200, duration: 1 }).I, 41 * 500 * Math.sqrt(Math.log(430 / 310)), 1e-12);
  rel(adiabaticShortCircuit({ material: 'lead', area: 500e-6, thetaInitialC: 80, thetaFinalC: 200, duration: 1 }).I, 11727, 1e-3);
  rel(adiabaticShortCircuit({ material: 'steel', area: 1000e-6, thetaInitialC: 60, thetaFinalC: 200, duration: 1 }).I, 78 * 1000 * Math.sqrt(Math.log(402 / 262)), 1e-12);
});

test('minimum area is the inverse of the adiabatic rating', () => {
  rel(adiabaticMinimumArea({ material: 'Cu', I: 90142, duration: 1, thetaInitialC: 90, thetaFinalC: 250 }), 630e-6, 1e-4);
  rel(adiabaticMinimumArea({ material: 'Cu', I: 40e3, duration: 0.5, thetaInitialC: 90, thetaFinalC: 250 }), (40e3 * Math.sqrt(0.5)) / (226 * 0.633111) * 1e-6, 1e-5);
});

test('non-adiabatic factors: conductor ε = sqrt(1 + 0.41·sqrt(1/630) + 0.12/630) = 1.00823 for 630 mm² Cu/XLPE, 1 s; sheath polynomial', () => {
  rel(nonAdiabaticFactorConductor({ area: 630e-6, duration: 1, material: 'Cu' }), 1.008228, 1e-5);
  assert.ok(nonAdiabaticFactorConductor({ area: 50e-6, duration: 5, material: 'Al' }) > nonAdiabaticFactorConductor({ area: 630e-6, duration: 1, material: 'Al' }));
  assert.throws(() => nonAdiabaticFactorConductor({ area: 1e-4, duration: 1, insulation: 'pvc' }), RangeError);
  // sheath: M = 0.7·2·sqrt(2.4e6/3.5)/(2·1.45e6·0.003) = 0.133254 s^-0.5; ε(1 s) = 1 + 0.61M − 0.069M² + 0.0043M³
  const s = nonAdiabaticFactorSheath({ thickness: 0.003, duration: 1, sigma1: 1.45e6, sigma2: 2.4e6, rho2: 3.5, sigma3: 2.4e6, rho3: 3.5 });
  rel(s.M, 0.133254, 1e-5);
  rel(s.epsilon, 1 + 0.61 * 0.133254 - 0.069 * 0.133254 ** 2 + 0.0043 * 0.133254 ** 3, 1e-6);
});

test('fault withstand: utilisation and the adiabatic temperature actually reached', () => {
  const w = faultWithstand({ Ifault: 45071, material: 'Cu', area: 630e-6, thetaInitialC: 90, thetaFinalC: 250, duration: 1 });
  rel(w.utilisation, 0.5, 1e-4); assert.equal(w.ok, true);
  // I²t is a quarter → ln ratio is a quarter: θ = 324.5·(484.5/324.5)^0.25 − 234.5 = 124.2 °C
  rel(w.thetaFinalAdiabaticC, 324.5 * (484.5 / 324.5) ** 0.25 - 234.5, 1e-4);
  assert.equal(faultWithstand({ Ifault: 100e3, material: 'Cu', area: 630e-6, thetaInitialC: 90, thetaFinalC: 250, duration: 1 }).ok, false);
});

test('IEC 60909 peak factor κ = 1.02 + 0.98·e^(−3R/X): R/X = 0.1 → 1.746; limits 2.0 and 1.02', () => {
  rel(peakFactor(0.1), 1.746, 1e-4);
  rel(peakFactor(0), 2.0, 1e-12);
  rel(peakFactor(0.07), 1.02 + 0.98 * Math.exp(-0.21), 1e-12);
  assert.ok(peakFactor(50) - 1.02 < 1e-12);
});

test('IEC 60909 initial short-circuit current by hand: Un = 33 kV, Zk = j1 Ω, c = 1.1 → 20.958 kA; peak = κ·√2·I″k', () => {
  const f = threePhaseFault({ Un: 33e3, Zk: { re: 0.1, im: 1 }, c: 1.1, f: 50, Tk: 1 });
  rel(f.Ik, (1.1 * 33e3) / (Math.sqrt(3) * Math.hypot(0.1, 1)), 1e-12);
  rel(threePhaseFault({ Un: 33e3, Zk: { re: 0, im: 1 } }).Ik, 20957.8, 1e-5);
  rel(f.kappa, 1.746, 1e-4);
  rel(f.ip, 1.746 * Math.SQRT2 * f.Ik, 1e-4);
  rel(f.Sk, Math.sqrt(3) * 33e3 * f.Ik, 1e-12);
  // minimum current uses c = 1.0
  rel(threePhaseFault({ Un: 33e3, Zk: { re: 0.1, im: 1 }, c: 1.0 }).Ik / f.Ik, 1 / 1.1, 1e-12);
});

test('heat factor m equals the normalised Joule integral of the decaying DC component i = √2·I″k·exp(−t/τ), 1/τ = −2f·ln(κ − 1)', () => {
  const kappa = 1.746; const f = 50; const Tk = 0.2;
  const tau = -1 / (2 * f * Math.log(kappa - 1));
  // closed form of (1/Tk)·∫ 2·exp(−2t/τ) dt and a Simpson quadrature of the same integrand
  const closed = (tau / Tk) * (1 - Math.exp((-2 * Tk) / tau));
  const n = 2000; let s = 0;
  for (let i = 0; i <= n; i++) { const t = (Tk * i) / n; s += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * 2 * Math.exp((-2 * t) / tau); }
  s *= Tk / n / 3 / Tk;
  const m = heatFactorM({ kappa, f, Tk });
  rel(m, closed, 1e-12); rel(m, s, 1e-9);
  // thermal equivalent current far from generators: Ith = I″k·sqrt(m + 1)
  const flt = threePhaseFault({ Un: 33e3, Zk: { re: 0.1, im: 1 }, Tk: 0.2 });
  rel(flt.Ith, flt.Ik * Math.sqrt(flt.m + 1), 1e-12);
  assert.ok(flt.m > 0.1 && flt.m < 0.2);
});

test('network feeder and transformer impedances (IEC 60909-0 §6.2, §6.3) by hand', () => {
  // 132 kV, S″kQ = 5000 MVA: ZQ = 1.1·132²/5000 = 3.83328 Ω; XQ = ZQ/sqrt(1.01) = 3.81426 Ω; RQ = 0.1·XQ
  const zq = gridImpedance({ Un: 132e3, Sk: 5000e6 });
  rel(Math.hypot(zq.re, zq.im), 3.83328, 1e-6); rel(zq.im, 3.81426, 1e-5); rel(zq.re, 0.381426, 1e-5);
  // transformer 100 MVA, 132 kV side, uk 12 %, Pk 400 kW: ZT = 0.12·174.24 = 20.9088 Ω, RT = 0.004·174.24
  const t = transformerImpedance({ Sr: 100e6, Ur: 132e3, uk: 0.12, Pk: 400e3 });
  const R = 0.004 * 174.24; const X = Math.sqrt(20.9088 ** 2 - R * R);
  const KT = (0.95 * 1.1) / (1 + 0.6 * (X / 174.24));
  rel(t.KT, KT, 1e-12); rel(t.Z.im, KT * X, 1e-12); rel(t.Z.re, KT * R, 1e-12);
  rel(t.KT, 0.9748, 1e-3);
  assert.equal(transformerImpedance({ Sr: 100e6, Ur: 132e3, uk: 0.12, correction: false }).KT, 1);
  // grid + transformer referred to the same voltage in series
  const z = seriesImpedance([zq, t.Z]);
  rel(z.im, zq.im + t.Z.im, 1e-12);
});

test('breaking factor μ: ratio 4, tmin 0.1 s → 0.62 + 0.72·e^(−1.28) = 0.8202; μ = 1 for I″kG/IrG ≤ 2', () => {
  rel(breakingFactorMu({ ratio: 4, tmin: 0.1 }), 0.820187, 1e-5);
  assert.equal(breakingFactorMu({ ratio: 1.5, tmin: 0.1 }), 1);
  rel(breakingFactorMu({ ratio: 6, tmin: 0.02 }), 0.84 + 0.26 * Math.exp(-1.56), 1e-12);
  const a = breakingFactorMu({ ratio: 6, tmin: 0.05 }); const b = breakingFactorMu({ ratio: 6, tmin: 0.1 });
  rel(breakingFactorMu({ ratio: 6, tmin: 0.075 }), 0.5 * (a + b), 1e-12);
});

test('busbar fault with grid, induction machines and converter-limited sources: contributions add', () => {
  const grid = { name: 'Grid via transformer', Z: { re: 0.2, im: 2.0 } };
  const zm = inductionMachineImpedance({ Ur: 33e3, Sr: 20e6, lockedRotorRatio: 6 });
  rel(Math.hypot(zm.re, zm.im), 33e3 ** 2 / (20e6 * 6), 1e-12);
  const In = 400e6 / (Math.sqrt(3) * 33e3);
  const out = busFault({ Un: 33e3, sources: [grid, { name: 'Motors', Z: zm }, { name: 'Wind turbines', I: 1.2 * In }] });
  const ig = (1.1 * 33e3) / (Math.sqrt(3) * Math.hypot(0.2, 2));
  rel(out.contributions[0].Ik, ig, 1e-12);
  rel(out.contributions[2].Ik, 1.2 * In, 1e-12);
  rel(out.Ik, ig + (1.1 * 33e3) / (Math.sqrt(3) * (33e3 ** 2 / 120e6)) + 1.2 * In, 1e-12);
  // synchronous generator correction KG = (Un/UrG)·cmax/(1 + x″d·sin φ)
  const gen = synchronousGeneratorImpedance({ Ur: 11e3, Sr: 50e6, xd2: 0.18, Un: 11e3, powerFactor: 0.8 });
  rel(gen.KG, 1.1 / (1 + 0.18 * 0.6), 1e-12);
  rel(gen.Z.im, gen.KG * 0.18 * (11e3 ** 2 / 50e6), 1e-12);
});

test('earth-fault current split between sheath and sea: reduction factor r = Rs/|Rs + Re + jX| by hand', () => {
  const Rs = 2e-4; const ds = 0.1; const f = 50; const rho = 0.25;
  const De = 658.87 * Math.sqrt(rho / f); const Re = Math.PI ** 2 * f * 1e-7;
  const X = 2 * Math.PI * f * 2e-7 * Math.log(De / 0.05);
  const out = earthFaultSplit({ Ifault: 20e3, Rs, ds, f, rhoEarth: rho });
  rel(out.reductionFactor, Rs / Math.hypot(Rs + Re, X), 1e-12);
  rel(out.Iearth, 20e3 * out.reductionFactor, 1e-12);
  rel(out.Isheath, (20e3 * Math.hypot(Re, X)) / Math.hypot(Rs + Re, X), 1e-12);
  // most of the fault current returns in the metallic path at power frequency
  assert.ok(out.Isheath > out.Iearth);
  const fixed = earthFaultSplit({ Ifault: 20e3, reductionFactor: 0.3 });
  rel(fixed.Iearth, 6000, 1e-12);
});

test('electromagnetic force by hand: ip = 50 kA, spacing 0.1 m → 2e-7·2.5e9/0.1 = 5000 N/m (two-phase), ×√3/2 (three-phase)', () => {
  rel(shortCircuitForce({ ip: 50e3, spacing: 0.1, fault: 'twoPhase' }), 5000, 1e-12);
  rel(shortCircuitForce({ ip: 50e3, spacing: 0.1, fault: 'threePhase' }), 5000 * Math.sqrt(3) / 2, 1e-12);
  rel(shortCircuitForce({ ip: 50e3, spacing: 0.1, k: 0.5 }), 2500, 1e-12);
});
