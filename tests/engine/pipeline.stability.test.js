import test from 'node:test';
import assert from 'node:assert/strict';
import * as F from '../../src/engine/pipeline/stability.js';
import { pipeWeights } from '../../src/engine/pipeline/pipeprops.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${a} vs ${b}`);

test('dispersion relation: satisfies w^2 = g k tanh(kd); deep and shallow limits', () => {
  const g = 9.80665;
  for (const [w, d] of [[0.5, 30], [1.2, 100], [0.3, 5], [2, 1000]]) {
    const k = F.waveNumber(w, d);
    close(g * k * Math.tanh(k * d), w * w, 1e-12);
  }
  close(F.waveNumber(2, 1000), 4 / g, 1e-10);                         // deep water k = w^2/g
  close(F.waveNumber(0.05, 2), 0.05 / Math.sqrt(g * 2), 2e-3);        // shallow water k = w/sqrt(gd)
});

test('JONSWAP: zeroth moment equals Hs^2/16 and gamma follows Tp/sqrt(Hs)', () => {
  for (const gamma of [1, 3.3, 5]) {
    const sea = { Hs: 4, Tp: 10, gamma };
    const wp = 2 * Math.PI / 10;
    let m0 = 0;
    const n = 20000; const a = 0.1 * wp; const b = 20 * wp; const h = (b - a) / n;
    for (let i = 0; i <= n; i++) m0 += (i === 0 || i === n ? 0.5 : 1) * F.jonswapSpectrum(a + i * h, sea) * h;
    close(m0, 16 / 16, 0.015);                                        // Hs^2/16 = 1.0 m2 (the 1 - 0.287 ln(gamma) fit is ~1 % accurate)
  }
  assert.equal(F.jonswapGamma(4, 6), 5);                              // phi = 3.0
  assert.equal(F.jonswapGamma(4, 12), 1);                             // phi = 6.0
  close(F.jonswapGamma(4, 8.6), Math.exp(5.75 - 1.15 * 4.3), 1e-12);  // phi = 4.3
});

test('near-bed spectral velocity: shallow-water limit Us = (Hs/2) sqrt(g/d); decays with depth', () => {
  // for kd -> 0: w/sinh(kd) -> sqrt(g/d), so M0 = (g/d) Hs^2/16 and Us = 2 sqrt(M0) = (Hs/2) sqrt(g/d)
  const r = F.nearBedSpectralVelocity({ Hs: 0.2, Tp: 25, depth: 2, gamma: 3.3 });
  close(r.Us, 0.1 * Math.sqrt(9.80665 / 2), 0.03);
  close(r.Tn, Math.sqrt(2 / 9.80665), 1e-12);
  const a = F.nearBedSpectralVelocity({ Hs: 6, Tp: 12, depth: 40 });
  const b = F.nearBedSpectralVelocity({ Hs: 6, Tp: 12, depth: 120 });
  assert.ok(a.Us > 3 * b.Us);
  assert.ok(b.Tu > a.Tu);                                             // deeper water filters the short waves
  // regular-wave cross-check: a very peaked spectrum behaves like a single component at Tp
  const pk = F.nearBedSpectralVelocity({ Hs: 6, Tp: 12, depth: 40, gamma: 5 });
  const w = 2 * Math.PI / 12; const k = F.waveNumber(w, 40);
  const single = 2 * Math.sqrt((w / Math.sinh(k * 40)) ** 2 * 36 / 16);
  assert.ok(Math.abs(pk.Us / single - 1) < 0.12, `${pk.Us} vs ${single}`);
});

test('design oscillation U* and T*: hand calculation (DNV-RP-F109 Eq. 3.15-3.17)', () => {
  // tau = 10800/10.8 = 1000: sqrt(2 ln 1000) = 3.71692 ; U* = 0.5*Us*(3.71692 + 0.5772/3.71692) = 1.93611 Us
  const r = F.designOscillation({ Us: 1, Tu: 10.8, gamma: 1, Tn: 1.08 });
  close(r.Ustar, 1.93611, 1e-5);
  // Tn/Tu = 0.1, kT = 1.25: T*/Tu = 1.25 - 5*0.25*0.1 = 1.125
  close(r.Tstar, 1.125 * 10.8, 1e-12);
  close(F.designOscillation({ Us: 1, Tu: 10, gamma: 3.3, Tn: 1 }).kT, 1.21, 1e-12);
  close(F.designOscillation({ Us: 1, Tu: 10, gamma: 5, Tn: 3 }).Tstar, 10, 1e-12);      // Tn/Tu > 0.2
  close(F.designOscillation({ Us: 2, Tu: 10.8, RD: 0.9 }).Ustar, 2 * 0.9 * 1.93611, 1e-5);
});

test('current over the pipe: log-law average (DNV-RP-F109 Eq. 3.3), hand calculation', () => {
  // z0 = 1e-5 m, D = 0.4 m, zr = 3 m:
  // [(1 + 2.5e-5) ln(40001) - 1]/ln(300001) = (1.000025*10.596660 - 1)/12.611541 = 0.760974
  const v = F.currentOverPipe({ Vref: 1, zRef: 3, z0: 1e-5, D: 0.4 });
  close(v, 0.760974, 2e-5);
  close(F.currentOverPipe({ Vref: 1, zRef: 3, z0: 1e-5, D: 0.4, angleDeg: 30 }), 0.5 * v, 1e-12);
  // independent check: numerical mean of the log profile over the pipe height
  let s = 0; const n = 200000;
  for (let i = 0; i < n; i++) { const z = (i + 0.5) * 0.4 / n; s += Math.log((z + 1e-5) / 1e-5) / Math.log((3 + 1e-5) / 1e-5); }
  close(v, s / n, 1e-4);
});

test('load reduction factors for penetration, trench and permeable seabed', () => {
  const a = F.loadReductionFactors({ zp: 0.08, D: 0.4 });                   // zp/D = 0.2
  close(a.rPenY, 0.72); close(a.rPenZ, 0.87);
  const b = F.loadReductionFactors({ zp: 0.24, D: 0.4, permeable: true });  // zp/D = 0.6
  close(b.rPenY, 0.3); close(b.rPenZ, 0.35); close(b.rTotZ, 0.35 * 0.7);
  const c = F.loadReductionFactors({ zp: 0, D: 0.4, trenchAngleDeg: 21 });
  close(c.rTrY, 1 - 0.18 * 2); close(c.rTrZ, 1 - 0.14 * 16 ** 0.43);
  assert.equal(F.loadReductionFactors({ zp: 0, D: 0.4 }).rTotY, 1);
});

test('passive resistance and initial penetration: hand calculations (DNV-RP-F109 3.4.6)', () => {
  // sand: kappa_s = 9000*0.4^2/144 = 10 ; FR/FC = (50 - 15)*0.1^1.25 = 35*0.0562341 = 1.96819
  const s = F.passiveResistance({ soil: 'sand', FC: 144, D: 0.4, zp: 0.04, gammaSoilSub: 9000 });
  close(s.kappa, 10, 1e-12);
  close(s.FR / 144, 1.96819, 1e-5);
  // the two sand branches meet at kappa_s = 26.7 (5*26.7 - 0.15*26.7^2 = 26.57, within 0.5 % of 26.7)
  const lo = F.passiveResistance({ soil: 'sand', FC: 9000 * 0.16 / 26.7, D: 0.4, zp: 0.04, gammaSoilSub: 9000 });
  close(lo.FR / lo.kappa / (9000 * 0.16 / 26.7), 0.1 ** 1.25, 6e-3);
  // clay: Gc = 5000/(0.4*18000) = 0.694444 ; kappa_c = 5000*0.4/2000 = 1
  // FR/FC = 4.1*1/0.694444^0.39 * 0.1^1.31 = 4.1*1.15283*0.0489779 = 0.23150
  const c = F.passiveResistance({ soil: 'clay', FC: 2000, D: 0.4, zp: 0.04, su: 5000, gammaSoilDry: 18000 });
  close(c.FR / 2000, 0.23150, 1e-4);
  assert.equal(F.passiveResistance({ soil: 'sand', FC: -5, D: 0.4, zp: 0.04, gammaSoilSub: 9000 }).FR, 0);
  // initial penetration, sand: 0.037*0.72^-0.67 = 0.037*1.24621 = 0.046110 -> 0.018444 m
  close(F.initialPenetration({ soil: 'sand', FC: 2000, D: 0.4, gammaSoilSub: 9000 }).zpi, 0.4 * 0.037 * 0.72 ** -0.67, 1e-12);
  close(F.initialPenetration({ soil: 'sand', FC: 2000, D: 0.4, gammaSoilSub: 9000 }).zpi, 0.018444, 1e-4);
  // clay: x = 0.694444^0.3 = 0.896365 ; 0.0071*x^3.2 + 0.062*x^0.7 = 0.0071*0.70463 + 0.062*0.92627 = 0.062432 -> 0.024973 m
  close(F.initialPenetration({ soil: 'clay', FC: 2000, D: 0.4, su: 5000, gammaSoilDry: 18000 }).zpi, 0.024973, 1e-4);
});

test('vertical stability: gamma_w b/(ws + b)', () => {
  const v = F.verticalStability({ ws: 500, b: 1300 });
  close(v.uc, 1.1 * 1300 / 1800, 1e-12);
  close(F.verticalStability({ ws: 130, b: 1300 }).uc, 1, 1e-12);       // SG = 1.1 exactly
});

test('Morison cycle: steady current and pure inertia limits', () => {
  const rho = 1025; const D = 0.5;
  const c = F.morisonCycle({ D, Ustar: 0, Vstar: 1.2, Tstar: 10 });
  close(c.FYmax, 0.5 * rho * D * 0.7 * 1.44, 1e-12);
  close(c.FZmax, 0.5 * rho * D * 0.9 * 1.44, 1e-12);
  // at 90 deg the velocity is zero and the force is pure inertia: rho*pi*D^2/4*CM*omega*U*
  const o = F.morisonCycle({ D, Ustar: 1, Vstar: 0, Tstar: 8, nPhase: 360 });
  close(Math.abs(o.FY[90]), rho * Math.PI / 4 * D * D * 3.29 * (2 * Math.PI / 8) * 1, 1e-9);
  assert.ok(Math.abs(o.FZ[90]) < 1e-9);
  close(o.K, 16, 1e-12);
});

test('absolute lateral stability: steady current reduces to ws = gamma (FD + mu FL)/mu', () => {
  const rho = 1025; const D = 0.5; const V = 1.2; const mu = 0.6;
  const FD = 0.5 * rho * D * 0.7 * V * V; const FL = 0.5 * rho * D * 0.9 * V * V;
  const wsReq = 1.1 * (FD + mu * FL) / mu;
  const r = F.absoluteLateralStability({ ws: wsReq, D, mu, Ustar: 0, Vstar: V, Tstar: 10, gammaSC: 1.1 });
  close(r.ucLateral, 1, 1e-12);
  assert.match(r.method, /morison/);
  close(F.requiredSubmergedWeight({ D, mu, Ustar: 0, Vstar: V, Tstar: 10, gammaSC: 1.1 }).ws, wsReq, 1e-9);
  // Fw scales the loads
  close(F.absoluteLateralStability({ ws: wsReq, D, mu, Ustar: 0, Vstar: V, Tstar: 10, gammaSC: 1.1, Fw: 1.2 }).ucLateral, 1.2, 1e-12);
});

test('absolute lateral stability with caller-supplied peak coefficients: hand calculation of Eq. 3.38/3.39', () => {
  // q = 0.5*1025*0.5*(1.0 + 0.5)^2 = 576.5625 N/m ; FY* = 1.2 q = 691.875 ; FZ* = 1.0 q = 576.5625
  // lateral: 1.3*(691.875 + 0.6*576.5625)/(0.6*3000) = 1.3*1037.8125/1800 = 0.74953 ; lift: 1.3*576.5625/3000 = 0.24984
  const r = F.absoluteLateralStability({ ws: 3000, D: 0.5, mu: 0.6, Ustar: 1.0, Vstar: 0.5, gammaSC: 1.3, method: 'coefficients', CY: 1.2, CZ: 1.0 });
  close(r.ucLateral, 0.74953, 1e-4);
  close(r.ucLift, 0.24984, 1e-4);
  assert.match(r.method, /coefficients/);
  // penetration lowers the loads and adds passive resistance
  const p = F.absoluteLateralStability({ ws: 3000, D: 0.5, mu: 0.6, Ustar: 1.0, Vstar: 0.5, gammaSC: 1.3, method: 'coefficients', CY: 1.2, CZ: 1.0,
    zp: 0.05, soil: 'sand', gammaSoilSub: 9000 });
  assert.ok(p.ucLateral < r.ucLateral && p.FR > 0);
  close(p.FY, 691.875 * 0.86, 1e-12);
  assert.throws(() => F.absoluteLateralStability({ ws: 3000, D: 0.5, mu: 0.6, Ustar: 1, Vstar: 0, gammaSC: 1.3, method: 'tables' }), RangeError);
  assert.throws(() => F.absoluteLateralStability({ ws: 3000, D: 0.5, mu: 0.6, Ustar: 1, Vstar: 0, Tstar: 10 }), /gammaSC/);
});

test('required concrete thickness: just satisfies the criteria and is minimal', () => {
  const pipe = { D: 0.3239, t: 0.0127, coatings: [{ thickness: 0.003, density: 1300 }], rhoContents: 0 };
  const stab = { mu: 0.6, Ustar: 1.0, Vstar: 0.4, Tstar: 11, gammaSC: 1.1 };
  const r = F.requiredConcreteThickness({ pipe, stability: stab });
  assert.ok(r.thickness > 0.01 && r.thickness < 0.2);
  assert.ok(r.check.uc <= 1 + 1e-9 && r.check.uc > 0.999);
  const thinner = pipeWeights({ ...pipe, coatings: [...pipe.coatings, { thickness: r.thickness - 0.002, density: 3040 }] });
  assert.ok(F.absoluteLateralStability({ ...stab, ws: thinner.submergedOperating, D: thinner.Dtotal }).uc > 1);
  assert.throws(() => F.requiredConcreteThickness({ pipe, stability: { ...stab, Ustar: 6 }, tMax: 0.05 }), RangeError);
});
