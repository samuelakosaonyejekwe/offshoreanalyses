import test from 'node:test';
import assert from 'node:assert/strict';
import * as m from '../../src/engine/geotech/soilmodels.js';

const deg = Math.PI / 180;
const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);
const mcc = { M: 1.2, lambda: 0.2, kappa: 0.04, N: 3.0 };

test('Modified Cam-Clay undrained triaxial (NC): reaches the critical-state line with su/p0 = (M/2) 2^(-Lambda) (Wood 1990)', () => {
  const r = m.mccTriaxial({ ...mcc, p0: 200e3, drainage: 'undrained', epsAMax: 0.3 });
  const ratio = (1.2 / 2) * 0.5 ** 0.8;
  close(m.mccUndrainedStrengthRatio(mcc), ratio, 1e-12);
  close(r.su / 200e3, ratio, 1e-3);
  close(r.eta, 1.2, 1e-4);
  close(r.pf, 200e3 * 0.5 ** 0.8, 1e-3);
  // constant volume: specific volume unchanged; positive excess pore pressure u = p0 + q/3 - p'
  close(r.vf, r.v[0], 1e-10);
  close(r.u[r.u.length - 1], 200e3 + r.qf / 3 - r.pf, 1e-12);
  assert.ok(r.u[r.u.length - 1] > 0);
});

test('Modified Cam-Clay undrained triaxial (OCR = 4): su/p0 = (M/2)(R/2)^Lambda; elastic path at constant p\' before yield; negative pore pressure', () => {
  const r = m.mccTriaxial({ ...mcc, p0: 100e3, ocr: 4, drainage: 'undrained', epsAMax: 0.4 });
  close(r.su / 100e3, 0.6 * 2 ** 0.8, 1e-3);
  close(r.eta, 1.2, 1e-4);
  close(r.p[5], 100e3, 1e-9);
  assert.ok(r.u[r.u.length - 1] < 0);
  // yield where the constant-p' path meets the initial ellipse: q = M sqrt(p (pc - p)) = 1.2 sqrt(100·300) kPa
  const qy = 1.2 * Math.sqrt(100e3 * 300e3);
  const i = r.p.findIndex((p) => Math.abs(p - 100e3) > 1e-4 * 100e3);
  assert.ok(Math.abs(r.q[i - 1] - qy) / qy < 0.03, `yield at q = ${r.q[i - 1]}`);
});

test('Modified Cam-Clay drained triaxial (NC): stress path q = 3 (p\' - p0), critical state at p\' = 3 p0/(3 - M), v = Gamma - lambda ln p\'', () => {
  const r = m.mccTriaxial({ ...mcc, p0: 200e3, drainage: 'drained', epsAMax: 0.8, nSteps: 800 });
  for (const i of [100, 400, 800]) close(r.q[i], 3 * (r.p[i] - 200e3), 1e-6);
  const pf = 3 * 200e3 / (3 - 1.2);
  close(r.pf, pf, 3e-3);
  close(r.eta, 1.2, 5e-3);
  const Gamma = 3.0 - (0.2 - 0.04) * Math.log(2);
  close(r.vf, Gamma - 0.2 * Math.log(r.pf / 1000), 1e-3);
  // initial state on the isotropic normal compression line
  close(r.v[0], 3.0 - 0.2 * Math.log(200), 1e-12);
  // volumetric strain consistent with the change in specific volume
  close(Math.exp(-r.epsV[800]) * r.v[0], r.vf, 2e-3);
});

test('Modified Cam-Clay integration converges with step refinement (undrained strength)', () => {
  const coarse = m.mccTriaxial({ ...mcc, p0: 150e3, ocr: 2, drainage: 'undrained', epsAMax: 0.3, nSteps: 100, subSteps: 2 });
  const fine = m.mccTriaxial({ ...mcc, p0: 150e3, ocr: 2, drainage: 'undrained', epsAMax: 0.3, nSteps: 800, subSteps: 10 });
  const exact = 150e3 * m.mccUndrainedStrengthRatio({ ...mcc, ocr: 2 });
  assert.ok(Math.abs(fine.su - exact) <= Math.abs(coarse.su - exact) + 1e-6 * exact);
  close(fine.su, exact, 1e-3);
});

test('Mohr-Coulomb drained triaxial compression: sigma1 = sigma3 Kp + 2 c sqrt(Kp); initial stiffness E', () => {
  const r = m.frictionalTriaxial({ E: 50e6, nu: 0.3, phi: 30 * deg, c: 10e3, p0: 100e3, drainage: 'drained' });
  close(r.qf + 100e3, 100e3 * 3 + 2 * 10e3 * Math.sqrt(3), 1e-9);
  close(r.q[1] / r.epsA[1], 50e6, 1e-6);
  // Poisson's ratio from the elastic volumetric strain: eps_v = (1 - 2 nu) eps_a
  close(r.epsV[1] / r.epsA[1], 0.4, 1e-6);
});

test('Mohr-Coulomb undrained (psi = 0): p\' constant after yield, q = M p0 + kq; dilation raises undrained strength', () => {
  const p = m.mohrCoulombTriaxialParams({ phi: 30 * deg, c: 10e3 });
  close(p.M, 1.2, 1e-12);
  close(p.kq, 6 * 10e3 * Math.cos(30 * deg) / 2.5, 1e-12);
  const r = m.frictionalTriaxial({ E: 50e6, nu: 0.3, phi: 30 * deg, c: 10e3, p0: 100e3, drainage: 'undrained' });
  close(r.pf, 100e3, 1e-9);
  close(r.qf, 1.2 * 100e3 + p.kq, 1e-9);
  const d = m.frictionalTriaxial({ E: 50e6, nu: 0.3, phi: 30 * deg, c: 10e3, psi: 10 * deg, p0: 100e3, drainage: 'undrained' });
  assert.ok(d.qf > r.qf && d.u[d.u.length - 1] < r.u[r.u.length - 1]);
  // drained test with dilation: plastic volumetric strain rate -Md per unit plastic shear strain
  const dd = m.frictionalTriaxial({ E: 50e6, nu: 0.3, phi: 30 * deg, psi: 10 * deg, p0: 100e3, drainage: 'drained', epsAMax: 0.05, nSteps: 500 });
  const Md = m.mohrCoulombTriaxialParams({ phi: 10 * deg }).M;
  const n = dd.epsA.length - 1;
  close((dd.epsV[n] - dd.epsV[n - 100]) / (dd.epsQ[n] - dd.epsQ[n - 100]), -Md, 1e-6);
});

test('Drucker-Prager cones matched to Mohr-Coulomb: compression cone reproduces the triaxial compression strength; plane-strain constants', () => {
  const dp = m.druckerPragerParams({ phi: 30 * deg, c: 10e3 });
  const mc = m.mohrCoulombTriaxialParams({ phi: 30 * deg, c: 10e3 });
  close(dp.M, mc.M, 1e-12);
  close(dp.kq, mc.kq, 1e-12);
  close(dp.alpha, 2 * 0.5 / (Math.sqrt(3) * 2.5), 1e-12);
  const ext = m.druckerPragerParams({ phi: 30 * deg, c: 10e3, match: 'extension' });
  close(ext.M, m.mohrCoulombTriaxialParams({ phi: 30 * deg, mode: 'extension' }).M, 1e-12);
  const ps = m.druckerPragerParams({ phi: 30 * deg, c: 10e3, match: 'plane-strain' });
  close(ps.alpha, Math.tan(30 * deg) / Math.sqrt(9 + 12 / 3), 1e-12);
  close(ps.k, 3 * 10e3 / Math.sqrt(13), 1e-12);
  const r = m.frictionalTriaxial({ model: 'drucker-prager', E: 50e6, nu: 0.3, phi: 30 * deg, c: 10e3, p0: 100e3, drainage: 'drained' });
  close(r.qf, 3 * (mc.M * 100e3 + mc.kq) / (3 - mc.M), 1e-9);
  // tension cut-off limits the tensile mean stress in an undrained dilative test
  const t = m.frictionalTriaxial({ E: 50e6, nu: 0.3, phi: 30 * deg, c: 10e3, p0: 0, drainage: 'undrained', tensionCutoff: 5e3 });
  assert.ok(Math.min(...t.p) >= -5e3 - 1e-6);
});

test('simple shear: drained strength c + sigma tan(phi) with dilation rate -tan(psi); undrained dilative hardening slope by hand', () => {
  const d = m.simpleShearMohrCoulomb({ G: 10e6, Moed: 30e6, phi: 30 * deg, c: 5e3, psi: 5 * deg, sigmaN0: 100e3, drainage: 'drained', gammaMax: 0.05 });
  const n = d.tau.length - 1;
  close(d.tau[n], 5e3 + 100e3 * Math.tan(30 * deg), 1e-12);
  const gy = d.tau[n] / 10e6;
  close(d.epsN[n], -Math.tan(5 * deg) * (0.05 - gy), 1e-9);
  const u = m.simpleShearMohrCoulomb({ G: 10e6, Moed: 30e6, phi: 30 * deg, c: 5e3, psi: 5 * deg, sigmaN0: 100e3, drainage: 'undrained', gammaMax: 0.05 });
  const tf = Math.tan(30 * deg), td = Math.tan(5 * deg);
  const slope = 10e6 * 30e6 * tf * td / (10e6 + 30e6 * tf * td);
  close((u.tau[n] - u.tau[n - 50]) / (u.gamma[n] - u.gamma[n - 50]), slope, 1e-9);
  assert.ok(u.u[n] < 0);
});

test('Duncan-Chang hyperbola: initial modulus, failure stress, asymptote and tangent modulus by hand', () => {
  const o = { K: 600, n: 0.5, Rf: 0.8, phi: 35 * deg, sigma3: 200e3 };
  const r = m.duncanChang({ ...o, epsA: 0.01, Kur: 1200 });
  const Ei = 600 * 101325 * Math.sqrt(200e3 / 101325);
  const qf = 2 * 200e3 * Math.sin(35 * deg) / (1 - Math.sin(35 * deg));
  close(r.Ei, Ei, 1e-12);
  close(r.qf, qf, 1e-12);
  close(r.q, 0.01 / (1 / Ei + 0.8 * 0.01 / qf), 1e-12);
  close(r.Et, Ei * (1 - 0.8 * r.q / qf) ** 2, 1e-12);
  close(r.Eur, 2 * Ei, 1e-12);
  // tangent modulus equals the derivative of the hyperbola
  const h = 1e-7;
  const q1 = m.duncanChang({ ...o, epsA: 0.01 + h }).q;
  const q0 = m.duncanChang({ ...o, epsA: 0.01 - h }).q;
  close((q1 - q0) / (2 * h), r.Et, 1e-6);
});

test('Janbu modulus and Hardin Gmax by hand', () => {
  close(m.janbuModulus({ m: 300, j: 0.5, sigma: 400e3 }), 300 * 100e3 * 2, 1e-12);
  close(m.gmaxHardin({ e: 0.7, sigmaM: 100e3 }), 625 / (0.3 + 0.7 * 0.49) * Math.sqrt(101325 * 100e3), 1e-12);
});

test('Hardin-Drnevich: G/Gmax = 0.5 and D = Dmax/2 at the reference strain', () => {
  const r = m.hardinDrnevich(1e-3, { gammaRef: 1e-3, Dmax: 0.28 });
  close(r.ratio, 0.5, 1e-14);
  close(r.damping, 0.14, 1e-14);
});

test('Darendeli (2001): reference strain 0.0352 % at PI = 0, OCR = 1, 1 atm; G/Gmax = 0.5 there; Masing damping 14.5 % for a = 1 hyperbola', () => {
  const r = m.darendeli(0.000352, { sigmaM: 101325 });
  close(r.gammaRef, 0.000352, 1e-12);
  close(r.ratio, 0.5, 1e-12);
  // PI = 30, OCR = 2, 4 atm: gamma_r = (0.0352 + 0.001·30·2^0.3246) 4^0.3483 %
  close(m.darendeli(1e-4, { sigmaM: 4 * 101325, PI: 30, ocr: 2 }).gammaRef, (0.0352 + 0.03 * 2 ** 0.3246) * 4 ** 0.3483 / 100, 1e-12);
  // small-strain damping: Dmin = 0.8005 % at 1 atm, 1 Hz
  close(m.darendeli(1e-9, { sigmaM: 101325 }).damping, 0.008005, 1e-4);
  // Masing damping of a pure hyperbola at gamma = gamma_r: (100/pi)(8 (1 - ln 2) - 2) = 14.48 %
  close((100 / Math.PI) * (8 * (1 - Math.log(2)) - 2), 14.48, 1e-3);
  // monotonic trends
  let prev = m.darendeli(1e-6, { sigmaM: 101325 });
  for (const g of [1e-5, 1e-4, 1e-3, 1e-2]) {
    const cur = m.darendeli(g, { sigmaM: 101325 });
    assert.ok(cur.ratio < prev.ratio && cur.damping > prev.damping);
    prev = cur;
  }
  assert.ok(m.darendeli(1e-3, { sigmaM: 101325, PI: 50 }).ratio > m.darendeli(1e-3, { sigmaM: 101325 }).ratio);
});

test('Ishibashi & Zhang (1993): K and m - m0 by hand at gamma = 1e-3, PI = 0, 100 kPa; plasticity makes the curve more linear', () => {
  const g = 1e-3;
  const K = 0.5 * (1 + Math.tanh(Math.log((0.000102 / g) ** 0.492)));
  const mm = 0.272 * (1 - Math.tanh(Math.log((0.000556 / g) ** 0.4)));
  const r = m.ishibashiZhang(g, { sigmaM: 100e3 });
  close(r.ratio, K * 100 ** mm, 1e-12);
  close(r.damping, 0.333 * (0.586 * r.ratio ** 2 - 1.547 * r.ratio + 1), 1e-12);
  assert.ok(m.ishibashiZhang(g, { sigmaM: 100e3, PI: 50 }).ratio > r.ratio);
  assert.ok(m.ishibashiZhang(1e-6, { sigmaM: 100e3 }).ratio > 0.97);
});

test('modulus-reduction dispatcher, tabulated curves and Idriss degradation index', () => {
  close(m.modulusReduction(1e-3, { model: 'hardin-drnevich', gammaRef: 1e-3 }).ratio, 0.5, 1e-14);
  const tab = { model: 'table', strain: [1e-5, 1e-4, 1e-3], ratio: [1, 0.8, 0.4], damping: [0.01, 0.03, 0.1] };
  close(m.modulusReduction(Math.sqrt(1e-9), tab).ratio, 0.9, 1e-12);
  close(m.modulusReduction(1, tab).damping, 0.1, 1e-14);
  assert.deepEqual(m.modulusReduction(0.01, { model: 'linear', damping: 0.05 }), { ratio: 1, damping: 0.05 });
  close(m.cyclicDegradationIndex(100, 0.1), 100 ** -0.1, 1e-14);
  assert.throws(() => m.modulusReduction(1e-3, { model: 'nope' }), RangeError);
});
