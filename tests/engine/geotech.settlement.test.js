import test from 'node:test';
import assert from 'node:assert/strict';
import * as s from '../../src/engine/geotech/settlement.js';

const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);

test('Boussinesq under the corner of a rectangle: published influence values (Newmark 1935; Das Table 10.x)', () => {
  close(s.boussinesqRectCorner(1, 1, 1), 0.1752, 3e-4);   // m = n = 1
  close(s.boussinesqRectCorner(2, 2, 1), 0.2325, 3e-4);   // m = n = 2
  close(s.boussinesqRectCorner(1, 2, 1), 0.1999, 3e-4);   // m = 1, n = 2
  close(s.boussinesqRectCorner(0.5, 0.5, 1), 0.0840, 6e-4); // m = n = 0.5
  assert.equal(s.boussinesqRectCorner(3, 4, 0), 0.25);
  // very large loaded area: quarter of the applied pressure
  close(s.boussinesqRectCorner(1e4, 1e4, 1), 0.25, 1e-6);
});

test('stress below a rectangle by superposition: centre = 4 corners of the quarter rectangle; symmetry outside', () => {
  close(s.stressBelowRect({ q: 100, B: 2, L: 2, z: 1 }), 400 * s.boussinesqRectCorner(1, 1, 1), 1e-12);
  // point at a corner
  close(s.stressBelowRect({ q: 100, B: 2, L: 4, z: 3, x: 1, y: 2 }), 100 * s.boussinesqRectCorner(2, 4, 3), 1e-12);
  // point outside the footing: add/subtract rectangles (hand superposition)
  const z = 2;
  const expect = 100 * (2 * s.boussinesqRectCorner(4, 1, z) - 2 * s.boussinesqRectCorner(2, 1, z));
  close(s.stressBelowRect({ q: 100, B: 2, L: 2, z, x: 3 }), expect, 1e-12);
});

test('circle centreline: 64.6 % of q at z = R; point load and 2:1 by hand', () => {
  close(s.stressBelowCircleCentre(1, 3, 3), 1 - 0.5 ** 1.5, 1e-14);
  close(s.stressPointLoad(1000, 0, 2), 3 * 1000 / (2 * Math.PI * 4), 1e-12);
  close(s.stress2to1({ q: 100, B: 2, L: 3, z: 2 }), 100 * 6 / 20, 1e-12);
  // integrating the point-load solution over a large circle reproduces the circle solution
  let sum = 0;
  const R = 2, z = 1.5, n = 4000;
  for (let i = 0; i < n; i++) { const r = (i + 0.5) * R / n; sum += s.stressPointLoad(1, r, z) * 2 * Math.PI * r * (R / n); }
  close(sum, s.stressBelowCircleCentre(1, R, z), 1e-6);
});

test('elastic settlement influence factors: square corner 0.561, centre 1.122, average 0.95; circle 1, 2/pi, pi/4', () => {
  close(s.influenceRectCorner(1), 0.561, 1e-3);
  const a = { q: 1, B: 1, L: 1, E: 1, nu: 0 };
  close(s.elasticSettlement({ ...a, point: 'corner' }).I, 0.561, 1e-3);
  close(s.elasticSettlement({ ...a, point: 'centre' }).I, 1.122, 1e-3);
  close(s.elasticSettlement({ ...a, point: 'average' }).I, 0.95, 6e-3);
  close(s.elasticSettlement({ ...a, L: 2, point: 'centre' }).I, 1.53, 5e-3);   // Timoshenko & Goodier: 1.53 for L/B = 2
  close(s.elasticSettlement({ ...a, L: 5, point: 'centre' }).I, 2.10, 5e-3);   // 2.10 for L/B = 5
  close(s.elasticSettlement({ ...a, shape: 'circle', point: 'rigid' }).I, Math.PI / 4, 1e-14);
  close(s.elasticSettlement({ ...a, shape: 'circle', point: 'average' }).I, 0.85, 2e-3);
  // hand: q = 100 kPa, B = 2 m circle, E = 20 MPa, nu = 0.3, rigid: s = 100e3·2·0.91·0.7854/20e6
  close(s.elasticSettlement({ q: 100e3, B: 2, E: 20e6, nu: 0.3, shape: 'circle', point: 'rigid' }).s, 100e3 * 2 * 0.91 * (Math.PI / 4) / 20e6, 1e-12);
});

test('Steinbrenner factors tend to the half-space factor for a very thick layer and vanish for H = 0', () => {
  const f = s.steinbrennerFactors(2, 1e5);
  close(f.F1, s.influenceRectCorner(2), 1e-4);
  assert.ok(f.F2 < 1e-3 * f.F1 * 100);
  assert.deepEqual(s.steinbrennerFactors(2, 0), { F1: 0, F2: 0 });
  // finite layer settles less than the half-space
  const half = s.elasticSettlement({ q: 1, B: 2, L: 4, E: 1, nu: 0.3, point: 'centre' }).s;
  const fin = s.elasticSettlement({ q: 1, B: 2, L: 4, E: 1, nu: 0.3, point: 'centre', H: 3 }).s;
  assert.ok(fin < half && fin > 0.3 * half);
});

test('rigid circular footing rotation (Borowicka): theta = 6 (1 - nu^2) M / (E D^3)', () => {
  const r = s.footingRotation({ M: 5e6, B: 10, E: 50e6, nu: 0.3, shape: 'circle' });
  close(r.theta, 6 * 0.91 * 5e6 / (50e6 * 1000), 1e-12);
});

test('Schmertmann (1978) strain-influence settlement: uniform Es, axisymmetric, hand integration of the Iz diagram', () => {
  const B = 2, dp = 150e3, sv0 = 20e3, gamma = 10e3, Es = 25e6;
  const izp = 0.5 + 0.1 * Math.sqrt(dp / (sv0 + gamma * B / 2));
  // area of the diagram: trapezoid 0.1 -> Izp over B/2, then triangle Izp -> 0 over 1.5 B
  const area = 0.5 * (0.1 + izp) * (B / 2) + 0.5 * izp * 1.5 * B;
  const C1 = 1 - 0.5 * sv0 / dp;
  const r = s.schmertmannSettlement({ B, dp, sigmaV0: sv0, gamma, layers: [{ zTop: 0, zBottom: 4, Es }], nSub: 400 });
  close(r.Izp, izp, 1e-12);
  close(r.C1, C1, 1e-12);
  close(r.s, C1 * dp * area / Es, 1e-5);
  // creep factor after 10 years: 1 + 0.2 log10(100) = 1.4
  close(s.schmertmannSettlement({ B, dp, sigmaV0: sv0, gamma, layers: [{ zTop: 0, zBottom: 4, Es }], timeYears: 10 }).C2, 1.4, 1e-12);
  // plane strain: diagram 0.2 -> Izp at B -> 0 at 4B
  const ps = s.schmertmannSettlement({ B, L: 40, dp, sigmaV0: sv0, gamma, layers: [{ zTop: 0, zBottom: 8, Es }], nSub: 800 });
  const izpP = 0.5 + 0.1 * Math.sqrt(dp / (sv0 + gamma * B));
  close(ps.s, C1 * dp * (0.5 * (0.2 + izpP) * B + 0.5 * izpP * 3 * B) / Es, 1e-5);
});

test('Janbu tangent-modulus strain: linear case j = 1 and logarithmic case j = 0', () => {
  close(s.janbuStrain({ sigma0: 100e3, sigma1: 150e3, m: 200, j: 1 }), 50e3 / (200 * 100e3), 1e-12);
  close(s.janbuStrain({ sigma0: 100e3, sigma1: 200e3, m: 20, j: 0 }), Math.log(2) / 20, 1e-12);
  close(s.janbuStrain({ sigma0: 100e3, sigma1: 400e3, m: 300, j: 0.5 }), (2 - 1) / (300 * 0.5), 1e-12);
});

test('oedometer strain: NC, OC below and across the preconsolidation pressure (hand)', () => {
  close(s.oedometerStrain({ sigma0: 100e3, dSigma: 100e3, Cc: 0.3, Cr: 0.05, e0: 1 }), 0.15 * Math.log10(2), 1e-12);
  close(s.oedometerStrain({ sigma0: 100e3, dSigma: 50e3, sigmaP: 200e3, Cc: 0.3, Cr: 0.05, e0: 1 }), 0.025 * Math.log10(1.5), 1e-12);
  close(s.oedometerStrain({ sigma0: 100e3, dSigma: 300e3, sigmaP: 200e3, Cc: 0.3, Cr: 0.05, e0: 1 }), 0.025 * Math.log10(2) + 0.15 * Math.log10(2), 1e-12);
});

test('layer-summation consolidation settlement: wide load on a thin NC layer equals the one-dimensional hand value', () => {
  // B = L = 1000 m so that the stress increment is q; layer 2 m thick at 9 to 11 m depth
  const r = s.consolidationSettlement({
    q: 50e3, B: 1000, L: 1000,
    layers: [{ thickness: 9, gamma: 8e3, Cc: 0, Cr: 0, e0: 1 }, { thickness: 2, gamma: 8e3, Cc: 0.4, Cr: 0.04, e0: 1.2 }], nSub: 1,
  });
  const s0 = 80e3;
  close(r.s, (0.4 / 2.2) * 2 * Math.log10((s0 + 50e3) / s0), 1e-5);
  // Boussinesq vs 2:1 are of the same order below a footing
  const a = s.consolidationSettlement({ q: 100e3, B: 3, L: 3, sigmaV0Base: 10e3, layers: [{ thickness: 6, gamma: 8e3, Cc: 0.3, Cr: 0.03, e0: 1 }] });
  const b2 = s.consolidationSettlement({ q: 100e3, B: 3, L: 3, sigmaV0Base: 10e3, stressMethod: '2:1', layers: [{ thickness: 6, gamma: 8e3, Cc: 0.3, Cr: 0.03, e0: 1 }] });
  assert.ok(a.s > 0 && Math.abs(a.s - b2.s) / a.s < 0.25);
});

test('Terzaghi consolidation: U = 50 % at Tv = 0.197 and 90 % at Tv = 0.848; inverse; approximations', () => {
  close(s.consolidationU(0.197), 0.5, 2e-3);
  close(s.consolidationU(0.848), 0.9, 2e-4);
  close(s.consolidationTv(0.5), 0.197, 2e-3);
  close(s.consolidationTv(0.9), 0.848, 2e-4);
  // classical approximations: Tv = (pi/4) U^2 (U < 0.6), Tv = 1.781 - 0.933 log10(100 - U%) (U > 0.6)
  close(s.consolidationTv(0.3), (Math.PI / 4) * 0.09, 1e-4);
  close(s.consolidationTv(0.8), 1.781 - 0.933 * Math.log10(20), 2e-3);
  for (const U of [0.01, 0.1, 0.45, 0.75, 0.99]) close(s.consolidationU(s.consolidationTv(U)), U, 1e-9);
  assert.equal(s.consolidationU(0), 0);
});

test('excess pore pressure isochrone integrates to 1 - U', () => {
  const Tv = 0.3;
  let sum = 0;
  const n = 2000;
  for (let i = 0; i < n; i++) sum += s.excessPorePressureRatio((i + 0.5) / n, Tv) / n;
  close(1 - sum, s.consolidationU(Tv), 1e-6);
});

test('Barron radial consolidation: F(n = 20) = 2.2539 by hand; Hansbo ideal-drain form; Carrillo combination', () => {
  const r = s.radialConsolidation({ ch: 1e-7, t: 3.0e7, de: 2, dw: 0.1 });
  close(r.F, (400 / 399) * Math.log(20) - 1199 / 1600, 1e-12);
  close(r.F, 2.2539, 1e-4);
  close(r.Ur, 1 - Math.exp(-8 * (1e-7 * 3e7 / 4) / r.F), 1e-12);
  const h = s.radialConsolidation({ ch: 1e-7, t: 3.0e7, de: 2, dw: 0.1, s: 2, khOverKs: 3 });
  close(h.F, Math.log(10) + 3 * Math.log(2) - 0.75, 1e-12);
  close(s.combinedConsolidation(0.3, 0.6), 0.72, 1e-12);
});

test('secondary compression and time-settlement curve', () => {
  close(s.secondaryCompression({ Calpha: 0.02, e0: 1, H: 5, t1: 1, t2: 100 }), 0.02 / 2 * 5 * 2, 1e-12);
  const c = s.timeSettlementCurve({ times: [0, 1e7, 1e12], si: 0.01, sc: 0.2, cv: 1e-7, Hdr: 3 });
  close(c.s[0], 0.01, 1e-12);
  close(c.s[1], 0.01 + 0.2 * s.consolidationU(1e-7 * 1e7 / 9), 1e-12);
  close(c.s[2], 0.21, 1e-9);
});

test('finite-volume consolidation reproduces Terzaghi U(Tv) for a single layer, a split layer and double drainage', () => {
  const k = 1e-9, mv = 5e-7, H = 4, gw = 9810;
  const cv = k / (mv * gw);
  const t50 = 0.197 * H * H / cv;
  const r = s.consolidationFD({ layers: [{ thickness: H, k, mv }], load: () => 100e3, tEnd: t50, nSteps: 800, nCells: 80, theta: 0.5 });
  close(r.U[r.U.length - 1], s.consolidationU(0.197), 2e-3);
  close(r.sFinal, mv * 100e3 * H, 1e-12);
  const r2 = s.consolidationFD({ layers: [{ thickness: 1, k, mv }, { thickness: 3, k, mv }], load: () => 100e3, tEnd: t50, nSteps: 800, nCells: 40, theta: 0.5 });
  close(r2.U[r2.U.length - 1], s.consolidationU(0.197), 2e-3);
  // double drainage is four times faster
  const r3 = s.consolidationFD({ layers: [{ thickness: H, k, mv }], load: () => 100e3, tEnd: t50 / 4, nSteps: 800, nCells: 80, bottom: 'drained', theta: 0.5 });
  close(r3.U[r3.U.length - 1], s.consolidationU(0.197), 2e-3);
});

test('finite-volume consolidation: a ramp load lags the instantaneous load; no settlement at the instant of loading', () => {
  const k = 1e-9, mv = 5e-7;
  const single = s.consolidationFD({ layers: [{ thickness: 4, k, mv }], load: () => 50e3, tEnd: 5e8, nSteps: 300 });
  const ramp = s.consolidationFD({ layers: [{ thickness: 4, k, mv }], load: (t) => 50e3 * Math.min(t / 2.5e8, 1), tEnd: 5e8, nSteps: 300 });
  const sInst = single.settlement[single.settlement.length - 1];
  const sRamp = ramp.settlement[ramp.settlement.length - 1];
  assert.ok(sRamp < sInst && sRamp > 0.5 * sInst);
  // undrained instant response: no settlement at t = 0 for an instantaneous load
  assert.equal(single.settlement[0], 0);
});

test('Biot column: incompressible constituents reduce to Terzaghi; compressible fluid gives the undrained loading efficiency', () => {
  const K = 20e6, G = 12e6, k = 1e-8, h = 5, q = 100e3;
  const Kv = K + 4 * G / 3;
  const t = 0.197 * h * h / (k * Kv / 9810);
  const r = s.biotColumn({ q, h, K, G, k, t });
  close(r.p0, q, 1e-12);
  close(r.U, 0.5, 2e-3);
  close(r.settlementFinal, q * h / Kv, 1e-12);
  // base pore pressure at Tv = 0.197, first two series terms written out:
  // (4/pi) exp(-pi^2 Tv/4) - (4/(3 pi)) exp(-9 pi^2 Tv/4) = 0.7831 - 0.0054 = 0.7777
  const Tv = 0.197;
  const series = (4 / Math.PI) * Math.exp(-(Math.PI ** 2) * Tv / 4) - (4 / (3 * Math.PI)) * Math.exp(-9 * (Math.PI ** 2) * Tv / 4);
  close(series, 0.7777, 2e-4);
  close(r.p / q, series, 1e-4);
  const M = 4e8;
  const c = s.biotColumn({ q, h, K, G, k, M, alpha: 0.9, t: 0 });
  close(c.p0, q * 0.9 * M / (Kv + 0.81 * M), 1e-12);
  close(c.settlementImmediate, q * h / (Kv + 0.81 * M), 1e-12);
  // Biot diffusion check against the finite-volume solver (storage 1/M* = 1/Kv for alpha = 1, M = inf)
  const fd = s.consolidationFD({ layers: [{ thickness: h, k, mv: 1 / Kv }], load: () => q, tEnd: t, nSteps: 600, nCells: 60, theta: 0.5 });
  close(fd.u[fd.u.length - 1] / q, r.p / q, 5e-3);
});
