import test from 'node:test';
import assert from 'node:assert/strict';
import * as B from '../../src/engine/pipeline/buckling.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${a} vs ${b}`);

const D = 0.3239; const t = 0.0127; const E = 207e9;
const EI = E * Math.PI / 64 * (D ** 4 - (D - 2 * t) ** 4);
const EA = E * Math.PI / 4 * (D ** 2 - (D - 2 * t) ** 2);
const w = 1200; const mu = 0.5;

test('Hobbs constants, mode 1 / vertical: re-derived from the beam-column solution', () => {
  // lift-off length L, y'''' + n^2 y'' = -w/EI with y = y' = y'' = 0 at the ends -> tan(nL/2) = nL/2
  let lo = Math.PI + 0.1; let hi = 1.5 * Math.PI - 1e-9;
  for (let i = 0; i < 100; i++) { const m = 0.5 * (lo + hi); if (Math.tan(m) - m < 0) lo = m; else hi = m; }
  const half = lo;                                        // nL/2 = 4.4934
  close(half, 4.4934, 1e-5);
  const k1 = (2 * half) ** 2;                             // P = k1 EI/L^2
  close(B.HOBBS_CONSTANTS[1].k1, k1, 1e-4);               // 80.76
  // amplitude y(0) = w L^4/EI [ (1 - 1/cos(nL/2))/k1^2 + 1/(8 k1) ] ; moment = w L^2 (1 - 1/cos(nL/2))/k1
  const c = Math.cos(half);
  close(B.HOBBS_CONSTANTS[1].k4, (1 - 1 / c) / k1 ** 2 + 1 / (8 * k1), 5e-4);     // 2.407e-3
  close(B.HOBBS_CONSTANTS[1].k5, (1 - 1 / c) / k1, 2e-4);                           // 0.06938
  const v = B.hobbsVertical({ L: 60, EI, EA, w, muAxial: mu });
  close(v.P, k1 * EI / 3600, 1e-4);
  close(v.amplitude, ((1 - 1 / c) / k1 ** 2 + 1 / (8 * k1)) * w * 60 ** 4 / EI, 1e-3);
  close(v.moment, (1 - 1 / c) / k1 * w * 3600, 2e-4);
  // lateral mode 1 has the same shape with phi*w as the lateral load, and k2 (lateral) = 4 k2 (vertical)
  close(B.HOBBS_CONSTANTS[1].k2, 4 * 1.597e-5, 5e-4);
});

test('Hobbs constants, mode infinity: re-derived in closed form', () => {
  // half-wave L with nL = 2 pi: y = q/(P n^2) (1 - cos nx) + q x (L - x)/(2P), q = phi*w
  // amplitude at L/2: q L^4/EI [ 2/(16 pi^4) + 1/(32 pi^2) ] ; moment 2 q/n^2 = q L^2/(2 pi^2)
  const c = B.HOBBS_CONSTANTS.inf;
  close(c.k4, 2 / (16 * Math.PI ** 4) + 1 / (32 * Math.PI ** 2), 2e-5);            // 4.4495e-3
  close(c.k5, 1 / (2 * Math.PI ** 2), 1e-4);                                          // 0.05066
  // axial shortening per unit length = (1/2L) int y'^2 dx = (q/P)^2 L^2 (2.5/(4 pi^2) + 1/12)/2 -> k2
  const k2 = (2.5 / (4 * Math.PI ** 2) + 1 / 12) / 2 / (16 * Math.PI ** 4);
  close(c.k2, k2, 2e-5);                                                              // 4.7050e-5
  // numerical confirmation of the shortening integral
  const L = 1; const n = 2 * Math.PI; let s = 0; const N = 20000;
  for (let i = 0; i < N; i++) { const x = (i + 0.5) / N; const yp = Math.sin(n * x) / n + (L - 2 * x) / 2; s += yp * yp / N; }
  close(s / 2 / (16 * Math.PI ** 4), c.k2, 2e-5);
  const r = B.hobbsLateral({ mode: 'inf', L: 40, EI, EA, w, muLateral: mu });
  close(r.P, 4 * Math.PI ** 2 * EI / 1600, 1e-12);
  close(r.P0 - r.P, k2 * EA * (mu * w / EI) ** 2 * 40 ** 6, 2e-5);
});

test('Hobbs mode infinity minimum: numerical minimum equals the closed form 2.29 (EA)^1/4 (EI)^1/4 fL^1/2', () => {
  const m = B.hobbsMinimum({ mode: 'inf', EI, EA, w, muLateral: mu });
  const c = B.infiniteModeForce({ EA, EI, fL: mu * w });
  close(c.constant, 2.289, 5e-4);
  close(m.P0min, c.Sinf, 1e-9);
  close(m.L, c.L, 1e-5);
  // stationarity: dP0/dL = 0 -> the friction term is one third of the buckle force term
  close(m.P0min, (4 / 3) * m.P, 1e-6);
  // it is a minimum
  for (const f of [0.9, 1.1]) assert.ok(B.hobbsLateral({ mode: 'inf', L: f * m.L, EI, EA, w, muLateral: mu }).P0 > m.P0min);
  assert.equal(B.infiniteModeForce({ EA, EI, fL: mu * w, S0: 2 * c.Sinf }).susceptible, true);
});

test('Hobbs lateral modes: minimum safe force ordering and critical temperature', () => {
  const p = [1, 2, 3, 4].map((mode) => B.hobbsMinimum({ mode, EI, EA, w, muLateral: mu, D }));
  // Hobbs (1984): the higher modes give slightly lower safe forces, mode 4 the lowest of the four
  assert.ok(p[3].P0min < p[2].P0min && p[2].P0min < p[1].P0min && p[1].P0min < p[0].P0min);
  assert.ok((p[0].P0min - p[3].P0min) / p[0].P0min < 0.1);
  for (const r of p) {
    close(r.deltaT, r.P0min / (EA * 11.7e-6), 1e-12);
    assert.ok(r.P0min > r.P && r.amplitude > 0);
    close(r.strain, r.moment * D / (2 * EI), 1e-12);
  }
  // pressure lowers the temperature needed: dT = (P0 - dp Ai (1 - 2 nu))/(EA alpha)
  const Ai = Math.PI / 4 * (D - 2 * t) ** 2;
  const q = B.hobbsMinimum({ mode: 3, EI, EA, w, muLateral: mu, deltaPi: 10e6, Ai });
  close(q.deltaT, (q.P0min - 10e6 * Ai * 0.4) / (EA * 11.7e-6), 1e-12);
  // explicit evaluation of the mode-2 formulas at L = 50 m
  const r = B.hobbsLateral({ mode: 2, L: 50, EI, EA, w, muLateral: mu });
  const P = 4 * Math.PI ** 2 * EI / 2500;
  close(r.P, P, 1e-12);
  close(r.P0, P + 1.0 * mu * w * 50 * (Math.sqrt(1 + 1.743e-4 * EA * mu * w * 50 ** 5 / EI ** 2) - 1), 1e-12);
  close(r.amplitude, 5.532e-3 * mu * w * 50 ** 4 / EI, 1e-12);
  close(r.moment, 0.1088 * mu * w * 2500, 1e-12);
  close(r.slipLength, (r.P0 - r.P) / (mu * w), 1e-12);
  // equal axial and lateral coefficients reproduce the single-coefficient form; lower axial friction lowers P0
  assert.ok(B.hobbsLateral({ mode: 2, L: 50, EI, EA, w, muLateral: mu, muAxial: 0.2 }).P0 < r.P0);
  const curve = B.hobbsCurve({ mode: 1, EI, EA, w, muLateral: mu, Lmin: 20, Lmax: 150, n: 50 });
  assert.equal(curve.L.length, 50);
  assert.ok(Math.min(...curve.P0) >= p[0].P0min * (1 - 1e-9));
});

test('Palmer et al. (1990) upheaval design curve: values and continuity at 4.49 and 8.06', () => {
  assert.equal(B.palmerPhiW(3), 0.0646);
  close(B.palmerPhiW(6), 5.68 / 36 - 88.35 / 1296, 1e-12);
  close(B.palmerPhiW(10), 9.6 / 100 - 343 / 1e4, 1e-12);
  // continuity (the published fits meet within about 1 %)
  close(B.palmerPhiW(4.49 - 1e-9), B.palmerPhiW(4.49 + 1e-9), 0.01);
  close(B.palmerPhiW(8.06 - 1e-9), B.palmerPhiW(8.06 + 1e-9), 0.01);
  // hand values at the joints: 5.68/20.1601 - 88.35/406.43 = 0.28174 - 0.21738 = 0.06436
  close(B.palmerPhiW(4.49), 0.06436, 1e-3);
  // 9.6/64.9636 - 343/4220.27 = 0.14777 - 0.08127 = 0.06650
  close(B.palmerPhiW(8.0600001), 0.06650, 1e-3);
});

test('upheaval download and cover: hand calculation', () => {
  // prop imperfection: touchdown distance (72 EI delta/w)^(1/4) each side
  const L = B.propImperfectionLength({ EI, delta: 0.3, w: 800 });
  close(L, 2 * (72 * EI * 0.3 / 800) ** 0.25, 1e-12);
  // verify the touchdown beam solution: y = a x^3 - w x^4/(24 EI), a = w l/(18 EI) gives y(l) = delta, y'(l) = 0
  const l = L / 2; const a = 800 * l / (18 * EI);
  close(a * l ** 3 - 800 * l ** 4 / (24 * EI), 0.3, 1e-12);
  assert.ok(Math.abs(3 * a * l * l - 800 * l ** 3 / (6 * EI)) < 1e-12);
  const P = 2e6;
  const r = B.upheavalDownload({ P, EI, delta: 0.3, wInstall: 800, wPipe: 1200, safetyFactor: 1.5, cover: { soil: 'sand', D: 0.33, gammaSub: 9000, f: 0.5 } });
  const phiL = L * Math.sqrt(P / EI);
  assert.ok(phiL > 8.06);
  const wReq = 1.5 * (9.6 / phiL ** 2 - 343 / phiL ** 4) * 0.3 * P * P / EI;
  close(r.wRequired, wReq, 1e-12);
  close(r.coverResistance, wReq - 1200, 1e-12);
  // cover H solves gamma' H D (1 + f H/D) = R
  close(9000 * r.coverHeight * 0.33 * (1 + 0.5 * r.coverHeight / 0.33), r.coverResistance, 1e-9);
  assert.equal(B.upheavalDownload({ P: 1e5, EI, delta: 0.1, L: 60, wPipe: 1200 }).stableWithoutCover, true);
});

test('buckle spacing: slip length, feed-in and spacing limits from axial friction', () => {
  const r = B.buckleSpacing({ P0: 2e6, Ppost: 0.5e6, f: 500, EA, PcritUnplanned: 1.5e6, PcritTrigger: 0.9e6 });
  close(r.slipLength, 3000, 1e-12);
  close(r.feedInIsolated, 1.5e6 ** 2 / (EA * 500), 1e-12);            // two sides of dP^2/(2 EA f)
  close(r.maxSpacing, 4000, 1e-12);
  close(r.minSpacing, 800, 1e-12);
  const s = B.buckleSpacing({ P0: 2e6, Ppost: 0.5e6, f: 500, EA, spacing: 2000 });
  assert.equal(s.interacting, true);
  close(s.forceAtAnchor, 0.5e6 + 500 * 1000, 1e-12);
  close(s.feedIn, (1.5e6 * 2000 - 500 * 2000 ** 2 / 4) / EA, 1e-12);
  // at spacing = 2 x slip length the interacting and isolated results coincide
  close(B.buckleSpacing({ P0: 2e6, Ppost: 0.5e6, f: 500, EA, spacing: 5999.999 }).feedIn, r.feedInIsolated, 1e-6);
});
