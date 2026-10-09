import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  takahashiImpulsive, godaPressure, caissonStability, sainflou, gravityWallStability, flotationCheck,
} from '../../src/engine/coastal/caisson.js';

const g = 9.80665;
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const abs = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b}`);

const base = { Hmax: 8, T: 11, h: 15, d: 9, hPrime: 11, hc: 4, hb: 16, B: 18, rho: 1030 };

test('Goda pressure coefficients: hand calculation from the published expressions', () => {
  const r = godaPressure(base);
  const k = (2 * Math.PI) / r.L;
  rel((2 * Math.PI / 11) ** 2, g * k * Math.tanh(k * 15), 1e-12, 'L solves the dispersion relation at depth h');
  const x = 2 * k * 15;
  rel(r.alpha1, 0.6 + 0.5 * (x / Math.sinh(x)) ** 2, 1e-12);
  rel(r.alpha2, Math.min(((16 - 9) / (3 * 16)) * (8 / 9) ** 2, (2 * 9) / 8), 1e-12);
  rel(r.alpha3, 1 - (11 / 15) * (1 - 1 / Math.cosh(k * 15)), 1e-12);
  rel(r.etaStar, 1.5 * 8, 1e-12); // 0.75 (1 + cos 0) Hmax
  const w = 1030 * g * 8;
  rel(r.p1, (r.alpha1 + r.alpha2) * w, 1e-12);
  rel(r.p2, r.p1 / Math.cosh(k * 15), 1e-12);
  rel(r.p3, r.alpha3 * r.p1, 1e-12);
  rel(r.p4, r.p1 * (1 - 4 / 12), 1e-12);
  rel(r.pu, r.alpha1 * r.alpha3 * w, 1e-12);
  rel(r.U, 0.5 * r.pu * 18, 1e-12);
  rel(r.Mu, (2 / 3) * r.U * 18, 1e-12);
  assert.equal(r.alphaI, 0);
});

test('Goda force and moment equal the numerical integral of the pressure diagram', () => {
  for (const o of [base, { ...base, hc: 20 }, { ...base, betaDeg: 30 }]) {
    const r = godaPressure(o);
    // p(z): linear p3 -> p1 from the base (z = −h') to still-water level, then linear to zero at η*,
    // truncated at the crest.
    const p = (z) => (z <= 0 ? r.p1 + ((r.p1 - r.p3) * z) / o.hPrime : Math.max(0, r.p1 * (1 - z / r.etaStar)));
    const top = Math.min(o.hc, r.etaStar);
    const n = 200000;
    const dz = (top + o.hPrime) / n;
    let F = 0;
    let M = 0;
    for (let i = 0; i < n; i++) {
      const z = -o.hPrime + (i + 0.5) * dz;
      F += p(z) * dz;
      M += p(z) * (z + o.hPrime) * dz;
    }
    rel(r.P, F, 1e-8, 'P');
    rel(r.Mp, M, 1e-8, 'Mp');
    rel(r.armP, M / F, 1e-8);
  }
  // Crest above η*: no truncation, p4 = 0.
  assert.equal(godaPressure({ ...base, hc: 20 }).p4, 0);
});

test('Goda limits: deep water α1 -> 0.6 and α3 -> 1 − h\'/h; shallow α1 -> 1.1; no mound α2 -> 0; obliquity reduces load', () => {
  const deep = godaPressure({ Hmax: 2, T: 4, h: 60, d: 60, hPrime: 20, hc: 3, hb: 60, B: 10 });
  rel(deep.alpha1, 0.6, 1e-6);
  rel(deep.alpha3, 1 - 20 / 60, 1e-6);
  assert.equal(deep.alpha2, 0);
  abs(deep.p2 / deep.p1, 0, 1e-6);
  const shallow = godaPressure({ Hmax: 0.5, T: 60, h: 2, d: 2, hPrime: 2, hc: 1, hb: 2, B: 5 });
  rel(shallow.alpha1, 1.1, 2e-3);
  rel(shallow.alpha3, 1, 2e-3); // long wave: uniform pressure over depth
  const ob = godaPressure({ ...base, betaDeg: 40 });
  const cb = Math.cos((40 * Math.PI) / 180);
  rel(ob.etaStar, 0.75 * (1 + cb) * 8, 1e-12);
  rel(ob.p1, 0.5 * (1 + cb) * (ob.alpha1 + ob.alpha2 * cb * cb) * 1030 * g * 8, 1e-12);
  assert.ok(ob.P < godaPressure(base).P);
  assert.throws(() => godaPressure({ ...base, d: 16 }), RangeError);
});

test('Takahashi impulsive coefficient: δ11 = δ22 = 0 gives αI1 = 1; αI0 capped at 2; enters as max(α2, αI)', () => {
  // BM/L = 0.12 and (h − d)/h = 0.6 → δ11 = δ22 = 0 → αI1 = cos 0 / cosh 0 = 1
  const t = takahashiImpulsive({ H: 6, d: 4, h: 10, BM: 12, L: 100 });
  rel(t.alphaI1, 1, 1e-12);
  rel(t.alphaI, 1.5, 1e-12);
  rel(takahashiImpulsive({ H: 10, d: 4, h: 10, BM: 12, L: 100 }).alphaI0, 2, 1e-12);
  // Hand calculation off the peak: BM/L = 0.05, (h − d)/h = 0.3
  const a = 0.05 - 0.12;
  const b = 0.3 - 0.6;
  const d11 = 0.93 * a + 0.36 * b; // −0.1731 → δ1 = 20 δ11
  const d22 = -0.36 * a + 0.93 * b; // −0.2538 → δ2 = 4.9 δ22
  const u = takahashiImpulsive({ H: 5, d: 7, h: 10, BM: 5, L: 100 });
  rel(u.alphaI1, Math.cos(4.9 * d22) / Math.cosh(20 * d11), 1e-12);
  // Low mound, no berm: negligible impulsive effect; high mound with a berm: αI governs.
  const low = godaPressure({ ...base, d: 14, hPrime: 14.5, BM: 0 });
  assert.ok(low.alphaI < 0.05);
  const high = godaPressure({ Hmax: 8, T: 11, h: 15, d: 6, hPrime: 8, hc: 4, hb: 16, B: 18, BM: 14 });
  assert.ok(high.alphaI > high.alpha2);
  rel(high.alphaStar, high.alphaI, 1e-12);
  const noBerm = godaPressure({ Hmax: 8, T: 11, h: 15, d: 6, hPrime: 8, hc: 4, hb: 16, B: 18 });
  assert.ok(high.P > noBerm.P);
});

test('caisson stability: safety factors by statics; heel pressure reproduces vertical and moment equilibrium', () => {
  const r = godaPressure(base);
  // Caisson 18 m wide, 15 m high (11 m submerged): W' = (2100·15 − 1030·11)·g·18 per metre
  const W = (2100 * 15 - 1030 * 11) * g * 18;
  const t = 9;
  const s = caissonStability({ P: r.P, Mp: r.Mp, U: r.U, Mu: r.Mu, W, t, B: 18, mu: 0.6 });
  rel(s.sfSliding, (0.6 * (W - r.U)) / r.P, 1e-12);
  rel(s.sfOverturning, (W * t - r.Mu) / r.Mp, 1e-12);
  rel(s.te, (W * t - r.Mu - r.Mp) / (W - r.U), 1e-12);
  // The linear contact-pressure block must carry We with its resultant at te from the heel.
  const b = s.contactWidth;
  const N = 0.5 * (s.pHeel + s.pToe) * b;
  const arm = (b * (s.pHeel + 2 * s.pToe)) / (3 * (s.pHeel + s.pToe));
  rel(N, s.We, 1e-12, 'vertical equilibrium');
  rel(arm, s.te, 1e-12, 'moment equilibrium');
  // Centred load: uniform pressure We/B. Load at B/3: triangular, 2We/B at the heel, zero at the toe.
  const c = caissonStability({ P: 1, Mp: 1e-9, U: 0, Mu: 0, W: 1000, t: 5, B: 10 });
  rel(c.pHeel, 100, 1e-9);
  rel(c.pToe, 100, 1e-9);
  const e = caissonStability({ P: 1, Mp: 1e-9, U: 0, Mu: 0, W: 1000, t: 10 / 3, B: 10 });
  rel(e.pHeel, 200, 1e-9);
  abs(e.pToe, 0, 1e-6);
  // Outside the middle third: partial contact, pe = 2We/(3te)
  const o = caissonStability({ P: 1, Mp: 1e-9, U: 0, Mu: 0, W: 1000, t: 2, B: 10 });
  rel(o.pHeel, 2000 / 6, 1e-9);
  rel(o.contactWidth, 6, 1e-9);
  assert.equal(o.fullContact, false);
});

test('Sainflou: h0 = (πH²/L) coth(kh), pressures, deep-water limit and numerical integration of the diagrams', () => {
  const H = 2;
  const T = 8;
  const h = 10;
  const rho = 1030;
  const r = sainflou({ H, T, h, rho });
  const k = (2 * Math.PI) / r.L;
  rel(r.h0, ((Math.PI * H * H) / r.L) / Math.tanh(k * h), 1e-12);
  rel(r.p2, (rho * g * H) / Math.cosh(k * h), 1e-12);
  rel(r.p1, ((r.p2 + rho * g * h) * (H + r.h0)) / (h + H + r.h0), 1e-12);
  const integrate = (p, a, b) => {
    const n = 100000;
    const dz = (b - a) / n;
    let F = 0;
    let M = 0;
    for (let i = 0; i < n; i++) {
      const z = a + (i + 0.5) * dz;
      F += p(z) * dz;
      M += p(z) * (z + h) * dz;
    }
    return [F, M];
  };
  const top = H + r.h0;
  const pc = (z) => (z <= 0 ? r.p1 + ((r.p1 - r.p2) * z) / h : r.p1 * (1 - z / top));
  const [Fc, Mc] = integrate(pc, -h, top);
  rel(r.crest.force, Fc, 1e-8);
  rel(r.crest.moment, Mc, 1e-8);
  const dt = H - r.h0;
  const ps = rho * g * dt;
  const pt = (z) => (z >= -dt ? -rho * g * z : ps + ((r.p2 - ps) * (-dt - z)) / (h - dt));
  const [Ft, Mt] = integrate(pt, -h, 0);
  rel(r.trough.force, Ft, 1e-8);
  rel(r.trough.moment, Mt, 1e-8);
  assert.ok(r.crest.force > r.trough.force);
  // Truncated wall: less force, and pTop is the pressure at the crest level.
  const tr = sainflou({ H, T, h, rho, hc: 1 });
  rel(tr.crest.pTop, r.p1 * (1 - 1 / top), 1e-12);
  assert.ok(tr.crest.force < r.crest.force);
  // Deep water: no pressure at the bed and h0 -> πH²/L.
  const deep = sainflou({ H: 1, T: 4, h: 80 });
  abs(deep.p2, 0, 1e-3);
  rel(deep.h0, Math.PI / deep.L, 1e-9);
});

test('gravity wall: rectangular block by hand, balanced water levels, drained base', () => {
  // Concrete block 4 m wide, 6 m high (ρ = 2400), water 5 m on one side and 2 m on the other, ρw = 1000.
  const rho = 1000;
  const W = 2400 * g * 4 * 6;
  const r = gravityWallStability({ W, xW: 2, B: 4, hHigh: 5, hLow: 2, mu: 0.5, rho });
  const Fh = 0.5 * rho * g * 25;
  const Fl = 0.5 * rho * g * 4;
  const U = rho * g * 3.5 * 4;
  rel(r.Fh, Fh, 1e-12);
  rel(r.Fl, Fl, 1e-12);
  rel(r.U, U, 1e-12);
  rel(r.armU, (4 * (10 + 2)) / (3 * 7), 1e-12); // centroid of the 5 m -> 2 m trapezoid from the low side
  rel(r.sfSliding, (0.5 * (W - U)) / (Fh - Fl), 1e-12);
  const mRes = W * 2 + (Fl * 2) / 3;
  const mOt = (Fh * 5) / 3 + U * r.armU;
  rel(r.sfOverturning, mRes / mOt, 1e-12);
  rel(r.resultantFromToe, (mRes - mOt) / (W - U), 1e-12);
  // Base pressures in the middle third follow N/B (1 ± 6e/B).
  assert.equal(r.middleThird, true);
  const N = W - U;
  rel(r.pMax, (N / 4) * (1 + (6 * Math.abs(r.eccentricity)) / 4), 1e-10);
  rel(r.pMin, (N / 4) * (1 - (6 * Math.abs(r.eccentricity)) / 4), 1e-10);
  // Equal water levels: no net thrust, uplift equals the buoyancy of the submerged part (Archimedes).
  const eq = gravityWallStability({ W, xW: 2, B: 4, hHigh: 3, hLow: 3, rho });
  assert.equal(eq.H, 0);
  assert.equal(eq.sfSliding, Infinity);
  rel(eq.U, rho * g * 3 * 4, 1e-12);
  rel(eq.sfFlotation, 2400 * 6 / (1000 * 3), 1e-12);
  // Halving the uplift (drained base) raises both safety factors.
  const dr = gravityWallStability({ W, xW: 2, B: 4, hHigh: 5, hLow: 2, mu: 0.5, rho, upliftFactor: 0.5 });
  assert.ok(dr.sfSliding > r.sfSliding && dr.sfOverturning > r.sfOverturning);
  assert.throws(() => gravityWallStability({ W, xW: 2, B: 4, hHigh: 2, hLow: 5 }), RangeError);
});

test('flotation: Archimedes buoyancy, safety factor and ballast to reach a target', () => {
  // Empty chamber 10 × 5 × 4 m fully below the water table, weight 1500 kN, fresh water.
  const f = flotationCheck({ weight: 1.5e6, displacedVolume: 200, rho: 1000 });
  rel(f.buoyancy, 1000 * g * 200, 1e-12);
  rel(f.sf, 1.5e6 / (1000 * g * 200), 1e-12);
  assert.equal(f.floats, true);
  const need = f.ballastForSf(1.1);
  rel(flotationCheck({ weight: 1.5e6 + need, displacedVolume: 200, rho: 1000 }).sf, 1.1, 1e-12);
  assert.equal(flotationCheck({ weight: 3e6, displacedVolume: 200, rho: 1000 }).floats, false);
  assert.equal(flotationCheck({ weight: 3e6, displacedVolume: 200, rho: 1000 }).ballastForSf(1.1), 0);
});
