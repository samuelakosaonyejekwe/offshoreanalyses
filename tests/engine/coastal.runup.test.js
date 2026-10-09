import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  iribarren, huntRunup, ROUGHNESS_FACTORS, obliqueFactor, bermFactor, roughnessFactorSurging, runupEurOtop,
  overtoppingSlope, overtoppingSteepSlope, overtoppingVertical, overtoppingRubbleMound, TOLERABLE_OVERTOPPING,
  classifyOvertopping, overtoppingVolumes, seaLevelRiseSensitivity,
} from '../../src/engine/coastal/runup.js';

const g = 9.80665;
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('Iribarren number and Hunt run-up, hand calculation', () => {
  // H = 2 m, T = 8 s, slope 1:4: L0 = g·64/(2π) = 99.88 m, ξ = 0.25/√(2/99.88) = 1.767
  const r = iribarren({ tanAlpha: 0.25, H: 2, T: 8 });
  rel(r.L0, (g * 64) / (2 * Math.PI), 1e-12);
  rel(r.xi, 0.25 / Math.sqrt(2 / r.L0), 1e-12);
  assert.ok(Math.abs(r.xi - 1.767) < 1e-3);
  assert.equal(r.breakerType, 'plunging');
  const h = huntRunup({ H: 2, T: 8, tanAlpha: 0.25 });
  rel(h.R, 2 * r.xi, 1e-12);
  assert.equal(h.valid, true);
  assert.equal(huntRunup({ H: 0.5, T: 12, tanAlpha: 0.5 }).valid, false);
  assert.throws(() => iribarren({ tanAlpha: 0, H: 1, T: 5 }), RangeError);
});

test('EurOtop (2018) Eq. 5.1/5.2 run-up: breaking branch, surging maximum, design values', () => {
  // Hm0 = 2 m, Tm−1,0 = 6 s, 1:4 smooth: ξ = 0.25/√(2/56.18) = 1.325
  const xi = 0.25 / Math.sqrt(2 / ((g * 36) / (2 * Math.PI)));
  const m = runupEurOtop({ Hm0: 2, Tm10: 6, tanAlpha: 0.25 });
  rel(m.xi, xi, 1e-12);
  rel(m.ratio, 1.65 * xi, 1e-12); // 2.19 < 4 − 1.5/√ξ = 2.70
  assert.equal(m.regime, 'breaking');
  rel(runupEurOtop({ Hm0: 2, Tm10: 6, tanAlpha: 0.25, design: true }).ratio, 1.75 * xi, 1e-12);
  // Steep smooth slope 1:1.5, long waves: the maximum governs.
  const s = runupEurOtop({ Hm0: 1, Tm10: 10, tanAlpha: 1 / 1.5 });
  rel(s.ratio, 4 - 1.5 / Math.sqrt(s.xi), 1e-12);
  assert.equal(s.regime, 'surging');
  rel(runupEurOtop({ Hm0: 1, Tm10: 10, tanAlpha: 1 / 1.5, design: true }).ratio, 1.07 * (4 - 1.5 / Math.sqrt(s.xi)), 1e-12);
  // Influence factors scale the breaking branch linearly.
  const f = runupEurOtop({ Hm0: 2, Tm10: 6, tanAlpha: 0.25, gammaF: 0.8, gammaBeta: 0.9, gammaB: 0.7 });
  rel(f.ratio, 1.65 * 0.8 * 0.9 * 0.7 * xi, 1e-12);
  // Rough slope in surging conditions uses γf,surging (Eq. 6.2).
  const rr = runupEurOtop({ Hm0: 1, Tm10: 10, tanAlpha: 0.5, gammaF: 0.4 });
  const gfs = 0.4 + ((rr.xi - 1.8) * 0.6) / 8.2;
  rel(roughnessFactorSurging({ gammaF: 0.4, xi: rr.xi }), gfs, 1e-12);
  rel(rr.ratio, gfs * (4 - 1.5 / Math.sqrt(rr.xi)), 1e-12);
  assert.equal(roughnessFactorSurging({ gammaF: 0.4, xi: 1.2 }), 0.4);
  assert.equal(roughnessFactorSurging({ gammaF: 0.4, xi: 12 }), 1);
});

test('influence factors: obliquity and berm', () => {
  rel(obliqueFactor({ betaDeg: 30, kind: 'runup' }), 1 - 0.0022 * 30, 1e-12);
  rel(obliqueFactor({ betaDeg: -30 }), 1 - 0.0033 * 30, 1e-12);
  rel(obliqueFactor({ betaDeg: 85, kind: 'rubble' }), 1 - 0.0063 * 80, 1e-12);
  // Berm at still-water level is most effective: rdb = 0, γb = 1 − B/Lberm, floored at 0.6.
  rel(bermFactor({ B: 4, Lberm: 20, db: 0, Hm0: 2, Ru2: 4 }), 0.8, 1e-12);
  assert.equal(bermFactor({ B: 15, Lberm: 20, db: 0, Hm0: 2, Ru2: 4 }), 0.6);
  // A berm deeper than 2 Hm0 or higher than Ru2% has no effect.
  assert.equal(bermFactor({ B: 4, Lberm: 20, db: 4.5, Hm0: 2, Ru2: 4 }), 1);
  assert.equal(bermFactor({ B: 4, Lberm: 20, db: -4.5, Hm0: 2, Ru2: 4 }), 1);
  rel(bermFactor({ B: 4, Lberm: 20, db: 2, Hm0: 2, Ru2: 4 }), 1 - 0.2 * (1 - 0.5), 1e-12);
  assert.equal(ROUGHNESS_FACTORS['rock-2layer-permeable'], 0.4);
  assert.equal(ROUGHNESS_FACTORS.tetrapod, 0.38);
});

test('EurOtop (2018) Eq. 5.10–5.13 overtopping on slopes: hand calculation of both branches', () => {
  const Hm0 = 2;
  const Tm10 = 6;
  const tanA = 1 / 3;
  const Rc = 3;
  const xi = tanA / Math.sqrt(Hm0 / ((g * 36) / (2 * Math.PI))); // 1.767
  const sc = Math.sqrt(g * 8);
  const qb = (0.023 / Math.sqrt(tanA)) * xi * Math.exp(-(((2.7 * Rc) / (xi * Hm0)) ** 1.3));
  const qm = 0.09 * Math.exp(-(((1.5 * Rc) / Hm0) ** 1.3));
  const r = overtoppingSlope({ Hm0, Tm10, tanAlpha: tanA, Rc });
  rel(r.qBreaking, qb * sc, 1e-12);
  rel(r.qMax, qm * sc, 1e-12);
  rel(r.q, Math.min(qb, qm) * sc, 1e-12);
  assert.equal(r.regime, 'breaking');
  // 2.7·3/(1.767·2) = 2.292, 2.292^1.3 = 2.94, 0.0398·1.767·e^−2.94 = 3.72e-3, × √(g·8) = 8.86 → 33 l/s per m
  assert.ok(Math.abs(r.q - 0.033) < 5e-4, `q = ${r.q}`);
  const d = overtoppingSlope({ Hm0, Tm10, tanAlpha: tanA, Rc, design: true });
  rel(d.qBreaking, (0.026 / Math.sqrt(tanA)) * xi * Math.exp(-(((2.5 * Rc) / (xi * Hm0)) ** 1.3)) * sc, 1e-12);
  rel(d.qMax, 0.1035 * Math.exp(-(((1.35 * Rc) / Hm0) ** 1.3)) * sc, 1e-12);
  assert.ok(d.q > r.q);
  // Zero freeboard: the exponentials vanish.
  const z = overtoppingSlope({ Hm0, Tm10: 12, tanAlpha: 0.5, Rc: 0 });
  rel(z.qDimensionless, 0.09, 1e-12);
  // Roughness and obliquity reduce the discharge; influence factors enter as Rc/(γ…).
  const rough = overtoppingSlope({ Hm0, Tm10, tanAlpha: tanA, Rc, gammaF: 0.55, gammaBeta: 0.9 });
  rel(rough.qBreaking, (0.023 / Math.sqrt(tanA)) * xi * Math.exp(-(((2.7 * Rc) / (xi * Hm0 * 0.55 * 0.9)) ** 1.3)) * sc, 1e-12);
  assert.ok(rough.q < r.q);
  assert.throws(() => overtoppingSlope({ Hm0, Tm10, tanAlpha: tanA, Rc: -1 }), RangeError);
});

test('steep-slope coefficients join the slope maximum at cot α = 2 and the vertical-wall formula at cot α = 0', () => {
  const a2 = overtoppingSteepSlope({ Hm0: 1.5, cotAlpha: 2, Rc: 1.2 });
  rel(a2.a, 0.09, 1e-12);
  rel(a2.b, 1.5, 1e-12);
  rel(a2.q, overtoppingSlope({ Hm0: 1.5, Tm10: 30, tanAlpha: 0.5, Rc: 1.2 }).qMax, 1e-12);
  const a0 = overtoppingSteepSlope({ Hm0: 1.5, cotAlpha: 0, Rc: 1.2 });
  assert.ok(Math.abs(a0.a - 0.047) < 5e-4, `a = ${a0.a}`); // 0.09 − 0.01·2^2.1 = 0.0471
  rel(a0.b, 2.35, 1e-12);
  rel(a0.q, overtoppingVertical({ Hm0: 1.5, Rc: 1.2 }).q, 3e-3);
});

test('EurOtop (2018) vertical walls: regime selection, hand values, continuity of the impulsive formula at Rc/Hm0 = 1.35', () => {
  const Hm0 = 1.5;
  const sc = Math.sqrt(g * Hm0 ** 3);
  // No foreshore
  const nf = overtoppingVertical({ Hm0, Rc: 2 });
  assert.equal(nf.regime, 'no-foreshore');
  rel(nf.q, 0.047 * Math.exp(-(((2.35 * 2) / Hm0) ** 1.3)) * sc, 1e-12);
  rel(overtoppingVertical({ Hm0, Rc: 2, design: true }).q, 0.054 * Math.exp(-(((2.12 * 2) / Hm0) ** 1.3)) * sc, 1e-12);
  // Foreshore, deep at the toe: h²/(Hm0 L) = 64/(1.5·56.2) = 0.76 > 0.23 → non-impulsive
  const ni = overtoppingVertical({ Hm0, Rc: 2, Tm10: 6, h: 8 });
  assert.equal(ni.regime, 'non-impulsive');
  rel(ni.q, 0.05 * Math.exp((-2.78 * 2) / Hm0) * sc, 1e-12);
  rel(overtoppingVertical({ Hm0, Rc: 2, Tm10: 6, h: 8, design: true }).q, 0.062 * Math.exp((-2.61 * 2) / Hm0) * sc, 1e-12);
  // Shallow toe: h = 2.5 m, Tm−1,0 = 10 s: 6.25/(1.5·156.1) = 0.027 → impulsive
  const L = (g * 100) / (2 * Math.PI);
  const im = overtoppingVertical({ Hm0, Rc: 1.5, Tm10: 10, h: 2.5 });
  assert.equal(im.regime, 'impulsive');
  rel(im.impulsiveness, 6.25 / (Hm0 * L), 1e-12);
  const f = Math.sqrt(Hm0 / (2.5 * (Hm0 / L)));
  rel(im.qDimensionless, 0.011 * f * Math.exp(-2.2), 1e-12);
  const hi = overtoppingVertical({ Hm0, Rc: 3, Tm10: 10, h: 2.5 });
  rel(hi.qDimensionless, (0.0014 * f) / 8, 1e-12);
  rel(overtoppingVertical({ Hm0, Rc: 3, Tm10: 10, h: 2.5, design: true }).qDimensionless, (0.002 * f) / 8, 1e-12);
  // The two impulsive branches meet at Rc/Hm0 = 1.35 (to the rounding of the published coefficients).
  const lo = overtoppingVertical({ Hm0, Rc: 1.35 * Hm0 * (1 - 1e-9), Tm10: 10, h: 2.5 }).q;
  const up = overtoppingVertical({ Hm0, Rc: 1.35 * Hm0 * (1 + 1e-9), Tm10: 10, h: 2.5 }).q;
  rel(lo, up, 0.02);
  // Impulsive conditions overtop more than non-impulsive ones for the same freeboard.
  assert.ok(im.q > overtoppingVertical({ Hm0, Rc: 1.5, Tm10: 6, h: 8 }).q);
});

test('rubble mound: roughness, obliquity and crest-berm reduction Cr = 3.06 exp(−1.5 Gc/Hm0)', () => {
  const base = overtoppingRubbleMound({ Hm0: 3, Tm10: 9, tanAlpha: 0.5, Rc: 4, gammaF: 0.4 });
  assert.equal(base.Cr, 1);
  rel(base.q, 0.09 * Math.exp(-(((1.5 * 4) / (3 * 0.4)) ** 1.3)) * Math.sqrt(g * 27), 1e-12);
  const wide = overtoppingRubbleMound({ Hm0: 3, Tm10: 9, tanAlpha: 0.5, Rc: 4, gammaF: 0.4, Gc: 6 });
  rel(wide.Cr, 3.06 * Math.exp(-3), 1e-12);
  rel(wide.q, base.q * wide.Cr, 1e-12);
  // Cr only reduces: for Gc < 0.75 Hm0 it is capped at 1.
  assert.equal(overtoppingRubbleMound({ Hm0: 3, Tm10: 9, tanAlpha: 0.5, Rc: 4, gammaF: 0.4, Gc: 2 }).Cr, 1);
});

test('tolerable-discharge classification', () => {
  const rows = classifyOvertopping({ q: 2e-3, Vmax: 0.5 });
  assert.equal(rows.length, TOLERABLE_OVERTOPPING.length);
  const byId = Object.fromEntries(rows.map((r) => [r.id, r]));
  assert.equal(byId['people-crest-Hm0-3'].ok, false);
  assert.equal(byId['people-crest-Hm0-1'].ok, true);
  assert.equal(byId['vehicles-crest-Hm0-3'].ok, true);
  rel(byId['people-crest-Hm0-2'].utilisationQ, 2, 1e-12);
  rel(byId['people-crest-Hm0-1'].utilisationV, 0.5 / 0.6, 1e-12);
  const custom = classifyOvertopping({ q: 2e-3, limits: [{ id: 'x', group: 'g', description: 'd', q: 1e-3, Vmax: null }] });
  assert.equal(custom[0].ok, false);
});

test('overtopping volumes: Weibull mean reproduces q, exceedance of Vmax is 1/Now, Pov from Rayleigh run-up', () => {
  const o = { q: 2e-3, Hm0: 2, Tm10: 6, Tm: 5, duration: 3 * 3600, Rc: 3, Ru2: 4.2 };
  const v = overtoppingVolumes(o);
  rel(v.Pov, Math.exp(-((Math.sqrt(-Math.log(0.02)) * 3 / 4.2) ** 2)), 1e-12);
  // A crest at the 2 % run-up level is overtopped by 2 % of the waves.
  rel(overtoppingVolumes({ ...o, Rc: 4.2 }).Pov, 0.02, 1e-12);
  rel(v.b, 0.73 + 55 * (2e-3 / (g * 2 * 6)) ** 0.8, 1e-12);
  rel(v.Nw, 2160, 1e-12);
  // Volume balance: number of overtopping waves × mean volume = q × duration.
  rel(v.Now * v.Vmean, o.q * o.duration, 1e-12);
  // Numerical mean of the Weibull distribution equals Vmean = a Γ(1 + 1/b).
  let mean = 0;
  const nq = 200000;
  for (let i = 0; i < nq; i++) mean += v.volumeAt(1 - (i + 0.5) / nq) / nq;
  rel(mean, v.Vmean, 2e-3);
  rel(v.exceedance(v.Vmax), 1 / v.Now, 1e-9);
  rel(v.exceedance(v.volumeAt(0.1)), 0.1, 1e-12);
  assert.ok(v.Vmax > 5 * v.Vmean);
  // Exponential distribution (b = 1): a = mean, Vmax = a ln(Now)
  const e = overtoppingVolumes({ ...o, b: 1 });
  rel(e.a, e.Vmean, 1e-10);
  rel(e.Vmax, e.Vmean * Math.log(e.Now), 1e-10);
});

test('sea-level-rise sensitivity: freeboard loss, depth-limited growth, crest raise that restores q', () => {
  const evaluate = (s) => overtoppingSlope({ Hm0: s.Hm0, Tm10: 6, tanAlpha: 1 / 3, Rc: s.Rc }).q;
  const r = seaLevelRiseSensitivity({ rises: [0, 0.25, 0.5, 1.0], Rc: 3, Hm0: 2, evaluate });
  rel(r[0].factor, 1, 1e-12);
  assert.equal(r[0].crestRaise, 0);
  for (let i = 1; i < r.length; i++) {
    assert.ok(r[i].factor > r[i - 1].factor);
    rel(r[i].Rc, 3 - r[i].rise, 1e-12);
    rel(r[i].crestRaise, r[i].rise, 1e-9); // without wave growth the crest must rise by exactly the SLR
  }
  const d = seaLevelRiseSensitivity({ rises: [0.5], Rc: 3, Hm0: 2, h: 4, depthLimited: true, evaluate });
  rel(d[0].Hm0, 2 * 4.5 / 4, 1e-12);
  assert.ok(d[0].factor > r[2].factor);
  assert.ok(d[0].crestRaise > 0.5);
  rel(evaluate({ Rc: 3 - 0.5 + d[0].crestRaise, Hm0: d[0].Hm0 }), evaluate({ Rc: 3, Hm0: 2 }), 1e-6);
});
