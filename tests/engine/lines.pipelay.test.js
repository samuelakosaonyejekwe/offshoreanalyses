import test from 'node:test';
import assert from 'node:assert/strict';
import * as P from '../../src/engine/lines/pipelay.js';
import { pipeProperties } from '../../src/engine/lines/riser.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const DEG = Math.PI / 180;
const D = 0.4064, t = 0.0191;
const pp = pipeProperties({ D, t, coatings: [{ thickness: 0.05, rho: 3040 }] });
const pipe = { EA: pp.EA, EI: pp.EI, w: pp.w };
const base = { pipe, D, hingeHeight: 160, stingerRadius: 120, stingerLength: 70 };

test('S-lay at high tension: pipe seated on the stinger has curvature 1/R, overbend strain = D/(2R) + T/EA, H = T − w·z, sagbend vs stiffened catenary', () => {
  const T = 1.5e6;
  const r = P.sLay({ ...base, tension: T });
  assert.ok(rel(r.topTension, T) < 1e-3);
  assert.ok(!r.bearsOnTip && r.departureAngleDeg < r.stingerTipAngleDeg);
  // seated nodes well inside the contact arc
  const rod = r.rod;
  for (const s of [15, 30, 45]) { const k = rod.s.findIndex((v) => v >= s); assert.ok(rel(-rod.curvature[k], 1 / 120) < 0.01, `s=${s}: ${rod.curvature[k]}`); }
  assert.ok(rel(r.overbendStrain, r.overbendStrainAnalytic) < 0.03, `${r.overbendStrain} vs ${r.overbendStrainAnalytic}`);
  assert.ok(rel(r.overbendStrainAnalytic, D / 240 + T / pp.EA) < 1e-12);
  // frictionless rollers: tension changes only with elevation
  assert.ok(rel(r.bottomTension, T - pp.w * 160) < 0.01);
  // sagbend curvature vs the boundary-layer stiffened catenary hanging from the lift-off elevation
  const zLift = 160 - 120 * (1 - Math.cos(r.departureAngleDeg * DEG));
  const sc = P.stiffenedCatenary({ EI: pp.EI, H: r.bottomTension, w: pp.w, h: zLift, topFree: false });
  assert.ok(rel(r.sagbendCurvature, sc.sagbendCurvature) < 0.02, `${r.sagbendCurvature} vs ${sc.sagbendCurvature}`);
  // global vertical equilibrium: rollers + seabed + end reactions carry the pipe weight
  const n = rod.s.length - 1; let Rv = 0;
  for (let k = 0; k <= n; k++) {
    const tr = 0.5 * ((k > 0 ? rod.s[k] - rod.s[k - 1] : 0) + (k < n ? rod.s[k + 1] - rod.s[k] : 0));
    const onStinger = rod.z[k] > 50, ang = onStinger ? Math.atan2(rod.z[k] - 40, rod.x[k]) : Math.PI / 2;
    Rv += rod.reaction[k] * tr * Math.sin(ang);
  }
  assert.ok(rel(Rv + rod.endA.Fz + rod.endB.Fz, pp.w * rod.s[n]) < 1e-4);
  // catenary lay-back from lift-off: a·asinh(S/a)
  const a = r.bottomTension / pp.w, S = Math.sqrt(zLift * (zLift + 2 * a));
  assert.ok(rel(r.tdpX - 120 * Math.sin(r.departureAngleDeg * DEG), a * Math.asinh(S / a)) < 0.03);
});

test('S-lay at low tension: the pipe bears on the last roller and the overbend strain exceeds D/(2R); minimum tension search meets the sagbend limit', () => {
  const lo = P.sLay({ ...base, tension: 6e5 });
  assert.ok(lo.bearsOnTip && lo.departureAngleDeg > lo.stingerTipAngleDeg);
  assert.ok(lo.overbendStrain > 1.5 * lo.overbendStrainAnalytic);
  const hi = P.sLay({ ...base, tension: 1.5e6 });
  assert.ok(lo.sagbendStrain > hi.sagbendStrain && lo.layBack < hi.layBack);
  const mt = P.minimumLayTension({ ...base, sagbendStrainLimit: 0.0009, tensionMin: 6e5, tensionMax: 1.5e6 });
  assert.ok(rel(mt.result.sagbendStrain, 0.0009) < 0.01, `${mt.result.sagbendStrain}`);
  assert.ok(mt.tension > 6e5 && mt.tension < 9e5);
});

test('J-lay: top tension from the catenary matched to the tower angle; sagbend vs stiffened catenary', () => {
  const j = P.jLay({ pipe, D, depth: 800, towerAngleDeg: 12 });
  const sa = Math.sin(12 * DEG), a = 800 * sa / (1 - sa);
  assert.ok(rel(j.topTension, pp.w * (a + 800)) < 3e-3);
  assert.ok(rel(j.bottomTension, pp.w * a) < 0.01);
  const sc = P.stiffenedCatenary({ EI: pp.EI, H: j.bottomTension, w: pp.w, h: 800, D, EA: pp.EA, topFree: false });
  assert.ok(rel(j.sagbendCurvature, sc.sagbendCurvature) < 0.05, `${j.sagbendCurvature} vs ${sc.sagbendCurvature}`);
  assert.ok(rel(j.layBack, a * Math.asinh(Math.sqrt(800 * (800 + 2 * a)) / a)) < 0.03);
  // stiffened catenary tends to the natural catenary curvature 1/a when EI → 0
  const soft = P.stiffenedCatenary({ EI: 1, H: 1e5, w: 1000, h: 500, topFree: false });
  assert.ok(rel(soft.sagbendCurvature, 1 / 100) < 1e-3);
});

test('reel-lay: strain geometry, thin-tube plastic moment limits My and (4/π)My, residual curvature after spring-back', () => {
  const r = P.reelLay({ D: 0.3239, t: 0.0191, reelRadius: 8.0, alignerRadius: 9.5, fy: 450e6 });
  assert.ok(rel(r.reelStrain, 0.16195 / (8 + 0.16195)) < 1e-12);
  assert.ok(rel(r.alignerStrain, 0.16195 / (9.5 + 0.16195)) < 1e-12);
  assert.ok(Number.isNaN(r.straightenerStrain));
  const rm = (0.3239 - 0.0191) / 2, My = Math.PI * rm * rm * 0.0191 * 450e6;
  // at first yield the moment is My and nothing remains after unloading (thin-wall I = π r³ t to within t/D terms)
  const Ry = 1 / (450e6 / (207e9 * rm)) - 0.16195;
  const y = P.reelLay({ D: 0.3239, t: 0.0191, reelRadius: Ry, fy: 450e6 });
  assert.ok(rel(y.reelMoment, My) < 1e-9);
  assert.ok(Math.abs(y.residualCurvature) < 0.01 * y.reelCurvature);
  // far beyond yield the moment tends to the fully plastic value (4/π) My
  assert.ok(r.reelMoment > My && r.reelMoment < 4 / Math.PI * My);
  assert.ok(rel(P.reelLay({ D: 0.3239, t: 0.0191, reelRadius: 0.5, fy: 450e6 }).reelMoment, 4 / Math.PI * My) < 0.01);
  assert.ok(r.residualCurvature > 0 && r.residualCurvature < r.reelCurvature);
  // ovalisation expressions: hand values
  assert.ok(rel(r.ovalityF101, 0.005 + 0.030 * (1 + 0.3239 / (120 * 0.0191)) * (2 * r.reelStrain * 0.3239 / 0.0191) ** 2) < 1e-12);
  assert.ok(rel(r.ovalityBrazier - 0.005, 2 * 0.91 * rm ** 4 * r.reelCurvature ** 2 / 0.0191 ** 2) < 1e-12);
});

test('cable lay over a chute: catenary identities, sidewall pressure T/R, slack; Zajac straight-line slack lay force balance', () => {
  const c = P.cableLay({ w: 300, h: 100, departureAngleDeg: 60, chuteRadius: 3, mbr: 2.4, swpLimit: 20e3, payoutSpeed: 0.515, vesselSpeed: 0.5 });
  assert.ok(rel(c.H, 300 * 100 * 0.5 / 0.5) < 1e-9);                 // H = w h cos α/(1 − cos α)
  assert.ok(rel(c.topTension, c.H + 300 * 100) < 1e-12);
  assert.ok(rel(c.departureAngleDeg, 60) < 1e-9);
  assert.ok(rel(c.sidewallPressure, c.topTension / 3) < 1e-12);
  assert.ok(rel(c.ucChuteRadius, 0.8) < 1e-12 && rel(c.ucTdpRadius, 2.4 / 100) < 1e-9 && rel(c.ucSidewall, c.topTension / 3 / 20e3) < 1e-12);
  assert.ok(rel(c.slack, 0.03) < 1e-9);
  assert.ok(rel(P.cableLay({ w: 300, h: 100, topTension: c.topTension }).layBack, c.layBack) < 1e-9);
  const z = P.cableSlackLay({ w: 300, d: 0.15, V: 1.5, h: 100, Cd: 1.2 });
  const al = z.angleDeg * DEG;
  assert.ok(rel(300 * Math.cos(al), 0.5 * 1025 * 1.2 * 0.15 * (1.5 * Math.sin(al)) ** 2) < 1e-9);
  assert.ok(rel(z.topTension, 3e4) < 1e-12 && rel(z.layBack, 100 / Math.tan(al)) < 1e-12);
});

test('abandonment and recovery, capstan equation e^(μθ) and shore pull-in: hand calculations', () => {
  const ar = P.abandonmentRecovery({ H: 2e5, wPipe: 1500, wWire: 250, depth: 400, headHeights: [0, 100, 400] });
  assert.ok(rel(ar[0].headTension, 2e5) < 1e-12 && rel(ar[0].vesselTension, 2e5 + 250 * 400) < 1e-12);
  assert.ok(rel(ar[1].headTension, 2e5 + 1500 * 100) < 1e-12);      // catenary: T = H + w z
  assert.ok(rel(ar[1].vesselTension, 2e5 + 1500 * 100 + 250 * 300) < 1e-12);
  assert.ok(rel(ar[2].vesselTension, 2e5 + 1500 * 400) < 1e-12);
  // capstan: μ = 0.3 over half a turn → e^(0.3π) = 2.5663
  assert.ok(Math.abs(P.capstan({ T1: 1000, mu: 0.3, theta: Math.PI }) / 1000 - 2.5663) < 1e-4);
  assert.equal(P.capstan({ T1: 1000, mu: 0, theta: 5 }), 1000);
  // pull-in: 500 m on seabed (w = 400 N/m, μ = 0.6) → 120 kN; 30° bend at μ = 0.2; 100 m on rollers up a 5° slope (w = 900, μ = 0.1)
  const sp = P.shorePullIn({ segments: [{ type: 'straight', length: 500, w: 400, mu: 0.6 }, { type: 'bend', mu: 0.2, angleDeg: 30 }, { type: 'straight', length: 100, w: 900, mu: 0.1, slopeDeg: 5 }], startFactor: 1.2 });
  const T1 = 120e3, T2 = T1 * Math.exp(0.2 * 30 * DEG), T3 = T2 + 900 * 100 * (0.1 * Math.cos(5 * DEG) + Math.sin(5 * DEG));
  assert.ok(rel(sp.tensions[0], T1) < 1e-12 && rel(sp.tensions[1], T2) < 1e-12 && rel(sp.winchTension, 1.2 * T3) < 1e-12);
  const pl = P.ploughTowForce({ Cw: 0.4, Wsub: 3e5, Cs: 10, gammaSub: 9000, depth: 1.5, Cd: 2e4, speed: 0.1 });
  assert.ok(rel(pl.force, 1.2e5 + 10 * 9000 * 3.375 + 2e4 * 0.1 * 2.25) < 1e-12);
});

test('quasi-dynamic tension bounds, limiting DAF and DNV-ST-F101 simplified laying criteria: hand calculations', () => {
  const q = P.layDynamicTension({ staticTension: 5e5, w: 1000, EA: 4e9, m: 300, suspendedLength: 400, topAngleDeg: 60, heaveAmp: 1.0, period: 8 });
  const at = Math.sin(60 * DEG), om = 2 * Math.PI / 8;
  assert.ok(rel(q.lower, 1000 + 300 * 400 * om * om * at) < 1e-12);
  assert.ok(rel(q.upper, 4e9 / 400 * at) < 1e-12);
  assert.ok(rel(q.dafUpper, 1 + q.upper / 5e5) < 1e-12 && q.minTensionUpper < 0);
  const ls = P.limitingSeaState({ staticTension: 5e5, maxTension: 9e5, minTension: 5e4, hsGrid: [1, 1.5, 2, 2.5, 3], dafOfHs: (hs) => 1 + 0.3 * hs });
  assert.ok(rel(ls.limitingDaf, 1.8) < 1e-12);       // min(9/5, 2 − 0.1)
  assert.equal(ls.allowableHs, 2.5);
  const cr = P.layCriteriaF101({ grade: 'X65', criterion: 1, overbendStrain: 0.002, D, t, fy: 448e6, sagbendTe: 5e5, sagbendMoment: 2e5, pe: 2e6 });
  assert.equal(cr.overbendLimit, 0.0025);
  assert.ok(rel(cr.ucOverbend, 0.8) < 1e-12);
  assert.equal(P.layCriteriaF101({ grade: 'X52', criterion: 2, overbendStrain: 0.002, D, t, fy: 358e6 }).overbendLimit, 0.0026);
  assert.equal(P.layCriteriaF101({ grade: 'X70', overbendStrain: 0.002, D, t, fy: 485e6 }).overbendLimit, 0.0027);
  const Di = D - 2 * t, A = Math.PI * (D * D - Di * Di) / 4, I = Math.PI * (D ** 4 - Di ** 4) / 64;
  const N = 5e5 - 2e6 * Math.PI * D * D / 4, sh = -2e6 * (D - t) / (2 * t);
  const eq = (sl) => Math.sqrt(sl * sl + sh * sh - sl * sh);
  const seq = Math.max(eq(N / A + 2e5 * D / (2 * I)), eq(N / A - 2e5 * D / (2 * I)));
  assert.ok(rel(cr.sagbendStress, seq) < 1e-12 && rel(cr.ucSagbend, seq / (0.87 * 448e6)) < 1e-12);
  assert.throws(() => P.layCriteriaF101({ grade: 'X80', overbendStrain: 0.002, D, t, fy: 1 }), /grade/);
});
