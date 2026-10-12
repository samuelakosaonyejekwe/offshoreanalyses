import test from 'node:test';
import assert from 'node:assert/strict';
import * as P from '../../src/engine/lines/pipelay.js';
import { pipeProperties } from '../../src/engine/lines/riser.js';
import { solveRodStatic } from '../../src/engine/lines/rod.js';

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

/* ------------------------------------------------------------------------------------------------------- */
/* Robustness of the S-lay solve: the result must not depend on the last bit of the arithmetic.             */

// 323.9 × 14.3 mm X65 pipe, 3 mm coating, 40 mm concrete, empty, in 120 m of water with the hinge 5 m above the surface
const layPipe = { EA: 2879103036.9407005, EI: 34569624.09166873, w: 1198.8017938116636 };
const layBase = { pipe: layPipe, D: 0.3239, hingeHeight: 125, stingerRadius: 150, stingerLength: 90, tension: 1e6, seabed: { kn: 1e6 }, nElem: 100 };
const LAY_KEYS = ['topTension', 'bottomTension', 'overbendStrain', 'sagbendStrain', 'sagbendMoment', 'departureAngleDeg', 'tdpX', 'layBack', 'maxRollerReaction'];
const maxRelDiff = (a, b) => Math.max(...LAY_KEYS.map((k) => rel(a[k], b[k])));
// deterministic pseudo-random numbers in (−0.5, 0.5)
const lcg = (seed) => () => { seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0; return seed / 4294967296 - 0.5; };
// a clean Newton solve from a good starting shape takes a handful of iterations; dozens mean it wandered
const CLEAN = 12;

test('S-lay, 120 m water depth, 1000 kN: converges in a few Newton iterations to T_top = T and H = T − w·z', () => {
  const r = P.sLay(layBase);
  assert.ok(r.rod.converged && r.rod.steps === 1 && r.rod.iterations <= CLEAN, `iterations ${r.rod.iterations}`);
  assert.ok(r.rod.residual < 1e-8 * 2e6);
  assert.ok(rel(r.topTension, 1000025.1876) < 1e-6, `${r.topTension}`);
  assert.ok(rel(r.bottomTension, 1e6 - layPipe.w * 125) < 2e-3, `${r.bottomTension}`);       // frictionless rollers
  assert.ok(!r.bearsOnTip && r.departureAngleDeg < r.stingerTipAngleDeg);
  // the starting shape already carries the lay tension and is seated on its supports
  const m = P.sLayRodModel(layBase), n = m.initial.x.length - 1;
  const strain0 = Math.hypot(m.initial.x[1] - m.initial.x[0], m.initial.z[1] - m.initial.z[0]) / (m.sections[0].length / m.sections[0].nElem) - 1;
  assert.ok(rel(strain0, 1e6 / layPipe.EA) < 0.05, `${strain0}`);
  assert.ok(rel(-m.initial.z[n], layPipe.w / 1e6) < 1e-9);                                    // seabed penetration w/kn
  const k = 10, d = Math.hypot(m.initial.x[k], m.initial.z[k] - (125 - 150));
  assert.ok(150 - d > 0.5 * (1e6 / 150) / 1e8 && 150 - d < 2 * (1e6 / 150) / 1e8, `${150 - d}`); // roller penetration ≈ (T/R)/k
});

test('S-lay is insensitive to perturbations of every input by ±1e-9 (relative)', () => {
  const ref = P.sLay(layBase);
  const vary = (f) => ({
    pipe: { EA: layPipe.EA * f(), EI: layPipe.EI * f(), w: layPipe.w * f() }, D: layBase.D * f(), hingeHeight: 125 * f(), stingerRadius: 150 * f(),
    stingerLength: 90 * f(), tension: 1e6 * f(), seabed: { kn: 1e6 * f() }, nElem: 100,
  });
  const cases = [vary(() => 1 + 1e-9), vary(() => 1 - 1e-9)];
  // one input at a time, each way
  for (let i = 0; i < 9; i++) for (const sg of [1, -1]) { let j = 0; cases.push(vary(() => (j++ === i ? 1 + sg * 1e-9 : 1))); }
  const rnd = lcg(2024);
  for (let t = 0; t < 10; t++) cases.push(vary(() => 1 + Math.sign(rnd()) * 1e-9));
  for (const c of cases) {
    const r = P.sLay(c);
    assert.ok(r.rod.converged && r.rod.iterations <= CLEAN, `iterations ${r.rod.iterations}`);
    assert.ok(maxRelDiff(r, ref) < 1e-6, `${maxRelDiff(r, ref)}`);
  }
});

test('S-lay is insensitive to 1e-10 m jitter of the starting shape and to last-bit noise in the Math functions', () => {
  const NOISY = ['sin', 'cos', 'tan', 'atan', 'atan2', 'asinh', 'hypot', 'exp'];
  // lift-off before the tip, and the pipe bearing on the stinger tip (low tension)
  for (const c of [layBase, { ...layBase, tension: 6e5 }, { ...layBase, tension: 7e5, nElem: 200 }]) {
    const ref = P.sLay(c);
    assert.ok(ref.rod.converged);
    const rnd = lcg(7);
    for (let t = 0; t < 8; t++) {
      const m = P.sLayRodModel(c);
      for (let i = 0; i < m.initial.x.length; i++) { m.initial.x[i] += 2e-10 * rnd(); m.initial.z[i] += 2e-10 * rnd(); }
      const rod = solveRodStatic(m), n = rod.s.length - 1;
      assert.ok(rod.converged && rod.iterations === ref.rod.iterations, `iterations ${rod.iterations} vs ${ref.rod.iterations}`);
      assert.ok(rel(rod.Te[0], ref.rod.Te[0]) < 1e-8 && rel(rod.Te[n], ref.rod.Te[n]) < 1e-8);
      let dz = 0;
      for (let i = 0; i <= n; i++) dz = Math.max(dz, Math.abs(rod.z[i] - ref.rod.z[i]), Math.abs(rod.x[i] - ref.rod.x[i]));
      assert.ok(dz < 1e-7, `shape differs by ${dz} m`);
    }
    // every transcendental result off by up to one unit in the last place, as between JavaScript engines
    const orig = NOISY.map((f) => Math[f]);
    try {
      for (let t = 0; t < 4; t++) {
        NOISY.forEach((f, i) => { Math[f] = (...a) => orig[i](...a) * (1 + 2.220446049250313e-16 * Math.round(2 * rnd())); });
        const r = P.sLay(c);
        assert.ok(r.rod.converged && r.rod.iterations === ref.rod.iterations, `iterations ${r.rod.iterations} vs ${ref.rod.iterations}`);
        assert.ok(maxRelDiff(r, ref) < 1e-6, `${maxRelDiff(r, ref)}`);
      }
    } finally { NOISY.forEach((f, i) => { Math[f] = orig[i]; }); }
  }
});

test('S-lay sweep, tension 600–2500 kN × 60–200 elements: every case converges cleanly; tension below w·z is rejected', () => {
  // All combinations are physically feasible: the tension exceeds the weight of the pipe hanging from the hinge,
  // w·z = 150 kN. Below about 710 kN the natural catenary would leave steeper than the stinger tip (34.4°), so the
  // pipe stays on the rollers to the last metres of the stinger or bears on its tip, which is still an equilibrium.
  // The only infeasible region is T ≤ w·z (no horizontal tension left), rejected as an invalid argument.
  const wz = layPipe.w * 125;
  for (const T of [600e3, 650e3, 700e3, 800e3, 1000e3, 1300e3, 1700e3, 2100e3, 2500e3]) {
    let prev = null;
    for (const nElem of [60, 80, 100, 130, 160, 200]) {
      const r = P.sLay({ ...layBase, tension: T, nElem });
      const tag = `T = ${T / 1e3} kN, nElem = ${nElem}`;
      assert.ok(r.rod.converged && r.rod.iterations <= 2 * CLEAN, `${tag}: iterations ${r.rod.iterations}`);
      assert.ok(rel(r.topTension, T) < 1e-3, `${tag}: top ${r.topTension}`);
      assert.ok(rel(r.bottomTension, T - wz) < 5e-3, `${tag}: bottom ${r.bottomTension}`);
      // bearsOnTip flags a last roller contact within 1.5 elements of the tip, so near the transition (650–720 kN)
      // it depends on the mesh; well below and well above it must not
      if (T <= 600e3) assert.ok(r.bearsOnTip, tag);
      if (T >= 800e3) assert.ok(!r.bearsOnTip && r.departureAngleDeg < r.stingerTipAngleDeg - 1, tag);
      assert.ok(r.sagbendStrain > 0 && r.layBack > 0 && r.suspendedLength > 125, tag);
      // mesh refinement changes the answer only a little (the touchdown point is reported at a node, so it is
      // resolved to one element length, up to about 2 % of the lay-back on the coarsest mesh)
      if (prev) assert.ok(rel(r.sagbendStrain, prev.sagbendStrain) < 0.02 && rel(r.tdpX, prev.tdpX) < 0.02, `${tag}: mesh sensitivity`);
      prev = r;
    }
  }
  // sagbend curvature falls and lay-back grows monotonically with tension (the total sagbend strain does not: its
  // axial part T/EA grows with tension and dominates above about 1000 kN here)
  const rs = [600e3, 800e3, 1000e3, 1500e3, 2500e3].map((T) => P.sLay({ ...layBase, tension: T }));
  for (let i = 1; i < rs.length; i++) assert.ok(rs[i].sagbendCurvature < rs[i - 1].sagbendCurvature && rs[i].layBack > rs[i - 1].layBack);
  assert.throws(() => P.sLay({ ...layBase, tension: 0.99 * wz }), /pipelay: argument tension/);
  assert.throws(() => P.sLay({ ...layBase, tension: wz }), /pipelay: argument tension/);
});

test('S-lay with the pipe bearing hard on the stinger tip (25° kink): converges, and the tension search is continuous', () => {
  for (const T of [6e5, 6.5e5, 7e5, 7.5e5, 8e5, 9e5]) {
    const r = P.sLay({ ...base, tension: T });
    assert.ok(r.rod.converged && r.bearsOnTip, `T = ${T}`);
    assert.ok(rel(r.topTension, T) < 2e-3 && rel(r.bottomTension, T - pp.w * 160) < 0.08, `T = ${T}: ${r.topTension}, ${r.bottomTension}`);
  }
  // neighbouring tensions give neighbouring answers (no jump between equilibria as a node passes the tip)
  const a = P.sLay({ ...base, tension: 7e5 }), b = P.sLay({ ...base, tension: 7e5 * (1 + 1e-7) });
  assert.ok(maxRelDiff(a, b) < 1e-4, `${maxRelDiff(a, b)}`);
});
