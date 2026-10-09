import test from 'node:test';
import assert from 'node:assert/strict';
import * as R from '../../src/engine/lines/riser.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const DEG = Math.PI / 180;
const pp = R.pipeProperties({ D: 0.2731, t: 0.0206, rhoInt: 800 });
const pipe = { EA: pp.EA, EI: pp.EI, w: pp.w, m: pp.m, ma: pp.ma, D: 0.2731, Cdn: 1.0 };

test('pipe properties: hand calculation for 273.1 × 20.6 mm steel pipe with 800 kg/m³ contents', () => {
  const Ae = Math.PI * 0.2731 ** 2 / 4, Ai = Math.PI * 0.2319 ** 2 / 4;
  assert.ok(rel(pp.As, Ae - Ai) < 1e-12);
  assert.ok(rel(pp.I, Math.PI * (0.2731 ** 4 - 0.2319 ** 4) / 64) < 1e-12);
  // apparent weight = (steel + contents) g − buoyancy of the outer area
  const w = (7850 * (Ae - Ai) + 800 * Ai - 1025 * Ae) * 9.80665;
  assert.ok(rel(pp.w, w) < 1e-12);
  // buoyancy modules: 0.9 m OD syntactic foam of 450 kg/m³ over the whole length
  const pb = R.pipeProperties({ D: 0.2731, t: 0.0206, rhoInt: 800, buoyancy: { Dout: 0.9, rho: 450 } });
  const Am = Math.PI * (0.9 ** 2 - 0.2731 ** 2) / 4;
  assert.ok(rel(pb.w, w + (450 - 1025) * Am * 9.80665) < 1e-12);
  assert.ok(pb.w < 0);
});

test('effective tension round trip and Lamé wall stress (closed-end cylinder: axial stress = Lamé constant A at Te = 0)', () => {
  const Ai = 0.04, Ae = 0.06;
  const Te = R.effectiveTension({ Tw: 1e6, pi: 2e7, pe: 1e7, Ai, Ae });
  assert.equal(Te, 1e6 - 2e7 * Ai + 1e7 * Ae);
  assert.ok(rel(R.trueWallTension({ Te, pi: 2e7, pe: 1e7, Ai, Ae }), 1e6) < 1e-12);
  const D = 0.4, t = 0.02, pi = 3e7, pe = 5e6, ro = 0.2, ri = 0.18;
  const s = R.wallStresses({ Te: 0, pi, pe, D, t });
  const A = (pi * ri * ri - pe * ro * ro) / (ro * ro - ri * ri);
  assert.ok(rel(s.axial, A) < 1e-12);
  assert.ok(rel(s.hoopInner, A + (pi - pe) * ro * ro / (ro * ro - ri * ri)) < 1e-12);
  assert.ok(rel(s.hoopInner + s.radialInner, s.hoopOuter + s.radialOuter) < 1e-12);   // σθ + σr constant
  // pure bending: outer-fibre stress M·ro/I
  const b = R.wallStresses({ Te: 0, M: 1e5, D, t });
  assert.ok(rel(b.vonMises, 1e5 * ro / (Math.PI * (ro ** 4 - ri ** 4) / 4)) < 1e-12);
});

test('analytic catenary riser: angle, top tension and span are mutually consistent (T = H + w h)', () => {
  const c = R.catenaryRiser({ h: 1000, w: 1000, hangoffAngleDeg: 12 });
  assert.ok(rel(c.topTension, c.H + 1000 * 1000) < 1e-12);
  assert.ok(rel(Math.tan(12 * DEG), c.H / c.V) < 1e-12);
  const c2 = R.catenaryRiser({ h: 1000, w: 1000, topTension: c.topTension });
  assert.ok(rel(c2.hangoffAngleDeg, 12) < 1e-9);
  // y = a (cosh(x/a) − 1) at the top
  assert.ok(rel(c.a * (Math.cosh(c.horizontalSpan / c.a) - 1), 1000) < 1e-9);
});

test('SCR with bending stiffness: global response matches the catenary, TDZ curvature follows the boundary layer (w/H)(1 − e^(−ξ/λ))', () => {
  const r = R.solveScr({ pipe, h: 1000, hangoffAngleDeg: 12, seabed: { kn: 1e7 }, nElem: 600 });
  const c = r.catenary;
  assert.ok(rel(r.topTension, c.topTension) < 2e-3);
  assert.ok(Math.abs(r.hangoffAngleDeg - 12) < 0.1);
  assert.ok(rel(r.bottomTension, c.H) < 5e-3);
  const bl = R.boundaryLayerTdp({ EI: pipe.EI, H: c.H, w: pipe.w });
  assert.ok(rel(bl.lambda, Math.sqrt(pipe.EI / c.H)) < 1e-12);
  // touchdown sits within ~λ of the catenary touchdown
  assert.ok(Math.abs(r.tdp.s - c.suspendedLength) < 1.5 * bl.lambda);
  const at = (s) => { let k = 0; while (r.s[k + 1] < s) k++; const f = (s - r.s[k]) / (r.s[k + 1] - r.s[k]); return r.curvature[k] + f * (r.curvature[k + 1] - r.curvature[k]); };
  for (const m of [1.5, 2, 3]) {
    const fe = at(r.tdp.s - m * bl.lambda), an = bl.catenaryCurvature * (1 - Math.exp(-m));
    assert.ok(rel(fe, an) < 0.08, `ξ = ${m}λ: ${fe} vs ${an}`);
  }
  // seabed reaction far behind the touchdown equals the apparent weight
  assert.ok(rel(r.reaction[r.reaction.length - 10], pipe.w) < 1e-3);
});

test('hang-off boundary layer: clamped and flex-joint end moments vs Δθ / (1/k + 1/sqrt(EI·T))', () => {
  const rc = R.solveScr({ pipe, h: 1000, hangoffAngleDeg: 12, topCondition: { type: 'clamped', angleDeg: 14 } });
  const an = R.boundaryLayerHangoff({ EI: pipe.EI, T: rc.topTension, deltaAngle: 2 * DEG });
  assert.ok(rel(Math.abs(rc.moment[0]), an.endMoment) < 0.03, `${rc.moment[0]} vs ${an.endMoment}`);
  const rf = R.solveScr({ pipe, h: 1000, hangoffAngleDeg: 12, topCondition: { type: 'flex', angleDeg: 14, kRot: 2e5 } });
  const af = R.boundaryLayerHangoff({ EI: pipe.EI, T: rf.topTension, deltaAngle: 2 * DEG, kRot: 2e5 });
  assert.ok(rel(Math.abs(rf.moment[0]), af.endMoment) < 0.03);
  assert.ok(rel(af.endMoment, 2 * DEG / (1 / 2e5 + 1 / Math.sqrt(pipe.EI * rf.topTension))) < 1e-12);
});

test('far/near vessel offsets: top tension vs the catenary of the same total length at the new span', () => {
  const base = R.solveScr({ pipe, h: 1000, hangoffAngleDeg: 12 });
  const L = base.s[base.s.length - 1], X0 = base.config.bottom.x, w = pipe.w, h = 1000;
  const cases = R.riserOffsetCases(base, [-40, 60]);
  for (const c of cases) {
    const X = X0 + c.offset;
    let lo = 1, hi = 1e5;      // L − S + a asinh(S/a) = X, S = sqrt(h(h + 2a)), solved by bisection
    for (let k = 0; k < 200; k++) { const a = 0.5 * (lo + hi), S = Math.sqrt(h * (h + 2 * a)); if (L - S + a * Math.asinh(S / a) > X) hi = a; else lo = a; }
    const a = 0.5 * (lo + hi);
    assert.ok(rel(c.topTension, w * (a + h)) < 4e-3, `offset ${c.offset}: ${c.topTension} vs ${w * (a + h)}`);
  }
  assert.ok(cases[1].topTension > base.topTension && cases[0].topTension < base.topTension);
  assert.ok(cases[1].tdpExcursion > 0 && cases[0].tdpExcursion < 0);
  const qs = R.quasiStaticHeave(base, { raoHeave: 0.8, waveAmp: 2.5 });
  assert.equal(qs.heaveAmp, 2);
  // raising the hang-off 2 m lifts roughly 2 m·(dS/dh) more pipe: tension amplitude ≈ dT/dh · 2 from the catenary
  const c1 = R.catenaryRiser({ h: 1002, w, H: base.catenary.H }), c0 = R.catenaryRiser({ h: 998, w, H: base.catenary.H });
  assert.ok(qs.tensionAmplitude > 0 && qs.tensionAmplitude < 3 * 0.5 * (c1.topTension - c0.topTension) + 0.2 * w * 1000);
});

test('lazy-wave riser: equilibrium (ΣFz = net weight, H constant), sag and hog bends, buoyancy-length solver hits the target clearance', () => {
  const o = { pipe, buoyant: { ...pipe, w: -1.2 * pipe.w }, upperLength: 900, buoyancyLength: 300, lowerLength: 700, top: { x: 0, z: 800 }, bottom: { x: 1100, z: 0 } };
  const r = R.solveWaveRiser(o);
  assert.ok(r.sag && r.hog && r.hog.z > r.sag.z && r.tdp.s > r.hog.s);
  let Rv = 0; const n = r.s.length - 1;
  for (let k = 0; k <= n; k++) { const tr = 0.5 * ((k > 0 ? r.s[k] - r.s[k - 1] : 0) + (k < n ? r.s[k + 1] - r.s[k] : 0)); Rv += r.reaction[k] * tr; }
  const net = pipe.w * (900 + 700) - 1.2 * pipe.w * 300;
  assert.ok(rel(r.endA.Fz + r.endB.Fz + Rv, net) < 1e-6);
  assert.ok(rel(-r.endA.Fx, r.endB.Fx) < 1e-6);
  assert.ok(rel(r.netBuoyancy, 1.2 * pipe.w * 300) < 1e-12);
  // tension at the sag bend equals the horizontal force (horizontal tangent)
  const kSag = r.s.findIndex((s) => s === r.sag.s);
  assert.ok(rel(r.Te[kSag], r.endB.Fx) < 5e-3);
  const sol = R.solveBuoyancyLength({ ...o, targetSagClearance: 150, minLength: 200, maxLength: 450 });
  assert.ok(Math.abs(sol.result.sagClearance - 150) < 0.5, `${sol.result.sagClearance}`);
  assert.ok(sol.buoyancyLength > 300 && sol.buoyancyLength < 450);
});

test('steep-wave variant: riser base clamped vertical carries a bending moment and no seabed is needed', () => {
  const r = R.solveWaveRiser({ pipe, buoyant: { ...pipe, w: -1.5 * pipe.w }, upperLength: 800, buoyancyLength: 250, lowerLength: 150,
    top: { x: 0, z: 800 }, bottom: { x: 500, z: 0 }, seabed: null, bottomCondition: { type: 'clamped', angleDeg: 0 } });
  assert.ok(Math.abs(r.bottomAngleDeg) < 1e-6);
  assert.ok(r.tdp === null && r.sag && r.hog);
  assert.ok(r.bottomTension > 0);
});

test('top-tensioned riser as a string with linearly varying tension: flex-joint angles = C/T, C = offset·w/ln(Tt/Tb)', () => {
  const L = 1000, w = 1500, TF = 1.3, off = 30;
  const r = R.solveTtr({ length: L, pipe: { EA: 5e9, EI: 1e3, w }, tensionFactor: TF, offset: off });
  const Tt = TF * w * L, Tb = Tt - w * L, C = off * w / Math.log(Tt / Tb);
  assert.ok(rel(r.bottomTension, Tb) < 5e-3);     // axial tension is Tb/cos(angle)
  assert.ok(rel(r.tensionFactor, TF) < 1e-3);
  assert.ok(rel(r.lowerAngleDeg, Math.atan(C / Tb) / DEG) < 0.02, `${r.lowerAngleDeg}`);
  assert.ok(rel(r.upperAngleDeg, Math.atan(C / Tt) / DEG) < 0.02, `${r.upperAngleDeg}`);
  // tensioner stiffness: converged tension satisfies T = T0 + k·stroke
  const k = 2e5;
  const r2 = R.solveTtr({ length: L, pipe: { EA: 5e9, EI: 2e8, w, D: 0.5, Cdn: 1 }, tensionFactor: TF, offset: off, tensionerStiffness: k, kFlexLower: 1e5, current: { U: 0.8 } });
  assert.ok(rel(r2.topTension, Tt + k * r2.stroke) < 1e-3);
  assert.ok(r2.stroke > 0 && r2.lowerAngleDeg > r.lowerAngleDeg);       // current pushes the riser further over
});

test('drilling riser offset sweep: offset at the 2° lower flex-joint limit vs the string closed form; drift-off closed form', () => {
  const L = 1000, w = 1500, TF = 1.3;
  const d = R.drillingRiserOffsets({ length: L, pipe: { EA: 5e9, EI: 1e3, w }, tensionFactor: TF, offsets: [0, 10, 20, 30, 40], nElem: 60 });
  const Tt = TF * w * L, Tb = Tt - w * L;
  const offAt2 = Math.tan(2 * DEG) * Tb * Math.log(Tt / Tb) / w;
  assert.ok(rel(d.offsetAtMeanLimit, offAt2) < 0.03, `${d.offsetAtMeanLimit} vs ${offAt2}`);
  assert.ok(Number.isNaN(d.offsetAtMaxLimit) || d.offsetAtMaxLimit > d.offsetAtMeanLimit);
  assert.ok(rel(d.rows[2].ucMean, Math.abs(d.rows[2].lowerAngleDeg) / 2) < 1e-12);
  const dr = R.driftOff({ mass: 6e7, force: 2e6, dragCoefficient: 5e5, limitOffset: 80, edsTime: 30, reactionTime: 20 });
  assert.equal(dr.terminalVelocity, 2); assert.equal(dr.timeConstant, 60);
  // x(t) = v∞ τ ln cosh(t/τ) = 80 m at t = τ·acosh(e^(80/120))
  assert.ok(rel(dr.timeToLimit, 60 * Math.acosh(Math.exp(80 / 120))) < 1e-12);
  assert.ok(rel(dr.redOffset, 120 * Math.log(Math.cosh((dr.timeToLimit - 30) / 60))) < 1e-9);
  assert.ok(dr.yellowOffset < dr.redOffset && dr.redOffset < 80);
  // early motion is F t²/(2M)
  assert.ok(rel(dr.x[1], 0.5 * (2e6 / 6e7) * dr.t[1] ** 2) < 0.01);
});

test('free-standing riser: base reactions balance can uplift, tower weight and jumper catenary forces', () => {
  const o = { length: 1200, pipe: { EA: 8e9, EI: 3e8, w: 2000 }, canNetUplift: 5e6, jumper: { length: 420, w: 900, vesselX: 250, vesselZ: 1350 } };
  const f = R.solveFreeStandingRiser(o);
  assert.ok(rel(-f.rod.endA.Fx, f.jumperH) < 1e-6);
  assert.ok(rel(-f.rod.endA.Fz, 5e6 - 2000 * 1200 - f.jumperTowerV) < 1e-6);
  assert.ok(rel(f.jumperTowerV + f.jumperVesselV, 900 * 420) < 1e-9);
  // jumper catenary: L² = h² + (2a sinh(X/2a))²
  const a = f.jumperH / 900, X = 250 - f.topOffset, h = 1350 - (1200 - f.topSetDown);
  assert.ok(rel(Math.hypot(h, 2 * a * Math.sinh(X / (2 * a))), 420) < 1e-5);
  // pinned base: moment equilibrium about the base, H·z_top + ∫w x ds = (uplift − V1)·x_top
  const r = f.rod; let mw = 0;
  for (let k = 0; k < r.s.length; k++) mw += 2000 * (1200 / 80) * (k === 0 || k === 80 ? 0.5 : 1) * r.x[k];
  assert.ok(rel(f.jumperH * r.z[80] + mw, (5e6 - f.jumperTowerV) * r.x[80]) < 1e-4);
});

test('stress along a riser and trench profile helper', () => {
  const r = R.solveScr({ pipe, h: 1000, hangoffAngleDeg: 12 });
  const st = R.stressAlong(r, { D: 0.2731, t: 0.0206, pIntRef: 2e7, zIntRef: 1000, rhoInt: 800, zSurface: 1000 });
  const k = 5, pi = 2e7 + 800 * 9.80665 * (1000 - r.z[k]), pe = 1025 * 9.80665 * (1000 - r.z[k]);
  assert.ok(rel(st.vonMises[k], R.wallStresses({ Te: r.Te[k], M: r.moment[k], pi, pe, D: 0.2731, t: 0.0206 }).vonMises) < 1e-12);
  assert.ok(st.max >= st.vonMises[k]);
  const tp = R.trenchProfile({ zSeabed: 0, xStart: 100, length: 50, depth: 1.2 });
  assert.equal(tp(90), 0); assert.ok(Math.abs(tp(125) + 1.2) < 1e-12); assert.equal(tp(151), 0);
  const rt = R.solveScr({ pipe, h: 1000, hangoffAngleDeg: 12, seabed: { kn: 1e6, profile: R.trenchProfile({ xStart: r.tdp.x - 20, length: 80, depth: 1 }) } });
  let zmin = 0; for (const z of rt.z) zmin = Math.min(zmin, z);
  assert.ok(zmin < -0.5 && zmin > -1.1, `${zmin}`);       // the pipe beds into the trench
  assert.throws(() => R.catenaryRiser({ h: -1, w: 1, H: 1 }), /argument h/);
});
