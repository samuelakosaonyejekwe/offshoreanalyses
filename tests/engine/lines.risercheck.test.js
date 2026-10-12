import test from 'node:test';
import assert from 'node:assert/strict';
import * as C from '../../src/engine/lines/risercheck.js';
import { solveRodStatic } from '../../src/engine/lines/rod.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const D = 0.3239, t = 0.0191, fy = 448e6, fu = 531e6;
const Di = D - 2 * t, A = Math.PI * (D * D - Di * Di) / 4, I = Math.PI * (D ** 4 - Di ** 4) / 64;

test('API RP 2RD von Mises: pure tension gives σe = T/A against Cf·(2/3)·σy; hand calculation with pressure and bending', () => {
  const a = C.vonMisesApi2rd({ T: 2e6, D, t, sigmaY: fy });
  assert.ok(rel(a.vonMises, 2e6 / A) < 1e-12);
  assert.ok(rel(a.allowable, fy * 2 / 3) < 1e-12);
  assert.ok(rel(C.vonMisesApi2rd({ T: 2e6, D, t, sigmaY: fy, Cf: 1.5 }).uc, a.uc / 1.5) < 1e-12);
  // hand calculation: pi = 30 MPa, pe = 10 MPa, T = 1 MN, M = 200 kN·m
  const pi = 30e6, pe = 10e6, T = 1e6, M = 2e5;
  const sr = -(pe * D + pi * Di) / (D + Di), sh = (pi - pe) * D / (2 * t) - pi;
  const vm = (sz) => Math.sqrt(0.5 * ((sr - sh) ** 2 + (sh - sz) ** 2 + (sz - sr) ** 2));
  const sb = M * (D - t) / (2 * I);
  const b = C.vonMisesApi2rd({ T, M, pi, pe, D, t, sigmaY: fy, Cf: 1.2 });
  assert.ok(rel(b.vonMises, Math.max(vm(T / A + sb), vm(T / A - sb))) < 1e-12);
  assert.ok(rel(b.uc, b.vonMises / (1.2 * fy * 2 / 3)) < 1e-12);
});

test('burst: DNV pb hand calculation; API RP 1111 thick- and thin-wall forms agree for thin pipe', () => {
  const pb = C.burstPressureDnv({ D, t, fy, fu });
  assert.ok(rel(pb, (2 / Math.sqrt(3)) * (2 * t / (D - t)) * Math.min(fy, fu / 1.15)) < 1e-12);
  const chk = C.burstCheckDnv({ pli: 35e6, pe: 5e6, D, t, fy, fu });
  assert.ok(rel(chk.uc, 30e6 * 1.15 * 1.14 / pb) < 1e-12);
  const thick = C.burstPressureApi1111({ D: 0.5, t: 0.01, S: fy, U: fu }), thin = C.burstPressureApi1111({ D: 0.5, t: 0.01, S: fy, U: fu, thin: true });
  assert.ok(rel(thick, thin) < 1e-3);                  // ln(D/Di) = 2t/(D − t) + O((t/D)³)
  assert.ok(rel(C.hoopStress({ pi: 20e6, pe: 0, D, t }).barlow, 20e6 * D / (2 * t)) < 1e-12);
});

test('collapse: the DNV closed form satisfies (pc − pel)(pc² − pp²) = pc·pel·pp·f0·D/t and lies below both limits', () => {
  for (const tt of [0.008, 0.0127, 0.0191, 0.03]) {
    for (const f0 of [0.005, 0.015, 0.03]) {
      const { pc, pel, pp } = C.collapsePressureDnv({ D, t: tt, fy, f0 });
      const lhs = (pc - pel) * (pc * pc - pp * pp), rhs = pc * pel * pp * f0 * D / tt;
      assert.ok(Math.abs(lhs - rhs) < 1e-9 * Math.abs(pel * pp * pp), `t=${tt} f0=${f0}`);
      assert.ok(pc > 0 && pc < Math.min(pel, pp));
    }
  }
  // elastic collapse hand value: 2E(t/D)³/(1 − ν²)
  const thin = C.collapsePressureDnv({ D: 1, t: 0.01, fy });
  assert.ok(rel(thin.pel, 2 * 207e9 * 1e-6 / 0.91) < 1e-12);
  // more ovality always reduces the collapse pressure
  assert.ok(C.collapsePressureDnv({ D, t, fy, f0: 0.03 }).pc < C.collapsePressureDnv({ D, t, fy, f0: 0.005 }).pc);
  const api = C.collapsePressureApi({ D, t, S: fy });
  assert.ok(rel(api.Pc, api.Py * api.Pe / Math.hypot(api.Py, api.Pe)) < 1e-12 && api.Pc < Math.min(api.Py, api.Pe));
  const cc = C.collapseCheckDnv({ D, t, fy, pe: 20e6 });
  assert.ok(rel(cc.uc, 20e6 * 1.14 * 1.15 / cc.pc) < 1e-12);
});

test('propagating buckling: ppr = 35 fy αfab (t/D)^2.5 and API 24 S (t/D)^2.4 hand values', () => {
  const p = C.propagatingBuckling({ D, t, fy, pe: 5e6 });
  assert.ok(rel(p.pprDnv, 35 * fy * (t / D) ** 2.5) < 1e-12);
  assert.ok(rel(p.ppApi, 24 * fy * (t / D) ** 2.4) < 1e-12);
  assert.ok(rel(p.uc, 5e6 * 1.15 * 1.14 / p.pprDnv) < 1e-12);
});

test('DNV LRFD combined loading: hand calculation of αc, Mk, Tk; utilisation is exactly 1 on the capacity surface', () => {
  const g = 1.14 * 1.15, dt = D / t;
  const beta = 0.4 * (60 - dt) / 45, ac = 1 - beta + beta * fu / fy;
  const Mk = fy * ac * (D - t) ** 2 * t, Tk = fy * ac * Math.PI * (D - t) * t;
  const r0 = C.combinedLoadingDnv({ D, t, fy, fu, MF: 1e5 });
  assert.ok(rel(r0.alphaC, ac) < 1e-12 && rel(r0.Mk, Mk) < 1e-12 && rel(r0.Tk, Tk) < 1e-12);
  // pure functional bending at Mk/(γF γSC γm)
  assert.ok(rel(C.combinedLoadingDnv({ D, t, fy, fu, MF: Mk / (1.1 * g) }).uc, 1) < 1e-12);
  // pure tension at Tk/sqrt(γSC γm) with environmental factor
  assert.ok(rel(C.combinedLoadingDnv({ D, t, fy, fu, TeE: Tk / Math.sqrt(g) / 1.3 }).uc, 1) < 1e-12);
  // load-effect factors: Md = 1.1 MF + 1.3 ME + 1.0 MA
  const lf = C.combinedLoadingDnv({ D, t, fy, fu, MF: 100, ME: 10, MA: 1, TeF: 1000, TeE: 100, TeA: 10 });
  assert.ok(rel(lf.Md, 110 + 13 + 1) < 1e-12 && rel(lf.Ted, 1100 + 130 + 10) < 1e-12);
  // internal overpressure: hand evaluation of Eq. 5.28
  const pld = 30e6, pe = 5e6, pb = C.burstPressureDnv({ D, t, fy, fu }), pr = (pld - pe) / pb;
  const qh = pr * 2 / Math.sqrt(3), b2 = (0.4 + qh) * (60 - dt) / 45, ac2 = 1 - b2 + b2 * fu / fy;
  const Mk2 = fy * ac2 * (D - t) ** 2 * t, Tk2 = fy * ac2 * Math.PI * (D - t) * t;
  const ri = C.combinedLoadingDnv({ D, t, fy, fu, MF: 2e5, ME: 1e5, TeF: 1e6, pld, pe });
  const Md = 1.1 * 2e5 + 1.3 * 1e5, Ted = 1.1e6;
  assert.equal(ri.mode, 'internal');
  assert.ok(rel(ri.uc, g * ((Md / Mk2) * Math.sqrt(1 - pr * pr) + (Ted / Tk2) ** 2) + pr * pr) < 1e-12);
  // external overpressure: hand evaluation of Eq. 5.30
  const re = C.combinedLoadingDnv({ D, t, fy, fu, MF: 2e5, TeF: 1e6, pe: 15e6, pmin: 1e6 });
  assert.equal(re.mode, 'external');
  assert.ok(rel(re.uc, g * g * ((2.2e5 / Mk + (1.1e6 / Tk) ** 2) ** 2 + (14e6 / re.pc) ** 2)) < 1e-12);
  assert.throws(() => C.combinedLoadingDnv({ D: 1, t: 0.01, fy, fu }), /D\/t/);
});

test('strain-based local buckling: εc = 0.78 (t/D − 0.01) αh^−1.5 hand value and pressure enhancement', () => {
  const r = C.localBucklingStrain({ strain: 0.004, D, t, alphaH: 0.92 });
  const ec = 0.78 * (t / D - 0.01) * 0.92 ** -1.5;
  assert.ok(rel(r.epsC, ec) < 1e-12 && rel(r.uc, 0.004 * 2.5 / ec) < 1e-12);
  const rp = C.localBucklingStrain({ strain: 0.004, D, t, alphaH: 0.92, pmin: 10e6, pe: 0, pb: 60e6 });
  assert.ok(rel(rp.epsC, ec * (1 + 5.75 * 10 / 60)) < 1e-12);
});

test('tapered stress joint: constant-curvature profile verified with the rod model (curvature uniform at κ0 along the joint)', () => {
  const T = 1.5e6, E = 207e9, dth = 6 * Math.PI / 180;
  const sj = C.stressJointTaper({ deltaAngle: dth, T, D, t, sigmaBendAllow: 250e6, nPoints: 21 });
  assert.ok(sj.length > 0 && sj.rootOD > D);
  assert.ok(rel(sj.curvature, 250e6 / (E * sj.rootOD / 2)) < 1e-9);
  assert.ok(rel(sj.length, dth / sj.curvature - sj.lambdaPipe) < 1e-9);
  const q = sj.length / sj.lambdaPipe;
  assert.ok(rel(sj.rootEI, E * I * (1 + q + 0.5 * q * q)) < 1e-12);
  // rod model: root clamped at Δθ to the pull direction, joint in 20 uniform pieces from root to tip, then pipe
  const sections = [];
  for (let k = 20; k >= 1; k--) sections.push({ length: sj.length / 20, EA: 1e11, EI: 0.5 * (sj.EI[k] + sj.EI[k - 1]), w: 0, nElem: 4 });
  sections.push({ length: 40 * sj.lambdaPipe, EA: 1e11, EI: E * I, w: 0, nElem: 200 });
  const r = solveRodStatic({ sections, endA: { x: 0, z: 0, angle: dth }, endB: { z: 0, Fx: T }, nSteps: 10 });
  for (const k of [4, 20, 40, 60, 76]) assert.ok(rel(Math.abs(r.curvature[k]), sj.curvature) < 0.05, `node ${k}: ${r.curvature[k]} vs ${sj.curvature}`);
  assert.ok(rel(Math.abs(r.moment[0]), sj.rootMoment) < 0.05);
});

test('flex joint, connector envelope (diamond envelope gives T/Tr + M/Mr) and fatigue hot-spot hand calculation', () => {
  const fj = C.flexJointCheck({ angleDeg: 6, allowableDeg: 10, kRot: 5e4 });
  assert.ok(rel(fj.uc, 0.6) < 1e-12 && rel(fj.moment, 5e4 * 6 * Math.PI / 180) < 1e-12);
  const env = [{ T: 5e6, M: 0 }, { T: 0, M: 2e6 }, { T: -5e6, M: 0 }, { T: 0, M: -2e6 }];
  const cu = C.connectorUtilisation({ T: 2e6, M: 0.5e6, envelope: env });
  assert.ok(rel(cu.uc, 2 / 5 + 0.5 / 2) < 1e-12);
  assert.ok(rel(cu.capacityT / cu.capacityM, 4) < 1e-12);
  // Δσ = SCF (ΔT/A + ΔM r/I); S–N: N = 1e12·S^-3 (S in MPa), 1e6 cycles
  const f = C.fatigueHotSpot({ dT: 2e5, dM: 5e4, D, t, scf: 1.3, cycles: 1e6, snCyclesToFailure: (S) => 1e12 * (S / 1e6) ** -3 });
  const S = 1.3 * (2e5 / A + 5e4 * (D / 2) / I);
  assert.ok(rel(f.stressRange, S) < 1e-12);
  assert.ok(rel(f.damage, 1e6 / (1e12 * (S / 1e6) ** -3)) < 1e-12);
});

test('DNV combined loading with net internal pressure at or above burst: finite utilisation and burstExceeded', () => {
  const pb = C.burstPressureDnv({ D, t, fy, fu });
  const g = 1.14 * 1.15;
  // 20 % above the burst resistance with some bending and tension: uc = γSC γm (Ted/Tk)² + 1.2² (no bending term)
  const over = C.combinedLoadingDnv({ D, t, fy, fu, pld: 1.2 * pb, MF: 1e5, TeF: 5e5 });
  assert.equal(over.burstExceeded, true); assert.equal(over.mode, 'internal');
  assert.ok(Number.isFinite(over.uc));
  assert.ok(rel(over.uc, g * (1.1 * 5e5 / over.Tk) ** 2 + 1.44) < 1e-12, `uc ${over.uc}`);
  for (const v of Object.values(over)) if (typeof v === 'number') assert.ok(Number.isFinite(v));
  // exactly at burst with no tension: uc = 1; the expression is continuous from below
  assert.ok(rel(C.combinedLoadingDnv({ D, t, fy, fu, pld: pb, MF: 1e5 }).uc, 1) < 1e-12);
  const below = C.combinedLoadingDnv({ D, t, fy, fu, pld: 0.999999 * pb, MF: 1e5 });
  assert.equal(below.burstExceeded, false);
  assert.ok(Math.abs(below.uc - 1) < 5e-3, `just below burst ${below.uc}`);
  // utilisation grows with the overpressure; external pressure reduces the net pressure
  assert.ok(C.combinedLoadingDnv({ D, t, fy, fu, pld: 2 * pb }).uc > over.uc);
  assert.equal(C.combinedLoadingDnv({ D, t, fy, fu, pld: 1.2 * pb, pe: 0.5 * pb }).burstExceeded, false);
  assert.equal(C.combinedLoadingDnv({ D, t, fy, fu, MF: 1e5 }).burstExceeded, false);
});
