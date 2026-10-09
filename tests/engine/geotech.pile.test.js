import test from 'node:test';
import assert from 'node:assert/strict';
import {
  solveLateralPile, pileHeadCurves, pileHeadStiffness, equivalentCantilever, ultimateLateralCapacity,
  criticalPileLength, steelUtilisation, solveAxialPile, axialLoadSettlement, hetenyiSemiInfinite, bromsClay,
  bromsSand, randolphLateral, randolphWrothAxial,
} from '../../src/engine/geotech/pile.js';
import { pyUser, pyLinear, zeroCurve } from '../../src/engine/geotech/py.js';
import { axialCapacity } from '../../src/engine/geotech/axialcapacity.js';
import { createProfile, exampleProfile } from '../../src/engine/geotech/soilprofile.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-12), `${a} vs ${b}`);

const EI = 1e9, k = 2e7; // beta = (k/4EI)^0.25 = 0.2659 1/m → beta L = 10.6 for L = 40 m
const linear = (extra = {}) => ({ L: 40, nElem: 320, EI, pyAt: () => k, ...extra });

test('Hetenyi semi-infinite beam on elastic foundation: y0 = 2H beta/k + 2M beta^2/k, theta0 = 2H beta^2/k + 4M beta^3/k', () => {
  const H = 1e5, M = 2e5;
  const b = Math.pow(k / (4 * EI), 0.25);
  const ref = hetenyiSemiInfinite({ H, M, k, EI, z: [0, 2, 5, 10] });
  close(ref.y0, (2 * H * b) / k + (2 * M * b * b) / k, 1e-12);
  close(ref.theta0, (2 * H * b * b) / k + (4 * M * b ** 3) / k, 1e-12);
  const r = solveLateralPile(linear({ head: { H, M } }));
  close(r.head.y, ref.y0, 1e-3);
  close(r.head.theta, ref.theta0, 1e-3);
  close(r.head.H, H, 1e-8); close(r.head.M, M, 1e-8);
  // distributions at z = 2, 5, 10 m (nodes every 0.125 m)
  for (const [j, zz] of [[1, 2], [2, 5], [3, 10]]) {
    const i = Math.round(zz / 0.125);
    assert.ok(Math.abs(r.y[i] - ref.y[j]) < 2e-3 * ref.y0);
    assert.ok(Math.abs(r.moment[i] - ref.moment[j]) < 3e-3 * r.maxMoment);
  }
  // H only: Mmax = 0.3224 H / beta at beta z = pi/4
  const h = solveLateralPile(linear({ head: { H } }));
  close(h.maxMoment, (0.3224 * H) / b, 2e-3);
  close(h.depthMaxMoment, Math.PI / 4 / b, 0.05);
  // soil reaction integrates to the applied shear (horizontal equilibrium)
  let sum = 0;
  for (let i = 0; i < h.z.length; i++) sum += h.p[i] * (i === 0 || i === h.z.length - 1 ? 0.0625 : 0.125);
  close(sum, H, 1e-6);
});

test('BNWF convergence with element count: observed order 2 (lumped nodal springs)', () => {
  const ref = hetenyiSemiInfinite({ H: 1e5, M: 2e5, k, EI });
  const err = [20, 40, 80, 160].map((n) => Math.abs(solveLateralPile(linear({ nElem: n, head: { H: 1e5, M: 2e5 } })).head.y / ref.y0 - 1));
  const orders = [1, 2, 3].map((i) => Math.log2(err[i - 1] / err[i]));
  for (let i = 1; i < 4; i++) assert.ok(err[i] < err[i - 1]);
  assert.ok(Math.abs(orders[2] - 2) < 0.1, `observed orders ${orders.map((o) => o.toFixed(2))}`);
  assert.ok(err[3] < 2e-3);
});

test('head conditions: fixed head y = H beta/k; prescribed displacement and springs are consistent', () => {
  const H = 1e5, b = Math.pow(k / (4 * EI), 0.25);
  const fx = solveLateralPile(linear({ head: { H, fixedRotation: true } }));
  assert.ok(Math.abs(fx.head.theta) < 1e-15);
  close(fx.head.y, (H * b) / k, 1e-3);
  close(fx.head.M, -H / (2 * b), 2e-3); // restraining moment
  const ks = hetenyiSemiInfinite({ k, EI });
  close(ks.KHH, 4 * EI * b ** 3, 1e-12); close(ks.KHM, -2 * EI * b * b, 1e-12); close(ks.KMM, 2 * EI * b, 1e-12);
  close(H / fx.head.y, ks.KHH, 1e-3);
  // prescribing the free-head response returns the free-head loads as reactions
  const free = solveLateralPile(linear({ head: { H, M: 3e5 } }));
  const pr = solveLateralPile(linear({ head: { y: free.head.y, theta: free.head.theta } }));
  close(pr.head.H, H, 1e-6); close(pr.head.M, 3e5, 1e-6);
  // a rotational head spring between free (0) and fixed (infinite) gives an intermediate deflection; H is shared
  const sp = solveLateralPile(linear({ head: { H, kRot: ks.KMM, kTrans: ks.KHH } }));
  assert.ok(sp.head.y < free.head.y && sp.head.y > 0);
  close(sp.head.H + ks.KHH * sp.head.y, H, 1e-6);
  assert.ok(Math.abs(sp.head.M + ks.KMM * sp.head.theta) < 1e-3);
});

test('cantilever checks: H L^3/3EI, Timoshenko shear term, and P-Delta amplification (tan kL - kL)/(k^3 EI)', () => {
  // 10 m pile without soil, clamped at the toe by very stiff base springs, loaded at the head
  const Ls = 10, D = 1, t = 0.02, E = 210e9;
  const A = (Math.PI * (1 - 0.96 ** 2)) / 4, I = (Math.PI * (1 - 0.96 ** 4)) / 64, H = 1e5;
  const base = { L: Ls, nElem: 40, D, t, E, pyAt: () => zeroCurve(), base: { shear: 1e16, moment: 1e16 }, head: { H } };
  const eb = solveLateralPile(base);
  const dEB = (H * Ls ** 3) / (3 * E * I);
  close(eb.head.y, dEB, 1e-6);
  close(eb.head.theta, (H * Ls ** 2) / (2 * E * I), 1e-6);
  close(eb.moment[eb.moment.length - 1], H * Ls, 1e-6); // clamping moment
  close(eb.shear[20], H, 1e-6);
  const G = E / 2.6, kappa = 0.5;
  const tim = solveLateralPile({ ...base, timoshenko: true, shearCoeff: kappa });
  close(tim.head.y - eb.head.y, (H * Ls) / (kappa * G * A), 1e-5);
  // axial compression P = 0.4 P_Euler (cantilever P_E = pi^2 EI / 4L^2)
  const P = (0.4 * Math.PI ** 2 * E * I) / (4 * Ls * Ls), kk = Math.sqrt(P / (E * I));
  const pd = solveLateralPile({ ...base, axialLoad: P });
  close(pd.head.y, (H * (Math.tan(kk * Ls) - kk * Ls)) / (kk ** 3 * E * I), 1e-4);
  assert.ok(pd.head.y > 1.5 * dEB);
});

const nlOpts = () => ({ L: 40, D: 2, t: 0.05, profile: exampleProfile('layeredSandClay'), head: { H: 2e6, M: 1e7 } });

test('pile-head stiffness matrix: symmetric, positive definite, secant reproduces the working load', () => {
  const r = pileHeadStiffness(nlOpts());
  for (const K of [r.tangent, r.secant]) {
    assert.ok(K.asymmetry < 1e-9); // flexibility cross terms from the two unit-load solves agree
    assert.ok(K.KHH > 0 && K.KMM > 0 && K.KHM < 0);
    assert.ok(K.KHH * K.KMM - K.KHM ** 2 > 0);
    close(K.KHH * K.fHH + K.KHM * K.fHM, 1, 1e-9);
  }
  const { y, theta, H, M } = r.working;
  close(r.secant.KHH * y + r.secant.KHM * theta, H, 1e-6);
  close(r.secant.KHM * y + r.secant.KMM * theta, M, 1e-6);
  assert.ok(r.tangent.KHH < r.secant.KHH); // softening soil
  // tangent matrix equals the finite-difference slope of the head response
  const dH = 2e3;
  const a = solveLateralPile({ ...nlOpts(), head: { H: 2e6 + dH, M: 1e7 } }), b = solveLateralPile({ ...nlOpts(), head: { H: 2e6 - dH, M: 1e7 } });
  close((a.head.y - b.head.y) / (2 * dH), r.tangent.fHH, 1e-3);
  close((a.head.theta - b.head.theta) / (2 * dH), r.tangent.fHM, 1e-3);
  // linear soil: both matrices equal Hetenyi
  const lin = pileHeadStiffness(linear({ head: { H: 1e5 } }));
  const ref = hetenyiSemiInfinite({ k, EI });
  close(lin.tangent.KHH, ref.KHH, 3e-3); close(lin.secant.KHM, ref.KHM, 3e-3); close(lin.tangent.KMM, ref.KMM, 3e-3);
});

test('equivalent cantilever lengths recover a real cantilever', () => {
  const L = 12, ei = 5e9;
  const K = { KHH: (12 * ei) / L ** 3, KHM: (-6 * ei) / L ** 2, KMM: (4 * ei) / L };
  const e = equivalentCantilever({ ...K, EI: ei });
  for (const v of [e.lateral, e.rotation, e.coupling, e.fixedHead]) close(v, L, 1e-9);
});

test('nonlinear response: equilibrium, mobilisation <= 1, head curves soften, cyclic softer than static', () => {
  const r = solveLateralPile(nlOpts());
  assert.equal(r.convergence.converged, true);
  assert.ok(r.convergence.residual < 1e-8);
  const h = r.z[1] - r.z[0];
  let sum = 0, mom = 0;
  for (let i = 0; i < r.z.length; i++) { const w = (i === 0 || i === r.z.length - 1 ? 0.5 : 1) * h; sum += r.p[i] * w; mom += r.p[i] * w * r.z[i]; }
  close(sum, 2e6, 1e-6); // horizontal equilibrium
  close(mom, -1e7, 1e-6); // moment equilibrium about the head: M + ∫ p z dz = 0
  assert.ok(r.maxMobilisation <= 1 + 1e-9);
  assert.ok(Math.abs(r.toe.y) < 1e-3 * r.head.y);
  const c = pileHeadCurves(nlOpts(), { factors: [0.25, 0.5, 0.75, 1] });
  const sec = c.H.map((H, i) => H / c.y[i]);
  for (let i = 1; i < 4; i++) assert.ok(sec[i] < sec[i - 1]);
  close(c.y[3], r.head.y, 1e-6);
  // cyclic degradation only acts beyond 3 yc: use a clay with small eps50 so that the head region is past it
  const softClay = createProfile([{ type: 'clay', top: 0, thickness: 60, gammaSub: 6000, suTop: 2000, suGrad: 1500, eps50: 0.005 }]);
  const soft = { L: 40, D: 2, t: 0.05, profile: softClay, head: { H: 2e6, M: 1e7 } };
  const st = solveLateralPile(soft), cy = solveLateralPile({ ...soft, loading: 'cyclic' });
  assert.ok(st.head.y > 3 * 2.5 * 0.005 * 2);
  assert.ok(cy.head.y > 1.05 * st.head.y, `${cy.head.y} vs ${st.head.y}`);
});

test('ultimate lateral capacity by load stepping reproduces Broms (1964) for a short rigid pile in clay', () => {
  const su = 50e3, D = 1, L = 8, e = 2;
  const br = bromsClay({ su, D, L, e });
  // hand check of the Broms equations: f = H/(9 su D), g = L − 1.5D − f, H (e + 1.5D + f/2) = 2.25 su D g²
  const f = br.Hu / (9 * su * D), g = L - 1.5 * D - f;
  close(br.Hu * (e + 1.5 * D + 0.5 * f), 2.25 * su * D * g * g, 1e-9);
  // same idealisation in the Winkler model: rigid-plastic 9 su D below 1.5 D, nothing above, very stiff pile
  const u = ultimateLateralCapacity(
    { L, nElem: 160, EI: 1e11, pyAt: (z) => (z < 1.5 * D ? zeroCurve() : pyUser({ y: [0, 1e-4], p: [0, 9 * su * D] })) },
    { eccentricity: e, yLimit: 10, nSteps: 60, tol: 1e-3 },
  );
  // approaching collapse the deflection grows without bound, so either stopping criterion marks the limit load
  assert.ok(['non-convergence', 'displacement limit'].includes(u.criterion));
  close(u.Hult, br.Hu, 0.02);
  // long-pile mode governs when the yield moment is low
  const long = bromsClay({ su, D, L: 30, e, My: 1e6 });
  assert.equal(long.mode, 'long');
  close(long.Hu * (e + 1.5 * D + long.Hu / (18 * su * D)), 1e6, 1e-9);
  const bs = bromsSand({ gammaSub: 10e3, phiDeg: 30, D: 1, L: 10, e: 0 });
  close(bs.Hu, 0.5 * 10e3 * 1 * 100 * 3, 1e-12); // Kp = 3: 0.5 gamma D L^3 Kp / L
});

test('critical pile length of an elastic pile is a few 1/beta; elastic continuum closed forms are self-consistent', () => {
  const b = Math.pow(k / (4 * EI), 0.25);
  const c = criticalPileLength({ EI, pyAt: () => k, head: { H: 1e5 }, nElem: 200 }, { Lmin: 4, Lmax: 40, nLengths: 37, tol: 0.02 });
  assert.ok(c.Lcrit * b > 2 && c.Lcrit * b < 4, `beta Lcrit = ${c.Lcrit * b}`);
  // Randolph (1981), homogeneous soil: lc = D (Ep/G*)^(2/7), rho_c = 1
  const D = 2, G0 = 30e6, nu = 0.3, ei = 2e11;
  const r = randolphLateral({ H: 1e6, M: 5e6, EI: ei, D, G0, nu });
  const Ep = ei / ((Math.PI * 1) / 4), Gs = G0 * (1 + 0.75 * nu);
  close(r.lc, D * Math.pow(Ep / Gs, 2 / 7), 1e-9);
  close(r.rhoc, 1, 1e-12);
  const c0 = Math.pow(Ep / Gs, 1 / 7) / Gs;
  close(r.y0, c0 * ((0.27 * 1e6) / (r.lc / 2) + (0.3 * 5e6) / (r.lc / 2) ** 2), 1e-9);
  close(r.KHH * r.y0 + r.KHM * r.theta0, 1e6, 1e-9);
  const g = randolphLateral({ EI: ei, D, G0: 0, Ggrad: 2e6, nu });
  close(g.rhoc, 0.5, 1e-9); // Gibson soil
});

test('steel utilisation: axial plus bending stress by hand', () => {
  const D = 2, t = 0.05, A = (Math.PI * (4 - 1.9 ** 2)) / 4, I = (Math.PI * (16 - 1.9 ** 4)) / 64;
  const u = steelUtilisation({ N: 10e6, M: 50e6, D, t, fy: 355e6, usageFactor: 0.6 });
  close(u.sigmaAxial, 10e6 / A, 1e-12); close(u.sigmaBending, 50e6 / I, 1e-12);
  close(u.utilisation, (10e6 / A + 50e6 / I) / (0.6 * 355e6), 1e-12);
  const r = solveLateralPile(nlOpts());
  close(Math.max(...r.sigmaBending), (r.maxMoment * D) / (2 * I), 1e-9);
});

const sandProfile = () => createProfile([{ type: 'sand', top: 0, thickness: 60, gammaSub: 10e3, phiDeg: 36, apiSand: 'dense sand' }]);

test('t–z analysis, rigid-pile limit: capacity equals the summed shaft and tip resistance', () => {
  const o = { profile: sandProfile(), L: 30, D: 1.5, t: 0.04, nElem: 120, E: 210e14 };
  const cap = axialCapacity({ profile: o.profile, L: 30, D: 1.5, t: 0.04 });
  const c = axialLoadSettlement(o, { wMax: 0.2, steps: 40 });
  close(c.capacity, c.QsUltimate + c.QpUltimate, 1e-6);
  close(c.QpUltimate, cap.Qp, 1e-12);
  close(c.QsUltimate, cap.Qfo, 2e-3); // nodal (trapezoidal) lumping versus Gauss integration
  close(c.capacity, cap.Qc, 2e-3);
  // tension: shaft only
  const tn = axialLoadSettlement({ ...o, mode: 'tension' }, { wMax: 0.05, steps: 20 });
  close(tn.capacity, tn.QsUltimate, 1e-6);
  assert.equal(tn.QpUltimate, 0);
  // rigid pile: full shaft mobilisation at 2.5 mm, tip follows the Q–z table (0.25 Qp at 0.002 D)
  const s = solveAxialPile({ ...o, headDisplacement: 0.003 });
  close(s.headLoad, c.QsUltimate + 0.25 * c.QpUltimate, 2e-5);
});

test('t–z analysis, elastic pile with linear springs: closed-form head stiffness', () => {
  // bar on uniform shaft springs ks [N/m per m] with a tip spring Kb:
  // K = EA mu (Kb + EA mu tanh(mu L)) / (EA mu + Kb tanh(mu L)), mu = sqrt(ks/EA)
  const D = 1, t = 0.03, L = 40, E = 210e9, EA = (E * Math.PI * (1 - 0.94 ** 2)) / 4;
  const kUnit = 40e6, ks = kUnit * Math.PI * D, Kb = 5e8; // kUnit: shaft stiffness per unit area [Pa/m]
  const mu = Math.sqrt(ks / EA), th = Math.tanh(mu * L);
  const Kref = (EA * mu * (Kb + EA * mu * th)) / (EA * mu + Kb * th);
  const err = [25, 50, 100, 200].map((n) => {
    const r = solveAxialPile({ L, D, t, E, nElem: n, tzAt: () => pyLinear(kUnit), tipCurve: pyLinear(Kb), headLoad: 1e6 });
    return Math.abs(1e6 / r.headDisplacement / Kref - 1);
  });
  assert.ok(err[3] < 1e-4, `errors ${err}`);
  assert.ok(Math.abs(Math.log2(err[1] / err[2]) - 2) < 0.2); // second order
  // axial force decays from the head load to the tip load
  const r = solveAxialPile({ L, D, t, E, nElem: 100, tzAt: () => pyLinear(kUnit), tipCurve: pyLinear(Kb), headLoad: 1e6 });
  close(r.tipLoad, Kb * r.tipDisplacement, 1e-12);
  assert.ok(r.force[0] > r.force[50] && r.force[50] > r.force[100]);
});

test('t–z analysis with API curves: load control agrees with displacement control; clay softens after peak', () => {
  const profile = exampleProfile('layeredSandClay'), o = { profile, L: 40, D: 2, t: 0.05 };
  const c = axialLoadSettlement(o);
  assert.ok(c.capacity < c.QsUltimate + c.QpUltimate); // clay friction is at residual when the tip is mobilised
  assert.ok(c.capacity > 0.9 * (c.QsUltimate + c.QpUltimate));
  const a = solveAxialPile({ ...o, headLoad: 0.5 * c.capacity });
  const b = solveAxialPile({ ...o, headDisplacement: a.headDisplacement });
  close(b.headLoad, 0.5 * c.capacity, 1e-6);
  assert.throws(() => solveAxialPile({ ...o, headLoad: 1.2 * (c.QsUltimate + c.QpUltimate) }), /capacity/);
  assert.throws(() => solveAxialPile({ ...o }), RangeError);
  // Randolph & Wroth closed form: positive, and stiffer for a stiffer soil
  const k1 = randolphWrothAxial({ L: 40, D: 2, EA: 6.4e10, GL: 50e6, nu: 0.3 }).stiffness;
  const k2 = randolphWrothAxial({ L: 40, D: 2, EA: 6.4e10, GL: 100e6, nu: 0.3 }).stiffness;
  assert.ok(k2 > k1 && k1 > 0);
  // rigid-pile limit of the same formula: K → G r0 [4/(1−nu) + 2 pi L/(zeta r0)]
  const rg = randolphWrothAxial({ L: 40, D: 2, EA: 1e20, GL: 50e6, nu: 0.3 });
  close(rg.stiffness, 50e6 * 1 * (4 / 0.7 + (2 * Math.PI * 40) / rg.zeta), 1e-6);
});
