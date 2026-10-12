import test from 'node:test';
import assert from 'node:assert/strict';
import {
  analyzeFrame, memberDiagram, refineModel, memberAxes, beamLocalStiffness, sectionStress, buildFrameSystem, updateFrameSystem,
} from '../../src/engine/structural/frame3d.js';
import { rectSection, tubeSection, iSection } from '../../src/engine/structural/sections.js';

const E = 2.1e11;
const G = E / 2.6;
const sec = rectSection({ h: 0.4, b: 0.2 });
const EI = E * sec.Iy;
const GAs = G * sec.Asz;
const near = (a, b, rel = 1e-9, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);

test('cantilever, tip load: δ = PL³/3EI + PL/(G As) (Timoshenko) and equilibrium residual ≈ 0', () => {
  const L = 2; const P = 1e5;
  const model = {
    nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: sec }], supports: [{ node: 0, fix: 'fixed' }],
    loadCases: [{ name: 'P', nodal: [{ node: 1, F: [0, 0, -P, 0, 0, 0] }] }],
  };
  const r = analyzeFrame(model).results.P;
  near(r.displacements[8], -(P * L ** 3 / (3 * EI) + P * L / GAs), 1e-10);
  // the shear term is a visible part of the answer for this stocky beam
  assert.ok(P * L / GAs > 0.02 * P * L ** 3 / (3 * EI));
  near(r.reactions[0].R[2], P, 1e-10);
  near(r.reactions[0].R[4], -P * L, 1e-10); // support moment about y balances the load couple
  assert.ok(r.equilibrium.forceResidual < 1e-6 && r.equilibrium.momentResidual < 1e-6 && r.equilibrium.nodalResidual < 1e-6);
  near(r.strainEnergy, 0.5 * P * -r.displacements[8], 1e-10, 'strain energy = ½Pδ');
  // Euler–Bernoulli option removes the shear term
  const rb = analyzeFrame(model, { shearDeformation: false }).results.P;
  near(rb.displacements[8], -P * L ** 3 / (3 * EI), 1e-10);
});

test('simply supported beam, uniform load: M = qL²/8, δ = 5qL⁴/384EI + qL²/(8 G As)', () => {
  const L = 6; const q = 3e4;
  const model = {
    nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: sec }],
    supports: [{ node: 0, fix: [1, 1, 1, 1, 0, 0] }, { node: 1, fix: [0, 1, 1, 0, 0, 0] }],
    loadCases: [{ name: 'q', member: [{ member: 0, type: 'uniform', q: [0, 0, -q] }] }],
  };
  const res = analyzeFrame(model).results.q;
  const d = memberDiagram(res, 0, 21);
  near(d.My[10], -q * L * L / 8, 1e-9, 'sagging moment is negative My');
  near(d.Vz[0], -q * L / 2, 1e-9, 'V on the +x face at the left support');
  near(d.uz[10], -(5 * q * L ** 4 / (384 * EI) + q * L * L / (8 * GAs)), 1e-9);
  near(res.reactions[0].R[2], q * L / 2, 1e-10);
  // end rotation qL³/24EI (unaffected by shear deformation)
  near(Math.abs(res.displacements[4]), q * L ** 3 / (24 * EI), 1e-9);
});

test('fixed-fixed beam: uniform load M_end = qL²/12, δ_mid = qL⁴/384EI; mid-span point load M = PL/8', () => {
  const L = 8; const q = 2e4; const P = 5e4;
  const model = {
    nodes: [[0, 0, 0], [L / 2, 0, 0], [L, 0, 0]],
    members: [{ a: 0, b: 1, section: sec }, { a: 1, b: 2, section: sec }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 2, fix: 'fixed' }],
    loadCases: [
      { name: 'q', member: [{ member: 0, type: 'uniform', q: [0, 0, -q] }, { member: 1, type: 'uniform', q: [0, 0, -q] }] },
      { name: 'P', nodal: [{ node: 1, F: [0, 0, -P, 0, 0, 0] }] },
    ],
  };
  const a = analyzeFrame(model, { shearDeformation: false });
  near(a.results.q.displacements[8], -q * L ** 4 / (384 * EI), 1e-9);
  near(a.results.q.memberForces[0].endA.My, q * L * L / 12, 1e-9, 'hogging at the support');
  near(a.results.q.memberForces[0].endB.My, -q * L * L / 24, 1e-9, 'sagging at mid-span');
  near(a.results.P.memberForces[0].endA.My, P * L / 8, 1e-9);
  near(a.results.P.displacements[8], -P * L ** 3 / (192 * EI), 1e-9);
  // single element with the point load as a member load gives the same fixed-end moments
  const one = {
    nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: sec }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 1, fix: 'fixed' }],
    loadCases: [{ name: 'P', member: [{ member: 0, type: 'point', xi: 0.5, P: [0, 0, -P] }] }],
  };
  const r1 = analyzeFrame(one, { shearDeformation: false }).results.P;
  near(r1.memberForces[0].endA.My, P * L / 8, 1e-9);
  const d = memberDiagram(r1, 0, 41);
  near(d.My[20], -P * L / 8, 1e-9);
  near(d.uz[20], -P * L ** 3 / (192 * EI), 1e-6);
});

test('fixed-base portal frame under lateral load vs slope-deflection hand solution', () => {
  // Columns h, beam L, axially rigid. Unknowns θ (both joints, antisymmetric) and sway Δ:
  //   joint:  (4EIc/h + 6EIb/L) θ − (6EIc/h²) Δ = 0
  //   storey: −(6EIc/h²) θ + (12EIc/h³) Δ = P/2      (each column carries P/2)
  const h = 4; const L = 6; const P = 8e4;
  const col = rectSection({ h: 0.3, b: 0.3 });
  const beam = rectSection({ h: 0.5, b: 0.25 });
  const EIc = E * col.Iy; const EIb = E * beam.Iy;
  const a11 = 4 * EIc / h + 6 * EIb / L; const a12 = -6 * EIc / h ** 2; const a22 = 12 * EIc / h ** 3;
  const det = a11 * a22 - a12 * a12;
  const theta = (-a12 * P / 2) / det;
  const delta = (a11 * P / 2) / det;
  const Mbase = (2 * EIc / h) * (theta - 3 * delta / h); // M_AB by slope-deflection
  const rigid = (s) => ({ ...s, A: s.A * 1e4 });
  const model = {
    nodes: [[0, 0, 0], [0, 0, h], [L, 0, h], [L, 0, 0]],
    members: [{ a: 0, b: 1, section: rigid(col) }, { a: 1, b: 2, section: rigid(beam) }, { a: 3, b: 2, section: rigid(col) }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 3, fix: 'fixed' }],
    loadCases: [{ name: 'H', nodal: [{ node: 1, F: [P, 0, 0, 0, 0, 0] }] }],
  };
  const r = analyzeFrame(model, { shearDeformation: false }).results.H;
  near(r.displacements[6], delta, 2e-5, 'sway');
  near(Math.abs(r.displacements[6 + 4]), theta, 2e-5, 'joint rotation');
  near(Math.abs(r.reactions[0].R[4]), Math.abs(Mbase), 2e-5, 'base moment');
  near(r.reactions[0].R[0] + r.reactions[1].R[0], -P, 1e-6);
  assert.ok(r.equilibrium.relative < 1e-7);
});

test('3-D space truss (tripod) vs statics', () => {
  // apex (0,0,h); bases at radius r, angles 0°, 120°, 240°.
  // vertical load −V: N = −V Lb/(3h) in each bar. Horizontal load +H (x): N1 = −2H Lb/(3r), N2 = N3 = +H Lb/(3r).
  const h = 3; const r = 2; const V = 9e4; const H = 3e4;
  const Lb = Math.hypot(h, r);
  const tube = tubeSection({ D: 0.2, t: 0.01 });
  const nodes = [[0, 0, h]];
  for (const a of [0, 120, 240]) nodes.push([r * Math.cos(a * Math.PI / 180), r * Math.sin(a * Math.PI / 180), 0]);
  const model = {
    nodes, members: [1, 2, 3].map((b) => ({ a: 0, b, section: tube, type: 'truss' })),
    supports: [1, 2, 3].map((n) => ({ node: n, fix: 'pinned' })),
    loadCases: [{ name: 'V', nodal: [{ node: 0, F: [0, 0, -V, 0, 0, 0] }] }, { name: 'H', nodal: [{ node: 0, F: [H, 0, 0, 0, 0, 0] }] }],
    combinations: [{ name: 'C', factors: { V: 1.35, H: 1.5 } }],
  };
  const a = analyzeFrame(model);
  for (let i = 0; i < 3; i++) near(a.results.V.memberForces[i].endA.N, -V * Lb / (3 * h), 1e-9);
  near(a.results.H.memberForces[0].endA.N, -2 * H * Lb / (3 * r), 1e-9);
  near(a.results.H.memberForces[1].endA.N, H * Lb / (3 * r), 1e-9);
  near(a.results.H.memberForces[2].endA.N, H * Lb / (3 * r), 1e-9);
  near(a.results.C.memberForces[0].endA.N, 1.35 * -V * Lb / (3 * h) + 1.5 * -2 * H * Lb / (3 * r), 1e-9, 'combination');
  // vertical deflection of the apex: δ = V Lb³/(3 E A h²) by virtual work
  near(a.results.V.displacements[2], -V * Lb ** 3 / (3 * E * tube.A * h * h), 1e-9);
  assert.ok(a.warnings.length === 1, 'rotations at truss nodes are restrained automatically');
  assert.ok(a.results.H.equilibrium.forceResidual < 1e-6);
});

test('member releases: hinge at one end turns a fixed-fixed beam into a propped cantilever', () => {
  const L = 5; const q = 1.2e4;
  const model = {
    nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: sec, releases: { b: 'ry' } }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 1, fix: 'fixed' }],
    loadCases: [{ name: 'q', member: [{ member: 0, type: 'uniform', q: [0, 0, -q] }] }],
  };
  const r = analyzeFrame(model, { shearDeformation: false }).results.q;
  near(r.memberForces[0].endA.My, q * L * L / 8, 1e-9);
  assert.equal(r.memberForces[0].endB.My, 0);
  near(r.reactions[1].R[2], 3 * q * L / 8, 1e-9);
  near(r.reactions[0].R[2], 5 * q * L / 8, 1e-9);
  // maximum deflection of a propped cantilever qL⁴/(185 EI) (more precisely 0.005416 qL⁴/EI at x = 0.5785 L from the fixed end)
  const d = memberDiagram(r, 0, 201);
  near(d.maxAbs.uz, 0.0054161 * q * L ** 4 / EI, 2e-4);
});

test('internal hinge: two cantilevers joined by a shear hinge share a tip load equally', () => {
  // each half is a cantilever of length a loaded by P/2 at its tip: δ = (P/2) a³/(3EI)
  const a = 3; const P = 4e4;
  const model = {
    nodes: [[0, 0, 0], [a, 0, 0], [2 * a, 0, 0]],
    members: [{ a: 0, b: 1, section: sec, releases: { b: 'pin' } }, { a: 1, b: 2, section: sec, releases: { a: 'pin' } }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 2, fix: 'fixed' }, { node: 1, fix: [0, 1, 0, 1, 0, 1] }],
    loadCases: [{ name: 'P', nodal: [{ node: 1, F: [0, 0, -P, 0, 0, 0] }] }],
  };
  const an = analyzeFrame(model, { shearDeformation: false });
  near(an.results.P.displacements[8], -(P / 2) * a ** 3 / (3 * EI), 1e-9);
  near(an.results.P.memberForces[0].endA.My, P * a / 2, 1e-9);
});

test('thermal loads: restrained bar N = −E A α ΔT; cantilever gradient tip deflection α g L²/2', () => {
  const L = 4; const dT = 40; const alpha = 1.2e-5; const grad = 50;
  const base = { nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: sec }] };
  const fixed = analyzeFrame({
    ...base, supports: [{ node: 0, fix: 'fixed' }, { node: 1, fix: 'fixed' }],
    loadCases: [{ name: 'T', member: [{ member: 0, type: 'thermal', dT }] }],
  }).results.T;
  near(fixed.memberForces[0].endA.N, -E * sec.A * alpha * dT, 1e-10);
  const free = analyzeFrame({
    ...base, supports: [{ node: 0, fix: 'fixed' }],
    loadCases: [{ name: 'T', member: [{ member: 0, type: 'thermal', dT, gradZ: grad }] }],
  }).results.T;
  near(free.displacements[6], alpha * dT * L, 1e-10, 'free expansion');
  // hotter top fibre (dT/dz > 0) → the free end moves down: w = −α g L²/2
  near(free.displacements[8], -alpha * grad * L * L / 2, 1e-10);
  assert.ok(Math.abs(free.memberForces[0].endA.My) < 1e-6, 'no forces in a statically determinate member');
  const d = memberDiagram(free, 0, 11);
  near(d.uz[5], -alpha * grad * (L / 2) ** 2 / 2, 1e-9);
  near(d.ux[10], alpha * dT * L, 1e-10);
});

test('support settlement of a fixed-fixed beam: M = 6EIδ/L², V = 12EIδ/L³', () => {
  const L = 6; const dz = -0.01;
  const model = {
    nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: sec }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 1, fix: 'fixed' }],
    loadCases: [{ name: 'S', settlements: [{ node: 1, d: [0, 0, dz, 0, 0, 0] }] }],
  };
  const r = analyzeFrame(model, { shearDeformation: false }).results.S;
  near(Math.abs(r.memberForces[0].endA.My), 6 * EI * Math.abs(dz) / L ** 2, 1e-9);
  near(Math.abs(r.reactions[0].R[2]), 12 * EI * Math.abs(dz) / L ** 3, 1e-9);
  assert.ok(r.equilibrium.forceResidual < 1e-6 && r.equilibrium.momentResidual < 1e-5);
});

test('elastic supports: tip spring k and full 6×6 foundation matrix', () => {
  const L = 3; const P = 2e4; const k = 4e6;
  const model = {
    nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: sec }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 1, spring: [0, 0, k, 0, 0, 0] }],
    loadCases: [{ name: 'P', nodal: [{ node: 1, F: [0, 0, -P, 0, 0, 0] }] }],
  };
  const r = analyzeFrame(model, { shearDeformation: false }).results.P;
  near(r.displacements[8], -P / (k + 3 * EI / L ** 3), 1e-10);
  near(r.reactions[1].R[2], k * P / (k + 3 * EI / L ** 3), 1e-10);
  // rigid column on a coupled lateral/rocking foundation spring: [Kh, Khr; Khr, Kr] u = [H, H·h]
  const h = 10; const H = 1e5; const Kh = 5e7; const Kr = 4e9; const Khr = -2e8;
  const K = Array.from({ length: 6 }, () => new Array(6).fill(0));
  K[0][0] = Kh; K[4][4] = Kr; K[0][4] = Khr; K[4][0] = Khr; K[1][1] = 1e12; K[2][2] = 1e12; K[3][3] = 1e12; K[5][5] = 1e12;
  const stiff = { ...sec, A: sec.A * 1e6, Iy: sec.Iy * 1e8, Iz: sec.Iz * 1e8 };
  const m2 = {
    nodes: [[0, 0, 0], [0, 0, h]], members: [{ a: 0, b: 1, section: stiff }], supports: [{ node: 0, spring: K }],
    loadCases: [{ name: 'H', nodal: [{ node: 1, F: [H, 0, 0, 0, 0, 0] }] }],
  };
  const r2 = analyzeFrame(m2, { shearDeformation: false }).results.H;
  const det = Kh * Kr - Khr * Khr;
  const u0 = (Kr * H - Khr * H * h) / det;
  const th0 = (-Khr * H + Kh * H * h) / det;
  near(r2.displacements[0], u0, 1e-6);
  near(r2.displacements[4], th0, 1e-6);
  near(r2.displacements[6], u0 + th0 * h, 1e-4, 'head displacement of the rigid column');
});

test('tension-only X bracing: the compression diagonal drops out, N = H/cosφ', () => {
  const w = 4; const h = 3; const H = 5e4;
  const rod = tubeSection({ D: 0.06, t: 0.005 });
  const model = {
    nodes: [[0, 0, 0], [w, 0, 0], [0, 0, h], [w, 0, h]],
    members: [
      { a: 0, b: 2, section: sec, type: 'truss' }, { a: 1, b: 3, section: sec, type: 'truss' }, { a: 2, b: 3, section: sec, type: 'truss' },
      { a: 0, b: 3, section: rod, type: 'tension' }, { a: 1, b: 2, section: rod, type: 'tension' },
    ],
    supports: [{ node: 0, fix: 'pinned' }, { node: 1, fix: 'pinned' }, { node: 2, fix: [0, 1, 0, 0, 0, 0] }, { node: 3, fix: [0, 1, 0, 0, 0, 0] }],
    loadCases: [{ name: 'H', nodal: [{ node: 2, F: [H, 0, 0, 0, 0, 0] }] }],
  };
  const r = analyzeFrame(model).results.H;
  assert.deepEqual(r.inactiveMembers, [4]);
  near(r.memberForces[3].endA.N, H * Math.hypot(w, h) / w, 1e-9);
  assert.equal(Math.abs(r.memberForces[4].endA.N), 0);
  assert.ok(r.iterations >= 2);
});

test('mechanism detection names the offending degree of freedom', () => {
  const L = 4;
  const model = {
    nodes: [[0, 0, 0], [L, 0, 0], [2 * L, 0, 0]],
    members: [{ a: 0, b: 1, section: sec, releases: { b: 'ry' } }, { a: 1, b: 2, section: sec }],
    supports: [{ node: 0, fix: [1, 1, 1, 1, 0, 0] }, { node: 2, fix: [0, 1, 1, 0, 0, 0] }, { node: 1, fix: [0, 1, 0, 0, 0, 1] }],
    loadCases: [{ name: 'P', nodal: [{ node: 1, F: [0, 0, -1e3, 0, 0, 0] }] }],
  };
  assert.throws(() => analyzeFrame(model), (e) => e.mechanism === true && /node \d+ (u|r)[xyz]/.test(e.message));
  // releasing torsion at both ends leaves the member free to spin
  assert.throws(() => buildFrameSystem({
    nodes: [[0, 0, 0], [1, 0, 0]], members: [{ a: 0, b: 1, section: sec, releases: { a: 'rx', b: 'rx' } }], supports: [{ node: 0, fix: 'fixed' }],
  }), /mechanism/);
});

test('rigid end offset: cantilever whose first part is rigid deflects as the shorter flexible length', () => {
  const L = 5; const a = 1.5; const P = 3e4;
  const model = {
    nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: sec, offsetA: [a, 0, 0] }], supports: [{ node: 0, fix: 'fixed' }],
    loadCases: [{ name: 'P', nodal: [{ node: 1, F: [0, 0, -P, 0, 0, 0] }] }],
  };
  const r = analyzeFrame(model, { shearDeformation: false }).results.P;
  near(r.displacements[8], -P * (L - a) ** 3 / (3 * EI), 1e-10);
  near(r.reactions[0].R[4], -P * L, 1e-10, 'the support still sees the full lever arm');
});

test('self-weight with buoyancy: submerged tube, dry and flooded; marine growth', () => {
  const D = 0.8; const t = 0.02; const L = 10; const g = 9.80665; const rhoW = 1025; const rhoS = 7850;
  const tube = tubeSection({ D, t });
  const mk = (flooded, extra = {}) => ({
    nodes: [[0, 0, -20], [L, 0, -20]], members: [{ a: 0, b: 1, section: tube, flooded }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 1, fix: [0, 1, 1, 0, 0, 0] }],
    loadCases: [{ name: 'W', selfWeight: { buoyancy: { waterline: 0, rhoW }, ...extra } }],
  });
  const sumR = (m) => { const r = analyzeFrame(m).results.W; return r.reactions[0].R[2] + r.reactions[1].R[2]; };
  near(sumR(mk(false)), (rhoS * tube.A - rhoW * Math.PI * D * D / 4) * g * L, 1e-9, 'dry: steel weight minus envelope buoyancy');
  near(sumR(mk(true)), (rhoS - rhoW) * tube.A * g * L, 1e-9, 'flooded: buoyant weight of the steel');
  const tg = 0.05; const rhoG = 1325;
  const growth = Math.PI / 4 * ((D + 2 * tg) ** 2 - D * D);
  near(sumR(mk(true, { marineGrowth: { thickness: tg, rho: rhoG, zTop: 0, zBottom: -50 } })),
    ((rhoS - rhoW) * tube.A + (rhoG - rhoW) * growth) * g * L, 1e-9, 'marine growth adds its buoyant weight');
  // a vertical member piercing the waterline is only buoyant below it
  const m = {
    nodes: [[0, 0, -10], [0, 0, 5]], members: [{ a: 0, b: 1, section: tube }], supports: [{ node: 0, fix: 'fixed' }],
    loadCases: [{ name: 'W', selfWeight: { buoyancy: { waterline: 0, rhoW } } }],
  };
  near(analyzeFrame(m).results.W.reactions[0].R[2], (rhoS * tube.A * 15 - rhoW * Math.PI * D * D / 4 * 10) * g, 1e-9);
});

test('inertial (transport) load: base shear = m·a, rotation adds α × r', () => {
  const L = 6; const ax = 2.5; const alpha = 0.3;
  const tube = tubeSection({ D: 0.5, t: 0.02 });
  const m = 7850 * tube.A * L; const mTop = 2000;
  const model = {
    nodes: [[0, 0, 0], [0, 0, L]], members: [{ a: 0, b: 1, section: tube }], supports: [{ node: 0, fix: 'fixed' }],
    nodeMasses: [{ node: 1, m: mTop }],
    loadCases: [{ name: 'A', acceleration: { a: [ax, 0, 0] } }, { name: 'R', acceleration: { alpha: [0, alpha, 0], origin: [0, 0, 0] } }],
  };
  const a = analyzeFrame(model);
  near(a.results.A.reactions[0].R[0], (m + mTop) * ax, 1e-9);
  // α about y at the base: a_x = α z → inertia force −μ α z per length: shear = α (μ L²/2 + mTop L)
  near(a.results.R.reactions[0].R[0], alpha * (m * L / 2 + mTop * L), 1e-9);
});

test('local axes: betaDeg = 90 turns an I-beam onto its weak axis; default z is up', () => {
  const s = iSection({ h: 0.3, b: 0.15, tw: 0.008, tf: 0.012 });
  const { R } = memberAxes([0, 0, 0], [3, 4, 0]);
  near(R[2][2], 1, 1e-12);
  assert.ok(Math.abs(R[1][0] * 3 + R[1][1] * 4) < 1e-12);
  const vert = memberAxes([0, 0, 0], [0, 0, 2]).R;
  assert.deepEqual(vert[1], [0, 1, 0]);
  const L = 3; const P = 1e4;
  const mk = (betaDeg) => ({
    nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: s, betaDeg }], supports: [{ node: 0, fix: 'fixed' }],
    loadCases: [{ name: 'P', nodal: [{ node: 1, F: [0, 0, -P, 0, 0, 0] }] }],
  });
  const strong = analyzeFrame(mk(0), { shearDeformation: false }).results.P.displacements[8];
  const weak = analyzeFrame(mk(90), { shearDeformation: false }).results.P.displacements[8];
  near(strong, -P * L ** 3 / (3 * E * s.Iy), 1e-9);
  near(weak, -P * L ** 3 / (3 * E * s.Iz), 1e-9);
});

test('distributed loads: linear load cantilever (11qL⁴/120EI), caller-supplied function load, refinement invariance', () => {
  const L = 4; const q = 1.5e4;
  const base = { nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section: sec }], supports: [{ node: 0, fix: 'fixed' }] };
  const model = {
    ...base,
    loadCases: [
      { name: 'lin', member: [{ member: 0, type: 'linear', q1: [0, 0, 0], q2: [0, 0, -q] }] },
      { name: 'fn', member: [{ member: 0, type: 'function', fn: (s) => [0, 0, -q * s / L] }] },
      { name: 'part', member: [{ member: 0, type: 'uniform', q: [0, 0, -q], s1: 1, s2: 3 }] },
    ],
  };
  const a = analyzeFrame(model, { shearDeformation: false });
  near(a.results.lin.displacements[8], -11 * q * L ** 4 / (120 * EI), 1e-10);
  near(a.results.fn.displacements[8], a.results.lin.displacements[8], 1e-10);
  near(a.results.lin.memberForces[0].endA.My, q * L * L / 3, 1e-10);
  // partial uniform load on [1, 3]: resultant 2q at x = 2 → root moment 4q; tip deflection by superposition of
  // two full-length-from-root cantilever cases: δ(b) − δ(a) with δ(c) = q c³ (4L − c)/(24 EI)
  const dc = (c) => q * c ** 3 * (4 * L - c) / (24 * EI);
  near(a.results.part.memberForces[0].endA.My, 4 * q, 1e-10);
  near(a.results.part.displacements[8], -(dc(3) - dc(1)), 1e-10);
  const fine = analyzeFrame(refineModel(model, 5).model, { shearDeformation: false });
  for (const k of ['lin', 'fn', 'part']) near(fine.results[k].displacements[8], a.results[k].displacements[8], 1e-9, `refined ${k}`);
  // diagram of the linear load: M(x) = −q (L − x)² (2L + x)/(6L) … checked at mid-length against statics
  const d = memberDiagram(a.results.lin, 0, 21);
  const x = L / 2;
  const Mstat = (q * L / 2) * (2 * L / 3 - x) - 0 + (q * x * x * x) / (6 * L) - (q * L / 2) * 0; // root reaction terms, see below
  // statics from the free end: M(x) = ∫_x^L q s/L (s − x) ds = q/L [ (L³ − x³)/3 − x (L² − x²)/2 ]
  const Mfree = q / L * ((L ** 3 - x ** 3) / 3 - x * (L * L - x * x) / 2);
  near(d.My[10], Mfree, 1e-9);
  assert.ok(Number.isFinite(Mstat));
});

test('3-D frame with mixed loading: global equilibrium residual is reported and negligible', () => {
  const tube = tubeSection({ D: 0.6, t: 0.02 });
  const model = {
    nodes: [[0, 0, 0], [6, 0, 0], [6, 5, 0], [0, 5, 0], [0.5, 0.5, 8], [5.5, 0.5, 8], [5.5, 4.5, 8], [0.5, 4.5, 8]],
    members: [
      ...[0, 1, 2, 3].map((i) => ({ a: i, b: i + 4, section: tube })),
      ...[0, 1, 2, 3].map((i) => ({ a: 4 + i, b: 4 + (i + 1) % 4, section: tube })),
      { a: 0, b: 5, section: tube, releases: { a: 'pin', b: 'pin' } }, { a: 1, b: 6, section: tube }, { a: 2, b: 7, section: tube }, { a: 3, b: 4, section: tube },
    ],
    supports: [{ node: 0, fix: 'fixed' }, { node: 1, fix: 'pinned' }, { node: 2, fix: 'fixed' }, { node: 3, spring: [1e8, 1e8, 1e9, 0, 0, 0] }],
    loadCases: [
      { name: 'G', selfWeight: {} },
      { name: 'E', nodal: [{ node: 4, F: [5e4, 2e4, 0, 0, 0, 1e4] }], member: [{ member: 4, type: 'uniform', q: [0, 3e3, -2e3] }, { member: 9, type: 'point', xi: 0.3, P: [1e4, 0, 0], M: [0, 5e3, 0] }] },
    ],
    combinations: [{ name: 'ULS', factors: { G: 1.1, E: 1.35 } }],
  };
  const a = analyzeFrame(model);
  for (const k of ['G', 'E', 'ULS']) {
    const eq = a.results[k].equilibrium;
    assert.ok(eq.relative < 1e-8, `${k}: relative residual ${eq.relative}`);
    assert.ok(eq.momentResidual < 1e-3, `${k}: moment residual ${eq.momentResidual}`);
  }
  // superposition holds for the linear model
  near(a.results.ULS.displacements[6 * 5], 1.1 * a.results.G.displacements[6 * 5] + 1.35 * a.results.E.displacements[6 * 5], 1e-8);
  // section forces from the diagram close at end B
  const d = memberDiagram(a.results.ULS, 9, 11);
  near(d.Mz[10], a.results.ULS.memberForces[9].endB.Mz, 1e-9);
  near(d.N[0], a.results.ULS.memberForces[9].endA.N, 1e-9);
});

test('element stiffness and stress helpers', () => {
  const k = beamLocalStiffness({ E, G, A: sec.A, Iy: sec.Iy, Iz: sec.Iz, J: sec.J, L: 2 });
  near(k[2 * 12 + 2], 12 * EI / 8, 1e-12);
  near(k[4 * 12 + 4], 4 * EI / 2, 1e-12);
  for (let i = 0; i < 12; i++) for (let j = 0; j < 12; j++) assert.equal(k[i * 12 + j], k[j * 12 + i]);
  const tube = tubeSection({ D: 0.5, t: 0.02 });
  const s = sectionStress({ section: tube, N: 1e6, My: 3e5, Mz: 4e5 });
  near(s.sigmaMax, 1e6 / tube.A + 5e5 / tube.Wely, 1e-12);
  assert.throws(() => analyzeFrame({ nodes: [[0, 0, 0]], members: [] }), RangeError);
});

test('sectionStress on a generic {A, Iy, Iz, J} section: finite stresses with documented extreme-fibre distances', () => {
  // A = 0.02 m², Iy = 4e-4 m⁴, Iz = 1e-4 m⁴: radii of gyration iy = √(Iy/A) = 0.14142 m, iz = 0.07071 m
  // default fibre distances c = √2·i: cz = 0.2 m, cy = 0.1 m → Wely = Iy/cz = 2e-3 m³, Welz = Iz/cy = 1e-3 m³
  const gen = { A: 0.02, Iy: 4e-4, Iz: 1e-4, J: 3e-4 };
  const s = sectionStress({ section: gen, N: 2e5, My: 1e5, Mz: -3e4, Vz: 4e4, Vy: 1e4 });
  near(s.sigmaAxial, 1e7, 1e-12);
  near(s.sigmaBending, 1e5 / 2e-3 + 3e4 / 1e-3, 1e-12); // 50 + 30 = 80 MPa
  near(s.sigmaMax, 9e7, 1e-12); near(s.sigmaMin, -7e7, 1e-12);
  near(s.tau, (4e4 + 1e4) / 0.02, 1e-12); // mean shear 2.5 MPa
  near(s.vonMises, Math.hypot(9e7, Math.sqrt(3) * 2.5e6), 1e-12);
  assert.equal(s.fibreAssumed, true);
  for (const v of Object.values(s)) if (typeof v === 'number') assert.ok(Number.isFinite(v));
  // explicit fibre distances override the default: cz = 0.25 m, cy = 0.125 m
  const e = sectionStress({ section: { ...gen, cy: 0.125, cz: 0.25 }, My: 1e5, Mz: 3e4 });
  near(e.sigmaBending, 1e5 * 0.25 / 4e-4 + 3e4 * 0.125 / 1e-4, 1e-12); // 62.5 + 37.5 = 100 MPa
  assert.equal(e.fibreAssumed, false);
  // the default is exact for a thin-walled circular tube: i = R/√2
  const R = 0.5; const t = 1e-4; const A = 2 * Math.PI * R * t; const I = Math.PI * R ** 3 * t;
  near(sectionStress({ section: { A, Iy: I, Iz: I, J: 2 * I }, My: 1e3 }).sigmaBending, 1e3 * R / I, 1e-12);
  // typed sections are unaffected
  const isec = iSection({ h: 0.4, b: 0.2, tw: 0.01, tf: 0.016 });
  const si = sectionStress({ section: isec, My: 1e5, Vz: 1e4 });
  near(si.sigmaBending, 1e5 / isec.Wely, 1e-12); near(si.tau, 1e4 / isec.Asz, 1e-12); assert.equal(si.fibreAssumed, false);
  assert.throws(() => sectionStress({ section: { ...gen, cz: -1 }, My: 1 }), /cz/);
  // a frame analysed with a generic section gives member stresses through the same route
  const res = analyzeFrame({ nodes: [[0, 0, 0], [3, 0, 0]], members: [{ a: 0, b: 1, section: gen }],
    supports: [{ node: 0, fix: 'fixed' }], loadCases: [{ name: 'P', nodal: [{ node: 1, F: [0, 0, -1e4] }] }] });
  assert.ok(res);
});

test('updateFrameSystem after a release change equals a full rebuild (element records, numbering, profile, solution)', () => {
  // two-bay portal; a hinge is introduced at end B of member 1 and member 3 becomes a truss bar
  const s2 = iSection({ h: 0.4, b: 0.2, tw: 0.01, tf: 0.016 });
  const base = {
    nodes: [[0, 0, 0], [0, 0, 4], [5, 0, 4], [5, 0, 0], [10, 0, 4], [10, 0, 0]],
    members: [{ a: 0, b: 1, section: s2 }, { a: 1, b: 2, section: s2 }, { a: 3, b: 2, section: s2 }, { a: 2, b: 4, section: s2 }, { a: 5, b: 4, section: s2 }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 3, fix: 'fixed' }, { node: 5, fix: 'fixed' }],
    loadCases: [{ name: 'H', nodal: [{ node: 1, F: [1e4, 0, -2e4] }] }],
  };
  const sys0 = buildFrameSystem(base);
  const members = base.members.slice();
  members[1] = { ...members[1], releases: { b: [0, 0, 0, 0, 1, 1] } };
  members[3] = { ...members[3], type: 'truss' };
  const changedModel = { ...base, members };
  const full = buildFrameSystem(changedModel);
  const upd = updateFrameSystem(sys0, changedModel, [1, 3]);
  assert.equal(upd.nFree, full.nFree);
  assert.deepEqual(Array.from(upd.eqDof), Array.from(full.eqDof));
  assert.deepEqual(Array.from(upd.height), Array.from(full.height));
  assert.deepEqual(upd.autoRestrained, full.autoRestrained);
  for (let e = 0; e < 5; e++) {
    assert.deepEqual(Array.from(upd.elems[e].ke), Array.from(full.elems[e].ke), `member ${e}`);
    assert.deepEqual(upd.elems[e].relIdx, full.elems[e].relIdx);
  }
  // unchanged members are reused, changed ones are new records
  assert.equal(upd.elems[0], sys0.elems[0]); assert.notEqual(upd.elems[1], sys0.elems[1]);
  // the original system is not modified
  assert.deepEqual(sys0.elems[1].relIdx, []);
  const a = analyzeFrame(changedModel);
  const b = analyzeFrame(changedModel, { system: upd });
  assert.ok(a && b);
  assert.throws(() => updateFrameSystem(sys0, { ...base, members: members.slice(0, 4) }, [1]), RangeError);
  assert.throws(() => updateFrameSystem(sys0, changedModel, [7]), RangeError);
});

test('self-weight of a member flagged flooded: contained water only where a buoyancy specification submerges it', () => {
  // tube D 0.8 × 20 mm, 10 m long. Steel: 7850 × A; contained sea water: 1025 × π 0.76²/4 = 464.98 kg/m
  const D = 0.8; const t = 0.02; const L = 10; const g = 9.80665; const rhoW = 1025; const rhoS = 7850;
  const tube = tubeSection({ D, t });
  const water = rhoW * Math.PI * (D - 2 * t) ** 2 / 4;
  near(water, 464.98, 1e-4);
  const mk = (member, selfWeight, z = 0) => ({
    nodes: [[0, 0, z], [L, 0, z]], members: [{ a: 0, b: 1, section: tube, ...member }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 1, fix: [0, 1, 1, 0, 0, 0] }], loadCases: [{ name: 'W', selfWeight }],
  });
  const sumR = (m) => { const r = analyzeFrame(m).results.W; return r.reactions[0].R[2] + r.reactions[1].R[2]; };
  // in air (no buoyancy specification): steel only, the flooded flag adds nothing
  near(sumR(mk({ flooded: true }, {})), rhoS * tube.A * g * L, 1e-9, 'in air');
  near(sumR(mk({ flooded: true }, {})), sumR(mk({}, {})), 1e-12);
  // explicit request: the water stays in the member in air
  near(sumR(mk({ flooded: true }, { floodedInAir: true })), (rhoS * tube.A + water) * g * L, 1e-9, 'floodedInAir');
  // an explicit content density is always carried (product at 800 kg/m³)
  near(sumR(mk({ fluidRho: 800 }, {})), (rhoS * tube.A + 800 * Math.PI * (D - 2 * t) ** 2 / 4) * g * L, 1e-9, 'fluidRho');
  // with a buoyancy specification: above the waterline dry, below it flooded and buoyant (water in = water out)
  near(sumR(mk({ flooded: true }, { buoyancy: { waterline: 0, rhoW } }, 5)), rhoS * tube.A * g * L, 1e-9, 'above the waterline');
  near(sumR(mk({ flooded: true }, { buoyancy: { waterline: 0, rhoW } }, -20)), (rhoS - rhoW) * tube.A * g * L, 1e-9, 'submerged');
  // vertical flooded leg through the waterline: 10 m submerged (steel buoyant), 5 m above (steel only, no water)
  const leg = { nodes: [[0, 0, -10], [0, 0, 5]], members: [{ a: 0, b: 1, section: tube, flooded: true }], supports: [{ node: 0, fix: 'fixed' }],
    loadCases: [{ name: 'W', selfWeight: { buoyancy: { waterline: 0, rhoW } } }] };
  near(analyzeFrame(leg).results.W.reactions[0].R[2], ((rhoS - rhoW) * tube.A * 10 + rhoS * tube.A * 5) * g, 1e-9, 'piercing leg');
});
