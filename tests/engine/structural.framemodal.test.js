import test from 'node:test';
import assert from 'node:assert/strict';
import { modalAnalysis, rayleighDamping, beamLocalMass, frameMass } from '../../src/engine/structural/framemodal.js';
import { buildFrameSystem } from '../../src/engine/structural/frame3d.js';
import { tubeSection } from '../../src/engine/structural/sections.js';

const E = 2.1e11; const rho = 7850;
const tube = tubeSection({ D: 0.3, t: 0.01 });
const EI = E * tube.Iy;
const m = rho * tube.A;
const L = 10;
const near = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

// horizontal beam along x, vibration in the x–z plane only
function beam(n, left, right, extra = {}) {
  const nodes = []; const members = []; const supports = [];
  for (let i = 0; i <= n; i++) {
    nodes.push([L * i / n, 0, 0]);
    if (i < n) members.push({ a: i, b: i + 1, section: tube, ...extra });
    const fix = [1, 1, 0, 1, 0, 1]; // axial, lateral and twist suppressed → pure bending modes
    const end = i === 0 ? left : (i === n ? right : null);
    if (end === 'pinned') fix[2] = 1;
    if (end === 'fixed') { fix[2] = 1; fix[4] = 1; }
    supports.push({ node: i, fix });
  }
  return { nodes, members, supports };
}
const fBeam = (lam) => lam ** 2 / (2 * Math.PI * L * L) * Math.sqrt(EI / m);

test('uniform beams: natural frequencies vs closed form (Blevins, Table 8-1)', () => {
  // λ_i: pinned–pinned iπ; cantilever 1.87510407, 4.69409113, 7.85475744; fixed–fixed 4.73004074, 7.85320462
  const pp = modalAnalysis(beam(16, 'pinned', 'pinned'), { nModes: 3, shearDeformation: false });
  [1, 2, 3].forEach((i, k) => near(pp.frequencies[k], fBeam(i * Math.PI), 2e-4, `pinned mode ${i}`));
  const cl = modalAnalysis(beam(16, 'fixed', 'free'), { nModes: 3, shearDeformation: false });
  [1.87510407, 4.69409113, 7.85475744].forEach((l, k) => near(cl.frequencies[k], fBeam(l), 2e-4, `cantilever mode ${k + 1}`));
  const ff = modalAnalysis(beam(16, 'fixed', 'fixed'), { nModes: 2, shearDeformation: false });
  [4.73004074, 7.85320462].forEach((l, k) => near(ff.frequencies[k], fBeam(l), 2e-4, `fixed mode ${k + 1}`));
  assert.ok(pp.converged && cl.converged && ff.converged);
});

test('consistent mass converges from above at fourth order; lumped (HRZ) mass converges', () => {
  const exact = fBeam(Math.PI);
  const e = [2, 4, 8].map((n) => modalAnalysis(beam(n, 'pinned', 'pinned'), { nModes: 1, shearDeformation: false }).frequencies[0] / exact - 1);
  assert.ok(e.every((v) => v > 0), 'upper bound');
  const order = Math.log2(e[1] / e[2]);
  assert.ok(order > 3.5 && order < 4.5, `observed order ${order}`);
  const lumped = modalAnalysis(beam(16, 'pinned', 'pinned'), { nModes: 1, lumped: true, shearDeformation: false }).frequencies[0];
  near(lumped, exact, 5e-3);
});

test('axial bar and tip-mass oscillator', () => {
  // fixed–free bar: f = (1/4L)√(E/ρ)
  const nodes = []; const members = []; const supports = [];
  const n = 12;
  for (let i = 0; i <= n; i++) {
    nodes.push([L * i / n, 0, 0]);
    if (i < n) members.push({ a: i, b: i + 1, section: tube });
    supports.push({ node: i, fix: [i === 0 ? 1 : 0, 1, 1, 1, 1, 1] });
  }
  const bar = modalAnalysis({ nodes, members, supports }, { nModes: 1 });
  near(bar.frequencies[0], Math.sqrt(E / rho) / (4 * L), 1e-3);
  // massless cantilever with a tip mass: f = (1/2π)√(3EI/(M L³))
  const M = 5000;
  const light = { E, rho: 1e-6 };
  const model = {
    nodes: [[0, 0, 0], [0, 0, L]], members: [{ a: 0, b: 1, section: tube, material: light }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 1, fix: [0, 1, 1, 1, 0, 1] }], nodeMasses: [{ node: 1, m: M, inertia: [0, 1e-3, 0] }],
  };
  const r = modalAnalysis(model, { nModes: 1, shearDeformation: false });
  near(r.frequencies[0], Math.sqrt(3 * EI / (M * L ** 3)) / (2 * Math.PI), 1e-5);
  near(r.effectiveMass[0].x, M, 1e-5, 'all the mass participates in the single sway mode');
});

test('effective modal masses sum to (nearly) the total mass; first cantilever mode carries 61.3 %', () => {
  const model = beam(12, 'fixed', 'free');
  const r = modalAnalysis(model, { nModes: 22, shearDeformation: false });
  near(r.massFraction[0].z, 0.6131, 2e-3, 'first-mode mass fraction of a uniform cantilever');
  const sum = r.effectiveMass.reduce((s, e) => s + e.z, 0);
  // the remainder is the share of the first element's mass that is carried directly by the support
  near(r.totalMass.z, m * L, 1e-9);
  assert.ok(sum > 0.96 * r.totalMass.z && sum <= r.totalMass.z * (1 + 1e-9), `sum of effective masses ${sum}`);
  // modes are M-orthonormal: φᵀ M φ = 1
  const sys = r.system;
  const x = new Float64Array(sys.nFree);
  for (let q = 0; q < sys.nFree; q++) x[q] = r.modes[0][sys.eqDof[q]];
  const y = r.mass.mulM(x);
  let s = 0;
  for (let q = 0; q < sys.nFree; q++) s += x[q] * y[q];
  near(s, 1, 1e-8);
});

test('contained water and hydrodynamic added mass lower the frequency by √(m/(m + Δm))', () => {
  const dry = modalAnalysis(beam(8, 'pinned', 'pinned'), { nModes: 1, shearDeformation: false }).frequencies[0];
  const rhoW = 1025; const D = tube.D; const di = D - 2 * tube.t;
  const wetModel = beam(8, 'pinned', 'pinned', { flooded: true });
  wetModel.nodes = wetModel.nodes.map((p) => [p[0], p[1], -30]);
  const wet = modalAnalysis(wetModel, { nModes: 1, shearDeformation: false, waterline: 0, rhoW, Ca: 1.0 });
  const dm = rhoW * Math.PI / 4 * (di * di + D * D);
  near(wet.frequencies[0] / dry, Math.sqrt(m / (m + dm)), 1e-6);
  near(wet.addedMass, rhoW * Math.PI / 4 * D * D * L, 1e-9);
  // above the waterline nothing changes
  const above = modalAnalysis(beam(8, 'pinned', 'pinned', { flooded: true }), { nModes: 1, shearDeformation: false, waterline: -5 });
  near(above.frequencies[0], dry, 1e-9);
});

test('Rayleigh damping: coefficients reproduce the target ratios; mass matrices carry the member mass', () => {
  const r = rayleighDamping({ f1: 0.5, f2: 3, zeta1: 0.02, zeta2: 0.03 });
  near(r.zetaAt(0.5), 0.02, 1e-12);
  near(r.zetaAt(3), 0.03, 1e-12);
  // equal ratios: α = 2ζ ω1 ω2/(ω1 + ω2), β = 2ζ/(ω1 + ω2)
  const q = rayleighDamping({ f1: 1, f2: 4, zeta1: 0.05 });
  near(q.alpha, 2 * 0.05 * (2 * Math.PI) * (8 * Math.PI) / (10 * Math.PI), 1e-12);
  near(q.beta, 2 * 0.05 / (10 * Math.PI), 1e-12);
  for (const lumped of [false, true]) {
    const M = beamLocalMass({ L: 3, mAxial: 40, mTrans: 55, lumped });
    let sx = 0; let sy = 0;
    for (const i of [0, 6]) for (const j of [0, 6]) sx += M[i * 12 + j];
    for (const i of [1, 7]) for (const j of [1, 7]) sy += M[i * 12 + j];
    near(sx, 120, 1e-12); near(sy, 165, 1e-12);
  }
  const sys = buildFrameSystem(beam(2, 'pinned', 'pinned'));
  near(frameMass(sys).structuralMass, m * L, 1e-12);
});
