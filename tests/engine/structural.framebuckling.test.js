import test from 'node:test';
import assert from 'node:assert/strict';
import {
  bucklingAnalysis, pDeltaAnalysis, beamGeometricStiffness, effectiveLength,
} from '../../src/engine/structural/framebuckling.js';
import { tubeSection } from '../../src/engine/structural/sections.js';

const E = 2.1e11;
const tube = tubeSection({ D: 0.3, t: 0.01 });
const EI = E * tube.Iy;
const L = 6;
const near = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

// vertical column of n elements, in-plane (x–z) behaviour only; unit compressive load at the top
function column(n, base, top) {
  const nodes = []; const members = []; const supports = [];
  for (let i = 0; i <= n; i++) {
    nodes.push([0, 0, L * i / n]);
    if (i < n) members.push({ a: i, b: i + 1, section: tube });
    // out-of-plane dofs and twist suppressed everywhere
    const fix = [0, 1, 0, 1, 0, 1];
    if (i === 0) { fix[0] = 1; fix[2] = 1; if (base === 'fixed') fix[4] = 1; }
    if (i === n) { if (top !== 'free') fix[0] = 1; if (top === 'fixed') fix[4] = 1; }
    supports.push({ node: i, fix });
  }
  return { nodes, members, supports, loadCases: [{ name: 'P', nodal: [{ node: n, F: [0, 0, -1, 0, 0, 0] }] }] };
}

const cases = [
  { name: 'pinned–pinned', base: 'pinned', top: 'pinned', exact: Math.PI ** 2 * EI / L ** 2, K: 1 },
  { name: 'fixed–free', base: 'fixed', top: 'free', exact: Math.PI ** 2 * EI / (4 * L ** 2), K: 2 },
  // fixed–pinned: tan(kL) = kL → kL = 4.493409458 → Pcr = 20.19073 EI/L²
  { name: 'fixed–pinned', base: 'fixed', top: 'pinned', exact: 4.493409458 ** 2 * EI / L ** 2, K: Math.PI / 4.493409458 },
  { name: 'fixed–fixed', base: 'fixed', top: 'fixed', exact: 4 * Math.PI ** 2 * EI / L ** 2, K: 0.5 },
];

for (const c of cases) {
  test(`Euler buckling, ${c.name}: converges to the closed form at fourth order`, () => {
    const err = [];
    let last;
    for (const n of [2, 4, 8]) {
      last = bucklingAnalysis(column(n, c.base, c.top), { loadCase: 'P', nModes: 1, shearDeformation: false });
      assert.ok(last.converged);
      err.push(last.factors[0] / c.exact - 1);
    }
    assert.ok(err.every((e) => e > -1e-9), 'consistent geometric stiffness gives an upper bound');
    assert.ok(Math.abs(err[2]) < 1e-3, `error with 8 elements ${err[2]}`);
    const order = Math.log2(err[1] / err[2]);
    assert.ok(order > 3.5 && order < 4.6, `observed order ${order}`);
    // effective length factor recovered for the whole column (sum over the n elements)
    const el = last.effectiveLengths[0];
    near(el.Ley / L, c.K, 5e-4, 'K factor');
  });
}

test('second and third pinned–pinned modes: 4 and 9 times the first', () => {
  const r = bucklingAnalysis(column(12, 'pinned', 'pinned'), { loadCase: 'P', nModes: 3, shearDeformation: false });
  near(r.factors[1] / r.factors[0], 4, 2e-3);
  near(r.factors[2] / r.factors[0], 9, 1e-2);
});

test('P-Δ / P-δ: cantilever beam-column, δ = (H/(k³EI))(tan kL − kL), M = (H/k) tan kL', () => {
  const P = 0.4 * Math.PI ** 2 * EI / (4 * L * L); // 40 % of the Euler load
  const H = 1e4;
  const m = column(8, 'fixed', 'free');
  m.loadCases = [{ name: 'PH', nodal: [{ node: 8, F: [H, 0, -P, 0, 0, 0] }] }];
  const r = pDeltaAnalysis(m, { loadCase: 'PH', shearDeformation: false });
  const k = Math.sqrt(P / EI);
  const delta = H / (k ** 3 * EI) * (Math.tan(k * L) - k * L);
  const Mbase = H / k * Math.tan(k * L);
  near(r.result.displacements[6 * 8], delta, 2e-4, 'tip deflection');
  near(Math.abs(r.result.memberForces[0].endA.My), Mbase, 5e-4, 'base moment');
  near(Math.abs(r.result.reactions[0].R[4]), Mbase, 5e-4, 'base reaction moment');
  // first-order values for reference: amplification ≈ 1/(1 − P/Pcr) = 1.667 (approximately)
  near(r.firstOrder.displacements[6 * 8], H * L ** 3 / (3 * EI), 1e-9);
  assert.ok(r.amplification > 1.6 && r.amplification < 1.7);
  assert.ok(r.converged && r.iterations <= 3);
});

test('sway frame: buckling factor of a pinned-base portal with a rigid beam, Pcr = π²EI/(2h)² per column', () => {
  const h = 5; const w = 6;
  const stiff = { ...tube, Iy: tube.Iy * 1e5, Iz: tube.Iz * 1e5, A: tube.A * 100 };
  const n = 6;
  const nodes = []; const members = []; const supports = [];
  for (let c = 0; c < 2; c++) {
    for (let i = 0; i <= n; i++) {
      nodes.push([c * w, 0, h * i / n]);
      const id = c * (n + 1) + i;
      // axially rigid columns, so that the beam cannot rotate by column shortening (hand-solution assumption)
      if (i < n) members.push({ a: id, b: id + 1, section: { ...tube, A: tube.A * 1e4 } });
      supports.push({ node: id, fix: i === 0 ? [1, 1, 1, 1, 0, 1] : [0, 1, 0, 1, 0, 1] });
    }
  }
  members.push({ a: n, b: 2 * n + 1, section: stiff });
  const model = {
    nodes, members, supports,
    loadCases: [{ name: 'P', nodal: [{ node: n, F: [0, 0, -1, 0, 0, 0] }, { node: 2 * n + 1, F: [0, 0, -1, 0, 0, 0] }] }],
  };
  const r = bucklingAnalysis(model, { loadCase: 'P', nModes: 1, shearDeformation: false });
  near(r.factors[0], Math.PI ** 2 * EI / (2 * h) ** 2, 2e-4);
});

test('geometric stiffness matrix and effective-length helper', () => {
  const k = beamGeometricStiffness({ N: -1e5, L: 2 });
  near(k[1 * 12 + 1], -1e5 * 36 / 60, 1e-12);
  near(k[5 * 12 + 5], -1e5 * 16 / 60, 1e-12);
  for (let i = 0; i < 12; i++) for (let j = 0; j < 12; j++) assert.equal(k[i * 12 + j], k[j * 12 + i]);
  const e = effectiveLength({ E, I: tube.Iy, Pcr: Math.PI ** 2 * EI / (0.7 * L) ** 2, L });
  near(e.K, 0.7, 1e-12);
  assert.throws(() => bucklingAnalysis(column(2, 'fixed', 'free'), {}), RangeError);
});

test('P-Δ swayAmplification isolates the lateral ratio when gravity shortening dominates the largest displacement', () => {
  // cantilever column at 40 % of its Euler load with a small lateral load: the largest displacement is the axial
  // shortening P L/(E A), which second-order effects do not change, so amplification ≈ 1; the lateral tip
  // displacement is amplified by 3 (tan kL − kL)/(kL)³ = 1.643 (≈ 1/(1 − 0.4) = 1.667)
  const P = 0.4 * Math.PI ** 2 * EI / (4 * L * L);
  const H = 1;
  const m = column(8, 'fixed', 'free');
  m.loadCases = [{ name: 'PH', nodal: [{ node: 8, F: [H, 0, -P, 0, 0, 0] }] }];
  const r = pDeltaAnalysis(m, { loadCase: 'PH', shearDeformation: false });
  const kL = Math.sqrt(P / EI) * L;
  const exact = 3 * (Math.tan(kL) - kL) / kL ** 3;
  near(r.swayAmplification, exact, 5e-4);
  assert.ok(exact > 1.6 && exact < 1.7);
  near(r.amplification, 1, 1e-3, 'largest displacement is the vertical shortening');
  near(r.firstOrder.maxDisplacement, P * L / (E * tube.A), 1e-6);
  assert.equal(r.sway.node, 8);
  near(r.sway.firstOrder, H * L ** 3 / (3 * EI), 1e-9);
  near(r.sway.secondOrder / r.sway.firstOrder, r.swayAmplification, 1e-12);
  // no lateral load at all: nothing to amplify
  const g = column(8, 'fixed', 'free');
  g.loadCases = [{ name: 'PH', nodal: [{ node: 8, F: [0, 0, -P, 0, 0, 0] }] }];
  assert.equal(pDeltaAnalysis(g, { loadCase: 'PH', shearDeformation: false }).swayAmplification, 1);
  assert.throws(() => pDeltaAnalysis(m, { loadCase: 'PH', verticalAxis: 'up' }), /verticalAxis/);
});

test('buckling with many slender members in tension: the positive factor is found although reversed-load modes lie closer to zero', () => {
  // a stocky fixed–free column (D 1.0 × 40 mm, 20 m, 1 MN compression) beside 12 slender ties (D 0.1 × 5 mm, 200 kN
  // tension), every member in 4 elements. Reversed, each tie would buckle at |λ| = π²EI/(4L²)/N = 0.0054, some 3000
  // times nearer to zero than the column factor λ = π²EI/(2L)² / 1 MN; there are 48 such modes (two planes, two per
  // tie resolved) against a starting subspace of 10 vectors.
  const col = tubeSection({ D: 1.0, t: 0.04 }); const tie = tubeSection({ D: 0.1, t: 0.005 });
  const nodes = []; const members = []; const supports = []; const nodal = [];
  const add = (x, section, F) => {
    const b = nodes.length;
    for (let i = 0; i <= 4; i++) { nodes.push([x, 0, 5 * i]); if (i < 4) members.push({ a: b + i, b: b + i + 1, section }); }
    supports.push({ node: b, fix: 'fixed' }); nodal.push({ node: b + 4, F: [0, 0, F, 0, 0, 0] });
  };
  add(0, col, -1e6);
  for (let k = 0; k < 12; k++) add(5 * (k + 1), tie, 2e5);
  const r = bucklingAnalysis({ nodes, members, supports, loadCases: [{ name: 'ref', nodal }] }, { loadCase: 'ref', nModes: 2, shearDeformation: false });
  const exact = Math.PI ** 2 * E * col.Iy / 40 ** 2 / 1e6;
  assert.ok(r.converged && r.iterations <= 15, `iterations ${r.iterations}`);
  near(r.factors[0], exact, 1e-3); near(r.factors[1], exact, 1e-3); // round column: two planes
  assert.ok(r.shift > 0 && r.shift < r.factors[0], `shift ${r.shift}`);
  // effective length of the column elements: K L = 2 × 20 m over four 5 m elements → K = 8
  const el = r.effectiveLengths.find((e) => e.member === 0);
  near(el.Ley, 40, 1e-3); near(el.Ky, 8, 1e-3);
  // all members in tension: nothing buckles
  const pull = { nodes, members, supports, loadCases: [{ name: 'ref', nodal: nodal.map((l) => ({ node: l.node, F: [0, 0, Math.abs(l.F[2]), 0, 0, 0] })) }] };
  assert.throws(() => bucklingAnalysis(pull, { loadCase: 'ref' }), /No positive buckling factor/);
});
