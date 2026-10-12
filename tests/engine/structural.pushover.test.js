import test from 'node:test';
import assert from 'node:assert/strict';
import {
  pushoverAnalysis, damagedStrength, progressiveCollapse, removeMembers, plasticInteraction, memberPlasticCapacity,
} from '../../src/engine/structural/pushover.js';
import { buildFrameSystem } from '../../src/engine/structural/frame3d.js';
import { tubeSection, rectSection } from '../../src/engine/structural/sections.js';
import { isoTubularStrengths } from '../../src/engine/structural/tubular.js';

const near = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const sec = rectSection({ h: 0.3, b: 0.15 });
const Mp = 250e3;
const big = { Np: 1e12, Mp, Pc: 1e12 }; // bending-only plasticity for the classical plastic-analysis benchmarks

test('fixed-fixed beam, uniform load: first hinges at w = 12 Mp/L², collapse at w = 16 Mp/L²', () => {
  const L = 8;
  const model = {
    nodes: [[0, 0, 0], [L / 2, 0, 0], [L, 0, 0]],
    members: [{ a: 0, b: 1, section: sec, capacity: big }, { a: 1, b: 2, section: sec, capacity: big }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 2, fix: 'fixed' }, { node: 1, fix: [0, 1, 0, 1, 0, 1] }],
    loadCases: [{ name: 'w', member: [{ member: 0, type: 'uniform', q: [0, 0, -1] }, { member: 1, type: 'uniform', q: [0, 0, -1] }] }],
  };
  const r = pushoverAnalysis(model, { pattern: 'w', pDelta: false, control: { node: 1, dof: 'uz' } });
  assert.equal(r.status, 'collapse');
  near(r.firstEventLoadFactor, 12 * Mp / L ** 2, 1e-8, 'first yield');
  near(r.ultimateLoadFactor, 16 * Mp / L ** 2, 1e-8, 'collapse load');
  near(r.redundancyFactor, 4 / 3, 1e-8);
  assert.deepEqual(r.events.map((e) => e.type).sort(), ['hinge-A', 'hinge-A', 'hinge-B', 'hinge-B']);
  // mid-span deflection when the last hinge forms: elastic part w1 L⁴/(384 EI) + simply supported increment 5 Δw L⁴/(384 EI)
  const EI = 2.1e11 * sec.Iy;
  const w1 = 12 * Mp / L ** 2; const dw = 4 * Mp / L ** 2;
  const r2 = pushoverAnalysis(model, { pattern: 'w', pDelta: false, control: { node: 1, dof: 'uz' }, shearDeformation: false });
  near(-r2.curve.displacement.at(-1), w1 * L ** 4 / (384 * EI) + 5 * dw * L ** 4 / (384 * EI), 1e-8);
  near(r2.ultimateBaseShear, 16 * Mp / L ** 2 * L, 1e-8, 'total load at collapse');
});

test('fixed-fixed beam, central point load: all three hinges form together at P = 8 Mp/L', () => {
  const L = 6;
  const model = {
    nodes: [[0, 0, 0], [L / 2, 0, 0], [L, 0, 0]],
    members: [{ a: 0, b: 1, section: sec, capacity: big }, { a: 1, b: 2, section: sec, capacity: big }],
    supports: [{ node: 0, fix: 'fixed' }, { node: 2, fix: 'fixed' }, { node: 1, fix: [0, 1, 0, 1, 0, 1] }],
    loadCases: [{ name: 'P', nodal: [{ node: 1, F: [0, 0, -1, 0, 0, 0] }] }],
  };
  const r = pushoverAnalysis(model, { pattern: 'P', pDelta: false });
  near(r.ultimateLoadFactor, 8 * Mp / L, 1e-8);
  near(r.RSR, 8 * Mp / L, 1e-8);
  assert.equal(r.status, 'collapse');
});

test('fixed-base portal frame vs plastic analysis: sway 4 Mp/h and combined mechanism 6 Mp = H h + V L/2', () => {
  const h = 4; const L = 8;
  const mk = (H, V) => ({
    nodes: [[0, 0, 0], [0, 0, h], [L / 2, 0, h], [L, 0, h], [L, 0, 0]],
    members: [[0, 1], [1, 2], [2, 3], [4, 3]].map(([a, b]) => ({ a, b, section: sec, capacity: big })),
    supports: [{ node: 0, fix: 'fixed' }, { node: 4, fix: 'fixed' }, ...[1, 2, 3].map((n) => ({ node: n, fix: [0, 1, 0, 1, 0, 1] }))],
    loadCases: [{ name: 'L', nodal: [{ node: 1, F: [H, 0, 0, 0, 0, 0] }, { node: 2, F: [0, 0, -V, 0, 0, 0] }] }],
  });
  // pure sway: H h = 4 Mp
  const sway = pushoverAnalysis(mk(1, 0), { pattern: 'L', pDelta: false });
  near(sway.ultimateLoadFactor, 4 * Mp / h, 1e-7);
  assert.equal(sway.events.length, 4);
  // beam mechanism alone: V L/2 = 4 Mp → V = 8 Mp/L
  const beam = pushoverAnalysis(mk(0, 1), { pattern: 'L', pDelta: false });
  near(beam.ultimateLoadFactor, 8 * Mp / L, 1e-7);
  // H = V = P with h = L/2 = 4: sway P = Mp, beam P = Mp, combined 8P = 6Mp → P = 0.75 Mp governs
  const comb = pushoverAnalysis(mk(1, 1), { pattern: 'L', pDelta: false, control: { node: 1, dof: 'ux' } });
  near(comb.ultimateLoadFactor, 0.75 * Mp, 1e-6);
  assert.ok(comb.firstEventLoadFactor < comb.ultimateLoadFactor);
  // the capacity curve is monotonic and ends at the collapse load
  for (let i = 1; i < comb.curve.loadFactor.length; i++) assert.ok(comb.curve.loadFactor[i] >= comb.curve.loadFactor[i - 1] - 1e-9);
  assert.ok(comb.curve.displacement.at(-1) > 0);
});

// three-bar truss: vertical bar (length h) and two bars at ±θ from the vertical, all meeting at the loaded node
const hT = 3; const theta = Math.PI / 4;
const rod = tubeSection({ D: 0.1, t: 0.005 });
function threeBar(capV, capI, Fz) {
  const x = hT * Math.tan(theta);
  return {
    nodes: [[0, 0, 0], [0, 0, hT], [-x, 0, hT], [x, 0, hT]],
    members: [
      { a: 0, b: 1, section: rod, type: 'truss', capacity: capV },
      { a: 0, b: 2, section: rod, type: 'truss', capacity: capI },
      { a: 0, b: 3, section: rod, type: 'truss', capacity: capI },
    ],
    supports: [{ node: 1, fix: 'pinned' }, { node: 2, fix: 'pinned' }, { node: 3, fix: 'pinned' }, { node: 0, fix: [0, 1, 0, 0, 0, 0] }],
    loadCases: [{ name: 'P', nodal: [{ node: 0, F: [0, 0, Fz, 0, 0, 0] }] }],
  };
}
const c = Math.cos(theta);

test('three-bar truss in tension: first yield Np(1 + 2cos³θ), collapse Np(1 + 2cosθ)', () => {
  const Np = 4e5;
  const cap = { Np, Mp: 1, Pc: 1e12 };
  const r = pushoverAnalysis(threeBar(cap, cap, -1), { pattern: 'P', pDelta: false, designFactor: 2e5 });
  near(r.firstEventLoadFactor, Np * (1 + 2 * c ** 3), 1e-9);
  near(r.ultimateLoadFactor, Np * (1 + 2 * c), 1e-9);
  near(r.RSR, Np * (1 + 2 * c) / 2e5, 1e-9);
  assert.deepEqual(r.events.map((e) => [e.member, e.type]), [[0, 'tension-yield'], [1, 'tension-yield'], [2, 'tension-yield']]);
  assert.equal(r.memberStates[0].axial, 'yielded');
  near(r.memberStates[1].N, Np, 1e-9);
});

test('three-bar truss in compression: buckling with post-buckling unloading and load shedding', () => {
  // elastic: Nv = P/(1 + 2cos³θ), Ni = P cos²θ/(1 + 2cos³θ). The vertical bar buckles at P1 = Pcv (1 + 2cos³θ),
  // then carries ρ Pcv; the inclined bars then buckle at P2 = ρ Pcv + 2 cosθ Pci (if P2 > P1).
  const rho = 0.4; const Pcv = 3e5;
  const strong = 5e5; // inclined bars can absorb the shed load
  const r = pushoverAnalysis(threeBar({ Np: 1e12, Mp: 1, Pc: Pcv }, { Np: 1e12, Mp: 1, Pc: strong }, 1), { pattern: 'P', pDelta: false, postBucklingRatio: rho });
  const P1 = Pcv * (1 + 2 * c ** 3);
  const P2 = rho * Pcv + 2 * c * strong;
  near(r.firstEventLoadFactor, P1, 1e-9);
  near(r.ultimateLoadFactor, P2, 1e-9);
  assert.equal(r.events[0].type, 'buckling');
  near(r.memberStates[0].N, -rho * Pcv, 1e-9);
  // weak inclined bars: the load shed by the vertical bar overloads them at constant load → P1 is the ultimate load
  const weak = 1.2e5;
  assert.ok(rho * Pcv + 2 * c * weak < P1);
  const r2 = pushoverAnalysis(threeBar({ Np: 1e12, Mp: 1, Pc: Pcv }, { Np: 1e12, Mp: 1, Pc: weak }, 1), { pattern: 'P', pDelta: false, postBucklingRatio: rho });
  // the inclined bars reach Pci first here? elastic Ni/Pci vs Nv/Pcv: cos²θ/weak = 4.17e-6 > 1/Pcv = 3.33e-6 → inclined buckle first
  const P1i = weak * (1 + 2 * c ** 3) / c ** 2;
  near(r2.firstEventLoadFactor, P1i, 1e-9);
  near(r2.ultimateLoadFactor, P1i, 1e-9, 'progressive failure at constant load');
  assert.equal(r2.status, 'collapse');
});

test('damaged strength ratio and member-removal screening of the three-bar truss', () => {
  const Np = 4e5;
  const cap = { Np, Mp: 1, Pc: 1e12 };
  const model = threeBar(cap, cap, -1);
  // without the vertical bar the two inclined bars carry 2 Np cosθ
  const d = damagedStrength(model, { pattern: 'P', pDelta: false, remove: [0], designFactor: 1e5 });
  near(d.damaged.ultimateLoadFactor, 2 * Np * c, 1e-9);
  near(d.damagedStrengthRatio, 2 * c / (1 + 2 * c), 1e-9);
  near(d.RSRdamaged, 2 * Np * c / 1e5, 1e-9);
  near(d.redundancyFactor, (1 + 2 * c) / (1 + 2 * c ** 3), 1e-9);
  // notional removal of each bar in turn under a 300 kN load with a dynamic amplification factor of 2
  const pc = progressiveCollapse({ ...model, loadCases: [{ name: 'P', nodal: [{ node: 0, F: [0, 0, -3e5, 0, 0, 0] }] }] }, { pattern: 'P', pDelta: false, daf: 2 });
  assert.equal(pc.cases.length, 3);
  near(pc.cases[0].ultimateLoadFactor, 2 * Np * c / 3e5, 1e-9);
  assert.equal(pc.cases[0].survives, false, 'λult = 1.886 < DAF = 2');
  assert.ok(!pc.allSurvive);
  const rm = removeMembers(model, [1]);
  assert.deepEqual(Array.from(rm.indexMap), [0, -1, 1]);
});

test('tube hinge interaction and default capacities; P-Δ lowers the collapse load of a cantilever column', () => {
  const tube = tubeSection({ D: 0.5, t: 0.02 });
  const fy = 355e6;
  const sys = buildFrameSystem({ nodes: [[0, 0, 0], [0, 0, 10]], members: [{ a: 0, b: 1, section: tube }], supports: [{ node: 0, fix: 'fixed' }] });
  const cap = memberPlasticCapacity(sys.elems[0]);
  near(cap.Np, tube.A * fy, 1e-12);
  near(cap.Mpy, tube.Wply * fy, 1e-12);
  near(cap.Pc, isoTubularStrengths({ D: 0.5, t: 0.02, fy, E: 2.1e11, L: 10, K: 1 }).Pc, 1e-12);
  assert.ok(Math.abs(plasticInteraction(cap, 0.5 * cap.Np, cap.Mpy * Math.cos(Math.PI / 4), 0)) < 1e-12);
  assert.ok(Math.abs(plasticInteraction(cap, 0, 0.6 * cap.Mpy, 0.8 * cap.Mpy)) < 1e-12, 'biaxial resultant moment');
  // cantilever column, vertical load V (held) and lateral load H (scaled): first-order collapse H1 = Mp cos(π/2·V/Np)/L
  const L = 10; const V = 0.3 * cap.Np;
  const n = 4;
  const nodes = []; const members = []; const supports = [{ node: 0, fix: 'fixed' }];
  for (let i = 0; i <= n; i++) {
    nodes.push([0, 0, L * i / n]);
    if (i < n) members.push({ a: i, b: i + 1, section: tube, K: 2, bucklingLength: L, capacity: { Pc: 1e12 } });
    if (i > 0) supports.push({ node: i, fix: [0, 1, 0, 1, 0, 1] });
  }
  const model = {
    nodes, members, supports,
    loadCases: [{ name: 'V', nodal: [{ node: n, F: [0, 0, -V, 0, 0, 0] }] }, { name: 'H', nodal: [{ node: n, F: [1, 0, 0, 0, 0, 0] }] }],
  };
  const H1 = cap.Mpy * Math.cos(Math.PI / 2 * 0.3) / L;
  const first = pushoverAnalysis(model, { pattern: 'H', base: 'V', pDelta: false, control: { node: n, dof: 'ux' } });
  near(first.ultimateLoadFactor, H1, 1e-6);
  const second = pushoverAnalysis(model, { pattern: 'H', base: 'V', pDelta: true, control: { node: n, dof: 'ux' } });
  // second order: H L + V Δ = Mp,red at the base → H2 = H1 − V Δ/L with the computed tip deflection
  const delta = second.curve.displacement.at(-1);
  assert.ok(second.ultimateLoadFactor < H1);
  near(second.ultimateLoadFactor, H1 - V * delta / L, 2e-2);
});

test('generic {A, Iy, Iz, J} section: default plastic capacities from the elastic modulus and a shape factor; pushover runs', () => {
  // A = 0.02 m², Iy = 4e-4, Iz = 1e-4 m⁴ → Wel = I/(√2 i) = 2e-3 and 1e-3 m³; shape factor 4/π (thin-walled tube)
  // fy = 355 MPa: Mpy = 4/π × 2e-3 × 355e6 = 904.0 kN·m, Mpz = 452.0 kN·m, Np = 7100 kN
  const gen = { A: 0.02, Iy: 4e-4, Iz: 1e-4, J: 3e-4 };
  const L = 4;
  const model = (section) => ({ nodes: [[0, 0, 0], [L, 0, 0]], members: [{ a: 0, b: 1, section, capacity: { Pc: 1e12 } }],
    supports: [{ node: 0, fix: 'fixed' }], loadCases: [{ name: 'P', nodal: [{ node: 1, F: [0, 0, -1e3] }] }] });
  const cap = memberPlasticCapacity(buildFrameSystem(model(gen)).elems[0]);
  const fy = cap.Np / gen.A;
  near(cap.Mpy, 4 / Math.PI * 2e-3 * fy, 1e-12); near(cap.Mpz, 4 / Math.PI * 1e-3 * fy, 1e-12);
  assert.equal(cap.assumedShape, true);
  // explicit fibre distance and shape factor: Mpy = 1.12 × Iy/0.25 × fy
  const c2 = memberPlasticCapacity(buildFrameSystem(model({ ...gen, cz: 0.25, cy: 0.1, shapeFactor: 1.12 })).elems[0]);
  near(c2.Mpy, 1.12 * 4e-4 / 0.25 * fy, 1e-12); assert.equal(c2.assumedShape, false);
  assert.throws(() => memberPlasticCapacity(buildFrameSystem(model({ ...gen, shapeFactor: 0.5 })).elems[0]), /shapeFactor/);
  // typed sections keep their own plastic modulus
  const tube = tubeSection({ D: 0.5, t: 0.02 });
  const ct = memberPlasticCapacity(buildFrameSystem(model(tube)).elems[0]);
  near(ct.Mpy, tube.Wply * fy, 1e-12); assert.equal(ct.assumedShape, false);
  // tip-loaded cantilever collapses when the root moment reaches Mp: λ × 1 kN × L = Mpy
  const r = pushoverAnalysis(model(gen), { pattern: 'P', pDelta: false });
  near(r.ultimateLoadFactor, cap.Mpy / (1e3 * L), 1e-6);
});

test('a load pattern without any load is rejected up front, also with maxStep', () => {
  const model = {
    nodes: [[0, 0, 0], [4, 0, 0]], members: [{ a: 0, b: 1, section: sec, capacity: big }], supports: [{ node: 0, fix: 'fixed' }],
    loadCases: [{ name: 'zero', nodal: [{ node: 1, F: [0, 0, 0] }] }, { name: 'P', nodal: [{ node: 1, F: [0, 0, -1e3] }] },
      { name: 'empty' }],
  };
  for (const pattern of ['zero', 'empty']) {
    assert.throws(() => pushoverAnalysis(model, { pattern, maxStep: 0.01 }), (e) => e instanceof RangeError && /applies no load/.test(e.message) && e.message.includes(pattern));
    assert.throws(() => pushoverAnalysis(model, { pattern, base: 'P' }), /applies no load/);
  }
  // the same model with a real pattern: root hinge at λ × 1 kN × 4 m = Mp
  near(pushoverAnalysis(model, { pattern: 'P', pDelta: false }).ultimateLoadFactor, Mp / 4e3, 1e-9);
});
