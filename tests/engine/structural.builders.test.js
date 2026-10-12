import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  jacketModel, monopileTowerModel, deckGridModel, wharfBentModel, dolphinModel, portalFrameModel, fromImportedFrame,
} from '../../src/engine/structural/builders.js';
import { analyzeFrame } from '../../src/engine/structural/frame3d.js';
import { modalAnalysis } from '../../src/engine/structural/framemodal.js';
import { tubeSection, iSectionByName } from '../../src/engine/structural/sections.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);
const G0 = 9.80665;

/* Every member end must be reachable from a supported node through members. */
function connected(model) {
  const adj = model.nodes.map(() => []);
  for (const m of model.members) { adj[m.a].push(m.b); adj[m.b].push(m.a); }
  const seen = new Uint8Array(model.nodes.length);
  const stack = model.supports.map((s) => s.node);
  for (const n of stack) seen[n] = 1;
  while (stack.length) { const n = stack.pop(); for (const k of adj[n]) if (!seen[k]) { seen[k] = 1; stack.push(k); } }
  return model.members.every((m) => seen[m.a] && seen[m.b]) && seen.every((v) => v === 1);
}

/* Static solve of the first load case + modal solve; returns the results after the generic checks. */
function solveAndCheck(model, { fMin, fMax, modal = {} } = {}) {
  assert.ok(connected(model), 'model is connected to its supports');
  for (const m of model.members) assert.notEqual(m.a, m.b);
  const r = analyzeFrame(model);
  for (const name of r.order) {
    const c = r.results[name];
    assert.ok(c.equilibrium.relative < 1e-9, `${name}: relative equilibrium residual ${c.equilibrium.relative}`);
    for (let k = 0; k < 3; k++) {
      const scale = Math.max(1, Math.abs(c.appliedSum[2]));
      assert.ok(Math.abs(c.reactionSum[k] + c.appliedSum[k]) < 1e-8 * scale, `${name}: ΣR + ΣF = 0 in direction ${k}`);
    }
    for (const v of c.displacements) assert.ok(Number.isFinite(v));
  }
  const mo = modalAnalysis(model, { nModes: 3, ...modal });
  assert.ok(mo.converged);
  if (fMin !== undefined) assert.ok(mo.frequencies[0] > fMin && mo.frequencies[0] < fMax, `first natural frequency ${mo.frequencies[0]} Hz outside ${fMin}–${fMax}`);
  return { static: r, modal: mo, first: r.results[r.order[0]] };
}

const leg = { D: 1.2, t: 0.04 }; const brace = { D: 0.6, t: 0.02 };
const jacketBase = { waterDepth: 50, topElevation: 15, topWidthX: 12, bays: 4, legSection: leg, braceSection: brace };

test('jacket: 4-leg and 3-leg, X / K / diagonal bracing — connected, stable, equilibrium, first frequency 1–4 Hz', () => {
  for (const legs of [4, 3]) {
    for (const bracing of ['X', 'K', 'diagonal']) {
      const m = jacketModel({ ...jacketBase, legs, bracing, deckMass: 500e3 });
      const { modal, first } = solveAndCheck(m, { fMin: 1, fMax: 4 });
      // gravity + buoyancy only: no horizontal resultant; the two sway modes of a symmetric X/K jacket coincide
      assert.ok(Math.abs(first.reactionSum[0]) < 1e-6 * first.reactionSum[2]);
      if (bracing !== 'diagonal') close(modal.frequencies[1], modal.frequencies[0], 1e-3);
      assert.equal(m.meta.legNodes.length, 5);
      assert.equal(m.meta.groups.leg.length, legs * 4);
      assert.equal(m.meta.groups.brace.length, bracing === 'diagonal' ? legs * 4 : legs * 8);
      assert.equal(m.supports.length, legs);
    }
  }
});

test('jacket: geometry follows the batter; deck mass enters self-weight; lateral load is carried to the supports', () => {
  const m = jacketModel({ ...jacketBase, deckMass: 500e3 });
  // H = 65 m, batter 1:8 → each leg moves 65/8 m outward in x and y: base width 12 + 2·65/8 = 28.25 m
  close(m.meta.baseWidth[0], 28.25, 1e-12);
  const p0 = m.nodes[m.meta.mudlineNodes[0]]; const p1 = m.nodes[m.meta.topNodes[0]];
  close(p0[0], -14.125, 1e-12); close(p0[2], -50, 1e-12); close(p1[0], -6, 1e-12); close(p1[2], 15, 1e-12);
  assert.deepEqual(m.meta.levels, [-50, -33.75, -17.5, -1.25, 15]);
  const noDeck = jacketModel({ ...jacketBase });
  const w1 = analyzeFrame(m).results.selfWeight.reactionSum[2];
  const w0 = analyzeFrame(noDeck).results.selfWeight.reactionSum[2];
  close(w1 - w0, 500e3 * G0, 1e-9, 'deck weight');
  // 1 MN at the top in x: ΣRx = −1 MN and the overturning moment about the mudline is 1 MN·65 m
  m.loadCases.push({ name: 'push', nodal: m.meta.topNodes.map((n) => ({ node: n, F: [250e3, 0, 0, 0, 0, 0] })) });
  const r = analyzeFrame(m, { cases: ['push'] }).results.push;
  close(r.reactionSum[0], -1e6, 1e-9);
  let My = 0;
  for (const rc of r.reactions) { const p = m.nodes[rc.node]; My += (p[2] + 50) * rc.R[0] - p[0] * rc.R[2] + rc.R[4]; }
  close(My, -65e6, 1e-8, 'moment equilibrium about the mudline y axis');
  assert.ok(r.equilibrium.relative < 1e-9);
  // vertical legs when batter = Infinity; base width overrides batter
  close(jacketModel({ ...jacketBase, batter: Infinity }).meta.baseWidth[0], 12, 1e-12);
  close(jacketModel({ ...jacketBase, baseWidthX: 30 }).meta.baseWidth[0], 30, 1e-12);
});

test('jacket: pile idealisations (pinned, springs, fixity stubs) all give stable models', () => {
  const pinned = solveAndCheck(jacketModel({ ...jacketBase, piles: { mode: 'pinned' } }), { fMin: 1, fMax: 8 });
  const fixed = solveAndCheck(jacketModel({ ...jacketBase }), { fMin: 1, fMax: 8 });
  const spr = solveAndCheck(jacketModel({ ...jacketBase, piles: { mode: 'springs', spring: [1e9, 1e9, 2e9, 1e10, 1e10, 1e10] } }), { fMin: 1, fMax: 8 });
  const fix = jacketModel({ ...jacketBase, piles: { mode: 'fixity', section: { D: 1.0, t: 0.04 } } });
  const stub = solveAndCheck(fix, { fMin: 1, fMax: 8 });
  // softer foundations lower the first frequency
  assert.ok(spr.modal.frequencies[0] < fixed.modal.frequencies[0]);
  assert.ok(stub.modal.frequencies[0] < fixed.modal.frequencies[0]);
  assert.ok(pinned.modal.frequencies[0] <= fixed.modal.frequencies[0] * (1 + 1e-9));
  // fixity stubs: 6 D below the mudline, fixed at the tip
  assert.equal(fix.meta.groups.pile.length, 4);
  for (const s of fix.supports) close(fix.nodes[s.node][2], -56, 1e-12);
  assert.throws(() => jacketModel({ ...jacketBase, legs: 5 }), RangeError);
  assert.throws(() => jacketModel({ ...jacketBase, piles: { mode: 'springs' } }), RangeError);
  assert.throws(() => jacketModel({ ...jacketBase, bracing: 'Z' }), RangeError);
});

test('monopile + tower: a prismatic stick clamped at the mudline reproduces the cantilever frequency 1.8751² √(EI/mL⁴)/2π', () => {
  // same tube for pile and tower, no RNA, dry (no waterline passed to the modal solver → no added mass)
  const D = 5; const t = 0.05; const L = 100;
  const m = monopileTowerModel({ waterDepth: 20, pileD: D, pileT: t, interfaceZ: 10, hubHeight: 80, towerBaseD: D, towerTopD: D, towerBaseT: t, foundation: { mode: 'fixed' }, nPile: 6, nTower: 14 });
  for (const mem of m.members) mem.flooded = false;
  const s = tubeSection({ D, t });
  const f = 1.87510407 ** 2 / (2 * Math.PI) * Math.sqrt(2.1e11 * s.Iy / (7850 * s.A * L ** 4));
  const mo = modalAnalysis(m, { nModes: 2, shearDeformation: false });
  close(mo.frequencies[0], f, 2e-3, 'Euler–Bernoulli cantilever (rotary inertia neglected in the closed form)');
  close(mo.frequencies[1], mo.frequencies[0], 1e-6, 'axisymmetric: two equal bending modes');
  close(mo.structuralMass, 7850 * s.A * L, 1e-9);
  // self-weight reaction of the dry stick in air: move the waterline out of the way
  m.loadCases = [{ name: 'dry', selfWeight: {} }];
  close(analyzeFrame(m).results.dry.reactionSum[2], 7850 * s.A * L * G0, 1e-9);
});

test('monopile + tower: tapered tower with RNA — stable, equilibrium, first frequency 0.15–0.4 Hz; foundation flexibility lowers it', () => {
  const base = { waterDepth: 30, pileD: 6, pileT: 0.06, interfaceZ: 10, hubHeight: 90, towerBaseD: 6, towerTopD: 3.87, towerBaseT: 0.027, towerTopT: 0.019, rnaMass: 350e3 };
  const fixity = monopileTowerModel(base);
  const a = solveAndCheck(fixity, { fMin: 0.15, fMax: 0.4, modal: { waterline: 0 } });
  const b = solveAndCheck(monopileTowerModel({ ...base, foundation: { mode: 'fixed' } }), { fMin: 0.15, fMax: 0.4, modal: { waterline: 0 } });
  const c = solveAndCheck(monopileTowerModel({ ...base, foundation: { mode: 'spring', spring: [8e8, 8e8, 5e9, 1e11, 1e11, 1e10] } }), { fMin: 0.1, fMax: 0.4, modal: { waterline: 0 } });
  assert.ok(a.modal.frequencies[0] < b.modal.frequencies[0]);
  assert.ok(c.modal.frequencies[0] < b.modal.frequencies[0]);
  // default fixity depth 4 D below the mudline; tower segments use the mean diameter of each segment
  close(fixity.nodes[fixity.supports[0].node][2], -30 - 24, 1e-12);
  close(fixity.nodes[fixity.meta.topNode][2], 90, 1e-12);
  const first = fixity.members[fixity.meta.towerMembers[0]].section; const last = fixity.members[fixity.meta.towerMembers.at(-1)].section;
  close(first.D, 6 - (6 - 3.87) / 24, 1e-12); close(last.D, 3.87 + (6 - 3.87) / 24, 1e-12);
  close(last.t, 0.019 + 0.008 / 24, 1e-12);
  assert.equal(fixity.meta.pileMembers.length, 9); assert.equal(fixity.meta.towerMembers.length, 12);
  assert.throws(() => monopileTowerModel({ ...base, interfaceZ: 100 }), RangeError);
  assert.throws(() => monopileTowerModel({ ...base, foundation: { mode: 'spring' } }), RangeError);
});

test('deck grid: tributary-area nodal loads sum to p·A, reactions balance, column options', () => {
  const I = iSectionByName('HEB 300');
  const col = { type: 'tube', D: 0.8, t: 0.03 };
  const counts = { corners: 4, perimeter: 10, all: 12 };
  let fPrev = 0;
  for (const columns of ['corners', 'perimeter', 'all']) {
    const m = deckGridModel({ lengthX: 30, widthY: 20, columnHeight: 6, girderSection: I, columnSection: col, columns, areaLoad: 10e3 });
    const r = solveAndCheck(m, { fMin: 0.5, fMax: 30 });
    close(r.static.results.deckLoad.appliedSum[2], -10e3 * 600, 1e-12, 'Σ nodal loads = p·A');
    close(r.static.results.deckLoad.reactionSum[2], 6e6, 1e-9);
    assert.equal(m.supports.length, counts[columns]);
    assert.equal(m.meta.columnNodes.length, counts[columns]);
    assert.ok(r.modal.frequencies[0] > fPrev, 'more columns → stiffer deck');
    fPrev = r.modal.frequencies[0];
  }
  // pinned column bases are softer in sway than fixed ones
  const fixed = modalAnalysis(deckGridModel({ lengthX: 30, widthY: 20, columnHeight: 6, girderSection: I, columnSection: col }), { nModes: 1 });
  const pinned = modalAnalysis(deckGridModel({ lengthX: 30, widthY: 20, columnHeight: 6, girderSection: I, columnSection: col, columnBase: 'pinned' }), { nModes: 1 });
  assert.ok(pinned.frequencies[0] < fixed.frequencies[0]);
  // symmetric load on a symmetric deck: equal corner reactions p·A/4 + self-weight/4
  const m = deckGridModel({ lengthX: 30, widthY: 20, columnHeight: 6, girderSection: I, columnSection: col, areaLoad: 10e3 });
  const rc = analyzeFrame(m).results.deckLoad.reactions;
  for (const r of rc) close(r.R[2], 1.5e6, 1e-8);
  assert.throws(() => deckGridModel({ lengthX: 30, widthY: 20, columnHeight: 6, girderSection: I, columnSection: col, columns: 'none' }), RangeError);
});

test('portal frame: sway stiffness of a fixed-base portal, δ = H h³/(12 E Ic)·(3k + 2)/(6k + 1), k = Ib h/(Ic L)', () => {
  // slope-deflection result for a single-bay portal with inextensible members; axial and shear flexibility suppressed
  const E = 2.1e11; const Ic = 2e-4; const Ib = 5e-4; const h = 5; const L = 8; const H = 100e3;
  const sec = (I) => ({ A: 1e3, Iy: I, Iz: I, J: I });
  const m = portalFrameModel({ span: L, height: h, columnSection: sec(Ic), rafterSection: sec(Ib), nRafter: 1 });
  m.loadCases = [{ name: 'sway', nodal: [{ node: m.meta.eavesNodes[0], F: [H, 0, 0, 0, 0, 0] }] }];
  const r = analyzeFrame(m, { shearDeformation: false }).results.sway;
  const k = Ib * h / (Ic * L);
  const delta = H * h ** 3 / (12 * E * Ic) * (3 * k + 2) / (6 * k + 1);
  close(r.displacements[6 * m.meta.eavesNodes[0]], delta, 1e-5);
  assert.ok(r.equilibrium.relative < 1e-9);
  // limits of the formula: rigid beam → H h³/(24 E Ic); no beam → H h³/(6 E Ic)
  close((3 * 1e9 + 2) / (6 * 1e9 + 1), 0.5, 1e-8);
});

test('portal frame: pitched multi-bay frame — connected, stable, apex geometry, plane-frame restraints', () => {
  const I = iSectionByName('HEB 300');
  const m = portalFrameModel({ span: 20, height: 6, rise: 2, bays: 2, columnSection: I, rafterSection: iSectionByName('IPE 400') });
  const r = solveAndCheck(m, { fMin: 1, fMax: 15 });
  assert.equal(m.meta.apexNodes.length, 2);
  assert.deepEqual(m.nodes[m.meta.apexNodes[0]], [10, 0, 8]);
  assert.deepEqual(m.nodes[m.meta.apexNodes[1]], [30, 0, 8]);
  assert.equal(m.meta.columnMembers.length, 3); assert.equal(m.meta.rafterMembers.length, 8);
  // no out-of-plane displacement anywhere
  const u = r.first.displacements;
  for (let n = 0; n < m.nodes.length; n++) assert.equal(u[6 * n + 1], 0);
  // pinned bases carry no moment
  const p = portalFrameModel({ span: 20, height: 6, columnSection: I, rafterSection: I, base: 'pinned' });
  const rp = solveAndCheck(p, { fMin: 0.5, fMax: 15 });
  for (const rc of rp.first.reactions.filter((x) => p.meta.baseNodes.includes(x.node))) assert.ok(Math.abs(rc.R[4]) < 1e-6);
});

test('wharf bent: vertical and raker piles — toes at seabed − fixity depth, raker slope, stable; rakers stiffen the bent', () => {
  const deck = { type: 'rect', h: 1, b: 1.2 };
  const args = { pileX: [0, 6, 12], deckElevation: 4, seabedZ: [-4, -8, -12], deckSection: deck, pileSection: { D: 0.9, t: 0.02 }, overhang: 1.5 };
  const plain = wharfBentModel(args);
  const rak = wharfBentModel({ ...args, rakers: [{ pile: 1, batter: 4, direction: -1 }] });
  const a = solveAndCheck(plain, { fMin: 0.3, fMax: 20 });
  const b = solveAndCheck(rak, { fMin: 0.3, fMax: 20 });
  assert.ok(b.modal.frequencies[0] > a.modal.frequencies[0]);
  // default fixity 5 D = 4.5 m below the local seabed
  assert.deepEqual(plain.meta.toeNodes.map((n) => plain.nodes[n][2]), [-8.5, -12.5, -16.5]);
  // raker from the middle pile head: toe at z = −12.5, x = 6 − (4 + 12.5)/4
  const toe = rak.nodes[rak.meta.toeNodes[3]];
  close(toe[0], 6 - 16.5 / 4, 1e-12); close(toe[2], -12.5, 1e-12);
  assert.equal(plain.meta.deckNodes.length, 5);
  // self-weight of the flooded (open-ended) piles is a net downward load
  assert.ok(a.first.appliedSum[2] < 0);
  // lateral load: horizontal reactions balance and the deck beam stays in the plane
  rak.loadCases.push({ name: 'berthing', nodal: [{ node: rak.meta.headNodes[0], F: [500e3, 0, 0, 0, 0, 0] }] });
  const r = analyzeFrame(rak, { cases: ['berthing'] }).results.berthing;
  let rx = 0;
  for (const rc of r.reactions) rx += rc.R[0];
  close(rx, -500e3, 1e-9);
  // pinned heads: stable with two or more piles, rejected for a single pile
  solveAndCheck(wharfBentModel({ ...args, pinnedHeads: true }), { fMin: 0.2, fMax: 20 });
  assert.throws(() => wharfBentModel({ ...args, pileX: [0], seabedZ: -5, pinnedHeads: true }), RangeError);
  assert.throws(() => wharfBentModel({ ...args, seabedZ: 5 }), RangeError);
  assert.throws(() => wharfBentModel({ ...args, pileX: [0, 0] }), RangeError);
});

test('dolphin: lateral stiffness of n vertical piles under a rigid cap, k = n·12EI/L³ (rotation of the cap prevented by pile axial stiffness to within a few %)', () => {
  // piles on a large radius: cap rotation ≈ 0 so each pile is fixed–guided, k = 12 E I / L³
  const pile = { D: 1.2, t: 0.03 };
  const m = dolphinModel({ nPiles: 4, radius: 6, capElevation: 5, seabedZ: -12, pileSection: pile, fixityDepth: 6, capMass: 200e3 });
  solveAndCheck(m, { fMin: 0.5, fMax: 10 });
  const s = tubeSection(pile);
  const L = 5 + 12 + 6;
  m.loadCases = [
    { name: 'x', nodal: [{ node: m.meta.capNode, F: [1e6, 0, 0, 0, 0, 0] }] },
    { name: 'y', nodal: [{ node: m.meta.capNode, F: [0, 1e6, 0, 0, 0, 0] }] },
  ];
  const r = analyzeFrame(m, { shearDeformation: false });
  const ux = r.results.x.displacements[6 * m.meta.capNode];
  const uy = r.results.y.displacements[6 * m.meta.capNode + 1];
  const kGuided = 4 * 12 * 2.1e11 * s.Iy / L ** 3;
  close(ux, uy, 1e-9, 'same stiffness in x and y');
  assert.ok(1e6 / ux < kGuided && 1e6 / ux > 0.93 * kGuided, `k = ${1e6 / ux}, fixed–guided bound ${kGuided}`);
  // lower bound: free cantilevers 3EI/L³ each
  assert.ok(1e6 / ux > 4 * 3 * 2.1e11 * s.Iy / L ** 3);
  // raked piles are stiffer laterally
  const raked = dolphinModel({ nPiles: 4, radius: 6, capElevation: 5, seabedZ: -12, pileSection: pile, fixityDepth: 6, batter: 5 });
  raked.loadCases = [m.loadCases[0]];
  raked.loadCases[0] = { name: 'x', nodal: [{ node: raked.meta.capNode, F: [1e6, 0, 0, 0, 0, 0] }] };
  assert.ok(analyzeFrame(raked).results.x.displacements[6 * raked.meta.capNode] < ux);
  close(Math.hypot(raked.nodes[raked.meta.toeNodes[0]][0], raked.nodes[raked.meta.toeNodes[0]][1]), 6 + 23 / 5, 1e-12);
  assert.throws(() => dolphinModel({ nPiles: 2, radius: 6, capElevation: 5, seabedZ: -12, pileSection: pile }), RangeError);
});

test('imported frame: conversion, default sections, automatic supports at the lowest level, degenerate members dropped', () => {
  const frame = {
    kind: 'frame', name: 'bent',
    nodes: new Float64Array([0, 0, 0, 0, 0, 5, 5, 0, 5, 5, 0, 0.001]),
    members: [{ a: 0, b: 1 }, { a: 1, b: 2, group: 'beam', section: { A: 0.01, Iy: 1e-4, Iz: 1e-4, J: 2e-4 } }, { a: 2, b: 3, section: { type: 'tube', D: 0.4, t: 0.01 } }, { a: 1, b: 1 }],
  };
  const m = fromImportedFrame(frame);
  assert.equal(m.members.length, 3);
  assert.deepEqual(m.members[0].section, { type: 'tube', D: 0.5, t: 0.02 });
  assert.equal(m.members[1].group, 'beam');
  assert.deepEqual(m.supports.map((s) => s.node), [0, 3]); // within 1e-3·height of the lowest node
  assert.equal(m.meta.unusedNodes, 0);
  solveAndCheck(m, { fMin: 1, fMax: 100 });
  const explicit = fromImportedFrame(frame, { supports: [{ node: 0, fix: 'fixed' }], supportType: 'pinned' });
  assert.equal(explicit.supports.length, 1);
  assert.equal(fromImportedFrame(frame, { supportType: 'pinned' }).supports[0].fix, 'pinned');
  assert.throws(() => fromImportedFrame({ kind: 'surface' }), RangeError);
  assert.throws(() => fromImportedFrame({ ...frame, members: [{ a: 0, b: 9 }] }), RangeError);
});
