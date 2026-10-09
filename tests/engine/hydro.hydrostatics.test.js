import test from 'node:test';
import assert from 'node:assert/strict';
import * as hs from '../../src/engine/hydro/hydrostatics.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b}`);
const barge = () => hs.boxMesh({ xMin: -50, xMax: 50, yMin: -10, yMax: 10, zMin: 0, zMax: 10 });

function sphereMesh(R, nLat, nLon, zc) {
  const pos = []; const idx = [];
  for (let i = 0; i <= nLat; i++) {
    const th = (Math.PI * i) / nLat;
    for (let j = 0; j < nLon; j++) {
      const ph = (2 * Math.PI * j) / nLon;
      pos.push(R * Math.sin(th) * Math.cos(ph), R * Math.sin(th) * Math.sin(ph), zc + R * Math.cos(th));
    }
  }
  for (let i = 0; i < nLat; i++) for (let j = 0; j < nLon; j++) {
    const a = i * nLon + j; const b = i * nLon + ((j + 1) % nLon); const c = a + nLon; const d = b + nLon;
    if (i > 0) idx.push(a, c, b);
    if (i < nLat - 1) idx.push(b, c, d);
  }
  return { positions: Float32Array.from(pos), indices: Uint32Array.from(idx) };
}

test('Box barge 100 × 20 at T = 4 m: V = LBT, KB = T/2, BM = B²/12T, BMl = L²/12T (closed form)', () => {
  const rho = 1025; const g = 9.80665;
  const h = hs.meshHydrostatics({ mesh: barge(), draft: 4, cog: [0, 0, 5], rho, g });
  close(h.volume, 8000, 1e-12); close(h.KB, 2, 1e-12);
  close(h.BMt, 400 / 48, 1e-12); close(h.BMl, 10000 / 48, 1e-12);
  close(h.waterplaneArea, 2000, 1e-12); close(h.wettedArea, 2000 + 2 * 4 * 120, 1e-12);
  close(h.GMt, 2 + 400 / 48 - 5, 1e-12); close(h.KMl, 2 + 10000 / 48, 1e-12);
  close(h.TPC, 1.025 * 2000 / 100, 1e-12);
  close(h.stiffness[2][2], rho * g * 2000, 1e-12);
  close(h.stiffness[3][3], rho * g * 8000 * h.GMt, 1e-10);
  close(h.stiffness[4][4], rho * g * 8000 * h.GMl, 1e-10);
  close(h.MCT, (rho * g * 8000 * h.GMl) / 100, 1e-10);
  near(h.cof[0], 0, 1e-9); near(h.stiffness[2][4], 0, 1e-3);
  const cf = hs.boxHydrostatics({ L: 100, B: 20, T: 4, KG: 5, rho, g });
  close(cf.GMt, h.GMt, 1e-12); close(cf.C44, h.stiffness[3][3], 1e-10); close(cf.KMt, h.KMt, 1e-12);
});

test('Orientation-independent: an inward-oriented mesh gives the same hydrostatics', () => {
  const m = barge();
  const flipped = { positions: m.positions, indices: Uint32Array.from(m.indices) };
  for (let i = 0; i < flipped.indices.length; i += 3) { const t = flipped.indices[i + 1]; flipped.indices[i + 1] = flipped.indices[i + 2]; flipped.indices[i + 2] = t; }
  const a = hs.meshHydrostatics({ mesh: m, draft: 3, heel: 0.2, trim: 0.03 });
  const b = hs.meshHydrostatics({ mesh: flipped, draft: 3, heel: 0.2, trim: 0.03 });
  close(b.volume, a.volume, 1e-12); close(b.Ixx, a.Ixx, 1e-12); close(b.cob[1], a.cob[1], 1e-12);
});

test('Heeled box (wall-sided range): immersed section is a trapezium — volume and centre of buoyancy by hand', () => {
  // heel phi about the keel centreline with the body origin T0 below the waterline:
  // section area = B·T0/cos(phi) (trapezium of mean depth T0/cos phi measured along the body z axis)
  const phi = 0.2; const T0 = 4;
  const h = hs.meshHydrostatics({ mesh: barge(), draft: T0, heel: phi });
  close(h.volume, (100 * 20 * T0) / Math.cos(phi), 1e-9);
  // body-frame centroid of the trapezium: mean depth d = T0/cos, y_b = -B² tan(phi)/(12 d), z_b = d/2 + B² tan²(phi)/(24 d)
  const d = T0 / Math.cos(phi);
  close(h.cobBody[1], (-400 * Math.tan(phi)) / (12 * d), 1e-9);
  close(h.cobBody[2], d / 2 + (400 * Math.tan(phi) ** 2) / (24 * d), 1e-9);
  // heeled waterplane breadth B/cos(phi)
  close(h.waterplaneArea, (100 * 20) / Math.cos(phi), 1e-9);
});

test('Half-submerged sphere mesh: V = 2/3 pi R³, centre of buoyancy 3R/8 below the waterplane, Aw = pi R²', () => {
  const R = 5;
  const h = hs.meshHydrostatics({ mesh: sphereMesh(R, 64, 128, R), draft: R });
  close(h.volume, (2 / 3) * Math.PI * R ** 3, 3e-3);
  close(-h.cob[2], (3 * R) / 8, 3e-3);
  close(h.waterplaneArea, Math.PI * R * R, 2e-3);
  close(h.Ixx, (Math.PI * R ** 4) / 4, 4e-3);
  close(h.wettedArea, 2 * Math.PI * R * R, 2e-3);
  // fully submerged: total volume and no waterplane
  const full = hs.clipIntegrals(sphereMesh(R, 64, 128, R), 3 * R, 0, 0);
  close(full.V, (4 / 3) * Math.PI * R ** 3, 3e-3); near(full.Aw, 0, 1e-6);
});

test('Vertical cylinder mesh vs closed form: V = pi D² T/4, KB = T/2, BM = D²/16T', () => {
  const D = 12; const T = 30;
  const h = hs.meshHydrostatics({ mesh: hs.cylinderMesh({ D, zMin: 0, zMax: 40, nTheta: 96 }), draft: T, cog: [0, 0, 12] });
  const cf = hs.verticalCylinderHydrostatics({ D, T, KG: 12 });
  close(h.volume, cf.volume, 1e-6); close(h.KB, cf.KB, 1e-6);
  close(h.BMt, D * D / (16 * T), 2e-3); close(h.GMt, cf.GM, 1e-3);
});

test('Four-column body: mesh hydrostatics match the parallel-axis closed form', () => {
  const cols = [[30, 30], [30, -30], [-30, 30], [-30, -30]];
  const mesh = hs.mergeMeshes([
    ...cols.map(([x, y]) => hs.cylinderMesh({ D: 10, zMin: 8, zMax: 40, x0: x, y0: y, nTheta: 96 })),
    hs.boxMesh({ xMin: -38, xMax: 38, yMin: 22, yMax: 38, zMin: 0, zMax: 8 }),
    hs.boxMesh({ xMin: -38, xMax: 38, yMin: -38, yMax: -22, zMin: 0, zMax: 8 }),
  ]);
  const h = hs.meshHydrostatics({ mesh, draft: 20, cog: [0, 0, 18] });
  const cf = hs.multiColumnHydrostatics({
    columns: cols.map(([x, y]) => ({ x, y, D: 10, zBottom: 8 })),
    submerged: [{ volume: 76 * 16 * 8, y: 30, z: 4 }, { volume: 76 * 16 * 8, y: -30, z: 4 }], T: 20, KG: 18,
  });
  close(h.volume, cf.volume, 1e-6); close(h.KB, cf.KB, 1e-6);
  close(h.BMt, cf.BMt, 1e-4); close(h.GMl, cf.GMl, 1e-3);
  // hand value: Aw = 4·pi·25, Ixx = 4 (pi 10⁴/64 + pi·25·900)
  close(cf.Ixx, 4 * ((Math.PI * 1e4) / 64 + Math.PI * 25 * 900), 1e-12);
});

test('Equilibrium: transverse and longitudinal weight offsets give tan(heel) = e/GM and tan(trim) = e/GMl', () => {
  const rho = 1025; const mass = rho * 8000;
  const GMt = 2 + 400 / 48 - 5; const GMl = 2 + 10000 / 48 - 5;
  const eq = hs.solveEquilibrium({ mesh: barge(), mass, cog: [1.0, 0.1, 5], rho });
  assert.ok(eq.converged);
  // weight to port (+y) heels port side down = negative heel; weight forward trims bow down = positive trim.
  // wall-sided: e = tan(phi) (GM + ½ BM tan² phi) -> correction of order 1e-4 here
  close(Math.tan(-eq.heel), 0.1 / GMt, 2e-3);
  close(Math.tan(eq.trim), 1.0 / GMl, 2e-3);
  near(hs.meshHydrostatics({ mesh: barge(), draft: eq.draft, heel: eq.heel, trim: eq.trim, rho }).displacement / mass, 1, 1e-8);
});

test('Equilibrium recovery: a prescribed floating position is recovered from mass and centre of gravity', () => {
  const mesh = barge(); const rho = 1025;
  const target = { draft: 4.2, heel: 0.12, trim: 0.02 };
  const h = hs.meshHydrostatics({ mesh, ...target, rho });
  // put G on the vertical through B, 0.5 m above it (GM stays positive: BM ≈ 8 m)
  const cog = hs.earthToBody([h.cob[0], h.cob[1], h.cob[2] + 0.5], target.draft, target.heel, target.trim);
  const eq = hs.solveEquilibrium({ mesh, mass: h.displacement, cog, rho });
  assert.ok(eq.converged);
  near(eq.draft, target.draft, 1e-6); near(eq.heel, target.heel, 1e-6); near(eq.trim, target.trim, 1e-6);
});

test('draftForVolume and a body that cannot float', () => {
  near(hs.draftForVolume({ mesh: barge(), volume: 5000 }), 2.5, 1e-9);
  assert.throws(() => hs.draftForVolume({ mesh: barge(), volume: 30000 }), RangeError);
});

test('Free-surface correction: rectangular tank i = l b³/12; mesh tank agrees', () => {
  const r = hs.freeSurfaceCorrection({ tanks: [{ l: 10, b: 8, rho: 1000 }], volume: 8000, rho: 1025 });
  close(r.perTank[0].iT, (10 * 512) / 12, 1e-12);
  close(r.GGt, (1000 * 10 * 512) / 12 / (1025 * 8000), 1e-12);
  const tank = hs.boxMesh({ xMin: 5, xMax: 15, yMin: 1, yMax: 9, zMin: 0, zMax: 6 });
  const m = hs.freeSurfaceCorrection({ tanks: [{ mesh: tank, fillLevel: 3, rho: 1000 }], volume: 8000 });
  close(m.GGt, r.GGt, 1e-9); close(m.perTank[0].iL, (8 * 1000) / 12, 1e-9);
  // a pressed-up (full) tank has no free surface
  near(hs.freeSurfaceCorrection({ tanks: [{ mesh: tank, fillLevel: 7 }], volume: 8000 }).GGt, 0, 1e-12);
});
