import test from 'node:test';
import assert from 'node:assert/strict';
import {
  pileLocalStiffness, analysePileGroup, pMultiplier, P_MULTIPLIER_TABLE_AASHTO, axialInteractionFactor,
  lateralInteractionFactor, groupStiffnessFromInteraction, converseLabarre, groupTorsionalCapacity,
  pileUtilisationTable,
} from '../../src/engine/geotech/pilegroup.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-9), `${a} vs ${b}`);
const stiff = { KHH: 2e8, KHM: -9e8, KMM: 8e9, Kaxial: 1.5e9, Ktorsion: 5e8 };

function grid(nx, ny, s, extra = {}) {
  const piles = [];
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) piles.push({ id: `${i}-${j}`, x: (i - (nx - 1) / 2) * s, y: (j - (ny - 1) / 2) * s, stiffness: stiff, ...extra });
  return piles;
}

test('vertical piles, pinned heads, under V and M: P = V/n ± M x / sum(x^2)', () => {
  const piles = grid(3, 3, 3, { stiffness: { ...stiff, head: 'pinned' } });
  const V = 18e6, My = 30e6, Mx = -12e6;
  const r = analysePileGroup({ piles, load: { Fz: -V, My, Mx } });
  const sx2 = piles.reduce((s, p) => s + p.x * p.x, 0), sy2 = piles.reduce((s, p) => s + p.y * p.y, 0);
  assert.equal(sx2, 54);
  r.piles.forEach((p, i) => {
    // +My (right-hand rule, z up) pushes the piles at +x down; +Mx lifts the piles at +y
    close(p.axial, V / 9 + (My * piles[i].x) / sx2 - (Mx * piles[i].y) / sy2, 1e-9);
    assert.ok(Math.abs(p.shear) < 1e-3 && Math.abs(p.moment) < 1e-3);
  });
  assert.ok(r.equilibriumError < 1e-12);
  close(r.displacement.uz, -V / (9 * stiff.Kaxial));
  close(r.displacement.ry, My / (stiff.Kaxial * sx2));
});

test('fixed-head vertical group under horizontal load: equal shears, global equilibrium, symmetric positive-definite K', () => {
  const piles = grid(2, 3, 4);
  const load = { Fx: 3e6, Fy: -1e6, Fz: -20e6, Mx: 4e6, My: 9e6, Mz: 2e6 };
  const r = analysePileGroup({ piles, load });
  assert.ok(r.equilibriumError < 1e-10);
  for (let a = 0; a < 6; a++) for (let b = 0; b < a; b++) assert.ok(Math.abs(r.K[a][b] - r.K[b][a]) <= 1e-9 * Math.abs(r.K[a][a]));
  // work of the loads is positive (positive-definite stiffness)
  const d = r.displacement;
  assert.ok(load.Fx * d.ux + load.Fy * d.uy + load.Fz * d.uz + load.Mx * d.rx + load.My * d.ry + load.Mz * d.rz > 0);
  // pure Fx on a doubly symmetric group: each pile takes Fx/n in shear
  const h = analysePileGroup({ piles, load: { Fx: 3e6 } });
  for (const p of h.piles) close(p.shear1, 3e6 / 6, 1e-9);
  // hand check of the group terms: Kxx = n KHH, Kzz = n Kaxial, Kzz-torsion = n Ktors + KHH sum(r^2)
  close(h.K[0][0], 6 * stiff.KHH); close(h.K[2][2], 6 * stiff.Kaxial);
  close(h.K[5][5], 6 * stiff.Ktorsion + stiff.KHH * piles.reduce((s, p) => s + p.x ** 2 + p.y ** 2, 0));
  // rotation stiffness about y: n KMM + Kaxial sum(x^2); coupling Fx–ry = −n KHM (sign for z up)
  close(h.K[4][4], 6 * stiff.KMM + stiff.Kaxial * piles.reduce((s, p) => s + p.x ** 2, 0));
  close(Math.abs(h.K[0][4]), 6 * Math.abs(stiff.KHM));
  // fixed head: the cap rotation is restrained by the axial couple, so the pile moment is close to −KHM/KHH × shear
  const single = analysePileGroup({ piles: [{ x: 0, y: 0, stiffness: stiff }], load: { Fx: 1e6 } });
  // single free-standing pile: free-head flexibility fHH = KMM / (KHH KMM − KHM²)
  close(single.displacement.ux, (1e6 * stiff.KMM) / (stiff.KHH * stiff.KMM - stiff.KHM ** 2), 1e-9);
  assert.ok(Math.abs(single.piles[0].moment) < 1e-3);
  // the head tilts in the direction of the load: rotation about +y for a load along +x
  assert.ok(single.displacement.ry > 0);
});

test('raked piles: symmetric A-frame carries horizontal load axially; equilibrium holds', () => {
  const b = 15, ax = { KHH: 1e3, KHM: 0, KMM: 1e3, Kaxial: 2e9, Ktorsion: 1e3 };
  // two piles raked outwards in the x–z plane with heads at the same point, plus two in the y–z plane
  const piles = [
    { x: 0, y: 0, batterDeg: b, azimuthDeg: 0, stiffness: ax }, { x: 0, y: 0, batterDeg: b, azimuthDeg: 180, stiffness: ax },
    { x: 0, y: 0, batterDeg: b, azimuthDeg: 90, stiffness: ax }, { x: 0, y: 0, batterDeg: b, azimuthDeg: 270, stiffness: ax },
  ];
  const Fx = 1e6;
  const r = analysePileGroup({ piles, load: { Fx } });
  assert.ok(r.equilibriumError < 1e-9);
  // truss solution: N = ± Fx / (2 sin b); the +x raked pile (toe towards +x) is compressed by a load along +x
  close(r.piles[0].axial, Fx / (2 * Math.sin((b * Math.PI) / 180)), 1e-4);
  close(r.piles[1].axial, -Fx / (2 * Math.sin((b * Math.PI) / 180)), 1e-4);
  assert.ok(Math.abs(r.piles[2].axial) < 1);
  // vertical load: N = V / (4 cos b)
  const v = analysePileGroup({ piles, load: { Fz: -8e6 } });
  for (const p of v.piles) close(p.axial, 8e6 / (4 * Math.cos((b * Math.PI) / 180)), 1e-5);
});

test('local stiffness: sign convention and pinned condensation', () => {
  const K = pileLocalStiffness(stiff);
  assert.equal(K[0][4], -stiff.KHM); assert.equal(K[1][3], stiff.KHM);
  const Kp = pileLocalStiffness({ ...stiff, head: 'pinned' });
  close(Kp[0][0], stiff.KHH - stiff.KHM ** 2 / stiff.KMM);
  assert.equal(Kp[3][3], 0);
  assert.throws(() => pileLocalStiffness({ ...stiff, KHM: 1 }), RangeError);
  assert.throws(() => analysePileGroup({ piles: [{ x: 0, y: 0, stiffness: { ...stiff, head: 'pinned', Ktorsion: 0 } }], load: { My: 1 } }), /mechanism/);
});

test('group reduction helpers: p-multipliers, interaction factors, Converse-Labarre, torsion, utilisation', () => {
  assert.deepEqual(P_MULTIPLIER_TABLE_AASHTO.rows, [[0.8, 0.4, 0.3], [1.0, 0.85, 0.7]]);
  close(pMultiplier({ spacingOverD: 3, row: 1 }), 0.8);
  close(pMultiplier({ spacingOverD: 4, row: 2 }), 0.625); // midway between 0.40 and 0.85
  close(pMultiplier({ spacingOverD: 5, row: 4 }), 0.7);
  close(pMultiplier({ spacingOverD: 2, row: 3 }), 0.3);
  close(pMultiplier({ spacingOverD: 4, row: 1, table: { spacingOverD: [3, 5], rows: [[0.5], [1]] } }), 0.75);
  // Randolph & Wroth: rm = 2.5 (1 − 0.3) 30 = 52.5 m; alpha = ln(52.5/3)/ln(52.5/0.5)
  const al = axialInteractionFactor({ s: 3, D: 1, L: 30, nu: 0.3 });
  close(al, Math.log(52.5 / 3) / Math.log(105), 1e-12);
  assert.equal(axialInteractionFactor({ s: 60, D: 1, L: 30, nu: 0.3 }), 0);
  // two identical piles under a rigid cap: K = 2 K1 / (1 + alpha)
  const g = groupStiffnessFromInteraction({ piles: [{ x: 0, y: 0 }, { x: 3, y: 0 }], K1: 1e9, alphaFn: (s) => axialInteractionFactor({ s, D: 1, L: 30, nu: 0.3 }) });
  close(g.Kgroup, 2e9 / (1 + al), 1e-12); close(g.efficiency, 1 / (1 + al), 1e-12);
  // 3 in a row: the centre pile carries less than the outer piles
  const g3 = groupStiffnessFromInteraction({ piles: [{ x: -3, y: 0 }, { x: 0, y: 0 }, { x: 3, y: 0 }], K1: 1e9, alphaFn: (s) => axialInteractionFactor({ s, D: 1, L: 30, nu: 0.3 }) });
  assert.ok(g3.loadShare[1] < g3.loadShare[0]); close(g3.loadShare[0], g3.loadShare[2], 1e-12);
  // Randolph lateral factor: 0.6 rho_c (Ep/Gc)^(1/7) (r0/s) (1 + cos² psi), with the close-spacing correction
  const EI = (30e9 * Math.PI) / 64; // solid 1 m pile, Ep = 30 GPa
  close(lateralInteractionFactor({ s: 8, D: 1, EI, Gc: 30e6, psiDeg: 90 }), 0.6 * Math.pow(1000, 1 / 7) * (0.5 / 8), 1e-9);
  const a0 = 0.6 * Math.pow(1000, 1 / 7) * (0.5 / 2) * 2;
  close(lateralInteractionFactor({ s: 2, D: 1, EI, Gc: 30e6, psiDeg: 0 }), 1 - 2 / Math.sqrt(27 * a0), 1e-9);
  // Converse–Labarre 3 × 3, D = 0.5 m, s = 1.5 m: theta = 18.435°, Eg = 1 − 18.435 × 12 / 810
  close(converseLabarre({ m: 3, n: 3, D: 0.5, s: 1.5 }), 1 - (18.4349488 * 12) / 810, 1e-8);
  // 4 piles on a 4 m square, Hu = 1 MN: T = 4 × 1 MN × 2.828 m (+ pile torsion)
  close(groupTorsionalCapacity({ piles: grid(2, 2, 4).map((p) => ({ ...p, Hu: 1e6, Tu: 1e5 })) }).Tult, 4e6 * Math.SQRT2 * 2 + 4e5, 1e-12);
  const r = analysePileGroup({ piles: grid(2, 2, 4), load: { Fz: -8e6, My: 40e6, Fx: 1e6 } });
  const tab = pileUtilisationTable({ piles: r.piles, capacity: { compression: 10e6, tension: 4e6, lateral: 1e6 } });
  assert.equal(tab.rows.length, 4);
  close(tab.maxUtilisation, Math.max(...tab.rows.map((x) => x.uc)));
  const tensionPile = tab.rows.find((x) => x.axial < 0);
  if (tensionPile) close(tensionPile.ucAxial, -tensionPile.axial / 4e6);
});
