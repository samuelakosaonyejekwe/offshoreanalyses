import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  boxMesh, boundaryFaces, selectFaces, faceSetArea, meshStats, voxelMesh, surfaceVolume, solveSolid, modalSolid,
  principalStresses, solidD, shape3d,
} from '../../src/engine/structural/solidfem.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);

const wobble = (lx, ly, lz, amp) => (x, y, z) => {
  const s = Math.sin((Math.PI * x) / lx) * Math.sin((Math.PI * y) / ly) * Math.sin((Math.PI * z) / lz);
  return [x + amp * lx * s, y - 0.8 * amp * ly * s, z + 0.6 * amp * lz * s];
};

/** Closed triangulated sphere by recursive subdivision of an octahedron (outward orientation). */
function sphereSurface(radius, levels, centre = [0, 0, 0]) {
  let verts = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
  let tris = [[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4], [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]];
  for (let l = 0; l < levels; l++) {
    const cache = new Map();
    const mid = (a, b) => {
      const key = a < b ? `${a}_${b}` : `${b}_${a}`;
      if (!cache.has(key)) {
        const m = [0, 1, 2].map((k) => verts[a][k] + verts[b][k]);
        const n = Math.hypot(...m);
        verts.push(m.map((v) => v / n));
        cache.set(key, verts.length - 1);
      }
      return cache.get(key);
    };
    const next = [];
    for (const [a, b, c] of tris) {
      const ab = mid(a, b); const bc = mid(b, c); const ca = mid(c, a);
      next.push([a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]);
    }
    tris = next;
  }
  return {
    kind: 'surface', name: 'sphere',
    positions: Float32Array.from(verts.flatMap((v) => v.map((x, k) => radius * x + centre[k]))),
    indices: Uint32Array.from(tris.flat()),
  };
}

function boxSurface(lx, ly, lz, o = [0, 0, 0]) {
  const p = [];
  for (const z of [0, lz]) for (const y of [0, ly]) for (const x of [0, lx]) p.push(x + o[0], y + o[1], z + o[2]);
  // vertex id = x + 2 y + 4 z
  const quads = [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]];
  const idx = [];
  for (const q of quads) idx.push(q[0], q[1], q[2], q[0], q[2], q[3]);
  return { kind: 'surface', name: 'box', positions: Float32Array.from(p), indices: Uint32Array.from(idx) };
}

test('shape functions: partition of unity and derivative sums', () => {
  for (const [type, n] of [['tet4', 4], ['tet10', 10], ['hex8', 8]]) {
    const N = new Float64Array(n); const a = new Float64Array(n); const b = new Float64Array(n); const c = new Float64Array(n);
    shape3d(type, 0.21, 0.17, 0.33, N, a, b, c);
    const sum = (v) => v.reduce((s, x) => s + x, 0);
    close(sum(N), 1, 1e-14);
    assert.ok(Math.abs(sum(a)) < 1e-14 && Math.abs(sum(b)) < 1e-14 && Math.abs(sum(c)) < 1e-14);
  }
});

test('patch test: hex8 (incompatible and full), tet4 and tet10 reproduce constant strain on distorted meshes', () => {
  const E = 200e9; const nu = 0.29;
  const A = [[1e-4, -2e-4, 0.5e-4], [1.5e-4, 2e-4, -1e-4], [-0.5e-4, 1e-4, 3e-4]];
  const disp = (x, y, z) => [0, 1, 2].map((i) => 1e-5 * (i + 1) + A[i][0] * x + A[i][1] * y + A[i][2] * z);
  const eps = [A[0][0], A[1][1], A[2][2], A[0][1] + A[1][0], A[1][2] + A[2][1], A[2][0] + A[0][2]];
  const D = solidD(E, nu);
  const sig = [0, 1, 2, 3, 4, 5].map((r) => eps.reduce((s, e, c) => s + D[r * 6 + c] * e, 0));
  for (const [type, hex] of [['hex8', 'incompatible'], ['hex8', 'full'], ['tet4'], ['tet10']]) {
    const model = boxMesh({ lx: 1, ly: 0.8, lz: 0.6, nx: 3, ny: 3, nz: 2, type, distort: wobble(1, 0.8, 0.6, 0.05) });
    const bn = boundaryFaces(model).nodes;
    const p = model.positions;
    const fix = Array.from(bn).map((n) => { const d = disp(p[3 * n], p[3 * n + 1], p[3 * n + 2]); return { nodes: [n], ux: d[0], uy: d[1], uz: d[2] }; });
    const r = solveSolid({ model, E, nu, fix, hex: hex || 'incompatible', solver: 'skyline' });
    let err = 0;
    for (let n = 0; n < p.length / 3; n++) {
      const d = disp(p[3 * n], p[3 * n + 1], p[3 * n + 2]);
      for (let k = 0; k < 3; k++) err = Math.max(err, Math.abs(r.u[3 * n + k] - d[k]));
    }
    assert.ok(err < 1e-12, `${type} ${hex || ''}: displacement error ${err}`);
    const names = ['sxx', 'syy', 'szz', 'sxy', 'syz', 'szx'];
    for (let n = 0; n < p.length / 3; n++) for (let c = 0; c < 6; c++) assert.ok(Math.abs(r.stress[names[c]][n] - sig[c]) < 1e-9 * E * 1e-4, `${type} ${names[c]}`);
    assert.ok(r.equilibriumResidual < 1e-9, `residual ${r.equilibriumResidual}`);
  }
});

test('cantilever with solid elements converges to Timoshenko beam theory', () => {
  const E = 210e9; const nu = 0.3; const L = 1; const b = 0.1; const h = 0.1; const P = 1000;
  const I = (b * h ** 3) / 12; const G = E / (2 * (1 + nu));
  const exact = (P * L ** 3) / (3 * E * I) + (P * L) / ((5 / 6) * G * b * h);
  const run = (type, nx, n, hex) => {
    const model = boxMesh({ lx: L, ly: b, lz: h, nx, ny: n, nz: n, type });
    const r = solveSolid({ model, E, nu, hex: hex || 'incompatible', fix: [{ set: 'xmin', ux: 0, uy: 0, uz: 0 }], faceLoads: [{ set: 'xmax', traction: [0, 0, -P / (b * h)] }], solver: 'skyline' });
    close(r.reactionSum.fz, P, 1e-8, 'reaction');
    const tip = model.sets.xmax.nodes;
    let w = 0;
    for (const n2 of tip) w += r.u[3 * n2 + 2];
    return { tip: -w / tip.length, r };
  };
  // Incompatible-mode hexahedra: within 2 % with a single element through the depth (the fully clamped root
  // restrains Poisson contraction, so the 3-D answer is slightly stiffer than beam theory).
  close(run('hex8', 10, 1).tip, exact, 0.02, 'hex8 incompatible, 10×1×1');
  close(run('hex8', 20, 2).tip, exact, 0.02, 'hex8 incompatible, 20×2×2');
  // The plain trilinear element locks in bending on the same coarse mesh.
  assert.ok(run('hex8', 10, 1, 'full').tip < 0.75 * exact);
  // Quadratic tetrahedra.
  close(run('tet10', 8, 1).tip, exact, 0.035, 'tet10 coarse');
  close(run('tet10', 12, 2).tip, exact, 0.02, 'tet10 refined');
  // Linear tetrahedra lock badly in bending: far too stiff on coarse meshes, converging monotonically from
  // below with an observed order that climbs towards the asymptotic 2 as the mesh is refined.
  const t1 = run('tet4', 8, 1).tip; const t2 = run('tet4', 16, 2).tip; const t3 = run('tet4', 32, 4).tip;
  assert.ok(t1 < t2 && t2 < t3 && t3 < exact, `tet4 sequence ${t1}, ${t2}, ${t3}`);
  const o1 = Math.log2((exact - t1) / (exact - t2)); const o2 = Math.log2((exact - t2) / (exact - t3));
  assert.ok(o1 > 0.4 && o2 > 1.0 && o2 > o1 && o2 < 2.2, `tet4 observed orders ${o1}, ${o2}`);
  // Bending stress at the root top fibre: σ = P L / Z (sampled one cell away from the clamp).
  const fine = run('hex8', 20, 2);
  const model = boxMesh({ lx: L, ly: b, lz: h, nx: 20, ny: 2, nz: 2 });
  let sx = 0;
  const p = model.positions;
  for (let n = 0; n < p.length / 3; n++) if (Math.abs(p[3 * n] - 0.25) < 1e-9 && Math.abs(p[3 * n + 2] - h) < 1e-9 && Math.abs(p[3 * n + 1] - b / 2) < 1e-9) sx = fine.r.stress.sxx[n];
  close(sx, (P * (L - 0.25) * (h / 2)) / I, 0.03, 'bending stress');
});

test('loads: face pressure, gravity, hydrostatic pressure, foundation springs — closed-form resultants', () => {
  const E = 30e9; const nu = 0.2;
  const model = boxMesh({ lx: 2, ly: 1, lz: 3, nx: 4, ny: 2, nz: 6 });
  // Uniaxial compression by pressure on top with frictionless supports: σzz = −p, uz(top) = −p L/E.
  const p = 2e6;
  const r1 = solveSolid({ model, E, nu, fix: [{ set: 'zmin', uz: 0 }, { set: 'xmin', ux: 0 }, { set: 'ymin', uy: 0 }], faceLoads: [{ set: 'zmax', pressure: p }], solver: 'skyline' });
  for (const n of model.sets.zmax.nodes) close(r1.u[3 * n + 2], (-p * 3) / E, 1e-9, 'shortening');
  close(r1.stress.szz[5], -p, 1e-9, 'σzz');
  close(r1.maxVonMises, p, 1e-9, 'von Mises');
  close(r1.reactionSum.fz, p * 2 * 1, 1e-10, 'reaction');
  close(r1.strainEnergy, ((p * p) / (2 * E)) * 6, 1e-9, 'strain energy = σ²V/2E');
  // Self-weight: total reaction = ρ g V.
  const r2 = solveSolid({ model, E, nu, fix: [{ set: 'zmin', ux: 0, uy: 0, uz: 0 }], gravity: { rho: 2400, g: [0, 0, -9.81] }, solver: 'skyline' });
  close(r2.reactionSum.fz, 2400 * 9.81 * 6, 1e-10, 'weight');
  // Hydrostatic pressure on the face x = 0 with the free surface at z = 2.4: F = ½ ρ g d² b, acting in +x.
  const r3 = solveSolid({ model, E, nu, fix: [{ set: 'zmin', ux: 0, uy: 0, uz: 0 }], hydrostatic: [{ set: 'xmin', rho: 1025, zSurface: 2.4 }], solver: 'skyline' });
  // The waterline cuts through a cell layer (2×2 Gauss on a kinked integrand), hence the looser tolerance.
  close(r3.appliedSum.fx, 0.5 * 1025 * 9.80665 * 2.4 * 2.4 * 1, 0.01, 'hydrostatic thrust');
  close(r3.reactionSum.fx, -r3.appliedSum.fx, 1e-9, 'equilibrium');
  // Block on a Winkler foundation under uniform pressure: settlement = p/k, rigid-body.
  const k = 5e7;
  const r4 = solveSolid({ model, E, nu, fix: [{ nodes: [0], ux: 0, uy: 0 }, { nodes: [4], uy: 0 }], foundation: [{ set: 'zmin', kz: k }], faceLoads: [{ set: 'zmax', pressure: p }], solver: 'skyline' });
  const bottom = model.sets.zmin.nodes;
  let mean = 0;
  for (const n of bottom) mean += r4.u[3 * n + 2] / bottom.length;
  close(mean, -p / k, 1e-3, 'foundation settlement');
  close(r4.reactionSum.fz, p * 2, 1e-9, 'spring reaction');
  assert.ok(r4.equilibriumResidual < 1e-9);
});

test('voxel mesher: exact for an aligned box, including rays through triangle diagonals and vertices', () => {
  // h = 0.25 puts cell-centre rays exactly on the face diagonals of the triangulated unit cube.
  const m = voxelMesh({ surface: boxSurface(1, 1, 1), h: 0.25 });
  assert.equal(m.stats.nCells, 64);
  assert.equal(m.stats.nNodes, 125);
  assert.equal(m.stats.voteDisagreement, 0);
  close(m.stats.volume, 1, 1e-12);
  close(faceSetArea(m, m.sets.zmin).area, 1, 1e-12, 'zmin set');
  close(faceSetArea(m, m.sets['+x']).area, 1, 1e-12, '+x set');
  close(faceSetArea(m, m.sets.boundary).area, 6, 1e-12, 'boundary');
  // Reversed triangle orientation does not matter for parity classification.
  const s = boxSurface(2, 1, 0.5, [3, -1, 7]);
  for (let t = 0; t < s.indices.length; t += 3) { const a = s.indices[t]; s.indices[t] = s.indices[t + 1]; s.indices[t + 1] = a; }
  const m2 = voxelMesh({ surface: s, h: 0.125 });
  assert.equal(m2.stats.nCells, 16 * 8 * 4);
  // Odd division: centres hit the cube's vertices' projections and edges nowhere, the box is still filled.
  const m3 = voxelMesh({ surface: boxSurface(1, 1, 1), targetCells: 1000 });
  assert.equal(m3.stats.nCells, 1000);
  const st = meshStats(m3);
  close(st.volume, 1, 1e-12);
  close(st.maxAspect, 1, 1e-9);
  assert.throws(() => voxelMesh({ surface: boxSurface(1, 1, 1), h: 1e-3 }), RangeError);
  assert.throws(() => voxelMesh({ surface: { positions: new Float32Array(9), indices: new Uint32Array([0, 1, 2]) }, h: 1 }), RangeError);
});

test('voxel mesher: sphere volume converges to 4πr³/3 as h → 0', () => {
  const r = 1;
  const surf = sphereSurface(r, 5, [0.013, -0.027, 0.041]);
  const sv = surfaceVolume(surf);
  const exact = (4 / 3) * Math.PI * r ** 3;
  close(sv.volume, exact, 2e-3, 'tessellated sphere volume');
  close(sv.area, 4 * Math.PI, 2e-3, 'tessellated sphere area');
  const errs = [];
  for (const h of [0.2, 0.1, 0.05, 0.025]) {
    const m = voxelMesh({ surface: surf, h });
    errs.push(Math.abs(m.stats.volume - sv.volume) / sv.volume);
    assert.equal(m.stats.removedCells, 0);
    assert.ok(m.stats.voteDisagreement < 1e-3, `vote disagreement ${m.stats.voteDisagreement}`);
  }
  // Centre-sampling error of a smooth body decays at least linearly in h (typically faster, with scatter).
  assert.ok(errs[0] < 0.08 && errs[1] < 0.02 && errs[2] < 0.008 && errs[3] < 0.003, `errors ${errs}`);
  // User plane set: equatorial cut does not exist on a full sphere, the top cap does.
  const m = voxelMesh({ surface: surf, h: 0.1, planes: [{ name: 'cap', point: [0, 0, 1], normal: [0, 0, 1], tol: 0.12 }] });
  assert.ok(m.sets.cap.quad4.length > 0);
  const capArea = faceSetArea(m, m.sets.cap);
  close(capArea.areaVector[2], capArea.area, 1e-12, 'cap faces all face +z');
});

test('voxel mesher: leaky and disconnected input is handled without failure', () => {
  // Remove two triangles of the box: the three-axis vote still fills the body.
  const s = boxSurface(1, 1, 1);
  const leaky = { positions: s.positions, indices: s.indices.slice(6) };
  const m = voxelMesh({ surface: leaky, h: 0.1 });
  assert.ok(m.stats.nCells > 0.85 * 1000 && m.stats.nCells <= 1000, `cells ${m.stats.nCells}`);
  assert.ok(m.stats.voteDisagreement > 0.05, 'the leak is reported');
  // Two separate bodies: only the largest is kept.
  const a = boxSurface(1, 1, 1); const b = boxSurface(0.4, 0.4, 0.4, [2, 0, 0]);
  const both = {
    positions: Float32Array.from([...a.positions, ...b.positions]),
    indices: Uint32Array.from([...a.indices, ...Array.from(b.indices, (i) => i + 8)]),
  };
  const m2 = voxelMesh({ surface: both, h: 0.1 });
  assert.equal(m2.stats.nCells, 1000);
  assert.equal(m2.stats.removedCells, 64);
});

test('voxel model solve: matrix-free CG agrees with the direct solver and satisfies equilibrium', () => {
  const E = 70e9; const nu = 0.33;
  const m = voxelMesh({ surface: boxSurface(1, 0.2, 0.2), h: 0.05 });
  const args = { model: m, E, nu, fix: [{ set: 'xmin', ux: 0, uy: 0, uz: 0 }], faceLoads: [{ set: 'xmax', traction: [0, 0, -1e5] }, { set: 'zmax', pressure: 2e4 }] };
  const d = solveSolid({ ...args, solver: 'skyline' });
  const it = solveSolid({ ...args, solver: 'ebe', tol: 1e-10 });
  assert.equal(it.solver, 'ebe');
  assert.ok(it.iterations > 0);
  let diff = 0;
  for (let i = 0; i < d.u.length; i++) diff = Math.max(diff, Math.abs(d.u[i] - it.u[i]));
  assert.ok(diff < 1e-6 * d.maxDisp, `max difference ${diff} vs ${d.maxDisp}`);
  close(it.reactionSum.fz, 1e5 * 0.04 + 2e4 * 0.2, 1e-7, 'reaction');
  assert.ok(it.equilibriumResidual < 1e-8, `residual ${it.equilibriumResidual}`);
  close(it.maxVonMises, d.maxVonMises, 1e-6);
  // CSR + preconditioned CG path
  const cg = solveSolid({ ...args, solver: 'pcg' });
  close(cg.maxDisp, d.maxDisp, 1e-6, 'pcg');
});

test('voxelised sphere under self-weight on a cut plane: mesh-convergent reaction and displacement', () => {
  const surf = sphereSurface(1, 4);
  const res = [];
  for (const h of [0.25, 0.125]) {
    const m = voxelMesh({ surface: surf, h });
    const r = solveSolid({ model: m, E: 10e9, nu: 0.25, fix: [{ set: 'zmin', ux: 0, uy: 0, uz: 0 }], gravity: { rho: 2000, g: [0, 0, -9.81] } });
    close(r.reactionSum.fz, 2000 * 9.81 * m.stats.volume, 1e-7, 'weight of the voxel body');
    assert.ok(r.equilibriumResidual < 1e-6);
    res.push(r.maxDisp);
  }
  assert.ok(res.every((v) => v > 0 && Number.isFinite(v)));
});

test('modal analysis: cantilever bending frequency f1 = (1.8751²/2π) √(EI/(ρ A L⁴))', () => {
  const E = 210e9; const nu = 0.3; const rho = 7850; const L = 1; const b = 0.05; const h = 0.02;
  const model = boxMesh({ lx: L, ly: b, lz: h, nx: 24, ny: 2, nz: 1 });
  const r = modalSolid({ model, E, nu, rho, fix: [{ set: 'xmin', ux: 0, uy: 0, uz: 0 }], nModes: 3 });
  const f = (I, A) => (1.8751 ** 2 / (2 * Math.PI)) * Math.sqrt((E * I) / (rho * A * L ** 4));
  close(r.totalMass, rho * L * b * h, 1e-12, 'mass');
  close(r.frequencies[0], f((b * h ** 3) / 12, b * h), 0.02, 'weak-axis bending');
  close(r.frequencies[1], f((h * b ** 3) / 12, b * h), 0.03, 'strong-axis bending');
  // Second mode in the weak plane: (4.6941/1.8751)² × f1
  close(r.frequencies[2], f((b * h ** 3) / 12, b * h) * (4.6941 / 1.8751) ** 2, 0.03, 'second weak-axis mode');
  // Tetrahedral mass matrices integrate to the element mass.
  const t = modalSolid({ model: boxMesh({ lx: 0.3, ly: 0.2, lz: 0.1, nx: 2, ny: 2, nz: 1, type: 'tet10' }), E, nu, rho, fix: [{ set: 'xmin', ux: 0, uy: 0, uz: 0 }], nModes: 1 });
  close(t.totalMass, rho * 0.3 * 0.2 * 0.1, 1e-12);
  assert.ok(t.frequencies[0] > 0);
});

test('principal stresses, selections and input validation', () => {
  const p = principalStresses({ sxx: 100, syy: -40, szz: 20, sxy: 30, syz: -10, szx: 15 });
  // Invariants: trace and von Mises from principals
  close(p.s1 + p.s2 + p.s3, 80, 1e-12);
  close(Math.sqrt(0.5 * ((p.s1 - p.s2) ** 2 + (p.s2 - p.s3) ** 2 + (p.s3 - p.s1) ** 2)), p.vonMises, 1e-12);
  const pure = principalStresses({ sxx: 0, syy: 0, szz: 0, sxy: 50 });
  close(pure.s1, 50, 1e-12); close(pure.s3, -50, 1e-12); close(pure.tresca, 100, 1e-12);
  const model = boxMesh({ lx: 1, ly: 1, lz: 1, nx: 2, ny: 2, nz: 2, type: 'tet4' });
  close(faceSetArea(model, boundaryFaces(model)).area, 6, 1e-12);
  const top = selectFaces(model, (f) => f.normal[2] > 0.9);
  close(faceSetArea(model, top).area, 1, 1e-12);
  assert.equal(meshStats(model).cellCounts.tet4, 48);
  assert.throws(() => solveSolid({ model, E: 1, nu: 0.3, fix: [] }), /no supports/);
  assert.throws(() => solveSolid({ model, E: 1, nu: 0.6, fix: [] }), RangeError);
  assert.throws(() => solveSolid({ model, E: 1, nu: 0.3, fix: [{ set: 'nope', ux: 0 }] }), RangeError);
  assert.throws(() => solveSolid({ model, E: 1, nu: 0.3, fix: [{ set: 'xmin', ux: 0, uy: 0, uz: 0 }], solver: 'ebe' }), RangeError);
  assert.throws(() => solveSolid({ model: { positions: model.positions, cells: [{ type: 'wedge6', conn: new Uint32Array(6) }] }, E: 1, nu: 0.3, fix: [] }), RangeError);
});
