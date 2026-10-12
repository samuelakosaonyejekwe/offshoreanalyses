import test from 'node:test';
import assert from 'node:assert/strict';
import { addedMassBem, wettedPanels } from '../../src/engine/hydro/bem.js';
import { sphere, tube, box, barge } from '../../src/mesh/generate.js';

const rho = 1000;
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
/** Double-precision copy of a generated surface moved by dz (keeps deep bodies free of float32 round-off). */
const sink = (m, dz) => { const p = Float64Array.from(m.positions); for (let i = 2; i < p.length; i += 3) p[i] += dz; return { ...m, positions: p }; };

test('wettedPanels: floating box is clipped at z = 0, lid dropped, normals into the fluid', () => {
  // box 4 × 2 × 3 centred at z = -0.5: draught 2, freeboard 1
  const b = box({ size: [4, 2, 3], center: [0, 0, -0.5], divisions: [4, 2, 3] });
  const w = wettedPanels(b);
  // wetted area = bottom + four sides to the waterline; displaced volume = L B T
  assert.ok(rel(w.area, 4 * 2 + 2 * (4 + 2) * 2) < 1e-6);
  assert.ok(rel(w.volume, 4 * 2 * 2) < 1e-6);
  assert.equal(w.flipped, false);
  let maxZ = -Infinity;
  for (let i = 2; i < w.vertices.length; i += 3) maxZ = Math.max(maxZ, w.vertices[i]);
  assert.equal(maxZ, 0);
  for (let p = 0; p < w.n; p++) {
    // a convex body centred on the z axis: outward normals point away from the centre line / keel centre
    const d = w.centroids[3 * p] * w.normals[3 * p] + w.centroids[3 * p + 1] * w.normals[3 * p + 1] + (w.centroids[3 * p + 2] + 1) * w.normals[3 * p + 2];
    assert.ok(d > 0);
    assert.ok(w.normals[3 * p + 2] <= 1e-12, 'no wetted panel faces upwards on a wall-sided box');
  }
  // waterline between mesh rows: the row of side triangles is cut
  const c = wettedPanels(box({ size: [4, 2, 3], center: [0, 0, -0.25], divisions: [4, 2, 3] }));
  assert.ok(c.clipped > 0);
  assert.ok(rel(c.volume, 4 * 2 * 1.75) < 1e-6);
  // a box floating flush (top face in the waterplane): the lid is not a wetted panel
  const flush = wettedPanels(box({ size: [2, 2, 1], center: [0, 0, -0.5] }));
  assert.ok(rel(flush.area, 4 + 4 * 2) < 1e-6);
});

test('wettedPanels: inward-oriented input is detected and reversed; invalid input throws', () => {
  const b = box({ size: [2, 2, 2], center: [0, 0, 0] });
  const idx = Uint32Array.from(b.indices);
  for (let f = 0; f < idx.length; f += 3) { const t = idx[f + 1]; idx[f + 1] = idx[f + 2]; idx[f + 2] = t; }
  const w = wettedPanels({ ...b, indices: idx });
  assert.equal(w.flipped, true);
  assert.ok(rel(w.volume, 4) < 1e-6);
  assert.throws(() => wettedPanels(box({ size: [1, 1, 1], center: [0, 0, 5] })), /no wetted panels/);
  assert.throws(() => wettedPanels({}), RangeError);
});

test('deeply submerged sphere: A11 = A22 = A33 = rho V / 2 (Lamb §92), second-order convergence', () => {
  // Faceted icospheres with 80, 320 and 1280 panels (edge halved each time), centre 60 radii below the surface so
  // the image contributes O((a/2d)³) ~ 1e-6. Reference: rho (2/3) pi a³ for the true sphere.
  const ref = rho * (2 / 3) * Math.PI;
  const err = [];
  for (const s of [1, 2, 3]) {
    const r = addedMassBem({ surface: sink(sphere({ radius: 1, subdivisions: s }), -60), rho, limit: 'infinite', cog: [0, 0, -60] });
    assert.equal(r.panels, 20 * 4 ** s);
    // the three translations agree to the icosahedral symmetry of the mesh and the tiny image effect
    assert.ok(rel(r.A[1][1], r.A[0][0]) < 1e-5 && rel(r.A[2][2], r.A[0][0]) < 1e-5);
    // no added inertia in rotation about the centre, no coupling
    assert.ok(Math.abs(r.A[3][3]) < 1e-3 * ref && Math.abs(r.A[4][4]) < 1e-3 * ref && Math.abs(r.A[5][5]) < 1e-3 * ref);
    assert.ok(Math.abs(r.A[0][1]) < 1e-9 * ref);
    err.push(rel(r.A[0][0], ref));
  }
  // observed order from successive halvings of the panel size: p = log2(e_h / e_{h/2})
  const p1 = Math.log2(err[0] / err[1]); const p2 = Math.log2(err[1] / err[2]);
  // observed on this mesh family: errors 7.8 %, 2.0 %, 0.51 %  ->  p = 1.95, 2.00
  assert.ok(err[2] < 0.006, `error at 1280 panels ${err[2]}`);
  assert.ok(p1 > 1.7 && p1 < 2.3 && p2 > 1.8 && p2 < 2.2, `observed order ${p1}, ${p2}`);
});

test('deep sphere: both limits and the unbounded fluid coincide; source formulation converges at first order', () => {
  const m = sink(sphere({ radius: 1, subdivisions: 2 }), -60);
  const cog = [0, 0, -60];
  const a = addedMassBem({ surface: m, rho, limit: 'infinite', cog });
  const b = addedMassBem({ surface: m, rho, limit: 'zero', cog });
  const c = addedMassBem({ surface: m, rho, limit: 'unbounded', cog });
  assert.ok(rel(a.A[2][2], c.A[2][2]) < 1e-5 && rel(b.A[2][2], c.A[2][2]) < 1e-5);
  // opposite-sign image lowers, equal-sign image raises the heave added mass of a submerged body
  assert.ok(a.A[2][2] < c.A[2][2] && c.A[2][2] < b.A[2][2]);
  // Hess–Smith source panels on the same meshes, measured against rho V_facets / 2 so that only the
  // integral-equation error remains: 15.8 %, 7.4 %, 3.6 %  ->  p = 1.1, 1.05
  const err = [];
  for (const s of [1, 2, 3]) {
    const r = addedMassBem({ surface: sink(sphere({ radius: 1, subdivisions: s }), -60), rho, cog, formulation: 'source' });
    err.push(rel(r.A[0][0], 0.5 * rho * r.volume));
  }
  const p1 = Math.log2(err[0] / err[1]); const p2 = Math.log2(err[1] / err[2]);
  assert.ok(p1 > 0.8 && p1 < 1.4 && p2 > 0.8 && p2 < 1.3, `observed order ${p1}, ${p2}`);
});

test('floating hemisphere: A33(inf) = A11(0) = rho pi a³ / 3 exactly; A33(0) from the Legendre series (0.8310)', () => {
  // omega -> inf: phi = 0 on z = 0, so phi is odd in z. The heaving hemisphere and its image then form a complete
  // sphere with d phi/dr = cos(theta) everywhere, i.e. a rigid sphere translating along z in unbounded fluid,
  // whose added mass rho (2/3) pi a³ is shared equally by the two halves:  A33(inf) = rho pi a³ / 3.
  // omega -> 0: d phi/dz = 0 on z = 0, phi even in z. Surge/sway of the double body is the same rigid sphere
  // translating horizontally:  A11(0) = A22(0) = rho pi a³ / 3.
  // Heave at omega -> 0: d phi/dr = -|mu| on the full sphere (mu = cos theta). Expanding in even Legendre
  // polynomials, -|mu| = sum c_n P_n with c_n = -(2n + 1) I_n, I_n = ∫₀¹ mu P_n dmu, the exterior potential is
  // phi = sum a (2n + 1) I_n / (n + 1) (a/r)^(n+1) P_n and
  //   A33(0) = -rho ∫ phi n_z dS = 2 pi rho a³ sum_{n even} (2n + 1) I_n² / (n + 1) = 0.8310 · rho (2/3) pi a³,
  // the value tabulated by Hulme, J. Fluid Mech. 121 (1982), for ka -> 0.
  const I = (n) => {
    const N = 4000; let s = 0; // Simpson on [0, 1]
    for (let i = 0; i <= N; i++) {
      const mu = i / N; let p0 = 1; let p1 = mu;
      for (let k = 2; k <= n; k++) { const pk = ((2 * k - 1) * mu * p1 - (k - 1) * p0) / k; p0 = p1; p1 = pk; }
      s += (i === 0 || i === N ? 1 : i % 2 ? 4 : 2) * mu * (n === 0 ? 1 : p1);
    }
    return s / (3 * N);
  };
  let series = 0;
  for (let n = 0; n <= 80; n += 2) series += (2 * n + 1) * I(n) ** 2 / (n + 1);
  const heave0 = 2 * Math.PI * rho * series;
  assert.ok(Math.abs(heave0 / (rho * (2 / 3) * Math.PI) - 0.8310) < 5e-4, 'series reproduces Hulme (1982)');

  const m = sphere({ radius: 1, subdivisions: 3 }); // 656 wetted panels after clipping at the waterplane
  const ref = rho * Math.PI / 3;
  const inf = addedMassBem({ surface: m, rho, limit: 'infinite' });
  const zero = addedMassBem({ surface: m, rho, limit: 'zero' });
  assert.ok(rel(inf.volume, (2 / 3) * Math.PI) < 0.01);
  assert.ok(rel(inf.A[2][2], ref) < 0.008, `A33(inf) ${inf.A[2][2] / ref}`);
  assert.ok(rel(zero.A[0][0], ref) < 0.008 && rel(zero.A[1][1], ref) < 0.008, `A11(0) ${zero.A[0][0] / ref}`);
  assert.ok(rel(zero.A[2][2], heave0) < 0.012, `A33(0) ${zero.A[2][2] / heave0}`);
  // the two limits differ: heave is larger under a rigid lid, surge smaller under phi = 0
  assert.ok(zero.A[2][2] > 1.5 * inf.A[2][2] && inf.A[0][0] < 0.7 * zero.A[0][0]);
});

test('long horizontal circular cylinder, deeply submerged: sway and heave approach rho pi r² L (Lamb §68)', () => {
  // strip value for an infinitely long cylinder; a finite cylinder lies below it because of the flow round the
  // ends, and approaches it as L/D grows (compare Lamb's k2 = 0.96 for a 10:1 spheroid).
  const run = (L, nd) => addedMassBem({ surface: sink(tube({ a: [-L / 2, 0, 0], b: [L / 2, 0, 0], radius: 1, segments: 20, heightDivisions: nd, capRings: 1 }), -80), rho, limit: 'zero', cog: [0, 0, -80] });
  const long = run(40, 30); const short = run(10, 12);
  const kLong = long.A[1][1] / (rho * Math.PI * 40); const kShort = short.A[1][1] / (rho * Math.PI * 10);
  assert.ok(rel(long.A[2][2], long.A[1][1]) < 1e-3, 'sway = heave for a body of revolution');
  assert.ok(kLong > 0.94 && kLong < 1.0, `L/D = 20: ${kLong}`);
  assert.ok(kShort < kLong && kShort > 0.8, `L/D = 5: ${kShort}`);
  // axial added mass of a slender body is small
  assert.ok(long.A[0][0] < 0.05 * long.A[1][1]);
});

test('added-mass matrix: symmetric to round-off, raw asymmetry reported, reference-point transfer, solvers agree', () => {
  const hull = barge({ length: 20, beam: 6, depth: 4, draft: 2, rakeBow: 3, elementSize: 1 });
  const r = addedMassBem({ surface: hull, rho, cog: [1, 0.5, -1], solver: 'lu' });
  for (let i = 0; i < 6; i++) {
    assert.ok(r.A[i][i] > 0);
    for (let j = 0; j < 6; j++) assert.equal(r.A[i][j], r.A[j][i]);
  }
  // positive definite (kinetic energy of the fluid): Cholesky succeeds
  const L = r.A.map((row) => row.slice());
  for (let i = 0; i < 6; i++) for (let j = 0; j <= i; j++) {
    let s = L[i][j];
    for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
    if (i === j) { assert.ok(s > 0, 'positive definite'); L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j];
  }
  assert.ok(r.symmetryError > 0 && r.symmetryError < 0.02, `raw asymmetry ${r.symmetryError}`);
  assert.ok(r.conditionNumber > 1 && r.conditionNumber < 100);
  assert.equal(r.solver, 'lu');
  const g = addedMassBem({ surface: hull, rho, cog: [1, 0.5, -1] });
  assert.equal(g.solver, 'gmres');
  assert.equal(g.conditionNumber, null);
  for (let i = 0; i < 6; i++) assert.ok(rel(g.A[i][i], r.A[i][i]) < 1e-8);
  // moving the reference point by d = (0, 0, 2) leaves the translations unchanged and transfers the moments:
  // with a sphere (no rotational added mass about its centre) A55 = A11 d², A15 = -A11 d, A24 = +A11 d
  const sp = sink(sphere({ radius: 1, subdivisions: 2 }), -60);
  const c = addedMassBem({ surface: sp, rho, cog: [0, 0, -58] });
  assert.ok(rel(c.A[4][4], 4 * c.A[0][0]) < 1e-3);
  assert.ok(rel(c.A[0][4], -2 * c.A[0][0]) < 1e-9 && rel(c.A[1][3], 2 * c.A[1][1]) < 1e-9);
});

test('finite depth (infinite-frequency limit): deep limit recovered, shallow water raises heave, formulations agree', () => {
  const hull = barge({ length: 20, beam: 6, depth: 4, draft: 2, elementSize: 1 });
  const deep = addedMassBem({ surface: hull, rho });
  const far = addedMassBem({ surface: hull, rho, depth: 2000 });
  assert.ok(rel(far.A[2][2], deep.A[2][2]) < 1e-5);
  const shallow = addedMassBem({ surface: hull, rho, depth: 3 });
  const mid = addedMassBem({ surface: hull, rho, depth: 6 });
  // 1 m under-keel clearance: the fluid is squeezed out sideways
  assert.ok(shallow.A[2][2] > mid.A[2][2] && mid.A[2][2] > deep.A[2][2] && shallow.A[2][2] > 1.5 * deep.A[2][2]);
  // truncation of the image series: doubling the retained pairs changes the result by < 0.1 %
  assert.ok(rel(addedMassBem({ surface: hull, rho, depth: 3, imagePairs: 24 }).A[2][2], shallow.A[2][2]) < 1e-3);
  // independent discretisation (source panels) gives the same amplification within mesh accuracy
  const s0 = addedMassBem({ surface: hull, rho, formulation: 'source' });
  const s3 = addedMassBem({ surface: hull, rho, depth: 3, formulation: 'source' });
  assert.ok(rel(s3.A[2][2] / s0.A[2][2], shallow.A[2][2] / deep.A[2][2]) < 0.05);
  assert.ok(rel(s0.A[2][2], deep.A[2][2]) < 0.05 && rel(s0.A[0][0], deep.A[0][0]) < 0.08);
});

test('addedMassBem: invalid arguments throw RangeError naming the argument', () => {
  const m = sphere({ radius: 1, subdivisions: 1 });
  assert.throws(() => addedMassBem({ surface: m, rho: -1 }), /rho/);
  assert.throws(() => addedMassBem({ surface: m, limit: 'high' }), /limit/);
  assert.throws(() => addedMassBem({ surface: m, cog: [0, 0] }), /cog/);
  assert.throws(() => addedMassBem({ surface: m, depth: 0.5 }), /depth must exceed the draught/);
  assert.throws(() => addedMassBem({ surface: m, limit: 'zero', depth: 50 }), /depth/);
  assert.throws(() => addedMassBem({ surface: m, formulation: 'dipole' }), /formulation/);
  assert.throws(() => addedMassBem({ surface: m, solver: 'qr' }), /solver/);
});
