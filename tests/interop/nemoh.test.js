import test from 'node:test';
import assert from 'node:assert/strict';
import { writeNemohCase, readNemohResults, writeWamitGdf, writeCapytaineScript, readCapytaineCsv, wettedMesh, parseNemohRadiation, parseNemohExcitation, parseNemohMatrixDat } from '../../src/interop/nemoh.js';
import { parseNemoh, parseGdf } from '../../src/geometry/importers/panels.js';
import { raoSolve, checkHydroDb } from '../../src/engine/hydro/hydrodb.js';
import { box, cylinder } from '../../src/mesh/generate.js';

const file = (r, path) => { const f = r.files.find((x) => x.path === path); assert.ok(f, `file ${path}`); return f.text; };
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

/** Unit normals and centroids of a triangle surface. */
function facets(m) {
  const out = [];
  for (let f = 0; f < m.indices.length; f += 3) {
    const [a, b, c] = [0, 1, 2].map((k) => [0, 1, 2].map((q) => m.positions[3 * m.indices[f + k] + q]));
    const u = b.map((v, k) => v - a[k]); const w = c.map((v, k) => v - a[k]);
    const n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
    const l = Math.hypot(...n);
    out.push({ n: n.map((v) => v / l), c: a.map((v, k) => (v + b[k] + c[k]) / 3), area: l / 2 });
  }
  return out;
}

// box 4 × 2 × 3 m, draught 2 m, mesh rows on the waterline: 4·2 bottom cells and 2·(4 + 2)·2 side cells, two triangles each
const BOX = box({ size: [4, 2, 3], center: [0, 0, -0.5], divisions: [4, 2, 3], name: 'barge box' });
const BOX_PANELS = 2 * (4 * 2 + 2 * (4 + 2) * 2);
// cylinder R = 1, draught 3 (z = -3 … +1): side 12 segments × 3 wetted rows × 2 triangles, bottom cap of 2 rings: 12·(2·2 − 1)
const CYL = cylinder({ radius: 1, height: 4, base: [0, 0, -3], segments: 12, heightDivisions: 4, capRings: 2, name: 'spar' });
const CYL_PANELS = 12 * 3 * 2 + 12 * 3;

test('NEMOH mesh of a box: panel and node counts, normals into the fluid, read back by the panel importer', () => {
  const r = writeNemohCase({ surface: BOX, omegas: [0.2, 0.4, 0.6, 0.8, 1.0], headingsDeg: [0, 45, 90], cog: [0, 0, -0.5], depth: 30, name: 'barge box' });
  assert.equal(r.mesh.nPanels, BOX_PANELS);
  // wetted nodes: lattice points on the bottom and sides with z <= 0
  assert.equal(r.mesh.nNodes, 5 * 3 + 2 * (4 + 2) * 2);
  assert.ok(rel(r.mesh.volume, 16) < 1e-6 && rel(r.mesh.area, 8 + 24) < 1e-6);
  assert.deepEqual(r.files.map((f) => f.path), ['Nemoh.cal', 'input_solver.txt', 'input.txt', 'ID.dat', 'mesh/barge_box.dat', 'results/']);
  const text = file(r, 'mesh/barge_box.dat');
  const lines = text.trim().split('\n');
  assert.equal(lines[0].trim().split(/\s+/).join(' '), '2 0');
  assert.equal(lines.length, 1 + r.mesh.nNodes + 1 + BOX_PANELS + 1);
  assert.deepEqual(lines[lines.length - 1].trim().split(/\s+/), ['0', '0', '0', '0']);
  // triangles are quadrilaterals with a repeated node
  const quad = lines[1 + r.mesh.nNodes + 1].trim().split(/\s+/);
  assert.equal(quad.length, 4);
  assert.equal(quad[3], quad[0]);

  const back = parseNemoh({ name: 'barge_box.dat', text }).models[0];
  assert.equal(back.indices.length / 3, BOX_PANELS);
  const fs = facets(back);
  let area = 0; let vol = 0;
  for (const f of fs) {
    area += f.area; vol += f.c[2] * f.n[2] * f.area;
    assert.ok(f.c[2] < 0, 'only wetted panels');
    // convex box: the normal points away from the keel centre, i.e. into the fluid
    assert.ok(f.n[0] * f.c[0] + f.n[1] * f.c[1] + f.n[2] * (f.c[2] + 1) > 0);
  }
  assert.ok(rel(area, 32) < 1e-6);
  assert.ok(rel(vol, 16) < 1e-6, 'divergence theorem with normals out of the hull');
});

test('NEMOH mesh and GDF of a vertical cylinder: counts, radial side normals, downward bottom normals', () => {
  const m = wettedMesh(CYL);
  assert.equal(m.nPanels, CYL_PANELS);
  // inscribed 12-gon: wetted volume = polygon area × draught
  assert.ok(rel(m.volume, 0.5 * 12 * Math.sin(Math.PI / 6) * 3) < 1e-5);
  let side = 0; let bottom = 0;
  for (let p = 0; p < m.nPanels; p++) {
    const n = m.normals.subarray(3 * p, 3 * p + 3); const c = m.centroids.subarray(3 * p, 3 * p + 3);
    if (Math.abs(n[2]) < 1e-6) { side++; const r = Math.hypot(c[0], c[1]); assert.ok((n[0] * c[0] + n[1] * c[1]) / r > 0.96, 'side normal points radially outwards'); } else { bottom++; assert.ok(n[2] < -0.999999 && Math.abs(c[2] + 3) < 1e-6); }
  }
  assert.equal(side, 72);
  assert.equal(bottom, 36);

  const g = writeWamitGdf({ surface: CYL, ulen: 2, name: 'spar' });
  assert.equal(g.files[0].path, 'spar.gdf');
  const lines = g.files[0].text.trim().split('\n');
  assert.deepEqual(lines[1].trim().split(/\s+/).slice(0, 2).map(Number), [2, 9.80665]);
  assert.deepEqual(lines[2].trim().split(/\s+/).slice(0, 2), ['0', '0']);
  assert.equal(Number(lines[3].trim().split(/\s+/)[0]), CYL_PANELS);
  assert.equal(lines.length, 4 + 4 * CYL_PANELS);
  const back = parseGdf({ name: 'spar.gdf', text: g.files[0].text });
  assert.equal(back.meta.npan, CYL_PANELS);
  assert.equal(back.meta.units.id ?? back.meta.units.name ?? 'm', 'm');
  const fs = facets(back.models[0]);
  assert.equal(fs.length, CYL_PANELS);
  let vol = 0;
  for (const f of fs) vol += f.c[2] * f.n[2] * f.area;
  assert.ok(rel(vol, m.volume) < 1e-6, 'GDF vertices are counter-clockwise seen from the fluid');

  // the same counts for the box through the GDF writer
  assert.equal(writeWamitGdf({ surface: BOX }).mesh.nPanels, BOX_PANELS);
});

test('Nemoh.cal: one body, six modes and six forces about the centre, frequency and heading ranges', () => {
  const r = writeNemohCase({ surface: BOX, rho: 1000, g: 9.81, depth: 30, omegas: [0.2, 0.4, 0.6, 0.8, 1.0], headingsDeg: [0, 45, 90], cog: [0.5, 0, -0.75], name: 'b' });
  const L = file(r, 'Nemoh.cal').split('\n').map((l) => l.split('!')[0].trim().split(/\s+/));
  assert.deepEqual([L[1][0], L[2][0], L[3][0]].map(Number), [1000, 9.81, 30]);
  assert.equal(L[6][0], '1');
  assert.equal(L[8][0], 'b.dat');
  assert.deepEqual(L[9].map(Number), [r.mesh.nNodes, r.mesh.nPanels]);
  assert.equal(L[10][0], '6');
  assert.deepEqual(L[11].map(Number), [1, 1, 0, 0, 0, 0, 0]);
  assert.deepEqual(L[14].map(Number), [2, 1, 0, 0, 0.5, 0, -0.75]);
  assert.deepEqual(L[16].map(Number), [2, 0, 0, 1, 0.5, 0, -0.75]);
  assert.equal(L[17][0], '6');
  assert.deepEqual(L[23].map(Number), [2, 0, 0, 1, 0.5, 0, -0.75]);
  assert.equal(L[24][0], '0');
  assert.deepEqual(L[26].map(Number), [1, 5, 0.2, 1]);
  assert.deepEqual(L[27].map(Number), [3, 0, 90]);
  // deep water is written as depth 0; non-uniform lists are replaced by a uniform range and reported
  const d = writeNemohCase({ surface: BOX, omegas: [0.2, 0.3, 1.0], headingsDeg: [0] });
  assert.equal(file(d, 'Nemoh.cal').split('\n')[3].split('!')[0].trim(), '0.');
  assert.ok(d.notes.some((n) => /not uniformly spaced/.test(n)));
  assert.ok(!r.notes.some((n) => /not uniformly spaced/.test(n)));
  assert.throws(() => writeNemohCase({ surface: BOX, omegas: [] }), /omegas/);
  assert.throws(() => writeNemohCase({ surface: BOX, omegas: [0.5], depth: 1 }), /depth must exceed the draught/);
  assert.throws(() => writeNemohCase({ surface: BOX, omegas: [0.5], cog: [0, 0] }), /cog/);
});

/** Result files in the layout NEMOH writes (Tecplot POINT zones; CA/CM blocks), from known matrices. */
function nemohFixture({ omega, A, B, headingsDeg, F, phase = 'rad' }) {
  const e = (v) => ` ${v.toExponential(6).padStart(14)}`;
  const rad = ['VARIABLES="w (rad/s)"'];
  for (let i = 1; i <= 6; i++) rad.push(`"A ${String(1).padStart(4)}${String(i).padStart(4)}" "B ${String(1).padStart(4)}${String(i).padStart(4)}"`);
  for (let j = 0; j < 6; j++) {
    rad.push(`Zone t="Motion of body    1 in DoF    ${j + 1}",I=${String(omega.length).padStart(5)},F=POINT`);
    omega.forEach((w, k) => { let l = e(w); for (let i = 0; i < 6; i++) l += e(A[k][i][j]) + e(B[k][i][j]); rad.push(l); });
  }
  const exc = ['VARIABLES="w (rad/s)"'];
  for (let i = 1; i <= 6; i++) exc.push(`"abs(F ${String(1).padStart(4)}${String(i).padStart(4)})" "angle(F ${String(1).padStart(4)}${String(i).padStart(4)})"`);
  headingsDeg.forEach((h, b) => {
    exc.push(`Zone t="Diffraction force - beta = ${h.toFixed(3).padStart(8)} deg",I=${String(omega.length).padStart(5)},F=POINT`);
    omega.forEach((w, k) => {
      const parts = [e(w)];
      // wrap the record over two lines, as list-directed Fortran output may do
      for (let i = 0; i < 6; i++) { const [re, im] = F[b][k][i]; parts.push(e(Math.hypot(re, im)) + e(Math.atan2(im, re) * (phase === 'deg' ? 180 / Math.PI : 1))); }
      exc.push(parts.slice(0, 4).join(''), parts.slice(4).join(''));
    });
  });
  const block = (M, label) => [` Nb of ${label} :   ${omega.length}`, ...omega.flatMap((w, k) => [e(label === 'periods' ? (2 * Math.PI) / w : w), ...M[k].map((row) => row.map(e).join(''))])].join('\n');
  return { rad: rad.join('\n'), exc: exc.join('\n'), cm: block(A, 'periods'), ca: block(B, 'periods'), cmW: block(A, 'frequencies') };
}

test('readNemohResults: radiation and excitation files → hydro database (conjugated to exp(+iωt)), usable by raoSolve', () => {
  const omega = [0.4, 0.8, 1.2];
  const headingsDeg = [0, 90, 180, 270];
  const A = omega.map((w, k) => Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => (i === j ? 1000 * (i + 1) * (1 + 0.1 * k) : 3 * (i + 1) + 0.5 * (j + 1) + k))));
  const B = omega.map((w, k) => Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => (i === j ? 50 * (i + 1) * w : 0.25 * (i - j) * (k + 1)))));
  const F = headingsDeg.map((h, b) => omega.map((w, k) => Array.from({ length: 6 }, (_, i) => [100 * (i + 1) * Math.cos(0.3 * i + b), -40 * (k + 1) * Math.sin(0.7 * i + b + 0.2)])));
  const fx = nemohFixture({ omega, A, B, headingsDeg, F });
  const cal = writeNemohCase({ surface: BOX, rho: 1000, g: 9.81, depth: 30, omegas: omega, headingsDeg: [0, 90, 180, 270] }).files[0];
  const db = readNemohResults({ files: [cal, { path: 'results/RadiationCoefficients.tec', text: fx.rad }, { path: 'results\\ExcitationForce.tec', text: fx.exc }] });
  assert.deepEqual([db.rho, db.g, db.depth], [1000, 9.81, 30]);
  assert.deepEqual(db.omega, omega);
  assert.deepEqual(db.headingsDeg, headingsDeg);
  assert.equal(db.periodicHeading, true);
  for (let k = 0; k < 3; k++) for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
    // zone j holds the force in mode i due to motion in mode j: no transposition
    assert.ok(Math.abs(db.A[k][i][j] - A[k][i][j]) <= 1e-6 * Math.abs(A[k][i][j]) + 1e-9);
    assert.ok(Math.abs(db.B[k][i][j] - B[k][i][j]) <= 1e-6 * Math.abs(B[k][i][j]) + 1e-9);
  }
  for (let b = 0; b < 4; b++) for (let k = 0; k < 3; k++) for (let i = 0; i < 6; i++) {
    assert.ok(Math.abs(db.Fre[b][k][i] - F[b][k][i][0]) < 1e-3);
    assert.ok(Math.abs(db.Fim[b][k][i] + F[b][k][i][1]) < 1e-3, 'imaginary part conjugated');
  }
  assert.match(db.source, /RadiationCoefficients\.tec, ExcitationForce\.tec/);
  // phase in degrees is detected; conjugation can be switched off; the parsers are exported on their own
  const deg = parseNemohExcitation(nemohFixture({ omega, A, B, headingsDeg, F, phase: 'deg' }).exc, { conjugate: false });
  assert.equal(deg.phaseUnit, 'deg');
  assert.ok(Math.abs(deg.Fim[1][2][3] - F[1][2][3][1]) < 1e-3);
  assert.equal(parseNemohRadiation(fx.rad).A[2][4][1], db.A[2][4][1]);

  // CM.dat / CA.dat route (per-block abscissa as period or as frequency), object-map input, explicit environment
  const db2 = readNemohResults({ files: { 'CM.dat': fx.cm, 'CA.dat': fx.ca }, rho: 1025, depth: null });
  assert.ok(db2.omega.every((w, k) => rel(w, omega[k]) < 1e-6));
  assert.ok(rel(db2.A[1][2][2], A[1][2][2]) < 1e-6 && rel(db2.B[2][5][5], B[2][5][5]) < 1e-6);
  assert.equal(db2.headingsDeg.length, 0);
  assert.ok(rel(parseNemohMatrixDat(fx.cmW).omega[2], 1.2) < 1e-6);
  assert.throws(() => readNemohResults({ files: { 'CM.dat': fx.cm } }), /CA\.dat/);
  assert.throws(() => readNemohResults({ files: [] }), /files/);

  // the database drives the frequency-domain motion solver
  const M = Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => (i === j ? 5000 : 0)));
  const C = Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => (i === j ? 2e4 : 0)));
  const rao = raoSolve({ db, M, C, headingDeg: 45 });
  assert.equal(rao.amp.length, 6);
  assert.ok(rao.amp.every((a) => a.every((v) => Number.isFinite(v) && v >= 0)));
  assert.equal(checkHydroDb(db).negativeDamping.length, 0);
  // an infinite-frequency added mass (e.g. from the panel solver) is carried along for the time-domain model
  const withInf = readNemohResults({ files: { 'RadiationCoefficients.tec': fx.rad }, Ainf: A[2] });
  assert.equal(withInf.Ainf[3][3], A[2][3][3]);
});

test('Capytaine: script names the mesh and the cases; the CSV it writes is read back into the database', () => {
  const r = writeCapytaineScript({ surface: CYL, rho: 1000, depth: 50, omegas: [0.5, 1.0], headingsDeg: [0, 30], cog: [0, 0, -1.5], name: 'spar' });
  assert.deepEqual(r.files.map((f) => f.path), ['run_capytaine.py', 'mesh/spar.dat']);
  const py = r.files[0].text;
  assert.match(py, /cpt\.load_mesh\("mesh\/spar\.dat", file_format="nemoh"\)/);
  assert.match(py, /OMEGAS = \[0\.5, 1\]/);
  assert.match(py, /HEADINGS_DEG = \[0, 30\]/);
  assert.match(py, /DEPTH = 50\n/);
  assert.match(py, /rotation_center=CENTRE/);
  assert.match(py, /CENTRE = \(0, 0, -1\.5\)/);
  assert.match(writeCapytaineScript({ surface: CYL, omegas: [1] }).files[0].text, /DEPTH = np\.inf/);
  assert.equal(parseNemoh({ name: 'spar.dat', text: r.files[1].text }).models[0].indices.length / 3, CYL_PANELS);
  assert.throws(() => writeCapytaineScript({ surface: CYL, omegas: [0, 1] }), /omegas/);

  // rows exactly as the script's format strings produce them (Python repr of floats)
  const rows = ['# solver=capytaine 2.2', '# rho=1000.0', '# g=9.80665', '# depth=50.0', 'kind,omega,headingDeg,i,j,re,im'];
  const a = (w, i, j) => 100 * i + j + w; const b = (w, i, j) => (i === j ? 10 * w * i : 0);
  for (const w of [0.5, 1.0]) {
    for (let j = 1; j <= 6; j++) for (let i = 1; i <= 6; i++) { rows.push(`A,${w},,${i},${j},${a(w, i, j)},`); rows.push(`B,${w},,${i},${j},${b(w, i, j)},`); }
    for (const h of [0, 30]) for (let i = 1; i <= 6; i++) rows.push(`F,${w},${h},${i},,${1e3 * i + h},${-2.5 * i * w}`);
  }
  for (let j = 1; j <= 6; j++) for (let i = 1; i <= 6; i++) rows.push(`AINF,,,${i},${j},${i === j ? 77 * i : 0},`);
  const db = readCapytaineCsv(rows.join('\n') + '\n');
  assert.deepEqual([db.rho, db.g, db.depth], [1000, 9.80665, 50]);
  assert.deepEqual(db.omega, [0.5, 1]);
  assert.deepEqual(db.headingsDeg, [0, 30]);
  assert.equal(db.A[1][2][4], a(1, 3, 5));
  assert.equal(db.B[0][3][3], b(0.5, 4, 4));
  assert.deepEqual([db.Fre[1][0][2], db.Fim[1][1][5]], [3030, -15]);
  assert.equal(db.Ainf[4][4], 385);
  assert.equal(db.source, 'capytaine 2.2');
  assert.equal(readCapytaineCsv(rows.join('\n').replace('# depth=50.0', '# depth=inf')).depth, null);
});
