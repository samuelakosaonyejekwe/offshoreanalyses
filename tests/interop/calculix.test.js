import test from 'node:test';
import assert from 'node:assert/strict';
import { writeCalculixCase, readCalculixFrd, readCalculixDat } from '../../src/interop/calculix.js';
import { parseAbaqus } from '../../src/geometry/importers/abaqus.js';

const byKind = (r, k) => r.models.find((m) => m.kind === k);
const STEEL = { name: 'S355', E: 2.1e11, nu: 0.3, rho: 7850, fy: 355e6 };

/** Cantilever along x: 5 nodes, 4 tubular members (D = 0.5 m, t = 20 mm), 8 m long. */
function cantilever() {
  const nodes = new Float64Array(15);
  for (let i = 0; i < 5; i++) nodes[3 * i] = 2 * i;
  const members = [];
  for (let i = 0; i < 4; i++) members.push({ a: i, b: i + 1, section: { type: 'tube', D: 0.5, t: 0.02 }, material: 'S355' });
  return { kind: 'frame', name: 'cantilever', nodes, members };
}

test('cantilever deck is re-read by the Abaqus importer to the same nodes, elements, section, material and sets', () => {
  const model = cantilever();
  const r = writeCalculixCase(model, { material: STEEL, supports: [{ nodes: [0], dofs: 'fixed' }], loads: [{ type: 'force', nodes: [4], F: [0, 0, -1e4] }], analysis: 'static' });
  assert.equal(r.files.length, 1);
  assert.equal(r.files[0].path, 'cantilever.inp');
  assert.deepEqual(r.counts, { nodes: 5, elements: 4, appNodes: 5 });
  const deck = r.files[0].text;
  const back = parseAbaqus({ name: 'cantilever.inp', text: deck });
  const f = byKind(back, 'frame');
  assert.equal(back.models.length, 1);
  assert.deepEqual(Array.from(f.nodes), Array.from(model.nodes));
  assert.deepEqual(f.members.map((m) => [m.a, m.b]), model.members.map((m) => [m.a, m.b]));
  for (const m of f.members) { assert.deepEqual(m.section, { type: 'tube', D: 0.5, t: 0.02 }); assert.equal(m.material, 'S355'); assert.deepEqual(m.orient, [0, 0, -1]); }
  assert.deepEqual(back.meta.elementTypes, { line2: 4 });
  assert.equal(back.meta.materials.S355.E, 2.1e11);
  assert.equal(back.meta.materials.S355.nu, 0.3);
  assert.equal(back.meta.materials.S355.rho, 7850);
  assert.deepEqual(Array.from(f.sets.nodes.SUP1), [0]);
  assert.deepEqual(Array.from(f.sets.nodes.NALL), [0, 1, 2, 3, 4]);
  assert.deepEqual(Array.from(f.sets.members.EALL), [0, 1, 2, 3]);
  // step content: fixed root (six degrees of freedom), tip load, linear static, outputs on the original nodes
  assert.match(deck, /\*BOUNDARY\nSUP1, 1, 1\nSUP1, 2, 2\nSUP1, 3, 3\nSUP1, 4, 4\nSUP1, 5, 5\nSUP1, 6, 6\n/);
  assert.match(deck, /\*STEP, INC=200\n\*STATIC\n\*CLOAD\n5, 3, -10000\n/);
  assert.match(deck, /\*BEAM SECTION, ELSET=BEAM1, MATERIAL=S355, SECTION=PIPE\n0\.25, 0\.02\n0, 0, -1\n/);
  assert.match(deck, /\*PLASTIC\n355000000, 0\.\n/);
  assert.match(deck, /\*NODE FILE, OUTPUT=2D\nU, RF\n\*EL FILE, OUTPUT=2D\nS, E, PEEQ\n/);
  assert.match(deck, /\*NODE PRINT, NSET=NSUP, TOTALS=YES\nRF\n/);
  assert.ok(deck.trimEnd().endsWith('*END STEP'));
});

test('analysis types: NLGEOM, buckling, frequency, dynamic; hardening; B32; distributed loads; vertical members', () => {
  const model = cantilever();
  const base = { material: { ...STEEL, hardening: { Et: 2.1e9 } }, supports: [{ nodes: [0] }], loads: [{ type: 'force', nodes: [4], F: [-1e5, 0, 0] }, { type: 'gravity' }, { type: 'line', q: [0, 0, -500] }] };
  const nl = writeCalculixCase(model, { ...base, analysis: 'static-nlgeom', steps: { increments: 50, initialIncrement: 0.05 } }).files[0].text;
  assert.match(nl, /\*STEP, NLGEOM, INC=50\n\*STATIC\n0\.05, 1\., 0\.000001, 0\.25\n/);
  // bilinear hardening: H = E Et / (E − Et) = 2.1e11·2.1e9 / 2.079e11 = 2.1212e9 Pa; at 20 % plastic strain 355e6 + 0.2 H
  const pl = /\*PLASTIC\n355000000, 0\.\n([0-9.e+]+), 0\.2\n/.exec(nl);
  assert.ok(pl && Math.abs(Number(pl[1]) - (355e6 + 0.2 * 2.1e11 * 2.1e9 / (2.1e11 - 2.1e9))) < 1);
  assert.match(nl, /\*DLOAD\nEALL, GRAV, 9\.80665, 0, 0, -1\n/);
  // line load lumped: q L / 2 = 500 N on end nodes of each 2 m member, 1000 N on interior nodes; tip also carries the point load in x
  assert.match(nl, /\n1, 3, -500\n/);
  assert.match(nl, /\n3, 3, -1000\n/);
  assert.match(nl, /\n5, 1, -100000\n/);

  const bk = writeCalculixCase(model, { ...base, analysis: 'buckle', steps: 3 });
  assert.match(bk.files[0].text, /\*STEP\n\*BUCKLE\n3\n/);
  assert.ok(!/\*PLASTIC/.test(bk.files[0].text) && bk.notes.some((n) => /yield stress is not written/.test(n)));
  const fr = writeCalculixCase(model, { material: STEEL, supports: [{ nodes: [0] }], analysis: 'frequency', steps: { nModes: 6 } }).files[0].text;
  assert.match(fr, /\*STEP\n\*FREQUENCY\n6\n\*NODE FILE, OUTPUT=2D\nU\n/);
  assert.ok(!/\*CLOAD/.test(fr));
  const dy = writeCalculixCase(model, { ...base, analysis: 'dynamic', steps: { timeStep: 0.01, endTime: 2, nlgeom: true } }).files[0].text;
  assert.match(dy, /\*STEP, NLGEOM, INC=210\n\*DYNAMIC, DIRECT\n0\.01, 2\n/);
  const tab = writeCalculixCase(model, { ...base, material: { ...STEEL, hardening: { table: [[355e6, 0], [470e6, 0.05], [510e6, 0.15]] } } }).files[0].text;
  assert.match(tab, /\*PLASTIC\n355000000, 0\n470000000, 0\.05\n510000000, 0\.15\n/);

  // B32: generated mid-side nodes; the importer splits each 3-node beam at its mid node
  const b32 = writeCalculixCase(model, { ...base, elementOrder: 2 });
  assert.deepEqual(b32.counts, { nodes: 9, elements: 4, appNodes: 5 });
  assert.match(b32.files[0].text, /\*ELEMENT, TYPE=B32, ELSET=EALL\n1, 1, 6, 2\n/);
  const f = byKind(parseAbaqus({ name: 'b32.inp', text: b32.files[0].text }), 'frame');
  assert.equal(f.members.length, 8);
  assert.equal(f.nodes.length / 3, 9);
  // the importer numbers nodes in its own order: compare the set of x coordinates and the member lengths
  assert.deepEqual(Array.from({ length: 9 }, (_, i) => f.nodes[3 * i]).sort((p, q) => p - q), [0, 1, 2, 3, 4, 5, 6, 7, 8]);
  for (const m of f.members) assert.ok(Math.abs(Math.abs(f.nodes[3 * m.b] - f.nodes[3 * m.a]) - 1) < 1e-12);

  // a vertical member cannot use the default (0, 0, -1) local 1-axis; generic sections become equivalent rectangles
  const tower = { kind: 'frame', name: 'tower', nodes: Float64Array.from([0, 0, 0, 0, 0, 10, 5, 0, 10]), members: [{ a: 0, b: 1, section: { A: 0.02, Iy: 4e-4, Iz: 1e-4, J: 3e-4 }, group: 'leg' }, { a: 1, b: 2, section: { type: 'tube', D: 0.3, t: 0.01 } }] };
  const tw = writeCalculixCase(tower, { material: STEEL, supports: [{ nodes: [0] }], loads: [{ type: 'force', nodes: [2], F: [0, 0, -1e3], M: [0, 50, 0] }] });
  // rectangle a × b with a b³ / 12 = Iy and b a³ / 12 = Iz:  b / a = 2, a⁴ = 12 Iy / 8  →  a = 0.15651, b = 0.31302
  const rect = /SECTION=RECT\n([0-9.]+), ([0-9.]+)\n0, 1, 0\n/.exec(tw.files[0].text);
  assert.ok(rect);
  const [a, b] = [Number(rect[1]), Number(rect[2])];
  assert.ok(Math.abs(a * b ** 3 / 12 - 4e-4) < 1e-9 && Math.abs(b * a ** 3 / 12 - 1e-4) < 1e-9);
  assert.ok(tw.notes.some((n) => /"leg".*no general beam section/.test(n)));
  assert.match(tw.files[0].text, /\n3, 3, -1000\n3, 5, 50\n/);
  const gen = writeCalculixCase(tower, { material: STEEL, supports: [{ nodes: [0] }], loads: [{ type: 'gravity' }], genericSection: 'general' });
  const leg = byKind(parseAbaqus({ name: 't.inp', text: gen.files[0].text }), 'frame').members[0].section;
  assert.deepEqual([leg.A, leg.Iy, leg.Iz, leg.J], [0.02, 4e-4, 1e-4, 3e-4]);
});

test('shell and solid models: S3 with thickness; C3D4/C3D10/C3D8/C3D20/C3D6, pyramids split, long records wrapped', () => {
  // two triangles forming a 1 × 1 plate
  const plate = { kind: 'surface', name: 'plate', positions: Float32Array.from([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0]), indices: Uint32Array.from([0, 1, 2, 0, 2, 3]), groups: [{ name: 'deck', start: 0, count: 6, thickness: 0.012 }] };
  const sh = writeCalculixCase(plate, { material: STEEL, supports: [{ nodes: [0, 3], dofs: 'pinned' }], loads: [{ type: 'pressure', p: 2500 }, { type: 'force', nodes: [1, 2], F: [0, 0, -10], share: true }] });
  const deck = sh.files[0].text;
  assert.match(deck, /\*ELEMENT, TYPE=S3, ELSET=EALL\n1, 1, 2, 3\n2, 1, 3, 4\n/);
  assert.match(deck, /\*SHELL SECTION, ELSET=SHELL1, MATERIAL=S355\n0\.012\n/);
  assert.match(deck, /\*BOUNDARY\nSUP1, 1, 1\nSUP1, 2, 2\nSUP1, 3, 3\n\*STEP/);
  assert.match(deck, /\*CLOAD\n2, 3, -5\n3, 3, -5\n\*DLOAD\nEALL, P, 2500\n/);
  const s = byKind(parseAbaqus({ name: 'p.inp', text: deck }), 'surface');
  assert.deepEqual(Array.from(s.indices), [0, 1, 2, 0, 2, 3]);
  assert.ok(Math.abs(s.groups[0].thickness - 0.012) < 1e-12);
  assert.throws(() => writeCalculixCase({ ...plate, groups: undefined }, { material: STEEL, supports: [{ nodes: [0] }], loads: [{ type: 'gravity' }] }), /thickness/);

  // one cell of each kind on disjoint node ranges (connectivity only matters here)
  const nNodes = 4 + 10 + 8 + 20 + 6 + 5;
  const positions = new Float32Array(3 * nNodes);
  for (let i = 0; i < nNodes; i++) { positions[3 * i] = i; positions[3 * i + 1] = (i * 7) % 5; positions[3 * i + 2] = (i * 3) % 4; }
  let o = 0;
  const conn = (n) => { const c = Uint32Array.from({ length: n }, (_, k) => o + k); o += n; return c; };
  const vol = { kind: 'volume', name: 'mix', positions, cells: [{ type: 'tet4', conn: conn(4) }, { type: 'tet10', conn: conn(10) }, { type: 'hex8', conn: conn(8) }, { type: 'hex20', conn: conn(20) }, { type: 'wedge6', conn: conn(6) }, { type: 'pyr5', conn: conn(5) }] };
  const v = writeCalculixCase(vol, { material: STEEL, supports: [{ nodes: [0, 1, 2], dofs: [1, 2, 3, 4] }], loads: [{ type: 'pressure', p: 1e5, elements: [2, 5], face: 1 }] });
  assert.equal(v.counts.elements, 7);
  assert.ok(v.notes.some((n) => /pyramid/.test(n)) && v.notes.some((n) => /rotational degrees of freedom/.test(n)));
  assert.match(v.files[0].text, /\*ELSET, ELSET=PRES1\n3, 6, 7\n/);
  assert.match(v.files[0].text, /\*DLOAD\nPRES1, P1, 100000\n/);
  assert.match(v.files[0].text, /\*SOLID SECTION, ELSET=EALL, MATERIAL=S355\n/);
  assert.match(v.files[0].text, /\*NODE FILE\nU, RF\n/);
  for (const l of v.files[0].text.split('\n')) assert.ok(l.length <= 132, 'CalculiX reads at most 132 characters per line');
  const back = parseAbaqus({ name: 'mix.inp', text: v.files[0].text });
  assert.deepEqual(back.meta.elementTypes, { tet4: 3, tet10: 1, hex8: 1, hex20: 1, wedge6: 1 });
  // the importer renumbers nodes, so connectivities are compared through the coordinates they point to
  const bv = byKind(back, 'volume');
  const coords = (pos, conn) => Array.from(conn).flatMap((i) => [pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]]);
  for (const [type, blk] of [['hex20', 3], ['tet10', 1], ['hex8', 2], ['wedge6', 4]]) assert.deepEqual(coords(bv.positions, bv.cells.find((c) => c.type === type).conn), coords(positions, vol.cells[blk].conn), type);
  const t4 = bv.cells.find((c) => c.type === 'tet4').conn;
  const py = vol.cells[5].conn;
  assert.deepEqual(coords(bv.positions, t4.subarray(4, 12)), coords(positions, [py[0], py[1], py[2], py[4], py[0], py[2], py[3], py[4]]), 'pyramid as two tetrahedra');
});

test('writeCalculixCase rejects inconsistent input', () => {
  const m = cantilever();
  const ok = { material: STEEL, supports: [{ nodes: [0] }], loads: [{ type: 'gravity' }] };
  assert.throws(() => writeCalculixCase({ kind: 'grid' }, ok), /model/);
  assert.throws(() => writeCalculixCase(m, { ...ok, analysis: 'modal' }), /analysis/);
  assert.throws(() => writeCalculixCase(m, { ...ok, material: { E: 2e11, nu: 0.5 } }), /material\.nu/);
  assert.throws(() => writeCalculixCase(m, { ...ok, material: { E: 2e11, nu: 0.3 } }), /material\.rho/);
  assert.throws(() => writeCalculixCase(m, { ...ok, supports: [] }), /supports/);
  assert.throws(() => writeCalculixCase(m, { ...ok, supports: [{ nodes: [9] }] }), /node index 9/);
  assert.throws(() => writeCalculixCase(m, { ...ok, supports: [{ nodes: 'BASE' }] }), /node set "BASE"/);
  assert.throws(() => writeCalculixCase(m, { ...ok, loads: [] }), /loads/);
  assert.throws(() => writeCalculixCase(m, { ...ok, loads: [{ type: 'pressure', p: 1 }] }), /pressure/);
  assert.throws(() => writeCalculixCase(m, { ...ok, analysis: 'dynamic' }), /timeStep/);
  assert.throws(() => writeCalculixCase(m, { ...ok, material: { ...STEEL, hardening: { Et: 3e11 } } }), /hardening\.Et/);
  // named node sets of an imported model are accepted
  const withSet = { ...m, sets: { nodes: { Base: Uint32Array.from([0]) } } };
  assert.match(writeCalculixCase(withSet, { ...ok, supports: [{ nodes: 'BASE' }] }).files[0].text, /\*NSET, NSET=SUP1\n1\n/);
});

/** One .frd data record in ccx's fixed format (1X, I2, I10, 6E12.5): negative values abut the previous field. */
const frdRow = (key, id, vals) => ` ${key}${String(id).padStart(10)}${vals.map((v) => v.toExponential(5).replace(/e([+-])(\d)$/, 'E$10$2').replace(/e([+-])(\d\d)$/, 'E$1$2').padStart(12)).join('')}`;

function frdFixture(nNodes, steps) {
  const L = ['    1Ccantilever', '    1UUSER', `    2C${String(nNodes + 2).padStart(30)}${' '.repeat(37)}1`];
  for (let i = 1; i <= nNodes + 2; i++) L.push(frdRow('-1', i, [2 * (i - 1), 0, 0]));
  L.push(' -3', `    3C${String(4).padStart(30)}${' '.repeat(37)}1`, ' -1         1   11    0    1', ' -2         1         2', ' -3');
  steps.forEach((s, k) => {
    const head = (name, nc) => { L.push(`    1PSTEP${String(k + 1).padStart(26)}${String(1).padStart(12)}${String(k + 1).padStart(12)}`, `  100CL  101${s.time.toExponential(5).toUpperCase().padStart(12)}${String(nNodes + 2).padStart(12)}${' '.repeat(20)} 0${String(k + 1).padStart(5)}${' '.repeat(10)} 1`, ` -4  ${name.padEnd(8)}${String(nc).padStart(5)}    1`); };
    head('DISP', 4);
    L.push(' -5  D1          1    2    1    0', ' -5  D2          1    2    2    0', ' -5  D3          1    2    3    0', ' -5  ALL         1    2    0    0    1ALL');
    for (let i = 1; i <= nNodes + 2; i++) L.push(frdRow('-1', i, s.disp(i)));
    L.push(' -3');
    if (s.stress) {
      head('STRESS', 6);
      for (const c of ['SXX', 'SYY', 'SZZ', 'SXY', 'SYZ', 'SZX']) L.push(` -5  ${c.padEnd(8)}    1    4    1    1`);
      for (let i = 1; i <= nNodes; i++) L.push(frdRow('-1', i, s.stress(i)));
      L.push(' -3');
    }
  });
  L.push(' 9999');
  return L.join('\n');
}

test('readCalculixFrd: DISP and STRESS blocks → node-aligned arrays, von Mises colours and deformed shape', () => {
  const model = cantilever();
  // step 1: half load; step 2: tip deflection w(x) = -1e-3 (x/8)² with a small axial part; stress only in step 2
  const disp = (f) => (i) => [-1.5e-5 * f * (i - 1), 0, -1e-3 * f * ((2 * (i - 1)) / 8) ** 2];
  const stress = (i) => [-1.2e8 * (5 - i) / 4, 0, 0, 2.5e6, 0, 0];
  const text = frdFixture(5, [{ time: 0.5, disp: disp(0.5) }, { time: 1, disp: disp(1), stress }]);
  // fixed-width records really abut: "…E+00-1.50000E-05"
  assert.ok(/\d-\d\.\d{5}E/.test(text));
  const r = readCalculixFrd(text, { nodeCount: 5, positions: model, scale: 100 });
  assert.equal(r.steps.length, 2);
  assert.deepEqual(r.steps.map((s) => s.time), [0.5, 1]);
  assert.deepEqual(r.steps[1].fields.DISP.components, ['D1', 'D2', 'D3']);
  assert.deepEqual(r.steps[1].fields.STRESS.components, ['SXX', 'SYY', 'SZZ', 'SXY', 'SYZ', 'SZX']);
  assert.equal(r.nodeCount, 5);
  assert.equal(r.displacement.length, 15);
  assert.ok(Math.abs(r.displacement[14] + 1e-3) < 1e-12 && Math.abs(r.displacement[12] + 6e-5) < 1e-12);
  // von Mises of a uniaxial stress with shear: sqrt(sxx² + 3 sxy²)
  assert.ok(Math.abs(r.nodeValues[0] - Math.hypot(1.2e8, Math.sqrt(3) * 2.5e6)) < 1);
  assert.ok(Math.abs(r.nodeValues[4] - Math.sqrt(3) * 2.5e6) < 1);
  assert.equal(r.valueLabel, 'von Mises stress');
  assert.equal(r.valueUnit, 'Pa');
  assert.equal(r.deformed.length, 5);
  assert.ok(Math.abs(r.deformed[4][2] + 0.1) < 1e-9 && Math.abs(r.deformed[4][0] - (8 - 6e-3)) < 1e-9);
  // an earlier step without stresses falls back to the displacement magnitude
  const s0 = readCalculixFrd(text, { nodeCount: 5, step: 0 });
  assert.equal(s0.valueLabel, 'Displacement magnitude');
  assert.equal(s0.stress, null);
  assert.equal(s0.deformed, null);
  assert.ok(Math.abs(s0.nodeValues[4] - Math.hypot(3e-5, 5e-4)) < 1e-12);
  // without nodeCount the largest node id in the file is used (here two extra nodes)
  assert.equal(readCalculixFrd(text).nodeCount, 7);
  assert.throws(() => readCalculixFrd('*HEADING\n'), /no result blocks/);
  assert.throws(() => readCalculixFrd(text, { step: 5 }), /step/);
});

test('readCalculixDat: displacement, reaction and integration-point stress tables', () => {
  const model = cantilever();
  const e = (v) => v.toExponential(6).toUpperCase().padStart(14);
  const L = ['', ' displacements (vx,vy,vz) for set NALL and time  0.1000000E+01', ''];
  for (let i = 1; i <= 5; i++) L.push(`${String(i).padStart(10)}${e(0)}${e(0)}${e(-2e-3 * ((i - 1) / 4) ** 2)}`);
  L.push('', ' forces (fx,fy,fz) for set NSUP and time  0.1000000E+01', '', `${String(1).padStart(10)}${e(0)}${e(0)}${e(1e4)}`);
  L.push('', ' total force (fx,fy,fz) for set NSUP and time  0.1000000E+01', '', `      ${e(0)}${e(0)}${e(1e4)}`);
  L.push('', ' stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set EALL and time  0.1000000E+01', '');
  for (let el = 1; el <= 4; el++) for (let ip = 1; ip <= 2; ip++) L.push(`${String(el).padStart(10)}${String(ip).padStart(4)}${e((ip === 1 ? 1 : -1) * 1e8 / el)}${e(0)}${e(0)}${e(0)}${e(0)}${e(0)}`);
  const r = readCalculixDat(L.join('\n'), { positions: model.nodes, scale: 10 });
  assert.deepEqual(r.tables.map((t) => [t.kind, t.set, t.rows.length]), [['displacements', 'NALL', 5], ['forces', 'NSUP', 1], ['total force', 'NSUP', 1], ['stresses', 'EALL', 8]]);
  assert.equal(r.nodeCount, 5);
  assert.equal(r.time, 1);
  assert.ok(Math.abs(r.displacement[14] + 2e-3) < 1e-12);
  assert.deepEqual(r.reactions, [{ node: 0, F: [0, 0, 1e4] }]);
  assert.deepEqual(r.totalForce, [0, 0, 1e4]);
  assert.deepEqual(r.elementStress.element, [0, 1, 2, 3]);
  assert.ok(Math.abs(r.elementStress.vonMises[1] - 5e7) < 1e-3);
  assert.equal(r.valueUnit, 'm');
  assert.ok(Math.abs(r.nodeValues[4] - 2e-3) < 1e-12);
  assert.ok(Math.abs(r.deformed[4][2] + 0.02) < 1e-12);
  assert.throws(() => readCalculixDat('nothing here'), /no result tables/);
});
