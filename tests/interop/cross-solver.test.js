// Code-to-code verification against files produced by the external solvers (see fixtures/README.md for the
// provenance of every file). The solvers are not needed to run these tests: the application's readers parse the
// stored outputs, and the built-in models must stay within the differences recorded in src/spec/cross-solver.json.

import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { SOLVER_EXPORTS, exportCase, readNemohResults, readCapytaineCsv, readCalculixDat, readCalculixFrd, writeCalculixCase, wettedMesh } from '../../src/interop/index.js';
import { addedMassBem } from '../../src/engine/hydro/bem.js';
import { buildFloater, floaterRAO, boxBargeRAO } from '../../src/engine/hydro/floater.js';
import { raoSolve, hydroAt } from '../../src/engine/hydro/hydrodb.js';
import { analyzeFrame } from '../../src/engine/structural/frame3d.js';
import { modalAnalysis } from '../../src/engine/structural/framemodal.js';
import { bucklingAnalysis } from '../../src/engine/structural/framebuckling.js';
import { tubeSection, rectSection } from '../../src/engine/structural/sections.js';
import { navierPlate } from '../../src/engine/structural/plates.js';
import { barge, cylinder, sphere } from '../../src/mesh/generate.js';

const fx = (p) => readFileSync(new URL(`./fixtures/${p}`, import.meta.url), 'utf8');
const SPEC = JSON.parse(readFileSync(new URL('../../src/spec/cross-solver.json', import.meta.url), 'utf8'));
const comparison = (id) => { const c = SPEC.comparisons.find((x) => x.id === id); assert.ok(c, `comparison ${id}`); return c; };
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

/**
 * Checks recomputed values against one recorded comparison: the reference read from the fixture must be the
 * recorded one, and the built-in model must not be further from it than recorded (plus `slack`, absolute on the
 * relative difference), so that a change which degrades the agreement fails here.
 */
function within(id, label, builtIn, reference, slack = 2e-3) {
  const c = comparison(id).cases.find((x) => x.label === label);
  assert.ok(c, `${id}: case "${label}"`);
  assert.ok(rel(reference, c.reference) < 1e-9, `${id} / ${label}: reference ${reference} is not the recorded ${c.reference}`);
  if (c.relDiff === null) return;
  const d = (builtIn - reference) / Math.abs(reference);
  assert.ok(Math.abs(d) <= Math.abs(c.relDiff) + slack, `${id} / ${label}: built-in differs by ${d}, recorded ${c.relDiff}`);
}

// --- the cases that were exported, solved and read back ---------------------------------------------------------
const RHO = 1025; const G = 9.80665;
const PERIODS = [4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 25];
const OM = PERIODS.map((T) => Number((2 * Math.PI / T).toPrecision(10)));
const HULLS = {
  barge: () => barge({ length: 20, beam: 8, depth: 4, draft: 2, divisions: [20, 8, 8], name: 'barge' }),
  cylinder: () => cylinder({ radius: 4, height: 44, base: [0, 0, -40], segments: 24, heightDivisions: 22, capRings: 4, name: 'cylinder' }),
  hemisphere: () => sphere({ radius: 5, subdivisions: 4, name: 'hemisphere' }),
};
const DOF = ['A11 surge', 'A22 sway', 'A33 heave', 'A44 roll', 'A55 pitch', 'A66 yaw'];
const capDb = (name) => readCapytaineCsv(fx(`capytaine/${name}_hydro.csv`));

const STEEL = { name: 'S355', E: 2.1e11, nu: 0.3, rho: 7850 };
const TUBE = { type: 'tube', D: 0.5, t: 0.02 };
function chain(points, nPer, section) {
  const xyz = [...points[0]]; const members = [];
  for (let s = 0; s + 1 < points.length; s++) {
    for (let k = 1; k <= nPer; k++) {
      const f = k / nPer;
      xyz.push(...[0, 1, 2].map((q) => points[s][q] + f * (points[s + 1][q] - points[s][q])));
      const n = xyz.length / 3 - 1;
      members.push({ a: n - 1, b: n, section });
    }
  }
  return { kind: 'frame', nodes: Float64Array.from(xyz), members };
}
const cantilever = (section = TUBE) => ({ ...chain([[0, 0, 0], [10, 0, 0]], 10, section), name: 'cantilever' });
const portal = () => ({ ...chain([[0, 0, 0], [0, 0, 5], [6, 0, 5], [6, 0, 0]], 5, TUBE), name: 'portal' });
function plate(a = 2, b = 1.5, nx = 12, ny = 10, t = 0.02) {
  const pos = []; const idx = []; const edge = [];
  const id = (i, j) => j * (nx + 1) + i;
  for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++) { pos.push(a * i / nx, b * j / ny, 0); if (i === 0 || j === 0 || i === nx || j === ny) edge.push(id(i, j)); }
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    if ((i + j) % 2 === 0) idx.push(id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j), id(i + 1, j + 1), id(i, j + 1));
    else idx.push(id(i, j), id(i + 1, j), id(i, j + 1), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
  }
  return { kind: 'surface', name: 'plate', positions: Float64Array.from(pos), indices: Uint32Array.from(idx), groups: [{ name: 'plate', start: 0, count: idx.length, thickness: t }], edge, centre: id(nx / 2, ny / 2), a, b, t };
}
/** The same frame as a model of the frame engine. */
function engineModel(m, section, fixedNodes, nodal) {
  const nodes = [];
  for (let i = 0; i < m.nodes.length / 3; i++) nodes.push([m.nodes[3 * i], m.nodes[3 * i + 1], m.nodes[3 * i + 2]]);
  return { nodes, materials: { steel: { E: STEEL.E, nu: STEEL.nu, rho: STEEL.rho } }, members: m.members.map((e) => ({ a: e.a, b: e.b, section, material: 'steel' })), supports: fixedNodes.map((n) => ({ node: n, fix: 'fixed' })), loadCases: [{ name: 'L', nodal }] };
}

test('cross-solver.json: every comparison is complete and consistent', () => {
  assert.deepEqual(SPEC.solvers.map((s) => s.name), ['Capytaine', 'NEMOH', 'CalculiX CrunchiX (ccx)']);
  assert.ok(SPEC.comparisons.length >= 26);
  for (const c of SPEC.comparisons) {
    for (const k of ['id', 'title', 'builtIn', 'reference', 'quantity', 'unit', 'validity', 'note']) assert.ok(typeof c[k] === 'string' && c[k].length > 0, `${c.id}.${k}`);
    assert.ok(['agrees', 'agrees within stated range', 'differs'].includes(c.verdict), c.id);
    assert.ok(c.cases.length > 0);
    for (const k of c.cases) {
      assert.ok(Number.isFinite(k.builtIn) && Number.isFinite(k.reference), `${c.id} / ${k.label}`);
      if (k.relDiff !== null) assert.ok(Math.abs(k.relDiff - (k.builtIn - k.reference) / Math.abs(k.reference)) < 1e-12, `${c.id} / ${k.label}`);
    }
  }
  assert.equal(new Set(SPEC.comparisons.map((c) => c.id)).size, SPEC.comparisons.length);
});

test('Capytaine 3.0.0 output: the CSV of the exported script is read to the recorded limits, and the phase convention is exp(+i omega t)', () => {
  const db = capDb('barge');
  assert.equal(db.source, 'capytaine 3.0.0');
  assert.equal(db.omega.length, 12);
  assert.deepEqual(db.headingsDeg, [0, 90, 180]);
  assert.equal(db.rho, 1025); assert.equal(db.depth, null);
  assert.ok(db.Ainf && db.A0, 'both frequency limits are in the file');
  // longest wave of the run, omega = 0.1 rad/s (wavelength 6.2 km), heading 0 = propagation towards +x
  const w = db.omega[0]; const A = db.A[0]; const B = db.B[0];
  assert.equal(w, 0.1);
  const F = { re: db.Fre[0][0], im: db.Fim[0][0] };
  // heave: hydrostatic limit rho g Awp in phase with the elevation (pressure at the keel, exp(-k T)), reduced by the
  // added-mass force …
  const k = w * w / G;
  assert.ok(rel(F.re[2], RHO * G * 160 * Math.exp(-2 * k) - w * w * A[2][2]) < 1e-3);
  // … and LEADING the elevation by the radiation-damping force + i omega B33 (it would lag under exp(-i omega t))
  assert.ok(F.im[2] > 0 && rel(F.im[2], w * B[2][2]) < 0.01, `Im F3 = ${F.im[2]}, omega B33 = ${w * B[2][2]}`);
  // surge: inertia force of the displaced and added mass in the horizontal particle acceleration + i omega² (at z ≈ 0)
  assert.ok(Math.abs(F.re[0]) < 1e-3 * F.im[0] && rel(F.im[0], w * w * (RHO * 320 + A[0][0])) < 5e-3);
  // pitch moment follows the wave slope: + i for waves towards +x, − i towards −x; beam sea rolls instead
  assert.ok(F.im[4] > 0 && db.Fim[2][0][4] < 0 && rel(-db.Fim[2][0][4], F.im[4]) < 1e-4);
  assert.ok(db.Fim[1][0][1] > 0 && Math.abs(db.Fim[1][0][0]) < 1e-3 * db.Fim[1][0][1]);
});

test('NEMOH v3 output of the exported barge: both result layouts are read, and agree with Capytaine through the other reader', () => {
  const files = ['Nemoh.cal', 'RadiationCoefficients.tec', 'ExcitationForce.tec', 'CM.dat', 'CA.dat'].map((p) => ({ path: p === 'Nemoh.cal' ? p : `results/${p}`, text: fx(`nemoh/${p}`) }));
  const n = SOLVER_EXPORTS.find((e) => e.id === 'nemoh').read(files);
  assert.equal(n.source, 'NEMOH (RadiationCoefficients.tec, ExcitationForce.tec)');
  assert.deepEqual(n.omega, [0.25, 0.5, 0.75, 1, 1.25]);
  assert.deepEqual(n.headingsDeg, [0, 90, 180]);
  assert.equal(n.rho, 1025); assert.equal(n.g, 9.80665); assert.equal(n.depth, null);
  // CM.dat / CA.dat carry the same coefficients (six significant digits)
  const m = readNemohResults({ files: files.filter((f) => !/Radiation/.test(f.path)) });
  assert.equal(m.source, 'NEMOH (CM.dat, CA.dat, ExcitationForce.tec)');
  assert.deepEqual(m.omega, n.omega);
  for (let w = 0; w < 5; w++) for (let i = 0; i < 6; i++) { assert.ok(rel(m.A[w][i][i], n.A[w][i][i]) < 1e-5); assert.ok(rel(m.B[w][i][i], n.B[w][i][i]) < 1e-5); }
  // … in the same row / column order (surge–pitch coupling, which is not exactly symmetric in the output)
  for (let w = 0; w < 5; w++) { assert.ok(rel(m.A[w][0][4], n.A[w][0][4]) < 1e-5); assert.ok(rel(m.A[w][4][0], n.A[w][4][0]) < 1e-5); }
  // the phase column is in radians and, once conjugated by the reader, follows the convention of the application
  assert.ok(n.Fim[0][0][2] > 0 && rel(n.Fim[0][0][2], n.omega[0] * n.B[0][2][2]) < 0.01);
  assert.ok(n.Fim[0][0][0] > 0 && n.Fim[2][0][0] < 0);

  const c = readCapytaineCsv(fx('capytaine/barge5_hydro.csv'));
  assert.deepEqual(c.omega, n.omega);
  const rec = comparison('nemoh-capytaine-barge');
  for (const k of rec.cases) {
    const mt = /^omega = ([\d.]+) rad\/s, (?:heading 0, )?(A|B|Re F|Im F)(\d)/.exec(k.label);
    const w = n.omega.indexOf(Number(mt[1])); const i = Number(mt[3]) - 1;
    const pick = (db) => (mt[2] === 'A' ? db.A[w][i][i] : mt[2] === 'B' ? db.B[w][i][i] : mt[2] === 'Re F' ? db.Fre[0][w][i] : db.Fim[0][w][i]);
    assert.ok(rel(pick(n), k.builtIn) < 1e-9 && rel(pick(c), k.reference) < 1e-9, k.label);
    if (k.relDiff !== null) assert.ok(Math.abs(k.relDiff) < 3e-3, `${k.label}: NEMOH and Capytaine differ by ${k.relDiff}`);
  }
  // the surge–pitch coupling is not symmetric on this mesh (11 %) and both solvers put the larger term in A51
  for (let w = 0; w < 5; w++) { assert.ok(n.A[w][4][0] > 1.05 * n.A[w][0][4] && c.A[w][4][0] > 1.05 * c.A[w][0][4]); assert.ok(rel(n.A[w][0][4], c.A[w][0][4]) < 2e-3 && rel(n.A[w][4][0], c.A[w][4][0]) < 2e-3); }
});

test('built-in panel method stays within the recorded differences from Capytaine (barge both limits, cylinder infinite)', () => {
  for (const [name, limit, key] of [['barge', 'infinite', 'Ainf'], ['barge', 'zero', 'A0'], ['cylinder', 'infinite', 'Ainf']]) {
    const ref = capDb(name)[key];
    const r = addedMassBem({ surface: HULLS[name](), rho: RHO, limit, cog: [0, 0, 0] });
    for (let i = 0; i < 6; i++) within(`added-mass-${limit}`, `${name} ${DOF[i]}`, r.A[i][i], ref[i][i]);
  }
  // the hemisphere references are stored as read (the 2592-panel solve itself is covered by hydro.bem.test.js)
  const h = capDb('hemisphere');
  for (const [limit, key] of [['infinite', 'Ainf'], ['zero', 'A0']]) for (let i = 0; i < 6; i++) {
    const k = comparison(`added-mass-${limit}`).cases.find((x) => x.label === `hemisphere ${DOF[i]}`);
    assert.ok(rel(h[key][i][i], k.reference) < 1e-9);
    assert.equal(k.relDiff === null, i >= 3, 'terms that vanish by symmetry carry no relative difference');
  }
  // exact value rho pi a³ / 3 for heave at infinite frequency: the built-in result recorded is the closer one
  const exact = RHO * Math.PI * 125 / 3;
  const k = comparison('added-mass-infinite').cases.find((x) => x.label === 'hemisphere A33 heave');
  assert.ok(rel(k.builtIn, exact) < 2e-3 && rel(k.reference, exact) > 0.01);
});

test('box panel model: barge RAOs stay within the recorded differences from Capytaine (three hulls)', () => {
  const hulls = [
    ['barge', 20, 8, 2, PERIODS], ['barge100', 100, 20, 5, [5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 25]], ['barge80', 80, 24, 3, [4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20]],
  ];
  for (const [name, L, B, T, periods] of hulls) {
    const db = capDb(name);
    const om = periods.map((t) => Number((2 * Math.PI / t).toPrecision(10)));
    const m = RHO * L * B * T; const z6 = () => Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
    // centre of gravity on the waterline, dry gyradii 0.35 B and 0.25 L in the model and in the reference
    const GMt = B * B / (12 * T) - T / 2; const GMl = L * L / (12 * T) - T / 2;
    const opt = { L, B, T, omega: om, GMt, kxxDry: 0.35 * B, rho: RHO, g: G };
    const head = boxBargeRAO({ ...opt, headingDeg: 180 }); const beam = boxBargeRAO({ ...opt, headingDeg: 90 });
    assert.equal(head.method, 'panel');
    const M = z6(); const C = z6(); const Bext = z6();
    [m, m, m, m * (0.35 * B) ** 2, m * (0.25 * L) ** 2, m * (0.25 * L) ** 2].forEach((v, i) => { M[i][i] = v; });
    C[2][2] = RHO * G * L * B; C[3][3] = m * G * GMt; C[4][4] = m * G * GMl;
    // the added roll damping of the model (5 % of critical of its own oscillator) goes into the reference as well
    Bext[3][3] = 2 * 0.05 * Math.sqrt((M[3][3] + beam.rollAddedInertia) * C[3][3]);
    const refHead = raoSolve({ db, M, C, Bext, omega: om, headingDeg: 180 }); const refBeam = raoSolve({ db, M, C, Bext, omega: om, headingDeg: 90 });
    periods.forEach((t, i) => {
      within(`${name}-rao-heave`, `T = ${t} s`, head.heave[i], refHead.amp[2][i]);
      within(`${name}-rao-pitch`, `T = ${t} s`, head.pitch[i], refHead.amp[4][i]);
      within(`${name}-rao-heaveBeam`, `T = ${t} s`, beam.heave[i], refBeam.amp[2][i]);
      within(`${name}-rao-roll`, `T = ${t} s`, beam.roll[i], refBeam.amp[3][i]);
    });
    // the range the model states: down to the wave whose length is the beam
    assert.ok(rel(beam.validPeriodMin.heave, Math.sqrt(2 * Math.PI * B / G)) < 1e-12 && beam.validPeriodMin.roll === beam.validPeriodMin.heave);
    // its parts against the solver's own coefficients at the longest period: added mass, and the roll inertia of
    // the free sway–roll pair against the same combination of the reference
    const h = hydroAt(db, om[om.length - 1], 90);
    assert.ok(rel(beam.heaveNaturalPeriod, 2 * Math.PI * Math.sqrt((m + hydroAt(db, 2 * Math.PI / beam.heaveNaturalPeriod, 90).A[2][2]) / C[2][2])) < 0.01, `${name}: heave natural period`);
    assert.ok(h.A[3][3] > 0 && beam.rollAddedInertia > 0.5 * h.A[3][3] && beam.rollAddedInertia < 1.1 * h.A[3][3], `${name}: roll added inertia`);
  }
  // the recorded agreement in numbers: heave and pitch within 5 % at every period in head seas for all three hulls
  for (const name of ['barge', 'barge100', 'barge80']) for (const key of ['heave', 'pitch']) for (const k of comparison(`${name}-rao-${key}`).cases) assert.ok(Math.abs(k.relDiff) <= 0.05, `${name} ${key} ${k.label}`);
  // beam seas: heave within 5 % except the wide barge at 4 s (wavelength equal to its beam); roll within 5 % from
  // 5 s (20 × 8 × 2), 7 s (100 × 20 × 5) and at every period (80 × 24 × 3), and never further than 12 %
  for (const k of comparison('barge-rao-heaveBeam').cases) assert.ok(Math.abs(k.relDiff) <= 0.05, k.label);
  for (const k of comparison('barge100-rao-heaveBeam').cases) assert.ok(Math.abs(k.relDiff) <= 0.05, k.label);
  for (const k of comparison('barge80-rao-heaveBeam').cases) assert.equal(Math.abs(k.relDiff) <= 0.05, k.label !== 'T = 4 s', k.label);
  for (const [name, from] of [['barge', 5], ['barge100', 7], ['barge80', 4]]) for (const k of comparison(`${name}-rao-roll`).cases) { assert.equal(Math.abs(k.relDiff) <= 0.05, Number(/\d+/.exec(k.label)[0]) >= from, `${name} roll ${k.label}`); assert.ok(Math.abs(k.relDiff) < 0.12, k.label); }
});

test('closed-form barge model (method closed-form): the differences from Capytaine that its validity limits state', () => {
  const db = capDb('barge');
  const m = RHO * 320; const z6 = () => Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  const GMt = 1 + 64 / 24 - 2; const GMl = 1 + 400 / 24 - 2;
  const head = boxBargeRAO({ L: 20, B: 8, T: 2, omega: OM, headingDeg: 180, method: 'closed-form', rho: RHO, g: G });
  const beam = boxBargeRAO({ L: 20, B: 8, T: 2, omega: OM, headingDeg: 90, GMt, kxxDry: 0.35 * 8, method: 'closed-form', rho: RHO, g: G });
  const M = z6(); const C = z6(); const Bext = z6();
  [m, m, m, m * (0.35 * 8) ** 2, m * 5 ** 2, m * 5 ** 2].forEach((v, i) => { M[i][i] = v; });
  C[2][2] = RHO * G * 160; C[3][3] = m * G * GMt; C[4][4] = m * G * GMl;
  Bext[3][3] = 2 * 0.05 * Math.sqrt((M[3][3] + beam.rollAddedInertia) * C[3][3]);
  const refHead = raoSolve({ db, M, C, Bext, omega: OM, headingDeg: 180 }); const refBeam = raoSolve({ db, M, C, Bext, omega: OM, headingDeg: 90 });
  const d = (a, b) => a / b - 1;
  // heave within 5 % from 8 s in head and in beam seas (twice the heave natural period of the model, the limit it
  // returns) and low below; pitch from 6 s; roll from 7 s (1.3 roll natural periods) and within 16 % below
  PERIODS.forEach((T, i) => {
    assert.equal(Math.abs(d(head.heave[i], refHead.amp[2][i])) <= 0.05, T >= 8, `heave ${T} s`);
    assert.equal(Math.abs(d(beam.heave[i], refBeam.amp[2][i])) <= 0.05, T >= 8, `heave, beam sea ${T} s`);
    if (T < 8) assert.ok(d(head.heave[i], refHead.amp[2][i]) < 0 && d(beam.heave[i], refBeam.amp[2][i]) < 0, `low at ${T} s`);
    assert.equal(Math.abs(d(head.pitch[i], refHead.amp[4][i])) <= 0.05, T >= 6, `pitch ${T} s`);
    assert.equal(Math.abs(d(beam.roll[i], refBeam.amp[3][i])) <= 0.05, T >= 7, `roll ${T} s`);
    assert.ok(Math.abs(d(beam.roll[i], refBeam.amp[3][i])) < 0.16, `roll ${T} s`);
  });
  assert.ok(beam.validPeriodMin.heave > 7 && beam.validPeriodMin.heave <= 8.03 && beam.validPeriodMin.roll > 6 && beam.validPeriodMin.roll <= 7);
  // the parts of the roll model against the solver's own coefficients: strip added inertia (no end effect) 6–11 %
  // above A44, Haskind radiation damping within 12 % of B44 from 5 s
  for (const w of OM) { const dd = beam.rollAddedInertia / hydroAt(db, w, 90).A[3][3] - 1; assert.ok(dd > 0.05 && dd < 0.11, `A44 at ${w} rad/s: ${dd}`); }
  PERIODS.forEach((T, i) => { if (T >= 5) assert.ok(rel(beam.rollRadiationDamping[i], hydroAt(db, OM[i], 90).B[3][3]) < 0.12, `B44 at T = ${T} s`); });
});

test('slender-member spar model stays within the recorded differences from Capytaine; excitation phases agree in sign', () => {
  const db = capDb('cylinder');
  const model = buildFloater({ members: [{ type: 'cylinder', p1: [0, 0, -40], p2: [0, 0, 4], D: 8, Ca: 1, end1: 'hemisphere', free1: true }], cog: [0, 0, -25], gyradii: [12, 12, 3], rho: RHO, g: G });
  // radiation damping from the Haskind relation, nothing added: the same physics as the reference
  const r = floaterRAO({ model, omega: OM.slice().reverse(), headingDeg: 0, viscous: false });
  const ref = raoSolve({ db, M: model.M, C: model.C, omega: OM, headingDeg: 0 });
  PERIODS.forEach((T, i) => {
    const q = OM.length - 1 - i; // floaterRAO was given ascending frequencies
    const h = hydroAt(db, OM[i], 0); const F = h.F;
    within('spar-force-surge', `T = ${T} s`, Math.hypot(r.force.re[0][q], r.force.im[0][q]), Math.hypot(F.re[0], F.im[0]));
    within('spar-force-heave', `T = ${T} s`, Math.hypot(r.force.re[2][q], r.force.im[2][q]), Math.hypot(F.re[2], F.im[2]));
    within('spar-rao-surge', `T = ${T} s`, r.amp[0][q], ref.amp[0][i]);
    within('spar-rao-heave', `T = ${T} s`, r.amp[2][q], ref.amp[2][i]);
    within('spar-rao-pitch', `T = ${T} s`, r.amp[4][q], ref.amp[4][i]);
    // the app's own long-wave model puts the surge force 90° ahead of the crest; the converted Capytaine force
    // has the same sign of phase at every period and tends to +90° in long waves
    assert.ok(r.force.im[0][q] > 0 && F.im[0] > 0, `surge force phase at T = ${T} s`);
    if (T >= 12) assert.ok(Math.abs(Math.atan2(F.im[0], F.re[0]) - Math.PI / 2) < 0.01);
    // the solver's own coefficients satisfy the Haskind relation of a body of revolution in deep water,
    // B11 = omega³ |X1|² / (4 rho g³) (also pitch) and B33 = omega³ |X3|² / (2 rho g³), to the accuracy of its mesh:
    // 2–3 % in surge and pitch, 6 % in heave (9 % at 7 s) …
    const c = OM[i] ** 3 / (RHO * G ** 3);
    for (const d of [0, 4]) assert.ok(rel(c * (F.re[d] ** 2 + F.im[d] ** 2) / 4, h.B[d][d]) < 0.03, `Haskind, mode ${d + 1}, T = ${T} s`);
    if (T >= 7) assert.ok(rel(c * (F.re[2] ** 2 + F.im[2] ** 2) / 2, h.B[2][2]) < 0.1, `Haskind, heave, T = ${T} s`);
    // … and the damping the model derives from its own excitation is within 12 % of them from 6 s, 8 % from 10 s
    if (T >= 6) for (const d of [0, 2, 4]) assert.ok(rel(r.radiationDamping[q][d][d], h.B[d][d]) < (T >= 10 ? 0.08 : 0.12), `B${d + 1}${d + 1} at T = ${T} s`);
  });
  // heave within 5 % at every period from 5 s, the resonance (13.1 s) included; surge and pitch within 5 % outside
  // 14–20 s, where the narrow pitch resonance (14.8 s against 14.6 s) decides
  for (const k of comparison('spar-rao-heave').cases) assert.equal(Math.abs(k.relDiff) <= 0.05, Number(/\d+/.exec(k.label)[0]) >= 5, k.label);
  for (const k of comparison('spar-force-heave').cases) assert.equal(Math.abs(k.relDiff) <= 0.05, Number(/\d+/.exec(k.label)[0]) >= 5, k.label);
  for (const id of ['spar-rao-surge', 'spar-rao-pitch']) for (const k of comparison(id).cases) { const T = Number(/\d+/.exec(k.label)[0]); if (T >= 5 && (T <= 12 || T >= 25)) assert.ok(Math.abs(k.relDiff) <= 0.05, `${id} ${k.label}`); assert.ok(Math.abs(k.relDiff) < 0.25, `${id} ${k.label}`); }
  // strip added mass against the solver at 14 s: surge within 4 %, pitch within 3 % with the free-end loss (15 % without)
  const h14 = hydroAt(db, OM[PERIODS.indexOf(14)], 0);
  assert.ok(rel(model.A[0][0], h14.A[0][0]) < 0.04 && rel(model.A[4][4], h14.A[4][4]) < 0.03 && rel(model.A[2][2], h14.A[2][2]) < 0.03);
  const plain = buildFloater({ members: [{ type: 'cylinder', p1: [0, 0, -40], p2: [0, 0, 4], D: 8, Ca: 1 }], cog: [0, 0, -25], gyradii: [12, 12, 3], rho: RHO, g: G });
  assert.ok(plain.A[4][4] / h14.A[4][4] > 1.13);
});

test('CalculiX 2.23 output: printed tables and the .frd file of decks written by the application are read back', () => {
  const st = readCalculixDat(fx('calculix/cantilever_static.dat'), { nodeCount: 11 });
  assert.deepEqual(st.tables.map((t) => `${t.kind}:${t.rows.length}`), ['displacements:21', 'forces:1', 'total force:1', 'stresses:160']);
  assert.equal(st.displacement[3 * 10 + 2], -0.01807179);
  assert.ok(rel(st.totalForce[2], 1e4) < 1e-6, 'reaction balances the 10 kN tip load');
  assert.deepEqual(st.reactions.map((r) => r.node), [0]);
  assert.equal(st.elementStress.element.length, 10);
  assert.equal(st.frequencies, null); assert.equal(st.bucklingFactors, null);
  // the same run through the .frd file (five significant digits there), via the registry's format detection
  const frd = SOLVER_EXPORTS.find((e) => e.id === 'calculix').read(fx('calculix/cantilever_static.frd'), { nodeCount: 11 });
  assert.deepEqual(Object.keys(frd.steps[0].fields), ['DISP', 'STRESS', 'TOSTRAIN', 'FORC', 'ERROR']);
  assert.ok(rel(frd.displacement[3 * 10 + 2], st.displacement[3 * 10 + 2]) < 1e-5);
  assert.equal(frd.valueLabel, 'von Mises stress');
  for (let i = 0; i < 11; i++) assert.ok(rel(frd.displacement[3 * i + 2] || 1e-30, st.displacement[3 * i + 2] || 1e-30) < 1e-4);

  // *FREQUENCY and *BUCKLE steps print no nodal tables: eigenvalues, frequencies and factors are read instead
  const fr = readCalculixDat(fx('calculix/cantilever_frequency.dat'));
  assert.deepEqual(fr.frequencies, [4.942547, 4.942547, 30.52148, 30.52148, 80.19153, 83.59951]);
  assert.ok(rel(Math.sqrt(fr.eigenvalues[0]) / (2 * Math.PI), fr.frequencies[0]) < 1e-6);
  assert.equal(fr.displacement, null);
  const bk = readCalculixDat(fx('calculix/cantilever_buckle.dat'));
  assert.deepEqual(bk.bucklingFactors, [4.531915, 4.531915, 39.94306, 39.94306]);
  assert.deepEqual(readCalculixDat(fx('calculix/portal_buckle.dat')).bucklingFactors, [18.26507, 43.81909, 51.59644, 152.1445]);
  // shell rows end with the name of the expanded element; the displacement table of the plate is read on the app's nodes
  const pl = readCalculixDat(fx('calculix/plate_static.displacements.dat'), { nodeCount: 143 });
  assert.equal(pl.tables[0].rows.length, 525);
  assert.equal(pl.displacement.length, 3 * 143);
  assert.equal(pl.displacement[3 * 71 + 2], -0.002167075);
});

test('frame statics, natural frequencies, buckling factors and plate bending stay within the recorded differences from CalculiX', () => {
  const dat = (j, n) => readCalculixDat(fx(`calculix/${j}.dat`), { nodeCount: n });
  const tube = tubeSection({ D: 0.5, t: 0.02 }); const rect = rectSection({ h: 0.4, b: 0.2 });
  const c = cantilever(); const p = portal();
  const tip = [{ node: 10, F: [0, 0, -1e4, 0, 0, 0] }];
  const cs = engineModel(c, tube, [0], tip);
  within('frame-static', 'cantilever tube 10 m, tip load 10 kN: tip deflection uz (B32R, PIPE)', analyzeFrame(cs).results.L.displacements[6 * 10 + 2], dat('cantilever_static', 11).displacement[32]);
  const rt = analyzeFrame(engineModel(c, rect, [0], tip)).results.L.displacements[6 * 10 + 2];
  within('frame-static', 'cantilever rectangle 0.2 × 0.4 m, tip load 10 kN: tip deflection uz (B31)', rt, dat('rect_static_b31', 11).displacement[32]);
  within('frame-static', 'cantilever rectangle 0.2 × 0.4 m, tip load 10 kN: tip deflection uz (B32)', rt, dat('rect_static_b32', 11).displacement[32]);
  // the depth h = 0.4 m of the rectangle is vertical in both codes: strong-axis value P L³ / (3 E I) = 14.88 mm, not 59.5 mm
  assert.ok(rel(-dat('rect_static_b31', 11).displacement[32], 1e4 * 1000 / (3 * 2.1e11 * 0.2 * 0.4 ** 3 / 12)) < 0.01);
  const ps = engineModel(p, tube, [0, 15], [{ node: 5, F: [1e5, 0, 0, 0, 0, 0] }]);
  within('frame-static', 'portal frame 6 × 5 m, 100 kN at the eaves: sway ux (B32R, PIPE)', analyzeFrame(ps).results.L.displacements[6 * 5], dat('portal_static', 16).displacement[15]);

  const fc = modalAnalysis(cs, { nModes: 6 }).frequencies; const fp = modalAnalysis(ps, { nModes: 6 }).frequencies;
  dat('cantilever_frequency', 11).frequencies.forEach((v, i) => within('frame-modal-cantilever', `cantilever tube 10 m, mode ${i + 1}`, fc[i], v));
  dat('portal_frequency', 16).frequencies.forEach((v, i) => within('frame-modal-portal', `portal frame 6 × 5 m, mode ${i + 1}`, fp[i], v));

  const bc = bucklingAnalysis(engineModel(c, tube, [0], [{ node: 10, F: [-1e6, 0, 0, 0, 0, 0] }]), { loadCase: 'L', nModes: 4 }).factors;
  const bp = bucklingAnalysis(engineModel(p, tube, [0, 15], [{ node: 5, F: [0, 0, -1e6, 0, 0, 0] }, { node: 10, F: [0, 0, -1e6, 0, 0, 0] }]), { loadCase: 'L', nModes: 4 }).factors;
  dat('cantilever_buckle', 11).bucklingFactors.forEach((v, i) => within('frame-buckling-cantilever', `cantilever tube 10 m, mode ${i + 1}`, bc[i], v));
  dat('portal_buckle', 16).bucklingFactors.forEach((v, i) => within('frame-buckling-portal', `portal frame 6 × 5 m, mode ${i + 1}`, bp[i], v));

  const w = -navierPlate({ a: 2, b: 1.5, t: 0.02, E: STEEL.E, nu: STEEL.nu, load: { type: 'uniform', q: 1e4 } }).w;
  const s6 = readCalculixDat(fx('calculix/plate_static.displacements.dat'), { nodeCount: 143 }).displacement[3 * 71 + 2];
  const s3 = readCalculixDat(fx('calculix/plate_static_s3.displacements.dat'), { nodeCount: 143 }).displacement[3 * 71 + 2];
  within('plate-bending', 'centre deflection, S6 shells (default of the writer)', w, s6);
  within('plate-bending', 'centre deflection, S3 shells (elementOrder 1)', w, s3);
  assert.ok(s6 < 0 && rel(w, s6) < 0.01, 'a positive pressure pushes against the element normal (+z here)');
  assert.ok(s3 / w < 0.15, 'S3 shells of CalculiX lock in bending');
});

test('decks carry what the CalculiX runs showed to be necessary', () => {
  const fix = [{ nodes: [0], dofs: 'fixed' }];
  const load = [{ type: 'force', nodes: [10], F: [0, 0, -1e4] }];
  // PIPE exists for B32R only: tubes are written as B32R with generated mid-side nodes, whatever the element order
  const t = writeCalculixCase(cantilever(), { material: STEEL, supports: fix, loads: load });
  assert.match(t.files[0].text, /\*ELEMENT, TYPE=B32R, ELSET=EALL\n1, 1, 12, 2\n/);
  assert.deepEqual(t.counts, { nodes: 21, elements: 10, appNodes: 11 });
  assert.ok(t.notes.some((n) => /B32R instead of B31/.test(n)));
  // default beam orientation = local y of the frame engine: horizontal for a member along x, so that h is vertical
  assert.match(t.files[0].text, /SECTION=PIPE\n0\.25, 0\.02\n0, 1, 0\n/);
  const r = writeCalculixCase(cantilever({ shape: 'rect', b: 0.2, h: 0.4 }), { material: STEEL, supports: fix, loads: load });
  assert.match(r.files[0].text, /\*ELEMENT, TYPE=B31, ELSET=EALL\n/);
  assert.match(r.files[0].text, /SECTION=RECT\n0\.2, 0\.4\n0, 1, 0\n/);
  const along = { kind: 'frame', nodes: Float64Array.from([0, 0, 0, 3, 4, 0]), members: [{ a: 0, b: 1, section: { shape: 'rect', b: 0.2, h: 0.4 } }] };
  assert.match(writeCalculixCase(along, { material: STEEL, supports: fix, loads: [{ type: 'force', nodes: [1], F: [0, 0, -1] }] }).files[0].text, /SECTION=RECT\n0\.2, 0\.4\n-0\.8, 0\.6, 0\n/);
  // a member's own orientation still wins
  along.members[0].orient = [0, 0, -1];
  assert.match(writeCalculixCase(along, { material: STEEL, supports: fix, loads: [{ type: 'force', nodes: [1], F: [0, 0, -1] }] }).files[0].text, /SECTION=RECT\n0\.2, 0\.4\n0, 0, -1\n/);

  // shells: S6 by default with shared mid-side nodes, supported edges held at their mid-side nodes, pressure sign reversed
  const pl = plate();
  const o = { material: STEEL, supports: [{ nodes: pl.edge, dofs: [3] }, { nodes: [0], dofs: [1, 2] }, { nodes: [12], dofs: [2] }], loads: [{ type: 'pressure', p: 1e4 }] };
  const s6 = writeCalculixCase(pl, o);
  assert.deepEqual(s6.counts, { nodes: 143 + (12 * 11 + 13 * 10 + 120), elements: 240, appNodes: 143 });
  assert.match(s6.files[0].text, /\*ELEMENT, TYPE=S6, ELSET=EALL\n1, 1, 2, 15, 144, 145, 146\n/);
  assert.match(s6.files[0].text, /\*DLOAD\nEALL, P, -10000\n/);
  const sup = /\*NSET, NSET=SUP1\n([\d, \n]+)\*NSET, NSET=SUP2\n1\n/.exec(s6.files[0].text)[1].split(/[,\s]+/).filter(Boolean).map(Number);
  assert.equal(sup.length, 44 + 44, '44 edge nodes and the 44 mid-side nodes of the boundary edges');
  const xyz = new Map(s6.files[0].text.slice(0, s6.files[0].text.indexOf('*ELEMENT')).split('\n').filter((l) => /^\d+, [-\d.e]+, [-\d.e]+, [-\d.e]+$/.test(l)).map((l) => { const v = l.split(', ').map(Number); return [v[0], v.slice(1)]; }));
  for (const id of sup) { const [x, y] = xyz.get(id); assert.ok(x === 0 || x === 2 || y === 0 || y === 1.5, `support node ${id} lies on the boundary`); }
  const s3 = writeCalculixCase(pl, { ...o, elementOrder: 1 });
  assert.match(s3.files[0].text, /\*ELEMENT, TYPE=S3, ELSET=EALL\n1, 1, 2, 15\n/);
  assert.ok(s3.notes.some((n) => /too stiff in bending/.test(n)));
});

test('exported hydrodynamic cases: NEMOH reads the mesh from the case folder; the Capytaine script survives a failing NetCDF export', () => {
  const o = { omegas: [0.25, 0.5, 0.75, 1, 1.25], headingsDeg: [0, 90, 180], cog: [0, 0, 0], name: 'barge' };
  const n = exportCase('nemoh', HULLS.barge(), o);
  assert.deepEqual(n.files.map((f) => f.path), ['Nemoh.cal', 'input_solver.txt', 'input.txt', 'ID.dat', 'barge.dat', 'results/']);
  // the Nemoh.cal that was run is the one the writer produces now
  assert.equal(n.files[0].text, fx('nemoh/Nemoh.cal'));
  assert.match(n.files[0].text, /\nbarge\.dat\t+! Name of mesh file\n413\t768\t/);
  const c = exportCase('capytaine', HULLS.barge(), o);
  const py = c.files[0].text;
  const csv = py.indexOf('with open("barge_hydro.csv", "w")'); const nc = py.indexOf('cpt.export_dataset("barge_hydro.nc", ds)');
  assert.ok(csv > 0 && nc > csv, 'the CSV is written before the NetCDF export is attempted');
  assert.match(py, /for kind, w in \(\("AINF", np\.inf\), \("A0", 0\.0\)\):/);
  assert.match(py, /except Exception as err:/);
});

test('WAMIT GDF files of the three hulls, as measured by an independent reader (Capytaine load_mesh)', () => {
  const seen = JSON.parse(fx('gdf/capytaine_gdf_reader.json'));
  for (const name of ['barge', 'cylinder', 'hemisphere']) {
    const m = wettedMesh(HULLS[name]());
    const g = exportCase('wamit-gdf', HULLS[name](), { name });
    const s = seen[name];
    assert.equal(g.mesh.nPanels, s.panels, name);
    assert.equal(Number(/^(\d+)\s+NPAN/m.exec(g.files[0].text)[1]), s.panels);
    assert.equal(m.nNodes, s.vertices, `${name}: the reader welds the repeated vertices to the node count of the app`);
    assert.ok(rel(m.area, s.area) < 1e-8 && rel(m.volume, s.volume_z) < 1e-8, name);
    // normals as the reader orients them (counter-clockwise vertex order) point into the fluid: ∫ z n_z dS = ∫ x n_x dS = V > 0
    assert.ok(s.volume_z > 0 && rel(s.volume_x, s.volume_z) < 1e-9);
    assert.ok(s.zmax === 0 && Math.abs(s.normal_sum[0]) < 1e-9 && Math.abs(s.normal_sum[1]) < 1e-9 && s.normal_sum[2] < 0);
  }
});
