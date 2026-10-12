import test from 'node:test';
import assert from 'node:assert/strict';
import { writeOpenfoamWaveTank, writeOpenfoamCurrentCase, readOpenfoamForces } from '../../src/interop/openfoam.js';
import { parseStl } from '../../src/geometry/importers/surfaces.js';
import { cylinder, sphere } from '../../src/mesh/generate.js';

const text = (r, path) => { const f = r.files.find((x) => x.path === path); assert.ok(f && typeof f.text === 'string', `file ${path}`); return f.text; };

/** Braces and parentheses balance and the FoamFile header names the object. */
function assertDictionary(path, body) {
  let b = 0; let p = 0;
  for (const ch of body) { if (ch === '{') b++; else if (ch === '}') b--; else if (ch === '(') p++; else if (ch === ')') p--; assert.ok(b >= 0 && p >= 0, `${path}: closing before opening`); }
  assert.equal(b, 0, `${path}: unbalanced braces`);
  assert.equal(p, 0, `${path}: unbalanced parentheses`);
  assert.match(body, new RegExp(`FoamFile\\s*\\{[^}]*object\\s+${path.split('/').pop().replace(/\.orig$/, '').replace('.', '\\.')};`), `${path}: FoamFile header`);
}

const PILE = cylinder({ radius: 0.2, height: 0.9, base: [0, 0, -0.55], segments: 24, heightDivisions: 6, capRings: 2, name: 'pile' });
const TANK = { surface: PILE, domain: { length: 6, width: 1.6, height: 0.7, waterDepth: 0.4 }, wave: { H: 0.06, T: 1.2, theory: 'stokes2' }, cellSize: 0.1, refinement: { structure: 1, box: 1, freeSurface: 1 }, endTime: 1.2, processors: 1 };

test('wave tank case: complete folder, consistent patches in every field, mesh and wave dictionaries', () => {
  const r = writeOpenfoamWaveTank(TANK);
  const paths = r.files.map((f) => f.path);
  for (const p of ['constant/triSurface/pile.stl', 'system/blockMeshDict', 'system/snappyHexMeshDict', 'system/controlDict', 'system/fvSchemes', 'system/fvSolution', 'system/decomposeParDict', 'system/setFieldsDict',
    'constant/transportProperties', 'constant/turbulenceProperties', 'constant/turbulenceProperties.kOmegaSST', 'constant/g', 'constant/waveProperties', '0/alpha.water', '0/alpha.water.orig', '0/U', '0/p_rgh', '0/k', '0/omega', '0/nut', 'Allrun', 'Allclean']) assert.ok(paths.includes(p), p);
  assert.equal(new Set(paths).size, paths.length);
  for (const f of r.files) if (f.text !== undefined && /^(0|constant|system)\//.test(f.path) && !f.path.endsWith('.stl')) assertDictionary(f.path.replace('turbulenceProperties.kOmegaSST', 'turbulenceProperties'), f.text);

  // tank: x from one third of the length up-wave of the pile, bottom at -waterDepth, 0.1 m cells
  assert.deepEqual(r.mesh.background, [60, 16, 7]);
  assert.deepEqual(r.mesh.bounds.min.map((v) => Number(v.toFixed(6))), [-2, -0.8, -0.4]);
  assert.deepEqual(r.mesh.bounds.max.map((v) => Number(v.toFixed(6))), [4, 0.8, 0.3]);
  const bm = text(r, 'system/blockMeshDict');
  assert.match(bm, /hex \(0 1 2 3 4 5 6 7\) \(60 16 7\) simpleGrading \(1 1 1\)/);
  assert.match(bm, /\(-2 -0\.8 -0\.4\)/);
  for (const p of ['inlet', 'outlet', 'sides', 'bottom', 'atmosphere']) assert.match(bm, new RegExp(`\\n    ${p}\\n    \\{\\n        type (patch|wall);`));
  // every field carries exactly the mesh patches plus the structure
  for (const f of ['alpha.water', 'U', 'p_rgh', 'k', 'omega', 'nut']) {
    const names = [...text(r, `0/${f}`).matchAll(/\n {4}(\w+)\n {4}\{/g)].map((m) => m[1]);
    assert.deepEqual(names.sort(), ['atmosphere', 'bottom', 'inlet', 'outlet', 'sides', 'structure'], f);
  }
  assert.match(text(r, '0/U'), /inlet\s*\{\s*type\s+waveVelocity;/);
  assert.match(text(r, '0/alpha.water'), /inlet\s*\{\s*type\s+waveAlpha;/);
  assert.match(text(r, '0/p_rgh'), /atmosphere\s*\{\s*type\s+totalPressure;/);
  assert.equal(text(r, '0/alpha.water'), text(r, '0/alpha.water.orig'));
  const wp = text(r, 'constant/waveProperties');
  assert.match(wp, /waveModel\s+StokesII;/);
  assert.match(wp, /waveHeight\s+0\.06;/);
  assert.match(wp, /wavePeriod\s+1\.2;/);
  assert.match(wp, /outlet\s*\{[^}]*waveModel\s+shallowWaterAbsorption;/);
  // water fills the tank below z = 0
  assert.match(text(r, 'system/setFieldsDict'), /box \(-3 -1\.8 -1\.4\) \(5 1\.8 0\);/);
  assert.match(text(r, 'constant/g'), /value\s+\(0 0 -9\.80665\);/);
  assert.match(text(r, 'constant/transportProperties'), /phases\s+\(water air\);/);
  assert.match(text(r, 'constant/turbulenceProperties'), /simulationType\s+laminar;/);
  assert.match(text(r, 'constant/turbulenceProperties.kOmegaSST'), /RASModel\s+kOmegaSST;/);
  const sn = text(r, 'system/snappyHexMeshDict');
  assert.match(sn, /pile\.stl\s*\{\s*type triSurfaceMesh;\s*name structure;/);
  assert.match(sn, /structure\s*\{\s*level \(1 1\);/);
  // box around the pile (radius 0.2 m, single-precision vertices) with four cells of clearance, clipped to the tank
  const rb = /refinementBox\s*\{\s*type searchableBox;\s*min \(([^)]*)\);\s*max \(([^)]*)\);/.exec(sn).slice(1).map((v) => v.split(' ').map(Number));
  [[-0.6, -0.6, -0.4], [0.6, 0.6, 0.3]].forEach((row, i) => row.forEach((v, k) => assert.ok(Math.abs(rb[i][k] - v) < 1e-6)));
  assert.match(sn, /freeSurface\s*\{\s*type searchableBox;\s*min \(-2 -0\.8 -0\.06\);\s*max \(4 0\.8 0\.06\);/);
  // the seed point of the mesher lies in the tank and outside the pile's bounding box
  const loc = /locationInMesh \(([-0-9.e]+) ([-0-9.e]+) ([-0-9.e]+)\);/.exec(sn).slice(1).map(Number);
  assert.ok(loc[0] > -2 && loc[0] < 4 && loc[1] > -0.8 && loc[1] < 0.8 && loc[2] > -0.4 && loc[2] < 0.3);
  assert.ok(Math.abs(loc[0]) > 0.25 || Math.abs(loc[1]) > 0.25);
  const cd = text(r, 'system/controlDict');
  assert.match(cd, /application\s+interFoam;/);
  assert.match(cd, /endTime\s+1\.2;/);
  assert.match(cd, /forces\s*\{\s*type\s+forces;\s*libs\s+\(forces\);\s*patches\s+\(structure\);\s*rhoInf\s+1025;\s*log/);
  assert.match(text(r, 'Allrun'), /runApplication blockMesh\nrunApplication snappyHexMesh -overwrite\nrunApplication checkMesh\ncp 0\/alpha\.water\.orig 0\/alpha\.water\nrunApplication setFields\nrunApplication \$\(getApplication\)\n/);
  assert.match(text(writeOpenfoamWaveTank({ ...TANK, processors: 4 }), 'Allrun'), /runApplication decomposePar\nrunParallel \$\(getApplication\)\nrunApplication reconstructPar\n/);
  assert.match(text(writeOpenfoamWaveTank({ ...TANK, processors: 4 }), 'system/decomposeParDict'), /numberOfSubdomains 4;/);
  assert.ok(r.notes.some((n) => /cells per wave height/.test(n)), 'coarse free-surface resolution is flagged');
});

test('wave tank: STL in both encodings is read back by the STL importer; theories, turbulence switch, current note', () => {
  const nT = PILE.indices.length / 3;
  const bin = writeOpenfoamWaveTank(TANK).files[0];
  assert.ok(bin.bytes instanceof Uint8Array && bin.bytes.length === 84 + 50 * nT);
  const b = parseStl({ name: 'pile.stl', bytes: bin.bytes }).models[0];
  assert.equal(b.indices.length / 3, nT);
  const asc = writeOpenfoamWaveTank({ ...TANK, stl: 'ascii' }).files[0];
  assert.ok(asc.text.startsWith('solid structure\n') && asc.text.trimEnd().endsWith('endsolid structure'));
  const a = parseStl({ name: 'pile.stl', text: asc.text }).models[0];
  assert.equal(a.indices.length / 3, nT);
  let zMin = Infinity; let zMax = -Infinity;
  for (let i = 2; i < a.positions.length; i += 3) { zMin = Math.min(zMin, a.positions[i] + (a.origin ? a.origin[2] : 0)); zMax = Math.max(zMax, a.positions[i] + (a.origin ? a.origin[2] : 0)); }
  assert.ok(Math.abs(zMin + 0.55) < 1e-6 && Math.abs(zMax - 0.35) < 1e-6);

  for (const [theory, model] of [['airy', 'StokesI'], ['stokes5', 'StokesV'], ['cnoidal', 'cnoidal'], [undefined, 'StokesII']]) {
    const w = text(writeOpenfoamWaveTank({ ...TANK, wave: { H: 0.06, T: 1.2, theory } }), 'constant/waveProperties');
    assert.match(w, new RegExp(`waveModel\\s+${model};`));
  }
  const k = writeOpenfoamWaveTank({ ...TANK, turbulence: 'kOmegaSST', current: 0.5 });
  assert.match(text(k, 'constant/turbulenceProperties'), /RASModel\s+kOmegaSST;/);
  assert.match(text(k, 'constant/turbulenceProperties.laminar'), /simulationType\s+laminar;/);
  assert.ok(k.notes.some((n) => /current = 0\.5 m\/s is NOT applied/.test(n)));

  assert.throws(() => writeOpenfoamWaveTank({ ...TANK, wave: { H: 0.06, T: 1.2, theory: 'jonswap' } }), /wave\.theory/);
  assert.throws(() => writeOpenfoamWaveTank({ ...TANK, wave: { H: 0.35, T: 1.2 } }), /breaking/);
  assert.throws(() => writeOpenfoamWaveTank({ ...TANK, domain: { ...TANK.domain, height: 0.3 } }), /domain\.height/);
  assert.throws(() => writeOpenfoamWaveTank({ ...TANK, domain: { ...TANK.domain, width: 0.3 } }), /structure reaches/);
  assert.throws(() => writeOpenfoamWaveTank({ ...TANK, cellSize: 0 }), /cellSize/);
  assert.throws(() => writeOpenfoamWaveTank({ ...TANK, surface: { kind: 'frame' } }), /surface/);
});

test('current case: simpleFoam set-up with inlet turbulence from intensity and length scale', () => {
  const ball = sphere({ radius: 0.5, subdivisions: 2, name: 'ball' });
  const r = writeOpenfoamCurrentCase({ surface: ball, domain: { length: 12, width: 6, height: 6 }, current: 1, cellSize: 0.5, iterations: 60, processors: 1, stl: 'ascii' });
  for (const f of r.files) if (/^(0|constant|system)\//.test(f.path) && !f.path.endsWith('.stl')) assertDictionary(f.path.replace('turbulenceProperties.laminar', 'turbulenceProperties'), f.text);
  assert.deepEqual(r.mesh.background, [24, 12, 12]);
  assert.match(text(r, 'system/controlDict'), /application\s+simpleFoam;/);
  assert.match(text(r, 'system/controlDict'), /rhoInf\s+1025;\s*rho\s+rhoInf;/);
  assert.match(text(r, '0/U'), /internalField\s+uniform \(1 0 0\);/);
  for (const f of ['U', 'p', 'k', 'omega', 'nut']) assert.deepEqual([...text(r, `0/${f}`).matchAll(/\n {4}(\w+)\n {4}\{/g)].map((m) => m[1]).sort(), ['inlet', 'outlet', 'sides', 'structure'], f);
  // k = 1.5 (I U)² = 1.5 · 0.05² = 0.00375 m²/s²; omega = sqrt(k) / (0.09^0.25 · 0.07 · D), D = 1 m (facet extent slightly less)
  assert.match(text(r, '0/k'), /internalField\s+uniform 0\.00375;/);
  const om = Number(/internalField\s+uniform ([0-9.e-]+);/.exec(text(r, '0/omega'))[1]);
  assert.ok(Math.abs(om - Math.sqrt(0.00375) / (0.09 ** 0.25 * 0.07 * 1)) / om < 0.06);
  assert.match(text(r, 'constant/transportProperties'), /transportModel\s+Newtonian;\s*nu\s+0\.00000105;/);
  assert.ok(!r.notes.some((n) => /Blockage/.test(n)));
  assert.ok(writeOpenfoamCurrentCase({ surface: ball, domain: { length: 12, width: 2, height: 2 }, current: 1, cellSize: 0.25 }).notes.some((n) => /Blockage/.test(n)));
  assert.throws(() => writeOpenfoamCurrentCase({ surface: ball, domain: { length: 12, width: 0.8, height: 6 }, current: 1, cellSize: 0.5 }), /strictly inside/);
  assert.throws(() => writeOpenfoamCurrentCase({ surface: ball, domain: { length: 12, width: 6, height: 6 }, current: 0, cellSize: 0.5 }), /current/);
});

// Lines taken from postProcessing/forces/0/ of an interFoam run of the tank case above with OpenFOAM v2406.
const FORCE_V2406 = `# Force
# CofR        : (0.000000e+00 0.000000e+00 0.000000e+00)
#
# Time        	total_x total_y total_z	pressure_x pressure_y pressure_z	viscous_x viscous_y viscous_z
0.00363636     -1.195098e-03 -3.552223e-05 -1.864935e+00 -1.195096e-03 -3.552191e-05 -1.864935e+00 -1.767730e-09 -3.236246e-10 8.109923e-08
0.00794613     -1.368761e-03 -3.287870e-05 -1.864935e+00 -1.368757e-03 -3.287817e-05 -1.864935e+00 -3.892857e-09 -5.319790e-10 1.852904e-07
1.176          -8.319950e-01 2.194867e-04 -1.860595e+00 -8.318737e-01 2.194893e-04 -1.860580e+00 -1.213385e-04 -2.554533e-09 -1.470987e-05
1.2            -6.987001e-01 1.576325e-04 -1.860707e+00 -6.985764e-01 1.576306e-04 -1.860699e+00 -1.236997e-04 1.969143e-09 -8.240153e-06
`;
const MOMENT_V2406 = `# Moment
# CofR        : (0.000000e+00 0.000000e+00 0.000000e+00)
#
# Time        	total_x total_y total_z	pressure_x pressure_y pressure_z	viscous_x viscous_y viscous_z
0.00363636     1.166993e-06 3.377585e-04 7.672442e-08 1.166974e-06 3.377583e-04 7.738504e-08 1.871091e-11 1.592087e-10 -6.606209e-10
1.176          3.776538e-05 1.458043e-01 1.274876e-07 3.776495e-05 1.457772e-01 8.612997e-08 4.329190e-10 2.713277e-05 4.135762e-08
1.2            3.083060e-05 1.144900e-01 1.972665e-07 3.082989e-05 1.144616e-01 1.537265e-07 7.038016e-10 2.840684e-05 4.354002e-08
`;
// Layout of forces.dat written by the forces function object up to OpenFOAM v2012 and by the Foundation releases.
const FORCES_BRACKETED = `# Forces
# CofR       : (1.000000e+00 0.000000e+00 -5.000000e-01)
# Time       forces(pressure viscous porous) moment(pressure viscous porous)
0.05	((1.2e+03 -3.0e+00 4.5e+02) (2.5e+01 1.0e-01 -2.0e+00) (0 0 0)) ((1.0e+01 -6.0e+02 3.0e+00) (5.0e-01 -1.5e+01 2.0e-01) (0 0 0))
0.1	((1.4e+03 -2.0e+00 4.0e+02) (2.7e+01 2.0e-01 -1.0e+00) (0 0 0)) ((1.2e+01 -7.0e+02 2.0e+00) (6.0e-01 -1.7e+01 1.0e-01) (0 0 0))
0.1	((1.5e+03 -2.0e+00 4.0e+02) (2.7e+01 2.0e-01 -1.0e+00) (0 0 0)) ((1.2e+01 -7.0e+02 2.0e+00) (6.0e-01 -1.7e+01 1.0e-01) (0 0 0))
`;

test('readOpenfoamForces: column layout (v2106+, force.dat + moment.dat) and bracketed layout (forces.dat)', () => {
  const r = readOpenfoamForces(FORCE_V2406, { momentText: MOMENT_V2406 });
  assert.equal(r.format, 'columns');
  assert.equal(r.n, 4);
  assert.deepEqual(r.cofr, [0, 0, 0]);
  assert.deepEqual(r.t, [0.00363636, 0.00794613, 1.176, 1.2]);
  assert.deepEqual([r.Fx[3], r.Fy[3], r.Fz[3]], [-0.6987001, 0.0001576325, -1.860707]);
  assert.deepEqual([r.pressure.Fx[2], r.viscous.Fx[2]], [-0.8318737, -0.0001213385]);
  // total = pressure + viscous to the precision of the file
  for (let i = 0; i < r.n; i++) assert.ok(Math.abs(r.Fx[i] - r.pressure.Fx[i] - r.viscous.Fx[i]) < 1e-6);
  // moment.dat has no record at t = 0.00794613: the gap is NaN, the rest is aligned on time
  assert.deepEqual([r.My[0], r.My[2], r.My[3]], [3.377585e-4, 0.1458043, 0.11449]);
  assert.ok(Number.isNaN(r.My[1]));
  const only = readOpenfoamForces(FORCE_V2406);
  assert.equal(only.Mx, null);

  const b = readOpenfoamForces(FORCES_BRACKETED);
  assert.equal(b.format, 'bracketed');
  assert.deepEqual(b.cofr, [1, 0, -0.5]);
  // a restart repeated t = 0.1: the later record is kept
  assert.deepEqual(b.t, [0.05, 0.1]);
  assert.deepEqual([b.Fx[0], b.Fy[0], b.Fz[0]], [1225, -2.9, 448]);
  assert.deepEqual([b.Mx[0], b.My[0], b.Mz[0]], [10.5, -615, 3.2]);
  assert.equal(b.Fx[1], 1527);
  assert.deepEqual([b.pressure.Fx[1], b.viscous.Fz[1]], [1500, -1]);
  // older files without the porous group: 2 × (pressure, viscous)
  const old = readOpenfoamForces('# Time forces(pressure viscous) moment(pressure viscous)\n2 ((10 0 1) (2 0 0)) ((0 5 0) (0 1 0))\n');
  assert.deepEqual([old.Fx[0], old.Fz[0], old.My[0]], [12, 1, 6]);
  assert.throws(() => readOpenfoamForces('# Force\n'), /no data rows/);
  assert.throws(() => readOpenfoamForces('1 (1 2 3) (4 5 6)\n'), /not a known layout/);
});
