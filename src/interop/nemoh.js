// Hand-off to frequency-domain potential-flow (radiation/diffraction) solvers and back.
//
//   writeNemohCase        NEMOH v3 case folder (Nemoh.cal, input_solver.txt, mesh/<name>.dat, results/ …)
//   readNemohResults      NEMOH results (RadiationCoefficients.tec or CA.dat/CM.dat, ExcitationForce.tec) → hydro database
//   writeWamitGdf         WAMIT low-order geometric data file (.gdf) of the same wetted mesh
//   writeCapytaineScript  Python script for Capytaine 2.x that solves the exported mesh and writes a CSV
//   readCapytaineCsv      that CSV → hydro database
//
// Targets: NEMOH v3.0 (Kurnia & Ducrozet, "NEMOH v3.0 User Manual", 2022; the load-case line of Nemoh.cal carries
// the frequency-type flag introduced in v3), WAMIT v7 low-order GDF (User Manual ch. 6.1), Capytaine 2.x
// (Ancellin & Dias, J. Open Source Softw. 4 (2019)).
//
// The wetted mesh is the hull clipped at z = 0 by wettedPanels() of src/engine/hydro/bem.js: flat triangles,
// vertices counter-clockwise seen from the fluid (normals into the fluid), written as quadrilaterals with a
// repeated vertex. No symmetry plane is exploited (symmetry flag 0).
//
// Conventions of the hydro database (src/engine/hydro/hydrodb.js): SI, exp(+i omega t), heading = propagation
// direction from +x. NEMOH and Capytaine use exp(-i omega t); the readers conjugate the excitation by default.

import { wettedPanels } from '../engine/hydro/bem.js';
import { createHydroDb, parseHydroCsv } from '../engine/hydro/hydrodb.js';

const G0 = 9.80665;
const num = (v) => (Number.isFinite(v) ? String(Number(v.toPrecision(12))) : 'NaN');
const fx = (v, w = 16) => v.toExponential(8).padStart(w);
const safeName = (s) => String(s || 'hull').replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/\.[A-Za-z0-9]+$/, '') || 'hull';

function needPos(v, name) { if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name} must be > 0`); return v; }
function needVec3(v, name) { if (!v || v.length !== 3 || ![0, 1, 2].every((k) => Number.isFinite(v[k]))) throw new RangeError(`${name} must be [x, y, z] in metres`); return v; }

/**
 * Wetted panel mesh of a hull with shared vertices: the hull is clipped at z = 0, coincident vertices are welded.
 * @param {object} surface closed (or already open-at-the-waterline) triangle surface, z up, waterline z = 0, m
 * @param {{weldTol?:number}} [opts] welding tolerance, m (default 1e-7 × largest coordinate magnitude, min 1e-9)
 * @returns {{nodes:Float64Array, panels:Uint32Array, nNodes:number, nPanels:number, normals:Float64Array,
 *   centroids:Float64Array, areas:Float64Array, area:number, volume:number, clipped:number, flipped:boolean}}
 *   nodes: 3 per node (m); panels: 3 node indices (0-based) per triangular panel, counter-clockwise seen from the
 *   fluid; unit normals into the fluid; areas m²; wetted area m²; displaced volume m³
 */
export function wettedMesh(surface, { weldTol } = {}) {
  const w = wettedPanels(surface);
  let big = 0;
  for (let i = 0; i < w.vertices.length; i++) big = Math.max(big, Math.abs(w.vertices[i]));
  const tol = weldTol === undefined ? Math.max(1e-9, 1e-7 * big) : weldTol;
  const map = new Map(); const nodes = []; const panels = new Uint32Array(3 * w.n);
  for (let k = 0; k < 3 * w.n; k++) {
    const x = w.vertices[3 * k]; const y = w.vertices[3 * k + 1]; const z = w.vertices[3 * k + 2];
    const key = `${Math.round(x / tol)},${Math.round(y / tol)},${Math.round(z / tol)}`;
    let id = map.get(key);
    if (id === undefined) { id = nodes.length / 3; nodes.push(x, y, z); map.set(key, id); }
    panels[k] = id;
  }
  return { nodes: Float64Array.from(nodes), panels, nNodes: nodes.length / 3, nPanels: w.n, normals: w.normals, centroids: w.centroids, areas: w.areas, area: w.area, volume: w.volume, clipped: w.clipped, flipped: w.flipped };
}

/** NEMOH mesh file ("2 isym", node table, "0 0. 0. 0.", panel table, "0 0 0 0"). */
function nemohMeshText(mesh) {
  const L = ['    2    0'];
  for (let i = 0; i < mesh.nNodes; i++) L.push(`${String(i + 1).padStart(8)} ${fx(mesh.nodes[3 * i])} ${fx(mesh.nodes[3 * i + 1])} ${fx(mesh.nodes[3 * i + 2])}`);
  L.push('       0   0.00000000e+0   0.00000000e+0   0.00000000e+0');
  for (let p = 0; p < mesh.nPanels; p++) {
    const a = mesh.panels[3 * p] + 1; const b = mesh.panels[3 * p + 1] + 1; const c = mesh.panels[3 * p + 2] + 1;
    L.push(`${String(a).padStart(8)}${String(b).padStart(8)}${String(c).padStart(8)}${String(a).padStart(8)}`);
  }
  L.push('       0       0       0       0');
  return L.join('\n') + '\n';
}

/** Uniform range (count, min, max) covering a list of values; `uniform` tells whether the list is reproduced. */
function uniformRange(values, name) {
  const v = Array.from(values || []);
  if (v.length === 0 || v.some((x) => !Number.isFinite(x))) throw new RangeError(`${name} must be a non-empty list of finite numbers`);
  const s = v.slice().sort((a, b) => a - b);
  const min = s[0]; const max = s[s.length - 1];
  let uniform = true;
  const step = s.length > 1 ? (max - min) / (s.length - 1) : 0;
  for (let i = 0; i < s.length; i++) if (Math.abs(s[i] - (min + i * step)) > 1e-6 * Math.max(Math.abs(max), Math.abs(step), 1e-12)) uniform = false;
  return { n: s.length, min, max, uniform };
}

/**
 * Writes a NEMOH v3.0 case: one rigid body, six degrees of freedom and six generalised forces about `cog`.
 *
 * Folder: Nemoh.cal, input_solver.txt (v3 solver settings), mesh/<name>.dat (solver-format mesh), results/,
 * plus ID.dat and input.txt, which only NEMOH v2 reads (v3 takes the case folder on the command line).
 * Run (v3):  preProc <case> ; hydrosCal <case> (optional) ; solver <case> ; postProc <case>.
 * NEMOH takes a uniform frequency range and a uniform heading range (count, min, max): lists that are not
 * uniformly spaced are replaced by the uniform range with the same count and limits, and a note says so.
 * Reference: NEMOH v3.0 User Manual (2022), "Description of input files"; Babarit & Delhommeau, EWTEC 2015.
 * Validity: consistently oriented hull surface, z up, still-water level z = 0; body free of the seabed.
 * @param {object} o
 * @param {object} o.surface closed hull surface (CONTRACT §4 `surface`), m
 * @param {number} [o.rho=1025] water density, kg/m³
 * @param {number} [o.g=9.80665] gravity, m/s²
 * @param {number|null} [o.depth=null] water depth, m (null or Infinity: deep water, written as 0)
 * @param {ArrayLike<number>} o.omegas wave frequencies, rad/s (> 0)
 * @param {ArrayLike<number>} [o.headingsDeg=[0]] wave propagation directions from +x, deg
 * @param {number[]} [o.cog=[0,0,0]] point about which rotations and moments are taken, m
 * @param {string} [o.name='hull'] body / mesh file name
 * @returns {{files:{path:string, text:string}[], notes:string[], mesh:{nNodes:number, nPanels:number,
 *   area:number, volume:number}}} files of the case folder (a path ending in '/' is an empty folder), remarks
 *   for the user, mesh counts with wetted area m² and displaced volume m³
 */
export function writeNemohCase({ surface, rho = 1025, g = G0, depth = null, omegas, headingsDeg = [0], cog = [0, 0, 0], name = 'hull' } = {}) {
  needPos(rho, 'rho'); needPos(g, 'g'); needVec3(cog, 'cog');
  const deep = depth === null || depth === undefined || depth === Infinity;
  if (!deep) needPos(depth, 'depth');
  const fr = uniformRange(omegas, 'omegas');
  if (!(fr.min > 0)) throw new RangeError('omegas must be > 0');
  const hr = uniformRange(headingsDeg, 'headingsDeg');
  const mesh = wettedMesh(surface);
  const nm = safeName(name);
  const notes = [];
  if (!fr.uniform) notes.push(`omegas are not uniformly spaced: NEMOH will compute ${fr.n} frequencies uniformly from ${num(fr.min)} to ${num(fr.max)} rad/s.`);
  if (!hr.uniform) notes.push(`headingsDeg are not uniformly spaced: NEMOH will compute ${hr.n} headings uniformly from ${num(hr.min)} to ${num(hr.max)} deg.`);
  if (mesh.flipped) notes.push('The surface was oriented inwards; panel normals were reversed to point into the fluid.');
  if (mesh.clipped) notes.push(`${mesh.clipped} triangle(s) were cut at the waterplane z = 0.`);
  notes.push(`${mesh.nPanels} triangular panels (written as quadrilaterals with a repeated node), ${mesh.nNodes} nodes, wetted area ${num(Number(mesh.area.toPrecision(6)))} m², displaced volume ${num(Number(mesh.volume.toPrecision(6)))} m³.`);
  notes.push('No symmetry is used (flag 0). Irregular-frequency removal, QTF and RAO post-processing are switched off; hydrostatics and the mass matrix stay in the application.');
  notes.push('NEMOH works with exp(-i omega t); readNemohResults converts the excitation to the exp(+i omega t) convention of the application.');
  let zmin = 0;
  for (let i = 2; i < mesh.nodes.length; i += 3) zmin = Math.min(zmin, mesh.nodes[i]);
  if (!deep && !(depth > -zmin)) throw new RangeError('depth must exceed the draught of the hull');
  const c = cog.map(num).join(' ');
  const cal = [
    '--- Environment ------------------------------------------------------------------------------------------------------------------',
    `${num(rho)}\t\t\t\t! RHO \t\t\t! KG/M**3 \t! Fluid specific volume`,
    `${num(g)}\t\t\t\t! G\t\t\t! M/S**2\t! Gravity`,
    `${deep ? '0.' : num(depth)}\t\t\t\t! DEPTH\t\t\t! M\t\t! Water depth (0 = infinite)`,
    '0.\t0.\t\t\t\t! XEFF YEFF\t\t! M\t\t! Wave measurement point',
    '--- Description of floating bodies -----------------------------------------------------------------------------------------------',
    '1\t\t\t\t! Number of bodies',
    '--- Body 1 -----------------------------------------------------------------------------------------------------------------------',
    `${nm}.dat\t\t\t! Name of mesh file`,
    `${mesh.nNodes}\t${mesh.nPanels}\t\t\t! Number of points and number of panels`,
    '6\t\t\t\t! Number of degrees of freedom',
    '1 1. 0. 0. 0. 0. 0.\t\t! Surge',
    '1 0. 1. 0. 0. 0. 0.\t\t! Sway',
    '1 0. 0. 1. 0. 0. 0.\t\t! Heave',
    `2 1. 0. 0. ${c}\t\t! Roll about a point`,
    `2 0. 1. 0. ${c}\t\t! Pitch about a point`,
    `2 0. 0. 1. ${c}\t\t! Yaw about a point`,
    '6\t\t\t\t! Number of resulting generalised forces',
    '1 1. 0. 0. 0. 0. 0.\t\t! Force in x direction',
    '1 0. 1. 0. 0. 0. 0.\t\t! Force in y direction',
    '1 0. 0. 1. 0. 0. 0.\t\t! Force in z direction',
    `2 1. 0. 0. ${c}\t\t! Moment force in x direction about a point`,
    `2 0. 1. 0. ${c}\t\t! Moment force in y direction about a point`,
    `2 0. 0. 1. ${c}\t\t! Moment force in z direction about a point`,
    '0\t\t\t\t! Number of lines of additional information',
    '--- Load cases to be solved -------------------------------------------------------------------------------------------------------',
    `1\t${fr.n}\t${num(fr.min)}\t${num(fr.max)}\t\t! Freq type 1,2,3=[rad/s,Hz,s], Number of wave frequencies/periods, Min, and Max`,
    `${hr.n}\t${num(hr.min)}\t${num(hr.max)}\t\t! Number of wave directions, Min and Max (degrees)`,
    '--- Post processing ---------------------------------------------------------------------------------------------------------------',
    '0\t0.1\t10.\t\t! IRF \t\t\t\t! IRF calculation (0 for no calculation), time step and duration',
    '0\t\t\t\t! Show pressure',
    '0\t0.\t180.\t\t! Kochin function \t\t! Number of directions of calculation (0 for no calculations), Min and Max (degrees)',
    '0\t50\t400.\t400.\t! Free surface elevation \t! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction',
    '0\t\t\t\t! Response Amplitude Operator (RAO), 0 no calculation, 1 calculated',
    '1\t\t\t\t! output freq type, 1,2,3=[rad/s,Hz,s]',
    '---QTF---',
    '0\t\t\t\t! QTF flag, 1 selected',
    '',
  ].join('\n');
  const solver = [
    '2\t\t\t\t! Gauss quadrature (GQ) surface integration, N^2 GQ Nodes, specify N(1,4)',
    '0.001\t\t\t\t! eps_zmin for determine minimum z of flow and source points of panel, zmin=eps_zmin*body_diameter',
    '1\t\t\t\t! 0 GAUSS ELIM.; 1 LU DECOMP.: 2 GMRES\t!Linear system solver',
    '10 1e-5 1000\t\t\t! Restart parameter, Relative Tolerance, max iter -> additional input for GMRES',
    '',
  ].join('\n');
  const v2input = [
    '--- Calculation parameters ------------------------------------------------------------------------------------',
    '1\t\t\t\t! Indiq_solver\t\t! - \t\t! Solver (0) Direct Gauss (1) GMRES (2) GMRES with FMM acceleration (2 not implemented yet)',
    '20\t\t\t\t! IRES\t\t\t! - \t\t! Restart parameter for GMRES',
    '5.E-07\t\t\t\t! TOL_GMRES\t\t! -\t\t! Stopping criterion for GMRES',
    '100\t\t\t\t! MAXIT\t\t\t! - \t\t! Maximum iterations for GMRES',
    '1\t\t\t\t! Sav_sol\t\t! -\t\t! Save solution for visualisation (0) No (1) Yes',
    '',
  ].join('\n');
  const files = [
    { path: 'Nemoh.cal', text: cal },
    { path: 'input_solver.txt', text: solver },
    { path: 'input.txt', text: v2input },
    { path: 'ID.dat', text: '1\n.\n' },
    { path: `mesh/${nm}.dat`, text: nemohMeshText(mesh) },
    { path: 'results/', text: '' },
  ];
  return { files, notes, mesh: { nNodes: mesh.nNodes, nPanels: mesh.nPanels, area: mesh.area, volume: mesh.volume } };
}

/**
 * WAMIT low-order geometric data file of the wetted hull: header, "ULEN GRAV", "ISX ISY", "NPAN", then four
 * vertices per panel (triangles repeat the last vertex), counter-clockwise when seen from the fluid.
 * Reference: WAMIT v7 User Manual, ch. 6.1 (low-order GDF).
 * Validity: as writeNemohCase; lengths are written in metres (GRAV = g in m/s²).
 * @param {object} o
 * @param {object} o.surface closed hull surface, m, z up, waterline z = 0
 * @param {number} [o.ulen=1] WAMIT characteristic length ULEN, m
 * @param {number} [o.g=9.80665] gravity, m/s²
 * @param {string} [o.name='hull'] name used in the header line and the file name
 * @returns {{files:{path:string, text:string}[], notes:string[], mesh:{nNodes:number, nPanels:number,
 *   area:number, volume:number}}}
 */
export function writeWamitGdf({ surface, ulen = 1, g = G0, name = 'hull' } = {}) {
  needPos(ulen, 'ulen'); needPos(g, 'g');
  const mesh = wettedMesh(surface);
  const nm = safeName(name);
  const L = [`${nm} - low-order GDF, wetted surface, metres`, `${num(ulen)} ${num(g)}  ULEN GRAV`, '0 0  ISX ISY', `${mesh.nPanels}  NPAN`];
  const v = (i) => `${fx(mesh.nodes[3 * i])} ${fx(mesh.nodes[3 * i + 1])} ${fx(mesh.nodes[3 * i + 2])}`;
  for (let p = 0; p < mesh.nPanels; p++) {
    const a = mesh.panels[3 * p]; const b = mesh.panels[3 * p + 1]; const c = mesh.panels[3 * p + 2];
    L.push(v(a), v(b), v(c), v(c));
  }
  const notes = [`${mesh.nPanels} panels, no symmetry planes (ISX = ISY = 0); normals point into the fluid.`];
  if (mesh.flipped) notes.push('The surface was oriented inwards; panel normals were reversed to point into the fluid.');
  notes.push('WAMIT output (.1, .3) is read back with hydroDbFromWamit of the hydro database module, using the same ULEN.');
  return { files: [{ path: `${nm}.gdf`, text: L.join('\n') + '\n' }], notes, mesh: { nNodes: mesh.nNodes, nPanels: mesh.nPanels, area: mesh.area, volume: mesh.volume } };
}

function fileMap(files) {
  const m = new Map();
  if (Array.isArray(files)) for (const f of files) m.set(String(f.path || f.name || ''), f.text ?? (f.bytes ? new TextDecoder().decode(f.bytes) : ''));
  else if (files && typeof files === 'object') for (const k of Object.keys(files)) m.set(k, typeof files[k] === 'string' ? files[k] : new TextDecoder().decode(files[k]));
  return m;
}
function findText(map, base) {
  const want = base.toLowerCase();
  for (const [k, v] of map) { const b = k.replace(/\\/g, '/').split('/').pop().toLowerCase(); if (b === want) return v; }
  return null;
}
const numbersOf = (line) => line.trim().split(/[\s,]+/).filter(Boolean).map(Number);

/** Zones of a Tecplot POINT file: [{ title, rows: number[][] }], tolerant of wrapped records. */
function tecZones(text, what) {
  const zones = [];
  let cur = null;
  for (const line of String(text).split(/\r?\n/)) {
    const zm = /^\s*zone\b(.*)$/i.exec(line);
    if (zm) {
      const t = /t\s*=\s*"([^"]*)"/i.exec(zm[1]); const im = /\bI\s*=\s*(\d+)/i.exec(zm[1]);
      cur = { title: t ? t[1] : '', count: im ? Number(im[1]) : 0, values: [] };
      zones.push(cur);
      continue;
    }
    if (!cur) continue;
    const v = numbersOf(line);
    if (v.length && v.every((x) => Number.isFinite(x))) for (const x of v) cur.values.push(x);
  }
  if (!zones.length) throw new RangeError(`${what}: no "Zone" records found`);
  for (const z of zones) {
    const n = z.count || 0;
    if (!n || z.values.length % n !== 0) throw new RangeError(`${what}: zone "${z.title}" does not hold I = ${n} complete rows`);
    const nc = z.values.length / n;
    z.rows = [];
    for (let i = 0; i < n; i++) z.rows.push(z.values.slice(i * nc, (i + 1) * nc));
  }
  return zones;
}

/**
 * Parses NEMOH `RadiationCoefficients.tec`: one zone per radiating degree of freedom j, rows
 * `omega  A_1j B_1j  A_2j B_2j … A_6j B_6j` (dimensional, SI).
 * Reference: NEMOH v3.0 User Manual, "Output files".
 * @param {string} text file content
 * @returns {{omega:number[], A:number[][][], B:number[][][]}} rad/s; [frequency][6][6] in kg … / N·s/m …
 */
export function parseNemohRadiation(text) {
  const zones = tecZones(text, 'RadiationCoefficients.tec');
  if (zones.length !== 6) throw new RangeError(`RadiationCoefficients.tec: expected 6 zones (one per degree of freedom of a single body), found ${zones.length}`);
  const omega = zones[0].rows.map((r) => r[0]);
  const A = omega.map(() => Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]));
  const B = omega.map(() => Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]));
  zones.forEach((z, j) => {
    if (z.rows.length !== omega.length || z.rows[0].length !== 13) throw new RangeError('RadiationCoefficients.tec: each zone must hold omega and six (A, B) pairs per row');
    z.rows.forEach((r, w) => { for (let i = 0; i < 6; i++) { A[w][i][j] = r[1 + 2 * i]; B[w][i][j] = r[2 + 2 * i]; } });
  });
  return { omega, A, B };
}

/**
 * Parses NEMOH `CA.dat` (radiation damping) or `CM.dat` (added mass): a header line with the number of
 * frequencies, then per frequency one line with the frequency (or period) followed by a 6×6 block.
 * @param {string} text file content
 * @param {{abscissa?:'auto'|'omega'|'period'}} [opts] meaning of the per-block value; 'auto' reads it from the
 *   header ("periods" → s, otherwise rad/s)
 * @returns {{omega:number[], M:number[][][]}} rad/s; [frequency][6][6], SI
 */
export function parseNemohMatrixDat(text, { abscissa = 'auto' } = {}) {
  const lines = String(text).split(/\r?\n/);
  const head = lines.shift() || '';
  const count = Number((/(\d+)\s*$/.exec(head.trim()) || [])[1]);
  const vals = [];
  for (const l of lines) for (const v of numbersOf(l)) if (Number.isFinite(v)) vals.push(v);
  if (!(count > 0) || vals.length < count * 37) throw new RangeError('NEMOH CA/CM file: expected a header with the number of frequencies and 37 numbers per frequency (single body, six modes)');
  const isPeriod = abscissa === 'period' || (abscissa === 'auto' && /period/i.test(head));
  const omega = []; const M = [];
  for (let k = 0; k < count; k++) {
    const o = 37 * k;
    omega.push(isPeriod ? (2 * Math.PI) / vals[o] : vals[o]);
    M.push(Array.from({ length: 6 }, (_, i) => vals.slice(o + 1 + 6 * i, o + 7 + 6 * i)));
  }
  return { omega, M };
}

/**
 * Parses NEMOH `ExcitationForce.tec`: one zone per wave direction ("… beta = <deg> deg"), rows
 * `omega  |F_1| arg(F_1) … |F_6| arg(F_6)` per unit wave amplitude.
 * @param {string} text file content
 * @param {{phaseUnit?:'auto'|'rad'|'deg', conjugate?:boolean}} [opts] unit of the phase column ('auto': degrees
 *   when any |phase| exceeds 2π + 0.01, radians otherwise); conjugate = true converts from NEMOH's
 *   exp(-i omega t) to exp(+i omega t)
 * @returns {{omega:number[], headingsDeg:number[], Fre:number[][][], Fim:number[][][], phaseUnit:string}}
 *   [heading][frequency][6], N/m and N·m/m
 */
export function parseNemohExcitation(text, { phaseUnit = 'auto', conjugate = true } = {}) {
  const zones = tecZones(text, 'ExcitationForce.tec');
  let unit = phaseUnit;
  if (unit === 'auto') {
    let big = 0;
    for (const z of zones) for (const r of z.rows) for (let i = 2; i < r.length; i += 2) big = Math.max(big, Math.abs(r[i]));
    unit = big > 2 * Math.PI + 0.01 ? 'deg' : 'rad';
  }
  if (unit !== 'rad' && unit !== 'deg') throw new RangeError("phaseUnit must be 'auto', 'rad' or 'deg'");
  const k = unit === 'deg' ? Math.PI / 180 : 1;
  const headingsDeg = []; const Fre = []; const Fim = [];
  const omega = zones[0].rows.map((r) => r[0]);
  zones.forEach((z, n) => {
    const m = /beta\s*=\s*([-+0-9.eE]+)\s*(rad|deg)?/i.exec(z.title);
    let h = m ? Number(m[1]) : n;
    if (m && m[2] && m[2].toLowerCase() === 'rad') h = (h * 180) / Math.PI;
    headingsDeg.push(h);
    if (z.rows.length !== omega.length || z.rows[0].length !== 13) throw new RangeError('ExcitationForce.tec: each zone must hold omega and six (modulus, phase) pairs per row');
    Fre.push(z.rows.map((r) => [0, 1, 2, 3, 4, 5].map((i) => r[1 + 2 * i] * Math.cos(k * r[2 + 2 * i]))));
    Fim.push(z.rows.map((r) => [0, 1, 2, 3, 4, 5].map((i) => (conjugate ? -1 : 1) * r[1 + 2 * i] * Math.sin(k * r[2 + 2 * i]))));
  });
  return { omega, headingsDeg, Fre, Fim, phaseUnit: unit };
}

/**
 * Builds a hydrodynamic database (createHydroDb of hydrodb.js) from the result files of a NEMOH run, ready for
 * raoSolve and the Cummins time-domain model. Radiation coefficients come from RadiationCoefficients.tec when
 * present, otherwise from CM.dat (added mass) and CA.dat (damping); excitation from ExcitationForce.tec.
 * rho, g and depth are read from Nemoh.cal when it is among the files, unless given.
 * NEMOH does not report the infinite-frequency added mass in these files: pass `Ainf` (for instance from
 * addedMassBem with limit 'infinite') when the time-domain model needs it.
 * Reference: NEMOH v3.0 User Manual, "Output files".
 * @param {object} o
 * @param {{path:string, text:string}[]|Object<string,string>} o.files result files (any folder prefix)
 * @param {number} [o.rho] kg/m³
 * @param {number} [o.g] m/s²
 * @param {number|null} [o.depth] m (null: deep)
 * @param {'auto'|'rad'|'deg'} [o.phaseUnit='auto'] unit of the excitation phase column
 * @param {boolean} [o.conjugate=true] convert the excitation to exp(+i omega t)
 * @param {'auto'|'omega'|'period'} [o.abscissa='auto'] abscissa of CA.dat / CM.dat
 * @param {number[][]} [o.Ainf] infinite-frequency added mass, 6×6 SI
 * @param {number[][]} [o.A0] zero-frequency added mass, 6×6 SI
 * @returns {object} database (createHydroDb) with an extra `source` string naming the files used
 */
export function readNemohResults({ files, rho, g, depth, phaseUnit = 'auto', conjugate = true, abscissa = 'auto', Ainf = null, A0 = null } = {}) {
  const map = fileMap(files);
  if (map.size === 0) throw new RangeError('files must hold the NEMOH result files');
  const cal = findText(map, 'Nemoh.cal');
  let env = {};
  if (cal) {
    const l = cal.split(/\r?\n/);
    const first = (s) => Number((s || '').trim().split(/[\s!]+/)[0]);
    env = { rho: first(l[1]), g: first(l[2]), depth: first(l[3]) };
  }
  const rho1 = rho ?? (env.rho > 0 ? env.rho : 1025);
  const g1 = g ?? (env.g > 0 ? env.g : G0);
  const depth1 = depth !== undefined ? depth : (env.depth > 0 ? env.depth : null);
  const used = [];
  let rad;
  const tec = findText(map, 'RadiationCoefficients.tec');
  if (tec) { rad = parseNemohRadiation(tec); used.push('RadiationCoefficients.tec'); } else {
    const cm = findText(map, 'CM.dat'); const ca = findText(map, 'CA.dat');
    if (!cm || !ca) throw new RangeError('files: RadiationCoefficients.tec, or both CM.dat and CA.dat, are required');
    const a = parseNemohMatrixDat(cm, { abscissa }); const b = parseNemohMatrixDat(ca, { abscissa });
    if (a.omega.length !== b.omega.length) throw new RangeError('CM.dat and CA.dat list different numbers of frequencies');
    rad = { omega: a.omega, A: a.M, B: b.M };
    used.push('CM.dat', 'CA.dat');
  }
  let exc = { headingsDeg: [], Fre: [], Fim: [] };
  const ft = findText(map, 'ExcitationForce.tec');
  if (ft) {
    exc = parseNemohExcitation(ft, { phaseUnit, conjugate });
    if (exc.omega.length !== rad.omega.length || exc.omega.some((w, i) => Math.abs(w - rad.omega[i]) > 1e-5 * Math.max(w, 1e-12))) throw new RangeError('ExcitationForce.tec and the radiation coefficients list different frequencies');
    used.push('ExcitationForce.tec');
  }
  const db = createHydroDb({ omega: rad.omega, A: rad.A, B: rad.B, headingsDeg: exc.headingsDeg, Fre: exc.Fre, Fim: exc.Fim, Ainf, A0, rho: rho1, g: g1, depth: depth1 });
  db.source = `NEMOH (${used.join(', ')})`;
  return db;
}

/**
 * Python script for Capytaine 2.x: loads the exported NEMOH-format mesh, solves radiation and diffraction for
 * the given frequencies and headings, saves the full dataset as NetCDF and writes `<name>_hydro.csv` in the CSV
 * schema of hydrodb.js (`kind,omega,headingDeg,i,j,re,im`), already converted to exp(+i omega t), together with
 * the infinite-frequency added mass (rows `AINF`).
 * Run:  pip install capytaine ; python run_capytaine.py   (in the folder holding the script and mesh/).
 * Reference: Capytaine documentation v2 (user manual: "Resolution", "Outputs"); Ancellin & Dias (2019).
 * Validity: as writeNemohCase. Capytaine removes nothing above z = 0 here: the mesh is already the wetted part.
 * @param {object} o same arguments as writeNemohCase (`surface`, `rho`, `g`, `depth`, `omegas`, `headingsDeg`,
 *   `cog`, `name`); frequencies and headings are used as listed (no uniform spacing needed)
 * @returns {{files:{path:string, text:string}[], notes:string[], mesh:{nNodes:number, nPanels:number,
 *   area:number, volume:number}}} run_capytaine.py and mesh/<name>.dat
 */
export function writeCapytaineScript({ surface, rho = 1025, g = G0, depth = null, omegas, headingsDeg = [0], cog = [0, 0, 0], name = 'hull' } = {}) {
  needPos(rho, 'rho'); needPos(g, 'g'); needVec3(cog, 'cog');
  const deep = depth === null || depth === undefined || depth === Infinity;
  if (!deep) needPos(depth, 'depth');
  const om = Array.from(omegas || []); const hd = Array.from(headingsDeg || []);
  if (!om.length || om.some((w) => !(w > 0) || !Number.isFinite(w))) throw new RangeError('omegas must be a non-empty list of frequencies > 0');
  if (!hd.length || hd.some((h) => !Number.isFinite(h))) throw new RangeError('headingsDeg must be a non-empty list of finite numbers');
  const mesh = wettedMesh(surface);
  const nm = safeName(name);
  const py = `#!/usr/bin/env python3
"""Radiation/diffraction of '${nm}' with Capytaine 2.x.

Reads  mesh/${nm}.dat  (NEMOH mesh format, wetted surface, z up, waterline z = 0).
Writes ${nm}_hydro.nc  (complete Capytaine dataset) and
       ${nm}_hydro.csv (kind,omega,headingDeg,i,j,re,im; SI; exp(+i omega t); heading = propagation direction).
"""
import numpy as np
import xarray as xr
import capytaine as cpt

RHO = ${num(rho)}
G = ${num(g)}
DEPTH = ${deep ? 'np.inf' : num(depth)}
OMEGAS = [${om.map(num).join(', ')}]
HEADINGS_DEG = [${hd.map(num).join(', ')}]
CENTRE = (${cog.map(num).join(', ')})
DOFS = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

mesh = cpt.load_mesh("mesh/${nm}.dat", file_format="nemoh")
body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=CENTRE), center_of_mass=CENTRE, name="${nm}")

test_matrix = xr.Dataset(coords={
    "omega": OMEGAS,
    "wave_direction": np.deg2rad(HEADINGS_DEG),
    "radiating_dof": list(body.dofs),
    "water_depth": [DEPTH],
    "rho": [RHO],
    "g": [G],
})
solver = cpt.BEMSolver()
ds = solver.fill_dataset(test_matrix, body)
cpt.io.xarray.separate_complex_values(ds).to_netcdf("${nm}_hydro.nc")

# only the optional dimensions of length one are dropped, never omega / wave_direction
ds = ds.squeeze(dim=[d for d in ("water_depth", "rho", "g") if d in ds.dims], drop=True)
rows = ["# solver=capytaine " + cpt.__version__, "# rho=%r" % RHO, "# g=%r" % G,
        "# depth=%s" % ("inf" if np.isinf(DEPTH) else repr(DEPTH)), "kind,omega,headingDeg,i,j,re,im"]
for w in OMEGAS:
    for j, rad in enumerate(DOFS):
        for i, inf in enumerate(DOFS):
            a = float(ds["added_mass"].sel(omega=w, radiating_dof=rad, influenced_dof=inf))
            b = float(ds["radiation_damping"].sel(omega=w, radiating_dof=rad, influenced_dof=inf))
            rows.append("A,%r,,%d,%d,%r," % (w, i + 1, j + 1, a))
            rows.append("B,%r,,%d,%d,%r," % (w, i + 1, j + 1, b))
    for h in HEADINGS_DEG:
        for i, inf in enumerate(DOFS):
            f = complex(ds["excitation_force"].sel(omega=w, influenced_dof=inf).sel(wave_direction=np.deg2rad(h), method="nearest"))
            # Capytaine uses exp(-i omega t); the conjugate is the exp(+i omega t) amplitude
            rows.append("F,%r,%r,%d,,%r,%r" % (w, h, i + 1, f.real, -f.imag))

# infinite-frequency added mass (constant term of the Cummins equation)
if np.isinf(DEPTH):
    for j, rad in enumerate(DOFS):
        res = solver.solve(cpt.RadiationProblem(body=body, radiating_dof=rad, omega=np.inf, rho=RHO, g=G), keep_details=False)
        for i, inf in enumerate(DOFS):
            rows.append("AINF,,,%d,%d,%r," % (i + 1, j + 1, res.added_masses[inf]))

with open("${nm}_hydro.csv", "w") as fh:
    fh.write("\\n".join(rows) + "\\n")
print("wrote ${nm}_hydro.nc and ${nm}_hydro.csv")
`;
  const notes = [
    `${mesh.nPanels} panels, ${mesh.nNodes} nodes; ${om.length} frequencies × ${hd.length} headings.`,
    'Requires Python 3 with capytaine 2.x (pip install capytaine). The CSV is read back with readCapytaineCsv; the NetCDF file is kept for other tools.',
    'The infinite-frequency added mass is requested for deep water only.',
  ];
  if (mesh.flipped) notes.push('The surface was oriented inwards; panel normals were reversed to point into the fluid.');
  return { files: [{ path: 'run_capytaine.py', text: py }, { path: `mesh/${nm}.dat`, text: nemohMeshText(mesh) }], notes, mesh: { nNodes: mesh.nNodes, nPanels: mesh.nPanels, area: mesh.area, volume: mesh.volume } };
}

/**
 * Reads the CSV written by the script of writeCapytaineScript (schema of parseHydroCsv; rho, g and depth are
 * taken from the `# key=value` comment lines unless given).
 * @param {string} text CSV content
 * @param {{rho?:number, g?:number, depth?:number|null}} [opts] overrides, SI
 * @returns {object} database (createHydroDb) with an extra `source` string
 */
export function readCapytaineCsv(text, opts = {}) {
  const meta = {};
  for (const line of String(text).split(/\r?\n/)) { const m = /^#\s*([A-Za-z]+)\s*=\s*(.+?)\s*$/.exec(line); if (m) meta[m[1].toLowerCase()] = m[2]; }
  const rho = opts.rho ?? (Number(meta.rho) > 0 ? Number(meta.rho) : 1025);
  const g = opts.g ?? (Number(meta.g) > 0 ? Number(meta.g) : G0);
  const d = Number(meta.depth);
  const depth = opts.depth !== undefined ? opts.depth : (Number.isFinite(d) && d > 0 ? d : null);
  const db = parseHydroCsv(text, { rho, g, depth });
  db.source = meta.solver || 'Capytaine CSV';
  return db;
}
