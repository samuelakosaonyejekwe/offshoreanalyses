// CFD mesh formats: SU2 native meshes (.su2) and OpenFOAM polyMesh directories.
//
// SU2: volume (3-D) or surface (2-D) elements plus the boundary markers. Marker faces of a 3-D mesh become a
// surface whose groups are the markers; marker edges of a 2-D mesh become a frame-kind model named
// "(boundary edges)" whose member groups are the markers (they are boundaries, not beams).
//
// OpenFOAM: `points`, `faces`, `owner` and `boundary` (ASCII, optionally gzipped). The boundary patches are
// returned as a triangulated surface with one group per patch. Cells are arbitrary polyhedra, which the
// volume model kind cannot hold, so the interior is not reconstructed — the cell count is reported only.

import { lineCursor, baseName, FeMesh, Grow, finish, surfaceModel, findFile, decodeText, pushPolygon, tally } from './common.js';
import { feResult } from './fe.js';
import { resolveUnits, scaleArray } from '../units.js';
import { gunzip } from '../zip.js';

const SU2 = { 3: ['line2', 2], 5: ['tri3', 3], 9: ['quad4', 4], 10: ['tet4', 4], 12: ['hex8', 8], 13: ['wedge6', 6, [0, 2, 1, 3, 5, 4]], 14: ['pyr5', 5] };

/** SU2 native mesh. */
export function parseSu2(input, opts = {}) {
  const cur = lineCursor(input);
  const name = baseName(input.name);
  const fe = new FeMesh(name);
  const warnings = [];
  const markers = [];
  let ndime = 0, eid = 0, tagNo = 0, sawKey = false;
  const kv = (l) => { const m = /^([A-Z_]+)\s*=\s*(.*)$/.exec(l); return m ? [m[1], m[2].trim()] : null; };
  const element = (l, tag) => {
    const p = l.split(/\s+/).map(Number);
    const def = SU2[p[0]];
    if (!def || p.length < def[1] + 1) { fe.skip(`su2 element type ${p[0]}`); return 0; }
    const nodes = p.slice(1, 1 + def[1]);
    fe.addElement(def[0], def[2] ? def[2].map((k) => nodes[k]) : nodes, ++eid, tag);
    return eid;
  };
  for (let line = cur.nextNonBlank(); line !== null; line = cur.nextNonBlank()) {
    if (line[0] === '%') continue;
    const k = kv(line);
    if (!k) continue;
    sawKey = true;
    if (k[0] === 'NDIME') ndime = Number(k[1]);
    else if (k[0] === 'NELEM') { const n = Number(k[1]); for (let i = 0; i < n; i++) { const l = cur.nextNonBlank(); if (l === null) break; element(l, 0); } } else if (k[0] === 'NPOIN') {
      const n = Number(k[1].split(/\s+/)[0]);
      for (let i = 0; i < n; i++) { const l = cur.nextNonBlank(); if (l === null) throw new Error(`${input.name || 'SU2'}: file ends after ${i} of ${n} points.`); const p = l.split(/\s+/); fe.addNode(i, Number(p[0]), Number(p[1]), ndime === 3 ? Number(p[2]) : 0); }
    } else if (k[0] === 'MARKER_TAG') {
      const tag = ++tagNo;
      const mk = { name: k[1], ids: [] };
      fe.setProp(tag, { name: k[1] });
      const ml = cur.nextNonBlank();
      const mv = ml && kv(ml);
      const n = mv && mv[0] === 'MARKER_ELEMS' ? Number(mv[1]) : 0;
      for (let i = 0; i < n; i++) { const l = cur.nextNonBlank(); if (l === null) break; const id = element(l, tag); if (id) mk.ids.push(id); }
      fe.addElemSet(k[1], mk.ids);
      markers.push({ name: k[1], elements: mk.ids.length });
    } else if (!['NMARK', 'MARKER_ELEMS', 'NZONE', 'IZONE', 'NCHUNK', 'NPERIODIC', 'AOA_OFFSET', 'AOS_OFFSET'].includes(k[0])) fe.skip(`keyword ${k[0]}`);
  }
  if (!sawKey || !ndime) throw new Error(`${input.name || 'SU2'}: not an SU2 mesh (expected NDIME=, NELEM= and NPOIN= keyword lines).`);
  const r = feResult(fe, 'su2', opts, { name: input.name || 'SU2', expected: 'expected an NPOIN= block', warnings, extra: { dimension: ndime, markers } });
  for (const m of r.models) {
    if (m.kind === 'frame' && ndime === 2) m.name = `${name} (boundary edges)`;
    if (m.kind === 'surface' && ndime === 3) m.name = `${name} (boundary markers)`;
  }
  return r;
}

function foamBody(text, what, fileName) {
  // strip comments, then the FoamFile dictionary
  let t = text.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/[^\n]*/g, '');
  const h = /FoamFile\s*\{([^}]*)\}/.exec(t);
  if (h) {
    if (/format\s+binary/.test(h[1])) throw new Error(`${fileName}: binary OpenFOAM files are not read. Convert the case with "foamFormatConvert" after setting writeFormat ascii in system/controlDict.`);
    t = t.slice(h.index + h[0].length);
  }
  const m = /(\d+)\s*\(/.exec(t);
  if (!m) throw new Error(`${fileName}: expected "${what}" as a list of the form "N ( … )".`);
  return { n: Number(m[1]), body: t.slice(m.index + m[0].length) };
}

/** OpenFOAM polyMesh directory (points, faces, owner, boundary). */
export async function parseOpenFoam(input, opts = {}) {
  const files = input.files;
  const get = async (nm) => {
    let f = findFile(files, `constant/polyMesh/${nm}`) || findFile(files, `polyMesh/${nm}`) || findFile(files, nm);
    if (!f) { const gz = findFile(files, `${nm}.gz`); if (gz) f = await gunzip(gz); }
    if (f && f[0] === 0x1f && f[1] === 0x8b) f = await gunzip(f);
    return f ? decodeText(f) : null;
  };
  if (!files || !files.size) throw new Error(`${input.name || 'polyMesh'}: an OpenFOAM mesh is a directory. Upload the files points, faces, owner and boundary from constant/polyMesh together (or as a ZIP of the folder).`);
  const [tp, tf, to, tb] = [await get('points'), await get('faces'), await get('owner'), await get('boundary')];
  const missing = [['points', tp], ['faces', tf]].filter((x) => !x[1]).map((x) => x[0]);
  if (missing.length) throw new Error(`${input.name || 'polyMesh'}: missing polyMesh file(s): ${missing.join(', ')} (expected points, faces, owner and boundary).`);
  const warnings = [];
  const P = foamBody(tp, 'points', 'points');
  const coords = new Float64Array(P.n * 3);
  {
    const re = /\(\s*([-+.\deE]+)\s+([-+.\deE]+)\s+([-+.\deE]+)\s*\)/g;
    let i = 0;
    for (let m = re.exec(P.body); m && i < P.n; m = re.exec(P.body), i++) { coords[3 * i] = +m[1]; coords[3 * i + 1] = +m[2]; coords[3 * i + 2] = +m[3]; }
    if (i < P.n) throw new Error(`points: the list declares ${P.n} points but only ${i} "(x y z)" entries were found.`);
  }
  const F = foamBody(tf, 'faces', 'faces');
  const faceStart = new Uint32Array(F.n + 1);
  const faceIds = new Grow(Uint32Array, F.n * 4);
  {
    const re = /(\d+)\s*\(([^)]*)\)/g;
    let i = 0;
    for (let m = re.exec(F.body); m && i < F.n; m = re.exec(F.body), i++) { for (const s of m[2].trim().split(/\s+/)) faceIds.push(Number(s)); faceStart[i + 1] = faceIds.n; }
    if (i < F.n) throw new Error(`faces: the list declares ${F.n} faces but only ${i} "n(i j k …)" entries were found.`);
  }
  let nCells = null;
  if (to) { const O = foamBody(to, 'owner', 'owner'); let mx = -1; const re = /\d+/g; let c = 0; for (let m = re.exec(O.body); m && c < O.n; m = re.exec(O.body), c++) { const v = Number(m[0]); if (v > mx) mx = v; } nCells = mx + 1; }
  const patches = [];
  if (tb) {
    const B = foamBody(tb, 'boundary', 'boundary');
    const re = /([A-Za-z_][\w.\-:]*)\s*\{([^}]*)\}/g;
    for (let m = re.exec(B.body); m; m = re.exec(B.body)) {
      const g = (k) => { const x = new RegExp(`\\b${k}\\s+([^;]+);`).exec(m[2]); return x ? x[1].trim() : null; };
      const nFaces = Number(g('nFaces')), startFace = Number(g('startFace'));
      if (nFaces >= 0 && startFace >= 0) patches.push({ name: m[1], type: g('type') || 'patch', nFaces, startFace });
    }
  } else warnings.push('No boundary file was supplied: every face (internal faces included) is returned as one group.');
  const units = resolveUnits(opts, 'm', 'm');
  scaleArray(coords, units.scale);
  const tris = new Grow(Uint32Array, 3072), groups = [], hist = {};
  let bad = 0;
  const emitFaces = (start, count) => {
    for (let f = start; f < start + count && f < F.n; f++) {
      const ids = Array.from(faceIds.a.subarray(faceStart[f], faceStart[f + 1]));
      if (ids.length < 3 || ids.some((k) => k >= P.n)) { bad++; continue; }
      tally(hist, ids.length === 3 ? 'tri3' : ids.length === 4 ? 'quad4' : 'polygon');
      pushPolygon(tris, ids, coords);
    }
  };
  if (patches.length) for (const pt of patches) { const s = tris.n; emitFaces(pt.startFace, pt.nFaces); if (tris.n > s) groups.push({ name: pt.name, start: s, count: tris.n - s, type: pt.type }); } else emitFaces(0, F.n);
  if (!tris.n) throw new Error(`${input.name || 'polyMesh'}: no boundary faces could be built (the boundary file lists ${patches.length} patch(es), the faces file ${F.n} face(s)).`);
  // keep only the boundary vertices
  const map = new Int32Array(P.n).fill(-1);
  const used = new Grow(Float64Array, 3072);
  const indices = tris.finish();
  for (let i = 0; i < indices.length; i++) { const k = indices[i]; if (map[k] < 0) { map[k] = used.n / 3; used.push3(coords[3 * k], coords[3 * k + 1], coords[3 * k + 2]); } indices[i] = map[k]; }
  if (bad) warnings.push(`${bad} face(s) with fewer than three or out-of-range point labels were skipped.`);
  warnings.push(`OpenFOAM cells are general polyhedra and are not reconstructed: only the boundary patches are returned${nCells !== null ? ` (the mesh has ${nCells} cells)` : ''}.`);
  const r = finish([surfaceModel(`${baseName(input.name)} (boundary)`, used.finish(), indices, groups.length ? groups : null)], warnings, {
    format: 'openfoam', units, crs: null, elementTypes: hist, skipped: { invalidFaces: bad, internalFaces: Math.max(0, F.n - patches.reduce((a, p) => a + p.nFaces, 0)), cells: nCells ?? 0 },
    points: P.n, faces: F.n, cells: nCells, patches,
  });
  return r;
}
