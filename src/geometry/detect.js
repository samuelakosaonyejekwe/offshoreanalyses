// Format detection for uploaded geometry: magic bytes first, then content sniffing of the first kilobytes,
// then the file extension. Extensions alone are not trusted where they are ambiguous in practice
// (.dat, .txt, .xml, .json, .msh, .grd, .asc, .stp, .nc, .inp, .csv, .zip).
//
// Returns format ids of registry.js. Two container ids are also returned and unpacked by the registry:
// 'zip' (multi-file upload) and 'gzip' (single compressed file).

import { extOf, decodeText } from './importers/common.js';
import { isZip, zipEntries } from './zip.js';

const HEAD_BYTES = 65536;

/** Extension → format id when the content gives no better answer. */
export const EXTENSION_MAP = Object.freeze({
  stl: 'stl', obj: 'obj', ply: 'ply', off: 'off', gltf: 'gltf', glb: 'glb', '3mf': '3mf',
  vtk: 'vtk', vtu: 'vtu', vtp: 'vtu', msh: 'msh', su2: 'su2',
  nc: 'netcdf', cdf: 'netcdf', nc4: 'netcdf4', e: 'exodus', exo: 'exodus', ex2: 'exodus', exii: 'exodus', g: 'exodus',
  h5: 'hdf5', hdf5: 'hdf5', he5: 'hdf5', med: 'med', rmed: 'med', cgns: 'cgns',
  bdf: 'nastran', nas: 'nastran', nastran: 'nastran', inp: 'abaqus', cdb: 'ansys-cdb', unv: 'unv',
  step: 'step', stp: 'step', p21: 'step', stpz: 'step', cis: 'cis2', ifc: 'ifc', ifczip: 'ifc', igs: 'iges', iges: 'iges',
  sdnf: 'sdnf', sdn: 'sdnf', dxf: 'dxf', svg: 'svg', gdf: 'gdf', nemoh: 'nemoh',
  csv: 'table', tsv: 'table', tab: 'table', xyz: 'table', txt: 'table', dat: 'table',
  json: 'table-json', xml: 'xml-model',
  geojson: 'geojson', geojsonl: 'geojson', ndjson: 'geojson', shp: 'shapefile', gml: 'gml', citygml: 'gml',
  kml: 'kml', kmz: 'kml', gpx: 'gpx', wkt: 'wkt',
  asc: 'esri-ascii', agr: 'esri-ascii', grd: 'surfer-grd', tif: 'geotiff', tiff: 'geotiff', gtiff: 'geotiff',
  las: 'las', laz: 'laz', e57: 'e57', pts: 'pts', ptx: 'ptx', xyzrgb: 'xyzrgb', xyzi: 'xyzrgb', xyzn: 'xyzrgb',
  gpkg: 'geopackage', sqlite: 'geopackage',
  x_t: 'parasolid', x_b: 'parasolid', xmt_txt: 'parasolid', xmt_bin: 'parasolid', xt: 'parasolid',
  sat: 'acis', sab: 'acis', asat: 'acis', asab: 'acis', dwg: 'dwg', dgn: 'native-cad',
  sldprt: 'native-cad', sldasm: 'native-cad', ipt: 'native-cad', iam: 'native-cad', prt: 'native-cad', asm: 'native-cad',
  catpart: 'native-cad', catproduct: 'native-cad', '3dm': 'native-cad', jt: 'native-cad', rvt: 'native-cad', skp: 'native-cad',
  fbx: 'native-cad', '3ds': 'native-cad', par: 'native-cad', psm: 'native-cad', f3d: 'native-cad', nwd: 'native-cad', rfa: 'native-cad',
  hst: 'hydrostar', cas: 'fluent', topojson: 'geojson',
});

const startsWith = (b, sig, off = 0) => {
  if (!b || b.length < off + sig.length) return false;
  for (let i = 0; i < sig.length; i++) if (b[off + i] !== sig[i]) return false;
  return true;
};
const ascii = (s) => Array.from(s, (c) => c.charCodeAt(0));

const SIG = {
  hdf5: [0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a],
  sqlite: ascii('SQLite format 3\u0000'),
  e57: ascii('ASTM-E57'),
  glb: ascii('glTF'),
  las: ascii('LASF'),
  binDxf: ascii('AutoCAD Binary DXF'),
};

function headText(bytes, text) {
  if (typeof text === 'string') return text.length > HEAD_BYTES ? text.slice(0, HEAD_BYTES) : text;
  if (!bytes) return '';
  return decodeText(bytes.length > HEAD_BYTES ? bytes.subarray(0, HEAD_BYTES) : bytes);
}

function looksBinary(bytes) {
  if (!bytes) return false;
  const n = Math.min(bytes.length, 2048);
  if (n >= 2 && ((bytes[0] === 0xff && bytes[1] === 0xfe) || (bytes[0] === 0xfe && bytes[1] === 0xff))) return false;
  let bad = 0;
  for (let i = 0; i < n; i++) { const c = bytes[i]; if (c === 0 || (c < 9) || (c > 13 && c < 32 && c !== 27)) bad++; }
  return bad > n * 0.02;
}

/** Local name of the root element of an XML text head, or ''. */
export function xmlRootName(head) {
  let i = 0;
  for (;;) {
    const lt = head.indexOf('<', i);
    if (lt < 0) return '';
    const c = head[lt + 1];
    if (c === '?') { const e = head.indexOf('?>', lt); if (e < 0) return ''; i = e + 2; continue; }
    if (c === '!') {
      if (head.startsWith('<!--', lt)) { const e = head.indexOf('-->', lt); if (e < 0) return ''; i = e + 3; continue; }
      const e = head.indexOf('>', lt); if (e < 0) return ''; i = e + 1; continue;
    }
    const m = /^<([A-Za-z_][\w.\-]*:)?([A-Za-z_][\w.\-]*)/.exec(head.slice(lt, lt + 200));
    return m ? m[2] : '';
  }
}

function numericRows(head, maxRows = 60) {
  const rows = [];
  let first = null;
  const lines = head.split(/\r\n|\n|\r/, maxRows + 8);
  if (head.length >= HEAD_BYTES) lines.pop();
  for (const l of lines) {
    const t = l.trim();
    if (!t || t[0] === '#' || t[0] === '!' || t.startsWith('//')) continue;
    const tok = t.split(/[\s,;|]+/).filter(Boolean);
    const nums = tok.map(Number);
    const allNum = nums.every((v) => Number.isFinite(v));
    if (first === null) first = { text: t, tokens: tok, numeric: allNum };
    if (allNum) rows.push(nums);
    else if (rows.length) rows.push(null);
    if (rows.length >= maxRows) break;
  }
  return { first, rows };
}

/** Airfoil coordinate file: an optional title line, then (almost) only x y pairs on a unit or percent chord. */
function sniffAirfoil(head) {
  const { first, rows } = numericRows(head, 80);
  if (!first) return false;
  const data = rows.filter(Boolean);
  if (data.length < 8) return false;
  let pairs = 0, inRange = 0;
  for (const r of data) {
    if (r.length === 2) { pairs++; if (r[0] >= -0.05 && r[0] <= 1.05 && Math.abs(r[1]) <= 0.6) inRange++; }
  }
  const lednicerCounts = data[0].length === 2 && data[0][0] > 2 && data[0][1] > 2;
  return pairs >= data.length - 1 && inRange >= pairs - (lednicerCounts ? 1 : 0) - 1 && inRange >= 8;
}

/** NEMOH mesh: "2 isym" header then "id x y z" node rows. */
function sniffNemoh(head) {
  const { first, rows } = numericRows(head, 12);
  if (!first || !first.numeric || first.tokens.length !== 2 || Number(first.tokens[0]) !== 2 || !(first.tokens[1] === '0' || first.tokens[1] === '1')) return false;
  const data = rows.filter(Boolean).slice(1);
  return data.length >= 1 && data.slice(0, 6).every((r) => r.length === 4 && Number.isInteger(r[0]));
}

/** Panel CSV: twelve numeric columns per row (four vertices) or x1…z4 headers. */
function sniffPanelCsv(head) {
  const { first, rows } = numericRows(head, 20);
  if (!first) return false;
  if (!first.numeric) {
    const h = first.tokens.map((s) => s.toLowerCase().replace(/[^a-z0-9]/g, ''));
    if (h.includes('x1') && h.includes('z4') && h.includes('y3')) return true;
  }
  const data = rows.filter(Boolean);
  return data.length >= 1 && data.every((r) => r.length === 12);
}

/** Point cloud with colour/intensity/normals: header-less rows of 4, 6, 7, 9 or 10 numbers. */
function sniffXyzRgb(head) {
  const { first, rows } = numericRows(head, 30);
  if (!first || !first.numeric) return false;
  const data = rows.filter(Boolean);
  if (data.length < 2) return false;
  const n = data[data.length > 2 ? 1 : 0].length;
  if (![6, 7, 9, 10].includes(n)) return false;
  let same = 0;
  for (const r of data) if (r.length === n) same++;
  if (same < data.length - 1) return false;
  // colours are integers 0–255 (or 0–1 floats) in the last three or the three after xyz[i]
  const col = n === 6 ? 3 : n === 7 ? 4 : n === 9 ? 3 : 4;
  const sample = data.filter((r) => r.length === n);
  const rgb = sample.every((r) => r[col] >= 0 && r[col] <= 255 && r[col + 1] >= 0 && r[col + 1] <= 255 && r[col + 2] >= 0 && r[col + 2] <= 255);
  const normals = sample.every((r) => Math.abs(Math.hypot(r[n - 3], r[n - 2], r[n - 1]) - 1) < 0.02);
  return rgb || normals;
}

function sniffText(head, ext) {
  const t = head.replace(/^﻿/, '');
  const lead = t.slice(0, 4096).trimStart();
  const up = lead.slice(0, 2048).toUpperCase();

  if (lead.startsWith('ISO-10303-21')) {
    const m = /FILE_SCHEMA\s*\(\s*\(([^;]*?)\)\s*\)\s*;/i.exec(t);
    const schema = m ? m[1].toUpperCase() : '';
    if (/IFC\d|IFC4|IFC2X/.test(schema)) return ['ifc', `ISO 10303-21 file with schema ${schema.replace(/['\s]/g, '')}`];
    if (/STRUCTURAL_FRAME_SCHEMA|CIS2|CIMSTEEL/.test(schema)) return ['cis2', 'ISO 10303-21 file with the CIS/2 structural frame schema'];
    return ['step', schema ? `ISO 10303-21 file with schema ${schema.replace(/['\s]/g, '')}` : 'ISO 10303-21 file'];
  }
  if (lead.startsWith('# vtk DataFile')) return ['vtk', 'legacy VTK header'];
  if (lead.startsWith('$MeshFormat') || lead.startsWith('$NOD') || lead.startsWith('$PostFormat')) return ['msh', 'Gmsh section header'];
  if (/^ply\s/.test(lead)) return ['ply', 'PLY magic'];
  if (/^(C|N|ST|CN|4|4N)*OFF\b/.test(lead)) return ['off', 'OFF header'];
  if (lead.startsWith('DSAA')) return ['surfer-grd', 'Surfer ASCII grid magic DSAA'];
  if (lead.startsWith('AutoCAD Binary DXF')) return ['dxf', 'binary DXF sentinel'];
  if (lead.startsWith('**ABCDEFGHIJKLMNOPQRSTUVWXYZ') || /^\*\*PART1;/m.test(lead)) return ['parasolid', 'Parasolid transmit-file header'];

  if (lead[0] === '<') {
    const root = xmlRootName(t);
    const r = root.toLowerCase();
    if (r === 'landxml') return ['landxml', 'root element LandXML'];
    if (r === 'kml') return ['kml', 'root element kml'];
    if (r === 'gpx') return ['gpx', 'root element gpx'];
    if (r === 'svg') return ['svg', 'root element svg'];
    if (r === 'vtkfile') return ['vtu', 'root element VTKFile'];
    if (r === 'citymodel' || r === 'featurecollection' || r === 'topocomplex' || /xmlns(:\w+)?\s*=\s*["']http:\/\/www\.opengis\.net\/(gml|citygml|infragml|landinfra)/i.test(t.slice(0, 8192))) return ['gml', `GML document (root element ${root})`];
    if (r === 'model' && /3dmanufacturing|3mf/i.test(t.slice(0, 2048))) return ['3mf', '3MF model part (unzipped)'];
    if (r === 'html') return [null, 'HTML page'];
    if (root) return ['xml-model', `XML document with root element ${root}`];
  }

  if (lead[0] === '{' || lead[0] === '[') {
    const h = t.slice(0, 16384);
    if (/"format"\s*:\s*"offshore-analyses-geometry"/.test(h)) return ['oa-json', 'application geometry JSON marker'];
    if (/"asset"\s*:/.test(h) && /"(scenes|meshes|accessors|buffers|nodes)"\s*:/.test(t)) return ['gltf', 'glTF JSON (asset + scene graph)'];
    if (/"type"\s*:\s*"(FeatureCollection|Feature|GeometryCollection|Point|MultiPoint|LineString|MultiLineString|Polygon|MultiPolygon|Topology)"/.test(h)) return ['geojson', 'GeoJSON type member'];
    return ['table-json', 'JSON document'];
  }

  if (/^(SRID\s*=\s*\d+\s*;)?\s*(POINT|LINESTRING|POLYGON|MULTIPOINT|MULTILINESTRING|MULTIPOLYGON|GEOMETRYCOLLECTION|TIN|POLYHEDRALSURFACE|TRIANGLE)\s*(ZM|Z|M)?\s*(\(|EMPTY)/i.test(lead)) return ['wkt', 'WKT geometry keyword'];
  if (/^\s*(ncols|nrows)\s+\d+/im.test(lead.slice(0, 400)) && /cellsize|dx\s|xll/i.test(lead.slice(0, 600))) return ['esri-ascii', 'ESRI ASCII grid header'];
  if (/^NDIME\s*=/m.test(lead) || (/^NPOIN\s*=/m.test(t) && /^NELEM\s*=/m.test(t))) return ['su2', 'SU2 NDIME/NPOIN/NELEM keywords'];
  if (/^\s*FoamFile\b/m.test(lead) || lead.startsWith('/*--------------------------------*- C++ -*')) return ['openfoam', 'OpenFOAM file header'];
  if (/^\(\s*[0-2]\s+["(]/.test(lead) || /^\(\s*10\s+\(/m.test(lead.slice(0, 600))) return ['fluent', 'Fluent mesh s-expressions'];
  if (/^\s{0,4}-1\s*\r?\n\s*\d{2,4}\s*\r?\n/.test(t.slice(0, 200)) && (ext === 'unv' || ext === '' || ext === 'dat' || ext === 'txt')) return ['unv', 'I-DEAS universal dataset delimiter'];

  // IGES: the flag letter in column 73 of the first records
  {
    const l1 = t.split(/\r?\n/, 3);
    if (l1[0] && l1[0].length >= 73 && l1[0][72] === 'S' && (ext === 'igs' || ext === 'iges' || (l1[1] && /[SG]/.test(l1[1][72] || '')))) return ['iges', 'IGES start-section records'];
    if (l1[0] && l1[0].length >= 73 && l1[0][72] === 'C') return ['iges', 'compressed-ASCII IGES flag record'];
  }

  const dxfLead = lead.replace(/\r/g, '');
  if (/^\s*0\s*\n\s*SECTION\s*\n/.test(dxfLead) || /^\s*999\s*\n.*\n\s*0\s*\n\s*SECTION/.test(dxfLead)) return ['dxf', 'DXF group-code pairs'];
  if (/^\s*(?:\*\*.*\r?\n\s*)*\*(HEADING|NODE|PART|PREPRINT|ELEMENT|ASSEMBLY|MATERIAL|INCLUDE|PARAMETER)\b/i.test(lead) || (ext === 'inp' && /^\*NODE\b/im.test(t))) return ['abaqus', 'Abaqus keyword lines'];
  if (/^\s*(NBLOCK|EBLOCK)\s*,/im.test(t) || /^\s*\/(PREP7|COM|BATCH|NOPR)\b/im.test(lead) || (ext === 'cdb' && /^\s*ET\s*,/im.test(t))) return ['ansys-cdb', 'ANSYS archive commands'];
  if (/^\s*BEGIN\s+BULK\b/im.test(t) || /^(GRID\*?|CBAR|CBEAM|CROD|CQUAD4|CTRIA3|CTETRA|CHEXA|PSHELL|MAT1)[\s,*]/m.test(t) || /^\s*(SOL\s+\d+|CEND)\s*$/im.test(lead)) return ['nastran', 'Nastran bulk-data cards'];
  if (/^\s*Packet\s+0*\d+\s*$/im.test(lead) || /^\s*Packet\s+00\b/im.test(lead)) return ['sdnf', 'SDNF packet header'];
  if (/^\s*solid\b/.test(lead) && /\bfacet\b|\bendsolid\b/.test(t)) return ['stl', 'ASCII STL keywords'];
  if (ext === 'obj' || (/^(v|vn|vt|f|o|g|mtllib|usemtl)\s/m.test(lead) && /^v\s+[-+.\d]/m.test(t) && /^(f|l|p)\s+\d/m.test(t))) return ['obj', 'Wavefront OBJ records'];
  if (ext === 'gdf') return ['gdf', 'extension .gdf'];
  if (ext === 'ptx') return ['ptx', 'extension .ptx'];
  if (ext === 'pts') return ['pts', 'extension .pts'];
  if (up.startsWith('HYDROSTAR') || ext === 'hst') return ['hydrostar', 'Hydrostar mesh'];
  if (sniffNemoh(t)) return ['nemoh', 'NEMOH mesh header "2 isym" and node rows'];
  if ((ext === 'dat' || ext === 'txt' || ext === '' || ext === 'csv' || ext === 'air' || ext === 'afl') && sniffAirfoil(t)) return ['airfoil', 'title line followed by x y pairs on a unit chord'];
  if ((ext === 'csv' || ext === 'txt' || ext === 'dat' || ext === 'tsv') && sniffPanelCsv(t)) return ['panel-csv', 'twelve coordinate columns per row'];
  if (['xyzrgb', 'xyzi', 'xyzn'].includes(ext) || (['txt', 'asc', 'xyz', 'pts', 'csv', 'dat', ''].includes(ext) && sniffXyzRgb(t))) return ['xyzrgb', 'header-less point rows with colour, intensity or normals'];
  return null;
}

function sniffBinary(bytes, ext) {
  if (startsWith(bytes, SIG.glb)) return ['glb', 'binary glTF magic'];
  if (startsWith(bytes, SIG.hdf5) || startsWith(bytes, SIG.hdf5, 512) || startsWith(bytes, SIG.hdf5, 1024) || startsWith(bytes, SIG.hdf5, 2048)) {
    if (ext === 'med' || ext === 'rmed') return ['med', 'HDF5 container (.med)'];
    if (ext === 'cgns') return ['cgns', 'HDF5 container (.cgns)'];
    if (['e', 'exo', 'ex2', 'exii', 'g'].includes(ext)) return ['exodus-hdf5', 'HDF5 container (Exodus on NetCDF-4)'];
    if (['nc', 'nc4', 'cdf'].includes(ext)) return ['netcdf4', 'HDF5 container (NetCDF-4)'];
    return ['hdf5', 'HDF5 signature'];
  }
  if (bytes.length >= 4 && bytes[0] === 0x43 && bytes[1] === 0x44 && bytes[2] === 0x46 && (bytes[3] === 1 || bytes[3] === 2 || bytes[3] === 5)) {
    return [['e', 'exo', 'ex2', 'exii', 'g'].includes(ext) ? 'exodus' : 'netcdf', `NetCDF classic magic CDF-${bytes[3]}`];
  }
  if (bytes.length >= 4 && ((bytes[0] === 0x49 && bytes[1] === 0x49 && (bytes[2] === 42 || bytes[2] === 43) && bytes[3] === 0) || (bytes[0] === 0x4d && bytes[1] === 0x4d && bytes[2] === 0 && (bytes[3] === 42 || bytes[3] === 43)))) return ['geotiff', 'TIFF byte-order mark and version'];
  if (startsWith(bytes, SIG.las)) return [ext === 'laz' || (bytes.length > 104 && (bytes[104] & 0x80)) ? 'laz' : 'las', 'LASF signature'];
  if (startsWith(bytes, SIG.e57)) return ['e57', 'ASTM-E57 signature'];
  if (startsWith(bytes, SIG.sqlite)) return ['geopackage', 'SQLite database header'];
  if (startsWith(bytes, ascii('DSBB'))) return ['surfer-grd', 'Surfer 6 binary grid magic DSBB'];
  if (startsWith(bytes, ascii('DSRR'))) return ['surfer-grd', 'Surfer 7 binary grid magic DSRR'];
  if (startsWith(bytes, SIG.binDxf)) return ['dxf', 'binary DXF sentinel'];
  if (startsWith(bytes, ascii('AC10')) && bytes.length > 6 && bytes[4] >= 0x30 && bytes[4] <= 0x39) return ['dwg', 'DWG version string'];
  if (bytes.length >= 100 && bytes[0] === 0 && bytes[1] === 0 && bytes[2] === 0x27 && bytes[3] === 0x0a) return ['shapefile', 'ESRI shapefile file code 9994'];
  if (startsWith(bytes, ascii('ply\n')) || startsWith(bytes, ascii('ply\r'))) return ['ply', 'PLY magic'];
  if (startsWith(bytes, ascii('# vtk DataFile'))) return ['vtk', 'legacy VTK header'];
  if (startsWith(bytes, ascii('$MeshFormat'))) return ['msh', 'Gmsh section header'];
  if (startsWith(bytes, ascii('<?xml')) || startsWith(bytes, ascii('<VTKFile'))) {
    const head = decodeText(bytes.subarray(0, Math.min(bytes.length, 4096)));
    if (/<VTKFile\b/.test(head)) return ['vtu', 'root element VTKFile'];
  }
  if (startsWith(bytes, ascii('PS\u0000\u0000'))) return ['parasolid', 'Parasolid binary transmit header'];
  if (bytes.length >= 84) {
    const n = bytes[80] | (bytes[81] << 8) | (bytes[82] << 16) | (bytes[83] * 16777216);
    if (84 + 50 * n === bytes.length) return ['stl', 'binary STL: length matches the triangle count'];
  }
  return null;
}

const POLYMESH = ['points', 'faces', 'owner'];
const base = (k) => { const s = k.replace(/\\/g, '/').toLowerCase(); return s.slice(s.lastIndexOf('/') + 1).replace(/\.gz$/, ''); };

/**
 * Decide what a set of file names is. Used for multi-file uploads and ZIP listings.
 * @param {Iterable<string>} names
 * @returns {{ id: string, main: string|null, reason: string }|null}
 */
export function detectFileSet(names) {
  const list = [...names].filter((n) => !/(^|\/)(__MACOSX\/|\.DS_Store$)/.test(n) && !n.endsWith('/'));
  const bases = list.map(base);
  if (POLYMESH.every((p) => bases.includes(p))) return { id: 'openfoam', main: list[bases.indexOf('points')], reason: 'polyMesh files points/faces/owner' };
  const byExt = (e) => list.find((n) => extOf(n) === e);
  if (list.some((n) => /(^|\/)3D\/[^/]+\.model$/i.test(n)) || (list.includes('[Content_Types].xml') && byExt('model'))) return { id: '3mf', main: null, reason: '3MF package parts' };
  if (byExt('shp')) return { id: 'shapefile', main: byExt('shp'), reason: 'contains a .shp' };
  if (byExt('gltf')) return { id: 'gltf', main: byExt('gltf'), reason: 'contains a .gltf' };
  if (byExt('kml')) return { id: 'kml', main: list.find((n) => base(n) === 'doc.kml') || byExt('kml'), reason: 'contains a .kml' };
  // companions that never stand alone
  const companions = new Set(['prj', 'dbf', 'shx', 'cpg', 'sbn', 'sbx', 'qpj', 'qix', 'bin', 'mtl', 'png', 'jpg', 'jpeg', 'aux', 'tfw', 'ovr', 'md', 'pdf']);
  const mains = list.filter((n) => !companions.has(extOf(n)));
  if (mains.length === 1) return { id: null, main: mains[0], reason: 'single data file with companions' };
  if (mains.length > 1) {
    const known = mains.filter((n) => EXTENSION_MAP[extOf(n)]);
    if (known.length >= 1) return { id: null, main: known[0], reason: known.length > 1 ? `several data files; using ${known[0]}` : 'one recognised data file', others: known.slice(1) };
  }
  return null;
}

/**
 * Detect the format of an importer input.
 * @param {{ name?: string, bytes?: Uint8Array, text?: string, files?: Map<string, Uint8Array> }} input
 * @returns {{ id: string, reason: string, main?: string }|null} null when nothing matches
 */
export function detectFormat(input) {
  const name = input.name || '';
  let ext = extOf(name);
  const bytes = input.bytes || null;

  if (input.files && input.files.size && !bytes && typeof input.text !== 'string') {
    const set = detectFileSet(input.files.keys());
    if (!set) return null;
    if (set.id) return set;
    const inner = detectFormat({ name: set.main, bytes: input.files.get(set.main) });
    return inner ? { ...inner, main: set.main, reason: `${inner.reason} (${set.reason})` } : null;
  }

  if (bytes && bytes.length >= 2 && bytes[0] === 0x1f && bytes[1] === 0x8b) return { id: 'gzip', reason: 'gzip magic' };
  if (bytes && isZip(bytes)) {
    let names = [];
    try { names = zipEntries(bytes).map((e) => e.name); } catch { return { id: 'zip', reason: 'ZIP signature' }; }
    if (ext === 'kmz') return { id: 'kml', reason: 'KMZ archive' };
    if (ext === '3mf') return { id: '3mf', reason: '3MF package' };
    const set = detectFileSet(names);
    if (set && (set.id === '3mf' || set.id === 'kml')) return { id: set.id, reason: set.reason };
    return { id: 'zip', reason: set ? `ZIP archive: ${set.reason}` : 'ZIP archive' };
  }

  if (bytes) {
    const b = sniffBinary(bytes, ext);
    if (b) return { id: b[0], reason: b[1] };
    if (looksBinary(bytes)) {
      const id = EXTENSION_MAP[ext];
      if (ext === 'stl') return { id: 'stl', reason: 'extension .stl (binary)' };
      if (id) return { id, reason: `extension .${ext} (binary content)` };
      return null;
    }
  }

  const head = headText(bytes, input.text);
  if (head) {
    const s = sniffText(head, ext);
    if (s) return s[0] ? { id: s[0], reason: s[1] } : null;
  }
  if (ext === 'grd') return { id: 'surfer-grd', reason: 'extension .grd' };
  const id = EXTENSION_MAP[ext];
  if (id) return { id, reason: `extension .${ext}` };
  if (head && /\d/.test(head)) {
    const { rows } = numericRows(head, 10);
    if (rows.filter(Boolean).length >= 2) return { id: 'table', reason: 'rows of numbers' };
  }
  return null;
}
