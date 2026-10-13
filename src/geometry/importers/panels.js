// Hydrodynamic panel meshes for potential-flow analysis: WAMIT low-order GDF, NEMOH mesh files and HydroStar
// (Bureau Veritas) .hst hull meshes.
//
// All describe the wetted body surface as quadrilateral panels (triangles as quads with a repeated vertex),
// possibly exploiting planes of symmetry. The result is a triangulated surface with coincident vertices
// welded; the original panels are kept as `quads` (four vertex indices per panel, a triangle repeating its
// last vertex). By default symmetric halves/quarters are mirrored into the full body, each mirror image a
// separate group with reversed winding so that normals keep pointing the same way relative to the body; pass
// `opts.expandSymmetry = false` to keep only what the file stores. `meta.symmetry` records what was declared.
//
// WAMIT: length unit inferred from GRAV (9.80665 → m, 32.174 → ft, 386.09 → in, 9806.65 → mm); ULEN is the
// non-dimensionalising length and does not scale the coordinates here. Higher-order GDF (B-spline patches,
// IGDEF ≠ 0 forms) is not read.
//
// HydroStar .hst: keyword lines around a node table (COORDINATES … ENDCOORDINATES, "id x y z" or "x y z") and
// a panel table (PANEL TYPE 0: "n1 n2 n3 [n4]"; PANEL TYPE 1: "id n1 n2 n3 [n4]"; some writers spell the
// keywords PANELS / ENDPANELS) closed by ENDFILE. Panels
// name their nodes by the identifiers of the node table, which need neither start at 1 nor be consecutive.
// SYMMETRY n (whole file) or SYMMETRY_BODY body n: 1 = the stored half is mirrored about y = 0, 2 = the stored
// quarter about y = 0 and x = 0. NUMPANEL body first last assigns panel numbers to the bodies of a multi-body
// file. Coordinates are metres. Other keywords (RHO, GRAVITY, REFLENGTH, ZONEDAMPING, tank and free-surface
// data …) are not geometry: they are passed over and listed in a warning. There is no public specification;
// the layout follows the files written by HydroStar's own mesher and the open readers and writers of the
// format (Capytaine, meshmagick, BEMRosetta, the NEMOH Gmsh converter), used as references only.

import { lineCursor, baseName, Grow, Welder, finish, surfaceModel, unitAdvice, tally } from './common.js';
import { resolveUnits, suggestUnits } from '../units.js';

/**
 * `parts` (optional) splits the panel list into named bodies, each `{ name, from, to, sym }` with its own
 * symmetry; every body × mirror image is then a group.
 */
function buildPanels(name, panels, sym, opts, units, format, extra, warnings, parts) {
  const s = units.scale;
  const expand = opts.expandSymmetry !== false;
  const imagesOf = (q) => {
    const im = [[1, 1]];
    if (expand && q.x) im.push([-1, 1]);
    if (expand && q.y) { im.push([1, -1]); if (q.x) im.push([-1, -1]); }
    return im;
  };
  const weld = new Welder(0, panels.length / 3);
  const tris = new Grow(Uint32Array, panels.length), quads = new Grow(Uint32Array, panels.length / 3);
  const groups = [];
  let nq = 0, nt = 0, bad = 0;
  const n = panels.length / 12;
  for (const part of parts || [{ name: null, from: 0, to: n, sym }]) for (const [mx, my] of imagesOf(part.sym)) {
    const start = tris.n;
    const flip = mx * my < 0;
    for (let p = part.from; p < part.to; p++) {
      const v = [0, 0, 0, 0];
      for (let k = 0; k < 4; k++) { const o = 12 * p + 3 * k; v[k] = weld.add(panels[o] * s * mx, panels[o + 1] * s * my, panels[o + 2] * s); }
      if (flip) { const t = v[1]; v[1] = v[3]; v[3] = t; }
      const u = v.filter((q, k) => v.indexOf(q) === k);
      if (u.length < 3) { bad++; continue; }
      if (u.length === 3) { tris.push(u[0]); tris.push(u[1]); tris.push(u[2]); quads.push(u[0]); quads.push(u[1]); quads.push(u[2]); quads.push(u[2]); nt++; } else { tris.push(v[0]); tris.push(v[1]); tris.push(v[2]); tris.push(v[0]); tris.push(v[2]); tris.push(v[3]); quads.push(v[0]); quads.push(v[1]); quads.push(v[2]); quads.push(v[3]); nq++; }
    }
    const image = mx === 1 && my === 1 ? 'as stored' : `mirror ${mx < 0 ? 'x' : ''}${my < 0 ? 'y' : ''}`;
    if (tris.n > start) groups.push({ name: part.name === null ? image : image === 'as stored' ? part.name : `${part.name} (${image})`, start, count: tris.n - start });
  }
  if (!tris.n) throw new Error(`${name}: no valid panels were found.`);
  if (bad) warnings.push(`${bad} degenerate panel(s) (fewer than three distinct vertices) were skipped.`);
  const coords = weld.finish();
  let above = 0;
  for (let i = 2; i < coords.length; i += 3) if (coords[i] > 1e-9) above++;
  if (above) warnings.push(`${above} vertex/vertices lie above z = 0: the mesh includes dry surface or is not referred to the still-water plane.`);
  const model = surfaceModel(name, coords, tris.finish(), groups.length > 1 ? groups : null, { quads: quads.finish() });
  const hist = {};
  if (nq) hist.quad4 = nq;
  if (nt) hist.tri3 = nt;
  const r = finish([model], warnings, { format, units, crs: null, elementTypes: hist, skipped: { degeneratePanels: bad }, symmetry: { x: !!sym.x, y: !!sym.y, expanded: expand && (sym.x || sym.y), note: 'x: mirror plane x = 0; y: mirror plane y = 0' }, panelsStored: n, ...extra });
  unitAdvice(units, r.meta.bbox, suggestUnits, r.warnings);
  return r;
}

/** WAMIT low-order geometric data file (.gdf). */
export function parseGdf(input, opts = {}) {
  const cur = lineCursor(input);
  const fname = input.name || 'GDF';
  const header = cur.next();
  const need = (what) => { const l = cur.nextNonBlank(); if (l === null) throw new Error(`${fname}: the file ends before ${what}.`); return l.split(/[\s,]+/).filter(Boolean); };
  if (header === null) throw new Error(`${fname}: empty file (expected a WAMIT GDF: header line, "ULEN GRAV", "ISX ISY", "NPAN", then panel vertices).`);
  const l2 = need('the "ULEN GRAV" line');
  const ulen = Number(l2[0]), grav = Number(l2[1]);
  const l3 = need('the "ISX ISY" line');
  const isx = Number(l3[0]), isy = Number(l3[1]);
  const l4 = need('the "NPAN" line');
  const npan = Number(l4[0]);
  if (!(ulen > 0) || !(grav > 0) || ![0, 1, -1].includes(isx) || ![0, 1, -1].includes(isy) || !(npan > 0)) throw new Error(`${fname}: not a WAMIT low-order GDF (expected line 2 "ULEN GRAV", line 3 "ISX ISY" with values 0/1, line 4 "NPAN").`);
  if (l4.length >= 2 && Number.isFinite(Number(l4[1])) && !/[a-df-zA-DF-Z]/.test(l4[1])) throw new Error(`${fname}: this is a higher-order GDF (line 4 holds "NPATCH IGDEF" = ${l4[0]} ${l4[1]}). Only low-order panel GDF files are read; export flat panels (ILOWHI = 0) from the mesher.`);
  const panels = new Grow(Float64Array, npan * 12);
  for (let l = cur.nextNonBlank(); l !== null && panels.n < npan * 12; l = cur.nextNonBlank()) for (const t of l.split(/[\s,]+/)) { if (t === '') continue; const v = Number(t); if (v === v && panels.n < npan * 12) panels.push(v); }
  const warnings = [];
  const got = Math.floor(panels.n / 12);
  if (got < npan) warnings.push(`NPAN declares ${npan} panels but the file holds ${got}.`);
  if (!got) throw new Error(`${fname}: no panel vertices follow the NPAN line (expected 12 numbers per panel).`);
  const hint = Math.abs(grav - 9.80665) < 0.05 ? 'm' : Math.abs(grav - 32.174) < 0.1 ? 'ft' : Math.abs(grav - 386.09) < 1 ? 'in' : Math.abs(grav - 9806.65) < 50 ? 'mm' : null;
  if (!hint) warnings.push(`GRAV = ${grav} does not match a standard gravity in m, ft, in or mm; the length unit could not be inferred.`);
  const units = resolveUnits(opts, hint, 'm');
  return buildPanels(baseName(input.name), panels.view().subarray(0, got * 12), { x: isx === 1, y: isy === 1 }, opts, units, 'gdf', { header: header.trim(), ulen, grav, npan }, warnings);
}

/** NEMOH mesh file: "2 isym" header, node table, panel connectivity. */
export function parseNemoh(input, opts = {}) {
  const cur = lineCursor(input);
  const fname = input.name || 'NEMOH mesh';
  const h = cur.nextNonBlank();
  const hp = h === null ? [] : h.split(/\s+/).map(Number);
  if (hp.length < 2 || hp[0] !== 2 || ![0, 1].includes(hp[1])) throw new Error(`${fname}: not a NEMOH mesh (expected the first line "2 0" or "2 1": the format flag and the xOz symmetry switch).`);
  const nodes = new Map();
  let l;
  for (l = cur.nextNonBlank(); l !== null; l = cur.nextNonBlank()) { const p = l.split(/\s+/).map(Number); if (p.length < 4) break; if (p[0] === 0) { l = cur.nextNonBlank(); break; } nodes.set(p[0], [p[1], p[2], p[3]]); }
  if (!nodes.size) throw new Error(`${fname}: no node rows ("id x y z") follow the header.`);
  const panels = new Grow(Float64Array, 1200);
  let missing = 0;
  for (; l !== null; l = cur.nextNonBlank()) {
    const p = l.split(/\s+/).map(Number);
    if (p.length < 4 || (p[0] === 0 && p[1] === 0 && p[2] === 0 && p[3] === 0)) break;
    const v = p.slice(0, 4).map((k) => nodes.get(k));
    if (v.some((q) => !q)) { missing++; continue; }
    for (const q of v) panels.push3(q[0], q[1], q[2]);
  }
  const warnings = [];
  if (missing) warnings.push(`${missing} panel(s) reference node numbers that are not in the node table and were skipped.`);
  if (!panels.n) throw new Error(`${fname}: no panel rows ("n1 n2 n3 n4") follow the node table.`);
  const units = resolveUnits(opts, 'm', 'm');
  // isym = 1: the stored half is mirrored about the xOz plane (y → −y)
  return buildPanels(baseName(input.name), panels.finish(), { x: false, y: hp[1] === 1 }, opts, units, 'nemoh', { nodes: nodes.size }, warnings);
}

/** HydroStar hull mesh (.hst). */
export function parseHst(input, opts = {}) {
  const cur = lineCursor(input);
  const fname = input.name || 'HST';
  const ids = new Map(); // node identifier → position in the node table
  const xyz = new Grow(Float64Array, 3 * 1024);
  const panels = new Grow(Float64Array, 12 * 1024), panelNo = [];
  const symOf = new Map(), ranges = [], other = {};
  const header = {};
  let fileSym = 0, nbody = null, sawCoords = false, sawPanels = false, ended = false, badNodes = 0, badPanels = 0, missing = 0, dupIds = 0, line;
  const toks = (l) => l.split(/[\s,;]+/).filter(Boolean);
  const block = (endKey, what, row) => {
    for (;;) {
      const l = cur.nextNonBlank();
      if (l === null) throw new Error(`${fname}: the file ends inside the ${what} (no ${endKey} line).`);
      const t = toks(l);
      if (t[0].toUpperCase() === endKey || t[0].toUpperCase() === `${endKey}S`) return;
      row(t.map(Number));
    }
  };
  while (!ended && (line = cur.nextNonBlank()) !== null) {
    const t = toks(line);
    const key = t[0].replace(/:$/, '').toUpperCase();
    if (key === 'COORDINATES') {
      sawCoords = true;
      block('ENDCOORDINATES', 'node table', (v) => {
        const withId = v.length >= 4, o = withId ? 1 : 0;
        if (v.length < 3 || !(v[o] === v[o] && v[o + 1] === v[o + 1] && v[o + 2] === v[o + 2]) || (withId && !Number.isInteger(v[0]))) { badNodes++; return; }
        const id = withId ? v[0] : xyz.n / 3 + 1;
        if (ids.has(id)) dupIds++;
        ids.set(id, xyz.n / 3);
        xyz.push3(v[o], v[o + 1], v[o + 2]);
      });
    } else if (key === 'PANEL' || key === 'PANELS') {
      const type = Number(t[t.length - 1]);
      if (t.length < 3 || t[1].toUpperCase() !== 'TYPE' || (type !== 0 && type !== 1)) throw new Error(`${fname}: "${line.trim()}" — only PANEL TYPE 0 (node numbers) and PANEL TYPE 1 (panel number, then node numbers) are read.`);
      if (!sawCoords) throw new Error(`${fname}: the panel table comes before the node table (COORDINATES … ENDCOORDINATES).`);
      sawPanels = true;
      block('ENDPANEL', 'panel table', (v) => {
        const nn = v.length - type;
        if (nn < 3 || nn > 4 || v.some((q) => !Number.isInteger(q))) { badPanels++; return; }
        const q = [];
        for (let k = type; k < v.length; k++) q.push(ids.get(v[k]));
        if (q.some((i) => i === undefined)) { missing++; return; }
        if (q.length === 3) q.push(q[2]);
        for (const i of q) panels.push3(xyz.a[3 * i], xyz.a[3 * i + 1], xyz.a[3 * i + 2]);
        panelNo.push(type === 1 ? v[0] : panelNo.length + 1);
      });
    } else if (key === 'ENDFILE') ended = true;
    else if (key === 'SYMMETRY') fileSym = Number(t[1]);
    else if (key === 'SYMMETRY_BODY') symOf.set(Number(t[1]), Number(t[2]));
    else if (key === 'NUMPANEL') ranges.push({ body: Number(t[1]), first: Number(t[2]), last: Number(t[3]) });
    else if (key === 'NBBODY' || key === 'NBODY') nbody = Number(t[1]);
    else if (key === 'PROJECT' || key === 'USER' || key === 'USERS') header[key === 'PROJECT' ? 'project' : 'user'] = line.trim().slice(t[0].length).trim();
    else if (/^[A-Z_]+$/.test(key)) tally(other, key);
    else if (!sawCoords && !sawPanels && !Object.keys(header).length && nbody === null) throw new Error(`${fname}: not a HydroStar mesh (expected keyword lines such as PROJECT, NBBODY, SYMMETRY_BODY, then COORDINATES … ENDCOORDINATES and PANEL TYPE 0 … ENDPANEL; found "${line.trim().slice(0, 40)}").`);
    else tally(other, 'line without a keyword');
  }
  if (!sawCoords) throw new Error(`${fname}: not a HydroStar mesh — no COORDINATES … ENDCOORDINATES node table was found.`);
  if (!xyz.n) throw new Error(`${fname}: the node table holds no "id x y z" rows.`);
  if (!sawPanels) throw new Error(`${fname}: the file holds ${xyz.n / 3} nodes but no panel table (PANEL TYPE 0 … ENDPANEL) — it is truncated or holds only points.`);
  if (!panels.n) throw new Error(`${fname}: no panel of the panel table could be read (${badPanels} malformed row(s), ${missing} naming nodes that are not in the node table).`);
  const warnings = [];
  if (!ended) warnings.push('The file has no ENDFILE line: it may be truncated.');
  if (badNodes) warnings.push(`${badNodes} row(s) of the node table without three numeric coordinates were skipped.`);
  if (dupIds) warnings.push(`${dupIds} node number(s) occur twice in the node table; panels use the last definition.`);
  if (badPanels) warnings.push(`${badPanels} row(s) of the panel table that are not three or four node numbers were skipped.`);
  if (missing) warnings.push(`${missing} panel(s) reference node numbers that are not in the node table and were skipped.`);
  const names = Object.keys(other);
  if (names.length) warnings.push(`Not geometry of the hull and not read: ${names.map((k) => `${k}${other[k] > 1 ? ` (×${other[k]})` : ''}`).join(', ')}.`);
  const symFlags = (v, what) => {
    if (![0, 1, 2].includes(v)) { warnings.push(`${what} = ${v} is not a known symmetry code (0, 1 or 2); no symmetry was applied.`); return { x: false, y: false }; }
    return { x: v === 2, y: v >= 1 };
  };
  // bodies: the NUMPANEL ranges where they cover panels that were read, otherwise one body
  const n = panelNo.length;
  let parts = null;
  const bodies = [];
  if (ranges.length > 1 || symOf.size > 1) {
    const pos = new Map(panelNo.map((v, i) => [v, i]));
    parts = [];
    let covered = 0;
    for (const r of ranges) {
      const a = pos.get(r.first), b = pos.get(r.last);
      if (a === undefined || b === undefined || b < a) { warnings.push(`NUMPANEL ${r.body} ${r.first} ${r.last} names panel numbers that are not in the panel table; body ${r.body} was not separated.`); parts = null; break; }
      parts.push({ name: `body ${r.body}`, from: a, to: b + 1, sym: symFlags(symOf.has(r.body) ? symOf.get(r.body) : fileSym, `SYMMETRY_BODY ${r.body}`) });
      bodies.push({ body: r.body, panels: b + 1 - a, symmetry: symOf.has(r.body) ? symOf.get(r.body) : fileSym });
      covered += b + 1 - a;
    }
    if (parts && (covered !== n || parts.some((p, i) => i && p.from < parts[i - 1].to))) { warnings.push('The NUMPANEL ranges do not cover the panel table exactly once; the bodies were not separated.'); parts = null; }
    if (!parts) bodies.length = 0;
  }
  const one = symOf.size ? [...symOf.values()][0] : fileSym;
  if (!parts && symOf.size > 1 && new Set(symOf.values()).size > 1) throw new Error(`${fname}: the bodies carry different symmetry codes (SYMMETRY_BODY) but the panels could not be assigned to them through NUMPANEL, so the mirror images cannot be built. Check the NUMPANEL lines.`);
  const sym = symFlags(one, symOf.size ? 'SYMMETRY_BODY' : 'SYMMETRY');
  const any = parts ? { x: parts.some((p) => p.sym.x), y: parts.some((p) => p.sym.y) } : sym;
  const units = resolveUnits(opts, 'm', 'm');
  return buildPanels(baseName(input.name), panels.finish(), any, opts, units, 'hydrostar', { ...header, nodes: xyz.n / 3, declaredBodies: nbody, bodies: bodies.length ? bodies : undefined, keywordsSkipped: other }, warnings, parts);
}
