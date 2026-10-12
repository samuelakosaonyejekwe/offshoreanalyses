// Hydrodynamic panel meshes for potential-flow analysis: WAMIT low-order GDF and NEMOH mesh files.
//
// Both describe the wetted body surface as quadrilateral panels (triangles as quads with a repeated vertex),
// possibly exploiting planes of symmetry. The result is a triangulated surface with coincident vertices
// welded; the original panels are kept as `quads` (four vertex indices per panel, a triangle repeating its
// last vertex). By default symmetric halves/quarters are mirrored into the full body, each mirror image a
// separate group with reversed winding so that normals keep pointing the same way relative to the body; pass
// `opts.expandSymmetry = false` to keep only what the file stores. `meta.symmetry` records what was declared.
//
// WAMIT: length unit inferred from GRAV (9.80665 → m, 32.174 → ft, 386.09 → in, 9806.65 → mm); ULEN is the
// non-dimensionalising length and does not scale the coordinates here. Higher-order GDF (B-spline patches,
// IGDEF ≠ 0 forms) is not read.

import { lineCursor, baseName, Grow, Welder, finish, surfaceModel, unitAdvice } from './common.js';
import { resolveUnits, suggestUnits } from '../units.js';

function buildPanels(name, panels, sym, opts, units, format, extra, warnings) {
  const s = units.scale;
  const expand = opts.expandSymmetry !== false;
  const images = [[1, 1]];
  if (expand && sym.x) images.push([-1, 1]);
  if (expand && sym.y) { images.push([1, -1]); if (sym.x) images.push([-1, -1]); }
  const weld = new Welder(0, panels.length / 3);
  const tris = new Grow(Uint32Array, panels.length), quads = new Grow(Uint32Array, panels.length / 3);
  const groups = [];
  let nq = 0, nt = 0, bad = 0;
  const n = panels.length / 12;
  for (const [mx, my] of images) {
    const start = tris.n;
    const flip = mx * my < 0;
    for (let p = 0; p < n; p++) {
      const v = [0, 0, 0, 0];
      for (let k = 0; k < 4; k++) { const o = 12 * p + 3 * k; v[k] = weld.add(panels[o] * s * mx, panels[o + 1] * s * my, panels[o + 2] * s); }
      if (flip) { const t = v[1]; v[1] = v[3]; v[3] = t; }
      const u = v.filter((q, k) => v.indexOf(q) === k);
      if (u.length < 3) { bad++; continue; }
      if (u.length === 3) { tris.push(u[0]); tris.push(u[1]); tris.push(u[2]); quads.push(u[0]); quads.push(u[1]); quads.push(u[2]); quads.push(u[2]); nt++; } else { tris.push(v[0]); tris.push(v[1]); tris.push(v[2]); tris.push(v[0]); tris.push(v[2]); tris.push(v[3]); quads.push(v[0]); quads.push(v[1]); quads.push(v[2]); quads.push(v[3]); nq++; }
    }
    if (tris.n > start) groups.push({ name: mx === 1 && my === 1 ? 'as stored' : `mirror ${mx < 0 ? 'x' : ''}${my < 0 ? 'y' : ''}`, start, count: tris.n - start });
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
