// IGES 5.x (fixed 80-column ASCII form), read without a geometry kernel.
//
// Entities read:
//   116 point                         → points
//   110 line                          → frame (end points welded; member.group = level)
//   100 circular arc, 106 copious data (forms 11–13, 63), 112 parametric spline curve,
//   126 rational B-spline curve, 102 composite curve of those        → polyline
//   128 rational B-spline surface (untrimmed, independent)            → surface (tessellated)
//   124 transformation matrices are applied; the global-section unit flag gives the length unit.
// When every curve in the file is closed and lies in one z-plane, the curves are also returned as a section.
//
// Entities that need a geometry kernel — trimmed and bounded surfaces (144, 143), surfaces of revolution and
// tabulated cylinders (120, 122), ruled, offset and parametric-spline surfaces (118, 140, 114), B-rep
// solids (186, 502–514), planes and CSG primitives — are counted and reported in `meta.needsAdapter`.
// Physically dependent entities (the basis surface and boundary curves of a trimmed surface) are not drawn
// on their own. Annotation, views and properties are skipped. The compressed-ASCII and binary forms of IGES
// are not read.

import { lineCursor, baseName, Grow, Welder, finish, surfaceModel, pointsModel, fortranFloat, tally, M4, cleanName, unitAdvice, closeLoop } from './common.js';
import { resolveUnits, suggestUnits, normaliseUnit } from '../units.js';
import { sampleCurve, tessellateSurface, arcPoints } from './nurbs.js';

const UNIT_FLAG = { 1: 'in', 2: 'mm', 4: 'ft', 5: 'mi', 6: 'm', 7: 'km', 8: 'mil', 9: 'um', 10: 'cm' };
const KERNEL = { 108: 'plane', 114: 'parametric spline surface', 118: 'ruled surface', 120: 'surface of revolution', 122: 'tabulated cylinder', 140: 'offset surface', 141: 'boundary', 143: 'bounded surface', 144: 'trimmed surface', 150: 'block', 152: 'wedge', 154: 'cylinder', 156: 'cone', 158: 'sphere', 160: 'torus', 162: 'solid of revolution', 164: 'solid of extrusion', 168: 'ellipsoid', 180: 'boolean tree', 186: 'manifold solid B-rep', 190: 'plane surface', 192: 'cylindrical surface', 194: 'conical surface', 196: 'spherical surface', 198: 'toroidal surface', 502: 'vertex list', 504: 'edge list', 508: 'loop', 510: 'face', 514: 'shell' };

/** Split free-format IGES data into fields, honouring Hollerith strings (nHtext). Stops at the record delimiter. */
export function igesFields(text, pd = ',', rd = ';', start = 0) {
  const out = [];
  let i = start;
  const n = text.length;
  for (;;) {
    while (i < n && text[i] === ' ') i++;
    if (i >= n) { out.push(''); return { fields: out, end: n }; }
    const m = /^(\d+)H/.exec(text.slice(i, i + 12));
    let val;
    if (m) { const len = Number(m[1]); val = text.slice(i + m[0].length, i + m[0].length + len); i += m[0].length + len; while (i < n && text[i] === ' ') i++; } else { let e = i; while (e < n && text[e] !== pd && text[e] !== rd) e++; val = text.slice(i, e).trim(); i = e; }
    out.push(val);
    if (i >= n || text[i] === rd) return { fields: out, end: i + 1 };
    i++;
  }
}

/** IGES. */
export function parseIges(input, opts = {}) {
  const cur = lineCursor(input);
  const fname = input.name || 'IGES';
  const name = baseName(input.name);
  let G = '';
  const D = [], Pl = [], S = [];
  let first = true;
  for (let l = cur.next(); l !== null; l = cur.next()) {
    if (l.length < 73) { if (!l.trim()) continue; l = l.padEnd(80); }
    const sec = l[72];
    if (first) { first = false; if (sec === 'C') throw new Error(`${fname}: compressed-ASCII IGES is not read. Export the file in the standard fixed-format ASCII form.`); if (sec !== 'S') throw new Error(`${fname}: not an IGES file (expected section letter "S" in column 73 of the first record).`); }
    if (sec === 'S') S.push(l.slice(0, 72).trimEnd()); else if (sec === 'G') G += l.slice(0, 72); else if (sec === 'D') D.push(l); else if (sec === 'P') Pl.push(l.slice(0, 64)); else if (sec === 'T') break;
  }
  if (!D.length || !Pl.length) throw new Error(`${fname}: the IGES file has no directory or parameter section (expected records flagged D and P in column 73).`);
  // global section: the first two fields define the delimiters themselves
  let pd = ',', rd = ';';
  {
    let m = /^\s*1H(.)/.exec(G);
    if (m) pd = m[1];
    const rest = G.slice(G.indexOf(pd, m ? m.index + m[0].length : 0) + 1);
    m = /^\s*1H(.)/.exec(rest);
    if (m) rd = m[1];
  }
  const glob = igesFields(G, pd, rd).fields;
  const unitFlag = parseInt(glob[13], 10);
  const unitName = unitFlag === 3 ? glob[14] : UNIT_FLAG[unitFlag];
  const units = resolveUnits(opts, unitName && normaliseUnit(unitName) ? unitName : null, 'mm');
  const s = units.scale;
  const warnings = [];
  if (unitFlag === 3 && !normaliseUnit(unitName)) warnings.push(`IGES unit name "${unitName}" was not recognised; ${units.unit} assumed.`);

  const nEnt = D.length >> 1;
  const de = (seq) => {
    const k = (seq - 1) >> 1;
    if (!(seq > 0) || k >= nEnt) return null;
    const a = D[2 * k], b = D[2 * k + 1];
    const f = (l, i) => l.slice(8 * i, 8 * i + 8).trim();
    const st = f(a, 8).padStart(8, '0');
    return { seq: 2 * k + 1, type: parseInt(f(a, 0), 10), pd: parseInt(f(a, 1), 10), level: parseInt(f(a, 4), 10) || 0, xform: parseInt(f(a, 6), 10) || 0, subord: Number(st.slice(2, 4)), use: Number(st.slice(4, 6)), count: parseInt(f(b, 3), 10) || 1, form: parseInt(f(b, 4), 10) || 0, label: f(b, 7) };
  };
  const params = (d) => {
    if (!(d.pd > 0) || d.pd > Pl.length) return null;
    let t = '';
    for (let i = d.pd - 1; i < Math.min(Pl.length, d.pd - 1 + d.count); i++) t += Pl[i];
    return igesFields(t, pd, rd).fields;
  };
  const N = (v) => { const x = v === undefined || v === '' ? 0 : fortranFloat(v); return x === x ? x : 0; };
  const xfCache = new Map();
  const xform = (seq, depth = 0) => {
    if (!seq) return null;
    if (xfCache.has(seq)) return xfCache.get(seq);
    const d = de(seq);
    const f = d && d.type === 124 ? params(d) : null;
    let M = null;
    if (f) {
      const v = f.slice(1, 13).map(N);
      M = M4.identity();
      M[0] = v[0]; M[4] = v[1]; M[8] = v[2]; M[12] = v[3];
      M[1] = v[4]; M[5] = v[5]; M[9] = v[6]; M[13] = v[7];
      M[2] = v[8]; M[6] = v[9]; M[10] = v[10]; M[14] = v[11];
      if (d.xform && depth < 8) { const P2 = xform(d.xform, depth + 1); if (P2) M = M4.mul(P2, M); }
    }
    xfCache.set(seq, M);
    return M;
  };
  const apply = (M, c) => { if (M) for (let i = 0; i < c.length; i += 3) M4.apply(M, c[i], c[i + 1], c[i + 2], c, i); if (s !== 1) for (let i = 0; i < c.length; i++) c[i] *= s; return c; };

  const hist = {}, skipped = {}, kernel = {};
  /** Curve of entity `seq` in definition space with its own transformation applied (model units). */
  const curve = (seq, depth = 0) => {
    const d = de(seq);
    if (!d || depth > 8) return null;
    const f = params(d);
    if (!f) return null;
    const M = xform(d.xform);
    const tf = (c) => { if (M) for (let i = 0; i < c.length; i += 3) M4.apply(M, c[i], c[i + 1], c[i + 2], c, i); return c; };
    switch (d.type) {
      case 110: return tf(Float64Array.from([N(f[1]), N(f[2]), N(f[3]), N(f[4]), N(f[5]), N(f[6])]));
      case 100: {
        const zt = N(f[1]), cx = N(f[2]), cy = N(f[3]);
        const a0 = Math.atan2(N(f[5]) - cy, N(f[4]) - cx), r = Math.hypot(N(f[4]) - cx, N(f[5]) - cy);
        let sw = Math.atan2(N(f[7]) - cy, N(f[6]) - cx) - a0;
        // counter-clockwise from start to end; coincident end points mean a full circle
        while (sw <= 1e-12) sw += 2 * Math.PI;
        return tf(arcPoints([cx, cy, zt], [1, 0, 0], [0, 1, 0], r, a0, sw));
      }
      case 106: {
        const ip = parseInt(f[1], 10), n = parseInt(f[2], 10);
        if (!(n > 0)) return null;
        const c = new Float64Array(n * 3);
        if (ip === 1) { const zt = N(f[3]); for (let i = 0; i < n; i++) { c[3 * i] = N(f[4 + 2 * i]); c[3 * i + 1] = N(f[5 + 2 * i]); c[3 * i + 2] = zt; } } else { const w = ip === 3 ? 6 : 3; for (let i = 0; i < n; i++) { c[3 * i] = N(f[3 + w * i]); c[3 * i + 1] = N(f[4 + w * i]); c[3 * i + 2] = N(f[5 + w * i]); } }
        return tf(c);
      }
      case 112: {
        const ndim = parseInt(f[3], 10), n = parseInt(f[4], 10);
        if (!(n > 0)) return null;
        const T = [];
        for (let i = 0; i <= n; i++) T.push(N(f[5 + i]));
        const base = 6 + n;
        const out = [];
        for (let i = 0; i < n; i++) {
          const q = (k) => N(f[base + 12 * i + k]);
          const h = T[i + 1] - T[i];
          for (let k = i === 0 ? 0 : 1; k <= 8; k++) { const u = (h * k) / 8; out.push(q(0) + u * (q(1) + u * (q(2) + u * q(3))), q(4) + u * (q(5) + u * (q(6) + u * q(7))), ndim >= 3 ? q(8) + u * (q(9) + u * (q(10) + u * q(11))) : q(8)); }
        }
        return tf(Float64Array.from(out));
      }
      case 126: {
        const K = parseInt(f[1], 10), Mdeg = parseInt(f[2], 10);
        if (!(K >= 1 && Mdeg >= 1)) return null;
        const nk = K + Mdeg + 2;
        const U = [], W = [];
        for (let i = 0; i < nk; i++) U.push(N(f[7 + i]));
        for (let i = 0; i <= K; i++) W.push(N(f[7 + nk + i]));
        const cp = new Float64Array(3 * (K + 1));
        for (let i = 0; i < cp.length; i++) cp[i] = N(f[7 + nk + K + 1 + i]);
        const o = 7 + nk + K + 1 + 3 * (K + 1);
        const pts = sampleCurve(Mdeg, U, cp, parseInt(f[5], 10) === 1 ? null : W, 8, N(f[o]), N(f[o + 1]));
        return pts ? tf(pts) : null;
      }
      case 102: {
        const n = parseInt(f[1], 10);
        const out = [];
        for (let i = 0; i < n; i++) {
          const c = curve(parseInt(f[2 + i], 10), depth + 1);
          if (!c) continue;
          // members are stored head to tail; reverse one that arrives backwards
          let rev = false;
          if (out.length) { const e = out.length; const ds = Math.hypot(out[e - 3] - c[0], out[e - 2] - c[1], out[e - 1] - c[2]), dr = Math.hypot(out[e - 3] - c[c.length - 3], out[e - 2] - c[c.length - 2], out[e - 1] - c[c.length - 1]); rev = dr < ds; }
          const m = c.length / 3;
          for (let q = 0; q < m; q++) { const k = 3 * (rev ? m - 1 - q : q); const e = out.length; if (e && Math.hypot(out[e - 3] - c[k], out[e - 2] - c[k + 1], out[e - 1] - c[k + 2]) < 1e-9) continue; out.push(c[k], c[k + 1], c[k + 2]); }
        }
        return out.length >= 6 ? tf(Float64Array.from(out)) : null;
      }
      default: return null;
    }
  };

  const segs = [], lines = [], pts = [];
  const sc = new Grow(Float64Array, 3072), si = new Grow(Uint32Array, 3072), groups = [];
  const perSpan = Math.max(2, Math.min(16, opts.tessellation || 6));
  for (let k = 0; k < nEnt; k++) {
    const d = de(2 * k + 1);
    if (!d || !(d.type > 0)) continue;
    if (KERNEL[d.type]) { tally(kernel, `${d.type} ${KERNEL[d.type]}`); continue; }
    if (d.type === 124) continue;
    const independent = d.subord === 0 || d.subord === 2;
    if (![100, 102, 106, 110, 112, 116, 126, 128].includes(d.type)) { tally(skipped, `entity ${d.type}`); continue; }
    if (!independent || (d.use !== 0 && d.use !== 6)) { tally(skipped, `${d.type} (dependent or non-geometry use)`); continue; }
    const label = d.label || `${d.type}/${d.seq}`;
    if (d.type === 116) { const f = params(d); if (f) { const c = apply(xform(d.xform), [N(f[1]), N(f[2]), N(f[3])]); pts.push(c[0], c[1], c[2]); tally(hist, '116 point'); } continue; }
    if (d.type === 128) {
      const f = params(d);
      if (!f) continue;
      const K1 = parseInt(f[1], 10), K2 = parseInt(f[2], 10), M1 = parseInt(f[3], 10), M2 = parseInt(f[4], 10);
      const n1 = K1 + M1 + 2, n2 = K2 + M2 + 2, nc = (K1 + 1) * (K2 + 1);
      if (!(K1 >= 1 && K2 >= 1 && M1 >= 1 && M2 >= 1)) { tally(skipped, '128 with an invalid definition'); continue; }
      let o = 10;
      const U = f.slice(o, o + n1).map(N); o += n1;
      const V = f.slice(o, o + n2).map(N); o += n2;
      const W = f.slice(o, o + nc).map(N); o += nc;
      const cp = Float64Array.from(f.slice(o, o + 3 * nc), N); o += 3 * nc;
      const t = tessellateSurface(M1, M2, U, V, cp, parseInt(f[7], 10) === 1 ? null : W, K1, K2, perSpan, [N(f[o]), N(f[o + 1]), N(f[o + 2]), N(f[o + 3])]);
      if (!t) { tally(skipped, '128 with an invalid definition'); continue; }
      apply(xform(d.xform), t.coords);
      const base = sc.n / 3, start = si.n;
      for (const v of t.coords) sc.push(v);
      for (const v of t.indices) si.push(base + v);
      groups.push({ name: label, start, count: si.n - start, level: d.level });
      tally(hist, '128 B-spline surface');
      continue;
    }
    if (d.type === 106 && ![11, 12, 13, 63].includes(d.form)) {
      if ([1, 2, 3].includes(d.form)) { const c = curve(d.seq); if (c) { if (s !== 1) for (let i = 0; i < c.length; i++) c[i] *= s; for (const v of c) pts.push(v); tally(hist, '106 data points', c.length / 3); } } else tally(skipped, `106 form ${d.form} (annotation)`);
      continue;
    }
    const c = curve(d.seq);
    if (!c || c.length < 6) { tally(skipped, `${d.type} with an invalid definition`); continue; }
    if (s !== 1) for (let i = 0; i < c.length; i++) c[i] *= s;
    if (d.type === 110) { segs.push({ c, level: d.level }); tally(hist, '110 line'); } else { lines.push({ c, name: label, level: d.level, type: d.type }); tally(hist, { 100: '100 circular arc', 102: '102 composite curve', 106: '106 copious data', 112: '112 parametric spline', 126: '126 B-spline curve' }[d.type]); }
  }

  const models = [];
  if (si.n) models.push(surfaceModel(name, sc.finish(), si.finish(), groups.length > 1 ? groups : null));
  if (segs.length) {
    let ext = 0;
    for (const sg of segs) for (const v of sg.c) ext = Math.max(ext, Math.abs(v));
    const weld = new Welder(Math.max(ext * 1e-8, 1e-12), segs.length * 2);
    const members = [];
    for (const sg of segs) { const a = weld.add(sg.c[0], sg.c[1], sg.c[2]), b = weld.add(sg.c[3], sg.c[4], sg.c[5]); if (a !== b) members.push({ a, b, group: `level ${sg.level}` }); }
    if (members.length) models.push({ kind: 'frame', name: models.length ? `${name} (lines)` : name, nodes: weld.finish(), members });
  }
  const closedFlat = lines.length > 0 && lines.every((l) => { const c = l.c, n = c.length; if (Math.hypot(c[0] - c[n - 3], c[1] - c[n - 2], c[2] - c[n - 1]) > 1e-9 * (1 + Math.abs(c[0]) + Math.abs(c[1]))) return false; for (let i = 2; i < n; i += 3) if (Math.abs(c[i] - lines[0].c[2]) > 1e-9 * (1 + Math.abs(c[i]))) return false; return true; });
  if (closedFlat && !segs.length) {
    const loops = [];
    for (const l of lines) { const two = new Float64Array((l.c.length / 3) * 2); for (let i = 0; i < l.c.length / 3; i++) { two[2 * i] = l.c[3 * i]; two[2 * i + 1] = l.c[3 * i + 1]; } const cl = closeLoop(two); if (cl) loops.push(cl); }
    if (loops.length) models.push({ kind: 'section', name: `${name} (profiles)`, loops, layers: lines.map((l, i) => ({ name: l.name, loops: [i] })) });
  }
  for (const l of lines) models.push({ kind: 'polyline', name: cleanName(l.name), points: l.c, level: l.level });
  if (pts.length) models.push(pointsModel(`${name} (points)`, Float64Array.from(pts)));

  const nKernel = Object.values(kernel).reduce((a, b) => a + b, 0);
  for (const k in skipped) warnings.push(`Skipped ${skipped[k]} × ${k}.`);
  if (nKernel && models.length) warnings.push(`${nKernel} entit${nKernel === 1 ? 'y needs' : 'ies need'} a geometry kernel and ${nKernel === 1 ? 'was' : 'were'} not read (${Object.keys(kernel).slice(0, 6).join(', ')}${Object.keys(kernel).length > 6 ? ', …' : ''}); only curves, points and untrimmed B-spline surfaces are returned.`);
  const meta = {
    format: 'iges', units, crs: null, elementTypes: hist, skipped: { ...skipped, ...(nKernel ? { 'entities needing a geometry kernel': nKernel } : {}) }, kernelEntities: kernel, needsAdapter: nKernel ? 'occt' : null,
    start: S.join(' ').trim().slice(0, 200), productId: glob[2] || null, sourceSystem: glob[4] || null, directoryEntries: nEnt, unitFlag: Number.isFinite(unitFlag) ? unitFlag : null,
  };
  if (!models.length) {
    if (opts.allowEmpty) return finish([], warnings, meta);
    if (nKernel) throw new Error(`${fname}: the file holds only entities that need a geometry kernel (${Object.keys(kernel).join(', ')}) — adapter "occt". Without it, export the model as STL / OBJ, or as IGES wireframe and untrimmed NURBS surfaces (entities 110, 126, 128).`);
    throw new Error(`${fname}: no geometry entities were found in the IGES directory (expected types 110, 100, 106, 116, 126 or 128).`);
  }
  const r = finish(models, warnings, meta);
  unitAdvice(units, r.meta.bbox, suggestUnits, r.warnings);
  return r;
}
