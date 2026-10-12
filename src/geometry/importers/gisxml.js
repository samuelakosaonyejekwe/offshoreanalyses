// XML-based civil and geospatial formats: LandXML 1.x, GML 2/3 (including CityGML and other application
// schemas that embed gml geometry), KML / KMZ and GPX.
//
// LandXML: TIN surfaces (Pnts / Faces, invisible faces dropped), breaklines, CgPoints, alignments (Line, Curve,
// Spiral, IrregularLine; Profile PVIs with parabolic vertical curves give the elevations), plan features and
// pipe networks (structures + circular pipes → frame). LandXML writes northing before easting; models are
// returned x = easting, y = northing. Spirals are drawn as a quadratic curve tangent at both ends (their
// stationing uses the stated length). Parcels, cross-sections, roadway design and survey data are not read.
//
// GML: every gml:Point, LineString / Curve, Polygon / Surface patch and Triangle found in the document,
// whatever the application schema. 3-D polygons inside solids and (multi/composite/triangulated) surfaces are
// triangulated into one surface per thematic element (WallSurface, RoofSurface, TINRelief …). srsName and
// srsDimension are honoured, including latitude-first axis order of URN-style geographic CRS names.
// Feature attributes, appearance, implicit geometries and xlink references are not read.

import { inputText, baseName, finish, findFile, decodeText, cleanName, FeMesh, tally, extOf, surfaceModel, pointsModel } from './common.js';
import { resolveUnits, normaliseUnit, epsgInfo } from '../units.js';
import { parseXml, kids, kid, attr, numAttr, desc, firstDesc, numberList, kidText } from '../xml.js';
import { isZip, readZip } from '../zip.js';
import { Features, crsFromName, patchesToMesh } from './vector.js';

// ---------------------------------------------------------------------------------------------------------
// LandXML
// ---------------------------------------------------------------------------------------------------------

/** Elevation along a profile defined by PVIs with optional symmetric vertical curves of length L. */
export function profileElevation(pvis, s) {
  const n = pvis.length;
  if (!n) return 0;
  if (n === 1) return pvis[0].z;
  let i = 0;
  while (i < n - 2 && s > pvis[i + 1].s) i++;
  const a = pvis[i], b = pvis[i + 1];
  const g = (b.z - a.z) / (b.s - a.s || 1);
  let z = a.z + g * (s - a.s);
  // within a vertical curve centred on a or b the grade changes linearly: add the parabolic offset
  for (const k of [i, i + 1]) {
    const p = pvis[k];
    if (!(p.L > 0) || k === 0 || k === n - 1) continue;
    const g1 = (p.z - pvis[k - 1].z) / (p.s - pvis[k - 1].s || 1), g2 = (pvis[k + 1].z - p.z) / (pvis[k + 1].s - p.s || 1);
    const s0 = p.s - p.L / 2, s1 = p.s + p.L / 2;
    if (s <= s0 || s >= s1) continue;
    const r = (g2 - g1) / (2 * p.L);
    z += k === i + 1 ? r * (s - s0) ** 2 : r * (s1 - s) ** 2;
  }
  return z;
}

/** LandXML. */
export function parseLandXml(input, opts = {}) {
  const root = parseXml(inputText(input));
  if (root.name !== 'LandXML') throw new Error(`${input.name || 'LandXML'}: not a LandXML document (expected root element <LandXML>, found <${root.tag}>).`);
  const name = baseName(input.name);
  const warnings = [];
  const hist = {}, skipped = {};
  const un = kid(root, 'Units');
  const ue = un && (kid(un, 'Metric') || kid(un, 'Imperial'));
  const unitHint = ue ? attr(ue, 'linearUnit', null) : null;
  const units = resolveUnits(opts, unitHint && normaliseUnit(unitHint) ? unitHint : null, 'm');
  const s = units.scale;
  const cs = kid(root, 'CoordinateSystem');
  let crs = null;
  if (cs) {
    const code = attr(cs, 'epsgCode');
    const info = code ? epsgInfo(Number(code)) : null;
    crs = { source: code ? `EPSG:${code}` : attr(cs, 'name', attr(cs, 'desc', null)), kind: info ? info.kind : 'projected', epsg: code ? Number(code) : null, name: info ? info.name : attr(cs, 'name', null), horizontalDatum: attr(cs, 'horizontalDatum', null), verticalDatum: attr(cs, 'verticalDatum', null), basis: 'file', units: 'm' };
  }
  const models = [];
  // "N E Z" → x = E, y = N
  const cg = new Map();
  const pt = (el) => {
    if (!el) return null;
    const ref = attr(el, 'pntRef');
    if (ref && cg.has(ref)) return cg.get(ref);
    const v = numberList(el.text);
    return v.length >= 2 ? [v[1] * s, v[0] * s, (v[2] || 0) * s] : null;
  };

  const cgPts = desc(root, 'CgPoint');
  if (cgPts.length) {
    const coords = [], names = [], codes = [];
    for (const c of cgPts) {
      const v = numberList(c.text);
      if (v.length < 2) { tally(skipped, 'CgPoint without coordinates'); continue; }
      const p = [v[1] * s, v[0] * s, (v[2] || 0) * s];
      coords.push(...p);
      const nm = attr(c, 'name', '');
      names.push(nm); codes.push(attr(c, 'code', attr(c, 'desc', '')));
      if (nm) cg.set(nm, p);
    }
    if (coords.length) {
      const fe = { kind: 'points', name: `${name} (CgPoints)`, attrs: { names } };
      if (codes.some((c) => c)) fe.attrs.codes = codes;
      models.push({ pts: Float64Array.from(coords), model: fe });
      tally(hist, 'point1', names.length);
    }
  }

  for (const surf of desc(root, 'Surface')) {
    const def = kid(surf, 'Definition');
    const sname = cleanName(attr(surf, 'name'), `${name} surface`);
    if (def) {
      const pn = kid(def, 'Pnts'), fc = kid(def, 'Faces');
      const ids = new Map();
      const P = pn ? kids(pn, 'P') : [];
      const coords = new Float64Array(P.length * 3);
      P.forEach((p, i) => { const v = numberList(p.text); coords[3 * i] = v[1] * s; coords[3 * i + 1] = v[0] * s; coords[3 * i + 2] = (v[2] || 0) * s; ids.set(attr(p, 'id', String(i + 1)), i); });
      const tri = [];
      let invisible = 0, bad = 0;
      for (const f of fc ? kids(fc, 'F') : []) {
        if (attr(f, 'i') === '1') { invisible++; continue; }
        const t = f.text.trim().split(/\s+/).map((k) => ids.get(k));
        if (t.length < 3 || t.some((k) => k === undefined)) { bad++; continue; }
        // LandXML does not fix the winding, and a TIN is a height field: orient every face normal upwards
        const [a, b, c] = t;
        const up = (coords[3 * b] - coords[3 * a]) * (coords[3 * c + 1] - coords[3 * a + 1]) - (coords[3 * b + 1] - coords[3 * a + 1]) * (coords[3 * c] - coords[3 * a]);
        if (up >= 0) tri.push(a, b, c); else tri.push(a, c, b);
      }
      if (invisible) tally(skipped, 'invisible TIN faces', invisible);
      if (bad) tally(skipped, 'TIN faces with unknown point ids', bad);
      if (tri.length) { models.push({ surface: { name: sname, coords, indices: Uint32Array.from(tri), surfType: attr(def, 'surfType', 'TIN') } }); tally(hist, 'tri3', tri.length / 3); } else if (P.length) warnings.push(`Surface "${sname}" has points but no visible faces and was skipped.`);
    }
    for (const bl of desc(surf, 'Breakline')) {
      const l = kid(bl, 'PntList3D');
      const v = l ? numberList(l.text) : [];
      if (v.length < 6) continue;
      const p = new Float64Array(v.length);
      for (let i = 0; i + 2 < v.length; i += 3) { p[i] = v[i + 1] * s; p[i + 1] = v[i] * s; p[i + 2] = v[i + 2] * s; }
      models.push({ polyline: { kind: 'polyline', name: cleanName(attr(bl, 'name'), `${sname} breakline`), points: p } });
      tally(hist, 'polyline');
    }
  }

  const coordGeom = (cgEl) => {
    const xy = [], sta = [];
    let st = 0;
    const push = (p, d) => {
      const n = xy.length / 3;
      if (n && Math.hypot(p[0] - xy[3 * n - 3], p[1] - xy[3 * n - 2]) < 1e-9) { st += d; sta[n - 1] = st; return; }
      st += d;
      xy.push(p[0], p[1], p[2] || 0); sta.push(st);
    };
    for (const e of cgEl.children) {
      const a = pt(kid(e, 'Start')), b = pt(kid(e, 'End'));
      if (e.name === 'Line' && a && b) { push(a, 0); push(b, Math.hypot(b[0] - a[0], b[1] - a[1])); } else if (e.name === 'Curve' && a && b) {
        const c = pt(kid(e, 'Center'));
        if (!c) { push(a, 0); push(b, Math.hypot(b[0] - a[0], b[1] - a[1])); continue; }
        const R = Math.hypot(a[0] - c[0], a[1] - c[1]);
        const a0 = Math.atan2(a[1] - c[1], a[0] - c[0]), a1 = Math.atan2(b[1] - c[1], b[0] - c[0]);
        // rot is the sense seen on the plan (x east, y north)
        const ccw = attr(e, 'rot', 'cw') === 'ccw';
        let sweep = a1 - a0;
        if (ccw) { while (sweep <= 0) sweep += 2 * Math.PI; } else while (sweep >= 0) sweep -= 2 * Math.PI;
        const n = Math.max(2, Math.ceil(Math.abs(sweep) / (Math.PI / 72)));
        push(a, 0);
        for (let k = 1; k <= n; k++) { const t = a0 + (sweep * k) / n; push(k === n ? b : [c[0] + R * Math.cos(t), c[1] + R * Math.sin(t), 0], (R * Math.abs(sweep)) / n); }
      } else if (e.name === 'Spiral' && a && b) {
        const pi = pt(kid(e, 'PI'));
        const L = numAttr(e, 'length', NaN) * s;
        push(a, 0);
        if (!pi) { push(b, L === L ? L : Math.hypot(b[0] - a[0], b[1] - a[1])); continue; }
        const n = 16, q = [];
        let len = 0;
        for (let k = 1; k <= n; k++) { const t = k / n, u = 1 - t; q.push([u * u * a[0] + 2 * u * t * pi[0] + t * t * b[0], u * u * a[1] + 2 * u * t * pi[1] + t * t * b[1], 0]); }
        let prev = a;
        for (const p of q) { len += Math.hypot(p[0] - prev[0], p[1] - prev[1]); prev = p; }
        prev = a;
        for (const p of q) { push(p, (Math.hypot(p[0] - prev[0], p[1] - prev[1]) * (L === L ? L : len)) / len); prev = p; }
        tally(skipped, 'spirals approximated by a quadratic curve');
      } else if (e.name === 'IrregularLine') {
        const l2 = kid(e, 'PntList2D'), l3 = kid(e, 'PntList3D');
        const v = numberList((l3 || l2 || { text: '' }).text), d = l3 ? 3 : 2;
        let prev = null;
        for (let i = 0; i + d - 1 < v.length; i += d) { const p = [v[i + 1] * s, v[i] * s, d === 3 ? v[i + 2] * s : 0]; push(p, prev ? Math.hypot(p[0] - prev[0], p[1] - prev[1]) : 0); prev = p; }
        if (!v.length && a && b) { push(a, 0); push(b, Math.hypot(b[0] - a[0], b[1] - a[1])); }
      } else tally(skipped, `CoordGeom element ${e.name}`);
    }
    return { xy, sta };
  };

  for (const al of desc(root, 'Alignment')) {
    const cgEl = kid(al, 'CoordGeom');
    if (!cgEl) { tally(skipped, 'Alignment without CoordGeom'); continue; }
    const { xy, sta } = coordGeom(cgEl);
    if (xy.length < 6) { tally(skipped, 'Alignment with fewer than two points'); continue; }
    const sta0 = numAttr(al, 'staStart', 0) * s;
    const pa = firstDesc(al, 'ProfAlign');
    const pvis = [];
    if (pa) for (const e of pa.children) { const v = numberList(e.text); if (v.length >= 2) pvis.push({ s: v[0] * s, z: v[1] * s, L: e.name === 'ParaCurve' || e.name === 'CircCurve' || e.name === 'UnsymParaCurve' ? (numAttr(e, 'length', 0) || numAttr(e, 'lengthIn', 0) + numAttr(e, 'lengthOut', 0)) * s : 0 }); }
    let st = sta.map((v) => v + sta0), P = xy;
    if (pvis.length >= 2) {
      // add stations where the profile bends so that the vertical shape is not lost between horizontal vertices
      const extra = [];
      for (const p of pvis) { extra.push(p.s); if (p.L > 0) for (let k = 0; k <= 8; k++) extra.push(p.s - p.L / 2 + (p.L * k) / 8); }
      const all = [...new Set(st.concat(extra.filter((v) => v > st[0] && v < st[st.length - 1])))].sort((x, y) => x - y);
      const out = [];
      let j = 0;
      for (const v of all) {
        while (j < st.length - 2 && st[j + 1] < v) j++;
        const d = st[j + 1] - st[j], f = d > 0 ? Math.min(1, Math.max(0, (v - st[j]) / d)) : 0;
        out.push(xy[3 * j] + f * (xy[3 * j + 3] - xy[3 * j]), xy[3 * j + 1] + f * (xy[3 * j + 4] - xy[3 * j + 1]), profileElevation(pvis, v));
      }
      st = all; P = out;
    }
    models.push({ polyline: { kind: 'polyline', name: cleanName(attr(al, 'name'), `${name} alignment`), points: Float64Array.from(P), attrs: { kp: Float64Array.from(st) }, hasProfile: pvis.length >= 2 } });
    tally(hist, 'polyline');
  }
  for (const pf of desc(root, 'PlanFeature')) {
    const cgEl = kid(pf, 'CoordGeom');
    if (!cgEl) continue;
    const { xy } = coordGeom(cgEl);
    if (xy.length >= 6) { models.push({ polyline: { kind: 'polyline', name: cleanName(attr(pf, 'name'), `${name} plan feature`), points: Float64Array.from(xy) } }); tally(hist, 'polyline'); }
  }

  for (const net of desc(root, 'PipeNetwork')) {
    const fe = new FeMesh(cleanName(attr(net, 'name'), `${name} pipe network`));
    const structs = new Map();
    let id = 0;
    for (const st of desc(net, 'Struct')) {
      const c = pt(kid(st, 'Center'));
      if (!c) continue;
      const inv = {};
      for (const iv of kids(st, 'Invert')) inv[attr(iv, 'refPipe', '')] = numAttr(iv, 'elev', NaN) * s;
      structs.set(attr(st, 'name', String(id + 1)), { c, inv, sump: numAttr(st, 'elevSump', NaN) * s, rim: numAttr(st, 'elevRim', NaN) * s });
    }
    let eid = 0, nid = 0, miss = 0;
    for (const pipe of desc(net, 'Pipe')) {
      const pn = attr(pipe, 'name', String(eid + 1));
      const A = structs.get(attr(pipe, 'refStart')), B = structs.get(attr(pipe, 'refEnd'));
      if (!A || !B) { miss++; continue; }
      const z = (S) => { const v = S.inv[pn]; return v === v && v !== undefined ? v : S.sump === S.sump ? S.sump : S.c[2]; };
      const a = ++nid, b = ++nid;
      fe.addNode(a, A.c[0], A.c[1], z(A)); fe.addNode(b, B.c[0], B.c[1], z(B));
      const cp = kid(pipe, 'CircPipe');
      const tag = ++eid;
      const rec = { name: pn };
      if (cp) { const D = numAttr(cp, 'diameter', 0) * s, t = numAttr(cp, 'thickness', 0) * s; if (D > 0 && t > 0) rec.section = { type: 'tube', D: D + 2 * t, t }; else if (D > 0) rec.innerDiameter = D; }
      fe.setProp(tag, rec);
      fe.addElement('line2', [a, b], eid, tag);
      if (rec.innerDiameter) fe.memberExtra.set(eid, { innerDiameter: rec.innerDiameter });
    }
    if (miss) tally(skipped, 'pipes with unknown structures', miss);
    if (fe.nodeCount) { const b = fe.build({ scale: 1 }); for (const m of b.models) models.push({ ready: m }); tally(hist, 'line2', eid); }
  }
  for (const nm of ['Parcels', 'Roadways', 'GradeModel', 'Survey', 'Monuments', 'Zones']) { const n = desc(root, nm).length; if (n) tally(skipped, `<${nm}>`, n); }

  const out = [];
  for (const m of models) {
    if (m.surface) out.push(surfaceModel(m.surface.name, m.surface.coords, m.surface.indices, null, { surfType: m.surface.surfType }));
    else if (m.model) out.push(pointsModel(m.model.name, m.pts, m.model.attrs));
    else if (m.polyline) out.push(m.polyline);
    else if (m.ready) out.push(m.ready);
  }
  if (!out.length) throw new Error(`${input.name || 'LandXML'}: no geometry was found (expected <Surfaces> with a TIN definition, <Alignments>, <CgPoints> or <PipeNetworks>).`);
  for (const k in skipped) warnings.push(`${skipped[k]} × ${k} not read.`);
  if (opts.vertical === 'depth') warnings.push('LandXML elevations are positive up by definition; vertical = "depth" was ignored.');
  return finish(out, warnings, { format: 'landxml', units, crs, elementTypes: hist, skipped, landxmlVersion: attr(root, 'version', null), application: (() => { const a = kid(root, 'Application'); return a ? attr(a, 'name', null) : null; })() });
}

// ---------------------------------------------------------------------------------------------------------
// GML / CityGML
// ---------------------------------------------------------------------------------------------------------

const SURFACE_PARENTS = new Set(['Solid', 'CompositeSurface', 'MultiSurface', 'TriangulatedSurface', 'Tin', 'Shell', 'CompositeSolid', 'MultiSolid', 'Surface', 'PolyhedralSurface']);
const GEOM = new Set(['Point', 'LineString', 'Curve', 'Polygon', 'PolygonPatch', 'Triangle', 'Rectangle', 'LinearRing', 'Envelope', 'Box']);

/** GML geometry embedded in any application schema (CityGML, InfraGML, WFS feature collections …). */
export function parseGml(input, opts = {}) {
  const text = inputText(input);
  const root = parseXml(text);
  const name = baseName(input.name);
  const F = new Features(name);
  const warnings = [];
  const city = /citygml/i.test(text.slice(0, 4096));
  const crsNames = new Set();
  let crs = null, swapped = false;
  const patches = new Map();
  const tuples = (el, ctx) => {
    // → flat xyz for the coordinates held by a posList / pos / coordinates / coord child
    let dim = ctx.dim, v = null;
    const pl = kid(el, 'posList');
    if (pl) { v = numberList(pl.text); dim = numAttr(pl, 'srsDimension', dim); } else {
      const pos = kids(el, 'pos').concat(kids(el, 'pointProperty').map((p) => firstDesc(p, 'pos')).filter(Boolean));
      if (pos.length) { const parts = pos.map((p) => numberList(p.text)); dim = numAttr(pos[0], 'srsDimension', dim) || parts[0].length; v = new Float64Array(parts.length * dim); parts.forEach((p, i) => v.set(p.subarray(0, dim), i * dim)); } else {
        const co = kid(el, 'coordinates');
        if (co) { const first = co.text.trim().split(/\s+/)[0] || ''; dim = first.split(attr(co, 'cs', ',')).length; v = numberList(co.text); } else {
          const cd = kids(el, 'coord');
          if (cd.length) { dim = kid(cd[0], 'Z') ? 3 : 2; v = new Float64Array(cd.length * dim); cd.forEach((c, i) => { v[i * dim] = +kidText(c, 'X', 0); v[i * dim + 1] = +kidText(c, 'Y', 0); if (dim === 3) v[i * dim + 2] = +kidText(c, 'Z', 0); }); }
        }
      }
    }
    if (!v || !v.length) return null;
    if (!dim) {
      const n = v.length;
      if (city) dim = 3;
      else if (n % 3 === 0 && n % 2 !== 0) dim = 3;
      else if (n % 2 === 0 && n % 3 !== 0) dim = 2;
      else dim = n >= 9 && v[0] === v[n - 3] && v[1] === v[n - 2] && v[2] === v[n - 1] && !(v[0] === v[n - 2] && v[1] === v[n - 1]) ? 3 : 2;
    }
    const n = Math.floor(v.length / dim);
    const out = new Float64Array(n * 3);
    const sw = ctx.swap;
    if (sw) swapped = true;
    for (let i = 0; i < n; i++) { out[3 * i] = v[i * dim + (sw ? 1 : 0)]; out[3 * i + 1] = v[i * dim + (sw ? 0 : 1)]; out[3 * i + 2] = dim >= 3 ? v[i * dim + 2] : 0; }
    return { c: out, dim };
  };
  const ring = (el, ctx) => {
    const lr = el.name === 'LinearRing' ? el : firstDesc(el, 'LinearRing');
    const t = lr ? tuples(lr, ctx) : null;
    if (!t) return null;
    const c = t.c, n = c.length;
    return { c: n >= 6 && c[0] === c[n - 3] && c[1] === c[n - 2] && c[2] === c[n - 1] ? c.subarray(0, n - 3) : c, dim: t.dim };
  };
  const walk = (el, ctx) => {
    const sn = attr(el, 'srsName');
    if (sn) {
      const c = crsFromName(sn);
      crsNames.add(sn);
      if (!crs && c) crs = c;
      const order = opts.axisOrder || (c && c.latFirst ? 'yx' : 'xy');
      ctx = { ...ctx, swap: order === 'yx', dim: numAttr(el, 'srsDimension', ctx.dim) };
    } else if (attr(el, 'srsDimension')) ctx = { ...ctx, dim: numAttr(el, 'srsDimension', ctx.dim) };
    const nm = el.name;
    if (!GEOM.has(nm)) {
      let next = ctx;
      if (SURFACE_PARENTS.has(nm)) next = { ...ctx, solid: true };
      // thematic CityGML / feature elements label the geometry found beneath them
      else if (/^[A-Z]/.test(nm) && el.tag.includes(':') && !el.tag.startsWith('gml:') && el.children.length) { const gn = kids(el, 'name')[0]; next = { ...ctx, label: nm, feature: gn ? gn.text.trim() : attr(el, 'id', ctx.feature) }; }
      for (const c of el.children) walk(c, next);
      return;
    }
    if (nm === 'Envelope' || nm === 'Box') return;
    if (nm === 'Point') { const t = tuples(el, ctx); if (t) F.addPoint(t.c[0], t.c[1], t.c[2], ctx.label ? { type: ctx.label } : null, ctx.feature); return; }
    if (nm === 'LineString' || nm === 'LinearRing') { const t = tuples(el, ctx); if (t) F.addLine(t.c, ctx.label ? { type: ctx.label } : null, ctx.feature); return; }
    if (nm === 'Curve') {
      const parts = [];
      for (const seg of desc(el, ['LineStringSegment', 'GeodesicString', 'Arc', 'ArcString', 'Circle', 'CubicSpline', 'Geodesic'])) { const t = tuples(seg, ctx); if (t) { parts.push(t.c); if (seg.name !== 'LineStringSegment' && seg.name !== 'GeodesicString' && seg.name !== 'Geodesic') tally(F.skipped, `curve segment ${seg.name} drawn through its control points`); } }
      if (parts.length) {
        const all = [];
        for (const p of parts) for (let i = 0; i < p.length; i += 3) { const n = all.length; if (n && all[n - 3] === p[i] && all[n - 2] === p[i + 1] && all[n - 1] === p[i + 2]) continue; all.push(p[i], p[i + 1], p[i + 2]); }
        F.addLine(all, ctx.label ? { type: ctx.label } : null, ctx.feature);
      }
      return;
    }
    // Polygon, PolygonPatch, Triangle, Rectangle
    const ext = kid(el, 'exterior') || kid(el, 'outerBoundaryIs');
    const outer = ext ? ring(ext, ctx) : null;
    if (!outer || outer.c.length < 9) { tally(F.skipped, 'polygon without exterior ring'); return; }
    const rings = [outer.c];
    for (const h of kids(el, 'interior').concat(kids(el, 'innerBoundaryIs'))) { const r = ring(h, ctx); if (r && r.c.length >= 9) rings.push(r.c); }
    if (outer.dim >= 3 && (ctx.solid || nm === 'Triangle' || nm === 'PolygonPatch')) {
      const key = ctx.label || 'surface';
      let list = patches.get(key);
      if (!list) patches.set(key, (list = []));
      list.push(rings);
    } else F.addPolygon(rings, ctx.label ? { type: ctx.label } : null, ctx.feature);
  };
  walk(root, { dim: 0, swap: false, solid: false, label: null, feature: undefined });
  for (const [label, list] of patches) { const m = patchesToMesh(list); if (m.indices.length) F.addMesh(m.coords, m.indices, patches.size > 1 || label !== 'surface' ? label : name); }
  if (crsNames.size > 1) warnings.push(`The document uses ${crsNames.size} different srsName values; all coordinates were treated as ${[...crsNames][0]}.`);
  if (crs && crs.kind === 'unknown') { warnings.push(`Coordinate system ${crs.source} is not one this application can classify; coordinates were read as projected metres.`); crs = { ...crs, kind: 'projected' }; }
  return F.build(opts, { format: 'gml', fileName: input.name, crs, warnIfDegrees: !crs, warnings, expected: `expected gml:Point, gml:LineString, gml:Polygon or surface patches inside <${root.tag}>`, extra: { rootElement: root.tag, axisOrder: swapped ? 'latitude/northing first in the file (swapped to x = east)' : 'x first', cityGml: city } });
}

// ---------------------------------------------------------------------------------------------------------
// KML / KMZ and GPX
// ---------------------------------------------------------------------------------------------------------

function kmlCoords(text) {
  const t = text.trim();
  if (!t) return new Float64Array(0);
  const tup = t.split(/\s+/);
  const out = new Float64Array(tup.length * 3);
  let n = 0;
  for (const s of tup) { const p = s.split(','); const x = Number(p[0]), y = Number(p[1]); if (!(x === x && y === y)) continue; out[3 * n] = x; out[3 * n + 1] = y; out[3 * n + 2] = Number(p[2]) || 0; n++; }
  return out.subarray(0, 3 * n);
}

/** KML, or KMZ (zipped KML). */
export async function parseKml(input, opts = {}) {
  let text;
  if (input.bytes && isZip(input.bytes)) {
    const files = await readZip(input.bytes, { filter: (n) => /\.kml$/i.test(n) });
    const doc = findFile(files, 'doc.kml') || [...files.values()][0];
    if (!doc) throw new Error(`${input.name || 'KMZ'}: the archive holds no .kml document.`);
    text = decodeText(doc);
  } else if (!input.bytes && typeof input.text !== 'string' && input.files) {
    let doc = findFile(input.files, 'doc.kml');
    if (!doc) for (const [k, v] of input.files) if (extOf(k) === 'kml') { doc = v; break; }
    if (!doc) throw new Error(`${input.name || 'KML'}: no .kml document was supplied.`);
    text = decodeText(doc);
  } else text = inputText(input);
  const root = parseXml(text);
  if (root.name !== 'kml' && !firstDesc(root, 'Placemark') && root.name !== 'Placemark') throw new Error(`${input.name || 'KML'}: not a KML document (expected root element <kml> with <Placemark> elements, found <${root.tag}>).`);
  const name = baseName(input.name);
  const F = new Features(name);
  const geom = (g, props, nm) => {
    switch (g.name) {
      case 'Point': { const c = kmlCoords(kidText(g, 'coordinates', '')); if (c.length) F.addPoint(c[0], c[1], c[2], props, nm); break; }
      case 'LineString': case 'LinearRing': F.addLine(kmlCoords(kidText(g, 'coordinates', '')), props, nm); break;
      case 'Polygon': {
        const rings = [];
        for (const bnd of [kid(g, 'outerBoundaryIs'), ...kids(g, 'innerBoundaryIs')]) { const lr = bnd && firstDesc(bnd, 'LinearRing'); if (!lr) continue; let c = kmlCoords(kidText(lr, 'coordinates', '')); const n = c.length; if (n >= 6 && c[0] === c[n - 3] && c[1] === c[n - 2]) c = c.subarray(0, n - 3); rings.push(c); }
        if (rings.length) F.addPolygon(rings, props, nm);
        break;
      }
      case 'MultiGeometry': case 'MultiTrack': for (const s of g.children) geom(s, props, nm); break;
      case 'Track': { const cs = kids(g, 'coord'); const c = new Float64Array(cs.length * 3); cs.forEach((e, i) => { const v = numberList(e.text); c[3 * i] = v[0]; c[3 * i + 1] = v[1]; c[3 * i + 2] = v[2] || 0; }); F.addLine(c, props, nm); break; }
      case 'Model': tally(F.skipped, 'Model (COLLADA reference)'); break;
      default: break;
    }
  };
  const marks = root.name === 'Placemark' ? [root] : desc(root, 'Placemark');
  for (const pm of marks) {
    const nm = kidText(pm, 'name', undefined);
    const props = {};
    const ed = kid(pm, 'ExtendedData');
    if (ed) { for (const d of desc(ed, 'Data')) props[attr(d, 'name', 'data')] = kidText(d, 'value', ''); for (const d of desc(ed, 'SimpleData')) props[attr(d, 'name', 'data')] = d.text.trim(); }
    for (const k in props) { const v = Number(props[k]); if (props[k] !== '' && v === v) props[k] = v; }
    const dsc = kidText(pm, 'description', '');
    if (dsc) props.description = dsc;
    for (const g of pm.children) geom(g, Object.keys(props).length ? props : null, nm);
  }
  const overlays = desc(root, 'GroundOverlay').length + desc(root, 'NetworkLink').length;
  const warnings = [];
  if (overlays) { tally(F.skipped, 'overlays and network links', overlays); }
  return F.build(opts, { format: 'kml', fileName: input.name, crs: { ...crsFromName('EPSG:4326'), basis: 'file' }, warnings, expected: 'expected <Placemark> elements with Point, LineString or Polygon geometry' });
}

/** GPX 1.0 / 1.1: waypoints → points, routes and track segments → polylines. */
export function parseGpx(input, opts = {}) {
  const root = parseXml(inputText(input));
  if (root.name !== 'gpx') throw new Error(`${input.name || 'GPX'}: not a GPX document (expected root element <gpx>, found <${root.tag}>).`);
  const name = baseName(input.name);
  const F = new Features(name);
  const ll = (e) => [numAttr(e, 'lon', NaN), numAttr(e, 'lat', NaN), Number(kidText(e, 'ele', 0)) || 0];
  for (const w of kids(root, 'wpt')) { const p = ll(w); const props = {}; for (const k of ['sym', 'type', 'desc', 'cmt', 'time']) { const v = kidText(w, k, ''); if (v) props[k] = v; } F.addPoint(p[0], p[1], p[2], props, kidText(w, 'name', undefined)); }
  const line = (pts, nm) => { const c = []; for (const e of pts) { const p = ll(e); if (p[0] === p[0] && p[1] === p[1]) c.push(p[0], p[1], p[2]); } F.addLine(c, null, nm); };
  for (const r of kids(root, 'rte')) line(kids(r, 'rtept'), kidText(r, 'name', undefined));
  for (const t of kids(root, 'trk')) { const segs = kids(t, 'trkseg'); const nm = kidText(t, 'name', undefined); segs.forEach((sg, i) => line(kids(sg, 'trkpt'), nm && segs.length > 1 ? `${nm} (${i + 1})` : nm)); }
  return F.build(opts, { format: 'gpx', fileName: input.name, crs: { ...crsFromName('EPSG:4326'), basis: 'file' }, expected: 'expected <wpt>, <rte> or <trk> elements', extra: { creator: attr(root, 'creator', null) } });
}
