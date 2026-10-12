// IFC (ISO 16739; IFC2x3, IFC4, IFC4.3) read without a geometry kernel, for structural idealisation.
//
// Linear elements — IfcBeam, IfcColumn, IfcMember, IfcPile (and their StandardCase forms) and the
// analysis-model IfcStructuralCurveMember — are returned as one frame: each element's axis is taken from its
// 'Axis' representation, or from the extrusion of its 'Body' representation, placed through the chain of
// IfcLocalPlacements. Kept per member: section (from the swept profile or the material profile set),
// material name, IFC class as group, GlobalId, Name, Tag and the spatial container (storey / facility part).
// IfcAlignment axis curves (polylines, indexed polycurves, composite curves and 4.3 curve segments on
// lines and circles) are returned as polylines.
//
// Not produced here: surfaces and solids of plates, slabs, walls, footings, proxies, fasteners and any other
// product; the registry routes full geometry to the web-ifc adapter. Those products are counted in
// `meta.otherProducts`. Member ends are welded within 0.1 mm only — physical-model axes generally stop short
// of the joints, so connectivity is not implied. Profile offsets within the cross-section (cardinal points)
// and boolean cut-backs do not move the axis.

import { inputBytes, baseName, Welder, finish, tally, M4, V3, cleanName, unitAdvice } from './common.js';
import { resolveUnits, suggestUnits, normaliseUnit, UNIT_SCALE } from '../units.js';
import { P21, Ref, Typed } from './step.js';
import { resolveSection, circleProps, rectProps, boxProps, iSectionProps } from './sections.js';

const LINEAR = ['IFCBEAM', 'IFCBEAMSTANDARDCASE', 'IFCCOLUMN', 'IFCCOLUMNSTANDARDCASE', 'IFCMEMBER', 'IFCMEMBERSTANDARDCASE', 'IFCPILE', 'IFCSTRUCTURALCURVEMEMBER', 'IFCSTRUCTURALCURVEMEMBERVARYING'];
const OTHER = ['IFCPLATE', 'IFCPLATESTANDARDCASE', 'IFCSLAB', 'IFCWALL', 'IFCWALLSTANDARDCASE', 'IFCFOOTING', 'IFCBUILDINGELEMENTPROXY', 'IFCDISCRETEACCESSORY', 'IFCMECHANICALFASTENER', 'IFCFASTENER', 'IFCELEMENTASSEMBLY', 'IFCRAILING', 'IFCSTAIR', 'IFCSTAIRFLIGHT', 'IFCROOF', 'IFCCOVERING', 'IFCFLOWSEGMENT', 'IFCPIPESEGMENT', 'IFCFLOWFITTING', 'IFCREINFORCINGBAR', 'IFCTENDON', 'IFCBUILDINGELEMENTPART', 'IFCCIVILELEMENT', 'IFCGEOGRAPHICELEMENT', 'IFCSTRUCTURALSURFACEMEMBER', 'IFCBEARING', 'IFCDEEPFOUNDATION', 'IFCCAISSONFOUNDATION', 'IFCMOORINGDEVICE', 'IFCNAVIGATIONELEMENT', 'IFCMARINEFACILITY', 'IFCPAVEMENT', 'IFCCOURSE', 'IFCEARTHWORKSFILL', 'IFCEARTHWORKSCUT'];
const PRETTY = { IFCBEAM: 'IfcBeam', IFCBEAMSTANDARDCASE: 'IfcBeam', IFCCOLUMN: 'IfcColumn', IFCCOLUMNSTANDARDCASE: 'IfcColumn', IFCMEMBER: 'IfcMember', IFCMEMBERSTANDARDCASE: 'IfcMember', IFCPILE: 'IfcPile', IFCSTRUCTURALCURVEMEMBER: 'IfcStructuralCurveMember', IFCSTRUCTURALCURVEMEMBERVARYING: 'IfcStructuralCurveMember' };

/** IFC SPF file → frame of linear elements (+ alignment polylines). */
export function parseIfc(input, opts = {}) {
  const fname = input.name || 'IFC';
  const p = new P21(inputBytes(input), fname);
  if (!p.header.schema.some((x) => /^IFC/i.test(x))) throw new Error(`${fname}: not an IFC file (FILE_SCHEMA is ${p.header.schema.join(', ') || 'missing'}; expected IFC2X3, IFC4 or IFC4X3).`);
  const name = baseName(input.name);
  const warnings = [];
  const hist = {}, skipped = {};
  const num = (v) => (v instanceof Typed ? Number(v.value) : typeof v === 'number' ? v : NaN);
  const str = (v) => (v instanceof Typed ? String(v.value ?? '') : typeof v === 'string' && !/^\.[A-Z0-9_]+\.$/.test(v) ? v : '');

  // --- length unit ---
  let lu = null;
  const SI = { '.KILO.': 1e3, '.HECTO.': 100, '.DECA.': 10, '.DECI.': 0.1, '.CENTI.': 0.01, '.MILLI.': 1e-3, '.MICRO.': 1e-6 };
  const siLen = (e) => (e && e.type === 'IFCSIUNIT' && e.args[1] === '.LENGTHUNIT.' ? { scale: SI[e.args[2]] || 1 } : null);
  for (const id of p.idsOf('IFCUNITASSIGNMENT')) {
    for (const r of p.get(id).args[0] || []) {
      const e = r instanceof Ref ? p.get(r.id) : null;
      if (!e) continue;
      const si = siLen(e);
      if (si) { lu = { name: { 1: 'm', 0.001: 'mm', 0.01: 'cm' }[si.scale] || 'custom', scale: si.scale }; break; }
      if (e.type === 'IFCCONVERSIONBASEDUNIT' && e.args[1] === '.LENGTHUNIT.') {
        const u = normaliseUnit(str(e.args[2]));
        const m = p.get(e.args[3]);
        const base = m ? siLen(p.get(m.args[1])) : null;
        if (m && base && num(m.args[0]) > 0) lu = { name: u || str(e.args[2]) || 'custom', scale: num(m.args[0]) * base.scale };
        else if (u) lu = { name: u, scale: UNIT_SCALE[u] };
        if (lu) break;
      }
    }
    if (lu) break;
  }
  const units = resolveUnits(opts, lu ? (normaliseUnit(lu.name) ? lu.name : null) : null, 'm');
  if (lu && !opts.units && !opts.scale) { units.scale = lu.scale; units.unit = normaliseUnit(lu.name) || lu.name; units.basis = 'file'; }
  if (!lu) warnings.push('No length unit was found in IfcUnitAssignment; metres assumed.');
  const s = units.scale;

  // --- placements ---
  const pt = (r) => { const e = r instanceof Ref ? p.get(r.id) : null; if (!e) return null; const c = e.args[e.args.length - 1]; return Array.isArray(c) ? [Number(c[0]) || 0, Number(c[1]) || 0, Number(c[2]) || 0] : null; };
  const axis3 = (r) => {
    const e = r instanceof Ref ? p.get(r.id) : null;
    if (!e) return M4.identity();
    const o = pt(e.args[0]) || [0, 0, 0];
    if (e.type === 'IFCAXIS2PLACEMENT2D') { const x = V3.norm(pt(e.args[1]) || [1, 0, 0]); return M4.fromAxes(o, x, [-x[1], x[0], 0], [0, 0, 1]); }
    const z = V3.norm(pt(e.args[1]) || [0, 0, 1]);
    let x = pt(e.args[2]) || (Math.abs(z[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]);
    const d = V3.dot(x, z);
    x = V3.norm([x[0] - d * z[0], x[1] - d * z[1], x[2] - d * z[2]]);
    return M4.fromAxes(o, x, V3.cross(z, x), z);
  };
  const plCache = new Map();
  let odd = 0;
  const placement = (r, depth = 0) => {
    if (!(r instanceof Ref) || depth > 64) return M4.identity();
    const c = plCache.get(r.id);
    if (c) return c;
    const e = p.get(r.id);
    let M;
    if (e && e.type === 'IFCLOCALPLACEMENT') M = M4.mul(placement(e.args[0], depth + 1), axis3(e.args[1]));
    else { M = M4.identity(); if (e) odd++; }
    plCache.set(r.id, M);
    return M;
  };

  // --- relations ---
  const materialOf = new Map(), typeOf = new Map(), containerOf = new Map();
  for (const id of p.idsOf('IFCRELASSOCIATESMATERIAL')) { const e = p.get(id); for (const r of e.args[4] || []) if (r instanceof Ref) materialOf.set(r.id, e.args[5]); }
  for (const id of p.idsOf('IFCRELDEFINESBYTYPE')) { const e = p.get(id); for (const r of e.args[4] || []) if (r instanceof Ref) typeOf.set(r.id, e.args[5]); }
  for (const id of p.idsOf('IFCRELCONTAINEDINSPATIALSTRUCTURE')) { const e = p.get(id); const st = p.get(e.args[5]); const nm = st ? cleanName(str(st.args[2]) || str(st.args[7]), st.type) : ''; for (const r of e.args[4] || []) if (r instanceof Ref) containerOf.set(r.id, nm); }
  const unresolved = new Set();
  const profile = (r) => {
    const e = r instanceof Ref ? p.get(r.id) : null;
    if (!e) return null;
    const a = e.args, nm = str(a[1]);
    const n = (k) => num(a[k]) * s;
    let sec = null;
    switch (e.type) {
      case 'IFCCIRCLEHOLLOWPROFILEDEF': if (n(3) > 0 && n(4) > 0) sec = { type: 'tube', D: 2 * n(3), t: n(4) }; break;
      case 'IFCCIRCLEPROFILEDEF': if (n(3) > 0) sec = { ...circleProps(2 * n(3)), shape: 'circle', D: 2 * n(3) }; break;
      case 'IFCRECTANGLEPROFILEDEF': case 'IFCROUNDEDRECTANGLEPROFILEDEF': if (n(3) > 0 && n(4) > 0) sec = { ...rectProps(n(3), n(4)), shape: 'rect', b: n(3), h: n(4) }; break;
      case 'IFCRECTANGLEHOLLOWPROFILEDEF': if (n(3) > 0 && n(4) > 0 && n(5) > 0) sec = { ...boxProps(n(3), n(4), n(5), n(5)), shape: 'box', b: n(3), h: n(4), tw: n(5), tf: n(5) }; break;
      case 'IFCISHAPEPROFILEDEF': if (n(3) > 0 && n(4) > 0) sec = { ...iSectionProps({ h: n(4), b: n(3), tw: n(5), tf: n(6), r: n(7) > 0 ? n(7) : 0 }), shape: 'I', h: n(4), b: n(3), tw: n(5), tf: n(6) }; break;
      default: break;
    }
    if (sec) { if (nm) sec.name = nm; return { section: sec }; }
    if (nm) { const res = resolveSection(nm, { dimUnit: 'mm' }); if (res.matched) return { section: { ...res.section, name: nm } }; unresolved.add(nm); return { sectionName: nm }; }
    tally(skipped, `profile ${e.type} without a name`);
    return null;
  };
  const material = (r, depth = 0) => {
    const e = r instanceof Ref ? p.get(r.id) : null;
    if (!e || depth > 6) return {};
    switch (e.type) {
      case 'IFCMATERIAL': return { material: str(e.args[0]) };
      case 'IFCMATERIALPROFILESETUSAGE': case 'IFCMATERIALPROFILESETUSAGETAPERING': case 'IFCMATERIALLAYERSETUSAGE': return material(e.args[0], depth + 1);
      case 'IFCMATERIALPROFILESET': { const f = (e.args[2] || [])[0]; return material(f, depth + 1); }
      case 'IFCMATERIALPROFILE': return { ...material(e.args[2], depth + 1), ...(profile(e.args[3]) || {}) };
      case 'IFCMATERIALLAYERSET': { const f = (e.args[0] || [])[0]; return material(f, depth + 1); }
      case 'IFCMATERIALLAYER': return material(e.args[0], depth + 1);
      case 'IFCMATERIALLIST': return material((e.args[0] || [])[0], depth + 1);
      case 'IFCMATERIALCONSTITUENTSET': return material((e.args[2] || [])[0], depth + 1);
      case 'IFCMATERIALCONSTITUENT': return material(e.args[2], depth + 1);
      default: return {};
    }
  };

  // --- curves in local coordinates → list of [x, y, z] ---
  let approxSegs = 0;
  const curve = (r, depth = 0) => {
    const e = r instanceof Ref ? p.get(r.id) : null;
    if (!e || depth > 8) return null;
    switch (e.type) {
      case 'IFCPOLYLINE': return (e.args[0] || []).map(pt).filter(Boolean);
      case 'IFCINDEXEDPOLYCURVE': { const l = p.get(e.args[0]); return l ? (l.args[0] || []).map((c) => [Number(c[0]) || 0, Number(c[1]) || 0, Number(c[2]) || 0]) : null; }
      case 'IFCTRIMMEDCURVE': {
        const b = p.get(e.args[0]);
        if (!b || b.type !== 'IFCLINE') { const c = curve(e.args[0], depth + 1); return c; }
        const o = pt(b.args[0]), v = p.get(b.args[1]);
        const d = v ? V3.norm(pt(v.args[0]) || [1, 0, 0]) : [1, 0, 0], mag = v ? num(v.args[1]) || 1 : 1;
        const at = (t) => { for (const x of t || []) if (x instanceof Ref) return pt(x); for (const x of t || []) { const q = num(x); if (q === q) return [o[0] + d[0] * mag * q, o[1] + d[1] * mag * q, o[2] + d[2] * mag * q]; } return null; };
        const a = at(e.args[1]), c = at(e.args[2]);
        return a && c ? [a, c] : null;
      }
      case 'IFCCOMPOSITECURVE': case 'IFCGRADIENTCURVE': case 'IFCSEGMENTEDREFERENCECURVE': {
        if (e.type !== 'IFCCOMPOSITECURVE' && e.args[2] instanceof Ref) { tally(skipped, `${e.type}: vertical profile not applied (horizontal axis only)`); return curve(e.args[2], depth + 1); }
        const out = [];
        const add = (q) => { const n = out.length; if (n && Math.hypot(out[n - 1][0] - q[0], out[n - 1][1] - q[1], out[n - 1][2] - q[2]) < 1e-9) return; out.push(q); };
        for (const sr of e.args[0] || []) {
          const sg = p.get(sr);
          if (!sg) continue;
          if (sg.type === 'IFCCOMPOSITECURVESEGMENT') { let c = curve(sg.args[2], depth + 1); if (c) { if (sg.args[1] === '.F.') c = c.slice().reverse(); c.forEach(add); } } else if (sg.type === 'IFCCURVESEGMENT') {
            const M = axis3(sg.args[1]);
            const par = p.get(sg.args[4]);
            const L = num(sg.args[3]);
            const o = [M[12], M[13], M[14]], x = [M[0], M[1], M[2]], y = [M[4], M[5], M[6]];
            add(o);
            if (par && par.type === 'IFCLINE') add([o[0] + x[0] * L, o[1] + x[1] * L, o[2] + x[2] * L]);
            else if (par && par.type === 'IFCCIRCLE' && num(par.args[1]) > 0) {
              const R = num(par.args[1]), turn = L >= 0 ? 1 : -1, ang = Math.abs(L) / R, n = Math.max(2, Math.ceil(ang / (Math.PI / 36)));
              for (let k = 1; k <= n; k++) { const t = (ang * k) / n, u = R * Math.sin(t), w = turn * R * (1 - Math.cos(t)); add([o[0] + u * x[0] + w * y[0], o[1] + u * x[1] + w * y[1], o[2] + u * x[2] + w * y[2]]); }
            } else approxSegs++;
          }
        }
        return out;
      }
      default: return null;
    }
  };

  // --- body → extrusion axis ---
  const extrusion = (r, M, depth = 0) => {
    const e = r instanceof Ref ? p.get(r.id) : null;
    if (!e || depth > 8) return null;
    switch (e.type) {
      case 'IFCEXTRUDEDAREASOLID': case 'IFCEXTRUDEDAREASOLIDTAPERED': {
        const T = M4.mul(M, axis3(e.args[1]));
        const d = V3.norm(pt(e.args[2]) || [0, 0, 1]), depthV = num(e.args[3]);
        if (!(depthV > 0)) return null;
        return { a: M4.apply(T, 0, 0, 0), b: M4.apply(T, d[0] * depthV, d[1] * depthV, d[2] * depthV), prof: e.args[0], y: M4.applyDir(T, 0, 1, 0) };
      }
      case 'IFCMAPPEDITEM': {
        const src = p.get(e.args[0]), tgt = p.get(e.args[1]);
        let T = M;
        if (tgt) {
          const o = pt(tgt.args[2]) || [0, 0, 0], sc = num(tgt.args[3]);
          const x = V3.norm(pt(tgt.args[0]) || [1, 0, 0]), yv = V3.norm(pt(tgt.args[1]) || [0, 1, 0]), z = tgt.args[4] instanceof Ref ? V3.norm(pt(tgt.args[4])) : V3.cross(x, yv);
          const k = sc === sc && sc > 0 ? sc : 1;
          T = M4.mul(M, M4.fromAxes(o, V3.scale(x, k), V3.scale(yv, k), V3.scale(z, k)));
        }
        const rep = src ? p.get(src.args[1]) : null;
        for (const it of rep ? rep.args[3] || [] : []) { const x = extrusion(it, T, depth + 1); if (x) return x; }
        return null;
      }
      case 'IFCBOOLEANCLIPPINGRESULT': case 'IFCBOOLEANRESULT': return extrusion(e.args[1], M, depth + 1);
      case 'IFCCSGSOLID': return extrusion(e.args[0], M, depth + 1);
      default: return null;
    }
  };

  const weld = new Welder(1e-4, 1024);
  const members = [];
  const classes = {};
  let noAxis = 0, zero = 0;
  for (const id of p.idsOf(...LINEAR)) {
    const e = p.get(id);
    const a = e.args;
    const M = placement(a[5]);
    const shape = p.get(a[6]);
    const reps = shape ? (shape.args[2] || []).map((r) => p.get(r)).filter(Boolean) : [];
    let pts = null, prof = null, orient = null;
    const body = reps.find((r) => /^body$/i.test(str(r.args[1]))) || reps.find((r) => /SweptSolid|Clipping|MappedRepresentation|AdvancedSweptSolid/i.test(str(r.args[2])));
    let ext = null;
    if (body) for (const it of body.args[3] || []) { ext = extrusion(it, M); if (ext) break; }
    const ax = reps.find((r) => /^axis$/i.test(str(r.args[1])));
    if (ax) { for (const it of ax.args[3] || []) { const c = curve(it); if (c && c.length >= 2) { pts = c.map((q) => M4.apply(M, q[0], q[1], q[2])); break; } } }
    if (!pts) {
      const topo = reps.find((r) => r.type === 'IFCTOPOLOGYREPRESENTATION');
      const edge = topo ? p.get((topo.args[3] || [])[0]) : null;
      if (edge && /^IFCEDGE/.test(edge.type)) { const v0 = p.get(edge.args[0]), v1 = p.get(edge.args[1]); const q0 = v0 && pt(v0.args[0]), q1 = v1 && pt(v1.args[0]); if (q0 && q1) pts = [M4.apply(M, ...q0), M4.apply(M, ...q1)]; }
    }
    if (!pts && ext) pts = [ext.a, ext.b];
    if (ext) { prof = ext.prof; orient = V3.norm(ext.y); }
    if (!pts) { noAxis++; continue; }
    const info = { ...material(materialOf.get(id) || materialOf.get((typeOf.get(id) || {}).id)) };
    if (prof && !info.section) Object.assign(info, profile(prof) || {});
    if (info.section && info.sectionName) delete info.sectionName;
    const cls = PRETTY[e.type];
    tally(classes, cls);
    for (let k = 0; k + 1 < pts.length; k++) {
      const na = weld.add(pts[k][0] * s, pts[k][1] * s, pts[k][2] * s), nb = weld.add(pts[k + 1][0] * s, pts[k + 1][1] * s, pts[k + 1][2] * s);
      if (na === nb) { zero++; continue; }
      const m = { a: na, b: nb, id, group: cls, globalId: str(a[0]), name: str(a[2]) };
      if (info.section) m.section = info.section; else if (info.sectionName) m.sectionName = info.sectionName;
      if (info.material) m.material = info.material;
      if (str(a[7])) m.tag = str(a[7]);
      if (str(a[4])) m.objectType = str(a[4]);
      if (containerOf.has(id)) m.container = containerOf.get(id);
      if (orient) m.orient = orient;
      members.push(m);
    }
  }
  const models = [];
  if (members.length) { models.push({ kind: 'frame', name, nodes: weld.finish(), members }); tally(hist, 'line2', members.length); }

  for (const id of p.idsOf('IFCALIGNMENT')) {
    const e = p.get(id);
    const M = placement(e.args[5]);
    const shape = p.get(e.args[6]);
    let done = false;
    for (const rr of shape ? shape.args[2] || [] : []) {
      const rep = p.get(rr);
      for (const it of rep ? rep.args[3] || [] : []) {
        const c = curve(it);
        if (!c || c.length < 2) continue;
        const pp = new Float64Array(c.length * 3);
        c.forEach((q, i) => { const w = M4.apply(M, q[0], q[1], q[2]); pp[3 * i] = w[0] * s; pp[3 * i + 1] = w[1] * s; pp[3 * i + 2] = w[2] * s; });
        const kp = new Float64Array(c.length);
        for (let i = 1; i < c.length; i++) kp[i] = kp[i - 1] + Math.hypot(pp[3 * i] - pp[3 * i - 3], pp[3 * i + 1] - pp[3 * i - 2], pp[3 * i + 2] - pp[3 * i - 1]);
        models.push({ kind: 'polyline', name: cleanName(str(e.args[2]), `${name} alignment`), points: pp, attrs: { kp }, globalId: str(e.args[0]) });
        tally(hist, 'polyline');
        done = true;
        break;
      }
      if (done) break;
    }
    if (!done) tally(skipped, 'IfcAlignment without a readable axis curve');
  }
  const h = p.histogram();
  const otherProducts = {};
  for (const t of OTHER) if (h[t]) otherProducts[t] = h[t];
  if (noAxis) { skipped.linearElementsWithoutAxis = noAxis; warnings.push(`${noAxis} beam/column/member/pile element(s) have neither an 'Axis' representation nor an extruded 'Body' and could not be idealised.`); }
  if (zero) skipped.zeroLengthMembers = zero;
  if (odd) warnings.push(`${odd} placement(s) are not IfcLocalPlacement (grid or linear placements) and were treated as the identity.`);
  if (approxSegs) warnings.push(`${approxSegs} alignment segment(s) on curves other than lines and circles (clothoids, parabolas …) are drawn as straight chords.`);
  if (unresolved.size) warnings.push(`Profile name(s) not resolved to properties (kept by name only): ${[...unresolved].slice(0, 12).join(', ')}${unresolved.size > 12 ? ' …' : ''}.`);
  const nOther = Object.values(otherProducts).reduce((x, y) => x + y, 0);
  if (nOther && models.length) warnings.push(`${nOther} other product(s) (${Object.keys(otherProducts).slice(0, 6).join(', ')}${Object.keys(otherProducts).length > 6 ? ', …' : ''}) carry surface/solid geometry that this reader does not produce; use the full-geometry IFC import for them.`);
  const proj = p.idsOf('IFCPROJECT')[0];
  const meta = {
    format: 'ifc', units, crs: null, elementTypes: hist, skipped, schema: p.header.schema, originatingSystem: p.header.originatingSystem, instances: p.count,
    project: proj ? str(p.get(proj).args[2]) : null, classes, otherProducts, needsAdapter: nOther ? 'web-ifc' : null,
    containers: [...new Set(containerOf.values())],
  };
  if (!models.length) {
    if (opts.allowEmpty) return finish([], warnings, meta);
    throw new Error(`${fname}: no beams, columns, members, piles, structural curve members or alignments with a readable axis were found${nOther ? ` (the file holds ${nOther} other products whose surface/solid geometry needs the full-geometry IFC reader, adapter "web-ifc")` : ''}.`);
  }
  const r = finish(models, warnings, meta);
  unitAdvice(units, r.meta.bbox, suggestUnits, r.warnings);
  return r;
}
