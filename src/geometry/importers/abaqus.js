// Abaqus input decks (.inp), flat or part/instance based.
//
// Read: *NODE, *ELEMENT (beams, trusses, pipes, shells, membranes, plane and axisymmetric continuum, solids),
// *NSET / *ELSET (lists, GENERATE, references to other sets), *BEAM SECTION (PIPE, BOX, I, RECT, CIRC),
// *BEAM GENERAL SECTION, *SHELL SECTION, *MEMBRANE SECTION, *SOLID SECTION, *MATERIAL with *ELASTIC and
// *DENSITY, *PART / *INSTANCE (translation and rotation) and *INCLUDE when the file was supplied.
// Everything else (steps, loads, boundary conditions, contact, springs, masses, connectors) is counted in
// `meta.skipped`.
//
// Section axes: Abaqus I11 (bending about the local 1-axis) is returned as `Iy`, I22 as `Iz`.

import { lineCursor, baseName, FeMesh, Grow, findFile, tally, M4, LineCursor } from './common.js';
import { feResult } from './fe.js';
import { resolveUnits } from '../units.js';
import { circleProps, rectProps, boxProps, iSectionProps } from './sections.js';

/** Abaqus element name → [canonical type, node count, canonical-from-native order]. */
export function abaqusType(name) {
  const t = name.toUpperCase();
  let m;
  if (/^(B2[1-3]|B3[13]|PIPE[23]1|T[23]D2|ELBOW31|FRAME[23]D|ITT[23]1)/.test(t)) return ['line2', 2];
  if (/^(B[23]2|PIPE[23]2|T[23]D3|ELBOW32)/.test(t)) return ['line3', 3, [0, 2, 1]];
  if (/^(S3|STRI3|M3D3|R3D3|SFM3D3|DS3|CPS3|CPE3|CPEG3|CAX3|DC2D3|AC2D3|CGAX3|SC6)/.test(t) && !/^SC6/.test(t)) return ['tri3', 3];
  if (/^(S6|STRI65|M3D6|DS6|CPS6|CPE6|CAX6|DC2D6|AC2D6|SFM3D6)/.test(t)) return ['tri6', 6];
  if (/^(S4|M3D4|R3D4|SFM3D4|DS4|CPS4|CPE4|CPEG4|CAX4|DC2D4|AC2D4|CGAX4)/.test(t)) return ['quad4', 4];
  if (/^(S8|M3D8|DS8|CPS8|CPE8|CAX8|DC2D8|AC2D8|SFM3D8)/.test(t)) return ['quad8', 8];
  if (/^(S9|M3D9)/.test(t)) return ['quad9', 9];
  if ((m = /^(?:C3D|DC3D|AC3D|DCC3D|Q3D)(\d+)/.exec(t))) {
    switch (Number(m[1])) {
      case 4: return ['tet4', 4]; case 10: return ['tet10', 10]; case 8: return ['hex8', 8]; case 20: return ['hex20', 20];
      case 27: return ['hex27', 27]; case 6: return ['wedge6', 6]; case 15: return ['wedge15', 15]; case 5: return ['pyr5', 5]; case 13: return ['pyr13', 13];
      default: return null;
    }
  }
  if (/^SC8/.test(t)) return ['hex8', 8];
  if (/^SC6/.test(t)) return ['wedge6', 6];
  return null;
}

function parseKeyword(text) {
  const p = text.slice(1).split(',');
  const name = p[0].trim().toUpperCase().replace(/\s+/g, ' ');
  const params = {};
  for (let i = 1; i < p.length; i++) {
    const s = p[i].trim();
    if (!s) continue;
    const e = s.indexOf('=');
    if (e < 0) params[s.toUpperCase().replace(/\s+/g, ' ')] = true;
    else params[s.slice(0, e).trim().toUpperCase().replace(/\s+/g, ' ')] = s.slice(e + 1).trim().replace(/^"|"$/g, '');
  }
  return { name, params };
}

const newScope = (name) => ({ name, ids: new Grow(Float64Array, 256), xyz: new Grow(Float64Array, 768), blocks: [], nsets: new Map(), elsets: new Map(), sections: [] });
const setKey = (s) => String(s).trim().toUpperCase();

function iGeneral(l, h, b1, b2, t1, t2, t3) {
  // unequal-flange I by parts, about the centroid; y is the axis parallel to the flanges (strong axis)
  const hw = h - t1 - t2;
  const parts = [[b1 * t1, t1 / 2, (b1 * t1 ** 3) / 12, (t1 * b1 ** 3) / 12], [hw * t3, t1 + hw / 2, (t3 * hw ** 3) / 12, (hw * t3 ** 3) / 12], [b2 * t2, h - t2 / 2, (b2 * t2 ** 3) / 12, (t2 * b2 ** 3) / 12]];
  const A = parts.reduce((a, p) => a + p[0], 0);
  const zc = parts.reduce((a, p) => a + p[0] * p[1], 0) / A;
  return { A, Iy: parts.reduce((a, p) => a + p[2] + p[0] * (p[1] - zc) ** 2, 0), Iz: parts.reduce((a, p) => a + p[3], 0), J: (b1 * t1 ** 3 + b2 * t2 ** 3 + (h - (t1 + t2) / 2) * t3 ** 3) / 3 };
}

function beamSection(kind, d, s) {
  const v = d.map((x) => x * s);
  switch (kind) {
    case 'PIPE': return v[0] > 0 && v[1] > 0 && v[1] < v[0] ? { type: 'tube', D: 2 * v[0], t: v[1] } : null;
    case 'CIRC': return v[0] > 0 ? { ...circleProps(2 * v[0]), shape: 'circle', D: 2 * v[0] } : null;
    case 'RECT': return v[0] > 0 && v[1] > 0 ? { ...rectProps(v[0], v[1]), shape: 'rect', b: v[0], h: v[1] } : null;
    case 'BOX': {
      if (!(v.length >= 3)) return null;
      const tw = ((v[2] || 0) + (v[4] || v[2] || 0)) / 2, tf = ((v[3] || v[2]) + (v[5] || v[3] || v[2])) / 2;
      return { ...boxProps(v[0], v[1], tw, tf), shape: 'box', b: v[0], h: v[1], tw, tf };
    }
    case 'I': {
      if (!(v.length >= 7)) return null;
      const [, h, b1, b2, t1, t2, t3] = v;
      if (Math.abs(b1 - b2) < 1e-9 * h && Math.abs(t1 - t2) < 1e-9 * h) return { ...iSectionProps({ h, b: b1, tw: t3, tf: t1 }), shape: 'I', h, b: b1, tw: t3, tf: t1 };
      return { ...iGeneral(v[0], h, b1, b2, t1, t2, t3), shape: 'I', h, b1, b2, t1, t2, tw: t3 };
    }
    default: return null;
  }
}

/** Abaqus input deck. */
export function parseAbaqus(input, opts = {}) {
  const name = baseName(input.name);
  const warnings = [];
  const units = resolveUnits(opts, null, 'm');
  const s = units.scale;
  const skipped = {};
  const materials = {};
  const global = newScope('');
  const parts = new Map();
  const instances = [];
  const asmSets = { n: [], e: [] };
  let scope = global, inAssembly = false, curMat = null, heading = '';
  const stack = [lineCursor(input)];
  let pushed = null;
  const next = () => {
    if (pushed !== null) { const l = pushed; pushed = null; return l; }
    while (stack.length) { const l = stack[stack.length - 1].next(); if (l !== null) return l; stack.pop(); }
    return null;
  };
  /** Data lines up to the next keyword; comment lines dropped. */
  const dataLine = () => {
    for (;;) {
      const l = next();
      if (l === null) return null;
      if (l.startsWith('**')) continue;
      if (l[0] === '*') { pushed = l; return null; }
      if (l.trim()) return l;
    }
  };
  const nums = (l) => l.split(',').map((x) => x.trim()).filter((x) => x !== '').map(Number);
  const addToSet = (map, key, vals) => { let a = map.get(key); if (!a) map.set(key, (a = [])); for (const v of vals) a.push(v); };
  const readSet = (kw, isNode) => {
    const nm = setKey(kw.params[isNode ? 'NSET' : 'ELSET'] || '');
    const items = [];
    for (let l = dataLine(); l !== null; l = dataLine()) {
      if (kw.params.GENERATE) { const [a, b, inc = 1] = nums(l); for (let i = a; i <= b && inc > 0; i += inc) items.push(i); } else for (const tk of l.split(',')) { const t = tk.trim(); if (!t) continue; const v = Number(t); items.push(v === v ? v : setKey(t)); }
    }
    if (!nm) return;
    if (inAssembly && scope === global && kw.params.INSTANCE) (isNode ? asmSets.n : asmSets.e).push({ name: nm, instance: setKey(kw.params.INSTANCE), items });
    else addToSet(isNode ? scope.nsets : scope.elsets, nm, items);
  };

  let any = false;
  for (let line = next(); line !== null; line = next()) {
    if (!line || line.startsWith('**') || line[0] !== '*') continue;
    let text = line;
    while (text.trimEnd().endsWith(',')) { const more = next(); if (more === null) break; text = text.trimEnd() + more; }
    const kw = parseKeyword(text);
    any = true;
    switch (kw.name) {
      case 'HEADING': { const l = dataLine(); if (l !== null) heading = l.trim(); break; }
      case 'NODE': {
        if (kw.params.INPUT) { tally(skipped, '*NODE, INPUT='); break; }
        const first = scope.ids.n;
        for (let l = dataLine(); l !== null; l = dataLine()) { const p = l.split(','); scope.ids.push(Number(p[0])); scope.xyz.push3(Number(p[1]) || 0, Number(p[2]) || 0, Number(p[3]) || 0); }
        if (kw.params.NSET) addToSet(scope.nsets, setKey(kw.params.NSET), scope.ids.view().subarray(first));
        break;
      }
      case 'ELEMENT': {
        const def = abaqusType(kw.params.TYPE || '');
        const elset = kw.params.ELSET ? setKey(kw.params.ELSET) : null;
        if (!def) {
          let n = 0;
          for (let l = dataLine(); l !== null; l = dataLine()) n++;
          tally(skipped, `element type ${kw.params.TYPE || '?'}`, n);
          break;
        }
        const blk = { type: def[0], n: def[1], order: def[2] || null, elset, ids: new Grow(Float64Array, 256), conn: new Grow(Float64Array, 256 * def[1]), src: kw.params.TYPE.toUpperCase() };
        let acc = [];
        for (let l = dataLine(); l !== null; l = dataLine()) {
          for (const v of nums(l)) acc.push(v);
          // long elements (C3D20 …) wrap over several data lines
          if (acc.length >= def[1] + 1) { blk.ids.push(acc[0]); for (let k = 0; k < def[1]; k++) blk.conn.push(acc[1 + (blk.order ? blk.order[k] : k)]); acc = []; }
        }
        if (acc.length) tally(skipped, 'incomplete element record');
        scope.blocks.push(blk);
        if (elset) addToSet(scope.elsets, elset, blk.ids.view());
        break;
      }
      case 'NSET': readSet(kw, true); break;
      case 'ELSET': readSet(kw, false); break;
      case 'BEAM SECTION': case 'BEAM GENERAL SECTION': case 'FRAME SECTION': {
        const kind = (kw.params.SECTION || 'RECT').toUpperCase();
        const l1 = dataLine(), l2 = l1 === null ? null : dataLine();
        for (let l = l2 === null ? null : dataLine(); l !== null; l = dataLine()); // remaining lines (integration points …)
        const d = l1 ? nums(l1) : [];
        let section;
        if (kind === 'GENERAL') section = d.length >= 5 ? { A: d[0] * s * s, Iy: d[1] * s ** 4, Iz: d[3] * s ** 4, J: d[4] * s ** 4 } : null;
        else section = beamSection(kind, d, s);
        const rec = { elset: setKey(kw.params.ELSET || ''), name: `${kind} ${kw.params.ELSET || ''}`.trim(), material: kw.params.MATERIAL || undefined, section };
        if (!section) warnings.push(`*${kw.name}, SECTION=${kind} on ${kw.params.ELSET || '?'}: this profile is kept by name only (no properties computed).`);
        const o = l2 ? nums(l2) : [];
        if (o.length >= 3 && kind !== 'GENERAL') rec.orient = o.slice(0, 3);
        scope.sections.push(rec);
        break;
      }
      case 'SHELL SECTION': case 'SHELL GENERAL SECTION': case 'MEMBRANE SECTION': {
        const l1 = dataLine();
        for (let l = l1 === null ? null : dataLine(); l !== null; l = dataLine());
        const d = l1 ? nums(l1) : [];
        scope.sections.push({ elset: setKey(kw.params.ELSET || ''), name: kw.params.ELSET || '', material: kw.params.MATERIAL || undefined, thickness: d[0] > 0 ? d[0] * s : undefined });
        break;
      }
      case 'SOLID SECTION': {
        for (let l = dataLine(); l !== null; l = dataLine());
        scope.sections.push({ elset: setKey(kw.params.ELSET || ''), name: kw.params.ELSET || '', material: kw.params.MATERIAL || undefined });
        break;
      }
      case 'MATERIAL': curMat = kw.params.NAME || `material ${Object.keys(materials).length + 1}`; materials[curMat] = { units: 'as written in the deck (not converted)' }; break;
      case 'ELASTIC': { const l = dataLine(); if (l && curMat && (!kw.params.TYPE || /^ISO/i.test(kw.params.TYPE))) { const d = nums(l); materials[curMat].E = d[0]; materials[curMat].nu = d[1] ?? 0; } for (let x = dataLine(); x !== null; x = dataLine()); break; }
      case 'DENSITY': { const l = dataLine(); if (l && curMat) materials[curMat].rho = nums(l)[0]; break; }
      case 'PART': scope = newScope(kw.params.NAME || `part ${parts.size + 1}`); parts.set(setKey(scope.name), scope); break;
      case 'END PART': scope = global; break;
      case 'ASSEMBLY': inAssembly = true; scope = global; break;
      case 'END ASSEMBLY': inAssembly = false; break;
      case 'INSTANCE': {
        const inst = { name: kw.params.NAME || `instance ${instances.length + 1}`, part: setKey(kw.params.PART || ''), M: M4.identity() };
        const l1 = dataLine(), l2 = l1 === null ? null : dataLine();
        const t = l1 ? nums(l1) : [];
        if (t.length >= 3) inst.M = M4.translation(t[0], t[1], t[2]);
        const r = l2 ? nums(l2) : [];
        if (r.length >= 7) {
          // rotation by r[6] degrees about the axis from point a (r[0..2]) to point b (r[3..5]), applied after the translation
          const R = M4.mul(M4.translation(r[0], r[1], r[2]), M4.mul(M4.rotation(r[3] - r[0], r[4] - r[1], r[5] - r[2], (r[6] * Math.PI) / 180), M4.translation(-r[0], -r[1], -r[2])));
          inst.M = M4.mul(R, inst.M);
        }
        instances.push(inst);
        // an instance may carry its own mesh instead of referring to a part
        if (!parts.has(inst.part)) { scope = newScope(inst.name); parts.set(setKey(inst.name), scope); inst.part = setKey(inst.name); }
        break;
      }
      case 'END INSTANCE': scope = global; break;
      case 'INCLUDE': {
        const f = kw.params.INPUT && findFile(input.files, kw.params.INPUT);
        if (f) stack.push(new LineCursor(f)); else warnings.push(`*INCLUDE, INPUT=${kw.params.INPUT || '?'}: the included file was not supplied and is missing from the model.`);
        break;
      }
      default: tally(skipped, `*${kw.name}`);
    }
  }
  if (!any) throw new Error(`${input.name || 'INP'}: no Abaqus keyword lines were found (expected lines starting with *, such as *NODE and *ELEMENT).`);

  // ---- emit scopes into one mesh --------------------------------------------------------------------
  const fe = new FeMesh(name);
  let nodeBase = 0, elemBase = 0, tagNo = 0;
  const instanceBase = new Map();
  const multi = (global.ids.n ? 1 : 0) + (instances.length || parts.size) > 1;
  const emit = (sc, M, label) => {
    const nb = nodeBase, eb = elemBase;
    const prefix = multi && label ? `${label}.` : '';
    let maxN = 0, maxE = 0;
    const ids = sc.ids.a, xyz = sc.xyz.a, p = [0, 0, 0];
    for (let i = 0; i < sc.ids.n; i++) {
      if (M) M4.apply(M, xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2], p); else { p[0] = xyz[3 * i]; p[1] = xyz[3 * i + 1]; p[2] = xyz[3 * i + 2]; }
      fe.addNode(ids[i] + nb, p[0], p[1], p[2]);
      if (ids[i] > maxN) maxN = ids[i];
    }
    const resolve = (map, key, seen = new Set()) => {
      const out = [];
      if (seen.has(key)) return out;
      seen.add(key);
      for (const it of map.get(key) || []) { if (typeof it === 'number') out.push(it); else for (const v of resolve(map, it, seen)) out.push(v); }
      return out;
    };
    // element → tag: a section assignment wins over the ELSET named on *ELEMENT
    const tagOf = new Map();
    const orientOf = new Map();
    for (const sec of sc.sections) {
      const tag = ++tagNo;
      const rec = { name: prefix + (sec.name || sec.elset) };
      if (sec.section) rec.section = sec.section;
      if (sec.material) rec.material = sec.material;
      if (sec.thickness !== undefined) rec.thickness = sec.thickness;
      fe.setProp(tag, rec);
      for (const id of resolve(sc.elsets, sec.elset)) { tagOf.set(id, tag); if (sec.orient) orientOf.set(id, sec.orient); }
    }
    const setTag = new Map();
    for (const b of sc.blocks) {
      let fallback = 0;
      if (b.elset) { fallback = setTag.get(b.elset); if (!fallback) { fallback = ++tagNo; setTag.set(b.elset, fallback); fe.setProp(fallback, { name: prefix + b.elset }); } }
      const nodes = new Array(b.n);
      for (let e = 0; e < b.ids.n; e++) {
        const id = b.ids.a[e];
        for (let k = 0; k < b.n; k++) nodes[k] = b.conn.a[e * b.n + k] + nb;
        fe.addElement(b.type, nodes, id + eb, tagOf.get(id) ?? fallback);
        const o = orientOf.get(id);
        if (o) fe.memberExtra.set(id + eb, { orient: o });
        if (id > maxE) maxE = id;
      }
    }
    for (const key of sc.nsets.keys()) fe.addNodeSet(prefix + key, resolve(sc.nsets, key).map((v) => v + nb));
    for (const key of sc.elsets.keys()) fe.addElemSet(prefix + key, resolve(sc.elsets, key).map((v) => v + eb));
    if (multi) { nodeBase += maxN; elemBase += maxE; }
    return { nb, eb };
  };
  if (instances.length) {
    for (const inst of instances) {
      const sc = parts.get(inst.part);
      if (!sc) { warnings.push(`*INSTANCE ${inst.name}: part "${inst.part}" is not defined in the file.`); continue; }
      instanceBase.set(setKey(inst.name), emit(sc, M4.isIdentity(inst.M) ? null : inst.M, inst.name));
    }
  } else for (const sc of parts.values()) emit(sc, null, sc.name);
  if (global.ids.n || global.blocks.length) emit(global, null, '');
  else { for (const [k, v] of global.nsets) if (v.every((x) => typeof x === 'number') && !multi) fe.addNodeSet(k, v); for (const [k, v] of global.elsets) if (v.every((x) => typeof x === 'number') && !multi) fe.addElemSet(k, v); }
  for (const st of asmSets.n) { const b = instanceBase.get(st.instance); if (b) fe.addNodeSet(st.name, st.items.filter((v) => typeof v === 'number').map((v) => v + b.nb)); }
  for (const st of asmSets.e) { const b = instanceBase.get(st.instance); if (b) fe.addElemSet(st.name, st.items.filter((v) => typeof v === 'number').map((v) => v + b.eb)); }

  const nSkip = Object.values(skipped).reduce((a, b) => a + b, 0);
  if (nSkip) warnings.push(`${nSkip} keyword block(s) / record(s) are not geometry this reader keeps (${Object.keys(skipped).slice(0, 10).join(', ')}${Object.keys(skipped).length > 10 ? ', …' : ''}).`);
  return feResult(fe, 'abaqus', opts, {
    name: input.name || 'INP', expected: 'expected a *NODE block', warnings, skipped: { keywords: skipped },
    extra: { heading, materials, parts: [...parts.values()].map((p) => p.name), instances: instances.map((i) => i.name) },
  });
}
