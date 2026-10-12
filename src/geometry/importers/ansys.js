// ANSYS Mechanical APDL archive (.cdb, written by CDWRITE): blocked NBLOCK / EBLOCK data in the SOLID layout,
// ET element-type table, CMBLOCK components, SECTYPE / SECDATA beam, pipe and shell sections, R real constants
// (shell thickness) and MP / MPDATA elastic constants.
//
// Not read: the unblocked N / EN command form, EBLOCK written without the SOLID key, element types outside the
// table below (contact, target, surface-effect, mass, spring, MPC …), loads and constraints. All are counted
// in `meta.skipped`.

import { lineCursor, baseName, FeMesh, fortranFloat, tally } from './common.js';
import { feResult } from './fe.js';
import { resolveUnits } from '../units.js';
import { circleProps, rectProps, boxProps, iSectionProps } from './sections.js';

// routine number → family
const FAMILY = {};
for (const n of [1, 3, 4, 8, 10, 16, 20, 23, 24, 33, 44, 54, 59, 180, 188, 288]) FAMILY[n] = 'line2';
for (const n of [189, 289, 290]) FAMILY[n] = 'line3';
for (const n of [28, 41, 43, 63, 131, 143, 157, 163, 181, 42, 55, 182, 25, 29, 75, 162]) FAMILY[n] = 'shell4';
for (const n of [91, 93, 99, 132, 150, 281, 77, 78, 82, 83, 88, 183, 223, 230, 121]) FAMILY[n] = 'shell8';
for (const n of [5, 45, 62, 64, 65, 69, 70, 96, 97, 164, 185, 190, 278]) FAMILY[n] = 'solid8';
for (const n of [90, 95, 117, 122, 147, 186, 226, 231, 279]) FAMILY[n] = 'solid20';
for (const n of [87, 92, 98, 123, 127, 148, 168, 187, 227, 232, 237]) FAMILY[n] = 'tet10';
for (const n of [285]) FAMILY[n] = 'tet4';

function parseFormat(l) {
  // "(3i9,6e21.13e3)" → { ni, wi, nf, wf }
  const m = /\((\d+)i(\d+)(?:,(\d+)e(\d+)\.\d+(?:e\d+)?)?\)/i.exec(l);
  if (!m) return null;
  return { ni: Number(m[1]), wi: Number(m[2]), nf: m[3] ? Number(m[3]) : 0, wf: m[4] ? Number(m[4]) : 0 };
}

const uniq = (a) => a.filter((v, i) => a.indexOf(v) === i);

/** ANSYS CDB archive. */
export function parseAnsysCdb(input, opts = {}) {
  const cur = lineCursor(input);
  const name = baseName(input.name);
  const fe = new FeMesh(name);
  const warnings = [];
  const units = resolveUnits(opts, null, 'm');
  const s = units.scale;
  const skipped = {};
  const et = new Map();
  const real = new Map();
  const materials = {};
  const sections = new Map();
  let secId = 0, any = false, unitsCmd = null;

  const addSolid8 = (n, id, tag) => {
    const [I, J, K, L, M, N, O, P] = n;
    if (K === L && M === N && N === O && O === P) fe.addElement('tet4', [I, J, K, M], id, tag);
    else if (M === N && N === O && O === P) fe.addElement('pyr5', [I, J, K, L, M], id, tag);
    else if (K === L && O === P) fe.addElement('wedge6', [I, J, K, M, N, O], id, tag);
    else fe.addElement('hex8', n.slice(0, 8), id, tag);
  };

  for (let line = cur.next(); line !== null; line = cur.next()) {
    const t = line.trim();
    if (!t || t[0] === '!') continue;
    const p = t.split(',').map((x) => x.trim());
    const cmd = p[0].toUpperCase();
    if (cmd === 'NBLOCK') {
      any = true;
      const fmt = parseFormat(cur.next() || '');
      if (!fmt) throw new Error(`${input.name || 'CDB'}: NBLOCK is not followed by a Fortran format line such as (3i9,6e21.13e3).`);
      for (let l = cur.next(); l !== null; l = cur.next()) {
        if (/^\s*-1\b/.test(l) || /^N,/i.test(l)) break;
        const id = parseInt(l.slice(0, fmt.wi), 10);
        if (!(id > 0)) continue;
        const o = fmt.ni * fmt.wi;
        const c = [0, 0, 0];
        // trailing zero coordinates are simply not written
        for (let k = 0; k < 3; k++) { const f = l.slice(o + k * fmt.wf, o + (k + 1) * fmt.wf); if (f.trim()) c[k] = fortranFloat(f); }
        fe.addNode(id, c[0], c[1], c[2]);
      }
    } else if (cmd === 'EBLOCK') {
      any = true;
      const solidKey = (p[2] || '').toUpperCase() === 'SOLID';
      const fmt = parseFormat(cur.next() || '');
      if (!fmt) throw new Error(`${input.name || 'CDB'}: EBLOCK is not followed by a Fortran format line such as (19i9).`);
      const ints = (l) => { const out = []; for (let k = 0; k * fmt.wi < l.length; k++) { const f = l.slice(k * fmt.wi, (k + 1) * fmt.wi).trim(); if (f === '') break; out.push(parseInt(f, 10)); } return out; };
      if (!solidKey) {
        let n = 0;
        for (let l = cur.next(); l !== null; l = cur.next()) { if (/^\s*-1\b/.test(l)) break; n++; }
        tally(skipped, 'EBLOCK without SOLID key (lines)', n);
        warnings.push('An EBLOCK written without the SOLID key was skipped. Re-write the archive with CDWRITE,DB (blocked format) so that element records carry the material, type, section and node count.');
        continue;
      }
      for (let l = cur.next(); l !== null; l = cur.next()) {
        if (/^\s*-1\b/.test(l)) break;
        const v = ints(l);
        if (v.length < 11) continue;
        const [mat, type, rc, sec, , , , , nn, , id] = v;
        let nodes = v.slice(11);
        while (nodes.length < nn) { const more = cur.next(); if (more === null) break; nodes = nodes.concat(ints(more)); }
        nodes = nodes.slice(0, nn);
        const routine = et.get(type);
        const fam = FAMILY[routine];
        if (!fam) { tally(skipped, routine ? `element routine ${routine}` : `element type ${type} (no ET)`); continue; }
        if (fam === 'line2' || fam === 'line3') {
          const tag = sec ? 1e6 + sec : rc ? 2e6 + rc : 3e6 + mat;
          if (!fe.props.has(tag)) fe.setProp(tag, { name: sec ? `section ${sec}` : rc ? `real ${rc}` : `material ${mat}`, material: `material ${mat}` });
          // BEAM188/PIPE288: I, J (+ orientation node K); BEAM189/PIPE289: I, J, mid K (+ orientation L)
          if (fam === 'line3' && nodes.length >= 3) { fe.addElement('line3', [nodes[0], nodes[1], nodes[2]], id, tag); if (nodes[3]) fe.memberExtra.set(id, { orientNode: nodes[3] }); } else if (nodes.length >= 2) { fe.addElement('line2', [nodes[0], nodes[1]], id, tag); if (nodes[2] && fam === 'line2' && [188, 288, 4, 44].includes(routine)) fe.memberExtra.set(id, { orientNode: nodes[2] }); }
        } else if (fam === 'shell4' || fam === 'shell8') {
          const tag = sec ? 1e6 + sec : rc ? 2e6 + rc : 3e6 + mat;
          if (!fe.props.has(tag)) fe.setProp(tag, { name: sec ? `section ${sec}` : rc ? `real ${rc}` : `material ${mat}`, material: `material ${mat}` });
          if (nodes.length >= 8 && fam === 'shell8') {
            if (nodes[2] === nodes[3]) fe.addElement('tri6', [nodes[0], nodes[1], nodes[2], nodes[4], nodes[5], nodes[7]], id, tag); else fe.addElement('quad8', nodes.slice(0, 8), id, tag);
          } else if (nodes.length >= 4) { if (nodes[2] === nodes[3]) fe.addElement('tri3', nodes.slice(0, 3), id, tag); else fe.addElement('quad4', nodes.slice(0, 4), id, tag); } else if (nodes.length === 3) fe.addElement('tri3', nodes, id, tag);
          else tally(skipped, 'incomplete shell record');
        } else {
          const tag = 3e6 + mat;
          if (!fe.props.has(tag)) fe.setProp(tag, { name: `material ${mat}`, material: `material ${mat}` });
          if (fam === 'tet4' && nodes.length >= 4) fe.addElement('tet4', nodes.slice(0, 4), id, tag);
          else if (fam === 'tet10' && nodes.length >= 10) fe.addElement('tet10', nodes.slice(0, 10), id, tag);
          else if (fam === 'solid20' && nodes.length >= 20) {
            const c = nodes.slice(0, 8);
            if (uniq(c).length === 8) fe.addElement('hex20', nodes.slice(0, 20), id, tag);
            else if (c[2] === c[3] && c[4] === c[5] && c[5] === c[6] && c[6] === c[7]) fe.addElement('tet10', [c[0], c[1], c[2], c[4], nodes[8], nodes[9], nodes[11], nodes[16], nodes[17], nodes[18]], id, tag);
            else { addSolid8(c, id, tag); tally(skipped, 'mid-side nodes of degenerate 20-node solids'); }
          } else if (nodes.length >= 8) addSolid8(nodes, id, tag);
          else if (nodes.length >= 4 && fam === 'tet10') fe.addElement('tet4', nodes.slice(0, 4), id, tag);
          else tally(skipped, 'incomplete solid record');
        }
      }
    } else if (cmd === 'ET') {
      const m = /(\d+)\s*$/.exec(p[2] || '');
      if (m) et.set(Number(p[1]), Number(m[1]));
    } else if (cmd === 'CMBLOCK') {
      const cname = p[1], isNode = /^NODE/i.test(p[2] || ''), count = Number((p[3] || '').split('!')[0]);
      const fmt = parseFormat(cur.next() || '') || { ni: 8, wi: 10 };
      const vals = [];
      while (vals.length < count) {
        const l = cur.next();
        if (l === null) break;
        for (let k = 0; k * fmt.wi < l.length && vals.length < count; k++) { const f = l.slice(k * fmt.wi, (k + 1) * fmt.wi).trim(); if (f) vals.push(parseInt(f, 10)); }
      }
      const ids = [];
      // a negative entry closes a range that starts at the previous entry
      for (let i = 0; i < vals.length; i++) { if (vals[i] < 0 && ids.length) { for (let v = ids[ids.length - 1] + 1; v <= -vals[i]; v++) ids.push(v); } else ids.push(vals[i]); }
      if (isNode) fe.addNodeSet(cname, ids); else fe.addElemSet(cname, ids);
    } else if (cmd === 'SECTYPE') {
      secId = Number(p[1]);
      sections.set(secId, { type: (p[2] || '').toUpperCase(), sub: (p[3] || '').toUpperCase(), name: p[4] || '' });
    } else if (cmd === 'SECDATA') {
      const sec = sections.get(secId);
      if (sec && !sec.data) sec.data = p.slice(1).map((x) => (x === '' ? 0 : fortranFloat(x)));
    } else if (cmd === 'R' || cmd === 'RLBLOCK') {
      if (cmd === 'R') real.set(Number(p[1]), p.slice(2).map((x) => (x === '' ? 0 : fortranFloat(x))));
      else tally(skipped, 'RLBLOCK');
    } else if (cmd === 'MP') {
      const lab = (p[1] || '').toUpperCase(), mat = `material ${Number(p[2])}`;
      (materials[mat] = materials[mat] || { units: 'as written in the archive (not converted)' })[lab] = fortranFloat(p[3] || '');
    } else if (cmd === 'MPDATA') {
      const lab = (p[3] || '').toUpperCase(), mat = `material ${Number(p[4])}`;
      (materials[mat] = materials[mat] || { units: 'as written in the archive (not converted)' })[lab] = fortranFloat(p[6] || '');
    } else if (cmd === '/UNITS') unitsCmd = (p[1] || '').toUpperCase();
    else if (/^(D|F|SFE|SF|BF|ACEL|CP|CE|N|EN|E|TYPE|MAT|REAL|SECNUM|ESYS|LOCAL|CSYS)$/.test(cmd)) tally(skipped, cmd);
  }
  if (!any) throw new Error(`${input.name || 'CDB'}: no NBLOCK or EBLOCK was found (expected a blocked archive written by CDWRITE).`);
  for (const m of Object.values(materials)) { if (m.EX !== undefined) m.E = m.EX; if (m.PRXY !== undefined || m.NUXY !== undefined) m.nu = m.PRXY ?? m.NUXY; if (m.DENS !== undefined) m.rho = m.DENS; }
  for (const [id, sec] of sections) {
    const tag = 1e6 + id;
    const d = (sec.data || []).map((v) => v * s);
    const rec = { name: sec.name || `${sec.type} ${sec.sub} ${id}`.replace(/\s+/g, ' ') };
    if (sec.type === 'PIPE' && d[0] > 0 && d[1] > 0) rec.section = { type: 'tube', D: d[0], t: d[1] };
    else if (sec.type === 'BEAM') {
      if (sec.sub === 'CTUBE' && d[1] > d[0] && d[0] > 0) rec.section = { type: 'tube', D: 2 * d[1], t: d[1] - d[0] };
      else if (sec.sub === 'CSOLID' && d[0] > 0) rec.section = { ...circleProps(2 * d[0]), shape: 'circle', D: 2 * d[0] };
      else if (sec.sub === 'RECT' && d[0] > 0 && d[1] > 0) rec.section = { ...rectProps(d[0], d[1]), shape: 'rect', b: d[0], h: d[1] };
      else if (sec.sub === 'HREC' && d.length >= 6) rec.section = { ...boxProps(d[0], d[1], (d[2] + d[3]) / 2, (d[4] + d[5]) / 2), shape: 'box', b: d[0], h: d[1], tw: (d[2] + d[3]) / 2, tf: (d[4] + d[5]) / 2 };
      else if (sec.sub === 'I' && d.length >= 6 && Math.abs(d[0] - d[1]) < 1e-9 * d[2] && Math.abs(d[3] - d[4]) < 1e-9 * d[2]) rec.section = { ...iSectionProps({ h: d[2], b: d[0], tw: d[5], tf: d[3] }), shape: 'I', h: d[2], b: d[0], tw: d[5], tf: d[3] };
      else warnings.push(`SECTYPE ${id} BEAM ${sec.sub}: this profile is kept by name only (no properties computed).`);
    } else if (sec.type === 'SHELL' && d[0] > 0) rec.thickness = d[0];
    fe.setProp(tag, rec);
  }
  for (const [id, r] of real) if (r[0] > 0) fe.setProp(2e6 + id, { thickness: r[0] * s });
  const hint = { SI: 'm', MKS: 'm', UMKS: 'um', CGS: 'cm', MPA: 'mm', BFT: 'ft', BIN: 'in' }[unitsCmd] || null;
  const nSkip = Object.values(skipped).reduce((a, b) => a + b, 0);
  if (nSkip) warnings.push(`${nSkip} record(s) are not geometry this reader keeps (${Object.keys(skipped).slice(0, 10).join(', ')}).`);
  return feResult(fe, 'ansys-cdb', opts, {
    name: input.name || 'CDB', expected: 'expected an NBLOCK', unitHint: hint, warnings, skipped: { records: skipped },
    extra: { materials, elementTypeTable: Object.fromEntries(et), sections: sections.size },
  });
}
