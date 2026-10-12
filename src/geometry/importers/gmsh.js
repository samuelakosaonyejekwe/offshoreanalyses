// Gmsh MSH, ASCII, file-format versions 2.x and 4.1.
//
// Physical groups are kept: every element is tagged with its (first) physical group, named from
// $PhysicalNames, and each physical group is also returned as an element set. Elements without a physical
// group fall back to their elementary entity tag. Binary MSH files and MSH 4.0 are rejected with the export
// setting to change; post-processing sections ($NodeData …) and periodic/partition data are skipped.

import { lineCursor, baseName, FeMesh, permute, tally } from './common.js';
import { feResult } from './fe.js';

// canonical position k ← native position ORDER[k] (see ELEM in common.js for the canonical ordering)
const HEX20 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15];
const HEX27 = HEX20.concat([20, 21, 22, 23, 24, 25, 26]);
const TET10 = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8];
export const GMSH_TYPES = {
  1: ['line2', 2], 2: ['tri3', 3], 3: ['quad4', 4], 4: ['tet4', 4], 5: ['hex8', 8], 6: ['wedge6', 6], 7: ['pyr5', 5],
  8: ['line3', 3], 9: ['tri6', 6], 10: ['quad9', 9], 11: ['tet10', 10, TET10], 12: ['hex27', 27, HEX27], 13: ['wedge18', 18],
  14: ['pyr14', 14], 15: ['point1', 1], 16: ['quad8', 8], 17: ['hex20', 20, HEX20], 18: ['wedge15', 15], 19: ['pyr13', 13],
};
const DIM_TAG = 1000000;

/** Gmsh MSH 2.x / 4.1 ASCII. */
export function parseGmsh(input, opts = {}) {
  const cur = lineCursor(input);
  const name = baseName(input.name);
  const fe = new FeMesh(name);
  const warnings = [];
  const skippedSections = {};
  const physNames = new Map(); // dim*DIM_TAG + tag → name
  const entityPhys = new Map(); // dim*DIM_TAG + entityTag → [physical tags]
  const physSets = new Map();
  let version = null;
  const need = (what) => { const l = cur.nextNonBlank(); if (l === null) throw new Error(`${input.name || 'MSH'}: file ends inside ${what}.`); return l; };
  const tagFor = (dim, phys, entity) => {
    const key = dim * DIM_TAG + (phys || entity);
    if (!fe.props.has(key)) fe.setProp(key, { name: phys ? physNames.get(dim * DIM_TAG + phys) || String(phys) : `entity ${entity}` });
    return key;
  };
  const addEl = (type, nodes, id, dim, physList, entity) => {
    const def = GMSH_TYPES[type];
    if (!def) { fe.skip(`gmsh element type ${type}`); return; }
    const [canon, , order] = def;
    if (canon === 'point1') { fe.skip('point element'); return; }
    fe.addElement(canon, order ? permute(nodes, order) : nodes, id, tagFor(dim, physList[0] || 0, entity));
    for (const p of physList) { if (!p) continue; const nm = physNames.get(dim * DIM_TAG + p) || `physical ${dim}D ${p}`; let s = physSets.get(nm); if (!s) physSets.set(nm, (s = [])); s.push(id); }
  };
  const DIM_OF = { line2: 1, line3: 1, tri3: 2, tri6: 2, quad4: 2, quad8: 2, quad9: 2, point1: 0 };

  let first = true;
  for (let line = cur.nextNonBlank(); line !== null; line = cur.nextNonBlank()) {
    if (line[0] !== '$') { if (first) throw new Error(`${input.name || 'MSH'}: not a Gmsh mesh (expected the first line to be $MeshFormat).`); continue; }
    first = false;
    const sec = line.slice(1).trim();
    if (sec === 'MeshFormat') {
      const p = need('$MeshFormat').split(/\s+/);
      version = Number(p[0]);
      if (Number(p[1]) !== 0) throw new Error(`${input.name || 'MSH'}: binary Gmsh files are not read. Re-save as ASCII (File → Export → "Version 4 ASCII" or "Version 2 ASCII", or Mesh.Binary = 0).`);
      if (!(version >= 2 && version < 3) && !(version >= 4.1 && version < 5)) throw new Error(`${input.name || 'MSH'}: Gmsh MSH version ${p[0]} is not read (supported: 2.x and 4.1). Re-save with Mesh.MshFileVersion = 4.1 or 2.2.`);
    } else if (sec === 'PhysicalNames') {
      const n = Number(need(sec));
      for (let i = 0; i < n; i++) { const m = /^(\d+)\s+(-?\d+)\s+"(.*)"\s*$/.exec(need(sec)); if (m) physNames.set(Number(m[1]) * DIM_TAG + Number(m[2]), m[3]); }
    } else if (sec === 'Entities' && version >= 4) {
      const c = need(sec).split(/\s+/).map(Number);
      for (let dim = 0; dim < 4; dim++) for (let i = 0; i < c[dim]; i++) {
        const p = need(sec).split(/\s+/).map(Number);
        const o = dim === 0 ? 4 : 7; // tag + point xyz, or tag + bounding box
        const np = p[o];
        if (np > 0) entityPhys.set(dim * DIM_TAG + p[0], p.slice(o + 1, o + 1 + np).map(Math.abs));
      }
    } else if (sec === 'Nodes' || sec === 'NOD') {
      if (version === null) version = 1;
      if (version >= 4) {
        const h = need(sec).split(/\s+/).map(Number);
        for (let b = 0; b < h[0]; b++) {
          const bh = need(sec).split(/\s+/).map(Number);
          const n = bh[3];
          const tags = new Array(n);
          for (let i = 0; i < n; i++) tags[i] = Number(need(sec));
          for (let i = 0; i < n; i++) { const p = need(sec).split(/\s+/); fe.addNode(tags[i], Number(p[0]), Number(p[1]), Number(p[2])); }
        }
      } else {
        const n = Number(need(sec));
        for (let i = 0; i < n; i++) { const p = need(sec).split(/\s+/); fe.addNode(Number(p[0]), Number(p[1]), Number(p[2]), Number(p[3])); }
      }
    } else if (sec === 'Elements' || sec === 'ELM') {
      if (version >= 4) {
        const h = need(sec).split(/\s+/).map(Number);
        for (let b = 0; b < h[0]; b++) {
          const [dim, ent, type, n] = need(sec).split(/\s+/).map(Number);
          const phys = entityPhys.get(dim * DIM_TAG + ent) || [];
          for (let i = 0; i < n; i++) { const p = need(sec).split(/\s+/).map(Number); addEl(type, p.slice(1), p[0], dim, phys, ent); }
        }
      } else {
        const n = Number(need(sec));
        for (let i = 0; i < n; i++) {
          const p = need(sec).split(/\s+/).map(Number);
          const type = p[1];
          let phys, ent, nodes;
          if (version >= 2) { const nt = p[2]; phys = nt > 0 ? p[3] : 0; ent = nt > 1 ? p[4] : 0; nodes = p.slice(3 + nt); } else { phys = p[2]; ent = p[3]; nodes = p.slice(5); }
          const def = GMSH_TYPES[type];
          const dim = def ? (DIM_OF[def[0]] ?? 3) : 3;
          addEl(type, nodes, p[0], dim, [phys], ent);
        }
      }
    } else if (sec.startsWith('End')) {
      // section terminator
    } else {
      // unknown or unused section: read past it
      tally(skippedSections, sec);
      const end = `$End${sec}`;
      for (let l = cur.next(); l !== null; l = cur.next()) if (l.trim() === end) break;
    }
  }
  if (version === null) throw new Error(`${input.name || 'MSH'}: not a Gmsh mesh (expected a $MeshFormat section).`);
  for (const [nm, ids] of physSets) fe.addElemSet(nm, ids);
  for (const k in skippedSections) if (!['Entities', 'PartitionedEntities', 'Periodic', 'GhostElements', 'Parametrizations'].includes(k)) warnings.push(`Section $${k} was skipped (not geometry).`);
  return feResult(fe, 'msh', opts, {
    name: input.name || 'MSH', expected: 'expected a $Nodes section', warnings, skipped: Object.keys(skippedSections).length ? { sections: skippedSections } : {},
    extra: { version, physicalGroups: [...physNames].map(([k, v]) => ({ dim: Math.floor(k / DIM_TAG), tag: k % DIM_TAG, name: v })) },
  });
}
