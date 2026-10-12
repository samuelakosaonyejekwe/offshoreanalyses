// I-DEAS universal file (.unv): datasets 2411 (nodes, double precision), 2412 (elements), 164 (units) and
// the group datasets 2467 / 2477 / 2435 / 2452 (same record layout). Every other dataset (results, coordinate
// systems, properties, legacy single-precision 15 / 781 nodes) is read past and counted in `meta.skipped`.

import { lineCursor, baseName, FeMesh, fortranFloat, tally } from './common.js';
import { feResult } from './fe.js';

// FE descriptor id → [canonical type, node count, canonical-from-native order]
const T = (type, n, order) => [type, n, order];
const TRI6 = [0, 2, 4, 1, 3, 5], QUAD8 = [0, 2, 4, 6, 1, 3, 5, 7];
const UNV_TYPES = {
  11: T('line2', 2), 21: T('line2', 2), 22: T('line2', 2), 23: T('line2', 2), 24: T('line3', 3, [0, 2, 1]),
  41: T('tri3', 3), 51: T('tri3', 3), 61: T('tri3', 3), 74: T('tri3', 3), 81: T('tri3', 3), 91: T('tri3', 3),
  42: T('tri6', 6, TRI6), 52: T('tri6', 6, TRI6), 62: T('tri6', 6, TRI6), 72: T('tri6', 6, TRI6), 82: T('tri6', 6, TRI6), 92: T('tri6', 6, TRI6),
  44: T('quad4', 4), 54: T('quad4', 4), 64: T('quad4', 4), 71: T('quad4', 4), 84: T('quad4', 4), 94: T('quad4', 4),
  45: T('quad8', 8, QUAD8), 55: T('quad8', 8, QUAD8), 65: T('quad8', 8, QUAD8), 75: T('quad8', 8, QUAD8), 85: T('quad8', 8, QUAD8), 95: T('quad8', 8, QUAD8),
  111: T('tet4', 4), 118: T('tet10', 10, [0, 2, 4, 9, 1, 3, 5, 6, 7, 8]),
  112: T('wedge6', 6), 113: T('wedge15', 15, [0, 2, 4, 9, 11, 13, 1, 3, 5, 6, 7, 8, 10, 12, 14]),
  115: T('hex8', 8), 116: T('hex20', 20, [0, 2, 4, 6, 12, 14, 16, 18, 1, 3, 5, 7, 13, 15, 17, 19, 8, 9, 10, 11]),
};
const BEAMS = new Set([11, 21, 22, 23, 24]);
const UNIT_CODE = { 1: 'm', 2: 'ft', 3: 'm', 4: 'ft', 5: 'mm', 6: 'cm', 7: 'in', 8: 'mm' };
const DELIM = /^\s*-1\s*$/;

/** I-DEAS universal file. */
export function parseUnv(input, opts = {}) {
  const cur = lineCursor(input);
  const name = baseName(input.name);
  const fe = new FeMesh(name);
  const warnings = [];
  const skipped = {};
  let unitHint = null, datasets = 0;
  const ints = (l) => l.trim().split(/\s+/).map(Number);
  for (;;) {
    // seek the opening delimiter of the next dataset
    let l = cur.next();
    while (l !== null && !DELIM.test(l)) l = cur.next();
    if (l === null) break;
    const idLine = cur.next();
    if (idLine === null) break;
    const ds = parseInt(idLine, 10);
    if (!(ds >= 0)) continue;
    datasets++;
    if (ds === 2411) {
      for (let h = cur.next(); h !== null && !DELIM.test(h); h = cur.next()) {
        const c = cur.next();
        if (c === null) break;
        const p = c.trim().split(/\s+/);
        fe.addNode(parseInt(h, 10), fortranFloat(p[0]), fortranFloat(p[1] || '0'), fortranFloat(p[2] || '0'));
      }
    } else if (ds === 2412) {
      for (let h = cur.next(); h !== null && !DELIM.test(h); h = cur.next()) {
        const r = ints(h);
        const [label, fed, phys, , , nn] = r;
        if (BEAMS.has(fed)) { const b = cur.next(); if (b === null) break; const on = ints(b)[0]; if (on) fe.memberExtra.set(label, { orientNode: on }); }
        let nodes = [];
        while (nodes.length < nn) { const c = cur.next(); if (c === null || DELIM.test(c)) break; nodes = nodes.concat(ints(c)); }
        const def = UNV_TYPES[fed];
        if (!def || nodes.length < def[1]) { tally(skipped, `element descriptor ${fed}`); continue; }
        if (!fe.props.has(phys)) fe.setProp(phys, { name: `property ${phys}` });
        fe.addElement(def[0], def[2] ? def[2].map((k) => nodes[k]) : nodes, label, phys);
      }
    } else if (ds === 2467 || ds === 2477 || ds === 2435 || ds === 2452) {
      for (let h = cur.next(); h !== null && !DELIM.test(h); h = cur.next()) {
        const r = ints(h);
        const count = r[7] || 0;
        const gname = (cur.next() || '').trim() || `group ${r[0]}`;
        const nodes = [], elems = [];
        let read = 0;
        while (read < count) {
          const c = cur.next();
          if (c === null || DELIM.test(c)) break;
          const v = ints(c);
          for (let i = 0; i + 1 < v.length && read < count; i += 4, read++) { if (v[i] === 7) nodes.push(v[i + 1]); else if (v[i] === 8) elems.push(v[i + 1]); else tally(skipped, `group entity type ${v[i]}`); }
        }
        if (nodes.length) fe.addNodeSet(gname, nodes);
        if (elems.length) fe.addElemSet(gname, elems);
      }
    } else if (ds === 164) {
      const h = cur.next();
      if (h !== null && !DELIM.test(h)) { const code = parseInt(h, 10); unitHint = UNIT_CODE[code] || null; if (!unitHint) warnings.push(`Units dataset 164 uses code ${code} (user defined); the length unit could not be determined from it.`); for (let x = cur.next(); x !== null && !DELIM.test(x); x = cur.next()); }
    } else {
      tally(skipped, `dataset ${ds}`);
      for (let x = cur.next(); x !== null && !DELIM.test(x); x = cur.next());
    }
  }
  if (!datasets) throw new Error(`${input.name || 'UNV'}: no universal-file datasets were found (expected blocks delimited by lines holding "    -1", such as dataset 2411 for nodes).`);
  const other = Object.keys(skipped).filter((k) => k.startsWith('dataset'));
  if (other.length) warnings.push(`Dataset(s) not read: ${other.map((k) => k.slice(8)).join(', ')}.`);
  if (skipped['dataset 15'] || skipped['dataset 781']) warnings.push('The file stores nodes in the legacy dataset 15 / 781; re-export with dataset 2411.');
  return feResult(fe, 'unv', opts, { name: input.name || 'UNV', expected: 'expected dataset 2411', unitHint, warnings, skipped, extra: { datasets } });
}
