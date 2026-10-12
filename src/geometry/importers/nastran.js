// Nastran bulk data (.bdf / .nas / .dat): small-field, large-field and free-field cards with continuations.
//
// Read: GRID (with CORD2R / CORD2C / CORD2S systems defined in the basic system), CBAR, CBEAM, CROD, CTUBE,
// CONROD, CTRIA3/CTRIAR, CTRIA6, CQUAD4/CQUADR, CQUAD8, CTETRA (4/10), CHEXA (8/20), CPENTA (6/15), CPYRAM (5/13),
// PBAR, PBARL, PBEAM, PBEAML, PROD, PTUBE, PSHELL, PSOLID, MAT1, INCLUDE (when the file was supplied with it).
// Everything else (loads, constraints, rigid elements, sets, case control) is counted in `meta.skipped`.
//
// Section axes: Nastran I1 is the bending inertia in element plane 1 (normally the strong axis). It is returned
// as `Iy` and I2 as `Iz`, matching the strong-axis-is-y convention of sections.js.

import { lineCursor, baseName, FeMesh, fortranFloat, findFile, tally, V3, LineCursor } from './common.js';
import { feResult } from './fe.js';
import { resolveUnits } from '../units.js';
import { circleProps, rectProps, boxProps, iSectionProps } from './sections.js';

function expandTabs(l) {
  if (l.indexOf('\t') < 0) return l;
  let out = '';
  for (const ch of l) { if (ch === '\t') out += ' '.repeat(8 - (out.length % 8)); else out += ch; }
  return out;
}

/** Split one physical line into { name, fields, cont } — `cont` is true for a continuation line. */
function splitCard(raw) {
  let l = expandTabs(raw);
  const d = l.indexOf('$');
  if (d >= 0) l = l.slice(0, d);
  if (!l.trim()) return null;
  if (l.indexOf(',') >= 0) {
    const p = l.split(',').map((s) => s.trim());
    const head = p[0];
    const cont = head === '' || head[0] === '+' || head[0] === '*';
    const large = head.endsWith('*') || head[0] === '*';
    // free field: any number of data fields; a trailing "+name" is the continuation mnemonic, not data
    const f = p.slice(1);
    if (f.length && /^[+*]/.test(f[f.length - 1])) f.pop();
    const w = large ? 4 : 8;
    while (f.length % w || !f.length) f.push('');
    return { name: cont ? '' : head.replace(/\*$/, '').toUpperCase(), fields: f, cont };
  }
  const head = l.slice(0, 8).trim();
  const cont = head === '' || head[0] === '+' || head[0] === '*';
  const large = head.endsWith('*') || head[0] === '*';
  const f = [];
  if (large) for (let i = 0; i < 4; i++) f.push(l.slice(8 + 16 * i, 24 + 16 * i).trim());
  else for (let i = 0; i < 8; i++) f.push(l.slice(8 + 8 * i, 16 + 8 * i).trim());
  return { name: cont ? '' : head.replace(/\*$/, '').toUpperCase(), fields: f, cont };
}

const int = (s) => (s === undefined || s === '' ? 0 : parseInt(s, 10) || 0);
const real = (s) => { const v = s === undefined ? NaN : fortranFloat(s); return v === v ? v : 0; };
const isReal = (s) => s !== undefined && s !== '' && /[.eEdD]/.test(s);

const HEX20 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15];

/** Nastran bulk data. */
export function parseNastran(input, opts = {}) {
  const name = baseName(input.name);
  const fe = new FeMesh(name);
  const warnings = [];
  const units = resolveUnits(opts, null, 'm');
  const s = units.scale;
  let skipped = {};
  const materials = {};
  const cords = new Map();
  const gridCp = [];
  const stack = [lineCursor(input)];
  let sawBulk = false, cards = 0;

  const solid = (c, counts, types, orders) => {
    const f = c.fields;
    const ids = [];
    for (let i = 2; i < f.length; i++) { if (f[i] === '') { if (ids.length >= counts[0]) break; continue; } ids.push(int(f[i])); }
    let k = counts.length - 1;
    while (k > 0 && ids.length < counts[k]) k--;
    if (ids.length < counts[0]) { fe.skip(`${c.name} with too few nodes`); return; }
    if (ids.length > counts[k] && k === 0 && counts.length > 1) warnings.push(`${c.name} ${f[0]}: some mid-side nodes are blank; read as the linear element.`);
    const nodes = ids.slice(0, counts[k]);
    fe.addElement(types[k], orders && orders[k] ? orders[k].map((q) => nodes[q]) : nodes, int(f[0]), int(f[1]) || int(f[0]));
  };

  const handle = (c) => {
    const f = c.fields;
    cards++;
    switch (c.name) {
      case 'GRID': {
        const id = int(f[0]), cp = int(f[1]);
        fe.addNode(id, real(f[2]), real(f[3]), real(f[4]));
        if (cp) gridCp.push(id, cp);
        break;
      }
      case 'CORD2R': case 'CORD2C': case 'CORD2S':
        cords.set(int(f[0]), { type: c.name[5], rid: int(f[1]), a: [real(f[2]), real(f[3]), real(f[4])], b: [real(f[5]), real(f[6]), real(f[7])], c: [real(f[8]), real(f[9]), real(f[10])] });
        break;
      case 'CBAR': case 'CBEAM': {
        const eid = int(f[0]);
        fe.addElement('line2', [int(f[2]), int(f[3])], eid, int(f[1]) || eid);
        const ex = {};
        if (isReal(f[4]) || f[5] || f[6]) ex.orient = [real(f[4]), real(f[5]), real(f[6])]; else if (f[4]) ex.orientNode = int(f[4]);
        if (f[8] || f[9]) ex.releases = { a: f[8] || '', b: f[9] || '' };
        const off = [10, 11, 12, 13, 14, 15].map((i) => real(f[i]) * s);
        if (off.some((v) => v !== 0)) ex.offsets = { a: off.slice(0, 3), b: off.slice(3) };
        if (Object.keys(ex).length) fe.memberExtra.set(eid, ex);
        break;
      }
      case 'CROD': case 'CTUBE': fe.addElement('line2', [int(f[2]), int(f[3])], int(f[0]), int(f[1]) || int(f[0])); break;
      case 'CONROD': {
        const eid = int(f[0]), tag = -eid;
        fe.addElement('line2', [int(f[1]), int(f[2])], eid, tag);
        fe.setProp(tag, { name: `CONROD ${eid}`, section: { A: real(f[4]) * s * s, Iy: 0, Iz: 0, J: real(f[5]) * s ** 4 }, material: `MAT1 ${int(f[3])}` });
        break;
      }
      case 'CTRIA3': case 'CTRIAR': fe.addElement('tri3', [int(f[2]), int(f[3]), int(f[4])], int(f[0]), int(f[1]) || int(f[0])); break;
      case 'CQUAD4': case 'CQUADR': fe.addElement('quad4', [int(f[2]), int(f[3]), int(f[4]), int(f[5])], int(f[0]), int(f[1]) || int(f[0])); break;
      case 'CTRIA6': solid({ name: c.name, fields: f.slice(0, 8) }, [3, 6], ['tri3', 'tri6']); break;
      case 'CQUAD8': solid({ name: c.name, fields: f.slice(0, 10) }, [4, 8], ['quad4', 'quad8']); break;
      case 'CTETRA': solid(c, [4, 10], ['tet4', 'tet10']); break;
      case 'CHEXA': solid(c, [8, 20], ['hex8', 'hex20'], [null, HEX20]); break;
      case 'CPENTA': solid(c, [6, 15], ['wedge6', 'wedge15']); break;
      case 'CPYRAM': solid(c, [5, 13], ['pyr5', 'pyr13']); break;
      case 'PBAR': case 'PBEAM': {
        const o = c.name === 'PBEAM' ? [2, 3, 4, 6] : [2, 3, 4, 5];
        fe.setProp(int(f[0]), { name: `${c.name} ${int(f[0])}`, material: `MAT1 ${int(f[1])}`, section: { A: real(f[o[0]]) * s * s, Iy: real(f[o[1]]) * s ** 4, Iz: real(f[o[2]]) * s ** 4, J: real(f[o[3]]) * s ** 4 } });
        break;
      }
      case 'PROD': fe.setProp(int(f[0]), { name: `PROD ${int(f[0])}`, material: `MAT1 ${int(f[1])}`, section: { A: real(f[2]) * s * s, Iy: 0, Iz: 0, J: real(f[3]) * s ** 4 } }); break;
      case 'PTUBE': {
        const D = real(f[2]) * s, t = real(f[3]) * s;
        // T = 0 or blank means a solid rod of diameter OD
        fe.setProp(int(f[0]), { name: `PTUBE ${int(f[0])}`, material: `MAT1 ${int(f[1])}`, section: t > 0 && 2 * t < D ? { type: 'tube', D, t } : { ...circleProps(D), shape: 'circle', D } });
        break;
      }
      case 'PBARL': case 'PBEAML': {
        const type = (f[3] || '').toUpperCase();
        const d = [];
        for (let i = 8; i < f.length && d.length < 8; i++) { if (f[i] === '' && d.length) break; if (f[i] !== '') d.push(real(f[i]) * s); }
        let section = null;
        if (type === 'TUBE' && d[0] > d[1] && d[1] > 0) section = { type: 'tube', D: 2 * d[0], t: d[0] - d[1] };
        else if (type === 'TUBE' && d[0] > 0) section = { ...circleProps(2 * d[0]), shape: 'circle', D: 2 * d[0] };
        else if (type === 'ROD' && d[0] > 0) section = { ...circleProps(2 * d[0]), shape: 'circle', D: 2 * d[0] };
        else if (type === 'BAR' && d[0] > 0 && d[1] > 0) section = { ...rectProps(d[0], d[1]), shape: 'rect', b: d[0], h: d[1] };
        else if (type === 'BOX' && d.length >= 4) section = { ...boxProps(d[0], d[1], d[3], d[2]), shape: 'box', b: d[0], h: d[1], tw: d[3], tf: d[2] };
        else if (type === 'I' && d.length >= 6 && Math.abs(d[1] - d[2]) < 1e-9 * d[1] && Math.abs(d[4] - d[5]) < 1e-9 * d[4]) section = { ...iSectionProps({ h: d[0], b: d[1], tw: d[3], tf: d[4] }), shape: 'I', h: d[0], b: d[1], tw: d[3], tf: d[4] };
        const rec = { name: `${c.name} ${int(f[0])} ${type}`, material: `MAT1 ${int(f[1])}` };
        if (section) rec.section = section; else warnings.push(`${c.name} ${int(f[0])}: cross-section type "${type}" is kept by name only (no properties computed).`);
        fe.setProp(int(f[0]), rec);
        break;
      }
      case 'PSHELL': fe.setProp(int(f[0]), { name: `PSHELL ${int(f[0])}`, material: `MAT1 ${int(f[1])}`, thickness: real(f[2]) * s }); break;
      case 'PSOLID': fe.setProp(int(f[0]), { name: `PSOLID ${int(f[0])}`, material: `MAT1 ${int(f[1])}` }); break;
      case 'PCOMP': case 'PCOMPG': fe.setProp(int(f[0]), { name: `${c.name} ${int(f[0])}` }); tally(skipped, c.name); break;
      case 'MAT1': {
        const E = real(f[1]), G = real(f[2]), nu = real(f[3]);
        materials[`MAT1 ${int(f[0])}`] = { E: E || (G && f[3] !== '' ? 2 * G * (1 + nu) : 0), G: G || (E && f[3] !== '' ? E / (2 * (1 + nu)) : 0), nu: f[3] !== '' ? nu : (E && G ? E / (2 * G) - 1 : 0), rho: real(f[4]), units: 'as written in the deck (not converted)' };
        break;
      }
      case 'ENDDATA': return false;
      default: tally(skipped, c.name || '(blank)');
    }
    return true;
  };

  let cur = null, more = true;
  const flush = () => { if (cur) { const c = cur; cur = null; return handle(c); } return true; };
  while (more && stack.length) {
    const lc = stack[stack.length - 1];
    const raw = lc.next();
    if (raw === null) { stack.pop(); continue; }
    if (raw[0] === '$') continue;
    const up = raw.trimStart().toUpperCase();
    if (!sawBulk && /^BEGIN\s+BULK/.test(up)) { sawBulk = true; cur = null; skipped = {}; cards = 0; continue; }
    if (up.startsWith('INCLUDE')) {
      more = flush();
      const m = /INCLUDE\s+['"]?([^'"]+)['"]?/i.exec(raw);
      const f = m && findFile(input.files, m[1].trim());
      if (f) stack.push(new LineCursor(f)); else warnings.push(`INCLUDE ${m ? m[1].trim() : ''}: the included file was not supplied and is missing from the model.`);
      continue;
    }
    const c = splitCard(raw);
    if (!c) continue;
    if (c.cont) { if (cur) for (const v of c.fields) cur.fields.push(v); continue; }
    more = flush();
    cur = c;
  }
  if (more) flush();

  // GRIDs given in local coordinate systems → basic system
  if (gridCp.length) {
    const frames = new Map();
    let unresolved = 0;
    for (const [cid, c] of cords) {
      if (c.rid !== 0) continue;
      const z = V3.norm(V3.sub(c.b, c.a));
      const y = V3.norm(V3.cross(z, V3.sub(c.c, c.a)));
      frames.set(cid, { o: c.a, x: V3.cross(y, z), y, z, type: c.type });
    }
    const xyz = fe.coords.a, D = Math.PI / 180;
    for (let i = 0; i < gridCp.length; i += 2) {
      const fr = frames.get(gridCp[i + 1]);
      if (!fr) { unresolved++; continue; }
      const k = 3 * fe.nodeIndex.get(gridCp[i]);
      let u = xyz[k], v = xyz[k + 1], w = xyz[k + 2];
      if (fr.type === 'C') { const r = u, th = v * D; u = r * Math.cos(th); v = r * Math.sin(th); } else if (fr.type === 'S') { const r = u, th = v * D, ph = w * D; u = r * Math.sin(th) * Math.cos(ph); v = r * Math.sin(th) * Math.sin(ph); w = r * Math.cos(th); }
      xyz[k] = fr.o[0] + u * fr.x[0] + v * fr.y[0] + w * fr.z[0];
      xyz[k + 1] = fr.o[1] + u * fr.x[1] + v * fr.y[1] + w * fr.z[1];
      xyz[k + 2] = fr.o[2] + u * fr.x[2] + v * fr.y[2] + w * fr.z[2];
    }
    if (unresolved) warnings.push(`${unresolved} GRID point(s) refer to a coordinate system that is not a CORD2R/CORD2C/CORD2S defined in the basic system; their coordinates were left untransformed.`);
  }
  if (!cards && !fe.nodeCount) throw new Error(`${input.name || 'BDF'}: no bulk-data cards were found (expected GRID and element cards such as CBAR, CQUAD4, CTETRA, optionally after BEGIN BULK).`);
  const nSkip = Object.values(skipped).reduce((a, b) => a + b, 0);
  if (nSkip) warnings.push(`${nSkip} card(s) of ${Object.keys(skipped).length} type(s) are not geometry this reader keeps (${Object.keys(skipped).slice(0, 10).join(', ')}${Object.keys(skipped).length > 10 ? ', …' : ''}).`);
  return feResult(fe, 'nastran', opts, {
    name: input.name || 'BDF', expected: 'expected GRID cards', warnings, skipped: { cards: skipped },
    extra: { materials, coordinateSystems: cords.size, properties: fe.props.size },
  });
}
