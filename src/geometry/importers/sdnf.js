// Steel Detailing Neutral File (SDNF) version 3.0.
//
// Read: Packet 00 (title), Packet 10 (linear members) and Packet 20 (plates).
//   members → frame: end coordinates, section size resolved to properties where the designation is known
//             (pipes / tubulars, box sections, IPE / HE A / HE B, the tabulated AISC W shapes — see
//             sections.js), grade as material, member type as group, plus cardinal point, rotation,
//             orientation vector, end offsets and releases carried on each member.
//   plates  → surface: each plate outline triangulated, one group per plate with its thickness.
// Packets 22 (holes), 30 (member loads), 40 (connection details), 50 (grids) and 60 (curved members) are
// counted in `meta.skipped`.
//
// The file is read as a token stream (numbers and quoted strings), so the reader does not depend on how a
// particular writer wraps its records.

import { lineCursor, baseName, Grow, Welder, finish, surfaceModel, tally, triangulate3D, unitAdvice } from './common.js';
import { resolveUnits, suggestUnits, normaliseUnit } from '../units.js';
import { resolveSection } from './sections.js';

function tokenize(line) {
  const out = [];
  const re = /"([^"]*)"|'([^']*)'|(\S+)/g;
  for (let m = re.exec(line); m; m = re.exec(line)) {
    if (m[1] !== undefined) out.push({ s: m[1] }); else if (m[2] !== undefined) out.push({ s: m[2] }); else { const v = Number(m[3].replace(/[dD](?=[-+]?\d)/, 'e')); out.push(v === v ? { n: v } : { s: m[3] }); }
  }
  return out;
}

/** SDNF 3.0. */
export function parseSdnf(input, opts = {}) {
  const cur = lineCursor(input);
  const name = baseName(input.name);
  const fname = input.name || 'SDNF';
  const packets = new Map();
  let tok = null;
  for (let line = cur.next(); line !== null; line = cur.next()) {
    const t = line.trim();
    if (!t || t[0] === '#') continue;
    const m = /^Packet\s+0*(\d+)\b/i.exec(t);
    if (m) { tok = []; const k = Number(m[1]); if (!packets.has(k)) packets.set(k, tok); else tok = packets.get(k); continue; }
    if (tok) for (const x of tokenize(t)) tok.push(x);
  }
  if (!packets.size) throw new Error(`${fname}: not an SDNF file (expected lines "Packet 00", "Packet 10" … introducing each data packet).`);
  const warnings = [];
  const skipped = {}, hist = {};
  const title = (packets.get(0) || []).filter((x) => x.s !== undefined).map((x) => x.s.trim()).filter(Boolean);

  const unitOf = (toks) => { const u = toks.find((x) => x.s !== undefined); return u ? u.s : null; };
  const p10 = packets.get(10) || [], p20 = packets.get(20) || [];
  const unitName = unitOf(p10) || unitOf(p20);
  const units = resolveUnits(opts, unitName && normaliseUnit(unitName) ? unitName : null, 'm');
  if (unitName && !normaliseUnit(unitName)) warnings.push(`Linear unit "${unitName}" was not recognised; ${units.unit} assumed.`);
  const s = units.scale;
  // section designations in an inch/foot file give dimensions in inches; metric files in millimetres
  const dimUnit = ['in', 'ft', 'us-ft'].includes(units.unit) ? 'in' : 'mm';
  const models = [];

  if (p10.length) {
    let i = 0;
    const isN = () => i < p10.length && p10[i].n !== undefined;
    const isS = () => i < p10.length && p10[i].s !== undefined;
    const nums = (max) => { const o = []; while (isN() && o.length < max) o.push(p10[i++].n); return o; };
    const strs = (max) => { const o = []; while (isS() && o.length < max) o.push(p10[i++].s.trim()); return o; };
    strs(1);
    const declared = nums(1)[0];
    const weld = new Welder(1e-6, 256);
    const members = [];
    const unresolved = new Set();
    let bad = 0, zero = 0;
    while (i < p10.length) {
      const head = nums(8);
      const s1 = strs(2);
      nums(4); // revision number
      const s2 = strs(2);
      const v = nums(32);
      if (!s1.length && !s2.length && !v.length) { i++; bad++; continue; }
      // some writers put all four strings together, others split them around the revision number
      const strsAll = s1.concat(s2);
      const [type, mark, size, grade] = strsAll.length >= 4 ? strsAll : strsAll.length === 3 ? [strsAll[0], '', strsAll[1], strsAll[2]] : [s1[0] || '', s1[1] || '', s2[0] || '', s2[1] || ''];
      if (v.length < 12) { bad++; continue; }
      const [rot, mirX, mirY, ox, oy, oz, x1, y1, z1, x2, y2, z2] = v;
      const a = weld.add(x1 * s, y1 * s, z1 * s), b = weld.add(x2 * s, y2 * s, z2 * s);
      if (a === b) { zero++; continue; }
      const m = { a, b, id: head[0], group: type || 'member' };
      const res = size ? resolveSection(size, { dimUnit }) : { matched: false };
      if (res.matched) m.section = { ...res.section, name: size }; else if (size) { m.sectionName = size; unresolved.add(size); }
      if (grade) m.material = grade;
      if (mark) m.mark = mark;
      m.cardinalPoint = head[1] ?? 0;
      m.rotationDeg = rot;
      if (mirX || mirY) m.mirror = { x: !!mirX, y: !!mirY };
      m.orient = [ox, oy, oz];
      if (v.length >= 20 && v.slice(14, 20).some((q) => q !== 0)) m.offsets = { a: v.slice(14, 17).map((q) => q * s), b: v.slice(17, 20).map((q) => q * s) };
      if (v.length >= 32 && v.slice(20, 32).some((q) => q !== 0)) m.releases = { a: v.slice(20, 26), b: v.slice(26, 32) };
      members.push(m);
    }
    if (members.length) { models.push({ kind: 'frame', name, nodes: weld.finish(), members }); tally(hist, 'line2', members.length); }
    if (declared !== undefined && declared !== members.length + zero + bad) warnings.push(`Packet 10 declares ${declared} member(s); ${members.length} were read.`);
    if (bad) { skipped.incompleteMembers = bad; warnings.push(`${bad} member record(s) in Packet 10 were incomplete and skipped.`); }
    if (zero) { skipped.zeroLengthMembers = zero; warnings.push(`${zero} zero-length member(s) were skipped.`); }
    if (unresolved.size) warnings.push(`Section size(s) not resolved to properties (kept by name only): ${[...unresolved].slice(0, 12).join(', ')}${unresolved.size > 12 ? ' …' : ''}.`);
  }

  if (p20.length) {
    let i = 0;
    const isN = () => i < p20.length && p20[i].n !== undefined;
    const isS = () => i < p20.length && p20[i].s !== undefined;
    const nums = (max) => { const o = []; while (isN() && o.length < max) o.push(p20[i++].n); return o; };
    const strs = (max) => { const o = []; while (isS() && o.length < max) o.push(p20[i++].s.trim()); return o; };
    const us = strs(2);
    nums(1);
    const tu = us[1] && normaliseUnit(us[1]) ? resolveUnits({}, us[1]).scale : s;
    const coords = new Grow(Float64Array, 3072), idx = new Grow(Uint32Array, 3072);
    const groups = [];
    let bad = 0;
    while (i < p20.length) {
      const head = nums(8);
      const st = strs(3);
      const v = nums(2);
      if (v.length < 2 || !(v[1] >= 3)) { bad++; if (!head.length && !st.length && !v.length) i++; continue; }
      const n = v[1];
      const c = nums(3 * n);
      if (c.length < 3 * n) { bad++; continue; }
      const flat = Float64Array.from(c, (q) => q * s);
      const t = triangulate3D(flat);
      if (!t.length) { bad++; continue; }
      const base = coords.n / 3, start = idx.n;
      for (const q of flat) coords.push(q);
      for (const k of t) idx.push(base + k);
      groups.push({ name: st[1] || `plate ${head[0] ?? groups.length + 1}`, start, count: idx.n - start, thickness: v[0] * tu, material: st[2] || undefined, type: st[0] || undefined, id: head[0] });
      tally(hist, 'plate');
    }
    if (groups.length) models.push(surfaceModel(models.length ? `${name} (plates)` : name, coords.finish(), idx.finish(), groups));
    if (bad) { skipped.incompletePlates = bad; warnings.push(`${bad} plate record(s) in Packet 20 were incomplete and skipped.`); }
  }
  const NAMES = { 22: 'holes', 30: 'member loads', 40: 'connection details', 50: 'grids', 60: 'curved members' };
  for (const k of packets.keys()) if (![0, 10, 20].includes(k) && packets.get(k).length) { skipped[`Packet ${k}${NAMES[k] ? ` (${NAMES[k]})` : ''}`] = 1; warnings.push(`Packet ${k}${NAMES[k] ? ` (${NAMES[k]})` : ''} is not read.`); }
  if (!models.length) throw new Error(`${fname}: no members (Packet 10) or plates (Packet 20) were found in the SDNF file.`);
  const r = finish(models, warnings, { format: 'sdnf', units, crs: null, elementTypes: hist, skipped, title, packets: [...packets.keys()].sort((a, b) => a - b) });
  unitAdvice(units, r.meta.bbox, suggestUnits, r.warnings);
  return r;
}
