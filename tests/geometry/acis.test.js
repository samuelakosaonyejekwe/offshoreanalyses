import test from 'node:test';
import assert from 'node:assert/strict';
import { parseAcis, readSat } from '../../src/geometry/importers/acis.js';
import { importFile, formatFor, GeometryFormatError } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { surfaceArea, signedVolume } from '../../src/mesh/measure.js';
import { enc, close, assertContract } from './fixtures/helpers.js';

// ---------------------------------------------------------------------------------------------------------
// Fixtures. They are written here, record by record, from the published structure of the SAT format, so that
// every expected area and volume is known in closed form. Files written by ACIS-based CAD systems are read
// in acis.real.test.js. `Sat` only numbers the records and adds the per-version record
// prefix; every field of every record is spelled out by the helpers below.
//   version < 700 : "-n type $attrib fields #"            (sequence numbers, counted strings "len text")
//   version ≥ 700 : "type $attrib -1 $history fields #"   (strings "@len text", "End-of-ACIS-data")
// ---------------------------------------------------------------------------------------------------------

const P = (i) => `$${i}`;
const NIL = '$-1';

class Sat {
  constructor(version, { mmPerUnit = 1000 } = {}) { this.version = version; this.mm = mmPerUnit; this.recs = []; }
  id() { this.recs.push(null); return this.recs.length - 1; }
  set(i, type, fields) { this.recs[i] = [type, fields]; return i; }
  add(type, fields) { return this.set(this.id(), type, fields); }
  text() {
    const modern = this.version >= 700;
    const str = (s) => `${modern ? '@' : ''}${s.length} ${s}`;
    let out = `${this.version} ${this.recs.length} 1 0\n${str('fixture')} ${str(modern ? 'ACIS 7.0 NT' : 'ACIS 4.0 NT')} ${str('Thu Jan 01 00:00:00 2004')}\n${this.mm} 9.9999999999999995e-007 1e-010\n`;
    this.recs.forEach(([type, fields], i) => { out += modern ? `${type} $-1 -1 $-1 ${fields.join(' ')} #\n` : `-${i} ${type} $-1 ${fields.join(' ')} #\n`; });
    return modern ? `${out}End-of-ACIS-data\n` : out;
  }
}

/** One body with one lump and one shell; faces are added with their surface record, sense and loops. */
function solid(version, opts = {}) {
  const s = new Sat(version, opts);
  const body = s.id(), lump = s.id(), shell = s.id();
  const faces = [], vinfo = new Map(), einfo = new Map(), uses = new Map(), memo = new Map();
  const b = {
    s,
    vertex(x, y, z) { const pt = s.add('point', [x, y, z]); const v = s.id(); vinfo.set(v, { pt, xyz: [x, y, z], edge: -1 }); return v; },
    edge(v0, v1, curve, t0, t1) {
      const e = s.id();
      einfo.set(e, { v0, v1, curve, t0, t1 });
      for (const v of [v0, v1]) if (vinfo.get(v).edge < 0) vinfo.get(v).edge = e;
      return e;
    },
    line(v0, v1) {
      const a = vinfo.get(v0).xyz, c = vinfo.get(v1).xyz, d = [c[0] - a[0], c[1] - a[1], c[2] - a[2]], l = Math.hypot(...d);
      return b.edge(v0, v1, s.add('straight-curve', [...a, ...d.map((x) => x / l), 'I', 'I']), 0, l);
    },
    /** Edge between two vertices of a polyhedron, created on first use: [edge, runs against the edge]. */
    between(v0, v1) {
      if (memo.has(`${v1},${v0}`)) return [memo.get(`${v1},${v0}`), true];
      if (!memo.has(`${v0},${v1}`)) memo.set(`${v0},${v1}`, b.line(v0, v1));
      return [memo.get(`${v0},${v1}`), false];
    },
    /** Circular arc from v0 to v1, counter-clockwise about `normal`; `major` points from the centre to parameter 0. */
    arc(v0, v1, centre, normal, major, t0, t1) { return b.edge(v0, v1, s.add('ellipse-curve', [...centre, ...normal, ...major, 1, 'I', 'I']), t0, t1); },
    /** Closed circle (or ellipse) through vertex v; runs counter-clockwise about `normal`. */
    circle(v, centre, normal, major, ratio = 1) { return b.edge(v, v, s.add('ellipse-curve', [...centre, ...normal, ...major, ratio, 'I', 'I']), 0, 2 * Math.PI); },
    face(surface, sense, loops) {
      const f = s.id();
      const loopIds = loops.map(() => s.id());
      loops.forEach((loop, li) => {
        const ce = loop.map(() => s.id());
        loop.forEach(([e, rev], k) => {
          if (!uses.has(e)) uses.set(e, []);
          uses.get(e).push(ce[k]);
          // next, previous, partner (filled in below), edge, sense, loop, pcurve
          s.set(ce[k], 'coedge', [P(ce[(k + 1) % ce.length]), P(ce[(k + ce.length - 1) % ce.length]), null, P(e), rev ? 'reversed' : 'forward', P(loopIds[li]), NIL]);
        });
        s.set(loopIds[li], 'loop', [li + 1 < loopIds.length ? P(loopIds[li + 1]) : NIL, P(ce[0]), P(f)]);
      });
      faces.push({ f, surface: s.add(surface[0], surface.slice(1)), sense, loop: loopIds.length ? loopIds[0] : -1 });
      return f;
    },
    plane: (root, normal, u = null) => ['plane-surface', ...root, ...normal, ...(u || (Math.abs(normal[2]) > 0.5 ? [1, 0, 0] : [0, 0, 1])), 'forward_v', 'I', 'I', 'I', 'I'],
    cone: (centre, axis, major, sin = 0, cos = 1, ratio = 1) => ['cone-surface', ...centre, ...axis, ...major, ratio, 'I', 'I', sin, cos, Math.hypot(...major), 'forward', 'I', 'I', 'I', 'I'],
    sphere: (centre, radius) => ['sphere-surface', ...centre, radius, 1, 0, 0, 0, 0, 1, 'forward_v', 'I', 'I', 'I', 'I'],
    torus: (centre, normal, major, minor, uv) => ['torus-surface', ...centre, ...normal, major, minor, ...uv, 'forward_v', 'I', 'I', 'I', 'I'],
    /** Planar face from vertex cycles (the first is the outer boundary), counter-clockwise about `normal`. */
    polygon(normal, cycles, extraLoops = []) {
      const loops = cycles.map((c) => c.map((v, k) => b.between(v, c[(k + 1) % c.length])));
      return b.face(b.plane(vinfo.get(cycles[0][0]).xyz, normal), 'forward', [...loops, ...extraLoops]);
    },
    text(transform = null) {
      const tr = transform ? s.add('transform', [...transform, 'no_rotate', 'no_reflect', 'no_shear']) : -1;
      s.set(body, 'body', [P(lump), NIL, tr < 0 ? NIL : P(tr)]);
      s.set(lump, 'lump', [NIL, P(shell), P(body)]);
      s.set(shell, 'shell', [NIL, NIL, faces.length ? P(faces[0].f) : NIL, NIL, P(lump)]);
      faces.forEach((q, i) => s.set(q.f, 'face', [i + 1 < faces.length ? P(faces[i + 1].f) : NIL, q.loop < 0 ? NIL : P(q.loop), P(shell), NIL, P(q.surface), q.sense, 'single']));
      for (const [e, q] of einfo) s.set(e, 'edge', [P(q.v0), q.t0, P(q.v1), q.t1, uses.has(e) ? P(uses.get(e)[0]) : NIL, P(q.curve), 'forward', ...(version >= 700 ? ['@7 unknown'] : [])]);
      for (const [v, q] of vinfo) s.set(v, 'vertex', [P(q.edge), P(q.pt)]);
      // the coedges of an edge point to each other in a ring
      for (const list of uses.values()) list.forEach((c, k) => { s.recs[c][1][2] = list.length > 1 ? P(list[(k + 1) % list.length]) : NIL; });
      return s.text();
    },
  };
  return b;
}

/** Rectangular block: 8 vertices, 12 straight edges, 6 planar faces with outward normals. */
function addBox(b, lo, hi, { top = [], bottom = [] } = {}) {
  const v = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]].map((c) => b.vertex(...c.map((t, k) => (t ? hi[k] : lo[k]))));
  b.polygon([0, 0, -1], [[v[0], v[3], v[2], v[1]]], bottom);
  b.polygon([0, 0, 1], [[v[4], v[5], v[6], v[7]]], top);
  b.polygon([0, -1, 0], [[v[0], v[1], v[5], v[4]]]);
  b.polygon([1, 0, 0], [[v[1], v[2], v[6], v[5]]]);
  b.polygon([0, 1, 0], [[v[2], v[3], v[7], v[6]]]);
  b.polygon([-1, 0, 0], [[v[3], v[0], v[4], v[7]]]);
}

const boxSat = (version, opts, transform) => { const b = solid(version, opts); addBox(b, [0, 0, 0], [2, 3, 4]); return b.text(transform); };

/** Cylinder of radius r between z0 and z1 about the z axis: side and two flat ends. */
function cylinderSat(version, r, z0, z1) {
  const b = solid(version);
  const vb = b.vertex(r, 0, z0), vt = b.vertex(r, 0, z1);
  const eb = b.circle(vb, [0, 0, z0], [0, 0, 1], [r, 0, 0]), et = b.circle(vt, [0, 0, z1], [0, 0, 1], [r, 0, 0]);
  b.face(b.cone([0, 0, z0], [0, 0, 1], [r, 0, 0]), 'forward', [[[eb, false]], [[et, true]]]);
  b.face(b.plane([0, 0, z1], [0, 0, 1]), 'forward', [[[et, false]]]);
  b.face(b.plane([0, 0, z0], [0, 0, -1]), 'forward', [[[eb, true]]]);
  return b.text();
}

/** Tubular member: outer and inner cylinders (the bore is a reversed face) and two annular ends. */
function tubeSat(version, ro, ri, len, transform) {
  const b = solid(version);
  const ring = (r, z) => { const v = b.vertex(r, 0, z); return b.circle(v, [0, 0, z], [0, 0, 1], [r, 0, 0]); };
  const ob = ring(ro, 0), ot = ring(ro, len), ib = ring(ri, 0), it = ring(ri, len);
  b.face(b.cone([0, 0, 0], [0, 0, 1], [ro, 0, 0]), 'forward', [[[ob, false]], [[ot, true]]]);
  b.face(b.cone([0, 0, 0], [0, 0, 1], [ri, 0, 0]), 'reversed', [[[ib, true]], [[it, false]]]);
  b.face(b.plane([0, 0, len], [0, 0, 1]), 'forward', [[[ot, false]], [[it, true]]]);
  b.face(b.plane([0, 0, 0], [0, 0, -1]), 'forward', [[[ob, true]], [[ib, false]]]);
  return b.text(transform);
}

const surf = (r) => r.models.find((m) => m.kind === 'surface');
const vol = (m) => signedVolume(m.positions, m.indices);
const area = (m) => surfaceArea(m.positions, m.indices);
const tight = (m) => { const q = inspect(m); assert.ok(q.watertight, `watertight: ${(q.reasons || q.issues.map((i) => i.message)).join('; ')}`); return q; };

// ---------------------------------------------------------------------------------------------------------

test('SAT header and records: both record layouts, counted strings, pointers', () => {
  for (const version of [400, 700]) {
    const { header, records, numbered } = readSat(boxSat(version));
    assert.equal(header.version, version);
    assert.equal(header.product, 'fixture');
    assert.equal(header.acis, version >= 700 ? 'ACIS 7.0 NT' : 'ACIS 4.0 NT');
    assert.equal(header.date, 'Thu Jan 01 00:00:00 2004');
    assert.equal(header.mmPerUnit, 1000);
    assert.equal(header.resabs, 9.9999999999999995e-7);
    assert.equal(header.records, records.length);
    assert.equal(numbered, version < 700);
    const count = (k) => records.filter((r) => r.kind === k).length;
    assert.deepEqual(['body', 'lump', 'shell', 'face', 'loop', 'coedge', 'edge', 'vertex', 'point'].map(count), [1, 1, 1, 6, 6, 24, 12, 8, 8]);
    assert.equal(records.filter((r) => r.type === 'plane-surface').length, 6);
    assert.equal(records.filter((r) => r.type === 'straight-curve').length, 12);
    assert.equal(records[0].type, 'body');
    assert.deepEqual(records[0].tok.slice(0, version >= 700 ? 4 : 2).map((t) => t.t), version >= 700 ? ['p', 'n', 'p', 'p'] : ['p', 'p']);
  }
  // the 7.0 edge record ends with a counted string
  const e = readSat(boxSat(700)).records.find((r) => r.kind === 'edge');
  assert.deepEqual(e.tok[e.tok.length - 1], { t: 's', v: 'unknown' });
});

test('SAT (a): box of six planar faces → closed surface with exact volume and area, in both record layouts', async () => {
  for (const version of [400, 700]) {
    const r = await importFile({ name: 'box.sat', bytes: enc(boxSat(version)) });
    assert.equal(r.format, 'acis');
    assert.equal(r.meta.detected.id, 'acis');
    assertContract(r, ['surface']);
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    assert.equal(m.name, 'box');
    assert.equal(m.positions.length / 3, 8);
    assert.equal(m.indices.length / 3, 12);
    assert.deepEqual(m.groups.map((g) => [g.name, g.count]), [1, 2, 3, 4, 5, 6].map((k) => [`face ${k}`, 6]));
    tight(m);
    close(vol(m), 24, 1e-12, 'volume');
    close(area(m), 52, 1e-12, 'area');
    assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [2, 3, 4] });
    assert.equal(r.meta.units.unit, 'm');
    assert.equal(r.meta.units.basis, 'file');
    assert.equal(r.meta.acisVersion, version);
    assert.equal(r.meta.closedBodies, 1);
    assert.deepEqual(r.meta.elementTypes, { plane: 6 });
  }
  assert.equal(formatFor('box.sat', enc(boxSat(700))).status, 'native');
  // detected from the content when the extension says nothing
  assert.equal((await importFile({ name: 'box.txt', bytes: enc(boxSat(400)) })).format, 'acis');
});

test('SAT: header unit and body transform are applied', () => {
  // millimetres per model unit = 1: a 2 × 3 × 4 mm block
  const mm = parseAcis({ name: 'b.sat', text: boxSat(700, { mmPerUnit: 1 }) });
  assert.equal(mm.meta.units.unit, 'mm');
  close(vol(surf(mm)) / 24e-9, 1, 1e-6);
  const inch = parseAcis({ name: 'b.sat', text: boxSat(400, { mmPerUnit: 25.4 }) });
  assert.equal(inch.meta.units.unit, 'in');
  close(inch.meta.bbox.max[2], 4 * 0.0254, 1e-7);
  // a unit the table does not name is kept as a scale
  const odd = parseAcis({ name: 'b.sat', text: boxSat(700, { mmPerUnit: 500 }) });
  assert.deepEqual(odd.meta.units, { unit: 'custom', scale: 0.5, basis: 'file' });
  close(odd.meta.bbox.max[1], 1.5, 1e-12);
  // the caller's unit wins over the header
  close(parseAcis({ name: 'b.sat', text: boxSat(700) }, { units: 'mm' }).meta.bbox.max[0], 0.002, 1e-9);
  // rotation by +90° about z (rows are the images of the axes), scale 2, then translation (10, 20, 30)
  const t = parseAcis({ name: 'b.sat', text: boxSat(700, {}, [0, 1, 0, -1, 0, 0, 0, 0, 1, 10, 20, 30, 2]) });
  assert.deepEqual(t.meta.bbox, { min: [4, 20, 30], max: [10, 24, 38] });
  close(vol(surf(t)), 8 * 24, 1e-9);
  tight(surf(t));
  // a mirroring transform keeps the normals outward
  const mir = parseAcis({ name: 'b.sat', text: boxSat(400, {}, [-1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]) });
  assert.deepEqual(mir.meta.bbox, { min: [-2, 0, 0], max: [0, 3, 4] });
  close(vol(surf(mir)), 24, 1e-12);
});

test('SAT (b): capped cylinder → closed; area and volume within the tessellation error and converging with deflection', () => {
  const R = 1, H = 3, V = Math.PI * R * R * H, A = 2 * Math.PI * R * H + 2 * Math.PI * R * R;
  let prevV = Infinity, prevA = Infinity;
  const segs = [];
  for (const deflection of [undefined, 0.004, 0.0005]) {
    const r = parseAcis({ name: 'cyl.sat', text: cylinderSat(700, R, 0, H) }, deflection ? { deflection } : {});
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    tight(m);
    assert.deepEqual(r.meta.elementTypes, { cylinder: 1, plane: 2 });
    // facets are inscribed: all vertices on the cylinder or inside the flat ends, none outside
    let n = 0;
    for (let i = 0; i < m.positions.length; i += 3) { const rr = Math.hypot(m.positions[i], m.positions[i + 1]); assert.ok(rr <= R + 1e-6); if (Math.abs(m.positions[i + 2]) < 1e-9 && Math.abs(rr - R) < 1e-6) n++; }
    segs.push(n);
    // an n-sided inscribed prism: V = (n/2) sin(2π/n) R² H
    close(vol(m), (n / 2) * Math.sin(2 * Math.PI / n) * R * R * H, 2e-6, 'prism volume');
    const eV = (V - vol(m)) / V, eA = (A - area(m)) / A;
    assert.ok(eV > 0 && eA > 0, 'inscribed: below the exact values');
    // sagitta of one facet relative to the radius bounds both errors
    const sag = 1 - Math.cos(Math.PI / n);
    assert.ok(eV < 2.1 * sag && eA < 2.1 * sag, `errors ${eV}, ${eA} against sagitta ${sag}`);
    if (deflection) assert.ok(sag * R <= deflection * (1 + 1e-9), `chord error ${sag * R} within the requested ${deflection}`);
    assert.ok(eV < prevV && eA < prevA, 'converging');
    prevV = eV; prevA = eA;
  }
  // default: 1 % of the diagonal would allow 33° per facet, so the 15° limit governs
  assert.equal(segs[0], 24);
  assert.ok(segs[1] > segs[0] && segs[2] > segs[1]);
  assert.ok(prevV < 1e-3);
  // the angular limit alone
  const fine = parseAcis({ name: 'cyl.sat', text: cylinderSat(400, R, 0, H) }, { angularStep: 5 });
  close(vol(surf(fine)), 36 * Math.sin(2 * Math.PI / 72) * R * R * H, 2e-6);
  assert.throws(() => parseAcis({ name: 'c.sat', text: cylinderSat(700, 1, 0, 1) }, { deflection: 0 }), /deflection: must be a positive length/);
});

test('SAT (c): plate with a circular hole — planar faces with two loops and the bore as a reversed cylinder', () => {
  for (const version of [400, 700]) {
    const b = solid(version);
    const r0 = 0.5, t = 0.25;
    const vb = b.vertex(2 + r0, 1.5, 0), vt = b.vertex(2 + r0, 1.5, t);
    const eb = b.circle(vb, [2, 1.5, 0], [0, 0, 1], [r0, 0, 0]), et = b.circle(vt, [2, 1.5, t], [0, 0, 1], [r0, 0, 0]);
    addBox(b, [0, 0, 0], [4, 3, t], { top: [[[et, true]]], bottom: [[[eb, false]]] });
    b.face(b.cone([2, 1.5, 0], [0, 0, 1], [r0, 0, 0]), 'reversed', [[[eb, true]], [[et, false]]]);
    const r = parseAcis({ name: 'plate.sat', text: b.text() });
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    const q = tight(m);
    assert.equal(m.groups.length, 7);
    assert.deepEqual(r.meta.elementTypes, { plane: 6, cylinder: 1 });
    // the hole is a 24-gon at the default tolerance (15° per facet)
    const hole = 12 * Math.sin(2 * Math.PI / 24) * r0 * r0;
    close(vol(m), (12 - hole) * t, 1e-6, 'volume');
    close(area(m), 2 * (12 - hole) + 2 * (4 + 3) * t + 24 * 2 * r0 * Math.sin(Math.PI / 24) * t, 1e-6, 'area');
    // one through-hole: genus 1
    assert.equal(q.stats.genus, 1);
  }
});

// Full sphere, written out in full: one face without loops on a sphere-surface (centre, radius, direction of
// u = 0 on the equator, pole direction).
const SPHERE_400 = `400 0 1 0
7 fixture 11 ACIS 4.0 NT 24 Thu Jan 01 00:00:00 2004
1000 9.9999999999999995e-007 1e-010
-0 body $-1 $1 $-1 $-1 #
-1 lump $-1 $-1 $2 $0 #
-2 shell $-1 $-1 $-1 $3 $-1 $1 #
-3 face $-1 $-1 $-1 $2 $-1 $4 forward single #
-4 sphere-surface $-1 1 -2 0.5 2 1 0 0 0 0 1 forward_v I I I I #
`;
const SPHERE_700 = `700 0 1 0
@7 fixture @11 ACIS 7.0 NT @24 Thu Jan 01 00:00:00 2004
1000 9.9999999999999995e-007 1e-010
body $-1 -1 $-1 $1 $-1 $-1 #
lump $-1 -1 $-1 $-1 $2 $0 #
shell $-1 -1 $-1 $-1 $-1 $3 $-1 $1 #
face $-1 -1 $-1 $-1 $-1 $2 $-1 $4 forward single #
sphere-surface $-1 -1 $-1 1 -2 0.5 2 1 0 0 0 0 1 forward_v I I I I #
End-of-ACIS-data
`;

test('SAT (d): sphere → closed surface on the sphere, volume and area converging', () => {
  const R = 2, V = (4 / 3) * Math.PI * R ** 3, A = 4 * Math.PI * R * R;
  for (const text of [SPHERE_400, SPHERE_700]) {
    let prev = Infinity;
    for (const [opts, lim] of [[{}, 0.03], [{ angularStep: 6 }, 0.006]]) {
      const r = parseAcis({ name: 'ball.sat', text }, opts);
      assert.deepEqual(r.warnings, []);
      const m = surf(r);
      tight(m);
      for (let i = 0; i < m.positions.length; i += 3) close(Math.hypot(m.positions[i] - 1, m.positions[i + 1] + 2, m.positions[i + 2] - 0.5), R, 2e-6, 'on the sphere');
      const eV = (V - vol(m)) / V, eA = (A - area(m)) / A;
      assert.ok(eV > 0 && eV < lim && eA > 0 && eA < lim, `errors ${eV}, ${eA}`);
      assert.ok(eV < prev);
      prev = eV;
      close(r.meta.bbox.max[2], 2.5, 1e-6);
      close(r.meta.bbox.min[1], -4, 0.02);
    }
  }
  // hemisphere: one circular edge shared by the spherical face and a flat face
  {
    const b = solid(700);
    const v = b.vertex(R, 0, 0);
    const e = b.circle(v, [0, 0, 0], [0, 0, 1], [R, 0, 0]);
    b.face(b.sphere([0, 0, 0], R), 'forward', [[[e, false]]]);
    b.face(b.plane([0, 0, 0], [0, 0, -1]), 'forward', [[[e, true]]]);
    const r = parseAcis({ name: 'dome.sat', text: b.text() }, { angularStep: 6 });
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    tight(m);
    close(vol(m), V / 2, 0.006 * V, 'hemisphere volume');
    assert.ok(r.meta.bbox.min[2] > -1e-9 && Math.abs(r.meta.bbox.max[2] - R) < 1e-6, 'upper half');
  }
});

test('SAT: sphere bounded through its pole, spherical cavity, end rings that start at different angles', () => {
  // one eighth of a ball: the spherical face is bounded by three arcs, two of which meet at the pole
  const R = 2, H = Math.PI / 2;
  for (const version of [400, 700]) {
    const b = solid(version);
    const O = b.vertex(0, 0, 0), A = b.vertex(R, 0, 0), B = b.vertex(0, R, 0), N = b.vertex(0, 0, R);
    const ab = b.arc(A, B, [0, 0, 0], [0, 0, 1], [R, 0, 0], 0, H), bn = b.arc(B, N, [0, 0, 0], [1, 0, 0], [0, R, 0], 0, H), na = b.arc(N, A, [0, 0, 0], [0, 1, 0], [0, 0, R], 0, H);
    const oa = b.line(O, A), ob = b.line(O, B), on = b.line(O, N);
    b.face(b.sphere([0, 0, 0], R), 'forward', [[[ab, false], [bn, false], [na, false]]]);
    b.face(b.plane([0, 0, 0], [0, 0, -1]), 'forward', [[[ob, false], [ab, true], [oa, true]]]);
    b.face(b.plane([0, 0, 0], [-1, 0, 0]), 'forward', [[[on, false], [bn, true], [ob, true]]]);
    b.face(b.plane([0, 0, 0], [0, -1, 0]), 'forward', [[[oa, false], [na, true], [on, true]]]);
    const r = parseAcis({ name: 'octant.sat', text: b.text() }, { angularStep: 6 });
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    tight(m);
    const V = Math.PI * R ** 3 / 6;
    assert.ok(vol(m) < V && vol(m) > 0.994 * V, `octant volume ${vol(m)} / ${V}`);
    assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [2, 2, 2] });
  }
  // a cavity is the same sphere seen from inside: a reversed face, or a negative radius
  for (const [radius, sense] of [[1, 'reversed'], [-1, 'forward']]) {
    const b = solid(700);
    addBox(b, [-2, -2, -2], [2, 2, 2]);
    b.face(b.sphere([0, 0, 0], radius), sense, []);
    const r = parseAcis({ name: 'cavity.sat', text: b.text() }, { angularStep: 6 });
    assert.deepEqual(r.warnings, []);
    const q = tight(surf(r));
    assert.equal(q.stats.components, 2);
    const V = 64 - (4 / 3) * Math.PI;
    assert.ok(vol(surf(r)) > V && vol(surf(r)) < V + 0.006 * 4.2, `volume ${vol(surf(r))} / ${V}`);
  }
  // the two rings of a cylinder need not start at the same angle: the seam then runs obliquely between them
  const b = solid(700);
  const ring = (z, a) => { const v = b.vertex(Math.cos(a), Math.sin(a), z); return b.circle(v, [0, 0, z], [0, 0, 1], [1, 0, 0]); };
  const e0 = ring(0, 0.3), e1 = ring(4, 2);
  b.face(b.cone([0, 0, 0], [0, 0, 1], [1, 0, 0]), 'forward', [[[e0, false]], [[e1, true]]]);
  b.face(b.plane([0, 0, 4], [0, 0, 1]), 'forward', [[[e1, false]]]);
  b.face(b.plane([0, 0, 0], [0, 0, -1]), 'forward', [[[e0, true]]]);
  const r = parseAcis({ name: 'skew.sat', text: b.text() });
  assert.deepEqual(r.warnings, []);
  tight(surf(r));
  assert.ok(vol(surf(r)) < 4 * Math.PI && vol(surf(r)) > 0.988 * 4 * Math.PI);
});

test('SAT (e): tubular member — two cylinders and annular ends; placed by a transform', () => {
  const ro = 0.6, ri = 0.55, L = 12;
  for (const version of [400, 700]) {
    const r = parseAcis({ name: 'brace.sat', text: tubeSat(version, ro, ri, L) });
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    const q = tight(m);
    assert.equal(m.groups.length, 4);
    assert.deepEqual(r.meta.elementTypes, { cylinder: 2, plane: 2 });
    const n = 24, poly = (rr) => (n / 2) * Math.sin(2 * Math.PI / n) * rr * rr, per = (rr) => 2 * n * rr * Math.sin(Math.PI / n);
    close(vol(m), (poly(ro) - poly(ri)) * L, 2e-5, 'steel volume');
    close(area(m), (per(ro) + per(ri)) * L + 2 * (poly(ro) - poly(ri)), 2e-4, 'area');
    assert.equal(q.stats.genus, 1);
    // every vertex on one of the two walls
    for (let i = 0; i < m.positions.length; i += 3) { const rr = Math.hypot(m.positions[i], m.positions[i + 1]); assert.ok(Math.abs(rr - ro) < 1e-6 || Math.abs(rr - ri) < 1e-6 || (rr < ro && rr > ri * Math.cos(Math.PI / n) - 1e-6)); }
  }
  // the same member laid along x: local z → global x, x → y, y → z, then moved to (5, 0, -20)
  const t = parseAcis({ name: 'brace.sat', text: tubeSat(700, ro, ri, L, [0, 1, 0, 0, 0, 1, 1, 0, 0, 5, 0, -20, 1]) });
  close(t.meta.bbox.min[0], 5, 1e-6); close(t.meta.bbox.max[0], 17, 1e-6);
  close(t.meta.bbox.max[1], ro, 1e-6); close(t.meta.bbox.min[2], -20 - ro, 1e-5);
  tight(surf(t));
  assert.ok(vol(surf(t)) > 0);
});

test('SAT: cones (apex and frustum), elliptical cylinder, torus and pipe elbow', () => {
  // cone of base radius 1 with its apex 2 above the base: radius = R + z·sin/cos, so sin/cos = −1/2
  const s = -1 / Math.sqrt(5), c = 2 / Math.sqrt(5);
  {
    const b = solid(700);
    const v = b.vertex(1, 0, 0), e = b.circle(v, [0, 0, 0], [0, 0, 1], [1, 0, 0]);
    b.face(b.cone([0, 0, 0], [0, 0, 1], [1, 0, 0], s, c), 'forward', [[[e, false]]]);
    b.face(b.plane([0, 0, 0], [0, 0, -1]), 'forward', [[[e, true]]]);
    const r = parseAcis({ name: 'cone.sat', text: b.text() });
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    tight(m);
    close(vol(m), 12 * Math.sin(Math.PI / 12) * 2 / 3, 1e-6, 'pyramid on a 24-gon');
    close(r.meta.bbox.max[2], 2, 1e-6);
    assert.deepEqual(r.meta.elementTypes, { cone: 1, plane: 1 });
  }
  {
    // frustum between z = 0 (radius 1) and z = 1 (radius 0.5)
    const b = solid(400);
    const v0 = b.vertex(1, 0, 0), v1 = b.vertex(0.5, 0, 1);
    const e0 = b.circle(v0, [0, 0, 0], [0, 0, 1], [1, 0, 0]), e1 = b.circle(v1, [0, 0, 1], [0, 0, 1], [0.5, 0, 0]);
    b.face(b.cone([0, 0, 0], [0, 0, 1], [1, 0, 0], s, c), 'forward', [[[e0, false]], [[e1, true]]]);
    b.face(b.plane([0, 0, 1], [0, 0, 1]), 'forward', [[[e1, false]]]);
    b.face(b.plane([0, 0, 0], [0, 0, -1]), 'forward', [[[e0, true]]]);
    const r = parseAcis({ name: 'frustum.sat', text: b.text() });
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    tight(m);
    close(vol(m), 12 * Math.sin(Math.PI / 12) * (1 + 0.5 + 0.25) / 3, 1e-6, 'frustum of 24-gons');
    for (let i = 0; i < m.positions.length; i += 3) assert.ok(Math.hypot(m.positions[i], m.positions[i + 1]) <= 1 - 0.5 * m.positions[i + 2] + 1e-6, 'inside the cone');
  }
  {
    // elliptical cylinder, semi-axes 2 and 1 (ratio 0.5), height 1
    const b = solid(700);
    const v0 = b.vertex(2, 0, 0), v1 = b.vertex(2, 0, 1);
    const e0 = b.circle(v0, [0, 0, 0], [0, 0, 1], [2, 0, 0], 0.5), e1 = b.circle(v1, [0, 0, 1], [0, 0, 1], [2, 0, 0], 0.5);
    b.face(b.cone([0, 0, 0], [0, 0, 1], [2, 0, 0], 0, 1, 0.5), 'forward', [[[e0, false]], [[e1, true]]]);
    b.face(b.plane([0, 0, 1], [0, 0, 1]), 'forward', [[[e1, false]]]);
    b.face(b.plane([0, 0, 0], [0, 0, -1]), 'forward', [[[e0, true]]]);
    const r = parseAcis({ name: 'oval.sat', text: b.text() }, { angularStep: 5 });
    const m = surf(r);
    tight(m);
    close(vol(m), 36 * Math.sin(2 * Math.PI / 72) * 2 * 1, 1e-5, 'affine image of a 72-gon prism');
    close(r.meta.bbox.max[1], 1, 1e-6);
  }
  {
    // full torus: one face without loops
    const b = solid(700);
    b.face(b.torus([0, 0, 0], [0, 0, 1], 3, 1, [1, 0, 0]), 'forward', []);
    const r = parseAcis({ name: 'ring.sat', text: b.text() }, { angularStep: 7.5 });
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    const q = tight(m);
    assert.equal(q.stats.genus, 1);
    const V = 2 * Math.PI ** 2 * 3, A = 4 * Math.PI ** 2 * 3;
    assert.ok(vol(m) < V && vol(m) > 0.985 * V, `torus volume ${vol(m)} / ${V}`);
    assert.ok(area(m) < A && area(m) > 0.99 * A);
  }
  {
    // 90° pipe elbow: quarter of a torus (bend radius 2, pipe radius 0.5) closed by two flat discs
    const R = 2, rr = 0.5;
    const b = solid(700);
    const va = b.vertex(R + rr, 0, 0), vb = b.vertex(0, R + rr, 0);
    const ea = b.circle(va, [R, 0, 0], [0, -1, 0], [rr, 0, 0]), eb = b.circle(vb, [0, R, 0], [1, 0, 0], [0, rr, 0]);
    b.face(b.torus([0, 0, 0], [0, 0, 1], R, rr, [1, 0, 0]), 'forward', [[[ea, true]], [[eb, false]]]);
    b.face(b.plane([R, 0, 0], [0, -1, 0]), 'forward', [[[ea, false]]]);
    b.face(b.plane([0, R, 0], [-1, 0, 0]), 'forward', [[[eb, true]]]);
    const r = parseAcis({ name: 'elbow.sat', text: b.text() }, { angularStep: 7.5 });
    assert.deepEqual(r.warnings, []);
    const m = surf(r);
    tight(m);
    const V = Math.PI ** 2 * R * rr * rr / 2;
    assert.ok(vol(m) < V && vol(m) > 0.985 * V, `elbow volume ${vol(m)} / ${V}`);
    assert.ok(r.meta.bbox.min[0] > -1e-6 && r.meta.bbox.min[1] > -1e-6, 'the quarter with x, y ≥ 0');
    assert.deepEqual(r.meta.elementTypes, { torus: 1, plane: 2 });
  }
});

test('SAT: spline-surface faces — a stored control net is trimmed by the face boundary, procedural ones are counted; nothing crashes', () => {
  // ruled patch (degree 1 × 1, 2 × 3 control points) standing on the edge x = 0 … 2 of a box top and leaning out to
  // y = 0.5 at mid-height; the face is the part of it inside a loop of six straight edges that lie on it. A second spline face has a
  // procedural definition without a stored B-spline.
  const net = '0 0 4 2 0 4 0 0.5 5 2 0.5 5 0 0 6 2 0 6';
  const mk = (version, second, loopPts = [[0, 0, 4], [2, 0, 4], [2, 0.5, 5], [2, 0, 6], [0, 0, 6], [0, 0.5, 5]]) => {
    const b = solid(version);
    addBox(b, [0, 0, 0], [2, 3, 4]);
    const v = loopPts.map((p) => b.vertex(...p));
    const loop = v.map((x, k) => b.between(x, v[(k + 1) % v.length]));
    const nubs = ['{', 'exactsur', ...(version >= 700 ? ['full'] : []), 'nubs', 1, 1, 'open', 'open', 'none', 'none', 2, 3, 0, 1, 1, 1, 0, 1, 1, 1, 2, 1, net, 0, '}'];
    b.face(['spline-surface', 'forward', ...nubs, 'I', 'I', 'I', 'I'], 'forward', [loop]);
    if (second === 'procedural') b.face(['spline-surface', 'forward', '{', 'sweepsur', 'none', 'F', 0, 'F', 1, '}', 'I', 'I', 'I', 'I'], 'forward', []);
    if (second === 'ref') b.face(['spline-surface', 'forward', '{', 'ref', 0, '}', 'I', 'I', 'I', 'I'], 'forward', []);
    return b.text();
  };
  const patchOf = (m, k = 6) => { const g = m.groups[k]; return { positions: m.positions, indices: m.indices.subarray(g.start, g.start + g.count) }; };
  for (const version of [400, 700]) {
    const r = parseAcis({ name: 'mixed.sat', text: mk(version, 'procedural') });
    assertContract(r, ['surface']);
    const m = surf(r);
    assert.equal(m.groups.length, 7, 'six planar faces and the patch');
    assert.ok(r.warnings.some((w) => /^1 of 8 faces could not be read \(1 × procedural spline surface without a stored B-spline \("sweepsur"\)\) — export those bodies as STEP/.test(w)), r.warnings.join(' | '));
    assert.equal(r.meta.skipped['face: procedural spline surface without a stored B-spline ("sweepsur")'], 1);
    assert.equal(r.meta.elementTypes.spline, 1);
    assert.deepEqual([r.meta.faces, r.meta.faceTotal, r.meta.facesSkipped], [7, 8, 1]);
    // the patch reaches y = 0.5 at mid-height and spans z = 4 … 6; its corners are the box-top vertices it shares
    const patch = patchOf(m);
    let ymax = 0, zmax = 0;
    for (const i of patch.indices) { ymax = Math.max(ymax, m.positions[3 * i + 1]); zmax = Math.max(zmax, m.positions[3 * i + 2]); }
    close(ymax, 0.5, 1e-6); close(zmax, 6, 1e-6);
    close(surfaceArea(patch.positions, patch.indices), 2 * 2 * Math.hypot(1, 0.5), 1e-5, 'two flat strips');
  }
  // trimming: a loop over the lower half of the patch only (its upper edge is the knot line z = 5, y = 0.5)
  const half = parseAcis({ name: 'half.sat', text: mk(700, null, [[0, 0, 4], [2, 0, 4], [2, 0.5, 5], [0, 0.5, 5]]) });
  close(surfaceArea(patchOf(surf(half)).positions, patchOf(surf(half)).indices), 2 * Math.hypot(1, 0.5), 1e-5, 'one strip');
  // "{ ref n }" stands for the n-th sub-type block of the file; a spline face without loops covers its surface
  const ref = parseAcis({ name: 'ref.sat', text: mk(700, 'ref') });
  assert.equal(surf(ref).groups.length, 8);
  close(surfaceArea(patchOf(surf(ref), 7).positions, patchOf(surf(ref), 7).indices), 2 * 2 * Math.hypot(1, 0.5), 1e-5, 'whole patch');
  // only free-form geometry: an error that says what to do, not an empty result
  const b = solid(700);
  b.face(['spline-surface', 'forward', '{', 'ref', 0, '}', 'I', 'I', 'I', 'I'], 'forward', []);
  assert.throws(() => parseAcis({ name: 'ff.sat', text: b.text() }), /ff\.sat: no geometry could be read from the SAT file\. 1 of 1 face could not be read \(1 × spline surface on a body without any vertex\) — export those bodies as STEP/);
  // a control net that does not pass by the face's boundary is not used
  const far = mk(700, null).replace(net, '0 0 40 2 0 40 0 0.5 50 2 0.5 50 0 0 60 2 0 60');
  const fr = parseAcis({ name: 'far.sat', text: far });
  assert.equal(surf(fr).groups.length, 6);
  assert.ok(fr.warnings.some((w) => /1 of 7 faces could not be read \(1 × spline surface: its boundary lies \d[\d.e+]* model units off the B-spline surface stored for it\)/.test(w)), fr.warnings.join(' | '));
});

// Wire body, written out in full (7.0 layout): a straight edge, a quarter-circle arc and an edge on an
// intcurve that stores its B-spline (exactcur; degree 2, knots 0 and 1 with multiplicity 2 as the format
// keeps them, three control points — a parabola from (2,1,0) through control point (3,1,0) to (3,2,0)).
const WIRE = (curve) => `700 0 1 0
@7 fixture @11 ACIS 7.0 NT @24 Thu Jan 01 00:00:00 2004
1 9.9999999999999995e-007 1e-010
body $-1 -1 $-1 $-1 $1 $-1 #
wire $-1 -1 $-1 $-1 $2 $0 $-1 out #
coedge $-1 -1 $-1 $3 $-1 $-1 $5 forward $1 $-1 #
coedge $-1 -1 $-1 $4 $2 $-1 $6 forward $1 $-1 #
coedge $-1 -1 $-1 $-1 $3 $-1 $7 forward $1 $-1 #
edge $-1 -1 $-1 $8 0 $9 1000 $2 $16 forward @7 unknown #
edge $-1 -1 $-1 $9 0 $10 1.5707963267948966 $3 $17 forward @7 unknown #
edge $-1 -1 $-1 $10 0 $11 1 $4 $18 forward @7 unknown #
vertex $-1 -1 $-1 $5 $12 #
vertex $-1 -1 $-1 $5 $13 #
vertex $-1 -1 $-1 $6 $14 #
vertex $-1 -1 $-1 $7 $15 #
point $-1 -1 $-1 0 0 0 #
point $-1 -1 $-1 1000 0 0 #
point $-1 -1 $-1 2000 1000 0 #
point $-1 -1 $-1 3000 2000 0 #
straight-curve $-1 -1 $-1 0 0 0 1 0 0 I I #
ellipse-curve $-1 -1 $-1 1000 1000 0 0 0 1 0 -1000 0 1 I I #
${curve}
End-of-ACIS-data
`;
const EXACT = `intcurve-curve $-1 -1 $-1 forward { exactcur full nubs 2 open 2
 0 2 1 2
 2000 1000 0
 3000 1000 0
 3000 2000 0
 0
 null_surface
 null_surface
 nullbs
 nullbs
 I I } I I #`;

test('SAT: a file with only wire bodies gives a frame; arcs and stored B-splines are followed, others become chords', () => {
  const r = parseAcis({ name: 'wire.sat', text: WIRE(EXACT) });
  assertContract(r, ['frame']);
  assert.deepEqual(r.warnings, []);
  assert.equal(r.meta.units.unit, 'mm');
  const f = r.models[0];
  assert.equal(f.name, 'wire');
  // 1 line + 6 arc chords (90° at 15°) + 8 spline chords
  assert.equal(f.members.length, 1 + 6 + 8);
  assert.equal(f.nodes.length / 3, 16);
  const node = (i) => [f.nodes[3 * i], f.nodes[3 * i + 1], f.nodes[3 * i + 2]];
  let onArc = 0, onParabola = 0;
  for (let i = 0; i < f.nodes.length / 3; i++) {
    const [x, y] = node(i);
    if (x > 1 - 1e-9 && x < 2 + 1e-9 && y < 1 + 1e-9 && Math.abs(Math.hypot(x - 1, y - 1) - 1) < 1e-9) onArc++;
    // parabola (2,1) → (3,2) with control point (3,1): x = 2 + 2t − t², y = 1 + t²
    if (x >= 2 - 1e-9) { const t = Math.sqrt(Math.max(0, y - 1)); if (Math.abs(x - (2 + 2 * t - t * t)) < 1e-9) onParabola++; }
  }
  assert.equal(onArc, 7);
  assert.equal(onParabola, 9);
  let len = 0;
  for (const m of f.members) { const a = node(m.a), b = node(m.b); len += Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]); }
  assert.ok(len > 1 + 12 * Math.sin(Math.PI / 24) - 1e-9 + Math.SQRT2 && len < 1 + Math.PI / 2 + 2, `length ${len}`);
  assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [3, 2, 0] });

  // the same edge on a procedural curve without stored data: a straight chord, and a warning that says so
  const c = parseAcis({ name: 'wire.sat', text: WIRE('intcurve-curve $-1 -1 $-1 forward { ref 0 } I I #') });
  assert.equal(c.models[0].members.length, 1 + 6 + 1);
  assert.ok(c.warnings.some((w) => /intcurve-curve without stored B-spline data was replaced by a straight chord/.test(w)), c.warnings.join(' | '));
  // a stored B-spline that does not reach the edge's vertices is not used
  const off = parseAcis({ name: 'wire.sat', text: WIRE(EXACT.replace('3000 2000 0\n 0', '3000 5000 0\n 0')) });
  assert.equal(off.models[0].members.length, 1 + 6 + 1);
  assert.ok(off.warnings.some((w) => /does not pass through its edge vertices/.test(w)));
  // an unknown curve type likewise
  const unk = parseAcis({ name: 'wire.sat', text: WIRE('helix-curve $-1 -1 $-1 0 0 0 0 0 1 I I #') });
  assert.ok(unk.warnings.some((w) => /curve of type "helix-curve" is not read; its edge was replaced by a straight chord/.test(w)));
});

test('SAT: unknown records are reported; damaged, truncated and binary files are rejected with the reason', async () => {
  const good = boxSat(700);
  const extra = good.replace('End-of-ACIS-data', 'widget-gadget $-1 -1 $-1 3 4 5 #\nwidget-gadget $-1 -1 $-1 6 #\nname_attrib-gen-attrib $-1 -1 $-1 $-1 $-1 $0 @4 leg1 #\nEnd-of-ACIS-data');
  const r = parseAcis({ name: 'x.sat', text: extra });
  assert.deepEqual(r.warnings, ['2 × record type "widget-gadget" not recognised and ignored.', 'The header announces 85 records but 88 were found.']);
  close(vol(surf(r)), 24, 1e-12);
  // a face on a surface type this reader does not know is left out and counted, the rest is kept
  const odd = parseAcis({ name: 'x.sat', text: good.replace(/plane-surface \$-1 -1 \$-1 0 0 0 0 0 -1/, 'meshsurf-surface $-1 -1 $-1 0 0 0 0 0 -1') });
  assert.equal(surf(odd).groups.length, 5);
  assert.ok(odd.warnings.some((w) => /^1 of 6 faces could not be read \(1 × surface type "meshsurf-surface"\) — export those bodies as STEP/.test(w)), odd.warnings.join(' | '));
  assert.equal(inspect(surf(odd)).watertight, false);
  // an analytic surface with missing values is skipped with its own message
  const bad = parseAcis({ name: 'x.sat', text: good.replace(/plane-surface \$-1 -1 \$-1 0 0 0 0 0 -1 1 0 0/, 'plane-surface $-1 -1 $-1 0 0 0') });
  assert.equal(surf(bad).groups.length, 5);
  assert.ok(bad.warnings.some((w) => /1 × plane-surface record with missing or unsupported values/.test(w)), bad.warnings.join(' | '));
  // a pointer past the end of the file
  const dangling = parseAcis({ name: 'x.sat', text: boxSat(400).replace(/\n-\d+ straight-curve[^\n]*#(?=\n-\d+ edge)/, '') });
  assert.ok(dangling.models.length >= 1);

  assert.throws(() => parseAcis({ name: 't.sat', text: good.slice(0, good.length >> 1) }), /t\.sat: the SAT file is truncated: record \d+ \(\w[\w-]*\) is not closed by "#"\./);
  assert.throws(() => parseAcis({ name: 'n.sat', text: 'solid cube\nfacet normal 0 0 1\n' }), /n\.sat: not an ACIS SAT text file \(expected a first line of four integers/);
  assert.throws(() => parseAcis({ name: 'h.sat', text: '700 0 1 0\n' }), /h\.sat: the SAT file holds a header but no entity records\./);
  // binary files go to the SAB reader (sab.test.js); the text reader says so when it is handed one
  assert.throws(() => parseAcis({ name: 'b.sab', bytes: enc('ACIS BinaryFile\u0004\0\0\0') }), /b\.sab: this is a binary ACIS file \(SAB\), which the SAT \(text\) reader does not take/);
  assert.equal(formatFor('part.sat').status, 'native');
});
