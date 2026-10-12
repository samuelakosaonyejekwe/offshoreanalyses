// Shared helpers for the geometry importers: byte/text access, growable typed arrays, line iteration over
// large inputs, bounding boxes, polygon triangulation, 4x4 transforms and a finite-element mesh collector that
// turns node/element tables from any solver deck into the model kinds of docs/CONTRACT.md §4.
//
// Everything here is environment-neutral (Node, browser main thread, Web Worker): no DOM, no fetch.

const UTF8 = new TextDecoder('utf-8');

/** Decode bytes to text. Honours UTF-8 / UTF-16 byte-order marks; otherwise UTF-8 (invalid bytes become U+FFFD). */
export function decodeText(bytes) {
  if (bytes.length >= 3 && bytes[0] === 0xef && bytes[1] === 0xbb && bytes[2] === 0xbf) return UTF8.decode(bytes.subarray(3));
  if (bytes.length >= 2 && bytes[0] === 0xff && bytes[1] === 0xfe) return new TextDecoder('utf-16le').decode(bytes.subarray(2));
  if (bytes.length >= 2 && bytes[0] === 0xfe && bytes[1] === 0xff) return new TextDecoder('utf-16be').decode(bytes.subarray(2));
  return UTF8.decode(bytes);
}

/** Text of an importer input ({ text } preferred, else decoded { bytes }). */
export function inputText(input) {
  if (typeof input.text === 'string') return input.text.charCodeAt(0) === 0xfeff ? input.text.slice(1) : input.text;
  if (input.bytes) return decodeText(input.bytes);
  throw new Error(`${input.name || 'input'}: no content supplied (expected text or bytes).`);
}

/** Bytes of an importer input ({ bytes } preferred, else UTF-8 encoded { text }). */
export function inputBytes(input) {
  if (input.bytes) return input.bytes;
  if (typeof input.text === 'string') return new TextEncoder().encode(input.text);
  throw new Error(`${input.name || 'input'}: no content supplied (expected bytes or text).`);
}

/** File name without directory and extension. */
export function baseName(name) {
  const s = String(name || 'model').replace(/\\/g, '/');
  const f = s.slice(s.lastIndexOf('/') + 1);
  const d = f.lastIndexOf('.');
  return d > 0 ? f.slice(0, d) : f || 'model';
}

/** Lower-case extension without the dot ('' when none). */
export function extOf(name) {
  const s = String(name || '').replace(/\\/g, '/');
  const f = s.slice(s.lastIndexOf('/') + 1);
  const d = f.lastIndexOf('.');
  return d >= 0 ? f.slice(d + 1).toLowerCase() : '';
}

/** Look a file up in a multi-file input by exact path, then by base file name, case-insensitively. */
export function findFile(files, wanted) {
  if (!files) return null;
  if (files.has(wanted)) return files.get(wanted);
  const w = String(wanted).replace(/\\/g, '/').toLowerCase();
  const wb = w.slice(w.lastIndexOf('/') + 1);
  let byBase = null;
  for (const [k, v] of files) {
    const kk = k.replace(/\\/g, '/').toLowerCase();
    if (kk === w || kk.endsWith('/' + w)) return v;
    if (!byBase && kk.slice(kk.lastIndexOf('/') + 1) === wb) byBase = v;
  }
  return byBase;
}

/** Growable typed array. `finish()` returns an exact-length copy. */
export class Grow {
  constructor(Type = Float64Array, capacity = 1024) {
    this.a = new Type(Math.max(4, capacity));
    this.n = 0;
  }
  ensure(extra) {
    if (this.n + extra <= this.a.length) return;
    let cap = this.a.length * 2;
    while (cap < this.n + extra) cap *= 2;
    const b = new this.a.constructor(cap);
    b.set(this.a.subarray(0, this.n));
    this.a = b;
  }
  push(v) {
    if (this.n === this.a.length) this.ensure(1);
    this.a[this.n++] = v;
  }
  push2(x, y) {
    this.ensure(2);
    this.a[this.n++] = x; this.a[this.n++] = y;
  }
  push3(x, y, z) {
    this.ensure(3);
    this.a[this.n++] = x; this.a[this.n++] = y; this.a[this.n++] = z;
  }
  get length() { return this.n; }
  view() { return this.a.subarray(0, this.n); }
  finish() { return this.a.slice(0, this.n); }
}

/**
 * Line iterator over a string or a byte buffer. Over bytes it never materialises the whole file as a string,
 * so multi-hundred-megabyte ASCII meshes can be read line by line. Handles LF, CRLF and lone CR.
 * `offset` is the position of the next unread byte/char (use it to switch to binary reads).
 */
export class LineCursor {
  constructor(source) {
    if (typeof source === 'string') { this.s = source; this.b = null; this.len = source.length; } else { this.b = source; this.s = null; this.len = source.length; }
    this.offset = 0;
    this.lineNo = 0;
    this.win = null; this.w0 = 0; this.w1 = 0; this.slowUntil = 0;
    if (this.b && this.len >= 3 && this.b[0] === 0xef && this.b[1] === 0xbb && this.b[2] === 0xbf) this.offset = 3;
    if (this.s && this.s.charCodeAt(0) === 0xfeff) this.offset = 1;
  }
  get done() { return this.offset >= this.len; }
  /** Next line without its terminator, or null at end of input. */
  next() {
    if (this.offset >= this.len) return null;
    this.lineNo++;
    const start = this.offset;
    if (this.s !== null) {
      let e = this.s.indexOf('\n', start);
      if (e < 0) e = this.len;
      this.offset = e + 1;
      if (e > start && this.s.charCodeAt(e - 1) === 13) e--;
      return this.s.slice(start, e);
    }
    const b = this.b;
    // Fast path: decode a 64 kB window once and cut lines out of it. Valid only while the window is pure
    // ASCII, so that character positions equal byte positions.
    if (start >= this.slowUntil) {
      if (!(start >= this.w0 && start < this.w1)) this.loadWindow(start);
      if (this.win !== null) {
        const rel = start - this.w0;
        let e = this.win.indexOf('\n', rel);
        if (e < 0 && this.w1 < this.len && rel > 0) { this.loadWindow(start); if (this.win !== null) e = this.win.indexOf('\n', 0); }
        if (this.win !== null && (e >= 0 || this.w1 >= this.len)) {
          const r0 = start - this.w0;
          if (e < 0) e = this.win.length;
          let line = this.win.slice(r0, e);
          let adv = e + 1;
          const cr = line.indexOf('\r');
          if (cr >= 0) {
            // CRLF ends here; a CR elsewhere is an old-style line break of its own
            if (cr === line.length - 1) line = line.slice(0, cr); else { line = line.slice(0, cr); adv = r0 + cr + 1; }
          }
          this.offset = this.w0 + adv;
          return line;
        }
        // a line longer than the window: scan it byte by byte below
      }
    }
    let e = start, hi = 0;
    while (e < this.len) {
      const c = b[e];
      if (c === 10 || c === 13) break;
      hi |= c;
      e++;
    }
    this.offset = e < this.len && b[e] === 13 && e + 1 < this.len && b[e + 1] === 10 ? e + 2 : e + 1;
    if (hi < 128 && e - start <= 8192) return String.fromCharCode.apply(null, b.subarray(start, e));
    return UTF8.decode(b.subarray(start, e));
  }
  loadWindow(start) {
    const end = Math.min(this.len, start + 65536);
    const text = UTF8.decode(this.b.subarray(start, end));
    if (text.length === end - start) { this.win = text; this.w0 = start; this.w1 = end; } else { this.win = null; this.w0 = 0; this.w1 = 0; this.slowUntil = end; }
  }
  /** Next line that is not blank (after trimming), trimmed; null at end. */
  nextNonBlank() {
    for (;;) {
      const l = this.next();
      if (l === null) return null;
      const t = l.trim();
      if (t) return t;
    }
  }
}

/** Line cursor over an importer input, preferring bytes (no whole-file string). */
export function lineCursor(input) {
  return new LineCursor(input.bytes ? input.bytes : inputText(input));
}

/** Split a line on whitespace and/or commas into numbers (NaN for unparsable tokens). */
export function parseNumbers(line) {
  const t = line.trim();
  if (!t) return [];
  const p = t.split(/[\s,;]+/);
  const out = new Array(p.length);
  for (let i = 0; i < p.length; i++) out[i] = Number(p[i]);
  return out;
}

/** Fortran-tolerant float: accepts 1.5D+3, 1.5+3 (exponent without E, Nastran style), blank → NaN. */
export function fortranFloat(s) {
  let t = s.trim();
  if (!t) return NaN;
  t = t.replace(/[dD]/, 'e');
  if (!/[eE]/.test(t)) {
    const m = /^([+-]?[\d.]+)([+-]\d+)$/.exec(t);
    if (m) t = m[1] + 'e' + m[2];
  }
  return Number(t);
}

// ---------------------------------------------------------------------------------------------------------
// Bounding boxes, precision-safe single-precision packing
// ---------------------------------------------------------------------------------------------------------

/** Axis-aligned bounds of a flat coordinate array. Non-finite values are ignored. */
export function bboxOf(arr, stride = 3, count = arr.length) {
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i + stride <= count; i += stride) {
    for (let k = 0; k < stride && k < 3; k++) {
      const v = arr[i + k];
      if (v === v && v !== Infinity && v !== -Infinity) { if (v < min[k]) min[k] = v; if (v > max[k]) max[k] = v; }
    }
  }
  for (let k = 0; k < 3; k++) if (min[k] === Infinity) { min[k] = 0; max[k] = 0; }
  return { min, max };
}

/**
 * Pack double-precision coordinates into the Float32Array required for `surface` and `points` models.
 * Single precision resolves only ~0.5 m at a UTM northing of 6 000 000 m, so when any |coordinate| exceeds
 * 10 km the data are stored relative to an integer-metre `origin` (returned) — world = positions + origin.
 */
export function packPositions(src, count = src.length) {
  const out = new Float32Array(count);
  let big = 0;
  for (let i = 0; i < count; i++) { const v = Math.abs(src[i]); if (v > big && v !== Infinity) big = v; }
  if (!(big > 1e4)) {
    for (let i = 0; i < count; i++) out[i] = src[i];
    return { positions: out };
  }
  const bb = bboxOf(src, 3, count);
  const o = [Math.round((bb.min[0] + bb.max[0]) / 2), Math.round((bb.min[1] + bb.max[1]) / 2), Math.round((bb.min[2] + bb.max[2]) / 2)];
  for (let i = 0; i + 2 < count; i += 3) { out[i] = src[i] - o[0]; out[i + 1] = src[i + 1] - o[1]; out[i + 2] = src[i + 2] - o[2]; }
  return { positions: out, origin: o };
}

/** Build a `surface` model from double-precision coordinates (applies packPositions). */
export function surfaceModel(name, coords, indices, groups, extra) {
  const p = packPositions(coords, coords.length);
  const m = { kind: 'surface', name, positions: p.positions, indices: indices instanceof Uint32Array ? indices : Uint32Array.from(indices) };
  if (p.origin) m.origin = p.origin;
  if (groups && groups.length) m.groups = groups;
  return extra ? Object.assign(m, extra) : m;
}

/** Build a `points` model from double-precision coordinates (applies packPositions). */
export function pointsModel(name, coords, attrs, extra) {
  const p = packPositions(coords, coords.length);
  const m = { kind: 'points', name, positions: p.positions };
  if (p.origin) m.origin = p.origin;
  if (attrs && Object.keys(attrs).length) m.attrs = attrs;
  return extra ? Object.assign(m, extra) : m;
}

/** Bounds of one model in world coordinates (origin applied). */
export function modelBounds(m) {
  let bb;
  const o = m.origin || [0, 0, 0];
  switch (m.kind) {
    case 'surface': case 'volume': case 'points': bb = bboxOf(m.positions); break;
    case 'frame': bb = bboxOf(m.nodes); break;
    case 'polyline': bb = bboxOf(m.points); break;
    case 'grid': {
      let lo = Infinity, hi = -Infinity;
      for (let i = 0; i < m.z.length; i++) { const v = m.z[i]; if (v === v) { if (v < lo) lo = v; if (v > hi) hi = v; } }
      if (lo === Infinity) { lo = 0; hi = 0; }
      return { min: [m.x0, m.y0, lo], max: [m.x0 + (m.nx - 1) * m.dx, m.y0 + (m.ny - 1) * m.dy, hi] };
    }
    case 'section': {
      const min = [Infinity, Infinity, 0], max = [-Infinity, -Infinity, 0];
      for (const l of m.loops || []) { const b = bboxOf(l, 2); for (let k = 0; k < 2; k++) { if (b.min[k] < min[k]) min[k] = b.min[k]; if (b.max[k] > max[k]) max[k] = b.max[k]; } }
      if (min[0] === Infinity) return { min: [0, 0, 0], max: [0, 0, 0] };
      return { min, max };
    }
    default: return { min: [0, 0, 0], max: [0, 0, 0] };
  }
  return { min: [bb.min[0] + o[0], bb.min[1] + o[1], bb.min[2] + o[2]], max: [bb.max[0] + o[0], bb.max[1] + o[1], bb.max[2] + o[2]] };
}

/** Entity counts of one model. */
export function modelCounts(m) {
  switch (m.kind) {
    case 'surface': return { vertices: m.positions.length / 3, triangles: m.indices.length / 3 };
    case 'volume': { let c = 0; for (const b of m.cells) c += b.count ?? b.conn.length / ELEM[b.type].n; return { vertices: m.positions.length / 3, cells: c }; }
    case 'frame': return { nodes: m.nodes.length / 3, members: m.members.length };
    case 'polyline': return { points: m.points.length / 3 };
    case 'grid': return { cells: m.nx * m.ny };
    case 'points': return { points: m.positions.length / 3 };
    case 'section': { let n = 0; for (const l of m.loops) n += l.length / 2; return { loops: m.loops.length, points: n }; }
    default: return {};
  }
}

/**
 * Standard importer return value. Fills `meta.counts` (totals by entity and by model kind), `meta.bbox`
 * (world coordinates, metres) and leaves everything the importer already put in `meta` untouched.
 */
export function finish(models, warnings = [], meta = {}) {
  const counts = { models: models.length, byKind: {} };
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (const m of models) {
    counts.byKind[m.kind] = (counts.byKind[m.kind] || 0) + 1;
    const c = modelCounts(m);
    for (const k in c) counts[k] = (counts[k] || 0) + c[k];
    const b = modelBounds(m);
    for (let k = 0; k < 3; k++) { if (b.min[k] < min[k]) min[k] = b.min[k]; if (b.max[k] > max[k]) max[k] = b.max[k]; }
  }
  if (!models.length) { min.fill(0); max.fill(0); }
  meta.counts = counts;
  meta.bbox = { min, max };
  return { models, warnings, meta };
}

/** Add to a histogram object. */
export function tally(hist, key, n = 1) {
  hist[key] = (hist[key] || 0) + n;
  return hist;
}

/** Collapse repeated warnings: call `w.add(text)` many times, `w.list()` yields "text (×N)". */
export class Warnings {
  constructor() { this.m = new Map(); }
  add(text) { this.m.set(text, (this.m.get(text) || 0) + 1); }
  list() { return [...this.m].map(([t, n]) => (n > 1 ? `${t} (×${n})` : t)); }
}

// ---------------------------------------------------------------------------------------------------------
// 4x4 transforms (column-major, as in glTF): p' = M · p
// ---------------------------------------------------------------------------------------------------------

export const M4 = {
  identity() { const m = new Float64Array(16); m[0] = m[5] = m[10] = m[15] = 1; return m; },
  /** a · b (apply b first, then a). */
  mul(a, b) {
    const o = new Float64Array(16);
    for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) {
      o[c * 4 + r] = a[r] * b[c * 4] + a[4 + r] * b[c * 4 + 1] + a[8 + r] * b[c * 4 + 2] + a[12 + r] * b[c * 4 + 3];
    }
    return o;
  },
  translation(x, y, z) { const m = M4.identity(); m[12] = x; m[13] = y; m[14] = z; return m; },
  scaling(x, y = x, z = x) { const m = M4.identity(); m[0] = x; m[5] = y; m[10] = z; return m; },
  /** Rotation by `angle` (rad) about a unit-normalised axis (right-hand rule). */
  rotation(ax, ay, az, angle) {
    const l = Math.hypot(ax, ay, az) || 1; ax /= l; ay /= l; az /= l;
    const c = Math.cos(angle), s = Math.sin(angle), t = 1 - c;
    const m = M4.identity();
    m[0] = t * ax * ax + c; m[1] = t * ax * ay + s * az; m[2] = t * ax * az - s * ay;
    m[4] = t * ax * ay - s * az; m[5] = t * ay * ay + c; m[6] = t * ay * az + s * ax;
    m[8] = t * ax * az + s * ay; m[9] = t * ay * az - s * ax; m[10] = t * az * az + c;
    return m;
  },
  /** From translation [x,y,z], unit quaternion [x,y,z,w] and scale [x,y,z]. */
  fromTRS(t = [0, 0, 0], q = [0, 0, 0, 1], s = [1, 1, 1]) {
    const [x, y, z, w] = q;
    const m = M4.identity();
    m[0] = (1 - 2 * (y * y + z * z)) * s[0]; m[1] = 2 * (x * y + z * w) * s[0]; m[2] = 2 * (x * z - y * w) * s[0];
    m[4] = 2 * (x * y - z * w) * s[1]; m[5] = (1 - 2 * (x * x + z * z)) * s[1]; m[6] = 2 * (y * z + x * w) * s[1];
    m[8] = 2 * (x * z + y * w) * s[2]; m[9] = 2 * (y * z - x * w) * s[2]; m[10] = (1 - 2 * (x * x + y * y)) * s[2];
    m[12] = t[0]; m[13] = t[1]; m[14] = t[2];
    return m;
  },
  /** Frame with the given origin and (not necessarily normalised) axis directions as columns. */
  fromAxes(o, x, y, z) {
    const m = M4.identity();
    m[0] = x[0]; m[1] = x[1]; m[2] = x[2];
    m[4] = y[0]; m[5] = y[1]; m[6] = y[2];
    m[8] = z[0]; m[9] = z[1]; m[10] = z[2];
    m[12] = o[0]; m[13] = o[1]; m[14] = o[2];
    return m;
  },
  /** Transform a point; writes into `out` (array-like) at `off` and returns it. */
  apply(m, x, y, z, out = [0, 0, 0], off = 0) {
    out[off] = m[0] * x + m[4] * y + m[8] * z + m[12];
    out[off + 1] = m[1] * x + m[5] * y + m[9] * z + m[13];
    out[off + 2] = m[2] * x + m[6] * y + m[10] * z + m[14];
    return out;
  },
  /** Transform a direction (no translation). */
  applyDir(m, x, y, z, out = [0, 0, 0]) {
    out[0] = m[0] * x + m[4] * y + m[8] * z;
    out[1] = m[1] * x + m[5] * y + m[9] * z;
    out[2] = m[2] * x + m[6] * y + m[10] * z;
    return out;
  },
  isIdentity(m) {
    for (let i = 0; i < 16; i++) if (Math.abs(m[i] - (i % 5 === 0 ? 1 : 0)) > 1e-12) return false;
    return true;
  },
};

export const V3 = {
  sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
  add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
  scale: (a, s) => [a[0] * s, a[1] * s, a[2] * s],
  dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
  cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
  len: (a) => Math.hypot(a[0], a[1], a[2]),
  norm(a) { const l = Math.hypot(a[0], a[1], a[2]); return l > 0 ? [a[0] / l, a[1] / l, a[2] / l] : [0, 0, 0]; },
};

// ---------------------------------------------------------------------------------------------------------
// Polygon triangulation: ear clipping with hole bridging (after Eberly, "Triangulation by Ear Clipping", and
// the linked-list formulation of Mapbox earcut, without its spatial hash — cost is O(n²) per ring).
// ---------------------------------------------------------------------------------------------------------

/**
 * Triangulate a polygon with holes.
 * @param {ArrayLike<number>} data flat coordinates, `dim` numbers per vertex (only the first two are used)
 * @param {number[]} [holeIndices] vertex index at which each hole ring starts
 * @param {number} [dim=2]
 * @returns {number[]} vertex indices, three per triangle (counter-clockwise in the x–y plane)
 */
export function earclip(data, holeIndices, dim = 2) {
  const hasHoles = holeIndices && holeIndices.length;
  const outerLen = hasHoles ? holeIndices[0] * dim : data.length;
  let outer = ringList(data, 0, outerLen, dim, true);
  const tris = [];
  if (!outer || outer.next === outer.prev) return tris;
  if (hasHoles) outer = bridgeHoles(data, holeIndices, outer, dim);
  clipEars(outer, tris, 0);
  return tris;
}

function ringList(data, start, end, dim, clockwise) {
  let last = null;
  let area = 0;
  for (let i = start, j = end - dim; i < end; i += dim) { area += (data[j] - data[i]) * (data[i + 1] + data[j + 1]); j = i; }
  if (clockwise === area > 0) for (let i = start; i < end; i += dim) last = insertNode(i / dim, data[i], data[i + 1], last);
  else for (let i = end - dim; i >= start; i -= dim) last = insertNode(i / dim, data[i], data[i + 1], last);
  if (last && ptEq(last, last.next)) { removeNode(last); last = last.next; }
  return last;
}

function insertNode(i, x, y, last) {
  const p = { i, x, y, prev: null, next: null, steiner: false };
  if (!last) { p.prev = p; p.next = p; } else { p.next = last.next; p.prev = last; last.next.prev = p; last.next = p; }
  return p;
}
function removeNode(p) { p.next.prev = p.prev; p.prev.next = p.next; }
const ptEq = (a, b) => a.x === b.x && a.y === b.y;
const triArea = (p, q, r) => (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
const sgn = (v) => (v > 0 ? 1 : v < 0 ? -1 : 0);

function filterPoints(start, end) {
  if (!start) return start;
  if (!end) end = start;
  let p = start, again;
  do {
    again = false;
    if (!p.steiner && (ptEq(p, p.next) || triArea(p.prev, p, p.next) === 0)) {
      removeNode(p);
      p = end = p.prev;
      if (p === p.next) break;
      again = true;
    } else p = p.next;
  } while (again || p !== end);
  return end;
}

function clipEars(ear, tris, pass) {
  if (!ear) return;
  let stop = ear;
  while (ear.prev !== ear.next) {
    const prev = ear.prev, next = ear.next;
    if (isEar(ear)) {
      tris.push(prev.i, ear.i, next.i);
      removeNode(ear);
      ear = next.next;
      stop = next.next;
      continue;
    }
    ear = next;
    if (ear === stop) {
      if (!pass) clipEars(filterPoints(ear), tris, 1);
      else if (pass === 1) { ear = cureLocalIntersections(filterPoints(ear), tris); clipEars(ear, tris, 2); } else if (pass === 2) splitAndClip(ear, tris);
      break;
    }
  }
}

function isEar(ear) {
  const a = ear.prev, b = ear, c = ear.next;
  if (triArea(a, b, c) >= 0) return false;
  let p = c.next;
  while (p !== a) {
    if (pointInTri(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) && triArea(p.prev, p, p.next) >= 0) return false;
    p = p.next;
  }
  return true;
}

function pointInTri(ax, ay, bx, by, cx, cy, px, py) {
  return (cx - px) * (ay - py) >= (ax - px) * (cy - py) && (ax - px) * (by - py) >= (bx - px) * (ay - py) && (bx - px) * (cy - py) >= (cx - px) * (by - py);
}

function cureLocalIntersections(start, tris) {
  let p = start;
  do {
    const a = p.prev, b = p.next.next;
    if (!ptEq(a, b) && segsIntersect(a, p, p.next, b) && locallyInside(a, b) && locallyInside(b, a)) {
      tris.push(a.i, p.i, b.i);
      removeNode(p); removeNode(p.next);
      p = start = b;
    }
    p = p.next;
  } while (p !== start);
  return filterPoints(p);
}

function splitAndClip(start, tris) {
  let a = start;
  do {
    let b = a.next.next;
    while (b !== a.prev) {
      if (a.i !== b.i && validDiagonal(a, b)) {
        let c = splitPolygon(a, b);
        a = filterPoints(a, a.next);
        c = filterPoints(c, c.next);
        clipEars(a, tris, 0);
        clipEars(c, tris, 0);
        return;
      }
      b = b.next;
    }
    a = a.next;
  } while (a !== start);
}

function bridgeHoles(data, holeIndices, outer, dim) {
  const queue = [];
  for (let i = 0; i < holeIndices.length; i++) {
    const start = holeIndices[i] * dim;
    const end = i < holeIndices.length - 1 ? holeIndices[i + 1] * dim : data.length;
    const list = ringList(data, start, end, dim, false);
    if (!list) continue;
    if (list === list.next) list.steiner = true;
    let p = list, left = list;
    do { if (p.x < left.x || (p.x === left.x && p.y < left.y)) left = p; p = p.next; } while (p !== list);
    queue.push(left);
  }
  queue.sort((a, b) => a.x - b.x);
  for (const h of queue) outer = bridgeHole(h, outer);
  return outer;
}

function bridgeHole(hole, outer) {
  const bridge = findBridge(hole, outer);
  if (!bridge) return outer;
  const rev = splitPolygon(bridge, hole);
  filterPoints(rev, rev.next);
  return filterPoints(bridge, bridge.next);
}

function findBridge(hole, outer) {
  let p = outer;
  const hx = hole.x, hy = hole.y;
  let qx = -Infinity, m = null;
  do {
    if (hy <= p.y && hy >= p.next.y && p.next.y !== p.y) {
      const x = p.x + ((hy - p.y) * (p.next.x - p.x)) / (p.next.y - p.y);
      if (x <= hx && x > qx) {
        qx = x;
        m = p.x < p.next.x ? p : p.next;
        if (x === hx) return m;
      }
    }
    p = p.next;
  } while (p !== outer);
  if (!m) return null;
  const stop = m, mx = m.x, my = m.y;
  let tanMin = Infinity;
  p = m;
  do {
    if (hx >= p.x && p.x >= mx && hx !== p.x && pointInTri(hy < my ? hx : qx, hy, mx, my, hy < my ? qx : hx, hy, p.x, p.y)) {
      const tan = Math.abs(hy - p.y) / (hx - p.x);
      if (locallyInside(p, hole) && (tan < tanMin || (tan === tanMin && (p.x > m.x || (p.x === m.x && triArea(m.prev, m, p.prev) < 0 && triArea(p.next, m, m.next) < 0))))) { m = p; tanMin = tan; }
    }
    p = p.next;
  } while (p !== stop);
  return m;
}

function segsIntersect(p1, q1, p2, q2) {
  const o1 = sgn(triArea(p1, q1, p2)), o2 = sgn(triArea(p1, q1, q2)), o3 = sgn(triArea(p2, q2, p1)), o4 = sgn(triArea(p2, q2, q1));
  if (o1 !== o2 && o3 !== o4) return true;
  if (o1 === 0 && onSeg(p1, p2, q1)) return true;
  if (o2 === 0 && onSeg(p1, q2, q1)) return true;
  if (o3 === 0 && onSeg(p2, p1, q2)) return true;
  if (o4 === 0 && onSeg(p2, q1, q2)) return true;
  return false;
}
const onSeg = (p, q, r) => q.x <= Math.max(p.x, r.x) && q.x >= Math.min(p.x, r.x) && q.y <= Math.max(p.y, r.y) && q.y >= Math.min(p.y, r.y);

function validDiagonal(a, b) {
  return a.next.i !== b.i && a.prev.i !== b.i && !crossesPolygon(a, b) &&
    ((locallyInside(a, b) && locallyInside(b, a) && middleInside(a, b) && (triArea(a.prev, a, b.prev) || triArea(a, b.prev, b))) ||
      (ptEq(a, b) && triArea(a.prev, a, a.next) > 0 && triArea(b.prev, b, b.next) > 0));
}
function crossesPolygon(a, b) {
  let p = a;
  do {
    if (p.i !== a.i && p.next.i !== a.i && p.i !== b.i && p.next.i !== b.i && segsIntersect(p, p.next, a, b)) return true;
    p = p.next;
  } while (p !== a);
  return false;
}
function locallyInside(a, b) {
  return triArea(a.prev, a, a.next) < 0 ? triArea(a, b, a.next) >= 0 && triArea(a, a.prev, b) >= 0 : triArea(a, b, a.prev) < 0 || triArea(a, a.next, b) < 0;
}
function middleInside(a, b) {
  let p = a, inside = false;
  const px = (a.x + b.x) / 2, py = (a.y + b.y) / 2;
  do {
    if (p.y > py !== p.next.y > py && p.next.y !== p.y && px < ((p.next.x - p.x) * (py - p.y)) / (p.next.y - p.y) + p.x) inside = !inside;
    p = p.next;
  } while (p !== a);
  return inside;
}
function splitPolygon(a, b) {
  const a2 = { i: a.i, x: a.x, y: a.y, prev: null, next: null, steiner: false };
  const b2 = { i: b.i, x: b.x, y: b.y, prev: null, next: null, steiner: false };
  const an = a.next, bp = b.prev;
  a.next = b; b.prev = a;
  a2.next = an; an.prev = a2;
  b2.next = a2; a2.prev = b2;
  bp.next = b2; b2.prev = bp;
  return b2;
}

/**
 * Triangulate a planar (or nearly planar) polygon given in 3-D. The rings are projected on the plane of
 * the Newell normal of the outer ring; triangles keep the winding of the outer ring.
 * @param {ArrayLike<number>} xyz flat x,y,z of all rings (outer first)
 * @param {number[]} [holeIndices] vertex index at which each hole starts
 * @returns {number[]} triangle vertex indices into the input vertices
 */
export function triangulate3D(xyz, holeIndices) {
  const n = xyz.length / 3;
  const outerEnd = holeIndices && holeIndices.length ? holeIndices[0] : n;
  if (outerEnd < 3) return [];
  if (outerEnd === 3 && n === 3) return [0, 1, 2];
  let nx = 0, ny = 0, nz = 0;
  for (let i = 0, j = outerEnd - 1; i < outerEnd; j = i++) {
    const ax = xyz[3 * j], ay = xyz[3 * j + 1], az = xyz[3 * j + 2], bx = xyz[3 * i], by = xyz[3 * i + 1], bz = xyz[3 * i + 2];
    nx += (ay - by) * (az + bz); ny += (az - bz) * (ax + bx); nz += (ax - bx) * (ay + by);
  }
  const l = Math.hypot(nx, ny, nz);
  if (!(l > 0)) return [];
  nx /= l; ny /= l; nz /= l;
  // In-plane basis (u, v) with u × v = n so that a ring that is counter-clockwise about n stays so in 2-D.
  let ux, uy, uz;
  if (Math.abs(nx) < 0.9) { ux = 0; uy = nz; uz = -ny; } else { ux = -nz; uy = 0; uz = nx; }
  const ul = Math.hypot(ux, uy, uz); ux /= ul; uy /= ul; uz /= ul;
  const vx = ny * uz - nz * uy, vy = nz * ux - nx * uz, vz = nx * uy - ny * ux;
  const flat = new Float64Array(2 * n);
  for (let i = 0; i < n; i++) {
    const x = xyz[3 * i], y = xyz[3 * i + 1], z = xyz[3 * i + 2];
    flat[2 * i] = x * ux + y * uy + z * uz;
    flat[2 * i + 1] = x * vx + y * vy + z * vz;
  }
  return earclip(flat, holeIndices, 2);
}

/** Signed area of a 2-D ring (positive when counter-clockwise). */
export function ringArea(flat, dim = 2) {
  let a = 0;
  for (let i = 0, j = flat.length - dim; i < flat.length; j = i, i += dim) a += flat[j] * flat[i + 1] - flat[i] * flat[j + 1];
  return a / 2;
}

// ---------------------------------------------------------------------------------------------------------
// Finite-element mesh collector
// ---------------------------------------------------------------------------------------------------------

/**
 * Canonical element types. Node ordering is the VTK / Abaqus ordering:
 *  line3  = end, end, mid
 *  tri6   = 3 corners, mids of edges 01, 12, 20
 *  quad8  = 4 corners, mids of edges 01, 12, 23, 30 (quad9 adds the centre)
 *  tet10  = 4 corners, mids of 01, 12, 20, 03, 13, 23
 *  hex20  = 8 corners (bottom face 0-3, top face 4-7), mids of bottom edges, top edges, vertical edges
 *  wedge6 = bottom triangle 0-2, top triangle 3-5;  pyr5 = base quad 0-3, apex 4
 * Importers must permute their native ordering into this one (see `permute`).
 */
export const ELEM = {
  point1: { n: 1, dim: 0 },
  line2: { n: 2, dim: 1 }, line3: { n: 3, dim: 1 },
  tri3: { n: 3, dim: 2 }, tri6: { n: 6, dim: 2 },
  quad4: { n: 4, dim: 2 }, quad8: { n: 8, dim: 2 }, quad9: { n: 9, dim: 2 },
  tet4: { n: 4, dim: 3 }, tet10: { n: 10, dim: 3 },
  hex8: { n: 8, dim: 3 }, hex20: { n: 20, dim: 3 }, hex27: { n: 27, dim: 3, reduce: 'hex20' },
  wedge6: { n: 6, dim: 3 }, wedge15: { n: 15, dim: 3, reduce: 'wedge6' }, wedge18: { n: 18, dim: 3, reduce: 'wedge6' },
  pyr5: { n: 5, dim: 3 }, pyr13: { n: 13, dim: 3, reduce: 'pyr5' }, pyr14: { n: 14, dim: 3, reduce: 'pyr5' },
};

const SURF_TRIS = {
  tri3: [0, 1, 2],
  tri6: [0, 3, 5, 3, 1, 4, 5, 4, 2, 3, 4, 5],
  quad4: [0, 1, 2, 0, 2, 3],
  quad8: [0, 4, 7, 4, 1, 5, 5, 2, 6, 6, 3, 7, 4, 5, 6, 4, 6, 7],
  quad9: [0, 4, 8, 4, 1, 8, 1, 5, 8, 5, 2, 8, 2, 6, 8, 6, 3, 8, 3, 7, 8, 7, 0, 8],
};

/** Reorder `nodes` by `order` (order[k] = position in the native list of canonical node k). */
export function permute(nodes, order) {
  const out = new Array(order.length);
  for (let i = 0; i < order.length; i++) out[i] = nodes[order[i]];
  return out;
}

/**
 * Collects nodes, elements, sets and property records from a finite-element file and builds §4 models:
 * line elements → `frame`, shell/face elements → `surface`, solid elements → `volume`.
 *
 * Extra (non-contract) fields written on the models so that nothing the solver needs is lost:
 *  - every model: `nodeIds` (Uint32Array, source node number of each model node) and
 *    `sets: { nodes: { name: Uint32Array(model node indices) }, ... }`
 *  - frame members: `id` (source element number), `group` (property/set name), `section`, `material`,
 *    plus anything registered through `memberExtra` (orientation vector, offsets, releases); frame
 *    `sets.members[name]` = member indices
 *  - surface: `elements: [{ type, conn, ids, tags }]` — the untriangulated shells (tri3/tri6/quad4/quad8/quad9)
 *    indexing the same `positions`; `groups` follow the property/physical tag; `sets.elements[name]` = indices
 *    into the concatenated `elements` blocks
 *  - volume: `cells: [{ type, conn, ids, tags }]`; `sets.cells[name]` = indices into the concatenated blocks;
 *    `positions` is a Float64Array (solver precision)
 */
export class FeMesh {
  constructor(name = 'mesh') {
    this.name = name;
    this.nodeIndex = new Map();
    this.nodeIdList = new Grow(Uint32Array, 1024);
    this.coords = new Grow(Float64Array, 3072);
    this.blocks = new Map();
    this.nodeSets = new Map();
    this.elemSets = new Map();
    this.props = new Map();
    this.memberExtra = new Map();
    this.skipped = {};
    this.nElems = 0;
  }

  get nodeCount() { return this.nodeIdList.n; }

  /** Add (or overwrite) a node. Coordinates in file units. */
  addNode(id, x, y, z) {
    const have = this.nodeIndex.get(id);
    if (have !== undefined) { const a = this.coords.a; a[3 * have] = x; a[3 * have + 1] = y; a[3 * have + 2] = z; return have; }
    const i = this.nodeIdList.n;
    this.nodeIndex.set(id, i);
    this.nodeIdList.push(id);
    this.coords.push3(x, y, z);
    return i;
  }

  /**
   * Add an element in canonical node order.
   * @param {string} type key of ELEM
   * @param {ArrayLike<number>} nodes source node numbers
   * @param {number} [id] source element number
   * @param {number} [tag] property / physical-group / block number used for grouping
   */
  addElement(type, nodes, id = 0, tag = 0) {
    const def = ELEM[type];
    if (!def) { tally(this.skipped, type); return; }
    let blk = this.blocks.get(type);
    if (!blk) {
      blk = { type, n: def.n, conn: new Grow(Uint32Array, 1024 * def.n), ids: new Grow(Uint32Array, 1024), tags: new Grow(Int32Array, 1024) };
      this.blocks.set(type, blk);
    }
    blk.conn.ensure(def.n);
    for (let k = 0; k < def.n; k++) blk.conn.a[blk.conn.n++] = nodes[k];
    blk.ids.push(id || ++this.nElems);
    if (id) this.nElems++;
    blk.tags.push(tag);
  }

  /** Record an unsupported element/record type for the skipped-records report. */
  skip(what, n = 1) { tally(this.skipped, what, n); }

  addNodeSet(name, ids) { appendSet(this.nodeSets, name, ids); }
  addElemSet(name, ids) { appendSet(this.elemSets, name, ids); }

  /**
   * Attach a property record to a tag: `{ name?, section?, material?, thickness?, kind? }`.
   * `section` follows the contract ({ type: 'tube', D, t } or { A, Iy, Iz, J }) in SI units — convert before
   * calling. Members and groups whose elements carry this tag pick it up.
   */
  setProp(tag, rec) { this.props.set(tag, Object.assign(this.props.get(tag) || {}, rec)); }

  /** Element-type histogram of everything added (canonical type names). */
  histogram() {
    const h = {};
    for (const b of this.blocks.values()) if (b.ids.n) h[b.type] = b.ids.n;
    return h;
  }

  /**
   * Build the models.
   * @param {object} [o]
   * @param {number} [o.scale=1] file units → metres
   * @param {string} [o.name]
   * @returns {{ models: object[], warnings: string[], elementTypes: object, skipped: object }}
   */
  build({ scale = 1, name = this.name } = {}) {
    const warnings = [];
    const nNodes = this.nodeIdList.n;
    const xyz = this.coords.a;
    const nodeIds = this.nodeIdList.a;
    const lines = [], shells = [], solids = [];
    for (const b of this.blocks.values()) {
      if (!b.ids.n) continue;
      const d = ELEM[b.type].dim;
      if (d === 1) lines.push(b); else if (d === 2) shells.push(b); else if (d === 3) solids.push(b); else tally(this.skipped, b.type, b.ids.n);
    }
    const kindsPresent = (lines.length ? 1 : 0) + (shells.length ? 1 : 0) + (solids.length ? 1 : 0);
    const label = (s) => (kindsPresent > 1 ? `${name} (${s})` : name);
    let missing = 0;
    const models = [];

    // Resolve a block to model-node indices, compacting the node table to the nodes actually used.
    const compactor = () => {
      const local = new Int32Array(nNodes).fill(-1);
      const order = new Grow(Uint32Array, 1024);
      return {
        local,
        map: (srcId) => {
          const gi = this.nodeIndex.get(srcId);
          if (gi === undefined) return -1;
          let li = local[gi];
          if (li < 0) { li = order.n; local[gi] = li; order.push(gi); }
          return li;
        },
        coords: (Type) => {
          const out = new Type(order.n * 3);
          for (let i = 0; i < order.n; i++) { const g = order.a[i] * 3; out[3 * i] = xyz[g] * scale; out[3 * i + 1] = xyz[g + 1] * scale; out[3 * i + 2] = xyz[g + 2] * scale; }
          return out;
        },
        ids: () => { const out = new Uint32Array(order.n); for (let i = 0; i < order.n; i++) out[i] = nodeIds[order.a[i]]; return out; },
      };
    };
    const nodeSetsFor = (local) => {
      const out = {};
      for (const [sname, ids] of this.nodeSets) {
        const g = new Grow(Uint32Array, ids.length || 4);
        for (const id of ids) { const gi = this.nodeIndex.get(id); if (gi !== undefined && local[gi] >= 0) g.push(local[gi]); }
        if (g.n) out[sname] = g.finish();
      }
      return out;
    };
    const elemSetsFor = (idToLocal) => {
      const out = {};
      for (const [sname, ids] of this.elemSets) {
        const g = new Grow(Uint32Array, ids.length || 4);
        for (const id of ids) { const li = idToLocal.get(id); if (li !== undefined) g.push(li); }
        if (g.n) out[sname] = g.finish();
      }
      return out;
    };
    const tagName = (tag) => { const p = this.props.get(tag); return p && p.name ? String(p.name) : String(tag); };

    if (lines.length) {
      const c = compactor();
      const members = [];
      const idToLocal = new Map();
      for (const b of lines) {
        const cn = b.conn.a;
        for (let e = 0; e < b.ids.n; e++) {
          const a = c.map(cn[e * b.n]), bb = c.map(cn[e * b.n + 1]);
          const mid = b.n === 3 ? c.map(cn[e * b.n + 2]) : -2;
          if (a < 0 || bb < 0 || mid === -1) { missing++; continue; }
          const tag = b.tags.a[e], id = b.ids.a[e];
          const p = this.props.get(tag);
          const base = { id, group: tagName(tag) };
          if (p && p.section) base.section = p.section;
          if (p && p.material !== undefined) base.material = p.material;
          const ex = this.memberExtra.get(id);
          if (ex) Object.assign(base, ex);
          if (!idToLocal.has(id)) idToLocal.set(id, members.length);
          if (mid >= 0) { members.push({ a, b: mid, ...base }); members.push({ a: mid, b: bb, ...base }); } else members.push({ a, b: bb, ...base });
        }
      }
      const m = { kind: 'frame', name: label('beams'), nodes: c.coords(Float64Array), members, nodeIds: c.ids() };
      const sets = { nodes: nodeSetsFor(c.local), members: elemSetsFor(idToLocal) };
      if (Object.keys(sets.nodes).length || Object.keys(sets.members).length) m.sets = sets;
      if (members.length) models.push(m);
    }

    if (shells.length) {
      const c = compactor();
      const buckets = new Map();
      const elements = [];
      const idToLocal = new Map();
      let ecount = 0;
      for (const b of shells) {
        const tri = SURF_TRIS[b.type];
        const cn = b.conn.a;
        const conn = new Grow(Uint32Array, b.conn.n), ids = new Grow(Uint32Array, b.ids.n), tags = new Grow(Int32Array, b.ids.n);
        const loc = new Array(b.n);
        for (let e = 0; e < b.ids.n; e++) {
          let ok = true;
          for (let k = 0; k < b.n; k++) { loc[k] = c.map(cn[e * b.n + k]); if (loc[k] < 0) ok = false; }
          if (!ok) { missing++; continue; }
          const tag = b.tags.a[e];
          let g = buckets.get(tag);
          if (!g) { g = new Grow(Uint32Array, 3072); buckets.set(tag, g); }
          g.ensure(tri.length);
          for (let k = 0; k < tri.length; k++) g.a[g.n++] = loc[tri[k]];
          for (let k = 0; k < b.n; k++) conn.push(loc[k]);
          ids.push(b.ids.a[e]); tags.push(tag);
          idToLocal.set(b.ids.a[e], ecount++);
        }
        if (ids.n) elements.push({ type: b.type, conn: conn.finish(), ids: ids.finish(), tags: tags.finish() });
      }
      let total = 0;
      for (const g of buckets.values()) total += g.n;
      const indices = new Uint32Array(total);
      const groups = [];
      let off = 0;
      for (const [tag, g] of buckets) {
        indices.set(g.view(), off);
        const grp = { name: tagName(tag), start: off, count: g.n };
        const p = this.props.get(tag);
        if (p && p.thickness !== undefined) grp.thickness = p.thickness;
        if (p && p.material !== undefined) grp.material = p.material;
        groups.push(grp);
        off += g.n;
      }
      if (total) {
        const m = surfaceModel(label('shells'), c.coords(Float64Array), indices, groups.length > 1 || (groups[0] && groups[0].name !== '0') ? groups : null, { elements, nodeIds: c.ids() });
        const sets = { nodes: nodeSetsFor(c.local), elements: elemSetsFor(idToLocal) };
        if (Object.keys(sets.nodes).length || Object.keys(sets.elements).length) m.sets = sets;
        models.push(m);
      }
    }

    if (solids.length) {
      const c = compactor();
      const cells = [];
      const idToLocal = new Map();
      let ccount = 0;
      const merged = new Map();
      for (const b of solids) {
        const red = ELEM[b.type].reduce;
        const outType = red || b.type;
        const on = ELEM[outType].n;
        if (red) warnings.push(`${b.ids.n} ${b.type} element(s) reduced to ${red}: the extra mid-side/centre nodes are not kept.`);
        let dst = merged.get(outType);
        if (!dst) { dst = { type: outType, conn: new Grow(Uint32Array, b.ids.n * on), ids: new Grow(Uint32Array, b.ids.n), tags: new Grow(Int32Array, b.ids.n) }; merged.set(outType, dst); }
        const cn = b.conn.a;
        const loc = new Array(on);
        for (let e = 0; e < b.ids.n; e++) {
          let ok = true;
          for (let k = 0; k < on; k++) { loc[k] = c.map(cn[e * b.n + k]); if (loc[k] < 0) ok = false; }
          if (!ok) { missing++; continue; }
          for (let k = 0; k < on; k++) dst.conn.push(loc[k]);
          dst.ids.push(b.ids.a[e]); dst.tags.push(b.tags.a[e]);
        }
      }
      for (const d of merged.values()) {
        if (!d.ids.n) continue;
        const ids = d.ids.finish();
        for (let i = 0; i < ids.length; i++) idToLocal.set(ids[i], ccount++);
        cells.push({ type: d.type, conn: d.conn.finish(), ids, tags: d.tags.finish() });
      }
      if (cells.length) {
        const m = { kind: 'volume', name: label('solids'), positions: c.coords(Float64Array), cells, nodeIds: c.ids() };
        const sets = { nodes: nodeSetsFor(c.local), cells: elemSetsFor(idToLocal) };
        const tagNames = {};
        for (const blk of cells) for (let i = 0; i < blk.tags.length; i++) { const t = blk.tags[i]; if (!(t in tagNames)) tagNames[t] = tagName(t); }
        m.tagNames = tagNames;
        if (Object.keys(sets.nodes).length || Object.keys(sets.cells).length) m.sets = sets;
        models.push(m);
      }
    }

    if (!models.length && nNodes) {
      const out = new Float64Array(nNodes * 3);
      for (let i = 0; i < nNodes * 3; i++) out[i] = xyz[i] * scale;
      models.push(pointsModel(`${name} (nodes)`, out, null, { nodeIds: nodeIds.slice(0, nNodes) }));
      warnings.push('No supported elements were found; the nodes are returned as a point set.');
    }
    if (missing) warnings.push(`${missing} element(s) reference undefined nodes and were skipped.`);
    for (const k in this.skipped) warnings.push(`Skipped ${this.skipped[k]} record(s) of unsupported type ${k}.`);
    return { models, warnings, elementTypes: this.histogram(), skipped: { ...this.skipped } };
  }
}

function appendSet(map, name, ids) {
  const cur = map.get(name);
  if (!cur) map.set(name, Array.from(ids)); else for (const v of ids) cur.push(v);
}

/** Tube section record in SI units (contract §4). */
export function tubeSection(D, t) { return { type: 'tube', D, t }; }

/**
 * Convert a closed chain of 2-D vertices into a Float64Array loop without a duplicated closing vertex.
 * Returns null when fewer than three distinct vertices remain.
 */
export function closeLoop(flat) {
  let n = flat.length;
  if (n >= 4 && Math.abs(flat[0] - flat[n - 2]) < 1e-12 && Math.abs(flat[1] - flat[n - 1]) < 1e-12) n -= 2;
  if (n < 6) return null;
  return Float64Array.from(n === flat.length ? flat : Array.prototype.slice.call(flat, 0, n));
}

// ---------------------------------------------------------------------------------------------------------
// Vertex welding, base64, importer metadata
// ---------------------------------------------------------------------------------------------------------

const F64 = new Float64Array(3);
const U32 = new Uint32Array(F64.buffer);

/**
 * Merge identical vertices while they are being read (triangle soups: STL, panel lists, 3DFACE).
 * `tol` = 0 welds bit-identical coordinates; a positive `tol` snaps to a lattice of that pitch first, which
 * merges points closer than about tol/2 (points straddling a lattice boundary may stay separate).
 */
export class Welder {
  constructor(tol = 0, capacity = 1024) {
    this.tol = tol;
    this.coords = new Grow(Float64Array, capacity * 3);
    this.size = 16;
    while (this.size < capacity * 2) this.size *= 2;
    this.table = new Int32Array(this.size).fill(-1);
    this.count = 0;
  }
  add(x, y, z) {
    const t = this.tol;
    // + 0 turns −0 into +0 so both hash alike
    if (t > 0) { F64[0] = Math.round(x / t) + 0; F64[1] = Math.round(y / t) + 0; F64[2] = Math.round(z / t) + 0; } else { F64[0] = x + 0; F64[1] = y + 0; F64[2] = z + 0; }
    let h = (Math.imul(U32[0] ^ U32[1], 0x9e3779b1) ^ Math.imul(U32[2] ^ U32[3], 0x85ebca6b) ^ Math.imul(U32[4] ^ U32[5], 0xc2b2ae35)) >>> 0;
    const kx = F64[0], ky = F64[1], kz = F64[2];
    const mask = this.size - 1;
    const c = this.coords.a;
    for (let i = h & mask; ; i = (i + 1) & mask) {
      const v = this.table[i];
      if (v < 0) {
        this.table[i] = this.count;
        this.coords.push3(x, y, z);
        if (++this.count * 2 > this.size) this.rehash();
        return this.count - 1;
      }
      const o = 3 * v;
      if (t > 0 ? (Math.round(c[o] / t) + 0 === kx && Math.round(c[o + 1] / t) + 0 === ky && Math.round(c[o + 2] / t) + 0 === kz) : (c[o] === x && c[o + 1] === y && c[o + 2] === z)) return v;
    }
  }
  rehash() {
    const pts = this.coords.a, n = this.count, t = this.tol;
    this.size *= 2;
    this.table = new Int32Array(this.size).fill(-1);
    const mask = this.size - 1;
    for (let v = 0; v < n; v++) {
      if (t > 0) { F64[0] = Math.round(pts[3 * v] / t) + 0; F64[1] = Math.round(pts[3 * v + 1] / t) + 0; F64[2] = Math.round(pts[3 * v + 2] / t) + 0; } else { F64[0] = pts[3 * v] + 0; F64[1] = pts[3 * v + 1] + 0; F64[2] = pts[3 * v + 2] + 0; }
      const h = (Math.imul(U32[0] ^ U32[1], 0x9e3779b1) ^ Math.imul(U32[2] ^ U32[3], 0x85ebca6b) ^ Math.imul(U32[4] ^ U32[5], 0xc2b2ae35)) >>> 0;
      let i = h & mask;
      while (this.table[i] >= 0) i = (i + 1) & mask;
      this.table[i] = v;
    }
  }
  /** Exact-length copy of the welded coordinates. */
  finish() { return this.coords.finish(); }
}

const B64 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
let B64_REV = null;

/** Decode base64 (standard or URL-safe alphabet; whitespace and padding tolerated) to bytes. */
export function base64Decode(s, start = 0, end = s.length) {
  if (!B64_REV) {
    B64_REV = new Int8Array(128).fill(-1);
    for (let i = 0; i < 64; i++) B64_REV[B64.charCodeAt(i)] = i;
    B64_REV[45] = 62; B64_REV[95] = 63;
  }
  const out = new Uint8Array(Math.floor(((end - start) * 3) / 4) + 3);
  let o = 0, acc = 0, bits = 0;
  for (let i = start; i < end; i++) {
    const c = s.charCodeAt(i);
    const v = c < 128 ? B64_REV[c] : -1;
    if (v < 0) continue;
    acc = (acc << 6) | v;
    bits += 6;
    if (bits >= 8) { bits -= 8; out[o++] = (acc >> bits) & 0xff; }
  }
  return out.subarray(0, o);
}

/** Encode bytes as standard base64 with padding. */
export function base64Encode(bytes) {
  const parts = [];
  const CH = 3 * 8192;
  for (let p = 0; p < bytes.length; p += CH) {
    const e = Math.min(bytes.length, p + CH);
    let s = '';
    let i = p;
    for (; i + 2 < e; i += 3) {
      const n = (bytes[i] << 16) | (bytes[i + 1] << 8) | bytes[i + 2];
      s += B64[n >> 18] + B64[(n >> 12) & 63] + B64[(n >> 6) & 63] + B64[n & 63];
    }
    if (i < e) {
      const n = (bytes[i] << 16) | ((i + 1 < e ? bytes[i + 1] : 0) << 8);
      s += B64[n >> 18] + B64[(n >> 12) & 63] + (i + 1 < e ? B64[(n >> 6) & 63] : '=') + '=';
    }
    parts.push(s);
  }
  return parts.join('');
}

/**
 * Advisory warning for files that declare no length unit: says what was assumed and, when the size is
 * implausible for an offshore structure in metres, what the unit probably is.
 * @param {object} units result of resolveUnits()
 * @param {{min:number[],max:number[]}} bbox bounds after scaling (metres)
 * @param {(b: object) => object|null} suggest units.suggestUnits
 */
export function unitAdvice(units, bbox, suggest, warnings) {
  if (units.basis !== 'assumed') return;
  const s = suggest({ min: bbox.min.map((v) => v / units.scale), max: bbox.max.map((v) => v / units.scale) });
  if (s) warnings.push(`The file declares no length unit; ${units.unit} was assumed, but ${s.reason}. Re-import with units = '${s.unit}' if so.`);
  else warnings.push(`The file declares no length unit; ${units.unit} was assumed.`);
}

/** Strip a model name down to something printable (control characters removed, trimmed). */
export function cleanName(s, fallback = 'model') {
  const t = String(s ?? '').replace(/[\u0000-\u001f]+/g, ' ').trim();
  return t || fallback;
}

/**
 * Split triangles that were collected per named bucket into one index array plus contiguous `groups`.
 * @param {Map<string, Grow>} buckets group name → Grow(Uint32Array) of vertex indices
 */
export function groupsFromBuckets(buckets) {
  let total = 0;
  for (const g of buckets.values()) total += g.n;
  const indices = new Uint32Array(total);
  const groups = [];
  let off = 0;
  for (const [name, g] of buckets) {
    if (!g.n) continue;
    indices.set(g.view(), off);
    groups.push({ name, start: off, count: g.n });
    off += g.n;
  }
  return { indices, groups };
}

/** Append a planar polygon (indices into a shared vertex table) to an index list as triangles. */
export function pushPolygon(out, idx, coords) {
  const n = idx.length;
  if (n < 3) return 0;
  if (n === 3) { out.push(idx[0]); out.push(idx[1]); out.push(idx[2]); return 1; }
  if (n === 4) {
    // split along the shorter diagonal so that warped quads stay well shaped
    const d = (a, b) => { const p = 3 * idx[a], q = 3 * idx[b]; return (coords[p] - coords[q]) ** 2 + (coords[p + 1] - coords[q + 1]) ** 2 + (coords[p + 2] - coords[q + 2]) ** 2; };
    if (d(0, 2) <= d(1, 3)) { out.push(idx[0]); out.push(idx[1]); out.push(idx[2]); out.push(idx[0]); out.push(idx[2]); out.push(idx[3]); } else { out.push(idx[0]); out.push(idx[1]); out.push(idx[3]); out.push(idx[1]); out.push(idx[2]); out.push(idx[3]); }
    return 2;
  }
  const flat = new Float64Array(3 * n);
  for (let i = 0; i < n; i++) { flat[3 * i] = coords[3 * idx[i]]; flat[3 * i + 1] = coords[3 * idx[i] + 1]; flat[3 * i + 2] = coords[3 * idx[i] + 2]; }
  let t = triangulate3D(flat);
  if (!t.length) { t = []; for (let i = 1; i + 1 < n; i++) t.push(0, i, i + 1); }
  for (let i = 0; i < t.length; i++) out.push(idx[t[i]]);
  return t.length / 3;
}
