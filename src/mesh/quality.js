/**
 * Element quality assessment for triangle, quadrilateral, tetrahedral, hexahedral, wedge and pyramid meshes,
 * size statistics, and the checks specific to hydrodynamic panel meshes.
 *
 * Metric definitions follow the Verdict library (Stimpson, Ernst, Knupp, Pébay & Thompson, "The Verdict
 * geometric quality library", Sandia report SAND2007-1751) unless stated:
 *   aspect ratio     tri: L_max·(L0 + L1 + L2) / (4√3·A)          = 1 for an equilateral triangle
 *                    quad, hex, wedge, pyramid: longest / shortest edge  = 1 for a square / cube
 *                    tet: circumradius / (3 · inradius)            = 1 for a regular tetrahedron
 *   min / max angle  interior corner angles of the faces [deg] (tet: dihedral angles, 70.53° when regular)
 *   skewness         equiangle skew max[(θmax − θe) / (180 − θe), (θe − θmin) / θe], θe = 60° (tri) or 90°
 *                    (quad, hex faces); tet: equivolume skew 1 − V / V_regular(same circumradius)
 *   scaled Jacobian  smallest corner Jacobian divided by the product of the corner edge lengths, normalised to
 *                    1 for the ideal element; ≤ 0 means the element is inverted or folded
 * The default limits are the common solver-default warning and error levels (aspect ratio 5 / 20, corner angle
 * 20°–120° / 10°–160° for triangles and 45°–135° / 20°–160° for quadrilaterals, skewness 0.5 / 0.85 / 0.95,
 * scaled Jacobian 0.5 / 0.2 / 0); they are a starting point, not a code requirement, and can be overridden.
 *
 * Panel meshes: a diffraction/radiation panel solution needs several panels per wavelength. With n panels per
 * wavelength (default 7) the shortest wave a mesh resolves is λ = n·h and the highest valid frequency follows
 * from the linear dispersion relation ω² = g·k·tanh(k·d).
 */
import { CELL_INFO, CELL_FACES, CELL_EDGES, cellVolume, bounds, orientConsistently, buildEdges, weldVertices } from './topology.js';

const DEG = 180 / Math.PI;
const SQ3 = Math.sqrt(3);

/** Default acceptance limits per element family: { warn, fail } pairs. */
export const QUALITY_LIMITS = {
  tri: { aspectRatio: [5, 20], minAngle: [20, 10], maxAngle: [120, 160], skewness: [0.5, 0.85], scaledJacobian: [0.5, 0.2] },
  quad: { aspectRatio: [5, 20], minAngle: [45, 20], maxAngle: [135, 160], skewness: [0.5, 0.85], scaledJacobian: [0.5, 0.2] },
  tet: { aspectRatio: [3, 10], minAngle: [20, 10], maxAngle: [130, 165], skewness: [0.8, 0.95], scaledJacobian: [0.3, 0.1] },
  hex: { aspectRatio: [10, 100], minAngle: [45, 20], maxAngle: [135, 160], skewness: [0.5, 0.85], scaledJacobian: [0.5, 0.2] },
};
QUALITY_LIMITS.wedge = QUALITY_LIMITS.hex;
QUALITY_LIMITS.pyr = QUALITY_LIMITS.hex;

const sub = (P, a, b) => [P[3 * b] - P[3 * a], P[3 * b + 1] - P[3 * a + 1], P[3 * b + 2] - P[3 * a + 2]];
const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
const len = (u) => Math.hypot(u[0], u[1], u[2]);
const angle = (u, v) => { const l = len(u) * len(v); return l > 0 ? Math.acos(Math.min(1, Math.max(-1, dot(u, v) / l))) * DEG : 0; };

/**
 * Quality of one triangle.
 * @returns {{ aspectRatio, minAngle, maxAngle, skewness, scaledJacobian, size, area }} angles [deg]; size = longest
 *   edge [m]; area [m²]. A zero-area triangle gives aspectRatio = Infinity, skewness = 1, scaledJacobian = 0.
 */
export function triQuality(P, a, b, c) {
  const e0 = sub(P, a, b), e1 = sub(P, b, c), e2 = sub(P, c, a);
  const l0 = len(e0), l1 = len(e1), l2 = len(e2);
  const A = len(cross(e0, e1)) / 2;
  const lmax = Math.max(l0, l1, l2);
  if (!(A > 0)) return { aspectRatio: Infinity, minAngle: 0, maxAngle: 180, skewness: 1, scaledJacobian: 0, size: lmax, area: 0 };
  const angs = [angle(e0, [-e2[0], -e2[1], -e2[2]]), angle(e1, [-e0[0], -e0[1], -e0[2]]), angle(e2, [-e1[0], -e1[1], -e1[2]])];
  const mn = Math.min(...angs), mx = Math.max(...angs);
  return {
    aspectRatio: lmax * (l0 + l1 + l2) / (4 * SQ3 * A), minAngle: mn, maxAngle: mx,
    skewness: Math.max((mx - 60) / 120, (60 - mn) / 60),
    scaledJacobian: 2 * A / Math.max(l0 * l1, l1 * l2, l2 * l0) * 2 / SQ3,
    size: lmax, area: A,
  };
}

/**
 * Quality of one quadrilateral (corners in order around the face).
 * @returns {{ aspectRatio, minAngle, maxAngle, skewness, scaledJacobian, warpage, size, area }} warpage = angle
 *   between the normals of the two triangles on either diagonal [deg], 0 for a flat panel
 */
export function quadQuality(P, a, b, c, d) {
  const v = [a, b, c, d];
  const e = [0, 1, 2, 3].map((k) => sub(P, v[k], v[(k + 1) % 4]));
  const l = e.map(len);
  const n = cross(sub(P, a, c), sub(P, b, d));
  const nl = len(n);
  const area = nl / 2;
  let mn = 180, mx = 0, sj = Infinity;
  for (let k = 0; k < 4; k++) {
    const u = e[k], w = e[(k + 3) % 4];
    const back = [-w[0], -w[1], -w[2]];
    const ang = angle(u, back);
    if (ang < mn) mn = ang; if (ang > mx) mx = ang;
    const j = nl > 0 && l[k] > 0 && l[(k + 3) % 4] > 0 ? dot(cross(w, u), n) / (nl * l[k] * l[(k + 3) % 4]) : 0;
    if (j < sj) sj = j;
  }
  const n1 = cross(e[0], e[1]), n2 = cross(e[2], e[3]);
  const lmin = Math.min(...l), lmax = Math.max(...l);
  return {
    aspectRatio: lmin > 0 ? lmax / lmin : Infinity, minAngle: mn, maxAngle: mx, skewness: Math.max((mx - 90) / 90, (90 - mn) / 90),
    scaledJacobian: sj, warpage: len(n1) > 0 && len(n2) > 0 ? angle(n1, n2) : 0, size: lmax, area,
  };
}

// edges meeting at each corner, ordered so that their triple product is positive for a valid cell
const CORNER_EDGES = {
  tet4: [[1, 2, 3], [2, 0, 3], [0, 1, 3], [0, 2, 1]],
  hex8: [[1, 3, 4], [2, 0, 5], [3, 1, 6], [0, 2, 7], [7, 5, 0], [4, 6, 1], [5, 7, 2], [6, 4, 3]],
  wedge6: [[1, 2, 3], [2, 0, 4], [0, 1, 5], [5, 4, 0], [3, 5, 1], [4, 3, 2]],
  pyr5: [[1, 3, 4], [2, 0, 4], [3, 1, 4], [0, 2, 4]],
};
// triple product of unit edges at a corner of the ideal element
const CORNER_IDEAL = { tet4: Math.SQRT1_2, hex8: 1, wedge6: SQ3 / 2, pyr5: Math.SQRT1_2 };

/**
 * Quality of one volume cell.
 * @param {ArrayLike<number>} P node coordinates
 * @param {ArrayLike<number>} conn connectivity array
 * @param {number} o offset of the cell in `conn`
 * @param {string} type cell type (quadratic cells are judged on their corner nodes)
 * @returns {{ aspectRatio, minAngle, maxAngle, skewness, scaledJacobian, size, volume }} size = longest edge [m];
 *   volume [m³] (negative when inverted). For pyramids the scaled Jacobian is normalised by the value of a
 *   pyramid with equilateral faces.
 */
export function cellQuality(P, conn, o, type) {
  const base = type === 'tet10' ? 'tet4' : type === 'hex20' ? 'hex8' : type;
  const edges = CELL_EDGES[base];
  if (!edges) throw new RangeError(`type: unsupported cell type '${type}'`);
  let lmin = Infinity, lmax = 0;
  for (const [i, j] of edges) { const l = len(sub(P, conn[o + i], conn[o + j])); if (l < lmin) lmin = l; if (l > lmax) lmax = l; }
  const volume = cellVolume(P, conn, o, base);
  let sj = Infinity;
  for (let k = 0; k < CORNER_EDGES[base].length; k++) {
    const [a, b, c] = CORNER_EDGES[base][k];
    const u = sub(P, conn[o + k], conn[o + a]), v = sub(P, conn[o + k], conn[o + b]), w = sub(P, conn[o + k], conn[o + c]);
    const l = len(u) * len(v) * len(w);
    const j = l > 0 ? dot(u, cross(v, w)) / l / CORNER_IDEAL[base] : 0;
    if (j < sj) sj = j;
  }
  let mn = 180, mx = 0, aspect, skew;
  if (base === 'tet4') {
    // dihedral angles from the outward face normals
    const faces = CELL_FACES.tet4.map((f) => cross(sub(P, conn[o + f[0]], conn[o + f[1]]), sub(P, conn[o + f[0]], conn[o + f[2]])));
    let areaSum = 0;
    for (const n of faces) areaSum += len(n) / 2;
    for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) { const d = 180 - angle(faces[i], faces[j]); if (d < mn) mn = d; if (d > mx) mx = d; }
    const a = sub(P, conn[o], conn[o + 1]), b = sub(P, conn[o], conn[o + 2]), c = sub(P, conn[o], conn[o + 3]);
    const det = 2 * dot(a, cross(b, c));
    const absV = Math.abs(volume);
    if (!(absV > 0) || det === 0) { aspect = Infinity; skew = 1; } else {
      const a2 = dot(a, a), b2 = dot(b, b), c2 = dot(c, c);
      const bc = cross(b, c), ca = cross(c, a), ab = cross(a, b);
      const R = len([a2 * bc[0] + b2 * ca[0] + c2 * ab[0], a2 * bc[1] + b2 * ca[1] + c2 * ab[1], a2 * bc[2] + b2 * ca[2] + c2 * ab[2]]) / Math.abs(det);
      const r = 3 * absV / areaSum;
      aspect = R / (3 * r);
      // regular tetrahedron inscribed in the same sphere: V = 8√3/27 · R³
      skew = Math.min(1, Math.max(0, 1 - absV / (8 * SQ3 / 27 * R * R * R)));
    }
  } else {
    let sk = 0;
    for (const f of CELL_FACES[base]) {
      const m = f.length;
      const te = m === 3 ? 60 : 90;
      for (let k = 0; k < m; k++) {
        const p = conn[o + f[k]], q = conn[o + f[(k + 1) % m]], r = conn[o + f[(k + m - 1) % m]];
        const ang = angle(sub(P, p, q), sub(P, p, r));
        if (ang < mn) mn = ang; if (ang > mx) mx = ang;
        sk = Math.max(sk, (ang - te) / (180 - te), (te - ang) / te);
      }
    }
    aspect = lmin > 0 ? lmax / lmin : Infinity;
    skew = sk;
  }
  return { aspectRatio: aspect, minAngle: mn, maxAngle: mx, skewness: skew, scaledJacobian: volume > 0 ? sj : Math.min(sj, 0), size: lmax, volume };
}

function stats(arr) {
  let n = 0, min = Infinity, max = -Infinity, sum = 0, sum2 = 0;
  for (let i = 0; i < arr.length; i++) { const v = arr[i]; if (!Number.isFinite(v)) continue; n++; sum += v; sum2 += v * v; if (v < min) min = v; if (v > max) max = v; }
  if (!n) return { min: NaN, max: NaN, mean: NaN, std: NaN };
  const mean = sum / n;
  return { min, max, mean, std: Math.sqrt(Math.max(0, sum2 / n - mean * mean)) };
}

function histogramOf(values, lo, hi, bins) {
  const edges = Array.from({ length: bins + 1 }, (_, i) => lo + (hi - lo) * i / bins);
  const counts = new Array(bins).fill(0);
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    if (Number.isNaN(v)) continue;
    counts[Math.min(bins - 1, Math.max(0, Math.floor((v - lo) / (hi - lo) * bins)))]++;
  }
  return { edges, counts };
}

const FAMILY = { tet4: 'tet', tet10: 'tet', hex8: 'hex', hex20: 'hex', wedge6: 'wedge', pyr5: 'pyr' };

/** Worst-first: 0 good, 1 warn, 2 fail for one element against a family's limits. */
function rate(q, lim) {
  let r = 0;
  const hi = (v, [w, f]) => { if (!(v <= f)) r = 2; else if (v > w && r < 1) r = 1; };
  const lo = (v, [w, f]) => { if (!(v >= f)) r = 2; else if (v < w && r < 1) r = 1; };
  hi(q.aspectRatio, lim.aspectRatio); lo(q.minAngle, lim.minAngle); hi(q.maxAngle, lim.maxAngle); hi(q.skewness, lim.skewness); lo(q.scaledJacobian, lim.scaledJacobian);
  return r;
}

/**
 * Wave-resolution limit of a panel mesh.
 * @param {number} panelSize characteristic panel length [m]
 * @param {{ panelsPerWavelength?: number, depth?: number, g?: number }} [opts] panels per wavelength (default 7;
 *   6–10 is the usual range); water depth [m] (default deep water); gravity [m/s²]
 * @returns {{ minWavelength: number, maxOmega: number, maxFrequencyHz: number, minPeriod: number }} shortest
 *   resolved wavelength [m], highest valid angular frequency [rad/s], frequency [Hz] and shortest period [s]
 */
export function panelFrequencyLimit(panelSize, { panelsPerWavelength = 7, depth = Infinity, g = 9.80665 } = {}) {
  if (!(panelSize > 0)) throw new RangeError('panelSize: must be > 0');
  if (!(panelsPerWavelength >= 1)) throw new RangeError('panelsPerWavelength: must be >= 1');
  if (!(depth > 0)) throw new RangeError('depth: must be > 0');
  const lambda = panelsPerWavelength * panelSize;
  const k = 2 * Math.PI / lambda;
  const omega = Math.sqrt(g * k * (Number.isFinite(depth) ? Math.tanh(k * depth) : 1));
  return { minWavelength: lambda, maxOmega: omega, maxFrequencyHz: omega / (2 * Math.PI), minPeriod: 2 * Math.PI / omega };
}

/**
 * Checks for a hydrodynamic panel mesh (surface model): normal direction, closure, panel size against the
 * shortest wave to be resolved, and panels above the waterline.
 * @param {object} model CONTRACT §4 surface model (z up, still water level at `waterline`)
 * @param {{ panelsPerWavelength?: number, depth?: number, g?: number, waterline?: number, minPeriod?: number }} [opts]
 *   see `panelFrequencyLimit`; waterline elevation [m] (default 0); shortest wave period of interest [s] — when
 *   given, the panel size needed for it is returned too
 * @returns {{ panels, panelSizeMax, panelSizeP95, panelSizeMean, minWavelength, maxOmega, maxFrequencyHz,
 *   minPeriod, normalsConsistent, normalsOutward, reversedPanels, closed, panelsAboveWaterline, wettedPanels,
 *   requiredPanelSize?, adequate?, notes: string[] }} sizes [m] (panel size = longest edge); the frequency limit
 *   is based on the 95th-percentile panel size so that a few large panels on flat areas do not govern, and
 *   considers wetted panels only when any exist
 */
export function panelCheck(model, { panelsPerWavelength = 7, depth = Infinity, g = 9.80665, waterline = 0, minPeriod } = {}) {
  if (!model || model.kind !== 'surface') throw new RangeError('model: a surface model is required');
  const P = model.positions, I = model.indices, F = (I.length / 3) | 0;
  if (!F) throw new RangeError('model: the surface has no panels');
  const sizes = [];
  let above = 0;
  for (let f = 0; f < F; f++) {
    const a = I[3 * f], b = I[3 * f + 1], c = I[3 * f + 2];
    const wet = Math.min(P[3 * a + 2], P[3 * b + 2], P[3 * c + 2]) < waterline - 1e-9;
    if (!wet) { above++; continue; }
    sizes.push(Math.max(len(sub(P, a, b)), len(sub(P, b, c)), len(sub(P, c, a))));
  }
  if (!sizes.length) for (let f = 0; f < F; f++) sizes.push(Math.max(len(sub(P, I[3 * f], I[3 * f + 1])), len(sub(P, I[3 * f + 1], I[3 * f + 2])), len(sub(P, I[3 * f + 2], I[3 * f]))));
  sizes.sort((x, y) => x - y);
  const p95 = sizes[Math.min(sizes.length - 1, Math.round(0.95 * (sizes.length - 1)))];
  const w = weldVertices(P, 1e-6 * bounds(P).diagonal);
  const idx = Uint32Array.from(I, (v) => w.map[v]);
  const ori = orientConsistently(w.positions, idx, w.count);
  const ed = buildEdges(idx, w.count);
  let open = 0;
  for (let e = 0; e < ed.edgeCount; e++) if (ed.edgeFaceStart[e + 1] - ed.edgeFaceStart[e] !== 2) open++;
  const lim = panelFrequencyLimit(p95, { panelsPerWavelength, depth, g });
  const notes = [];
  const closed = open === 0;
  if (!ori.orientable) notes.push('Panel normals cannot be made consistent; the surface is not a valid panel model.');
  else if (ori.flipped) notes.push(closed ? `${ori.flipped} panel(s) have normals pointing into the body; panel codes need all normals pointing into the fluid.` : `${ori.flipped} panel(s) are wound against their neighbours.`);
  if (!closed) notes.push('The surface is open. That is expected for a wetted-hull panel model cut at the waterline, but hydrostatics from an open mesh need the waterplane closed.');
  if (above && sizes.length < F) notes.push(`${above} panel(s) lie entirely on or above the waterline and do not belong in a diffraction model.`);
  const out = {
    panels: F, panelSizeMax: sizes[sizes.length - 1], panelSizeP95: p95, panelSizeMean: sizes.reduce((s, v) => s + v, 0) / sizes.length,
    ...lim, normalsConsistent: ori.orientable && (ori.flipped === 0 || ori.flipped === F), normalsOutward: closed && ori.orientable && ori.flipped === 0,
    reversedPanels: ori.flipped, closed, panelsAboveWaterline: above, wettedPanels: F - above, notes,
  };
  if (minPeriod > 0) {
    // wavelength at the shortest period from the dispersion relation (fixed-point on k)
    const om = 2 * Math.PI / minPeriod;
    let k = om * om / g;
    if (Number.isFinite(depth)) for (let i = 0; i < 100; i++) k = om * om / (g * Math.tanh(k * depth));
    out.requiredPanelSize = 2 * Math.PI / k / panelsPerWavelength;
    out.adequate = p95 <= out.requiredPanelSize * (1 + 1e-9);
    notes.push(out.adequate
      ? `The mesh resolves waves down to ${lim.minPeriod.toFixed(2)} s; the shortest period of interest is ${minPeriod} s.`
      : `Panels are too large for ${minPeriod} s waves: ${p95.toPrecision(3)} m against ${out.requiredPanelSize.toPrecision(3)} m needed for ${panelsPerWavelength} panels per wavelength. Results are valid only above ${lim.minPeriod.toFixed(2)} s.`);
  }
  return out;
}

const METRICS = ['aspectRatio', 'minAngle', 'maxAngle', 'skewness', 'scaledJacobian', 'size'];

/**
 * Assess the element quality of a model.
 *
 * surface → triangles; volume → every cell; frame → member lengths; polyline → segment lengths; grid → cell
 * shape from the spacing. Points and sections have no elements and return an empty assessment.
 *
 * @param {object} model CONTRACT §4 model; a surface may carry `quads` (Uint32Array, four ids per face) which
 *   are then assessed as quadrilaterals alongside the triangles
 * @param {{ limits?: object, bins?: number, panel?: boolean | object }} [opts] limit overrides merged into
 *   QUALITY_LIMITS per family ({ tri: { minAngle: [25, 15] } }); histogram bins (default 10); `panel` adds the
 *   panel-mesh check for surfaces (pass options for `panelCheck`, or false to skip; default on)
 * @returns {{ perElement: { aspectRatio, minAngle, maxAngle, skewness, scaledJacobian, size: Float32Array,
 *   rating: Uint8Array }, histogram: { metric: string, edges: number[], counts: number[] },
 *   histograms: object, summary: object, grade: 'A'|'B'|'C'|'D'|'-', diagnosis: string[], panel?: object }}
 *   perElement arrays are in element order (surface faces, then quads; volume cells in running block order);
 *   rating 0 good, 1 marginal, 2 poor. The headline histogram is of the scaled Jacobian (1 = ideal).
 *   summary: { elements, byType, metric ranges { min, max, mean, std }, size: { min, max, mean, std, ratio },
 *   marginal, poor, inverted, worst: { element, metric, value } }.
 *   grade: A all elements good; B no poor elements and at most 5 % marginal; C at most 1 % poor and none
 *   inverted; D otherwise.
 */
export function assess(model, opts = {}) {
  if (!model || typeof model !== 'object') throw new RangeError('model: a geometry model is required');
  const bins = opts.bins ?? 10;
  const limits = {};
  for (const fam of Object.keys(QUALITY_LIMITS)) limits[fam] = { ...QUALITY_LIMITS[fam], ...((opts.limits && opts.limits[fam]) || {}) };
  const rows = [];
  const fams = [];
  const byType = {};
  const add = (q, fam, type) => { rows.push(q); fams.push(fam); byType[type] = (byType[type] || 0) + 1; };
  const diagnosis = [];
  if (model.kind === 'surface') {
    const P = model.positions, I = model.indices;
    for (let f = 0; f < I.length; f += 3) add(triQuality(P, I[f], I[f + 1], I[f + 2]), 'tri', 'tri3');
    if (model.quads) for (let f = 0; f < model.quads.length; f += 4) add(quadQuality(P, model.quads[f], model.quads[f + 1], model.quads[f + 2], model.quads[f + 3]), 'quad', 'quad4');
  } else if (model.kind === 'volume') {
    for (const b of model.cells || []) {
      const info = CELL_INFO[b.type];
      if (!info) throw new RangeError(`cells.type: unsupported cell type '${b.type}'`);
      for (let o = 0; o + info.nodes <= b.conn.length; o += info.nodes) add(cellQuality(model.positions, b.conn, o, b.type), FAMILY[b.type], b.type);
    }
  } else if (model.kind === 'frame' || model.kind === 'polyline') {
    const N = model.kind === 'frame' ? model.nodes : model.points;
    const segs = model.kind === 'frame' ? model.members.map((m) => [m.a, m.b]) : Array.from({ length: Math.max(0, N.length / 3 - 1) }, (_, i) => [i, i + 1]);
    const size = Float32Array.from(segs, ([a, b]) => len(sub(N, a, b)));
    const s = stats(size);
    const ratio = s.min > 0 ? s.max / s.min : Infinity;
    const what = model.kind === 'frame' ? 'member' : 'segment';
    if (!size.length) diagnosis.push(`The ${model.kind} has no ${what}s.`);
    else {
      diagnosis.push(`${size.length} ${what}s, ${s.min.toPrecision(3)}–${s.max.toPrecision(3)} m long (mean ${s.mean.toPrecision(3)} m).`);
      if (!(s.min > 0)) diagnosis.push(`At least one ${what} has zero length; repair the model.`);
      else if (ratio > 20) diagnosis.push(`The longest ${what} is ${ratio.toFixed(0)} times the shortest. For beam elements this spread harms the conditioning of the stiffness matrix; subdivide the long ${what}s or merge the short ones.`);
    }
    const grade = !size.length ? '-' : !(s.min > 0) ? 'D' : ratio > 100 ? 'C' : ratio > 20 ? 'B' : 'A';
    const lo = Number.isFinite(s.min) ? s.min : 0, hi = s.max > lo ? s.max : lo + 1;
    const histogram = { metric: 'size', ...histogramOf(size, lo, hi, bins) };
    return { perElement: { size, rating: new Uint8Array(size.length) }, histogram, histograms: { size: histogram }, summary: { elements: size.length, byType: { [what]: size.length }, size: { ...s, ratio }, marginal: 0, poor: 0, inverted: 0, worst: null }, grade, diagnosis };
  } else if (model.kind === 'grid') {
    const ax = Math.abs(model.dx), ay = Math.abs(model.dy);
    const ar = Math.max(ax, ay) / Math.min(ax, ay);
    diagnosis.push(`Regular grid of ${model.nx} × ${model.ny} nodes at ${ax} × ${ay} m spacing${ar > 1.5 ? `; the cells are stretched ${ar.toFixed(1)} : 1, so slopes are resolved unequally in the two directions` : ''}.`);
    const histogram = { metric: 'size', edges: [Math.min(ax, ay), Math.max(ax, ay)], counts: [Math.max(0, (model.nx - 1) * (model.ny - 1))] };
    return { perElement: {}, histogram, histograms: { size: histogram }, summary: { elements: Math.max(0, (model.nx - 1) * (model.ny - 1)), byType: { cell: Math.max(0, (model.nx - 1) * (model.ny - 1)) }, aspectRatio: { min: ar, max: ar, mean: ar, std: 0 }, size: { min: Math.min(ax, ay), max: Math.max(ax, ay), mean: (ax + ay) / 2, std: Math.abs(ax - ay) / 2, ratio: ar }, marginal: 0, poor: 0, inverted: 0, worst: null }, grade: ar > 5 ? 'C' : ar > 2 ? 'B' : 'A', diagnosis };
  } else if (model.kind === 'points' || model.kind === 'section') {
    return { perElement: {}, histogram: { metric: 'size', edges: [], counts: [] }, histograms: {}, summary: { elements: 0, byType: {}, marginal: 0, poor: 0, inverted: 0, worst: null }, grade: '-', diagnosis: [`A ${model.kind} model has no elements to assess${model.kind === 'section' ? '; triangulate it first' : ''}.`] };
  } else throw new RangeError(`model.kind: '${model.kind}' is not recognised`);

  const n = rows.length;
  const perElement = { rating: new Uint8Array(n) };
  for (const m of METRICS) perElement[m] = new Float32Array(n);
  let marginal = 0, poor = 0, inverted = 0, worst = null;
  const counts = { aspectRatio: 0, minAngle: 0, maxAngle: 0, skewness: 0, scaledJacobian: 0 };
  for (let i = 0; i < n; i++) {
    const q = rows[i], lim = limits[fams[i]];
    for (const m of METRICS) perElement[m][i] = q[m];
    const r = rate(q, lim);
    perElement.rating[i] = r;
    if (r === 1) marginal++; else if (r === 2) poor++;
    if (!(q.scaledJacobian > 0)) inverted++;
    if (!worst || q.scaledJacobian < worst.value) worst = { element: i, metric: 'scaledJacobian', value: q.scaledJacobian };
    if (r === 2) {
      if (!(q.aspectRatio <= lim.aspectRatio[1])) counts.aspectRatio++;
      if (!(q.minAngle >= lim.minAngle[1])) counts.minAngle++;
      if (!(q.maxAngle <= lim.maxAngle[1])) counts.maxAngle++;
      if (!(q.skewness <= lim.skewness[1])) counts.skewness++;
      if (!(q.scaledJacobian >= lim.scaledJacobian[1])) counts.scaledJacobian++;
    }
  }
  const summary = { elements: n, byType, marginal, poor, inverted, worst };
  for (const m of METRICS) summary[m] = stats(perElement[m]);
  summary.size.ratio = summary.size.min > 0 ? summary.size.max / summary.size.min : Infinity;
  const histograms = {
    scaledJacobian: { metric: 'scaledJacobian', ...histogramOf(perElement.scaledJacobian, 0, 1, bins) },
    skewness: { metric: 'skewness', ...histogramOf(perElement.skewness, 0, 1, bins) },
    minAngle: { metric: 'minAngle', ...histogramOf(perElement.minAngle, 0, 90, bins) },
    aspectRatio: { metric: 'aspectRatio', ...histogramOf(perElement.aspectRatio, 1, 11, bins) },
  };
  let grade;
  if (!n) grade = '-';
  else if (!poor && !marginal) grade = 'A';
  else if (!poor && marginal <= 0.05 * n) grade = 'B';
  else if (!inverted && poor <= 0.01 * n) grade = 'C';
  else grade = 'D';
  const pct = (k) => `${(100 * k / n).toPrecision(k / n < 0.1 ? 2 : 3)} %`;
  if (!n) diagnosis.push('The mesh has no elements.');
  else {
    diagnosis.push(`${n} elements: ${n - marginal - poor} good, ${marginal} marginal (${pct(marginal)}), ${poor} poor (${pct(poor)}).`);
    if (inverted) diagnosis.push(`${inverted} element(s) are inverted or have zero size (scaled Jacobian ≤ 0). A solver will reject them or return nonsense; repair or re-mesh before use.`);
    if (counts.minAngle) diagnosis.push(`${counts.minAngle} element(s) have a corner angle below the limit (smallest ${summary.minAngle.min.toFixed(1)}°): slivers stiffen the element artificially and degrade stress accuracy. Re-mesh with a smaller growth rate or collapse the short edges.`);
    if (counts.maxAngle) diagnosis.push(`${counts.maxAngle} element(s) have a corner angle above the limit (largest ${summary.maxAngle.max.toFixed(1)}°): nearly flat elements interpolate gradients badly.`);
    if (counts.aspectRatio) diagnosis.push(`${counts.aspectRatio} element(s) are too stretched (aspect ratio up to ${summary.aspectRatio.max.toPrecision(3)}). Stretching is acceptable only where the solution varies in one direction, as in a boundary layer.`);
    if (counts.skewness) diagnosis.push(`${counts.skewness} element(s) are too skewed (up to ${summary.skewness.max.toFixed(2)}).`);
    if (counts.scaledJacobian && !inverted) diagnosis.push(`${counts.scaledJacobian} element(s) have a low scaled Jacobian (down to ${summary.scaledJacobian.min.toFixed(2)}): strongly distorted corners.`);
    if (Number.isFinite(summary.size.ratio) && summary.size.ratio > 50) diagnosis.push(`Element size varies by a factor of ${summary.size.ratio.toFixed(0)} (${summary.size.min.toPrecision(3)}–${summary.size.max.toPrecision(3)} m). Check that the transition is gradual; abrupt size changes reflect waves and concentrate error.`);
    if (grade === 'A') diagnosis.push('All elements are within the recommended limits.');
  }
  const out = { perElement, histogram: histograms.scaledJacobian, histograms, summary, grade, diagnosis };
  if (model.kind === 'surface' && opts.panel !== false && n) {
    try { out.panel = panelCheck(model, typeof opts.panel === 'object' ? opts.panel : {}); } catch { /* panel data are optional */ }
  }
  return out;
}
