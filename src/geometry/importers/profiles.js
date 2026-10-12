// Two-dimensional profile formats: SVG vector drawings (cross-sections) and airfoil coordinate files
// (Selig and Lednicer layouts).
//
// SVG → section: every closed sub-path of <path>, plus <rect>, <circle>, <ellipse> and <polygon>, becomes a
// loop; curves are flattened; `transform` attributes are applied; y is flipped so that the section is y-up.
// Loops are grouped into `layers` by the nearest named <g> (id or inkscape:label). Open paths, <line> and
// <polyline> are returned as polylines. SVG is a drawing format: the physical size comes from the width /
// height attributes when they carry real units (mm, cm, in, pt, pc); otherwise one user unit is taken as one
// CSS pixel (1/96 in) and a warning says so — set `opts.units` to state what one user unit is.
// Not read: <use>, <text>, <image>, clipping, masks, patterns, stroke widths (outlines are centre-lines).
//
// Airfoils → section with one closed loop running from the trailing edge over the upper surface to the
// leading edge and back along the lower surface (counter-clockwise for x towards the trailing edge, y up).
// Coordinates are fractions of chord; percent-chord files are normalised. `opts.chord` (m) scales the loop.

import { inputText, baseName, finish, M4, tally, cleanName, closeLoop, ringArea } from './common.js';
import { normaliseUnit, UNIT_SCALE } from '../units.js';
import { parseXml, attr } from '../xml.js';
import { arcPoints } from './nurbs.js';

const CSS = { mm: 1e-3, cm: 1e-2, in: 0.0254, pt: 0.0254 / 72, pc: 0.0254 / 6, px: 0.0254 / 96, m: 1, q: 2.5e-4 };

function cssLength(v) {
  const m = /^\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*([a-zA-Z%]*)\s*$/.exec(v || '');
  if (!m) return null;
  return { value: Number(m[1]), unit: m[2].toLowerCase() };
}

function parseTransform(t) {
  let M = M4.identity();
  const re = /(matrix|translate|scale|rotate|skewX|skewY)\s*\(([^)]*)\)/g;
  for (let m = re.exec(t || ''); m; m = re.exec(t)) {
    const a = m[2].split(/[\s,]+/).filter(Boolean).map(Number);
    let T = M4.identity();
    switch (m[1]) {
      case 'matrix': T[0] = a[0]; T[1] = a[1]; T[4] = a[2]; T[5] = a[3]; T[12] = a[4]; T[13] = a[5]; break;
      case 'translate': T = M4.translation(a[0] || 0, a[1] || 0, 0); break;
      case 'scale': T = M4.scaling(a[0], a[1] ?? a[0], 1); break;
      case 'rotate': { const R = M4.rotation(0, 0, 1, ((a[0] || 0) * Math.PI) / 180); T = a.length >= 3 ? M4.mul(M4.translation(a[1], a[2], 0), M4.mul(R, M4.translation(-a[1], -a[2], 0))) : R; break; }
      case 'skewX': T[4] = Math.tan(((a[0] || 0) * Math.PI) / 180); break;
      case 'skewY': T[1] = Math.tan(((a[0] || 0) * Math.PI) / 180); break;
      default: break;
    }
    M = M4.mul(M, T);
  }
  return M;
}

/** Flatten an SVG path `d` string into sub-paths: `[{ pts: number[] (flat x, y), closed }]`. */
export function flattenPath(d) {
  const tok = d.match(/[a-zA-Z]|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/g) || [];
  const subs = [];
  let cur = null, i = 0, cmd = '', x = 0, y = 0, sx = 0, sy = 0, px = 0, py = 0, prev = '';
  const num = () => Number(tok[i++]);
  const flag = () => {
    // arc flags may be written without separators ("a5 5 0 01 10 0")
    const t = tok[i];
    if (t.length > 1 && /^[01]/.test(t)) { tok[i] = t.slice(1); return Number(t[0]); }
    i++;
    return Number(t);
  };
  const start = (nx, ny) => { cur = { pts: [nx, ny], closed: false }; subs.push(cur); };
  const line = (nx, ny) => { if (!cur) start(x, y); cur.pts.push(nx, ny); };
  const cubic = (x1, y1, x2, y2, x3, y3) => { const n = 24; for (let k = 1; k <= n; k++) { const t = k / n, u = 1 - t; line(u * u * u * x + 3 * u * u * t * x1 + 3 * u * t * t * x2 + t * t * t * x3, u * u * u * y + 3 * u * u * t * y1 + 3 * u * t * t * y2 + t * t * t * y3); } };
  const quad = (x1, y1, x2, y2) => { const n = 16; for (let k = 1; k <= n; k++) { const t = k / n, u = 1 - t; line(u * u * x + 2 * u * t * x1 + t * t * x2, u * u * y + 2 * u * t * y1 + t * t * y2); } };
  const arc = (rx, ry, phiDeg, fa, fs, x2, y2) => {
    // endpoint → centre parameterisation, SVG 1.1 implementation notes F.6.5
    rx = Math.abs(rx); ry = Math.abs(ry);
    if (!rx || !ry || (x === x2 && y === y2)) { line(x2, y2); return; }
    const phi = (phiDeg * Math.PI) / 180, c = Math.cos(phi), s = Math.sin(phi);
    const dx = (x - x2) / 2, dy = (y - y2) / 2;
    const x1p = c * dx + s * dy, y1p = -s * dx + c * dy;
    const lam = (x1p * x1p) / (rx * rx) + (y1p * y1p) / (ry * ry);
    if (lam > 1) { rx *= Math.sqrt(lam); ry *= Math.sqrt(lam); }
    const den = rx * rx * y1p * y1p + ry * ry * x1p * x1p;
    const k = (fa !== fs ? 1 : -1) * Math.sqrt(Math.max(0, (rx * rx * ry * ry - den) / den));
    const cxp = (k * rx * y1p) / ry, cyp = (-k * ry * x1p) / rx;
    const cx = c * cxp - s * cyp + (x + x2) / 2, cy = s * cxp + c * cyp + (y + y2) / 2;
    const a0 = Math.atan2((y1p - cyp) / ry, (x1p - cxp) / rx);
    let sw = Math.atan2((-y1p - cyp) / ry, (-x1p - cxp) / rx) - a0;
    if (fs && sw < 0) sw += 2 * Math.PI; else if (!fs && sw > 0) sw -= 2 * Math.PI;
    const n = Math.max(2, Math.ceil(Math.abs(sw) / (Math.PI / 36)));
    for (let q = 1; q <= n; q++) { const t = a0 + (sw * q) / n; const ex = rx * Math.cos(t), ey = ry * Math.sin(t); if (q === n) line(x2, y2); else line(cx + c * ex - s * ey, cy + s * ex + c * ey); }
  };
  while (i < tok.length) {
    if (/[a-zA-Z]/.test(tok[i])) cmd = tok[i++];
    else if (cmd === 'M') cmd = 'L'; else if (cmd === 'm') cmd = 'l';
    else if (!cmd) { i++; continue; }
    const rel = cmd === cmd.toLowerCase();
    const C = cmd.toUpperCase();
    const ox = rel ? x : 0, oy = rel ? y : 0;
    switch (C) {
      case 'M': x = ox + num(); y = oy + num(); sx = x; sy = y; start(x, y); cmd = rel ? 'm' : 'M'; break;
      case 'L': { const nx = ox + num(), ny = oy + num(); line(nx, ny); x = nx; y = ny; break; }
      case 'H': { const nx = ox + num(); line(nx, y); x = nx; break; }
      case 'V': { const ny = oy + num(); line(x, ny); y = ny; break; }
      case 'C': { const x1 = ox + num(), y1 = oy + num(), x2 = ox + num(), y2 = oy + num(), x3 = ox + num(), y3 = oy + num(); cubic(x1, y1, x2, y2, x3, y3); px = x2; py = y2; x = x3; y = y3; break; }
      case 'S': { const x1 = /[CS]/i.test(prev) ? 2 * x - px : x, y1 = /[CS]/i.test(prev) ? 2 * y - py : y; const x2 = ox + num(), y2 = oy + num(), x3 = ox + num(), y3 = oy + num(); cubic(x1, y1, x2, y2, x3, y3); px = x2; py = y2; x = x3; y = y3; break; }
      case 'Q': { const x1 = ox + num(), y1 = oy + num(), x2 = ox + num(), y2 = oy + num(); quad(x1, y1, x2, y2); px = x1; py = y1; x = x2; y = y2; break; }
      case 'T': { const x1 = /[QT]/i.test(prev) ? 2 * x - px : x, y1 = /[QT]/i.test(prev) ? 2 * y - py : y; const x2 = ox + num(), y2 = oy + num(); quad(x1, y1, x2, y2); px = x1; py = y1; x = x2; y = y2; break; }
      case 'A': { const rx = num(), ry = num(), ph = num(), fa = flag(), fs = flag(), x2 = ox + num(), y2 = oy + num(); arc(rx, ry, ph, fa, fs, x2, y2); x = x2; y = y2; break; }
      case 'Z': if (cur) cur.closed = true; x = sx; y = sy; cur = null; break;
      default: i++;
    }
    prev = C;
    if (C === 'M') prev = 'M';
  }
  return subs;
}

/** SVG drawing → section. */
export function parseSvg(input, opts = {}) {
  const root = parseXml(inputText(input));
  if (root.name !== 'svg') throw new Error(`${input.name || 'SVG'}: not an SVG document (expected root element <svg>, found <${root.tag}>).`);
  const name = baseName(input.name);
  const warnings = [];
  const hist = {}, skipped = {};
  // physical size of one user unit
  const vb = (attr(root, 'viewBox', '') || '').split(/[\s,]+/).filter(Boolean).map(Number);
  const w = cssLength(attr(root, 'width')), h = cssLength(attr(root, 'height'));
  let unit, basis, vbScale = 1;
  if (vb.length === 4 && w && vb[2] > 0 && CSS[w.unit || 'px']) vbScale = w.value / vb[2];
  if (opts.units || opts.scale) {
    const u = opts.units ? normaliseUnit(opts.units) : null;
    if (opts.units && !u) throw new RangeError(`units: unknown length unit "${opts.units}".`);
    unit = (u ? UNIT_SCALE[u] : opts.scale); basis = 'option'; vbScale = 1;
  } else if (w && w.unit && w.unit !== 'px' && w.unit !== '%' && CSS[w.unit]) { unit = CSS[w.unit] * vbScale; basis = 'file'; } else {
    unit = CSS.px * vbScale; basis = 'assumed';
    warnings.push('The SVG states no physical size (width/height in mm, cm or in), so one user unit was taken as one CSS pixel (1/96 inch). Set units to what one drawing unit represents.');
  }
  const loops = [], loopLayer = [], polylines = [];
  const height = vb.length === 4 ? vb[3] + 2 * vb[1] : h ? h.value : 0;
  const out = (M, pts, closed, layer, kind) => {
    const n = pts.length / 2;
    const c = new Float64Array(n * 2);
    for (let i = 0; i < n; i++) { const q = M4.apply(M, pts[2 * i], pts[2 * i + 1], 0); c[2 * i] = q[0] * unit; c[2 * i + 1] = (height - q[1]) * unit; }
    if (closed) { const l = closeLoop(c); if (l && Math.abs(ringArea(l)) > 0) { loops.push(l); loopLayer.push(layer); tally(hist, kind); return; } tally(skipped, 'degenerate closed path'); return; }
    if (n >= 2) { const p3 = new Float64Array(n * 3); for (let i = 0; i < n; i++) { p3[3 * i] = c[2 * i]; p3[3 * i + 1] = c[2 * i + 1]; } polylines.push({ c: p3, layer }); tally(hist, `${kind} (open)`); }
  };
  const N = (el, k, d = 0) => { const l = cssLength(attr(el, k)); return l ? l.value : d; };
  const walk = (el, M, layer) => {
    const T = attr(el, 'transform') ? M4.mul(M, parseTransform(attr(el, 'transform'))) : M;
    switch (el.name) {
      case 'svg': case 'g': case 'a': case 'switch': {
        const lab = attr(el, 'label') || attr(el, 'id');
        const nl = el.name === 'g' && lab ? cleanName(lab) : layer;
        for (const c of el.children) walk(c, T, nl);
        break;
      }
      case 'path': for (const sp of flattenPath(attr(el, 'd', ''))) out(T, sp.pts, sp.closed, layer, 'path'); break;
      case 'rect': { const x = N(el, 'x'), y = N(el, 'y'), ww = N(el, 'width'), hh = N(el, 'height'); if (ww > 0 && hh > 0) out(T, [x, y, x + ww, y, x + ww, y + hh, x, y + hh], true, layer, 'rect'); break; }
      case 'circle': case 'ellipse': {
        const cx = N(el, 'cx'), cy = N(el, 'cy'), rx = el.name === 'circle' ? N(el, 'r') : N(el, 'rx'), ry = el.name === 'circle' ? rx : N(el, 'ry');
        if (!(rx > 0 && ry > 0)) break;
        const a = arcPoints([cx, cy, 0], [rx, 0, 0], [0, ry, 0], 1, 0, 2 * Math.PI);
        const p = [];
        for (let i = 0; i < a.length - 3; i += 3) p.push(a[i], a[i + 1]);
        out(T, p, true, layer, el.name);
        break;
      }
      case 'polygon': case 'polyline': { const p = (attr(el, 'points', '') || '').split(/[\s,]+/).filter(Boolean).map(Number); if (p.length >= 4) out(T, p.slice(0, p.length - (p.length % 2)), el.name === 'polygon', layer, el.name); break; }
      case 'line': out(T, [N(el, 'x1'), N(el, 'y1'), N(el, 'x2'), N(el, 'y2')], false, layer, 'line'); break;
      case 'defs': case 'symbol': case 'clipPath': case 'mask': case 'pattern': case 'marker': case 'metadata': case 'style': case 'title': case 'desc': case 'namedview': case 'linearGradient': case 'radialGradient': case 'filter': break;
      default: tally(skipped, `<${el.name}>`);
    }
  };
  walk(root, M4.identity(), 'drawing');
  const models = [];
  if (loops.length) {
    const byLayer = new Map();
    loopLayer.forEach((l, i) => { let a = byLayer.get(l); if (!a) byLayer.set(l, (a = [])); a.push(i); });
    models.push({ kind: 'section', name, loops, layers: [...byLayer].map(([nm, idx]) => ({ name: nm, loops: idx })) });
  }
  polylines.forEach((p, i) => models.push({ kind: 'polyline', name: `${p.layer} ${i + 1}`, points: p.c, layer: p.layer }));
  if (!models.length) throw new Error(`${input.name || 'SVG'}: no outlines were found (expected <path>, <rect>, <circle>, <ellipse>, <polygon> or <polyline> elements).`);
  for (const k in skipped) warnings.push(`Skipped ${skipped[k]} × ${k}.`);
  if (polylines.length) warnings.push(`${polylines.length} open path(s) cannot bound a section and are returned as polylines.`);
  return finish(models, warnings, { format: 'svg', units: { unit: 'user unit', scale: unit, basis }, crs: null, elementTypes: hist, skipped, viewBox: vb.length === 4 ? vb : null, yAxis: 'flipped to y-up' });
}

/** Airfoil coordinate file (Selig or Lednicer layout). */
export function parseAirfoil(input, opts = {}) {
  const text = inputText(input);
  const fname = input.name || 'airfoil';
  const lines = text.split(/\r\n|\n|\r/);
  let title = '';
  const blocks = [[]];
  let seenData = false;
  for (const raw of lines) {
    const t = raw.trim();
    if (!t) { if (seenData && blocks[blocks.length - 1].length) blocks.push([]); continue; }
    if (t[0] === '#' || t[0] === '!') continue;
    const p = t.split(/[\s,;]+/).filter(Boolean).map(Number);
    if (p.length >= 2 && p[0] === p[0] && p[1] === p[1]) { blocks[blocks.length - 1].push([p[0], p[1]]); seenData = true; } else if (!seenData && !title) title = t;
  }
  let pts = blocks.flat();
  if (pts.length < 5) throw new Error(`${fname}: not an airfoil coordinate file (expected a name line followed by at least five "x y" coordinate pairs).`);
  // Lednicer: a leading "N_upper N_lower" pair of point counts, then the two surfaces each from the leading edge
  let layout = 'selig';
  const [c0, c1] = pts[0];
  if (c0 > 1.5 && c1 > 1.5 && Number.isInteger(Math.round(c0)) && Math.abs(c0 - Math.round(c0)) < 1e-9 && Math.round(c0) + Math.round(c1) === pts.length - 1) {
    const nu = Math.round(c0);
    const upper = pts.slice(1, 1 + nu), lower = pts.slice(1 + nu);
    pts = upper.slice().reverse().concat(lower.slice(1));
    layout = 'lednicer';
  } else if (blocks.length >= 2 && blocks[0].length >= 2 && blocks[1].length >= 2 && blocks[0][0][0] < blocks[0][blocks[0].length - 1][0] && blocks[1][0][0] < blocks[1][blocks[1].length - 1][0]) {
    // two blocks both running leading edge → trailing edge without the count line
    pts = blocks[0].slice().reverse().concat(blocks[1].slice(1));
    layout = 'lednicer';
  }
  let xmin = Infinity, xmax = -Infinity;
  for (const q of pts) { if (q[0] < xmin) xmin = q[0]; if (q[0] > xmax) xmax = q[0]; }
  const span = xmax - xmin;
  if (!(span > 0)) throw new Error(`${fname}: the coordinates have no chordwise extent.`);
  const warnings = [];
  // files in percent of chord (0 … 100) are normalised; anything else is taken as given
  const norm = span > 50 && span < 150 ? 100 : 1;
  if (norm !== 1) warnings.push('Coordinates are in percent of chord and were divided by 100.');
  if (Math.abs(span / norm - 1) > 0.02) warnings.push(`The chordwise extent is ${(span / norm).toPrecision(4)} rather than 1; coordinates are kept as given (not re-normalised).`);
  const chord = opts.chord > 0 ? opts.chord : 1;
  const flat = [];
  for (const q of pts) flat.push((q[0] / norm) * chord, (q[1] / norm) * chord);
  // orientation: trailing edge → upper surface → leading edge → lower surface is counter-clockwise
  let loop = closeLoop(flat);
  if (!loop) throw new Error(`${fname}: fewer than three distinct points.`);
  if (ringArea(loop) < 0) { const r = new Float64Array(loop.length); for (let i = 0; i < loop.length; i += 2) { r[i] = loop[loop.length - 2 - i]; r[i + 1] = loop[loop.length - 1 - i]; } loop = r; warnings.push('The points ran clockwise (lower surface first); the loop was reversed to the upper-surface-first convention.'); }
  // thickness and camber from the two surfaces, sampled at the upper-surface stations
  let tmax = 0, cmax = 0, xt = 0;
  const n = loop.length / 2;
  let le = 0;
  for (let i = 1; i < n; i++) if (loop[2 * i] < loop[2 * le]) le = i;
  const yAt = (from, to, x) => { for (let i = from; i < to; i++) { const x0 = loop[2 * i], x1 = loop[2 * ((i + 1) % n)]; if ((x0 - x) * (x1 - x) <= 0 && x0 !== x1) return loop[2 * i + 1] + ((x - x0) / (x1 - x0)) * (loop[2 * ((i + 1) % n) + 1] - loop[2 * i + 1]); } return null; };
  for (let i = 1; i < le; i++) { const x = loop[2 * i], yu = loop[2 * i + 1], yl = yAt(le, n, x); if (yl === null) continue; const t = yu - yl; if (t > tmax) { tmax = t; xt = x; } const c = Math.abs((yu + yl) / 2); if (c > cmax) cmax = c; }
  const model = { kind: 'section', name: cleanName(title, baseName(input.name)), loops: [loop], chord };
  return finish([model], warnings, {
    format: 'airfoil', units: { unit: opts.chord > 0 ? 'm' : 'chord', scale: chord, basis: opts.chord > 0 ? 'option' : 'file' }, crs: null, elementTypes: { point: n }, skipped: {},
    airfoil: { name: title || baseName(input.name), layout, points: n, leadingEdgeIndex: le, maxThickness: tmax / chord, maxThicknessAt: (xt - loop[2 * le]) / chord, maxCamber: cmax / chord, trailingEdgeGap: Math.abs(loop[1] - loop[2 * (n - 1) + 1]) / chord },
  });
}
