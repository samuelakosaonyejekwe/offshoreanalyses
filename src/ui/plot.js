// Shared SVG plot kit. Renders the plot objects produced by analysis modules (docs/CONTRACT.md §2.2).
// One y-axis per chart, categorical colours in fixed slot order, crosshair + tooltip on every form.
import { h, fmt, download, toCsv, icon } from './dom.js';

const SERIES = ['--s1', '--s2', '--s3', '--s4', '--s5', '--s6', '--s7', '--s8'];
const SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'];
const W = 640, H = 360;

const cssVar = (name) => `var(${name})`;
const color = (i) => cssVar(SERIES[i % SERIES.length]);

function niceTicks(lo, hi, n = 6) {
  if (!(hi > lo)) { const d = Math.abs(lo) || 1; lo -= d * 0.5; hi += d * 0.5; }
  const span = hi - lo;
  const step0 = 10 ** Math.floor(Math.log10(span / n));
  const err = span / n / step0;
  const step = step0 * (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1);
  const a = Math.floor(lo / step) * step, b = Math.ceil(hi / step) * step;
  const t = [];
  for (let v = a; v <= b + step * 1e-9; v += step) t.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return { lo: a, hi: b, ticks: t };
}
function logTicks(lo, hi) {
  const a = Math.floor(Math.log10(lo)), b = Math.ceil(Math.log10(hi));
  const t = [];
  for (let e = a; e <= b; e++) t.push(10 ** e);
  return { lo: 10 ** a, hi: 10 ** b, ticks: t };
}
function extent(arrs) {
  let lo = Infinity, hi = -Infinity;
  for (const a of arrs) for (const v of a) if (Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
  return lo === Infinity ? [0, 1] : [lo, hi];
}
function axis(ext, opts = {}, pad = 0) {
  let [lo, hi] = ext;
  if (opts.log) { lo = Math.max(lo, 1e-300); return { ...logTicks(lo, Math.max(hi, lo * 10)), log: true }; }
  if (opts.min != null) lo = opts.min;
  if (opts.max != null) hi = opts.max;
  if (opts.zero || (lo > 0 && lo < 0.25 * hi && opts.zero !== false)) lo = Math.min(lo, 0);
  const sp = (hi - lo) * pad;
  return niceTicks(lo - (opts.min != null || lo === 0 ? 0 : sp), hi + (opts.max != null ? 0 : sp));
}
const scale = (ax, a, b) => ax.log
  ? (v) => a + (Math.log10(Math.max(v, 1e-300)) - Math.log10(ax.lo)) / (Math.log10(ax.hi) - Math.log10(ax.lo)) * (b - a)
  : (v) => a + (v - ax.lo) / (ax.hi - ax.lo) * (b - a);
const axLabel = (a) => (a?.label ? a.label + (a.unit && a.unit !== '-' ? ` [${a.unit}]` : '') : '');

let tipEl;
function tip(evt, html) {
  if (!tipEl) { tipEl = h('div.tip'); document.body.append(tipEl); }
  if (!html) { tipEl.style.display = 'none'; return; }
  tipEl.innerHTML = html; tipEl.style.display = 'block';
  const r = tipEl.getBoundingClientRect();
  const x = Math.min(window.innerWidth - r.width - 8, Math.max(8, evt.clientX + 14));
  const y = evt.clientY - r.height - 12 < 8 ? evt.clientY + 18 : evt.clientY - r.height - 12;
  tipEl.style.left = `${x}px`; tipEl.style.top = `${y}px`;
}
const esc = (s) => String(s).replace(/[&<>]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));

function frame(p, xa, ya, m) {
  const g = h('g');
  const x = scale(xa, m.l, W - m.r), y = scale(ya, H - m.b, m.t);
  for (const t of ya.ticks) {
    g.append(h('line', { x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), stroke: 'var(--grid)', 'stroke-width': 1 }));
    g.append(h('text', { x: m.l - 8, y: y(t) + 4, 'text-anchor': 'end', 'font-size': 11, fill: 'var(--ink-3)' }, fmt(t, 3)));
  }
  for (const t of xa.ticks) {
    g.append(h('line', { x1: x(t), x2: x(t), y1: H - m.b, y2: H - m.b + 4, stroke: 'var(--line)' }));
    g.append(h('text', { x: x(t), y: H - m.b + 17, 'text-anchor': 'middle', 'font-size': 11, fill: 'var(--ink-3)' }, fmt(t, 3)));
  }
  g.append(h('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: 'var(--line)' }));
  if (ya.lo < 0 && ya.hi > 0 && !ya.log) g.append(h('line', { x1: m.l, x2: W - m.r, y1: y(0), y2: y(0), stroke: 'var(--ink-3)', 'stroke-width': 1 }));
  g.append(h('text', { x: (m.l + W - m.r) / 2, y: H - 6, 'text-anchor': 'middle', 'font-size': 12, fill: 'var(--ink-2)' }, axLabel(p.x)));
  g.append(h('text', { x: 14, y: (m.t + H - m.b) / 2, 'text-anchor': 'middle', 'font-size': 12, fill: 'var(--ink-2)', transform: `rotate(-90 14 ${(m.t + H - m.b) / 2})` }, axLabel(p.y)));
  return { g, x, y };
}

function xyChart(p, svg) {
  const series = (p.series || []).filter((s) => s && s.x && s.y && s.x.length);
  const m = { l: 62, r: 18, t: 14, b: 44 };
  const xa = axis(extent(series.map((s) => s.x)), p.x || {}, 0);
  const ya = axis(extent(series.map((s) => s.y)), p.y || {}, 0.06);
  if (p.y?.reverse) [ya.lo, ya.hi] = [ya.hi, ya.lo];
  const { g, x, y } = frame(p, xa, ya, m);
  svg.append(g);
  series.forEach((s, i) => {
    const c = s.color || color(i);
    const pts = [];
    for (let k = 0; k < s.x.length; k++) if (Number.isFinite(s.x[k]) && Number.isFinite(s.y[k])) pts.push(`${x(s.x[k]).toFixed(1)},${y(s.y[k]).toFixed(1)}`);
    const markers = p.type === 'scatter' || s.markers;
    if (p.type === 'area' && pts.length) {
      const y0 = y(ya.log ? ya.lo : Math.max(Math.min(0, ya.hi), Math.min(ya.lo, ya.hi)));
      svg.append(h('polygon', { points: `${x(s.x[0]).toFixed(1)},${y0} ${pts.join(' ')} ${x(s.x[s.x.length - 1]).toFixed(1)},${y0}`, fill: c, opacity: 0.18 }));
    }
    if (p.type !== 'scatter' || s.line) svg.append(h('polyline', { points: pts.join(' '), fill: 'none', stroke: c, 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round', 'stroke-dasharray': s.dash ? '6 5' : null }));
    if (markers) for (const pt of pts.length > 600 ? pts.filter((_, k) => k % Math.ceil(pts.length / 600) === 0) : pts) {
      const [cx, cy] = pt.split(',');
      svg.append(h('circle', { cx, cy, r: 4, fill: c, stroke: 'var(--surface)', 'stroke-width': 2 }));
    }
    if (series.length > 1 && series.length <= 4 && pts.length && p.type !== 'scatter') {
      const [lx, ly] = pts[pts.length - 1].split(',').map(Number);
      svg.append(h('text', { x: Math.min(lx + 4, W - 4), y: ly - 6, 'font-size': 10.5, fill: 'var(--ink-2)', 'text-anchor': 'end' }, s.name || ''));
    }
  });
  // crosshair + tooltip: nearest x across all series
  const cross = h('line', { y1: m.t, y2: H - m.b, stroke: 'var(--ink-3)', 'stroke-width': 1, 'stroke-dasharray': '3 3', visibility: 'hidden' });
  const hit = h('rect', { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: 'transparent' });
  const onMove = (e) => {
    const r = svg.getBoundingClientRect();
    const px = (e.clientX - r.left) / r.width * W, py = (e.clientY - r.top) / r.height * H;
    let best = null;
    for (let i = 0; i < series.length; i++) {
      const s = series[i];
      for (let k = 0; k < s.x.length; k++) {
        const dx = x(s.x[k]) - px, dy = p.type === 'scatter' ? y(s.y[k]) - py : 0;
        const d = dx * dx + dy * dy;
        if (!best || d < best.d) best = { d, xv: s.x[k], px: x(s.x[k]) };
      }
    }
    if (!best) return;
    cross.setAttribute('x1', best.px); cross.setAttribute('x2', best.px); cross.setAttribute('visibility', 'visible');
    const rows = [];
    series.forEach((s, i) => {
      let kk = 0, dd = Infinity;
      for (let k = 0; k < s.x.length; k++) { const d = Math.abs(x(s.x[k]) - best.px); if (d < dd) { dd = d; kk = k; } }
      if (dd < 12) rows.push(`<span style="display:inline-block;width:8px;height:8px;border-radius:2px;background:${s.color || color(i)};margin-right:6px"></span>${esc(s.name || 'value')}: <b>${fmt(s.y[kk])}</b>`);
    });
    tip(e, `${esc(p.x?.label || 'x')}: <b>${fmt(best.xv)}</b><br>${rows.join('<br>')}`);
  };
  hit.addEventListener('pointermove', onMove);
  hit.addEventListener('pointerleave', () => { cross.setAttribute('visibility', 'hidden'); tip(); });
  svg.append(cross, hit);
  return series.map((s, i) => ({ name: s.name, color: s.color || color(i), dash: s.dash }));
}

function barChart(p, svg) {
  const cats = p.categories || [];
  const series = p.series || [];
  const m = { l: 62, r: 14, t: 14, b: cats.some((c) => String(c).length > 8) ? 78 : 44 };
  const ya = axis(extent(series.map((s) => s.y)), { ...(p.y || {}), zero: true }, 0.05);
  const xa = { lo: 0, hi: cats.length, ticks: [] };
  const { g, y } = frame({ ...p, x: {} }, xa, ya, m);
  svg.append(g);
  const bw = (W - m.l - m.r) / Math.max(1, cats.length);
  const inner = Math.min(46, (bw * 0.78) / Math.max(1, series.length));
  const rot = m.b > 60;
  cats.forEach((c, i) => {
    const cx = m.l + bw * (i + 0.5);
    svg.append(h('text', { x: cx, y: H - m.b + 16, 'text-anchor': rot ? 'end' : 'middle', 'font-size': 11, fill: 'var(--ink-3)', transform: rot ? `rotate(-35 ${cx} ${H - m.b + 16})` : null }, String(c).slice(0, 26)));
    series.forEach((s, j) => {
      const v = s.y[i];
      if (!Number.isFinite(v)) return;
      const x0 = cx - (inner * series.length) / 2 + inner * j + 1;
      const y0 = y(Math.max(0, ya.lo)), y1 = y(v);
      const col = s.colors?.[i] || s.color || color(j);
      const r = h('rect', { x: x0, y: Math.min(y0, y1), width: inner - 2, height: Math.max(1, Math.abs(y1 - y0)), rx: 3, fill: col });
      r.addEventListener('pointermove', (e) => tip(e, `${esc(c)}<br>${esc(s.name || p.y?.label || 'value')}: <b>${fmt(v)}</b>`));
      r.addEventListener('pointerleave', () => tip());
      svg.append(r);
    });
  });
  if (p.x?.label) svg.append(h('text', { x: (m.l + W - m.r) / 2, y: H - 6, 'text-anchor': 'middle', 'font-size': 12, fill: 'var(--ink-2)' }, axLabel(p.x)));
  if (p.limit != null) svg.append(h('line', { x1: m.l, x2: W - m.r, y1: y(p.limit), y2: y(p.limit), stroke: 'var(--critical)', 'stroke-width': 1.5, 'stroke-dasharray': '6 4' }));
  return series.length > 1 ? series.map((s, j) => ({ name: s.name, color: s.color || color(j), box: true })) : [];
}

function histChart(p, svg) {
  const vals = Array.from(p.values || []).filter(Number.isFinite);
  const [lo, hi] = extent([vals]);
  const n = p.bins || Math.max(8, Math.min(40, Math.round(Math.sqrt(vals.length))));
  const wBin = (hi - lo) / n || 1;
  const counts = new Array(n).fill(0);
  for (const v of vals) counts[Math.min(n - 1, Math.floor((v - lo) / wBin))]++;
  const m = { l: 62, r: 14, t: 14, b: 44 };
  const xa = niceTicks(lo, hi), ya = niceTicks(0, Math.max(...counts, 1));
  const { g, x, y } = frame({ x: p.x, y: p.y || { label: 'Count' } }, xa, ya, m);
  svg.append(g);
  counts.forEach((c, i) => {
    const x0 = x(lo + i * wBin), x1 = x(lo + (i + 1) * wBin);
    const r = h('rect', { x: x0 + 1, y: y(c), width: Math.max(1, x1 - x0 - 2), height: y(0) - y(c), rx: 2, fill: color(0) });
    r.addEventListener('pointermove', (e) => tip(e, `${fmt(lo + i * wBin)} – ${fmt(lo + (i + 1) * wBin)}<br>count: <b>${c}</b> (${fmt(100 * c / vals.length, 3)} %)`));
    r.addEventListener('pointerleave', () => tip());
    svg.append(r);
  });
  for (const mk of p.markers || []) {
    svg.append(h('line', { x1: x(mk.value), x2: x(mk.value), y1: m.t, y2: H - m.b, stroke: 'var(--ink)', 'stroke-width': 1.5, 'stroke-dasharray': '5 4' }));
    svg.append(h('text', { x: x(mk.value) + 4, y: m.t + 10, 'font-size': 10.5, fill: 'var(--ink-2)' }, mk.label));
  }
  return [];
}

function heatChart(p, svg) {
  const xs = Array.from(p.xs || []), ys = Array.from(p.ys || []), z = p.z || [];
  const m = { l: 62, r: 74, t: 14, b: 44 };
  const [zlo, zhi] = extent(z);
  const diverging = zlo < 0 && zhi > 0 && p.diverging !== false;
  const zmax = Math.max(Math.abs(zlo), Math.abs(zhi));
  const xa = { lo: xs[0], hi: xs[xs.length - 1], ticks: niceTicks(xs[0], xs[xs.length - 1], 6).ticks.filter((t) => t >= Math.min(xs[0], xs[xs.length - 1]) && t <= Math.max(xs[0], xs[xs.length - 1])) };
  const ya = { lo: ys[0], hi: ys[ys.length - 1], ticks: niceTicks(ys[0], ys[ys.length - 1], 6).ticks.filter((t) => t >= Math.min(ys[0], ys[ys.length - 1]) && t <= Math.max(ys[0], ys[ys.length - 1])) };
  if (xa.lo === xa.hi) xa.hi = xa.lo + 1;
  if (ya.lo === ya.hi) ya.hi = ya.lo + 1;
  const x = scale(xa, m.l, W - m.r), y = scale(ya, H - m.b, m.t);
  const cw = (W - m.l - m.r) / Math.max(1, xs.length), ch = (H - m.t - m.b) / Math.max(1, ys.length);
  const col = (v) => {
    if (!Number.isFinite(v)) return 'var(--surface-2)';
    if (diverging) {
      const t = Math.min(1, Math.abs(v) / (zmax || 1));
      return `color-mix(in srgb, ${v >= 0 ? '#d03b3b' : '#2a78d6'} ${(t * 100).toFixed(0)}%, var(--div-mid))`;
    }
    const t = (v - zlo) / ((zhi - zlo) || 1) * (SEQ.length - 1);
    const i = Math.min(SEQ.length - 2, Math.floor(t));
    return `color-mix(in srgb, ${SEQ[i + 1]} ${((t - i) * 100).toFixed(0)}%, ${SEQ[i]})`;
  };
  const xStep = xs.length > 1 ? (x(xs[xs.length - 1]) - x(xs[0])) / (xs.length - 1) : cw;
  const yStep = ys.length > 1 ? (y(ys[ys.length - 1]) - y(ys[0])) / (ys.length - 1) : -ch;
  const g = h('g');
  for (let j = 0; j < ys.length; j++) for (let i = 0; i < xs.length; i++) {
    const v = z[j]?.[i];
    g.append(h('rect', { x: m.l + (i * (W - m.l - m.r)) / xs.length, y: (yStep < 0 ? H - m.b - (j + 1) * ch : m.t + j * ch), width: cw + 0.6, height: ch + 0.6, fill: col(v), 'data-i': i, 'data-j': j }));
  }
  void xStep;
  g.addEventListener('pointermove', (e) => {
    const t = e.target; if (!t.dataset?.i) return;
    const i = +t.dataset.i, j = +t.dataset.j;
    tip(e, `${esc(p.x?.label || 'x')}: <b>${fmt(xs[i])}</b><br>${esc(p.y?.label || 'y')}: <b>${fmt(ys[j])}</b><br>${esc(p.zLabel || 'value')}: <b>${fmt(z[j][i])}</b>${p.zUnit ? ' ' + esc(p.zUnit) : ''}`);
  });
  g.addEventListener('pointerleave', () => tip());
  svg.append(g);
  if (p.type === 'contour') drawContours(svg, xs, ys, z, zlo, zhi, m, cw, ch, yStep < 0);
  const xt = (i) => m.l + ((i + 0.5) * (W - m.l - m.r)) / xs.length;
  const stepX = Math.ceil(xs.length / 8), stepY = Math.ceil(ys.length / 8);
  for (let i = 0; i < xs.length; i += stepX) svg.append(h('text', { x: xt(i), y: H - m.b + 16, 'text-anchor': 'middle', 'font-size': 11, fill: 'var(--ink-3)' }, fmt(xs[i], 3)));
  for (let j = 0; j < ys.length; j += stepY) svg.append(h('text', { x: m.l - 8, y: (yStep < 0 ? H - m.b - (j + 0.5) * ch : m.t + (j + 0.5) * ch) + 4, 'text-anchor': 'end', 'font-size': 11, fill: 'var(--ink-3)' }, fmt(ys[j], 3)));
  svg.append(h('text', { x: (m.l + W - m.r) / 2, y: H - 6, 'text-anchor': 'middle', 'font-size': 12, fill: 'var(--ink-2)' }, axLabel(p.x)));
  svg.append(h('text', { x: 14, y: (m.t + H - m.b) / 2, 'text-anchor': 'middle', 'font-size': 12, fill: 'var(--ink-2)', transform: `rotate(-90 14 ${(m.t + H - m.b) / 2})` }, axLabel(p.y)));
  // colour bar
  const bh = H - m.t - m.b, n = 24;
  for (let k = 0; k < n; k++) {
    const v = diverging ? -zmax + (2 * zmax * (k + 0.5)) / n : zlo + ((zhi - zlo) * (k + 0.5)) / n;
    svg.append(h('rect', { x: W - m.r + 16, y: H - m.b - ((k + 1) * bh) / n, width: 12, height: bh / n + 0.5, fill: col(v) }));
  }
  const top = diverging ? zmax : zhi, bot = diverging ? -zmax : zlo;
  svg.append(h('text', { x: W - m.r + 32, y: m.t + 9, 'font-size': 10.5, fill: 'var(--ink-3)' }, fmt(top, 3)));
  svg.append(h('text', { x: W - m.r + 32, y: H - m.b, 'font-size': 10.5, fill: 'var(--ink-3)' }, fmt(bot, 3)));
  svg.append(h('text', { x: W - 6, y: (m.t + H - m.b) / 2, 'font-size': 11, fill: 'var(--ink-2)', 'text-anchor': 'middle', transform: `rotate(-90 ${W - 6} ${(m.t + H - m.b) / 2})` }, (p.zLabel || '') + (p.zUnit ? ` [${p.zUnit}]` : '')));
  return [];
}

// Marching-squares iso-lines over the cell-centre lattice.
function drawContours(svg, xs, ys, z, zlo, zhi, m, cw, ch, yUp) {
  const levels = niceTicks(zlo, zhi, 7).ticks.filter((t) => t > zlo && t < zhi);
  const px = (i) => m.l + (i + 0.5) * cw;
  const py = (j) => (yUp ? H - m.b - (j + 0.5) * ch : m.t + (j + 0.5) * ch);
  for (const L of levels) {
    let d = '';
    for (let j = 0; j < ys.length - 1; j++) for (let i = 0; i < xs.length - 1; i++) {
      const a = z[j][i], b = z[j][i + 1], c = z[j + 1][i + 1], e = z[j + 1][i];
      if (![a, b, c, e].every(Number.isFinite)) continue;
      const pts = [];
      const cut = (v1, v2, x1, y1, x2, y2) => { if ((v1 < L) !== (v2 < L)) { const t = (L - v1) / (v2 - v1); pts.push([x1 + t * (x2 - x1), y1 + t * (y2 - y1)]); } };
      cut(a, b, px(i), py(j), px(i + 1), py(j));
      cut(b, c, px(i + 1), py(j), px(i + 1), py(j + 1));
      cut(e, c, px(i), py(j + 1), px(i + 1), py(j + 1));
      cut(a, e, px(i), py(j), px(i), py(j + 1));
      for (let k = 0; k + 1 < pts.length; k += 2) d += `M${pts[k][0].toFixed(1)},${pts[k][1].toFixed(1)}L${pts[k + 1][0].toFixed(1)},${pts[k + 1][1].toFixed(1)}`;
    }
    if (d) svg.append(h('path', { d, fill: 'none', stroke: 'var(--ink)', 'stroke-width': 0.8, opacity: 0.55, 'pointer-events': 'none' }));
  }
}

function polarChart(p, svg) {
  const series = p.series || [];
  const cx = W / 2, cy = H / 2 + 4, R = H / 2 - 30;
  const rmax = niceTicks(0, extent(series.map((s) => s.r))[1], 4);
  const rr = (v) => (v / rmax.hi) * R;
  for (const t of rmax.ticks.slice(1)) {
    svg.append(h('circle', { cx, cy, r: rr(t), fill: 'none', stroke: 'var(--grid)' }));
    svg.append(h('text', { x: cx + 3, y: cy - rr(t) - 2, 'font-size': 10, fill: 'var(--ink-3)' }, fmt(t, 3)));
  }
  // compass convention: 0° at top (north / bow), clockwise positive
  const ang = (deg) => ((deg - 90) * Math.PI) / 180;
  for (let d = 0; d < 360; d += 30) {
    svg.append(h('line', { x1: cx, y1: cy, x2: cx + R * Math.cos(ang(d)), y2: cy + R * Math.sin(ang(d)), stroke: 'var(--grid)' }));
    svg.append(h('text', { x: cx + (R + 14) * Math.cos(ang(d)), y: cy + (R + 14) * Math.sin(ang(d)) + 4, 'text-anchor': 'middle', 'font-size': 10.5, fill: 'var(--ink-3)' }, `${d}°`));
  }
  series.forEach((s, i) => {
    const pts = s.thetaDeg.map((t, k) => `${(cx + rr(s.r[k]) * Math.cos(ang(t))).toFixed(1)},${(cy + rr(s.r[k]) * Math.sin(ang(t))).toFixed(1)}`);
    const closed = Math.abs(((s.thetaDeg[s.thetaDeg.length - 1] - s.thetaDeg[0]) % 360)) < 1e-6 || s.closed;
    svg.append(h(closed || s.fill ? 'polygon' : 'polyline', { points: pts.join(' '), fill: s.fill ? color(i) : 'none', 'fill-opacity': 0.15, stroke: color(i), 'stroke-width': 2, 'stroke-linejoin': 'round' }));
    s.thetaDeg.forEach((t, k) => {
      const c = h('circle', { cx: cx + rr(s.r[k]) * Math.cos(ang(t)), cy: cy + rr(s.r[k]) * Math.sin(ang(t)), r: 7, fill: 'transparent' });
      c.addEventListener('pointermove', (e) => tip(e, `${esc(s.name || '')}<br>${fmt(t)}°: <b>${fmt(s.r[k])}</b>${p.rUnit ? ' ' + esc(p.rUnit) : ''}`));
      c.addEventListener('pointerleave', () => tip());
      svg.append(c);
    });
  });
  return series.map((s, i) => ({ name: s.name, color: color(i) }));
}

function tornadoChart(p, svg) {
  const items = (p.items || []).slice().sort((a, b) => Math.abs(b.high - b.low) - Math.abs(a.high - a.low)).slice(0, 14);
  const m = { l: 170, r: 20, t: 14, b: 40 };
  const xa = niceTicks(...extent([items.map((i) => i.low), items.map((i) => i.high), [p.base]]));
  const x = scale(xa, m.l, W - m.r);
  const bh = (H - m.t - m.b) / Math.max(1, items.length);
  for (const t of xa.ticks) {
    svg.append(h('line', { x1: x(t), x2: x(t), y1: m.t, y2: H - m.b, stroke: 'var(--grid)' }));
    svg.append(h('text', { x: x(t), y: H - m.b + 16, 'text-anchor': 'middle', 'font-size': 11, fill: 'var(--ink-3)' }, fmt(t, 3)));
  }
  items.forEach((it, i) => {
    const yy = m.t + i * bh + 3, hh = Math.max(4, Math.min(22, bh - 6));
    const bar = (v, c, lab) => {
      const r = h('rect', { x: Math.min(x(p.base), x(v)), y: yy, width: Math.max(1, Math.abs(x(v) - x(p.base))), height: hh, rx: 3, fill: c });
      r.addEventListener('pointermove', (e) => tip(e, `${esc(it.name)} ${lab}<br>${esc(p.x?.label || 'result')}: <b>${fmt(v)}</b>`));
      r.addEventListener('pointerleave', () => tip());
      svg.append(r);
    };
    bar(it.low, color(0), 'low'); bar(it.high, color(1), 'high');
    svg.append(h('text', { x: m.l - 8, y: yy + hh / 2 + 4, 'text-anchor': 'end', 'font-size': 11, fill: 'var(--ink-2)' }, String(it.name).slice(0, 26)));
  });
  svg.append(h('line', { x1: x(p.base), x2: x(p.base), y1: m.t, y2: H - m.b, stroke: 'var(--ink)', 'stroke-width': 1.5 }));
  svg.append(h('text', { x: (m.l + W - m.r) / 2, y: H - 6, 'text-anchor': 'middle', 'font-size': 12, fill: 'var(--ink-2)' }, axLabel(p.x)));
  return [{ name: 'Input at low value', color: color(0), box: true }, { name: 'Input at high value', color: color(1), box: true }];
}

function plotToCsv(p) {
  if (p.series && p.categories) return toCsv([{ key: 'c', label: p.x?.label || 'category' }, ...p.series.map((s, i) => ({ key: `s${i}`, label: s.name || `series ${i + 1}` }))], p.categories.map((c, k) => Object.fromEntries([['c', c], ...p.series.map((s, i) => [`s${i}`, s.y[k]])])));
  if (p.series) {
    const rows = [];
    p.series.forEach((s) => { const xs = s.x || s.thetaDeg || []; const ys = s.y || s.r || []; for (let k = 0; k < xs.length; k++) rows.push({ s: s.name || '', x: xs[k], y: ys[k] }); });
    return toCsv([{ key: 's', label: 'series' }, { key: 'x', label: p.x?.label || 'x', unit: p.x?.unit }, { key: 'y', label: p.y?.label || 'y', unit: p.y?.unit }], rows);
  }
  if (p.z) return [['', ...p.xs].join(','), ...p.ys.map((yv, j) => [yv, ...p.z[j]].join(','))].join('\n');
  if (p.values) return Array.from(p.values).join('\n');
  if (p.items) return toCsv([{ key: 'name', label: 'input' }, { key: 'low', label: 'low' }, { key: 'high', label: 'high' }], p.items);
  return '';
}

function svgToPng(svg, name) {
  const clone = svg.cloneNode(true);
  const cs = getComputedStyle(document.documentElement);
  // resolve CSS variables so the exported image is self-contained
  let xml = new XMLSerializer().serializeToString(clone).replace(/var\((--[\w-]+)\)/g, (_, v) => cs.getPropertyValue(v).trim() || '#888');
  xml = xml.replace(/color-mix\([^)]*\)\)?/g, '#3987e5');
  const img = new Image();
  img.onload = () => {
    const c = document.createElement('canvas'); c.width = W * 2; c.height = H * 2;
    const ctx = c.getContext('2d');
    ctx.fillStyle = cs.getPropertyValue('--surface').trim() || '#fff'; ctx.fillRect(0, 0, c.width, c.height);
    ctx.drawImage(img, 0, 0, c.width, c.height);
    c.toBlob((b) => download(`${name}.png`, b));
  };
  img.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(xml)}`;
}

/** Render one plot object into a figure element. `model3d` plots are handled by the caller (viewer). */
export function renderPlot(p) {
  const svg = h('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': p.title || 'chart', 'font-family': 'system-ui, sans-serif' });
  let legend = [];
  try {
    if (p.type === 'bar') legend = barChart(p, svg);
    else if (p.type === 'hist') legend = histChart(p, svg);
    else if (p.type === 'heatmap' || p.type === 'contour') legend = heatChart(p, svg);
    else if (p.type === 'polar') legend = polarChart(p, svg);
    else if (p.type === 'tornado') legend = tornadoChart(p, svg);
    else legend = xyChart(p, svg);
  } catch (err) {
    svg.append(h('text', { x: W / 2, y: H / 2, 'text-anchor': 'middle', fill: 'var(--ink-3)', 'font-size': 13 }, `Chart could not be drawn: ${err.message}`));
  }
  const name = (p.title || 'chart').replace(/[^\w-]+/g, '_').slice(0, 60);
  return h('div.figure',
    h('div.row.spread', h('h3', p.title || ''),
      h('div.tools',
        h('button.btn.small', { type: 'button', title: 'Download the plotted data as CSV', onclick: () => download(`${name}.csv`, plotToCsv(p), 'text/csv') }, icon('download'), 'CSV'),
        h('button.btn.small', { type: 'button', title: 'Download the chart as a PNG image', onclick: () => svgToPng(svg, name) }, 'PNG'))),
    p.subtitle ? h('small', p.subtitle) : null,
    svg,
    legend.length > 1 || (legend.length === 1 && legend[0].box) ? h('div.legend', legend.map((l) => h('span', h(l.box ? 'i.box' : 'i', { style: { background: l.color, opacity: l.dash ? 0.6 : 1 } }), l.name || ''))) : null);
}
