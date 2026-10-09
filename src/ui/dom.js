// Small DOM toolkit shared by every view: element factory, icons, number formatting, toast, download.

const SVG_NS = 'http://www.w3.org/2000/svg';
const SVG_TAGS = new Set(['svg', 'g', 'path', 'line', 'rect', 'circle', 'text', 'polyline', 'polygon', 'defs', 'clipPath', 'title', 'tspan', 'pattern', 'ellipse']);

/** h('div.card#id', { onclick, title }, child, 'text', [more]) */
export function h(spec, attrs, ...kids) {
  const m = /^([a-zA-Z0-9]+)?((?:[.#][\w-]+)*)$/.exec(spec) || [];
  const tag = m[1] || 'div';
  const el = SVG_TAGS.has(tag) ? document.createElementNS(SVG_NS, tag) : document.createElement(tag);
  for (const part of (m[2] || '').match(/[.#][\w-]+/g) || []) {
    if (part[0] === '.') el.classList.add(part.slice(1)); else el.id = part.slice(1);
  }
  if (attrs != null && (typeof attrs !== 'object' || attrs instanceof Node || Array.isArray(attrs))) { kids.unshift(attrs); attrs = null; }
  if (attrs) {
    for (const [k, v] of Object.entries(attrs)) {
      if (v == null || v === false) continue;
      if (k.startsWith('on') && typeof v === 'function') el.addEventListener(k.slice(2), v);
      else if (k === 'class') el.classList.add(...String(v).split(/\s+/).filter(Boolean));
      else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
      else if (k === 'html') el.innerHTML = v;
      else if (k === 'value' || k === 'checked' || k === 'disabled' || k === 'selected' || k === 'open') el[k] = v;
      else el.setAttribute(k, v === true ? '' : v);
    }
  }
  append(el, kids);
  return el;
}

function append(el, kids) {
  for (const k of kids) {
    if (k == null || k === false) continue;
    if (Array.isArray(k)) append(el, k);
    else el.append(k instanceof Node ? k : document.createTextNode(String(k)));
  }
}

export const clear = (el) => { while (el.firstChild) el.removeChild(el.firstChild); return el; };

const ICONS = {
  back: 'M15 5l-7 7 7 7', fwd: 'M9 5l7 7-7 7', menu: 'M4 7h16M4 12h16M4 17h16', home: 'M4 11l8-7 8 7v8a1 1 0 0 1-1 1h-4v-6h-6v6H5a1 1 0 0 1-1-1z',
  globe: 'M12 3a9 9 0 1 0 0 18 9 9 0 0 0 0-18zM3 12h18M12 3c3 3 3 15 0 18M12 3c-3 3-3 15 0 18', upload: 'M12 16V4m0 0l-4 4m4-4l4 4M5 20h14',
  cube: 'M12 3l8 4.5v9L12 21l-8-4.5v-9zM12 12l8-4.5M12 12v9M12 12L4 7.5', layers: 'M12 4l9 5-9 5-9-5zM3 14l9 5 9-5', flask: 'M9 3h6M10 3v6l-5 9a2 2 0 0 0 2 3h10a2 2 0 0 0 2-3l-5-9V3',
  check: 'M5 12.5l4.5 4.5L19 7.5', bulb: 'M9 18h6M10 21h4M12 3a6 6 0 0 0-4 10.5c.7.7 1 1.5 1 2.5h6c0-1 .3-1.8 1-2.5A6 6 0 0 0 12 3z',
  link: 'M10 14a4 4 0 0 0 5.7 0l3-3a4 4 0 0 0-5.7-5.7l-1 1M14 10a4 4 0 0 0-5.7 0l-3 3a4 4 0 0 0 5.7 5.7l1-1', install: 'M12 4v11m0 0l-4-4m4 4l4-4M5 20h14',
  sun: 'M12 8a4 4 0 1 0 0 8 4 4 0 0 0 0-8zM12 2v2M12 20v2M2 12h2M20 12h2M5 5l1.5 1.5M17.5 17.5L19 19M5 19l1.5-1.5M17.5 6.5L19 5', play: 'M7 5l12 7-12 7z',
  info: 'M12 3a9 9 0 1 0 0 18 9 9 0 0 0 0-18zM12 11v6M12 7.5v.5', warn: 'M12 4l9 16H3zM12 10v5M12 17.5v.5', grid: 'M4 4h7v7H4zM13 4h7v7h-7zM4 13h7v7H4zM13 13h7v7h-7z',
  download: 'M12 4v11m0 0l-4-4m4 4l4-4M5 20h14', search: 'M11 4a7 7 0 1 0 0 14 7 7 0 0 0 0-14zM20 20l-4-4', coin: 'M12 3a9 9 0 1 0 0 18 9 9 0 0 0 0-18zM14.5 9.5c-.5-1-1.5-1.5-2.5-1.5-1.5 0-2.5.8-2.5 2s1 1.7 2.5 2 2.5.8 2.5 2-1 2-2.5 2c-1 0-2-.5-2.5-1.5M12 6.5V8M12 16v1.5',
  help: 'M12 3a9 9 0 1 0 0 18 9 9 0 0 0 0-18zM9.5 9.5a2.5 2.5 0 1 1 3.5 2.3c-.6.3-1 .8-1 1.7M12 16.5v.5', book: 'M5 4h10a3 3 0 0 1 3 3v13H8a3 3 0 0 1-3-3zM5 17a3 3 0 0 1 3-3h10',
  leaf: 'M5 19c0-9 5-14 15-14 0 10-5 15-14 15zM5 19l8-8', trash: 'M5 7h14M10 7V4h4v3M7 7l1 13h8l1-13', refresh: 'M4 12a8 8 0 0 1 14-5l2 2M20 5v4h-4M20 12a8 8 0 0 1-14 5l-2-2M4 19v-4h4',
  chart: 'M4 20V4M4 20h16M8 16v-5M12 16V8M16 16v-8', settings: 'M12 9a3 3 0 1 0 0 6 3 3 0 0 0 0-6zM4 12h2M18 12h2M12 4v2M12 18v2M6.3 6.3l1.4 1.4M16.3 16.3l1.4 1.4M6.3 17.7l1.4-1.4M16.3 7.7l1.4-1.4',
};
export function icon(name) {
  return h('svg', { viewBox: '0 0 24 24', fill: 'none', stroke: 'currentColor', 'stroke-width': 2, 'stroke-linecap': 'round', 'stroke-linejoin': 'round', 'aria-hidden': 'true' },
    h('path', { d: ICONS[name] || ICONS.info }));
}

/** Engineering number format: 4 significant figures, thin-space thousands, exponent beyond 1e-4..1e7. */
export function fmt(v, sig = 4) {
  if (v == null || v === '') return '–';
  if (typeof v !== 'number') return String(v);
  if (!Number.isFinite(v)) return Number.isNaN(v) ? 'NaN' : (v > 0 ? '∞' : '−∞');
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a < 1e-4 || a >= 1e7) {
    const [m, e] = v.toExponential(Math.max(0, sig - 1)).split('e');
    return `${Number(m)}×10${sup(Number(e))}`;
  }
  const digits = Math.max(0, sig - 1 - Math.floor(Math.log10(a)));
  const s = Number(v.toFixed(Math.min(digits, 12))).toString();
  const [i, f] = s.split('.');
  const ii = i.replace(/\B(?=(\d{3})+(?!\d))/g, ' ');
  return (f ? `${ii}.${f}` : ii).replace('-', '−');
}
const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const sup = (n) => String(n).split('').map((c) => SUP[c] || c).join('');

let toastTimer;
export function toast(msg, ms = 3200) {
  document.querySelector('.toast')?.remove();
  const t = h('div.toast', { role: 'status' }, msg);
  document.body.append(t);
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => t.remove(), ms);
}

export function download(name, content, type = 'text/plain') {
  const blob = content instanceof Blob ? content : new Blob([content], { type });
  const url = URL.createObjectURL(blob);
  const a = h('a', { href: url, download: name });
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

export function toCsv(columns, rows) {
  const esc = (v) => { const s = v == null ? '' : String(v); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
  return [columns.map((c) => esc(c.unit ? `${c.label} [${c.unit}]` : c.label)).join(','),
    ...rows.map((r) => columns.map((c) => esc(r[c.key])).join(','))].join('\n');
}

export function tabs(items, onChange, initial = 0) {
  const bar = h('div.tabs', { role: 'tablist' });
  const body = h('div');
  const show = (i) => {
    [...bar.children].forEach((b, j) => b.setAttribute('aria-selected', String(i === j)));
    clear(body);
    const out = items[i].render();
    if (out instanceof Promise) out.then((n) => { if (bar.children[i].getAttribute('aria-selected') === 'true') { clear(body); body.append(n); } });
    else body.append(out);
    onChange?.(i);
  };
  items.forEach((it, i) => bar.append(h('button', { role: 'tab', type: 'button', onclick: () => show(i) }, it.label)));
  show(Math.min(initial, items.length - 1));
  return h('div', bar, body);
}

export function table(columns, rows, { max = 400 } = {}) {
  const shown = rows.slice(0, max);
  return h('div.tablewrap', { style: { maxHeight: '420px', overflowY: 'auto' } }, h('table.t',
    h('thead', h('tr', columns.map((c) => h('th', { class: c.left ? 'l' : '' }, c.label, c.unit ? h('small', ` [${c.unit}]`) : null)))),
    h('tbody', shown.map((r) => h('tr', columns.map((c) => {
      const v = r[c.key];
      return h('td', { class: c.left || typeof v === 'string' ? 'l' : '' }, v instanceof Node ? v : fmt(v, c.sig));
    })))),
    rows.length > max ? h('tfoot', h('tr', h('td', { colspan: columns.length, class: 'l muted' }, `Showing ${max} of ${rows.length} rows — export CSV for all.`))) : null));
}

export const debounce = (fn, ms = 250) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };
