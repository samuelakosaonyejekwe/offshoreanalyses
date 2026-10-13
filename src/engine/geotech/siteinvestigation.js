// From site-investigation records to design soil values: reading a CPT table as exported by common tools
// (delimited text with a header row, GEF, AGS), and reducing a CPT interpretation or a layered borehole
// description to a strength line, a friction angle and a unit weight. The correlations themselves are in cpt.js.
// SI units throughout (m, Pa, N/m³) unless a name says otherwise.

import { cptToLayers } from './cpt.js';

/** Pressure units accepted for cone resistance, sleeve friction and pore pressure, as factors to Pa. */
export const PRESSURE_UNITS = { Pa: 1, kPa: 1e3, MPa: 1e6, bar: 1e5, 'kgf/cm²': 98066.5, tsf: 95760.5, psi: 6894.757 };
/** Length units accepted for depth, as factors to m. */
export const LENGTH_UNITS = { m: 1, cm: 0.01, mm: 0.001, ft: 0.3048 };

const UNIT_ALIASES = [
  [/^mpa$/i, 'MPa'], [/^kpa$|^kn\/m(2|²|\^2)$/i, 'kPa'], [/^pa$/i, 'Pa'], [/^bar$/i, 'bar'], [/^kg(f)?\/cm(2|²|\^2)$/i, 'kgf/cm²'], [/^tsf$|^ton\/ft(2|²)$/i, 'tsf'], [/^psi$/i, 'psi'],
  [/^m$|^met(re|er)s?$/i, 'm'], [/^cm$/i, 'cm'], [/^mm$/i, 'mm'], [/^ft$|^feet$/i, 'ft'],
];
const unitOf = (s) => { const t = String(s || '').trim(); for (const [re, u] of UNIT_ALIASES) if (re.test(t)) return u; return null; };
// column names used by spreadsheets, CPT software and the AGS 3 / AGS 4 data formats
const COLUMN_NAMES = {
  z: /^(depth|z|dept|penetration( length| depth)?|sondeerlengte|scp[tg]_dpth|test depth|d)$/i,
  qt: /^(qt|corrected cone( resistance)?|scpt_qt)$/i,
  qc: /^(qc|q_c|cone|cone resistance|tip|tip resistance|conusweerstand|scp[tg]_res|cone end resistance)$/i,
  fs: /^(fs|f_s|sleeve|sleeve friction|local friction|friction|plaatselijke wrijving|scp[tg]_fres)$/i,
  u2: /^(u2|u_2|u|pore pressure|pwp|pwp2|waterspanning u2|scp[tg]_pwp2|porewater pressure)$/i,
};
const num = (s) => { const t = String(s ?? '').trim().replace(/^"|"$/g, '').replace(/\s/g, ''); if (t === '') return NaN; const v = Number(/^-?\d+,\d+$/.test(t) ? t.replace(',', '.') : t); return Number.isFinite(v) ? v : NaN; };

function splitLine(line, sep) {
  if (sep === ' ') return line.trim().split(/\s+/);
  const out = []; let cur = '', quoted = false;
  for (const ch of line) { if (ch === '"') quoted = !quoted; else if (ch === sep && !quoted) { out.push(cur); cur = ''; } else cur += ch; }
  out.push(cur);
  return out.map((c) => c.trim());
}

/** GEF (Dutch CPT exchange format): columns are identified by quantity number in #COLUMNINFO. */
function parseGef(text) {
  const head = text.slice(0, text.search(/#EOH\s*=?/i)), body = text.slice(text.search(/#EOH\s*=?/i)).replace(/^#EOH\s*=?\s*/i, '');
  const get = (k) => [...head.matchAll(new RegExp(`^#${k}\\s*=\\s*(.*)$`, 'gim'))].map((m) => m[1].trim());
  const sep = (get('COLUMNSEPARATOR')[0] || ' ')[0] || ' ', rec = get('RECORDSEPARATOR')[0]?.[0];
  const voids = Object.fromEntries(get('COLUMNVOID').map((s) => s.split(',').map((x) => x.trim())).map(([c, v]) => [Number(c) - 1, Number(v)]));
  // quantity numbers of the GEF-CPT report: 1 penetration length, 2 cone resistance, 3 local friction, 6 pore pressure u2, 11 corrected depth
  const want = { 1: 'z', 11: 'zc', 2: 'qc', 3: 'fs', 6: 'u2' }, cols = {}, units = {};
  for (const s of get('COLUMNINFO')) { const [c, unit, , q] = s.split(',').map((x) => x.trim()); const key = want[Number(q)]; if (key) { cols[key] = Number(c) - 1; units[key] = unitOf(unit); } }
  if (cols.zc != null) { cols.z = cols.zc; units.z = units.zc; }
  const lines = (rec ? body.split(rec) : body.split(/\r?\n/)).map((l) => l.replace(/[\r\n]/g, ' ').trim()).filter(Boolean);
  const rows = lines.map((l) => splitLine(l, sep === ' ' ? ' ' : sep).map(num)).map((r) => r.map((v, i) => (voids[i] != null && v === voids[i] ? NaN : v)));
  return { cols, units, rows, format: 'GEF' };
}

/** AGS 3 / AGS 4: the static cone group (SCPT), with HEADING, UNIT and DATA rows of quoted values. */
function parseAgs(text) {
  const lines = text.split(/\r?\n/).map((l) => splitLine(l, ','));
  const start = lines.findIndex((l) => /^"?\*?\*?(GROUP)?"?$/i.test(l[0]) ? /SCPT/i.test(l[1] || '') : /^\*\*SCPT/i.test(l[0]));
  const block = lines.slice(start < 0 ? 0 : start + 1);
  const end = block.findIndex((l) => /^\*?\*?GROUP$/i.test(l[0]) || /^\*\*[A-Z]/.test(l[0]));
  const part = end < 0 ? block : block.slice(0, end);
  const heading = part.find((l) => /^HEADING$/i.test(l[0]) || l.some((c) => /SCP[TG]_DPTH/i.test(c)));
  if (!heading) return null;
  const names = heading.map((c) => c.replace(/^\*/, '')), unitRow = part.find((l) => /^(UNIT|<UNITS>)$/i.test(l[0]));
  const cols = {}, units = {};
  names.forEach((n, i) => { for (const [k, re] of Object.entries(COLUMN_NAMES)) if (cols[k] == null && re.test(n)) { cols[k] = i; units[k] = unitRow ? unitOf(unitRow[i]) : null; } });
  const rows = part.filter((l) => /^DATA$/i.test(l[0]) || (l !== heading && l !== unitRow && Number.isFinite(num(l[cols.z])))).map((l) => l.map(num));
  return { cols, units, rows, format: 'AGS' };
}

function parseDelimited(text) {
  const lines = text.split(/\r?\n/).map((l) => l.trim()).filter((l) => l && !/^[#%!]/.test(l));
  if (!lines.length) return null;
  const probe = lines.slice(0, 30).join('\n');
  const sep = ['\t', ';', ','].map((s) => [s, probe.split(s).length]).sort((a, b) => b[1] - a[1]).find(([, n]) => n > lines.slice(0, 30).length)?.[0] || ' ';
  // a decimal comma with comma-separated columns cannot be told apart: semicolon or tab files may use it
  const cells = lines.map((l) => splitLine(l, sep));
  const first = cells.findIndex((r) => r.filter((c) => Number.isFinite(num(c))).length >= 3 && Number.isFinite(num(r[0])));
  const cols = {}, units = {};
  const header = cells.slice(0, Math.max(first, 0));
  for (const row of header) row.forEach((c, i) => {
    const m = /^(.*?)\s*[[(]\s*([^\])]+)\s*[\])]\s*$/.exec(c), name = (m ? m[1] : c).replace(/^"|"$/g, '').trim();
    for (const [k, re] of Object.entries(COLUMN_NAMES)) if (cols[k] == null && re.test(name)) { cols[k] = i; if (m) units[k] = unitOf(m[2]); }
    if (!m && unitOf(c) && Object.values(cols).includes(i)) { const k = Object.keys(cols).find((x) => cols[x] === i); units[k] ||= unitOf(c); }
  });
  // a units row under the names ("m", "MPa", "kPa")
  for (const row of header) if (row.every((c) => !c || unitOf(c))) for (const [k, i] of Object.entries(cols)) units[k] ||= unitOf(row[i]);
  if (cols.z == null && !header.length) { cols.z = 0; cols.qc = 1; cols.fs = 2; if (cells[first]?.length > 3) cols.u2 = 3; }
  return { cols, units, rows: cells.slice(Math.max(first, 0)).map((r) => r.map(num)), format: header.length ? 'table with header' : 'table without header (columns taken as depth, cone resistance, sleeve friction, pore pressure)' };
}

/**
 * Read a CPT table. Recognises GEF files, the SCPT group of AGS files, and delimited text (comma, semicolon, tab
 * or blanks) with a header naming the columns — depth, qc or qt, fs and optionally u2 — or with no header at all,
 * in which case the columns are taken in that order. Units written in the header ("qc [MPa]", a units row, GEF
 * column info) are used; the `units` argument fills in the rest.
 * Rows that are not numbers, or that lack depth, cone resistance or sleeve friction, are skipped and counted;
 * rows are sorted by depth and repeated depths dropped.
 * @param {string} text file or pasted content
 * @param {{ depth?: string, qc?: string, fs?: string, u2?: string }} [units] fallback units (defaults m, MPa, kPa, kPa)
 * @returns {{ z: number[], qc: number[], fs: number[], u2: number[]|null, qcIsQt: boolean, format: string,
 *   units: object, skipped: number, warnings: string[] }} z [m]; qc, fs, u2 [Pa]
 */
export function parseCptTable(text, units = {}) {
  const src = String(text || '');
  const p = /#EOH/i.test(src) && /#COLUMNINFO/i.test(src) ? parseGef(src) : (/SCP[TG]_DPTH/i.test(src) && parseAgs(src)) || parseDelimited(src);
  if (!p) throw new RangeError('text: no data rows found');
  const { cols } = p, warnings = [];
  const qKey = cols.qt != null ? 'qt' : 'qc';
  if (cols.z == null || cols[qKey] == null || cols.fs == null) throw new RangeError(`text: could not find the ${cols.z == null ? 'depth' : cols[qKey] == null ? 'cone resistance (qc or qt)' : 'sleeve friction (fs)'} column — name the columns depth, qc, fs, u2 in the first row`);
  const used = { depth: p.units.z || units.depth || 'm', qc: p.units[qKey] || units.qc || 'MPa', fs: p.units.fs || units.fs || 'kPa', u2: p.units.u2 || units.u2 || 'kPa' };
  for (const [k, table] of [['depth', LENGTH_UNITS], ['qc', PRESSURE_UNITS], ['fs', PRESSURE_UNITS], ['u2', PRESSURE_UNITS]]) if (!table[used[k]]) throw new RangeError(`units: unknown ${k} unit "${used[k]}"`);
  const fz = LENGTH_UNITS[used.depth], fq = PRESSURE_UNITS[used.qc], ff = PRESSURE_UNITS[used.fs], fu = PRESSURE_UNITS[used.u2];
  let rows = [], skipped = 0;
  for (const r of p.rows) {
    const z = r[cols.z] * fz, q = r[cols[qKey]] * fq, f = r[cols.fs] * ff, u = cols.u2 != null ? r[cols.u2] * fu : NaN;
    if (!Number.isFinite(z) || !Number.isFinite(q) || !Number.isFinite(f) || !(z > 0) || !(q > 0) || f < 0) { skipped++; continue; }
    rows.push([z, q, f, u]);
  }
  rows.sort((a, b) => a[0] - b[0]);
  rows = rows.filter((r, i) => i === 0 || r[0] > rows[i - 1][0] + 1e-9);
  if (rows.length < 2) throw new RangeError('text: fewer than two usable readings (depth > 0, cone resistance > 0)');
  const hasU = cols.u2 != null && rows.every((r) => Number.isFinite(r[3]));
  if (cols.u2 != null && !hasU) warnings.push('The pore-pressure column has gaps, so it is not used: cone resistance is taken as measured.');
  if (skipped) warnings.push(`${skipped} row${skipped > 1 ? 's were' : ' was'} skipped (not numbers, depth or cone resistance not above zero, or a missing value).`);
  const qMed = rows.map((r) => r[1]).sort((a, b) => a - b)[Math.floor(rows.length / 2)];
  if (qMed < 2e4 || qMed > 1.2e8) warnings.push(`The median cone resistance reads ${(qMed / 1e6).toPrecision(3)} MPa, outside what soils give (about 0.05 to 100 MPa): check the units chosen.`);
  return { z: rows.map((r) => r[0]), qc: rows.map((r) => r[1]), fs: rows.map((r) => r[2]), u2: hasU ? rows.map((r) => r[3]) : null, qcIsQt: qKey === 'qt', format: p.format, units: used, skipped, warnings };
}

/**
 * Least-squares strength line su = su0 + k z through points below the mudline. A fitted line that would give a
 * negative strength at the mudline is replaced by the best line through the origin (su0 = 0).
 * @param {number[]} z depths [m]
 * @param {number[]} su strengths [Pa]
 * @returns {{ su0: number, suGrad: number, n: number, r2: number, throughOrigin: boolean }} [Pa], [Pa/m]
 */
export function fitStrengthLine(z, su) {
  const n = z.length;
  if (n < 2 || su.length !== n) throw new RangeError('z: at least two points of equal length are needed');
  const zm = z.reduce((s, v) => s + v, 0) / n, sm = su.reduce((s, v) => s + v, 0) / n;
  let sxx = 0, sxy = 0, syy = 0, szz = 0, szs = 0;
  for (let i = 0; i < n; i++) { sxx += (z[i] - zm) ** 2; sxy += (z[i] - zm) * (su[i] - sm); syy += (su[i] - sm) ** 2; szz += z[i] * z[i]; szs += z[i] * su[i]; }
  let k = sxx > 0 ? sxy / sxx : 0, su0 = sm - k * zm, throughOrigin = false;
  if (su0 < 0) { k = szs / szz; su0 = 0; throughOrigin = true; }
  if (k < 0) { k = 0; su0 = sm; throughOrigin = false; }
  let sse = 0; for (let i = 0; i < n; i++) sse += (su[i] - su0 - k * z[i]) ** 2;
  return { su0, suGrad: k, n, r2: syy > 0 ? Math.max(0, 1 - sse / syy) : 1, throughOrigin };
}

const quantile = (a, p) => { const s = a.slice().sort((x, y) => x - y); const h = (s.length - 1) * p, i = Math.floor(h); return s[i] + (h - i) * ((s[i + 1] ?? s[i]) - s[i]); };

/**
 * Reduce a CPT sounding to the design values the analyses share: a strength line for the cohesive soil, a
 * cautious friction angle for the sand, and a submerged unit weight.
 *   su0, suGrad — least-squares line through the undrained strengths su = (qt − σv0)/Nkt of every reading
 *                 classed as clay or silt (Ic > 2.6); the same line with Nkt = 10 and 20 gives the stated range;
 *   phiDeg      — the lower quartile of the friction angles of the readings classed as sand (a cautious
 *                 characteristic value: three quarters of the readings are stronger);
 *   gammaSub    — mean submerged unit weight of the readings of the soil that makes up most of the sounding.
 * Reference: correlations of cpt.js (Robertson 2009; Lunne, Robertson & Powell 1997; Kulhawy & Mayne 1990;
 * Baldi et al. 1986 / Jamiolkowski et al. 2003; Robertson & Cabal 2010). Validity: as cptToLayers; a strength
 * line needs at least three cohesive readings, a friction angle at least three sand readings.
 * @param {{ z: number[], qc: number[], fs: number[], u2?: number[]|null, qcIsQt?: boolean, a?: number,
 *   Nkt?: number, gammaW?: number, phiMethod?: string, drMethod?: string }} o as cptToLayers; qcIsQt true when
 *   the table already holds the corrected resistance qt; Nkt default 15
 * @returns {{ su0: number|null, suGrad: number|null, suRange: object|null, suFit: object|null, phiDeg: number|null,
 *   phiRange: number[]|null, Dr: number|null, gammaSub: number, gammaSubClay: number|null, gammaSubSand: number|null,
 *   share: { clay: number, silt: number, sand: number }, depth: number, layers: object[], points: object }}
 */
export function cptDesignValues({ z, qc, fs, u2 = null, qcIsQt = false, a = 0.8, Nkt = 15, gammaW = 10052, phiMethod = 'kulhawy-mayne', drMethod = 'jamiolkowski' }) {
  const { layers, points: P } = cptToLayers({ z, qc, fs, u2: qcIsQt ? undefined : u2 || undefined, a, Nkt, gammaW, phiMethod, drMethod });
  const idx = (g) => P.group.map((x, i) => (g.includes(x) ? i : -1)).filter((i) => i >= 0);
  const coh = idx(['clay', 'silt']), sand = idx(['sand']);
  // each reading stands for the interval up to its depth
  const dz = z.map((v, i) => v - (i ? z[i - 1] : 0)), total = z[z.length - 1];
  const share = Object.fromEntries(['clay', 'silt', 'sand'].map((g) => [g, idx([g]).reduce((s, i) => s + dz[i], 0) / total]));
  const wmean = (ii, arr) => ii.reduce((s, i) => s + arr[i] * dz[i], 0) / ii.reduce((s, i) => s + dz[i], 0);
  const out = { su0: null, suGrad: null, suRange: null, suFit: null, phiDeg: null, phiRange: null, Dr: null, gammaSubClay: null, gammaSubSand: null, share, depth: total, layers, points: P };
  if (coh.length >= 3) {
    const fit = fitStrengthLine(coh.map((i) => z[i]), coh.map((i) => P.su[i]));
    Object.assign(out, { su0: fit.su0, suGrad: fit.suGrad, suFit: fit, gammaSubClay: wmean(coh, P.gamma) - gammaW,
      suRange: { NktLow: 10, NktHigh: 20, su0: [fit.su0 * Nkt / 20, fit.su0 * Nkt / 10], suGrad: [fit.suGrad * Nkt / 20, fit.suGrad * Nkt / 10] } });
  }
  if (sand.length >= 3) {
    const phi = sand.map((i) => Math.min(Math.max(P.phiDeg[i], 20), 45));
    Object.assign(out, { phiDeg: quantile(phi, 0.25), phiRange: [quantile(phi, 0), quantile(phi, 0.5), quantile(phi, 1)], Dr: Math.min(1, Math.max(0, wmean(sand, P.Dr))), gammaSubSand: wmean(sand, P.gamma) - gammaW });
  }
  out.gammaSub = share.sand > 0.5 && out.gammaSubSand != null ? out.gammaSubSand : out.gammaSubClay ?? out.gammaSubSand ?? wmean(z.map((_, i) => i), P.gamma) - gammaW;
  return out;
}

/**
 * Reduce a layered borehole description to the same design values. Each layer: { top, bottom [m], soil 'clay' |
 * 'sand', suTop, suBottom [Pa] (clay), phiDeg (sand), gammaSub [N/m³] }.
 *   su0, suGrad — least-squares line through the strength profile of the clay layers (linear between suTop and
 *                 suBottom), sampled every quarter metre so that thick layers count for more;
 *   phiDeg      — the lowest friction angle among the sand layers (cautious);
 *   gammaSub    — thickness-weighted mean of the soil that makes up most of the profile.
 * @returns {{ su0: number|null, suGrad: number|null, suFit: object|null, phiDeg: number|null, gammaSub: number|null,
 *   share: { clay: number, sand: number }, depth: number }}
 */
export function boreholeDesignValues(layers) {
  if (!Array.isArray(layers) || !layers.length) throw new RangeError('layers: at least one layer is needed');
  let depth = 0, tClay = 0, tSand = 0, gClay = 0, gSand = 0, phi = Infinity;
  const zs = [], sus = [];
  for (const [i, l] of layers.entries()) {
    const t = l.bottom - l.top;
    if (!(l.top >= 0) || !(t > 0)) throw new RangeError(`layers[${i}]: bottom must lie below top, and top at or below the mudline`);
    if (l.soil !== 'clay' && l.soil !== 'sand') throw new RangeError(`layers[${i}].soil: must be 'clay' or 'sand'`);
    depth = Math.max(depth, l.bottom);
    if (l.soil === 'clay') {
      if (!(l.suTop >= 0) || !(l.suBottom >= 0)) throw new RangeError(`layers[${i}]: a clay layer needs suTop and suBottom >= 0`);
      const w = Math.max(1, Math.round(t / 0.25));
      for (let k = 0; k <= w; k++) { zs.push(l.top + (t * k) / w); sus.push(l.suTop + ((l.suBottom - l.suTop) * k) / w); }
      tClay += t; if (l.gammaSub > 0) gClay += l.gammaSub * t;
    } else {
      if (!(l.phiDeg > 0 && l.phiDeg < 50)) throw new RangeError(`layers[${i}].phiDeg: a sand layer needs a friction angle between 0 and 50 degrees`);
      phi = Math.min(phi, l.phiDeg); tSand += t; if (l.gammaSub > 0) gSand += l.gammaSub * t;
    }
  }
  const fit = zs.length ? fitStrengthLine(zs, sus) : null;
  const gC = tClay && gClay ? gClay / tClay : null, gS = tSand && gSand ? gSand / tSand : null;
  return { su0: fit ? fit.su0 : null, suGrad: fit ? fit.suGrad : null, suFit: fit, phiDeg: Number.isFinite(phi) ? phi : null, gammaSub: tSand > tClay ? gS ?? gC : gC ?? gS, share: { clay: tClay / (tClay + tSand), sand: tSand / (tClay + tSand) }, depth };
}
