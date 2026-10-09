// Container, parsers and solvers for imported frequency-dependent hydrodynamic coefficients
// (radiation/diffraction results): added mass A(omega), damping B(omega), excitation F(omega, beta).
// Conventions: SI, exp(+i omega t), unit wave amplitude, heading = propagation direction from +x (deg),
// phases are leads relative to the incident crest at the origin — the conventions of WAMIT output.
//
// JSON schema (parseHydroJson):
//   { "rho": 1025, "g": 9.80665, "depth": null | m,
//     "omega": [nw],                    rad/s
//     "headingsDeg": [nb],
//     "A": [nw][6][6], "B": [nw][6][6], kg, kg·m, kg·m² / N·s/m, N·s, N·m·s
//     "Fre": [nb][nw][6], "Fim": [nb][nw][6],   N/m (i <= 3), N·m/m (i >= 4)
//     "Ainf": [6][6] (optional), "A0": [6][6] (optional) }
// CSV schema (parseHydroCsv), one coefficient per row, header required:
//   kind,omega,headingDeg,i,j,re,im
//   kind = A | B (value in `re`, i and j 1-based) | F (complex, i 1-based, j empty) | AINF | A0 (omega empty)

import { G, dispersion, csolve, needPos } from './_wavehelp.js';

const mat6 = () => Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
const copy6 = (m) => m.map((r) => Array.from(r));

/**
 * Builds a validated hydrodynamic database (frequencies sorted ascending, headings ascending).
 * Reference: data container; conventions of the WAMIT v7 User Manual (numeric output files).
 * @param {object} o
 * @param {ArrayLike<number>} o.omega rad/s
 * @param {number[][][]} o.A added mass per frequency, 6×6, kg / kg·m / kg·m²
 * @param {number[][][]} o.B radiation damping per frequency, 6×6, N·s/m …
 * @param {number[]} [o.headingsDeg] deg
 * @param {number[][][]} [o.Fre] excitation, real part [heading][frequency][6], N/m and N·m/m
 * @param {number[][][]} [o.Fim] imaginary part
 * @param {number[][]} [o.Ainf] infinite-frequency added mass
 * @param {number[][]} [o.A0] zero-frequency added mass
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @param {number|null} [o.depth=null] water depth, m (null = deep)
 * @returns {{omega:number[], headingsDeg:number[], A:number[][][], B:number[][][], Fre:number[][][], Fim:number[][][],
 *   Ainf:number[][]|null, A0:number[][]|null, rho:number, g:number, depth:number|null, periodicHeading:boolean}}
 */
export function createHydroDb({ omega, A, B, headingsDeg = [], Fre = [], Fim = [], Ainf = null, A0 = null, rho = 1025, g = G, depth = null }) {
  if (!omega || omega.length === 0) throw new RangeError('omega must be a non-empty array');
  if (!A || !B || A.length !== omega.length || B.length !== omega.length) throw new RangeError('A and B must have one 6×6 matrix per frequency');
  if (Fre.length !== headingsDeg.length || Fim.length !== headingsDeg.length) throw new RangeError('Fre and Fim must have one entry per heading');
  const order = Array.from(omega, (w, i) => i).sort((p, q) => omega[p] - omega[q]);
  const horder = headingsDeg.map((h, i) => i).sort((p, q) => headingsDeg[p] - headingsDeg[q]);
  const om = order.map((i) => {
    if (!(omega[i] > 0)) throw new RangeError('omega must be > 0');
    return omega[i];
  });
  for (let i = 1; i < om.length; i++) if (om[i] === om[i - 1]) throw new RangeError('omega contains duplicates');
  const hd = horder.map((i) => headingsDeg[i]);
  const chk = (m) => { if (!m || m.length !== 6 || m[0].length !== 6) throw new RangeError('coefficient matrices must be 6×6'); return copy6(m); };
  const db = {
    omega: om, headingsDeg: hd, A: order.map((i) => chk(A[i])), B: order.map((i) => chk(B[i])),
    Fre: horder.map((h) => order.map((i) => Array.from(Fre[h][i]))), Fim: horder.map((h) => order.map((i) => Array.from(Fim[h][i]))),
    Ainf: Ainf ? chk(Ainf) : null, A0: A0 ? chk(A0) : null, rho, g, depth, periodicHeading: false,
  };
  if (hd.length >= 3) {
    let maxGap = 0;
    for (let i = 1; i < hd.length; i++) maxGap = Math.max(maxGap, hd[i] - hd[i - 1]);
    db.periodicHeading = hd[0] + 360 - hd[hd.length - 1] <= 1.5 * maxGap;
  }
  return db;
}

function numericRows(text) {
  const rows = [];
  for (const line of String(text).split(/\r?\n/)) {
    const tok = line.trim().split(/[\s,]+/).filter(Boolean);
    if (tok.length === 0) continue;
    const nums = tok.map(Number);
    if (nums.some((v) => !Number.isFinite(v))) continue; // header or comment line
    rows.push(nums);
  }
  return rows;
}

const perToOmega = (per, periodIs) => (periodIs === 'omega' ? per : (2 * Math.PI) / per);

/**
 * Parses a WAMIT `.1` file (added mass and damping): lines `PER I J Abar [Bbar]`, dimensionalised as
 * A_ij = rho L^k Abar_ij, B_ij = rho L^k omega Bbar_ij with k = 3 (i, j <= 3), 4 (one of them > 3), 5 (both > 3).
 * PER = 0 denotes zero period (infinite frequency), PER < 0 infinite period (zero frequency); those lines carry
 * only the added mass.
 * Reference: WAMIT v7 User Manual (numeric output files; definitions of the non-dimensional output).
 * @param {string} text file content
 * @param {{rho?:number, L?:number, periodIs?:string}} [opts] density kg/m³, WAMIT length scale ULEN m (default 1),
 *   'period' (default, first column in s) or 'omega' (first column in rad/s)
 * @returns {{omega:number[], A:number[][][], B:number[][][], Ainf:number[][]|null, A0:number[][]|null}} SI
 */
export function parseWamit1(text, { rho = 1025, L = 1, periodIs = 'period' } = {}) {
  needPos(rho, 'rho'); needPos(L, 'L');
  const groups = new Map();
  let Ainf = null; let A0 = null;
  for (const r of numericRows(text)) {
    if (r.length < 4) continue;
    const [per, i, j, a, b] = r;
    if (!(i >= 1 && i <= 6 && j >= 1 && j <= 6)) throw new RangeError('WAMIT .1: mode indices must be 1..6');
    const k = 3 + (i > 3 ? 1 : 0) + (j > 3 ? 1 : 0);
    const scale = rho * L ** k;
    if (per <= 0) {
      if (per === 0) { Ainf = Ainf || mat6(); Ainf[i - 1][j - 1] = scale * a; } else { A0 = A0 || mat6(); A0[i - 1][j - 1] = scale * a; }
      continue;
    }
    const w = perToOmega(per, periodIs);
    if (!groups.has(per)) groups.set(per, { w, A: mat6(), B: mat6() });
    const gr = groups.get(per);
    gr.A[i - 1][j - 1] = scale * a;
    gr.B[i - 1][j - 1] = scale * w * (b || 0);
  }
  if (groups.size === 0) throw new RangeError('WAMIT .1: no finite-period records found');
  const list = Array.from(groups.values()).sort((p, q) => p.w - q.w);
  return { omega: list.map((x) => x.w), A: list.map((x) => x.A), B: list.map((x) => x.B), Ainf, A0 };
}

function parseWamitComplex(text, scaleOf, periodIs) {
  const per = new Map(); const heads = new Set();
  for (const r of numericRows(text)) {
    if (r.length < 7) continue;
    const [p, beta, i, , , re, im] = r;
    if (!(i >= 1 && i <= 6)) throw new RangeError('WAMIT: mode index must be 1..6');
    if (!(p > 0)) continue;
    heads.add(beta);
    if (!per.has(p)) per.set(p, new Map());
    const m = per.get(p);
    if (!m.has(beta)) m.set(beta, { re: [0, 0, 0, 0, 0, 0], im: [0, 0, 0, 0, 0, 0] });
    const s = scaleOf(i);
    m.get(beta).re[i - 1] = s * re; m.get(beta).im[i - 1] = s * im;
  }
  if (per.size === 0) throw new RangeError('WAMIT: no records found');
  const periods = Array.from(per.keys()).sort((a, b) => perToOmega(a, periodIs) - perToOmega(b, periodIs));
  const headingsDeg = Array.from(heads).sort((a, b) => a - b);
  const zero = () => ({ re: [0, 0, 0, 0, 0, 0], im: [0, 0, 0, 0, 0, 0] });
  return {
    omega: periods.map((p) => perToOmega(p, periodIs)), headingsDeg,
    re: headingsDeg.map((h) => periods.map((p) => (per.get(p).get(h) || zero()).re)),
    im: headingsDeg.map((h) => periods.map((p) => (per.get(p).get(h) || zero()).im)),
  };
}

/**
 * Parses a WAMIT `.3` (diffraction) or `.2` (Haskind) exciting-force file: lines
 * `PER BETA I Mod Pha Re Im`, dimensionalised as X_i = rho g L^m Xbar_i per unit wave amplitude, m = 2 for
 * forces (i <= 3) and 3 for moments.
 * Reference: WAMIT v7 User Manual (numeric output files; definitions of the non-dimensional output).
 * @param {string} text
 * @param {{rho?:number, g?:number, L?:number, periodIs?:string}} [opts]
 * @returns {{omega:number[], headingsDeg:number[], Fre:number[][][], Fim:number[][][]}} [heading][frequency][6], N/m, N·m/m
 */
export function parseWamit3(text, { rho = 1025, g = G, L = 1, periodIs = 'period' } = {}) {
  needPos(rho, 'rho'); needPos(L, 'L');
  const r = parseWamitComplex(text, (i) => rho * g * L ** (i > 3 ? 3 : 2), periodIs);
  return { omega: r.omega, headingsDeg: r.headingsDeg, Fre: r.re, Fim: r.im };
}

/**
 * Parses a WAMIT `.4` file (response amplitude operators): lines `PER BETA I Mod Pha Re Im`, with
 * xi_i = xibar_i / L^n per unit wave amplitude, n = 0 for translations and 1 for rotations.
 * Reference: WAMIT v7 User Manual (numeric output files; definitions of the non-dimensional output).
 * @param {string} text
 * @param {{L?:number, periodIs?:string}} [opts]
 * @returns {{omega:number[], headingsDeg:number[], re:number[][][], im:number[][][]}} [heading][frequency][6], m/m, rad/m
 */
export function parseWamit4(text, { L = 1, periodIs = 'period' } = {}) {
  needPos(L, 'L');
  return parseWamitComplex(text, (i) => (i > 3 ? 1 / L : 1), periodIs);
}

/**
 * Builds a database from WAMIT `.1` and `.3` (or `.2`) file contents. The two files must list the same periods.
 * @param {{text1:string, text3?:string, rho?:number, g?:number, L?:number, depth?:number|null, periodIs?:string}} o
 * @returns {object} database (createHydroDb)
 */
export function hydroDbFromWamit({ text1, text3, rho = 1025, g = G, L = 1, depth = null, periodIs = 'period' }) {
  const c = parseWamit1(text1, { rho, L, periodIs });
  let f = { headingsDeg: [], Fre: [], Fim: [] };
  if (text3) {
    f = parseWamit3(text3, { rho, g, L, periodIs });
    if (f.omega.length !== c.omega.length || f.omega.some((w, i) => Math.abs(w - c.omega[i]) > 1e-6 * w)) throw new RangeError('WAMIT .1 and .3 files list different periods');
  }
  return createHydroDb({ omega: c.omega, A: c.A, B: c.B, headingsDeg: f.headingsDeg, Fre: f.Fre, Fim: f.Fim, Ainf: c.Ainf, A0: c.A0, rho, g, depth });
}

/**
 * Parses the JSON schema documented at the top of this file.
 * @param {string|object} json text or parsed object
 * @returns {object} database (createHydroDb)
 */
export function parseHydroJson(json) {
  const o = typeof json === 'string' ? JSON.parse(json) : json;
  if (!o || !o.omega || !o.A || !o.B) throw new RangeError('hydro JSON needs omega, A and B');
  return createHydroDb({
    omega: o.omega, A: o.A, B: o.B, headingsDeg: o.headingsDeg || [], Fre: o.Fre || [], Fim: o.Fim || [],
    Ainf: o.Ainf || null, A0: o.A0 || null, rho: o.rho === undefined ? 1025 : o.rho, g: o.g === undefined ? G : o.g,
    depth: o.depth === undefined ? null : o.depth,
  });
}

/**
 * Parses the CSV schema documented at the top of this file (header `kind,omega,headingDeg,i,j,re,im`).
 * @param {string} text
 * @param {{rho?:number, g?:number, depth?:number|null}} [opts]
 * @returns {object} database (createHydroDb)
 */
export function parseHydroCsv(text, { rho = 1025, g = G, depth = null } = {}) {
  const lines = String(text).split(/\r?\n/).map((l) => l.trim()).filter((l) => l && !l.startsWith('#'));
  if (lines.length < 2) throw new RangeError('hydro CSV is empty');
  const head = lines[0].split(',').map((s) => s.trim());
  const col = (name) => { const i = head.indexOf(name); if (i < 0) throw new RangeError(`hydro CSV: missing column '${name}'`); return i; };
  const ck = col('kind'); const cw = col('omega'); const ch = col('headingDeg'); const ci = col('i'); const cj = col('j'); const cr = col('re'); const cim = col('im');
  const freq = new Map(); const heads = new Set();
  let Ainf = null; let A0 = null;
  const rec = (w) => { if (!freq.has(w)) freq.set(w, { A: mat6(), B: mat6(), F: new Map() }); return freq.get(w); };
  for (let n = 1; n < lines.length; n++) {
    const t = lines[n].split(',').map((s) => s.trim());
    const kind = t[ck].toUpperCase();
    const i = Number(t[ci]) - 1; const j = Number(t[cj]) - 1;
    const re = Number(t[cr]); const im = t[cim] === '' || t[cim] === undefined ? 0 : Number(t[cim]);
    if (!(i >= 0 && i < 6) || !Number.isFinite(re)) throw new RangeError(`hydro CSV line ${n + 1}: bad index or value`);
    if (kind === 'AINF') { Ainf = Ainf || mat6(); Ainf[i][j] = re; continue; }
    if (kind === 'A0') { A0 = A0 || mat6(); A0[i][j] = re; continue; }
    const w = Number(t[cw]);
    if (!(w > 0)) throw new RangeError(`hydro CSV line ${n + 1}: omega must be > 0`);
    if (kind === 'A') rec(w).A[i][j] = re;
    else if (kind === 'B') rec(w).B[i][j] = re;
    else if (kind === 'F') {
      const h = Number(t[ch]);
      heads.add(h);
      const F = rec(w).F;
      if (!F.has(h)) F.set(h, { re: [0, 0, 0, 0, 0, 0], im: [0, 0, 0, 0, 0, 0] });
      F.get(h).re[i] = re; F.get(h).im[i] = im;
    } else throw new RangeError(`hydro CSV line ${n + 1}: unknown kind '${t[ck]}'`);
  }
  const omega = Array.from(freq.keys()).sort((a, b) => a - b);
  const headingsDeg = Array.from(heads).sort((a, b) => a - b);
  const zero = { re: [0, 0, 0, 0, 0, 0], im: [0, 0, 0, 0, 0, 0] };
  return createHydroDb({
    omega, A: omega.map((w) => freq.get(w).A), B: omega.map((w) => freq.get(w).B), headingsDeg,
    Fre: headingsDeg.map((h) => omega.map((w) => (freq.get(w).F.get(h) || zero).re)),
    Fim: headingsDeg.map((h) => omega.map((w) => (freq.get(w).F.get(h) || zero).im)),
    Ainf, A0, rho, g, depth,
  });
}

function bracket(xs, x) {
  const n = xs.length;
  if (n === 1 || x <= xs[0]) return [0, 0, 0];
  if (x >= xs[n - 1]) return [n - 1, n - 1, 0];
  let lo = 0; let hi = n - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid; }
  return [lo, hi, (x - xs[lo]) / (xs[hi] - xs[lo])];
}

/**
 * Coefficients at an arbitrary frequency and heading: linear interpolation in omega (clamped to the tabulated
 * range) and, for the excitation, linear interpolation of real and imaginary parts in heading (periodic when
 * the tabulated headings span the full circle, otherwise clamped).
 * @param {object} db database
 * @param {number} omega rad/s
 * @param {number} [headingDeg] deg (required for excitation)
 * @returns {{A:number[][], B:number[][], F:{re:number[], im:number[]}|null}} SI
 */
export function hydroAt(db, omega, headingDeg) {
  const [i0, i1, f] = bracket(db.omega, omega);
  const A = mat6(); const B = mat6();
  for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
    A[i][j] = db.A[i0][i][j] + f * (db.A[i1][i][j] - db.A[i0][i][j]);
    B[i][j] = db.B[i0][i][j] + f * (db.B[i1][i][j] - db.B[i0][i][j]);
  }
  let F = null;
  if (headingDeg !== undefined && db.headingsDeg.length > 0) {
    const hs = db.headingsDeg; const nh = hs.length;
    let h0; let h1; let fh;
    if (db.periodicHeading) {
      let h = headingDeg;
      while (h < hs[0]) h += 360;
      while (h >= hs[0] + 360) h -= 360;
      if (h > hs[nh - 1]) { h0 = nh - 1; h1 = 0; fh = (h - hs[nh - 1]) / (hs[0] + 360 - hs[nh - 1]); } else [h0, h1, fh] = bracket(hs, h);
    } else [h0, h1, fh] = bracket(hs, headingDeg);
    F = { re: new Array(6), im: new Array(6) };
    for (let i = 0; i < 6; i++) {
      for (const part of ['re', 'im']) {
        const src = part === 're' ? db.Fre : db.Fim;
        const a = src[h0][i0][i] + f * (src[h0][i1][i] - src[h0][i0][i]);
        const b = src[h1][i0][i] + f * (src[h1][i1][i] - src[h1][i0][i]);
        F[part][i] = a + fh * (b - a);
      }
    }
  }
  return { A, B, F };
}

/**
 * Response amplitude operators from imported coefficients:
 * [-omega² (M + A(omega)) + i omega (B(omega) + Bext) + C + Kext] xi = F(omega, beta).
 * Reference: linear frequency-domain equation of motion, Newman (1977) ch. 6; DNV-RP-C205 (2019) §7.
 * Linear potential theory; viscous damping must be supplied through Bext.
 * @param {object} o
 * @param {object} o.db database with excitation
 * @param {number[][]} o.M rigid-body mass matrix 6×6 about the database origin, kg …
 * @param {number[][]} o.C hydrostatic stiffness 6×6, N/m …
 * @param {number[][]} [o.Kext] external stiffness 6×6
 * @param {number[][]} [o.Bext] additional linear damping 6×6
 * @param {ArrayLike<number>} [o.omega] rad/s (default: the tabulated frequencies)
 * @param {number} o.headingDeg deg
 * @returns {{omega:Float64Array, period:Float64Array, headingDeg:number, depth:number, g:number, re:Float64Array[],
 *   im:Float64Array[], amp:Float64Array[], phase:Float64Array[]}} per DOF over omega, m/m and rad/m, phase rad
 */
export function raoSolve({ db, M, C, Kext, Bext, omega, headingDeg }) {
  if (!db || db.headingsDeg.length === 0) throw new RangeError('db must contain excitation forces');
  if (!Number.isFinite(headingDeg)) throw new RangeError('headingDeg must be finite');
  const om = omega || db.omega;
  const nw = om.length;
  const mk = () => Array.from({ length: 6 }, () => new Float64Array(nw));
  const re = mk(); const im = mk(); const amp = mk(); const phase = mk();
  const period = new Float64Array(nw);
  for (let iw = 0; iw < nw; iw++) {
    const w = om[iw];
    const h = hydroAt(db, w, headingDeg);
    const Ar = mat6(); const Ai = mat6();
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
      Ar[i][j] = C[i][j] + (Kext ? Kext[i][j] : 0) - w * w * (M[i][j] + h.A[i][j]);
      Ai[i][j] = w * (h.B[i][j] + (Bext ? Bext[i][j] : 0));
    }
    const x = csolve(Ar, Ai, h.F.re.slice(), h.F.im.slice());
    for (let i = 0; i < 6; i++) {
      re[i][iw] = x.re[i]; im[i][iw] = x.im[i];
      amp[i][iw] = Math.hypot(x.re[i], x.im[i]); phase[i][iw] = Math.atan2(x.im[i], x.re[i]);
    }
    period[iw] = (2 * Math.PI) / w;
  }
  return { omega: Float64Array.from(om), period, headingDeg, depth: db.depth === null ? Infinity : db.depth, g: db.g, re, im, amp, phase };
}

/**
 * Consistency checks of a hydrodynamic database:
 *  - non-negative diagonal damping, B_ii(omega) >= 0 (energy radiation);
 *  - symmetry of A and B (A_ij = A_ji), reported as the largest |X_ij - X_ji| relative to sqrt(|X_ii X_jj|);
 *  - Haskind/energy relation between damping and excitation, when the headings span the full circle:
 *      B_ii = k / (8 pi rho g Cg) ∫₀^2π |X_i(beta)|² d beta,  Cg = (omega / 2k)(1 + 2kd / sinh 2kd).
 * Reference: Newman (1977) Marine Hydrodynamics ch. 6 (Haskind relations, damping–excitation identity).
 * The Haskind check assumes zero forward speed and needs a heading spacing fine enough to integrate |X|².
 * @param {object} db database
 * @param {{tol?:number}} [opts] relative tolerance for flags (default 0.05)
 * @returns {{ok:boolean, negativeDamping:{omega:number, dof:number, value:number}[], asymmetryA:number,
 *   asymmetryB:number, haskind:{dof:number, maxRelError:number}[]|null}}
 */
export function checkHydroDb(db, { tol = 0.05 } = {}) {
  const negativeDamping = [];
  let asymmetryA = 0; let asymmetryB = 0;
  const nw = db.omega.length;
  const bScale = [0, 0, 0, 0, 0, 0];
  for (let iw = 0; iw < nw; iw++) for (let i = 0; i < 6; i++) bScale[i] = Math.max(bScale[i], Math.abs(db.B[iw][i][i]));
  for (let iw = 0; iw < nw; iw++) {
    for (let i = 0; i < 6; i++) {
      if (db.B[iw][i][i] < -1e-6 * bScale[i]) negativeDamping.push({ omega: db.omega[iw], dof: i + 1, value: db.B[iw][i][i] });
      for (let j = i + 1; j < 6; j++) {
        const sa = Math.sqrt(Math.abs(db.A[iw][i][i] * db.A[iw][j][j]));
        const sb = Math.sqrt(Math.abs(db.B[iw][i][i] * db.B[iw][j][j]));
        if (sa > 0) asymmetryA = Math.max(asymmetryA, Math.abs(db.A[iw][i][j] - db.A[iw][j][i]) / sa);
        if (sb > 0) asymmetryB = Math.max(asymmetryB, Math.abs(db.B[iw][i][j] - db.B[iw][j][i]) / sb);
      }
    }
  }
  let haskind = null;
  if (db.periodicHeading) {
    const hs = db.headingsDeg; const nh = hs.length;
    haskind = [];
    for (let i = 0; i < 6; i++) {
      if (!(bScale[i] > 0)) continue;
      let maxErr = 0;
      for (let iw = 0; iw < nw; iw++) {
        const w = db.omega[iw];
        const d = db.depth === null ? Infinity : db.depth;
        const k = dispersion(w, d, db.g);
        const cg = Number.isFinite(d) && 2 * k * d < 50 ? (w / (2 * k)) * (1 + (2 * k * d) / Math.sinh(2 * k * d)) : w / (2 * k);
        let integral = 0;
        for (let h = 0; h < nh; h++) {
          const hn = (h + 1) % nh;
          const dh = ((h === nh - 1 ? hs[0] + 360 - hs[h] : hs[hn] - hs[h]) * Math.PI) / 180;
          const x0 = db.Fre[h][iw][i] ** 2 + db.Fim[h][iw][i] ** 2;
          const x1 = db.Fre[hn][iw][i] ** 2 + db.Fim[hn][iw][i] ** 2;
          integral += 0.5 * (x0 + x1) * dh;
        }
        const bH = (k / (8 * Math.PI * db.rho * db.g * cg)) * integral;
        maxErr = Math.max(maxErr, Math.abs(bH - db.B[iw][i][i]) / bScale[i]);
      }
      haskind.push({ dof: i + 1, maxRelError: maxErr });
    }
  }
  const ok = negativeDamping.length === 0 && asymmetryA <= tol && asymmetryB <= tol && (haskind === null || haskind.every((h) => h.maxRelError <= tol));
  return { ok, negativeDamping, asymmetryA, asymmetryB, haskind };
}
