/**
 * Airfoil section data: tabulated polars, 360-degree extrapolation, rotational (3-D) stall-delay corrections,
 * generic parametric polars, section geometry from coordinate files and cavitation screening.
 *
 * Conventions: tables are stored against angle of attack in DEGREES (`alphaDeg`), evaluators take RADIANS.
 * Coefficients are non-dimensional on chord and the local relative dynamic pressure. Cm is about the quarter chord,
 * nose-up positive.
 *
 * Reynolds number: a polar is valid for one chord Reynolds number Re = rho W c / mu. Measured multi-megawatt blade
 * sections operate at Re = 3e6 - 15e6; tidal blades at 1e6 - 1e7. This module does not scale coefficients with Re:
 * supply polars measured or computed at the operating Re (several polars can be blended with `blendPolars`).
 */

const D2R = Math.PI / 180;
const R2D = 180 / Math.PI;

function need(cond, name, msg) {
  if (!cond) throw new RangeError(`${name}: ${msg}`);
}

function wrapPi(a) {
  if (a > Math.PI || a < -Math.PI) {
    a = a - 2 * Math.PI * Math.floor((a + Math.PI) / (2 * Math.PI));
  }
  return a;
}

/**
 * Chord Reynolds number.
 * @param {number} W relative speed [m/s]
 * @param {number} chord chord [m]
 * @param {number} nu kinematic viscosity [m2/s] (air 15 degC: 1.46e-5; sea water 10 degC: 1.35e-6)
 * @returns {number} Re [-]
 */
export function reynoldsNumber(W, chord, nu) {
  need(nu > 0, 'nu', 'must be > 0');
  return Math.abs(W) * chord / nu;
}

/**
 * Build a polar object with fast evaluators from tabulated data.
 * Linear interpolation in the table; outside the tabulated range the end values are held (extrapolate the polar to
 * +-180 deg first with `extrapolateViterna` for use in BEM).
 * @param {object} o
 * @param {string} [o.name]
 * @param {number[]} o.alphaDeg angle of attack, strictly ascending [deg]
 * @param {number[]} o.cl lift coefficient [-]
 * @param {number[]} o.cd drag coefficient [-]
 * @param {number[]} [o.cm] quarter-chord moment coefficient [-] (zeros when omitted)
 * @param {number} [o.thickness] relative thickness t/c [-]
 * @param {number} [o.Re] Reynolds number of the data [-]
 * @param {string} [o.label] provenance label (e.g. 'measured', 'generic parametric polar ...')
 * @returns {{name:string,label:string,thickness:(number|undefined),Re:(number|undefined),alphaDeg:number[],cl:number[],
 *   cd:number[],cm:number[],clAt:function(number):number,cdAt:function(number):number,cmAt:function(number):number}}
 *   polar; `clAt/cdAt/cmAt(alpha)` take alpha in rad (wrapped to +-pi) and return coefficients [-]
 */
export function createPolar(o) {
  need(o && Array.isArray(o.alphaDeg) || (o && ArrayBuffer.isView(o.alphaDeg)), 'alphaDeg', 'array required');
  const n = o.alphaDeg.length;
  need(n >= 2, 'alphaDeg', 'at least two rows required');
  need(o.cl && o.cl.length === n && o.cd && o.cd.length === n, 'cl/cd', 'must match alphaDeg length');
  const a = new Float64Array(n);
  const cl = new Float64Array(n);
  const cd = new Float64Array(n);
  const cm = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    a[i] = o.alphaDeg[i] * D2R;
    cl[i] = o.cl[i];
    cd[i] = o.cd[i];
    cm[i] = o.cm ? o.cm[i] : 0;
    need(Number.isFinite(a[i]) && Number.isFinite(cl[i]) && Number.isFinite(cd[i]), 'polar', `non-finite value in row ${i}`);
    if (i > 0) need(a[i] > a[i - 1], 'alphaDeg', 'must be strictly ascending');
  }
  let hint = 0;
  // interval search with a moving hint: successive BEM evaluations are close in alpha
  function locate(x) {
    if (x <= a[0]) { hint = 0; return 0; }
    if (x >= a[n - 1]) { hint = n - 2; return 1; }
    let i = hint;
    if (x < a[i] || x > a[i + 1]) {
      let lo = 0;
      let hi = n - 1;
      while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (a[mid] <= x) lo = mid; else hi = mid;
      }
      i = lo;
      hint = i;
    }
    return (x - a[i]) / (a[i + 1] - a[i]);
  }
  const at = (arr) => (alpha) => {
    const s = locate(wrapPi(alpha));
    const i = hint;
    return arr[i] + s * (arr[i + 1] - arr[i]);
  };
  return {
    name: o.name || 'polar',
    label: o.label || 'user-supplied polar',
    thickness: o.thickness,
    Re: o.Re,
    alphaDeg: Array.from(o.alphaDeg),
    cl: Array.from(cl),
    cd: Array.from(cd),
    cm: Array.from(cm),
    clAt: at(cl),
    cdAt: at(cd),
    cmAt: at(cm),
  };
}

function numericRows(text) {
  const rows = [];
  const lines = String(text).split(/\r?\n/);
  for (let li = 0; li < lines.length; li++) {
    const raw = lines[li].split('!')[0].trim();
    if (!raw) { rows.push(null); continue; }
    const tok = raw.split(/[\s,;]+/).filter((s) => s.length);
    const nums = [];
    for (const t of tok) {
      if (!/^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$/.test(t)) break;
      nums.push(Number(t.replace(/[dD]/, 'e')));
    }
    rows.push(nums.length === tok.length && nums.length > 0 ? nums : null);
  }
  return { rows, lines };
}

function longestTable(rows, minCols) {
  let best = [0, 0];
  let start = -1;
  for (let i = 0; i <= rows.length; i++) {
    const ok = i < rows.length && rows[i] && rows[i].length >= minCols &&
      (start < 0 || i === start || rows[i][0] > rows[i - 1][0]);
    if (ok) {
      if (start < 0) start = i;
    } else {
      if (start >= 0 && i - start > best[1] - best[0]) best = [start, i];
      start = (i < rows.length && rows[i] && rows[i].length >= minCols) ? i : -1;
    }
  }
  return best;
}

/**
 * Parse an AeroDyn-style airfoil table (AeroDyn v13 or v15 text files): free header lines followed by a block of
 * rows `alpha[deg] Cl Cd [Cm]`. The longest block of rows with ascending first column is taken as the table; a
 * header line mentioning "Reynolds" or the keyword `Re` is read as the Reynolds number in millions (AeroDyn
 * convention). Other header parameters (stall angle, Beddoes-Leishman constants) are ignored.
 * @param {string} text file contents
 * @param {object} [opts]
 * @param {string} [opts.name]
 * @param {number} [opts.thickness] relative thickness t/c [-]
 * @returns {object} polar (see createPolar)
 */
export function parseAeroDynPolar(text, opts = {}) {
  const { rows, lines } = numericRows(text);
  // numericRows rejects lines with trailing words, so header lines of the form "value  Name  description" are skipped
  const [s, e] = longestTable(rows, 3);
  need(e - s >= 2, 'text', 'no table of alpha, Cl, Cd rows found');
  let Re;
  for (let i = 0; i < s; i++) {
    const m = /^\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+.*\b(Reynolds|Re)\b/i.exec(lines[i]);
    if (m) { Re = Number(m[1]) * 1e6; break; }
  }
  const tab = rows.slice(s, e);
  return createPolar({
    name: opts.name || (lines[0] || 'AeroDyn polar').trim().slice(0, 60),
    alphaDeg: tab.map((r) => r[0]), cl: tab.map((r) => r[1]), cd: tab.map((r) => r[2]),
    cm: tab[0].length >= 4 ? tab.map((r) => (r.length >= 4 ? r[3] : 0)) : undefined,
    thickness: opts.thickness, Re, label: 'imported AeroDyn-style table',
  });
}

/**
 * Parse a simple delimited polar: rows `alpha[deg], Cl, Cd[, Cm]`, separated by comma, semicolon or white space,
 * with an optional header row. Rows are sorted by alpha.
 * @param {string} text file contents
 * @param {object} [opts]
 * @param {string} [opts.name]
 * @param {number} [opts.thickness] t/c [-]
 * @param {number} [opts.Re] Reynolds number [-]
 * @returns {object} polar (see createPolar)
 */
export function parseCsvPolar(text, opts = {}) {
  const { rows } = numericRows(text);
  const tab = rows.filter((r) => r && r.length >= 3).sort((p, q) => p[0] - q[0]);
  need(tab.length >= 2, 'text', 'no rows of alpha, Cl, Cd found');
  return createPolar({
    name: opts.name || 'CSV polar',
    alphaDeg: tab.map((r) => r[0]), cl: tab.map((r) => r[1]), cd: tab.map((r) => r[2]),
    cm: tab.every((r) => r.length >= 4) ? tab.map((r) => r[3]) : undefined,
    thickness: opts.thickness, Re: opts.Re, label: 'imported CSV table',
  });
}

function lin(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid] <= x) lo = mid; else hi = mid;
  }
  return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo]);
}

/**
 * Parse an airfoil coordinate file (geometry only, no aerodynamic data are derived).
 * Selig format: name line, then x y pairs running from the trailing edge over the upper surface to the leading
 * edge and back along the lower surface. Lednicer format (a line with the two point counts, then upper and lower
 * surfaces each from leading to trailing edge) is detected and converted.
 * Thickness is measured perpendicular to the chord line (the "British" convention), camber as the mid-line between
 * upper and lower surface, both sampled at 201 cosine-spaced chord stations. The leading-edge radius is a
 * least-squares osculating circle through the points within `leWindow` of the leading edge:
 * (x - r)^2 + y^2 = r^2  ->  r = sum x (x^2 + y^2) / (2 sum x^2); its accuracy depends on the point density there.
 * @param {string} text file contents
 * @param {object} [opts]
 * @param {number} [opts.leWindow=0.005] chord fraction used for the leading-edge circle fit [-]
 * @returns {{name:string,x:number[],yUpper:number[],yLower:number[],thickness:number,thicknessPos:number,
 *   camber:number,camberPos:number,leRadius:number,teThickness:number,area:number,nPoints:number}}
 *   all lengths as fractions of chord [-]; area as fraction of chord^2
 */
export function parseSeligCoordinates(text, opts = {}) {
  const leWindow = opts.leWindow ?? 0.005;
  const lines = String(text).split(/\r?\n/);
  const { rows } = numericRows(text);
  let name = 'airfoil';
  for (let i = 0; i < lines.length; i++) {
    if (lines[i].trim()) { if (!rows[i]) name = lines[i].trim(); break; }
  }
  let pts = rows.filter((r) => r && r.length >= 2).map((r) => [r[0], r[1]]);
  need(pts.length >= 7, 'text', 'fewer than 7 coordinate pairs found');
  if (pts[0][0] > 1.5 && pts[0][1] > 1.5) {
    // Lednicer: counts line, upper LE->TE, lower LE->TE
    const nu = Math.round(pts[0][0]);
    const nl = Math.round(pts[0][1]);
    need(pts.length >= 1 + nu + nl, 'text', 'Lednicer point counts exceed the data');
    const up = pts.slice(1, 1 + nu).reverse();
    const lo = pts.slice(1 + nu, 1 + nu + nl);
    pts = up.concat(lo.slice(1));
  }
  let iLe = 0;
  let xMin = Infinity;
  let xMax = -Infinity;
  for (let i = 0; i < pts.length; i++) {
    if (pts[i][0] < xMin) { xMin = pts[i][0]; iLe = i; }
    if (pts[i][0] > xMax) xMax = pts[i][0];
  }
  const c = xMax - xMin;
  need(c > 0, 'text', 'degenerate chord');
  const yLe = pts[iLe][1];
  const norm = pts.map((p) => [(p[0] - xMin) / c, (p[1] - yLe) / c]);
  let a = norm.slice(0, iLe + 1).reverse();
  let b = norm.slice(iLe);
  need(a.length >= 3 && b.length >= 3, 'text', 'cannot split into upper and lower surfaces');
  const meanY = (s) => s.reduce((t, p) => t + p[1], 0) / s.length;
  if (meanY(a) < meanY(b)) { const t = a; a = b; b = t; }
  const mono = (s) => {
    const xs = [s[0][0]];
    const ys = [s[0][1]];
    for (let i = 1; i < s.length; i++) {
      if (s[i][0] > xs[xs.length - 1] + 1e-12) { xs.push(s[i][0]); ys.push(s[i][1]); }
    }
    return { xs, ys };
  };
  const U = mono(a);
  const L = mono(b);
  const m = 201;
  const x = new Array(m);
  const yu = new Array(m);
  const yl = new Array(m);
  let t = 0; let tx = 0; let cam = 0; let camx = 0; let area = 0;
  for (let i = 0; i < m; i++) {
    x[i] = 0.5 * (1 - Math.cos(Math.PI * i / (m - 1)));
    yu[i] = lin(U.xs, U.ys, x[i]);
    yl[i] = lin(L.xs, L.ys, x[i]);
    const th = yu[i] - yl[i];
    const cm = 0.5 * (yu[i] + yl[i]);
    if (th > t) { t = th; tx = x[i]; }
    if (Math.abs(cm) > Math.abs(cam)) { cam = cm; camx = x[i]; }
    if (i > 0) area += 0.5 * (th + (yu[i - 1] - yl[i - 1])) * (x[i] - x[i - 1]);
  }
  // remove the chord-line offset of the trailing edge from the camber (chord line joins LE and TE mid-point)
  const yTe = 0.5 * (yu[m - 1] + yl[m - 1]);
  cam = 0; camx = 0;
  for (let i = 0; i < m; i++) {
    const cm = 0.5 * (yu[i] + yl[i]) - yTe * x[i];
    if (Math.abs(cm) > Math.abs(cam)) { cam = cm; camx = x[i]; }
  }
  let sxx = 0; let sxr = 0; let cnt = 0;
  let win = leWindow;
  while (cnt < 2 && win < 0.1) {
    sxx = 0; sxr = 0; cnt = 0;
    for (const p of norm) {
      if (p[0] > 1e-9 && p[0] <= win) { sxx += p[0] * p[0]; sxr += p[0] * (p[0] * p[0] + p[1] * p[1]); cnt++; }
    }
    if (cnt < 2) win *= 2;
  }
  return {
    name, x, yUpper: yu, yLower: yl,
    thickness: t, thicknessPos: tx, camber: cam, camberPos: camx,
    leRadius: cnt >= 1 ? sxr / (2 * sxx) : NaN,
    teThickness: yu[m - 1] - yl[m - 1], area, nPoints: pts.length,
  };
}

/**
 * Maximum (flat-plate-like) drag coefficient at 90 deg as a function of blade aspect ratio.
 * Cd,max = 1.11 + 0.018 AR for AR <= 50, 2.01 above.
 * Reference: Viterna, L.A. and Corrigan, R.D. (1982), "Fixed pitch rotor performance of large horizontal axis wind
 * turbines", NASA CP-2230.
 * @param {number} aspectRatio blade aspect ratio (tip radius / representative chord) [-]
 * @returns {number} Cd,max [-]
 */
export function viternaCdMax(aspectRatio) {
  need(aspectRatio > 0, 'aspectRatio', 'must be > 0');
  return aspectRatio <= 50 ? 1.11 + 0.018 * aspectRatio : 2.01;
}

/**
 * Extrapolate a polar to the full +-180 deg range.
 *
 * Positive post-stall range alpha_s..90 deg: Viterna-Corrigan (1982) flat-plate model matched at the stall point
 * (alpha_s, Cl_s, Cd_s):
 *   Cd = B1 sin^2(alpha) + B2 cos(alpha),            B1 = Cd,max,  B2 = (Cd_s - Cd,max sin^2 alpha_s) / cos alpha_s
 *   Cl = A1 sin(2 alpha) + A2 cos^2(alpha)/sin(alpha), A1 = B1/2,  A2 = (Cl_s - Cd,max sin alpha_s cos alpha_s) sin alpha_s / cos^2 alpha_s
 * Beyond 90 deg and on the negative side the same functions are reflected, with the lift of the reflected
 * branches multiplied by `backFactor` (0.7 is the value used by the NREL AirfoilPrep tool to represent the poorer
 * lift of an airfoil in reversed/inverted flow; it is a convention, not part of the 1982 paper), lift taken linearly
 * to zero at +-180 deg, and a linear bridge between -alpha_s and the first tabulated angle.
 * Cm outside the table (assumption, documented): centre of pressure moves linearly with |alpha| from its position at
 * the matching point to mid-chord at 90 deg (plate normal to the flow) and to 0.75 c at 180 deg (reversed flow).
 * @param {object} polar polar from createPolar covering at least the attached range and positive stall
 * @param {object} [o]
 * @param {number} [o.aspectRatio=10] blade aspect ratio for Cd,max [-]
 * @param {number} [o.cdMax] override of Cd,max [-]; default max(table Cd, viternaCdMax(aspectRatio))
 * @param {number} [o.alphaMatchDeg] matching (stall) angle [deg]; default the last tabulated angle
 * @param {number} [o.backFactor=0.7] lift factor for the reflected branches [-] (0.6 - 1.0)
 * @param {number} [o.cdMin=0.001] lower bound on Cd [-]
 * @param {number} [o.stepDeg=2] spacing of the generated points [deg]
 * @returns {object} polar covering -180..180 deg
 */
export function extrapolateViterna(polar, o = {}) {
  const AR = o.aspectRatio ?? 10;
  const back = o.backFactor ?? 0.7;
  const cdFloor = o.cdMin ?? 0.001;
  const step = o.stepDeg ?? 2;
  need(step > 0 && step <= 10, 'stepDeg', 'must be in (0, 10]');
  need(back > 0 && back <= 1, 'backFactor', 'must be in (0, 1]');
  let aDeg = polar.alphaDeg.slice();
  let cl = polar.cl.slice();
  let cd = polar.cd.slice();
  let cm = polar.cm.slice();
  if (o.alphaMatchDeg !== undefined) {
    need(o.alphaMatchDeg > aDeg[0] && o.alphaMatchDeg <= aDeg[aDeg.length - 1], 'alphaMatchDeg', 'must lie inside the table');
    const am = o.alphaMatchDeg;
    const clm = polar.clAt(am * D2R);
    const cdm = polar.cdAt(am * D2R);
    const cmm = polar.cmAt(am * D2R);
    let k = aDeg.findIndex((v) => v >= am - 1e-9);
    aDeg = aDeg.slice(0, k); cl = cl.slice(0, k); cd = cd.slice(0, k); cm = cm.slice(0, k);
    aDeg.push(am); cl.push(clm); cd.push(cdm); cm.push(cmm);
  }
  const n = aDeg.length;
  const ah = aDeg[n - 1] * D2R;
  need(ah > 0 && ah < Math.PI / 2, 'polar', 'last tabulated (matching) angle must be between 0 and 90 deg');
  const clH = cl[n - 1];
  const cdH = cd[n - 1];
  const al = aDeg[0] * D2R;
  const clL = cl[0];
  const cdL = cd[0];
  const cdMax = o.cdMax ?? Math.max(Math.max(...cd), viternaCdMax(AR));
  const sh = Math.sin(ah);
  const ch = Math.cos(ah);
  const A2 = (clH - cdMax * sh * ch) * sh / (ch * ch);
  const B2 = (cdH - cdMax * sh * sh) / ch;
  const vCl = (x) => { const xx = Math.max(x, 1e-4); return 0.5 * cdMax * Math.sin(2 * xx) + A2 * Math.cos(xx) ** 2 / Math.sin(xx); };
  const vCd = (x) => { const xx = Math.max(x, 1e-4); return cdMax * Math.sin(xx) ** 2 + B2 * Math.cos(xx); };
  // centre-of-pressure model for Cm outside the table
  const cnH = clH * ch + cdH * sh;
  let xcpH = Math.abs(cnH) > 1e-6 ? 0.25 - cm[n - 1] / cnH : 0.25;
  xcpH = Math.min(0.5, Math.max(0.1, xcpH));
  const xcp = (absA) => (absA <= ah ? xcpH
    : absA <= Math.PI / 2 ? xcpH + (0.5 - xcpH) * (absA - ah) / (Math.PI / 2 - ah)
      : 0.5 + 0.25 * (absA - Math.PI / 2) / (Math.PI / 2));
  const outA = []; const outCl = []; const outCd = []; const outCm = [];
  const push = (aRad, l, d) => {
    const dd = Math.max(d, cdFloor);
    outA.push(aRad * R2D); outCl.push(l); outCd.push(dd);
    const cn = l * Math.cos(aRad) + dd * Math.sin(aRad);
    outCm.push(-cn * (xcp(Math.abs(aRad)) - 0.25));
  };
  const pi = Math.PI;
  const dS = step * D2R;
  const grid = (from, to) => { // open interval points on a regular grid
    const g = [];
    const nSeg = Math.max(1, Math.ceil((to - from) / dS - 1e-9));
    for (let i = 1; i < nSeg; i++) g.push(from + (to - from) * i / nSeg);
    return g;
  };
  // -180 .. -180 + ah : linear lift
  push(-pi, 0, vCd(0));
  for (const x of grid(-pi, -pi + ah)) push(x, (x + pi) / ah * clH * back, vCd(x + pi));
  // -180 + ah .. -90
  push(-pi + ah, vCl(ah) * back, vCd(ah));
  for (const x of grid(-pi + ah, -pi / 2)) push(x, vCl(x + pi) * back, vCd(x + pi));
  push(-pi / 2, 0, vCd(pi / 2));
  // -90 .. -ah
  const lowEnd = Math.min(-ah, al);
  for (const x of grid(-pi / 2, -ah)) if (x < al - 1e-9) push(x, -vCl(-x) * back, vCd(-x));
  if (al > -ah + 1e-9) {
    push(-ah, -clH * back, cdH);
    for (const x of grid(-ah, al)) {
      const s = (x + ah) / (al + ah);
      push(x, -clH * back + s * (clL + clH * back), cdH + s * (cdL - cdH));
    }
  }
  void lowEnd;
  for (let i = 0; i < n; i++) { outA.push(aDeg[i]); outCl.push(cl[i]); outCd.push(Math.max(cd[i], cdFloor)); outCm.push(cm[i]); }
  for (const x of grid(ah, pi / 2)) push(x, vCl(x), vCd(x));
  push(pi / 2, 0, vCd(pi / 2));
  for (const x of grid(pi / 2, pi - ah)) push(x, -vCl(pi - x) * back, vCd(pi - x));
  push(pi - ah, -vCl(ah) * back, vCd(ah));
  for (const x of grid(pi - ah, pi)) push(x, (x - pi) / ah * clH * back, vCd(pi - x));
  push(pi, 0, vCd(0));
  // guard against duplicates created by a table that already extends beyond -ah
  const A = []; const CL = []; const CD = []; const CM = [];
  const idx = outA.map((v, i) => i).sort((p, q) => outA[p] - outA[q]);
  for (const i of idx) {
    if (A.length && outA[i] - A[A.length - 1] < 1e-9) continue;
    A.push(outA[i]); CL.push(outCl[i]); CD.push(outCd[i]); CM.push(outCm[i]);
  }
  const res = createPolar({
    name: polar.name, alphaDeg: A, cl: CL, cd: CD, cm: CM, thickness: polar.thickness, Re: polar.Re,
    label: `${polar.label}; Viterna-Corrigan 360-deg extrapolation (Cd,max = ${cdMax.toFixed(3)})`,
  });
  res.cdMax = cdMax;
  res.viterna = { A1: 0.5 * cdMax, A2, B1: cdMax, B2, alphaMatchDeg: ah * R2D };
  return res;
}

/**
 * Zero-lift angle and attached-flow lift-curve slope from a table (least squares over alpha0 - 2 deg .. alpha0 +
 * `spanDeg`).
 * @param {object} polar polar from createPolar
 * @param {object} [o]
 * @param {number} [o.spanDeg=6] width of the fit window above the zero-lift angle [deg]
 * @returns {{alpha0:number, slope:number}} zero-lift angle [rad], lift-curve slope [1/rad]
 */
export function polarLinearFit(polar, o = {}) {
  const span = o.spanDeg ?? 6;
  const a = polar.alphaDeg;
  const cl = polar.cl;
  let best = null;
  for (let i = 0; i < a.length - 1; i++) {
    if (a[i] < -30 || a[i + 1] > 30) continue;
    if (cl[i] === 0 || cl[i] * cl[i + 1] < 0 || cl[i + 1] === 0) {
      if (cl[i + 1] <= cl[i]) continue;
      const z = a[i] - cl[i] * (a[i + 1] - a[i]) / (cl[i + 1] - cl[i]);
      if (best === null || Math.abs(z) < Math.abs(best)) best = z;
    }
  }
  need(best !== null, 'polar', 'no zero-lift crossing with positive slope found within +-30 deg');
  let sx = 0; let sy = 0; let sxx = 0; let sxy = 0; let m = 0;
  const nS = 25;
  for (let i = 0; i <= nS; i++) {
    const x = (best - 2 + (span + 2) * i / nS) * D2R;
    const y = polar.clAt(x);
    sx += x; sy += y; sxx += x * x; sxy += x * y; m++;
  }
  const slope = (m * sxy - sx * sy) / (m * sxx - sx * sx);
  return { alpha0: (sx * slope - sy) / (m * slope) , slope };
}

/**
 * Rotational stall delay, Du-Selig lift correction with Eggers drag adjustment.
 *   f_L = 1/(2 pi) [ 1.6 (c/r)/0.1267 * (a - (c/r)^(d R/(Lambda r))) / (b + (c/r)^(d R/(Lambda r))) - 1 ]
 *   Cl,3D = Cl,2D + f_L (Cl,p - Cl,2D),   Cl,p = m (alpha - alpha0),   Lambda = Omega R / sqrt(V^2 + (Omega R)^2)
 * with a = b = d = 1 (Du, Z. and Selig, M.S. (1998), "A 3-D stall-delay model for horizontal axis wind turbine
 * performance prediction", AIAA-98-0021). The drag change follows Eggers, Chaney and Kashyap (2003) as implemented in
 * NREL AirfoilPrep: dCd = dCl (sin alpha - k cos alpha)/(cos alpha + k sin alpha), k = 0.12. Above `alphaMaxCorrDeg`
 * the correction is faded with ((90 - alpha)/(90 - alphaMaxCorr))^2 (AirfoilPrep convention). Only alpha0 <= alpha <=
 * 90 deg is modified. Valid for inboard attached/lightly stalled sections of horizontal-axis rotors; use c/r < 1.
 * @param {object} polar 2-D polar (createPolar)
 * @param {object} o
 * @param {number} o.cOverR local chord / local radius c/r [-]
 * @param {number} o.rOverR local radius / tip radius r/R [-]
 * @param {number} o.tsr tip-speed ratio Omega R / V at the design/operating point [-]
 * @param {number} [o.a=1] empirical constant a [-]
 * @param {number} [o.b=1] empirical constant b [-]
 * @param {number} [o.d=1] empirical constant d [-]
 * @param {number|'fit'} [o.liftSlope=2*Math.PI] m [1/rad]; 'fit' uses the slope of the table (AirfoilPrep practice)
 * @param {number} [o.eggersK=0.12] Eggers drag constant k [-]; 0 leaves the drag unchanged
 * @param {number} [o.alphaMaxCorrDeg=30] angle above which the correction is faded out [deg]
 * @returns {object} corrected polar
 */
export function stallDelayDuSelig(polar, o) {
  need(o && o.cOverR > 0 && o.cOverR < 1.5, 'cOverR', 'must be in (0, 1.5)');
  need(o.rOverR > 0 && o.rOverR <= 1, 'rOverR', 'must be in (0, 1]');
  need(o.tsr > 0, 'tsr', 'must be > 0');
  const a = o.a ?? 1; const b = o.b ?? 1; const d = o.d ?? 1;
  const k = o.eggersK ?? 0.12;
  const aMax = (o.alphaMaxCorrDeg ?? 30) * D2R;
  const fit = polarLinearFit(polar);
  const m = o.liftSlope === 'fit' ? fit.slope : (o.liftSlope ?? 2 * Math.PI);
  const lam = o.tsr / Math.sqrt(1 + o.tsr * o.tsr);
  const ex = d / (lam * o.rOverR);
  const p = Math.pow(o.cOverR, ex);
  const fL = Math.max(0, (1 / m) * (1.6 * o.cOverR / 0.1267 * (a - p) / (b + p) - 1));
  const cl = polar.cl.slice();
  const cd = polar.cd.slice();
  for (let i = 0; i < cl.length; i++) {
    const al = polar.alphaDeg[i] * D2R;
    if (al < fit.alpha0 || al > Math.PI / 2) continue;
    const adj = al > aMax ? ((Math.PI / 2 - al) / (Math.PI / 2 - aMax)) ** 2 : 1;
    const dCl = fL * (m * (al - fit.alpha0) - cl[i]) * adj;
    if (dCl <= 0) continue;
    cl[i] += dCl;
    cd[i] += dCl * (Math.sin(al) - k * Math.cos(al)) / (Math.cos(al) + k * Math.sin(al));
    if (cd[i] < 1e-4) cd[i] = 1e-4;
  }
  const res = createPolar({ ...polar, cl, cd, label: `${polar.label}; Du-Selig/Eggers rotational correction (c/r = ${o.cOverR.toFixed(3)})` });
  res.stallDelayFactor = fL;
  return res;
}

/**
 * Rotational stall delay, Snel's correction: Cl,3D = Cl,2D + k (c/r)^2 (Cl,p - Cl,2D), k = 3, Cl,p = m (alpha - alpha0).
 * Reference: Snel, H., Houwink, R. and Bosschers, J. (1994), "Sectional prediction of lift coefficients on rotating
 * wind turbine blades in stall", ECN-C-93-052. Drag is not modified. Applied for alpha0 <= alpha <= `alphaEndDeg`,
 * faded linearly to zero between `alphaFadeDeg` and `alphaEndDeg` (fade angles are a common implementation choice).
 * @param {object} polar 2-D polar
 * @param {object} o
 * @param {number} o.cOverR c/r [-]
 * @param {number} [o.k=3] Snel constant [-]
 * @param {number} [o.liftSlope=2*Math.PI] m [1/rad]
 * @param {number} [o.alphaFadeDeg=30] start of fade-out [deg]
 * @param {number} [o.alphaEndDeg=50] end of fade-out [deg]
 * @returns {object} corrected polar
 */
export function stallDelaySnel(polar, o) {
  need(o && o.cOverR > 0 && o.cOverR < 1.5, 'cOverR', 'must be in (0, 1.5)');
  const k = o.k ?? 3;
  const m = o.liftSlope ?? 2 * Math.PI;
  const f0 = (o.alphaFadeDeg ?? 30) * D2R;
  const f1 = (o.alphaEndDeg ?? 50) * D2R;
  const fit = polarLinearFit(polar);
  const g = Math.min(1, k * o.cOverR * o.cOverR);
  const cl = polar.cl.slice();
  for (let i = 0; i < cl.length; i++) {
    const al = polar.alphaDeg[i] * D2R;
    if (al < fit.alpha0 || al > f1) continue;
    const fade = al <= f0 ? 1 : (f1 - al) / (f1 - f0);
    const dCl = g * (m * (al - fit.alpha0) - cl[i]) * fade;
    if (dCl > 0) cl[i] += dCl;
  }
  return createPolar({ ...polar, cl, label: `${polar.label}; Snel 3 (c/r)^2 rotational correction` });
}

/**
 * Linear blend of two polars on the union of their angle grids: result = (1 - w) A + w B.
 * @param {object} pA polar A
 * @param {object} pB polar B
 * @param {number} w weight of B, 0..1 [-]
 * @returns {object} blended polar (thickness blended when both have one)
 */
export function blendPolars(pA, pB, w) {
  need(w >= 0 && w <= 1, 'w', 'must be in [0, 1]');
  const set = new Set([...pA.alphaDeg, ...pB.alphaDeg].map((v) => Math.round(v * 1e6) / 1e6));
  const aDeg = [...set].sort((p, q) => p - q);
  const mix = (fa, fb) => aDeg.map((v) => (1 - w) * fa(v * D2R) + w * fb(v * D2R));
  return createPolar({
    name: `${pA.name}/${pB.name} blend`,
    alphaDeg: aDeg, cl: mix(pA.clAt, pB.clAt), cd: mix(pA.cdAt, pB.cdAt), cm: mix(pA.cmAt, pB.cmAt),
    thickness: pA.thickness !== undefined && pB.thickness !== undefined ? (1 - w) * pA.thickness + w * pB.thickness : undefined,
    Re: pA.Re, label: `blend of [${pA.label}] and [${pB.label}]`,
  });
}

/**
 * Interpolate a polar for a given relative thickness from a family of polars that each carry `thickness`.
 * Linear in thickness between the two bracketing family members; clamped outside the family range.
 * @param {object[]} polars family, each with `.thickness` [-]
 * @param {number} thickness target t/c [-]
 * @returns {object} polar
 */
export function polarByThickness(polars, thickness) {
  need(Array.isArray(polars) && polars.length >= 1, 'polars', 'non-empty array required');
  const fam = polars.slice().sort((p, q) => p.thickness - q.thickness);
  for (const p of fam) need(p.thickness > 0, 'polars', 'every polar needs a thickness');
  if (thickness <= fam[0].thickness) return fam[0];
  if (thickness >= fam[fam.length - 1].thickness) return fam[fam.length - 1];
  let i = 0;
  while (fam[i + 1].thickness < thickness) i++;
  const w = (thickness - fam[i].thickness) / (fam[i + 1].thickness - fam[i].thickness);
  if (w < 1e-9) return fam[i];
  if (w > 1 - 1e-9) return fam[i + 1];
  return blendPolars(fam[i], fam[i + 1], w);
}

/**
 * Potential-flow lift-curve slope with the first-order thickness correction of a Joukowski section,
 * dCl/dalpha = 2 pi (1 + 0.77 t/c), multiplied by a viscous efficiency `eta`.
 * Reference: Abbott and von Doenhoff, Theory of Wing Sections (thin-airfoil theory 2 pi; Joukowski thickness factor
 * 1 + 0.77 t/c). Real sections reach 85-95 % of the potential value.
 * @param {number} thickness t/c [-]
 * @param {number} [eta=0.9] viscous efficiency [-] (0.85 - 0.95)
 * @returns {number} lift-curve slope [1/rad]
 */
export function liftSlopeThickness(thickness, eta = 0.9) {
  need(thickness >= 0 && thickness <= 1, 'thickness', 'must be in [0, 1]');
  return 2 * Math.PI * (1 + 0.77 * thickness) * eta;
}

/**
 * Editable default parameters of the generic parametric polar as a function of relative thickness.
 * GENERIC values of representative order of magnitude for smooth wind-turbine sections at Re ~ 6e6; they are NOT
 * the data of any DU, NACA, FFA or other named airfoil. Rows are interpolated linearly in thickness.
 * Columns: thickness t/c [-], zero-lift angle [deg], stall angle [deg], Cl,max [-], minimum Cd [-],
 * drag-bucket curvature k in Cd = Cd,min + k (Cl,lin - Cl,minDrag)^2 [-].
 */
export const GENERIC_SECTION_DEFAULTS = [
  { thickness: 0.18, alpha0Deg: -4.0, alphaStallDeg: 13.0, clMax: 1.45, cdMin: 0.0060, kBucket: 0.0050 },
  { thickness: 0.21, alpha0Deg: -4.0, alphaStallDeg: 12.5, clMax: 1.40, cdMin: 0.0065, kBucket: 0.0055 },
  { thickness: 0.25, alpha0Deg: -3.5, alphaStallDeg: 12.0, clMax: 1.35, cdMin: 0.0075, kBucket: 0.0060 },
  { thickness: 0.30, alpha0Deg: -3.0, alphaStallDeg: 11.0, clMax: 1.30, cdMin: 0.0090, kBucket: 0.0080 },
  { thickness: 0.35, alpha0Deg: -3.0, alphaStallDeg: 10.0, clMax: 1.25, cdMin: 0.0110, kBucket: 0.0110 },
  { thickness: 0.40, alpha0Deg: -3.0, alphaStallDeg: 9.0, clMax: 1.15, cdMin: 0.0140, kBucket: 0.0150 },
];

function genericDefaults(t) {
  const T = GENERIC_SECTION_DEFAULTS;
  const xs = T.map((r) => r.thickness);
  const g = (key) => lin(xs, T.map((r) => r[key]), t);
  return { alpha0Deg: g('alpha0Deg'), alphaStallDeg: g('alphaStallDeg'), clMax: g('clMax'), cdMin: g('cdMin'), kBucket: g('kBucket') };
}

/**
 * Generic parametric polar — replace with measured/computed polars for design.
 *
 * Model (all parameters are arguments; defaults from GENERIC_SECTION_DEFAULTS by thickness):
 *  - attached flow: Cl = m (alpha - alpha0), m = liftSlopeThickness(t, eta);
 *  - stall rounding: a parabola tangent to the attached line that peaks with zero slope at (alphaStall, clMax);
 *    this needs m (alphaStall - alpha0)/2 <= clMax <= m (alphaStall - alpha0), otherwise RangeError;
 *  - negative stall: mirror image about alpha0 scaled by `negStallFactor`;
 *  - drag: Cd = cdMin + kBucket (m (alpha - alpha0) - clMinDrag)^2 between the two stall angles (parabolic bucket);
 *  - beyond stall: Viterna-Corrigan with aspect-ratio-dependent Cd,max (extrapolateViterna);
 *  - moment: thin-airfoil value of a parabolic camber line with the same zero-lift angle, Cm = pi alpha0 / 2
 *    (alpha0 in rad), in attached flow; centre-of-pressure model of extrapolateViterna beyond stall.
 * @param {object} o
 * @param {number} o.thickness relative thickness t/c [-], 0.05 - 0.6
 * @param {number} [o.alpha0Deg] zero-lift angle [deg]
 * @param {number} [o.alphaStallDeg] positive stall angle [deg]
 * @param {number} [o.clMax] maximum lift coefficient [-]
 * @param {number} [o.cdMin] minimum drag coefficient [-]
 * @param {number} [o.kBucket] drag-bucket curvature [-]
 * @param {number} [o.clMinDrag] lift coefficient of minimum drag [-]; default the lift at zero incidence, -m alpha0
 * @param {number} [o.eta=0.9] viscous lift-slope efficiency [-]
 * @param {number} [o.liftSlope] override of m [1/rad]
 * @param {number} [o.negStallFactor=1] ratio of negative to positive stall margin [-]
 * @param {number} [o.aspectRatio=17] blade aspect ratio for Viterna Cd,max [-]
 * @param {number} [o.backFactor=0.7] see extrapolateViterna
 * @param {string} [o.name]
 * @returns {object} polar covering -180..180 deg, labelled as generic
 */
export function genericPolar(o) {
  need(o && o.thickness >= 0.05 && o.thickness <= 0.6, 'thickness', 'must be in [0.05, 0.6]');
  const def = genericDefaults(o.thickness);
  const a0 = (o.alpha0Deg ?? def.alpha0Deg) * D2R;
  const as = (o.alphaStallDeg ?? def.alphaStallDeg) * D2R;
  const clMax = o.clMax ?? def.clMax;
  const cdMin = o.cdMin ?? def.cdMin;
  const kB = o.kBucket ?? def.kBucket;
  const m = o.liftSlope ?? liftSlopeThickness(o.thickness, o.eta ?? 0.9);
  const nf = o.negStallFactor ?? 1;
  need(as > a0, 'alphaStallDeg', 'must exceed alpha0Deg');
  need(cdMin > 0, 'cdMin', 'must be > 0');
  need(kB >= 0, 'kBucket', 'must be >= 0');
  need(nf > 0 && nf <= 1.5, 'negStallFactor', 'must be in (0, 1.5]');
  const span = as - a0;
  need(clMax <= m * span * (1 + 1e-9) && clMax >= 0.5 * m * span, 'clMax',
    `must lie between ${(0.5 * m * span).toFixed(3)} and ${(m * span).toFixed(3)} for the given slope and stall angle`);
  const dRound = 2 * (span - clMax / m);       // width of the rounding parabola before stall [rad]
  const kPar = dRound > 1e-12 ? m / (2 * dRound) : 0;
  const clMinDrag = o.clMinDrag ?? -m * a0;
  const clPos = (x) => (x <= as - dRound ? m * (x - a0) : clMax - kPar * (x - as) ** 2);
  const asN = a0 - nf * span;
  const aDeg = []; const cl = []; const cd = []; const cm = [];
  const nPts = Math.max(24, Math.ceil((as - asN) * R2D / 0.5));
  for (let i = 0; i <= nPts; i++) {
    const x = asN + (as - asN) * i / nPts;
    aDeg.push(x * R2D);
    cl.push(x >= a0 ? clPos(x) : -nf * clPos(a0 + (a0 - x) / nf));
    cd.push(cdMin + kB * (m * (x - a0) - clMinDrag) ** 2);
    cm.push(0.5 * Math.PI * a0);
  }
  const base = createPolar({
    name: o.name || `generic-${Math.round(o.thickness * 100)}`, alphaDeg: aDeg, cl, cd, cm, thickness: o.thickness,
    label: 'generic parametric polar — replace with measured/computed polars for design',
  });
  const full = extrapolateViterna(base, { aspectRatio: o.aspectRatio ?? 17, backFactor: o.backFactor ?? 0.7, stepDeg: 2 });
  full.generic = true;
  full.params = { thickness: o.thickness, alpha0Deg: a0 * R2D, alphaStallDeg: as * R2D, clMax, cdMin, kBucket: kB, clMinDrag, liftSlope: m };
  return full;
}

/**
 * Circular-cylinder section: Cl = 0, constant Cd. The drag coefficient is an argument because it depends on Reynolds
 * number and roughness (smooth cylinder: about 1.2 sub-critical, 0.3 - 0.7 super/trans-critical). The NREL 5-MW
 * definition uses 0.50 (Cylinder1) and 0.35 (Cylinder2) for its root sections.
 * @param {object} [o]
 * @param {number} [o.cd=0.5] drag coefficient [-]
 * @param {string} [o.name='cylinder']
 * @returns {object} polar covering -180..180 deg
 */
export function cylinderPolar(o = {}) {
  const cd = o.cd ?? 0.5;
  need(cd > 0, 'cd', 'must be > 0');
  return createPolar({
    name: o.name || 'cylinder', alphaDeg: [-180, 0, 180], cl: [0, 0, 0], cd: [cd, cd, cd], cm: [0, 0, 0], thickness: 1,
    label: `circular cylinder, Cl = 0, Cd = ${cd}`,
  });
}

/**
 * Screening estimate of the minimum pressure coefficient on a section near its design incidence.
 * Linearised superposition of surface-velocity increments (Abbott and von Doenhoff, ch. 4):
 *   v/V = 1 + kT (t/c) + Cl/4,    Cp,min = 1 - (v/V)^2
 * where kT = 1 is the exact potential-flow value for an ellipse at zero incidence and Cl/4 is the velocity increment
 * of a uniformly loaded (a = 1) mean line. Leading-edge suction peaks away from the shock-free-entry incidence are
 * NOT captured (they can be several times larger): use for first screening of blade-tip cavitation margin only.
 * @param {object} o
 * @param {number} o.cl section lift coefficient [-]
 * @param {number} o.thickness t/c [-]
 * @param {number} [o.kT=1] thickness velocity factor [-] (1.0 ellipse; 1.1 - 1.3 for NACA-type sections)
 * @returns {number} Cp,min [-] (negative)
 */
export function cpMinEstimate(o) {
  need(o && o.thickness >= 0 && o.thickness < 1, 'thickness', 'must be in [0, 1)');
  need(Number.isFinite(o.cl), 'cl', 'must be finite');
  const v = 1 + (o.kT ?? 1) * o.thickness + Math.abs(o.cl) / 4;
  return 1 - v * v;
}

/**
 * Cavitation check for a blade section: sigma = (pAtm + rho g h - pv) / (0.5 rho W^2); inception when sigma < -Cp,min.
 * @param {object} o
 * @param {number} o.cpMin minimum pressure coefficient [-]
 * @param {number} o.W relative speed at the section [m/s]
 * @param {number} o.depth submergence of the section below the free surface [m]
 * @param {number} [o.rho=1025] water density [kg/m3]
 * @param {number} [o.pAtm=101325] atmospheric pressure [Pa]
 * @param {number} [o.pVapour=1700] vapour pressure [Pa] (sea water about 1.2 kPa at 10 degC, 1.7 kPa at 15 degC)
 * @param {number} [o.g=9.80665] gravity [m/s2]
 * @returns {{sigma:number, margin:number, cavitates:boolean}} cavitation number [-], sigma + Cp,min [-], flag
 */
export function cavitationCheck(o) {
  need(o && o.W > 0, 'W', 'must be > 0');
  need(o.depth >= 0, 'depth', 'must be >= 0');
  const rho = o.rho ?? 1025;
  const sigma = ((o.pAtm ?? 101325) + rho * (o.g ?? 9.80665) * o.depth - (o.pVapour ?? 1700)) / (0.5 * rho * o.W * o.W);
  const margin = sigma + o.cpMin;
  return { sigma, margin, cavitates: margin < 0 };
}

/**
 * Resolve lift and drag into rotor-plane components for inflow angle phi (between the relative velocity and the
 * rotor plane): Cn = Cl cos(phi) + Cd sin(phi) (normal to the rotor plane), Ct = Cl sin(phi) - Cd cos(phi)
 * (tangential, driving positive).
 * @param {number} cl lift coefficient [-]
 * @param {number} cd drag coefficient [-]
 * @param {number} phi inflow angle [rad]
 * @returns {{cn:number, ct:number}} normal and tangential coefficients [-]
 */
export function normalTangential(cl, cd, phi) {
  const s = Math.sin(phi);
  const c = Math.cos(phi);
  return { cn: cl * c + cd * s, ct: cl * s - cd * c };
}
