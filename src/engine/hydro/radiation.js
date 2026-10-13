// Frequency-dependent radiation and wave excitation of a surface-piercing hull in deep water at zero speed:
// added mass A(omega), wave damping B(omega) and exciting forces by a three-dimensional panel method with the
// free-surface (wave) Green function.
//
// Method: flat quadrilateral panels on the wetted hull, constant potential per panel, collocation at the panel
// centres, Green's third identity (the 'potential' formulation of bem.js) with
//
//   G(x, xi) = 1/r + 1/r' + 2 K [ F(X, Y) - i pi exp(Y) J0(X) ],   K = omega²/g, X = K R, Y = K (z + zeta) <= 0,
//   F(X, Y) = PV ∫₀^∞ exp(u Y) J0(u X) / (u - 1) du,
//
// r' the distance to the image of the source point in z = 0, R the horizontal distance (time factor exp(+i omega t),
// outgoing waves). The Rankine parts 1/r + 1/r' are integrated in closed form over every panel with the routines
// of bem.js; the wave part is smooth away from the free surface and is taken at the panel centre times the area.
// F follows from the function of one variable F(X, 0) by integrating its differential equation in Y,
//
//   F(X, 0) = -(pi/2) [H0(X) + Y0(X)],   dF/dY = F + 1/d,   d = sqrt(X² + Y²),
//   F(X, Y) = exp(Y) [ F(X, 0) - ∫₀^|Y| exp(t) / sqrt(X² + t²) dt ]
//
// (H0 Struve function, Y0 Bessel function of the second kind; Newman 1985, whose F is 1/d + 2 F here). The
// logarithmic singularity at the origin, -exp(Y) ln(d + |Y|), is split off and the smooth remainder and its X
// derivative are tabulated (X to 100, |Y| to 16, on a grid uniform in sqrt(X) and sqrt(|Y|)) and interpolated with cubic Lagrange polynomials; dF/dY
// needs no table. Beyond |Y| = 16 the wave part is below 2e-7 of the image term and the expansion
// F ~ -sum n! P_n(|Y|/d) / d^(n+1) is used.
//
// Two planes of symmetry (x = 0 and y = 0) are used when the hull has them: only a quarter of the hull carries
// unknowns and the six modes fall into four classes that are solved separately.
//
// Exciting forces come from the Haskind relation, X_i = i omega rho ∫ (phi_0 n_i - phi_i d phi_0/dn) dS with the
// radiation potentials phi_i and the undisturbed wave phi_0, so no diffraction problem is solved and any heading
// costs one integral over the hull.
//
// Irregular frequencies (resonances of the fictitious flow inside the hull, where the integral equation loses
// uniqueness) are suppressed by null-field equations C phi = d at points inside the hull (Schenck 1968), added
// to the surface equations K phi = b as (K + C^H C) phi = b + C^H d: the true solution satisfies both, and the
// added term removes the null space of K at an irregular frequency.
//
// References: Wehausen & Laitone, "Surface waves", Handbuch der Physik IX (1960) §13 (Green function);
// Newman, "Algorithms for the free-surface Green function", J. Eng. Math. 19 (1985) 57–67; Abramowitz & Stegun,
// Handbook of Mathematical Functions (1964) ch. 9 (Bessel functions: ascending series and Hankel expansions) and
// ch. 12 (Struve functions: series, expansion of H0 - Y0); Newman, Marine Hydrodynamics (1977) ch. 6 (radiation
// problem, Haskind relations); Schenck, "Improved integral formulation for acoustic radiation
// problems", J. Acoust. Soc. Am. 44 (1968) 41–58 (interior null-field points); Morino & Kuo, AIAA J. 12 (1974).
//
// Validity: deep water, zero forward speed, linear potential flow (no viscous damping). Constant panels need
// about eight per wavelength; boxRadiation chooses its mesh for wavelengths down to the beam of the box.

import { panelGeometry, panelIntegral } from './bem.js';
import { G } from './_wavehelp.js';

const EULER = 0.5772156649015329;
const X_MAX = 100; const HQ = 0.02; const NQ = Math.round(Math.sqrt(X_MAX) / HQ) + 1; // X = q², q uniform
const S_MAX = 4; const HS = 0.05; const NS = Math.round(S_MAX / HS) + 1; // Y = -s², |Y| <= 16
const SERIES_MAX = 20;

// Bessel functions J0, J1, Y0, Y1 and the Struve function H0 with dH0/dx, by ascending series below x = 20 and
// by the Hankel expansions and the expansion of H0 - Y0 above (Abramowitz & Stegun ch. 9 and ch. 12). Used only
// to build the tables; the tabulated F is checked against direct quadrature of its defining integral in the tests.
function besselStruve(x, out) {
  if (x < SERIES_MAX) {
    const q = 0.25 * x * x; const lg = Math.log(0.5 * x) + EULER;
    let term = 1; let j0 = 1; let dj0 = 0; let s = 0; let ds = 0; let h = 0;
    for (let k = 1; k < 80; k++) {
      // term = (-q)^k / (k!)²; d term/dx = 2 k term / x
      term *= -q / (k * k); h += 1 / k;
      j0 += term; dj0 += (2 * k * term) / x;
      s -= h * term; ds -= (h * 2 * k * term) / x;
      if (Math.abs(term) < 1e-17 && k > q) break;
    }
    const y0 = (2 / Math.PI) * (lg * j0 + s);
    const dy0 = (2 / Math.PI) * (j0 / x + lg * dj0 + ds);
    // H0 = (2/pi) sum (-1)^k x^(2k+1) / ((2k+1)!!)²
    let t = x; let h0 = x; let dh0 = 1;
    for (let k = 1; k < 80; k++) {
      t *= (-x * x) / ((2 * k + 1) * (2 * k + 1));
      h0 += t; dh0 += ((2 * k + 1) * t) / x;
      if (Math.abs(t) < 1e-17 && k > x) break;
    }
    out.j0 = j0; out.j1 = -dj0; out.y0 = y0; out.y1 = -dy0; out.h0 = (2 / Math.PI) * h0; out.dh0 = (2 / Math.PI) * dh0;
    return out;
  }
  const amp = Math.sqrt(2 / (Math.PI * x));
  for (let nu = 0; nu < 2; nu++) {
    // a_k = prod_{j<=k} (4 nu² - (2j - 1)²) / (k! 8^k);  P = sum (-1)^k a_2k / x^2k,  Q = sum (-1)^k a_2k+1 / x^(2k+1)
    let a = 1; let P = 1; let Q = 0; let last = Infinity;
    for (let k = 1; k < 60; k++) {
      a *= (4 * nu * nu - (2 * k - 1) * (2 * k - 1)) / (8 * k * x);
      if (Math.abs(a) > last) break; // asymptotic series: stop at the smallest term
      last = Math.abs(a);
      const sg = (k >> 1) % 2 === 0 ? 1 : -1;
      if (k % 2 === 0) P += sg * a; else Q += sg * a;
      if (last < 1e-17) break;
    }
    const chi = x - (0.5 * nu + 0.25) * Math.PI;
    const jj = amp * (P * Math.cos(chi) - Q * Math.sin(chi)); const yy = amp * (P * Math.sin(chi) + Q * Math.cos(chi));
    if (nu === 0) { out.j0 = jj; out.y0 = yy; } else { out.j1 = jj; out.y1 = yy; }
  }
  // H0 - Y0 = (2/pi) sum (-1)^k ((2k-1)!!)² / x^(2k+1)
  let t = 1 / x; let d = t; let dd = -t / x; let last = Infinity;
  for (let k = 1; k < 60; k++) {
    t *= (-(2 * k - 1) * (2 * k - 1)) / (x * x);
    if (Math.abs(t) > last) break;
    last = Math.abs(t);
    d += t; dd -= ((2 * k + 1) * t) / x;
    if (last < 1e-18) break;
  }
  out.h0 = out.y0 + (2 / Math.PI) * d; out.dh0 = -out.y1 + (2 / Math.PI) * dd;
  return out;
}

// 8-point Gauss–Legendre rule on (-1, 1), computed by Newton iteration on P_8
const GL = (() => {
  const n = 8; const x = new Float64Array(n); const w = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let z = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5)); let dp = 1;
    for (let it = 0; it < 50; it++) {
      let p0 = 1; let p1 = z;
      for (let k = 2; k <= n; k++) { const p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
      dp = (n * (z * p1 - p0)) / (z * z - 1);
      const dz = p1 / dp; z -= dz;
      if (Math.abs(dz) < 1e-15) break;
    }
    x[i] = z; w[i] = 2 / ((1 - z * z) * dp * dp);
  }
  return { x, w };
})();

let TAB = null;

// Tables of the wave part of the Green function: T(X, Y) = F + exp(Y) ln(d + |Y|) and dT/dX on a grid uniform
// in q = sqrt(X) and s = sqrt(|Y|) (in these variables the remainder is smooth at the origin too), and J0, J1 on
// the X grid. Rows are independent and are filled up to the largest X asked for.
function tables(iqNeed) {
  if (!TAB) TAB = { T: new Float64Array(NQ * NS), TX: new Float64Array(NQ * NS), J0: new Float64Array(NQ), J1: new Float64Array(NQ), rows: 0 };
  if (iqNeed < TAB.rows) return TAB;
  const { T, TX, J0, J1 } = TAB;
  const bs = {};
  const to = Math.min(NQ, iqNeed + 26);
  for (let ix = TAB.rows; ix < to; ix++) {
    const X = (ix * HQ) ** 2;
    let c; let dc; // c = F(X, 0) + ln X
    if (ix === 0) { c = Math.LN2 - EULER; dc = -1; J0[0] = 1; J1[0] = 0; } else {
      besselStruve(X, bs);
      c = -0.5 * Math.PI * (bs.h0 + bs.y0) + Math.log(X);
      dc = -0.5 * Math.PI * (bs.dh0 - bs.y1) + 1 / X;
      J0[ix] = bs.j0; J1[ix] = bs.j1;
    }
    // I2 = ∫₀^a (e^t - 1)/sqrt(X² + t²) dt and I3 = X ∫₀^a (e^t - 1)/(X² + t²)^(3/2) dt, accumulated over the Y
    // grid; with t = X sinh v the integrands are (e^t - 1) and (e^t - 1)/(X cosh² v)
    let i2 = 0; let i3 = 0;
    T[ix * NS] = c; TX[ix * NS] = dc;
    for (let is = 1; is < NS; is++) {
      const a0 = ((is - 1) * HS) ** 2; const a1 = (is * HS) ** 2;
      if (ix === 0) {
        const hm = 0.5 * (a1 - a0); const mid = 0.5 * (a1 + a0);
        for (let q = 0; q < 8; q++) { const t = mid + hm * GL.x[q]; i2 += (GL.w[q] * hm * Math.expm1(t)) / t; }
        i3 = 1;
      } else {
        const v0 = Math.asinh(a0 / X); const v1 = Math.asinh(a1 / X);
        const hm = 0.5 * (v1 - v0); const mid = 0.5 * (v1 + v0);
        for (let q = 0; q < 8; q++) {
          const v = mid + hm * GL.x[q]; const e = Math.expm1(X * Math.sinh(v)); const ch = Math.cosh(v);
          i2 += GL.w[q] * hm * e; i3 += (GL.w[q] * hm * e) / (X * ch * ch);
        }
      }
      const ey = Math.exp(-a1);
      T[ix * NS + is] = ey * (c - i2); TX[ix * NS + is] = ey * (dc + i3);
    }
  }
  TAB.rows = to;
  return TAB;
}

const _wx = new Float64Array(4); const _ws = new Float64Array(4);
// cubic Lagrange weights for a point at fraction u of the interval between nodes 1 and 2 of four equidistant nodes
function lagrange(u, w) {
  const um = u - 1; const up = u + 1; const u2 = u - 2;
  w[0] = (-u * um * u2) / 6; w[1] = (up * um * u2) / 2; w[2] = (-up * u * u2) / 2; w[3] = (up * u * um) / 6;
}

/**
 * Wave part of the deep-water free-surface Green function, G = 1/r + 1/r' + 2 K [F - i pi exp(Y) J0(X)].
 * Returns F = PV ∫₀^∞ exp(uY) J0(uX)/(u - 1) du with its derivatives, and the Bessel functions J0(X), J1(X).
 * Reference: Newman (1985) J. Eng. Math. 19, 57–67; Wehausen & Laitone (1960) eq. 13.17.
 * Validity: 0 <= X <= 100, Y <= 0; interpolation error below 1e-4 (absolute) away from the origin, where F
 * behaves as -exp(Y) (ln((d + |Y|)/2) + gamma). At X = Y = 0 the function is singular (returns ±Infinity).
 * @param {number} X K times the horizontal distance, -
 * @param {number} Y K (z + zeta), <= 0
 * @param {{F:number, FX:number, FY:number, J0:number, J1:number}} [out] receives the result (reused when given)
 * @returns {{F:number, FX:number, FY:number, J0:number, J1:number}} F, dF/dX, dF/dY, J0(X), J1(X), all -
 */
export function waveGreenFunction(X, Y, out = {}) {
  if (!(X >= 0) || !(Y <= 0)) throw new RangeError('waveGreenFunction: X must be >= 0 and Y <= 0');
  if (X > X_MAX) throw new RangeError(`waveGreenFunction: X = ${X} is beyond the tabulated range (${X_MAX})`);
  const a = -Y; const d = Math.sqrt(X * X + a * a);
  const qx = Math.sqrt(X) / HQ;
  let ix = Math.floor(qx) - 1;
  if (ix < 0) ix = 0; else if (ix > NQ - 4) ix = NQ - 4;
  const tb = TAB !== null && ix + 4 <= TAB.rows ? TAB : tables(ix + 4);
  lagrange(qx - ix - 1, _wx);
  out.J0 = _wx[0] * tb.J0[ix] + _wx[1] * tb.J0[ix + 1] + _wx[2] * tb.J0[ix + 2] + _wx[3] * tb.J0[ix + 3];
  out.J1 = _wx[0] * tb.J1[ix] + _wx[1] * tb.J1[ix + 1] + _wx[2] * tb.J1[ix + 2] + _wx[3] * tb.J1[ix + 3];
  const s = Math.sqrt(a);
  if (s >= S_MAX) {
    // F ~ -(1/d + a/d³ + (2a² - X²)/d⁵ + (6a³ - 9aX²)/d⁷): the first terms of -sum n! P_n(a/d)/d^(n+1)
    const d2 = d * d; const d3 = d2 * d; const d5 = d3 * d2; const d7 = d5 * d2;
    const p2 = 2 * a * a - X * X; const p3 = 6 * a * a * a - 9 * a * X * X;
    out.F = -(1 / d + a / d3 + p2 / d5 + p3 / d7);
    out.FX = X / d3 + (3 * a * X) / d5 + (2 * X) / d5 + (5 * X * p2) / d7 + (18 * a * X) / d7 + (7 * X * p3) / (d7 * d2);
    out.FY = out.F + 1 / d;
    return out;
  }
  let is = Math.floor(s / HS) - 1;
  if (is < 0) is = 0; else if (is > NS - 4) is = NS - 4;
  lagrange(s / HS - is - 1, _ws);
  let t = 0; let tx = 0;
  for (let p = 0; p < 4; p++) {
    const o = (ix + p) * NS + is;
    t += _wx[p] * (_ws[0] * tb.T[o] + _ws[1] * tb.T[o + 1] + _ws[2] * tb.T[o + 2] + _ws[3] * tb.T[o + 3]);
    tx += _wx[p] * (_ws[0] * tb.TX[o] + _ws[1] * tb.TX[o + 1] + _ws[2] * tb.TX[o + 2] + _ws[3] * tb.TX[o + 3]);
  }
  const ey = Math.exp(Y);
  out.F = t - ey * Math.log(d + a);
  out.FX = d > 0 ? tx - (ey * X) / (d * (d + a)) : 0;
  out.FY = out.F + 1 / d;
  return out;
}

// mirror copies of a field point in the planes x = 0 and y = 0: [sign x, sign y]
const MIRRORS = [[1, 1], [-1, 1], [1, -1], [-1, -1]];
// parity of the six radiation potentials under x -> -x and y -> -y (surge … yaw)
const PARITY = [[-1, 1], [1, -1], [1, 1], [1, -1], [-1, 1], [-1, -1]];
const G2 = 1 / Math.sqrt(3);

/**
 * Prepares a set of quadrilateral hull panels for radiationSolve: centroids, normals, areas, the two triangles
 * of every quadrilateral (for the closed-form Rankine integrals of bem.js) and a 2 × 2 Gauss rule per panel.
 * Reference: geometry only.
 * @param {ArrayLike<number>} quads 12 numbers per panel: four vertices in order, counter-clockwise when seen
 *   from the fluid (the normal (v1 - v0) × (v3 - v0) points into the fluid), m; a triangle repeats its last
 *   vertex. With symmetry only the part of the hull in x >= 0 and/or y >= 0 is given.
 * @returns {object} panel set: { n, quads, centroids, normals, areas, tri, triGeo, triStart, gauss }
 */
export function hullPanels(quads) {
  if (!quads || quads.length < 12 || quads.length % 12 !== 0) throw new RangeError('quads must hold four vertices (12 numbers) per panel');
  const n = quads.length / 12;
  const centroids = new Float64Array(3 * n); const normals = new Float64Array(3 * n); const areas = new Float64Array(n);
  const tv = []; const tn = []; const triStart = new Int32Array(n + 1);
  const gauss = new Float64Array(16 * n); // per panel 4 points: x, y, z, weight (m²)
  for (let p = 0; p < n; p++) {
    const q = 12 * p;
    let ax = 0; let ay = 0; let az = 0; let cx = 0; let cy = 0; let cz = 0; let at = 0;
    const tris = [];
    for (const [i1, i2] of [[1, 2], [2, 3]]) {
      const ux = quads[q + 3 * i1] - quads[q]; const uy = quads[q + 3 * i1 + 1] - quads[q + 1]; const uz = quads[q + 3 * i1 + 2] - quads[q + 2];
      const vx = quads[q + 3 * i2] - quads[q]; const vy = quads[q + 3 * i2 + 1] - quads[q + 1]; const vz = quads[q + 3 * i2 + 2] - quads[q + 2];
      const nx = uy * vz - uz * vy; const ny = uz * vx - ux * vz; const nz = ux * vy - uy * vx;
      const a2 = Math.hypot(nx, ny, nz);
      ax += nx; ay += ny; az += nz;
      if (!(a2 > 0)) continue;
      at += a2 / 2;
      cx += (a2 / 2) * (quads[q] + quads[q + 3 * i1] + quads[q + 3 * i2]) / 3;
      cy += (a2 / 2) * (quads[q + 1] + quads[q + 3 * i1 + 1] + quads[q + 3 * i2 + 1]) / 3;
      cz += (a2 / 2) * (quads[q + 2] + quads[q + 3 * i1 + 2] + quads[q + 3 * i2 + 2]) / 3;
      tris.push([i1, i2, nx / a2, ny / a2, nz / a2]);
    }
    const al = Math.hypot(ax, ay, az);
    if (!(at > 0) || !(al > 0)) throw new RangeError(`panel ${p} has no area`);
    normals[3 * p] = ax / al; normals[3 * p + 1] = ay / al; normals[3 * p + 2] = az / al;
    centroids[3 * p] = cx / at; centroids[3 * p + 1] = cy / at; centroids[3 * p + 2] = cz / at; areas[p] = at;
    if (centroids[3 * p + 2] >= 0) throw new RangeError('panels must lie below the still-water plane z = 0');
    for (const [i1, i2, nx, ny, nz] of tris) {
      tv.push(quads[q], quads[q + 1], quads[q + 2]);
      for (const i of [i1, i2]) tv.push(quads[q + 3 * i], quads[q + 3 * i + 1], quads[q + 3 * i + 2]);
      tn.push(nx, ny, nz);
    }
    triStart[p + 1] = tn.length / 3;
    // bilinear map of the unit square, Jacobian |dP/du × dP/dv|
    let k = 0;
    for (const u of [-G2, G2]) for (const v of [-G2, G2]) {
      const N = [(1 - u) * (1 - v) / 4, (1 + u) * (1 - v) / 4, (1 + u) * (1 + v) / 4, (1 - u) * (1 + v) / 4];
      const dU = [-(1 - v) / 4, (1 - v) / 4, (1 + v) / 4, -(1 + v) / 4]; const dV = [-(1 - u) / 4, -(1 + u) / 4, (1 + u) / 4, (1 - u) / 4];
      const P = [0, 0, 0]; const Pu = [0, 0, 0]; const Pv = [0, 0, 0];
      for (let c = 0; c < 4; c++) for (let d = 0; d < 3; d++) { const val = quads[q + 3 * c + d]; P[d] += N[c] * val; Pu[d] += dU[c] * val; Pv[d] += dV[c] * val; }
      const jac = Math.hypot(Pu[1] * Pv[2] - Pu[2] * Pv[1], Pu[2] * Pv[0] - Pu[0] * Pv[2], Pu[0] * Pv[1] - Pu[1] * Pv[0]);
      gauss.set([P[0], P[1], P[2], jac], 16 * p + 4 * k); k++;
    }
  }
  const tri = { n: tn.length / 3, vertices: Float64Array.from(tv), normals: Float64Array.from(tn) };
  return { n, quads: Float64Array.from(quads), centroids, normals, areas, tri, triGeo: panelGeometry(tri), triStart, gauss };
}

// Complex Gaussian elimination with partial pivoting; a: n×n interleaved (re, im) row-major, b: nrhs columns
// stored as b[r][2 i]. Overwrites both; the solution is left in b.
function csolveDense(a, n, b) {
  const n2 = 2 * n; const nb = b.length;
  for (let c = 0; c < n; c++) {
    let p = c; let best = -1;
    for (let r = c; r < n; r++) { const re = a[r * n2 + 2 * c]; const im = a[r * n2 + 2 * c + 1]; const m = re * re + im * im; if (m > best) { best = m; p = r; } }
    if (!(best > 0)) throw new Error('radiationSolve: the influence matrix is singular (duplicate or degenerate panels?)');
    if (p !== c) {
      for (let j = 2 * c; j < n2; j++) { const t = a[c * n2 + j]; a[c * n2 + j] = a[p * n2 + j]; a[p * n2 + j] = t; }
      for (let q = 0; q < nb; q++) { const col = b[q]; let t = col[2 * c]; col[2 * c] = col[2 * p]; col[2 * p] = t; t = col[2 * c + 1]; col[2 * c + 1] = col[2 * p + 1]; col[2 * p + 1] = t; }
    }
    const pr = a[c * n2 + 2 * c]; const pi = a[c * n2 + 2 * c + 1]; const den = pr * pr + pi * pi;
    const rc = c * n2;
    for (let r = c + 1; r < n; r++) {
      const rr = r * n2;
      const ar = a[rr + 2 * c]; const ai = a[rr + 2 * c + 1];
      if (ar === 0 && ai === 0) continue;
      const fr = (ar * pr + ai * pi) / den; const fi = (ai * pr - ar * pi) / den;
      for (let j = 2 * c + 2; j < n2; j += 2) {
        const xr = a[rc + j]; const xi = a[rc + j + 1];
        a[rr + j] -= fr * xr - fi * xi; a[rr + j + 1] -= fr * xi + fi * xr;
      }
      for (let q = 0; q < nb; q++) { const col = b[q]; const xr = col[2 * c]; const xi = col[2 * c + 1]; col[2 * r] -= fr * xr - fi * xi; col[2 * r + 1] -= fr * xi + fi * xr; }
    }
  }
  for (const col of b) {
    for (let r = n - 1; r >= 0; r--) {
      let sr = col[2 * r]; let si = col[2 * r + 1];
      const rr = r * n2;
      for (let j = r + 1; j < n; j++) { const xr = col[2 * j]; const xi = col[2 * j + 1]; sr -= a[rr + 2 * j] * xr - a[rr + 2 * j + 1] * xi; si -= a[rr + 2 * j] * xi + a[rr + 2 * j + 1] * xr; }
      const pr = a[rr + 2 * r]; const pi = a[rr + 2 * r + 1]; const den = pr * pr + pi * pi;
      col[2 * r] = (sr * pr + si * pi) / den; col[2 * r + 1] = (si * pr - sr * pi) / den;
    }
  }
}

/**
 * Frequency-independent part of the panel problem: Rankine influence matrices (1/r and its image in z = 0,
 * closed form) and the geometry the wave part needs, for the field points of every mirror copy.
 * Reference: see the file header. Internal to radiationSolve; exported for reuse over many frequencies.
 * @param {object} o
 * @param {object} o.panels from hullPanels
 * @param {{x?:boolean, y?:boolean}} [o.symmetry] planes of symmetry x = 0 and/or y = 0 (the panels then cover
 *   only x >= 0 and/or y >= 0)
 * @param {number[]} [o.ref=[0,0,0]] reference point of rotations and moments, m (on a plane of symmetry)
 * @param {number[][]} [o.interiorPoints] points inside the hull, below z = 0 and off the planes of symmetry,
 *   for the null-field equations that remove irregular frequencies, m
 * @returns {object} setup for radiationSolve
 */
export function radiationSetup({ panels, symmetry = {}, ref = [0, 0, 0], interiorPoints = [] }) {
  if (!panels || !panels.tri) throw new RangeError('panels must come from hullPanels');
  const symX = !!symmetry.x; const symY = !!symmetry.y;
  if ((symX && ref[0] !== 0) || (symY && ref[1] !== 0)) throw new RangeError('ref must lie on the planes of symmetry');
  const mirrors = MIRRORS.filter(([sx, sy]) => (symX || sx === 1) && (symY || sy === 1));
  const nm = mirrors.length; const n = panels.n; const np = interiorPoints.length; const nr = n + np;
  const C = panels.centroids; const N = panels.normals; const Ar = panels.areas;
  for (let i = 0; i < n; i++) if ((symX && C[3 * i] <= 0) || (symY && C[3 * i + 1] <= 0)) throw new RangeError('with symmetry the panels must lie in x > 0 / y > 0');
  // S: ∫ (1/r + 1/r') dS, D: minus the normal-dipole integral (signed solid angle), rows = field points
  const S = []; const D = []; const R = []; const dRn = [];
  const zs = new Float64Array(nr * n);
  const r = new Float64Array(5);
  for (let m = 0; m < nm; m++) {
    const Sm = new Float64Array(nr * n); const Dm = new Float64Array(nr * n); const Rm = new Float64Array(nr * n); const dm = new Float64Array(nr * n);
    const [sx, sy] = mirrors[m];
    for (let i = 0; i < nr; i++) {
      const src = i < n ? C : interiorPoints[i - n]; const o = i < n ? 3 * i : 0;
      const px = sx * src[o]; const py = sy * src[o + 1]; const pz = src[o + 2];
      if (!(pz < 0)) throw new RangeError('interior points must lie below z = 0');
      for (let j = 0; j < n; j++) {
        let pot = 0; let ker = 0;
        const self = m === 0 && i === j;
        for (let t = panels.triStart[j]; t < panels.triStart[j + 1]; t++) {
          panelIntegral(panels.tri, panels.triGeo, t, px, py, pz, self, r); pot += r[0]; ker += r[4];
          panelIntegral(panels.tri, panels.triGeo, t, px, py, -pz, false, r); pot += r[0]; ker += r[4];
        }
        Sm[i * n + j] = pot; Dm[i * n + j] = ker;
        const dx = C[3 * j] - px; const dy = C[3 * j + 1] - py; const rr = Math.hypot(dx, dy);
        Rm[i * n + j] = rr; dm[i * n + j] = rr > 0 ? (dx * N[3 * j] + dy * N[3 * j + 1]) / rr : 0;
        if (m === 0) zs[i * n + j] = pz + C[3 * j + 2];
      }
    }
    S.push(Sm); D.push(Dm); R.push(Rm); dRn.push(dm);
  }
  let rMax = 0;
  for (const Rm of R) for (let i = 0; i < Rm.length; i++) if (Rm[i] > rMax) rMax = Rm[i];
  // generalised normals at the panel centres
  const gn = [];
  for (let k = 0; k < 6; k++) gn.push(new Float64Array(n));
  for (let i = 0; i < n; i++) {
    const rx = C[3 * i] - ref[0]; const ry = C[3 * i + 1] - ref[1]; const rz = C[3 * i + 2] - ref[2];
    const nx = N[3 * i]; const ny = N[3 * i + 1]; const nz = N[3 * i + 2];
    gn[0][i] = nx; gn[1][i] = ny; gn[2][i] = nz;
    gn[3][i] = ry * nz - rz * ny; gn[4][i] = rz * nx - rx * nz; gn[5][i] = rx * ny - ry * nx;
  }
  // classes of modes that share a system matrix
  const classes = [];
  if (symX || symY) {
    for (const [ex, ey] of MIRRORS) {
      if ((!symX && ex === -1) || (!symY && ey === -1)) continue;
      const modes = [];
      for (let k = 0; k < 6; k++) if ((!symX || PARITY[k][0] === ex) && (!symY || PARITY[k][1] === ey)) modes.push(k);
      classes.push({ signs: mirrors.map(([sx, sy]) => (sx === -1 ? ex : 1) * (sy === -1 ? ey : 1)), modes });
    }
  } else classes.push({ signs: [1], modes: [0, 1, 2, 3, 4, 5] });
  // Rankine part of every class: system matrix 2 pi I + sum sign D, and right-hand sides -(sum sign S) n_k
  for (const cl of classes) {
    cl.K0 = new Float64Array(nr * n); cl.rhs0 = cl.modes.map(() => new Float64Array(nr));
    for (let i = 0; i < nr; i++) for (let j = 0; j < n; j++) {
      const q = i * n + j;
      let d = i === j ? 2 * Math.PI : 0; let sv = 0;
      for (let m = 0; m < nm; m++) { d += cl.signs[m] * D[m][q]; sv += cl.signs[m] * S[m][q]; }
      cl.K0[q] = d;
      for (let c = 0; c < cl.modes.length; c++) cl.rhs0[c][i] -= sv * gn[cl.modes[c]][j];
    }
  }
  return { panels, mirrors, classes, n, np, nr, R, dRn, zs, gn, ref: Array.from(ref), rMax, areas: Ar };
}

const _g = { F: 0, FX: 0, FY: 0, J0: 0, J1: 0 };

/**
 * Radiation problem at one frequency: added mass, wave damping and the radiation potentials on the panels.
 * A_jk = -rho Re ∫ phi_k n_j dS, B_jk = rho omega Im ∫ phi_k n_j dS for unit velocity in mode k (exp(+i omega t)).
 * The raw matrices are not exactly symmetric (discretisation error); the symmetric parts are returned and
 * `symmetryError` reports the largest asymmetry of the added mass that was removed.
 * Reference: Newman (1977) §6.15–6.19; see the file header.
 * Validity: K rMax <= 100 with K = omega²/g and rMax the largest horizontal distance between two panels.
 * @param {object} o
 * @param {object} o.setup from radiationSetup
 * @param {number} o.omega rad/s
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @param {number[]} [o.modes] modes to solve (0 surge … 5 yaw; default all)
 * @param {boolean} [o.nullField=true] use the interior points of the setup at this frequency
 * @returns {{omega:number, A:number[][], B:number[][], phi:(Float64Array|null)[], symmetryError:number}} A in
 *   kg, kg·m, kg·m²; B in N·s/m, N·s, N·m·s (about `ref`); phi per mode: complex potential per unit velocity on
 *   the panels (re, im interleaved), m²/s per m/s or per rad/s
 */
export function radiationSolve({ setup, omega, rho = 1025, g = G, modes, nullField = true }) {
  if (!(omega > 0)) throw new RangeError('omega must be > 0');
  const { n, mirrors, classes, R, dRn, zs, gn, areas, panels } = setup;
  const K = (omega * omega) / g;
  if (K * setup.rMax > X_MAX) throw new RangeError('radiationSolve: frequency too high for the size of the hull (K × largest panel distance > 100)');
  const nr = nullField ? setup.nr : n; const nm = mirrors.length;
  const N = panels.normals;
  const want = modes ? new Set(modes) : null;
  // classes to solve, each with its complex system matrix (rows = field points) and right-hand sides
  const act = [];
  for (const cl of classes) {
    const idx = []; cl.modes.forEach((k, c) => { if (!want || want.has(k)) idx.push(c); });
    if (idx.length === 0) continue;
    // the system matrix is work space kept with the class; the right-hand sides become the returned potentials
    if (!cl.work) cl.work = new Float64Array(2 * setup.nr * n);
    const kc = cl.work;
    for (let q = 0; q < nr * n; q++) { kc[2 * q] = cl.K0[q]; kc[2 * q + 1] = 0; }
    const rhs = idx.map((c) => { const b = new Float64Array(2 * nr); for (let i = 0; i < nr; i++) b[2 * i] = cl.rhs0[c][i]; return b; });
    act.push({ cl, kc, rhs, gns: idx.map((c) => gn[cl.modes[c]]), ms: idx.map((c) => cl.modes[c]) });
  }
  // wave part: S += A_j Gw, K += -A_j dGw/dn at the source panel. Gw depends on the pair of points only through
  // R and z + zeta, so on the hull it is evaluated once for (i, j) and (j, i).
  const twoK = 2 * K; const twoK2 = 2 * K * K; const na = act.length;
  const add = (i, j, q, m, F, FX, FY, eJ0, eJ1) => {
    const a = areas[j]; const nz = N[3 * j + 2]; const dn = dRn[m][q];
    const sr = a * twoK * F; const si = -a * twoK * eJ0;
    const kr = -a * twoK2 * (FX * dn + FY * nz); const ki = -a * twoK2 * (eJ1 * dn - eJ0 * nz);
    for (let c = 0; c < na; c++) {
      const ac = act[c]; const sg = ac.cl.signs[m];
      ac.kc[2 * q] += sg * kr; ac.kc[2 * q + 1] += sg * ki;
      for (let r = 0; r < ac.rhs.length; r++) { const v = sg * ac.gns[r][j]; ac.rhs[r][2 * i] -= sr * v; ac.rhs[r][2 * i + 1] -= si * v; }
    }
  };
  for (let m = 0; m < nm; m++) {
    const Rm = R[m];
    for (let i = 0; i < nr; i++) {
      for (let j = i < n ? i : 0; j < n; j++) {
        const q = i * n + j;
        const Y = K * zs[q];
        waveGreenFunction(K * Rm[q], Y, _g);
        const e = Math.PI * Math.exp(Y);
        add(i, j, q, m, _g.F, _g.FX, _g.FY, e * _g.J0, e * _g.J1);
        if (i < n && j !== i) add(j, i, j * n + i, m, _g.F, _g.FX, _g.FY, e * _g.J0, e * _g.J1);
      }
    }
  }
  const phi = [null, null, null, null, null, null];
  for (const { kc, rhs, ms } of act) {
    if (nr > n) {
      // null-field equations C phi = d added to the surface equations with unit weight: (K + C^H C) phi = b + C^H d
      for (let p = n; p < nr; p++) {
        const co = 2 * p * n;
        for (let i = 0; i < n; i++) {
          const ar = kc[co + 2 * i]; const ai = -kc[co + 2 * i + 1];
          const io = 2 * i * n;
          for (let q = 0; q < 2 * n; q += 2) { const br = kc[co + q]; const bi = kc[co + q + 1]; kc[io + q] += ar * br - ai * bi; kc[io + q + 1] += ar * bi + ai * br; }
          for (const col of rhs) { const br = col[2 * p]; const bi = col[2 * p + 1]; col[2 * i] += ar * br - ai * bi; col[2 * i + 1] += ar * bi + ai * br; }
        }
      }
    }
    csolveDense(kc, n, rhs);
    for (let c = 0; c < ms.length; c++) phi[ms[c]] = nr === n ? rhs[c] : rhs[c].slice(0, 2 * n);
  }
  const rawA = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]); const rawB = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  for (const cl of classes) {
    for (const k of cl.modes) {
      if (!phi[k]) continue;
      for (const j of cl.modes) {
        let re = 0; let im = 0;
        for (let i = 0; i < n; i++) { const w = gn[j][i] * areas[i]; re += phi[k][2 * i] * w; im += phi[k][2 * i + 1] * w; }
        rawA[j][k] = -rho * nm * re; rawB[j][k] = rho * omega * nm * im;
      }
    }
  }
  const A = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]); const B = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  let symmetryError = 0;
  for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
    const both = phi[i] && phi[j];
    A[i][j] = both ? 0.5 * (rawA[i][j] + rawA[j][i]) : rawA[i][j]; B[i][j] = both ? 0.5 * (rawB[i][j] + rawB[j][i]) : rawB[i][j];
    if (both && j > i) { const sc = Math.sqrt(Math.abs(rawA[i][i] * rawA[j][j])); if (sc > 0) symmetryError = Math.max(symmetryError, Math.abs(rawA[i][j] - rawA[j][i]) / sc); }
  }
  return { omega, A, B, phi, symmetryError };
}

/**
 * Wave exciting forces from the Haskind relation, X_i = i omega rho ∫ (phi_0 n_i - phi_i d phi_0/dn) dS, with the
 * undisturbed deep-water wave phi_0 = (i g/omega) exp(K z - i K (x cos beta + y sin beta)) of unit amplitude,
 * integrated with a 2 × 2 Gauss rule on every panel and its mirror copies.
 * Reference: Newman (1962) J. Ship Research 6(3), 10–17; Newman (1977) §6.18–6.19.
 * @param {object} o
 * @param {object} o.setup from radiationSetup
 * @param {(ArrayLike<number>|null)[]} o.phi radiation potentials per mode from radiationSolve
 * @param {number} o.omega wave frequency, rad/s
 * @param {number} [o.headingDeg=0] direction the waves travel towards, from +x towards +y, deg
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{re:number[], im:number[]}} force per unit wave amplitude, N/m and N·m/m about `ref`; a positive
 *   phase is a lead on the crest at the origin; modes without a potential are 0
 */
export function haskindExcitation({ setup, phi, omega, headingDeg = 0, rho = 1025, g = G }) {
  const { panels, mirrors, n, ref } = setup;
  const K = (omega * omega) / g;
  const cb = Math.cos((headingDeg * Math.PI) / 180); const sb = Math.sin((headingDeg * Math.PI) / 180);
  const N = panels.normals; const gp = panels.gauss;
  // I_i = ∫ (phi_0 n_i - phi_i d phi_0/dn) dS / (i g/omega)
  const Ir = [0, 0, 0, 0, 0, 0]; const Ii = [0, 0, 0, 0, 0, 0];
  const nk = [0, 0, 0, 0, 0, 0];
  for (let p = 0; p < n; p++) {
    for (const [sx, sy] of mirrors) {
      const nx = sx * N[3 * p]; const ny = sy * N[3 * p + 1]; const nz = N[3 * p + 2];
      // d phi_0/dn = phi_0 (K nz - i K (nx cb + ny sb))
      const dr = K * nz; const di = -K * (nx * cb + ny * sb);
      for (let q = 0; q < 4; q++) {
        const o = 16 * p + 4 * q;
        const x = sx * gp[o]; const y = sy * gp[o + 1]; const z = gp[o + 2]; const w = gp[o + 3];
        const amp = w * Math.exp(K * z); const ph = -K * (x * cb + y * sb);
        const er = amp * Math.cos(ph); const ei = amp * Math.sin(ph);
        const rx = x - ref[0]; const ry = y - ref[1]; const rz = z - ref[2];
        nk[0] = nx; nk[1] = ny; nk[2] = nz; nk[3] = ry * nz - rz * ny; nk[4] = rz * nx - rx * nz; nk[5] = rx * ny - ry * nx;
        for (let k = 0; k < 6; k++) {
          if (!phi[k]) continue;
          // phi_k on the mirror copy = parity × phi_k on the panel
          const par = (sx === -1 ? PARITY[k][0] : 1) * (sy === -1 ? PARITY[k][1] : 1);
          const fr = par * phi[k][2 * p]; const fi = par * phi[k][2 * p + 1];
          const tr = nk[k] - (fr * dr - fi * di); const ti = -(fr * di + fi * dr);
          Ir[k] += er * tr - ei * ti; Ii[k] += er * ti + ei * tr;
        }
      }
    }
  }
  // X = i omega rho (i g/omega) I = -rho g I
  return { re: Ir.map((v) => -rho * g * v), im: Ii.map((v) => -rho * g * v) };
}

/**
 * Wave damping from the energy radiated to infinity (Haskind relation, deep water),
 *   B_ij = omega³ / (4 pi rho g³) ∫₀^{2π} Re{ X_i(beta) X_j*(beta) } d beta,
 * with the exciting forces of haskindExcitation. The matrix is symmetric and positive semi-definite by
 * construction, which the pressure integral rho omega Im ∫ phi_j n_i dS of radiationSolve is not: for roll about
 * the waterline, a small difference of the contributions of the sides and the bottom, that integral converges
 * slowly with the mesh, while the far-field form agrees with it within 1–2 % in the other modes.
 * Midpoint rule over the headings (over one quadrant when the hull has both planes of symmetry).
 * Reference: Newman (1962) J. Ship Research 6(3), 10–17; Newman (1977) §6.19.
 * @param {object} o
 * @param {object} o.setup from radiationSetup
 * @param {(ArrayLike<number>|null)[]} o.phi radiation potentials per mode from radiationSolve
 * @param {number} o.omega rad/s
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @param {number} [o.nHeadings] headings in the integration range (default: enough for the phase variation
 *   over the hull, 0.6 K rMax + 6 per quadrant)
 * @returns {number[][]} 6×6, N·s/m, N·s, N·m·s about `ref`; rows and columns of modes without a potential are 0
 */
export function haskindDamping({ setup, phi, omega, rho = 1025, g = G, nHeadings }) {
  const K = (omega * omega) / g;
  const quarter = setup.mirrors.length === 4;
  const span = quarter ? 90 : 360;
  const n = nHeadings === undefined ? (quarter ? 1 : 4) * (Math.ceil(0.6 * K * setup.rMax) + 6) : nHeadings;
  if (!(n >= 2)) throw new RangeError('nHeadings must be >= 2');
  const B = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  for (let h = 0; h < n; h++) {
    const x = haskindExcitation({ setup, phi, omega, headingDeg: (span * (h + 0.5)) / n, rho, g });
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) B[i][j] += x.re[i] * x.re[j] + x.im[i] * x.im[j];
  }
  const c = ((omega ** 3) / (4 * Math.PI * rho * g ** 3)) * ((2 * Math.PI) / n);
  for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
    // over a quadrant only the products that are even in both wave-number components survive the full circle
    const even = !quarter || (PARITY[i][0] === PARITY[j][0] && PARITY[i][1] === PARITY[j][1]);
    B[i][j] = even ? c * B[i][j] : 0;
  }
  return B;
}

/**
 * Wetted panels of a rectangular box L × B × T (quarter x >= 0, y >= 0, for radiationSetup with both planes of
 * symmetry): bottom nx × ny, side nx × nz, end ny × nz, with the panel edges drawn towards the keel edges and
 * the corners, where the potential varies most, by the mapping s -> (1 - c) s + c sin(pi s/2).
 * Reference: geometry only.
 * @param {{L:number, B:number, T:number, nx:number, ny:number, nz:number, grading?:number}} o m; panel counts
 *   over the half length, half beam and draught; grading c in 0..1 (0 = uniform, 1 = cosine spacing, the default)
 * @returns {Float64Array} 12 numbers per panel for hullPanels
 */
export function boxQuarterPanels({ L, B, T, nx, ny, nz, grading = 1 }) {
  if (!(L > 0 && B > 0 && T > 0)) throw new RangeError('L, B and T must be > 0');
  if (!(nx >= 1 && ny >= 1 && nz >= 1)) throw new RangeError('nx, ny and nz must be >= 1');
  const map = (s) => (1 - grading) * s + grading * Math.sin(0.5 * Math.PI * s);
  const xs = Array.from({ length: nx + 1 }, (_, i) => 0.5 * L * map(i / nx));
  const ys = Array.from({ length: ny + 1 }, (_, i) => 0.5 * B * map(i / ny));
  const zsv = Array.from({ length: nz + 1 }, (_, i) => -T * map(i / nz));
  const q = [];
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) q.push(xs[i], ys[j], -T, xs[i], ys[j + 1], -T, xs[i + 1], ys[j + 1], -T, xs[i + 1], ys[j], -T);
  for (let i = 0; i < nx; i++) for (let k = 0; k < nz; k++) q.push(xs[i], 0.5 * B, zsv[k + 1], xs[i], 0.5 * B, zsv[k], xs[i + 1], 0.5 * B, zsv[k], xs[i + 1], 0.5 * B, zsv[k + 1]);
  for (let j = 0; j < ny; j++) for (let k = 0; k < nz; k++) q.push(0.5 * L, ys[j], zsv[k + 1], 0.5 * L, ys[j + 1], zsv[k + 1], 0.5 * L, ys[j + 1], zsv[k], 0.5 * L, ys[j], zsv[k]);
  return Float64Array.from(q);
}

const _boxCache = new Map();
const BOX_CACHE_MAX = 8;
const _lw = new Float64Array(4);

/**
 * Hydrodynamic model of a rectangular box L × B × T floating in deep water at zero speed: radiation problems in
 * surge, sway, heave, roll and pitch solved by radiationSolve on a grid of frequencies and kept (potentials
 * included) for interpolation, so that coefficients and exciting forces at any frequency and heading cost one
 * pass over the panels. Results are cached per geometry.
 * Mesh: a quarter of the hull with both planes of symmetry; ny panels over the half beam, nx over the half
 * length (at most 14; longer than wide on long hulls), nz over the draught, drawn towards
 * the edges (boxQuarterPanels). Frequencies: `nOmega` uniform steps up to the wavelength `lambdaMin` (default:
 * the beam; raised if K × hull diagonal would leave the tabulated range of the Green function).
 * Damping: from the radiated energy (haskindDamping), symmetric and positive semi-definite.
 * Irregular frequencies: the lowest one of a box is known in closed form, K = kappa coth(kappa T) with
 * kappa = pi sqrt(1/L² + 1/B²) (first mode of the interior problem with phi = 0 on the hull); from 0.9 times that
 * frequency the null-field equations at four interior points are added.
 * Reference: see the file header. The interior mode is sin(pi (x + L/2)/L) sin(pi (y + B/2)/B) sinh(kappa (z + T)),
 * which satisfies Laplace's equation, vanishes on the sides and the bottom and meets d phi/dz = K phi on z = 0.
 * Validity: deep water, zero speed. Measured against Capytaine 3.0.0 in docs/engine/hydro.md.
 * @param {object} o
 * @param {number} o.L length, m
 * @param {number} o.B beam, m
 * @param {number} o.T draught, m
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @param {{nx:number, ny:number, nz:number, grading?:number}} [o.mesh] panel counts (default: chosen as above)
 * @param {number} [o.nOmega=14] number of solved frequencies
 * @param {number} [o.lambdaMin] shortest wavelength solved, m (default B)
 * @param {number[]} [o.modes=[0,1,2,3,4]] modes to solve (0 surge, 1 sway, 2 heave, 3 roll, 4 pitch); the others
 *   are left at zero. A model that already holds all five is re-used.
 * @returns {object} { L, B, T, rho, g, setup, panels (number of unknowns), modes, omega (Float64Array, solved grid),
 *   A, B (per grid frequency 6×6 about the waterline centre; yaw not solved), phi, omegaMax, irregularOmega }
 *   for boxCoefficients and boxExcitation
 */
export function boxHydrodynamics({ L, B, T, rho = 1025, g = G, mesh, nOmega = 14, lambdaMin, modes = [0, 1, 2, 3, 4] }) {
  if (!(L > 0 && B > 0 && T > 0)) throw new RangeError('L, B and T must be > 0');
  if (!(nOmega >= 4)) throw new RangeError('nOmega must be >= 4');
  const ms = Array.from(new Set(modes)).sort();
  if (ms.length === 0 || ms.some((k) => !(k >= 0 && k <= 4))) throw new RangeError('modes must be a non-empty subset of 0..4');
  const base = [L, B, T, rho, g, mesh ? [mesh.nx, mesh.ny, mesh.nz, mesh.grading].join('/') : '', nOmega, lambdaMin || ''].join('|');
  const key = `${base}|${ms.join('')}`;
  const hit = _boxCache.get(key) || _boxCache.get(`${base}|01234`);
  if (hit) return hit;
  let m = mesh;
  if (!m) {
    const ny = 5;
    m = { nx: Math.min(14, Math.max(4, Math.round(1.5 * ny * Math.sqrt(L / B)))), ny, nz: Math.min(6, Math.max(2, Math.round((1.2 * ny * T) / (0.5 * B)))) };
  }
  const panels = hullPanels(boxQuarterPanels({ L, B, T, nx: m.nx, ny: m.ny, nz: m.nz, grading: m.grading === undefined ? 1 : m.grading }));
  const interiorPoints = [[0.29, 0.37, 0.41], [0.67, 0.71, 0.23], [0.83, 0.19, 0.61], [0.13, 0.87, 0.77]].map(([a, b, c]) => [0.5 * L * a, 0.5 * B * b, -T * c]);
  const setup = radiationSetup({ panels, symmetry: { x: true, y: true }, interiorPoints });
  const kappa = Math.PI * Math.hypot(1 / L, 1 / B);
  const irregularOmega = Math.sqrt((g * kappa) / Math.tanh(kappa * T));
  const Kmax = Math.min((2 * Math.PI) / (lambdaMin || B), (0.98 * X_MAX) / setup.rMax);
  const omegaMax = Math.sqrt(g * Kmax);
  const omega = new Float64Array(nOmega); const A = []; const Bd = []; const phi = [];
  for (let i = 0; i < nOmega; i++) {
    omega[i] = (omegaMax * (i + 1)) / nOmega;
    const r = radiationSolve({ setup, omega: omega[i], rho, g, modes: ms, nullField: omega[i] >= 0.9 * irregularOmega });
    A.push(r.A); Bd.push(haskindDamping({ setup, phi: r.phi, omega: omega[i], rho, g })); phi.push(r.phi);
  }
  const model = { L, B, T, rho, g, setup, panels: panels.n, modes: ms, omega, A, B: Bd, phi, omegaMax, irregularOmega };
  if (_boxCache.size >= BOX_CACHE_MAX) _boxCache.delete(_boxCache.keys().next().value);
  _boxCache.set(key, model);
  return model;
}

// grid interval and cubic Lagrange weights at omega (clamped to the solved range)
function boxWeights(model, omega) {
  const n = model.omega.length; const h = model.omega[0];
  const u = Math.min(Math.max(omega / h - 1, 0), n - 1);
  let i0 = Math.floor(u) - 1;
  if (i0 < 0) i0 = 0; else if (i0 > n - 4) i0 = n - 4;
  lagrange(u - i0 - 1, _lw);
  return i0;
}

/**
 * Added mass and wave damping of a box model at any frequency: cubic interpolation between the solved
 * frequencies; below the first and above the last solved frequency the end values are held.
 * Reference: interpolation only; see boxHydrodynamics.
 * @param {object} model from boxHydrodynamics
 * @param {number} omega rad/s
 * @returns {{A:number[][], B:number[][]}} 6×6 about the waterline centre, SI (yaw row and column are 0)
 */
export function boxCoefficients(model, omega) {
  if (!(omega > 0)) throw new RangeError('omega must be > 0');
  const i0 = boxWeights(model, omega);
  const A = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]); const B = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  for (let p = 0; p < 4; p++) {
    const Ap = model.A[i0 + p]; const Bp = model.B[i0 + p]; const w = _lw[p];
    for (let i = 0; i < 5; i++) for (let j = 0; j < 5; j++) { A[i][j] += w * Ap[i][j]; B[i][j] += w * Bp[i][j]; }
  }
  return { A, B };
}

/**
 * Wave exciting forces on a box model at any frequency and heading from the Haskind relation, with the radiation
 * potentials interpolated between the solved frequencies and the undisturbed wave taken at the frequency asked
 * for (so the phase variation of the wave along the hull is exact).
 * Reference: see haskindExcitation.
 * @param {object} model from boxHydrodynamics
 * @param {number} omega wave frequency, rad/s
 * @param {number} [headingDeg=0] direction the waves travel towards, from +x towards +y, deg
 * @returns {{re:number[], im:number[]}} N/m and N·m/m about the waterline centre (yaw is 0)
 */
export function boxExcitation(model, omega, headingDeg = 0) {
  if (!(omega > 0)) throw new RangeError('omega must be > 0');
  const i0 = boxWeights(model, omega);
  const n2 = 2 * model.setup.n;
  const phi = [null, null, null, null, null, null];
  for (let k = 0; k < 5; k++) {
    if (!model.phi[i0][k]) continue;
    const f = new Float64Array(n2);
    for (let p = 0; p < 4; p++) { const src = model.phi[i0 + p][k]; const w = _lw[p]; for (let i = 0; i < n2; i++) f[i] += w * src[i]; }
    phi[k] = f;
  }
  return haskindExcitation({ setup: model.setup, phi, omega, headingDeg, rho: model.rho, g: model.g });
}
