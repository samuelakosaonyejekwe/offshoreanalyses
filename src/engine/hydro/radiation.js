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
// Cost. The Rankine matrices are computed once per hull (closed form within 2.5 panel diameters, the panel's
// 2 × 2 Gauss rule beyond: 3e-5 of the coefficients) and reused at every frequency and for every class. At a
// frequency the wave part is assembled from data that do not depend on it (square roots, logarithms and
// reciprocals of the pair distances are kept with the setup), and the system is solved as K_r + i U V^T: K_r
// real, and U V^T the outgoing-wave part exp(K(z + zeta)) J0(K R), which separates into a short series of
// products of functions of the field point and of the source point (Graf's addition theorem; 5–15 terms per
// class). One real elimination with U as extra right-hand sides and a small complex system replace the complex
// elimination; frequencies with null-field rows keep the complex one. The same series gives the heading integral
// of the Haskind damping in closed form. Both are checked against the complex elimination and against the
// heading quadrature in the tests.
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
// distance, in panel diameters, beyond which the Rankine integrals over a panel are taken with its Gauss rule
const FAR = 2.5;

/**
 * Prepares a set of quadrilateral hull panels for radiationSolve: centroids, normals, areas, the two triangles
 * of every quadrilateral (for the closed-form Rankine integrals of bem.js) and a 2 × 2 Gauss rule per panel.
 * Reference: geometry only.
 * @param {ArrayLike<number>} quads 12 numbers per panel: four vertices in order, counter-clockwise when seen
 *   from the fluid (the normal (v1 - v0) × (v3 - v0) points into the fluid), m; a triangle repeats its last
 *   vertex. With symmetry only the part of the hull in x >= 0 and/or y >= 0 is given.
 * @returns {object} panel set: { n, quads, centroids, normals, areas, diam, tri, triGeo, triStart, gauss }
 */
export function hullPanels(quads) {
  if (!quads || quads.length < 12 || quads.length % 12 !== 0) throw new RangeError('quads must hold four vertices (12 numbers) per panel');
  const n = quads.length / 12;
  const centroids = new Float64Array(3 * n); const normals = new Float64Array(3 * n); const areas = new Float64Array(n); const diam = new Float64Array(n);
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
    for (let c = 0; c < 4; c++) diam[p] = Math.max(diam[p], 2 * Math.hypot(quads[q + 3 * c] - centroids[3 * p], quads[q + 3 * c + 1] - centroids[3 * p + 1], quads[q + 3 * c + 2] - centroids[3 * p + 2]));
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
  return { n, quads: Float64Array.from(quads), centroids, normals, areas, diam, tri, triGeo: panelGeometry(tri), triStart, gauss };
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

const _r5 = new Float64Array(5);

// One row of the frequency-independent matrices: field point pt (a panel centre or an interior point, possibly
// mirrored) against every panel; `self` is the index of the panel the point is the centre of, or -1. Kept as a
// function of its own (called once per row, with the point in a typed array) so that it is compiled early.
function setupRow(row, i, self, pt) {
  const { panels, n } = row;
  const C = panels.centroids; const N = panels.normals; const gp = panels.gauss; const diam = panels.diam; const triStart = panels.triStart;
  const px = pt[0]; const py = pt[1]; const pz = pt[2]; const r = _r5;
  const Sm = row.S; const Dm = row.D; const Rm = row.R; const dm = row.dRn; const sq = row.sq; const lg = row.lg; const gd = row.gd; const iv = row.iv;
  for (let j = 0; j < n; j++) {
    let pot = 0; let ker = 0;
    const cx = C[3 * j] - px; const cy = C[3 * j + 1] - py; const far2 = FAR * FAR * diam[j] * diam[j];
    const nx = N[3 * j]; const ny = N[3 * j + 1]; const nz = N[3 * j + 2];
    for (let im = 0; im < 2; im++) {
      const fz = im === 0 ? pz : -pz; const cz = C[3 * j + 2] - fz;
      if (cx * cx + cy * cy + cz * cz > far2) {
        // far from the panel: 2 × 2 Gauss rule for ∫ 1/r dS and for the signed solid angle -∫ n·(p - xi)/r³ dS
        for (let q = 16 * j; q < 16 * j + 16; q += 4) {
          const dx = px - gp[q]; const dy = py - gp[q + 1]; const dz = fz - gp[q + 2];
          const r2 = dx * dx + dy * dy + dz * dz; const w = gp[q + 3] / Math.sqrt(r2);
          pot += w; ker -= (w * (dx * nx + dy * ny + dz * nz)) / r2;
        }
      } else {
        for (let t = triStart[j]; t < triStart[j + 1]; t++) { panelIntegral(panels.tri, panels.triGeo, t, px, py, fz, self === j && im === 0, r); pot += r[0]; ker += r[4]; }
      }
    }
    const q = i * n + j;
    Sm[q] = pot; Dm[q] = ker;
    const rr = Math.sqrt(cx * cx + cy * cy);
    Rm[q] = rr; dm[q] = rr > 0 ? (cx * nx + cy * ny) / rr : 0;
    const a = -(pz + C[3 * j + 2]); const d = Math.sqrt(rr * rr + a * a);
    if (row.zs) { row.zs[q] = -a; row.sqA[q] = Math.sqrt(a); }
    sq[q] = Math.sqrt(rr); lg[q] = Math.log(d + a); gd[q] = rr / (d * (d + a)); iv[q] = 1 / d;
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
  const C = panels.centroids; const N = panels.normals; const Ar = panels.areas; const gp = panels.gauss;
  for (let i = 0; i < n; i++) if ((symX && C[3 * i] <= 0) || (symY && C[3 * i + 1] <= 0)) throw new RangeError('with symmetry the panels must lie in x > 0 / y > 0');
  // S: ∫ (1/r + 1/r') dS, D: minus the normal-dipole integral (signed solid angle), rows = field points
  const S = []; const D = []; const R = []; const dRn = [];
  // per pair and mirror, for the wave part: sqrt(R), ln(d + |z + zeta|), R/(d (d + |z + zeta|)) and 1/d with
  // d = sqrt(R² + (z + zeta)²); none of them depends on the frequency
  const sqR = []; const lgD = []; const gD = []; const invD = [];
  const zs = new Float64Array(nr * n); const sqA = new Float64Array(nr * n);
  const pt = new Float64Array(3);
  for (let m = 0; m < nm; m++) {
    const row = { panels, n, S: new Float64Array(nr * n), D: new Float64Array(nr * n), R: new Float64Array(nr * n), dRn: new Float64Array(nr * n),
      sq: new Float64Array(nr * n), lg: new Float64Array(nr * n), gd: new Float64Array(nr * n), iv: new Float64Array(nr * n), zs: m === 0 ? zs : null, sqA: m === 0 ? sqA : null };
    const [sx, sy] = mirrors[m];
    for (let i = 0; i < nr; i++) {
      const src = i < n ? C : interiorPoints[i - n]; const o = i < n ? 3 * i : 0;
      pt[0] = sx * src[o]; pt[1] = sy * src[o + 1]; pt[2] = src[o + 2];
      if (!(pt[2] < 0)) throw new RangeError('interior points must lie below z = 0');
      setupRow(row, i, m === 0 && i < n ? i : -1, pt);
    }
    S.push(row.S); D.push(row.D); R.push(row.R); dRn.push(row.dRn); sqR.push(row.sq); lgD.push(row.lg); gD.push(row.gd); invD.push(row.iv);
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
  let rhoMax = 0;
  for (let i = 0; i < panels.gauss.length; i += 4) rhoMax = Math.max(rhoMax, Math.hypot(panels.gauss[i], panels.gauss[i + 1]));
  for (let i = 0; i < n; i++) rhoMax = Math.max(rhoMax, Math.hypot(C[3 * i], C[3 * i + 1]));
  const zRow = new Float64Array(nr);
  for (let i = 0; i < nr; i++) zRow[i] = i < n ? C[3 * i + 2] : interiorPoints[i - n][2];
  return { panels, mirrors, classes, n, np, nr, R, dRn, zs, sqA, sqR, lgD, gD, invD, zRow, gn, ref: Array.from(ref), rMax, rhoMax, areas: Ar };
}

// Bessel functions J_0 … J_P at z by downward recurrence (Miller's algorithm), normalised with
// J_0 + 2 (J_2 + J_4 + …) = 1 (Abramowitz & Stegun ch. 9, generating function and backward recurrence).
function besselJs(z, P, out) {
  if (z < 1e-9) { out.fill(0, 0, P + 1); out[0] = 1; return; }
  const top = Math.max(P, Math.ceil(z));
  const M = 2 * Math.ceil(0.5 * (top + 12 + Math.sqrt(48 * top + 60)));
  let jp = 0; let j = 1e-280; let sum = 0;
  for (let k = M; k > 0; k--) {
    const jm = ((2 * k) / z) * j - jp;
    jp = j; j = jm;
    if (k - 1 <= P) out[k - 1] = j;
    if (((k - 1) & 1) === 0) sum += k === 1 ? j : 2 * j;
    if (Math.abs(j) > 1e250) { j *= 1e-250; jp *= 1e-250; sum *= 1e-250; for (let q = k - 1; q <= P; q++) out[q] *= 1e-250; }
  }
  const inv = 1 / sum;
  for (let q = 0; q <= P; q++) out[q] *= inv;
}

// Separable form of the outgoing-wave (imaginary) part of the Green function: with polar coordinates (r, theta)
// of the field point and (rho, theta') of the source point about the vertical axis through the origin,
//   exp(K (z + zeta)) J0(K R) = sum_p eps_p h_p(x) · h_p(xi),  h = exp(K z) J_p(K r) {cos, sin}(p theta),
// eps_0 = 1, eps_p = 2 (Graf's addition theorem, Abramowitz & Stegun ch. 9). For a hull with both planes of
// symmetry each class of modes keeps one family: cos with p even (heave), cos with p odd (surge, pitch), sin with
// p odd (sway, roll), sin with p even (yaw). Fills h and the normal derivative dh/dn (times `scale`) of the
// family of class signs (ex, ey) for np points (`per` points per panel normal): out.h[i * nb + b],
// out.d[i * nb + b]; returns the number of orders and the first one.
function waveFamily(K, P, ex, ey, np, xyz, stride, normals, per, out) {
  const cosine = ey === 1; const first = cosine ? (ex === 1 ? 0 : 1) : (ex === 1 ? 1 : 2);
  const nb = first > P ? 0 : Math.floor((P - first) / 2) + 1;
  if (!out.h || out.h.length < np * nb) { out.h = new Float64Array(np * nb); out.d = new Float64Array(np * nb); }
  if (!out.J || out.J.length < P + 3) out.J = new Float64Array(P + 3);
  const J = out.J; const h = out.h; const d = out.d;
  for (let i = 0; i < np; i++) {
    const x = xyz[stride * i]; const y = xyz[stride * i + 1]; const z = xyz[stride * i + 2];
    const no = 3 * (per === 1 ? i : Math.floor(i / per));
    const nx = normals[no]; const ny = normals[no + 1]; const nz = normals[no + 2];
    const rho = Math.hypot(x, y); const ct = x / rho; const st = y / rho;
    const nr = nx * ct + ny * st; const nt = -nx * st + ny * ct;
    besselJs(K * rho, P + 1, J);
    const e = Math.exp(K * z);
    // cos(p theta), sin(p theta) by the angle-addition recurrence
    let cp = 1; let sp = 0;
    for (let q = 0; q < first; q++) { const c2 = cp * ct - sp * st; sp = sp * ct + cp * st; cp = c2; }
    const c2t = ct * ct - st * st; const s2t = 2 * st * ct;
    for (let b = 0; b < nb; b++) {
      const pp = first + 2 * b;
      const Jp = J[pp]; const dJ = pp === 0 ? -J[1] : 0.5 * (J[pp - 1] - J[pp + 1]);
      const tr = cosine ? cp : sp; const dtr = cosine ? -sp : cp; // d/d(p theta)
      h[i * nb + b] = e * Jp * tr;
      d[i * nb + b] = e * (K * nz * Jp * tr + K * nr * dJ * tr + (nt * pp * Jp * dtr) / rho);
      const c3 = cp * c2t - sp * s2t; sp = sp * c2t + cp * s2t; cp = c3;
    }
  }
  return { nb, first };
}

// Gaussian elimination with partial pivoting on a real n × nc matrix (row-major) whose columns n … nc - 1 are
// right-hand sides; they are overwritten by the solutions.
function rsolveAugmented(a, n, nc) {
  for (let c = 0; c < n; c++) {
    let p = c; let best = Math.abs(a[c * nc + c]);
    for (let r = c + 1; r < n; r++) { const v = Math.abs(a[r * nc + c]); if (v > best) { best = v; p = r; } }
    if (!(best > 0)) throw new Error('radiationSolve: the influence matrix is singular (duplicate or degenerate panels?)');
    const rc = c * nc;
    if (p !== c) { const rp = p * nc; for (let j = c; j < nc; j++) { const t = a[rc + j]; a[rc + j] = a[rp + j]; a[rp + j] = t; } }
    const inv = 1 / a[rc + c];
    for (let r = c + 1; r < n; r++) {
      const rr = r * nc; const f = a[rr + c] * inv;
      if (f === 0) continue;
      let j = c + 1;
      for (; j + 3 < nc; j += 4) { a[rr + j] -= f * a[rc + j]; a[rr + j + 1] -= f * a[rc + j + 1]; a[rr + j + 2] -= f * a[rc + j + 2]; a[rr + j + 3] -= f * a[rc + j + 3]; }
      for (; j < nc; j++) a[rr + j] -= f * a[rc + j];
    }
  }
  for (let r = n - 1; r >= 0; r--) {
    const rr = r * nc; const inv = 1 / a[rr + r];
    for (let k = n; k < nc; k++) a[rr + k] *= inv;
    for (let q = 0; q < r; q++) { const f = a[q * nc + r]; if (f !== 0) { const qq = q * nc; for (let k = n; k < nc; k++) a[qq + k] -= f * a[rr + k]; } }
  }
}

const _fam = {}; const _famG = {};
// order at which the series of waveFamily is cut: the terms fall off faster than exponentially beyond p = z
const familyOrder = (z) => Math.ceil(z + 2 * Math.cbrt(z) + 3);

const _g = { F: 0, FX: 0, FY: 0, J0: 0, J1: 0 };
const _t4 = new Float64Array(4);

// Smooth parts of the tables at grid coordinates qx = sqrt(X)/HQ and qs = sqrt(|Y|)/HS (qs below the last node):
// reads them from _t4[0], _t4[1] and writes T, dT/dX, J0, J1 to _t4 (numbers passed as arguments between
// functions are allocated on the heap by the engine; a typed array is not). The same interpolation as waveGreenFunction, without its logarithm and
// exponential, which radiationSolve takes from the frequency-independent pair data of the setup.
function tableAt(tb) {
  const qx = _t4[0]; const qs = _t4[1];
  let ix = Math.floor(qx) - 1;
  if (ix < 0) ix = 0; else if (ix > NQ - 4) ix = NQ - 4;
  let u = qx - ix - 1;
  const x0 = (-u * (u - 1) * (u - 2)) / 6; const x1 = ((u + 1) * (u - 1) * (u - 2)) / 2; const x2 = (-(u + 1) * u * (u - 2)) / 2; const x3 = ((u + 1) * u * (u - 1)) / 6;
  let is = Math.floor(qs) - 1;
  if (is < 0) is = 0; else if (is > NS - 4) is = NS - 4;
  u = qs - is - 1;
  const s0 = (-u * (u - 1) * (u - 2)) / 6; const s1 = ((u + 1) * (u - 1) * (u - 2)) / 2; const s2 = (-(u + 1) * u * (u - 2)) / 2; const s3 = ((u + 1) * u * (u - 1)) / 6;
  const T = tb.T; const TX = tb.TX;
  let o = ix * NS + is;
  let t = x0 * (s0 * T[o] + s1 * T[o + 1] + s2 * T[o + 2] + s3 * T[o + 3]); let tx = x0 * (s0 * TX[o] + s1 * TX[o + 1] + s2 * TX[o + 2] + s3 * TX[o + 3]);
  o += NS; t += x1 * (s0 * T[o] + s1 * T[o + 1] + s2 * T[o + 2] + s3 * T[o + 3]); tx += x1 * (s0 * TX[o] + s1 * TX[o + 1] + s2 * TX[o + 2] + s3 * TX[o + 3]);
  o += NS; t += x2 * (s0 * T[o] + s1 * T[o + 1] + s2 * T[o + 2] + s3 * T[o + 3]); tx += x2 * (s0 * TX[o] + s1 * TX[o + 1] + s2 * TX[o + 2] + s3 * TX[o + 3]);
  o += NS; t += x3 * (s0 * T[o] + s1 * T[o + 1] + s2 * T[o + 2] + s3 * T[o + 3]); tx += x3 * (s0 * TX[o] + s1 * TX[o + 1] + s2 * TX[o + 2] + s3 * TX[o + 3]);
  _t4[0] = t; _t4[1] = tx;
  _t4[2] = x0 * tb.J0[ix] + x1 * tb.J0[ix + 1] + x2 * tb.J0[ix + 2] + x3 * tb.J0[ix + 3];
  _t4[3] = x0 * tb.J1[ix] + x1 * tb.J1[ix + 1] + x2 * tb.J1[ix + 2] + x3 * tb.J1[ix + 3];
}

// Wave part of the influence matrices for mirror copy m at wave number K: wS = A_j Gw and wK = -A_j dGw/dn at
// the source panel (real parts; imaginary parts wSi, wKi unless `realOnly`). Gw depends on the pair of points
// only through R and z + zeta, so on the hull it is evaluated once for (i, j) and (j, i).
function waveMatrices(setup, m, K, nr, realOnly, tb) {
  const { n, areas, zs, sqA } = setup; const N = setup.panels.normals;
  const twoK = 2 * K; const twoK2 = 2 * K * K;
  const sK = Math.sqrt(K); const logK = Math.log(K); const invK = 1 / K; const qxs = sK / HQ; const qss = sK / HS;
  const ez = setup.ez; const t4 = _t4; const qsMax = NS - 1; const splitSolve = realOnly;
  const sq = setup.sqR[m]; const lg = setup.lgD[m]; const gd = setup.gD[m]; const iv = setup.invD[m]; const dn = setup.dRn[m]; const Rm = setup.R[m];
  const WS = setup.wS[m]; const WK = setup.wK[m]; const WSi = splitSolve ? null : setup.wSi[m]; const WKi = splitSolve ? null : setup.wKi[m];
  for (let i = 0; i < nr; i++) {
    const both = i < n; const ezi = ez[i]; const ai = both ? areas[i] : 0; const nzi = both ? N[3 * i + 2] : 0;
    for (let j = both ? i : 0; j < n; j++) {
      const q = i * n + j;
      const ey = ezi * ez[j]; const qs = qss * sqA[q];
      let F; let FX; let FY;
      if (qs < qsMax) {
        t4[0] = qxs * sq[q]; t4[1] = qs;
        tableAt(tb);
        F = t4[0] - ey * (logK + lg[q]); FX = t4[1] - ey * gd[q] * invK; FY = F + iv[q] * invK;
      } else {
        waveGreenFunction(K * Rm[q], K * zs[q], _g);
        F = _g.F; FX = _g.FX; FY = _g.FY; t4[2] = _g.J0; t4[3] = _g.J1;
      }
      const aj = areas[j]; const nzj = N[3 * j + 2];
      WS[q] = aj * twoK * F; WK[q] = -aj * twoK2 * (FX * dn[q] + FY * nzj);
      const two = both && j !== i; const q2 = j * n + i;
      if (two) { WS[q2] = ai * twoK * F; WK[q2] = -ai * twoK2 * (FX * dn[q2] + FY * nzi); }
      if (splitSolve) continue;
      const e0 = Math.PI * ey * t4[2]; const e1 = Math.PI * ey * t4[3];
      WSi[q] = -aj * twoK * e0; WKi[q] = -aj * twoK2 * (e1 * dn[q] - e0 * nzj);
      if (two) { WSi[q2] = -ai * twoK * e0; WKi[q2] = -ai * twoK2 * (e1 * dn[q2] - e0 * nzi); }
    }
  }
}

// Signed sums of the mirror matrices of one class: system matrix (added to the Rankine part already in ac.kc)
// and right-hand sides -S n_k.
function combineClass(setup, ac, nr, splitSolve) {
  const n = setup.n; const nm = setup.mirrors.length;
  const kc = ac.kc; const sg = ac.cl.signs; const nrh = ac.rhs.length;
  const g0 = ac.gns[0]; const g1 = nrh > 1 ? ac.gns[1] : g0; const b0 = ac.rhs[0]; const b1 = nrh > 1 ? ac.rhs[1] : null;
  for (let part = 0; part < (splitSolve ? 1 : 2); part++) {
    const WKs = part === 0 ? setup.wK : setup.wKi; const WSs = part === 0 ? setup.wS : setup.wSi;
    if (nm === 4) {
      const s1 = sg[1]; const s2 = sg[2]; const s3 = sg[3];
      const K0 = WKs[0]; const K1 = WKs[1]; const K2 = WKs[2]; const K3 = WKs[3]; const S0 = WSs[0]; const S1 = WSs[1]; const S2 = WSs[2]; const S3 = WSs[3];
      for (let i = 0; i < nr; i++) {
        let a0 = 0; let a1 = 0;
        for (let j = 0; j < n; j++) {
          const q = i * n + j;
          kc[2 * q + part] += K0[q] + s1 * K1[q] + s2 * K2[q] + s3 * K3[q];
          const sv = S0[q] + s1 * S1[q] + s2 * S2[q] + s3 * S3[q];
          a0 += sv * g0[j]; a1 += sv * g1[j];
        }
        b0[2 * i + part] -= a0;
        if (b1) b1[2 * i + part] -= a1;
      }
    } else {
      for (let m = 0; m < nm; m++) {
        const w = sg[m]; const WK = WKs[m]; const WS = WSs[m];
        for (let i = 0; i < nr; i++) {
          const ro = i * n;
          for (let j = 0; j < n; j++) kc[2 * (ro + j) + part] += w * WK[ro + j];
          for (let r = 0; r < nrh; r++) { const gv = ac.gns[r]; let sum = 0; for (let j = 0; j < n; j++) sum += WS[ro + j] * gv[j]; ac.rhs[r][2 * i + part] -= w * sum; }
        }
      }
    }
  }
}

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
 * @param {boolean} [o.split=true] solve through the real part of the system and the separable wave part (used
 *   with both planes of symmetry when no null-field rows are active); false forces the complex elimination
 * @returns {{omega:number, A:number[][], B:number[][], phi:(Float64Array|null)[], symmetryError:number}} A in
 *   kg, kg·m, kg·m²; B in N·s/m, N·s, N·m·s (about `ref`); phi per mode: complex potential per unit velocity on
 *   the panels (re, im interleaved), m²/s per m/s or per rad/s
 */
export function radiationSolve({ setup, omega, rho = 1025, g = G, modes, nullField = true, split = true }) {
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
  const splitSolve = split && nm === 4 && nr === n; // then only the real parts are assembled (see below)
  // wave part: one matrix pair per mirror copy (waveMatrices), then the signed sums of every class
  const tb = tables(Math.min(NQ - 1, Math.ceil(Math.sqrt(K * setup.rMax) / HQ) + 4));
  if (!setup.ez) { setup.ez = new Float64Array(setup.nr); setup.wS = mirrors.map(() => new Float64Array(setup.nr * n)); setup.wK = mirrors.map(() => new Float64Array(setup.nr * n)); }
  if (!splitSolve && !setup.wSi) { setup.wSi = mirrors.map(() => new Float64Array(setup.nr * n)); setup.wKi = mirrors.map(() => new Float64Array(setup.nr * n)); }
  for (let i = 0; i < nr; i++) setup.ez[i] = Math.exp(K * setup.zRow[i]);
  for (let m = 0; m < nm; m++) waveMatrices(setup, m, K, nr, splitSolve, tb);
  for (const ac of act) combineClass(setup, ac, nr, splitSolve);
  const phi = [null, null, null, null, null, null];
  // Without null-field rows and with both planes of symmetry the system is K_r + i U V^T with K_r real and
  // U V^T the separable wave part of low rank (waveFamily): one real elimination with U and the real right-hand
  // sides as columns, then a small complex system. phi = y + i Z c, K_r y = b_r, K_r Z = U,
  // (I + i V^T Z) c = t - V^T y, where the imaginary part of the right-hand side is U t.
  const P = splitSolve ? familyOrder(K * setup.rhoMax) : 0;
  for (const { kc, rhs, ms, gns, cl } of act) {
    if (splitSolve) {
      const { nb } = waveFamily(K, P, cl.signs[1], cl.signs[2], n, panels.centroids, 3, N, 1, _fam);
      const h = _fam.h; const dh = _fam.d; const nmd = ms.length; const nc = n + nb + nmd;
      if (!cl.aug || cl.aug.length < n * nc) cl.aug = new Float64Array(n * (n + nb + 8 + nmd));
      const aug = cl.aug; const c2 = 2 * Math.PI * K;
      for (let i = 0; i < n; i++) {
        const ro = i * nc;
        for (let j = 0; j < n; j++) aug[ro + j] = kc[2 * (i * n + j)];
        for (let b = 0; b < nb; b++) aug[ro + n + b] = (cl.signs[1] === 1 && cl.signs[2] === 1 && b === 0 ? c2 : 2 * c2) * h[i * nb + b];
        for (let c = 0; c < nmd; c++) aug[ro + n + nb + c] = rhs[c][2 * i];
      }
      rsolveAugmented(aug, n, nc);
      // V^T Z and the right-hand sides t - V^T y of the small system
      const sys = new Float64Array(2 * nb * nb); const cols = ms.map(() => new Float64Array(2 * nb));
      for (let j = 0; j < n; j++) {
        const ro = j * nc; const a4 = 4 * areas[j];
        for (let b = 0; b < nb; b++) {
          const v = a4 * dh[j * nb + b]; const hv = a4 * h[j * nb + b];
          for (let q = 0; q < nb; q++) sys[2 * (b * nb + q) + 1] += v * aug[ro + n + q];
          for (let c = 0; c < nmd; c++) cols[c][2 * b] += hv * gns[c][j] - v * aug[ro + n + nb + c];
        }
      }
      for (let b = 0; b < nb; b++) sys[2 * (b * nb + b)] = 1;
      if (nb > 0) csolveDense(sys, nb, cols);
      for (let c = 0; c < nmd; c++) {
        const f = new Float64Array(2 * n); const cc = cols[c];
        for (let i = 0; i < n; i++) {
          const ro = i * nc; let re = aug[ro + n + nb + c]; let im = 0;
          for (let b = 0; b < nb; b++) { const zb = aug[ro + n + b]; re -= zb * cc[2 * b + 1]; im += zb * cc[2 * b]; }
          f[2 * i] = re; f[2 * i + 1] = im;
        }
        phi[ms[c]] = f;
      }
      continue;
    }
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
  const N = panels.normals; const gp = panels.gauss; const nm = mirrors.length;
  // I_i = ∫ (phi_0 n_i - phi_i d phi_0/dn) dS / (i g/omega)
  const Ir = [0, 0, 0, 0, 0, 0]; const Ii = [0, 0, 0, 0, 0, 0];
  const ks = []; const par = [];
  for (let k = 0; k < 6; k++) if (phi[k]) { ks.push(k); par.push(mirrors.map(([sx, sy]) => (sx === -1 ? PARITY[k][0] : 1) * (sy === -1 ? PARITY[k][1] : 1))); }
  const nk = ks.length;
  for (let p = 0; p < n; p++) {
    const n0x = N[3 * p]; const n0y = N[3 * p + 1]; const nz = N[3 * p + 2];
    for (let q = 0; q < 4; q++) {
      const o = 16 * p + 4 * q;
      const x0 = gp[o]; const y0 = gp[o + 1]; const z = gp[o + 2];
      // exp(-i K (x cb + y sb)) of the mirror copies from one cosine and sine per coordinate
      const amp = gp[o + 3] * Math.exp(K * z);
      const cx = Math.cos(K * x0 * cb); const sx0 = Math.sin(K * x0 * cb); const cy = Math.cos(K * y0 * sb); const sy0 = Math.sin(K * y0 * sb);
      const rz = z - ref[2];
      for (let m = 0; m < nm; m++) {
        const sx = mirrors[m][0]; const sy = mirrors[m][1];
        const ax = -sx * sx0; const ay = -sy * sy0; // imaginary parts of the two factors
        const er = amp * (cx * cy - ax * ay); const ei = amp * (cx * ay + ax * cy);
        const nx = sx * n0x; const ny = sy * n0y;
        // d phi_0/dn = phi_0 (K nz - i K (nx cb + ny sb))
        const dr = K * nz; const di = -K * (nx * cb + ny * sb);
        const rx = sx * x0 - ref[0]; const ry = sy * y0 - ref[1];
        for (let c = 0; c < nk; c++) {
          const k = ks[c];
          const nkv = k === 0 ? nx : k === 1 ? ny : k === 2 ? nz : k === 3 ? ry * nz - rz * ny : k === 4 ? rz * nx - rx * nz : rx * ny - ry * nx;
          // phi_k on the mirror copy = parity × phi_k on the panel
          const pr = par[c][m]; const fr = pr * phi[k][2 * p]; const fi = pr * phi[k][2 * p + 1];
          const tr = nkv - (fr * dr - fi * di); const ti = -(fr * di + fi * dr);
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
 * With both planes of symmetry the heading integral is taken in closed form (Fourier–Bessel series of the
 * undisturbed wave); otherwise, or when `nHeadings` is given, by the midpoint rule over the headings.
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
  if (quarter && nHeadings === undefined) {
    // both planes of symmetry: the heading integral in closed form. With exp(-i K rho cos(theta - beta)) =
    // sum (-i)^p J_p(K rho) exp(i p (theta - beta)) the forces are Fourier series in beta, and by Parseval
    // B_ij = rho omega³/(2 g) sum_p eps_p Re{C_ip C_jp*}, C_ip = ∫ (n_i - phi_i d/dn) h_p dS (h_p of waveFamily).
    const { panels, classes, n, ref } = setup; const gp = panels.gauss; const N = panels.normals;
    const B = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
    const P = familyOrder(K * setup.rhoMax);
    for (const cl of classes) {
      const ms = cl.modes.filter((k) => phi[k]);
      if (ms.length === 0) continue;
      const { nb, first } = waveFamily(K, P, cl.signs[1], cl.signs[2], 4 * n, gp, 4, N, 4, _famG);
      const h = _famG.h; const dh = _famG.d;
      const Cr = ms.map(() => new Float64Array(nb)); const Ci = ms.map(() => new Float64Array(nb));
      for (let p = 0; p < n; p++) {
        const nx = N[3 * p]; const ny = N[3 * p + 1]; const nz = N[3 * p + 2];
        for (let q = 0; q < 4; q++) {
          const o = 16 * p + 4 * q; const w = 4 * gp[o + 3];
          const rx = gp[o] - ref[0]; const ry = gp[o + 1] - ref[1]; const rz = gp[o + 2] - ref[2];
          const row = (4 * p + q) * nb;
          for (let c = 0; c < ms.length; c++) {
            const k = ms[c];
            const nkv = k === 0 ? nx : k === 1 ? ny : k === 2 ? nz : k === 3 ? ry * nz - rz * ny : k === 4 ? rz * nx - rx * nz : rx * ny - ry * nx;
            const fr = phi[k][2 * p]; const fi = phi[k][2 * p + 1]; const cr = Cr[c]; const ci = Ci[c];
            for (let b = 0; b < nb; b++) { const dv = w * dh[row + b]; cr[b] += w * nkv * h[row + b] - fr * dv; ci[b] -= fi * dv; }
          }
        }
      }
      const f = (rho * omega ** 3) / (2 * g);
      for (let a = 0; a < ms.length; a++) for (let c = 0; c < ms.length; c++) {
        let sum = 0;
        for (let b = 0; b < nb; b++) sum += (first === 0 && b === 0 ? 1 : 2) * (Cr[a][b] * Cr[c][b] + Ci[a][b] * Ci[c][b]);
        B[ms[a]][ms[c]] = f * sum;
      }
    }
    return B;
  }
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

const _boxCache = new Map(); const _unitCache = new Map();
const BOX_CACHE_MAX = 12;
const _lw = new Float64Array(4);

/**
 * Hydrodynamic model of a rectangular box L × B × T floating in deep water at zero speed: radiation problems in
 * surge, sway, heave, roll and pitch solved by radiationSolve on a grid of frequencies and kept (potentials
 * included) for interpolation, so that coefficients and exciting forces at any frequency and heading cost one
 * pass over the panels. The box is solved in non-dimensional form (unit beam, density and gravity; Froude
 * similarity) and kept by its ratios L/B and T/B, so hulls that differ only in scale, density or gravity share
 * one solution; the dimensional model is cached per hull.
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
 * @returns {object} { L, B, T, rho, g, panels (number of unknowns), modes, omega (Float64Array, solved grid, rad/s),
 *   omegaMax, irregularOmega, unit } for boxCoefficients and boxExcitation; `unit` is the non-dimensional
 *   solution { lb, tb, mesh, setup, omega, A, B (per grid frequency 6×6 about the waterline centre; yaw not
 *   solved), phi, omegaMax, irregularOmega } in units of B, rho B³ and sqrt(B/g)
 */
export function boxHydrodynamics({ L, B, T, rho = 1025, g = G, mesh, nOmega = 14, lambdaMin, modes = [0, 1, 2, 3, 4] }) {
  if (!(L > 0 && B > 0 && T > 0)) throw new RangeError('L, B and T must be > 0');
  if (!(nOmega >= 4)) throw new RangeError('nOmega must be >= 4');
  const ms = Array.from(new Set(modes)).sort();
  if (ms.length === 0 || ms.some((k) => !(k >= 0 && k <= 4))) throw new RangeError('modes must be a non-empty subset of 0..4');
  const tail = [mesh ? [mesh.nx, mesh.ny, mesh.nz, mesh.grading].join('/') : '', nOmega, lambdaMin ? (lambdaMin / B).toPrecision(12) : ''].join('|');
  const dimKey = [L, B, T, rho, g, tail].join('|');
  const hit = _boxCache.get(`${dimKey}|${ms.join('')}`) || _boxCache.get(`${dimKey}|01234`);
  if (hit) return hit;
  // Froude similarity: the box is solved with unit beam, density and gravity, and kept by its ratios L/B and T/B
  const lb = Number((L / B).toPrecision(12)); const tb = Number((T / B).toPrecision(12));
  const unitKey = [lb, tb, tail].join('|');
  let unit = _unitCache.get(`${unitKey}|${ms.join('')}`) || _unitCache.get(`${unitKey}|01234`);
  if (!unit) {
    let m = mesh;
    if (!m) {
      const ny = 5;
      m = { nx: Math.min(14, Math.max(4, Math.round(1.5 * ny * Math.sqrt(lb)))), ny, nz: Math.min(6, Math.max(2, Math.round((1.2 * ny * tb) / 0.5))) };
    }
    const panels = hullPanels(boxQuarterPanels({ L: lb, B: 1, T: tb, nx: m.nx, ny: m.ny, nz: m.nz, grading: m.grading === undefined ? 1 : m.grading }));
    const interiorPoints = [[0.29, 0.37, 0.41], [0.67, 0.71, 0.23], [0.83, 0.19, 0.61], [0.13, 0.87, 0.77]].map(([a, b, c]) => [0.5 * lb * a, 0.5 * b, -tb * c]);
    const setup = radiationSetup({ panels, symmetry: { x: true, y: true }, interiorPoints });
    const kappa = Math.PI * Math.hypot(1 / lb, 1);
    const irregularOmega = Math.sqrt(kappa / Math.tanh(kappa * tb));
    const Kmax = Math.min((2 * Math.PI * B) / (lambdaMin || B), (0.98 * X_MAX) / setup.rMax);
    const omegaMax = Math.sqrt(Kmax);
    const omega = new Float64Array(nOmega); const A = []; const Bd = []; const phi = [];
    for (let i = 0; i < nOmega; i++) {
      omega[i] = (omegaMax * (i + 1)) / nOmega;
      const r = radiationSolve({ setup, omega: omega[i], rho: 1, g: 1, modes: ms, nullField: omega[i] >= 0.9 * irregularOmega });
      A.push(r.A); Bd.push(haskindDamping({ setup, phi: r.phi, omega: omega[i], rho: 1, g: 1 })); phi.push(r.phi);
    }
    unit = { lb, tb, mesh: { nx: m.nx, ny: m.ny, nz: m.nz }, setup, panels: panels.n, modes: ms, omega, A, B: Bd, phi, omegaMax, irregularOmega };
    if (_unitCache.size >= BOX_CACHE_MAX) _unitCache.delete(_unitCache.keys().next().value);
    _unitCache.set(`${unitKey}|${ms.join('')}`, unit);
  }
  const ws = Math.sqrt(g / B);
  const model = { L, B, T, rho, g, unit, panels: unit.panels, modes: unit.modes, omega: Float64Array.from(unit.omega, (w) => w * ws), omegaMax: unit.omegaMax * ws, irregularOmega: unit.irregularOmega * ws };
  if (_boxCache.size >= 4 * BOX_CACHE_MAX) _boxCache.delete(_boxCache.keys().next().value);
  _boxCache.set(`${dimKey}|${ms.join('')}`, model);
  return model;
}

// grid interval and cubic Lagrange weights at the non-dimensional frequency (clamped to the solved range)
function boxWeights(unit, omega) {
  const n = unit.omega.length; const h = unit.omega[0];
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
  const { unit, B: beam, rho, g } = model;
  const i0 = boxWeights(unit, omega * Math.sqrt(beam / g));
  const A = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]); const B = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  for (let p = 0; p < 4; p++) {
    const Ap = unit.A[i0 + p]; const Bp = unit.B[i0 + p]; const w = _lw[p];
    for (let i = 0; i < 5; i++) for (let j = 0; j < 5; j++) { A[i][j] += w * Ap[i][j]; B[i][j] += w * Bp[i][j]; }
  }
  // back to dimensions: rho B³ per mass, one more power of B per rotational index, sqrt(g/B) for the damping
  const a0 = rho * beam ** 3; const wsc = Math.sqrt(g / beam);
  for (let i = 0; i < 5; i++) for (let j = 0; j < 5; j++) { const sc = a0 * (i > 2 ? beam : 1) * (j > 2 ? beam : 1); A[i][j] *= sc; B[i][j] *= sc * wsc; }
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
  const { unit, B: beam, rho, g } = model;
  const wu = omega * Math.sqrt(beam / g);
  const i0 = boxWeights(unit, wu);
  const n2 = 2 * unit.setup.n;
  const phi = [null, null, null, null, null, null];
  for (let k = 0; k < 5; k++) {
    if (!unit.phi[i0][k]) continue;
    const f = new Float64Array(n2);
    for (let p = 0; p < 4; p++) { const src = unit.phi[i0 + p][k]; const w = _lw[p]; for (let i = 0; i < n2; i++) f[i] += w * src[i]; }
    phi[k] = f;
  }
  const x = haskindExcitation({ setup: unit.setup, phi, omega: wu, headingDeg, rho: 1, g: 1 });
  const f0 = rho * g * beam * beam;
  for (let k = 0; k < 6; k++) { const sc = f0 * (k > 2 ? beam : 1); x.re[k] *= sc; x.im[k] *= sc; }
  return x;
}
