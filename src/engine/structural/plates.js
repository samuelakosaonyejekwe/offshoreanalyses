/**
 * Plates and stiffened panels: classical bending solutions, buckling, effective width, ultimate strength,
 * grillages and hull-girder bending.
 *
 * Conventions: rectangular plate a × b in the x–y plane (x from 0 to a, y from 0 to b), thickness t, lateral
 * pressure positive in +z, deflection w positive in +z. Bending moments per unit width are positive when they
 * produce tension on the face opposite the load (sagging). Compressive membrane stresses are entered as
 * positive numbers in the buckling functions.
 *
 * References: Timoshenko & Woinowsky-Krieger, "Theory of Plates and Shells" (2nd ed.); Timoshenko & Gere,
 * "Theory of Elastic Stability" (2nd ed.) ch. 9; EN 1993-1-5:2006 §4.4, §4.5, §5.3, §10; EN 1993-1-1:2005
 * §6.3.1; Faulkner (1975) J. Ship Research 19(1); Hughes & Paik, "Ship Structural Analysis and Design" (2010);
 * Smith (1977) PRADS (progressive collapse of hull girders).
 */
import { solveDense } from '../math/linalg.js';
import { eigGenSym } from '../math/eigen.js';
import { gaussLegendre } from '../math/quad.js';
import { feSystem } from './fem2d.js';

function pos(v, name, fn) {
  if (!(typeof v === 'number' && v > 0 && Number.isFinite(v))) throw new RangeError(`${fn}: ${name} must be a positive number`);
}
function nonNeg(v, name, fn) {
  if (!(typeof v === 'number' && v >= 0 && Number.isFinite(v))) throw new RangeError(`${fn}: ${name} must be a non-negative number`);
}
function poisson(nu, fn) {
  if (!(typeof nu === 'number' && nu >= 0 && nu < 0.5)) throw new RangeError(`${fn}: nu must lie in [0, 0.5)`);
}

/**
 * Flexural rigidity of an isotropic plate, D = E t³ / (12 (1 − ν²)).
 * @param {object} o
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} o.t thickness [m]
 * @returns {number} D [N·m]
 */
export function plateRigidity({ E, nu, t }) {
  pos(E, 'E', 'plateRigidity'); poisson(nu, 'plateRigidity'); pos(t, 't', 'plateRigidity');
  return (E * t * t * t) / (12 * (1 - nu * nu));
}

/* ------------------------------------------------------------------------------------------------------------
 * Navier solution (all edges simply supported)
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Kirchhoff bending of a simply supported rectangular plate by Navier's double sine series
 * (Timoshenko & Woinowsky-Krieger §28–§30): w = Σ Σ q_mn / (D π⁴ (m²/a² + n²/b²)²) sin(mπx/a) sin(nπy/b).
 * Valid for thin plates (t ≲ min(a, b)/10) with small deflections (w ≲ t/2).
 * @param {object} o
 * @param {number} o.a plate length in x [m]
 * @param {number} o.b plate length in y [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {object} o.load { type: 'uniform', q [Pa] } | { type: 'patch', P [N] total, x0, y0 [m] centre, u, v [m]
 *   patch sizes in x and y } | { type: 'point', P [N], x0, y0 [m] }
 * @param {number} [o.x] evaluation point x [m] (default a/2)
 * @param {number} [o.y] evaluation point y [m] (default b/2)
 * @param {number} [o.nTerms] highest harmonic in each direction (default 99; point loads converge slowly in
 *   the moments, which are singular under the load)
 * @returns {{w: number, mx: number, my: number, mxy: number, qx: number, qy: number, sigmaX: number,
 *   sigmaY: number, D: number, alpha: number|null, betaX: number|null, betaY: number|null}}
 *   w [m]; mx, my, mxy [N·m/m] (mx bends the x-direction fibres; mxy = −D(1 − ν) w,xy); qx, qy shear forces
 *   [N/m]; extreme-fibre bending stresses [Pa]; for uniform load the classical coefficients
 *   alpha = w D / (q s⁴), betaX = mx / (q s²), betaY = my / (q s²) with s = min(a, b)
 */
export function navierPlate(o) {
  const fn = 'navierPlate';
  const { a, b, t, E, nu, load, nTerms = 99 } = o;
  pos(a, 'a', fn); pos(b, 'b', fn);
  const D = plateRigidity({ E, nu, t });
  if (!load || !['uniform', 'patch', 'point'].includes(load.type)) throw new RangeError(`${fn}: load.type must be 'uniform', 'patch' or 'point'`);
  if (!(Number.isInteger(nTerms) && nTerms >= 1)) throw new RangeError(`${fn}: nTerms must be a positive integer`);
  const x = o.x === undefined ? a / 2 : o.x;
  const y = o.y === undefined ? b / 2 : o.y;
  if (!(x >= 0 && x <= a)) throw new RangeError(`${fn}: x must lie within the plate`);
  if (!(y >= 0 && y <= b)) throw new RangeError(`${fn}: y must lie within the plate`);
  if (load.type !== 'uniform') {
    if (!(load.x0 >= 0 && load.x0 <= a && load.y0 >= 0 && load.y0 <= b)) throw new RangeError(`${fn}: load position must lie within the plate`);
    if (!Number.isFinite(load.P)) throw new RangeError(`${fn}: load.P must be a number`);
    if (load.type === 'patch') {
      pos(load.u, 'load.u', fn); pos(load.v, 'load.v', fn);
      if (load.x0 - load.u / 2 < -1e-12 * a || load.x0 + load.u / 2 > a * (1 + 1e-12) || load.y0 - load.v / 2 < -1e-12 * b || load.y0 + load.v / 2 > b * (1 + 1e-12)) {
        throw new RangeError(`${fn}: patch must lie within the plate`);
      }
    }
  } else if (!Number.isFinite(load.q)) throw new RangeError(`${fn}: load.q must be a number`);
  const pi = Math.PI;
  let w = 0; let mx = 0; let my = 0; let mxy = 0; let qx = 0; let qy = 0;
  const step = load.type === 'uniform' ? 2 : 1;
  for (let m = 1; m <= nTerms; m += step) {
    const am = (m * pi) / a;
    const sx = Math.sin(am * x);
    const cx = Math.cos(am * x);
    for (let n = 1; n <= nTerms; n += step) {
      const bn = (n * pi) / b;
      let qmn;
      if (load.type === 'uniform') qmn = (16 * load.q) / (pi * pi * m * n);
      else if (load.type === 'point') qmn = ((4 * load.P) / (a * b)) * Math.sin(am * load.x0) * Math.sin(bn * load.y0);
      else {
        qmn = ((16 * load.P) / (load.u * load.v * pi * pi * m * n)) * Math.sin(am * load.x0) * Math.sin(bn * load.y0)
          * Math.sin((am * load.u) / 2) * Math.sin((bn * load.v) / 2);
      }
      if (qmn === 0) continue;
      const k2 = am * am + bn * bn;
      const amn = qmn / (D * k2 * k2);
      const sy = Math.sin(bn * y);
      const cy = Math.cos(bn * y);
      w += amn * sx * sy;
      mx += D * amn * (am * am + nu * bn * bn) * sx * sy;
      my += D * amn * (nu * am * am + bn * bn) * sx * sy;
      mxy -= D * (1 - nu) * amn * am * bn * cx * cy;
      qx += D * amn * am * k2 * cx * sy;
      qy += D * amn * bn * k2 * sx * cy;
    }
  }
  const s = Math.min(a, b);
  const uni = load.type === 'uniform' && load.q !== 0;
  return {
    w, mx, my, mxy, qx, qy, sigmaX: (6 * mx) / (t * t), sigmaY: (6 * my) / (t * t), D,
    alpha: uni ? (w * D) / (load.q * s ** 4) : null,
    betaX: uni ? mx / (load.q * s * s) : null,
    betaY: uni ? my / (load.q * s * s) : null,
  };
}

/* ------------------------------------------------------------------------------------------------------------
 * Rayleigh–Ritz solution for general edge conditions
 * ---------------------------------------------------------------------------------------------------------- */

function polyMul(p, q) {
  const r = new Array(p.length + q.length - 1).fill(0);
  for (let i = 0; i < p.length; i++) for (let j = 0; j < q.length; j++) r[i + j] += p[i] * q[j];
  return r;
}
function polyDer(p) {
  const r = [];
  for (let i = 1; i < p.length; i++) r.push(i * p[i]);
  return r.length ? r : [0];
}
function polyVal(p, x) {
  let s = 0;
  for (let i = p.length - 1; i >= 0; i--) s = s * x + p[i];
  return s;
}

const EDGE_POWER = { F: 0, S: 1, C: 2 };

/**
 * One-dimensional Ritz basis on ξ ∈ [−1, 1]: φ_i = (1 + ξ)^p (1 − ξ)^q P_i(ξ), with p, q = 0 (free), 1 (simply
 * supported) or 2 (clamped) so that every function satisfies the geometric edge conditions, and P_i the
 * Legendre polynomials (ascending-power coefficient arrays). Returns basis values/derivatives at the Gauss
 * points and the five integral matrices needed by the plate energy.
 */
function ritzBasis(e0, e1, n) {
  const p0 = EDGE_POWER[e0];
  const p1 = EDGE_POWER[e1];
  let w = [1];
  for (let k = 0; k < p0; k++) w = polyMul(w, [1, 1]);
  for (let k = 0; k < p1; k++) w = polyMul(w, [1, -1]);
  const leg = [[1], [0, 1]];
  for (let k = 2; k < n; k++) {
    // (k) P_k = (2k − 1) ξ P_{k−1} − (k − 1) P_{k−2}
    const a = polyMul([0, 1], leg[k - 1]).map((v) => (v * (2 * k - 1)) / k);
    const b = leg[k - 2].map((v) => (v * (k - 1)) / k);
    leg.push(a.map((v, i) => v - (b[i] || 0)));
  }
  const f = []; const d1 = []; const d2 = [];
  for (let i = 0; i < n; i++) {
    const p = polyMul(w, leg[i]);
    f.push(p); d1.push(polyDer(p)); d2.push(polyDer(polyDer(p)));
  }
  const g = gaussLegendre(Math.min(64, n + 6));
  const ng = g.x.length;
  const V = []; const V1 = []; const V2 = [];
  for (let i = 0; i < n; i++) {
    V.push(g.x.map((x) => polyVal(f[i], x)));
    V1.push(g.x.map((x) => polyVal(d1[i], x)));
    V2.push(g.x.map((x) => polyVal(d2[i], x)));
  }
  const integ = (A, B) => {
    const M = [];
    for (let i = 0; i < n; i++) {
      const row = new Array(n);
      for (let j = 0; j < n; j++) {
        let s = 0;
        for (let k = 0; k < ng; k++) s += g.w[k] * A[i][k] * B[j][k];
        row[j] = s;
      }
      M.push(row);
    }
    return M;
  };
  return { n, f, d1, d2, I00: integ(V, V), I11: integ(V1, V1), I22: integ(V2, V2), I02: integ(V, V2), I01: integ(V, V1), g };
}

function ritzSetup(o, fn) {
  const { a, b, t, E, nu, edges = {}, nBasis = 8 } = o;
  pos(a, 'a', fn); pos(b, 'b', fn);
  const D = plateRigidity({ E, nu, t });
  if (!(Number.isInteger(nBasis) && nBasis >= 2 && nBasis <= 14)) throw new RangeError(`${fn}: nBasis must be an integer in 2..14`);
  const ed = { x0: edges.x0 || 'S', x1: edges.x1 || 'S', y0: edges.y0 || 'S', y1: edges.y1 || 'S' };
  for (const k of Object.keys(ed)) if (!(ed[k] in EDGE_POWER)) throw new RangeError(`${fn}: edges.${k} must be 'S', 'C' or 'F'`);
  const X = ritzBasis(ed.x0, ed.x1, nBasis);
  const Y = ritzBasis(ed.y0, ed.y1, nBasis);
  const n = nBasis;
  const N = n * n;
  const sx = 2 / a; const sy = 2 / b; const J = (a * b) / 4;
  const K = [];
  for (let r = 0; r < N; r++) K.push(new Array(N).fill(0));
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const r = i * n + j;
    for (let k = 0; k < n; k++) for (let l = 0; l < n; l++) {
      const c = k * n + l;
      K[r][c] = D * J * (sx ** 4 * X.I22[i][k] * Y.I00[j][l] + sy ** 4 * X.I00[i][k] * Y.I22[j][l]
        + nu * sx * sx * sy * sy * (X.I02[k][i] * Y.I02[j][l] + X.I02[i][k] * Y.I02[l][j])
        + 2 * (1 - nu) * sx * sx * sy * sy * X.I11[i][k] * Y.I11[j][l]);
    }
  }
  return { a, b, t, nu, D, X, Y, n, N, sx, sy, J, K, ed };
}

function ritzEval(S, c, x, y) {
  const xi = (2 * x) / S.a - 1;
  const et = (2 * y) / S.b - 1;
  const { n, X, Y, sx, sy, D, nu } = S;
  const fx = []; const fx1 = []; const fx2 = []; const fy = []; const fy1 = []; const fy2 = [];
  for (let i = 0; i < n; i++) {
    fx.push(polyVal(X.f[i], xi)); fx1.push(polyVal(X.d1[i], xi)); fx2.push(polyVal(X.d2[i], xi));
    fy.push(polyVal(Y.f[i], et)); fy1.push(polyVal(Y.d1[i], et)); fy2.push(polyVal(Y.d2[i], et));
  }
  let w = 0; let wxx = 0; let wyy = 0; let wxy = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const cc = c[i * n + j];
    w += cc * fx[i] * fy[j];
    wxx += cc * sx * sx * fx2[i] * fy[j];
    wyy += cc * sy * sy * fx[i] * fy2[j];
    wxy += cc * sx * sy * fx1[i] * fy1[j];
  }
  return { w, mx: -D * (wxx + nu * wyy), my: -D * (wyy + nu * wxx), mxy: -D * (1 - nu) * wxy };
}

/**
 * Kirchhoff bending of a rectangular plate with any combination of simply supported (S), clamped (C) and free
 * (F) edges by the Rayleigh–Ritz method with polynomial trial functions that satisfy the geometric edge
 * conditions (boundary-weighted Legendre polynomials). Converges rapidly for smooth loads; replaces the
 * Lévy/series tables for clamped and mixed edges. Validity as for `navierPlate`. Moments on clamped edges
 * converge more slowly than deflections (third digit with nBasis ≈ 10).
 * @param {object} o
 * @param {number} o.a plate length in x [m]
 * @param {number} o.b plate length in y [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {{x0?: string, x1?: string, y0?: string, y1?: string}} [o.edges] edge condition 'S' | 'C' | 'F' on
 *   x = 0, x = a, y = 0, y = b (default all 'S')
 * @param {object} o.load as in `navierPlate` (uniform, patch or point)
 * @param {number} [o.x] evaluation point x [m] (default a/2)
 * @param {number} [o.y] evaluation point y [m] (default b/2)
 * @param {number} [o.nBasis] trial functions per direction, 2..14 (default 8)
 * @param {number} [o.grid] when given, also return the deflection sampled on a (grid + 1)² lattice
 * @returns {{w: number, mx: number, my: number, mxy: number, sigmaX: number, sigmaY: number, wMax: number,
 *   wGrid?: number[][],
 *   edgeMoments: {x0: number, x1: number, y0: number, y1: number}, D: number, strainEnergy: number}}
 *   values at (x, y): w [m], moments [N·m/m], bending stresses [Pa]; wMax: largest |w| on a 21 × 21 sampling
 *   grid [m]; edgeMoments: normal bending moment at each edge mid-point [N·m/m] (negative = hogging at a
 *   clamped edge); strainEnergy [J]; wGrid[j][i] = w at (a i/grid, b j/grid) [m]
 */
export function ritzPlate(o) {
  const fn = 'ritzPlate';
  const S = ritzSetup(o, fn);
  const { load } = o;
  const { a, b, n, N, X, Y } = S;
  if (!load || !['uniform', 'patch', 'point'].includes(load.type)) throw new RangeError(`${fn}: load.type must be 'uniform', 'patch' or 'point'`);
  const F = new Array(N).fill(0);
  if (load.type === 'point') {
    if (!(load.x0 >= 0 && load.x0 <= a && load.y0 >= 0 && load.y0 <= b)) throw new RangeError(`${fn}: load position must lie within the plate`);
    const xi = (2 * load.x0) / a - 1;
    const et = (2 * load.y0) / b - 1;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) F[i * n + j] = load.P * polyVal(X.f[i], xi) * polyVal(Y.f[j], et);
  } else {
    let x0 = 0; let x1 = a; let y0 = 0; let y1 = b; let q = load.q;
    if (load.type === 'patch') {
      pos(load.u, 'load.u', fn); pos(load.v, 'load.v', fn);
      x0 = load.x0 - load.u / 2; x1 = load.x0 + load.u / 2; y0 = load.y0 - load.v / 2; y1 = load.y0 + load.v / 2;
      if (x0 < -1e-12 * a || x1 > a * (1 + 1e-12) || y0 < -1e-12 * b || y1 > b * (1 + 1e-12)) throw new RangeError(`${fn}: patch must lie within the plate`);
      q = load.P / (load.u * load.v);
    }
    if (!Number.isFinite(q)) throw new RangeError(`${fn}: load magnitude must be a number`);
    const g = X.g;
    const line = (B, lo, hi, len) => {
      // ∫ φ_i over [lo, hi] in physical length units, by Gauss–Legendre on the sub-interval.
      const out = [];
      const c = (lo + hi) / len - 1;
      const hw = (hi - lo) / len;
      for (let i = 0; i < n; i++) {
        let s = 0;
        for (let k = 0; k < g.x.length; k++) s += g.w[k] * polyVal(B.f[i], c + hw * g.x[k]);
        out.push(s * hw * (len / 2));
      }
      return out;
    };
    const ix = line(X, x0, x1, a);
    const iy = line(Y, y0, y1, b);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) F[i * n + j] = q * ix[i] * iy[j];
  }
  let c;
  try {
    c = solveDense(S.K, F);
  } catch {
    throw new Error('Plate is not sufficiently supported for bending (rigid-body mode with these edge conditions).');
  }
  let energy = 0;
  for (let i = 0; i < N; i++) energy += 0.5 * c[i] * F[i];
  if (!Number.isFinite(energy)) throw new Error('Plate is not sufficiently supported for bending (rigid-body mode with these edge conditions).');
  const x = o.x === undefined ? a / 2 : o.x;
  const y = o.y === undefined ? b / 2 : o.y;
  const at = ritzEval(S, c, x, y);
  let wMax = 0;
  for (let i = 0; i <= 20; i++) for (let j = 0; j <= 20; j++) wMax = Math.max(wMax, Math.abs(ritzEval(S, c, (a * i) / 20, (b * j) / 20).w));
  return {
    ...at, sigmaX: (6 * at.mx) / (S.t * S.t), sigmaY: (6 * at.my) / (S.t * S.t), wMax,
    edgeMoments: { x0: ritzEval(S, c, 0, b / 2).mx, x1: ritzEval(S, c, a, b / 2).mx, y0: ritzEval(S, c, a / 2, 0).my, y1: ritzEval(S, c, a / 2, b).my },
    D: S.D, strainEnergy: energy,
    ...(o.grid ? { wGrid: Array.from({ length: o.grid + 1 }, (_, j) => Array.from({ length: o.grid + 1 }, (__, i) => ritzEval(S, c, (a * i) / o.grid, (b * j) / o.grid).w)) } : {}),
  };
}

/**
 * Bending coefficients of a uniformly loaded rectangular plate with all edges clamped, computed by the Ritz
 * solution (not tabulated): w_max = alpha q s⁴ / D, edge and centre moments = beta q s² with s = min(a, b).
 * @param {object} o
 * @param {number} o.aspect a / b [-], 0.2..5
 * @param {number} [o.nu] Poisson's ratio (default 0.3; the deflection and edge-moment coefficients do not depend
 *   on it, the centre moments do)
 * @param {number} [o.nBasis] trial functions per direction (default 10)
 * @returns {{alpha: number, betaEdgeLong: number, betaEdgeShort: number, betaCentreShortSpan: number,
 *   betaCentreLongSpan: number}} coefficients [-]: centre deflection; moment at the middle of the long edges
 *   (largest, hogging, negative) and of the short edges; centre moments for the short-span and long-span
 *   directions
 */
export function clampedPlateCoefficients({ aspect, nu = 0.3, nBasis = 10 }) {
  const fn = 'clampedPlateCoefficients';
  if (!(aspect >= 0.2 && aspect <= 5)) throw new RangeError(`${fn}: aspect must lie in 0.2..5`);
  const a = aspect; const b = 1;
  const r = ritzPlate({ a, b, t: 1, E: 12 * (1 - nu * nu), nu, edges: { x0: 'C', x1: 'C', y0: 'C', y1: 'C' }, load: { type: 'uniform', q: 1 }, nBasis });
  const s = Math.min(a, b);
  // With E = 12(1 − ν²) and t = 1 the rigidity is exactly 1.
  const xShort = a <= b;
  return {
    alpha: r.w / s ** 4,
    betaEdgeLong: (xShort ? r.edgeMoments.x0 : r.edgeMoments.y0) / (s * s),
    betaEdgeShort: (xShort ? r.edgeMoments.y0 : r.edgeMoments.x0) / (s * s),
    betaCentreShortSpan: (xShort ? r.mx : r.my) / (s * s),
    betaCentreLongSpan: (xShort ? r.my : r.mx) / (s * s),
  };
}

/**
 * Elastic buckling of a rectangular plate with general edge conditions under uniform Nx, Ny, Nxy by the
 * Rayleigh–Ritz method (same trial functions as `ritzPlate`). The result is an upper bound that converges
 * quickly with nBasis; shear buckling of long plates needs nBasis ≥ 10–12.
 * @param {object} o
 * @param {number} o.a plate length in x [m]
 * @param {number} o.b plate length in y [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {object} [o.edges] edge conditions as in `ritzPlate`
 * @param {number} [o.sigmaX] reference compressive stress in x [Pa] (compression positive)
 * @param {number} [o.sigmaY] reference compressive stress in y [Pa]
 * @param {number} [o.tau] reference shear stress [Pa]
 * @param {number} [o.nBasis] trial functions per direction (default 10)
 * @returns {{factor: number, sigmaXcr: number, sigmaYcr: number, tauCr: number, k: number, sigmaE: number}}
 *   factor: critical multiplier on the reference stresses [-]; critical stresses [Pa]; k: buckling coefficient
 *   referred to the largest reference stress, σ_cr = k σE with σE = π² D / (b² t) [Pa]
 */
export function ritzPlateBuckling(o) {
  const fn = 'ritzPlateBuckling';
  const { sigmaX = 0, sigmaY = 0, tau = 0 } = o;
  const S = ritzSetup({ nBasis: 10, ...o }, fn);
  if (sigmaX === 0 && sigmaY === 0 && tau === 0) throw new RangeError(`${fn}: give at least one non-zero reference stress`);
  const { n, N, X, Y, sx, sy, J, t } = S;
  const G = [];
  for (let r = 0; r < N; r++) G.push(new Array(N).fill(0));
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const r = i * n + j;
    for (let k = 0; k < n; k++) for (let l = 0; l < n; l++) {
      G[r][k * n + l] = t * J * (sigmaX * sx * sx * X.I11[i][k] * Y.I00[j][l] + sigmaY * sy * sy * X.I00[i][k] * Y.I11[j][l]
        + tau * sx * sy * (X.I01[k][i] * Y.I01[j][l] + X.I01[i][k] * Y.I01[l][j]));
    }
  }
  // G φ = μ K φ with K positive definite; the largest μ gives the smallest positive load factor 1/μ.
  let eg;
  try {
    eg = eigGenSym(G, S.K, {});
  } catch {
    throw new Error('Plate is not sufficiently supported for a buckling analysis with these edge conditions.');
  }
  let mu = -Infinity;
  for (const v of eg.values) if (v > mu) mu = v;
  if (!(mu > 0)) throw new Error('The reference stresses are not compressive: no buckling.');
  const factor = 1 / mu;
  const sigmaE = (Math.PI * Math.PI * S.D) / (S.b * S.b * t);
  const ref = Math.max(Math.abs(sigmaX), Math.abs(sigmaY), Math.abs(tau));
  return { factor, sigmaXcr: factor * sigmaX, sigmaYcr: factor * sigmaY, tauCr: factor * tau, k: (factor * ref) / sigmaE, sigmaE };
}

/**
 * Natural frequencies of a rectangular Kirchhoff plate with general edge conditions (Rayleigh–Ritz).
 * @param {object} o
 * @param {number} o.a plate length in x [m]
 * @param {number} o.b plate length in y [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} o.rho density [kg/m³]
 * @param {number} [o.addedMass] non-structural or hydrodynamic added mass per unit area [kg/m²]
 * @param {object} [o.edges] edge conditions as in `ritzPlate`
 * @param {number} [o.nModes] number of modes returned (default 4)
 * @param {number} [o.nBasis] trial functions per direction (default 8)
 * @returns {{frequencies: number[], omega: number[]}} natural frequencies [Hz] and [rad/s], ascending
 *   (rigid-body modes of insufficiently supported plates appear as ≈ 0)
 */
export function ritzPlateFrequencies(o) {
  const fn = 'ritzPlateFrequencies';
  const { rho, addedMass = 0, nModes = 4 } = o;
  pos(rho, 'rho', fn); nonNeg(addedMass, 'addedMass', fn);
  const S = ritzSetup(o, fn);
  const { n, N, X, Y, J, t } = S;
  const M = [];
  for (let r = 0; r < N; r++) M.push(new Array(N).fill(0));
  const m = rho * t + addedMass;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) for (let l = 0; l < n; l++) {
    M[i * n + j][k * n + l] = m * J * X.I00[i][k] * Y.I00[j][l];
  }
  const eg = eigGenSym(S.K, M, { nModes: Math.min(N, nModes) });
  const omega = eg.values.slice(0, nModes).map((v) => Math.sqrt(Math.max(0, v)));
  return { frequencies: omega.map((w) => w / (2 * Math.PI)), omega };
}

/* ------------------------------------------------------------------------------------------------------------
 * Large deflection
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Large-deflection (membrane-stiffened) response of a simply supported rectangular plate under uniform
 * pressure: one-term Galerkin solution of the von Kármán equations with w = w0 sin(πx/a) sin(πy/b)
 * (Timoshenko & Woinowsky-Krieger ch. 13; Chia, "Nonlinear Analysis of Plates", 1980, §2). The Airy stress
 * function is F = Px y²/2 + Py x²/2 + (E w0²/32) [(a/b)² cos(2πx/a) + (b/a)² cos(2πy/b)], and the centre
 * deflection follows from
 *   D π⁴ (1/a² + 1/b²)² w0 + E t π⁴ w0³ [c_edge + (1/a⁴ + 1/b⁴)/16] = 16 q / π²,
 * with c_edge = (1/a⁴ + 1/b⁴ + 2ν/(a² b²)) / (8 (1 − ν²)) for edges held against pull-in (immovable in the
 * mean, remaining straight) and c_edge = 0 for edges free to pull in.
 * Approximate: the one-term shape over-predicts the small-deflection limit by about 2.5 % for a square plate
 * and the accuracy degrades for w0/t ≳ 2–3 and for aspect ratios above about 3.
 * @param {object} o
 * @param {number} o.a plate length in x [m]
 * @param {number} o.b plate length in y [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.nu Poisson's ratio [-]
 * @param {number} o.q uniform pressure [Pa] (>= 0)
 * @param {'immovable'|'movable'} [o.edges] in-plane edge restraint (default 'immovable')
 * @returns {{w0: number, w0Linear: number, ratio: number, membraneX: number, membraneY: number, bendingX: number,
 *   bendingY: number, sigmaXmax: number, sigmaYmax: number, vonMisesCentre: number}}
 *   w0: centre deflection [m]; w0Linear: one-term small-deflection value [m]; ratio = w0/t [-]; membrane and
 *   bending stresses at the plate centre [Pa]; total tensile-face stresses and von Mises stress at the centre [Pa]
 */
export function largeDeflectionPlate({ a, b, t, E, nu, q, edges = 'immovable' }) {
  const fn = 'largeDeflectionPlate';
  pos(a, 'a', fn); pos(b, 'b', fn); nonNeg(q, 'q', fn);
  const D = plateRigidity({ E, nu, t });
  if (edges !== 'immovable' && edges !== 'movable') throw new RangeError(`${fn}: edges must be 'immovable' or 'movable'`);
  const pi = Math.PI;
  const ia2 = 1 / (a * a); const ib2 = 1 / (b * b);
  const c1 = D * pi ** 4 * (ia2 + ib2) ** 2;
  const cEdge = edges === 'immovable' ? (ia2 * ia2 + ib2 * ib2 + 2 * nu * ia2 * ib2) / (8 * (1 - nu * nu)) : 0;
  const c3 = E * t * pi ** 4 * (cEdge + (ia2 * ia2 + ib2 * ib2) / 16);
  const rhs = (16 * q) / (pi * pi);
  const wl = rhs / c1;
  // Monotonic cubic: Newton from the linear solution converges unconditionally.
  let w = wl;
  for (let it = 0; it < 100; it++) {
    const f = c1 * w + c3 * w ** 3 - rhs;
    const dw = f / (c1 + 3 * c3 * w * w);
    w -= dw;
    if (Math.abs(dw) <= 1e-14 * Math.abs(w)) break;
  }
  const k = (E * w * w * pi * pi) / 8;
  const px = edges === 'immovable' ? (k / (1 - nu * nu)) * (ia2 + nu * ib2) : 0;
  const py = edges === 'immovable' ? (k / (1 - nu * nu)) * (ib2 + nu * ia2) : 0;
  const membraneX = px + k * ia2;
  const membraneY = py + k * ib2;
  const bendingX = (6 * D * w * pi * pi * (ia2 + nu * ib2)) / (t * t);
  const bendingY = (6 * D * w * pi * pi * (ib2 + nu * ia2)) / (t * t);
  const sx = membraneX + bendingX;
  const sy = membraneY + bendingY;
  return {
    w0: w, w0Linear: wl, ratio: w / t, membraneX, membraneY, bendingX, bendingY, sigmaXmax: sx, sigmaYmax: sy,
    vonMisesCentre: Math.sqrt(sx * sx + sy * sy - sx * sy),
  };
}

/* ------------------------------------------------------------------------------------------------------------
 * Buckling coefficients, slenderness, effective width, ultimate strength
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Elastic buckling coefficient k of a rectangular plate (σ_cr = k π² E / (12 (1 − ν²)) (t/b)², b = loaded-edge
 * width for compression in x).
 *  - simply supported, uniaxial compression in x: k = min over m of (m b/a + a/(m b))² (Timoshenko & Gere §9.2)
 *  - simply supported, biaxial σy = ratio·σx: k = min over m, n of (m²/α² + n²)² / (m²/α² + ratio n²), α = a/b
 *    (Timoshenko & Gere §9.3)
 *  - simply supported, shear: k = 5.34 + 4/α² (α ≥ 1), 4 + 5.34/α² (α < 1) (EN 1993-1-5 Annex A.3)
 *  - clamped edges (all four): computed by `ritzPlateBuckling` (no closed form)
 * @param {object} o
 * @param {number} o.aspect α = a / b [-] (0.1..20)
 * @param {'uniaxial'|'biaxial'|'shear'} [o.loading] default 'uniaxial'
 * @param {'ss'|'clamped'} [o.edges] default 'ss'
 * @param {number} [o.ratio] σy/σx for biaxial loading (negative = transverse tension)
 * @param {number} [o.nu] Poisson's ratio for the clamped Ritz computation (default 0.3; k is independent of ν
 *   for these edge conditions)
 * @returns {{k: number, m: number|null, n: number|null}} coefficient [-] and the half-wave numbers of the
 *   critical mode where they are defined
 */
export function plateBucklingCoefficient({ aspect, loading = 'uniaxial', edges = 'ss', ratio = 0, nu = 0.3 }) {
  const fn = 'plateBucklingCoefficient';
  if (!(aspect >= 0.1 && aspect <= 20)) throw new RangeError(`${fn}: aspect must lie in 0.1..20`);
  if (!['uniaxial', 'biaxial', 'shear'].includes(loading)) throw new RangeError(`${fn}: loading must be 'uniaxial', 'biaxial' or 'shear'`);
  if (edges !== 'ss' && edges !== 'clamped') throw new RangeError(`${fn}: edges must be 'ss' or 'clamped'`);
  if (edges === 'clamped') {
    const c = { x0: 'C', x1: 'C', y0: 'C', y1: 'C' };
    const base = { a: aspect, b: 1, t: 1, E: 1, nu, edges: c, nBasis: aspect > 3 || aspect < 1 / 3 ? 14 : 12 };
    const r = loading === 'shear' ? ritzPlateBuckling({ ...base, tau: 1 })
      : ritzPlateBuckling({ ...base, sigmaX: 1, sigmaY: loading === 'biaxial' ? ratio : 0 });
    return { k: r.sigmaE ? (loading === 'shear' ? r.tauCr : r.sigmaXcr) / r.sigmaE : NaN, m: null, n: null };
  }
  if (loading === 'shear') {
    return { k: aspect >= 1 ? 5.34 + 4 / (aspect * aspect) : 4 + 5.34 / (aspect * aspect), m: null, n: null };
  }
  const beta = loading === 'biaxial' ? ratio : 0;
  let best = Infinity; let bm = null; let bn = null;
  const mMax = Math.max(4, Math.ceil(3 * aspect) + 3);
  const nMax = loading === 'biaxial' ? Math.max(4, Math.ceil(3 / aspect) + 3) : 1;
  for (let m = 1; m <= mMax; m++) for (let n = 1; n <= nMax; n++) {
    const p = (m * m) / (aspect * aspect);
    const den = p + beta * n * n;
    if (den <= 0) continue;
    const k = (p + n * n) ** 2 / den;
    if (k < best) { best = k; bm = m; bn = n; }
  }
  if (!Number.isFinite(best)) throw new Error('The stress combination is not compressive: no buckling.');
  return { k: best, m: bm, n: bn };
}

/**
 * Elastic critical stress and slenderness parameters of a plate panel.
 * σ_E = π² E / (12 (1 − ν²)) (t/b)², σ_cr = k σ_E, λp = √(fy/σ_cr) (EN 1993-1-5 §4.4(2)),
 * β = (b/t) √(fy/E) (Faulkner).
 * @param {object} o
 * @param {number} o.b plate width (loaded edge) [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} o.fy yield stress [Pa]
 * @param {number} [o.k] buckling coefficient (default 4)
 * @returns {{sigmaE: number, sigmaCr: number, lambdaP: number, beta: number}} [Pa, Pa, -, -]
 */
export function plateSlenderness({ b, t, E, nu = 0.3, fy, k = 4 }) {
  const fn = 'plateSlenderness';
  pos(b, 'b', fn); pos(t, 't', fn); pos(E, 'E', fn); poisson(nu, fn); pos(fy, 'fy', fn); pos(k, 'k', fn);
  const sigmaE = ((Math.PI * Math.PI * E) / (12 * (1 - nu * nu))) * (t / b) ** 2;
  const sigmaCr = k * sigmaE;
  return { sigmaE, sigmaCr, lambdaP: Math.sqrt(fy / sigmaCr), beta: (b / t) * Math.sqrt(fy / E) };
}

/**
 * Effective-width ratio ρ = b_eff / b of a compressed plate.
 *  - 'vonKarman': ρ = 1/λp (λp > 1) (von Kármán, Sechler & Donnell 1932)
 *  - 'winter': ρ = (1 − 0.22/λp)/λp (λp > 0.673) (Winter 1947)
 *  - 'en1993-internal': ρ = (λp − 0.055 (3 + ψ))/λp² ≤ 1 for λp > 0.5 + √(0.085 − 0.055 ψ) (EN 1993-1-5 Eq. 4.2)
 *  - 'en1993-outstand': ρ = (λp − 0.188)/λp² ≤ 1 for λp > 0.748 (EN 1993-1-5 Eq. 4.3)
 *  - 'faulkner': ρ = 2/β − 1/β² for β ≥ 1 (Faulkner 1975), β = (b/t)√(fy/E)
 * @param {object} o
 * @param {'vonKarman'|'winter'|'en1993-internal'|'en1993-outstand'|'faulkner'} o.method formula
 * @param {number} [o.lambdaP] plate slenderness √(fy/σcr) [-] (all methods except 'faulkner')
 * @param {number} [o.beta] Faulkner slenderness (b/t)√(fy/E) [-] ('faulkner')
 * @param {number} [o.psi] stress ratio σ2/σ1 across the width, −3 < ψ ≤ 1 ('en1993-internal'; default 1)
 * @returns {number} ρ [-] in (0, 1]
 */
export function effectiveWidthRatio({ method, lambdaP, beta, psi = 1 }) {
  const fn = 'effectiveWidthRatio';
  if (method === 'faulkner') {
    pos(beta, 'beta', fn);
    return beta <= 1 ? 1 : 2 / beta - 1 / (beta * beta);
  }
  pos(lambdaP, 'lambdaP', fn);
  switch (method) {
    case 'vonKarman': return lambdaP <= 1 ? 1 : 1 / lambdaP;
    case 'winter': return lambdaP <= 0.673 ? 1 : (1 - 0.22 / lambdaP) / lambdaP;
    case 'en1993-internal': {
      if (!(psi > -3 && psi <= 1)) throw new RangeError(`${fn}: psi must lie in (-3, 1]`);
      if (lambdaP <= 0.5 + Math.sqrt(0.085 - 0.055 * psi)) return 1;
      return Math.min(1, (lambdaP - 0.055 * (3 + psi)) / (lambdaP * lambdaP));
    }
    case 'en1993-outstand': return lambdaP <= 0.748 ? 1 : Math.min(1, (lambdaP - 0.188) / (lambdaP * lambdaP));
    default: throw new RangeError(`${fn}: unknown method ${method}`);
  }
}

/**
 * Column buckling reduction factor χ (Perry–Robertson form of EN 1993-1-1 §6.3.1.2):
 * Φ = 0.5 [1 + α (λ − 0.2) + λ²], χ = 1 / (Φ + √(Φ² − λ²)) ≤ 1. Imperfection factors α (Table 6.1):
 * a0 0.13, a 0.21, b 0.34, c 0.49, d 0.76.
 * @param {object} o
 * @param {number} o.lambda relative slenderness √(fy/σ_E) [-]
 * @param {'a0'|'a'|'b'|'c'|'d'} [o.curve] buckling curve (default 'b')
 * @param {number} [o.alpha] imperfection factor overriding the curve [-]
 * @returns {number} χ [-]
 */
export function columnReduction({ lambda, curve = 'b', alpha }) {
  const fn = 'columnReduction';
  nonNeg(lambda, 'lambda', fn);
  const table = { a0: 0.13, a: 0.21, b: 0.34, c: 0.49, d: 0.76 };
  const al = alpha === undefined ? table[curve] : alpha;
  if (!(al >= 0)) throw new RangeError(`${fn}: unknown curve or negative alpha`);
  if (lambda <= 0.2) return 1;
  const phi = 0.5 * (1 + al * (lambda - 0.2) + lambda * lambda);
  return Math.min(1, 1 / (phi + Math.sqrt(phi * phi - lambda * lambda)));
}

/**
 * Ultimate strength of an unstiffened rectangular plate panel (all edges simply supported) under in-plane
 * compression in x and/or shear, per EN 1993-1-5:
 *  - plate-type reduction ρ from §4.4 (internal element, ψ = 1) or Faulkner's formula;
 *  - column-type behaviour of short wide panels by §4.5.3–4.5.4: σ_cr,c = π² E t² / (12 (1 − ν²) a²), χc on
 *    curve a (α = 0.21), ρc = (ρ − χc) ξ (2 − ξ) + χc with ξ = σ_cr,p/σ_cr,c − 1 limited to 0..1;
 *  - shear: χw from Table 5.1 with λw = 0.76 √(fy/τ_cr) (η = 1 by default, non-rigid end post unless stated);
 *  - interaction by the reduced-stress format of §10 Eq. (10.5) for uniaxial compression with shear.
 * Partial factors are not applied (characteristic resistance).
 * @param {object} o
 * @param {number} o.a panel length in the direction of compression [m]
 * @param {number} o.b panel width (loaded edge) [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} o.fy yield stress [Pa]
 * @param {'en1993'|'faulkner'} [o.method] plate reduction formula (default 'en1993')
 * @param {boolean} [o.rigidEndPost] shear: rigid end post (default false)
 * @param {number} [o.eta] shear: η of EN 1993-1-5 §5.1 (default 1.0)
 * @param {number} [o.sigmaX] applied compressive stress [Pa] for the utilisation (default 0)
 * @param {number} [o.tau] applied shear stress [Pa] for the utilisation (default 0)
 * @returns {{k: number, sigmaCr: number, lambdaP: number, rho: number, rhoC: number, sigmaU: number,
 *   kTau: number, tauCr: number, lambdaW: number, chiW: number, tauU: number, utilisation: number}}
 *   sigmaU = ρc fy and tauU = χw fy/√3 [Pa]; utilisation = √((σx/σU)² + (τ/τU)²) [-]
 */
export function plateUltimateStrength(o) {
  const fn = 'plateUltimateStrength';
  const { a, b, t, E, nu = 0.3, fy, method = 'en1993', rigidEndPost = false, eta = 1, sigmaX = 0, tau = 0 } = o;
  pos(a, 'a', fn); pos(b, 'b', fn);
  if (method !== 'en1993' && method !== 'faulkner') throw new RangeError(`${fn}: method must be 'en1993' or 'faulkner'`);
  const aspect = Math.min(20, Math.max(0.1, a / b));
  const k = plateBucklingCoefficient({ aspect, loading: 'uniaxial' }).k;
  const s = plateSlenderness({ b, t, E, nu, fy, k });
  const rho = method === 'faulkner' ? effectiveWidthRatio({ method: 'faulkner', beta: s.beta }) : effectiveWidthRatio({ method: 'en1993-internal', lambdaP: s.lambdaP, psi: 1 });
  const sigmaCrC = (Math.PI * Math.PI * E * t * t) / (12 * (1 - nu * nu) * a * a);
  const chiC = columnReduction({ lambda: Math.sqrt(fy / sigmaCrC), curve: 'a' });
  const xi = Math.min(1, Math.max(0, s.sigmaCr / sigmaCrC - 1));
  const rhoC = (rho - chiC) * xi * (2 - xi) + chiC;
  const kTau = plateBucklingCoefficient({ aspect, loading: 'shear' }).k;
  const tauCr = kTau * s.sigmaE;
  const lambdaW = 0.76 * Math.sqrt(fy / tauCr);
  let chiW;
  if (lambdaW < 0.83 / eta) chiW = eta;
  else if (!rigidEndPost || lambdaW < 1.08) chiW = 0.83 / lambdaW;
  else chiW = 1.37 / (0.7 + lambdaW);
  const sigmaU = rhoC * fy;
  const tauU = (chiW * fy) / Math.sqrt(3);
  return {
    k, sigmaCr: s.sigmaCr, lambdaP: s.lambdaP, rho, rhoC, sigmaU, kTau, tauCr, lambdaW, chiW, tauU,
    utilisation: Math.sqrt((sigmaX / sigmaU) ** 2 + (tau / tauU) ** 2),
  };
}

/* ------------------------------------------------------------------------------------------------------------
 * Stiffened panels
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Section properties of a stiffener (flat bar, angle or tee) with an attached plate flange.
 * @param {object} o
 * @param {number} o.be attached (effective) plate width [m]
 * @param {number} o.tp plate thickness [m]
 * @param {number} o.hw web height [m]
 * @param {number} o.tw web thickness [m]
 * @param {number} [o.bf] flange width [m] (0 for a flat bar)
 * @param {number} [o.tf] flange thickness [m]
 * @returns {{A: number, As: number, zNA: number, I: number, r: number, Zplate: number, Zflange: number}}
 *   total and stiffener-only area [m²]; neutral axis above the plate mid-plane [m]; second moment about the
 *   neutral axis [m⁴]; radius of gyration [m]; section moduli at the plate outer fibre and stiffener tip [m³]
 */
export function stiffenerWithPlate({ be, tp, hw, tw, bf = 0, tf = 0 }) {
  const fn = 'stiffenerWithPlate';
  nonNeg(be, 'be', fn); pos(tp, 'tp', fn); pos(hw, 'hw', fn); pos(tw, 'tw', fn); nonNeg(bf, 'bf', fn); nonNeg(tf, 'tf', fn);
  const parts = [
    { A: be * tp, z: 0, i: (be * tp ** 3) / 12 },
    { A: hw * tw, z: tp / 2 + hw / 2, i: (tw * hw ** 3) / 12 },
    { A: bf * tf, z: tp / 2 + hw + tf / 2, i: (bf * tf ** 3) / 12 },
  ];
  let A = 0; let Az = 0;
  for (const p of parts) { A += p.A; Az += p.A * p.z; }
  const zNA = Az / A;
  let I = 0;
  for (const p of parts) I += p.i + p.A * (p.z - zNA) ** 2;
  return { A, As: hw * tw + bf * tf, zNA, I, r: Math.sqrt(I / A), Zplate: I / (zNA + tp / 2), Zflange: I / (tp / 2 + hw + tf - zNA) };
}

/**
 * Column (beam-column) buckling of a stiffener with its effective plate flange — Perry–Robertson approach as
 * used in ship and offshore panel checks (Hughes & Paik ch. 14; EN 1993-1-5 §4.5.3 column-type behaviour):
 * the plate flange is reduced to its effective width, the Euler stress of the stiffener/plate column is
 * σ_E = π² E I_e / (A_e L²), and the ultimate average stress follows from χ(λ) with λ = √(fy/σ_E).
 * The load capacity is referred to the full area: σ_u = χ fy A_e / A.
 * @param {object} o
 * @param {number} o.L stiffener span between transverse frames [m]
 * @param {number} o.s stiffener spacing [m]
 * @param {number} o.tp plate thickness [m]
 * @param {number} o.hw web height [m]
 * @param {number} o.tw web thickness [m]
 * @param {number} [o.bf] flange width [m]
 * @param {number} [o.tf] flange thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} o.fy yield stress [Pa]
 * @param {'en1993'|'faulkner'|'full'} [o.effectiveWidth] plate effectiveness model (default 'en1993')
 * @param {'a0'|'a'|'b'|'c'|'d'} [o.curve] column curve (default 'b'; EN 1993-1-5 uses αe = α + 0.09/(i/e))
 * @param {number} [o.endFixity] effective-length factor on L (default 1, simply supported at frames)
 * @returns {{rho: number, be: number, sigmaE: number, lambda: number, chi: number, sigmaU: number, Pu: number,
 *   section: object, sectionFull: object}} effective width ratio [-] and width [m]; Euler stress [Pa];
 *   slenderness [-]; reduction [-]; ultimate average stress on the gross area [Pa]; ultimate axial force of one
 *   stiffener + plate [N]; effective and gross section properties
 */
export function stiffenerColumnBuckling(o) {
  const fn = 'stiffenerColumnBuckling';
  const { L, s, tp, hw, tw, bf = 0, tf = 0, E, nu = 0.3, fy, effectiveWidth = 'en1993', curve = 'b', endFixity = 1 } = o;
  pos(L, 'L', fn); pos(s, 's', fn); pos(E, 'E', fn); pos(fy, 'fy', fn); pos(endFixity, 'endFixity', fn);
  const sl = plateSlenderness({ b: s, t: tp, E, nu, fy, k: 4 });
  let rho;
  if (effectiveWidth === 'full') rho = 1;
  else if (effectiveWidth === 'faulkner') rho = effectiveWidthRatio({ method: 'faulkner', beta: sl.beta });
  else if (effectiveWidth === 'en1993') rho = effectiveWidthRatio({ method: 'en1993-internal', lambdaP: sl.lambdaP });
  else throw new RangeError(`${fn}: effectiveWidth must be 'en1993', 'faulkner' or 'full'`);
  const be = rho * s;
  const section = stiffenerWithPlate({ be, tp, hw, tw, bf, tf });
  const sectionFull = stiffenerWithPlate({ be: s, tp, hw, tw, bf, tf });
  const Le = endFixity * L;
  const sigmaE = (Math.PI * Math.PI * E * section.I) / (section.A * Le * Le);
  const lambda = Math.sqrt(fy / sigmaE);
  const chi = columnReduction({ lambda, curve });
  const Pu = chi * fy * section.A;
  return { rho, be, sigmaE, lambda, chi, sigmaU: Pu / sectionFull.A, Pu, section, sectionFull };
}

/**
 * Elastic tripping (lateral–torsional buckling about the line of attachment) of a stiffener, classical
 * thin-walled open-section result with the stiffener rotating about its toe and no rotational restraint from
 * the plating (conservative): σ_T = (G It + π² E Iω m² / L²) / Ip, lowest for m = 1 half-wave (Hughes & Paik
 * §13.1; Timoshenko & Gere §5.5 torsional buckling with an enforced axis of rotation).
 * Section constants about the toe: flat bar It = hw tw³/3, Ip = hw³ tw/3, Iω = hw³ tw³/36;
 * tee/angle (flange treated as symmetric about the web for Iω): It = (hw tw³ + bf tf³)/3,
 * Ip = hw³ tw/3 + hw² bf tf + bf³ tf/12, Iω = tf bf³ hw²/12.
 * @param {object} o
 * @param {number} o.L span between tripping supports [m]
 * @param {number} o.hw web height [m]
 * @param {number} o.tw web thickness [m]
 * @param {number} [o.bf] flange width [m]
 * @param {number} [o.tf] flange thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} [o.fy] yield stress [Pa]; when given, the inelastic value is returned using the
 *   Johnson–Ostenfeld correction σ = fy (1 − fy/(4 σ_T)) for σ_T > fy/2
 * @returns {{sigmaT: number, sigmaTinelastic: number|null, It: number, Ip: number, Iw: number}} elastic and
 *   plasticity-corrected tripping stress [Pa]; torsion constant [m⁴], polar moment [m⁴], warping constant [m⁶]
 */
export function stiffenerTripping({ L, hw, tw, bf = 0, tf = 0, E, nu = 0.3, fy }) {
  const fn = 'stiffenerTripping';
  pos(L, 'L', fn); pos(hw, 'hw', fn); pos(tw, 'tw', fn); nonNeg(bf, 'bf', fn); nonNeg(tf, 'tf', fn); pos(E, 'E', fn); poisson(nu, fn);
  const G = E / (2 * (1 + nu));
  const flat = !(bf > 0 && tf > 0);
  const It = flat ? (hw * tw ** 3) / 3 : (hw * tw ** 3 + bf * tf ** 3) / 3;
  const Ip = flat ? (hw ** 3 * tw) / 3 : (hw ** 3 * tw) / 3 + hw * hw * bf * tf + (bf ** 3 * tf) / 12;
  const Iw = flat ? (hw ** 3 * tw ** 3) / 36 : (tf * bf ** 3 * hw * hw) / 12;
  const sigmaT = (G * It + (Math.PI * Math.PI * E * Iw) / (L * L)) / Ip;
  let inel = null;
  if (fy !== undefined) { pos(fy, 'fy', fn); inel = sigmaT <= fy / 2 ? sigmaT : fy * (1 - fy / (4 * sigmaT)); }
  return { sigmaT, sigmaTinelastic: inel, It, Ip, Iw };
}

/**
 * Simply supported orthotropic plate (Huber's equation Dx w,xxxx + 2H w,xxyy + Dy w,yyyy = q) by Navier series —
 * the smeared-stiffness model of a grillage or stiffened panel (Timoshenko & Woinowsky-Krieger §86–§87).
 * Typical smeared rigidities: Dx = E Ix / sx (girders along x, spacing sx measured along y), Dy = E Iy / sy,
 * H between the plate value and √(Dx Dy) depending on the torsional stiffness of the members
 * (H = 0 for torsionally weak open stiffeners without plate action; H = √(Dx Dy) for the isotropic limit).
 * Also returns the elastic buckling load for uniform compression Nx:
 * Nx,cr = (π²/b²) min over m of [Dx (m b/a)² + 2H + Dy (a/(m b))²].
 * @param {object} o
 * @param {number} o.a panel length in x [m]
 * @param {number} o.b panel length in y [m]
 * @param {number} o.Dx bending rigidity for curvature in x [N·m]
 * @param {number} o.Dy bending rigidity for curvature in y [N·m]
 * @param {number} o.H effective torsional rigidity [N·m]
 * @param {number} [o.q] uniform pressure [Pa] (default 0)
 * @param {number} [o.nTerms] highest harmonic (default 49)
 * @returns {{wCentre: number, mxCentre: number, myCentre: number, NxCr: number, mCr: number}}
 *   centre deflection [m]; centre moments per unit width, Poisson coupling neglected (mx = −Dx w,xx) [N·m/m];
 *   critical compressive force per unit width [N/m] and its half-wave number
 */
export function orthotropicPlate({ a, b, Dx, Dy, H, q = 0, nTerms = 49 }) {
  const fn = 'orthotropicPlate';
  pos(a, 'a', fn); pos(b, 'b', fn); pos(Dx, 'Dx', fn); pos(Dy, 'Dy', fn); nonNeg(H, 'H', fn);
  const pi = Math.PI;
  let w = 0; let mx = 0; let my = 0;
  for (let m = 1; m <= nTerms; m += 2) for (let n = 1; n <= nTerms; n += 2) {
    const am = (m * pi) / a; const bn = (n * pi) / b;
    const amn = ((16 * q) / (pi * pi * m * n)) / (Dx * am ** 4 + 2 * H * am * am * bn * bn + Dy * bn ** 4);
    const sgn = Math.sin((m * pi) / 2) * Math.sin((n * pi) / 2);
    w += amn * sgn; mx += Dx * amn * am * am * sgn; my += Dy * amn * bn * bn * sgn;
  }
  let NxCr = Infinity; let mCr = 1;
  for (let m = 1; m <= Math.max(6, Math.ceil(4 * (a / b) * Math.pow(Dy / Dx, 0.25)) + 2); m++) {
    const r = (m * b) / a;
    const v = ((pi * pi) / (b * b)) * (Dx * r * r + 2 * H + Dy / (r * r));
    if (v < NxCr) { NxCr = v; mCr = m; }
  }
  return { wCentre: w, mxCentre: mx, myCentre: my, NxCr, mCr };
}

/**
 * Grillage analysis by the beam-grid (direct stiffness) method: a rectangular grid of girders along x and
 * transverse beams along y, rigidly connected at their intersections, loaded normal to its plane. Each node
 * has a deflection and two rotations; members carry bending (Euler–Bernoulli) and St-Venant torsion.
 * The four outer edges are support lines (the members end on them).
 * @param {object} o
 * @param {number} o.a grid length in x [m]
 * @param {number} o.b grid length in y [m]
 * @param {number} o.nGirders number of girders running along x (equally spaced over b, edges excluded)
 * @param {number} o.nBeams number of beams running along y (equally spaced over a, edges excluded)
 * @param {number} o.EIx bending stiffness of one girder (along x) [N·m²]
 * @param {number} o.EIy bending stiffness of one beam (along y) [N·m²]
 * @param {number} [o.GJx] torsional stiffness of one girder [N·m²] (default 0)
 * @param {number} [o.GJy] torsional stiffness of one beam [N·m²] (default 0)
 * @param {number} [o.q] uniform pressure [Pa] lumped to the intersections by tributary area
 * @param {Array<{i: number, j: number, P: number}>} [o.pointLoads] loads [N] at intersection (beam i = 1..nBeams,
 *   girder j = 1..nGirders)
 * @param {'simple'|'clamped'} [o.support] end condition of every member on the edges (default 'simple')
 * @returns {{wMax: number, w: number[][], girderMomentMax: number, beamMomentMax: number, reactionSum: number,
 *   appliedSum: number}} deflections [m] (w[j][i] at girder j, beam i, edges included); largest bending moments
 *   in the girders and beams [N·m]; total support reaction and applied load [N]
 */
export function grillageBeamGrid(o) {
  const fn = 'grillageBeamGrid';
  const { a, b, nGirders, nBeams, EIx, EIy, GJx = 0, GJy = 0, q = 0, pointLoads = [], support = 'simple' } = o;
  pos(a, 'a', fn); pos(b, 'b', fn); pos(EIx, 'EIx', fn); pos(EIy, 'EIy', fn); nonNeg(GJx, 'GJx', fn); nonNeg(GJy, 'GJy', fn);
  if (!(Number.isInteger(nGirders) && nGirders >= 1)) throw new RangeError(`${fn}: nGirders must be a positive integer`);
  if (!(Number.isInteger(nBeams) && nBeams >= 1)) throw new RangeError(`${fn}: nBeams must be a positive integer`);
  if (support !== 'simple' && support !== 'clamped') throw new RangeError(`${fn}: support must be 'simple' or 'clamped'`);
  const nx = nBeams + 2; const ny = nGirders + 2;
  const dx = a / (nBeams + 1); const dy = b / (nGirders + 1);
  const id = (i, j) => j * nx + i;
  const nn = nx * ny;
  // Members: girders along x at j = 1..nGirders; beams along y at i = 1..nBeams.
  const mem = [];
  for (let j = 1; j <= nGirders; j++) for (let i = 0; i < nx - 1; i++) mem.push({ n1: id(i, j), n2: id(i + 1, j), L: dx, c: 1, s: 0, EI: EIx, GJ: GJx, girder: true });
  for (let i = 1; i <= nBeams; i++) for (let j = 0; j < ny - 1; j++) mem.push({ n1: id(i, j), n2: id(i, j + 1), L: dy, c: 0, s: 1, EI: EIy, GJ: GJy, girder: false });
  const ke = new Float64Array(36);
  const kl = new Float64Array(36);
  const T = new Float64Array(36);
  const localK = (m) => {
    // Local dofs per node: w, twist about the member axis, slope β = dw/ds.
    kl.fill(0);
    const { L, EI, GJ } = m;
    const k1 = (12 * EI) / L ** 3; const k2 = (6 * EI) / L ** 2; const k3 = (4 * EI) / L; const k4 = (2 * EI) / L;
    const set = (r, c, v) => { kl[r * 6 + c] = v; kl[c * 6 + r] = v; };
    set(0, 0, k1); set(0, 2, k2); set(0, 3, -k1); set(0, 5, k2);
    set(2, 2, k3); set(2, 3, -k2); set(2, 5, k4);
    set(3, 3, k1); set(3, 5, -k2); set(5, 5, k3);
    const kt = GJ / L;
    set(1, 1, kt); set(1, 4, -kt); set(4, 4, kt);
    // Global rotations (θx, θy), right-handed: twist = c θx + s θy; slope dw/ds = s θx − c θy.
    T.fill(0);
    for (let n = 0; n < 2; n++) {
      const o3 = 3 * n;
      T[o3 * 6 + o3] = 1;
      T[(o3 + 1) * 6 + o3 + 1] = m.c; T[(o3 + 1) * 6 + o3 + 2] = m.s;
      T[(o3 + 2) * 6 + o3 + 1] = m.s; T[(o3 + 2) * 6 + o3 + 2] = -m.c;
    }
  };
  const elemMatrix = (e) => {
    localK(mem[e]);
    for (let r = 0; r < 6; r++) for (let c = 0; c < 6; c++) {
      let sum = 0;
      for (let p = 0; p < 6; p++) {
        const tp = T[p * 6 + r];
        if (tp === 0) continue;
        for (let q2 = 0; q2 < 6; q2++) sum += tp * kl[p * 6 + q2] * T[q2 * 6 + c];
      }
      ke[r * 6 + c] = sum;
    }
    return ke;
  };
  const pair = new Int32Array(2);
  const elemNodes = (e) => { pair[0] = mem[e].n1; pair[1] = mem[e].n2; return pair; };
  const fixed = new Float64Array(3 * nn).fill(NaN);
  const used = new Uint8Array(nn);
  for (const m of mem) { used[m.n1] = 1; used[m.n2] = 1; }
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const n = id(i, j);
    const edge = i === 0 || i === nx - 1 || j === 0 || j === ny - 1;
    if (!used[n]) { fixed[3 * n] = 0; fixed[3 * n + 1] = 0; fixed[3 * n + 2] = 0; continue; }
    if (edge) {
      fixed[3 * n] = 0;
      // The member ending here is restrained against twist by the support line; clamped also fixes its slope.
      const alongX = j !== 0 && j !== ny - 1;
      if (alongX) { fixed[3 * n + 1] = 0; if (support === 'clamped') fixed[3 * n + 2] = 0; } else { fixed[3 * n + 2] = 0; if (support === 'clamped') fixed[3 * n + 1] = 0; }
    }
  }
  const f = new Float64Array(3 * nn);
  for (let j = 1; j <= nGirders; j++) for (let i = 1; i <= nBeams; i++) f[3 * id(i, j)] += q * dx * dy;
  for (const p of pointLoads) {
    if (!(p.i >= 1 && p.i <= nBeams && p.j >= 1 && p.j <= nGirders)) throw new RangeError(`${fn}: pointLoads index out of range`);
    f[3 * id(p.i, p.j)] += p.P;
  }
  const sys = feSystem({ nNodes: nn, dpn: 3, nElem: mem.length, elemNodes, elemMatrix, fixed, solver: 'skyline' });
  const u = sys.solve(f);
  const ku = sys.mulK(u);
  let reactionSum = 0; let appliedSum = 0; let wMax = 0;
  for (let n = 0; n < nn; n++) {
    if (sys.isFixed[3 * n]) reactionSum += ku[3 * n] - f[3 * n];
    appliedSum += f[3 * n];
    wMax = Math.max(wMax, Math.abs(u[3 * n]));
  }
  let gm = 0; let bm = 0;
  const ul = new Float64Array(6);
  for (let e = 0; e < mem.length; e++) {
    const m = mem[e];
    localK(m);
    const ug = [u[3 * m.n1], u[3 * m.n1 + 1], u[3 * m.n1 + 2], u[3 * m.n2], u[3 * m.n2 + 1], u[3 * m.n2 + 2]];
    for (let r = 0; r < 6; r++) { let s2 = 0; for (let c = 0; c < 6; c++) s2 += T[r * 6 + c] * ug[c]; ul[r] = s2; }
    for (const r of [2, 5]) {
      let mom = 0;
      for (let c = 0; c < 6; c++) mom += kl[r * 6 + c] * ul[c];
      if (m.girder) gm = Math.max(gm, Math.abs(mom)); else bm = Math.max(bm, Math.abs(mom));
    }
  }
  const w = [];
  for (let j = 0; j < ny; j++) { const row = []; for (let i = 0; i < nx; i++) row.push(u[3 * id(i, j)]); w.push(row); }
  return { wMax, w, girderMomentMax: gm, beamMomentMax: bm, reactionSum, appliedSum };
}

/* ------------------------------------------------------------------------------------------------------------
 * Hull girder
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Elastic hull-girder (or any thin-walled built-up beam) section properties from lumped longitudinal elements.
 * @param {object} o
 * @param {Array<{area: number, z: number, i0?: number, E?: number}>} o.elements longitudinally effective
 *   elements: area [m²], height of centroid above baseline z [m], own second moment i0 [m⁴] (default 0),
 *   Young's modulus [Pa] for transformed sections (default: all equal)
 * @returns {{A: number, zNA: number, I: number, zTop: number, zBottom: number, Ztop: number, Zbottom: number}}
 *   area [m²] (transformed to the first element's modulus when moduli differ), neutral axis height [m], second
 *   moment about the neutral axis [m⁴], extreme fibre heights [m], section moduli at top (deck) and bottom
 *   (keel) [m³]
 */
export function hullGirderSection({ elements }) {
  const fn = 'hullGirderSection';
  if (!Array.isArray(elements) || elements.length === 0) throw new RangeError(`${fn}: elements must be a non-empty array`);
  const E0 = elements[0].E || 1;
  let A = 0; let Az = 0; let zTop = -Infinity; let zBottom = Infinity;
  for (const e of elements) {
    pos(e.area, 'elements[].area', fn);
    if (!Number.isFinite(e.z)) throw new RangeError(`${fn}: elements[].z must be a number`);
    const n = (e.E || E0) / E0;
    A += n * e.area; Az += n * e.area * e.z;
    if (e.z > zTop) zTop = e.z;
    if (e.z < zBottom) zBottom = e.z;
  }
  const zNA = Az / A;
  let I = 0;
  for (const e of elements) { const n = (e.E || E0) / E0; I += n * ((e.i0 || 0) + e.area * (e.z - zNA) ** 2); }
  return { A, zNA, I, zTop, zBottom, Ztop: I / (zTop - zNA), Zbottom: I / (zNA - zBottom) };
}

/**
 * Ultimate hull-girder bending moment by a Smith-type progressive-collapse (incremental curvature) analysis.
 * The cross-section is divided into elements acting independently; plane sections remain plane. Each element
 * follows an elastic–perfectly-plastic curve in tension (yield at fy) and, in compression, is elastic up to
 * its ultimate compressive strength sigmaU (plate/stiffener buckling strength, supplied per element, e.g. from
 * `stiffenerColumnBuckling`) and then holds that stress (no post-ultimate unloading — a simplification that
 * makes the peak moment an upper estimate compared with full Smith load-shortening curves). For every
 * curvature the neutral axis is found from axial-force equilibrium.
 * @param {object} o
 * @param {Array<{area: number, z: number, fy: number, sigmaU?: number}>} o.elements area [m²], height z [m],
 *   yield stress [Pa], ultimate compressive stress [Pa] (default fy: "hard corner")
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nSteps] curvature increments per bending sense (default 200)
 * @param {number} [o.maxStrainFactor] largest extreme-fibre strain as a multiple of the yield strain (default 5)
 * @returns {{sagging: {Mu: number, curvature: number[], moment: number[]}, hogging: {Mu: number,
 *   curvature: number[], moment: number[]}, Mp: number, MyieldFirst: number}} ultimate moments [N·m] (sagging:
 *   top in compression, positive values) with their moment–curvature curves [1/m, N·m]; fully plastic moment
 *   ignoring buckling [N·m]; first-yield moment of the elastic section [N·m]
 */
export function hullGirderUltimate({ elements, E, nSteps = 200, maxStrainFactor = 5 }) {
  const fn = 'hullGirderUltimate';
  pos(E, 'E', fn);
  const sec = hullGirderSection({ elements });
  let fyMax = 0; let fyMin = Infinity;
  for (const e of elements) {
    pos(e.fy, 'elements[].fy', fn);
    if (e.sigmaU !== undefined && !(e.sigmaU > 0 && e.sigmaU <= e.fy * (1 + 1e-12))) throw new RangeError(`${fn}: elements[].sigmaU must lie in (0, fy]`);
    fyMax = Math.max(fyMax, e.fy); fyMin = Math.min(fyMin, e.fy);
  }
  const depth = sec.zTop - sec.zBottom;
  if (!(depth > 0)) throw new RangeError(`${fn}: elements must span a finite depth`);
  const stress = (e, eps) => {
    // eps > 0 tension.
    const s = E * eps;
    if (s >= 0) return Math.min(s, e.fy);
    return Math.max(s, -(e.sigmaU === undefined ? e.fy : e.sigmaU));
  };
  const run = (sense) => {
    // sense +1: sagging (compression at the top). Strain at z: eps = -sense · κ · (z − z0).
    const curvature = [0]; const moment = [0];
    const kMax = (maxStrainFactor * fyMax) / E / (0.5 * depth);
    let Mu = 0;
    for (let s = 1; s <= nSteps; s++) {
      const kappa = (kMax * s) / nSteps;
      const force = (z0) => { let F = 0; for (const e of elements) F += stress(e, -sense * kappa * (e.z - z0)) * e.area; return F; };
      // Axial force is monotonic in z0; bracket on an extended depth and bisect.
      let lo = sec.zBottom - depth; let hi = sec.zTop + depth;
      let flo = force(lo);
      for (let it = 0; it < 80; it++) {
        const mid = 0.5 * (lo + hi);
        const fm = force(mid);
        if ((fm > 0) === (flo > 0)) { lo = mid; flo = fm; } else hi = mid;
      }
      const z0 = 0.5 * (lo + hi);
      let M = 0;
      for (const e of elements) M += stress(e, -sense * kappa * (e.z - z0)) * e.area * (e.z - z0);
      M *= -sense;
      curvature.push(kappa); moment.push(M);
      if (M > Mu) Mu = M;
    }
    return { Mu, curvature, moment };
  };
  // Plastic moment: equal-force axis with every element at yield.
  const sorted = elements.slice().sort((p, q2) => p.z - q2.z);
  let tot = 0;
  for (const e of sorted) tot += e.area * e.fy;
  let acc = 0; let zp = sorted[0].z;
  for (const e of sorted) { acc += e.area * e.fy; if (acc >= tot / 2) { zp = e.z; break; } }
  let Mp = 0;
  for (const e of elements) Mp += e.area * e.fy * Math.abs(e.z - zp);
  return { sagging: run(1), hogging: run(-1), Mp, MyieldFirst: fyMin * Math.min(sec.Ztop, sec.Zbottom) };
}

/* ------------------------------------------------------------------------------------------------------------
 * Deck plate under wheel / patch load
 * ---------------------------------------------------------------------------------------------------------- */

/**
 * Elastic response of a deck plate panel between stiffeners under a rectangular wheel or patch load.
 * The panel is analysed as a Kirchhoff plate with all edges simply supported (Navier; conservative for
 * deflection and field moment) or all edges clamped (Ritz; continuous plating under a single loaded panel).
 * Reports the largest bending stress found at the load centre and at the edge mid-points and its utilisation
 * against yield; plastic (yield-line/membrane) reserve is not included.
 * @param {object} o
 * @param {number} o.a panel length (between transverse frames) [m]
 * @param {number} o.b panel width (stiffener spacing) [m]
 * @param {number} o.t plate thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} o.P total patch load [N]
 * @param {number} o.u patch length in x [m]
 * @param {number} o.v patch length in y [m]
 * @param {number} [o.x0] patch centre x [m] (default a/2)
 * @param {number} [o.y0] patch centre y [m] (default b/2)
 * @param {'ss'|'clamped'} [o.edges] edge condition (default 'ss')
 * @param {number} [o.fy] yield stress [Pa] for the utilisation
 * @returns {{w: number, mx: number, my: number, edgeMomentMax: number, sigmaMax: number, vonMises: number,
 *   utilisation: number|null}} deflection under the load centre [m]; moments there [N·m/m]; largest edge
 *   moment magnitude [N·m/m]; largest surface bending stress [Pa]; von Mises surface stress under the load
 *   [Pa]; sigma/fy [-]
 */
export function deckPatchLoad(o) {
  const fn = 'deckPatchLoad';
  const { a, b, t, E, nu = 0.3, P, u, v, edges = 'ss', fy } = o;
  pos(a, 'a', fn); pos(b, 'b', fn); pos(t, 't', fn);
  if (edges !== 'ss' && edges !== 'clamped') throw new RangeError(`${fn}: edges must be 'ss' or 'clamped'`);
  const x0 = o.x0 === undefined ? a / 2 : o.x0;
  const y0 = o.y0 === undefined ? b / 2 : o.y0;
  const load = { type: 'patch', P, x0, y0, u, v };
  let r; let edgeM = 0;
  if (edges === 'ss') r = navierPlate({ a, b, t, E, nu, load, x: x0, y: y0, nTerms: 199 });
  else {
    r = ritzPlate({ a, b, t, E, nu, load, x: x0, y: y0, edges: { x0: 'C', x1: 'C', y0: 'C', y1: 'C' }, nBasis: 12 });
    edgeM = Math.max(...Object.values(r.edgeMoments).map(Math.abs));
  }
  const s6 = 6 / (t * t);
  const sx = s6 * r.mx; const sy = s6 * r.my;
  const sigmaMax = Math.max(Math.abs(sx), Math.abs(sy), s6 * edgeM);
  const vm = Math.sqrt(sx * sx + sy * sy - sx * sy);
  let util = null;
  if (fy !== undefined) { pos(fy, 'fy', fn); util = Math.max(vm, s6 * edgeM * Math.sqrt(1 - nu + nu * nu)) / fy; }
  return { w: r.w, mx: r.mx, my: r.my, edgeMomentMax: edgeM, sigmaMax, vonMises: vm, utilisation: util };
}
