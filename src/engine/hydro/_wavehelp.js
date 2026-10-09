// Private helpers shared by the hydro engines: linear (Airy) wave theory, JONSWAP spectrum, spectral moments,
// small complex linear solves and a box-bounded least-squares solver. Not part of the public engine API.
//
// Conventions used throughout src/engine/hydro:
//   * z is up, z = 0 at the still-water level, the seabed is at z = -depth.
//   * A regular wave of amplitude a travelling towards heading beta (measured from +x towards +y) has
//     elevation eta = a cos(k (x cos beta + y sin beta) - omega t).
//   * Complex amplitudes use the exp(+i omega t) convention: q(t) = Re{ q_hat exp(i omega t) }, so the
//     incident elevation amplitude at (x, y) is a exp(-i k (x cos beta + y sin beta)). A positive phase is a lead
//     relative to the wave crest passing the origin.

/** Standard gravity, m/s². */
export const G = 9.80665;

/**
 * Wave number from the linear dispersion relation omega² = g k tanh(k d).
 * Reference: Airy theory, e.g. DNV-RP-C205 (2019) §3.2. Newton iteration, relative tolerance 1e-14.
 * @param {number} omega angular frequency, rad/s (> 0)
 * @param {number} [depth=Infinity] water depth, m (> 0)
 * @param {number} [g=9.80665] gravity, m/s²
 * @returns {number} wave number k, rad/m
 */
export function dispersion(omega, depth = Infinity, g = G) {
  if (!(omega > 0)) throw new RangeError('omega must be > 0');
  if (!(depth > 0)) throw new RangeError('depth must be > 0');
  const k0 = (omega * omega) / g;
  if (!Number.isFinite(depth) || k0 * depth > 40) return k0;
  let k = k0 / Math.sqrt(Math.tanh(k0 * depth));
  for (let i = 0; i < 60; i++) {
    const th = Math.tanh(k * depth);
    const f = g * k * th - omega * omega;
    const df = g * th + g * k * depth * (1 - th * th);
    const dk = f / df;
    k -= dk;
    if (Math.abs(dk) < 1e-14 * k) break;
  }
  return k;
}

/**
 * Depth-decay ratios of linear theory, evaluated without overflow in deep water.
 * @param {number} k wave number, rad/m
 * @param {number} z elevation, m (<= 0 normally)
 * @param {number} depth water depth, m (may be Infinity)
 * @param {{chsh:number, shsh:number, chch:number}} out receives cosh k(z+d)/sinh kd, sinh k(z+d)/sinh kd,
 *   cosh k(z+d)/cosh kd (all dimensionless)
 * @returns {{chsh:number, shsh:number, chch:number}} out
 */
export function depthRatios(k, z, depth, out) {
  const e1 = Math.exp(k * z);
  if (!Number.isFinite(depth)) {
    out.chsh = e1; out.shsh = e1; out.chch = e1;
    return out;
  }
  const e2 = Math.exp(-k * (z + 2 * depth));
  const e3 = Math.exp(-2 * k * depth);
  out.chsh = (e1 + e2) / (1 - e3);
  out.shsh = (e1 - e2) / (1 - e3);
  out.chch = (e1 + e2) / (1 + e3);
  return out;
}

/**
 * Pre-computes the constants of a regular Airy wave.
 * @param {{H:number, T:number, depth?:number, headingDeg?:number, phaseDeg?:number}} wave height m, period s,
 *   depth m, propagation heading deg from +x, phase lag deg
 * @param {number} [g] gravity, m/s²
 * @returns {{a:number, omega:number, k:number, depth:number, cb:number, sb:number, phase:number, g:number}}
 */
export function prepWave(wave, g = G) {
  if (!wave) throw new RangeError('wave is required');
  const { H, T, depth = Infinity, headingDeg = 0, phaseDeg = 0 } = wave;
  if (!(H >= 0)) throw new RangeError('wave.H must be >= 0');
  if (!(T > 0)) throw new RangeError('wave.T must be > 0');
  const omega = (2 * Math.PI) / T;
  const b = (headingDeg * Math.PI) / 180;
  return {
    a: H / 2, omega, k: dispersion(omega, depth, g), depth, g,
    cb: Math.cos(b), sb: Math.sin(b), phase: (phaseDeg * Math.PI) / 180,
  };
}

const _r = { chsh: 0, shsh: 0, chch: 0 };

/**
 * Free-surface elevation of a prepared regular wave.
 * @returns {number} eta, m
 */
export function airyEta(w, x, y, t) {
  return w.a * Math.cos(w.k * (x * w.cb + y * w.sb) - w.omega * t - w.phase);
}

/**
 * Airy particle kinematics at a point. `zEff` is the elevation at which the depth functions are evaluated
 * (the caller applies any stretching).
 * @param {object} w prepared wave (prepWave)
 * @param {number} x m
 * @param {number} y m
 * @param {number} zEff m
 * @param {number} t s
 * @param {{ux:number,uy:number,uz:number,ax:number,ay:number,az:number,eta:number,pd:number}} out velocities m/s,
 *   accelerations m/s², elevation m, dynamic pressure per unit density pd (m²/s², multiply by rho for Pa)
 */
export function airyPoint(w, x, y, zEff, t, out) {
  const psi = w.k * (x * w.cb + y * w.sb) - w.omega * t - w.phase;
  const c = Math.cos(psi);
  const s = Math.sin(psi);
  depthRatios(w.k, zEff, w.depth, _r);
  const uh = w.a * w.omega * _r.chsh * c;
  const ah = w.a * w.omega * w.omega * _r.chsh * s;
  out.ux = uh * w.cb; out.uy = uh * w.sb;
  out.uz = w.a * w.omega * _r.shsh * s;
  out.ax = ah * w.cb; out.ay = ah * w.sb;
  out.az = -w.a * w.omega * w.omega * _r.shsh * c;
  out.eta = w.a * c;
  out.pd = w.g * w.a * _r.chch * c;
  return out;
}

/**
 * JONSWAP spectral density, DNV-RP-C205 (2019) §3.5.5:
 * S(w) = A_g (5/16) Hs² wp⁴ w⁻⁵ exp(-1.25 (w/wp)⁻⁴) gamma^exp(-0.5 ((w - wp)/(sigma wp))²),
 * A_g = 1 - 0.287 ln(gamma), sigma = 0.07 (w <= wp) or 0.09 (w > wp). gamma = 1 gives Pierson–Moskowitz.
 * @param {number} omega rad/s
 * @param {{Hs:number, Tp:number, gamma?:number}} sea Hs m, Tp s, peak enhancement (1..7)
 * @returns {number} S, m²·s/rad
 */
export function jonswap(omega, { Hs, Tp, gamma = 3.3 }) {
  if (!(Hs >= 0)) throw new RangeError('Hs must be >= 0');
  if (!(Tp > 0)) throw new RangeError('Tp must be > 0');
  if (!(gamma >= 1)) throw new RangeError('gamma must be >= 1');
  if (!(omega > 0)) return 0;
  const wp = (2 * Math.PI) / Tp;
  const r = omega / wp;
  const sigma = omega <= wp ? 0.07 : 0.09;
  const pm = (5 / 16) * Hs * Hs * wp ** 4 * omega ** -5 * Math.exp(-1.25 * r ** -4);
  const peak = gamma ** Math.exp(-0.5 * ((r - 1) / sigma) ** 2);
  return (1 - 0.287 * Math.log(gamma)) * pm * peak;
}

/**
 * Spectral moments m_n = ∫ w^n S(w) dw by the trapezoidal rule.
 * @param {ArrayLike<number>} omega rad/s, ascending
 * @param {ArrayLike<number>} S spectral density, unit²·s/rad
 * @returns {{m0:number, m1:number, m2:number, m4:number}} unit²·(rad/s)^n
 */
export function spectrumMoments(omega, S) {
  let m0 = 0; let m1 = 0; let m2 = 0; let m4 = 0;
  for (let i = 1; i < omega.length; i++) {
    const w0 = omega[i - 1]; const w1 = omega[i];
    const h = 0.5 * (w1 - w0);
    m0 += h * (S[i - 1] + S[i]);
    m1 += h * (S[i - 1] * w0 + S[i] * w1);
    m2 += h * (S[i - 1] * w0 * w0 + S[i] * w1 * w1);
    m4 += h * (S[i - 1] * w0 ** 4 + S[i] * w1 ** 4);
  }
  return { m0, m1, m2, m4 };
}

/**
 * Solves the complex system (Ar + i Ai) x = (br + i bi) by Gaussian elimination with partial pivoting.
 * Inputs are overwritten. Intended for the 6×6 systems of rigid-body response.
 * @param {number[][]} Ar real part, n×n
 * @param {number[][]} Ai imaginary part, n×n
 * @param {number[]} br real right-hand side
 * @param {number[]} bi imaginary right-hand side
 * @returns {{re:number[], im:number[]}}
 */
export function csolve(Ar, Ai, br, bi) {
  const n = br.length;
  for (let c = 0; c < n; c++) {
    let p = c; let best = -1;
    for (let r = c; r < n; r++) {
      const m = Ar[r][c] * Ar[r][c] + Ai[r][c] * Ai[r][c];
      if (m > best) { best = m; p = r; }
    }
    if (!(best > 0)) throw new Error('Singular dynamic stiffness matrix (unrestrained degree of freedom at zero frequency?)');
    if (p !== c) {
      let t = Ar[p]; Ar[p] = Ar[c]; Ar[c] = t;
      t = Ai[p]; Ai[p] = Ai[c]; Ai[c] = t;
      let s = br[p]; br[p] = br[c]; br[c] = s;
      s = bi[p]; bi[p] = bi[c]; bi[c] = s;
    }
    const pr = Ar[c][c]; const pi = Ai[c][c];
    const den = pr * pr + pi * pi;
    for (let r = c + 1; r < n; r++) {
      const ar = Ar[r][c]; const ai = Ai[r][c];
      if (ar === 0 && ai === 0) continue;
      const fr = (ar * pr + ai * pi) / den;
      const fi = (ai * pr - ar * pi) / den;
      for (let j = c; j < n; j++) {
        Ar[r][j] -= fr * Ar[c][j] - fi * Ai[c][j];
        Ai[r][j] -= fr * Ai[c][j] + fi * Ar[c][j];
      }
      br[r] -= fr * br[c] - fi * bi[c];
      bi[r] -= fr * bi[c] + fi * br[c];
    }
  }
  const re = new Array(n).fill(0);
  const im = new Array(n).fill(0);
  for (let r = n - 1; r >= 0; r--) {
    let sr = br[r]; let si = bi[r];
    for (let j = r + 1; j < n; j++) {
      sr -= Ar[r][j] * re[j] - Ai[r][j] * im[j];
      si -= Ar[r][j] * im[j] + Ai[r][j] * re[j];
    }
    const pr = Ar[r][r]; const pi = Ai[r][r];
    const den = pr * pr + pi * pi;
    re[r] = (sr * pr + si * pi) / den;
    im[r] = (si * pr - sr * pi) / den;
  }
  return { re, im };
}

/**
 * Box-bounded linear least squares: minimise ½|A x - b|² + ½ lambda |x|² subject to lo <= x <= hi, by cyclic
 * projected coordinate descent on the normal equations (convergent for a positive semi-definite Hessian with a
 * positive diagonal). Intended for small problems (tens of unknowns).
 * @param {number[][]} A m×n
 * @param {number[]} b length m
 * @param {number[]} lo lower bounds, length n
 * @param {number[]} hi upper bounds, length n
 * @param {{lambda?:number, maxIter?:number, tol?:number, x0?:number[]}} [opts]
 * @returns {{x:number[], residual:number, iterations:number, converged:boolean}} residual = |A x - b|
 */
export function boundedLeastSquares(A, b, lo, hi, { lambda = 0, maxIter = 20000, tol = 1e-12, x0 } = {}) {
  const m = A.length;
  const n = lo.length;
  const Gm = [];
  const h = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    Gm.push(new Array(n).fill(0));
    for (let j = 0; j <= i; j++) {
      let s = 0;
      for (let r = 0; r < m; r++) s += A[r][i] * A[r][j];
      Gm[i][j] = s; Gm[j][i] = s;
    }
    Gm[i][i] += lambda;
    for (let r = 0; r < m; r++) h[i] += A[r][i] * b[r];
  }
  const x = new Array(n);
  for (let i = 0; i < n; i++) {
    if (lo[i] > hi[i]) throw new RangeError('lower bound exceeds upper bound');
    x[i] = Math.min(hi[i], Math.max(lo[i], x0 ? x0[i] : 0));
  }
  let scale = 0;
  for (let i = 0; i < n; i++) scale = Math.max(scale, Math.abs(hi[i]) < Infinity ? Math.abs(hi[i]) : 0, Math.abs(lo[i]) < Infinity ? Math.abs(lo[i]) : 0);
  if (scale === 0) scale = 1;
  let it = 0; let converged = false;
  for (; it < maxIter; it++) {
    let change = 0;
    for (let j = 0; j < n; j++) {
      if (!(Gm[j][j] > 0)) continue;
      let s = h[j];
      for (let k2 = 0; k2 < n; k2++) if (k2 !== j) s -= Gm[j][k2] * x[k2];
      const xn = Math.min(hi[j], Math.max(lo[j], s / Gm[j][j]));
      change = Math.max(change, Math.abs(xn - x[j]));
      x[j] = xn;
    }
    if (change <= tol * scale) { converged = true; it++; break; }
  }
  let res = 0;
  for (let r = 0; r < m; r++) {
    let s = -b[r];
    for (let j = 0; j < n; j++) s += A[r][j] * x[j];
    res += s * s;
  }
  return { x, residual: Math.sqrt(res), iterations: it, converged };
}

/** Throws RangeError unless v is a finite number satisfying v > 0. */
export function needPos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name} must be a finite number > 0`);
}

/** Throws RangeError unless v is a length-3 numeric vector. */
export function needVec3(v, name) {
  if (!v || v.length !== 3 || !Number.isFinite(v[0] + v[1] + v[2])) throw new RangeError(`${name} must be a finite [x, y, z] vector`);
}
