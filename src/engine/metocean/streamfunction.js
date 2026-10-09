// Fourier approximation ("stream function") theory for steady nonlinear waves on a horizontal bed.
//
// Formulation: Rienecker, M.M. & Fenton, J.D. (1981) "A Fourier approximation method for steady water waves",
// J. Fluid Mech. 104, 119–137; Fenton, J.D. (1988) "The numerical solution of steady water wave problems",
// Computers & Geosciences 14(3), 357–368. Equivalent in content to Dean's (1965) stream-function theory.
//
// In the frame moving with the wave the flow is steady, with stream function (y measured up from the bed)
//   ψ(x, y) = B0 y + Σ_{j=1..N} B_j sinh(j k y)/cosh(j k d) cos(j k x).
// The unknowns are the surface ordinates η_m at N+1 points over half a wavelength, B0..B_N, k, c, the volume
// flux Q and the Bernoulli constant R. They satisfy ψ = −Q and ½(u² + w²) + g η = R at every collocation point,
// together with the mean-depth, wave-height, period (or length) and mean-current conditions. The nonlinear system
// is solved by Newton's method with an analytical Jacobian, approaching the target height in steps.

import { solveDense } from '../math/linalg.js';

const G = 9.80665;
const TWO_PI = 2 * Math.PI;

function linearK(T, d, U, g) {
  // linear dispersion with uniform current, used only for scaling and the starting guess
  const om = TWO_PI / T;
  let k = (om * om) / g;
  k /= Math.sqrt(Math.tanh(k * d));
  for (let i = 0; i < 60; i++) {
    const th = Math.tanh(k * d);
    const sig = Math.sqrt(g * k * th);
    const f = sig + k * U - om;
    const cg = (0.5 * sig / k) * (1 + (k * d * (1 - th * th)) / th);
    const dk = f / (cg + U);
    let kn = k - dk;
    if (!(kn > 0)) kn = 0.5 * k;
    if (Math.abs(kn - k) < 1e-14 * k) { k = kn; break; }
    k = kn;
  }
  return k;
}

// Residual vector F and Jacobian J of the Rienecker–Fenton system in units where g = 1.
function assemble(z, p, F, J) {
  const { N, D, H, cosT, sinT, tc, ts } = p;
  const M = 2 * N + 6;
  const iB = N + 1; const ik = 2 * N + 2; const ic = ik + 1; const iQ = ik + 2; const iR = ik + 3;
  const k = z[ik]; const c = z[ic]; const Q = z[iQ]; const R = z[iR];
  for (let i = 0; i < M; i++) J[i].fill(0);
  for (let m = 0; m <= N; m++) {
    const eta = z[m];
    let psi = z[iB] * eta; let u = z[iB]; let v = 0;
    let psiK = 0; let uK = 0; let vK = 0; let uE = 0; let vE = 0;
    const rk = J[m]; const rd = J[N + 1 + m];
    for (let j = 1; j <= N; j++) {
      const jk = j * k;
      const e = Math.exp(jk * (eta - D));
      const r = Math.exp(-2 * jk * eta);
      const q = Math.exp(-2 * jk * D);
      const S = (e * (1 - r)) / (1 + q);
      const C = (e * (1 + r)) / (1 + q);
      const th = (1 - q) / (1 + q);
      const dS = j * eta * C - j * D * S * th;
      const dC = j * eta * S - j * D * C * th;
      const cj = cosT[m * (N + 1) + j]; const sj = sinT[m * (N + 1) + j];
      const Bj = z[iB + j];
      psi += Bj * S * cj;
      u += jk * Bj * C * cj;
      v += jk * Bj * S * sj;
      psiK += Bj * cj * dS;
      uK += j * Bj * cj * (C + k * dC);
      vK += j * Bj * sj * (S + k * dS);
      uE += jk * jk * Bj * S * cj;
      vE += jk * jk * Bj * C * sj;
      rk[iB + j] = S * cj;
      tc[j] = jk * C * cj;
      ts[j] = jk * S * sj;
    }
    F[m] = psi + Q;
    rk[m] = u; rk[iB] = eta; rk[ik] = psiK; rk[iQ] = 1;
    F[N + 1 + m] = 0.5 * (u * u + v * v) + eta - R;
    rd[m] = u * uE + v * vE + 1; rd[iB] = u; rd[ik] = u * uK + v * vK; rd[iR] = -1;
    for (let j = 1; j <= N; j++) rd[iB + j] = u * tc[j] + v * ts[j];
  }
  // mean depth (trapezoidal rule over the half wavelength is spectrally accurate for the periodic profile)
  let mean = 0.5 * (z[0] + z[N]);
  for (let m = 1; m < N; m++) mean += z[m];
  let row = 2 * N + 2;
  F[row] = mean / N - D;
  J[row][0] = 0.5 / N; J[row][N] = 0.5 / N;
  for (let m = 1; m < N; m++) J[row][m] = 1 / N;
  // wave height
  row++;
  F[row] = z[0] - z[N] - H;
  J[row][0] = 1; J[row][N] = -1;
  // period or wavelength
  row++;
  if (p.kFixed !== null) { F[row] = k - p.kFixed; J[row][ik] = 1; }
  else { F[row] = k * c * p.T - TWO_PI; J[row][ik] = c * p.T; J[row][ic] = k * p.T; }
  // mean current: Eulerian c + B0 = cE (B0 = −mean speed in the wave frame), or Stokes c − Q/d = cS
  row++;
  if (p.stokes) { F[row] = c - Q / D - p.cur; J[row][ic] = 1; J[row][iQ] = -1 / D; }
  else { F[row] = c + z[iB] - p.cur; J[row][ic] = 1; J[row][iB] = 1; }
}

function newton(z, p, tol, maxIter) {
  const M = z.length;
  const F = new Array(M).fill(0);
  const J = Array.from({ length: M }, () => new Array(M).fill(0));
  let res = Infinity;
  let it = 0;
  const N = p.N;
  for (; it < maxIter; it++) {
    assemble(z, p, F, J);
    res = 0;
    for (let i = 0; i < M; i++) { const a = Math.abs(F[i]); if (!(a <= res)) res = a; }
    if (!Number.isFinite(res)) return { ok: false, res, it };
    if (res < tol) return { ok: true, res, it };
    const rhs = F.map((v) => -v);
    let dz;
    try { dz = solveDense(J.map((r) => r.slice()), rhs); } catch (e) { return { ok: false, res, it }; }
    // keep the surface above the bed and the wavenumber positive
    let lam = 1;
    for (let tries = 0; tries < 30; tries++) {
      let bad = !(z[2 * N + 2] + lam * dz[2 * N + 2] > 0);
      for (let m = 0; m <= N && !bad; m++) if (!(z[m] + lam * dz[m] > 0)) bad = true;
      if (!bad) break;
      lam *= 0.5;
    }
    for (let i = 0; i < M; i++) {
      if (!Number.isFinite(dz[i])) return { ok: false, res, it };
      z[i] += lam * dz[i];
    }
  }
  assemble(z, p, F, J);
  res = 0;
  for (let i = 0; i < M; i++) { const a = Math.abs(F[i]); if (!(a <= res)) res = a; }
  return { ok: res < tol, res, it };
}

// Physical admissibility of a converged solution: beyond the breaking limit the algebraic system still has
// roots, but they have a non-monotonic surface or particles overtaking the crest. Those are rejected.
function admissible(z, p) {
  const { N, D, cosT } = p;
  const iB = N + 1; const k = z[2 * N + 2];
  const tolEta = 1e-9 * p.H;
  for (let m = 0; m <= N; m++) {
    if (m > 0 && z[m] > z[m - 1] + tolEta) return false;
    const eta = z[m];
    let u = z[iB];
    for (let j = 1; j <= N; j++) {
      const jk = j * k;
      const C = (Math.exp(jk * (eta - D)) * (1 + Math.exp(-2 * jk * eta))) / (1 + Math.exp(-2 * jk * D));
      u += jk * z[iB + j] * C * cosT[m * (N + 1) + j];
    }
    // in the wave frame the fluid must everywhere move backwards relative to the crest
    if (!(u < 0)) return false;
  }
  return true;
}

function tryHeightSteps(p0, Htarget, nSteps, tol) {
  const { N, D } = p0;
  const M = 2 * N + 6;
  const iB = N + 1; const ik = 2 * N + 2;
  const p = { ...p0, H: Htarget / nSteps };
  // linear starting solution at the first (small) height
  const k0 = p.kStart;
  const sig = Math.sqrt(k0 * Math.tanh(k0 * D));
  const cRel = sig / k0;
  let z = new Array(M).fill(0);
  for (let m = 0; m <= N; m++) z[m] = D + 0.5 * p.H * p.cosT[m * (N + 1) + 1];
  z[iB] = -cRel;
  z[iB + 1] = (0.5 * p.H) / sig;
  z[ik] = k0;
  z[ik + 1] = cRel + p.cur; // phase speed (to first order both current definitions coincide)
  z[ik + 2] = cRel * D;
  z[ik + 3] = 0.5 * cRel * cRel + D;
  let prev = null;
  let iterations = 0;
  let res = 0;
  for (let s = 1; s <= nSteps; s++) {
    p.H = (Htarget * s) / nSteps;
    let guess;
    if (s === 1) guess = z;
    else if (prev === null) { guess = z.map((v, i) => (i <= N ? D + ((v - D) * s) / (s - 1) : v)); }
    else guess = z.map((v, i) => 2 * v - prev[i]);
    const out = newton(guess, p, tol, 40);
    iterations += out.it;
    if (!out.ok) return null;
    prev = s === 1 ? null : z;
    z = guess;
    res = out.res;
  }
  if (!admissible(z, p)) return null;
  return { z, iterations, res };
}

/**
 * Solve a steady nonlinear wave by the Fourier (stream-function) method of order N.
 * Reference: Rienecker & Fenton (1981) J. Fluid Mech. 104; Fenton (1988) Computers & Geosciences 14(3);
 * Dean (1965) J. Geophys. Res. 70(18). Recommended by API RP 2A-WSD §5.3.3 and ISO 19901-1 §A.8.4 for steep
 * and shallow-water regular waves.
 * Validity: finite depth, non-breaking waves. The highest waves need N ≥ 20 and converge slowly because the
 * crest tends to a 120° corner; above roughly 97 % of the limiting height the iteration may fail to converge
 * (check report.dynamicResidual for the truncation error). Deep water is obtained with d ≥ about 1.5 L.
 * Throws Error when Newton's method does not converge or converges to a non-physical root (surface not monotonic
 * between crest and trough, or particle speed at the surface reaching the wave speed), i.e. at or beyond breaking.
 * @param {object} o
 * @param {number} o.H wave height, m
 * @param {number} o.d still-water (mean) depth, m
 * @param {number} [o.T] wave period in the fixed frame, s (give T or L)
 * @param {number} [o.L] wavelength, m
 * @param {number} [o.order=10] number of Fourier terms N, 1..30
 * @param {number} [o.currentEulerian=0] Eulerian mean current in the wave direction, m/s
 * @param {number} [o.currentStokes] depth-averaged mass-transport (Stokes drift) velocity, m/s; when given it
 *   replaces currentEulerian (use 0 for a closed flume)
 * @param {number} [o.heightSteps] number of height steps; default chosen from the relative height and doubled
 *   automatically on failure
 * @param {number} [o.tol=1e-12] convergence tolerance on the dimensionless residuals
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{H:number,d:number,g:number,order:number,k:number,L:number,T:number,c:number,omega:number,
 *   currentEulerian:number,currentStokes:number,B0:number,B:number[],Q:number,R:number,etaNodes:number[],
 *   etaCoef:number[],crest:number,trough:number,
 *   report:{converged:boolean,iterations:number,heightSteps:number,residual:number,kinematicResidual:number,
 *   dynamicResidual:number}}}
 *   k rad/m; L m; T s; c m/s; omega rad/s; B0 m/s (minus the mean speed in the wave frame); B[j] m²/s
 *   (B[0] unused = 0); Q volume flux per unit width in the wave frame, m²/s; R Bernoulli constant with the bed as
 *   datum, m²/s²; etaNodes surface elevation above still water at phases mπ/N, m; etaCoef cosine-series
 *   coefficients of the surface, m; crest, trough elevations above still water, m. Report: residual is the largest
 *   dimensionless equation residual at the collocation points; kinematicResidual = max |ψ + Q|/Q and
 *   dynamicResidual = max |½(u²+w²) + g y − R|/(g H) evaluated on the interpolated surface midway between
 *   collocation points (a measure of truncation error).
 */
export function streamFunctionWave({ H, d, T, L, order = 10, currentEulerian = 0, currentStokes, heightSteps,
  tol = 1e-12, g = G }) {
  if (!(H > 0) || !Number.isFinite(H)) throw new RangeError('H must be a positive finite number');
  if (!(d > 0) || !Number.isFinite(d)) throw new RangeError('d must be a positive finite number');
  if (!Number.isInteger(order) || order < 1 || order > 30) throw new RangeError('order must be an integer from 1 to 30');
  const hasL = L !== undefined && L !== null;
  if (hasL) { if (!(L > 0) || !Number.isFinite(L)) throw new RangeError('L must be a positive finite number'); }
  else if (!(T > 0) || !Number.isFinite(T)) throw new RangeError('T must be a positive finite number (or give L)');
  const stokes = currentStokes !== undefined && currentStokes !== null;
  const cur = stokes ? currentStokes : currentEulerian;
  if (!Number.isFinite(cur)) throw new RangeError('current must be a finite number');
  const N = order;
  // length scale 1/k (linear estimate) keeps every unknown of order one in deep and shallow water alike
  const kRef = hasL ? TWO_PI / L : linearK(T, d, cur, g);
  const ell = 1 / kRef;
  const vel = Math.sqrt(g * ell);
  const cosT = new Float64Array((N + 1) * (N + 1));
  const sinT = new Float64Array((N + 1) * (N + 1));
  for (let m = 0; m <= N; m++) {
    for (let j = 0; j <= N; j++) {
      cosT[m * (N + 1) + j] = Math.cos((j * m * Math.PI) / N);
      sinT[m * (N + 1) + j] = Math.sin((j * m * Math.PI) / N);
    }
  }
  const p0 = {
    N, D: d / ell, H: 0, T: hasL ? 0 : (T * vel) / ell, kFixed: hasL ? 1 : null, kStart: 1, stokes, cur: cur / vel,
    cosT, sinT, tc: new Float64Array(N + 1), ts: new Float64Array(N + 1),
  };
  const Hn = H / ell;
  // relative height against a Miche-type limit decides the initial number of height steps
  const Hb = Math.min(0.142 * TWO_PI * Math.tanh(p0.D), 0.8 * p0.D);
  let steps = heightSteps ?? Math.max(1, Math.ceil(1 + 10 * (Hn / Hb) ** 2));
  if (!Number.isInteger(steps) || steps < 1) throw new RangeError('heightSteps must be a positive integer');
  let sol = null;
  for (let attempt = 0; attempt < 5 && sol === null; attempt++) {
    sol = tryHeightSteps(p0, Hn, steps, tol);
    if (sol === null) steps *= 2;
  }
  if (sol === null) {
    throw new Error(`stream-function wave has no admissible converged solution (H = ${H} m, d = ${d} m, order ${N}): `
      + 'the wave is at or beyond the breaking limit, or the order is too low');
  }
  const z = sol.z;
  const iB = N + 1; const ik = 2 * N + 2;
  const k = z[ik] / ell;
  const c = z[ik + 1] * vel;
  const B = new Array(N + 1).fill(0);
  for (let j = 1; j <= N; j++) B[j] = z[iB + j] * ell * vel;
  const etaNodes = new Array(N + 1);
  for (let m = 0; m <= N; m++) etaNodes[m] = z[m] * ell - d;
  // discrete cosine transform of the surface ordinates
  const etaCoef = new Array(N + 1).fill(0);
  for (let j = 0; j <= N; j++) {
    let s = 0.5 * (etaNodes[0] + etaNodes[N] * Math.cos(j * Math.PI));
    for (let m = 1; m < N; m++) s += etaNodes[m] * cosT[m * (N + 1) + j];
    etaCoef[j] = ((j === 0 || j === N ? 1 : 2) * s) / N;
  }
  const out = {
    H, d, g, order: N, k, L: TWO_PI / k, T: TWO_PI / (k * c), c, omega: k * c,
    currentEulerian: c + z[iB] * vel, currentStokes: c - (z[ik + 2] * ell * vel) / d,
    B0: z[iB] * vel, B, Q: z[ik + 2] * ell * vel, R: z[ik + 3] * g * ell,
    etaNodes, etaCoef, crest: etaNodes[0], trough: etaNodes[N],
    report: { converged: true, iterations: sol.iterations, heightSteps: steps, residual: sol.res,
      kinematicResidual: 0, dynamicResidual: 0 },
  };
  // truncation-error measure between collocation points
  let rk = 0; let rd = 0;
  for (let m = 0; m < N; m++) {
    const theta = ((m + 0.5) * Math.PI) / N;
    let eta = 0;
    for (let j = 0; j <= N; j++) eta += etaCoef[j] * Math.cos(j * theta);
    const f = field(out, theta, eta + d);
    rk = Math.max(rk, Math.abs(f.psi + out.Q) / out.Q);
    rd = Math.max(rd, Math.abs(0.5 * (f.u * f.u + f.w * f.w) + g * (eta + d) - out.R) / (g * H));
  }
  out.report.kinematicResidual = rk;
  out.report.dynamicResidual = rd;
  return out;
}

// Wave-frame stream function, velocity and velocity gradients at phase theta and height y above the bed.
function field(sol, theta, y) {
  const { k, d, B, order } = sol;
  let psi = sol.B0 * y; let u = sol.B0; let w = 0; let ux = 0; let uy = 0;
  for (let j = 1; j <= order; j++) {
    const jk = j * k;
    const e = Math.exp(jk * (y - d));
    const r = Math.exp(-2 * jk * y);
    const q = Math.exp(-2 * jk * d);
    const S = (e * (1 - r)) / (1 + q);
    const C = (e * (1 + r)) / (1 + q);
    const cj = Math.cos(j * theta); const sj = Math.sin(j * theta);
    const b = B[j];
    psi += b * S * cj;
    u += jk * b * C * cj;
    w += jk * b * S * sj;
    ux -= jk * jk * b * C * sj;
    uy += jk * jk * b * S * cj;
  }
  return { psi, u, w, ux, uy };
}

/**
 * Free-surface elevation of a solved stream-function wave at phase position x − c t.
 * The cosine series through the collocation ordinates is refined by Newton iteration onto the streamline ψ = −Q.
 * Reference: Fenton (1988) Eq. 3.5 and §3.
 * @param {object} sol result of `streamFunctionWave`
 * @param {number} [x=0] horizontal position, m (crest at x = c t)
 * @param {number} [t=0] time, s
 * @returns {number} surface elevation above still water, m
 */
export function streamFunctionEta(sol, x = 0, t = 0) {
  const theta = sol.k * (x - sol.c * t);
  let eta = 0;
  for (let j = 0; j <= sol.order; j++) eta += sol.etaCoef[j] * Math.cos(j * theta);
  let y = eta + sol.d;
  for (let i = 0; i < 6; i++) {
    const f = field(sol, theta, y);
    const dy = (f.psi + sol.Q) / f.u;
    if (!Number.isFinite(dy) || Math.abs(dy) > 0.25 * sol.H) break;
    y -= dy;
    if (Math.abs(dy) < 1e-13 * sol.d) break;
  }
  return y - sol.d;
}

/**
 * Surface profile of a solved stream-function wave sampled over one wavelength.
 * @param {object} sol result of `streamFunctionWave`
 * @param {number} [n=73] number of points from x = −L/2 to +L/2 (crest at x = 0)
 * @returns {{x:Float64Array,eta:Float64Array}} x m; eta m above still water
 */
export function streamFunctionProfile(sol, n = 73) {
  if (!Number.isInteger(n) || n < 2) throw new RangeError('n must be an integer >= 2');
  const x = new Float64Array(n); const eta = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    x[i] = sol.L * (i / (n - 1) - 0.5);
    eta[i] = streamFunctionEta(sol, x[i], 0);
  }
  return { x, eta };
}

/**
 * Kinematics and pressure of a solved stream-function wave at a point in the fixed frame.
 *   u = c + ∂ψ/∂y,  w = −∂ψ/∂x  (wave frame translated at c);  local accelerations ∂/∂t = −c ∂/∂x;
 *   particle accelerations Du/Dt = (u − c) u_x + w u_z;  p = ρ [R − g y − ½((u − c)² + w²)].
 * Reference: Rienecker & Fenton (1981) Eqs. 6–9; Fenton (1988) §2.
 * Validity: −d ≤ z ≤ instantaneous surface; above the surface zeros are returned with wet = false.
 * @param {object} sol result of `streamFunctionWave`
 * @param {object} p
 * @param {number} [p.x=0] m (crest at x = c t)
 * @param {number} [p.z=0] m above still-water level
 * @param {number} [p.t=0] s
 * @param {number} [p.rho=1025] kg/m³
 * @returns {{eta:number,u:number,w:number,ax:number,az:number,axTotal:number,azTotal:number,p:number,pDyn:number,
 *   wet:boolean}} eta m; u, w m/s (u includes the mean current); ax, az local accelerations m/s²; axTotal, azTotal
 *   particle accelerations m/s²; p gauge pressure Pa; pDyn = p + ρ g z Pa.
 */
export function streamFunctionKinematics(sol, { x = 0, z = 0, t = 0, rho = 1025 } = {}) {
  const { d, c, g } = sol;
  if (z < -d * (1 + 1e-12)) throw new RangeError('z is below the seabed');
  const eta = streamFunctionEta(sol, x, t);
  if (z > eta + 1e-9 * sol.H) return { eta, u: 0, w: 0, ax: 0, az: 0, axTotal: 0, azTotal: 0, p: 0, pDyn: 0, wet: false };
  const theta = sol.k * (x - c * t);
  const y = z + d;
  const f = field(sol, theta, y);
  const pr = rho * (sol.R - g * y - 0.5 * (f.u * f.u + f.w * f.w));
  return {
    eta, u: f.u + c, w: f.w,
    ax: -c * f.ux, az: -c * f.uy,
    axTotal: f.u * f.ux + f.w * f.uy, azTotal: f.u * f.uy - f.w * f.ux,
    p: pr, pDyn: pr + rho * g * z, wet: true,
  };
}
