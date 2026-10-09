// Dynamic positioning: three-degree-of-freedom (surge, sway, yaw) thrust allocation, capability polars,
// closed-loop station-keeping simulation and LQR design.
// Vessel frame: x forward, y to port, yaw positive counter-clockwise seen from above. A thruster at (x, y)
// delivering thrust T in direction alpha (deg from +x towards +y, direction of the force on the vessel) gives
//   Fx = T cos(alpha), Fy = T sin(alpha), Mz = x Fy - y Fx.

import { inverse, det, pinv, leastSquares, matMul, transpose } from '../math/linalg.js';
import { makeRng } from '../math/rng.js';
import { needPos } from './_wavehelp.js';

const DEG = Math.PI / 180;
const wrap360 = (a) => ((a % 360) + 360) % 360;
const wrapPi = (a) => Math.atan2(Math.sin(a), Math.cos(a));

function checkThrusters(thrusters) {
  if (!Array.isArray(thrusters) || thrusters.length === 0) throw new RangeError('thrusters must be a non-empty array');
  return thrusters.map((t, i) => {
    if (!Number.isFinite(t.x) || !Number.isFinite(t.y)) throw new RangeError(`thruster ${i}: x and y must be finite`);
    if (!(t.Tmax > 0)) throw new RangeError(`thruster ${i}: Tmax must be > 0`);
    const type = t.type || 'azimuth';
    if (type !== 'azimuth' && type !== 'tunnel' && type !== 'fixed') throw new RangeError(`thruster ${i}: type must be 'azimuth', 'tunnel' or 'fixed'`);
    if (type === 'fixed' && !Number.isFinite(t.angleDeg)) throw new RangeError(`thruster ${i}: a fixed thruster needs angleDeg`);
    const angle = type === 'tunnel' ? 90 : t.angleDeg;
    const Tmin = t.Tmin !== undefined ? t.Tmin : type === 'azimuth' ? 0 : -t.Tmax;
    return { name: t.name || `T${i + 1}`, x: t.x, y: t.y, type, Tmax: t.Tmax, Tmin, angle, forbidden: t.forbidden || [] };
  });
}

function inSector(a, sector) {
  const span = wrap360(sector[1] - sector[0]);
  const rel = wrap360(a - sector[0]);
  return rel > 1e-9 && rel < span - 1e-9;
}

/**
 * Static thrust allocation for azimuth, tunnel and fixed-direction thrusters.
 *  - 'pinv': weighted minimum-norm solution u = W⁻¹ Bᵀ (B W⁻¹ Bᵀ)⁺ tau with W = diag(1/Tmax²), i.e. the
 *    allocation minimising the sum of squared utilisations; limits are reported but not enforced.
 *  - 'bounded' (default): redistributed pseudo-inverse / active-set iteration. The minimum-norm solution is
 *    computed for the free thrusters; the worst limit violation is then fixed at its bound (an azimuth thruster
 *    at Tmax in its current direction, a tunnel/fixed thruster at its limit, an azimuth thruster pointing into
 *    a forbidden sector constrained to the nearest sector edge) and the remaining demand is re-allocated,
 *    until no limits are violated or no free thrusters remain. If the demand cannot be met, the largest
 *    achievable fraction of it (same direction, by bisection) is allocated and `feasible` is false.
 * Reference: Fossen (2011) Handbook of Marine Craft Hydrodynamics and Motion Control §12.3 (control
 * allocation, generalised inverse and constrained allocation); Sørdalen (1997) Optimal thrust allocation for
 * marine vessels, Control Eng. Practice 5(9). The bounded method is feasible-seeking, not a certified global
 * optimum; thrust losses and thruster interaction are not modelled.
 * @param {object} o
 * @param {{name?:string, x:number, y:number, type?:string, Tmax:number, Tmin?:number, angleDeg?:number,
 *   forbidden?:number[][]}[]} o.thrusters position m; type 'azimuth' (default) | 'tunnel' | 'fixed'; Tmax N;
 *   Tmin N (reverse thrust, default -Tmax for tunnel/fixed, 0 for azimuth); angleDeg for 'fixed'; forbidden
 *   azimuth sectors [fromDeg, toDeg] counter-clockwise
 * @param {{Fx:number, Fy:number, Mz:number}} o.demand N, N, N·m (vessel frame)
 * @param {string} [o.method='bounded']
 * @param {number[]} [o.disabled] indices of failed thrusters
 * @returns {{thrusters:{name:string, T:number, angleDeg:number, Tx:number, Ty:number, utilisation:number}[],
 *   achieved:{Fx:number, Fy:number, Mz:number}, residual:{Fx:number, Fy:number, Mz:number}, feasible:boolean,
 *   maxUtilisation:number, demandFraction:number, iterations:number}} thrust N (signed for tunnel/fixed), direction
 *   deg, components N; demandFraction = delivered fraction of the demand (1 when feasible)
 */
export function allocateThrust({ thrusters, demand, method = 'bounded', disabled = [] }) {
  const th = checkThrusters(thrusters);
  if (!demand || !Number.isFinite(demand.Fx + demand.Fy + demand.Mz)) throw new RangeError('demand must have finite Fx, Fy, Mz');
  if (method !== 'bounded' && method !== 'pinv') throw new RangeError("method must be 'bounded' or 'pinv'");
  const full = allocateCore(th, demand, method, disabled);
  full.demandFraction = 1;
  if (method === 'pinv' || full.feasible) return full;
  // demand beyond the capability: deliver the largest achievable fraction of it, in the demanded direction
  let lo = 0; let hi = 1; let best = null;
  for (let i = 0; i < 30 && hi - lo > 1e-4; i++) {
    const f = 0.5 * (lo + hi);
    const r = allocateCore(th, { Fx: f * demand.Fx, Fy: f * demand.Fy, Mz: f * demand.Mz }, method, disabled);
    if (r.feasible) { lo = f; best = r; } else hi = f;
  }
  if (!best) best = allocateCore(th, { Fx: 0, Fy: 0, Mz: 0 }, method, disabled);
  best.residual = { Fx: demand.Fx - best.achieved.Fx, Fy: demand.Fy - best.achieved.Fy, Mz: demand.Mz - best.achieved.Mz };
  best.feasible = false; best.demandFraction = lo; best.iterations += full.iterations;
  return best;
}

function allocateCore(th, demand, method, disabled) {
  const n = th.length;
  let Lc = 0; let Tsum = 0;
  for (const t of th) { Lc = Math.max(Lc, Math.hypot(t.x, t.y)); Tsum += t.Tmax; }
  if (!(Lc > 0)) Lc = 1;
  // state: mode 2 = free azimuth (Tx, Ty), 1 = free scalar along `ang`, 0 = fixed value
  const st = th.map((t, i) => {
    if (disabled.includes(i)) return { mode: 0, Tx: 0, Ty: 0 };
    if (t.type === 'azimuth') return { mode: 2, Tx: 0, Ty: 0, triedOther: false };
    return { mode: 1, ang: t.angle, lo: t.Tmin, hi: t.Tmax, Tx: 0, Ty: 0 };
  });
  const tau = [demand.Fx, demand.Fy, demand.Mz / Lc];
  let iterations = 0;
  const maxIt = method === 'pinv' ? 1 : 8 * n + 10;
  for (; iterations < maxIt; iterations++) {
    const rem = tau.slice();
    const cols = []; const wInv = []; const owner = [];
    for (let i = 0; i < n; i++) {
      const t = th[i]; const s = st[i];
      if (s.mode === 0) {
        rem[0] -= s.Tx; rem[1] -= s.Ty; rem[2] -= (t.x * s.Ty - t.y * s.Tx) / Lc;
      } else if (s.mode === 2) {
        cols.push([1, 0, -t.y / Lc], [0, 1, t.x / Lc]);
        wInv.push(t.Tmax * t.Tmax, t.Tmax * t.Tmax); owner.push(i, i);
      } else {
        const c = Math.cos(s.ang * DEG); const sn = Math.sin(s.ang * DEG);
        cols.push([c, sn, (t.x * sn - t.y * c) / Lc]);
        wInv.push(t.Tmax * t.Tmax); owner.push(i);
      }
    }
    if (cols.length === 0) break;
    const Gm = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let k = 0; k < cols.length; k++) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) Gm[a][b] += cols[k][a] * wInv[k] * cols[k][b];
    const Gp = pinv(Gm);
    const lam = [0, 1, 2].map((a) => Gp[a][0] * rem[0] + Gp[a][1] * rem[1] + Gp[a][2] * rem[2]);
    let k = 0;
    for (let i = 0; i < n; i++) {
      const s = st[i];
      if (s.mode === 2) {
        s.Tx = wInv[k] * (cols[k][0] * lam[0] + cols[k][1] * lam[1] + cols[k][2] * lam[2]);
        s.Ty = wInv[k + 1] * (cols[k + 1][0] * lam[0] + cols[k + 1][1] * lam[1] + cols[k + 1][2] * lam[2]);
        k += 2;
      } else if (s.mode === 1) {
        s.u = wInv[k] * (cols[k][0] * lam[0] + cols[k][1] * lam[1] + cols[k][2] * lam[2]);
        s.Tx = s.u * Math.cos(s.ang * DEG); s.Ty = s.u * Math.sin(s.ang * DEG);
        k += 1;
      }
    }
    if (method === 'pinv') { iterations = 1; break; }
    // forbidden sectors first, then the worst limit violation
    let acted = false;
    for (let i = 0; i < n && !acted; i++) {
      const s = st[i]; const t = th[i];
      if (s.mode !== 2) continue;
      const mag = Math.hypot(s.Tx, s.Ty);
      if (mag < 1e-9 * t.Tmax) continue;
      const a = wrap360(Math.atan2(s.Ty, s.Tx) / DEG);
      for (const sec of t.forbidden) {
        if (!inSector(a, sec)) continue;
        const d0 = wrap360(a - sec[0]); const d1 = wrap360(sec[1] - a);
        s.mode = 1; s.ang = d0 <= d1 ? sec[0] : sec[1]; s.other = d0 <= d1 ? sec[1] : sec[0];
        s.lo = 0; s.hi = t.Tmax; acted = true;
        break;
      }
    }
    if (acted) continue;
    let worst = -1; let ratio = 1 + 1e-9;
    for (let i = 0; i < n; i++) {
      const s = st[i]; const t = th[i];
      let r = 0;
      if (s.mode === 2) r = Math.hypot(s.Tx, s.Ty) / t.Tmax;
      else if (s.mode === 1) {
        if (s.u > s.hi) r = s.hi > 0 ? s.u / s.hi : Infinity;
        else if (s.u < s.lo) r = s.lo < 0 ? s.u / s.lo : Infinity;
      }
      if (r > ratio) { ratio = r; worst = i; }
    }
    if (worst < 0) break;
    const s = st[worst]; const t = th[worst];
    if (s.mode === 2) {
      const mag = Math.hypot(s.Tx, s.Ty);
      s.Tx *= t.Tmax / mag; s.Ty *= t.Tmax / mag; s.mode = 0;
    } else if (s.u < s.lo && s.other !== undefined && !s.triedOther) {
      // pushing against the wrong edge of a forbidden sector: try the other edge once
      s.triedOther = true; const o2 = s.ang; s.ang = s.other; s.other = o2;
    } else {
      const u = Math.min(s.hi, Math.max(s.lo, s.u));
      s.Tx = u * Math.cos(s.ang * DEG); s.Ty = u * Math.sin(s.ang * DEG); s.u = u; s.mode = 0; s.fixedAng = s.ang;
    }
  }
  let Fx = 0; let Fy = 0; let Mz = 0; let maxU = 0; let limitsOk = true;
  const out = th.map((t, i) => {
    const s = st[i];
    Fx += s.Tx; Fy += s.Ty; Mz += t.x * s.Ty - t.y * s.Tx;
    let T; let angleDeg; let util;
    if (t.type === 'azimuth') {
      T = Math.hypot(s.Tx, s.Ty);
      angleDeg = T > 0 ? wrap360(Math.atan2(s.Ty, s.Tx) / DEG) : 0;
      util = T / t.Tmax;
      if (T > 1e-9 * t.Tmax && t.forbidden.some((sec) => inSector(angleDeg, sec))) limitsOk = false;
    } else {
      T = s.Tx * Math.cos(t.angle * DEG) + s.Ty * Math.sin(t.angle * DEG);
      angleDeg = t.angle;
      util = T >= 0 ? T / t.Tmax : t.Tmin < 0 ? T / t.Tmin : Infinity;
    }
    if (util > 1 + 1e-6) limitsOk = false;
    maxU = Math.max(maxU, util);
    return { name: t.name, T, angleDeg, Tx: s.Tx, Ty: s.Ty, utilisation: util };
  });
  const residual = { Fx: demand.Fx - Fx, Fy: demand.Fy - Fy, Mz: demand.Mz - Mz };
  const resN = Math.hypot(residual.Fx, residual.Fy, residual.Mz / Lc);
  return {
    thrusters: out, achieved: { Fx, Fy, Mz }, residual, feasible: limitsOk && resN <= 1e-6 * Tsum, maxUtilisation: maxU, iterations,
  };
}

/**
 * DP capability polar: for each environment direction, the largest wind speed at which the thrusters can
 * balance the mean environmental force (bisection on wind speed with the bounded allocation), for the intact
 * system and for the worst single failure (or the supplied failure groups).
 * The caller supplies the collinear wind–wave–current relation through `envForce(windSpeed, directionDeg)`,
 * which returns the total mean force in the vessel frame for that wind speed and direction.
 * Reference: IMCA M 140 Rev. 1 (2000) Specification for DP capability plots; DNV-ST-0111 (2021) Assessment of
 * station keeping capability of DP vessels. Static balance only: no dynamic
 * allowance unless `maxUtilisation` < 1 is specified (e.g. 0.8 to reserve 20 % thrust for dynamics).
 * The environmental force is assumed to grow monotonically with wind speed.
 * @param {object} o
 * @param {object[]} o.thrusters as allocateThrust
 * @param {function(number, number): {Fx:number, Fy:number, Mz:number}} o.envForce (windSpeed m/s, directionDeg) -> N, N, N·m
 * @param {number[]} [o.directionsDeg] default 0…350 step 10
 * @param {number} [o.vMax=50] upper bound of the search, m/s
 * @param {'single'|'none'|number[][]} [o.failures='single'] failure cases: each thruster in turn, none, or groups of indices
 * @param {number} [o.maxUtilisation=1] allowed fraction of Tmax
 * @returns {{directionDeg:number[], intact:number[], failure:number[]|null, failureCase:(number[]|null)[]|null}}
 *   limiting wind speeds m/s (equal to vMax where the limit is not reached) and the governing failed thruster
 *   indices per direction
 */
export function dpCapability({ thrusters, envForce, directionsDeg, vMax = 50, failures = 'single', maxUtilisation = 1 }) {
  if (typeof envForce !== 'function') throw new RangeError('envForce must be a function');
  needPos(vMax, 'vMax');
  if (!(maxUtilisation > 0 && maxUtilisation <= 1)) throw new RangeError('maxUtilisation must be within (0, 1]');
  const th = thrusters.map((t) => ({ ...t, Tmax: t.Tmax * maxUtilisation, Tmin: t.Tmin === undefined ? undefined : t.Tmin * maxUtilisation }));
  const dirs = directionsDeg || Array.from({ length: 36 }, (_, i) => 10 * i);
  const cases = failures === 'single' ? th.map((_, i) => [i]) : failures === 'none' ? [] : failures;
  const limit = (dir, disabled) => {
    const ok = (v) => allocateThrust({ thrusters: th, demand: envForce(v, dir), disabled }).feasible;
    if (ok(vMax)) return vMax;
    if (!ok(0)) return 0;
    let lo = 0; let hi = vMax;
    for (let i = 0; i < 40 && hi - lo > 1e-3; i++) {
      const mid = 0.5 * (lo + hi);
      if (ok(mid)) lo = mid; else hi = mid;
    }
    return lo;
  };
  const out = { directionDeg: dirs.slice(), intact: [], failure: cases.length ? [] : null, failureCase: cases.length ? [] : null };
  for (const d of dirs) {
    out.intact.push(limit(d, []));
    if (cases.length) {
      let worst = Infinity; let wc = null;
      for (const c of cases) { const v = limit(d, c); if (v < worst) { worst = v; wc = c; } }
      out.failure.push(worst); out.failureCase.push(wc);
    }
  }
  return out;
}

/**
 * Solves the continuous-time algebraic Riccati equation Aᵀ P + P A - P B R⁻¹ Bᵀ P + Q = 0 and returns the LQR
 * gain K = R⁻¹ Bᵀ P (control u = -K x). Method: matrix sign function of the Hamiltonian
 * H = [A, -B R⁻¹ Bᵀ; -Q, -Aᵀ] by the Newton iteration Z ← ½ (c Z + (c Z)⁻¹) with determinant scaling
 * c = |det Z|^(-1/2n); with sign(H) = [W11 W12; W21 W22], P is the least-squares solution of
 * [W12; W22 + I] P = -[W11 + I; W21].
 * Reference: Roberts (1980) Linear model reduction and solution of the algebraic Riccati equation by use of
 * the sign function, Int. J. Control 32(4); Anderson & Moore (1990) Optimal Control: Linear Quadratic Methods.
 * Requires (A, B) stabilisable, (A, Q) detectable, R positive definite. Intended for small systems (n <= ~20).
 * @param {{A:number[][], B:number[][], Q:number[][], R:number[][]}} o state n×n, input n×m, weights n×n and m×m
 * @returns {{K:number[][], P:number[][], residual:number, iterations:number}} gain m×n, Riccati solution n×n,
 *   max-norm of the Riccati residual relative to max|Q|
 */
export function lqr({ A, B, Q, R }) {
  const n = A.length;
  if (!n || A[0].length !== n || B.length !== n || Q.length !== n) throw new RangeError('A, B, Q dimensions are inconsistent');
  const m = B[0].length;
  if (R.length !== m) throw new RangeError('R must be m×m');
  const Rinv = inverse(R);
  const Bt = transpose(B);
  const Gm = matMul(matMul(B, Rinv), Bt);
  let Z = [];
  for (let i = 0; i < 2 * n; i++) Z.push(new Array(2 * n).fill(0));
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    Z[i][j] = A[i][j]; Z[i][n + j] = -Gm[i][j]; Z[n + i][j] = -Q[i][j]; Z[n + i][n + j] = -A[j][i];
  }
  let iterations = 0;
  for (; iterations < 200; iterations++) {
    const d = Math.abs(det(Z));
    if (!(d > 0) || !Number.isFinite(d)) throw new Error('LQR: Hamiltonian is singular (system not stabilisable/detectable?)');
    const c = d ** (-1 / (2 * n));
    const Zi = inverse(Z);
    let diff = 0; let nrm = 0;
    const Zn = Z.map((row, i) => row.map((v, j) => {
      const nv = 0.5 * (c * v + Zi[i][j] / c);
      diff = Math.max(diff, Math.abs(nv - v)); nrm = Math.max(nrm, Math.abs(nv));
      return nv;
    }));
    Z = Zn;
    if (diff <= 1e-13 * nrm) { iterations++; break; }
  }
  const lhs = []; const rhs = [];
  for (let i = 0; i < n; i++) { lhs.push(Z[i].slice(n)); rhs.push(Z[i].slice(0, n).map((v, j) => -(v + (i === j ? 1 : 0)))); }
  for (let i = 0; i < n; i++) { lhs.push(Z[n + i].slice(n).map((v, j) => v + (i === j ? 1 : 0))); rhs.push(Z[n + i].slice(0, n).map((v) => -v)); }
  const P = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let j = 0; j < n; j++) {
    const sol = leastSquares(lhs, rhs.map((r) => r[j])).x;
    for (let i = 0; i < n; i++) P[i][j] = sol[i];
  }
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const s = 0.5 * (P[i][j] + P[j][i]); P[i][j] = s; P[j][i] = s; }
  const K = matMul(matMul(Rinv, Bt), P);
  const AtP = matMul(transpose(A), P); const PA = matMul(P, A); const PGP = matMul(matMul(P, Gm), P);
  let residual = 0; let qmax = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    residual = Math.max(residual, Math.abs(AtP[i][j] + PA[i][j] - PGP[i][j] + Q[i][j]));
    qmax = Math.max(qmax, Math.abs(Q[i][j]));
  }
  return { K, P, residual: qmax > 0 ? residual / qmax : residual, iterations };
}

/**
 * LQR position/velocity feedback gains for the low-speed DP model eta' = nu, M nu' + D nu = tau
 * (state [x, y, psi, u, v, r] in the vessel-parallel frame): tau = -Kp e - Kd nu.
 * Reference: Fossen (2011) ch. 13 (linear quadratic optimal control); lqr() above.
 * @param {{M:number[][], D:number[][], Q:number[][], R:number[][]}} o 3×3 inertia (incl. added mass) and linear
 *   damping; Q 6×6 state weights; R 3×3 control weights
 * @returns {{Kp:number[][], Kd:number[][], K:number[][], residual:number}} N/m … (3×3 each)
 */
export function dpLqrGains({ M, D, Q, R }) {
  const Mi = inverse(M);
  const MiD = matMul(Mi, D);
  const A = []; const B = [];
  for (let i = 0; i < 6; i++) { A.push(new Array(6).fill(0)); B.push([0, 0, 0]); }
  for (let i = 0; i < 3; i++) {
    A[i][3 + i] = 1;
    for (let j = 0; j < 3; j++) { A[3 + i][3 + j] = -MiD[i][j]; B[3 + i][j] = Mi[i][j]; }
  }
  const r = lqr({ A, B, Q, R });
  return { Kp: r.K.map((row) => row.slice(0, 3)), Kd: r.K.map((row) => row.slice(3)), K: r.K, residual: r.residual };
}

/**
 * PID gains by pole placement on the decoupled diagonal model m x'' + d x' = tau:
 *   Kp = m wn², Kd = 2 zeta wn m - d, Ki = (wn/10) Kp.
 * Reference: Fossen (2011) §12.2 (PID control and acceleration feedback, pole-placement algorithm).
 * @param {{M:number[][], D:number[][], omegaN:number[], zeta?:number[]}} o natural frequencies rad/s per DOF
 *   (typical DP bandwidth: periods of 60–150 s), damping ratios (default 0.7)
 * @returns {{Kp:number[], Ki:number[], Kd:number[]}}
 */
export function dpPidGains({ M, D, omegaN, zeta = [0.7, 0.7, 0.7] }) {
  const Kp = []; const Ki = []; const Kd = [];
  for (let i = 0; i < 3; i++) {
    needPos(omegaN[i], 'omegaN');
    Kp.push(M[i][i] * omegaN[i] ** 2);
    Kd.push(2 * zeta[i] * omegaN[i] * M[i][i] - D[i][i]);
    Ki.push((omegaN[i] / 10) * Kp[i]);
  }
  return { Kp, Ki, Kd };
}

function steadyKalmanGain(Phi, H, Qd, Rd) {
  const n = Phi.length; const m = H.length;
  let P = Qd.map((r, i) => r.map((v, j) => (i === j ? 100 * (v || 1e-6) : 0)));
  const Ht = transpose(H); const PhiT = transpose(Phi);
  let K = null;
  for (let it = 0; it < 5000; it++) {
    const Pp = matMul(matMul(Phi, P), PhiT);
    for (let i = 0; i < n; i++) Pp[i][i] += Qd[i][i];
    const S = matMul(matMul(H, Pp), Ht);
    for (let i = 0; i < m; i++) S[i][i] += Rd[i];
    const Kn = matMul(matMul(Pp, Ht), inverse(S));
    const KH = matMul(Kn, H);
    P = matMul(KH.map((r, i) => r.map((v, j) => (i === j ? 1 : 0) - v)), Pp);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const s = 0.5 * (P[i][j] + P[j][i]); P[i][j] = s; P[j][i] = s; }
    let diff = 0; let nrm = 0;
    if (K) for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) { diff = Math.max(diff, Math.abs(Kn[i][j] - K[i][j])); nrm = Math.max(nrm, Math.abs(Kn[i][j])); }
    K = Kn;
    if (it > 10 && diff <= 1e-10 * nrm) break;
  }
  return K;
}

/**
 * Closed-loop DP station-keeping simulation in three degrees of freedom.
 * Plant (low speed, Coriolis and nonlinear damping neglected): eta' = R(psi) nu, M nu' + D nu = tau_thr +
 * R(psi)ᵀ tau_env(t). Measurements = position/heading + first-order wave motion + sensor noise. A steady-state
 * Kalman filter (states: position, velocity, slowly varying bias force, and a second-order wave-frequency
 * model per axis) separates low-frequency motion from wave-frequency motion (wave filtering); a PID controller
 * acts on the filtered low-frequency estimates; the commanded force is distributed by `allocateThrust`
 * (bounded), and each thruster follows its command with thrust and azimuth rate limits (thrust is reduced by
 * the cosine of the azimuth error while slewing).
 * Reference: Fossen (2011) ch. 11 (wave filtering, Kalman filter design for DP), ch. 12; Balchen,
 * Jenssen & Sælid (1976) Dynamic positioning using Kalman filtering and optimal control theory.
 * Screening-level simulation: thrust losses, thruster–hull interaction, power limits and second-order wave
 * forces must be included in `envForce` by the caller if required.
 * @param {object} o
 * @param {number[][]} o.M 3×3 inertia incl. added mass (kg, kg·m, kg·m²)
 * @param {number[][]} o.D 3×3 linear damping (N·s/m …)
 * @param {object[]} o.thrusters as allocateThrust
 * @param {{x:number, y:number, psiDeg:number}} [o.setpoint] earth-fixed position m and heading deg
 * @param {{x?:number, y?:number, psiDeg?:number}} [o.initial] initial position (default: at the setpoint)
 * @param {{Kp:number[], Ki:number[], Kd:number[]}} o.gains PID gains per DOF (e.g. from dpPidGains)
 * @param {function(number): {Fx:number, Fy:number, Mz:number}} [o.envForce] (t) -> mean + slowly varying force in the
 *   EARTH-fixed frame, N, N, N·m
 * @param {{sigma:number[], period:number}} [o.waveMotion] first-order motion in the measurements: standard
 *   deviations [m, m, rad] and central period s
 * @param {number[]} [o.noiseStd=[0,0,0]] sensor noise standard deviations [m, m, rad]
 * @param {{wavePeriod?:number, waveDamping?:number, R?:number[], qBias?:number[], qVel?:number[], qWave?:number[]}} [o.observer]
 *   filter tuning: wave model period s (default waveMotion.period or 9) and relative damping (0.1); measurement
 *   variances; process-noise intensities. Defaults are starting points derived from the inputs.
 * @param {number} [o.thrustRate=0.1] thrust slew rate, fraction of Tmax per second
 * @param {number} [o.azimuthRateDeg=12] azimuth slew rate, deg/s
 * @param {number} [o.dt=0.2] s
 * @param {number} o.duration s
 * @param {number} [o.statsStart] start of the statistics window, s (default 25 % of the duration)
 * @param {number} [o.seed=1]
 * @returns {{t:Float64Array, x:Float64Array, y:Float64Array, psiDeg:Float64Array, footprint:{maxRadius:number,
 *   r95:number, stdX:number, stdY:number, stdPsiDeg:number, maxPsiDevDeg:number, meanX:number, meanY:number},
 *   utilisation:{name:string, mean:number, max:number}[], saturatedFraction:number}} decimated histories (<= 2000
 *   points) of the low-frequency position m and heading deg; footprint relative to the setpoint (m, deg);
 *   thrust utilisation per thruster; fraction of steps with an infeasible allocation
 */
export function simulateDp({ M, D, thrusters, setpoint = { x: 0, y: 0, psiDeg: 0 }, initial = {}, gains, envForce, waveMotion, noiseStd = [0, 0, 0], observer = {}, thrustRate = 0.1, azimuthRateDeg = 12, dt = 0.2, duration, statsStart, seed = 1 }) {
  needPos(dt, 'dt'); needPos(duration, 'duration'); needPos(thrustRate, 'thrustRate'); needPos(azimuthRateDeg, 'azimuthRateDeg');
  if (!M || M.length !== 3 || !D || D.length !== 3) throw new RangeError('M and D must be 3×3');
  if (!gains || !gains.Kp || !gains.Kd) throw new RangeError('gains must provide Kp and Kd (and optionally Ki)');
  const th = checkThrusters(thrusters);
  const nT = th.length;
  const Mi = inverse(M);
  const MiD = matMul(Mi, D);
  const Ki = gains.Ki || [0, 0, 0];
  const rng = makeRng(seed);
  let Fref = 0; let Lc = 0;
  for (const t of th) { Fref += t.Tmax; Lc = Math.max(Lc, Math.hypot(t.x, t.y)); }
  if (!(Lc > 0)) Lc = 1;
  // --- observer ---------------------------------------------------------------------------------------
  const Tw = observer.wavePeriod || (waveMotion && waveMotion.period) || 9;
  const lamW = observer.waveDamping === undefined ? 0.1 : observer.waveDamping;
  const w0 = (2 * Math.PI) / Tw;
  const nS = 15;
  const Ac = Array.from({ length: nS }, () => new Array(nS).fill(0));
  for (let i = 0; i < 3; i++) {
    Ac[i][3 + i] = 1;
    for (let j = 0; j < 3; j++) { Ac[3 + i][3 + j] = -MiD[i][j]; Ac[3 + i][6 + j] = Mi[i][j]; }
    Ac[9 + i][12 + i] = 1;
    Ac[12 + i][9 + i] = -w0 * w0; Ac[12 + i][12 + i] = -2 * lamW * w0;
  }
  const Adt = Ac.map((r) => r.map((v) => v * dt));
  const A2 = matMul(Adt, Adt); const A3 = matMul(A2, Adt);
  const Phi = Adt.map((r, i) => r.map((v, j) => (i === j ? 1 : 0) + v + A2[i][j] / 2 + A3[i][j] / 6));
  // input matrix: Gamma = (I dt + A dt²/2 + A² dt³/6) Bc, Bc = [0; M⁻¹; 0; 0; 0]
  const Gam = Array.from({ length: nS }, () => [0, 0, 0]);
  for (let i = 0; i < nS; i++) for (let j = 0; j < 3; j++) {
    let s = 0;
    for (let k = 0; k < 3; k++) s += ((i === 3 + k ? 1 : 0) + Adt[i][3 + k] / 2 + A2[i][3 + k] / 6) * Mi[k][j];
    Gam[i][j] = s * dt;
  }
  const H = Array.from({ length: 3 }, (_, i) => { const r = new Array(nS).fill(0); r[i] = 1; r[12 + i] = 1; return r; });
  const small = [0.05, 0.05, 0.1 * DEG];
  const sigW = [0, 1, 2].map((i) => Math.max(waveMotion ? waveMotion.sigma[i] : 0, small[i]));
  const fScale = [Fref, Fref, Fref * Lc];
  const Rd = observer.R || [0, 1, 2].map((i) => Math.max(noiseStd[i], small[i]) ** 2);
  const qBias = observer.qBias || fScale.map((f) => (0.01 * f) ** 2);
  const qVel = observer.qVel || [0, 1, 2].map((i) => ((0.01 * fScale[i]) / M[i][i]) ** 2);
  const qWave = observer.qWave || sigW.map((s) => s * s * 4 * lamW * w0);
  const Qd = Array.from({ length: nS }, () => new Array(nS).fill(0));
  for (let i = 0; i < 3; i++) { Qd[3 + i][3 + i] = qVel[i] * dt; Qd[6 + i][6 + i] = qBias[i] * dt; Qd[12 + i][12 + i] = qWave[i] * dt; }
  const Kk = steadyKalmanGain(Phi, H, Qd, Rd);
  // --- states -----------------------------------------------------------------------------------------
  const psiD = setpoint.psiDeg * DEG;
  const eta = [initial.x === undefined ? setpoint.x : initial.x, initial.y === undefined ? setpoint.y : initial.y,
    (initial.psiDeg === undefined ? setpoint.psiDeg : initial.psiDeg) * DEG];
  const nu = [0, 0, 0];
  const cd = Math.cos(psiD); const sd = Math.sin(psiD);
  const toSp = (e) => [cd * (e[0] - setpoint.x) + sd * (e[1] - setpoint.y), -sd * (e[0] - setpoint.x) + cd * (e[1] - setpoint.y), wrapPi(e[2] - psiD)];
  const xh = new Array(nS).fill(0);
  const e0 = toSp(eta);
  for (let i = 0; i < 3; i++) xh[i] = e0[i];
  const integ = [0, 0, 0];
  const thr = th.map((t) => ({ T: 0, ang: t.type === 'azimuth' ? 0 : t.angle }));
  const wfW = [0.8, 0.9, 1.0, 1.1, 1.2].map((f) => f * w0);
  const wfPh = [0, 1, 2].map(() => wfW.map(() => 2 * Math.PI * rng.next()));
  const wfAmp = Math.sqrt(2 / wfW.length);
  const nSteps = Math.round(duration / dt);
  const stride = Math.max(1, Math.ceil((nSteps + 1) / 2000));
  const nOut = Math.floor(nSteps / stride) + 1;
  const tO = new Float64Array(nOut); const xO = new Float64Array(nOut); const yO = new Float64Array(nOut); const pO = new Float64Array(nOut);
  const t0s = statsStart === undefined ? 0.25 * duration : statsStart;
  let nSt = 0; let sx = 0; let sy = 0; let sxx = 0; let syy = 0; let sp = 0; let spp = 0; let maxR = 0; let maxP = 0;
  const radii = [];
  const uMean = new Float64Array(nT); const uMax = new Float64Array(nT);
  let nSat = 0;
  const tauB = [0, 0, 0];
  const deriv = (e, v, tb, te, out) => {
    const c = Math.cos(e[2]); const s = Math.sin(e[2]);
    out[0] = c * v[0] - s * v[1]; out[1] = s * v[0] + c * v[1]; out[2] = v[2];
    const f0 = tb[0] + c * te[0] + s * te[1]; const f1 = tb[1] - s * te[0] + c * te[1]; const f2 = tb[2] + te[2];
    for (let i = 0; i < 3; i++) {
      out[3 + i] = Mi[i][0] * f0 + Mi[i][1] * f1 + Mi[i][2] * f2 - (MiD[i][0] * v[0] + MiD[i][1] * v[1] + MiD[i][2] * v[2]);
    }
  };
  const k1 = new Array(6); const k2 = new Array(6); const k3 = new Array(6); const k4 = new Array(6);
  for (let n = 0; n <= nSteps; n++) {
    const t = n * dt;
    // measurement in the setpoint-parallel frame
    const e = toSp(eta);
    const y = [0, 1, 2].map((i) => {
      let wf = 0;
      if (waveMotion) for (let k = 0; k < wfW.length; k++) wf += wfAmp * waveMotion.sigma[i] * Math.sin(wfW[k] * t + wfPh[i][k]);
      return e[i] + wf + noiseStd[i] * rng.normal();
    });
    // filter update (prediction was made at the end of the previous step)
    const inno = [0, 1, 2].map((i) => { let d = y[i] - xh[i] - xh[12 + i]; if (i === 2) d = wrapPi(d); return d; });
    for (let i = 0; i < nS; i++) xh[i] += Kk[i][0] * inno[0] + Kk[i][1] * inno[1] + Kk[i][2] * inno[2];
    // PID on the low-frequency estimates
    const cmd = [0, 1, 2].map((i) => -gains.Kp[i] * xh[i] - gains.Kd[i] * xh[3 + i] - Ki[i] * integ[i]);
    const ce = Math.cos(xh[2]); const se = Math.sin(xh[2]);
    const demand = { Fx: ce * cmd[0] + se * cmd[1], Fy: -se * cmd[0] + ce * cmd[1], Mz: cmd[2] };
    const alloc = allocateThrust({ thrusters, demand });
    if (alloc.feasible) for (let i = 0; i < 3; i++) integ[i] += xh[i] * dt; else nSat++;
    // thruster dynamics
    tauB[0] = 0; tauB[1] = 0; tauB[2] = 0;
    for (let i = 0; i < nT; i++) {
      const tc = alloc.thrusters[i]; const s = thr[i]; const tt = th[i];
      let target = tc.T;
      if (tt.type === 'azimuth') {
        if (tc.T > 1e-6 * tt.Tmax) {
          const da = wrapPi((tc.angleDeg - s.ang) * DEG) / DEG;
          const lim = azimuthRateDeg * dt;
          s.ang += Math.max(-lim, Math.min(lim, da));
          target = tc.T * Math.max(0, Math.cos((tc.angleDeg - s.ang) * DEG));
        }
      }
      const dT = thrustRate * tt.Tmax * dt;
      s.T += Math.max(-dT, Math.min(dT, target - s.T));
      const ca = Math.cos(s.ang * DEG); const sa = Math.sin(s.ang * DEG);
      tauB[0] += s.T * ca; tauB[1] += s.T * sa; tauB[2] += s.T * (tt.x * sa - tt.y * ca);
      const u = Math.abs(s.T) / (s.T >= 0 ? tt.Tmax : Math.abs(tt.Tmin) || tt.Tmax);
      if (t >= t0s) { uMean[i] += u; uMax[i] = Math.max(uMax[i], u); }
    }
    // outputs and statistics
    if (n % stride === 0) { const k = n / stride; tO[k] = t; xO[k] = eta[0]; yO[k] = eta[1]; pO[k] = eta[2] / DEG; }
    if (t >= t0s) {
      const dx = eta[0] - setpoint.x; const dy = eta[1] - setpoint.y; const dp = wrapPi(eta[2] - psiD) / DEG;
      const r = Math.hypot(dx, dy);
      nSt++; sx += dx; sy += dy; sxx += dx * dx; syy += dy * dy; sp += dp; spp += dp * dp;
      maxR = Math.max(maxR, r); maxP = Math.max(maxP, Math.abs(dp)); radii.push(r);
    }
    if (n === nSteps) break;
    // filter prediction with the applied thrust expressed in the setpoint-parallel frame
    const tauS = [ce * tauB[0] - se * tauB[1], se * tauB[0] + ce * tauB[1], tauB[2]];
    const xp = new Array(nS);
    for (let i = 0; i < nS; i++) {
      let s = Gam[i][0] * tauS[0] + Gam[i][1] * tauS[1] + Gam[i][2] * tauS[2];
      const row = Phi[i];
      for (let j = 0; j < nS; j++) s += row[j] * xh[j];
      xp[i] = s;
    }
    for (let i = 0; i < nS; i++) xh[i] = xp[i];
    // plant step (RK4, thrust and environment held over the step)
    const fe = envForce ? envForce(t) : null;
    const te = fe ? [fe.Fx, fe.Fy, fe.Mz] : [0, 0, 0];
    const ev = [0, 0, 0]; const vv = [0, 0, 0];
    deriv(eta, nu, tauB, te, k1);
    for (let i = 0; i < 3; i++) { ev[i] = eta[i] + 0.5 * dt * k1[i]; vv[i] = nu[i] + 0.5 * dt * k1[3 + i]; }
    deriv(ev, vv, tauB, te, k2);
    for (let i = 0; i < 3; i++) { ev[i] = eta[i] + 0.5 * dt * k2[i]; vv[i] = nu[i] + 0.5 * dt * k2[3 + i]; }
    deriv(ev, vv, tauB, te, k3);
    for (let i = 0; i < 3; i++) { ev[i] = eta[i] + dt * k3[i]; vv[i] = nu[i] + dt * k3[3 + i]; }
    deriv(ev, vv, tauB, te, k4);
    for (let i = 0; i < 3; i++) {
      eta[i] += (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
      nu[i] += (dt / 6) * (k1[3 + i] + 2 * k2[3 + i] + 2 * k3[3 + i] + k4[3 + i]);
    }
  }
  radii.sort((a, b) => a - b);
  const sd2 = (s, ss) => Math.sqrt(Math.max(ss / nSt - (s / nSt) ** 2, 0));
  return {
    t: tO, x: xO, y: yO, psiDeg: pO,
    footprint: {
      maxRadius: maxR, r95: radii.length ? radii[Math.min(radii.length - 1, Math.floor(0.95 * radii.length))] : 0,
      stdX: sd2(sx, sxx), stdY: sd2(sy, syy), stdPsiDeg: sd2(sp, spp), maxPsiDevDeg: maxP, meanX: sx / nSt, meanY: sy / nSt,
    },
    utilisation: th.map((t, i) => ({ name: t.name, mean: uMean[i] / Math.max(nSt, 1), max: uMax[i] })),
    saturatedFraction: nSat / (nSteps + 1),
  };
}
