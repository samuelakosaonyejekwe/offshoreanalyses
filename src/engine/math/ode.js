/**
 * Time integration of ordinary differential equations and of second-order structural dynamics systems.
 *
 * First-order systems are written y' = f(t, y); t in seconds (or any consistent independent variable), y in the
 * caller's units. Structural systems are M ü + C u̇ + K u = f(t) in consistent SI units (kg, N·s/m, N/m, N, m).
 *
 * References: Hairer, Nørsett & Wanner, "Solving Ordinary Differential Equations I", 2nd ed. (HNW): classical RK4
 * §II.1, Dormand–Prince 5(4) pair §II.5 (Table 5.2), step-size control §II.4, dense output §II.6;
 * Gustafsson, ACM TOMS 17 (1991) 533–554 (PI step-size control);
 * Hairer & Wanner, "Solving Ordinary Differential Equations II" (implicit Euler, BDF2, simplified Newton §IV.8);
 * Newmark, J. Eng. Mech. Div. ASCE 85 (1959) 67–94; Hilber, Hughes & Taylor, Earthquake Eng. Struct. Dyn. 5 (1977)
 * 283–292 (HHT-α); Hughes, "The Finite Element Method" (1987) §9.
 */

import { luFactor, luSolve } from './linalg.js';

const EPS = 2.220446049250313e-16;

function isVec(v) {
  return Array.isArray(v) || ArrayBuffer.isView(v);
}

function checkOdeArgs(name, f, t0, y0, t1) {
  if (typeof f !== 'function') throw new RangeError(`${name}: f must be a function (t, y) -> dy`);
  if (!Number.isFinite(t0) || !Number.isFinite(t1)) throw new RangeError(`${name}: t0 and t1 must be finite`);
  if (!isVec(y0) || y0.length < 1) throw new RangeError(`${name}: y0 must be a non-empty array`);
}

function checkSteps(name, nSteps) {
  if (!Number.isInteger(nSteps) || nSteps < 1) throw new RangeError(`${name}: nSteps must be an integer >= 1`);
}

/**
 * Classical fourth-order Runge–Kutta with fixed step h = (t1 − t0) / nSteps (HNW §II.1).
 * Validity: non-stiff problems; global error O(h⁴). The state passed to f is a scratch array: copy it if kept.
 * @param {(t: number, y: Float64Array) => ArrayLike<number>} f right-hand side dy/dt
 * @param {number} t0 start time [s]
 * @param {ArrayLike<number>} y0 initial state
 * @param {number} t1 end time [s] (may be < t0)
 * @param {number} nSteps number of steps (>= 1)
 * @param {(t: number, y: Float64Array, k: number) => (void | boolean)} [onStep] called after each step k = 1..nSteps
 *   with the new time and state; returning false stops the integration (output truncated at that step)
 * @returns {{ t: Float64Array, y: Float64Array[] }} t[k] times [s]; y[k] state at t[k] (k = 0..nSteps)
 */
export function rk4(f, t0, y0, t1, nSteps, onStep) {
  checkOdeArgs('rk4', f, t0, y0, t1);
  checkSteps('rk4', nSteps);
  const n = y0.length;
  const h = (t1 - t0) / nSteps;
  let t = new Float64Array(nSteps + 1);
  let ys = new Array(nSteps + 1);
  let y = Float64Array.from(y0);
  const tmp = new Float64Array(n);
  t[0] = t0;
  ys[0] = y;
  for (let k = 0; k < nSteps; k++) {
    const tk = t0 + k * h;
    const k1 = f(tk, y);
    // f may return the same buffer on every call: accumulate the weighted sum as we go
    const acc = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      acc[i] = k1[i];
      tmp[i] = y[i] + 0.5 * h * k1[i];
    }
    const k2 = f(tk + 0.5 * h, tmp);
    for (let i = 0; i < n; i++) {
      acc[i] += 2 * k2[i];
      tmp[i] = y[i] + 0.5 * h * k2[i];
    }
    const k3 = f(tk + 0.5 * h, tmp);
    for (let i = 0; i < n; i++) {
      acc[i] += 2 * k3[i];
      tmp[i] = y[i] + h * k3[i];
    }
    const k4 = f(tk + h, tmp);
    const yn = acc;
    for (let i = 0; i < n; i++) yn[i] = y[i] + (h / 6) * (acc[i] + k4[i]);
    y = yn;
    t[k + 1] = k === nSteps - 1 ? t1 : t0 + (k + 1) * h;
    ys[k + 1] = y;
    if (onStep && onStep(t[k + 1], y, k + 1) === false) {
      t = t.slice(0, k + 2);
      ys = ys.slice(0, k + 2);
      break;
    }
  }
  return { t, y: ys };
}

// Dormand–Prince 5(4) coefficients (HNW Table II.5.2)
const A21 = 1 / 5;
const A31 = 3 / 40;
const A32 = 9 / 40;
const A41 = 44 / 45;
const A42 = -56 / 15;
const A43 = 32 / 9;
const A51 = 19372 / 6561;
const A52 = -25360 / 2187;
const A53 = 64448 / 6561;
const A54 = -212 / 729;
const A61 = 9017 / 3168;
const A62 = -355 / 33;
const A63 = 46732 / 5247;
const A64 = 49 / 176;
const A65 = -5103 / 18656;
const B1 = 35 / 384;
const B3 = 500 / 1113;
const B4 = 125 / 192;
const B5 = -2187 / 6784;
const B6 = 11 / 84;
// error coefficients: fifth-order minus embedded fourth-order weights
const E1 = 71 / 57600;
const E3 = -71 / 16695;
const E4 = 71 / 1920;
const E5 = -17253 / 339200;
const E6 = 22 / 525;
const E7 = -1 / 40;
// continuous extension (Shampine; HNW §II.6, as used in DOPRI5)
const D1 = -12715105075 / 11282082432;
const D3 = 87487479700 / 32700410799;
const D4 = -10690763975 / 1880347072;
const D5 = 701980252875 / 199316789632;
const D6 = -1453857185 / 822651844;
const D7 = 69997945 / 29380423;

/**
 * Adaptive Dormand–Prince 5(4) explicit Runge–Kutta integrator (the DOPRI5 scheme) with local extrapolation,
 * FSAL, PI step-size control (Gustafsson 1991; β = 0.04) and fourth-order dense output (HNW §II.4–II.6).
 * Error control per step: RMS over components of err_i / (atol + rtol·max(|y_i|, |y_i,new|)) <= 1.
 * Validity: non-stiff or mildly stiff problems (stiff problems force tiny steps: use bdf2). Throws Error when
 * maxSteps is exceeded or the step size underflows.
 * @param {(t: number, y: Float64Array) => ArrayLike<number>} f right-hand side dy/dt
 * @param {number} t0 start time [s]
 * @param {ArrayLike<number>} y0 initial state
 * @param {number} t1 end time [s] (may be < t0)
 * @param {{ rtol?: number, atol?: number | ArrayLike<number>, hMax?: number, hInit?: number, maxSteps?: number,
 *   dense?: ArrayLike<number> }} [opts] rtol: relative tolerance (default 1e-6); atol: absolute tolerance, scalar
 *   or per component, in state units (default 1e-9); hMax: largest step [s]; hInit: first step [s] (automatic if
 *   omitted); maxSteps: cap on attempted steps (default 1e6); dense: output times [s], monotonic from t0 towards
 *   t1 and inside [t0, t1] — when given, the result is reported at exactly these times by dense output
 * @returns {{ t: Float64Array, y: Float64Array[], steps: number, rejected: number, nfev: number }} t: output times
 *   [s] (accepted step ends including t0, or the dense times); y[k]: state at t[k]; steps: accepted steps;
 *   rejected: rejected steps; nfev: number of f evaluations
 */
export function rk45(f, t0, y0, t1, { rtol = 1e-6, atol = 1e-9, hMax, hInit, maxSteps = 1e6, dense } = {}) {
  checkOdeArgs('rk45', f, t0, y0, t1);
  if (!(rtol >= 0)) throw new RangeError('rk45: rtol must be >= 0');
  const n = y0.length;
  const atolVec = isVec(atol);
  if (atolVec && atol.length !== n) throw new RangeError('rk45: atol array must have the length of y0');
  if (!atolVec && !(atol >= 0)) throw new RangeError('rk45: atol must be >= 0');
  if (!atolVec && atol === 0 && rtol === 0) throw new RangeError('rk45: rtol and atol must not both be zero');
  const span = t1 - t0;
  const dir = span >= 0 ? 1 : -1;
  const hmax = hMax !== undefined ? Math.abs(hMax) : Math.abs(span);
  if (hMax !== undefined && !(hmax > 0)) throw new RangeError('rk45: hMax must be > 0');
  if (dense !== undefined) {
    if (!isVec(dense)) throw new RangeError('rk45: dense must be an array of output times');
    for (let i = 0; i < dense.length; i++) {
      const td = dense[i];
      if (!((td - t0) * dir >= 0 && (t1 - td) * dir >= 0)) throw new RangeError(`rk45: dense[${i}] is outside [t0, t1]`);
      if (i > 0 && (td - dense[i - 1]) * dir < 0) throw new RangeError('rk45: dense times must be monotonic from t0 towards t1');
    }
  }
  const sc = (i, a, b) => (atolVec ? atol[i] : atol) + rtol * Math.max(Math.abs(a), Math.abs(b));

  let y = Float64Array.from(y0);
  let ynew = new Float64Array(n);
  const ytmp = new Float64Array(n);
  let k1 = new Float64Array(n);
  const k2 = new Float64Array(n);
  const k3 = new Float64Array(n);
  const k4 = new Float64Array(n);
  const k5 = new Float64Array(n);
  const k6 = new Float64Array(n);
  let k7 = new Float64Array(n);
  let nfev = 0;
  const evalF = (t, yy, out) => {
    const d = f(t, yy);
    nfev++;
    if (!d || d.length !== n) throw new RangeError('rk45: f must return an array of the same length as y');
    for (let i = 0; i < n; i++) out[i] = d[i];
  };

  const tOut = [];
  const yOut = [];
  let iDense = 0;
  if (dense) {
    while (iDense < dense.length && dense[iDense] === t0) {
      tOut.push(t0);
      yOut.push(Float64Array.from(y));
      iDense++;
    }
  } else {
    tOut.push(t0);
    yOut.push(Float64Array.from(y));
  }
  if (span === 0) return { t: Float64Array.from(tOut), y: yOut, steps: 0, rejected: 0, nfev: 0 };

  let t = t0;
  evalF(t, y, k1);
  // starting step (HNW §II.4, algorithm "hinit")
  let h;
  if (hInit !== undefined) {
    if (!(Math.abs(hInit) > 0)) throw new RangeError('rk45: hInit must be non-zero');
    h = Math.min(Math.abs(hInit), hmax, Math.abs(span));
  } else {
    let d0 = 0;
    let d1 = 0;
    for (let i = 0; i < n; i++) {
      const s = sc(i, y[i], y[i]) || 1e-300;
      d0 += (y[i] / s) * (y[i] / s);
      d1 += (k1[i] / s) * (k1[i] / s);
    }
    d0 = Math.sqrt(d0 / n);
    d1 = Math.sqrt(d1 / n);
    let h0 = d0 < 1e-5 || d1 < 1e-5 ? 1e-6 : 0.01 * (d0 / d1);
    h0 = Math.min(h0, hmax, Math.abs(span));
    for (let i = 0; i < n; i++) ytmp[i] = y[i] + dir * h0 * k1[i];
    evalF(t + dir * h0, ytmp, k2);
    let d2 = 0;
    for (let i = 0; i < n; i++) {
      const s = sc(i, y[i], y[i]) || 1e-300;
      const v = (k2[i] - k1[i]) / s;
      d2 += v * v;
    }
    d2 = Math.sqrt(d2 / n) / h0;
    const dm = Math.max(d1, d2);
    const h1 = dm <= 1e-15 ? Math.max(1e-6, h0 * 1e-3) : Math.pow(0.01 / dm, 0.2);
    h = Math.min(100 * h0, h1, hmax, Math.abs(span));
  }
  h *= dir;

  const beta = 0.04;
  const expo1 = 0.2 - beta * 0.75;
  const safe = 0.9;
  const facc1 = 5; // maximum step reduction factor
  const facc2 = 0.1; // 1 / maximum step growth factor
  let facold = 1e-4;
  let steps = 0;
  let rejected = 0;
  let lastRejected = false;
  let done = false;
  while (!done) {
    if (steps + rejected >= maxSteps) throw new Error(`rk45: maxSteps (${maxSteps}) exceeded at t = ${t}`);
    if (Math.abs(h) <= 16 * EPS * Math.max(Math.abs(t), Math.abs(t1))) {
      throw new Error(`rk45: step size underflow at t = ${t} (problem too stiff or singular)`);
    }
    let last = false;
    if ((t + 1.01 * h - t1) * dir > 0) {
      h = t1 - t;
      last = true;
    }
    for (let i = 0; i < n; i++) ytmp[i] = y[i] + h * A21 * k1[i];
    evalF(t + 0.2 * h, ytmp, k2);
    for (let i = 0; i < n; i++) ytmp[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);
    evalF(t + 0.3 * h, ytmp, k3);
    for (let i = 0; i < n; i++) ytmp[i] = y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
    evalF(t + 0.8 * h, ytmp, k4);
    for (let i = 0; i < n; i++) ytmp[i] = y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
    evalF(t + (8 / 9) * h, ytmp, k5);
    for (let i = 0; i < n; i++) {
      ytmp[i] = y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i]);
    }
    const tn = last ? t1 : t + h;
    evalF(tn, ytmp, k6);
    for (let i = 0; i < n; i++) {
      ynew[i] = y[i] + h * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i] + B6 * k6[i]);
    }
    evalF(tn, ynew, k7);
    let err = 0;
    for (let i = 0; i < n; i++) {
      const e = h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i] + E7 * k7[i]);
      const v = e / (sc(i, y[i], ynew[i]) || 1e-300);
      err += v * v;
    }
    err = Math.sqrt(err / n);
    if (!Number.isFinite(err)) {
      // overflow inside the step: retry with a much smaller step
      rejected++;
      lastRejected = true;
      h *= 0.1;
      continue;
    }
    const fac11 = Math.pow(err, expo1);
    if (err <= 1) {
      steps++;
      let fac = fac11 / Math.pow(facold, beta);
      fac = Math.max(facc2, Math.min(facc1, fac / safe));
      let hnew = h / fac;
      facold = Math.max(err, 1e-4);
      if (dense) {
        // all requested times in (t, tn]
        while (iDense < dense.length && (dense[iDense] - tn) * dir <= 0) {
          const td = dense[iDense++];
          const out = new Float64Array(n);
          if (td === tn) {
            out.set(ynew);
          } else {
            const th = (td - t) / h;
            const th1 = 1 - th;
            for (let i = 0; i < n; i++) {
              const ydiff = ynew[i] - y[i];
              const bspl = h * k1[i] - ydiff;
              const rc4 = ydiff - h * k7[i] - bspl;
              const rc5 = h * (D1 * k1[i] + D3 * k3[i] + D4 * k4[i] + D5 * k5[i] + D6 * k6[i] + D7 * k7[i]);
              out[i] = y[i] + th * (ydiff + th1 * (bspl + th * (rc4 + th1 * rc5)));
            }
          }
          tOut.push(td);
          yOut.push(out);
        }
      } else {
        tOut.push(tn);
        yOut.push(Float64Array.from(ynew));
      }
      t = tn;
      const sw = y;
      y = ynew;
      ynew = sw;
      const kw = k1;
      k1 = k7;
      k7 = kw;
      if (last) {
        done = true;
      } else {
        if (Math.abs(hnew) > hmax) hnew = dir * hmax;
        if (lastRejected) hnew = dir * Math.min(Math.abs(hnew), Math.abs(h));
        lastRejected = false;
        h = hnew;
      }
    } else {
      rejected++;
      lastRejected = true;
      h /= Math.min(facc1, fac11 / safe);
    }
  }
  if (dense) {
    // guard against round-off leaving trailing requested times (== t1) unserved
    while (iDense < dense.length) {
      tOut.push(dense[iDense++]);
      yOut.push(Float64Array.from(y));
    }
  }
  return { t: Float64Array.from(tOut), y: yOut, steps, rejected, nfev };
}

/** Finite-difference Jacobian of f(t, ·) at y (forward differences), returned as rows. */
function fdJacobian(f, t, y, f0, n) {
  const J = new Array(n);
  for (let i = 0; i < n; i++) J[i] = new Float64Array(n);
  const yp = Float64Array.from(y);
  const sq = Math.sqrt(EPS);
  for (let j = 0; j < n; j++) {
    const h = sq * Math.max(Math.abs(y[j]), 1);
    yp[j] = y[j] + h;
    const fp = f(t, yp);
    for (let i = 0; i < n; i++) J[i][j] = (fp[i] - f0[i]) / h;
    yp[j] = y[j];
  }
  return J;
}

/**
 * Solve y − g·h·f(t, y) = rhs for y by simplified Newton iteration (Jacobian frozen within the step and refreshed
 * if convergence is slow).
 */
function implicitStage(name, f, jac, t, yGuess, rhs, gh, tol, n, step) {
  const y = Float64Array.from(yGuess);
  let lu = null;
  const res = new Array(n);
  for (let it = 0; it < 60; it++) {
    const fr = f(t, y);
    if (!fr || fr.length !== n) throw new RangeError(`${name}: f must return an array of the same length as y`);
    const fy = Float64Array.from(fr);
    if (it % 6 === 0) {
      const J = jac ? jac(t, y) : fdJacobian(f, t, y, fy, n);
      const A = new Array(n);
      for (let i = 0; i < n; i++) {
        const row = new Float64Array(n);
        for (let j = 0; j < n; j++) row[j] = -gh * J[i][j];
        row[i] += 1;
        A[i] = row;
      }
      lu = luFactor(A);
    }
    let rmax = 0;
    for (let i = 0; i < n; i++) {
      res[i] = y[i] - gh * fy[i] - rhs[i];
      if (!Number.isFinite(res[i])) throw new Error(`${name}: non-finite residual at step ${step}`);
      rmax = Math.max(rmax, Math.abs(res[i]));
    }
    const dy = luSolve(lu, res);
    let dmax = 0;
    let ymax = 0;
    for (let i = 0; i < n; i++) {
      y[i] -= dy[i];
      dmax = Math.max(dmax, Math.abs(dy[i]));
      ymax = Math.max(ymax, Math.abs(y[i]));
    }
    if (dmax <= tol * (1 + ymax)) return y;
  }
  throw new Error(`${name}: Newton iteration did not converge at step ${step} (t = ${t})`);
}

function implicitFixedStep(name, order, f, t0, y0, t1, nSteps, opts) {
  checkOdeArgs(name, f, t0, y0, t1);
  checkSteps(name, nSteps);
  const { jac, tol = 1e-10 } = opts || {};
  if (!(tol > 0)) throw new RangeError(`${name}: tol must be > 0`);
  const n = y0.length;
  const h = (t1 - t0) / nSteps;
  const t = new Float64Array(nSteps + 1);
  const ys = new Array(nSteps + 1);
  t[0] = t0;
  ys[0] = Float64Array.from(y0);
  const rhs = new Float64Array(n);
  for (let k = 0; k < nSteps; k++) {
    const tn = k === nSteps - 1 ? t1 : t0 + (k + 1) * h;
    const yk = ys[k];
    if (order === 1 || k === 0) {
      ys[k + 1] = implicitStage(name, f, jac, tn, yk, yk, h, tol, n, k + 1);
    } else {
      const ym = ys[k - 1];
      for (let i = 0; i < n; i++) rhs[i] = (4 * yk[i] - ym[i]) / 3;
      ys[k + 1] = implicitStage(name, f, jac, tn, yk, rhs, (2 / 3) * h, tol, n, k + 1);
    }
    t[k + 1] = tn;
  }
  return { t, y: ys };
}

/**
 * Implicit (backward) Euler with fixed step: y_{k+1} = y_k + h f(t_{k+1}, y_{k+1}). First order, L-stable.
 * The nonlinear stage equation is solved by simplified Newton iteration with an LU-factorised iteration matrix
 * I − h·J (J analytic if `jac` is given, forward finite differences otherwise).
 * Validity: stiff or non-stiff systems of moderate size (dense Jacobian). Throws Error if Newton fails.
 * @param {(t: number, y: Float64Array) => ArrayLike<number>} f right-hand side dy/dt
 * @param {number} t0 start time [s]
 * @param {ArrayLike<number>} y0 initial state
 * @param {number} t1 end time [s]
 * @param {number} nSteps number of steps (>= 1)
 * @param {{ jac?: (t: number, y: Float64Array) => ArrayLike<number>[], tol?: number }} [opts] jac: Jacobian
 *   ∂f/∂y as rows; tol: Newton tolerance on the max-norm of the update relative to 1 + ‖y‖∞ (default 1e-10)
 * @returns {{ t: Float64Array, y: Float64Array[] }} t[k] times [s]; y[k] state at t[k] (k = 0..nSteps)
 */
export function implicitEuler(f, t0, y0, t1, nSteps, opts) {
  return implicitFixedStep('implicitEuler', 1, f, t0, y0, t1, nSteps, opts);
}

/**
 * Two-step backward differentiation formula with fixed step:
 * y_{k+1} = (4 y_k − y_{k−1})/3 + (2h/3) f(t_{k+1}, y_{k+1}). Second order, L-stable (HW II §V.1).
 * The first step is a backward Euler step (its O(h²) local error does not reduce the global order).
 * Validity: stiff systems of moderate size (dense Jacobian). Throws Error if Newton fails.
 * @param {(t: number, y: Float64Array) => ArrayLike<number>} f right-hand side dy/dt
 * @param {number} t0 start time [s]
 * @param {ArrayLike<number>} y0 initial state
 * @param {number} t1 end time [s]
 * @param {number} nSteps number of steps (>= 1)
 * @param {{ jac?: (t: number, y: Float64Array) => ArrayLike<number>[], tol?: number }} [opts] as for implicitEuler
 * @returns {{ t: Float64Array, y: Float64Array[] }} t[k] times [s]; y[k] state at t[k] (k = 0..nSteps)
 */
export function bdf2(f, t0, y0, t1, nSteps, opts) {
  return implicitFixedStep('bdf2', 2, f, t0, y0, t1, nSteps, opts);
}

// ---------------------------------------------------------------------------------------------------------------
// Structural dynamics
// ---------------------------------------------------------------------------------------------------------------

/** Dense matrix (rows, or a scalar for SDOF) -> flat row-major Float64Array of order n; null/undefined -> zeros. */
function flatMatrix(A, n, name, optional) {
  const out = new Float64Array(n * n);
  if (A === undefined || A === null) {
    if (!optional) throw new RangeError(`${name} is required`);
    return out;
  }
  if (typeof A === 'number') {
    if (n !== 1) throw new RangeError(`${name} must be an n×n matrix`);
    out[0] = A;
    return out;
  }
  if (!Array.isArray(A) || A.length !== n) throw new RangeError(`${name} must be an n×n matrix (n = ${n})`);
  for (let i = 0; i < n; i++) {
    if (!isVec(A[i]) || A[i].length !== n) throw new RangeError(`${name} must be an n×n matrix (row ${i})`);
    for (let j = 0; j < n; j++) {
      const v = A[i][j];
      if (!Number.isFinite(v)) throw new RangeError(`${name} contains a non-finite entry at (${i}, ${j})`);
      out[i * n + j] = v;
    }
  }
  return out;
}

function toVec(v, n, name) {
  if (typeof v === 'number') {
    if (n !== 1) throw new RangeError(`${name} must be an array of length ${n}`);
    return Float64Array.of(v);
  }
  if (!isVec(v) || v.length !== n) throw new RangeError(`${name} must be an array of length ${n}`);
  return Float64Array.from(v);
}

function stateSize(u0, name) {
  if (typeof u0 === 'number') return 1;
  if (!isVec(u0) || u0.length < 1) throw new RangeError(`${name}: u0 must be a number or a non-empty array`);
  return u0.length;
}

/** out = A x for a flat row-major matrix. */
function mulFlat(A, x, out, n) {
  for (let i = 0; i < n; i++) {
    let s = 0;
    const r = i * n;
    for (let j = 0; j < n; j++) s += A[r + j] * x[j];
    out[i] = s;
  }
}

function rowsFromFlat(A, n) {
  const R = new Array(n);
  for (let i = 0; i < n; i++) R[i] = A.subarray(i * n, (i + 1) * n);
  return R;
}

function loadInto(fn, t, out, n, name) {
  if (!fn) {
    out.fill(0);
    return;
  }
  const v = fn(t);
  if (typeof v === 'number') {
    if (n !== 1) throw new RangeError(`${name} must return an array of length ${n}`);
    out[0] = v;
    return;
  }
  if (!v || v.length !== n) throw new RangeError(`${name} must return an array of length ${n}`);
  for (let i = 0; i < n; i++) out[i] = v[i];
}

function newmarkParams(name, beta, gamma, alphaHHT) {
  if (!(alphaHHT <= 0 && alphaHHT >= -1 / 3 - 1e-12)) throw new RangeError(`${name}: alphaHHT must be in [-1/3, 0]`);
  // HHT-α: unconditional stability and second-order accuracy require γ = ½ − α, β = (1 − α)²/4
  const g = gamma !== undefined ? gamma : 0.5 - alphaHHT;
  const b = beta !== undefined ? beta : 0.25 * (1 - alphaHHT) * (1 - alphaHHT);
  if (!(b >= 0)) throw new RangeError(`${name}: beta must be >= 0`);
  if (!(g >= 0)) throw new RangeError(`${name}: gamma must be >= 0`);
  return [b, g];
}

/**
 * Newmark-β / HHT-α direct time integration of the linear system M ü + C u̇ + K u = f(t).
 * Scheme (Hughes 1987 §9.1–9.3), solved for the new acceleration with one LU factorisation of
 * S = M + (1+α)γΔt C + (1+α)βΔt² K:
 *   u_{n+1} = u_n + Δt v_n + Δt²[(½ − β) a_n + β a_{n+1}],  v_{n+1} = v_n + Δt[(1 − γ) a_n + γ a_{n+1}],
 *   M a_{n+1} + (1+α) C v_{n+1} − α C v_n + (1+α) K u_{n+1} − α K u_n = (1+α) f_{n+1} − α f_n.
 * Defaults β = ¼, γ = ½ (average acceleration: unconditionally stable, no numerical damping, period elongation
 * ≈ (ωΔt)²/12). When alphaHHT ≠ 0 and beta/gamma are not given, γ = ½ − α and β = (1 − α)²/4 are used.
 * The initial acceleration follows from equilibrium at t = 0, so M must be non-singular.
 * Validity: linear systems with dense matrices of moderate order (O(n³) factorisation, O(n²) per step);
 * −1/3 <= alphaHHT <= 0; explicit central difference (β = 0) is conditionally stable (Δt < 2/ω_max).
 * @param {{ M: number[][] | number, C?: number[][] | number, K: number[][] | number,
 *   f?: (t: number) => (ArrayLike<number> | number), u0: ArrayLike<number> | number, v0: ArrayLike<number> | number,
 *   dt: number, nSteps: number, beta?: number, gamma?: number, alphaHHT?: number }} p
 *   M mass [kg], C damping [N·s/m] (default 0), K stiffness [N/m] (n×n rows, or scalars for SDOF); f(t) load
 *   vector [N] (default 0); u0 [m], v0 [m/s] initial conditions; dt time step [s]; nSteps number of steps;
 *   beta, gamma Newmark parameters [-]; alphaHHT HHT-α parameter [-]
 * @returns {{ t: Float64Array, u: Float64Array[], v: Float64Array[], a: Float64Array[] }} t[k] [s]; u[k] [m],
 *   v[k] [m/s], a[k] [m/s²] at t[k], k = 0..nSteps
 */
export function newmark({ M, C, K, f, u0, v0, dt, nSteps, beta, gamma, alphaHHT = 0 } = {}) {
  const n = stateSize(u0, 'newmark');
  if (!(dt > 0) || !Number.isFinite(dt)) throw new RangeError('newmark: dt must be > 0');
  checkSteps('newmark', nSteps);
  if (f !== undefined && f !== null && typeof f !== 'function') throw new RangeError('newmark: f must be a function of t');
  const [b, g] = newmarkParams('newmark', beta, gamma, alphaHHT);
  const al = alphaHHT;
  const Mf = flatMatrix(M, n, 'newmark: M');
  const Cf = flatMatrix(C, n, 'newmark: C', true);
  const Kf = flatMatrix(K, n, 'newmark: K');
  const t = new Float64Array(nSteps + 1);
  const U = new Array(nSteps + 1);
  const V = new Array(nSteps + 1);
  const A = new Array(nSteps + 1);
  let u = toVec(u0, n, 'newmark: u0');
  let v = toVec(v0 === undefined ? (n === 1 ? 0 : new Float64Array(n)) : v0, n, 'newmark: v0');
  const fn = new Float64Array(n);
  const fn1 = new Float64Array(n);
  const w1 = new Float64Array(n);
  const w2 = new Float64Array(n);
  const rhs = new Float64Array(n);
  loadInto(f, 0, fn, n, 'newmark: f(t)');
  // initial acceleration from equilibrium
  mulFlat(Cf, v, w1, n);
  mulFlat(Kf, u, w2, n);
  for (let i = 0; i < n; i++) rhs[i] = fn[i] - w1[i] - w2[i];
  let a;
  try {
    a = Float64Array.from(luSolve(luFactor(rowsFromFlat(Mf, n)), rhs));
  } catch (err) {
    throw new RangeError('newmark: M is singular (the initial acceleration cannot be computed)');
  }
  const S = new Array(n);
  const cC = (1 + al) * g * dt;
  const cK = (1 + al) * b * dt * dt;
  for (let i = 0; i < n; i++) {
    const row = new Float64Array(n);
    for (let j = 0; j < n; j++) row[j] = Mf[i * n + j] + cC * Cf[i * n + j] + cK * Kf[i * n + j];
    S[i] = row;
  }
  let lu;
  try {
    lu = luFactor(S);
  } catch (err) {
    throw new RangeError('newmark: effective matrix M + (1+α)γΔt C + (1+α)βΔt² K is singular');
  }
  U[0] = u;
  V[0] = v;
  A[0] = a;
  const up = new Float64Array(n);
  const vp = new Float64Array(n);
  for (let k = 0; k < nSteps; k++) {
    const tn1 = (k + 1) * dt;
    loadInto(f, tn1, fn1, n, 'newmark: f(t)');
    for (let i = 0; i < n; i++) {
      const upi = u[i] + dt * v[i] + (0.5 - b) * dt * dt * a[i];
      const vpi = v[i] + (1 - g) * dt * a[i];
      // arguments of K and C in the α-weighted equilibrium, without the a_{n+1} contribution
      up[i] = (1 + al) * upi - al * u[i];
      vp[i] = (1 + al) * vpi - al * v[i];
      w1[i] = upi;
      w2[i] = vpi;
    }
    const upred = Float64Array.from(w1);
    const vpred = Float64Array.from(w2);
    mulFlat(Cf, vp, w1, n);
    mulFlat(Kf, up, w2, n);
    for (let i = 0; i < n; i++) rhs[i] = (1 + al) * fn1[i] - al * fn[i] - w1[i] - w2[i];
    const an = Float64Array.from(luSolve(lu, rhs));
    for (let i = 0; i < n; i++) {
      upred[i] += b * dt * dt * an[i];
      vpred[i] += g * dt * an[i];
    }
    u = upred;
    v = vpred;
    a = an;
    fn.set(fn1);
    t[k + 1] = tn1;
    U[k + 1] = u;
    V[k + 1] = v;
    A[k + 1] = a;
  }
  return { t, u: U, v: V, a: A };
}

/**
 * Implicit Newmark-β integration of the nonlinear system M ü + C u̇ + f_int(u, u̇) = f_ext(t) with full
 * Newton–Raphson equilibrium iterations in every step (Hughes 1987 §9.4; Bathe, "Finite Element Procedures"
 * §9.5.2). Unknown per step: a_{n+1}; iteration matrix S = M + γΔt (C + C_t) + βΔt² K_t with the tangent
 * returned by fint. Converged when ‖R‖₂ <= tol·max(‖f_ext‖, ‖f_int‖, ‖M a‖) or the displacement correction is
 * below 1e-3·tol·‖u‖₂.
 * Validity: smooth or piecewise-smooth internal force with a consistent tangent; β > 0. Throws Error naming the
 * step when the iteration does not converge in maxIter iterations.
 * @param {{ M: number[][] | number, C?: number[][] | number,
 *   fint: (u: Float64Array, v: Float64Array) => { f: ArrayLike<number> | number, Kt: number[][] | number, Ct?: number[][] | number },
 *   fext?: (t: number) => (ArrayLike<number> | number), u0: ArrayLike<number> | number, v0: ArrayLike<number> | number,
 *   dt: number, nSteps: number, beta?: number, gamma?: number, tol?: number, maxIter?: number }} p
 *   M mass [kg]; C linear viscous damping [N·s/m] (default 0); fint returns the internal (restoring) force f [N],
 *   its tangent stiffness Kt = ∂f/∂u [N/m] and optionally Ct = ∂f/∂u̇ [N·s/m]; fext external load [N] (default
 *   0); u0 [m], v0 [m/s]; dt [s]; nSteps; beta (default ¼), gamma (default ½); tol relative convergence
 *   tolerance (default 1e-8); maxIter Newton iterations per step (default 25)
 * @returns {{ t: Float64Array, u: Float64Array[], v: Float64Array[], a: Float64Array[], iterations: number,
 *   iterationsPerStep: Int32Array }} time histories as for newmark; iterations: total Newton iterations;
 *   iterationsPerStep[k]: iterations used in step k+1
 */
export function newmarkNonlinear({ M, C, fint, fext, u0, v0, dt, nSteps, beta = 0.25, gamma = 0.5, tol = 1e-8, maxIter = 25 } = {}) {
  const n = stateSize(u0, 'newmarkNonlinear');
  if (typeof fint !== 'function') throw new RangeError('newmarkNonlinear: fint must be a function (u, v) -> { f, Kt }');
  if (fext !== undefined && fext !== null && typeof fext !== 'function') throw new RangeError('newmarkNonlinear: fext must be a function of t');
  if (!(dt > 0) || !Number.isFinite(dt)) throw new RangeError('newmarkNonlinear: dt must be > 0');
  checkSteps('newmarkNonlinear', nSteps);
  if (!(beta > 0)) throw new RangeError('newmarkNonlinear: beta must be > 0');
  if (!(gamma >= 0)) throw new RangeError('newmarkNonlinear: gamma must be >= 0');
  if (!(tol > 0)) throw new RangeError('newmarkNonlinear: tol must be > 0');
  const Mf = flatMatrix(M, n, 'newmarkNonlinear: M');
  const Cf = flatMatrix(C, n, 'newmarkNonlinear: C', true);
  const t = new Float64Array(nSteps + 1);
  const U = new Array(nSteps + 1);
  const V = new Array(nSteps + 1);
  const A = new Array(nSteps + 1);
  const iterationsPerStep = new Int32Array(nSteps);
  let u = toVec(u0, n, 'newmarkNonlinear: u0');
  let v = toVec(v0 === undefined ? (n === 1 ? 0 : new Float64Array(n)) : v0, n, 'newmarkNonlinear: v0');
  const fe = new Float64Array(n);
  const fi = new Float64Array(n);
  const w1 = new Float64Array(n);
  const w2 = new Float64Array(n);
  const R = new Float64Array(n);
  const norm2 = (x) => {
    let s = 0;
    for (let i = 0; i < n; i++) s += x[i] * x[i];
    return Math.sqrt(s);
  };
  const internal = (uu, vv) => {
    const r = fint(uu, vv);
    if (!r || r.f === undefined || r.Kt === undefined) throw new RangeError('newmarkNonlinear: fint must return { f, Kt }');
    loadInto(() => r.f, 0, fi, n, 'newmarkNonlinear: fint().f');
    return r;
  };
  loadInto(fext, 0, fe, n, 'newmarkNonlinear: fext(t)');
  internal(u, v);
  mulFlat(Cf, v, w1, n);
  for (let i = 0; i < n; i++) R[i] = fe[i] - w1[i] - fi[i];
  let a;
  try {
    a = Float64Array.from(luSolve(luFactor(rowsFromFlat(Mf, n)), R));
  } catch (err) {
    throw new RangeError('newmarkNonlinear: M is singular (the initial acceleration cannot be computed)');
  }
  U[0] = u;
  V[0] = v;
  A[0] = a;
  let total = 0;
  const bdt2 = beta * dt * dt;
  const gdt = gamma * dt;
  for (let k = 0; k < nSteps; k++) {
    const tn1 = (k + 1) * dt;
    loadInto(fext, tn1, fe, n, 'newmarkNonlinear: fext(t)');
    const un = new Float64Array(n);
    const vn = new Float64Array(n);
    const an = Float64Array.from(a); // start from the previous acceleration
    const upred = new Float64Array(n);
    const vpred = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      upred[i] = u[i] + dt * v[i] + (0.5 - beta) * dt * dt * a[i];
      vpred[i] = v[i] + (1 - gamma) * dt * a[i];
      un[i] = upred[i] + bdt2 * an[i];
      vn[i] = vpred[i] + gdt * an[i];
    }
    let converged = false;
    let it = 0;
    for (;;) {
      const r = internal(un, vn);
      mulFlat(Mf, an, w1, n);
      mulFlat(Cf, vn, w2, n);
      for (let i = 0; i < n; i++) R[i] = fe[i] - w1[i] - w2[i] - fi[i];
      const rn = norm2(R);
      if (!Number.isFinite(rn)) throw new Error(`newmarkNonlinear: non-finite residual at step ${k + 1} (t = ${tn1} s)`);
      const ref = Math.max(norm2(fe), norm2(fi), norm2(w1));
      if (rn <= tol * ref) {
        converged = true;
        break;
      }
      if (it >= maxIter) break;
      it++;
      const Kt = flatMatrix(r.Kt, n, 'newmarkNonlinear: fint().Kt');
      const Ct = r.Ct !== undefined && r.Ct !== null ? flatMatrix(r.Ct, n, 'newmarkNonlinear: fint().Ct') : null;
      const S = new Array(n);
      for (let i = 0; i < n; i++) {
        const row = new Float64Array(n);
        for (let j = 0; j < n; j++) {
          const q = i * n + j;
          row[j] = Mf[q] + gdt * (Cf[q] + (Ct ? Ct[q] : 0)) + bdt2 * Kt[q];
        }
        S[i] = row;
      }
      let da;
      try {
        da = luSolve(luFactor(S), R);
      } catch (err) {
        throw new Error(`newmarkNonlinear: singular iteration matrix at step ${k + 1} (t = ${tn1} s)`);
      }
      let du2 = 0;
      let u2 = 0;
      for (let i = 0; i < n; i++) {
        an[i] += da[i];
        un[i] = upred[i] + bdt2 * an[i];
        vn[i] = vpred[i] + gdt * an[i];
        du2 += bdt2 * da[i] * (bdt2 * da[i]);
        u2 += un[i] * un[i];
      }
      if (Math.sqrt(du2) <= tol * Math.sqrt(u2) * 1e-3) {
        // correction far below the tolerance: accept without a further residual evaluation being decisive
        converged = true;
        break;
      }
    }
    if (!converged) {
      throw new Error(`newmarkNonlinear: equilibrium iteration did not converge at step ${k + 1} (t = ${tn1} s)`);
    }
    total += it;
    iterationsPerStep[k] = it;
    u = un;
    v = vn;
    a = an;
    t[k + 1] = tn1;
    U[k + 1] = u;
    V[k + 1] = v;
    A[k + 1] = a;
  }
  return { t, u: U, v: V, a: A, iterations: total, iterationsPerStep };
}
