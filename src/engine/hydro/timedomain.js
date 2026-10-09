// Time-domain rigid-body simulation with radiation memory (Cummins equation), wave excitation by superposition,
// free-decay simulation and decay-test analysis.
//   (M + A_inf) x'' + ∫₀ᵗ K(t - tau) x'(tau) d tau + B_lin x' + C x = F_wave(t) + F_ext(t, x, x')

import { makeRng } from '../math/rng.js';
import { needPos } from './_wavehelp.js';

const isMatrix = (v) => Array.isArray(v) || (v && typeof v.length === 'number' && typeof v[0] !== 'number');

/**
 * Retardation (impulse-response) function from radiation damping by the cosine transform
 *   K(t) = (2/pi) ∫₀^∞ B(omega) cos(omega t) d omega,
 * integrated exactly for a piecewise-linear B(omega) (Filon-type, no aliasing at large t). B is taken as zero
 * outside the tabulated range, and is extended linearly to B(0) = 0 below the first frequency.
 * Reference: Cummins (1962) The impulse response function and ship motions, Schiffstechnik 9; Ogilvie (1964).
 * The damping table must extend to frequencies where B has decayed to a small fraction of its peak.
 * @param {object} o
 * @param {ArrayLike<number>} o.omega rad/s, ascending
 * @param {ArrayLike<number> | number[][][]} o.B damping per frequency: scalars (N·s/m) or n×n matrices
 * @param {number} o.dt lag step, s
 * @param {number} o.tMax kernel length, s
 * @returns {{dt:number, t:Float64Array, n:number, K:Float64Array[]}} K[lag] is the n×n kernel stored row-major
 *   (n = 1 for scalar input), N/m (per unit of B·rad/s)
 */
export function retardationKernel({ omega, B, dt, tMax }) {
  needPos(dt, 'dt'); needPos(tMax, 'tMax');
  const nw = omega.length;
  if (nw < 2 || B.length !== nw) throw new RangeError('omega and B must have equal length >= 2');
  const scalar = !isMatrix(B[0]);
  const n = scalar ? 1 : B[0].length;
  const nn = n * n;
  const get = (iw, e) => (scalar ? B[iw] : B[iw][(e / n) | 0][e % n]);
  const w = [0]; const idx = [-1];
  for (let i = 0; i < nw; i++) {
    if (!(omega[i] > w[w.length - 1])) throw new RangeError('omega must be positive and ascending');
    w.push(omega[i]); idx.push(i);
  }
  const nLag = Math.floor(tMax / dt + 1e-9) + 1;
  const t = new Float64Array(nLag);
  const K = [];
  for (let l = 0; l < nLag; l++) {
    const tt = l * dt;
    t[l] = tt;
    const row = new Float64Array(nn);
    for (let s = 1; s < w.length; s++) {
      const w0 = w[s - 1]; const w1 = w[s];
      const h = w1 - w0;
      let c0; let c1; // weights of B at w0 and w1 in ∫ B cos(w t) dw over the segment
      if (tt * h < 1e-4) {
        const cm = Math.cos(0.5 * (w0 + w1) * tt);
        c0 = 0.5 * h * cm; c1 = 0.5 * h * cm;
      } else {
        const s0 = Math.sin(w0 * tt); const s1 = Math.sin(w1 * tt);
        const k0 = Math.cos(w0 * tt); const k1 = Math.cos(w1 * tt);
        // ∫ (b0 + m (w - w0)) cos(w t) dw = [b sin(w t)/t] + m (cos(w1 t) - cos(w0 t))/t²
        const q = (k1 - k0) / (tt * tt * h);
        c0 = -s0 / tt - q; c1 = s1 / tt + q;
      }
      for (let e = 0; e < nn; e++) {
        const b0 = idx[s - 1] < 0 ? 0 : get(idx[s - 1], e);
        row[e] += c0 * b0 + c1 * get(idx[s], e);
      }
    }
    for (let e = 0; e < nn; e++) row[e] *= 2 / Math.PI;
    K.push(row);
  }
  return { dt, t, n, K };
}

/**
 * Infinite-frequency added mass from A(omega) and the retardation kernel (Ogilvie relation):
 *   A_inf = A(omega) + (1/omega) ∫₀^∞ K(t) sin(omega t) dt, averaged over the supplied frequencies.
 * Reference: Ogilvie (1964) Recent progress toward the understanding and prediction of ship motions.
 * The kernel must be long enough to have decayed; the scatter of the per-frequency estimates is returned.
 * @param {object} o
 * @param {ArrayLike<number>} o.omega rad/s
 * @param {ArrayLike<number> | number[][][]} o.A added mass per frequency (scalars or n×n matrices)
 * @param {{dt:number, n:number, K:Float64Array[]}} o.kernel from retardationKernel
 * @returns {{Ainf:number|number[][], spread:number}} kg …; spread = largest (max - min) of the estimates relative
 *   to the mean diagonal magnitude
 */
export function infiniteFrequencyAddedMass({ omega, A, kernel }) {
  const nw = omega.length;
  const scalar = !isMatrix(A[0]);
  const n = kernel.n; const nn = n * n;
  const sum = new Float64Array(nn); const mn = new Float64Array(nn).fill(Infinity); const mx = new Float64Array(nn).fill(-Infinity);
  const nLag = kernel.K.length;
  for (let iw = 0; iw < nw; iw++) {
    const w = omega[iw];
    for (let e = 0; e < nn; e++) {
      let integ = 0;
      for (let l = 1; l < nLag; l++) {
        integ += 0.5 * kernel.dt * (kernel.K[l - 1][e] * Math.sin(w * (l - 1) * kernel.dt) + kernel.K[l][e] * Math.sin(w * l * kernel.dt));
      }
      const a = (scalar ? A[iw] : A[iw][(e / n) | 0][e % n]) + integ / w;
      sum[e] += a; mn[e] = Math.min(mn[e], a); mx[e] = Math.max(mx[e], a);
    }
  }
  let scale = 0; let spread = 0;
  for (let i = 0; i < n; i++) scale = Math.max(scale, Math.abs(sum[i * n + i] / nw));
  for (let e = 0; e < nn; e++) spread = Math.max(spread, mx[e] - mn[e]);
  if (scalar) return { Ainf: sum[0] / nw, spread: scale > 0 ? spread / scale : 0 };
  const Ainf = [];
  for (let i = 0; i < n; i++) { Ainf.push([]); for (let j = 0; j < n; j++) Ainf[i].push(sum[i * n + j] / nw); }
  return { Ainf, spread: scale > 0 ? spread / scale : 0 };
}

/**
 * Regular wave components representing a spectrum: a_i = sqrt(2 S(omega_i) d omega_i), with seeded uniform
 * random phases (deterministic for a seed).
 * Reference: DNV-RP-C205 (2019) §3 (linear random wave model).
 * The record repeats after 2 pi / d omega when the frequencies are equally spaced.
 * @param {{omega:ArrayLike<number>, S:ArrayLike<number>, seed?:number}} o rad/s, m²·s/rad
 * @returns {{omega:Float64Array, amp:Float64Array, phase:Float64Array}} rad/s, m, rad
 */
export function waveComponents({ omega, S, seed = 1 }) {
  const n = omega.length;
  if (S.length !== n || n < 2) throw new RangeError('omega and S must have equal length >= 2');
  const rng = makeRng(seed);
  const amp = new Float64Array(n); const phase = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const dw = i === 0 ? omega[1] - omega[0] : i === n - 1 ? omega[n - 1] - omega[n - 2] : 0.5 * (omega[i + 1] - omega[i - 1]);
    amp[i] = Math.sqrt(2 * Math.max(S[i], 0) * dw);
    phase[i] = 2 * Math.PI * rng.next();
  }
  return { omega: Float64Array.from(omega), amp, phase };
}

/**
 * First-order wave excitation in the time domain by superposition:
 *   F_j(t) = sum_i a_i Re{ X_j(omega_i) exp(i (omega_i t + phase_i)) }.
 * Reference: linear superposition, DNV-RP-C205 (2019) §7.
 * @param {object} o
 * @param {{omega:ArrayLike<number>, amp:ArrayLike<number>, phase:ArrayLike<number>}} o.components rad/s, m, rad
 * @param {ArrayLike<number>[]} o.Fre excitation per unit amplitude, real part [dof][component], N/m …
 * @param {ArrayLike<number>[]} o.Fim imaginary part
 * @returns {function(number): number[]} force vector as a function of time (s), N …
 */
export function waveExcitationForce({ components, Fre, Fim }) {
  const nc = components.omega.length;
  const nd = Fre.length;
  for (let j = 0; j < nd; j++) if (Fre[j].length !== nc || Fim[j].length !== nc) throw new RangeError('Fre/Fim must have one value per component');
  return (t) => {
    const f = new Array(nd).fill(0);
    for (let i = 0; i < nc; i++) {
      const th = components.omega[i] * t + components.phase[i];
      const c = components.amp[i] * Math.cos(th); const s = components.amp[i] * Math.sin(th);
      for (let j = 0; j < nd; j++) f[j] += Fre[j][i] * c - Fim[j][i] * s;
    }
    return f;
  };
}

function toMatrix(v, n, name, optional) {
  if (v === undefined || v === null) {
    if (!optional) throw new RangeError(`${name} is required`);
    return Array.from({ length: n }, () => new Array(n).fill(0));
  }
  if (typeof v === 'number') { if (n !== 1) throw new RangeError(`${name} must be an n×n matrix`); return [[v]]; }
  if (v.length !== n || v[0].length !== n) throw new RangeError(`${name} must be ${n}×${n}`);
  return v.map((r) => Array.from(r));
}

function invert(Min, n) {
  const a = Min.map((r, i) => [...r, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(a[r][c]) > Math.abs(a[p][c])) p = r;
    if (!(Math.abs(a[p][c]) > 0)) throw new RangeError('M + Ainf is singular');
    [a[c], a[p]] = [a[p], a[c]];
    const d = a[c][c];
    for (let j = 0; j < 2 * n; j++) a[c][j] /= d;
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const f = a[r][c];
      if (f !== 0) for (let j = 0; j < 2 * n; j++) a[r][j] -= f * a[c][j];
    }
  }
  return a.map((r) => r.slice(n));
}

/**
 * Integrates the Cummins equation with a fixed step:
 *   (M + Ainf) x'' + ∫₀ᵗ K(t - tau) x'(tau) d tau + Blin x' + C x = force(t, x, x').
 * The convolution uses the trapezoidal rule on the stored velocity history (kernel lag step = dt); time
 * stepping is a second-order predictor–corrector (explicit predictor, trapezoidal corrector with the memory
 * term re-evaluated). Without a kernel this is the constant-coefficient model
 *   (M + A) x'' + B x' + C x = force, with A passed as Ainf and B as Blin.
 * `force` may be any nonlinear function of time and state (moorings, viscous drag, fenders, wave excitation).
 * Reference: Cummins (1962); Ogilvie (1964); DNV-RP-C205 (2019) §7 (time-domain analysis).
 * Requires dt well below the shortest natural period (about T/40 or finer) and equal to the kernel step.
 * @param {object} o
 * @param {number|number[][]} o.M mass matrix, kg …
 * @param {number|number[][]} [o.Ainf] infinite-frequency (or constant) added mass
 * @param {number|number[][]} [o.Blin] constant linear damping
 * @param {number|number[][]} o.C stiffness
 * @param {{dt:number, n:number, K:Float64Array[]}} [o.kernel] from retardationKernel
 * @param {function(number, number[], number[]): number[]} [o.force] (t, x, v) -> force vector
 * @param {number[]} [o.x0] initial displacement
 * @param {number[]} [o.v0] initial velocity
 * @param {number} o.dt time step, s
 * @param {number} o.nSteps number of steps
 * @returns {{t:Float64Array, x:Float64Array[], v:Float64Array[]}} x[dof][step], v[dof][step]
 */
export function simulateCummins({ M, Ainf, Blin, C, kernel, force, x0, v0, dt, nSteps }) {
  needPos(dt, 'dt');
  if (!(nSteps >= 1)) throw new RangeError('nSteps must be >= 1');
  const n = typeof M === 'number' ? 1 : M.length;
  const Mm = toMatrix(M, n, 'M'); const Am = toMatrix(Ainf, n, 'Ainf', true);
  const Bm = toMatrix(Blin, n, 'Blin', true); const Cm = toMatrix(C, n, 'C');
  if (kernel) {
    if (kernel.n !== n) throw new RangeError('kernel dimension does not match M');
    if (Math.abs(kernel.dt - dt) > 1e-9 * dt) throw new RangeError('dt must equal the kernel lag step');
  }
  const Mt = Mm.map((r, i) => r.map((v, j) => v + Am[i][j]));
  const Minv = invert(Mt, n);
  const nLag = kernel ? kernel.K.length : 0;
  const t = new Float64Array(nSteps + 1);
  const X = Array.from({ length: n }, () => new Float64Array(nSteps + 1));
  const V = Array.from({ length: n }, () => new Float64Array(nSteps + 1));
  let x = (x0 || new Array(n).fill(0)).slice(); let v = (v0 || new Array(n).fill(0)).slice();
  if (x.length !== n || v.length !== n) throw new RangeError('x0 and v0 must have one value per degree of freedom');
  const hist = new Float64Array(n); // memory contribution from past velocities (lags >= 1)
  const mu = new Float64Array(n);
  const rhs = new Float64Array(n);
  // Memory force at step m for a trial current velocity vc: dt [ ½ K0 vc + sum_{l=1}^{m} w_l K_l v_{m-l} ]
  const past = (m) => {
    hist.fill(0);
    if (!kernel) return;
    const lmax = Math.min(m, nLag - 1);
    for (let l = 1; l <= lmax; l++) {
      const wgt = l === m ? 0.5 : 1;
      const Kl = kernel.K[l];
      for (let i = 0; i < n; i++) {
        let s = 0;
        for (let j = 0; j < n; j++) s += Kl[i * n + j] * V[j][m - l];
        hist[i] += wgt * s;
      }
    }
  };
  const accel = (tt, xc, vc, m) => {
    const f = force ? force(tt, xc, vc) : null;
    for (let i = 0; i < n; i++) {
      let s = hist[i];
      if (kernel && m > 0) for (let j = 0; j < n; j++) s += 0.5 * kernel.K[0][i * n + j] * vc[j];
      mu[i] = kernel ? s * dt : 0;
      let r = (f ? f[i] : 0) - mu[i];
      for (let j = 0; j < n; j++) r -= Bm[i][j] * vc[j] + Cm[i][j] * xc[j];
      rhs[i] = r;
    }
    const a = new Array(n);
    for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += Minv[i][j] * rhs[j]; a[i] = s; }
    return a;
  };
  for (let i = 0; i < n; i++) { X[i][0] = x[i]; V[i][0] = v[i]; }
  past(0);
  let a = accel(0, x, v, 0);
  for (let m = 0; m < nSteps; m++) {
    const tn = (m + 1) * dt;
    const xp = x.map((xv, i) => xv + dt * v[i] + 0.5 * dt * dt * a[i]);
    const vp = v.map((vv, i) => vv + dt * a[i]);
    past(m + 1);
    let an = accel(tn, xp, vp, m + 1);
    let xn = xp; let vn = vp;
    for (let c = 0; c < 2; c++) {
      vn = v.map((vv, i) => vv + 0.5 * dt * (a[i] + an[i]));
      xn = x.map((xv, i) => xv + dt * v[i] + 0.25 * dt * dt * (a[i] + an[i]));
      an = accel(tn, xn, vn, m + 1);
    }
    x = xn; v = vn; a = an;
    t[m + 1] = tn;
    for (let i = 0; i < n; i++) { X[i][m + 1] = x[i]; V[i][m + 1] = v[i]; }
  }
  return { t, x: X, v: V };
}

/**
 * Free-decay simulation of one degree of freedom with linear and quadratic damping:
 *   m x'' + [memory] + b1 x' + b2 |x'| x' + c x = 0,  x(0) = x0.
 * Reference: Cummins equation (see simulateCummins); Faltinsen (1990) ch. 7 (free-decay tests).
 * @param {object} o
 * @param {number} o.mass total inertia M + A (or M + Ainf with a kernel), kg or kg·m²
 * @param {number} o.stiffness N/m or N·m/rad
 * @param {number} [o.b1=0] linear damping, N·s/m
 * @param {number} [o.b2=0] quadratic damping, N·s²/m²
 * @param {object} [o.kernel] scalar retardation kernel
 * @param {number} o.x0 initial offset
 * @param {number} o.dt s
 * @param {number} o.duration s
 * @returns {{t:Float64Array, x:Float64Array, v:Float64Array}}
 */
export function freeDecay({ mass, stiffness, b1 = 0, b2 = 0, kernel, x0, dt, duration }) {
  needPos(mass, 'mass'); needPos(stiffness, 'stiffness'); needPos(duration, 'duration');
  const r = simulateCummins({
    M: mass, Blin: b1, C: stiffness, kernel, x0: [x0], dt, nSteps: Math.round(duration / dt),
    force: b2 ? (t, x, v) => [-b2 * Math.abs(v[0]) * v[0]] : undefined,
  });
  return { t: r.t, x: r.x[0], v: r.v[0] };
}

/**
 * Decay-test analysis. Successive extremes X_n (alternating sign, absolute values) of a free-decay record of
 *   x'' + p1 x' + p2 |x'| x' + p3 x = 0
 * satisfy, by energy balance over one cycle with equivalent linearisation,
 *   (2/Tm) ln(X_{n-1}/X_{n+1}) = p1 + (16/3) p2 X_n / Tm,
 * so a straight-line fit of the left side against (16/3) X_n / Tm gives p1 (intercept) and p2 (slope).
 * Also returns the mean logarithmic decrement per cycle and the equivalent damping ratio.
 * Reference: Faltinsen (1990) Sea Loads on Ships and Offshore Structures ch. 7 (damping from free-decay tests);
 * the relation is re-derived here from the energy balance stated above.
 * Needs at least four extremes; valid for light damping (ratio below about 0.2).
 * @param {object} o
 * @param {ArrayLike<number>} o.t s
 * @param {ArrayLike<number>} o.x response (zero mean, or the mean is removed when `removeMean`)
 * @param {boolean} [o.removeMean=false]
 * @param {number} [o.inertia] total inertia M + A to convert p1, p2 into dimensional b1, b2
 * @returns {{period:number, peaks:{t:number, x:number}[], logDecrement:number, dampingRatio:number, p1:number,
 *   p2:number, b1:number|null, b2:number|null, r2:number}} damped period s; log decrement per cycle; equivalent
 *   linear damping ratio from the mean decrement; p1 1/s, p2 1/unit; b1 = p1·inertia, b2 = p2·inertia; r² of the fit
 */
export function decayAnalysis({ t, x, removeMean = false, inertia }) {
  const n = x.length;
  if (n < 8 || t.length !== n) throw new RangeError('t and x must have equal length >= 8');
  let mean = 0;
  if (removeMean) { for (let i = 0; i < n; i++) mean += x[i]; mean /= n; }
  const peaks = [];
  if (Math.abs(x[0] - mean) > Math.abs(x[1] - mean)) peaks.push({ t: t[0], x: x[0] - mean });
  for (let i = 1; i < n - 1; i++) {
    const y0 = x[i - 1] - mean; const y1 = x[i] - mean; const y2 = x[i + 1] - mean;
    if ((y1 > y0 && y1 >= y2 && y1 > 0) || (y1 < y0 && y1 <= y2 && y1 < 0)) {
      // parabolic refinement of the extreme
      const den = y0 - 2 * y1 + y2;
      const dx = den !== 0 ? (0.5 * (y0 - y2)) / den : 0;
      const h = t[i + 1] - t[i];
      peaks.push({ t: t[i] + dx * h, x: y1 - 0.25 * (y0 - y2) * dx });
    }
  }
  if (peaks.length < 4) throw new RangeError('record contains fewer than four extremes');
  const np = peaks.length;
  const period = (2 * (peaks[np - 1].t - peaks[0].t)) / (np - 1);
  const xs = []; const ys = [];
  let dsum = 0;
  for (let i = 1; i < np - 1; i++) {
    const Tm = peaks[i + 1].t - peaks[i - 1].t;
    const ratio = Math.abs(peaks[i - 1].x) / Math.abs(peaks[i + 1].x);
    if (!(ratio > 0) || !Number.isFinite(ratio)) continue;
    const ld = Math.log(ratio);
    dsum += ld;
    xs.push(((16 / 3) * Math.abs(peaks[i].x)) / Tm);
    ys.push((2 / Tm) * ld);
  }
  const m = xs.length;
  const logDecrement = dsum / m;
  let sx = 0; let sy = 0; let sxx = 0; let sxy = 0; let syy = 0;
  for (let i = 0; i < m; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; syy += ys[i] * ys[i]; }
  const det = m * sxx - sx * sx;
  let p1; let p2; let r2;
  if (m >= 2 && det > 1e-12 * m * sxx) {
    p2 = (m * sxy - sx * sy) / det; p1 = (sy - p2 * sx) / m;
    const sst = syy - (sy * sy) / m;
    let sse = 0;
    for (let i = 0; i < m; i++) sse += (ys[i] - p1 - p2 * xs[i]) ** 2;
    r2 = sst > 0 ? 1 - sse / sst : 1;
  } else { p1 = sy / m; p2 = 0; r2 = 1; }
  return {
    period, peaks, logDecrement, dampingRatio: logDecrement / Math.sqrt(4 * Math.PI * Math.PI + logDecrement * logDecrement),
    p1, p2, b1: inertia ? p1 * inertia : null, b2: inertia ? p2 * inertia : null, r2,
  };
}
