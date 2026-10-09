// One-dimensional nonlinear shallow-water solver: finite volumes, HLL flux, MUSCL reconstruction,
// hydrostatic reconstruction (well-balanced, positivity preserving), wetting/drying, semi-implicit Manning
// friction, SSP-RK2 time stepping. SI units.

const G = 9.80665;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

function minmod(a, b) {
  if (a * b <= 0) return 0;
  return Math.abs(a) < Math.abs(b) ? a : b;
}
function vanLeer(a, b) {
  const p = a * b;
  if (p <= 0) return 0;
  return (2 * p) / (a + b);
}
function noLimiter() {
  return 0;
}

function seriesFn(spec, key, name) {
  const v = spec[key];
  if (typeof v === 'function') return v;
  if (typeof v === 'number') return () => v;
  if (v && v.t && v.values) {
    const ts = v.t;
    const ys = v.values;
    need(ts.length === ys.length && ts.length > 0, `${name}.${key} series must have equal-length t and values`);
    let j = 0;
    return (t) => {
      const n = ts.length;
      if (t <= ts[0]) return ys[0];
      if (t >= ts[n - 1]) return ys[n - 1];
      if (ts[j] > t) j = 0;
      while (j < n - 2 && ts[j + 1] < t) j++;
      return ys[j] + ((ys[j + 1] - ys[j]) * (t - ts[j])) / (ts[j + 1] - ts[j]);
    };
  }
  throw new RangeError(`${name}.${key} must be a number, a function of time or { t, values }`);
}

/**
 * Solves the one-dimensional nonlinear shallow-water (Saint-Venant) equations
 *   ∂h/∂t + ∂(hu)/∂x = 0,   ∂(hu)/∂t + ∂(hu² + ½ g h²)/∂x = −g h ∂zb/∂x − g n² u|u| / h^(1/3)
 * on a uniform grid with a Godunov-type finite-volume scheme:
 *   - HLL approximate Riemann solver with Einfeldt wave-speed estimates and dry-bed speeds (Toro 2001);
 *   - MUSCL reconstruction of h, free surface and u with a minmod or van Leer limiter (first order next to
 *     dry cells);
 *   - hydrostatic reconstruction of Audusse et al. (2004), which preserves lake-at-rest exactly, including
 *     at wet/dry fronts, and keeps depths non-negative;
 *   - Manning friction treated point-implicitly; two-stage SSP Runge–Kutta (Shu & Osher 1988); CFL-limited
 *     adaptive time step.
 * Boundary types (each side): 'wall' (reflective), 'transmissive' (zero gradient), 'periodic' (both sides),
 * 'stage' (prescribed free-surface level, velocity from the outgoing Riemann invariant), 'discharge'
 * (prescribed unit discharge, depth from the outgoing invariant), 'wave' (incident wave train entering
 * through an absorbing–generating characteristic boundary).
 * Reference: Toro (2001) "Shock-Capturing Methods for Free-Surface Shallow Flows"; Audusse, Bouchut,
 * Bristeau, Klein & Perthame (2004) SIAM J. Sci. Comput. 25(6); Harten, Lax & van Leer (1983).
 * Validity: hydrostatic, depth-averaged flow (long waves, bores, run-up, dam-break); no dispersion, so
 * short waves steepen into bores; characteristic boundaries assume subcritical flow at the boundary.
 * @param {object} o
 * @param {number} o.dx cell size [m]
 * @param {ArrayLike<number>} o.zb bed elevation at the n cell centres [m]
 * @param {ArrayLike<number>} [o.h0] initial depth per cell [m] (give h0 or eta0)
 * @param {number|ArrayLike<number>} [o.eta0] initial free-surface elevation [m]; depth = max(eta0 − zb, 0)
 * @param {ArrayLike<number>} [o.u0] initial velocity per cell [m/s], default 0
 * @param {number} o.tEnd end time [s]
 * @param {number} [o.x0] coordinate of the left edge of the domain [m], default 0
 * @param {number} [o.cfl] Courant number [-], default 0.45 (must be <= 0.5)
 * @param {number|ArrayLike<number>} [o.manning] Manning n [s/m^(1/3)], default 0
 * @param {'minmod'|'vanleer'|'none'} [o.limiter] slope limiter ('none' = first order), default 'minmod'
 * @param {number} [o.hDry] depth below which a cell is treated as dry [m], default 1e-6
 * @param {object} [o.left] boundary { type, ... }, default { type: 'wall' }
 *   stage: { type:'stage', eta } ; discharge: { type:'discharge', q } (q positive in +x) [m²/s];
 *   wave: { type:'wave', eta, depth } with eta(t) the incident surface elevation [m] and depth the
 *   still-water depth at the boundary [m] — or { amplitude [m], period [s], phase [rad] } for a sine train.
 *   eta and q may be numbers, functions of time [s] or series { t: number[], values: number[] }.
 * @param {object} [o.right] as left
 * @param {number[]} [o.outputTimes] times at which snapshots are stored [s]
 * @param {'right'|'left'|null} [o.landSide] side on which the shoreline is tracked, default 'right'
 * @param {number} [o.hShore] depth defining the wet front for shoreline and inundation output [m], default 1e-3
 * @param {number} [o.g] gravity [m/s²]
 * @param {number} [o.maxSteps] safety limit on the number of time steps, default 5e6
 * @returns {{x:Float64Array,zb:Float64Array,t:number,h:Float64Array,u:Float64Array,eta:Float64Array,
 *   snapshots:Array<{t:number,h:Float64Array,u:Float64Array}>,
 *   shoreline:{t:Float64Array,x:Float64Array,z:Float64Array,maxRunup:number,maxRunupTime:number},
 *   envelope:{hMax:Float64Array,etaMax:Float64Array,uMax:Float64Array,inundationLimit:number},
 *   volume:{initial:number,final:number,boundaryNet:number,clipped:number,error:number,relativeError:number},
 *   steps:number,dtMin:number,dtMax:number}}
 *   x cell centres [m]; final h [m], u [m/s], eta = h + zb [m]; shoreline history: position of the wet
 *   front [m] and its free-surface elevation [m], maxRunup = highest shoreline elevation [m];
 *   envelope: maxima over time per cell (etaMax is NaN where never wet), inundationLimit = furthest
 *   landward wet position [m]; volume per unit width [m²]: boundaryNet = time-integrated inflow minus
 *   outflow, clipped = volume added by positivity clipping, error = final − initial − boundaryNet − clipped
 */
export function solveSwe1d({
  dx, zb, h0, eta0, u0, tEnd, x0 = 0, cfl = 0.45, manning = 0, limiter = 'minmod', hDry = 1e-6,
  left = { type: 'wall' }, right = { type: 'wall' }, outputTimes = [], landSide = 'right', hShore = 1e-3,
  g = G, maxSteps = 5e6,
}) {
  need(Number.isFinite(dx) && dx > 0, 'dx must be positive');
  need(zb && zb.length >= 4, 'zb must have at least 4 cells');
  need(Number.isFinite(tEnd) && tEnd >= 0, 'tEnd must be >= 0');
  need(cfl > 0 && cfl <= 0.5, 'cfl must lie in (0, 0.5]');
  need(hDry > 0, 'hDry must be positive');
  need(h0 !== undefined || eta0 !== undefined, 'h0 or eta0 must be given');
  const lim = { minmod, vanleer: vanLeer, none: noLimiter }[limiter];
  need(lim !== undefined, "limiter must be 'minmod', 'vanleer' or 'none'");
  const n = zb.length;
  const N = n + 4;
  const iL = 2;
  const iR = n + 1;
  const types = ['wall', 'transmissive', 'periodic', 'stage', 'discharge', 'wave'];
  need(types.includes(left.type), `left.type must be one of ${types.join(', ')}`);
  need(types.includes(right.type), `right.type must be one of ${types.join(', ')}`);
  need((left.type === 'periodic') === (right.type === 'periodic'), 'periodic boundaries must be set on both sides');

  const z = new Float64Array(N);
  const h = new Float64Array(N);
  const q = new Float64Array(N);
  const nMan = new Float64Array(N);
  for (let i = 0; i < n; i++) {
    need(Number.isFinite(zb[i]), 'zb must be finite');
    z[i + 2] = zb[i];
    const hi = h0 !== undefined ? h0[i] : Math.max((typeof eta0 === 'number' ? eta0 : eta0[i]) - zb[i], 0);
    need(Number.isFinite(hi) && hi >= 0, 'initial depth must be >= 0');
    h[i + 2] = hi;
    q[i + 2] = u0 !== undefined && hi > hDry ? hi * u0[i] : 0;
    const nm = typeof manning === 'number' ? manning : manning[i];
    need(Number.isFinite(nm) && nm >= 0, 'manning must be >= 0');
    nMan[i + 2] = nm;
  }
  let hasFriction = false;
  for (let i = iL; i <= iR; i++) if (nMan[i] > 0) hasFriction = true;
  // Ghost bed levels: mirrored at walls (so that reflections are symmetric), wrapped when periodic,
  // otherwise extended flat.
  const setGhostBed = (type, g1, g2, in1, in2, w1, w2) => {
    if (type === 'wall') {
      z[g1] = z[in1];
      z[g2] = z[in2];
    } else if (type === 'periodic') {
      z[g1] = z[w1];
      z[g2] = z[w2];
    } else {
      z[g1] = z[in1];
      z[g2] = z[in1];
    }
  };
  setGhostBed(left.type, 1, 0, iL, iL + 1, iR, iR - 1);
  setGhostBed(right.type, iR + 1, iR + 2, iR, iR - 1, iL, iL + 1);

  const prep = (b, name, edge) => {
    const r = { type: b.type };
    if (b.type === 'stage') r.eta = seriesFn(b, 'eta', name);
    else if (b.type === 'discharge') r.q = seriesFn(b, 'q', name);
    else if (b.type === 'wave') {
      if (b.eta !== undefined) r.eta = seriesFn(b, 'eta', name);
      else {
        need(Number.isFinite(b.amplitude) && b.period > 0, `${name} wave needs eta or amplitude and period`);
        const om = (2 * Math.PI) / b.period;
        const ph = b.phase || 0;
        const a = b.amplitude;
        r.eta = (t) => a * Math.sin(om * t + ph);
      }
      r.depth = b.depth !== undefined ? b.depth : h[edge];
      need(r.depth > 0, `${name}.depth must be positive`);
      r.level = b.level !== undefined ? b.level : z[edge] + r.depth;
    }
    return r;
  };
  const bl = prep(left, 'left', iL);
  const br = prep(right, 'right', iR);

  // sgn = +1 for the left boundary (outgoing invariant u − 2c), −1 for the right (u + 2c).
  const applyOne = (b, hh, qq, t, g1, g2, in1, in2, w1, w2, sgn) => {
    switch (b.type) {
      case 'wall':
        hh[g1] = hh[in1];
        qq[g1] = -qq[in1];
        hh[g2] = hh[in2];
        qq[g2] = -qq[in2];
        return;
      case 'transmissive':
        hh[g1] = hh[in1];
        qq[g1] = qq[in1];
        hh[g2] = hh[in1];
        qq[g2] = qq[in1];
        return;
      case 'periodic':
        hh[g1] = hh[w1];
        qq[g1] = qq[w1];
        hh[g2] = hh[w2];
        qq[g2] = qq[w2];
        return;
      default:
    }
    const hi = hh[in1];
    const ui = hi > hDry ? qq[in1] / hi : 0;
    const ci = Math.sqrt(g * hi);
    const Rout = ui - sgn * 2 * ci; // invariant carried from the interior toward the boundary
    let hb;
    let ub;
    if (b.type === 'stage') {
      hb = Math.max(b.eta(t) - z[in1], 0);
      ub = hi > hDry ? Rout + sgn * 2 * Math.sqrt(g * hb) : 0;
    } else if (b.type === 'discharge') {
      const qb = b.q(t);
      hb = Math.max(hi, hDry * 10);
      // Solve qb/hb − sgn·2√(g hb) = Rout by Newton iteration starting from the interior depth.
      for (let it = 0; it < 30; it++) {
        const f = qb / hb - sgn * 2 * Math.sqrt(g * hb) - Rout;
        const df = -qb / (hb * hb) - sgn * Math.sqrt(g / hb);
        let hn = hb - f / df;
        if (!(hn > 0)) hn = 0.5 * hb;
        const done = Math.abs(hn - hb) < 1e-13 * hb;
        hb = hn;
        if (done) break;
      }
      ub = qb / hb;
    } else {
      // Incident simple wave: u_i = ±2 (c_i − c_0); its incoming invariant combined with the outgoing one.
      const c0 = Math.sqrt(g * b.depth);
      const cInc = Math.sqrt(g * Math.max(b.depth + b.eta(t), 0));
      const Rin = sgn * (4 * cInc - 2 * c0);
      const cb = (sgn * (Rin - Rout)) / 4;
      hb = cb > 0 ? (cb * cb) / g : 0;
      ub = 0.5 * (Rin + Rout);
    }
    hh[g1] = hb;
    qq[g1] = hb * ub;
    hh[g2] = hb;
    qq[g2] = hb * ub;
  };
  const applyBc = (hh, qq, t) => {
    applyOne(bl, hh, qq, t, 1, 0, iL, iL + 1, iR, iR - 1, 1);
    applyOne(br, hh, qq, t, iR + 1, iR + 2, iR, iR - 1, iL, iL + 1, -1);
  };

  const u = new Float64Array(N);
  const hm = new Float64Array(N);
  const hp = new Float64Array(N);
  const zm = new Float64Array(N);
  const zp = new Float64Array(N);
  const um = new Float64Array(N);
  const up = new Float64Array(N);
  const Fh = new Float64Array(N + 1);
  const FqL = new Float64Array(N + 1);
  const FqR = new Float64Array(N + 1);
  const dh = new Float64Array(N);
  const dq = new Float64Array(N);
  const h1 = new Float64Array(N);
  const q1 = new Float64Array(N);
  const halfG = 0.5 * g;

  const rhs = (hh, qq, t) => {
    applyBc(hh, qq, t);
    for (let i = 0; i < N; i++) u[i] = hh[i] > hDry ? qq[i] / hh[i] : 0;
    for (let i = 0; i < N; i++) {
      let sh = 0;
      let se = 0;
      let su = 0;
      if (i > 0 && i < N - 1 && hh[i] > hDry && hh[i - 1] > hDry && hh[i + 1] > hDry) {
        sh = 0.5 * lim(hh[i] - hh[i - 1], hh[i + 1] - hh[i]);
        se = 0.5 * lim(hh[i] + z[i] - hh[i - 1] - z[i - 1], hh[i + 1] + z[i + 1] - hh[i] - z[i]);
        su = 0.5 * lim(u[i] - u[i - 1], u[i + 1] - u[i]);
      }
      const e = hh[i] + z[i];
      hm[i] = hh[i] - sh;
      hp[i] = hh[i] + sh;
      zm[i] = e - se - hm[i];
      zp[i] = e + se - hp[i];
      um[i] = u[i] - su;
      up[i] = u[i] + su;
    }
    for (let f = iL; f <= iR + 1; f++) {
      const a = f - 1;
      const zs = zp[a] > zm[f] ? zp[a] : zm[f];
      let hL = hp[a] + zp[a] - zs;
      let hR = hm[f] + zm[f] - zs;
      if (hL < 0) hL = 0;
      if (hR < 0) hR = 0;
      let fh = 0;
      let fq = 0;
      if (hL > 0 || hR > 0) {
        const uL = up[a];
        const uR = um[f];
        const cL = Math.sqrt(g * hL);
        const cR = Math.sqrt(g * hR);
        let sL;
        let sR;
        if (hL <= 0) {
          sL = uR - 2 * cR;
          sR = uR + cR;
        } else if (hR <= 0) {
          sL = uL - cL;
          sR = uL + 2 * cL;
        } else {
          const us = 0.5 * (uL + uR) + cL - cR;
          const cs = 0.5 * (cL + cR) + 0.25 * (uL - uR);
          sL = Math.min(uL - cL, us - cs);
          sR = Math.max(uR + cR, us + cs);
        }
        const qL = hL * uL;
        const qR = hR * uR;
        if (sL >= 0) {
          fh = qL;
          fq = qL * uL + halfG * hL * hL;
        } else if (sR <= 0) {
          fh = qR;
          fq = qR * uR + halfG * hR * hR;
        } else {
          const inv = 1 / (sR - sL);
          fh = (sR * qL - sL * qR + sL * sR * (hR - hL)) * inv;
          fq = (sR * (qL * uL + halfG * hL * hL) - sL * (qR * uR + halfG * hR * hR) + sL * sR * (qR - qL)) * inv;
        }
      }
      Fh[f] = fh;
      FqL[f] = fq + halfG * (hp[a] * hp[a] - hL * hL);
      FqR[f] = fq + halfG * (hm[f] * hm[f] - hR * hR);
    }
    const idx = 1 / dx;
    for (let i = iL; i <= iR; i++) {
      dh[i] = -(Fh[i + 1] - Fh[i]) * idx;
      dq[i] = -(FqL[i + 1] - FqR[i] + halfG * (hp[i] + hm[i]) * (zp[i] - zm[i])) * idx;
    }
  };

  let clipped = 0;
  const stage = (hIn, qIn, hOut, qOut, dt) => {
    for (let i = iL; i <= iR; i++) {
      let hn = hIn[i] + dt * dh[i];
      let qn = qIn[i] + dt * dq[i];
      if (hn < 0) {
        clipped += -hn * dx;
        hn = 0;
      }
      if (hn <= hDry) qn = 0;
      else if (hasFriction && nMan[i] > 0) {
        qn /= 1 + (dt * g * nMan[i] * nMan[i] * Math.abs(qn)) / Math.pow(hn, 7 / 3);
      }
      hOut[i] = hn;
      qOut[i] = qn;
    }
  };

  const x = new Float64Array(n);
  for (let i = 0; i < n; i++) x[i] = x0 + (i + 0.5) * dx;
  const hMax = new Float64Array(n);
  const etaMax = new Float64Array(n).fill(Number.NaN);
  const uMax = new Float64Array(n);
  const shT = [];
  const shX = [];
  const shZ = [];
  let inundationLimit = Number.NaN;
  let maxRunup = Number.NEGATIVE_INFINITY;
  let maxRunupTime = Number.NaN;
  const record = (t) => {
    for (let i = 0; i < n; i++) {
      const hi = h[i + 2];
      if (hi > hMax[i]) hMax[i] = hi;
      if (hi > hShore) {
        const e = hi + z[i + 2];
        if (!(etaMax[i] >= e)) etaMax[i] = e;
        const ua = Math.abs(q[i + 2] / hi);
        if (ua > uMax[i]) uMax[i] = ua;
      }
    }
    if (landSide === 'right' || landSide === 'left') {
      // Walk from the sea side to the first dry cell: the wet front connected to the sea.
      let last = -1;
      if (landSide === 'right') {
        for (let i = 0; i < n && h[i + 2] > hShore; i++) last = i;
      } else {
        for (let i = n - 1; i >= 0 && h[i + 2] > hShore; i--) last = i;
      }
      if (last >= 0) {
        const xs = landSide === 'right' ? x[last] + 0.5 * dx : x[last] - 0.5 * dx;
        const zs = h[last + 2] + z[last + 2];
        shT.push(t);
        shX.push(xs);
        shZ.push(zs);
        if (zs > maxRunup) {
          maxRunup = zs;
          maxRunupTime = t;
        }
        if (!(landSide === 'right' ? inundationLimit >= xs : inundationLimit <= xs)) inundationLimit = xs;
      }
    }
  };
  const totalVolume = () => {
    let s = 0;
    for (let i = iL; i <= iR; i++) s += h[i];
    return s * dx;
  };
  const snapshots = [];
  const snap = (t) => {
    const hs = new Float64Array(n);
    const us = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      hs[i] = h[i + 2];
      us[i] = h[i + 2] > hDry ? q[i + 2] / h[i + 2] : 0;
    }
    snapshots.push({ t, h: hs, u: us });
  };

  const outs = Array.from(outputTimes).filter((v) => v >= 0 && v <= tEnd).sort((a, b) => a - b);
  let io = 0;
  const vol0 = totalVolume();
  let boundaryNet = 0;
  let t = 0;
  let steps = 0;
  let dtMin = Number.POSITIVE_INFINITY;
  let dtMax = 0;
  record(0);
  while (io < outs.length && outs[io] <= 0) {
    snap(0);
    io++;
  }
  while (t < tEnd) {
    need(steps < maxSteps, 'maxSteps exceeded');
    applyBc(h, q, t);
    let smax = 0;
    for (let i = 1; i <= iR + 1; i++) {
      if (h[i] > hDry) {
        const s = Math.abs(q[i] / h[i]) + Math.sqrt(g * h[i]);
        if (s > smax) smax = s;
      }
    }
    let target = tEnd;
    if (io < outs.length && outs[io] < target) target = outs[io];
    let dt = smax > 0 ? (cfl * dx) / smax : target - t;
    let hit = false;
    if (t + dt >= target - 1e-12 * Math.max(target, 1)) {
      dt = target - t;
      hit = true;
    }
    need(Number.isFinite(dt) && dt > 0, 'time step collapsed (non-finite state)');
    rhs(h, q, t);
    let fIn = Fh[iL] - Fh[iR + 1];
    stage(h, q, h1, q1, dt);
    rhs(h1, q1, t + dt);
    fIn += Fh[iL] - Fh[iR + 1];
    stage(h1, q1, h1, q1, dt);
    for (let i = iL; i <= iR; i++) {
      h[i] = 0.5 * (h[i] + h1[i]);
      q[i] = h[i] > hDry ? 0.5 * (q[i] + q1[i]) : 0;
    }
    boundaryNet += 0.5 * dt * fIn;
    t = hit ? target : t + dt;
    steps++;
    if (dt < dtMin) dtMin = dt;
    if (dt > dtMax) dtMax = dt;
    record(t);
    while (io < outs.length && outs[io] <= t) {
      snap(t);
      io++;
    }
  }
  // Clipping is counted once per stage, but each stage contributes with weight ½ to the final update.
  clipped *= 0.5;
  const hOut = new Float64Array(n);
  const uOut = new Float64Array(n);
  const eOut = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    hOut[i] = h[i + 2];
    uOut[i] = h[i + 2] > hDry ? q[i + 2] / h[i + 2] : 0;
    eOut[i] = h[i + 2] + z[i + 2];
  }
  const volF = totalVolume();
  const error = volF - vol0 - boundaryNet - clipped;
  return {
    x, zb: Float64Array.from(zb), t, h: hOut, u: uOut, eta: eOut, snapshots,
    shoreline: { t: Float64Array.from(shT), x: Float64Array.from(shX), z: Float64Array.from(shZ), maxRunup, maxRunupTime },
    envelope: { hMax, etaMax, uMax, inundationLimit },
    volume: { initial: vol0, final: volF, boundaryNet, clipped, error, relativeError: vol0 > 0 ? error / vol0 : error },
    steps, dtMin, dtMax,
  };
}

/**
 * Analytic dam-break solutions on a horizontal frictionless bed with the dam at x = 0, released at t = 0.
 * Wet bed (Stoker 1957): left rarefaction, constant state (hm, um) and a right-moving bore, with hm from
 *   2 (√(g hL) − √(g hm)) = (hm − hR) √[ g (hm + hR) / (2 hm hR) ].
 * Dry bed (Ritter 1892): h = (2 √(g hL) − x/t)² / (9 g), u = (2/3)(√(g hL) + x/t) for −√(g hL) t <= x <= 2 √(g hL) t.
 * Reference: Stoker (1957) "Water Waves" §10.8; Ritter (1892) Z. Ver. Dtsch. Ing. 36.
 * @param {object} o
 * @param {number} o.hL upstream depth [m]
 * @param {number} [o.hR] downstream depth [m], default 0 (dry bed)
 * @param {number} o.x position relative to the dam [m]
 * @param {number} o.t time after release [s] (> 0)
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{h:number, u:number, hm:number, um:number, shockSpeed:number}} depth [m], velocity [m/s],
 *   middle state [m], [m/s] and bore speed [m/s] (front speed 2√(g hL) on a dry bed)
 */
export function damBreakAnalytic({ hL, hR = 0, x, t, g = G }) {
  need(Number.isFinite(hL) && hL > 0, 'hL must be positive');
  need(Number.isFinite(hR) && hR >= 0 && hR < hL, 'hR must lie in [0, hL)');
  need(Number.isFinite(t) && t > 0, 't must be positive');
  need(Number.isFinite(x), 'x must be finite');
  const cL = Math.sqrt(g * hL);
  const xi = x / t;
  if (hR === 0) {
    if (xi <= -cL) return { h: hL, u: 0, hm: 0, um: 2 * cL, shockSpeed: 2 * cL };
    if (xi >= 2 * cL) return { h: 0, u: 0, hm: 0, um: 2 * cL, shockSpeed: 2 * cL };
    return { h: (2 * cL - xi) ** 2 / (9 * g), u: (2 / 3) * (cL + xi), hm: 0, um: 2 * cL, shockSpeed: 2 * cL };
  }
  const f = (hm) => 2 * (cL - Math.sqrt(g * hm)) - (hm - hR) * Math.sqrt((g * (hm + hR)) / (2 * hm * hR));
  let lo = hR;
  let hi = hL;
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (f(mid) > 0) lo = mid;
    else hi = mid;
  }
  const hm = 0.5 * (lo + hi);
  const cm = Math.sqrt(g * hm);
  const um = 2 * (cL - cm);
  const S = (hm * um) / (hm - hR);
  let h;
  let u;
  if (xi <= -cL) {
    h = hL;
    u = 0;
  } else if (xi <= um - cm) {
    h = (2 * cL - xi) ** 2 / (9 * g);
    u = (2 / 3) * (cL + xi);
  } else if (xi <= S) {
    h = hm;
    u = um;
  } else {
    h = hR;
    u = 0;
  }
  return { h, u, hm, um, shockSpeed: S };
}
