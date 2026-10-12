/**
 * Blast response of structural members by the equivalent single-degree-of-freedom method.
 *
 * Reference: Biggs, "Introduction to Structural Dynamics" (1964) ch. 5 (also UFC 3-340-02 ch. 3 and FABIG
 * Technical Note 4). The distributed member is replaced by a mass–spring system whose displacement is the
 * deflection at the point of maximum deflection:
 *     K_LM m ü + R(u) = F(t),   K_LM = K_M / K_L,
 * with the load factor K_L = ∫ p φ dx / ∫ p dx and the mass factor K_M = ∫ m φ² dx / ∫ m dx evaluated for an assumed
 * deflected shape φ(x) normalised to 1 at the reference point. The shapes used here are the static deflected
 * shape under the same load distribution (elastic range) and the rigid-plastic hinge mechanism (plastic range);
 * the factors are integrated numerically from those shapes rather than tabulated. Dynamic support reactions
 * follow from the dynamic equilibrium of the member with the inertia forces distributed as φ(x).
 *
 * Resistance: elastic–perfectly-plastic with elastic unloading. SI units (m, kg, s, N, Pa); rotations in radians
 * unless the name ends in `Deg`.
 */
import { gaussLegendre } from '../math/quad.js';
import { brent } from '../math/roots.js';

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/* ------------------------------------------------------------------ deflected shapes and transformation factors */

/*
 * Assumed deflected shapes φ(ξ), ξ = x/L, normalised to 1 at the reference point (midspan, or the tip of a
 * cantilever whose fixed end is at ξ = 0). Elastic shapes are the Euler–Bernoulli static deflection curves:
 *   simple,     uniform:  w ∝ ξ − 2ξ³ + ξ⁴            (5/16 at midspan)
 *   simple,     point:    w ∝ 3ξ − 4ξ³                (ξ ≤ ½, symmetric)
 *   fixed,      uniform:  w ∝ ξ²(1 − ξ)²              (1/16 at midspan)
 *   fixed,      point:    w ∝ ξ²(3 − 4ξ)              (ξ ≤ ½, symmetric; ¼ at midspan)
 *   cantilever, uniform:  w ∝ ξ⁴ − 4ξ³ + 6ξ²          (3 at the tip)
 *   cantilever, point:    w ∝ 3ξ² − ξ³                (2 at the tip)
 * Plastic shapes are rigid links between hinges: 1 − |2ξ − 1| (hinge at midspan, and at the supports of a fixed
 * beam) and ξ (hinge at the root of a cantilever).
 */
const sym = (f) => (x) => f(x <= 0.5 ? x : 1 - x);
const SHAPES = {
  simple: {
    uniform: (x) => 16 / 5 * (x - 2 * x ** 3 + x ** 4),
    point: sym((x) => 3 * x - 4 * x ** 3),
    plastic: (x) => 1 - Math.abs(2 * x - 1),
  },
  fixed: {
    uniform: (x) => 16 * x * x * (1 - x) * (1 - x),
    point: sym((x) => 4 * x * x * (3 - 4 * x)),
    plastic: (x) => 1 - Math.abs(2 * x - 1),
  },
  cantilever: {
    uniform: (x) => (x ** 4 - 4 * x ** 3 + 6 * x * x) / 3,
    point: (x) => (3 * x * x - x ** 3) / 2,
    plastic: (x) => x,
  },
};

function checkCase(support, load, range) {
  if (!SHAPES[support]) throw new RangeError('support must be "simple", "fixed" or "cantilever"');
  if (load !== 'uniform' && load !== 'point') throw new RangeError('load must be "uniform" or "point"');
  if (range !== undefined && range !== 'elastic' && range !== 'plastic') throw new RangeError('range must be "elastic" or "plastic"');
}

/**
 * Load, mass and load–mass transformation factors of a prismatic beam with uniform mass, and the coefficients of
 * the dynamic support reaction V = a·R + b·F, derived from the assumed deflected shape (Biggs 1964 §5.3–5.5).
 *   K_L = ∫φ dξ (uniform load) or φ(load point) = 1 (point load at midspan / at the tip of a cantilever);
 *   K_M = ∫φ² dξ; K_LM = K_M/K_L.
 * Reactions: moments about the centroid x̄ of the inertia forces (distributed as φ) of the half span between a support
 * and midspan (whole span for a cantilever): V x̄ = ΣM + F' (x̄ − x_F), with ΣM the sum of the bending moments at the ends
 * of that free body expressed through the resistance (ΣM = R L/8 uniform, R L/4 point; cantilever R L/2, R L), F' the load
 * on the free body (F/2, or F for the cantilever) and x_F its line of action.
 * Published values for comparison (Biggs Tables 5.1–5.3): simple/uniform elastic K_L 0.64, K_M 0.50, V = 0.39R + 0.11F;
 * plastic 0.50, 0.33, 0.38R + 0.12F; fixed/uniform elastic 0.53, 0.41, 0.36R + 0.14F.
 * @param {object} o
 * @param {'simple'|'fixed'|'cantilever'} o.support end conditions (both ends, or fixed–free)
 * @param {'uniform'|'point'} [o.load='uniform'] uniformly distributed load, or a point load at midspan (at the tip for a cantilever)
 * @param {'elastic'|'plastic'} [o.range='elastic'] response range
 * @returns {{KL:number, KM:number, KLM:number, reaction:{a:number, b:number}, inertiaCentroid:number, shape:function(number):number}}
 *   factors [-]; reaction coefficients of V = a·R + b·F at each support [-]; x̄/L [-]; φ(ξ)
 */
export function beamTransformationFactors({ support, load = 'uniform', range = 'elastic' }) {
  checkCase(support, load, range);
  const phi = range === 'plastic' ? SHAPES[support].plastic : SHAPES[support][load];
  const g = gaussLegendre(16);
  // integrate over [0, ½] and [½, 1] separately: the plastic and point-load shapes have a kink at midspan
  let i1 = 0; let i2 = 0; let h0 = 0; let h1 = 0; let c0 = 0; let c1 = 0;
  for (const [a, b] of [[0, 0.5], [0.5, 1]]) {
    for (let q = 0; q < g.x.length; q++) {
      const x = 0.5 * (a + b) + 0.5 * (b - a) * g.x[q]; const w = 0.5 * (b - a) * g.w[q];
      const f = phi(x);
      i1 += w * f; i2 += w * f * f;
      c0 += w * f; c1 += w * f * x;
      if (a === 0) { h0 += w * f; h1 += w * f * x; }
    }
  }
  const KL = load === 'uniform' ? i1 : 1;
  const KM = i2;
  let a; let b; let xbar;
  if (support === 'cantilever') {
    xbar = c1 / c0;
    const sumM = load === 'uniform' ? 0.5 : 1; const xF = load === 'uniform' ? 0.5 : 1;
    a = sumM / xbar; b = 1 - xF / xbar;
  } else {
    xbar = h1 / h0;
    const sumM = load === 'uniform' ? 0.125 : 0.25; const xF = load === 'uniform' ? 0.25 : 0.5;
    a = sumM / xbar; b = 0.5 * (1 - xF / xbar);
  }
  return { KL, KM, KLM: KM / KL, reaction: { a, b }, inertiaCentroid: xbar, shape: phi };
}

/**
 * Equivalent SDOF properties of a prismatic beam for blast response: elastic stiffness, ultimate resistance from
 * the plastic moments, the equivalent bilinear (elastic–perfectly-plastic) resistance function and the
 * transformation factors.
 * Stiffness (Euler–Bernoulli, load F total): simple 384EI/(5L³) uniform, 48EI/L³ point; fixed 384EI/L³ uniform, 192EI/L³ point;
 * cantilever 8EI/L³ uniform, 3EI/L³ point. Ultimate resistance by the hinge mechanism: simple 8Mp/L, 4Mp/L; fixed
 * 8(Mps + Mpm)/L, 4(Mps + Mpm)/L; cantilever 2Mps/L, Mps/L. A fixed beam yields in two stages (first hinges where the elastic
 * moment first reaches the plastic moment, then the remaining simply supported or double-cantilever system); the
 * tri-linear curve is replaced by the bilinear one of equal area up to the displacement at which the mechanism forms
 * (Biggs §5.5; for Mps = Mpm under uniform load this gives the published k_E = 307 EI/L³).
 * @param {object} o
 * @param {'simple'|'fixed'|'cantilever'} o.support
 * @param {'uniform'|'point'} [o.load='uniform']
 * @param {number} o.L span [m]
 * @param {number} o.EI bending stiffness [N·m²]
 * @param {number} o.Mp plastic moment at midspan [N·m] (dynamic strength; for a cantilever the moment at the root)
 * @param {number} [o.MpSupport=o.Mp] plastic moment at the supports of a fixed beam [N·m]
 * @param {number} o.mass total mass of the beam including attached mass spread along it [kg]
 * @returns {{k:number, kElastic:number, Rm:number, Re:number, ue:number, uMechanism:number, KLM:{elastic:number, plastic:number},
 *   factors:{elastic:object, plastic:object}, period:number, stages:Array<{k:number, R:number, u:number}>}}
 *   k equivalent bilinear stiffness [N/m], first-stage stiffness [N/m], ultimate resistance Rm and first-yield resistance Re [N],
 *   equivalent yield displacement ue = Rm/k [m], displacement at mechanism formation [m], natural period 2π√(K_LM,e m/k) [s],
 *   stages of the multi-linear resistance (stiffness, resistance and displacement at the end of each stage)
 */
export function equivalentSdof({ support, load = 'uniform', L, EI, Mp, MpSupport, mass }) {
  checkCase(support, load);
  pos('L', L); pos('EI', EI); pos('Mp', Mp); pos('mass', mass);
  const Ms = MpSupport === undefined ? Mp : pos('MpSupport', MpSupport);
  const u = load === 'uniform';
  const c = EI / L ** 3;
  const stages = [];
  let Rm;
  if (support === 'simple') {
    Rm = (u ? 8 : 4) * Mp / L;
    stages.push({ k: (u ? 384 / 5 : 48) * c, R: Rm });
  } else if (support === 'cantilever') {
    Rm = (u ? 2 : 1) * Mp / L;
    stages.push({ k: (u ? 8 : 3) * c, R: Rm });
  } else {
    Rm = (u ? 8 : 4) * (Ms + Mp) / L;
    // elastic moments: uniform R L/12 at the supports and R L/24 at midspan; point R L/8 at both
    const rSup = (u ? 12 : 8) * Ms / L; const rMid = (u ? 24 : 8) * Mp / L;
    const Re = Math.min(rSup, rMid);
    stages.push({ k: (u ? 384 : 192) * c, R: Math.min(Re, Rm) });
    if (Re < Rm * (1 - 1e-12)) {
      // support hinges first → simply supported remainder; midspan hinge first → two cantilevers of length L/2
      const k2 = rSup <= rMid ? (u ? 384 / 5 : 48) * c : (u ? 128 : 48) * c;
      stages.push({ k: k2, R: Rm });
    }
  }
  let area = 0; let uPrev = 0; let rPrev = 0;
  for (const s of stages) {
    s.u = uPrev + (s.R - rPrev) / s.k;
    area += 0.5 * (s.R + rPrev) * (s.u - uPrev);
    uPrev = s.u; rPrev = s.R;
  }
  const uMechanism = uPrev;
  const ue = 2 * (uMechanism - area / Rm);
  const k = Rm / ue;
  const fe = beamTransformationFactors({ support, load, range: 'elastic' });
  const fp = beamTransformationFactors({ support, load, range: 'plastic' });
  return {
    k, kElastic: stages[0].k, Rm, Re: stages[0].R, ue, uMechanism, KLM: { elastic: fe.KLM, plastic: fp.KLM },
    factors: { elastic: fe, plastic: fp }, period: 2 * Math.PI * Math.sqrt(fe.KLM * mass / k), stages,
  };
}

/* ------------------------------------------------------------------ load histories */

/**
 * Build a force (or pressure) history function from a description.
 *   { type: 'triangular', peak, duration, riseTime = 0 }  linear rise to the peak then linear decay to zero at `duration`
 *   { type: 'rectangular', peak, duration }
 *   { type: 'table', t: [...], F: [...] }                  piecewise linear, zero outside the table
 *   a function of time is passed through (give `duration` separately where it is needed).
 * @param {object|function(number):number} load load description
 * @returns {{f:function(number):number, duration:number, peak:number, impulse:number}} F(t), end of loading [s], largest |F|, ∫F dt
 */
export function loadHistory(load) {
  if (typeof load === 'function') return { f: load, duration: NaN, peak: NaN, impulse: NaN };
  if (!load || typeof load !== 'object') throw new RangeError('load must be an object or a function');
  if (load.type === 'triangular') {
    const P = load.peak; const td = pos('load.duration', load.duration); const tr = load.riseTime === undefined ? 0 : nonNeg('load.riseTime', load.riseTime);
    if (!Number.isFinite(P)) throw new RangeError('load.peak must be finite');
    if (tr >= td) throw new RangeError('load.riseTime must be less than load.duration');
    return {
      f: (t) => (t < 0 || t >= td ? 0 : (t < tr ? P * t / tr : P * (td - t) / (td - tr))), duration: td, peak: Math.abs(P), impulse: 0.5 * P * td,
    };
  }
  if (load.type === 'rectangular') {
    const P = load.peak; const td = pos('load.duration', load.duration);
    if (!Number.isFinite(P)) throw new RangeError('load.peak must be finite');
    return { f: (t) => (t < 0 || t >= td ? 0 : P), duration: td, peak: Math.abs(P), impulse: P * td };
  }
  if (load.type === 'table') {
    const ts = load.t; const fs = load.F;
    if (!ts || !fs || ts.length !== fs.length || ts.length < 2) throw new RangeError('load.t and load.F must have the same length (>= 2)');
    for (let i = 1; i < ts.length; i++) if (!(ts[i] > ts[i - 1])) throw new RangeError('load.t must be strictly ascending');
    let peak = 0; let imp = 0;
    for (let i = 0; i < ts.length; i++) { peak = Math.max(peak, Math.abs(fs[i])); if (i > 0) imp += 0.5 * (fs[i] + fs[i - 1]) * (ts[i] - ts[i - 1]); }
    let j = 0;
    const f = (t) => {
      if (t < ts[0] || t > ts[ts.length - 1]) return 0;
      if (t < ts[j]) j = 0;
      while (j < ts.length - 2 && t > ts[j + 1]) j++;
      return fs[j] + (fs[j + 1] - fs[j]) * (t - ts[j]) / (ts[j + 1] - ts[j]);
    };
    return { f, duration: ts[ts.length - 1], peak, impulse: imp };
  }
  throw new RangeError('load.type must be "triangular", "rectangular" or "table"');
}

/* ------------------------------------------------------------------ SDOF integration */

/*
 * Velocity-Verlet integration of  me(state) ü + c u̇ + R(u) = F(t)  with an elastic–perfectly-plastic spring.
 * The plastic offset `up` carries the history: R = k (u − up) limited to [−Rneg, Rm]. `onStep` receives each state.
 * The external load is sampled at the middle of each step (midpoint rule for its impulse), so that a jump in F(t)
 * falling on a step boundary — the end of a rectangular pulse — is integrated exactly.
 * Returns the extreme values; stops at tEnd, or (stopAtPeak) at the first velocity reversal after the load has ended
 * or once the load can no longer increase the response.
 */
function integrate({ m, k, Rm, Rneg, klmE, klmP, c, f, u0, v0, dtOf, tEnd, stopAtPeak, tLoadEnd, onStep }) {
  let u = u0; let v = v0; let up = 0; let t = 0;
  let R = k * u; let yielding = false;
  const resist = () => {
    R = k * (u - up); yielding = false;
    if (R > Rm) { up = u - Rm / k; R = Rm; yielding = true; } else if (R < -Rneg) { up = u + Rneg / k; R = -Rneg; yielding = true; }
  };
  resist();
  let me = (yielding ? klmP : klmE) * m;
  let F = f(0);
  let uMax = u; let tMax = 0; let uMin = u; let tMin = 0; let rMax = Math.abs(R); let steps = 0;
  if (onStep) onStep(t, u, v, R, F);
  while (t < tEnd) {
    const dt = Math.min(dtOf(t), tEnd - t);
    if (!(dt > 0)) break;
    const Fm = f(t + 0.5 * dt);
    const vh = v + 0.5 * dt * (Fm - R - c * v) / me;
    u += dt * vh;
    t += dt;
    resist();
    me = (yielding ? klmP : klmE) * m;
    F = f(t);
    const vPrev = v;
    v = (vh + 0.5 * dt * (Fm - R) / me) / (1 + 0.5 * dt * c / me);
    steps++;
    if (u > uMax) { uMax = u; tMax = t; }
    if (u < uMin) { uMin = u; tMin = t; }
    if (Math.abs(R) > rMax) rMax = Math.abs(R);
    if (onStep) onStep(t, u, v, R, F);
    if (stopAtPeak && vPrev > 0 && v <= 0 && t >= tLoadEnd) break;
  }
  return { uMax, tMax, uMin, tMin, up, rMax, steps, tLast: t, vLast: v };
}

/**
 * Response of an equivalent SDOF system with elastic–perfectly-plastic resistance to a load history:
 * K_LM m ü + c u̇ + R(u) = F(t), R = k(u − u_p) limited to ±Rm, elastic unloading. Explicit velocity-Verlet integration
 * (second order, time step ≤ T/400 and ≤ load duration/400 by default). Reference: Biggs (1964) §2.7, §5.
 * @param {object} o
 * @param {number} o.mass actual mass [kg]
 * @param {number} o.k stiffness [N/m]
 * @param {number} [o.Rm=Infinity] ultimate resistance [N] (Infinity = elastic)
 * @param {number} [o.RmRebound=o.Rm] resistance in the rebound direction [N]
 * @param {number|{elastic:number, plastic:number}} [o.klm=1] load–mass factor [-] (the plastic value applies while the spring is yielding)
 * @param {number} [o.zeta=0] damping ratio of the elastic system [-] (blast response is normally computed undamped)
 * @param {object|function(number):number} o.load force history: see loadHistory [N]
 * @param {number} [o.dt] time step [s]
 * @param {number} [o.tEnd] end of the analysis [s] (default: load duration + 3 natural periods, extended until the first peak has passed)
 * @param {number} [o.u0=0] initial displacement [m]
 * @param {number} [o.v0=0] initial velocity [m/s] (an ideal impulse I is applied as v0 = I/(K_LM m))
 * @param {{a:number, b:number}} [o.reaction] coefficients of the dynamic reaction V = a·R + b·F
 * @returns {{t:Float64Array, u:Float64Array, v:Float64Array, R:Float64Array, F:Float64Array, V:(Float64Array|undefined), uMax:number, tMax:number,
 *   uMin:number, ue:number, ductility:number, plasticDisplacement:number, dlf:number, period:number, omega:number, vMax:(number|undefined),
 *   peakLoad:number}} histories (time [s], displacement [m], velocity [m/s], resistance [N], load [N], reaction [N]); maximum and minimum
 *   displacement [m], time of the maximum [s], yield displacement Rm/k [m], ductility uMax/ue [-], permanent offset [m], dynamic load factor
 *   uMax k / Fpeak [-], elastic period [s] and circular frequency [rad/s], largest |V| [N], largest |F| [N]
 */
export function sdofBlastResponse({ mass, k, Rm = Infinity, RmRebound, klm = 1, zeta = 0, load, dt, tEnd, u0 = 0, v0 = 0, reaction }) {
  pos('mass', mass); pos('k', k); nonNeg('zeta', zeta);
  if (!(Rm > 0)) throw new RangeError('Rm must be positive');
  const Rneg = RmRebound === undefined ? Rm : RmRebound;
  if (!(Rneg > 0)) throw new RangeError('RmRebound must be positive');
  const klmE = typeof klm === 'number' ? klm : klm.elastic; const klmP = typeof klm === 'number' ? klm : klm.plastic;
  pos('klm (elastic)', klmE); pos('klm (plastic)', klmP);
  const lh = loadHistory(load);
  const omega = Math.sqrt(k / (klmE * mass));
  const T = 2 * Math.PI / omega;
  const td = Number.isFinite(lh.duration) ? lh.duration : 0;
  // default step: at most T/400 and td/400, and a whole number of steps over the load duration
  const h = dt === undefined ? (td > 0 ? td / Math.ceil(td / Math.min(T / 400, td / 400) - 1e-9) : T / 400) : pos('dt', dt);
  const c = 2 * zeta * omega * klmE * mass;
  const ts = []; const us = []; const vs = []; const Rs = []; const Fs = [];
  let peakLoad = 0;
  const onStep = (t, u, v, R, F) => { ts.push(t); us.push(u); vs.push(v); Rs.push(R); Fs.push(F); if (Math.abs(F) > peakLoad) peakLoad = Math.abs(F); };
  const base = { m: mass, k, Rm, Rneg, klmE, klmP, c, f: lh.f, u0, v0, dtOf: () => h, onStep };
  let end = tEnd === undefined ? td + 3 * T : pos('tEnd', tEnd);
  let r = integrate({ ...base, tEnd: end });
  if (tEnd === undefined) {
    // a strongly plastic response can peak later than three elastic periods after the load: extend until the peak has passed
    let guard = 0;
    while (r.tMax >= end - 2 * h && guard++ < 8) {
      end *= 2; ts.length = 0; us.length = 0; vs.length = 0; Rs.length = 0; Fs.length = 0;
      r = integrate({ ...base, tEnd: end });
    }
  }
  const ue = Rm / k;
  const out = {
    t: Float64Array.from(ts), u: Float64Array.from(us), v: Float64Array.from(vs), R: Float64Array.from(Rs), F: Float64Array.from(Fs),
    uMax: r.uMax, tMax: r.tMax, uMin: r.uMin, ue, ductility: Number.isFinite(ue) ? r.uMax / ue : 0, plasticDisplacement: r.up,
    dlf: peakLoad > 0 ? r.uMax * k / peakLoad : NaN, period: T, omega, peakLoad,
  };
  if (reaction) {
    const V = new Float64Array(ts.length); let vMax = 0;
    for (let i = 0; i < ts.length; i++) { V[i] = reaction.a * Rs[i] + reaction.b * Fs[i]; if (Math.abs(V[i]) > vMax) vMax = Math.abs(V[i]); }
    out.V = V; out.vMax = vMax;
  }
  return out;
}

/* Peak displacement only (no stored history), with a time step refined during a short load. */
function peakOnly({ mass, k, Rm, klmE, klmP, f, td, T }) {
  const hLoad = Math.min(T / 300, td / 60); const hFree = T / 300;
  const r = integrate({
    m: mass, k, Rm, Rneg: Rm, klmE, klmP, c: 0, f, u0: 0, v0: 0, dtOf: (t) => (t < td ? Math.min(hLoad, td - t + 1e-15 * td) : hFree),
    tEnd: td + 400 * T, stopAtPeak: true, tLoadEnd: 0,
  });
  return r.uMax;
}

/**
 * Strain energy absorbed by the elastic–perfectly-plastic spring at displacement u: ½ k u² (u ≤ ue), Rm (u − ue/2) beyond.
 * @param {object} o
 * @param {number} o.k stiffness [N/m]
 * @param {number} [o.Rm=Infinity] ultimate resistance [N]
 * @param {number} o.u displacement [m]
 * @returns {number} energy [J]
 */
export function strainEnergyEpp({ k, Rm = Infinity, u }) {
  pos('k', k); nonNeg('u', u);
  const ue = Rm / k;
  return u <= ue ? 0.5 * k * u * u : Rm * (u - 0.5 * ue);
}

/**
 * Pressure–impulse (iso-damage) curve of an elastic–perfectly-plastic SDOF system: for each load duration the peak
 * load that produces exactly the target maximum displacement is found by root finding on the time integration.
 * Asymptotes from energy balance (Baker et al., Explosion Hazards and Evaluation, 1983, ch. 4): impulsive
 * I = √(2 K_LM m U(u)), quasi-static (suddenly applied, long duration) F = U(u)/u, U = strain energy at the target.
 * The load must not increase after its start (triangular with zero rise time, or rectangular) so that the first
 * displacement peak is the maximum.
 * @param {object} o
 * @param {number} o.mass mass [kg]
 * @param {number} o.k stiffness [N/m]
 * @param {number} [o.Rm=Infinity] ultimate resistance [N]
 * @param {number|{elastic:number, plastic:number}} [o.klm=1] load–mass factor [-]
 * @param {number} [o.ductility] target ductility μ = u/ue [-] (needs a finite Rm), or
 * @param {number} [o.uTarget] target maximum displacement [m]
 * @param {'triangular'|'rectangular'} [o.shape='triangular'] pulse shape
 * @param {number} [o.nPoints=30] number of durations, logarithmically spaced
 * @param {[number, number]} [o.durationRatio=[0.02, 50]] range of td/T [-]
 * @returns {{points:Array<{td:number, F:number, I:number}>, impulsiveAsymptote:number, quasiStaticAsymptote:number, uTarget:number, period:number}}
 *   durations [s], peak loads [N], impulses [N·s]; asymptotes I [N·s] and F [N]. Divide F and I by the loaded area for pressure units.
 */
export function pressureImpulseCurve({ mass, k, Rm = Infinity, klm = 1, ductility, uTarget, shape = 'triangular', nPoints = 30, durationRatio = [0.02, 50] }) {
  pos('mass', mass); pos('k', k);
  if (!(Rm > 0)) throw new RangeError('Rm must be positive');
  if (shape !== 'triangular' && shape !== 'rectangular') throw new RangeError('shape must be "triangular" or "rectangular"');
  if (!Number.isInteger(nPoints) || nPoints < 2) throw new RangeError('nPoints must be an integer >= 2');
  const klmE = typeof klm === 'number' ? klm : klm.elastic; const klmP = typeof klm === 'number' ? klm : klm.plastic;
  let ut;
  if (uTarget !== undefined) ut = pos('uTarget', uTarget);
  else {
    pos('ductility', ductility);
    if (!Number.isFinite(Rm)) throw new RangeError('ductility needs a finite Rm; use uTarget for an elastic system');
    ut = ductility * Rm / k;
  }
  const [r0, r1] = durationRatio;
  if (!(r0 > 0 && r1 > r0)) throw new RangeError('durationRatio must be [min, max] with 0 < min < max');
  const T = 2 * Math.PI * Math.sqrt(klmE * mass / k);
  const U = strainEnergyEpp({ k, Rm, u: ut });
  const yielding = ut > Rm / k;
  const iAsym = Math.sqrt(2 * (yielding ? klmP : klmE) * mass * U);
  const fAsym = U / ut;
  const points = [];
  for (let i = 0; i < nPoints; i++) {
    const td = T * r0 * (r1 / r0) ** (i / (nPoints - 1));
    const g = (lnF) => {
      const P = Math.exp(lnF);
      const f = shape === 'triangular' ? (t) => (t >= td ? 0 : P * (1 - t / td)) : (t) => (t >= td ? 0 : P);
      return peakOnly({ mass, k, Rm, klmE, klmP, f, td, T }) - ut;
    };
    // the response increases monotonically with the peak load; the quasi-static asymptote bounds it from below
    let lo = Math.log(0.5 * fAsym); let hi = Math.log(2 * fAsym);
    let guard = 0;
    while (g(hi) < 0) { lo = hi; hi += Math.log(4); if (++guard > 60) throw new Error('pressureImpulseCurve: could not bracket the peak load'); }
    const F = Math.exp(brent(g, lo, hi, { tol: 1e-10 }));
    points.push({ td, F, I: (shape === 'triangular' ? 0.5 : 1) * F * td });
  }
  return { points, impulsiveAsymptote: iAsym, quasiStaticAsymptote: fAsym, uTarget: ut, period: T };
}

/**
 * Support rotation of a member from its maximum deflection, the usual blast damage measure together with the
 * ductility (UFC 3-340-02; FABIG Technical Note 4): θ = atan(u/(L/2)) for members supported at both ends, atan(u/L) for
 * a cantilever. Acceptance limits are project criteria.
 * @param {object} o
 * @param {number} o.deflection maximum deflection [m]
 * @param {number} o.L span [m]
 * @param {'simple'|'fixed'|'cantilever'} [o.support='simple']
 * @returns {{rotation:number, rotationDeg:number}} support rotation [rad] and [deg]
 */
export function supportRotation({ deflection, L, support = 'simple' }) {
  nonNeg('deflection', deflection); pos('L', L);
  if (!SHAPES[support]) throw new RangeError('support must be "simple", "fixed" or "cantilever"');
  const th = Math.atan(deflection / (support === 'cantilever' ? L : L / 2));
  return { rotation: th, rotationDeg: th * 180 / Math.PI };
}

/**
 * Blast response of a prismatic beam: equivalent SDOF from equivalentSdof, time integration, ductility, support
 * rotation and dynamic reactions V = a·R + b·F (coefficients of the elastic shape while the response is elastic, of the
 * plastic mechanism once the resistance has reached Rm).
 * @param {object} o
 * @param {'simple'|'fixed'|'cantilever'} o.support
 * @param {'uniform'|'point'} [o.load='uniform'] load distribution
 * @param {number} o.L span [m]
 * @param {number} o.EI bending stiffness [N·m²]
 * @param {number} o.Mp dynamic plastic moment at midspan (root of a cantilever) [N·m]
 * @param {number} [o.MpSupport=o.Mp] dynamic plastic moment at the supports of a fixed beam [N·m]
 * @param {number} o.mass total mass [kg]
 * @param {object|function(number):number} o.history history of the total load on the beam [N]: see loadHistory (pressure × loaded width × L)
 * @param {number} [o.dt] time step [s]
 * @param {number} [o.tEnd] end time [s]
 * @returns {{sdof:object, response:object, ductility:number, rotationDeg:number, uMax:number, reactionMax:number, dlf:number}} equivalent system,
 *   sdofBlastResponse result, ductility [-], support rotation [deg], maximum deflection [m], largest dynamic reaction per support [N], DLF [-]
 */
export function beamBlastResponse({ support, load = 'uniform', L, EI, Mp, MpSupport, mass, history, dt, tEnd }) {
  const sd = equivalentSdof({ support, load, L, EI, Mp, MpSupport, mass });
  const res = sdofBlastResponse({ mass, k: sd.k, Rm: sd.Rm, klm: sd.KLM, load: history, dt, tEnd });
  const re = sd.factors.elastic.reaction; const rp = sd.factors.plastic.reaction;
  let reactionMax = 0;
  const V = new Float64Array(res.t.length);
  for (let i = 0; i < V.length; i++) {
    const c = Math.abs(res.R[i]) >= sd.Rm * (1 - 1e-9) ? rp : re;
    V[i] = c.a * res.R[i] + c.b * res.F[i];
    if (Math.abs(V[i]) > reactionMax) reactionMax = Math.abs(V[i]);
  }
  res.V = V; res.vMax = reactionMax;
  return {
    sdof: sd, response: res, ductility: res.ductility, rotationDeg: supportRotation({ deflection: Math.max(res.uMax, 0), L, support }).rotationDeg,
    uMax: res.uMax, reactionMax, dlf: res.dlf,
  };
}

/* ------------------------------------------------------------------ strain rate and drag */

/**
 * Dynamic increase factor of the yield (flow) stress by the Cowper–Symonds relation: DIF = 1 + (ε̇/D)^(1/q)
 * (Cowper & Symonds 1957, Brown University Tech. Rep. 28). The constants are material data and must be supplied;
 * values commonly quoted for the yield stress (Jones, Structural Impact, 1989, Table 8.1): mild steel D = 40.4 1/s, q = 5;
 * high-tensile steel D = 3200 1/s, q = 5; stainless steel 304 D = 100 1/s, q = 10; aluminium alloy D = 6500 1/s, q = 4.
 * @param {object} o
 * @param {number} o.strainRate ε̇ [1/s] (>= 0)
 * @param {number} o.D reference strain rate [1/s]
 * @param {number} o.q exponent [-]
 * @returns {number} DIF [-] (>= 1)
 */
export function dynamicIncreaseFactor({ strainRate, D, q }) {
  nonNeg('strainRate', strainRate); pos('D', D); pos('q', q);
  return 1 + (strainRate / D) ** (1 / q);
}

/**
 * Mean strain rate of a member reaching yield at time tE after the start of loading: ε̇ = fy/(E tE) (the
 * average-rate form used in UFC 3-340-02 ch. 5, ε̇ = fds/(Es tE)); use with dynamicIncreaseFactor and iterate on the dynamic yield strength if required.
 * @param {object} o
 * @param {number} o.fy (dynamic) yield strength [Pa]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} o.timeToYield time to reach yield [s]
 * @returns {number} strain rate [1/s]
 */
export function meanStrainRate({ fy, E, timeToYield }) {
  pos('fy', fy); pos('E', E); pos('timeToYield', timeToYield);
  return fy / (E * timeToYield);
}

/**
 * Drag load of the blast wind on pipework, equipment and open structure: F = Cd · ½ρu² · A, with the dynamic
 * pressure q = ½ρu² of the gas flow behind the front (Morison-type drag term; inertia and pressure-difference
 * terms are not included). Appropriate for obstacles small enough not to disturb the overpressure field; larger items are
 * loaded by the pressure difference across them. The size limit between the two regimes is a project criterion.
 * @param {object} o
 * @param {number} [o.dynamicPressure] q = ½ρu² [Pa]; or give rho and velocity
 * @param {number} [o.rho] gas density [kg/m³]
 * @param {number} [o.velocity] gas velocity [m/s]
 * @param {number} o.Cd drag coefficient [-] (shape and Reynolds-number dependent; ≈ 1.2 for a smooth circular cylinder in subcritical flow)
 * @param {number} [o.area] projected area normal to the flow [m²]; or give D and length
 * @param {number} [o.D] diameter / width [m]
 * @param {number} [o.length=1] length [m]
 * @returns {{force:number, forcePerLength:(number|undefined), dynamicPressure:number}} drag force [N], per unit length [N/m], q [Pa]
 */
export function blastDragLoad({ dynamicPressure, rho, velocity, Cd, area, D, length = 1 }) {
  pos('Cd', Cd);
  let q = dynamicPressure;
  if (q === undefined) { pos('rho', rho); if (!Number.isFinite(velocity)) throw new RangeError('velocity must be finite'); q = 0.5 * rho * velocity * velocity; } else nonNeg('dynamicPressure', q);
  if (area !== undefined) return { force: Cd * q * nonNeg('area', area), forcePerLength: undefined, dynamicPressure: q };
  pos('D', D); pos('length', length);
  return { force: Cd * q * D * length, forcePerLength: Cd * q * D, dynamicPressure: q };
}
