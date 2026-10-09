/**
 * Transient, cyclic and emergency thermal rating of cables, and a 2-D steady conduction solver for trench
 * cross-sections.
 *
 * References
 *  - IEC 60853-2:1989 Calculation of the cyclic and emergency current rating of cables — Part 2 (thermal ladder
 *    network, Van Wormer coefficient, two-loop equivalent network, exponential-integral soil response and
 *    conductor-to-surface attainment factor).
 *  - IEC 60853-1:1985 (cyclic rating factor M from the six hourly loss-load ordinates before the peak).
 *  - Van Wormer (1955), AIEE Trans. 74, "An improved approximate technique for calculating cable temperature
 *    transients".
 *  - Carslaw & Jaeger, Conduction of Heat in Solids, §10.3 (continuous line source).
 *  - IEC 60949 (adiabatic heating, used by the adiabatic check).
 *
 * The cable is a lumped ladder solved exactly by eigen-decomposition; the surroundings respond as a line source
 * with its image; the two are coupled with the attainment factor of IEC 60853-2 §4.2.3. Losses are superposed
 * linearly (resistance held at a stated temperature, with an optional temperature feedback pass).
 * SI units; temperatures in °C where the name ends in C, time in seconds.
 */

import { eigSym } from '../math/eigen.js';
import { metalProperties } from './conductors.js';

const TWO_PI = 2 * Math.PI;
const EULER = 0.5772156649015329;

function pos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a positive finite number`);
}

/**
 * Exponential integral E1(x) = ∫ₓ^∞ e^(−t)/t dt = −Ei(−x), x > 0.
 * Power series for x ≤ 1, modified Lentz continued fraction for x > 1 (Abramowitz & Stegun 5.1.11, 5.1.22).
 * @param {number} x argument (> 0)
 * @returns {number} E1(x)
 */
export function expIntegralE1(x) {
  if (!(x > 0)) throw new RangeError('x: must be positive');
  if (x > 700) return 0;
  if (x <= 1) {
    let sum = 0; let term = 1;
    for (let k = 1; k < 60; k++) {
      term *= -x / k;
      const d = -term / k;
      sum += d;
      if (Math.abs(d) < 1e-17 * Math.abs(sum)) break;
    }
    return -EULER - Math.log(x) + sum;
  }
  const tiny = 1e-300;
  let b = x + 1; let c = 1 / tiny; let d = 1 / b; let h = d;
  for (let i = 1; i < 200; i++) {
    const a = -i * i;
    b += 2;
    d = 1 / (a * d + b);
    c = b + a / c;
    const del = c * d;
    h *= del;
    if (Math.abs(del - 1) < 1e-16) break;
  }
  return h * Math.exp(-x);
}

/**
 * Van Wormer coefficient: fraction of the insulation (or covering) thermal capacitance placed at the inner node.
 * Long durations: p = 1/(2·ln(Di/dc)) − 1/((Di/dc)² − 1); short durations (< about one third of the cable time
 * constant): p* = 1/ln(Di/dc) − 1/(Di/dc − 1).
 * Reference: IEC 60853-2 §4.2.2; Van Wormer (1955).
 * @param {object} o
 * @param {number} o.outerDiameter outer diameter of the layer, m
 * @param {number} o.innerDiameter inner diameter of the layer, m
 * @param {boolean} [o.shortDuration=false]
 * @returns {number} p (0–0.5)
 */
export function vanWormer({ outerDiameter, innerDiameter, shortDuration = false }) {
  pos(outerDiameter, 'outerDiameter'); pos(innerDiameter, 'innerDiameter');
  if (outerDiameter <= innerDiameter) throw new RangeError('outerDiameter: must exceed innerDiameter');
  const r = outerDiameter / innerDiameter;
  return shortDuration ? 1 / Math.log(r) - 1 / (r - 1) : 1 / (2 * Math.log(r)) - 1 / (r * r - 1);
}

/**
 * Thermal ladder network of a cable (per core) from the geometry returned by ampacity.cableGeometry.
 * Nodes: conductor (Qc + p·Qi); sheath ((1 − p)·Qi + Qs + Qjacket + p′·Qbedding/n); armour ((1 − p′)·Qbedding/n +
 * Qarmour/n + p″·Qserving/n). Resistances between nodes: T1, n·T2, n·T3 (per-core equivalent: the stages beyond
 * the cores carry the heat of n cores). The part (1 − p″) of the outer covering capacitance is at the surface
 * node, whose temperature is imposed by the surroundings, and so does not appear (IEC 60853-2 practice).
 * Stages with zero resistance are merged with the next node.
 * Reference: IEC 60853-2 §4.2.
 * @param {object} g geometry from cableGeometry
 * @param {object} [o]
 * @param {boolean} [o.shortDuration=false] use the short-duration Van Wormer coefficient for the insulation
 * @returns {{Q:number[], T:number[], roles:string[][], n:number}} node capacitances Q (J/(K·m)), resistance T[i]
 *          from node i towards the surface (K·m/W), roles of each node, cores n
 */
export function cableLadder(g, { shortDuration = false } = {}) {
  if (!g || !g.layers) throw new RangeError('g: geometry object from cableGeometry required');
  const n = g.cores;
  const cap = (name) => { const l = g.layers.find((x) => x.name === name); return l ? l.heatCapacity * l.area : 0; };
  const lay = (name) => g.layers.find((x) => x.name === name);
  const ins = lay('insulation');
  const p = vanWormer({ outerDiameter: ins.outerDiameter, innerDiameter: ins.innerDiameter, shortDuration });
  const Qi = cap('insulation');
  const vw = (l) => (l ? vanWormer({ outerDiameter: l.outerDiameter, innerDiameter: l.innerDiameter }) : 0);
  const hasArmour = !!lay('armour');
  const bed = lay('bedding'); const serv = lay('serving');
  const nodes = [{ Q: cap('conductor') + p * Qi, T: g.T1, roles: ['conductor'] }];
  if (hasArmour) {
    const pb = vw(bed); const ps = vw(serv);
    nodes.push({ Q: (1 - p) * Qi + cap('sheath') + cap('coreJacket') + (pb * cap('bedding')) / n, T: n * g.T2, roles: ['sheath'] });
    nodes.push({ Q: ((1 - pb) * cap('bedding') + cap('armour') + ps * cap('serving')) / n, T: n * g.T3, roles: ['armour'] });
  } else {
    // Covering(s) in series form T3; their capacitance is split with the Van Wormer coefficient of the whole.
    const first = lay('coreJacket') || bed || serv;
    let pc = 0; let Qcov = 0;
    if (first) {
      const last = serv || bed || first;
      pc = vanWormer({ outerDiameter: last.outerDiameter, innerDiameter: first.innerDiameter });
      Qcov = cap('coreJacket') + (cap('bedding') + cap('serving')) / n;
    }
    nodes.push({ Q: (1 - p) * Qi + cap('sheath') + pc * Qcov, T: n * g.T3, roles: ['sheath', 'armour'] });
  }
  // Merge nodes separated by zero resistance.
  const merged = [];
  for (const nd of nodes) {
    const prev = merged[merged.length - 1];
    if (prev && !(prev.T > 0)) { prev.Q += nd.Q; prev.T = nd.T; prev.roles.push(...nd.roles); } else merged.push({ ...nd, roles: [...nd.roles] });
  }
  while (merged.length > 1 && !(merged[merged.length - 1].T > 0)) {
    const last = merged.pop();
    merged[merged.length - 1].Q += last.Q; merged[merged.length - 1].roles.push(...last.roles);
  }
  return { Q: merged.map((x) => x.Q), T: merged.map((x) => x.T), roles: merged.map((x) => x.roles), n };
}

/**
 * Exact step response of a thermal ladder: node i has capacitance Q[i] and is connected to node i+1 (or to the
 * reference after the last node) through T[i]; heat w[i] (W/m) is injected at node i from t = 0.
 * Solved by symmetrised eigen-decomposition: θ_i(t) = Σ_k a_ik·(1 − e^(−t/τk)).
 * @param {object} o
 * @param {number[]} o.Q capacitances, J/(K·m) (> 0)
 * @param {number[]} o.T resistances, K·m/W (> 0)
 * @param {number[]} o.w injected heat per node, W/m
 * @returns {{timeConstants:number[], amplitudes:number[][], steady:number[], rise:(t:number, node?:number)=>number}}
 *          τk (s, ascending), a_ik (K), steady rises (K) and evaluator of the rise of a node at time t (s)
 */
export function ladderStepResponse({ Q, T, w }) {
  const m = Q.length;
  if (m < 1 || T.length !== m || w.length !== m) throw new RangeError('Q, T, w: must be arrays of equal length ≥ 1');
  for (let i = 0; i < m; i++) { pos(Q[i], `Q[${i}]`); pos(T[i], `T[${i}]`); }
  // Conductance matrix G (symmetric tridiagonal), then S = D^(-1/2)·G·D^(-1/2) with D = diag(Q).
  const G = Array.from({ length: m }, () => new Array(m).fill(0));
  for (let i = 0; i < m; i++) {
    const gi = 1 / T[i];
    G[i][i] += gi;
    if (i + 1 < m) { G[i + 1][i + 1] += gi; G[i][i + 1] -= gi; G[i + 1][i] -= gi; }
  }
  const sq = Q.map(Math.sqrt);
  const S = G.map((row, i) => row.map((v, j) => v / (sq[i] * sq[j])));
  let values; let vectors;
  if (m === 1) { values = [S[0][0]]; vectors = [[1]]; } else ({ values, vectors } = eigSym(S));
  // θ(t) = Σ_k v_k·(v_kᵀ·D^(-1/2)·w)/λ_k·(1 − e^(−λ_k t)) in the scaled variables y = D^(1/2)·θ.
  const amplitudes = Array.from({ length: m }, () => new Array(m).fill(0));
  for (let k = 0; k < m; k++) {
    let proj = 0;
    for (let j = 0; j < m; j++) proj += (vectors[k][j] * w[j]) / sq[j];
    for (let i = 0; i < m; i++) amplitudes[i][k] = ((vectors[k][i] / sq[i]) * proj) / values[k];
  }
  const order = values.map((v, k) => k).sort((a, b) => values[b] - values[a]);
  const timeConstants = order.map((k) => 1 / values[k]);
  const amps = amplitudes.map((row) => order.map((k) => row[k]));
  const steady = amps.map((row) => row.reduce((s, a) => s + a, 0));
  const rise = (t, node = 0) => {
    if (t <= 0) return 0;
    let s = 0;
    for (let k = 0; k < m; k++) s += amps[node][k] * (1 - Math.exp(-t / timeConstants[k]));
    return s;
  };
  return { timeConstants, amplitudes: amps, steady, rise };
}

/**
 * Two-loop equivalent network response (closed form).
 * M0 = 0.5·(QA·(TA + TB) + QB·TB), N0 = QA·TA·QB·TB, a, b = (M0 ± sqrt(M0² − N0))/N0,
 * Ta = (1/QA − b·(TA + TB))/(a − b), Tb = TA + TB − Ta, θc(t) = Wc·[Ta·(1 − e^(−a·t)) + Tb·(1 − e^(−b·t))].
 * Reference: IEC 60853-2 §4.2.3 (response of the conductor to a step of conductor loss Wc).
 * @param {object} o
 * @param {number} o.TA @param {number} o.TB resistances of the first and second loop, K·m/W
 * @param {number} o.QA @param {number} o.QB capacitances of the first and second loop, J/(K·m)
 * @returns {{a:number, b:number, Ta:number, Tb:number, rise:(t:number)=>number}} a, b in 1/s; Ta, Tb in K·m/W;
 *          rise(t) is the conductor rise per unit conductor loss, K/(W/m)
 */
export function twoLoopResponse({ TA, TB, QA, QB }) {
  pos(TA, 'TA'); pos(TB, 'TB'); pos(QA, 'QA'); pos(QB, 'QB');
  const M0 = 0.5 * (QA * (TA + TB) + QB * TB);
  const N0 = QA * TA * QB * TB;
  const root = Math.sqrt(M0 * M0 - N0);
  const a = (M0 + root) / N0;
  const b = (M0 - root) / N0;
  const Ta = (1 / QA - b * (TA + TB)) / (a - b);
  const Tb = TA + TB - Ta;
  return { a, b, Ta, Tb, rise: (t) => (t <= 0 ? 0 : Ta * (1 - Math.exp(-a * t)) + Tb * (1 - Math.exp(-b * t))) };
}

/**
 * Temperature rise of the surface of a buried cable after a step of heat W (line source and its image):
 * θe(t) = ρT·W/(4π)·[E1(De²/(16·δ·t)) − E1(L²/(δ·t))] + Σ_k ρT·Wk/(4π)·[E1(dpk²/(4δt)) − E1(d′pk²/(4δt))].
 * For t → ∞ this tends to ρT·W/(2π)·ln(4L/De), i.e. the steady T4 with ln(2u) in place of ln(u + sqrt(u² − 1)).
 * Reference: IEC 60853-2 §4.2.4 (single cable) and §4.2.5 (groups); Carslaw & Jaeger §10.3.
 * @param {object} o
 * @param {number} o.t time after the step, s
 * @param {number} o.rhoSoil soil thermal resistivity, K·m/W
 * @param {number} o.diffusivity soil thermal diffusivity δ, m²/s (typically 0.3e-6 … 0.8e-6; 0.5e-6 for moist soil)
 * @param {number} o.depth depth of the cable axis L, m
 * @param {number} o.De cable external diameter, m
 * @param {{d:number, dImage:number, weight?:number}[]} [o.others] other equally loaded cables: distance from the
 *        rated cable to cable k and to its image (m); weight = Wk/W (default 1)
 * @returns {number} rise per unit heat, K/(W/m)
 */
export function soilStepResponse({ t, rhoSoil, diffusivity, depth, De, others = [] }) {
  pos(rhoSoil, 'rhoSoil'); pos(diffusivity, 'diffusivity'); pos(depth, 'depth'); pos(De, 'De');
  if (!(t > 0)) return 0;
  const dt = diffusivity * t;
  let s = expIntegralE1((De * De) / (16 * dt)) - expIntegralE1((depth * depth) / dt);
  for (const o of others) {
    s += (o.weight ?? 1) * (expIntegralE1((o.d * o.d) / (4 * dt)) - expIntegralE1((o.dImage * o.dImage) / (4 * dt)));
  }
  return (rhoSoil / (4 * Math.PI)) * s;
}

/**
 * Complete transient model of a cable in its surroundings.
 * Conductor rise after a step of conductor loss Wc (per core) is
 *   θ(t) = θcable(t) + α(t)·Wtot·Rext(t),  α(t) = θcable(t)/θcable(∞),  Wtot = n·Wc·(1 + λ1 + λ2),
 * where θcable is the exact ladder response (sheath and armour losses injected at their nodes) and Rext(t) is the
 * soil step response (buried) or a constant external resistance (water, air: no external heat storage).
 * Reference: IEC 60853-2 §4.2–4.4.
 * @param {object} o
 * @param {object} o.geometry geometry from ampacity.cableGeometry
 * @param {number} [o.lambda1=0] @param {number} [o.lambda2=0] sheath and armour loss factors
 * @param {object} o.external { type: 'buried', rhoSoil, diffusivity, depth, others? } or { type: 'resistance', T4 }
 * @param {boolean} [o.shortDuration=false] short-duration Van Wormer split
 * @returns {{step:(t:number)=>number, steady:number, cableSteady:number, externalSteady:number,
 *            timeConstants:number[], ladder:object, lossRatio:number, dielectricSteady:number}}
 *          step(t) and steady in K per (W/m of conductor loss per core); timeConstants of the cable in s;
 *          dielectricSteady = rise per (W/m of dielectric loss per core) in steady state
 */
export function transientModel({ geometry, lambda1 = 0, lambda2 = 0, external, shortDuration = false }) {
  if (!external) throw new RangeError('external: required');
  const g = geometry;
  const ladder = cableLadder(g, { shortDuration });
  const m = ladder.Q.length;
  const w = new Array(m).fill(0);
  w[0] = 1;
  const iS = ladder.roles.findIndex((r) => r.includes('sheath'));
  const iA = ladder.roles.findIndex((r) => r.includes('armour'));
  if (iS >= 0) w[iS] += lambda1; else w[m - 1] += lambda1;
  if (iA >= 0) w[iA] += lambda2; else w[m - 1] += lambda2;
  const lad = ladderStepResponse({ Q: ladder.Q, T: ladder.T, w });
  const cableSteady = lad.steady[0];
  const lossRatio = g.cores * (1 + lambda1 + lambda2);
  let ext; let externalSteadyUnit;
  if (external.type === 'buried') {
    pos(external.rhoSoil, 'external.rhoSoil'); pos(external.diffusivity, 'external.diffusivity'); pos(external.depth, 'external.depth');
    const base = { rhoSoil: external.rhoSoil, diffusivity: external.diffusivity, depth: external.depth, De: g.De, others: external.others ?? [] };
    ext = (t) => soilStepResponse({ ...base, t });
    externalSteadyUnit = (external.rhoSoil / TWO_PI) * (Math.log((4 * external.depth) / g.De)
      + (external.others ?? []).reduce((s, o) => s + (o.weight ?? 1) * Math.log(o.dImage / o.d), 0));
  } else if (external.type === 'resistance') {
    if (!(external.T4 >= 0)) throw new RangeError('external.T4: must be ≥ 0');
    ext = () => external.T4;
    externalSteadyUnit = external.T4;
  } else throw new RangeError('external.type: must be buried or resistance');
  const step = (t) => {
    if (!(t > 0)) return 0;
    const c = lad.rise(t, 0);
    return c + (c / cableSteady) * lossRatio * ext(t);
  };
  const dielectricSteady = 0.5 * g.T1 + g.cores * (g.T2 + g.T3) + g.cores * externalSteadyUnit;
  return { step, steady: cableSteady + lossRatio * externalSteadyUnit, cableSteady, externalSteady: lossRatio * externalSteadyUnit,
    timeConstants: lad.timeConstants, ladder, lossRatio, dielectricSteady };
}

/**
 * Time for the step response to reach a fraction of its steady value (default 1 − 1/e: an equivalent overall
 * thermal time constant, which for buried cables is dominated by the soil).
 * @param {object} o
 * @param {{step:Function, steady:number}} o.model model from transientModel
 * @param {number} [o.fraction=0.6321] fraction of the steady rise (0–1)
 * @returns {number} time, s
 */
export function responseTime({ model, fraction = 1 - Math.exp(-1) }) {
  if (!(fraction > 0 && fraction < 1)) throw new RangeError('fraction: must be in (0, 1)');
  const target = fraction * model.steady;
  let lo = 1; let hi = 3600;
  while (model.step(hi) < target) { hi *= 4; if (hi > 1e14) throw new RangeError('model: response does not reach the requested fraction'); }
  while (model.step(lo) > target && lo > 1e-6) lo /= 4;
  for (let i = 0; i < 100; i++) { const mid = Math.sqrt(lo * hi); if (model.step(mid) < target) lo = mid; else hi = mid; }
  return Math.sqrt(lo * hi);
}

/**
 * Conductor temperature history for a current time series (dynamic rating), by superposition of step responses.
 * The current is piecewise constant: current[j] applies from j·dt to (j+1)·dt; the returned temperature k is at
 * the end of step k (time (k+1)·dt). Conductor loss per core is I²·R with R at `thetaRefC`; when `feedback` is
 * true the loss of each step is re-evaluated with R at the computed temperature (three passes), which removes
 * most of the conservatism at part load. Dielectric loss is treated as permanently applied.
 * Reference: IEC 60853-2 §4 (superposition principle); cost O(N²), N ≤ 30000 steps.
 * @param {object} o
 * @param {object} o.model model from transientModel
 * @param {number[]|Float64Array} o.current current of each step, A
 * @param {number} o.dt step length, s
 * @param {number} o.R conductor AC resistance at thetaRefC, Ω/m
 * @param {number} [o.thetaRefC=90] temperature at which R is given, °C
 * @param {number} [o.alpha20=0.00393] temperature coefficient for the feedback pass, 1/K
 * @param {boolean} [o.feedback=false] re-evaluate losses with temperature-dependent resistance
 * @param {number} [o.Wd=0] dielectric loss per core, W/m
 * @param {number} o.ambientC ambient temperature, °C
 * @param {number} [o.initialCurrent=0] current carried for a long time before t = 0 (steady preload), A
 * @returns {{conductorC:Float64Array, time:Float64Array, maxC:number, tMax:number, meanC:number}} temperatures at
 *          the end of each step (°C), their times (s), maximum and its time
 */
export function temperatureHistory({ model, current, dt, R, thetaRefC = 90, alpha20 = 0.00393, feedback = false, Wd = 0,
  ambientC, initialCurrent = 0 }) {
  const N = current.length;
  if (!(N >= 1)) throw new RangeError('current: must be a non-empty array');
  if (N > 30000) throw new RangeError('current: more than 30000 steps; resample the series to a longer step');
  pos(dt, 'dt'); pos(R, 'R');
  if (!Number.isFinite(ambientC)) throw new RangeError('ambientC: must be finite');
  const S = new Float64Array(N + 1);
  for (let k = 1; k <= N; k++) S[k] = model.step(k * dt);
  const base = ambientC + Wd * model.dielectricSteady;
  const W = new Float64Array(N);
  const theta = new Float64Array(N);
  const dW = new Float64Array(N);
  const rAt = (thC) => R * (1 + alpha20 * (thC - 20)) / (1 + alpha20 * (thetaRefC - 20));
  let W0 = initialCurrent * initialCurrent * R;
  const passes = feedback ? 3 : 1;
  for (let pass = 0; pass < passes; pass++) {
    for (let j = 0; j < N; j++) {
      // Resistance at the mean of the temperatures bounding the step (previous pass).
      const r = pass === 0 ? R : rAt(0.5 * (theta[j] + (j > 0 ? theta[j - 1] : theta[0])));
      W[j] = current[j] * current[j] * r;
    }
    if (pass > 0 && initialCurrent > 0) {
      // Preload steady state with its own temperature-dependent resistance.
      let th = base + W0 * model.steady;
      for (let i = 0; i < 50; i++) th = base + initialCurrent * initialCurrent * rAt(th) * model.steady;
      W0 = (th - base) / model.steady;
    }
    dW[0] = W[0] - W0;
    for (let j = 1; j < N; j++) dW[j] = W[j] - W[j - 1];
    for (let k = 0; k < N; k++) {
      let s = 0;
      for (let j = 0; j <= k; j++) s += dW[j] * S[k + 1 - j];
      theta[k] = base + W0 * model.steady + s;
    }
  }
  let maxC = -Infinity; let iMax = 0; let sum = 0;
  const time = new Float64Array(N);
  for (let k = 0; k < N; k++) { time[k] = (k + 1) * dt; sum += theta[k]; if (theta[k] > maxC) { maxC = theta[k]; iMax = k; } }
  return { conductorC: theta, time, maxC, tMax: time[iMax], meanC: sum / N };
}

/**
 * Statistics of a conductor temperature history against a limit: probability of exceedance (fraction of time
 * steps above the limit), hours above the limit and the margin of the hottest step. Used to judge how much of the
 * thermal inertia of a buried cable can be credited for an intermittent (wind) load.
 * @param {object} o
 * @param {ArrayLike<number>} o.conductorC temperature history, °C
 * @param {number} o.limitC limit temperature, °C
 * @param {number} [o.dt] step length, s (for durationAbove)
 * @returns {{probability:number, durationAbove:number, maxC:number, marginK:number, p99C:number}}
 *          probability (0–1), duration above the limit (s; NaN without dt), maximum (°C), limit − max (K) and the
 *          99th percentile temperature (°C)
 */
export function exceedanceStatistics({ conductorC, limitC, dt }) {
  const N = conductorC.length;
  if (!(N >= 1)) throw new RangeError('conductorC: must be a non-empty array');
  if (!Number.isFinite(limitC)) throw new RangeError('limitC: must be finite');
  let above = 0; let maxC = -Infinity;
  for (let i = 0; i < N; i++) { if (conductorC[i] > limitC) above++; if (conductorC[i] > maxC) maxC = conductorC[i]; }
  const sorted = Float64Array.from(conductorC).sort();
  const p99C = sorted[Math.min(N - 1, Math.floor(0.99 * (N - 1) + 0.5))];
  return { probability: above / N, durationAbove: dt ? above * dt : NaN, maxC, marginK: limitC - maxC, p99C };
}

/**
 * Cyclic rating factor M for a repeating load cycle, directly from the numerical model: the periodic-steady
 * conductor rise for the cycle with unit peak is found by superposition over `nCycles` repetitions and
 * M = sqrt(θsteady(peak current) / θpeak(cycle)), so that M × (steady-state rating) is the permissible peak current.
 * Resistance variation with temperature and dielectric loss are ignored, as in IEC 60853.
 * @param {object} o
 * @param {object} o.model model from transientModel
 * @param {number[]} o.cycle load cycle as current / peak current per equal step (e.g. 24 hourly values)
 * @param {number} [o.dt=3600] step length, s
 * @param {number} [o.nCycles=30] repetitions used to reach the periodic state
 * @returns {{M:number, lossLoadFactor:number, peakRiseRatio:number}} M ≥ 1, loss-load factor μ of the cycle and
 *          θpeak/θsteady
 */
export function cyclicRatingFactor({ model, cycle, dt = 3600, nCycles = 30 }) {
  const n = cycle.length;
  if (!(n >= 2)) throw new RangeError('cycle: at least two steps required');
  let peak = 0;
  for (const c of cycle) { if (!(c >= 0)) throw new RangeError('cycle: values must be ≥ 0'); if (c > peak) peak = c; }
  if (!(peak > 0)) throw new RangeError('cycle: peak must be positive');
  const y = cycle.map((c) => (c / peak) * (c / peak));
  const mu = y.reduce((s, v) => s + v, 0) / n;
  // Zero-mean part of the loss cycle superposed on the steady response to the mean loss: converges to the periodic
  // state as soon as the slope of the step response has died away over one cycle.
  const N = n * nCycles;
  const S = new Float64Array(N + 1);
  for (let k = 1; k <= N; k++) S[k] = model.step(k * dt);
  let worst = -Infinity;
  for (let k = N - n; k < N; k++) {
    let s = 0;
    for (let j = 0; j <= k; j++) {
      const dW = (y[j % n] - mu) - (j === 0 ? 0 : (y[(j - 1) % n] - mu));
      s += dW * S[k + 1 - j];
    }
    const th = mu * model.steady + s;
    if (th > worst) worst = th;
  }
  const ratio = worst / model.steady;
  return { M: 1 / Math.sqrt(ratio), lossLoadFactor: mu, peakRiseRatio: ratio };
}

/**
 * Cyclic rating factor by the IEC 60853 formula
 *   M = 1/sqrt( Σ_{i=0..5} Y_i·[θR(i+1)/θR(∞) − θR(i)/θR(∞)] + μ·[1 − θR(6)/θR(∞)] ),
 * where Y_i is the loss-load ordinate (I/Imax)² of the hour ending i hours before the time of maximum conductor
 * temperature (Y0 is the peak hour), μ the loss-load factor of the 24 h cycle and θR(i) the conductor rise i hours
 * after a step of rated current. Reference: IEC 60853-1 §3 / IEC 60853-2 §5.
 * @param {object} o
 * @param {object} o.model model from transientModel
 * @param {number[]} o.Y six loss-load ordinates Y0…Y5 (0–1), Y0 normally 1
 * @param {number} o.mu loss-load factor of the cycle (0–1)
 * @returns {number} M
 */
export function cyclicRatingFactorIEC({ model, Y, mu }) {
  if (!Array.isArray(Y) || Y.length !== 6) throw new RangeError('Y: six ordinates Y0…Y5 required');
  if (!(mu >= 0 && mu <= 1)) throw new RangeError('mu: must be in [0, 1]');
  const th = [0];
  for (let i = 1; i <= 6; i++) th.push(model.step(i * 3600) / model.steady);
  let s = mu * (1 - th[6]);
  for (let i = 0; i < 6; i++) s += Y[i] * (th[i + 1] - th[i]);
  return 1 / Math.sqrt(s);
}

/**
 * Emergency rating: the current that may be carried for a given time, starting from a steady preload, without
 * the conductor exceeding its emergency temperature.
 *   θ(t) = θamb + θd + I1²·R·S(∞) + (I2² − I1²)·R·S(t)  ⇒  I2 = sqrt(I1² + (θmax − θ1)/(R·S(t))).
 * R is held at its value for θmax for both the preload and the overload (conservative; IEC 60853-2 §6 gives a
 * correction for resistance variation that is not applied here).
 * @param {object} o
 * @param {object} o.model model from transientModel
 * @param {number} o.R conductor AC resistance at thetaMaxC, Ω/m
 * @param {number} o.preloadCurrent steady current before the emergency I1, A
 * @param {number} o.duration duration of the emergency load, s
 * @param {number} o.thetaMaxC permissible conductor temperature at the end of the emergency, °C
 * @param {number} o.ambientC ambient temperature, °C
 * @param {number} [o.Wd=0] dielectric loss per core, W/m
 * @returns {{I:number, preloadC:number, steadyRating:number, ratio:number}} emergency current (A), conductor
 *          temperature before the emergency (°C), continuous rating on the same basis (A) and I/steadyRating
 */
export function emergencyRating({ model, R, preloadCurrent, duration, thetaMaxC, ambientC, Wd = 0 }) {
  pos(R, 'R'); pos(duration, 'duration');
  if (!(preloadCurrent >= 0)) throw new RangeError('preloadCurrent: must be ≥ 0');
  const base = ambientC + Wd * model.dielectricSteady;
  const preloadC = base + preloadCurrent * preloadCurrent * R * model.steady;
  if (preloadC > thetaMaxC) throw new RangeError('preloadCurrent: preload temperature already exceeds thetaMaxC');
  const I = Math.sqrt(preloadCurrent * preloadCurrent + (thetaMaxC - preloadC) / (R * model.step(duration)));
  const steadyRating = Math.sqrt(Math.max(0, thetaMaxC - base) / (R * model.steady));
  return { I, preloadC, steadyRating, ratio: I / steadyRating };
}

/**
 * Adiabatic conductor temperature after a short-circuit: θf = (θi + β)·exp(I²·t/(K²·S²)) − β.
 * Hook for checking a transient thermal state against a following fault (IEC 60949, adiabatic method).
 * @param {object} o
 * @param {number} o.I short-circuit current (rms), A
 * @param {number} o.duration s
 * @param {number} o.area conductor cross-section, m²
 * @param {string} o.material metal key ('Cu', 'Al', 'lead', 'steel')
 * @param {number} o.thetaInitialC temperature at the start of the fault, °C
 * @returns {number} final temperature, °C
 */
export function adiabaticFinalTemperature({ I, duration, area, material, thetaInitialC }) {
  pos(area, 'area');
  if (!(I >= 0)) throw new RangeError('I: must be ≥ 0');
  if (!(duration >= 0)) throw new RangeError('duration: must be ≥ 0');
  const m = metalProperties(material);
  if (!Number.isFinite(m.K)) throw new RangeError('material: no IEC 60949 constants for this metal');
  const S = area * 1e6;
  return (thetaInitialC + m.beta) * Math.exp((I * I * duration) / (m.K * m.K * S * S)) - m.beta;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* 2-D steady conduction in a trench cross-section                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

function gradedLines(lo, hi, fineLo, fineHi, h, ratio) {
  // Cell edges: uniform spacing h on [fineLo, fineHi], growing geometrically by `ratio` towards lo and hi.
  const nFine = Math.max(1, Math.round((fineHi - fineLo) / h));
  const hf = (fineHi - fineLo) / nFine;
  const edges = [];
  for (let i = 0; i <= nFine; i++) edges.push(fineLo + i * hf);
  const grow = (from, to, dir) => {
    const total = Math.abs(to - from);
    if (total < 1e-12) return [];
    // Number of cells of the geometric series that covers the distance, then exact rescaling to the boundary.
    let nCells = 0; let acc = 0; let s = hf;
    while (acc < total - 1e-12 && nCells < 100000) { s *= ratio; acc += s; nCells++; }
    const scale = total / acc;
    const out = [];
    let x = from; let step = hf;
    for (let i = 0; i < nCells; i++) { step *= ratio; x += dir * step * scale; out.push(x); }
    out[out.length - 1] = to;
    return out;
  };
  const left = grow(fineLo, lo, -1).reverse();
  const right = grow(fineHi, hi, 1);
  return Float64Array.from([...left, ...edges, ...right]);
}

/**
 * Steady 2-D heat conduction in a vertical cross-section of the seabed/ground with several cables and regions of
 * different thermal resistivity (cell-centred finite volumes on a graded tensor-product grid, harmonic-mean face
 * conductances, conjugate gradients with a modified incomplete Cholesky preconditioner).
 * The grid is uniform (spacing h) over the box enclosing the cables and backfill regions plus a margin, and grows
 * geometrically (ratio) to the far boundaries. Cable heat is distributed over the cells covered by each cable in
 * proportion to the covered area and rescaled so the total is exact. Refinement level r halves h and takes the
 * 2^r-th root of the growth ratio, i.e. a uniform refinement of the whole grid, for grid-convergence studies.
 * Coordinates: x horizontal (m), z depth below the surface (m, positive down).
 * Top boundary: isothermal at ambient, or convective with coefficient h (W/(m²·K)); sides and bottom: adiabatic
 * or isothermal at ambient. The domain should extend at least ~20 burial depths so the truncation is negligible
 * against the semi-infinite image solution.
 * Reference: Patankar, Numerical Heat Transfer and Fluid Flow (1980), ch. 4; used to cross-check IEC 60287-2-1 T4.
 * @param {object} o
 * @param {{x:number, depth:number, diameter:number, W:number, rhoT?:number}[]} o.cables cable axis position (m),
 *        diameter (m), heat output (W/m) and optional internal thermal resistivity (K·m/W; default = local soil,
 *        which reproduces the line-source superposition and converges at second order; a small value such as
 *        0.01 approximates an isothermal cable surface, with a stair-step boundary that is first-order accurate —
 *        about 2 % on T4 with 8 cells per diameter)
 * @param {number} o.rhoSoil native soil resistivity, K·m/W
 * @param {{x0:number, x1:number, z0:number, z1:number, rhoT:number}[]} [o.regions] rectangular regions (backfill,
 *        rock berm, dry zone ...); later entries overwrite earlier ones
 * @param {number} [o.ambientC=0] surface / far-field temperature, °C
 * @param {{type:'isothermal'|'convective', h?:number}} [o.surface] top boundary (default isothermal)
 * @param {'adiabatic'|'isothermal'} [o.sides='adiabatic'] @param {'adiabatic'|'isothermal'} [o.bottom='adiabatic']
 * @param {number} [o.width] domain width, m (default 40 × deepest cable) @param {number} [o.height] domain depth, m
 *        (default 20 × deepest cable)
 * @param {number} [o.h] fine cell size, m (default smallest cable diameter / 8)
 * @param {number} [o.ratio=1.15] geometric growth ratio of the far-field cells
 * @param {number} [o.refinement=0] uniform refinement level
 * @param {number} [o.tol=1e-10] relative residual tolerance
 * @param {{x:number, z:number}[]} [o.probes] points where the temperature is reported
 * @returns {{nx:number, nz:number, xEdges:Float64Array, zEdges:Float64Array, T:Float64Array, cables:{meanC:number,
 *            surfaceMeanC:number, maxC:number, T4:number}[], probes:number[], iterations:number, residual:number,
 *            heatBalanceError:number, hFine:number}}
 *          T is row-major (index iz·nx + ix, °C); cables[k].T4 = (surfaceMeanC − ambient)/W_k (K·m/W, meaningful for
 *          a single cable or equally loaded identical cables); heatBalanceError = |boundary outflow − ΣW|/ΣW
 */
export function solveTrench2D({ cables, rhoSoil, regions = [], ambientC = 0, surface = { type: 'isothermal' },
  sides = 'adiabatic', bottom = 'adiabatic', width, height, h, ratio = 1.15, refinement = 0, tol = 1e-10, probes = [] }) {
  if (!Array.isArray(cables) || cables.length === 0) throw new RangeError('cables: must be a non-empty array');
  pos(rhoSoil, 'rhoSoil');
  if (!(ratio > 1)) throw new RangeError('ratio: must be > 1');
  let zMax = 0; let dMin = Infinity; let bx0 = Infinity; let bx1 = -Infinity; let bz1 = 0; let bz0 = Infinity;
  for (const c of cables) {
    pos(c.diameter, 'cables[].diameter'); pos(c.depth, 'cables[].depth');
    if (c.depth < c.diameter / 2) throw new RangeError('cables[].depth: cable must be fully below the surface');
    zMax = Math.max(zMax, c.depth); dMin = Math.min(dMin, c.diameter);
    bx0 = Math.min(bx0, c.x - c.diameter); bx1 = Math.max(bx1, c.x + c.diameter);
    bz0 = Math.min(bz0, c.depth - c.diameter); bz1 = Math.max(bz1, c.depth + c.diameter);
  }
  for (const r of regions) {
    bx0 = Math.min(bx0, r.x0); bx1 = Math.max(bx1, r.x1); bz0 = Math.min(bz0, r.z0); bz1 = Math.max(bz1, r.z1);
  }
  const xc = 0.5 * (bx0 + bx1);
  const W = width ?? 40 * zMax;
  const H = height ?? 20 * zMax;
  if (H <= bz1) throw new RangeError('height: domain must extend below the cables and regions');
  if (W <= bx1 - bx0) throw new RangeError('width: domain must be wider than the cables and regions');
  const scale = 2 ** refinement;
  const hf = (h ?? dMin / 8) / scale;
  const rr = Math.pow(ratio, 1 / scale);
  const xEdges = gradedLines(xc - W / 2, xc + W / 2, bx0, bx1, hf, rr);
  const zEdges = gradedLines(0, H, Math.max(0, bz0), bz1, hf, rr);
  const nx = xEdges.length - 1; const nz = zEdges.length - 1;
  const n = nx * nz;
  if (n > 4e6) throw new RangeError('h/refinement: grid exceeds 4 million cells');
  const dx = new Float64Array(nx); const dz = new Float64Array(nz); const xm = new Float64Array(nx); const zm = new Float64Array(nz);
  for (let i = 0; i < nx; i++) { dx[i] = xEdges[i + 1] - xEdges[i]; xm[i] = 0.5 * (xEdges[i + 1] + xEdges[i]); }
  for (let j = 0; j < nz; j++) { dz[j] = zEdges[j + 1] - zEdges[j]; zm[j] = 0.5 * (zEdges[j + 1] + zEdges[j]); }

  // Conductivity field and sources.
  const k = new Float64Array(n).fill(1 / rhoSoil);
  for (const r of regions) {
    pos(r.rhoT, 'regions[].rhoT');
    for (let j = 0; j < nz; j++) {
      if (zm[j] < r.z0 || zm[j] > r.z1) continue;
      for (let i = 0; i < nx; i++) if (xm[i] >= r.x0 && xm[i] <= r.x1) k[j * nx + i] = 1 / r.rhoT;
    }
  }
  const q = new Float64Array(n);
  const cover = [];
  const NS = 8;
  let totalW = 0;
  for (const c of cables) {
    const rad = c.diameter / 2; const r2 = rad * rad;
    const cells = []; let area = 0;
    for (let j = 0; j < nz; j++) {
      if (zEdges[j + 1] < c.depth - rad || zEdges[j] > c.depth + rad) continue;
      for (let i = 0; i < nx; i++) {
        if (xEdges[i + 1] < c.x - rad || xEdges[i] > c.x + rad) continue;
        let inside = 0;
        for (let a = 0; a < NS; a++) {
          const xs = xEdges[i] + ((a + 0.5) / NS) * dx[i] - c.x;
          for (let b = 0; b < NS; b++) {
            const zs = zEdges[j] + ((b + 0.5) / NS) * dz[j] - c.depth;
            if (xs * xs + zs * zs <= r2) inside++;
          }
        }
        if (inside > 0) {
          const a = (inside / (NS * NS)) * dx[i] * dz[j];
          cells.push([j * nx + i, a, inside]); area += a;
        }
      }
    }
    if (cells.length === 0) throw new RangeError('h: grid too coarse to resolve a cable');
    for (const [idx, a, inside] of cells) {
      q[idx] += (c.W * a) / area;
      if (c.rhoT !== undefined && inside * 2 >= NS * NS) k[idx] = 1 / c.rhoT;
    }
    cover.push({ cells, area });
    totalW += c.W;
  }

  // Face conductances: east (i → i+1) and south (j → j+1), plus boundary conductances on the diagonal.
  const gE = new Float64Array(n); const gS = new Float64Array(n); const diag = new Float64Array(n); const gB = new Float64Array(n);
  for (let j = 0; j < nz; j++) {
    for (let i = 0; i < nx; i++) {
      const p = j * nx + i;
      if (i + 1 < nx) gE[p] = dz[j] / (dx[i] / (2 * k[p]) + dx[i + 1] / (2 * k[p + 1]));
      if (j + 1 < nz) gS[p] = dx[i] / (dz[j] / (2 * k[p]) + dz[j + 1] / (2 * k[p + nx]));
    }
  }
  const convective = surface && surface.type === 'convective';
  if (convective) pos(surface.h, 'surface.h');
  for (let i = 0; i < nx; i++) {
    gB[i] += dx[i] / (dz[0] / (2 * k[i]) + (convective ? 1 / surface.h : 0));
    if (bottom === 'isothermal') { const p = (nz - 1) * nx + i; gB[p] += dx[i] / (dz[nz - 1] / (2 * k[p])); }
  }
  if (sides === 'isothermal') {
    for (let j = 0; j < nz; j++) {
      const p0 = j * nx; const p1 = j * nx + nx - 1;
      gB[p0] += dz[j] / (dx[0] / (2 * k[p0])); gB[p1] += dz[j] / (dx[nx - 1] / (2 * k[p1]));
    }
  }
  for (let j = 0; j < nz; j++) {
    for (let i = 0; i < nx; i++) {
      const p = j * nx + i;
      let d = gB[p] + gE[p] + gS[p];
      if (i > 0) d += gE[p - 1];
      if (j > 0) d += gS[p - nx];
      diag[p] = d;
    }
  }
  const apply = (x, y) => {
    for (let j = 0; j < nz; j++) {
      const row = j * nx;
      for (let i = 0; i < nx; i++) {
        const p = row + i;
        let s = diag[p] * x[p];
        if (i > 0) s -= gE[p - 1] * x[p - 1];
        if (i + 1 < nx) s -= gE[p] * x[p + 1];
        if (j > 0) s -= gS[p - nx] * x[p - nx];
        if (j + 1 < nz) s -= gS[p] * x[p + nx];
        y[p] = s;
      }
    }
  };
  // Modified incomplete Cholesky MIC(0) preconditioner for the 5-point stencil (Gustafsson 1978): the dropped
  // fill-in is lumped onto the diagonal with relaxation 0.97, which keeps the iteration count nearly independent
  // of the grid grading.
  const pc = new Float64Array(n);
  for (let j = 0; j < nz; j++) {
    for (let i = 0; i < nx; i++) {
      const p = j * nx + i;
      let e = diag[p];
      if (i > 0) { const w = gE[p - 1] * pc[p - 1]; e -= w * w + 0.97 * gE[p - 1] * gS[p - 1] * pc[p - 1] * pc[p - 1]; }
      if (j > 0) { const w = gS[p - nx] * pc[p - nx]; e -= w * w + 0.97 * gS[p - nx] * gE[p - nx] * pc[p - nx] * pc[p - nx]; }
      if (e < 0.25 * diag[p]) e = diag[p];
      pc[p] = 1 / Math.sqrt(e);
    }
  }
  const precond = (rv, zv) => {
    for (let j = 0; j < nz; j++) {
      for (let i = 0; i < nx; i++) {
        const p = j * nx + i;
        let t = rv[p];
        if (i > 0) t += gE[p - 1] * pc[p - 1] * zv[p - 1];
        if (j > 0) t += gS[p - nx] * pc[p - nx] * zv[p - nx];
        zv[p] = t * pc[p];
      }
    }
    for (let j = nz - 1; j >= 0; j--) {
      for (let i = nx - 1; i >= 0; i--) {
        const p = j * nx + i;
        let t = zv[p];
        if (i + 1 < nx) t += gE[p] * pc[p] * zv[p + 1];
        if (j + 1 < nz) t += gS[p] * pc[p] * zv[p + nx];
        zv[p] = t * pc[p];
      }
    }
  };
  const x = new Float64Array(n); const r = Float64Array.from(q); const z = new Float64Array(n);
  const pv = new Float64Array(n); const Ap = new Float64Array(n);
  let rz = 0; let bnorm = 0;
  precond(r, z);
  for (let i = 0; i < n; i++) { pv[i] = z[i]; rz += r[i] * z[i]; bnorm += q[i] * q[i]; }
  bnorm = Math.sqrt(bnorm) || 1;
  let iterations = 0; let res = 1;
  const maxIter = 20 * (nx + nz) + 2000;
  for (; iterations < maxIter; iterations++) {
    apply(pv, Ap);
    let pAp = 0;
    for (let i = 0; i < n; i++) pAp += pv[i] * Ap[i];
    const alpha = rz / pAp;
    let rn = 0;
    for (let i = 0; i < n; i++) { x[i] += alpha * pv[i]; r[i] -= alpha * Ap[i]; rn += r[i] * r[i]; }
    res = Math.sqrt(rn) / bnorm;
    if (res < tol) { iterations++; break; }
    precond(r, z);
    let rzNew = 0;
    for (let i = 0; i < n; i++) rzNew += r[i] * z[i];
    const beta = rzNew / rz;
    rz = rzNew;
    for (let i = 0; i < n; i++) pv[i] = z[i] + beta * pv[i];
  }
  if (!(res < tol * 1e3)) throw new Error(`solveTrench2D: conjugate gradients did not converge (residual ${res.toExponential(2)})`);

  let outflow = 0;
  for (let i = 0; i < n; i++) outflow += gB[i] * x[i];
  const T = new Float64Array(n);
  for (let i = 0; i < n; i++) T[i] = x[i] + ambientC;

  const locate = (mids, v) => {
    // Index i such that mids[i] ≤ v < mids[i+1], clamped.
    let lo = 0; let hi = mids.length - 1;
    if (v <= mids[0]) return 0;
    if (v >= mids[hi]) return hi - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (mids[mid] <= v) lo = mid; else hi = mid; }
    return lo;
  };
  const sample = (px, pz) => {
    const i = locate(xm, px); const j = locate(zm, pz);
    const tx = Math.min(1, Math.max(0, (px - xm[i]) / (xm[i + 1] - xm[i])));
    const tz = Math.min(1, Math.max(0, (pz - zm[j]) / (zm[j + 1] - zm[j])));
    const a = T[j * nx + i]; const b = T[j * nx + i + 1]; const c = T[(j + 1) * nx + i]; const d = T[(j + 1) * nx + i + 1];
    return (1 - tz) * ((1 - tx) * a + tx * b) + tz * ((1 - tx) * c + tx * d);
  };
  const cableOut = cables.map((c, ci) => {
    let mean = 0; let maxC = -Infinity;
    for (const [idx, a] of cover[ci].cells) { mean += T[idx] * a; if (T[idx] > maxC) maxC = T[idx]; }
    mean /= cover[ci].area;
    let surf = 0; const NP = 72;
    for (let m = 0; m < NP; m++) {
      const ang = (TWO_PI * (m + 0.5)) / NP;
      surf += sample(c.x + 0.5 * c.diameter * Math.cos(ang), c.depth + 0.5 * c.diameter * Math.sin(ang));
    }
    surf /= NP;
    return { meanC: mean, surfaceMeanC: surf, maxC, T4: c.W !== 0 ? (surf - ambientC) / c.W : NaN };
  });
  return { nx, nz, xEdges, zEdges, T, cables: cableOut, probes: probes.map((p) => sample(p.x, p.z)), iterations, residual: res,
    heatBalanceError: totalW !== 0 ? Math.abs(outflow - totalW) / Math.abs(totalW) : 0, hFine: hf };
}
