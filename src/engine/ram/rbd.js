// Reliability block diagrams.
//
// Block format:
//   { type: 'component', id?, life: { dist: 'exponential', rate } | { dist: 'weibull', shape, scale, loc? } |
//                                   { dist: 'lognormal', mu, sigma } | { dist: 'constant', reliability } |
//                                   { R: (t) => number },
//                         availability? | mtbf?, mttr? | lambda?, mu? }       (for availabilityOf)
//   { type: 'series' | 'parallel', children: [block, ...] }
//   { type: 'koon', k, children: [...] }                    works when at least k children work
//   { type: 'standby', children: [...], switchReliability? } cold standby, units used in the order listed
//   { type: 'network', source, sink, edges: [{ from, to, block, id? }], directed? }   e.g. a bridge
// Time and rate units must be consistent (e.g. h and 1/h). Blocks are statistically independent.

import { erfc } from '../math/special.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

/**
 * Survival probability R(t) of a component life distribution.
 * @param {object} life see the file header
 * @param {number} t time [time unit]
 * @returns {number} reliability [-]
 */
export function lifeReliability(life, t) {
  need(life && typeof life === 'object', 'life must be an object');
  need(t >= 0, 't must be non-negative');
  if (typeof life.R === 'function') return life.R(t);
  switch (life.dist) {
    case 'exponential':
      need(life.rate >= 0, 'life.rate must be non-negative');
      return Math.exp(-life.rate * t);
    case 'weibull': {
      need(life.shape > 0 && life.scale > 0, 'life.shape and life.scale must be positive');
      const x = t - (life.loc || 0);
      return x <= 0 ? 1 : Math.exp(-((x / life.scale) ** life.shape));
    }
    case 'lognormal':
      need(life.sigma > 0, 'life.sigma must be positive');
      return t === 0 ? 1 : 0.5 * erfc((Math.log(t) - life.mu) / (life.sigma * Math.SQRT2));
    case 'constant':
      need(life.reliability >= 0 && life.reliability <= 1, 'life.reliability must be in [0, 1]');
      return life.reliability;
    default:
      throw new RangeError(`life.dist: unknown distribution "${life.dist}"`);
  }
}

function kOfN(k, probs) {
  // Distribution of the number of working items among independent, non-identical items (dynamic programme).
  const n = probs.length;
  const d = new Float64Array(n + 1);
  d[0] = 1;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j >= 1; j--) d[j] = d[j] * (1 - probs[i]) + d[j - 1] * probs[i];
    d[0] *= 1 - probs[i];
  }
  let s = 0;
  for (let j = k; j <= n; j++) s += d[j];
  return Math.min(1, s);
}

/**
 * Minimal path sets of a two-terminal network by depth-first enumeration of simple paths.
 * @param {object} o
 * @param {{from:(string|number), to:(string|number)}[]} o.edges
 * @param {string|number} o.source
 * @param {string|number} o.sink
 * @param {boolean} [o.directed=false]
 * @returns {number[][]} paths as arrays of edge indices
 */
export function minimalPaths({ edges, source, sink, directed = false }) {
  need(Array.isArray(edges) && edges.length > 0, 'edges must be a non-empty array');
  const adj = new Map();
  const add = (a, b, i) => {
    if (!adj.has(a)) adj.set(a, []);
    adj.get(a).push([b, i]);
  };
  edges.forEach((e, i) => {
    add(e.from, e.to, i);
    if (!directed) add(e.to, e.from, i);
  });
  const paths = [];
  const visited = new Set([source]);
  const stack = [];
  const dfs = (node) => {
    if (node === sink) { paths.push(stack.slice()); return; }
    for (const [next, i] of adj.get(node) || []) {
      if (visited.has(next)) continue;
      visited.add(next);
      stack.push(i);
      dfs(next);
      stack.pop();
      visited.delete(next);
    }
  };
  dfs(source);
  need(paths.length <= 24, 'network has more than 24 minimal paths; simplify it or split it into modules');
  return paths;
}

function networkProbability(block, probOf) {
  const paths = minimalPaths(block);
  if (paths.length === 0) return 0;
  // Edges that share an id are the same physical component.
  const key = block.edges.map((e, i) => (e.id !== undefined ? `id:${e.id}` : `e:${i}`));
  const uniq = [...new Set(key)];
  const pos = new Map(uniq.map((k, i) => [k, i]));
  const p = new Array(uniq.length);
  block.edges.forEach((e, i) => { p[pos.get(key[i])] = probOf(e.block); });
  const masks = paths.map((path) => path.reduce((m, i) => m | (1 << pos.get(key[i])), 0));
  need(uniq.length <= 30, 'network is limited to 30 distinct components');
  const n = masks.length;
  let total = 0;
  for (let s = 1; s < 1 << n; s++) {
    let u = 0;
    let bits = 0;
    for (let i = 0; i < n; i++) if (s & (1 << i)) { u |= masks[i]; bits++; }
    let prod = 1;
    for (let j = 0; j < uniq.length; j++) if (u & (1 << j)) prod *= p[j];
    total += bits % 2 === 1 ? prod : -prod;
  }
  return Math.min(1, Math.max(0, total));
}

function structure(block, probOf, standbyOf) {
  need(block && typeof block === 'object', 'block must be an object');
  switch (block.type) {
    case 'series': return block.children.reduce((a, c) => a * probOf(c), 1);
    case 'parallel': return 1 - block.children.reduce((a, c) => a * (1 - probOf(c)), 1);
    case 'koon':
      need(Number.isInteger(block.k) && block.k >= 0 && block.k <= block.children.length, 'koon block needs integer k in [0, n]');
      return kOfN(block.k, block.children.map(probOf));
    case 'network': return networkProbability(block, probOf);
    case 'standby': return standbyOf(block);
    default: throw new RangeError(`block.type: unknown block type "${block.type}"`);
  }
}

/**
 * System reliability R_sys(t) of a block diagram.
 *  series: Π Rᵢ; parallel: 1 − Π(1 − Rᵢ); k-out-of-n with non-identical items by dynamic programming;
 *  cold standby by numerical convolution, R(t) = R₁(t) + p_sw·∫₀ᵗ f₁(x)·R_rest(t − x) dx, the standby unit
 *  not ageing until switched in (perfect switching p_sw = 1 or imperfect p_sw < 1 per switch-over);
 *  network: exact union of the minimal path sets by inclusion–exclusion.
 * Reference: Rausand & Høyland, System Reliability Theory (2nd ed.); IEC 61078.
 * @param {object} block see the file header
 * @param {number} t mission time [time unit]
 * @param {object} [o]
 * @param {number} [o.standbySteps=400] convolution steps for standby blocks [-]
 * @returns {number} reliability [-]
 */
export function systemReliability(block, t, { standbySteps = 400 } = {}) {
  need(t >= 0, 't must be non-negative');
  const rel = (b, time) => {
    if (b.type === 'component') return lifeReliability(b.life, time);
    return structure(b, (c) => rel(c, time), (sb) => standby(sb, time));
  };
  const standby = (sb, time) => {
    need(Array.isArray(sb.children) && sb.children.length > 0, 'standby block needs children');
    const ps = sb.switchReliability === undefined ? 1 : sb.switchReliability;
    need(ps >= 0 && ps <= 1, 'switchReliability must be in [0, 1]');
    if (time === 0) return 1;
    const n = standbySteps;
    const h = time / n;
    const kids = sb.children;
    let rest = new Float64Array(n + 1);
    for (let k = 0; k <= n; k++) rest[k] = rel(kids[kids.length - 1], k * h);
    for (let u = kids.length - 2; u >= 0; u--) {
      const Ru = new Float64Array(n + 1);
      for (let k = 0; k <= n; k++) Ru[k] = rel(kids[u], k * h);
      const next = new Float64Array(n + 1);
      for (let k = 0; k <= n; k++) {
        let s = 0;
        for (let j = 0; j < k; j++) s += (Ru[j] - Ru[j + 1]) * 0.5 * (rest[k - j] + rest[k - j - 1]);
        next[k] = Ru[k] + ps * s;
      }
      rest = next;
    }
    return rest[n];
  };
  return Math.min(1, Math.max(0, rel(block, t)));
}

/**
 * Reliability function over a set of times.
 * @param {object} block
 * @param {ArrayLike<number>} times [time unit]
 * @returns {Float64Array} R_sys at each time [-]
 */
export function reliabilityCurve(block, times) {
  const out = new Float64Array(times.length);
  for (let i = 0; i < times.length; i++) out[i] = systemReliability(block, times[i]);
  return out;
}

/**
 * Mean time to failure MTTF = ∫₀^∞ R(t) dt by Simpson quadrature on geometrically growing panels.
 * Reference: Rausand & Høyland
 * @param {object|function(number):number} blockOrR block diagram or a reliability function R(t)
 * @param {object} [o]
 * @param {number} [o.scale=1] initial guess of the time scale [time unit]
 * @param {number} [o.tail=1e-10] reliability at which the integration stops [-]
 * @returns {number} MTTF [time unit]; Infinity when R does not decay
 */
export function mttf(blockOrR, { scale = 1, tail = 1e-10 } = {}) {
  need(scale > 0, 'scale must be positive');
  const R = typeof blockOrR === 'function' ? blockOrR : (t) => systemReliability(blockOrR, t, { standbySteps: 120 });
  let tMax = scale;
  let k = 0;
  while (R(tMax) > tail) {
    tMax *= 2;
    if (++k > 200) return Infinity;
  }
  while (tMax > 1e-300 && R(tMax / 2) <= tail) tMax /= 2;
  const panels = 40;
  let total = 0;
  let b = tMax;
  for (let p = 0; p < panels; p++) {
    const a = p === panels - 1 ? 0 : b / 2;
    const n = 16;
    const h = (b - a) / n;
    let s = R(a) + R(b);
    for (let i = 1; i < n; i++) s += (i % 2 === 1 ? 4 : 2) * R(a + i * h);
    total += (s * h) / 3;
    b = a;
  }
  return total;
}

/**
 * Steady-state availability of a repairable block diagram with independent repair of every component.
 * Component availability: `availability`, or mtbf/(mtbf + mttr), or μ/(λ + μ).
 * Standby blocks are treated as active parallel (upper-bound approximation, exact when switching is perfect and
 * standby failure rates equal active rates); use ram/markov.js for exact standby or shared-crew models.
 * Reference: Rausand & Høyland; IEC 61078.
 * @param {object} block
 * @returns {number} availability [-]
 */
export function availabilityOf(block) {
  const av = (b) => {
    if (b.type === 'component') {
      if (b.availability !== undefined) {
        need(b.availability >= 0 && b.availability <= 1, 'availability must be in [0, 1]');
        return b.availability;
      }
      if (b.mtbf !== undefined) {
        need(b.mtbf > 0 && b.mttr >= 0, 'mtbf must be positive and mttr non-negative');
        return b.mtbf / (b.mtbf + b.mttr);
      }
      need(b.lambda >= 0 && b.mu > 0, 'component needs availability, mtbf/mttr or lambda/mu');
      return b.mu / (b.lambda + b.mu);
    }
    return structure(b, av, (sb) => 1 - sb.children.reduce((a, c) => a * (1 - av(c)), 1));
  };
  return av(block);
}

/**
 * Series system of repairable components: availability, failure frequency, mean up and down time.
 * A = Π Aᵢ (independent repair); failure frequency ω = A·Σλᵢ; MUT = 1/Σλᵢ; MDT = (1 − A)/ω.
 * Reference: Rausand & Høyland (frequency–duration relations).
 * @param {{lambda:number, mu:number}[]} components failure and repair rates [1/time]
 * @returns {{availability:number, failureFrequency:number, meanUpTime:number, meanDownTime:number}}
 */
export function seriesRepairable(components) {
  need(Array.isArray(components) && components.length > 0, 'components must be a non-empty array');
  let A = 1;
  let lam = 0;
  for (const c of components) {
    need(c.lambda >= 0 && c.mu > 0, 'each component needs lambda >= 0 and mu > 0');
    A *= c.mu / (c.lambda + c.mu);
    lam += c.lambda;
  }
  const w = A * lam;
  return { availability: A, failureFrequency: w, meanUpTime: lam > 0 ? 1 / lam : Infinity, meanDownTime: w > 0 ? (1 - A) / w : 0 };
}

/**
 * Parallel system of repairable components (system down only when all are down), independent repair:
 * U = Π Uᵢ, failure frequency ω = U·Σμᵢ, MDT = 1/Σμᵢ, MUT = (1 − U)/ω.
 * Reference: Rausand & Høyland
 * @param {{lambda:number, mu:number}[]} components [1/time]
 * @returns {{availability:number, failureFrequency:number, meanUpTime:number, meanDownTime:number}}
 */
export function parallelRepairable(components) {
  need(Array.isArray(components) && components.length > 0, 'components must be a non-empty array');
  let U = 1;
  let mu = 0;
  for (const c of components) {
    need(c.lambda >= 0 && c.mu > 0, 'each component needs lambda >= 0 and mu > 0');
    U *= c.lambda / (c.lambda + c.mu);
    mu += c.mu;
  }
  const w = U * mu;
  return { availability: 1 - U, failureFrequency: w, meanUpTime: w > 0 ? (1 - U) / w : Infinity, meanDownTime: 1 / mu };
}
