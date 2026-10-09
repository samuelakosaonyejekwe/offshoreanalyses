// Fault trees (exact quantification by binary decision diagram, minimal cut sets by MOCUS), importance measures,
// event trees and FMEA/FMECA tables.
//
// Tree format:
//   { type: 'or' | 'and', children: [node, ...] }
//   { type: 'vote', k, children: [...] }              true when at least k children are true
//   { type: 'not', child: node }                      (makes the tree non-coherent: no cut sets)
//   { type: 'basic', id, ...probability model }       repeated events share an id
// Probability model of a basic event (on the node or in the `events` map keyed by id):
//   { p }                               fixed probability
//   { lambda, t }                       non-repairable: 1 − exp(−λt)
//   { lambda, mu }                      repairable, steady state: λ/(λ + μ)
//   { lambda, mu, t }                   repairable at time t: λ/(λ + μ)·(1 − exp(−(λ + μ)t))
//   { lambda, testInterval }            periodically tested, mean unavailability: 1 − (1 − exp(−λτ))/(λτ)
// Rates and times must use consistent units (e.g. 1/h and h).

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

/**
 * Probability (unavailability / unreliability) of a basic event from its model, see the file header.
 * Reference: NUREG-0492 Fault Tree Handbook; IEC 61025; Rausand & Høyland, System Reliability Theory (2nd ed.).
 * @param {object} e basic-event model
 * @returns {number} probability [-]
 */
export function basicEventProbability(e) {
  need(e && typeof e === 'object', 'basic event model must be an object');
  if (e.p !== undefined) {
    need(e.p >= 0 && e.p <= 1, 'p must be in [0, 1]');
    return e.p;
  }
  need(e.lambda >= 0, 'lambda must be non-negative');
  if (e.mu !== undefined) {
    need(e.mu >= 0, 'mu must be non-negative');
    const s = e.lambda + e.mu;
    if (s === 0) return 0;
    const q = e.lambda / s;
    return e.t === undefined ? q : q * (1 - Math.exp(-s * e.t));
  }
  if (e.testInterval !== undefined) {
    need(e.testInterval > 0, 'testInterval must be positive');
    const x = e.lambda * e.testInterval;
    return x < 1e-8 ? x / 2 : 1 - (1 - Math.exp(-x)) / x;
  }
  need(e.t >= 0, 'basic event needs p, or lambda with t, mu or testInterval');
  return 1 - Math.exp(-e.lambda * e.t);
}

function collect(tree, events) {
  const ids = [];
  const index = new Map();
  const model = [];
  let coherent = true;
  const visit = (n) => {
    need(n && typeof n === 'object', 'tree node must be an object');
    if (n.type === 'basic') {
      need(n.id !== undefined, 'basic event needs an id');
      if (!index.has(n.id)) {
        index.set(n.id, ids.length);
        ids.push(n.id);
        model.push(null);
      }
      const own = n.p !== undefined || n.lambda !== undefined;
      const i = index.get(n.id);
      if (own) model[i] = n;
      else if (model[i] === null && events && events[n.id]) model[i] = events[n.id];
      return;
    }
    if (n.type === 'not') {
      coherent = false;
      visit(n.child);
      return;
    }
    need(n.type === 'and' || n.type === 'or' || n.type === 'vote', `unknown gate type "${n.type}"`);
    need(Array.isArray(n.children) && n.children.length > 0, `gate "${n.type}" needs children`);
    if (n.type === 'vote') need(Number.isInteger(n.k) && n.k >= 0 && n.k <= n.children.length, 'vote gate needs integer k in [0, n]');
    n.children.forEach(visit);
  };
  visit(tree);
  const p = model.map((m, i) => {
    need(m !== null, `basic event "${ids[i]}" has no probability model`);
    return basicEventProbability(m);
  });
  return { ids, index, p, coherent };
}

// Reduced ordered binary decision diagram with terminals 0 (false) and 1 (true).
function makeBdd() {
  const vr = [Infinity, Infinity];
  const lo = [0, 1];
  const hi = [0, 1];
  const unique = new Map();
  const memo = new Map();
  const mk = (v, l, h) => {
    if (l === h) return l;
    const key = `${v},${l},${h}`;
    let id = unique.get(key);
    if (id === undefined) {
      id = vr.length;
      vr.push(v); lo.push(l); hi.push(h);
      unique.set(key, id);
    }
    return id;
  };
  const apply = (op, f, g) => {
    if (op === 0) { // and
      if (f === 0 || g === 0) return 0;
      if (f === 1) return g;
      if (g === 1) return f;
    } else {
      if (f === 1 || g === 1) return 1;
      if (f === 0) return g;
      if (g === 0) return f;
    }
    if (f === g) return f;
    const a = f < g ? f : g;
    const b = f < g ? g : f;
    const key = `${op},${a},${b}`;
    let r = memo.get(key);
    if (r !== undefined) return r;
    const v = Math.min(vr[a], vr[b]);
    const al = vr[a] === v ? lo[a] : a;
    const ah = vr[a] === v ? hi[a] : a;
    const bl = vr[b] === v ? lo[b] : b;
    const bh = vr[b] === v ? hi[b] : b;
    r = mk(v, apply(op, al, bl), apply(op, ah, bh));
    memo.set(key, r);
    return r;
  };
  const not = (f) => {
    if (f < 2) return 1 - f;
    const key = `n,${f}`;
    let r = memo.get(key);
    if (r === undefined) {
      r = mk(vr[f], not(lo[f]), not(hi[f]));
      memo.set(key, r);
    }
    return r;
  };
  const variable = (v) => mk(v, 0, 1);
  const prob = (f, p) => {
    const cache = new Map();
    const go = (x) => {
      if (x < 2) return x;
      let r = cache.get(x);
      if (r === undefined) {
        const pv = p[vr[x]];
        r = pv * go(hi[x]) + (1 - pv) * go(lo[x]);
        cache.set(x, r);
      }
      return r;
    };
    return go(f);
  };
  return { apply, not, variable, prob, size: () => vr.length };
}

function buildBdd(bdd, node, index) {
  switch (node.type) {
    case 'basic': return bdd.variable(index.get(node.id));
    case 'not': return bdd.not(buildBdd(bdd, node.child, index));
    case 'and': return node.children.reduce((acc, c) => bdd.apply(0, acc, buildBdd(bdd, c, index)), 1);
    case 'or': return node.children.reduce((acc, c) => bdd.apply(1, acc, buildBdd(bdd, c, index)), 0);
    default: {
      const kids = node.children.map((c) => buildBdd(bdd, c, index));
      const n = kids.length;
      const memo = new Map();
      // atLeast(k, i): at least k of children i..n-1 are true.
      const atLeast = (k, i) => {
        if (k <= 0) return 1;
        if (n - i < k) return 0;
        const key = k * (n + 1) + i;
        let r = memo.get(key);
        if (r === undefined) {
          const yes = bdd.apply(0, kids[i], atLeast(k - 1, i + 1));
          const no = bdd.apply(0, bdd.not(kids[i]), atLeast(k, i + 1));
          r = bdd.apply(1, yes, no);
          memo.set(key, r);
        }
        return r;
      };
      return atLeast(node.k, 0);
    }
  }
}

function minimise(sets) {
  sets.sort((a, b) => a.length - b.length);
  const out = [];
  const seen = new Set();
  for (const s of sets) {
    const key = s.join(',');
    if (seen.has(key)) continue;
    let absorbed = false;
    for (const m of out) {
      if (m.length > s.length) break;
      let i = 0;
      for (let j = 0; j < s.length && i < m.length; j++) if (s[j] === m[i]) i++;
      if (i === m.length) { absorbed = true; break; }
    }
    if (!absorbed) { out.push(s); seen.add(key); }
  }
  return out;
}

function unionSorted(a, b) {
  const out = [];
  let i = 0;
  let j = 0;
  while (i < a.length || j < b.length) {
    if (j >= b.length || (i < a.length && a[i] < b[j])) out.push(a[i++]);
    else if (i >= a.length || b[j] < a[i]) out.push(b[j++]);
    else { out.push(a[i]); i++; j++; }
  }
  return out;
}

function mocus(node, index, limit) {
  if (node.type === 'basic') return [[index.get(node.id)]];
  need(node.type !== 'not', 'minimal cut sets are defined for coherent trees only (no NOT gates)');
  const kids = node.children.map((c) => mocus(c, index, limit));
  const andSets = (lists) => {
    let acc = [[]];
    for (const l of lists) {
      const next = [];
      for (const a of acc) for (const b of l) next.push(unionSorted(a, b));
      need(next.length <= 50 * limit, 'maxCutSets: cut-set expansion exceeds the limit');
      acc = minimise(next);
    }
    return acc;
  };
  let res;
  if (node.type === 'or') res = minimise([].concat(...kids));
  else if (node.type === 'and') res = andSets(kids);
  else {
    const n = kids.length;
    if (node.k === 0) return [[]];
    const all = [];
    const pick = [];
    const rec = (start) => {
      if (pick.length === node.k) {
        for (const s of andSets(pick.map((i) => kids[i]))) all.push(s);
        return;
      }
      for (let i = start; i <= n - (node.k - pick.length); i++) {
        pick.push(i);
        rec(i + 1);
        pick.pop();
      }
    };
    rec(0);
    res = minimise(all);
  }
  need(res.length <= limit, 'maxCutSets: number of cut sets exceeds the limit');
  return res;
}

/**
 * Minimal cut sets of a coherent fault tree by top-down expansion with absorption (MOCUS).
 * Reference: Fussell & Vesely (1972) Trans. ANS 15, 262–263; NUREG-0492.
 * @param {object} tree fault tree (no NOT gates)
 * @param {object} [o]
 * @param {number} [o.maxCutSets=20000] safety limit on the expansion [-]
 * @returns {string[][]} minimal cut sets as arrays of basic-event ids, ordered by size
 */
export function minimalCutSets(tree, { maxCutSets = 20000 } = {}) {
  const ids = [];
  const index = new Map();
  const visit = (n) => {
    if (n.type === 'basic') {
      if (!index.has(n.id)) { index.set(n.id, ids.length); ids.push(n.id); }
    } else if (n.type === 'not') visit(n.child);
    else {
      need(Array.isArray(n.children), `gate "${n.type}" needs children`);
      n.children.forEach(visit);
    }
  };
  visit(tree);
  return mocus(tree, index, maxCutSets).map((s) => s.map((i) => ids[i]));
}

/**
 * Probability of the union of cut sets by inclusion–exclusion (exact for independent basic events), for
 * cross-checking; cost grows as 2^n, so it is limited to 20 cut sets.
 * @param {string[][]} cutSets cut sets (arrays of ids)
 * @param {Object<string, number>} p probability of each basic event [-]
 * @returns {number} probability [-]
 */
export function inclusionExclusion(cutSets, p) {
  need(Array.isArray(cutSets), 'cutSets must be an array');
  need(cutSets.length <= 20, 'inclusion-exclusion is limited to 20 cut sets');
  const n = cutSets.length;
  let total = 0;
  for (let mask = 1; mask < 1 << n; mask++) {
    const u = new Set();
    let bits = 0;
    for (let i = 0; i < n; i++) if (mask & (1 << i)) { bits++; for (const id of cutSets[i]) u.add(id); }
    let prod = 1;
    for (const id of u) {
      need(p[id] !== undefined, `p has no entry for "${id}"`);
      prod *= p[id];
    }
    total += bits % 2 === 1 ? prod : -prod;
  }
  return total;
}

/**
 * Quantify a fault tree: exact top-event probability (binary decision diagram, valid with repeated events,
 * voting and NOT gates), minimal cut sets, rare-event and min-cut upper bounds, and importance measures.
 *   Birnbaum          I_B(i) = P(top | i failed) − P(top | i working)
 *   Criticality       I_C(i) = I_B(i)·p_i/P(top)
 *   Fussell–Vesely    I_FV(i) = P(union of minimal cut sets containing i)/P(top)   (exact, by BDD)
 *   Risk achievement worth RAW = P(top | i failed)/P(top); risk reduction worth RRW = P(top)/P(top | i working)
 * Reference: Rauzy (1993) Reliability Engineering & System Safety 40, 203–211 (BDD); NUREG-0492; IEC 61025;
 * Rausand & Høyland (importance).
 * Validity: statistically independent basic events (model common causes as explicit shared events).
 * @param {object} tree fault tree, see the file header
 * @param {object} [o]
 * @param {Object<string, object>} [o.events] probability models keyed by basic-event id
 * @param {number} [o.maxCutSets=20000]
 * @param {boolean} [o.cutSets=true] compute cut sets and the measures that need them
 * @returns {{probability:number, rareEvent:(number|null), minCutUpperBound:(number|null), cutSets:(object[]|null),
 *   importance: {id:string, p:number, birnbaum:number, criticality:number, fussellVesely:(number|null), raw:number,
 *   rrw:number}[], basicEvents: {id:string, p:number}[], coherent:boolean}} cutSets entries are
 *   { events: string[], probability, order }, sorted by probability
 */
export function evaluateFaultTree(tree, { events, maxCutSets = 20000, cutSets = true } = {}) {
  const { ids, index, p, coherent } = collect(tree, events);
  const bdd = makeBdd();
  const root = buildBdd(bdd, tree, index);
  const top = bdd.prob(root, p);
  let sets = null;
  let rare = null;
  let mcub = null;
  let csOut = null;
  if (cutSets && coherent) {
    sets = mocus(tree, index, maxCutSets);
    rare = 0;
    let surv = 1;
    csOut = sets.map((s) => {
      const q = s.reduce((a, i) => a * p[i], 1);
      rare += q;
      surv *= 1 - q;
      return { events: s.map((i) => ids[i]), probability: q, order: s.length };
    }).sort((a, b) => b.probability - a.probability);
    mcub = 1 - surv;
  }
  const importance = ids.map((id, i) => {
    const q = p.slice();
    q[i] = 1;
    const p1 = bdd.prob(root, q);
    q[i] = 0;
    const p0 = bdd.prob(root, q);
    let fv = null;
    if (sets) {
      let f = 0;
      for (const s of sets) {
        if (!s.includes(i)) continue;
        f = bdd.apply(1, f, s.reduce((acc, v) => bdd.apply(0, acc, bdd.variable(v)), 1));
      }
      fv = top > 0 ? bdd.prob(f, p) / top : 0;
    }
    return {
      id, p: p[i], birnbaum: p1 - p0, criticality: top > 0 ? ((p1 - p0) * p[i]) / top : 0, fussellVesely: fv,
      raw: top > 0 ? p1 / top : Infinity, rrw: p0 > 0 ? top / p0 : Infinity,
    };
  });
  return { probability: top, rareEvent: rare, minCutUpperBound: mcub, cutSets: csOut, importance, basicEvents: ids.map((id, i) => ({ id, p: p[i] })), coherent };
}

/**
 * Event-tree evaluation. Each node asks one question (barrier demand); each branch carries a conditional
 * probability and leads to another node or to an outcome with an optional consequence value.
 *   node:   { name, branches: [ { label, probability, next: node } | { label, probability, outcome, consequence } ] }
 * Branch probabilities at a node must sum to 1 (one may be omitted and is then taken as the complement).
 * Reference: IEC 62502; NUREG/CR-2300 (PRA Procedures Guide).
 * @param {object} o
 * @param {number} o.initiatingFrequency frequency of the initiating event [1/yr]
 * @param {object} o.tree root node
 * @returns {{sequences: {path:string[], probability:number, frequency:number, outcome:string, consequence:number,
 *   risk:number}[], outcomes: Object<string, {frequency:number, probability:number, risk:number}>, totalRisk:number}}
 *   frequency [1/yr]; risk = frequency × consequence [consequence unit/yr]
 */
export function evaluateEventTree({ initiatingFrequency, tree }) {
  need(initiatingFrequency >= 0, 'initiatingFrequency must be non-negative');
  const sequences = [];
  const walk = (node, prob, path) => {
    need(node && Array.isArray(node.branches) && node.branches.length > 0, 'event-tree node needs branches');
    let known = 0;
    let missing = 0;
    for (const b of node.branches) {
      if (b.probability === undefined) missing++;
      else {
        need(b.probability >= 0 && b.probability <= 1, `branch probability at "${node.name}" must be in [0, 1]`);
        known += b.probability;
      }
    }
    need(missing <= 1, `at most one branch probability may be omitted at "${node.name}"`);
    need(missing === 1 ? known <= 1 + 1e-9 : Math.abs(known - 1) < 1e-9, `branch probabilities at "${node.name}" must sum to 1`);
    for (const b of node.branches) {
      const pb = b.probability === undefined ? Math.max(0, 1 - known) : b.probability;
      const next = path.concat(b.label === undefined ? '' : String(b.label));
      if (b.next) walk(b.next, prob * pb, next);
      else {
        const consequence = b.consequence === undefined ? 0 : b.consequence;
        const frequency = initiatingFrequency * prob * pb;
        sequences.push({ path: next, probability: prob * pb, frequency, outcome: b.outcome === undefined ? next.join('/') : b.outcome, consequence, risk: frequency * consequence });
      }
    }
  };
  walk(tree, 1, []);
  const outcomes = {};
  let totalRisk = 0;
  for (const s of sequences) {
    const o = outcomes[s.outcome] || (outcomes[s.outcome] = { frequency: 0, probability: 0, risk: 0 });
    o.frequency += s.frequency;
    o.probability += s.probability;
    o.risk += s.risk;
    totalRisk += s.risk;
  }
  return { sequences, outcomes, totalRisk };
}

/**
 * FMEA / FMECA table helper.
 *  Risk priority number RPN = severity × occurrence × detection (each ranked 1–10) — IEC 60812 / AIAG practice.
 *  Failure-mode criticality number C_m = β·α·λ_p·t (conditional probability of the effect, failure-mode ratio,
 *  part failure rate, operating time) — MIL-STD-1629A Task 102; item criticality is the sum over its modes.
 * Rows may carry either set of fields or both.
 * @param {{item:string, mode:string, severity?:number, occurrence?:number, detection?:number, failureRate?:number,
 *   modeRatio?:number, effectProbability?:number, operatingTime?:number}[]} rows failureRate [1/h] with
 *   operatingTime [h] (or any consistent pair)
 * @returns {{rows: object[], byRpn: object[], byCriticality: object[], itemCriticality: Object<string, number>}}
 *   rows gain `rpn`, `criticality`, `rpnRank` and `criticalityRank` (1 = worst; null when not computable)
 */
export function fmeca(rows) {
  need(Array.isArray(rows), 'rows must be an array');
  const out = rows.map((r, i) => {
    let rpn = null;
    if (r.severity !== undefined && r.occurrence !== undefined && r.detection !== undefined) {
      for (const k of ['severity', 'occurrence', 'detection']) need(r[k] >= 1 && r[k] <= 10, `rows[${i}].${k} must be in 1..10`);
      rpn = r.severity * r.occurrence * r.detection;
    }
    let criticality = null;
    if (r.failureRate !== undefined && r.operatingTime !== undefined) {
      const alpha = r.modeRatio === undefined ? 1 : r.modeRatio;
      const beta = r.effectProbability === undefined ? 1 : r.effectProbability;
      need(r.failureRate >= 0 && r.operatingTime >= 0, `rows[${i}] failureRate and operatingTime must be non-negative`);
      need(alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1, `rows[${i}] modeRatio and effectProbability must be in [0, 1]`);
      criticality = beta * alpha * r.failureRate * r.operatingTime;
    }
    return { ...r, rpn, criticality, rpnRank: null, criticalityRank: null };
  });
  const byRpn = out.filter((r) => r.rpn !== null).sort((a, b) => b.rpn - a.rpn || (b.severity - a.severity));
  byRpn.forEach((r, i) => { r.rpnRank = i + 1; });
  const byCriticality = out.filter((r) => r.criticality !== null).sort((a, b) => b.criticality - a.criticality);
  byCriticality.forEach((r, i) => { r.criticalityRank = i + 1; });
  const itemCriticality = {};
  for (const r of byCriticality) itemCriticality[r.item] = (itemCriticality[r.item] || 0) + r.criticality;
  return { rows: out, byRpn, byCriticality, itemCriticality };
}
