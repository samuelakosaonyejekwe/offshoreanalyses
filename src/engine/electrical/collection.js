/**
 * Offshore wind-farm electrical collection and export system: layouts, cable sizing, radial load flow, energy
 * losses, reliability and option comparison.
 *
 * References
 *  - Esau & Williams, "On teleprocessing system design, part II", IBM Syst. J. 5 (1966) 142: capacitated
 *    minimum-spanning-tree heuristic.
 *  - Kersting, Distribution System Modeling and Analysis (2012) ch. 10; Shirmohammadi et al., IEEE Trans. Power
 *    Syst. 3 (1988) 753: backward/forward sweep load flow for radial networks.
 *  - IEC 60287 (cable ratings through ampacity.js); IEC 60502-2 (8.0 mm XLPE insulation for 18/30 (36) kV).
 *
 * The cable catalogue built here is GENERIC: three-core copper XLPE wire-armoured array cables with plausible
 * layer thicknesses, rated by the IEC 60287 engine for the stated installation. It is a planning aid — replace
 * dimensions, ratings and all costs with manufacturer data for design. Default costs are indicative placeholders.
 * SI units; complex phasors { re, im }; voltages line-to-line unless stated.
 */

import { cableGeometry, rateCable } from './ampacity.js';
import { inductance, transformer as transformerModel, compensatedLink, hvdcLink } from './circuit.js';
import { cx, cadd, csub, cmul, cdiv, cconj, cabs, cscale } from './complex.js';

function pos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a positive finite number`);
}
function nonNeg(v, name) {
  if (!(v >= 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a non-negative finite number`);
}

/** Standard conductor sections used for the generic array-cable catalogue, m². */
export const ARRAY_CABLE_SIZES = Object.freeze([95, 120, 150, 185, 240, 300, 400, 500, 630, 800].map((a) => a * 1e-6));

/**
 * Generic three-core array-cable construction for ampacity.cableGeometry.
 * Layers: stranded compacted conductor, 0.8 mm conductor screen, XLPE insulation (8.0 mm at 33 kV per
 * IEC 60502-2; 9.0 mm assumed at 66 kV — designs range roughly 9–12 mm), 1.0 mm insulation screen and tapes,
 * individual copper screen of 16 mm² (33 kV) or 25 mm² (66 kV) per core, 2 mm PE core jacket, 2 mm bedding,
 * one layer of 5 mm galvanised steel wires, 4 mm yarn serving. All of these are assumptions of a generic design.
 * @param {object} o
 * @param {number} o.voltage nominal system voltage, V (33e3 or 66e3; other values use the nearer design)
 * @param {number} o.area conductor cross-section, m² @param {'Cu'|'Al'} [o.material='Cu']
 * @returns {object} cable description
 */
export function genericArrayCable({ voltage, area, material = 'Cu' }) {
  pos(voltage, 'voltage'); pos(area, 'area');
  const hv = voltage > 50e3;
  const screenArea = hv ? 25e-6 : 16e-6;
  return {
    cores: 3,
    conductor: { material, area, fillFactor: 0.9, construction: 'roundStranded' },
    conductorScreenThickness: 0.8e-3,
    insulation: { material: 'xlpe', thickness: hv ? 9.0e-3 : 8.0e-3, tanDelta: hv ? 0.001 : 0.004 },
    insulationScreenThickness: 1.0e-3,
    // The wire/tape screen is represented as an equivalent tube with the stated copper cross-section.
    sheath: { material: 'Cu', thickness: 0.3e-3, R20: 1.7241e-8 / screenArea },
    coreJacket: { thickness: 2.0e-3, rhoT: 3.5 },
    bedding: { thickness: 2.0e-3, rhoT: 6.0 },
    armour: { material: 'steel', wireDiameter: 5.0e-3 },
    serving: { thickness: 4.0e-3, rhoT: 6.0 },
  };
}

/**
 * Catalogue of generic array cables for a voltage level: current rating for the given installation (IEC 60287
 * through ampacity.rateCable), effective AC resistance at maximum temperature R·(1 + λ1 + λ2) (so that load-flow
 * losses include screen and armour losses), reactance, capacitance and cost.
 * Costs: `cost(area, voltage)` returns the supplied-and-installed cost per metre; the default is an indicative
 * placeholder (100 + 0.40 per mm² of copper, ×1.15 at 66 kV, currency units per metre) to be replaced.
 * @param {object} o
 * @param {number} o.voltage V @param {number} [o.f=50] Hz
 * @param {number[]} [o.areas=ARRAY_CABLE_SIZES] m² @param {'Cu'|'Al'} [o.material='Cu']
 * @param {object} [o.installation] ampacity installation (default buried 1.0 m, 0.7 K·m/W, 15 °C seabed)
 * @param {number} [o.thetaMaxC=90] @param {number} [o.armourLossReduction=1]
 * @param {(area:number, voltage:number)=>number} [o.cost] cost per metre
 * @returns {{area:number, voltage:number, I:number, r:number, x:number, c:number, l:number, De:number, cost:number,
 *            generic:true}[]} A, Ω/m, Ω/m, F/m, H/m, m, currency/m
 */
export function genericCatalogue({ voltage, f = 50, areas = ARRAY_CABLE_SIZES, material = 'Cu', installation, thetaMaxC = 90,
  armourLossReduction = 1, cost }) {
  pos(voltage, 'voltage'); pos(f, 'f');
  const inst = installation ?? { type: 'buried', depth: 1.0, rhoSoil: 0.7, ambientC: 15 };
  const costFn = cost ?? ((a, u) => (100 + 0.40 * a * 1e6) * (u > 50e3 ? 1.15 : 1));
  return Array.from(areas, (area) => {
    const g = cableGeometry(genericArrayCable({ voltage, area, material }));
    const r = rateCable({ cable: g, installation: inst, system: { f, U: voltage, thetaMaxC, armourLossReduction } });
    const l = inductance({ dc: g.dc, gmd: g.s });
    return { area, voltage, I: r.I, r: r.R * (1 + r.lambda1 + r.lambda2), x: 2 * Math.PI * f * l, c: g.C, l, De: g.De,
      cost: costFn(area, voltage), generic: true };
  });
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Layout                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Cable length of a route segment between two structures: plan distance × seabed slack factor plus the cable in
 * the J-tubes / hang-offs at each end.
 * @param {object} o
 * @param {{x:number, y:number}} o.a @param {{x:number, y:number}} o.b plan coordinates, m
 * @param {number} [o.slackFactor=1.03] ratio of laid length to straight plan distance (seabed slope, micro-routing)
 * @param {number} [o.endLength=0] cable length at each end (J-tube rise + deck routing + burial transition), m
 * @returns {number} cable length, m
 */
export function routeLength({ a, b, slackFactor = 1.03, endLength = 0 }) {
  if (!(slackFactor >= 1)) throw new RangeError('slackFactor: must be ≥ 1');
  nonNeg(endLength, 'endLength');
  return Math.hypot(a.x - b.x, a.y - b.y) * slackFactor + 2 * endLength;
}

function segCross(p1, p2, p3, p4) {
  // Proper intersection of segments p1p2 and p3p4 (shared end points do not count).
  const o = (a, b, c) => (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  const d1 = o(p3, p4, p1); const d2 = o(p3, p4, p2); const d3 = o(p1, p2, p3); const d4 = o(p1, p2, p4);
  return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}

/**
 * Network of a collection system from a parent list. Node 0 is the offshore substation; node i ≥ 1 is turbine
 * i − 1. Each turbine has one upstream branch towards its parent.
 * @param {object} o
 * @param {{x:number, y:number}} o.substation @param {{x:number, y:number}[]} o.turbines m
 * @param {number[]} o.parent parent node of each turbine (0 = substation, k = turbine k − 1)
 * @param {number} [o.slackFactor=1.03] @param {number} [o.endLength=0] m
 * @returns {{nodes:{x:number, y:number}[], parent:Int32Array, length:Float64Array, downstream:Int32Array,
 *            order:Int32Array, totalLength:number}} parent/length/downstream are indexed by node (entry 0 unused);
 *          downstream = number of turbines fed through the node's upstream branch; order = nodes sorted from the
 *          substation outwards
 */
export function buildNetwork({ substation, turbines, parent, slackFactor = 1.03, endLength = 0 }) {
  const n = turbines.length;
  if (!Array.isArray(parent) && !ArrayBuffer.isView(parent)) throw new RangeError('parent: must be an array');
  if (parent.length !== n) throw new RangeError('parent: one entry per turbine required');
  const nodes = [substation, ...turbines];
  const par = new Int32Array(n + 1); par[0] = -1;
  const length = new Float64Array(n + 1);
  const children = Array.from({ length: n + 1 }, () => []);
  for (let i = 1; i <= n; i++) {
    const p = parent[i - 1];
    if (!(Number.isInteger(p) && p >= 0 && p <= n && p !== i)) throw new RangeError(`parent[${i - 1}]: invalid node index`);
    par[i] = p; children[p].push(i);
    length[i] = routeLength({ a: nodes[i], b: nodes[p], slackFactor, endLength });
  }
  const order = new Int32Array(n + 1); let cnt = 0;
  const queue = [0];
  while (queue.length) { const u = queue.shift(); order[cnt++] = u; for (const c of children[u]) queue.push(c); }
  if (cnt !== n + 1) throw new RangeError('parent: network is not a tree rooted at the substation (cycle or disconnected)');
  const downstream = new Int32Array(n + 1);
  for (let k = n; k >= 1; k--) { const u = order[k]; downstream[u] += 1; downstream[par[u]] += downstream[u]; }
  let totalLength = 0;
  for (let i = 1; i <= n; i++) totalLength += length[i];
  return { nodes, parent: par, length, downstream, order, totalLength };
}

/**
 * Parent list for radial strings given as ordered lists of turbine indices from the substation outwards.
 * @param {object} o
 * @param {number} o.nTurbines @param {number[][]} o.strings each string lists turbine indices (0-based) in
 *        connection order starting with the turbine nearest the substation
 * @returns {number[]} parent node per turbine (see buildNetwork)
 */
export function radialStrings({ nTurbines, strings }) {
  const parent = new Array(nTurbines).fill(-1);
  for (const s of strings) {
    for (let k = 0; k < s.length; k++) {
      const t = s[k];
      if (!(Number.isInteger(t) && t >= 0 && t < nTurbines)) throw new RangeError('strings: turbine index out of range');
      if (parent[t] !== -1) throw new RangeError(`strings: turbine ${t} appears more than once`);
      parent[t] = k === 0 ? 0 : s[k - 1] + 1;
    }
  }
  if (parent.includes(-1)) throw new RangeError('strings: every turbine must belong to a string');
  return parent;
}

/**
 * Ring (loop) layout: pairs of radial strings whose far ends are joined by a normally open link, so that either
 * string can be back-fed after a cable fault. Returns the normal (open-ring) radial network plus the closing
 * links; every cable of a ring must then be sized for the turbines of both strings (ringDownstream).
 * @param {object} o
 * @param {{x:number, y:number}} o.substation @param {{x:number, y:number}[]} o.turbines
 * @param {number[][]} o.strings radial strings as for radialStrings
 * @param {[number, number][]} o.rings pairs of string indices closed into a ring
 * @param {number} [o.slackFactor=1.03] @param {number} [o.endLength=0]
 * @returns {{network:object, links:{a:number, b:number, length:number}[], ringDownstream:Int32Array,
 *            totalLength:number}} links join the last turbines (node indices) of each pair; ringDownstream is the
 *          sizing turbine count per node branch; totalLength includes the links (m)
 */
export function ringLayout({ substation, turbines, strings, rings, slackFactor = 1.03, endLength = 0 }) {
  const parent = radialStrings({ nTurbines: turbines.length, strings });
  const network = buildNetwork({ substation, turbines, parent, slackFactor, endLength });
  const ringDownstream = Int32Array.from(network.downstream);
  const links = []; let extra = 0;
  for (const [i, j] of rings) {
    if (!strings[i] || !strings[j] || i === j) throw new RangeError('rings: invalid pair of string indices');
    const a = strings[i][strings[i].length - 1] + 1; const b = strings[j][strings[j].length - 1] + 1;
    const length = routeLength({ a: network.nodes[a], b: network.nodes[b], slackFactor, endLength });
    links.push({ a, b, length }); extra += length;
    const total = strings[i].length + strings[j].length;
    for (const t of [...strings[i], ...strings[j]]) ringDownstream[t + 1] = total;
  }
  return { network, links, ringDownstream, totalLength: network.totalLength + extra };
}

/**
 * Capacitated minimum-spanning-tree layout by the Esau–Williams heuristic. Starting from a star (every turbine
 * on its own feeder), the link i–j with the largest saving s = gate(component of i) − d(i, j) is accepted
 * repeatedly, provided the merged subtree has at most `maxPerString` turbines and (optionally) the new cable does
 * not cross an existing one. Branched strings can result; set `noBranching` to keep simple daisy chains
 * (connections only between the ends of chains).
 * Reference: Esau & Williams (1966). A heuristic: typically within a few per cent of the optimum cable length,
 * with no optimality guarantee; cable cost differences between sections are not considered.
 * @param {object} o
 * @param {{x:number, y:number}} o.substation @param {{x:number, y:number}[]} o.turbines m
 * @param {number} o.maxPerString largest number of turbines on one feeder
 * @param {boolean} [o.noCrossing=true] @param {boolean} [o.noBranching=false]
 * @returns {{parent:number[], strings:number[][], totalDistance:number, crossings:number}} parent list for
 *          buildNetwork, the turbines (0-based) of each feeder, total straight-line length (m) and the number of
 *          cable crossings remaining (0 when noCrossing succeeded)
 */
export function esauWilliams({ substation, turbines, maxPerString, noCrossing = true, noBranching = false }) {
  const n = turbines.length;
  if (!(n >= 1)) throw new RangeError('turbines: must be a non-empty array');
  if (!(Number.isInteger(maxPerString) && maxPerString >= 1)) throw new RangeError('maxPerString: must be a positive integer');
  const pts = [substation, ...turbines];
  const d = (i, j) => Math.hypot(pts[i].x - pts[j].x, pts[i].y - pts[j].y);
  // comp[i]: component id; gate[c]: node of component c connected to the root; adjacency of accepted edges.
  const comp = new Int32Array(n + 1); const size = new Int32Array(n + 1); const gate = new Int32Array(n + 1);
  const adj = Array.from({ length: n + 1 }, () => []);
  const edges = []; // accepted non-root edges [i, j]
  for (let i = 1; i <= n; i++) { comp[i] = i; size[i] = 1; gate[i] = i; }
  const crossesAny = (i, j) => {
    for (const [a, b] of edges) if (segCross(pts[i], pts[j], pts[a], pts[b])) return true;
    // Root links of the other components.
    for (let c = 1; c <= n; c++) {
      if (size[c] === 0 || c === comp[i]) continue; // the gate of i's component is removed by the merge
      if (segCross(pts[i], pts[j], pts[0], pts[gate[c]])) return true;
    }
    return false;
  };
  for (;;) {
    let best = 0; let bi = -1; let bj = -1;
    for (let i = 1; i <= n; i++) {
      const ci = comp[i];
      const gi = d(0, gate[ci]);
      // Without branching, i must be an end of its chain (its gate link disappears in the merge).
      if (noBranching && adj[i].length > 1) continue;
      for (let j = 1; j <= n; j++) {
        const cj = comp[j];
        if (ci === cj || size[ci] + size[cj] > maxPerString) continue;
        const s = gi - d(i, j);
        if (s <= best + 1e-12) continue;
        if (noBranching) {
          // j must be the free (far) end of the other chain.
          const endJ = size[cj] === 1 || (adj[j].length === 1 && j !== gate[cj]);
          if (!endJ) continue;
        }
        if (noCrossing && crossesAny(i, j)) continue;
        best = s; bi = i; bj = j;
      }
    }
    if (bi < 0) break;
    const ci = comp[bi]; const cj = comp[bj];
    edges.push([bi, bj]); adj[bi].push(bj); adj[bj].push(bi);
    for (let k = 1; k <= n; k++) if (comp[k] === ci) comp[k] = cj;
    size[cj] += size[ci]; size[ci] = 0;
  }
  // Orient each component from its gate.
  const parent = new Array(n).fill(-1);
  const strings = [];
  for (let c = 1; c <= n; c++) {
    if (size[c] === 0) continue;
    const g = gate[c];
    parent[g - 1] = 0;
    const members = [g - 1];
    const stack = [g]; const seen = new Set([g]);
    while (stack.length) {
      const u = stack.pop();
      for (const v of adj[u]) if (!seen.has(v)) { seen.add(v); parent[v - 1] = u; members.push(v - 1); stack.push(v); }
    }
    strings.push(members);
  }
  let totalDistance = 0; const all = [];
  for (let i = 1; i <= n; i++) { totalDistance += d(i, parent[i - 1]); all.push([i, parent[i - 1]]); }
  let crossings = 0;
  for (let a = 0; a < all.length; a++) for (let b = a + 1; b < all.length; b++) {
    if (segCross(pts[all[a][0]], pts[all[a][1]], pts[all[b][0]], pts[all[b][1]])) crossings++;
  }
  return { parent, strings, totalDistance, crossings };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Sizing and load flow                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Select a cable section for every branch from the number of turbines it carries:
 * I = n·Pturbine/(√3·U·cosφ); the smallest catalogue entry with rating ≥ I is chosen. `maxTypes` limits the number
 * of different sections (largest needed plus the best smaller ones by length saved), as is usual for procurement.
 * @param {object} o
 * @param {object} o.network network from buildNetwork @param {object[]} o.catalogue from genericCatalogue or own data
 * @param {number} o.turbinePower rated power per turbine, W @param {number} o.voltage V
 * @param {number} [o.powerFactor=0.95] @param {Int32Array|number[]} [o.downstream] override of the turbine count
 *        per branch (e.g. ringDownstream) @param {number} [o.maxTypes] maximum number of distinct sections
 * @returns {{type:Int32Array, current:Float64Array, utilisation:Float64Array, quantities:{area:number, length:number,
 *            cost:number}[], totalCost:number, totalLength:number}} type = catalogue index per node branch
 *          (entry 0 unused); quantities per section (m, currency)
 */
export function sizeCables({ network, catalogue, turbinePower, voltage, powerFactor = 0.95, downstream, maxTypes }) {
  pos(turbinePower, 'turbinePower'); pos(voltage, 'voltage');
  if (!(powerFactor > 0 && powerFactor <= 1)) throw new RangeError('powerFactor: must be in (0, 1]');
  if (!Array.isArray(catalogue) || catalogue.length === 0) throw new RangeError('catalogue: must be a non-empty array');
  const cat = catalogue.map((c, i) => ({ ...c, index: i })).sort((a, b) => a.I - b.I);
  const n = network.parent.length - 1;
  const ds = downstream ?? network.downstream;
  const current = new Float64Array(n + 1);
  const need = new Int32Array(n + 1);
  for (let i = 1; i <= n; i++) {
    current[i] = (ds[i] * turbinePower) / (Math.sqrt(3) * voltage * powerFactor);
    const k = cat.findIndex((c) => c.I >= current[i]);
    if (k < 0) throw new RangeError(`catalogue: no cable can carry ${current[i].toFixed(0)} A (${ds[i]} turbines); reduce the string size`);
    need[i] = k;
  }
  let allowed = [...new Set(need.subarray(1))].sort((a, b) => a - b);
  if (maxTypes !== undefined && allowed.length > maxTypes) {
    if (!(maxTypes >= 1)) throw new RangeError('maxTypes: must be ≥ 1');
    // Keep the largest; drop the section whose removal costs least until the limit is met.
    const costWith = (set) => {
      let c = 0;
      for (let i = 1; i <= n; i++) c += cat[set.find((k) => k >= need[i])].cost * network.length[i];
      return c;
    };
    while (allowed.length > maxTypes) {
      let bestSet = null; let bestCost = Infinity;
      for (let r = 0; r < allowed.length - 1; r++) {
        const set = allowed.filter((_, idx) => idx !== r);
        const c = costWith(set);
        if (c < bestCost) { bestCost = c; bestSet = set; }
      }
      allowed = bestSet;
    }
  }
  const type = new Int32Array(n + 1); const utilisation = new Float64Array(n + 1);
  const q = new Map(); let totalCost = 0; let totalLength = 0;
  for (let i = 1; i <= n; i++) {
    const k = allowed.find((kk) => kk >= need[i]);
    const c = cat[k];
    type[i] = c.index; utilisation[i] = current[i] / c.I;
    const row = q.get(c.area) ?? { area: c.area, length: 0, cost: 0 };
    row.length += network.length[i]; row.cost += c.cost * network.length[i];
    q.set(c.area, row); totalCost += c.cost * network.length[i]; totalLength += network.length[i];
  }
  return { type, current, utilisation, quantities: [...q.values()].sort((a, b) => a.area - b.area), totalCost, totalLength };
}

/**
 * Backward/forward-sweep AC load flow of a radial collection network (balanced, per phase).
 * Each branch is a nominal-π section (series r + jx, half the capacitance at each end); turbines inject
 * S = P + jQ at constant power; the substation busbar (node 0) is the slack at the given voltage.
 * Backward sweep: branch current (towards the substation) J_b = I_inj(c) − jωC·V_c/2 terms + Σ children;
 * forward sweep: V_c = V_p + Z_b·J_b. Iterated until the largest voltage change is below `tol`.
 * Reference: Shirmohammadi et al. (1988); Kersting ch. 10.
 * @param {object} o
 * @param {object} o.network from buildNetwork
 * @param {{r:number, x:number, c:number}[]} o.catalogue per-length parameters (Ω/m, Ω/m, F/m)
 * @param {Int32Array|number[]} o.type catalogue index per node branch (from sizeCables)
 * @param {number} o.voltage slack line-to-line voltage, V @param {number} [o.f=50] Hz
 * @param {number|number[]} o.P active power injected per turbine, W (scalar or one per turbine)
 * @param {number|number[]} [o.Q=0] reactive power injected per turbine, var (positive = exported/over-excited)
 * @param {number} [o.tol=1e-10] relative voltage tolerance @param {number} [o.maxIter=100]
 * @returns {{V:{re:number, im:number}[], U:Float64Array, I:Float64Array, branchLoss:Float64Array, losses:number,
 *            Psubstation:number, Qsubstation:number, Pgenerated:number, lossFraction:number, Umax:number,
 *            iterations:number}} node phase voltages (V), line-to-line magnitudes (V), largest branch current of
 *          each node branch (A), three-phase loss per branch and total (W), power delivered to the substation
 *          busbar (W, var; Q > 0 = reactive power flowing into the busbar), generated power, loss fraction
 */
export function radialLoadFlow({ network, catalogue, type, voltage, f = 50, P, Q = 0, tol = 1e-10, maxIter = 100 }) {
  pos(voltage, 'voltage'); pos(f, 'f');
  const n = network.parent.length - 1;
  const w = 2 * Math.PI * f;
  const E = voltage / Math.sqrt(3);
  const Pn = (i) => (typeof P === 'number' ? P : P[i - 1]);
  const Qn = (i) => (typeof Q === 'number' ? Q : Q[i - 1]);
  const Z = new Array(n + 1); const Bh = new Float64Array(n + 1); // half susceptance of each node branch
  const Bnode = new Float64Array(n + 1);
  for (let i = 1; i <= n; i++) {
    const c = catalogue[type[i]];
    if (!c) throw new RangeError(`type[${i}]: no such catalogue entry`);
    const L = network.length[i];
    Z[i] = cx(c.r * L, c.x * L);
    Bh[i] = 0.5 * w * (c.c ?? 0) * L;
    Bnode[i] += Bh[i]; Bnode[network.parent[i]] += Bh[i];
  }
  const V = Array.from({ length: n + 1 }, () => cx(E, 0));
  const J = Array.from({ length: n + 1 }, () => cx(0, 0));
  let iterations = 0; let converged = false;
  for (; iterations < maxIter && !converged; iterations++) {
    // Backward sweep (leaves first).
    for (let i = 0; i <= n; i++) J[i] = cx(0, 0);
    for (let k = n; k >= 1; k--) {
      const u = network.order[k];
      // Injection of the turbine minus the charging current of the half-capacitances lumped at the node.
      const inj = cconj(cdiv(cx(Pn(u) / 3, Qn(u) / 3), V[u]));
      const ich = cmul(cx(0, Bnode[u]), V[u]);
      J[u] = cadd(J[u], csub(inj, ich));
      J[network.parent[u]] = cadd(J[network.parent[u]], J[u]);
    }
    // Forward sweep.
    let dv = 0;
    for (let k = 1; k <= n; k++) {
      const u = network.order[k];
      const vNew = cadd(V[network.parent[u]], cmul(Z[u], J[u]));
      dv = Math.max(dv, cabs(csub(vNew, V[u])));
      V[u] = vNew;
    }
    converged = dv < tol * E;
  }
  if (!converged) throw new Error('radialLoadFlow: no convergence (voltage collapse or inconsistent data)');
  const U = new Float64Array(n + 1); const I = new Float64Array(n + 1); const branchLoss = new Float64Array(n + 1);
  let losses = 0; let Pgenerated = 0; let Umax = 0;
  for (let i = 0; i <= n; i++) { U[i] = cabs(V[i]) * Math.sqrt(3); if (U[i] > Umax) Umax = U[i]; }
  for (let i = 1; i <= n; i++) {
    const Iser = cabs(J[i]);
    // Current at the two cable ends includes the local half charging current.
    const iChild = cabs(cadd(J[i], cmul(cx(0, Bh[i]), V[i])));
    const iParent = cabs(csub(J[i], cmul(cx(0, Bh[i]), V[network.parent[i]])));
    I[i] = Math.max(Iser, iChild, iParent);
    branchLoss[i] = 3 * Iser * Iser * Z[i].re;
    losses += branchLoss[i]; Pgenerated += Pn(i);
  }
  // Power into the substation busbar: sum of branch currents arriving minus the charging of the halves at node 0.
  const I0 = csub(J[0], cmul(cx(0, Bnode[0]), V[0]));
  const S0 = cscale(cmul(V[0], cconj(I0)), 3);
  return { V, U, I, branchLoss, losses, Psubstation: S0.re, Qsubstation: S0.im, Pgenerated,
    lossFraction: Pgenerated > 0 ? losses / Pgenerated : 0, Umax, iterations };
}

/**
 * Annual energy loss of the collection system from a power-duration distribution: the load flow is solved at
 * each production level and the losses are weighted by the probability of that level.
 * The distribution comes from the caller (e.g. Weibull wind speed × turbine power curve, binned).
 * @param {object} o arguments of radialLoadFlow without P, plus
 * @param {number} o.turbinePower rated power per turbine, W
 * @param {{fraction:number, probability:number}[]} o.distribution production level (fraction of rated) and its
 *        probability; probabilities should sum to 1 (they are normalised)
 * @param {number} [o.hours=8760] @param {number} [o.powerFactorTan=0] Q/P of the turbines
 * @returns {{energyLoss:number, energyGenerated:number, lossFraction:number, lossAtRated:number, levels:object[]}}
 *          Wh per period, –, W
 */
export function annualEnergyLoss({ turbinePower, distribution, hours = 8760, powerFactorTan = 0, ...lf }) {
  pos(turbinePower, 'turbinePower');
  if (!Array.isArray(distribution) || distribution.length === 0) throw new RangeError('distribution: must be a non-empty array');
  let sp = 0;
  for (const d of distribution) { nonNeg(d.probability, 'distribution[].probability'); nonNeg(d.fraction, 'distribution[].fraction'); sp += d.probability; }
  pos(sp, 'distribution: total probability');
  let eLoss = 0; let eGen = 0; const levels = [];
  for (const d of distribution) {
    const p = d.fraction * turbinePower;
    const r = p > 0 ? radialLoadFlow({ ...lf, P: p, Q: p * powerFactorTan }) : radialLoadFlow({ ...lf, P: 0, Q: 0 });
    const wgt = d.probability / sp;
    eLoss += wgt * r.losses * hours; eGen += wgt * r.Pgenerated * hours;
    levels.push({ fraction: d.fraction, probability: wgt, losses: r.losses, Qsubstation: r.Qsubstation, Umax: r.Umax });
  }
  const rated = radialLoadFlow({ ...lf, P: turbinePower, Q: turbinePower * powerFactorTan });
  return { energyLoss: eLoss, energyGenerated: eGen, lossFraction: eGen > 0 ? eLoss / eGen : 0, lossAtRated: rated.losses, levels };
}

/**
 * Expected energy not delivered because of array-cable failures (N−1, one fault at a time).
 * Radial: a fault on a branch disconnects every turbine downstream of it for the repair time,
 * E = Σ_b λ·L_b·MTTR·n_b·Pturbine·CF. Ring: after the switching time the turbines are back-fed, so only the
 * switching time counts (cables must be sized for it, see ringLayout).
 * @param {object} o
 * @param {object} o.network @param {number} o.turbinePower W @param {number} o.capacityFactor mean production / rated
 * @param {number} o.failureRate cable failures per metre per year (typical array-cable statistics are of the order
 *        0.001–0.01 per km per year; project data) @param {number} o.mttr mean time to repair, hours
 * @param {boolean} [o.ring=false] @param {number} [o.switchingTime=4] hours to reconfigure a ring
 * @returns {{energyNotDelivered:number, expectedFailures:number, unavailability:number}} Wh/year, failures/year and
 *          energy not delivered as a fraction of the annual production
 */
export function energyNotDelivered({ network, turbinePower, capacityFactor, failureRate, mttr, ring = false, switchingTime = 4 }) {
  pos(turbinePower, 'turbinePower'); nonNeg(failureRate, 'failureRate'); nonNeg(mttr, 'mttr'); nonNeg(switchingTime, 'switchingTime');
  if (!(capacityFactor >= 0 && capacityFactor <= 1)) throw new RangeError('capacityFactor: must be in [0, 1]');
  const n = network.parent.length - 1;
  const outage = ring ? Math.min(switchingTime, mttr) : mttr;
  let e = 0; let nf = 0;
  for (let i = 1; i <= n; i++) {
    const lam = failureRate * network.length[i];
    nf += lam;
    e += lam * outage * network.downstream[i] * turbinePower * capacityFactor;
  }
  const annual = n * turbinePower * capacityFactor * 8760;
  return { energyNotDelivered: e, expectedFailures: nf, unavailability: annual > 0 ? e / annual : 0 };
}

/**
 * Cable quantities by section for the cost model (kilometres per conductor size).
 * @param {{area:number, length:number}[]} quantities from sizeCables
 * @returns {{sizeMm2:number, km:number}[]}
 */
export function cableQuantities(quantities) {
  return quantities.map((q) => ({ sizeMm2: Math.round(q.area * 1e6), km: q.length / 1000 }));
}

/**
 * Complete collection-system design for one voltage level: layout (given parent list or Esau–Williams), sizing,
 * load flow at rated power, annual losses and reliability; one row of a comparison table.
 * @param {object} o
 * @param {{x:number, y:number}} o.substation @param {{x:number, y:number}[]} o.turbines
 * @param {number} o.turbinePower W @param {number} o.voltage V @param {number} [o.f=50]
 * @param {number[]} [o.parent] layout; omitted → Esau–Williams with the largest string the biggest cable allows
 * @param {object[]} [o.catalogue] cable catalogue (default generic for the voltage)
 * @param {number} [o.maxPerString] @param {number} [o.powerFactor=0.95] @param {number} [o.maxTypes=2]
 * @param {{fraction:number, probability:number}[]} [o.distribution] @param {number} [o.slackFactor=1.03]
 * @param {number} [o.endLength=0] @param {number} [o.failureRate] per m per year @param {number} [o.mttr] hours
 * @param {number} [o.capacityFactor]
 * @returns {object} { voltage, maxPerString, nStrings, parent, network, sizing, loadFlow, annual, reliability,
 *          cableKm: [{ sizeMm2, km }], totalCost }
 */
export function designCollection({ substation, turbines, turbinePower, voltage, f = 50, parent, catalogue, maxPerString,
  powerFactor = 0.95, maxTypes = 2, distribution, slackFactor = 1.03, endLength = 0, failureRate, mttr, capacityFactor }) {
  const cat = catalogue ?? genericCatalogue({ voltage, f });
  const Imax = Math.max(...cat.map((c) => c.I));
  const perTurbine = turbinePower / (Math.sqrt(3) * voltage * powerFactor);
  const cap = Math.floor(Imax / perTurbine + 1e-9);
  if (cap < 1) throw new RangeError('catalogue: the largest cable cannot carry one turbine at this voltage');
  const mps = Math.min(maxPerString ?? cap, cap);
  const par = parent ?? esauWilliams({ substation, turbines, maxPerString: mps }).parent;
  const network = buildNetwork({ substation, turbines, parent: par, slackFactor, endLength });
  const sizing = sizeCables({ network, catalogue: cat, turbinePower, voltage, powerFactor, maxTypes });
  const lf = { network, catalogue: cat, type: sizing.type, voltage, f };
  const loadFlow = radialLoadFlow({ ...lf, P: turbinePower });
  const annual = distribution ? annualEnergyLoss({ ...lf, turbinePower, distribution }) : null;
  const reliability = failureRate !== undefined && mttr !== undefined && capacityFactor !== undefined
    ? energyNotDelivered({ network, turbinePower, capacityFactor, failureRate, mttr }) : null;
  let nStrings = 0;
  for (let i = 1; i < network.parent.length; i++) if (network.parent[i] === 0) nStrings++;
  return { voltage, maxPerString: mps, nStrings, parent: par, network, sizing, loadFlow, annual, reliability,
    cableKm: cableQuantities(sizing.quantities), totalCost: sizing.totalCost };
}

/**
 * Comparison table of collection voltages (e.g. 33 kV against 66 kV): cable length, cost, strings, losses.
 * @param {object} o arguments of designCollection without `voltage`, `parent` and `catalogue`, plus
 * @param {number[]} [o.voltages=[33e3, 66e3]] @param {(voltage:number)=>object[]} [o.catalogueFor] own catalogues
 * @returns {{voltage:number, nStrings:number, maxPerString:number, cableLength:number, cableCost:number,
 *            lossAtRated:number, lossFractionAtRated:number, annualLossFraction:?number, Umax:number}[]}
 */
export function compareCollectionVoltages({ voltages = [33e3, 66e3], catalogueFor, ...rest }) {
  return voltages.map((voltage) => {
    const d = designCollection({ ...rest, voltage, catalogue: catalogueFor ? catalogueFor(voltage) : undefined });
    return { voltage, nStrings: d.nStrings, maxPerString: d.maxPerString, cableLength: d.sizing.totalLength, cableCost: d.totalCost,
      lossAtRated: d.loadFlow.losses, lossFractionAtRated: d.loadFlow.lossFraction,
      annualLossFraction: d.annual ? d.annual.lossFraction : null, Umax: d.loadFlow.Umax };
  });
}

/**
 * Offshore substation transformer and export link for a given collected power: transformer losses and
 * regulation (circuit.transformer), then either an HVAC cable with shunt compensation (circuit.compensatedLink)
 * or an HVDC link (circuit.hvdcLink). Rows of an AC/DC comparison table.
 * @param {object} o
 * @param {number} o.P collected active power at the substation busbar, W @param {number} [o.Q=0] var
 * @param {{Sr:number, uk:number, Pk:number, P0:number}} [o.transformer] per transformer; `nTransformers` in parallel
 * @param {number} [o.nTransformers=1]
 * @param {object} [o.ac] { r, l, c, f, length, U, compensation, nCircuits } export cable constants per circuit
 * @param {object} [o.dc] arguments of hvdcLink (configuration, U, r, length, converterLossFraction); the pole
 *        current follows from the power
 * @returns {{transformerLoss:number, options:{name:string, Pdelivered:number, losses:number, efficiency:number,
 *            ImaxCable:?number, detail:object}[]}} W
 */
export function exportSystem({ P, Q = 0, transformer, nTransformers = 1, ac, dc }) {
  pos(P, 'P');
  let transformerLoss = 0; let Pafter = P;
  if (transformer) {
    pos(nTransformers, 'nTransformers');
    const s = Math.hypot(P, Q) / nTransformers;
    const t = transformerModel({ ...transformer, loading: s / transformer.Sr, powerFactor: Math.min(1, P / Math.hypot(P, Q)) });
    transformerLoss = nTransformers * t.losses;
    Pafter = P - transformerLoss;
  }
  const options = [];
  if (ac) {
    const nc = ac.nCircuits ?? 1;
    const link = compensatedLink({ r: ac.r, l: ac.l, c: ac.c, g: ac.g ?? 0, f: ac.f ?? 50, length: ac.length, U: ac.U,
      P: Pafter / nc, Q: (ac.Q ?? 0) / nc, compensation: ac.compensation ?? {} });
    options.push({ name: 'HVAC', Pdelivered: link.Pr * nc, losses: P - link.Pr * nc, efficiency: (link.Pr * nc) / P,
      ImaxCable: link.ImaxCable, detail: link });
  }
  if (dc) {
    const upp = (dc.configuration === 'symmetricMonopole' || dc.configuration === 'bipole') ? 2 * dc.U : dc.U;
    const eta = dc.converterLossFraction ?? 0;
    const I = (Pafter * (1 - eta)) / upp;
    const link = hvdcLink({ ...dc, I });
    options.push({ name: 'HVDC', Pdelivered: link.Pout, losses: P - link.Pout, efficiency: link.Pout / P, ImaxCable: I, detail: link });
  }
  return { transformerLoss, options };
}
