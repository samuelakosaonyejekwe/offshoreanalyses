// Array layout tools for wind, tidal and wave farms: regular/staggered/boundary layouts inside a polygon,
// constraint checks, yield-maximising layout optimisation, spacing/orientation sweeps, capacity density and
// inter-array cable length (minimum spanning tree and a capacitated-string heuristic).
//
// Coordinates are plane metres (x east, y north). Polygons are arrays of [x, y] vertices, not closed.

import { makeRng } from '../math/rng.js';
import { farmAep } from './aep.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function finitePos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}
function checkPolygon(poly, name = 'polygon') {
  need(Array.isArray(poly) && poly.length >= 3, `${name} must have at least three vertices`);
}

/**
 * Signed area of a polygon by the shoelace formula (positive for counter-clockwise vertices).
 * @param {number[][]} polygon vertices [x, y] (m)
 * @returns {number} signed area (m²)
 */
export function polygonArea(polygon) {
  checkPolygon(polygon);
  let a = 0;
  for (let i = 0, n = polygon.length; i < n; i++) {
    const p = polygon[i];
    const q = polygon[(i + 1) % n];
    a += p[0] * q[1] - q[0] * p[1];
  }
  return a / 2;
}

/**
 * Point-in-polygon test by ray casting (boundary points count as inside within rounding).
 * @param {number[]} pt [x, y] (m)
 * @param {number[][]} polygon vertices (m)
 * @returns {boolean}
 */
export function pointInPolygon(pt, polygon) {
  let inside = false;
  const [x, y] = pt;
  for (let i = 0, n = polygon.length, j = n - 1; i < n; j = i++) {
    const [xi, yi] = polygon[i];
    const [xj, yj] = polygon[j];
    if (yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

function distToSegment(px, py, ax, ay, bx, by) {
  const dx = bx - ax;
  const dy = by - ay;
  const L2 = dx * dx + dy * dy;
  const t = L2 > 0 ? Math.max(0, Math.min(1, ((px - ax) * dx + (py - ay) * dy) / L2)) : 0;
  return Math.hypot(px - ax - t * dx, py - ay - t * dy);
}

/**
 * Shortest distance from a point to the boundary of a polygon.
 * @param {number[]} pt [x, y] (m)
 * @param {number[][]} polygon vertices (m)
 * @returns {number} distance (m)
 */
export function distanceToBoundary(pt, polygon) {
  let d = Infinity;
  for (let i = 0, n = polygon.length; i < n; i++) {
    const a = polygon[i];
    const b = polygon[(i + 1) % n];
    d = Math.min(d, distToSegment(pt[0], pt[1], a[0], a[1], b[0], b[1]));
  }
  return d;
}

/**
 * Constraint check of a layout: inside the site polygon (with a setback), outside exclusion zones and with the
 * minimum centre spacing respected.
 * @param {object} o
 * @param {number[][]} o.points turbine positions [x, y] (m)
 * @param {number[][]} [o.polygon] site boundary (m)
 * @param {number[][][]} [o.exclusions] exclusion polygons (m)
 * @param {number} [o.minSpacing=0] minimum centre-to-centre distance (m)
 * @param {number} [o.setback=0] minimum distance to the site boundary (m)
 * @returns {{feasible:boolean, outside:number, excluded:number, spacingViolations:number, minDistance:number,
 *   penalty:number}} counts of violations, smallest pair distance (m), total violation length (m)
 */
export function checkLayout({ points, polygon, exclusions = [], minSpacing = 0, setback = 0 }) {
  need(Array.isArray(points), 'points must be an array');
  let outside = 0;
  let excluded = 0;
  let viol = 0;
  let minD = Infinity;
  let penalty = 0;
  for (let i = 0; i < points.length; i++) {
    const p = points[i];
    if (polygon) {
      const inside = pointInPolygon(p, polygon);
      const d = setback > 0 || !inside ? distanceToBoundary(p, polygon) : Infinity;
      if (!inside) {
        outside++;
        penalty += d + setback;
      } else if (d < setback) {
        outside++;
        penalty += setback - d;
      }
    }
    for (const ex of exclusions) {
      if (pointInPolygon(p, ex)) {
        excluded++;
        penalty += distanceToBoundary(p, ex);
      }
    }
    for (let j = i + 1; j < points.length; j++) {
      const d = Math.hypot(p[0] - points[j][0], p[1] - points[j][1]);
      if (d < minD) minD = d;
      if (d < minSpacing - 1e-9) {
        viol++;
        penalty += minSpacing - d;
      }
    }
  }
  return { feasible: outside + excluded + viol === 0, outside, excluded, spacingViolations: viol, minDistance: minD, penalty };
}

/**
 * Regular or staggered grid layout clipped to a polygon. Rows run along the direction `rowAngleDeg` (measured
 * from +x, counter-clockwise); positions along a row are `sx·D` apart and rows `sy·D` apart; with `stagger`
 * every other row is shifted by half the in-row spacing.
 * To place rows perpendicular to a prevailing wind coming from compass direction θ, use rowAngleDeg = −θ
 * (see `rowAngleForWind`).
 * @param {object} o
 * @param {number[][]} o.polygon site boundary (m)
 * @param {number} o.D rotor diameter (m)
 * @param {number} o.sx in-row spacing (rotor diameters)
 * @param {number} o.sy row spacing (rotor diameters)
 * @param {number} [o.rowAngleDeg=0] row direction (deg from +x, counter-clockwise)
 * @param {boolean} [o.stagger=false]
 * @param {number[][][]} [o.exclusions] exclusion polygons (m)
 * @param {number} [o.setback=0] minimum distance to the boundary (m)
 * @param {number[]} [o.origin] grid origin [x, y] (m); default polygon centroid of vertices
 * @param {number} [o.maxCount=Infinity] keep at most this many positions (those nearest the origin are dropped last)
 * @returns {{points:number[][], count:number}}
 */
export function gridLayout({ polygon, D, sx, sy, rowAngleDeg = 0, stagger = false, exclusions = [], setback = 0, origin, maxCount = Infinity }) {
  checkPolygon(polygon);
  finitePos(D, 'D');
  finitePos(sx, 'sx');
  finitePos(sy, 'sy');
  let cx = 0;
  let cy = 0;
  for (const p of polygon) {
    cx += p[0] / polygon.length;
    cy += p[1] / polygon.length;
  }
  if (origin) [cx, cy] = origin;
  let rMax = 0;
  for (const p of polygon) rMax = Math.max(rMax, Math.hypot(p[0] - cx, p[1] - cy));
  const a = (rowAngleDeg * Math.PI) / 180;
  const ex = [Math.cos(a), Math.sin(a)];
  const ey = [-Math.sin(a), Math.cos(a)];
  const dxs = sx * D;
  const dys = sy * D;
  const ni = Math.ceil(rMax / dxs) + 1;
  const nj = Math.ceil(rMax / dys) + 1;
  need((2 * ni + 1) * (2 * nj + 1) < 4e6, 'grid is too fine for the site extent');
  const points = [];
  for (let j = -nj; j <= nj; j++) {
    const shift = stagger && Math.abs(j) % 2 === 1 ? 0.5 * dxs : 0;
    for (let i = -ni; i <= ni; i++) {
      const u = i * dxs + shift;
      const v = j * dys;
      const p = [cx + u * ex[0] + v * ey[0], cy + u * ex[1] + v * ey[1]];
      if (!pointInPolygon(p, polygon)) continue;
      if (setback > 0 && distanceToBoundary(p, polygon) < setback) continue;
      if (exclusions.some((e) => pointInPolygon(p, e))) continue;
      points.push(p);
    }
  }
  if (points.length > maxCount) {
    points.sort((p, q) => Math.hypot(p[0] - cx, p[1] - cy) - Math.hypot(q[0] - cx, q[1] - cy));
    points.length = maxCount;
  }
  return { points, count: points.length };
}

/**
 * Row direction (deg from +x, counter-clockwise) that puts rows perpendicular to a wind coming from compass
 * direction `wdDeg` (clockwise from north).
 * @param {number} wdDeg wind direction, from (deg)
 * @returns {number} row angle (deg)
 */
export function rowAngleForWind(wdDeg) {
  need(Number.isFinite(wdDeg), 'wdDeg must be finite');
  return -wdDeg;
}

/**
 * Boundary-conforming layout: `nBoundary` positions equally spaced along the perimeter (pulled inward by the
 * setback toward the vertex centroid), optionally followed by a grid infill of the interior.
 * @param {object} o
 * @param {number[][]} o.polygon site boundary (m)
 * @param {number} o.nBoundary number of perimeter positions
 * @param {number} [o.setback=0] inward offset toward the centroid (m)
 * @param {object} [o.infill] options for `gridLayout` (D, sx, sy, …) used for the interior; positions closer than
 *   `minSpacing` to a perimeter position are dropped
 * @param {number} [o.minSpacing=0] minimum distance between positions (m)
 * @returns {{points:number[][], boundaryCount:number, count:number}}
 */
export function boundaryLayout({ polygon, nBoundary, setback = 0, infill, minSpacing = 0 }) {
  checkPolygon(polygon);
  need(Number.isInteger(nBoundary) && nBoundary >= 1, 'nBoundary must be an integer >= 1');
  const n = polygon.length;
  const seg = [];
  let per = 0;
  let cx = 0;
  let cy = 0;
  for (let i = 0; i < n; i++) {
    const a = polygon[i];
    const b = polygon[(i + 1) % n];
    const L = Math.hypot(b[0] - a[0], b[1] - a[1]);
    seg.push(L);
    per += L;
    cx += a[0] / n;
    cy += a[1] / n;
  }
  const points = [];
  for (let k = 0; k < nBoundary; k++) {
    let s = (k * per) / nBoundary;
    let i = 0;
    while (i < n - 1 && s > seg[i]) s -= seg[i++];
    const a = polygon[i];
    const b = polygon[(i + 1) % n];
    const t = seg[i] > 0 ? s / seg[i] : 0;
    let x = a[0] + t * (b[0] - a[0]);
    let y = a[1] + t * (b[1] - a[1]);
    const d = Math.hypot(cx - x, cy - y);
    if (setback > 0 && d > 0) {
      x += ((cx - x) / d) * setback;
      y += ((cy - y) / d) * setback;
    }
    points.push([x, y]);
  }
  const boundaryCount = points.length;
  if (infill) {
    const inner = gridLayout({ ...infill, polygon }).points;
    for (const p of inner) {
      if (points.every((q) => Math.hypot(p[0] - q[0], p[1] - q[1]) >= minSpacing)) points.push(p);
    }
  }
  return { points, boundaryCount, count: points.length };
}

/**
 * Build an annual-energy objective for layout optimisation from a turbine definition and a wind rose.
 * @param {object} o
 * @param {object} o.turbine turbine template without position: { D, hubHeight, curve } or { D, hubHeight, power, ct }
 * @param {object[]} o.rose wind rose (see `farmAep`)
 * @param {object} [o.aep] further options for `farmAep` (du, uMin, uMax, wake, …)
 * @returns {(points:number[][])=>number} annual energy (Wh) of a layout
 */
export function makeAepObjective({ turbine, rose, aep = {} }) {
  need(turbine && Number.isFinite(turbine.D), 'turbine must define D');
  return (points) => farmAep({ ...aep, rose, turbines: points.map((p) => ({ ...turbine, x: p[0], y: p[1] })) }).aepWh;
}

/**
 * Layout optimisation: maximise an objective (annual energy) over turbine coordinates subject to boundary,
 * exclusion and minimum-spacing constraints. Infeasible candidates are rejected (random and pattern search) or
 * penalised (particle swarm). All methods are deterministic for a given seed.
 *  - 'random': seeded random search — move one turbine by a random step (shrinking), keep if feasible and better.
 *  - 'pattern': greedy pattern search — each turbine tries ±step along 8 compass directions; the step halves
 *    when no move improves (Hooke–Jeeves-type coordinate polling).
 *  - 'pso': particle swarm on the 2N coordinates with a penalty on constraint violation; the initial layout is
 *    one of the particles, so the result is never worse than the start.
 * @param {object} o
 * @param {number[][]} o.initial feasible starting layout (m)
 * @param {(points:number[][])=>number} o.objective value to maximise (e.g. from `makeAepObjective`)
 * @param {number[][]} o.polygon site boundary (m)
 * @param {number} o.minSpacing minimum centre spacing (m)
 * @param {number[][][]} [o.exclusions]
 * @param {number} [o.setback=0] (m)
 * @param {'random'|'pattern'|'pso'} [o.method='pattern']
 * @param {number} [o.maxEvals=300] objective-evaluation budget
 * @param {number} [o.step] initial move length (m); default 2 × minSpacing
 * @param {number} [o.minStep] smallest step (m); default step/64
 * @param {number} [o.seed=1]
 * @param {number} [o.particles=12] swarm size (pso)
 * @returns {{points:number[][], value:number, initialValue:number, improvement:number,
 *   history:{evals:number, value:number}[], evals:number, feasible:boolean}} improvement = value/initialValue − 1
 */
export function optimiseLayout({
  initial,
  objective,
  polygon,
  minSpacing,
  exclusions = [],
  setback = 0,
  method = 'pattern',
  maxEvals = 300,
  step,
  minStep,
  seed = 1,
  particles = 12,
}) {
  need(Array.isArray(initial) && initial.length > 0, 'initial must be a non-empty layout');
  need(typeof objective === 'function', 'objective must be a function');
  checkPolygon(polygon);
  need(minSpacing >= 0, 'minSpacing must be >= 0');
  need(['random', 'pattern', 'pso'].includes(method), "method must be 'random', 'pattern' or 'pso'");
  const cons = (pts) => checkLayout({ points: pts, polygon, exclusions, minSpacing, setback });
  need(cons(initial).feasible, 'initial layout violates the constraints');
  const rng = makeRng(seed);
  const n = initial.length;
  let best = initial.map((p) => [p[0], p[1]]);
  let bestVal = objective(best);
  const initialValue = bestVal;
  let evals = 1;
  const history = [{ evals, value: bestVal }];
  let h = step ?? (minSpacing > 0 ? 2 * minSpacing : Math.sqrt(Math.abs(polygonArea(polygon))) / 10);
  const hMin = minStep ?? h / 64;

  const tryMove = (i, dx, dy) => {
    const cand = best.map((p) => [p[0], p[1]]);
    cand[i][0] += dx;
    cand[i][1] += dy;
    if (!cons(cand).feasible) return false;
    const v = objective(cand);
    evals++;
    if (v > bestVal) {
      best = cand;
      bestVal = v;
      history.push({ evals, value: v });
      return true;
    }
    return false;
  };

  if (method === 'random') {
    let fails = 0;
    let guard = 0;
    while (evals < maxEvals && guard++ < 50 * maxEvals) {
      const i = rng.int(n);
      const ang = 2 * Math.PI * rng.next();
      const r = h * rng.next();
      if (tryMove(i, r * Math.cos(ang), r * Math.sin(ang))) fails = 0;
      else if (++fails > 4 * n && h > hMin) {
        h /= 2;
        fails = 0;
      }
    }
  } else if (method === 'pattern') {
    const dirs = [[1, 0], [-1, 0], [0, 1], [0, -1], [0.7071, 0.7071], [-0.7071, 0.7071], [0.7071, -0.7071], [-0.7071, -0.7071]];
    while (evals < maxEvals && h >= hMin) {
      let improved = false;
      for (let i = 0; i < n && evals < maxEvals; i++) {
        for (const d of dirs) {
          if (evals >= maxEvals) break;
          if (tryMove(i, d[0] * h, d[1] * h)) {
            improved = true;
            break;
          }
        }
      }
      if (!improved) h /= 2;
    }
  } else {
    let xMin = Infinity;
    let xMax = -Infinity;
    let yMin = Infinity;
    let yMax = -Infinity;
    for (const p of polygon) {
      xMin = Math.min(xMin, p[0]);
      xMax = Math.max(xMax, p[0]);
      yMin = Math.min(yMin, p[1]);
      yMax = Math.max(yMax, p[1]);
    }
    const scale = Math.abs(initialValue) || 1;
    const span = Math.max(xMax - xMin, yMax - yMin);
    const fit = (pts) => {
      const c = cons(pts);
      evals++;
      // one full spacing of violation costs about 100 % of the initial yield
      const v = objective(pts) - (scale * c.penalty) / Math.max(minSpacing, span / 100);
      return { v, feasible: c.feasible };
    };
    const unpack = (z) => Array.from({ length: n }, (_, i) => [z[2 * i], z[2 * i + 1]]);
    const P = [];
    for (let k = 0; k < particles; k++) {
      const z = new Float64Array(2 * n);
      const vel = new Float64Array(2 * n);
      for (let i = 0; i < n; i++) {
        const j = k === 0 ? 0 : 1;
        z[2 * i] = initial[i][0] + j * h * (rng.next() - 0.5);
        z[2 * i + 1] = initial[i][1] + j * h * (rng.next() - 0.5);
        vel[2 * i] = 0.2 * h * (rng.next() - 0.5);
        vel[2 * i + 1] = 0.2 * h * (rng.next() - 0.5);
      }
      const f = k === 0 ? { v: initialValue, feasible: true } : fit(unpack(z));
      P.push({ z, vel, pz: z.slice(), pv: f.v });
      if (f.feasible && f.v > bestVal) {
        bestVal = f.v;
        best = unpack(z);
        history.push({ evals, value: bestVal });
      }
    }
    let gz = P.reduce((a, b) => (b.pv > a.pv ? b : a)).pz.slice();
    let gv = Math.max(...P.map((p) => p.pv));
    while (evals < maxEvals) {
      for (const p of P) {
        if (evals >= maxEvals) break;
        for (let d = 0; d < 2 * n; d++) {
          p.vel[d] = 0.7 * p.vel[d] + 1.5 * rng.next() * (p.pz[d] - p.z[d]) + 1.5 * rng.next() * (gz[d] - p.z[d]);
          p.z[d] += p.vel[d];
        }
        const pts = unpack(p.z);
        const f = fit(pts);
        if (f.v > p.pv) {
          p.pv = f.v;
          p.pz = p.z.slice();
        }
        if (f.v > gv) {
          gv = f.v;
          gz = p.z.slice();
        }
        if (f.feasible && f.v > bestVal) {
          bestVal = f.v;
          best = pts;
          history.push({ evals, value: bestVal });
        }
      }
    }
  }
  return {
    points: best,
    value: bestVal,
    initialValue,
    improvement: initialValue !== 0 ? bestVal / initialValue - 1 : 0,
    history,
    evals,
    feasible: cons(best).feasible,
  };
}

/**
 * Sweep of grid spacing: for each (sx, sy) pair build a grid layout and evaluate the objective.
 * @param {object} o
 * @param {object} o.grid options for `gridLayout` except sx, sy
 * @param {number[]} o.sxValues in-row spacings (rotor diameters)
 * @param {number[]} [o.syValues] row spacings (rotor diameters); default equal to sxValues pairwise
 * @param {(points:number[][])=>number} o.objective layout value (e.g. annual energy, Wh)
 * @returns {{sx:number, sy:number, count:number, value:number, valuePerTurbine:number}[]}
 */
export function spacingSweep({ grid, sxValues, syValues, objective }) {
  need(Array.isArray(sxValues) && sxValues.length > 0, 'sxValues must be a non-empty array');
  need(typeof objective === 'function', 'objective must be a function');
  const out = [];
  const pairs = syValues ? sxValues.flatMap((a) => syValues.map((b) => [a, b])) : sxValues.map((a) => [a, a]);
  for (const [sx, sy] of pairs) {
    const lay = gridLayout({ ...grid, sx, sy });
    const value = lay.count ? objective(lay.points) : 0;
    out.push({ sx, sy, count: lay.count, value, valuePerTurbine: lay.count ? value / lay.count : 0 });
  }
  return out;
}

/**
 * Sweep of row orientation for a grid layout.
 * @param {object} o
 * @param {object} o.grid options for `gridLayout` except rowAngleDeg
 * @param {number[]} o.anglesDeg row angles (deg from +x, counter-clockwise)
 * @param {(points:number[][])=>number} o.objective layout value
 * @returns {{results:{rowAngleDeg:number, count:number, value:number, valuePerTurbine:number}[],
 *   best:{rowAngleDeg:number, count:number, value:number, valuePerTurbine:number}}} best by value per turbine
 */
export function orientationSweep({ grid, anglesDeg, objective }) {
  need(Array.isArray(anglesDeg) && anglesDeg.length > 0, 'anglesDeg must be a non-empty array');
  need(typeof objective === 'function', 'objective must be a function');
  const results = anglesDeg.map((rowAngleDeg) => {
    const lay = gridLayout({ ...grid, rowAngleDeg });
    const value = lay.count ? objective(lay.points) : 0;
    return { rowAngleDeg, count: lay.count, value, valuePerTurbine: lay.count ? value / lay.count : 0 };
  });
  const best = results.reduce((a, b) => (b.valuePerTurbine > a.valuePerTurbine ? b : a));
  return { results, best };
}

function convexHull(points) {
  const pts = points.map((p) => [p[0], p[1]]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  if (pts.length < 3) return pts;
  const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const lower = [];
  for (const p of pts) {
    while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) lower.pop();
    lower.push(p);
  }
  const upper = [];
  for (let i = pts.length - 1; i >= 0; i--) {
    const p = pts[i];
    while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) upper.pop();
    upper.push(p);
  }
  return lower.slice(0, -1).concat(upper.slice(0, -1));
}

/**
 * Capacity density of an array. Numerically W/m² equals MW/km². Two areas are reported: the site polygon (when
 * given) and the convex hull of the turbine positions.
 * @param {object} o
 * @param {number[][]} o.points positions (m)
 * @param {number} o.ratedPower rated power per turbine (W)
 * @param {number[][]} [o.polygon] site boundary (m)
 * @returns {{installed:number, hullArea:number, hullDensity:number, siteArea:number, siteDensity:number}}
 *   installed capacity (W), areas (m²), densities (W/m² ≡ MW/km²)
 */
export function capacityDensity({ points, ratedPower, polygon }) {
  need(Array.isArray(points) && points.length > 0, 'points must be a non-empty array');
  finitePos(ratedPower, 'ratedPower');
  const installed = points.length * ratedPower;
  const hull = convexHull(points);
  const hullArea = hull.length >= 3 ? Math.abs(polygonArea(hull)) : 0;
  const siteArea = polygon ? Math.abs(polygonArea(polygon)) : NaN;
  return {
    installed,
    hullArea,
    hullDensity: hullArea > 0 ? installed / hullArea : NaN,
    siteArea,
    siteDensity: siteArea > 0 ? installed / siteArea : NaN,
  };
}

/**
 * Minimum spanning tree over a set of points by Prim's algorithm (dense O(n²) form, Euclidean lengths).
 * @param {number[][]} points [x, y] (m)
 * @param {number} [root=0] index of the starting node
 * @returns {{edges:{from:number, to:number, length:number}[], length:number}} tree edges (parent → child), total (m)
 */
export function minimumSpanningTree(points, root = 0) {
  need(Array.isArray(points) && points.length > 0, 'points must be a non-empty array');
  const n = points.length;
  need(Number.isInteger(root) && root >= 0 && root < n, 'root must be a valid point index');
  const inTree = new Uint8Array(n);
  const dist = new Float64Array(n).fill(Infinity);
  const parent = new Int32Array(n).fill(-1);
  dist[root] = 0;
  const edges = [];
  let total = 0;
  for (let it = 0; it < n; it++) {
    let u = -1;
    for (let i = 0; i < n; i++) if (!inTree[i] && (u < 0 || dist[i] < dist[u])) u = i;
    inTree[u] = 1;
    if (parent[u] >= 0) {
      edges.push({ from: parent[u], to: u, length: dist[u] });
      total += dist[u];
    }
    for (let v = 0; v < n; v++) {
      if (inTree[v]) continue;
      const d = Math.hypot(points[u][0] - points[v][0], points[u][1] - points[v][1]);
      if (d < dist[v]) {
        dist[v] = d;
        parent[v] = u;
      }
    }
  }
  return { edges, length: total };
}

/**
 * Inter-array cable estimate: radial strings from a substation with at most `maxPerString` turbines per string,
 * by the Esau–Williams savings heuristic for the capacitated minimum spanning tree (Esau & Williams 1966).
 * Starting from a star (every turbine connected to the substation), the connection with the largest saving
 * (own gate length − link length) is accepted while the merged string respects the capacity. Without a capacity
 * limit the unconstrained minimum spanning tree (Prim) is returned. Lengths are straight-line plan lengths;
 * apply `routeFactor` for burial slack, J-tube/riser length and routing.
 * @param {object} o
 * @param {number[][]} o.turbines turbine positions [x, y] (m)
 * @param {number[]} o.substation substation position [x, y] (m)
 * @param {number} [o.maxPerString=Infinity] maximum turbines per string
 * @param {number} [o.routeFactor=1] multiplier on plan length (-)
 * @returns {{edges:{from:number, to:number, length:number}[], length:number, strings:number[][],
 *   mstLength:number}} edges use turbine indices, −1 = substation; total length (m) incl. routeFactor; strings as
 *   lists of turbine indices; unconstrained MST plan length (m) for comparison
 */
export function cableLayout({ turbines, substation, maxPerString = Infinity, routeFactor = 1 }) {
  need(Array.isArray(turbines) && turbines.length > 0, 'turbines must be a non-empty array');
  need(Array.isArray(substation) && substation.length >= 2, 'substation must be [x, y]');
  need(maxPerString >= 1, 'maxPerString must be >= 1');
  finitePos(routeFactor, 'routeFactor');
  const n = turbines.length;
  const all = [substation, ...turbines];
  const mst = minimumSpanningTree(all, 0);
  const d = (i, j) => Math.hypot(turbines[i][0] - turbines[j][0], turbines[i][1] - turbines[j][1]);
  const d0 = (i) => Math.hypot(turbines[i][0] - substation[0], turbines[i][1] - substation[1]);
  let edges;
  if (maxPerString >= n) {
    edges = mst.edges.map((e) => ({ from: e.from - 1, to: e.to - 1, length: e.length }));
  } else {
    const comp = Int32Array.from({ length: n }, (_, i) => i);
    const size = new Int32Array(n).fill(1);
    const gate = Float64Array.from({ length: n }, (_, i) => d0(i)); // gate length of each component (by id)
    const links = [];
    for (;;) {
      let bestSave = 1e-12;
      let bi = -1;
      let bj = -1;
      for (let i = 0; i < n; i++) {
        const ci = comp[i];
        for (let j = 0; j < n; j++) {
          const cj = comp[j];
          if (ci === cj || size[ci] + size[cj] > maxPerString) continue;
          const save = gate[ci] - d(i, j);
          if (save > bestSave) {
            bestSave = save;
            bi = i;
            bj = j;
          }
        }
      }
      if (bi < 0) break;
      const ci = comp[bi];
      const cj = comp[bj];
      links.push({ from: bj, to: bi, length: d(bi, bj) });
      size[cj] += size[ci];
      for (let i = 0; i < n; i++) if (comp[i] === ci) comp[i] = cj;
    }
    edges = links;
    const gateNode = new Map();
    for (let i = 0; i < n; i++) {
      const c = comp[i];
      if (!gateNode.has(c) || d0(i) < d0(gateNode.get(c))) gateNode.set(c, i);
    }
    // the component gate is its node nearest the substation (never longer than the gate tracked during merging)
    for (const g of gateNode.values()) edges.push({ from: -1, to: g, length: d0(g) });
    const strings = new Map();
    for (let i = 0; i < n; i++) {
      if (!strings.has(comp[i])) strings.set(comp[i], []);
      strings.get(comp[i]).push(i);
    }
    let len = 0;
    for (const e of edges) len += e.length;
    return { edges, length: len * routeFactor, strings: [...strings.values()], mstLength: mst.length };
  }
  // strings of the unconstrained tree: one per substation branch
  const adj = Array.from({ length: n }, () => []);
  const roots = [];
  for (const e of edges) {
    if (e.from < 0) roots.push(e.to);
    else adj[e.from].push(e.to);
  }
  const strings = roots.map((r) => {
    const list = [];
    const stack = [r];
    while (stack.length) {
      const u = stack.pop();
      list.push(u);
      for (const v of adj[u]) stack.push(v);
    }
    return list;
  });
  return { edges, length: mst.length * routeFactor, strings, mstLength: mst.length };
}
