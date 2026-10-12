import test from 'node:test';
import assert from 'node:assert/strict';
import {
  polygonArea, pointInPolygon, distanceToBoundary, checkLayout, gridLayout, rowAngleForWind, boundaryLayout,
  makeAepObjective, optimiseLayout, spacingSweep, orientationSweep, capacityDensity, minimumSpanningTree,
  cableLayout,
} from '../../src/engine/aero/layout.js';
import { cubicPowerCurve } from '../../src/engine/aero/aep.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);
const square = (L) => [[0, 0], [L, 0], [L, L], [0, L]];

test('polygon geometry: area, containment and boundary distance', () => {
  close(polygonArea(square(2000)), 4e6, 1e-12);
  close(polygonArea(square(2000).reverse()), -4e6, 1e-12);
  // L-shaped polygon: 3×3 square minus a 2×2 corner = 5
  const L = [[0, 0], [3, 0], [3, 1], [1, 1], [1, 3], [0, 3]];
  close(polygonArea(L), 5, 1e-12);
  assert.equal(pointInPolygon([0.5, 2.5], L), true);
  assert.equal(pointInPolygon([2, 2], L), false);
  close(distanceToBoundary([500, 1000], square(2000)), 500, 1e-12);
  close(distanceToBoundary([2, 2], L), 1, 1e-12);
});

test('grid layout: count, spacing, stagger, exclusion zone and setback', () => {
  const g = gridLayout({ polygon: square(2000), D: 100, sx: 5, sy: 5, origin: [250, 250] });
  // nodes at 250, 750, 1250, 1750 in each direction
  assert.equal(g.count, 16);
  assert.equal(checkLayout({ points: g.points, polygon: square(2000), minSpacing: 500 }).feasible, true);
  close(checkLayout({ points: g.points }).minDistance, 500, 1e-9);
  const st = gridLayout({ polygon: square(2000), D: 100, sx: 5, sy: 5, origin: [250, 250], stagger: true });
  // staggered rows: nearest neighbours across rows are √(250² + 500²) apart
  const dmin = checkLayout({ points: st.points }).minDistance;
  close(dmin, 500, 1e-9);
  assert.ok(st.points.some((p) => Math.abs(p[0] - 500) < 1e-9));
  const ex = gridLayout({ polygon: square(2000), D: 100, sx: 5, sy: 5, origin: [250, 250], exclusions: [[[600, 600], [1400, 600], [1400, 1400], [600, 1400]]] });
  assert.equal(ex.count, 12);
  const sb = gridLayout({ polygon: square(2000), D: 100, sx: 5, sy: 5, origin: [250, 250], setback: 300 });
  assert.equal(sb.count, 4);
  const rot = gridLayout({ polygon: square(2000), D: 100, sx: 5, sy: 5, rowAngleDeg: 30 });
  close(checkLayout({ points: rot.points }).minDistance, 500, 1e-9);
  assert.equal(gridLayout({ polygon: square(2000), D: 100, sx: 5, sy: 5, origin: [250, 250], maxCount: 5 }).count, 5);
  assert.equal(rowAngleForWind(270), -270);
});

test('constraint check counts each kind of violation', () => {
  const r = checkLayout({
    points: [[100, 100], [150, 100], [2500, 100], [1000, 1000]],
    polygon: square(2000),
    exclusions: [[[900, 900], [1100, 900], [1100, 1100], [900, 1100]]],
    minSpacing: 200,
  });
  assert.equal(r.outside, 1);
  assert.equal(r.excluded, 1);
  assert.equal(r.spacingViolations, 1);
  assert.equal(r.feasible, false);
  close(r.minDistance, 50, 1e-12);
  // penalty: 500 m outside + 100 m inside the exclusion + 150 m of spacing shortfall
  close(r.penalty, 750, 1e-12);
});

test('boundary-conforming layout places positions on the perimeter at equal arc spacing', () => {
  const b = boundaryLayout({ polygon: square(1000), nBoundary: 8 });
  assert.equal(b.count, 8);
  for (const p of b.points) assert.ok(distanceToBoundary(p, square(1000)) < 1e-9);
  close(Math.hypot(b.points[1][0] - b.points[0][0], b.points[1][1] - b.points[0][1]), 500, 1e-9);
  const withInfill = boundaryLayout({ polygon: square(1000), nBoundary: 8, infill: { D: 100, sx: 4, sy: 4 }, minSpacing: 300 });
  assert.ok(withInfill.count > 8);
  assert.equal(checkLayout({ points: withInfill.points, minSpacing: 300 }).feasible, true);
});

test('capacity density: four 10 MW turbines on the corners of a 1 km square → 40 MW/km²', () => {
  const r = capacityDensity({ points: square(1000), ratedPower: 10e6, polygon: square(2000) });
  close(r.hullArea, 1e6, 1e-12);
  close(r.hullDensity, 40, 1e-12);
  close(r.siteDensity, 10, 1e-12);
});

test('minimum spanning tree (Prim): unit square has length 3; collinear points chain', () => {
  close(minimumSpanningTree(square(1)).length, 3, 1e-12);
  const line = [[0, 0], [5, 0], [1, 0], [3, 0]];
  const t = minimumSpanningTree(line);
  close(t.length, 5, 1e-12);
  assert.equal(t.edges.length, 3);
  // star-shaped set: MST connects every outer point to the centre
  const star = [[0, 0], [10, 0], [-10, 0], [0, 10], [0, -10]];
  close(minimumSpanningTree(star, 2).length, 40, 1e-12);
});

test('inter-array cabling: capacity-constrained strings never beat the MST and respect the string limit', () => {
  const turbines = [];
  for (let i = 0; i < 4; i++) for (let j = 0; j < 3; j++) turbines.push([800 * (i + 1), 800 * j]);
  const sub = [0, 800];
  const free = cableLayout({ turbines, substation: sub });
  close(free.length, free.mstLength, 1e-12);
  assert.equal(free.edges.length, turbines.length);
  for (const cap of [1, 2, 4, 6]) {
    const c = cableLayout({ turbines, substation: sub, maxPerString: cap });
    assert.ok(c.length >= free.mstLength - 1e-9, `cap ${cap}`);
    assert.ok(Math.max(...c.strings.map((s) => s.length)) <= cap);
    assert.equal(c.strings.reduce((s, x) => s + x.length, 0), turbines.length);
    assert.equal(c.edges.length, turbines.length);
    assert.equal(c.edges.filter((e) => e.from === -1).length, c.strings.length);
  }
  // one turbine per string is the star: Σ distance to substation
  const star = cableLayout({ turbines, substation: sub, maxPerString: 1 });
  close(star.length, turbines.reduce((s, p) => s + Math.hypot(p[0] - sub[0], p[1] - sub[1]), 0), 1e-12);
  // looser limits can only shorten the heuristic's result here, and the route factor scales it
  assert.ok(cableLayout({ turbines, substation: sub, maxPerString: 6 }).length <= star.length);
  close(cableLayout({ turbines, substation: sub, maxPerString: 3, routeFactor: 1.1 }).length, 1.1 * cableLayout({ turbines, substation: sub, maxPerString: 3 }).length, 1e-12);
});

const curve = cubicPowerCurve({ ratedPower: 3e6, D: 100, Cp: 0.45, step: 1 });
const northerly = [{ dirDeg: 0, frequency: 1, A: 9, k: 2 }];
const objective = makeAepObjective({ turbine: { D: 100, hubHeight: 90, curve }, rose: northerly, aep: { du: 2, uMin: 4, uMax: 24, wake: { model: 'jensen', turbulenceModel: 'none' } } });

test('layout optimisation improves the yield of a waked pair and keeps the constraints (all three methods)', () => {
  const polygon = square(1500);
  const initial = [[750, 300], [750, 1100]]; // aligned with the only wind direction
  const noWake = objective([[100, 700], [1400, 700]]);
  const start = objective(initial);
  assert.ok(start < noWake);
  for (const method of ['pattern', 'random', 'pso']) {
    const r = optimiseLayout({ initial, objective, polygon, minSpacing: 300, method, maxEvals: 120, seed: 7 });
    assert.ok(r.value > start * 1.02, `${method}: ${r.value / start}`);
    assert.ok(r.value <= noWake * (1 + 1e-9));
    assert.equal(r.feasible, true);
    assert.equal(checkLayout({ points: r.points, polygon, minSpacing: 300 }).feasible, true);
    for (let i = 1; i < r.history.length; i++) assert.ok(r.history[i].value >= r.history[i - 1].value);
    close(r.improvement, r.value / r.initialValue - 1, 1e-12);
    assert.ok(r.evals <= 120 + 12);
    // deterministic for a seed
    const again = optimiseLayout({ initial, objective, polygon, minSpacing: 300, method, maxEvals: 120, seed: 7 });
    assert.deepEqual(again.points, r.points);
  }
  assert.throws(() => optimiseLayout({ initial: [[750, 300], [750, 400]], objective, polygon, minSpacing: 300 }), RangeError);
});

test('spacing and orientation sweeps: wider spacing and rows across the wind give more energy per turbine', () => {
  const polygon = square(3000);
  const sp = spacingSweep({ grid: { polygon, D: 100, origin: [1500, 1500] }, sxValues: [4, 8], objective });
  assert.ok(sp[1].valuePerTurbine > sp[0].valuePerTurbine);
  assert.ok(sp[0].count > sp[1].count);
  // tight in-row spacing (3D) and wide row spacing (10D): rows must run east–west for a northerly wind
  const or = orientationSweep({ grid: { polygon, D: 100, sx: 3, sy: 10, origin: [1500, 1500] }, anglesDeg: [0, 90], objective });
  assert.equal(or.best.rowAngleDeg, 0);
  assert.ok(or.results[0].valuePerTurbine > or.results[1].valuePerTurbine);
  assert.equal(rowAngleForWind(0), -0);
});

test('capacity density of a single row or a single turbine: zero hull area and density with a note, never a non-finite value', () => {
  const row = [[0, 0], [1000, 0], [2000, 0], [3000, 0]];
  const r = capacityDensity({ points: row, ratedPower: 10e6 });
  assert.equal(r.installed, 40e6); assert.equal(r.hullArea, 0); assert.equal(r.hullDensity, 0);
  assert.equal(r.degenerate, true); assert.match(r.note, /one line/);
  // an oblique row (round-off sliver in the hull area) and a single turbine behave the same
  const ob = capacityDensity({ points: [[0, 0], [700.1, 300.3], [1400.2, 600.6], [2100.3, 900.9]], ratedPower: 10e6 });
  assert.equal(ob.hullDensity, 0); assert.equal(ob.degenerate, true);
  assert.equal(capacityDensity({ points: [[5, 5]], ratedPower: 10e6 }).hullDensity, 0);
  // the site polygon still gives a density for the row: 40 MW on 3 km × 0.5 km = 26.67 W/m²
  const s = capacityDensity({ points: row, ratedPower: 10e6, polygon: [[0, -250], [3000, -250], [3000, 250], [0, 250]] });
  close(s.siteDensity, 40e6 / 1.5e6, 1e-12); assert.equal(s.hullDensity, 0);
  // a proper layout: 4 turbines on a 1 km square → 40 W/m² (40 MW/km²), no note
  const q = capacityDensity({ points: square(1000), ratedPower: 10e6 });
  close(q.hullDensity, 40, 1e-12); assert.equal(q.degenerate, false); assert.equal(q.note, null);
});
