import test from 'node:test';
import assert from 'node:assert/strict';
import { solveCatenary } from '../../src/engine/lines/catenary.js';
import { solveMultiSegment, compileLine } from '../../src/engine/lines/multiseg.js';

const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);

test('single segment reproduces catenary.js (suspended, grounded with friction, slack)', () => {
  for (const p of [
    { X: 800, Z: 300, L: 900, w: 1200, EA: 4e8, CB: 0 },
    { X: 830, Z: 100, L: 850, w: 1000, EA: 5e8, CB: 0.4 },
    { X: 830, Z: 100, L: 850, w: 1000, EA: 5e8, CB: 5 },
    { X: 50, Z: 100, L: 400, w: 500, EA: 1e9, CB: 0 },
  ]) {
    const a = solveCatenary(p);
    const b = solveMultiSegment({ X: p.X, Z: p.Z, segments: [{ L: p.L, w: p.w, EA: p.EA, CB: p.CB }] });
    assert.equal(a.slack, b.slack);
    assert.ok(Math.abs(a.H - b.H) <= 1e-8 * a.T && Math.abs(a.V - b.V) <= 1e-8 * a.T);
    assert.ok(Math.abs(a.anchorH - b.anchorH) <= 1e-8 * a.T && Math.abs(a.anchorV - b.anchorV) <= 1e-8 * a.T);
    assert.ok(Math.abs(a.groundedLength - b.groundedLength) < 1e-6);
    assert.ok(Math.abs(a.stretchedLength - b.stretchedLength) < 1e-6);
    if (!a.slack) {
      for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
        assert.ok(Math.abs(a.stiffness[i][j] - b.stiffness[i][j]) < 1e-6 * Math.abs(a.stiffness[0][0]));
      }
    }
  }
});

test('a uniform line cut into three segments gives the single-line solution (touchdown inside a segment)', () => {
  const p = { X: 830, Z: 100, L: 850, w: 1000, EA: 5e8, CB: 0.4 };
  const a = solveCatenary(p);
  for (const cuts of [[100, 300, 450], [600, 200, 50], [a.groundedLength, 200, 850 - 200 - a.groundedLength]]) {
    const b = solveMultiSegment({ X: p.X, Z: p.Z, CB: 0.4, segments: cuts.map((L) => ({ L, w: p.w, EA: p.EA })), nProfile: 10 });
    assert.ok(rel(b.H, a.H) < 1e-8 && rel(b.V, a.V) < 1e-8 && rel(b.anchorH, a.anchorH) < 1e-8);
    assert.ok(Math.abs(b.groundedLength - a.groundedLength) < 1e-5);
    const m = b.profile.x.length;
    assert.ok(Math.abs(b.profile.x[m - 1] - p.X) < 1e-7 && Math.abs(b.profile.z[m - 1] - p.Z) < 1e-7);
    assert.ok(Math.abs(b.profile.x[0]) < 1e-6 && Math.abs(b.profile.z[0]) < 1e-9);
  }
});

// Independent integration of the cable ODE from the anchor, segment by segment, with jumps at point loads.
function integrate(sol, segments, P, nSub = 2000) {
  let x = 0, z = 0, V = sol.anchorV;
  const H = sol.H;
  segments.forEach((sg, i) => {
    if (i > 0) V += P[i] ?? 0;
    const h = sg.L / nSub;
    const f = (s) => { const v = V + sg.w * s, T = Math.hypot(H, v), e = 1 + T / sg.EA; return [e * H / T, e * v / T]; };
    for (let k = 0; k < nSub; k++) {
      const a = f(k * h), b = f((k + 0.5) * h), c = f((k + 1) * h);
      x += h * (a[0] + 4 * b[0] + c[0]) / 6; z += h * (a[1] + 4 * b[1] + c[1]) / 6;
    }
    V += sg.w * sg.L;
  });
  return { x, z, V };
}

test('chain-wire-chain with clump weight and buoy: junction force equilibrium and ODE integration', () => {
  const segments = [
    { name: 'ground chain', L: 300, w: 2000, EA: 8e8 },
    { name: 'wire', L: 500, w: 400, EA: 6e8 },
    { name: 'top chain', L: 150, w: 2000, EA: 8e8 },
  ];
  const points = [{ junction: 1, weight: 60e3 }, { junction: 2, weight: -40e3 }];
  const X = 820, Z = 350;
  const sol = solveMultiSegment({ X, Z, segments, points, seabed: false, nProfile: 4000 });
  // jump of the vertical tension component at each junction equals the point load
  for (const j of sol.junctions) assert.ok(Math.abs(j.Vabove - j.Vbelow - j.weight) < 1e-9 * sol.T);
  // the same from the computed geometry: V = H dz/dx just above and just below the junction
  const pr = sol.profile, n = 4001;
  const slope = (i, j) => (pr.z[j] - pr.z[i]) / (pr.x[j] - pr.x[i]);
  [[1, 60e3], [2, -40e3]].forEach(([k, Pk]) => {
    const below = slope(k * n - 2, k * n - 1), above = slope(k * n, k * n + 1);
    assert.ok(Math.abs(sol.H * (above - below) - Pk) < 2e-3 * Math.abs(Pk));
  });
  // global vertical equilibrium
  const W = 300 * 2000 + 500 * 400 + 150 * 2000 + 60e3 - 40e3;
  assert.ok(rel(sol.V - sol.anchorV, W) < 1e-10);
  // independent numerical integration from the anchor lands on the fairlead
  const end = integrate(sol, segments, { 1: 60e3, 2: -40e3 });
  assert.ok(Math.abs(end.x - X) < 1e-6 && Math.abs(end.z - Z) < 1e-6 && rel(end.V, sol.V) < 1e-10);
  // segment end tensions are continuous in H and consistent with the reported angles
  for (const s of sol.segments) assert.ok(rel(s.Ttop * Math.cos(s.angleTop), sol.H) < 1e-10);
});

test('grounded chain-polyester-chain: touchdown in the lowest segment, friction reduces the anchor load', () => {
  const segments = [
    { L: 400, w: 3000, EA: 9e8, CB: 0.05 },
    { L: 1200, w: 60, EA: 1.5e8 },
    { L: 100, w: 3000, EA: 9e8 },
  ];
  const sol = solveMultiSegment({ X: 1500, Z: 600, segments });
  const Lg = sol.groundedLength;
  assert.ok(Lg > 0 && Lg < 400);
  assert.ok(rel(sol.anchorH, sol.H - 0.05 * 3000 * Lg) < 1e-10);
  assert.equal(sol.anchorV, 0);
  assert.ok(rel(sol.V, 3000 * (400 - Lg) + 60 * 1200 + 3000 * 100) < 1e-10);
  assert.ok(sol.anchorH > 0);
  assert.equal(sol.sagClearance, null);
});

test('lazy-wave configuration: buoyant section forms a hog bend, sag bend stays clear of the seabed', () => {
  const segments = [
    { L: 250, w: 800, EA: 5e8 },
    { L: 120, w: -1600, EA: 5e8 },
    { L: 300, w: 800, EA: 5e8 },
  ];
  const X = 420, Z = 250;
  const sol = solveMultiSegment({ X, Z, segments, nProfile: 200 });
  const pr = sol.profile;
  // a local maximum (hog) inside the buoyant section and a local minimum (sag) in the upper section
  let hog = -Infinity, sag = Infinity;
  for (let i = 0; i < pr.z.length; i++) {
    if (pr.seg[i] === 1) hog = Math.max(hog, pr.z[i]);
    if (pr.seg[i] === 2) sag = Math.min(sag, pr.z[i]);
  }
  assert.ok(hog > sag + 1);
  assert.ok(sol.sagClearance > 0 && Math.abs(sol.sagClearance - sag) < 0.05);
  assert.ok(sol.groundedLength > 0 && sol.anchorV === 0);
  // ODE integration from touchdown: start at the end of the grounded length with V = 0
  const Lg = sol.groundedLength;
  const end = integrate({ H: sol.H, anchorV: 0 }, [{ ...segments[0], L: 250 - Lg }, segments[1], segments[2]], {});
  assert.ok(Math.abs(end.x + sol.touchdownX - X) < 1e-6 && Math.abs(end.z - Z) < 1e-6);
});

test('tangent stiffness matches finite differences; compiled line gives identical results', () => {
  const def = {
    segments: [{ L: 400, w: 3000, EA: 9e8, CB: 0.3 }, { L: 1200, w: 60, EA: 1.5e8 }, { L: 100, w: 3000, EA: 9e8 }],
    points: [{ junction: 2, weight: 20e3 }],
  };
  const line = compileLine(def);
  const f = (X, Z) => solveMultiSegment({ X, Z, line });
  const s0 = f(1500, 600), d = 1e-3;
  assert.equal(s0.H, solveMultiSegment({ X: 1500, Z: 600, ...def }).H);
  const fd = [
    [(f(1500 + d, 600).H - f(1500 - d, 600).H) / (2 * d), (f(1500, 600 + d).H - f(1500, 600 - d).H) / (2 * d)],
    [(f(1500 + d, 600).V - f(1500 - d, 600).V) / (2 * d), (f(1500, 600 + d).V - f(1500, 600 - d).V) / (2 * d)],
  ];
  for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
    assert.ok(Math.abs(fd[i][j] - s0.stiffness[i][j]) < 1e-5 * Math.abs(s0.stiffness[0][0]));
  }
});

test('invalid definitions throw RangeError', () => {
  assert.throws(() => solveMultiSegment({ X: 10, Z: 5, segments: [] }), /segments/);
  assert.throws(() => solveMultiSegment({ X: 10, Z: 5, segments: [{ L: 10, w: 1, EA: 1e6 }], points: [{ junction: 1, weight: 1 }] }), /junction/);
});
