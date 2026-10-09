import test from 'node:test';
import assert from 'node:assert/strict';
import {
  inextensibleCatenary, inextensibleCatenaryShape, inextensibleHorizontalStiffness, solveCatenary, catenaryProfile,
  lengthForFairleadTension, lengthForFairleadAngle, spanForFairleadTension, minLengthNoUplift,
} from '../../src/engine/lines/catenary.js';

const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);

test('inextensible closed form: s = sqrt(z(z + 2H/w)), x = (H/w) acosh(1 + w z/H), T = H + w z', () => {
  // Hand calculation: H = 500 kN, w = 1000 N/m, z = 100 m -> a = 500 m
  const r = inextensibleCatenary({ H: 5e5, w: 1000, z: 100 });
  assert.ok(rel(r.s, Math.sqrt(100 * (100 + 1000))) < 1e-14); // 331.662 m
  assert.ok(rel(r.s, 331.66247903554) < 1e-12);
  assert.ok(rel(r.x, 500 * Math.acosh(1.2)) < 1e-14);
  assert.ok(rel(r.x, 311.18125) < 1e-6);
  assert.equal(r.T, 6e5);
  // shape functions are consistent with the end values
  const sh = inextensibleCatenaryShape({ H: 5e5, w: 1000, x: r.x });
  assert.ok(rel(sh.z, 100) < 1e-12 && rel(sh.s, r.s) < 1e-12 && rel(sh.T, r.T) < 1e-12);
});

test('elastic solver reproduces the inextensible closed form when EA -> infinity (partly grounded line)', () => {
  const H = 5e5, w = 1000, z = 100, Lb = 150;
  const c = inextensibleCatenary({ H, w, z });
  const sol = solveCatenary({ X: c.x + Lb, Z: z, L: c.s + Lb, w, EA: 1e16, CB: 0 });
  assert.ok(rel(sol.H, H) < 1e-7);
  assert.ok(rel(sol.V, w * c.s) < 1e-7);
  assert.ok(rel(sol.T, H + w * z) < 1e-7);
  assert.ok(rel(sol.groundedLength, Lb) < 1e-7);
  assert.ok(rel(sol.touchdownX, Lb) < 1e-7);
  assert.ok(Math.abs(sol.anchorV) < 1e-6 && rel(sol.anchorH, H) < 1e-7);
  // Faltinsen (1990) Eq. 8.26 horizontal stiffness
  assert.ok(rel(sol.stiffness[0][0], inextensibleHorizontalStiffness({ H, w, z })) < 1e-5);
});

// Independent check: integrate the governing ODE of the extensible cable in unstretched arc length,
// dx/ds = (1 + T/EA) H/T, dz/ds = (1 + T/EA) V(s)/T, V(s) = VA + w s, with classical RK4.
function integrateCable({ H, VA, w, EA, L, n }) {
  const f = (s) => {
    const V = VA + w * s, T = Math.hypot(H, V), e = 1 + T / EA;
    return [e * H / T, e * V / T];
  };
  let x = 0, z = 0;
  const h = L / n;
  for (let i = 0; i < n; i++) {
    const s = i * h;
    const k1 = f(s), k2 = f(s + h / 2), k4 = f(s + h);
    x += h * (k1[0] + 4 * k2[0] + k4[0]) / 6; // RK4 for a pure quadrature right-hand side
    z += h * (k1[1] + 4 * k2[1] + k4[1]) / 6;
  }
  return { x, z };
}

test('elastic catenary (fully suspended) matches numerical integration of the cable ODE', () => {
  const p = { X: 800, Z: 300, L: 900, w: 1200, EA: 4e8 };
  const sol = solveCatenary({ ...p, seabed: false, nProfile: 5 });
  assert.equal(sol.groundedLength, 0);
  const end = integrateCable({ H: sol.H, VA: sol.anchorV, w: p.w, EA: p.EA, L: p.L, n: 4000 });
  assert.ok(Math.abs(end.x - p.X) < 1e-7 && Math.abs(end.z - p.Z) < 1e-7);
  // vertical equilibrium of the whole line and stretched length = integral of (1 + T/EA) ds
  assert.ok(rel(sol.V - sol.anchorV, p.w * p.L) < 1e-12);
  const mid = integrateCable({ H: sol.H, VA: sol.anchorV, w: p.w, EA: p.EA, L: p.L / 2, n: 2000 });
  assert.ok(Math.abs(sol.profile.x[2] - mid.x) < 1e-7 && Math.abs(sol.profile.z[2] - mid.z) < 1e-7);
  let ls = 0;
  const n = 4000;
  for (let i = 0; i < n; i++) {
    const s = (i + 0.5) * p.L / n;
    ls += (1 + Math.hypot(sol.H, sol.anchorV + p.w * s) / p.EA) * p.L / n;
  }
  assert.ok(rel(sol.stretchedLength, ls) < 1e-8);
});

test('seabed friction: anchor tension = H - CB w Lb, and zero once friction absorbs all of H', () => {
  const p = { X: 830, Z: 100, L: 850, w: 1000, EA: 5e8 };
  const s1 = solveCatenary({ ...p, CB: 0.2 });
  assert.ok(s1.groundedLength > 100 && s1.H > 2e5);
  assert.ok(rel(s1.anchorH, s1.H - 0.2 * p.w * s1.groundedLength) < 1e-12);
  assert.ok(s1.anchorH > 0);
  const s2 = solveCatenary({ ...p, CB: 5 });
  assert.ok(s2.H - 5 * p.w * s2.groundedLength < 0);
  assert.equal(s2.anchorH, 0);
  // grounded part: suspended weight is carried by the fairlead only
  assert.ok(rel(s1.V, p.w * s1.suspendedLength) < 1e-12);
  // profile: grounded nodes at z = 0, tension rises at CB w per metre on the seabed
  const pr = catenaryProfile(s1, 400);
  const i = 20;
  assert.ok(pr.s[i] < s1.groundedLength && pr.z[i] === 0);
  assert.ok(rel((pr.T[i] - pr.T[i - 1]) / (pr.s[i] - pr.s[i - 1]), 0.2 * p.w) < 1e-9);
  assert.ok(Math.abs(pr.x[399] - p.X) < 1e-7 && Math.abs(pr.z[399] - p.Z) < 1e-7);
});

test('tangent stiffness equals the finite-difference derivative of (H, V) with respect to (X, Z)', () => {
  for (const p of [
    { X: 780, Z: 100, L: 850, w: 1000, EA: 5e8, CB: 0.5 },
    { X: 800, Z: 300, L: 900, w: 1200, EA: 4e8, seabed: false },
    { X: 700, Z: 150, L: 700, w: 300, EA: 2e8 },
  ]) {
    const s0 = solveCatenary(p), d = 1e-3;
    const xp = solveCatenary({ ...p, X: p.X + d }), xm = solveCatenary({ ...p, X: p.X - d });
    const zp = solveCatenary({ ...p, Z: p.Z + d }), zm = solveCatenary({ ...p, Z: p.Z - d });
    const fd = [[(xp.H - xm.H) / (2 * d), (zp.H - zm.H) / (2 * d)], [(xp.V - xm.V) / (2 * d), (zp.V - zm.V) / (2 * d)]];
    const scale = Math.abs(s0.stiffness[0][0]) + Math.abs(s0.stiffness[1][1]);
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) assert.ok(Math.abs(fd[i][j] - s0.stiffness[i][j]) < 1e-5 * scale);
  }
});

test('taut, vertical, slack and weightless limits', () => {
  // weightless taut bar: T = EA (chord - L)/L
  const a = solveCatenary({ X: 300, Z: 400, L: 499, w: 0, EA: 1e8 });
  assert.ok(rel(a.T, 1e8 / 499) < 1e-12 && rel(a.H / a.T, 0.6) < 1e-12);
  // nearly weightless taut line tends to the same value
  const b = solveCatenary({ X: 300, Z: 400, L: 499, w: 1e-3, EA: 1e8, seabed: false });
  assert.ok(rel(b.T, a.T) < 1e-5);
  // vertical taut line: Z = L + (V L - w L^2/2)/EA
  const c = solveCatenary({ X: 0, Z: 101, L: 100, w: 500, EA: 1e7 });
  assert.ok(rel(c.V, 1e7 * 0.01 + 0.5 * 500 * 100) < 1e-9 && c.H < 1e-3);
  // slack line: hangs vertically, V = weight of the hanging length, remainder on the seabed, H = 0
  const d = solveCatenary({ X: 50, Z: 100, L: 400, w: 500, EA: 1e12 });
  assert.equal(d.slack, true);
  assert.equal(d.H, 0);
  assert.ok(rel(d.V, 500 * 100) < 1e-6 && rel(d.groundedLength, 300) < 1e-6);
});

test('buoyant line (w < 0) is the mirror image of the heavy line', () => {
  const heavy = solveCatenary({ X: 500, Z: -200, L: 600, w: 400, EA: 3e8, seabed: false, nProfile: 7 });
  const buoy = solveCatenary({ X: 500, Z: 200, L: 600, w: -400, EA: 3e8, nProfile: 7 });
  assert.ok(rel(buoy.H, heavy.H) < 1e-9 && rel(buoy.V, -heavy.V) < 1e-9);
  assert.ok(rel(buoy.anchorV, -heavy.anchorV) < 1e-9);
  for (let i = 0; i < 7; i++) assert.ok(Math.abs(buoy.profile.z[i] + heavy.profile.z[i]) < 1e-9);
  // buoyant line arches above the chord
  assert.ok(buoy.profile.z[3] > 100);
});

test('sloping seabed: global force equilibrium of the line (frictionless and with friction)', () => {
  for (const CB of [0, 0.3]) {
    const beta = 5 * Math.PI / 180, w = 1000, L = 900;
    const sol = solveCatenary({ X: 860, Z: 180, L, w, EA: 6e8, CB, slope: beta });
    const Lb = sol.groundedLength;
    assert.ok(Lb > 50 && Lb < L);
    const N = w * Lb * Math.cos(beta); // seabed normal reaction
    const F = CB * N; // fully mobilised friction, acting down-slope on the line
    // horizontal: H - TA cos(b) - N sin(b) - F cos(b) = 0 ; vertical: V - TA sin(b) + N cos(b) - F sin(b) - w L = 0
    assert.ok(Math.abs(sol.H - sol.anchorT * Math.cos(beta) - N * Math.sin(beta) - F * Math.cos(beta)) < 1e-6 * sol.H);
    assert.ok(Math.abs(sol.V - sol.anchorT * Math.sin(beta) + N * Math.cos(beta) - F * Math.sin(beta) - w * L) < 1e-6 * w * L);
    // touchdown tangent equals the slope, grounded nodes lie on the slope
    const pr = catenaryProfile(sol, 300);
    assert.ok(Math.abs(pr.z[10] - pr.x[10] * Math.tan(beta)) < 1e-9);
    assert.ok(Math.abs(sol.touchdownZ - sol.touchdownX * Math.tan(beta)) < 1e-9);
  }
  // zero slope reproduces the horizontal-seabed result
  const a = solveCatenary({ X: 800, Z: 180, L: 900, w: 1000, EA: 6e8, CB: 0.3, slope: 0 });
  const b = solveCatenary({ X: 800, Z: 180, L: 900, w: 1000, EA: 6e8, CB: 0.3 });
  assert.equal(a.H, b.H);
});

test('inverse problems: length / span for target tension and angle; minimum length for no uplift', () => {
  const base = { X: 780, Z: 100, w: 1000, EA: 5e8, CB: 0 };
  const r1 = lengthForFairleadTension({ ...base, T: 6e5 });
  assert.ok(rel(solveCatenary({ ...base, L: r1.L }).T, 6e5) < 1e-8);
  const r2 = spanForFairleadTension({ Z: 100, L: 850, w: 1000, EA: 5e8, T: 6e5 });
  assert.ok(rel(solveCatenary({ ...base, X: r2.X, L: 850 }).T, 6e5) < 1e-8);
  const ang = 40 * Math.PI / 180;
  const r3 = lengthForFairleadAngle({ ...base, angle: ang });
  assert.ok(rel(solveCatenary({ ...base, L: r3.L }).angle, ang) < 1e-8);
  // inextensible closed form: L = z sqrt(2T/(w z) - 1); T = 600 kN, w = 1 kN/m, z = 100 m -> 331.662 m
  const m0 = minLengthNoUplift({ Z: 100, w: 1000, T: 6e5 });
  assert.ok(rel(m0.L, 100 * Math.sqrt(11)) < 1e-12);
  // elastic: a line of exactly that length is tangent to the seabed at the anchor at the design tension
  const m1 = minLengthNoUplift({ Z: 100, w: 1000, T: 6e5, EA: 5e8 });
  const chk = solveCatenary({ X: m1.X, Z: 100, L: m1.L, w: 1000, EA: 5e8 });
  assert.ok(rel(chk.T, 6e5) < 1e-7 && Math.abs(chk.anchorV) < 1e-3 * 6e5 * 1e-3);
  assert.ok(m1.L < m0.L);
});

test('invalid input throws RangeError naming the argument', () => {
  assert.throws(() => solveCatenary({ X: -1, Z: 10, L: 20, w: 1, EA: 1 }), /X/);
  assert.throws(() => solveCatenary({ X: 10, Z: 10, L: 0, w: 1, EA: 1 }), /L/);
  assert.throws(() => inextensibleCatenary({ H: 0, w: 1, z: 1 }), RangeError);
  assert.throws(() => minLengthNoUplift({ Z: 100, w: 1000, T: 5e4 }), /T/);
});
