import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { waveGreenFunction, hullPanels, radiationSetup, radiationSolve, haskindExcitation, haskindDamping, boxQuarterPanels, boxHydrodynamics, boxCoefficients, boxExcitation } from '../../src/engine/hydro/radiation.js';
import { addedMassBem } from '../../src/engine/hydro/bem.js';
import { box } from '../../src/mesh/generate.js';
import { readCapytaineCsv } from '../../src/interop/index.js';
import { hydroAt } from '../../src/engine/hydro/hydrodb.js';

const rho = 1025; const g = 9.80665;
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const fixture = (n) => readCapytaineCsv(fs.readFileSync(new URL(`../interop/fixtures/capytaine/${n}_hydro.csv`, import.meta.url), 'utf8'));

// Bessel function J0 from its integral (1/pi) ∫ cos(x sin t) dt, midpoint rule (spectrally accurate)
const J0 = (x) => { const n = 2 * Math.ceil(Math.abs(x)) + 40; let s = 0; for (let i = 0; i < n; i++) s += Math.cos(x * Math.sin((Math.PI * (i + 0.5)) / n)); return s / n; };
// principal value of ∫₀^∞ e^{uY} J0(uX)/(u - 1) du: the pole is removed on (0, 2), composite Simpson rule
function principalValue(X, Y) {
  const f = (u) => Math.exp(u * Y) * J0(u * X); const f1 = f(1);
  const simpson = (fn, a, b, n) => { const h = (b - a) / n; let s = fn(a) + fn(b); for (let i = 1; i < n; i++) s += fn(a + i * h) * (i % 2 ? 4 : 2); return (s * h) / 3; };
  const near = (u) => (Math.abs(u - 1) < 1e-9 ? (f(1 + 1e-5) - f(1 - 1e-5)) / 2e-5 : (f(u) - f1) / (u - 1));
  return simpson(near, 0, 2, 400) + simpson((u) => f(u) / (u - 1), 2, 2 + 30 / Math.abs(Y), 3000);
}

test('Wave Green function: F equals the principal-value integral, dF/dX its derivative, J0 the Bessel integral', () => {
  for (const [X, Y] of [[0.3, -0.4], [1, -0.5], [2.5, -1], [7.3, -0.6], [15.2, -2], [0, -1]]) {
    const r = { ...waveGreenFunction(X, Y) };
    const ref = principalValue(X, Y);
    assert.ok(Math.abs(r.F - ref) < 2e-5, `F(${X}, ${Y}) = ${r.F}, quadrature ${ref}`);
    assert.ok(Math.abs(r.J0 - J0(X)) < 1e-6, `J0(${X})`);
    if (X > 0) {
      const h = 1e-3; const dref = (principalValue(X + h, Y) - principalValue(X - h, Y)) / (2 * h);
      assert.ok(Math.abs(r.FX - dref) < 2e-4, `dF/dX(${X}, ${Y}) = ${r.FX}, quadrature ${dref}`);
    }
  }
});

test('Wave Green function: differential equation in Y, logarithmic singularity, far field, deep limit, range', () => {
  // dF/dY = F + 1/d (what makes G satisfy the free-surface condition), checked against a finite difference
  for (const [X, Y] of [[0.8, -0.3], [3, -1.2], [12, -0.2], [40, -3]]) {
    const r = { ...waveGreenFunction(X, Y) }; const h = 1e-4;
    const fd = (waveGreenFunction(X, Y + h).F - waveGreenFunction(X, Y - h).F) / (2 * h);
    assert.ok(Math.abs(r.FY - (r.F + 1 / Math.hypot(X, Y))) < 1e-12);
    assert.ok(Math.abs(fd - r.FY) < 2e-4, `dF/dY at (${X}, ${Y}): ${fd} vs ${r.FY}`);
  }
  // near the origin F -> -exp(Y) (ln((d + |Y|)/2) + gamma)
  const s = waveGreenFunction(0.003, -0.004);
  assert.ok(Math.abs(s.F + Math.exp(-0.004) * (Math.log((0.005 + 0.004) / 2) + 0.5772156649)) < 0.01);
  // far from the source on the free surface F -> -pi exp(Y) Y0(X): amplitude pi sqrt(2/(pi X)) exp(Y)
  let lo = Infinity; let hi = -Infinity;
  for (let X = 60; X < 66.3; X += 0.05) { const v = waveGreenFunction(X, -0.5).F; lo = Math.min(lo, v); hi = Math.max(hi, v); }
  const amp = 0.5 * (hi - lo);
  assert.ok(rel(amp, Math.PI * Math.sqrt(2 / (Math.PI * 63)) * Math.exp(-0.5)) < 0.03, `far-field amplitude ${amp}`);
  // deep below the surface only the image term is left: F -> -1/d (1 + |Y|/d² + …), continuous across |Y| = 16
  assert.ok(rel(waveGreenFunction(3, -15.99).F, waveGreenFunction(3, -16.01).F) < 2e-3);
  assert.ok(rel(waveGreenFunction(0, -40).F, -(1 / 40 + 1 / 1600 + 2 / 64000 + 6 / 2560000)) < 1e-9);
  assert.throws(() => waveGreenFunction(1, 0.1), /Y <= 0/);
  assert.throws(() => waveGreenFunction(101, -1), /tabulated range/);
});

// quarter (x, y >= 0) of a hemisphere of radius R below z = 0: nt rings from the pole, np sectors
function hemisphereQuarter(R, nt, np) {
  const P = (i, j) => { const th = (0.5 * Math.PI * i) / nt; const ph = (0.5 * Math.PI * j) / np; return [R * Math.sin(th) * Math.cos(ph), R * Math.sin(th) * Math.sin(ph), -R * Math.cos(th)]; };
  const q = [];
  for (let i = 0; i < nt; i++) for (let j = 0; j < np; j++) {
    // counter-clockwise seen from outside; the ring at the pole gives triangles (repeated last vertex)
    const v = i === 0 ? [P(0, 0), P(1, j + 1), P(1, j), P(1, j)] : [P(i, j), P(i, j + 1), P(i + 1, j + 1), P(i + 1, j)];
    for (const p of v) q.push(...p);
  }
  return Float64Array.from(q);
}

test('Floating hemisphere: zero-frequency added mass (Hulme 1982), Capytaine at Ka = 0.13, energy relation, symmetry', () => {
  const R = 5; const V = (2 / 3) * Math.PI * R ** 3;
  const panels = hullPanels(hemisphereQuarter(R, 9, 9));
  // normals of the panel set point into the fluid: ∫ z n_z dS = V
  let vol = 0;
  for (let i = 0; i < panels.n; i++) vol += 4 * panels.centroids[3 * i + 2] * panels.normals[3 * i + 2] * panels.areas[i];
  assert.ok(rel(vol, V) < 0.02, `volume ${vol}`);
  const setup = radiationSetup({ panels, symmetry: { x: true, y: true } });
  // Hulme (1982) J. Fluid Mech. 121, as quoted by Storti & D'Elía (2001) arXiv:math/0110302 eq. (38):
  // A33/(rho V) = 0.830951 and A11/(rho V) = 0.5 for Ka -> 0
  const low = radiationSolve({ setup, omega: Math.sqrt((g * 1e-4) / R), rho, g });
  assert.ok(rel(low.A[2][2] / (rho * V), 0.830951) < 0.015, `heave ${low.A[2][2] / (rho * V)}`);
  assert.ok(rel(low.A[0][0] / (rho * V), 0.5) < 0.015, `surge ${low.A[0][0] / (rho * V)}`);
  // equal up to the switch between closed-form and Gauss panel integrals at mirrored distances
  assert.ok(rel(low.A[1][1], low.A[0][0]) < 1e-5, 'sway = surge');
  // Capytaine 3.0.0 on a 2592-panel sphere mesh (tests/interop/fixtures), omega = 0.5 rad/s
  const db = fixture('hemisphere'); const w = 0.5; const ref = hydroAt(db, w, 0);
  const r = radiationSolve({ setup, omega: w, rho, g });
  for (const [i, tolA, tolB] of [[0, 0.02, 0.04], [2, 0.02, 0.03]]) {
    assert.ok(rel(r.A[i][i], ref.A[i][i]) < tolA, `A${i + 1}${i + 1} ${r.A[i][i]} vs ${ref.A[i][i]}`);
    assert.ok(rel(r.B[i][i], ref.B[i][i]) < tolB, `B${i + 1}${i + 1} ${r.B[i][i]} vs ${ref.B[i][i]}`);
  }
  const X = haskindExcitation({ setup, phi: r.phi, omega: w, headingDeg: 0, rho, g });
  for (const i of [0, 2]) {
    assert.ok(rel(Math.hypot(X.re[i], X.im[i]), Math.hypot(ref.F.re[i], ref.F.im[i])) < 0.015, `|X${i + 1}|`);
    assert.ok(Math.abs(Math.atan2(X.im[i], X.re[i]) - Math.atan2(ref.F.im[i], ref.F.re[i])) < 0.01, `phase X${i + 1}`);
  }
  // energy radiated = work done by the damping (Newman 1977 §6.19), body of revolution in deep water:
  // B33 = omega³ |X3|² / (2 rho g³), B11 = omega³ |X1|² / (4 rho g³)
  for (const wk of [0.5, 1.2, 1.8]) {
    const q = radiationSolve({ setup, omega: wk, rho, g, modes: [0, 2] });
    const x = haskindExcitation({ setup, phi: q.phi, omega: wk, headingDeg: 0, rho, g });
    assert.ok(q.B[2][2] > 0 && q.B[0][0] > 0);
    assert.ok(rel(q.B[2][2], (wk ** 3 * (x.re[2] ** 2 + x.im[2] ** 2)) / (2 * rho * g ** 3)) < 0.02, `heave energy balance at ${wk}`);
    assert.ok(rel(q.B[0][0], (wk ** 3 * (x.re[0] ** 2 + x.im[0] ** 2)) / (4 * rho * g ** 3)) < 0.02, `surge energy balance at ${wk}`);
  }
  // the same body without using its symmetry (full mesh) gives the same coefficients
  const coarse = hemisphereQuarter(R, 4, 4); const full = [];
  for (const [sx, sy] of [[1, 1], [-1, 1], [1, -1], [-1, -1]]) {
    for (let p = 0; p < coarse.length; p += 12) {
      const v = [0, 1, 2, 3].map((c) => [sx * coarse[p + 3 * c], sy * coarse[p + 3 * c + 1], coarse[p + 3 * c + 2]]);
      if (sx * sy < 0) v.reverse(); // a single reflection reverses the orientation
      for (const c of v) full.push(...c);
    }
  }
  const sq = radiationSolve({ setup: radiationSetup({ panels: hullPanels(coarse), symmetry: { x: true, y: true } }), omega: 1.2, rho, g });
  const sf = radiationSolve({ setup: radiationSetup({ panels: hullPanels(full) }), omega: 1.2, rho, g });
  for (const i of [0, 1, 2]) { assert.ok(rel(sf.A[i][i], sq.A[i][i]) < 1e-6, `A${i}${i} full vs quarter`); assert.ok(rel(sf.B[i][i], sq.B[i][i]) < 1e-6); }
  assert.ok(Math.abs(sf.A[0][4] - sq.A[0][4]) < 1e-6 * sq.A[0][0] * R);
});

test('Box 20 × 8 × 2 m against Capytaine: added mass, damping and exciting forces from 4 to 25 s', () => {
  const db = fixture('barge');
  const m = boxHydrodynamics({ L: 20, B: 8, T: 2, rho, g });
  let worstA = 0; let worstB = 0; let worstX = 0;
  // damping is compared on the scale of its largest value over the periods (it vanishes in long waves)
  const Bmax = [0, 1, 2, 3, 4].map((i) => Math.max(...db.B.map((b) => b[i][i])));
  for (const T of [4, 5, 6, 8, 10, 14, 25]) {
    const w = Number(((2 * Math.PI) / T).toPrecision(10));
    const { A, B } = boxCoefficients(m, w); const ref = hydroAt(db, w, 180); const refB = hydroAt(db, w, 90);
    for (const i of [0, 1, 2, 4]) { worstA = Math.max(worstA, rel(A[i][i], ref.A[i][i])); worstB = Math.max(worstB, Math.abs(B[i][i] - ref.B[i][i]) / Bmax[i]); }
    // roll about the waterline centre is a small difference of side and bottom pressures: wider band
    assert.ok(rel(A[3][3], ref.A[3][3]) < 0.04, `A44 at ${T} s`);
    const xh = boxExcitation(m, w, 180); const xb = boxExcitation(m, w, 90);
    for (const [x, r, list] of [[xh, ref, [0, 2, 4]], [xb, refB, [1, 2]]]) for (const i of list) {
      worstX = Math.max(worstX, rel(Math.hypot(x.re[i], x.im[i]), Math.hypot(r.F.re[i], r.F.im[i])));
      const dph = Math.atan2(x.im[i], x.re[i]) - Math.atan2(r.F.im[i], r.F.re[i]);
      assert.ok(Math.abs(Math.atan2(Math.sin(dph), Math.cos(dph))) < 0.03, `phase of X${i + 1} at ${T} s`);
    }
    assert.ok(rel(Math.hypot(xb.re[3], xb.im[3]), Math.hypot(refB.F.re[3], refB.F.im[3])) < 0.06, `roll moment at ${T} s`);
  }
  assert.ok(worstA < 0.03, `added mass ${worstA}`);
  assert.ok(worstB < 0.05, `damping ${worstB}`);
  assert.ok(worstX < 0.03, `excitation ${worstX}`);
});

test('Box: symmetric coefficient matrices, positive damping, Haskind energy balance over all headings, limits', () => {
  const L = 20; const B = 8; const T = 2;
  const setup = radiationSetup({ panels: hullPanels(boxQuarterPanels({ L, B, T, nx: 8, ny: 4, nz: 2, grading: 1 })), symmetry: { x: true, y: true } });
  for (const w of [0.6, 1.2, 1.9]) {
    const r = radiationSolve({ setup, omega: w, rho, g });
    // raw asymmetry of the collocation solution (removed in the returned matrices)
    assert.ok(r.symmetryError < 0.04, `asymmetry ${r.symmetryError} at ${w}`);
    for (let i = 0; i < 6; i++) { assert.ok(r.B[i][i] > 0 && r.A[i][i] > 0, `diagonal ${i} at ${w}`); for (let j = 0; j < 6; j++) assert.equal(r.A[i][j], r.A[j][i]); }
    // modes of different symmetry do not couple
    for (const [i, j] of [[0, 1], [0, 2], [1, 2], [2, 3], [2, 4], [3, 4], [0, 3], [1, 4]]) assert.equal(r.A[i][j], 0);
    // B_ii = omega³/(4 pi rho g³) ∫ |X_i(beta)|² d beta (deep water; Newman 1962)
    const n = 48; const sum = [0, 0, 0, 0, 0, 0];
    for (let k = 0; k < n; k++) { const x = haskindExcitation({ setup, phi: r.phi, omega: w, headingDeg: (360 * k) / n, rho, g }); for (let i = 0; i < 6; i++) sum[i] += x.re[i] ** 2 + x.im[i] ** 2; }
    // roll is left out: about the waterline its pressure integral is a small difference of side and bottom parts
    for (const i of [0, 1, 2, 4, 5]) assert.ok(rel(r.B[i][i], ((w ** 3 / (4 * Math.PI * rho * g ** 3)) * sum[i] * 2 * Math.PI) / n) < 0.05, `energy balance mode ${i} at ${w}`);
    // haskindDamping is that integral over a quadrant: same values, symmetric, positive semi-definite sway–roll block
    const Bh = haskindDamping({ setup, phi: r.phi, omega: w, rho, g });
    for (const i of [0, 1, 2, 3, 4, 5]) assert.ok(rel(Bh[i][i], ((w ** 3 / (4 * Math.PI * rho * g ** 3)) * sum[i] * 2 * Math.PI) / n) < 2e-3, `quadrant rule mode ${i}`);
    assert.ok(Bh[1][1] * Bh[3][3] >= Bh[1][3] ** 2 && Bh[0][0] * Bh[4][4] >= Bh[0][4] ** 2 && Bh[1][3] === Bh[3][1] && Bh[0][1] === 0);
    assert.ok(Math.abs(Bh[0][4] - r.B[0][4]) < 0.05 * Math.sqrt(Bh[0][0] * Bh[4][4]), 'surge–pitch coupling');
  }
  // long waves: heave force -> rho g L B, pitch moment -> wave slope × rho g ∫x² dA less the Smith correction
  const lw = radiationSolve({ setup, omega: 0.05, rho, g });
  const x = haskindExcitation({ setup, phi: lw.phi, omega: 0.05, headingDeg: 180, rho, g });
  assert.ok(rel(Math.hypot(x.re[2], x.im[2]), rho * g * L * B) < 2e-3);
  // surge and sway added mass tend to the rigid-lid (omega -> 0) values of the Rankine solver in bem.js
  const a0 = addedMassBem({ surface: box({ size: [L, B, 2 * T], center: [0, 0, 0], divisions: [16, 8, 8] }), rho, limit: 'zero' }).A;
  assert.ok(rel(lw.A[0][0], a0[0][0]) < 0.04, `surge ${lw.A[0][0]} vs ${a0[0][0]}`);
  assert.ok(rel(lw.A[1][1], a0[1][1]) < 0.04, `sway ${lw.A[1][1]} vs ${a0[1][1]}`);
});

test('Box: convergence with the number of panels', () => {
  const L = 20; const B = 8; const T = 2; const w = 1.2;
  const solve = (nx, ny, nz) => radiationSolve({ setup: radiationSetup({ panels: hullPanels(boxQuarterPanels({ L, B, T, nx, ny, nz, grading: 1 })), symmetry: { x: true, y: true } }), omega: w, rho, g, modes: [2, 4] });
  const c = solve(6, 2, 2); const m = solve(9, 3, 3); const f = solve(12, 4, 4);
  // heave and pitch added mass and damping: the change from one mesh to the next shrinks, and is small at the end
  for (const [Q, i] of [['A', 2], ['B', 2], ['A', 4], ['B', 4]]) {
    const d1 = Math.abs(m[Q][i][i] - c[Q][i][i]); const d2 = Math.abs(f[Q][i][i] - m[Q][i][i]);
    assert.ok(d2 < d1, `${Q}${i + 1}${i + 1}: ${d1} then ${d2}`);
    assert.ok(d2 / Math.abs(f[Q][i][i]) < 0.02, `${Q}${i + 1}${i + 1} last change ${d2 / Math.abs(f[Q][i][i])}`);
  }
});

test('Irregular frequency of a box: the null-field equations remove the spurious peak in the heave damping', () => {
  const L = 20; const B = 8; const T = 2;
  const kappa = Math.PI * Math.hypot(1 / L, 1 / B); const wIrr = Math.sqrt((g * kappa) / Math.tanh(kappa * T));
  const m = boxHydrodynamics({ L, B, T, rho, g });
  assert.ok(rel(m.irregularOmega, wIrr) < 1e-12 && wIrr > 2.4 && wIrr < 2.5);
  const ip = [[0.29, 0.37, 0.41], [0.67, 0.71, 0.23], [0.83, 0.19, 0.61], [0.13, 0.87, 0.77]].map(([a, b, c]) => [0.5 * L * a, 0.5 * B * b, -T * c]);
  const setup = radiationSetup({ panels: hullPanels(boxQuarterPanels({ L, B, T, nx: 8, ny: 4, nz: 2, grading: 1 })), symmetry: { x: true, y: true }, interiorPoints: ip });
  const ws = []; for (let w = 2.3; w <= 2.62; w += 0.02) ws.push(w);
  const jump = (nullField) => { const b = ws.map((w) => radiationSolve({ setup, omega: w, rho, g, modes: [2], nullField }).B[2][2]); let j = 0; for (let i = 1; i < b.length; i++) j = Math.max(j, Math.abs(b[i] - b[i - 1])); return { j, b }; };
  const plain = jump(false); const fixed = jump(true);
  // without them B33 swings through the irregular frequency; with them it falls smoothly and stays positive
  assert.ok(plain.j > 5 * fixed.j, `largest step ${plain.j} (plain) vs ${fixed.j}`);
  assert.ok(fixed.b.every((v, i) => v > 0 && (i === 0 || v < fixed.b[i - 1])), 'monotonic, positive');
});

test('Box model: interpolation between solved frequencies, cache, and processor time', () => {
  const L = 60; const B = 16; const T = 4;
  const c0 = process.cpuUsage();
  const m = boxHydrodynamics({ L, B, T, rho, g });
  const cu = process.cpuUsage(c0);
  const ms = (cu.user + cu.system) / 1000;   // processor time: independent of how busy the machine is
  // a fraction of a second on a laptop; generous bound so that a loaded machine does not fail spuriously
  assert.ok(ms < 15000, `box model took ${ms.toFixed(0)} ms`);
  assert.equal(boxHydrodynamics({ L, B, T, rho, g }), m, 'cached per geometry');
  assert.ok(rel((2 * Math.PI * g) / m.omegaMax ** 2, B) < 1e-12, 'solved up to the wavelength equal to the beam');
  // Froude similarity: a hull of the same proportions at another scale, density and gravity is the same solution
  const big = boxHydrodynamics({ L: 2 * L, B: 2 * B, T: 2 * T, rho: 1000, g: 3 });
  assert.equal(big.unit, m.unit, 'one non-dimensional solution per L/B and T/B');
  { const w = 0.6; const w2 = w * Math.sqrt(3 / g) / Math.SQRT2; const a = boxCoefficients(m, w); const b = boxCoefficients(big, w2);
    const xa = boxExcitation(m, w, 130); const xb = boxExcitation(big, w2, 130);
    assert.ok(rel(b.A[2][2], a.A[2][2] * (1000 / rho) * 8) < 1e-12 && rel(b.A[4][4], a.A[4][4] * (1000 / rho) * 32) < 1e-12 && rel(b.A[0][4], a.A[0][4] * (1000 / rho) * 16) < 1e-12);
    assert.ok(rel(b.B[2][2], a.B[2][2] * (1000 / rho) * 8 * Math.sqrt(3 / g) / Math.SQRT2) < 1e-12);
    assert.ok(rel(xb.re[2], xa.re[2] * (1000 * 3) / (rho * g) * 4) < 1e-12 && rel(xb.im[4], xa.im[4] * (1000 * 3) / (rho * g) * 8) < 1e-12); }
  // halfway between two solved frequencies against a direct solution there (in the units of the solution)
  const u = m.unit; const sa = (k) => rho * B ** (k > 2 ? 4 : 3); const sw = Math.sqrt(g / B);
  const scaleA = [0, 1, 2, 3, 4].map((k) => Math.max(...u.A.map((a) => a[k][k])));
  for (const i of [3, 8, 12]) {
    const wu = 0.5 * (u.omega[i] + u.omega[i + 1]); const w = wu * sw;
    const direct = radiationSolve({ setup: u.setup, omega: wu, rho: 1, g: 1, modes: [0, 1, 2, 3, 4], nullField: wu >= 0.9 * u.irregularOmega });
    const { A, B: D } = boxCoefficients(m, w);
    direct.B = haskindDamping({ setup: u.setup, phi: direct.phi, omega: wu, rho: 1, g: 1 });
    for (const k of [0, 1, 2, 3, 4]) { const s2 = sa(k) * B ** (k > 2 ? 1 : 0); assert.ok(Math.abs(A[k][k] / s2 - direct.A[k][k]) < 0.02 * scaleA[k], `A${k} at ${w}`); assert.ok(Math.abs(D[k][k] / (s2 * sw) - direct.B[k][k]) < 0.03 * Math.max(...u.B.map((b) => b[k][k])), `B${k} at ${w}`); }
    const xi = boxExcitation(m, w, 150); const xd = haskindExcitation({ setup: u.setup, phi: direct.phi, omega: wu, headingDeg: 150, rho: 1, g: 1 });
    const fs = (k) => rho * g * B ** (k > 2 ? 3 : 2);
    const scale = [0, 1, 2, 3, 4].map((k) => Math.max(...[60, 150].map((h) => { const x = boxExcitation(m, m.omega[2], h); return Math.hypot(x.re[k], x.im[k]); })));
    for (const k of [0, 1, 2, 3, 4]) assert.ok(Math.hypot(xi.re[k] - fs(k) * xd.re[k], xi.im[k] - fs(k) * xd.im[k]) < 0.03 * scale[k], `X${k} at ${w}`);
  }
  // outside the solved range the end values are held
  for (const k of [0, 2, 4]) {
    assert.ok(rel(boxCoefficients(m, 0.5 * m.omega[0]).A[k][k], boxCoefficients(m, m.omega[0]).A[k][k]) < 1e-9);
    assert.ok(rel(boxCoefficients(m, 2 * m.omegaMax).B[k][k], boxCoefficients(m, m.omegaMax).B[k][k]) < 1e-9);
  }
  assert.throws(() => boxHydrodynamics({ L: 0, B, T }), /must be > 0/);
});

test('Box model: processor time of a new hull once the code is warm (regression guard)', () => {
  // Processor time, as in tests/suites.test.js, and measured against a loop of 10⁶ multiply–adds timed alongside,
  // so that the bound means "the solver got heavier", not "the machine was busy". Recorded on a loaded machine
  // (reference loop 6–7 ms): 75–90 loop units per new hull of 127 panels in the steady state (the version before
  // the split solution: 160–190); an unloaded desktop runs the loop in about 1 ms. The first hulls of a process
  // are not the measure: they run partly in code the engine has not yet optimised.
  const cpu = () => { const u = process.cpuUsage(); return (u.user + u.system) / 1000; };
  const loop = () => { const a = new Float64Array(1e6).fill(1.1); let acc = 0; const t0 = cpu(); for (let k = 0; k < 5; k++) for (let i = 0; i < a.length; i++) acc += a[i] * 1.0001; return (cpu() - t0) / 5 + 0 * acc; };
  for (const L of [87, 89, 91, 93]) boxHydrodynamics({ L, B: 20, T: 5, rho, g }); // warm-up
  let best = Infinity; let ref = Infinity;
  for (const L of [95, 97, 99, 101, 103]) {
    ref = Math.min(ref, loop());
    const t0 = cpu(); const m = boxHydrodynamics({ L, B: 20, T: 5, rho, g }); const t = cpu() - t0;
    ref = Math.min(ref, loop());
    assert.ok(m.panels >= 100 && m.panels <= 130);
    best = Math.min(best, t);
  }
  assert.ok(best < 260 * ref, `new hull: ${best.toFixed(0)} ms of processor time = ${(best / ref).toFixed(0)} reference loops (${ref.toFixed(1)} ms each)`);
  assert.ok(best < 8000, `new hull: ${best.toFixed(0)} ms of processor time`);
  // a hull that differs only in scale is not solved again
  const t0 = cpu(); boxHydrodynamics({ L: 202, B: 40, T: 10, rho, g }); const ts = cpu() - t0;
  assert.ok(ts < 0.05 * best + 5, `scaled hull took ${ts.toFixed(1)} ms`);
});

test('Split solution (real system plus separable wave part) equals the complex elimination; closed-form heading integral', () => {
  const setup = radiationSetup({ panels: hullPanels(boxQuarterPanels({ L: 20, B: 8, T: 2, nx: 8, ny: 4, nz: 2 })), symmetry: { x: true, y: true } });
  for (const w of [0.4, 1.3, 2.2]) {
    const a = radiationSolve({ setup, omega: w, rho, g, split: false }); const b = radiationSolve({ setup, omega: w, rho, g });
    // the two differ only by the interpolation of J0 and J1 in the tabulated path (the split path sums the series)
    for (const k of [0, 1, 2, 3, 4, 5]) {
      let d = 0; let sc = 0;
      for (let i = 0; i < a.phi[k].length; i++) { d = Math.max(d, Math.abs(a.phi[k][i] - b.phi[k][i])); sc = Math.max(sc, Math.abs(a.phi[k][i])); }
      assert.ok(d < 2e-5 * sc, `mode ${k} at ${w}: ${d / sc}`);
      assert.ok(rel(b.A[k][k], a.A[k][k]) < 1e-5, `A${k}${k} at ${w}`);
    }
    // damping: Fourier–Bessel sum against the midpoint rule over 60 headings of a quadrant
    const B1 = haskindDamping({ setup, phi: b.phi, omega: w, rho, g }); const B2 = haskindDamping({ setup, phi: b.phi, omega: w, rho, g, nHeadings: 60 });
    for (const [i, j] of [[0, 0], [1, 1], [2, 2], [3, 3], [4, 4], [5, 5], [1, 3], [0, 4]]) assert.ok(Math.abs(B1[i][j] - B2[i][j]) < 1e-9 * Math.sqrt(B2[i][i] * B2[j][j]), `B${i}${j} at ${w}`);
  }
});
