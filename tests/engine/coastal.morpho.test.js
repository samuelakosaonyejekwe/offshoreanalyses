import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  exnerEvolve, deanParameter, deanProfile, bruunRule, depthOfClosure, depthOfClosureOuter,
} from '../../src/engine/coastal/morpho.js';

const g = 9.80665;
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

// With qs = (1 − p) c zb the Exner equation is linear advection at speed c: after one trip round a periodic
// domain the exact solution is the initial bed.
function advect(n, scheme, limiter, cfl = 0.5) {
  const L = 100;
  const dx = L / n;
  const p = 0.4;
  const z0 = Float64Array.from({ length: n }, (_, i) => Math.exp(-((((i + 0.5) * dx - 50) / 8) ** 2)));
  const q = new Float64Array(n);
  const r = exnerEvolve({
    dx, zb0: z0, porosity: p, tEnd: 100, scheme, limiter, cfl, left: 'periodic', right: 'periodic',
    transport: (z) => {
      for (let i = 0; i < n; i++) q[i] = (1 - p) * 1.0 * z[i];
      return q;
    },
  });
  let e = 0;
  for (let i = 0; i < n; i++) e += Math.abs(r.zb[i] - z0[i]) * dx;
  return { e, r };
}

test('Exner with linear transport: Lax–Wendroff converges at second order, upwind at first order', () => {
  const lw = [100, 200, 400].map((n) => advect(n, 'laxWendroff', 'none').e);
  assert.ok(Math.log2(lw[0] / lw[1]) > 1.9 && Math.log2(lw[1] / lw[2]) > 1.95, `LW errors ${lw}`);
  const vl = [100, 200, 400].map((n) => advect(n, 'laxWendroff', 'vanleer').e);
  assert.ok(Math.log2(vl[0] / vl[1]) > 1.8 && Math.log2(vl[1] / vl[2]) > 1.8, `van Leer errors ${vl}`);
  const up = [200, 400, 800].map((n) => advect(n, 'upwind', 'none').e);
  const o = Math.log2(up[1] / up[2]);
  assert.ok(o > 0.7 && o < 1.1, `upwind order ${o}`);
  assert.ok(vl[2] < 0.02 * up[1]);
});

test('Exner: exact translation at Courant number 1, monotone limited solution, mass conserved to round-off', () => {
  const one = advect(50, 'upwind', 'none', 1.0);
  assert.ok(one.e < 1e-10, `error at unit Courant number ${one.e}`);
  rel(one.r.courantMax, 1, 1e-9);
  for (const lim of ['minmod', 'vanleer', 'superbee']) {
    const { r } = advect(100, 'laxWendroff', lim);
    assert.ok(Math.abs(r.mass.relativeError) < 1e-12);
    assert.equal(r.mass.boundaryNet, 0);
    // TVD: no new extrema
    assert.ok(Math.max(...r.zb) <= 1 + 1e-12 && Math.min(...r.zb) >= -1e-12, lim);
  }
});

test('Exner with nonlinear transport qs = a u³, u = Q/(η − zb): hump migrates downstream, mass balance with open boundaries', () => {
  const n = 200;
  const dx = 5;
  const Q = 10; // unit discharge [m²/s]
  const eta = 10;
  const aCoef = 1e-4;
  const p = 0.4;
  const z0 = Float64Array.from({ length: n }, (_, i) => 2 * Math.exp(-((((i + 0.5) * dx - 300) / 60) ** 2)));
  const q = new Float64Array(n);
  const transport = (z) => {
    for (let i = 0; i < n; i++) q[i] = aCoef * (Q / (eta - z[i])) ** 3;
    return q;
  };
  const centroid = (z) => {
    let m = 0;
    let mx = 0;
    for (let i = 0; i < n; i++) {
      m += z[i];
      mx += z[i] * (i + 0.5) * dx;
    }
    return mx / m;
  };
  const tEnd = 3.0e5;
  const r = exnerEvolve({ dx, zb0: z0, transport, tEnd, porosity: p });
  assert.ok(r.courantMax <= 0.8 + 1e-9);
  assert.ok(Math.abs(r.mass.error) < 1e-9 * r.mass.initial, `mass error ${r.mass.error}`);
  // Zero-gradient boundaries on a flat far field: inflow equals outflow, so the hump volume is unchanged.
  rel(r.mass.final, r.mass.initial, 1e-6);
  // Bed celerity c = (1/(1−p)) dqs/dz = 3 a Q³ /((1−p)(η − z)⁴): between its flat-bed and crest values.
  const cLow = (3 * aCoef * Q ** 3) / ((1 - p) * eta ** 4);
  const cHigh = (3 * aCoef * Q ** 3) / ((1 - p) * (eta - 2) ** 4);
  const shift = centroid(r.zb) - centroid(z0);
  assert.ok(shift > cLow * tEnd && shift < cHigh * tEnd, `shift ${shift} vs [${cLow * tEnd}, ${cHigh * tEnd}]`);
  // The crest travels faster than the foot, so the lee face steepens (a migrating dune).
  let fMax = 0;
  let bMax = 0;
  for (let i = 1; i < n; i++) {
    const s = (r.zb[i] - r.zb[i - 1]) / dx;
    fMax = Math.max(fMax, -s);
    bMax = Math.max(bMax, s);
  }
  assert.ok(fMax > 1.5 * bMax);
  assert.ok(Math.max(...r.zb) <= 2 + 1e-9);
});

test('Exner with a prescribed upstream sediment supply: bed volume grows by the supply surplus', () => {
  const n = 50;
  const dx = 2;
  const q = new Float64Array(n).fill(1e-4);
  const r = exnerEvolve({
    dx, zb0: new Float64Array(n), transport: () => q, tEnd: 1000, dt: 10, porosity: 0.35, left: { q: 3e-4 }, right: 'zeroGradient',
  });
  rel(r.mass.boundaryNet, (3e-4 - 1e-4) * 1000, 1e-12);
  rel(r.mass.final - r.mass.initial, 0.2, 1e-10);
  // Uniform transport capacity: all the surplus deposits in the first cell.
  rel(r.zb[0], 0.2 / (0.65 * dx), 1e-10);
  assert.equal(r.zb[10], 0);
  assert.equal(r.steps, 100);
  assert.throws(() => exnerEvolve({ dx, zb0: new Float64Array(n), transport: () => q, tEnd: 1, left: 'periodic' }), RangeError);
});

test('Dean equilibrium profile: h = A y^(2/3), A = 0.067 ws^0.44 (ws in cm/s)', () => {
  // ws = 2.6 cm/s (about 0.2 mm sand): A = 0.067·2.6^0.44 = 0.102 m^(1/3)
  const A = deanParameter(0.026);
  rel(A, 0.067 * 2.6 ** 0.44, 1e-12);
  assert.ok(Math.abs(A - 0.102) < 1e-3);
  const p = deanProfile({ y: [0, 100, 500], A: 0.1 });
  assert.equal(p.h[0], 0);
  rel(p.h[1], 0.1 * 100 ** (2 / 3), 1e-12); // 2.154 m
  rel(p.h[2] / p.h[1], 5 ** (2 / 3), 1e-12);
  rel(p.widthToDepth(p.h[2]), 500, 1e-12);
  rel(deanProfile({ y: 100, A: 0.1 }).h, 2.1544, 1e-4);
});

test('Bruun rule: R = S L/(B + h*) and its sediment balance on a translated Dean profile', () => {
  // S = 0.5 m, L = 1000 m, B = 2 m, h* = 8 m: R = 0.5·1000/10 = 50 m
  const r = bruunRule({ S: 0.5, hClosure: 8, B: 2, L: 1000 });
  rel(r.retreat, 50, 1e-12);
  rel(r.ratio, 100, 1e-12);
  const d = bruunRule({ S: 0.5, hClosure: 8, B: 2, A: 0.1 });
  rel(d.L, 80 ** 1.5, 1e-12);
  rel(d.retreat, (0.5 * 80 ** 1.5) / 10, 1e-12);
  // Balance: eroded volume R (B + h*) equals the volume needed to raise the active profile, S L.
  rel(d.retreat * (2 + 8), 0.5 * d.L, 1e-12);
  // Numerical check of that balance: translate the profile z(y) (berm of height B, Dean below) up by S and
  // landward by R and integrate the difference out to the closure depth.
  const A = 0.1;
  const B = 2;
  const S = 0.01;
  const hc = 8;
  const L = (hc / A) ** 1.5;
  const R = bruunRule({ S, hClosure: hc, B, L }).retreat;
  const z = (y) => (y < 0 ? B : -A * y ** (2 / 3));
  const n = 2000000;
  const y0 = -5;
  const y1 = L;
  const dy = (y1 - y0) / n;
  let dV = 0;
  for (let i = 0; i < n; i++) {
    const y = y0 + (i + 0.5) * dy;
    dV += (z(y + R) + S - z(y)) * dy;
  }
  // Net change vanishes to first order in S (residual is the O(S·R) edge term at the closure depth).
  assert.ok(Math.abs(dV) < 0.02 * S * L, `net volume change ${dV} vs S L = ${S * L}`);
});

test('depth of closure: Hallermeier and Birkemeier by hand', () => {
  // He = 3 m, Te = 9 s: 2.28·3 − 68.5·9/(g·81) = 6.84 − 0.776 = 6.06 m
  const h = depthOfClosure({ He: 3, Te: 9 });
  rel(h, 2.28 * 3 - (68.5 * 9) / (g * 81), 1e-12);
  assert.ok(Math.abs(h - 6.06) < 0.01);
  rel(depthOfClosure({ He: 3, Te: 9, method: 'birkemeier' }), 1.75 * 3 - (57.9 * 9) / (g * 81), 1e-12);
  assert.ok(depthOfClosure({ He: 3, Te: 9, method: 'birkemeier' }) < h);
  // Outer limit: (1.0 − 0.3·0.6)·7·√(g/(5000·0.25e-3)) = 0.82·7·2.801 = 16.1 m
  const o = depthOfClosureOuter({ HsMean: 1.0, sigmaHs: 0.6, TsMean: 7, d50: 0.25e-3 });
  rel(o, 0.82 * 7 * Math.sqrt(g / 1.25), 1e-12);
  assert.ok(o > h);
  assert.throws(() => depthOfClosure({ He: 3, Te: 9, method: 'x' }), RangeError);
});

test('Exner with spatially varying transport on a level bed: bounded step, bounded scour, grid-convergent (local celerity)', () => {
  // qs = a u³ with u = Q(x)/(η − zb): the unit discharge rises by 50 % through a contraction centred at x = 500 m.
  // The bed starts level, so the jump ratio Δqs/Δzb between neighbours is infinite; the local celerity
  // 3 qs/((1 − p)(η − zb)) is finite everywhere. Equilibrium scour in the throat makes u uniform:
  // η − zb = 5 × 1.5 → zb = −2.5 m, which bounds the erosion.
  const p = 0.4; const aC = 1e-4; const eta = 5; const L = 1000;
  const run = (n, extra = {}) => {
    const dx = L / n;
    const Qx = Float64Array.from({ length: n }, (_, i) => 5 * (1 + 0.5 * Math.exp(-((((i + 0.5) * dx - 500) / 100) ** 2))));
    return exnerEvolve({ dx, zb0: new Float64Array(n), porosity: p, tEnd: 2e6, maxSteps: 20000,
      transport: (z) => Float64Array.from(z, (zi, i) => aC * (Qx[i] / (eta - zi)) ** 3), ...extra });
  };
  const a = run(50); const b = run(100);
  for (const r of [a, b]) {
    assert.ok(r.courantMax <= 0.8 + 1e-9, `Courant ${r.courantMax}`);
    assert.ok(Math.min(...r.zb) > -2.5 && Math.min(...r.zb) < -2.0, `scour ${Math.min(...r.zb)}`);
    assert.ok(Math.abs(r.mass.error) < 1e-9);
    assert.ok(r.dtMin > 10 && r.dtMax < 1e5, `steps ${r.dtMin}…${r.dtMax}`);
  }
  rel(Math.min(...b.zb), Math.min(...a.zb), 0.01);
  // the raw jump ratio collapses the step on the same problem
  assert.throws(() => run(50, { celerity: 'jump' }), /maxSteps/);
  // a supplied celerity function reproduces the local estimate
  const dx = L / 50;
  const Qx = Float64Array.from({ length: 50 }, (_, i) => 5 * (1 + 0.5 * Math.exp(-((((i + 0.5) * dx - 500) / 100) ** 2))));
  const c = run(50, { celerity: (z) => Float64Array.from(z, (zi, i) => (3 * aC * Qx[i] ** 3) / ((1 - p) * (eta - zi) ** 4)) });
  rel(Math.min(...c.zb), Math.min(...a.zb), 1e-4);
});

test('Exner bed-slope term with uniform transport is linear diffusion: a Gaussian mound spreads with σ² = σ0² + 2 D t', () => {
  // qs = q0 everywhere → ∂zb/∂t = D ∂²zb/∂x², D = β q0/(1 − p) = 2 × 1e-4 / 0.6 = 3.333e-4 m²/s.
  // σ0 = 10 m, t = 3e5 s: σ² = 100 + 200 = 300 → peak ratio σ0/σ = 1/√3 = 0.5774.
  const n = 200; const dx = 1; const q = new Float64Array(n).fill(1e-4);
  const z0 = Float64Array.from({ length: n }, (_, i) => Math.exp(-0.5 * (((i + 0.5) * dx - 100) / 10) ** 2));
  const r = exnerEvolve({ dx, zb0: z0, transport: () => q, tEnd: 3e5, porosity: 0.4, bedSlopeBeta: 2 });
  rel(Math.max(...r.zb), 1 / Math.sqrt(3), 2e-3);
  let m = 0; let v = 0;
  for (let i = 0; i < n; i++) { m += r.zb[i]; v += r.zb[i] * ((i + 0.5) * dx - 100) ** 2; }
  rel(v / m, 300, 2e-3);
  assert.ok(Math.abs(r.mass.relativeError) < 1e-12);
  // explicit diffusion limit Δt <= cfl Δx² (1 − p)/(2 β q0) = 0.8 × 0.6 / 4e-4 = 1200 s
  assert.ok(r.dtMax <= 1200 * (1 + 1e-9), `dtMax ${r.dtMax}`);
  assert.throws(() => exnerEvolve({ dx, zb0: z0, transport: () => q, tEnd: 1, bedSlopeBeta: -1 }), /bedSlopeBeta/);
});

test('Exner adaptive step with zero celerity: supply onto a level bed is stepped by the bed-change limit, dtMax is honoured', () => {
  // uniform capacity 1e-4, supply 3e-4 m²/s for 1000 s, p = 0.35, dx = 2 m: first cell rises 0.2/(0.65 × 2) = 0.15385 m
  const n = 50; const dx = 2; const q = new Float64Array(n).fill(1e-4);
  const base = { dx, zb0: new Float64Array(n), transport: () => q, tEnd: 1000, porosity: 0.35, left: { q: 3e-4 } };
  const r = exnerEvolve(base);
  rel(r.zb[0], 0.2 / (0.65 * dx), 1e-10);
  // default limit cfl × 0.001 Δx = 1.6 mm per step → 0.15385 / 0.0016 = 96.2 → 97 steps
  assert.equal(r.steps, 97);
  const one = exnerEvolve({ ...base, maxBedChange: 1 });
  assert.equal(one.steps, 1); rel(one.zb[0], 0.2 / (0.65 * dx), 1e-10);
  const capped = exnerEvolve({ ...base, maxBedChange: 1, dtMax: 50 });
  assert.equal(capped.steps, 20); assert.ok(capped.dtMax <= 50 + 1e-9);
  assert.throws(() => exnerEvolve({ ...base, maxBedChange: 0 }), /maxBedChange/);
  assert.throws(() => exnerEvolve({ ...base, celerity: 'roe' }), /celerity/);
});
