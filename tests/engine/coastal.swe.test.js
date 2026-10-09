import { test } from 'node:test';
import assert from 'node:assert/strict';
import { solveSwe1d, damBreakAnalytic } from '../../src/engine/coastal/swe.js';

const g = 9.80665;
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const open = { left: { type: 'transmissive' }, right: { type: 'transmissive' } };

test('Stoker analytic solution satisfies the Rankine–Hugoniot conditions across the bore', () => {
  const s = damBreakAnalytic({ hL: 1, hR: 0.1, x: 0, t: 1 });
  // mass: S (hm − hR) = hm um ; momentum: S hm um = hm um² + g/2 (hm² − hR²)
  rel(s.shockSpeed * (s.hm - 0.1), s.hm * s.um, 1e-12);
  rel(s.shockSpeed * s.hm * s.um, s.hm * s.um ** 2 + 0.5 * g * (s.hm ** 2 - 0.01), 1e-10);
  rel(s.um, 2 * (Math.sqrt(g) - Math.sqrt(g * s.hm)), 1e-12); // Riemann invariant through the rarefaction
  // At the dam site the rarefaction is critical: h = 4/9 hL, u = 2/3 √(g hL) (also Ritter's dry-bed values)
  rel(s.h, 4 / 9, 1e-12);
  rel(s.u, (2 / 3) * Math.sqrt(g), 1e-12);
  const r = damBreakAnalytic({ hL: 1, x: 0, t: 1 });
  rel(r.h, 4 / 9, 1e-12);
  rel(r.shockSpeed, 2 * Math.sqrt(g), 1e-12);
});

function damBreak(n, hR, tEnd, limiter = 'minmod') {
  const L = 20;
  const dx = L / n;
  const h0 = Float64Array.from({ length: n }, (_, i) => ((i + 0.5) * dx < 10 ? 1 : hR));
  const s = solveSwe1d({ dx, zb: new Float64Array(n), h0, tEnd, limiter, ...open });
  let e = 0;
  for (let i = 0; i < n; i++) e += Math.abs(s.h[i] - damBreakAnalytic({ hL: 1, hR, x: s.x[i] - 10, t: tEnd }).h) * dx;
  return { s, e };
}

test('wet-bed dam break converges to the Stoker (1957) solution; plateau depth, velocity and bore position', () => {
  const t = 1.5;
  const a = damBreak(200, 0.1, t);
  const b = damBreak(400, 0.1, t);
  const c = damBreak(800, 0.1, t);
  assert.ok(c.e < 0.015, `L1 error ${c.e}`);
  // Discontinuous solution: first-order convergence in L1 is the best possible.
  assert.ok(a.e / b.e > 1.6 && b.e / c.e > 1.6, `ratios ${a.e / b.e}, ${b.e / c.e}`);
  const ex = damBreakAnalytic({ hL: 1, hR: 0.1, x: 0, t });
  // Middle of the plateau between the rarefaction tail and the bore
  const xm = 10 + 0.5 * (ex.um - Math.sqrt(g * ex.hm) + ex.shockSpeed) * t;
  const im = Math.round(xm / (20 / 800) - 0.5);
  rel(c.s.h[im], ex.hm, 2e-3, 'plateau depth');
  rel(c.s.u[im], ex.um, 2e-3, 'plateau velocity');
  // Bore front: last cell with h above the mid-jump level
  let ifr = 0;
  for (let i = 0; i < 800; i++) if (c.s.h[i] > 0.5 * (ex.hm + 0.1)) ifr = i;
  assert.ok(Math.abs(c.s.x[ifr] - (10 + ex.shockSpeed * t)) < 3 * (20 / 800), 'bore position within three cells');
  assert.ok(Math.abs(c.s.volume.relativeError) < 1e-12);
  // No new extrema (TVD): depths stay within the initial bounds.
  assert.ok(Math.max(...c.s.h) <= 1 + 1e-9 && Math.min(...c.s.h) >= 0.1 - 1e-9);
});

test('dry-bed dam break converges to the Ritter (1892) solution and stays non-negative', () => {
  const t = 1.0;
  const a = damBreak(400, 0, t);
  const b = damBreak(1600, 0, t);
  assert.ok(b.e < 0.01, `L1 error ${b.e}`);
  assert.ok(b.e < 0.6 * a.e);
  assert.ok(Math.min(...b.s.h) >= 0);
  assert.equal(b.s.volume.clipped, 0);
  assert.ok(Math.abs(b.s.volume.relativeError) < 1e-12);
  // Critical flow at the dam site and the parabolic profile inside the fan
  const dx = 20 / 1600;
  const at = (x) => Math.round(x / dx - 0.5);
  rel(0.5 * (b.s.h[at(10) - 1] + b.s.h[at(10)]), 4 / 9, 5e-3);
  for (const xr of [-2, -1, 1, 2, 3]) {
    rel(b.s.h[at(10 + xr)], (2 * Math.sqrt(g) - xr / t) ** 2 / (9 * g), 1e-2, `x=${xr}`);
    rel(b.s.u[at(10 + xr)], (2 / 3) * (Math.sqrt(g) + xr / t), 1e-2, `u at x=${xr}`);
  }
  // The numerical front (1 mm depth) trails the theoretical tip 2√(g hL) t, where the depth vanishes
  // quadratically, but is well beyond the 1 mm contour of the exact solution minus a few cells.
  const tip = 10 + 2 * Math.sqrt(g) * t;
  const xs = b.s.shoreline.x[b.s.shoreline.x.length - 1];
  assert.ok(xs < tip + dx && xs > tip - 0.8, `front ${xs} vs tip ${tip}`);
});

test('lake at rest over an uneven bed with islands and a dry beach is preserved to round-off', () => {
  const n = 240;
  const dx = 0.1;
  const zb = Float64Array.from({ length: n }, (_, i) => {
    const x = (i + 0.5) * dx;
    return 0.8 * Math.exp(-((x - 7) ** 2) / 2) + 1.5 * Math.exp(-((x - 14) ** 2)) + 0.2 * Math.sin(x) + (x > 20 ? 0.6 * (x - 20) : 0);
  });
  assert.ok(Math.max(...zb) > 1.2, 'the bed pierces the surface');
  for (const limiter of ['minmod', 'vanleer', 'none']) {
    for (const bc of [{ type: 'wall' }, { type: 'transmissive' }]) {
      const s = solveSwe1d({ dx, zb, eta0: 1.0, tEnd: 5, limiter, left: bc, right: { type: 'wall' }, manning: 0.02 });
      let umax = 0;
      let emax = 0;
      for (let i = 0; i < n; i++) {
        umax = Math.max(umax, Math.abs(s.u[i]));
        if (s.h[i] > 0) emax = Math.max(emax, Math.abs(s.eta[i] - 1));
        else assert.ok(zb[i] >= 1);
      }
      assert.ok(umax < 1e-12, `${limiter}: spurious velocity ${umax}`);
      assert.ok(emax < 1e-13, `${limiter}: surface error ${emax}`);
      assert.ok(s.steps > 100);
    }
  }
});

// Right-running simple wave on a flat bed: u = 2 (c − c0) everywhere, so each depth value travels at
// u + c = 3c − 2c0 and the exact solution follows from the characteristics until breaking.
function simpleWaveError(n, limiter) {
  const L = 40;
  const dx = L / n;
  const c0 = Math.sqrt(g);
  const hi = (x) => 1 + 0.1 * Math.exp(-(((x - 12) / 3) ** 2));
  const h0 = new Float64Array(n);
  const u0 = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    h0[i] = hi((i + 0.5) * dx);
    u0[i] = 2 * (Math.sqrt(g * h0[i]) - c0);
  }
  const t = 2;
  const s = solveSwe1d({ dx, zb: new Float64Array(n), h0, u0, tEnd: t, limiter, ...open });
  let e = 0;
  for (let i = 0; i < n; i++) {
    const x = s.x[i];
    let lo = x - 10 * t;
    let up = x;
    for (let k = 0; k < 80; k++) {
      const m = 0.5 * (lo + up);
      if (m + (3 * Math.sqrt(g * hi(m)) - 2 * c0) * t - x > 0) up = m;
      else lo = m;
    }
    e += Math.abs(s.h[i] - hi(0.5 * (lo + up))) * dx;
  }
  return e;
}

test('second-order convergence on a smooth simple wave (exact solution by characteristics)', () => {
  const e = [100, 200, 400, 800].map((n) => simpleWaveError(n, 'vanleer'));
  const orders = e.slice(1).map((v, i) => Math.log2(e[i] / v));
  assert.ok(orders[1] > 1.85 && orders[2] > 1.9, `van Leer orders ${orders}`);
  assert.ok(e[3] < 4e-4);
  const m = [200, 400, 800].map((n) => simpleWaveError(n, 'minmod'));
  // minmod clips the smooth crest, which costs a fraction of an order in L1
  assert.ok(Math.log2(m[0] / m[1]) > 1.6 && Math.log2(m[1] / m[2]) > 1.7, `minmod errors ${m}`);
  const f = [200, 400].map((n) => simpleWaveError(n, 'none'));
  const o1 = Math.log2(f[0] / f[1]);
  assert.ok(o1 > 0.7 && o1 < 1.3, `first-order scheme order ${o1}`);
  assert.ok(e[2] < 0.1 * f[1]);
});

test('Manning friction: uniform flow decays as u0 / (1 + g n² u0 t / h^(4/3))', () => {
  const n = 16;
  const nm = 0.03;
  const h = 2;
  const u0 = 1.5;
  const per = { type: 'periodic' };
  const s = solveSwe1d({
    dx: 1, zb: new Float64Array(n), h0: new Float64Array(n).fill(h), u0: new Float64Array(n).fill(u0), tEnd: 200,
    manning: nm, left: per, right: per, outputTimes: [50, 200],
  });
  const exact = (t) => u0 / (1 + (g * nm * nm * u0 * t) / h ** (4 / 3));
  assert.equal(s.snapshots.length, 2);
  assert.equal(s.snapshots[0].t, 50);
  rel(s.snapshots[0].u[5], exact(50), 5e-4);
  rel(s.u[5], exact(200), 5e-4);
  rel(s.h[5], h, 1e-12);
});

test('run-up of a pulse on a plane beach: volume conserved to round-off, shoreline history and envelope', () => {
  const n = 400;
  const dx = 0.25;
  const zb = Float64Array.from({ length: n }, (_, i) => {
    const x = (i + 0.5) * dx;
    return x < 50 ? 0 : (x - 50) * 0.05;
  });
  const eta0 = Float64Array.from({ length: n }, (_, i) => 1 + 0.3 * Math.exp(-((((i + 0.5) * dx - 20) / 4) ** 2)));
  const s = solveSwe1d({ dx, zb, eta0, tEnd: 40, manning: 0.01 });
  assert.ok(Math.abs(s.volume.relativeError) < 1e-12, `volume error ${s.volume.relativeError}`);
  assert.equal(s.volume.boundaryNet, 0);
  assert.equal(s.volume.clipped, 0);
  assert.ok(Math.min(...s.h) >= 0);
  // Still-water shoreline at x = 70 m (z = 1). The pulse must run up above it and recede.
  rel(s.shoreline.x[0], 70, 0.01);
  assert.ok(s.shoreline.maxRunup > 1.1 && s.shoreline.maxRunup < 1.6, `run-up ${s.shoreline.maxRunup}`);
  assert.ok(s.envelope.inundationLimit > 72 && s.envelope.inundationLimit < 82);
  assert.ok(s.shoreline.maxRunupTime > 5 && s.shoreline.maxRunupTime < 40);
  assert.equal(s.shoreline.t.length, s.steps + 1);
  // Envelope consistency
  for (let i = 0; i < n; i++) assert.ok(s.envelope.hMax[i] >= s.h[i]);
  const iLim = Math.round(s.envelope.inundationLimit / dx) + 2;
  assert.ok(Number.isNaN(s.envelope.etaMax[iLim]) && s.envelope.etaMax[iLim - 4] > 1);
});

test('discharge boundary: inflow volume balance and the bore it drives satisfy the jump conditions', () => {
  const n = 400;
  const dx = 0.25;
  const q = 0.5;
  const t = 15;
  const s = solveSwe1d({
    dx, zb: new Float64Array(n), eta0: 1, tEnd: t, left: { type: 'discharge', q }, right: { type: 'transmissive' },
  });
  // The discharge is imposed through the ghost state, so the numerical boundary flux matches it closely
  // but not identically; the conservation diagnostic uses the actual flux and closes to round-off.
  rel(s.volume.boundaryNet, q * t, 5e-4);
  assert.ok(Math.abs(s.volume.error) < 1e-10);
  // Exact: constant state (hm, um = q/hm) behind a bore moving into still water of depth 1:
  // S = q/(hm − 1) and q²/hm + g hm²/2 − g/2 = S q
  let lo = 1.0001;
  let hi = 2;
  for (let k = 0; k < 100; k++) {
    const hm = 0.5 * (lo + hi);
    const f = (q * q) / hm + 0.5 * g * (hm * hm - 1) - (q * q) / (hm - 1);
    if (f > 0) hi = hm;
    else lo = hm;
  }
  const hm = 0.5 * (lo + hi);
  rel(s.h[40], hm, 2e-3, 'depth behind the bore');
  rel(s.h[40] * s.u[40], q, 2e-3, 'discharge behind the bore');
  let ifr = 0;
  for (let i = 0; i < n; i++) if (s.h[i] > 0.5 * (hm + 1)) ifr = i;
  assert.ok(Math.abs(s.x[ifr] - (q / (hm - 1)) * t) < 4 * dx, 'bore position');
});

test('stage boundary: a slow rise of the boundary level fills a closed basin to that level', () => {
  const n = 40;
  const dx = 1;
  // Smooth 400 s ramp from 1.0 m to 1.5 m, long compared with the 42 s seiche period of the basin.
  const eta = (t) => (t < 400 ? 1.25 - 0.25 * Math.cos((Math.PI * t) / 400) : 1.5);
  const s = solveSwe1d({ dx, zb: new Float64Array(n), eta0: 1, tEnd: 600, manning: 0.05, left: { type: 'stage', eta }, right: { type: 'wall' } });
  for (const i of [0, 20, 39]) rel(s.h[i], 1.5, 2e-3, `cell ${i}`);
  rel(s.volume.boundaryNet, 0.5 * n * dx, 5e-3);
  assert.ok(Math.abs(s.volume.error) < 1e-10);
  // A tabulated series gives the same result as the equivalent function.
  const ts = Array.from({ length: 61 }, (_, i) => i * 10);
  const series = { t: ts, values: ts.map(eta) };
  const s2 = solveSwe1d({ dx, zb: new Float64Array(n), eta0: 1, tEnd: 600, manning: 0.05, left: { type: 'stage', eta: series }, right: { type: 'wall' } });
  rel(s2.h[20], s.h[20], 1e-3);
});

test('incident-wave boundary generates the prescribed long wave and absorbs outgoing waves', () => {
  const n = 400;
  const dx = 0.25;
  const zb = new Float64Array(n);
  // Absorption: an initial hump leaves through both characteristic boundaries without reflection.
  const hump = Float64Array.from({ length: n }, (_, i) => 1 + 0.05 * Math.exp(-((((i + 0.5) * dx - 50) / 4) ** 2)));
  const quiet = { type: 'wave', eta: 0, depth: 1 };
  const a = solveSwe1d({ dx, zb, eta0: hump, tEnd: 60, left: quiet, right: quiet });
  let dev = 0;
  for (let i = 0; i < n; i++) dev = Math.max(dev, Math.abs(a.eta[i] - 1));
  assert.ok(dev < 1e-6, `residual ${dev}`);
  rel(a.volume.final, 100, 1e-9);
  // Generation: a small-amplitude sine train of period 8 s travels at √(g h) with the prescribed amplitude.
  const amp = 0.02;
  const b = solveSwe1d({ dx, zb, eta0: 1, tEnd: 28, left: { type: 'wave', amplitude: amp, period: 8, depth: 1 }, right: quiet });
  let mx = 0;
  let mn = 0;
  for (let i = 40; i < 240; i++) {
    mx = Math.max(mx, b.envelope.etaMax[i] - 1);
    mn = Math.min(mn, b.eta[i] - 1);
  }
  rel(mx, amp, 0.05, 'crest amplitude');
  assert.ok(mn < -0.8 * amp);
  // Wavelength = √(g h) T ≈ 25 m: surface elevations half a wavelength apart are in antiphase.
  const Lw = Math.sqrt(g) * 8;
  const i1 = 120;
  const i2 = i1 + Math.round(Lw / 2 / dx);
  assert.ok(Math.abs(b.eta[i1] - 1 + (b.eta[i2] - 1)) < 0.2 * amp, 'antiphase half a wavelength apart');
  assert.ok(Math.abs(b.volume.error) < 1e-10);
});

test('input validation', () => {
  const zb = new Float64Array(10);
  assert.throws(() => solveSwe1d({ dx: 0, zb, eta0: 1, tEnd: 1 }), RangeError);
  assert.throws(() => solveSwe1d({ dx: 1, zb, tEnd: 1 }), RangeError);
  assert.throws(() => solveSwe1d({ dx: 1, zb, eta0: 1, tEnd: 1, cfl: 0.9 }), RangeError);
  assert.throws(() => solveSwe1d({ dx: 1, zb, eta0: 1, tEnd: 1, left: { type: 'periodic' } }), RangeError);
  assert.throws(() => solveSwe1d({ dx: 1, zb, eta0: 1, tEnd: 1, left: { type: 'sponge' } }), RangeError);
  assert.throws(() => damBreakAnalytic({ hL: 1, hR: 2, x: 0, t: 1 }), RangeError);
});
