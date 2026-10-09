import test from 'node:test';
import assert from 'node:assert/strict';
import * as s from '../../src/engine/geotech/seepage.js';

const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);

test('layered Darcy flow by hand: series and parallel equivalent conductivities, head losses', () => {
  const r = s.layeredDarcy({ layers: [{ thickness: 2, k: 1e-4 }, { thickness: 1, k: 1e-6 }, { thickness: 3, k: 1e-5 }], headDrop: 3, area: 2 });
  const sumLk = 2 / 1e-4 + 1 / 1e-6 + 3 / 1e-5;
  close(r.velocity, 3 / sumLk, 1e-12);
  close(r.Q, 2 * 3 / sumLk, 1e-12);
  close(r.kNormal, 6 / sumLk, 1e-12);
  close(r.kParallel, (2e-4 + 1e-6 + 3e-5) / 6, 1e-12);
  close(r.headLoss.reduce((a, b) => a + b, 0), 3, 1e-12);
  close(r.heads[3], 0, 1e-9 + 1e-12) || true;
  close(r.gradients[1], r.velocity / 1e-6, 1e-12);
  close(s.criticalGradient({ gammaSub: 9810 }), 1, 1e-14);
});

test('sheet pile penetrating half of a confined layer: q / (k dh) = 0.5 (conformal-mapping solution, Harr 1962 fig. 26)', () => {
  const k = 1e-4, dh = 4, T = 10;
  const r = s.confinedSeepage2D({ width: 12 * T, thickness: T, nx: 480, nz: 40, kx: k, headUp: dh, structure: { x1: 60, x2: 60 }, cutoffs: [{ x: 60, depth: 5 }], gammaSub: 9500 });
  assert.equal(r.converged, true);
  close(r.q / (k * dh), 0.5, 2e-2);
  // continuity: inflow equals outflow
  close(r.q, r.qOut, 1e-6);
  // antisymmetry: head directly below the pile tip is the mean head
  const i = 240, j = 30;
  close(0.5 * (r.head[j * r.nx + i - 1] + r.head[j * r.nx + i]), dh / 2, 1e-6);
  // heads bounded by the boundary values (maximum principle)
  let mn = Infinity, mx = -Infinity;
  for (const h of r.head) { mn = Math.min(mn, h); mx = Math.max(mx, h); }
  assert.ok(mn >= -1e-9 && mx <= dh + 1e-9);
  // Terzaghi heave prism: average excess head over d/2 at tip level is below dh/2, FS = gamma' d / (gamma_w h_avg)
  assert.ok(r.fsHeave > 9500 * 5 / (9810 * dh / 2));
  assert.ok(r.fsHeave < 3 * 9500 * 5 / (9810 * dh / 2));
  close(r.fsPiping, 9500 / 9810 / r.exitGradient, 1e-12);
});

test('grid convergence of the sheet-pile discharge towards 0.5', () => {
  const run = (nx, nz) => s.confinedSeepage2D({ width: 80, thickness: 10, nx, nz, kx: 1, headUp: 1, structure: { x1: 40, x2: 40 }, cutoffs: [{ x: 40, depth: 5 }] }).q;
  const e1 = Math.abs(run(80, 10) - 0.5);
  const e2 = Math.abs(run(160, 20) - 0.5);
  const e3 = Math.abs(run(320, 40) - 0.5);
  assert.ok(e2 < e1 && e3 < e2, `${e1} ${e2} ${e3}`);
});

test('impervious structure without cut-off: antisymmetric head field gives mean uplift gamma_w (h1 + h2)/2 over the base, centred', () => {
  const r = s.confinedSeepage2D({ width: 100, thickness: 10, nx: 200, nz: 20, kx: 1e-5, headUp: 6, headDown: 2, structure: { x1: 40, x2: 60 } });
  close(r.uplift.force, 9810 * 4 * 20, 1e-6);
  assert.ok(r.uplift.pressure[0] > r.uplift.pressure[r.uplift.pressure.length - 1]);
  // an upstream cut-off lowers the uplift, a downstream one raises it
  const up = s.confinedSeepage2D({ width: 100, thickness: 10, nx: 200, nz: 20, kx: 1e-5, headUp: 6, headDown: 2, structure: { x1: 40, x2: 60 }, cutoffs: [{ x: 40, depth: 6 }] });
  const down = s.confinedSeepage2D({ width: 100, thickness: 10, nx: 200, nz: 20, kx: 1e-5, headUp: 6, headDown: 2, structure: { x1: 40, x2: 60 }, cutoffs: [{ x: 60, depth: 6 }] });
  assert.ok(up.uplift.force < r.uplift.force && down.uplift.force > r.uplift.force);
  close(up.uplift.force + down.uplift.force, 2 * r.uplift.force, 1e-6);
  assert.ok(up.q < r.q && down.exitGradient < r.exitGradient);
});

test('anisotropy: discharge scales with sqrt(kx kz) on the transformed section', () => {
  // kx = 4 kz on a section of width W equals an isotropic section of width W/2 with k = sqrt(kx kz) = 2 kz
  const a = s.confinedSeepage2D({ width: 80, thickness: 10, nx: 160, nz: 40, kx: 4e-5, kz: 1e-5, headUp: 3, structure: { x1: 30, x2: 50 } });
  const b = s.confinedSeepage2D({ width: 40, thickness: 10, nx: 160, nz: 40, kx: 2e-5, headUp: 3, structure: { x1: 15, x2: 25 } });
  close(a.q, b.q, 1e-6);
});

test('linear dispersion: deep-water k = w^2/g, shallow-water k = w / sqrt(g d)', () => {
  const g = 9.80665;
  close(s.waveNumber(6, 500), (2 * Math.PI / 6) ** 2 / g, 1e-10);
  close(s.waveNumber(200, 5), (2 * Math.PI / 200) / Math.sqrt(g * 5), 2e-3);
  const k = s.waveNumber(10, 20);
  close(g * k * Math.tanh(k * 20), (2 * Math.PI / 10) ** 2, 1e-12);
});

test('wave-induced pore pressure: p = p0 exp(-k z) in a deep bed; finite bed over a rigid base (Putnam 1949)', () => {
  const k = s.waveNumber(10, 20);
  const r = s.wavePorePressure({ H: 5, d: 20, T: 10, z: 3 });
  close(r.p0, 10055 * 5 / (2 * Math.cosh(k * 20)), 1e-12);
  close(r.p, r.p0 * Math.exp(-3 * k), 1e-12);
  close(r.dpdz, k * r.p, 1e-12);
  const f = s.wavePorePressure({ H: 5, d: 20, T: 10, z: 3, thickness: 8 });
  close(f.ratio, Math.cosh(k * 5) / Math.cosh(k * 8), 1e-12);
  // no flow through the impervious base; thick bed tends to the exponential solution
  assert.ok(s.wavePorePressure({ H: 5, d: 20, T: 10, z: 8, thickness: 8 }).dpdz === 0);
  close(s.wavePorePressure({ H: 5, d: 20, T: 10, z: 3, thickness: 2000 }).ratio, r.ratio, 1e-9);
});

test('buried pipe wave uplift: twice the undisturbed gradient times the pipe volume; core attenuation by hand', () => {
  const o = { H: 6, d: 15, T: 9, z: 1.5, D: 0.8 };
  const w = s.wavePorePressure(o);
  const r = s.buriedPipeWaveUplift(o);
  close(r.F, 2 * (Math.PI * 0.64 / 4) * w.dpdz, 1e-12);
  // the factor 2: integrate the surface pressure of an impervious cylinder in a uniform gradient G,
  // p = 2 G a cos(theta) -> force per length = integral of p cos(theta) a dtheta = 2 G pi a^2
  const G = 1000, a = 0.4;
  let F = 0;
  const n = 100000;
  for (let i = 0; i < n; i++) { const th = (i + 0.5) * 2 * Math.PI / n; F += 2 * G * a * Math.cos(th) * Math.cos(th) * a * (2 * Math.PI / n); }
  close(F, 2 * G * Math.PI * a * a, 1e-9);
  close(s.coreAttenuation({ p0: 20e3, x: 10, delta: 0.5, Lprime: 60 }), 20e3 * Math.exp(-0.5 * 2 * Math.PI * 10 / 60), 1e-12);
  const b = s.burcharthDelta({ n: 0.36, Lp: 80, Hs: 4, b: 30 });
  close(b.delta, 0.014 * 0.6 * 6400 / 120, 1e-12);
  close(b.Lprime, 80 / Math.sqrt(1.4), 1e-12);
});
