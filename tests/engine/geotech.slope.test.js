import test from 'node:test';
import assert from 'node:assert/strict';
import * as s from '../../src/engine/geotech/slope.js';

const deg = Math.PI / 180;
const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);

test('infinite slope: dry sand tan(phi)/tan(beta); seepage parallel to the slope (gamma\'/gamma_sat) tan(phi)/tan(beta)', () => {
  close(s.infiniteSlope({ beta: 20 * deg, z: 3, gamma: 18e3, phi: 32 * deg }).fs, Math.tan(32 * deg) / Math.tan(20 * deg), 1e-12);
  const w = s.infiniteSlope({ beta: 20 * deg, z: 3, gamma: 20e3, phi: 32 * deg, hw: 3 });
  close(w.fs, ((20e3 - 9810) / 20e3) * Math.tan(32 * deg) / Math.tan(20 * deg), 1e-12);
  // ru form is identical for ru = gamma_w cos^2(beta) / gamma
  close(s.infiniteSlope({ beta: 20 * deg, z: 3, gamma: 20e3, phi: 32 * deg, ru: 9810 * Math.cos(20 * deg) ** 2 / 20e3 }).fs, w.fs, 1e-12);
});

test('infinite slope: c-phi hand calculation, undrained submarine slope, pseudo-static load', () => {
  const b = 25 * deg, z = 4, g = 19e3;
  const r = s.infiniteSlope({ beta: b, z, gamma: g, c: 5e3, phi: 28 * deg });
  close(r.fs, (5e3 + g * z * Math.cos(b) ** 2 * Math.tan(28 * deg)) / (g * z * Math.sin(b) * Math.cos(b)), 1e-12);
  const u = s.infiniteSlope({ beta: 5 * deg, z: 10, gamma: 6e3, mode: 'undrained', su: 12e3 });
  close(u.fs, 12e3 / (6e3 * 10 * Math.sin(5 * deg) * Math.cos(5 * deg)), 1e-12);
  // seismic, undrained, buoyant driving weight but total inertia
  const q = s.infiniteSlope({ beta: 5 * deg, z: 10, gamma: 6e3, gammaTotal: 16e3, mode: 'undrained', su: 12e3, kh: 0.05 });
  close(q.fs, 12e3 / (6e3 * 10 * Math.sin(5 * deg) * Math.cos(5 * deg) + 0.05 * 16e3 * 10 * Math.cos(5 * deg) ** 2), 1e-12);
  // dry cohesionless slope: yield coefficient kh = tan(phi - beta) gives FS = 1
  close(s.infiniteSlope({ beta: 20 * deg, z: 2, gamma: 18e3, phi: 35 * deg, kh: Math.tan(15 * deg) }).fs, 1, 1e-12);
});

test('Taylor (1937) stability number for phi = 0, beta = 60 deg: c/(gamma H FS) = 0.191 from the critical circle', () => {
  const H = 10, beta = 60 * deg, c = 20e3, gamma = 18e3;
  const run = H / Math.tan(beta);
  const surface = [[-40, 0], [0, 0], [run, -H], [run + 40, -H]];
  const r = s.criticalCircle({ surface, soil: { gamma, c, phi: 0 } });
  close(c / (gamma * H * r.fs), 0.191, 1.5e-2);
  // for phi = 0 the ordinary and Bishop methods coincide (moment equilibrium only)
  close(r.fsBishop, r.fsOrdinary, 1e-9);
  // the critical circle is a toe circle for beta > 53 deg
  assert.ok(Math.abs(r.xExit - run) < 0.05 * H, `exit at ${r.xExit}`);
});

test('Taylor (1937): vertical cut in phi = 0 soil, stability number 0.261', () => {
  const H = 6, c = 30e3, gamma = 19e3;
  const surface = [[-30, 0], [0, 0], [0.001, -H], [30, -H]];
  const r = s.criticalCircle({ surface, soil: { gamma, c, phi: 0 }, nSlices: 80 });
  close(c / (gamma * H * r.fs), 0.261, 2e-2);
});

test('Bishop simplified on a long shallow circle in dry cohesionless soil tends to the infinite-slope value tan(phi)/tan(beta)', () => {
  const beta = Math.atan(0.5);
  const Lh = 200;
  const surface = [[-50, 0], [0, 0], [Lh, -Lh * 0.5], [Lh + 50, -Lh * 0.5]];
  // circle through two points of the slope face 100 m apart with a sagitta of 0.5 m
  const x1 = 50, x2 = 150;
  const chord = (x2 - x1) / Math.cos(beta);
  const sag = 0.5;
  const R = chord * chord / (8 * sag) + sag / 2;
  const mx = 100, my = -50;
  const xc = mx + (R - sag) * Math.sin(beta);
  const yc = my + (R - sag) * Math.cos(beta);
  const r = s.circleStability({ surface, circle: { xc, yc, R }, soil: { gamma: 18e3, c: 0, phi: 30 * deg }, nSlices: 200 });
  const inf = Math.tan(30 * deg) / 0.5;
  close(r.fsBishop, inf, 5e-3);
  close(r.fsOrdinary, inf, 5e-3);
  assert.ok(r.fsBishop >= inf * (1 - 1e-6));
  assert.equal(r.direction, 1);
});

test('Bishop exceeds the ordinary method for frictional soil; pore pressure and seismic load reduce FS', () => {
  const H = 8;
  const surface = [[-30, 0], [0, 0], [16, -H], [50, -H]];
  const soil = { gamma: 19e3, c: 8e3, phi: 25 * deg };
  const circle = { xc: 12, yc: 10, R: 19 };
  const r = s.circleStability({ surface, circle, soil });
  assert.ok(r.fsBishop > r.fsOrdinary);
  assert.ok(s.circleStability({ surface, circle, soil, ru: 0.3 }).fsBishop < r.fsBishop);
  assert.ok(s.circleStability({ surface, circle, soil, kh: 0.1 }).fsBishop < r.fsBishop);
  // slice-count convergence
  close(s.circleStability({ surface, circle, soil, nSlices: 400 }).fsBishop, s.circleStability({ surface, circle, soil, nSlices: 200 }).fsBishop, 1e-4);
  // a circle that misses the ground returns null
  assert.equal(s.circleStability({ surface, circle: { xc: 0, yc: 50, R: 5 }, soil }), null);
  // critical circle is no higher than an arbitrary trial circle; weak layer lowers it further
  const crit = s.criticalCircle({ surface, soil });
  assert.ok(crit.fs <= r.fsBishop + 1e-9);
  const layered = s.criticalCircle({ surface, soil: [{ yTop: 100, gamma: 19e3, c: 8e3, phi: 25 * deg }, { yTop: -9, gamma: 17e3, c: 4e3, phi: 0 }] });
  assert.ok(layered.fs < crit.fs);
});

test('trench side slope: sand slope returns the infinite-slope value as the lower bound', () => {
  const r = s.trenchSideSlope({ depth: 3, slopeAngle: Math.atan(1 / 3), soil: { gamma: 9e3, c: 0, phi: 32 * deg } });
  close(r.fsInfinite, Math.tan(32 * deg) * 3, 1e-12);
  assert.ok(r.fs >= r.fsInfinite * 0.999 && r.fs < r.fsInfinite * 1.1, `fs ${r.fs}`);
});

test('Newmark sliding block: rectangular pulse d = 0.5 (A - ky) t0^2 A / ky (Newmark 1965)', () => {
  const A = 4, ky = 1, t0 = 0.5, dt = 0.001;
  const acc = [];
  for (let i = 0; i <= 4000; i++) acc.push(i * dt < t0 ? A : 0);
  const r = s.newmarkDisplacement({ acc, dt, ky });
  close(r.displacement, 0.5 * (A - ky) * t0 * t0 * A / ky, 5e-3);
  close(r.maxVelocity, (A - ky) * t0, 5e-3);
  close(r.slipTime, t0 * A / ky, 5e-3);
  // no slip below the yield acceleration
  assert.equal(s.newmarkDisplacement({ acc: acc.map((v) => v * 0.2), dt, ky }).displacement, 0);
});

test('Newmark sliding block: half-sine pulse against an independent fine integration; polarity option', () => {
  const A = 5, ky = 1.5, T = 1, dt = 0.01;
  const acc = [];
  for (let i = 0; i <= 300; i++) { const t = i * dt; acc.push(t < T / 2 ? A * Math.sin(2 * Math.PI * t / T) : 0); }
  // reference: explicit integration with a very small step
  let v = 0, d = 0;
  const h = 1e-6;
  for (let t = 0; t < 3; t += h) {
    const a = t < T / 2 ? A * Math.sin(2 * Math.PI * t / T) : 0;
    if (v > 0 || a > ky) { v += (a - ky) * h; if (v < 0) v = 0; d += v * h; }
  }
  close(s.newmarkDisplacement({ acc, dt, ky }).displacement, d, 5e-3);
  const neg = acc.map((x) => -x);
  assert.equal(s.newmarkDisplacement({ acc: neg, dt, ky }).displacement, 0);
  close(s.newmarkDisplacement({ acc: neg, dt, ky, bothPolarities: true }).displacement, d, 5e-3);
  close(s.yieldAcceleration(1.3, 20 * deg), 0.3 * 9.80665 * Math.sin(20 * deg), 1e-12);
});

test('submarine slope under waves: reduces to the gravity-only undrained infinite slope without waves; wave shear by hand', () => {
  const b = 3 * deg, k = 2 * Math.PI / 150;
  const calm = s.submarineSlopeWave({ beta: b, H: 0, k, d: 30, gammaSub: 6e3, su0: 2e3, suGrad: 1.2e3 });
  // FS(z) = (su0 + g z)/(gamma' z sin cos) decreases with depth: minimum at zMax
  const zMax = 150;
  close(calm.fsMin, (2e3 + 1.2e3 * zMax) / (6e3 * zMax * Math.sin(b) * Math.cos(b)), 1e-9);
  const storm = s.submarineSlopeWave({ beta: b, H: 10, k, d: 30, gammaSub: 6e3, su0: 2e3, suGrad: 1.2e3 });
  close(storm.p0, 10055 * 10 / (2 * Math.cosh(k * 30)), 1e-12);
  const i = 49;
  const z = storm.z[i];
  close(storm.fs[i], (2e3 + 1.2e3 * z) / (6e3 * z * Math.sin(b) * Math.cos(b) + storm.p0 * k * z * Math.exp(-k * z)), 1e-12);
  assert.ok(storm.fsMin < calm.fsMin);
});
