import test from 'node:test';
import assert from 'node:assert/strict';
import * as dr from '../../src/engine/hydro/drift.js';
import { jonswap, dispersion } from '../../src/engine/hydro/_wavehelp.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);
const rho = 1025; const g = 9.80665;

test('Maruo: full reflection in deep water gives ½ rho g a² per unit length; finite-depth factor 1 + 2kh/sinh 2kh', () => {
  close(dr.maruoDrift({ a: 2, R: 1 }), 0.5 * rho * g * 4, 1e-12);
  close(dr.maruoDrift({ a: 2, R: 0.5 }), 0.5 * rho * g * 4 * 0.25, 1e-12);
  const kh = dispersion((2 * Math.PI) / 12, 15) * 15;
  close(dr.maruoDrift({ a: 1, R: 1, T: 12, depth: 15 }), 0.5 * rho * g * (1 + (2 * kh) / Math.sinh(2 * kh)), 1e-12);
  // shallow-water limit of the factor is 2
  close(dr.maruoDrift({ a: 1, R: 1, T: 200, depth: 2 }) / (0.5 * rho * g), 2, 2e-3);
});

test('Vertical cylinder, short-wave limit: numerical integral of ½ rho g a² cos³(theta) r equals (2/3) rho g a² r', () => {
  const r = 6; const a = 1.5; const n = 20000;
  let s = 0;
  for (let i = 0; i < n; i++) { const th = -Math.PI / 2 + ((i + 0.5) * Math.PI) / n; s += 0.5 * rho * g * a * a * Math.cos(th) ** 3 * r * (Math.PI / n); }
  const c = dr.cylinderDriftShortWave({ a, D: 2 * r });
  close(c.force, s, 1e-7); close(c.coefficient, (2 / 3) * rho * g * r, 1e-12);
});

test('Mean drift in irregular waves: constant coefficient D0 gives 2 D0 m0 = D0 Hs²/8', () => {
  const omega = Array.from({ length: 4000 }, (_, i) => 0.05 + (i * 5.95) / 3999);
  const S = omega.map((w) => jonswap(w, { Hs: 5, Tp: 11, gamma: 1 }));
  close(dr.meanDriftForce({ omega, S, D: omega.map(() => 3e4) }), (3e4 * 25) / 8, 3e-3);
});

test('Newman slow-drift spectrum: box spectrum S0 on [w1, w2] with constant D0 gives 8 S0² D0² (w2 - w1 - mu)', () => {
  const omega = Array.from({ length: 6001 }, (_, i) => 0.2 + i * 2e-4);   // 0.2 … 1.4
  const S0 = 2.5; const D0 = 1e4;
  const S = omega.map((w) => (w >= 0.5 && w <= 1.0 ? S0 : 0));
  const mu = [0, 0.05, 0.2, 0.6];
  const SF = dr.newmanDriftSpectrum({ omega, S, D: omega.map(() => D0), mu });
  mu.forEach((m, i) => {
    const exact = 8 * S0 * S0 * D0 * D0 * Math.max(0.5 - m, 0);
    assert.ok(Math.abs(SF[i] - exact) <= 2e-3 * 8 * S0 * S0 * D0 * D0 * 0.5, `mu ${m}: ${SF[i]} vs ${exact}`);
  });
});

test('Slow-drift response: white-noise force spectrum gives sigma² = pi S_F / (2 c k); extreme = mean + sigma sqrt(2 ln N)', () => {
  const m = 1.2e8; const k = 2.5e5; const zeta = 0.08;
  const mu = Array.from({ length: 4001 }, (_, i) => 1e-5 + i * 2.5e-4);    // to 1 rad/s, wn = 0.0456
  const S0 = 4e10;
  const r = dr.slowDriftResponse({ mass: m, stiffness: k, dampingRatio: zeta, mu, SF: mu.map(() => S0), meanForce: 5e5, duration: 10800 });
  const c = 2 * zeta * Math.sqrt(k * m);
  close(r.sigmaWhiteNoise, Math.sqrt((Math.PI * S0) / (2 * c * k)), 1e-12);
  close(r.sigma, r.sigmaWhiteNoise, 1e-2);
  close(r.naturalPeriod, 2 * Math.PI * Math.sqrt(m / k), 1e-12);
  close(r.meanOffset, 2, 1e-12);
  close(r.extreme, 2 + r.sigma * Math.sqrt(2 * Math.log(10800 / r.naturalPeriod)), 1e-12);
});

test('OCIMF-form and projected-area forces (hand calculations)', () => {
  const w = dr.ocimfForce({ kind: 'wind', V: 30, angleDeg: 60, Cx: 0.4, Cy: (a) => 0.01 * a, Cxy: { angleDeg: [0, 90, 180], value: [0, 0.1, 0] }, AT: 800, AL: 3000, L: 250, rho: 1.25 });
  const q = 0.5 * 1.25 * 900;
  close(w.Fx, q * 0.4 * 800, 1e-12); close(w.Fy, q * 0.6 * 3000, 1e-12); close(w.Mz, q * (0.1 * 60 / 90) * 3000 * 250, 1e-12);
  const c = dr.ocimfForce({ kind: 'current', V: 1.5, angleDeg: 90, Cx: 0.05, Cy: 0.6, Cxy: 0.02, L: 250, T: 15 });
  close(c.Fy, 0.5 * 1025 * 2.25 * 0.6 * 250 * 15, 1e-12); close(c.Mz, 0.5 * 1025 * 2.25 * 0.02 * 250 * 250 * 15, 1e-12);
  const p = dr.projectedAreaForce({ V: 20, directionDeg: 30, Afront: 500, Aside: 2000, Cfront: 1.0, Cside: 1.1, xSide: 10, rho: 1.226 });
  close(p.Fx, 0.5 * 1.226 * 400 * 500 * Math.cos(Math.PI / 6) ** 2, 1e-12);
  close(p.Fy, 0.5 * 1.226 * 400 * 1.1 * 2000 * 0.25, 1e-12); close(p.Mz, p.Fy * 10, 1e-12);
});

test('Total mean environmental force versus heading: collinear loads are largest beam-on; end-on value by hand', () => {
  const wind = (rel) => dr.projectedAreaForce({ V: 25, directionDeg: rel, Afront: 600, Aside: 2500, Cfront: 1, Cside: 1, rho: 1.226 });
  const cur = (rel) => dr.projectedAreaForce({ V: 1, directionDeg: rel, Afront: 300, Aside: 3000, Cfront: 1, Cside: 1, rho: 1025 });
  const r = dr.meanEnvironmentalForce({ headingsDeg: Array.from({ length: 36 }, (_, i) => 10 * i), loads: [{ name: 'wind', directionDeg: 90, force: wind }, { name: 'current', directionDeg: 90, force: cur }] });
  assert.ok([0, 180].includes(r.worst.headingDeg));     // environment towards 90°, heading 0° -> beam-on
  close(r.worst.F, 0.5 * 1.226 * 625 * 2500 + 0.5 * 1025 * 3000, 1e-9);
  // end-on (heading 90°): only the frontal areas are loaded
  close(r.F[9], 0.5 * 1.226 * 625 * 600 + 0.5 * 1025 * 300, 1e-9);
  assert.ok(r.best.F <= r.F[9] && Math.abs(r.Fy[9]) < 1e-6 * r.F[9]);
});
