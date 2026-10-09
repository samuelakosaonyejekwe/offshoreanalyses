import test from 'node:test';
import assert from 'node:assert/strict';
import * as W from '../../src/engine/metocean/waves.js';

const g = W.G;
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('dispersion: deep-water wavelength L = g T² / 2π', () => {
  for (const T of [4, 8, 12, 16]) {
    const L0 = (g * T * T) / (2 * Math.PI);
    assert.ok(rel(W.wavelength(T, Infinity), L0) < 1e-14);
    assert.ok(rel(W.wavelength(T, 5000), L0) < 1e-12);
  }
});

test('dispersion: shallow-water celerity c = sqrt(g d)', () => {
  const d = 5; const T = 900;
  // kd = 0.005: c/sqrt(gd) = sqrt(tanh(kd)/kd) = 1 - (kd)²/6 to 4e-6
  assert.ok(rel(W.celerity(T, d), Math.sqrt(g * d)) < 1e-5);
  assert.ok(rel(W.groupVelocity(T, d), Math.sqrt(g * d)) < 3e-5);
});

test('dispersion: residual of ω² = g k tanh(kd) below 1e-12 over a wide range', () => {
  for (const T of [1, 3, 6, 10, 15, 25, 100, 1000]) {
    for (const d of [0.5, 5, 30, 150, 1000, 4000]) {
      const k = W.wavenumber({ T, d });
      const om = (2 * Math.PI) / T;
      assert.ok(Math.abs(g * k * Math.tanh(k * d) - om * om) / (om * om) < 1e-12, `T=${T} d=${d}`);
    }
  }
});

test('dispersion with current: Doppler-shifted relation satisfied, following current lengthens the wave', () => {
  const T = 8; const d = 20;
  const k0 = W.wavenumber({ T, d });
  for (const U of [-1.5, -0.5, 0.5, 2]) {
    const k = W.wavenumber({ T, d, U });
    const om = (2 * Math.PI) / T;
    assert.ok(Math.abs(g * k * Math.tanh(k * d) - (om - k * U) ** 2) / (om * om) < 1e-12);
    assert.ok(om - k * U > 0);
    assert.equal(k < k0, U > 0);
  }
});

test('dispersion with current: deep-water blocking at U = -g T / (8π) (cg_rel = -U at k = g/4U²)', () => {
  const T = 4;
  const Ublock = -(g * T) / (8 * Math.PI); // -1.561 m/s
  assert.ok(W.wavenumber({ T, d: 500, U: 0.98 * Ublock }) > 0);
  assert.throws(() => W.wavenumber({ T, d: 500, U: 1.02 * Ublock }), RangeError);
  // just inside the limit the solution approaches the stopping wavenumber 4 ω²/g
  const k = W.wavenumber({ T, d: 5000, U: 0.999999 * Ublock });
  assert.ok(rel(k, (4 * (2 * Math.PI / T) ** 2) / g) < 5e-3);
});

test('wave properties: deep-water energy flux P = ρ g² H² T / (32π); n = 1/2 deep, 1 shallow', () => {
  const H = 2; const T = 9; const rho = 1025;
  const p = W.waveProperties({ T, d: Infinity, H, rho });
  assert.ok(rel(p.P, (rho * g * g * H * H * T) / (32 * Math.PI)) < 1e-13);
  assert.equal(p.regime, 'deep');
  assert.ok(Math.abs(p.n - 0.5) < 1e-15);
  assert.ok(rel(W.energyFlux({ H, T, d: Infinity, rho }), p.P) < 1e-14);
  const s = W.waveProperties({ T: 600, d: 3, H: 0.1 });
  assert.equal(s.regime, 'shallow');
  assert.ok(Math.abs(s.n - 1) < 1e-3);
  assert.ok(rel(W.energyDensity(2, 1000, 10), 5000) < 1e-15);
  assert.ok(rel(W.ursellNumber(2, 100, 10), 20) < 1e-15);
  assert.ok(rel(W.steepness(2, 100), 0.02) < 1e-15);
});

test('Airy kinematics satisfy continuity, irrotationality (Laplace) and the bed condition numerically', () => {
  const o = { H: 2, T: 8, d: 25 };
  const h = 1e-3;
  for (const [x, z, t] of [[3, -2, 1], [40, -10, 3.3], [-17, -24, 6], [5, -0.5, 0]]) {
    const f = (dx, dz, dt = 0) => W.airyKinematics({ ...o, x: x + dx, z: z + dz, t: t + dt, clipAboveSurface: false });
    const ux = (f(h, 0).u - f(-h, 0).u) / (2 * h);
    const wz = (f(0, h).w - f(0, -h).w) / (2 * h);
    const uz = (f(0, h).u - f(0, -h).u) / (2 * h);
    const wx = (f(h, 0).w - f(-h, 0).w) / (2 * h);
    const scale = Math.abs(f(0, 0).u) + Math.abs(f(0, 0).w);
    assert.ok(Math.abs(ux + wz) < 1e-7 * scale);
    assert.ok(Math.abs(uz - wx) < 1e-7 * scale);
    // local accelerations are the time derivatives of the velocities
    const c0 = f(0, 0);
    assert.ok(Math.abs((f(0, 0, h).u - f(0, 0, -h).u) / (2 * h) - c0.ax) < 1e-6);
    assert.ok(Math.abs((f(0, 0, h).w - f(0, 0, -h).w) / (2 * h) - c0.az) < 1e-6);
  }
  assert.ok(Math.abs(W.airyKinematics({ ...o, x: 7, z: -25, t: 2 }).w) < 1e-15);
});

test('Airy kinematics satisfy the linearised free-surface conditions at z = 0', () => {
  const o = { H: 1.5, T: 7, d: 12, rho: 1025 };
  const h = 1e-4;
  for (const [x, t] of [[0, 0], [11, 2.2], [-30, 5.1]]) {
    const f = (dt) => W.airyKinematics({ ...o, x, z: 0, t: t + dt, clipAboveSurface: false });
    const c0 = f(0);
    // kinematic: ∂η/∂t = w ; dynamic: p_dyn = ρ g η (i.e. ∂φ/∂t + g η = 0)
    assert.ok(Math.abs((f(h).eta - f(-h).eta) / (2 * h) - c0.w) < 1e-7);
    assert.ok(Math.abs(c0.pDyn - 1025 * g * c0.eta) < 1e-9 * 1025 * g);
    // combined condition ∂²φ/∂t² + g ∂φ/∂z = 0  <=>  ∂u/∂t... expressed with velocities: ∂w/∂t = -g ∂η/∂x·(-1)
    const ex = (W.airyKinematics({ ...o, x: x + h, z: 0, t }).eta - W.airyKinematics({ ...o, x: x - h, z: 0, t }).eta) / (2 * h);
    assert.ok(Math.abs(c0.ax + g * ex) < 1e-6); // linear x-momentum at the surface: ∂u/∂t = -g ∂η/∂x
  }
});

test('Airy: bed orbital velocity amplitude π H / (T sinh kd) and deep-water surface value π H / T', () => {
  const H = 3; const T = 10; const d = 20;
  const k = W.wavenumber({ T, d });
  assert.ok(rel(W.airyKinematics({ H, T, d, z: -d }).u, (Math.PI * H) / (T * Math.sinh(k * d))) < 1e-13);
  assert.ok(rel(W.airyKinematics({ H, T, d: Infinity, z: 0, clipAboveSurface: false }).u, (Math.PI * H) / T) < 1e-14);
  assert.ok(rel(W.airyKinematics({ H, T, d: Infinity, z: -10 }).u, ((Math.PI * H) / T) * Math.exp((-10 * 4 * Math.PI ** 2) / (g * T * T))) < 1e-13);
});

test('Airy stretching options above still-water level', () => {
  const o = { H: 4, T: 9, d: 30 };
  const k = W.wavenumber(o);
  const u0 = W.airyKinematics({ ...o, z: 0 }).u; // crest phase, SWL value
  const eta = 2;
  // Wheeler: the SWL value is moved to the crest, the bed value stays at the bed
  assert.ok(rel(W.airyKinematics({ ...o, z: eta, stretching: 'wheeler' }).u, u0) < 1e-13);
  assert.ok(rel(W.airyKinematics({ ...o, z: -30, stretching: 'wheeler' }).u, W.airyKinematics({ ...o, z: -30 }).u) < 1e-13);
  // Wheeler mid-column: z' = (z - η) d / (d + η)
  const zp = ((-5 - eta) * 30) / 32;
  assert.ok(rel(W.airyKinematics({ ...o, z: -5, stretching: 'wheeler' }).u, W.airyKinematics({ ...o, z: zp }).u) < 1e-13);
  // vertical: constant above SWL
  assert.ok(rel(W.airyKinematics({ ...o, z: 1.5, stretching: 'vertical' }).u, u0) < 1e-14);
  // linear: u(0) + z du/dz(0) = u0 (1 + z k tanh kd)
  assert.ok(rel(W.airyKinematics({ ...o, z: 1.5, stretching: 'linear' }).u, u0 * (1 + 1.5 * k * Math.tanh(k * 30))) < 1e-13);
  // plain extrapolation grows as cosh
  assert.ok(rel(W.airyKinematics({ ...o, z: 1.5, stretching: 'none' }).u, (u0 * Math.cosh(k * 31.5)) / Math.cosh(k * 30)) < 1e-13);
  // above the instantaneous surface: dry
  const dry = W.airyKinematics({ ...o, z: 2.5 });
  assert.equal(dry.wet, false); assert.equal(dry.u, 0);
  assert.throws(() => W.airyKinematics({ ...o, z: -31 }), RangeError);
  assert.throws(() => W.airyKinematics({ ...o, z: 0, stretching: 'bogus' }), RangeError);
});

test('Stokes 2nd order: deep-water surface η = a cos θ + ½ k a² cos 2θ; first-order part equals Airy', () => {
  const H = 3; const T = 8;
  const k = (4 * Math.PI ** 2) / (g * T * T);
  const a = H / 2;
  const s = W.stokes2Kinematics({ H, T, d: Infinity, x: 0, z: 0, t: 0 });
  assert.ok(rel(s.eta, a + 0.5 * k * a * a) < 1e-13);
  const tr = W.stokes2Kinematics({ H, T, d: Infinity, x: Math.PI / k, z: -3, t: 0 });
  assert.ok(rel(tr.eta, -a + 0.5 * k * a * a) < 1e-13);
  assert.ok(rel(s.eta - tr.eta, H) < 1e-13);
  // deep water: no second-order velocity potential, Stokes drift a² ω k e^{2kz}
  const airy = W.airyKinematics({ H, T, d: Infinity, x: 5, z: -3, t: 1 });
  const st = W.stokes2Kinematics({ H, T, d: Infinity, x: 5, z: -3, t: 1 });
  assert.ok(rel(st.u, airy.u) < 1e-13);
  assert.ok(rel(st.stokesDrift, a * a * ((2 * Math.PI) / T) * k * Math.exp(-6 * k)) < 1e-13);
});

test('Stokes 2nd order: finite-depth hand calculation and consistency with the Bernoulli equation', () => {
  const H = 2; const T = 8; const d = 10; const rho = 1025;
  const k = W.wavenumber({ T, d }); const om = (2 * Math.PI) / T; const a = H / 2;
  const sh = Math.sinh(k * d); const ch = Math.cosh(k * d);
  const s = W.stokes2Kinematics({ H, T, d, x: 0, z: -4, t: 0, rho });
  assert.ok(rel(s.eta, a + (k * a * a * ch * (2 + Math.cosh(2 * k * d))) / (4 * sh ** 3)) < 1e-13);
  assert.ok(rel(s.u, (a * om * Math.cosh(6 * k)) / sh + (0.75 * a * a * om * k * Math.cosh(12 * k)) / sh ** 4) < 1e-13);
  // continuity and irrotationality of the full second-order field
  const h = 1e-3;
  const f = (dx, dz) => W.stokes2Kinematics({ H, T, d, x: 13 + dx, z: -4 + dz, t: 2 });
  assert.ok(Math.abs((f(h, 0).u - f(-h, 0).u + f(0, h).w - f(0, -h).w) / (2 * h)) < 1e-7);
  assert.ok(Math.abs((f(0, h).u - f(0, -h).u - f(h, 0).w + f(-h, 0).w) / (2 * h)) < 1e-7);
  // second-order dynamic pressure equals -ρ ∂φ2/∂t - ½ρ(u1² + w1²) + ½ρ(u1² + w1²)|bed-mean constant:
  // check the oscillating cos 2θ part against the difference of two phases a quarter second-harmonic period apart
  const pA = W.stokes2Kinematics({ H, T, d, x: 0, z: -4, t: 0, rho });
  const pB = W.stokes2Kinematics({ H, T, d, x: Math.PI / (2 * k), z: -4, t: 0, rho }); // θ = π/2: cos θ = 0, cos 2θ = -1
  const pC = W.stokes2Kinematics({ H, T, d, x: Math.PI / k, z: -4, t: 0, rho }); // θ = π
  const cos2Amp = 0.5 * (0.5 * (pA.pDyn + pC.pDyn) - pB.pDyn);
  const expected = rho * (0.75 * a * a * om * om * Math.cosh(12 * k) / sh ** 4 - 0.25 * a * a * om * om / sh ** 2);
  assert.ok(rel(cos2Amp, expected) < 1e-12);
});

test('Stokes 5th order (Fenton 1985): free-surface residuals decay as the sixth power of steepness', () => {
  const resid = (kd, eps) => {
    const d = 20; const k = kd / d; const H = (2 * eps) / k; const L = (2 * Math.PI) / k;
    const s = W.stokes5Wave({ H, d, L });
    let bmin = Infinity; let bmax = -Infinity; let kin = 0;
    for (let i = 0; i <= 40; i++) {
      const x = (L * i) / 80;
      const e = W.stokes5Kinematics(s, { x, z: -1 }).eta;
      const q = W.stokes5Kinematics(s, { x, z: e });
      const h = 1e-5 * L;
      const ex = (W.stokes5Kinematics(s, { x: x + h, z: -1 }).eta - W.stokes5Kinematics(s, { x: x - h, z: -1 }).eta) / (2 * h);
      const b = 0.5 * ((q.u - s.c) ** 2 + q.w ** 2) + g * e;
      bmin = Math.min(bmin, b); bmax = Math.max(bmax, b);
      kin = Math.max(kin, Math.abs(q.w - (q.u - s.c) * ex));
    }
    return { dyn: ((bmax - bmin) * k) / g, kin: kin / s.c };
  };
  for (const kd of [0.8, 1.5, 3, 40]) {
    const a = resid(kd, 0.04); const b = resid(kd, 0.08);
    const pd = Math.log2(b.dyn / a.dyn); const pk = Math.log2(b.kin / a.kin);
    assert.ok(pd > 5.7 && pd < 6.4, `dynamic order ${pd} at kd=${kd}`);
    assert.ok(pk > 5.7 && pk < 6.4, `kinematic order ${pk} at kd=${kd}`);
    assert.ok(a.dyn < 1e-5 && a.kin < 1e-5);
  }
});

test('Stokes 5th order: deep-water celerity c = sqrt(g/k)(1 + ε²/2 + ε⁴/8) and surface coefficients', () => {
  const L = 150; const k = (2 * Math.PI) / L; const eps = 0.25; const H = (2 * eps) / k;
  const s = W.stokes5Wave({ H, d: Infinity, L });
  assert.ok(rel(s.c, Math.sqrt(g / k) * (1 + eps ** 2 / 2 + eps ** 4 / 8)) < 1e-13);
  // kη = ε cos θ + ½ε² cos 2θ + ⅜ε³(cos 3θ − cos θ) + ⅓ε⁴(cos 2θ + cos 4θ) + ε⁵(…): crest value
  const B53 = 99 / 128; const B55 = 125 / 384;
  const crest = eps + eps ** 2 / 2 + (eps ** 4 * 2) / 3 + 0 * eps ** 3 + eps ** 5 * (-(B53 + B55) + B53 + B55);
  assert.ok(rel(s.crest * k, crest) < 1e-13);
  assert.ok(rel(s.crest - s.trough, H) < 1e-13);
  const c = W.stokes5Coefficients(Infinity);
  assert.deepEqual([c.C0, c.C2, c.C4, c.B[2][2], c.B[3][1], c.B[4][2], c.B[4][4]], [1, 0.5, 0.125, 0.5, -0.375, 1 / 3, 1 / 3]);
});

test('Stokes 5th order: period closure with Eulerian and Stokes currents; small-wave limit is Airy', () => {
  const s = W.stokes5Wave({ H: 5, d: 30, T: 10, currentEulerian: 0.8 });
  assert.ok(rel(s.T, 10) < 1e-12);
  assert.ok(rel(s.currentEulerian, 0.8) < 1e-12);
  assert.ok(s.currentStokes > 0.8);
  const st = W.stokes5Wave({ H: 5, d: 30, T: 10, currentStokes: 0 });
  assert.ok(rel(st.T, 10) < 1e-12);
  assert.ok(Math.abs(st.currentStokes) < 1e-12);
  assert.ok(st.currentEulerian < 0); // return flow balancing the Stokes drift
  // time-mean of u at a point below the trough equals the Eulerian current
  let m = 0; const n = 64;
  for (let i = 0; i < n; i++) m += W.stokes5Kinematics(s, { x: 0, z: -10, t: (i * s.T) / n }).u;
  assert.ok(Math.abs(m / n - 0.8) < 1e-12);
  const small = W.stokes5Wave({ H: 0.002, d: 30, T: 10 });
  assert.ok(rel(small.L, W.wavelength(10, 30)) < 1e-7);
  const a = W.airyKinematics({ H: 0.002, T: 10, d: 30, x: 3, z: -7, t: 1 });
  const f = W.stokes5Kinematics(small, { x: 3, z: -7, t: 1 });
  // differences are of relative order ε = kH/2 = 4.6e-5 (second harmonic)
  assert.ok(rel(f.u, a.u) < 1e-4 && rel(f.az, a.az) < 1e-4 && rel(f.pDyn, a.pDyn) < 1e-4);
});

test('elliptic integrals and Jacobi functions: tabulated values and identities (Abramowitz & Stegun 17.3, 16.9)', () => {
  const { K, E } = W.ellipticKE(0.5);
  assert.ok(Math.abs(K - 1.8540746773013719) < 1e-14);
  assert.ok(Math.abs(E - 1.3506438810476755) < 1e-14);
  assert.ok(Math.abs(W.ellipticKE(0).K - Math.PI / 2) < 1e-15);
  // Legendre relation E K' + E' K − K K' = π/2
  for (const m of [0.1, 0.37, 0.9, 0.999]) {
    const a = W.ellipticKE(m); const b = W.ellipticKE(1 - m);
    assert.ok(Math.abs(a.E * b.K + b.E * a.K - a.K * b.K - Math.PI / 2) < 1e-13);
    for (const u of [0.2, 1.1, 3.7, -2.4]) {
      const j = W.jacobiElliptic(u, m);
      assert.ok(Math.abs(j.sn ** 2 + j.cn ** 2 - 1) < 1e-14);
      assert.ok(Math.abs(j.dn ** 2 + m * j.sn ** 2 - 1) < 1e-14);
      // d(sn)/du = cn dn
      const h = 1e-5;
      const ds = (W.jacobiElliptic(u + h, m).sn - W.jacobiElliptic(u - h, m).sn) / (2 * h);
      assert.ok(Math.abs(ds - j.cn * j.dn) < 1e-9);
    }
    const q = W.jacobiElliptic(a.K, m);
    assert.ok(Math.abs(q.sn - 1) < 1e-13 && Math.abs(q.cn) < 1e-7);
  }
  const s = W.jacobiElliptic(0.7, 1);
  assert.ok(Math.abs(s.sn - Math.tanh(0.7)) < 1e-15 && Math.abs(s.cn - 1 / Math.cosh(0.7)) < 1e-15);
  assert.ok(Math.abs(W.jacobiElliptic(0.7, 0).sn - Math.sin(0.7)) < 1e-15);
});

test('cnoidal wave: closure relations, zero-mean surface, height, continuity and the solitary limit', () => {
  const H = 1; const d = 5; const T = 15;
  const c = W.cnoidalWave({ H, d, T });
  assert.ok(rel(c.T, T) < 1e-12);
  assert.ok(rel(c.L, c.c * T) < 1e-12);
  assert.ok(rel(c.ursell, (16 * c.m * c.K ** 2) / 3) < 1e-12); // H L²/d³ = 16 m K²/3
  assert.ok(rel(c.crest - c.trough, H) < 1e-13);
  // mean elevation over a wavelength is zero
  let mean = 0; const n = 2000;
  for (let i = 0; i < n; i++) mean += W.cnoidalKinematics(c, { x: (c.L * (i + 0.5)) / n, z: -1 }).eta;
  assert.ok(Math.abs(mean / n) < 1e-9);
  assert.ok(rel(W.cnoidalKinematics(c, { x: 0, z: -1 }).eta, c.crest) < 1e-13);
  assert.ok(rel(W.cnoidalKinematics(c, { x: c.L / 2, z: -1 }).eta, c.trough) < 1e-10);
  // continuity u_x + w_z = 0 and ax = ∂u/∂t
  const h = 1e-3;
  const f = (dx, dz, dt = 0) => W.cnoidalKinematics(c, { x: 9 + dx, z: -2 + dz, t: 0.4 + dt });
  assert.ok(Math.abs((f(h, 0).u - f(-h, 0).u + f(0, h).w - f(0, -h).w) / (2 * h)) < 1e-8);
  assert.ok(Math.abs((f(0, 0, h).u - f(0, 0, -h).u) / (2 * h) - f(0, 0).ax) < 1e-6);
  assert.ok(Math.abs((f(0, 0, h).w - f(0, 0, -h).w) / (2 * h) - f(0, 0).az) < 1e-6);
  // solving by wavelength returns the same wave
  assert.ok(rel(W.cnoidalWave({ H, d, L: c.L }).m1, c.m1) < 1e-9);
  // long-period limit: first-order solitary wave, c = sqrt(gd)(1 + H/2d), crest -> H, profile -> H sech²
  const s = W.cnoidalWave({ H, d, T: 400 });
  assert.ok(rel(s.c, Math.sqrt(g * d) * (1 + H / (2 * d))) < 2e-2);
  assert.ok(s.crest > 0.97 * H);
  const beta = Math.sqrt((3 * H) / (4 * d ** 3));
  assert.ok(Math.abs(W.cnoidalKinematics(s, { x: 6, z: -1 }).eta - s.trough - H / Math.cosh(beta * 6) ** 2) < 1e-6);
  assert.throws(() => W.cnoidalWave({ H: 1, d: 50, T: 5 }), RangeError);
});

test('cnoidal wave: small-m limit recovers the long-wave dispersion c = sqrt(gd)(1 − (kd)²/6)', () => {
  const d = 10; const L = 150; const H = 0.02; // Ursell = 0.45
  const c = W.cnoidalWave({ H, d, L });
  const kd = (2 * Math.PI * d) / L;
  assert.ok(rel(c.c, Math.sqrt(g * d) * (1 - kd ** 2 / 6)) < 2e-3);
});

test('solitary wave: celerity sqrt(g(d+H)), volume 4 d sqrt(dH/3) by quadrature, continuity', () => {
  const H = 1.2; const d = 6;
  const s0 = W.solitaryWave({ H, d, x: 0, z: 0 });
  assert.ok(rel(s0.c, Math.sqrt(g * (d + H))) < 1e-15);
  assert.ok(rel(s0.eta, H) < 1e-15);
  assert.ok(rel(s0.u, (s0.c * H) / (d + H)) < 1e-15);
  let vol = 0; const n = 20000; const X = 400;
  for (let i = 0; i < n; i++) vol += (W.solitaryWave({ H, d, x: -X + (2 * X * (i + 0.5)) / n, z: -d }).eta * 2 * X) / n;
  assert.ok(rel(vol, s0.volume) < 1e-9);
  assert.ok(rel(s0.volume, 4 * d * Math.sqrt((d * H) / 3)) < 1e-15);
  const h = 1e-3;
  const f = (dx, dz) => W.solitaryWave({ H, d, x: 4 + dx, z: -2 + dz, t: 0.3 });
  assert.ok(Math.abs((f(h, 0).u - f(-h, 0).u + f(0, h).w - f(0, -h).w) / (2 * h)) < 1e-8);
});

test('breaking limits: Miche 0.142 L in deep water and 0.142·2π d ≈ 0.89 d in shallow water; McCowan 0.78 d; Goda', () => {
  assert.ok(rel(W.breakingHeightMiche(100, Infinity), 14.2) < 1e-15);
  assert.ok(rel(W.breakingHeightMiche(1e5, 10), 0.142 * 2 * Math.PI * 10) < 1e-6);
  assert.ok(rel(W.breakingHeightMcCowan(10), 7.8) < 1e-15);
  // Goda: deep water -> A L0; shallow flat bed -> A·1.5π d = 0.80 d
  const T = 10; const L0 = (g * T * T) / (2 * Math.PI);
  assert.ok(rel(W.breakingHeightGoda({ T, d: 5000 }), 0.17 * L0) < 1e-12);
  assert.ok(rel(W.breakingHeightGoda({ T: 30, d: 1 }), 0.17 * 1.5 * Math.PI * 1) < 2e-3);
  // hand calculation with slope 1/20: exponent = 1.5 π d/L0 (1 + 15·(0.05)^(4/3))
  const ex = ((1.5 * Math.PI * 8) / L0) * (1 + 15 * 0.05 ** (4 / 3));
  assert.ok(rel(W.breakingHeightGoda({ T, d: 8, slope: 0.05 }), 0.17 * L0 * (1 - Math.exp(-ex))) < 1e-13);
  const b = W.breakingLimits({ T: 10, d: 8 });
  assert.equal(b.governing, Math.min(b.miche, b.mccowan));
});

test('wave-theory selector returns the expected regions (approximate chart boundaries)', () => {
  assert.equal(W.selectWaveTheory({ H: 0.2, T: 8, d: 100 }).theory, 'airy');
  assert.equal(W.selectWaveTheory({ H: 8, T: 10, d: 100 }).theory, 'stokes5');
  assert.equal(W.selectWaveTheory({ H: 1.5, T: 16, d: 5 }).theory, 'cnoidal');
  assert.equal(W.selectWaveTheory({ H: 3.7, T: 16, d: 5 }).theory, 'streamFunction');
  assert.equal(W.selectWaveTheory({ H: 6, T: 16, d: 5 }).theory, 'breaking');
  const r = W.selectWaveTheory({ H: 21, T: 10, d: 100 });
  assert.equal(r.theory, 'streamFunction');
  assert.ok(r.streamFunctionOrder >= 9 && r.approximate === true);
  assert.ok(rel(r.HgT2, 21 / (g * 100)) < 1e-14);
});

test('shoaling and refraction: Snell law, Ks = sqrt(c0 / 2cg) from deep water, Green\'s law in shallow water', () => {
  const T = 10; const d = 6;
  const r = W.shoalingRefraction({ T, d, theta0: (30 * Math.PI) / 180, H0: 2 });
  const c0 = (g * T) / (2 * Math.PI);
  assert.ok(rel(Math.sin(r.theta) / r.c, 0.5 / c0) < 1e-13);
  assert.ok(rel(r.Ks, Math.sqrt(c0 / (2 * r.cg))) < 1e-13);
  assert.ok(rel(r.Kr, Math.sqrt(Math.cos(Math.PI / 6) / Math.cos(r.theta))) < 1e-14);
  assert.ok(rel(r.H, 2 * r.Ks * r.Kr) < 1e-14);
  assert.ok(Math.abs(W.shoalingRefraction({ T, d }).Kr - 1) < 1e-15);
  // between two shallow depths H ∝ d^(-1/4)
  const g1 = W.shoalingRefraction({ T: 300, d: 2, d0: 8 });
  assert.ok(rel(g1.Ks, 4 ** 0.25) < 1e-3);
  // minimum shoaling coefficient 0.913 at kd = 1.2 (Dean & Dalrymple Fig. 4.9)
  let mn = Infinity;
  for (let dd = 5; dd < 80; dd += 0.05) mn = Math.min(mn, W.shoalingRefraction({ T, d: dd }).Ks);
  assert.ok(Math.abs(mn - 0.913) < 1e-3);
});
