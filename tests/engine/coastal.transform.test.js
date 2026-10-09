import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  linearWave, shoalingCoefficient, refractionSnell, transformProfile, breakingMiche, breakingMcCowan,
  breakingWeggel, breakingGoda, godaSurfZoneHeights, battjesStiveGamma, battjesJanssenQb, energyFluxBalance,
  fresnel, diffractionSemiInfinite, diffractionGap, reflectionSeeligAhrens, waveCurrent,
} from '../../src/engine/coastal/transform.js';

const g = 9.80665;
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const abs = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b}`);

test('linear dispersion: residual, deep-water and shallow-water limits', () => {
  for (const [T, h] of [[8, 10], [12, 3], [5, 200], [15, 0.5]]) {
    const w = linearWave({ T, h });
    rel(g * w.k * Math.tanh(w.k * h), (2 * Math.PI / T) ** 2, 1e-12);
  }
  const deep = linearWave({ T: 6, h: 500 });
  rel(deep.L, (g * 36) / (2 * Math.PI), 1e-10);
  rel(deep.n, 0.5, 1e-10);
  const shallow = linearWave({ T: 60, h: 1 });
  rel(shallow.C, Math.sqrt(g * 1), 2e-3);
  rel(shallow.n, 1, 5e-3);
});

test('shoaling coefficient: minimum 0.913 at kh = 1.20 (h/L0 = 0.159) and Green\'s law in shallow water', () => {
  const T = 10;
  const L0 = (g * T * T) / (2 * Math.PI);
  let min = 9;
  let at = 0;
  for (let r = 0.05; r < 0.4; r += 0.0005) {
    const ks = shoalingCoefficient({ T, h: r * L0 });
    if (ks < min) {
      min = ks;
      at = r;
    }
  }
  abs(min, 0.9129, 2e-4);
  abs(at, (1.1997 * Math.tanh(1.1997)) / (2 * Math.PI), 2e-3);
  // Shallow water: Cg -> sqrt(gh), so Ks -> sqrt(C0 / (2 sqrt(gh))) and H ∝ h^(-1/4).
  const h = 0.002 * L0;
  rel(shoalingCoefficient({ T, h }), Math.sqrt(L0 / T / (2 * Math.sqrt(g * h))), 5e-3);
  rel(shoalingCoefficient({ T, h: 500 * L0 / 100 }), 1, 1e-9);
});

test('Snell refraction: sin θ = tanh(kh) sin θ0 and Kr = sqrt(cos θ0 / cos θ)', () => {
  const T = 9;
  const h = 6;
  const w = linearWave({ T, h });
  const r = refractionSnell({ T, h, theta0Deg: 40 });
  const th = Math.asin(Math.tanh(w.k * h) * Math.sin((40 * Math.PI) / 180));
  rel(r.theta, th, 1e-12);
  rel(r.Kr, Math.sqrt(Math.cos((40 * Math.PI) / 180) / Math.cos(th)), 1e-12);
  assert.equal(refractionSnell({ T, h, theta0Deg: 0 }).Kr, 1);
  assert.throws(() => refractionSnell({ T, h, theta0Deg: 95 }), RangeError);
});

test('profile transformation conserves shoreward energy flux E Cg cos θ', () => {
  const depths = [200, 50, 20, 10, 5, 3];
  const T = 10;
  const p = transformProfile({ H0: 2, T, theta0Deg: 30, depths });
  const C0 = (g * T) / (2 * Math.PI);
  const f0 = 2 * 2 * 0.5 * C0 * Math.cos(Math.PI / 6);
  depths.forEach((h, i) => {
    const w = linearWave({ T, h });
    rel(p.H[i] ** 2 * w.Cg * Math.cos((p.thetaDeg[i] * Math.PI) / 180), f0, 1e-10);
  });
  const capped = transformProfile({ H0: 2, T, depths: [10, 1], gammaB: 0.78 });
  assert.equal(capped.broken[0], 0);
  assert.equal(capped.broken[1], 1);
  rel(capped.H[1], 0.78, 1e-12);
});

test('breaking criteria: Miche, McCowan, Weggel and Goda limits and hand calculations', () => {
  const T = 8;
  const L0 = (g * T * T) / (2 * Math.PI);
  rel(breakingMiche({ T, h: 2000 }).Hmax, 0.142 * L0, 1e-9); // deep water: H/L = 0.142
  rel(breakingMiche({ T: 30, h: 0.5 }).Hmax / 0.5, 0.142 * 2 * Math.PI, 2e-3); // shallow: 0.89 h
  rel(breakingMcCowan({ h: 2 }), 1.56, 1e-12);
  // Weggel on a flat bed reduces to Hb/hb = 0.78 (a = 0, b = 1.56/2).
  const w0 = breakingWeggel({ h: 3, T, slope: 0 });
  rel(w0.index, 0.78, 1e-12);
  // Hand calculation, m = 0.05: a = 43.75(1 − e^−0.95) = 26.83, b = 1.56/(1 + e^−0.975) = 1.132
  const w = breakingWeggel({ h: 3, T, slope: 0.05 });
  abs(w.a, 26.83, 0.01);
  abs(w.b, 1.132, 0.001);
  rel(w.Hb, (w.b * 3) / (1 + (w.a * 3) / (g * 64)), 1e-12);
  rel(w.Hb / 3, w.b - (w.a * w.Hb) / (g * 64), 1e-12); // satisfies the implicit published form
  // Goda: shallow flat-bed limit Hb/hb -> 0.17 · 1.5π = 0.801; deep limit Hb -> 0.17 L0.
  rel(breakingGoda({ h: 0.001 * L0, T, slope: 0 }).index, 0.17 * 1.5 * Math.PI, 3e-3);
  rel(breakingGoda({ h: 50 * L0, T, slope: 0 }).Hb, 0.17 * L0, 1e-9);
  assert.ok(breakingGoda({ h: 3, T, slope: 0.1 }).Hb > breakingGoda({ h: 3, T, slope: 0.01 }).Hb);
});

test('Goda surf-zone heights: offshore branch, breaker-limited branch by hand, monotonic decay', () => {
  const T = 10;
  const L0 = (g * T * T) / (2 * Math.PI);
  const H0 = 4;
  const off = godaSurfZoneHeights({ H0, T, h: 0.3 * L0, slope: 0.02 });
  rel(off.Hs, off.Ks * H0, 1e-12);
  rel(off.Hmax, 1.8 * off.Ks * H0, 1e-12);
  assert.equal(off.surfZone, false);
  // Hand calculation at h = 4 m, slope 1/50, H0'/L0 = 4/156.13 = 0.02562
  const s = H0 / L0;
  const b0 = 0.028 * s ** -0.38 * Math.exp(20 * 0.02 ** 1.5);
  const b1 = 0.52 * Math.exp(4.2 * 0.02);
  const bm = Math.max(0.92, 0.32 * s ** -0.29 * Math.exp(2.4 * 0.02));
  const r = godaSurfZoneHeights({ H0, T, h: 4, slope: 0.02 });
  rel(r.Hs, Math.min(b0 * H0 + b1 * 4, bm * H0), 1e-12);
  abs(r.Hs, 2.74, 0.02);
  assert.equal(r.surfZone, true);
  assert.ok(r.Hmax > r.Hs && r.Hmax < 1.8 * r.Hs);
  let prev = Infinity;
  for (let h = 6; h >= 0.5; h -= 0.5) {
    const v = godaSurfZoneHeights({ H0, T, h, slope: 0.02 }).Hs;
    assert.ok(v <= prev + 1e-12);
    prev = v;
  }
});

test('Battjes–Janssen Qb satisfies (1 − Qb)/(−ln Qb) = (Hrms/Hm)² and its limits', () => {
  for (const b of [0.06, 0.2, 0.35, 0.5, 0.7, 0.9, 0.99, 0.9999]) {
    const q = battjesJanssenQb(b);
    assert.ok(q > 0 && q < 1);
    rel((1 - q) / -Math.log(q), b * b, 1e-9, `b=${b}`);
  }
  assert.equal(battjesJanssenQb(0), 0);
  assert.equal(battjesJanssenQb(1), 1);
  assert.equal(battjesJanssenQb(1.3), 1);
  rel(battjesStiveGamma(0.02), 0.5 + 0.4 * Math.tanh(0.66), 1e-12);
});

test('energy-flux balance without dissipation reproduces linear shoaling/refraction and analytic set-down', () => {
  const n = 301;
  const x = Float64Array.from({ length: n }, (_, i) => i * 5);
  const depth = x.map((v) => 20 - 0.012 * v); // 20 m to 2 m
  const Tp = 9;
  const r = energyFluxBalance({ x, depth, Hrms0: 0.4, Tp, theta0Deg: 25, alpha: 0, fw: 0 });
  assert.equal(r.nWet, n);
  const w0 = linearWave({ T: Tp, h: r.h[0] });
  const sn = Math.sin((25 * Math.PI) / 180) / w0.C;
  for (const i of [50, 150, 300]) {
    const w = linearWave({ T: Tp, h: r.h[i] });
    const th = Math.asin(sn * w.C);
    const H = 0.4 * Math.sqrt((w0.Cg * Math.cos((25 * Math.PI) / 180)) / (w.Cg * Math.cos(th)));
    rel(r.Hrms[i], H, 1e-9, 'Hrms');
    rel(r.thetaDeg[i], (th * 180) / Math.PI, 1e-9, 'theta');
  }
  // Normal incidence: set-down η = −H² k / (8 sinh 2kh) (Longuet-Higgins & Stewart 1964).
  const r0 = energyFluxBalance({ x, depth, Hrms0: 0.4, Tp, alpha: 0 });
  for (const i of [100, 200, 300]) {
    const w = linearWave({ T: Tp, h: r0.h[i] });
    rel(r0.eta[i], (-(r0.Hrms[i] ** 2) * w.k) / (8 * Math.sinh(2 * w.k * r0.h[i])), 2e-3, 'set-down');
  }
  assert.ok(r0.eta[300] < r0.eta[0] && r0.eta[300] < 0);
});

test('energy-flux balance with Battjes–Janssen breaking: flux budget closes, Hrms <= Hm, set-up at the shore', () => {
  const n = 801;
  const x = Float64Array.from({ length: n }, (_, i) => i * 1);
  const depth = x.map((v) => 8 - 0.0101 * v); // crosses the still-water line near x = 792
  const Tp = 8;
  const rho = 1025;
  const r = energyFluxBalance({ x, depth, Hrms0: 1.0, Tp, gamma: 0.8, alpha: 1, fw: 0.01, rho });
  assert.ok(r.nWet > 780, 'march reaches the swash region');
  let diss = 0;
  for (let i = 1; i < r.nWet; i++) diss += 0.5 * (r.Db[i] + r.Df[i] + r.Db[i - 1] + r.Df[i - 1]) * (x[i] - x[i - 1]);
  rel(r.energyFlux[0] - r.energyFlux[r.nWet - 1], diss, 2e-3, 'flux budget');
  for (let i = 0; i < r.nWet; i++) {
    const w = linearWave({ T: Tp, h: r.h[i] });
    const Hm = (0.88 / w.k) * Math.tanh((0.8 * w.k * r.h[i]) / 0.88);
    assert.ok(r.Hrms[i] <= Hm * (1 + 1e-9));
    assert.ok(r.Qb[i] >= 0 && r.Qb[i] <= 1);
  }
  // Set-down seaward of breaking, set-up shoreward; in the saturated inner surf zone Hrms/h stays below γ.
  assert.ok(Math.min(...r.eta.subarray(0, r.nWet)) < 0);
  assert.ok(r.setupMax > 0.05 && r.setupMax < 0.3, `setup ${r.setupMax}`);
  const k = r.nWet - 20;
  assert.ok(r.Hrms[k] / r.h[k] < 0.8 && r.Hrms[k] / r.h[k] > 0.3);
  // Momentum balance: dSxx/dx + ρ g h dη/dx = 0 integrated over the profile.
  let lhs = 0;
  for (let i = 1; i < r.nWet; i++) lhs += rho * g * 0.5 * (r.h[i] + r.h[i - 1]) * (r.eta[i] - r.eta[i - 1]);
  rel(lhs, r.Sxx[0] - r.Sxx[r.nWet - 1], 1e-6, 'momentum');
});

test('bottom-friction dissipation equals the Rayleigh average of ½ ρ fw <|u|³> (first principles)', () => {
  const Tp = 10;
  const h = 5;
  const H = 0.5;
  const fw = 0.02;
  const rho = 1025;
  const r = energyFluxBalance({ x: [0, 1], depth: [h, h], Hrms0: H, Tp, alpha: 0, fw, rho, setup: false });
  const w = linearWave({ T: Tp, h });
  // <|cos|³> = 4/(3π); E[H³] for a Rayleigh distribution with rms Hrms is (3√π/4) Hrms³, integrated numerically.
  let m3 = 0;
  const dH = 1e-4;
  for (let v = dH / 2; v < 6 * H; v += dH) m3 += v ** 3 * ((2 * v) / (H * H)) * Math.exp(-((v / H) ** 2)) * dH;
  const ub3 = (Math.PI / (Tp * Math.sinh(w.k * h))) ** 3 * m3;
  rel(r.Df[0], 0.5 * rho * fw * (4 / (3 * Math.PI)) * ub3, 1e-6);
});

test('Fresnel integrals against Abramowitz & Stegun Table 7.7 and the defining derivative', () => {
  const tab = [[0.5, 0.4923442, 0.0647324], [1, 0.7798934, 0.4382591], [2, 0.4882534, 0.3434157], [3, 0.6057208, 0.4963130], [5, 0.5636312, 0.4991914]];
  for (const [x, C, S] of tab) {
    const f = fresnel(x);
    abs(f.C, C, 2e-7, `C(${x})`);
    abs(f.S, S, 2e-7, `S(${x})`);
    const m = fresnel(-x);
    assert.equal(m.C, -f.C);
    assert.equal(m.S, -f.S);
  }
  // dC/dx = cos(πx²/2) across the series/asymptotic switch at x = 3
  for (const x of [2.5, 3.0, 3.5, 7.3]) {
    const d = 1e-5;
    abs((fresnel(x + d).C - fresnel(x - d).C) / (2 * d), Math.cos((Math.PI * x * x) / 2), 1e-6);
    abs((fresnel(x + d).S - fresnel(x - d).S) / (2 * d), Math.sin((Math.PI * x * x) / 2), 1e-6);
  }
  abs(fresnel(60).C, 0.5, 6e-3);
});

test('semi-infinite breakwater: Kd = 1/2 on the shadow line, -> 0 in the shadow, -> 1 in the open', () => {
  const L = 10;
  for (const t0 of [45, 90, 135]) {
    const a = (t0 * Math.PI) / 180;
    const r = 80 * L;
    abs(diffractionSemiInfinite({ x: r * Math.cos(a), y: r * Math.sin(a), L, theta0Deg: t0 }).Kd, 0.5, 0.01);
    const shadow = a / 4;
    assert.ok(diffractionSemiInfinite({ x: r * Math.cos(shadow), y: r * Math.sin(shadow), L, theta0Deg: t0 }).Kd < 0.06);
    const open = a + (Math.PI - a) / 2;
    abs(diffractionSemiInfinite({ x: r * Math.cos(open), y: r * Math.sin(open), L, theta0Deg: t0 }).Kd, 1, 0.06);
  }
  // On the exposed face (θ -> 2π) with normal incidence there is full reflection: standing wave, |F| up to 2.
  let mx = 0;
  for (let y = -30; y < -20; y += 0.05) mx = Math.max(mx, diffractionSemiInfinite({ x: 400, y, L }).Kd);
  abs(mx, 2, 0.05);
});

test('Sommerfeld solution satisfies the Helmholtz equation and the no-flux condition on the breakwater', () => {
  const L = 10;
  const k = (2 * Math.PI) / L;
  const d = 0.01;
  const F = (x, y) => diffractionSemiInfinite({ x, y, L, theta0Deg: 60 });
  for (const [x, y] of [[12, 7], [-9, 14], [30, 3], [5, -20], [-15, -4]]) {
    const c = F(x, y);
    for (const part of ['re', 'im']) {
      const lap = (F(x + d, y)[part] + F(x - d, y)[part] + F(x, y + d)[part] + F(x, y - d)[part] - 4 * c[part]) / (d * d);
      abs(lap + k * k * c[part], 0, 2e-4, `Helmholtz ${part} at ${x},${y}`);
    }
  }
  // ∂F/∂y = 0 on both faces of the wall (x > 0, y = ±0)
  for (const x of [3, 11, 40]) {
    for (const side of [1, -1]) {
      for (const part of ['re', 'im']) {
        const a = F(x, side * 1e-9)[part];
        const b = F(x, side * 1e-3)[part];
        const c = F(x, side * 2e-3)[part];
        abs((-3 * a + 4 * b - c) / 2e-3, 0, 2e-4, `wall ${x} ${side}`);
      }
    }
  }
  // Far from the tip on the open side the field is the undisturbed incident wave e^(−ik r cos(θ−θ0)).
  const p = F(-300, 500);
  const ph = -k * (-300 * Math.cos(Math.PI / 3) + 500 * Math.sin(Math.PI / 3));
  abs(p.re, Math.cos(ph), 0.03);
  abs(p.im, Math.sin(ph), 0.03);
});

test('breakwater gap: symmetric for normal incidence, tends to the incident wave in a wide gap, decays behind a narrow one', () => {
  const L = 10;
  const a = diffractionGap({ x: 17, y: 60, L, B: 40 });
  const b = diffractionGap({ x: -17, y: 60, L, B: 40 });
  rel(a.Kd, b.Kd, 1e-12);
  abs(diffractionGap({ x: 0, y: 30, L, B: 2000 }).Kd, 1, 0.02);
  const near = diffractionGap({ x: 0, y: 2 * L, L, B: 2 * L }).Kd;
  const far = diffractionGap({ x: 0, y: 20 * L, L, B: 2 * L }).Kd;
  assert.ok(far < near && far > 0.1 && far < 0.5, `far ${far}, near ${near}`);
  // Deep in the lee of an arm the gap solution reduces to the single-tip solution of the nearer tip.
  const g1 = diffractionGap({ x: 300, y: 30, L, B: 100 }).Kd;
  const s1 = diffractionSemiInfinite({ x: 250, y: 30, L }).Kd;
  abs(g1, s1, 0.02);
});

test('Seelig & Ahrens reflection: Kr = a ξ²/(b + ξ²)', () => {
  rel(reflectionSeeligAhrens({ xi: 3 }), 9 / 14.5, 1e-12);
  rel(reflectionSeeligAhrens({ xi: 3, a: 0.6, b: 6.6 }), (0.6 * 9) / 15.6, 1e-12);
  assert.ok(reflectionSeeligAhrens({ xi: 50 }) < 1);
  assert.throws(() => reflectionSeeligAhrens({ xi: -1 }), RangeError);
});

test('wave–current interaction: deep-water closed forms (Longuet-Higgins & Stewart 1961) and blocking at U = −C0/4', () => {
  const T = 8;
  const C0 = (g * T) / (2 * Math.PI);
  const h = 3000;
  for (const U of [-2.5, -1, 0.5, 2]) {
    const r = waveCurrent({ T, h, U, H0: 1.5 });
    const C = 0.5 * C0 * (1 + Math.sqrt(1 + (4 * U) / C0)); // relative phase speed
    rel(r.Cr, C, 1e-9, 'C');
    rel(r.heightRatio, C0 / Math.sqrt(C * (C + 2 * U)), 1e-9, 'H/H0');
    rel(r.H, 1.5 * r.heightRatio, 1e-12);
    rel((2 * Math.PI / T - r.k * U) ** 2, g * r.k * Math.tanh(r.k * h), 1e-10, 'Doppler dispersion');
  }
  const still = waveCurrent({ T, h: 12, U: 0 });
  rel(still.heightRatio, 1, 1e-12);
  rel(still.k, linearWave({ T, h: 12 }).k, 1e-12);
  rel(waveCurrent({ T, h, U: 0 }).Ublock, -C0 / 4, 1e-9);
  assert.equal(waveCurrent({ T, h, U: -C0 / 4 - 0.01 }).blocked, true);
  assert.equal(waveCurrent({ T, h, U: -C0 / 4 + 0.01 }).blocked, false);
  // Opposing currents steepen and shorten the waves, following currents do the opposite (finite depth).
  const opp = waveCurrent({ T, h: 10, U: -1 });
  const fol = waveCurrent({ T, h: 10, U: 1 });
  assert.ok(opp.heightRatio > 1 && opp.L < still.L && fol.heightRatio < 1 && fol.L > still.L);
  rel(opp.Cga, opp.Cgr - 1, 1e-12);
});
