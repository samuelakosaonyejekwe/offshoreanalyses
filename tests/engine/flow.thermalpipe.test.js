import test from 'node:test';
import assert from 'node:assert/strict';
import {
  nusseltInternal, internalFilmCoefficient, externalFilmCoefficient, burialFilmCoefficient, equivalentSoilThickness,
  overallU, segmentTemperature, temperatureProfile, insulationThickness, cooldownLumped, radialTransient,
  cooldownMultiCapacity, timeToTemperature, pipeInPipeAnnulus, heatingPower, buriedThermalResistance,
} from '../../src/engine/flow/thermalpipe.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('internal Nusselt correlations — hand calculations at Re = 1e5, Pr = 0.7', () => {
  // Dittus–Boelter cooling: 0.023·1e5^0.8·0.7^0.3 = 0.023·10000·0.8985 = 206.7
  assert.ok(rel(nusseltInternal({ Re: 1e5, Pr: 0.7, method: 'dittusBoelter' }), 0.023 * 1e4 * 0.7 ** 0.3) < 1e-12);
  assert.ok(rel(nusseltInternal({ Re: 1e5, Pr: 0.7, method: 'dittusBoelter', heating: true }), 199.4) < 1e-3);
  // Gnielinski: f = (0.79 ln 1e5 − 1.64)^−2 = 0.017993; Nu = (f/8)(Re − 1000)Pr/[1 + 12.7√(f/8)(Pr^(2/3) − 1)] = 178.6
  const f = (0.79 * Math.log(1e5) - 1.64) ** -2;
  const nu = (f / 8) * 99000 * 0.7 / (1 + 12.7 * Math.sqrt(f / 8) * (0.7 ** (2 / 3) - 1));
  assert.ok(rel(nusseltInternal({ Re: 1e5, Pr: 0.7, method: 'gnielinski' }), nu) < 1e-12);
  assert.ok(rel(nu, 178.6) < 1e-3);
  assert.equal(nusseltInternal({ Re: 1500, Pr: 50 }), 3.66);
  // Sieder–Tate: 1.86 (Re Pr D/L)^(1/3) (μ/μw)^0.14 = 1.86·(1000·100/50)^(1/3)·1 = 23.43
  assert.ok(rel(nusseltInternal({ Re: 1000, Pr: 100, method: 'siederTate', LoverD: 50 }), 1.86 * 2000 ** (1 / 3)) < 1e-12);
  // 'auto' is continuous across the transition band
  assert.ok(rel(nusseltInternal({ Re: 2999.999, Pr: 5 }), nusseltInternal({ Re: 3000, Pr: 5 })) < 1e-5);
  const h = internalFilmCoefficient({ v: 2, D: 0.2, rho: 800, mu: 2e-3, k: 0.13, cp: 2000 });
  assert.ok(rel(h.Re, 1.6e5) < 1e-12 && rel(h.Pr, 2000 * 2e-3 / 0.13) < 1e-12 && rel(h.h, h.Nu * 0.13 / 0.2) < 1e-12);
});

test('Churchill–Bernstein cross-flow — hand calculation (Re = 1e4, Pr = 7)', () => {
  const Re = 1e4;
  const Pr = 7;
  const nu = 0.3 + 0.62 * 100 * Pr ** (1 / 3) / (1 + (0.4 / Pr) ** (2 / 3)) ** 0.25 * (1 + (Re / 282000) ** 0.625) ** 0.8;
  // seawater-like: ρ = 1000, μ = 1e-3, D = 0.1 → v = 0.1 m/s; k = 0.6, cp = 4200 → Pr = 7
  const r = externalFilmCoefficient({ v: 0.1, D: 0.1, rho: 1000, mu: 1e-3, k: 0.6, cp: 4200 });
  assert.ok(rel(r.Nu, nu) < 1e-12);
  assert.ok(nu > 120 && nu < 135); // 0.3 + 118.6/1.037·1.097
  assert.ok(rel(r.h, nu * 6) < 1e-12);
});

test('overall U of a three-layer pipe — hand calculation', () => {
  // Di = 0.3 m; steel 20 mm (k = 45), insulation 50 mm (k = 0.2), concrete 40 mm (k = 1.5); hi = 500, ho = 1000
  const u = overallU({ Di: 0.3, layers: [{ t: 0.02, k: 45 }, { t: 0.05, k: 0.2 }, { t: 0.04, k: 1.5 }], hi: 500, ho: 1000 });
  const R = 1 / (500 * Math.PI * 0.3) + Math.log(0.34 / 0.3) / (2 * Math.PI * 45) + Math.log(0.44 / 0.34) / (2 * Math.PI * 0.2)
    + Math.log(0.52 / 0.44) / (2 * Math.PI * 1.5) + 1 / (1000 * Math.PI * 0.52);
  assert.ok(rel(u.resistance, R) < 1e-12);
  assert.ok(rel(u.U, 1 / (R * Math.PI * 0.3)) < 1e-12);
  assert.ok(rel(u.Do, 0.52) < 1e-12);
  // insulation dominates: U ≈ 2 k/(Di ln(0.44/0.34)) = 5.17 W/m²K upper bound
  assert.ok(u.U < 2 * 0.2 / (0.3 * Math.log(0.44 / 0.34)) && u.U > 4.5);
  // referring U to the outer diameter keeps UA unchanged
  const uo = overallU({ Di: 0.3, layers: [{ t: 0.02, k: 45 }, { t: 0.05, k: 0.2 }, { t: 0.04, k: 1.5 }], hi: 500, ho: 1000, Dref: 0.52 });
  assert.ok(rel(uo.U * 0.52, u.U * 0.3) < 1e-12);
  assert.throws(() => overallU({ Di: 0.3 }), RangeError);
});

test('burial: cosh⁻¹ shape factor, deep-burial logarithmic limit, film correction and partial burial', () => {
  const D = 0.5;
  const b = burialFilmCoefficient({ D, H: 1.5, kSoil: 1.2 });
  assert.ok(rel(b.h, 2 * 1.2 / (D * Math.acosh(6))) < 1e-12);
  assert.ok(rel(b.resistance, Math.acosh(6) / (2 * Math.PI * 1.2)) < 1e-12);
  // deep burial: cosh⁻¹(2H/D) → ln(4H/D)
  const deep = burialFilmCoefficient({ D, H: 50, kSoil: 1.2 });
  assert.ok(rel(deep.resistance, Math.log(400) / (2 * Math.PI * 1.2)) < 1e-4);
  // finite seabed film coefficient = extra cover k/h
  const f = burialFilmCoefficient({ D, H: 1.5, kSoil: 1.2, hSea: 60 });
  assert.ok(rel(f.h, 2 * 1.2 / (D * Math.acosh(2 * (1.5 + 0.02) / D))) < 1e-12);
  // partial burial interpolates monotonically between exposed and just-buried
  const hE = 400;
  const half = burialFilmCoefficient({ D, H: 0, kSoil: 1.2, hSea: 60, hExposed: hE });
  const just = burialFilmCoefficient({ D, H: D / 2, kSoil: 1.2, hSea: 60, hExposed: hE });
  const exp = burialFilmCoefficient({ D, H: -D / 2, kSoil: 1.2, hSea: 60, hExposed: hE });
  assert.equal(exp.h, hE);
  assert.equal(half.exposedFraction, 0.5);
  assert.ok(rel(half.h, 0.5 * (hE + just.h)) < 1e-12);
  assert.ok(just.h < half.h && half.h < exp.h);
  // equivalent annulus reproduces the steady resistance; IEC 60287 T4 is the same expression
  const t = equivalentSoilThickness({ D, H: 1.5 });
  assert.ok(rel(Math.log((D + 2 * t) / D) / (2 * Math.PI * 1.2), b.resistance) < 1e-12);
  const t4 = buriedThermalResistance({ depth: 1.5, De: D, rhoT: 1 / 1.2, dT: 30 });
  assert.ok(rel(t4.T4, b.resistance) < 1e-12);
  assert.ok(rel(t4.q, 30 / b.resistance) < 1e-12);
});

test('steady temperature: exponential profile by hand, Joule–Thomson asymptote, comparison with RK4 integration', () => {
  // ṁ cp/(U π D) = 50·2500/(3·π·0.3) = 44 210 m
  const lam = 50 * 2500 / (3 * Math.PI * 0.3);
  const s = segmentTemperature({ Tin: 350, Tamb: 278, L: 20000, mDot: 50, cp: 2500, U: 3, D: 0.3 });
  assert.ok(rel(s.lambda, lam) < 1e-12);
  assert.ok(rel(s.Tout, 278 + 72 * Math.exp(-20000 / lam)) < 1e-12);
  // energy balance: heat lost = ṁ cp (Tin − Tout)
  assert.ok(rel(s.heatLoss, 50 * 2500 * (350 - s.Tout)) < 1e-10);
  // with JT cooling and an uphill slope the far-field temperature is T_amb + λ(μJT dP/dx − g sinθ/cp)
  const args = { Tin: 320, Tamb: 278, mDot: 30, cp: 2400, U: 2, D: 0.4, muJT: 4e-6, dPdx: -40, dzdx: 0.01 };
  const lam2 = 30 * 2400 / (2 * Math.PI * 0.4);
  const far = segmentTemperature({ ...args, L: 1e7 });
  assert.ok(rel(far.Tout, 278 + lam2 * (4e-6 * -40 - 9.80665 * 0.01 / 2400)) < 1e-10);
  assert.ok(far.Tout < 278);
  // RK4 integration of dT/dx = −(T − Tamb)/λ + μJT dP/dx − g/cp dz/dx over 15 km
  const rhs = (T) => -(T - 278) / lam2 + 4e-6 * -40 - 9.80665 * 0.01 / 2400;
  let T = 320;
  const n = 1500;
  const h = 15000 / n;
  for (let i = 0; i < n; i++) {
    const k1 = rhs(T);
    const k2 = rhs(T + 0.5 * h * k1);
    const k3 = rhs(T + 0.5 * h * k2);
    const k4 = rhs(T + h * k3);
    T += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  }
  assert.ok(Math.abs(segmentTemperature({ ...args, L: 15000 }).Tout - T) < 1e-9);
  // marching over several profile intervals equals one analytic segment for uniform conditions
  const p = temperatureProfile({ Tin: 350, mDot: 50, cp: 2500, U: 3, D: 0.3, Tamb: 278,
    profile: [{ x: 0 }, { x: 3000 }, { x: 11000 }, { x: 20000 }] });
  assert.ok(Math.abs(p.Tarrival - s.Tout) < 1e-10);
  assert.ok(rel(p.heatLoss, s.heatLoss) < 1e-10);
  // adiabatic segment
  assert.equal(segmentTemperature({ Tin: 350, Tamb: 278, L: 1000, mDot: 50, cp: 2500, U: 0, D: 0.3 }).Tout, 350);
});

test('insulation sizing meets a target U-value and a target arrival temperature', () => {
  const layers = [{ t: 0.015, k: 45 }, { t: 0, k: 0.17 }, { t: 0.04, k: 1.8 }];
  const a = insulationThickness({ Di: 0.25, layers, index: 1, hi: 800, ho: 500, Utarget: 2.5 });
  assert.ok(Math.abs(a.U - 2.5) < 1e-8);
  const check = overallU({ Di: 0.25, layers: layers.map((l, i) => (i === 1 ? { ...l, t: a.t } : l)), hi: 800, ho: 500 });
  assert.ok(Math.abs(check.U - 2.5) < 1e-8);
  const b = insulationThickness({ Di: 0.25, layers, index: 1, hi: 800, ho: 500, Tarrival: 318, Tin: 350, Tamb: 278,
    L: 15000, mDot: 25, cp: 2300 });
  const arr = segmentTemperature({ Tin: 350, Tamb: 278, L: 15000, mDot: 25, cp: 2300, U: b.U, D: 0.25 }).Tout;
  assert.ok(Math.abs(arr - 318) < 1e-6);
  assert.ok(b.t > 0.005 && b.t < 0.3);
  // a bare pipe that already satisfies the target needs no insulation; an impossible target throws
  assert.equal(insulationThickness({ Di: 0.25, layers, index: 1, hi: 800, ho: 500, Utarget: 1e4 }).t, 0);
  assert.throws(() => insulationThickness({ Di: 0.25, layers, index: 1, hi: 800, ho: 500, Utarget: 0.01 }), RangeError);
});

test('lumped cooldown: analytic exponential and time to the hydrate temperature', () => {
  // U = 3 W/m²K on D = 0.3; heat capacity 150 kJ/(m·K) → τ = 150e3/(3π·0.3) = 53 052 s
  const c = cooldownLumped({ T0: 340, Tamb: 277, U: 3, D: 0.3, heatCapacity: 150e3, Tthreshold: 293, t: [0, 53051.65] });
  assert.ok(rel(c.tau, 150e3 / (3 * Math.PI * 0.3)) < 1e-12);
  assert.ok(rel(c.timeTo, c.tau * Math.log(63 / 16)) < 1e-12);
  assert.ok(Math.abs(c.T[1] - (277 + 63 / Math.E)) < 1e-3);
  assert.equal(cooldownLumped({ T0: 340, Tamb: 277, U: 3, D: 0.3, heatCapacity: 1e5, Tthreshold: 270 }).timeTo, Infinity);
});

const wall = [{ t: 0.02, k: 45, rho: 7850, cp: 480 }, { t: 0.06, k: 0.2, rho: 900, cp: 1800 }];

test('radial transient: exact energy conservation and steady state equal to the multilayer U-value', () => {
  const o = { Di: 0.3, contents: { rho: 800, cp: 2000 }, layers: wall, hi: 300, ho: 800, Tamb: 277, T0: 340 };
  const r = radialTransient({ ...o, tEnd: 24 * 3600, dt: 120 });
  assert.ok(rel(r.energyStored0 - r.energyStored, r.energyLost) < 1e-10);
  const u = overallU({ Di: 0.3, layers: wall, hi: 300, ho: 800 });
  assert.ok(rel(r.Usteady, u.U) < 1e-12);
  // initial steady profile: heat flow is uniform and equals UA·ΔT
  assert.ok(rel(r.qLoss[0], u.UA * 63) < 1e-9);
  assert.ok(rel(r.qFluid[0], u.UA * 63) < 1e-9);
  // contents cool monotonically towards ambient; inner wall is colder than the contents
  for (let i = 1; i < r.t.length; i++) assert.ok(r.Tcontents[i] < r.Tcontents[i - 1] && r.Tcontents[i] > 277);
  assert.ok(r.TwallInner[10] < r.Tcontents[10]);
  // warm-up from ambient with the fluid held at 340 K converges to the same steady loss; the heat supplied equals
  // the heat stored plus the heat lost
  const w = radialTransient({ ...o, mode: 'warmup', tEnd: 5 * 24 * 3600, dt: 300 });
  assert.ok(rel(w.qLoss[w.qLoss.length - 1], u.UA * 63) < 1e-4);
  assert.ok(w.qLoss[0] === 0 && w.qFluid[0] > w.qFluid[w.qFluid.length - 1]);
  assert.ok(rel(w.energyIn, w.energyStored - w.energyStored0 + w.energyLost) < 1e-10);
});

test('radial transient reduces to the lumped exponential when the wall is thin and highly conductive', () => {
  // one thin steel layer, well-mixed contents, all resistance in the outer film: τ = ΣC/(ho π Do)
  const Di = 0.2;
  const layer = { t: 0.005, k: 5000, rho: 7850, cp: 480 };
  const ho = 15;
  const Cc = 1000 * 4000 * Math.PI * Di * Di / 4;
  const Cw = 7850 * 480 * Math.PI * (0.105 ** 2 - 0.1 ** 2);
  const tau = (Cc + Cw) / (ho * Math.PI * 0.21);
  const r = radialTransient({ Di, contents: { rho: 1000, cp: 4000 }, layers: [layer], ho, Tamb: 280, T0: 350,
    initial: 'uniform', tEnd: tau, dt: tau / 4000, theta: 0.5 });
  assert.ok(Math.abs(r.Tcontents[r.Tcontents.length - 1] - (280 + 70 / Math.E)) < 0.01);
  const lump = cooldownLumped({ T0: 350, Tamb: 280, U: ho, D: 0.21, heatCapacity: Cc + Cw, Tthreshold: 320 });
  assert.ok(rel(timeToTemperature({ t: r.t, T: r.Tcontents, threshold: 320 }), lump.timeTo) < 2e-3);
});

test('radial transient: first-order (θ = 1) and second-order (θ = ½) convergence in time step; grid convergence', () => {
  const o = { Di: 0.3, contents: { rho: 800, cp: 2000 }, layers: wall, hi: 300, ho: 800, Tamb: 277, T0: 340, tEnd: 14400 };
  const end = (dt, theta, cells = 8) => {
    const r = radialTransient({ ...o, dt, theta, cellsPerLayer: cells });
    return r.Tcontents[r.Tcontents.length - 1];
  };
  const ref = end(7.5, 0.5);
  const e1 = Math.abs(end(480, 1) - ref);
  const e2 = Math.abs(end(240, 1) - ref);
  assert.ok(e1 / e2 > 1.8 && e1 / e2 < 2.2, `backward Euler ratio ${e1 / e2}`);
  const c1 = Math.abs(end(960, 0.5) - ref);
  const c2 = Math.abs(end(480, 0.5) - ref);
  assert.ok(c1 / c2 > 3.5 && c1 / c2 < 4.5, `Crank–Nicolson ratio ${c1 / c2}`);
  // spatial refinement: second-order convergence of the finite-volume discretisation
  const g = (cells) => end(30, 0.5, cells);
  const fine = g(64);
  const r48 = Math.abs(g(4) - fine) / Math.abs(g(8) - fine);
  assert.ok(r48 > 3 && r48 < 5.5, `grid ratio ${r48}`);
  // the lumped multi-capacity model (one node per layer) is within 1 K of the refined solution here
  const mc = cooldownMultiCapacity({ ...o, dt: 30 });
  assert.ok(Math.abs(mc.Tcontents[mc.Tcontents.length - 1] - fine) < 1.0);
});

test('time-to-temperature interpolation', () => {
  assert.equal(timeToTemperature({ t: [0, 10, 20], T: [50, 40, 30], threshold: 35 }), 15);
  assert.equal(timeToTemperature({ t: [0, 10, 20], T: [50, 40, 30], threshold: 60 }), 0);
  assert.equal(timeToTemperature({ t: [0, 10, 20], T: [50, 40, 30], threshold: 10 }), Infinity);
  assert.equal(timeToTemperature({ t: [0, 10, 20], T: [10, 20, 40], threshold: 30 }), 15); // warm-up
});

test('pipe-in-pipe annulus: black-body radiation, conduction limit and convective enhancement', () => {
  const gas = { k: 0.026, rho: 1.2, mu: 1.8e-5, cp: 1005 };
  // black surfaces: q = σ π D1 (T1⁴ − T2⁴)
  const b = pipeInPipeAnnulus({ D1: 0.2, D2: 0.3, T1: 350, T2: 280, eps1: 1, eps2: 1, gas });
  assert.ok(rel(b.qRad, 5.670374419e-8 * Math.PI * 0.2 * (350 ** 4 - 280 ** 4)) < 1e-12);
  // grey surfaces: divide by 1/ε1 + (D1/D2)(1/ε2 − 1)
  const g = pipeInPipeAnnulus({ D1: 0.2, D2: 0.3, T1: 350, T2: 280, eps1: 0.3, eps2: 0.6, gas });
  assert.ok(rel(g.qRad, b.qRad / (1 / 0.3 + (0.2 / 0.3) * (1 / 0.6 - 1))) < 1e-12);
  assert.ok(g.kEff > gas.k); // natural convection active in a 50 mm gap (Ra_c ≫ 100)
  assert.ok(rel(g.q, g.qRad + g.qCond) < 1e-12);
  // narrow gap: pure conduction 2πkΔT/ln(D2/D1)
  const n = pipeInPipeAnnulus({ D1: 0.2, D2: 0.204, T1: 350, T2: 280, gas });
  assert.equal(n.kEff, gas.k);
  assert.ok(rel(n.qCond, 2 * Math.PI * 0.026 * 70 / Math.log(1.02)) < 1e-12);
  assert.ok(rel(n.kEquivalent * 2 * Math.PI * 70 / Math.log(1.02), n.q) < 1e-12);
});

test('heating power demand — hand calculation', () => {
  // hold 25 °C above ambient with U = 4 W/m²K on D = 0.3 m: q′ = 4π·0.3·25 = 94.2 W/m; efficiency 0.7; 10 km
  const p = heatingPower({ U: 4, D: 0.3, L: 10000, Ttarget: 302, Tamb: 277, efficiency: 0.7, heatCapacity: 120e3, warmupTime: 48 * 3600 });
  assert.ok(rel(p.holdPerLength, 4 * Math.PI * 0.3 * 25 / 0.7) < 1e-12);
  assert.ok(rel(p.holdPower, p.holdPerLength * 1e4) < 1e-12);
  // warm-up in 48 h: sensible 120e3·25/172800 = 17.4 W/m plus mean loss at half the temperature rise
  assert.ok(rel(p.warmupPerLength, (120e3 * 25 / 172800 + 4 * Math.PI * 0.3 * 12.5) / 0.7) < 1e-12);
});
