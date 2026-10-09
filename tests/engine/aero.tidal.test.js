import test from 'node:test';
import assert from 'node:assert/strict';
import {
  TIDAL_SPEEDS_DEG_PER_HOUR, constituentPeriod, synthesizeTide, springNeap, velocityDistribution,
  currentProfile, rotorAveragedVelocity, actuatorRotor, rotorFromSolver, tidalPowerCurve, cavitationCheck,
  channelBlockage, channelBlockageLimit, blockageCorrection, tidalArrayFlow, bidirectionalYield, tidalEnergy,
  extremeThrust, supportLoads, kineticPowerDensity, fluxResource, tidalRangeEbb,
} from '../../src/engine/aero/tidal.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);

test('constituent periods: M2 = 12.42 h, S2 = 12 h, K1 = 23.93 h, M4 = M2/2', () => {
  close(constituentPeriod('M2').periodHours, 12.4206012, 1e-7);
  close(constituentPeriod('S2').periodHours, 12, 1e-12);
  close(constituentPeriod('K1').periodHours, 23.9344697, 1e-7);
  close(constituentPeriod('O1').periodHours, 25.8193417, 1e-7);
  close(constituentPeriod('M4').periodHours, constituentPeriod('M2').periodHours / 2, 1e-12);
  // astronomical identities among the speeds: K1 + O1 = M2, K1 + P1 = S2, K2 = 2 K1
  const s = TIDAL_SPEEDS_DEG_PER_HOUR;
  close(s.K1 + s.O1, s.M2, 1e-8);
  close(s.K1 + s.P1, s.S2, 1e-8);
  close(2 * s.K1, s.K2, 1e-8);
  assert.throws(() => constituentPeriod('X9'), RangeError);
});

test('harmonic synthesis: M2 zero-crossing period and spring–neap cycle of 14.77 days', () => {
  const m2 = synthesizeTide({ constituents: [{ name: 'M2', amplitude: 2, phaseDeg: 90 }], duration: 5 * 86400, dt: 60 });
  // upward zero crossings by linear interpolation
  const zc = [];
  for (let i = 1; i < m2.u.length; i++) {
    if (m2.u[i - 1] < 0 && m2.u[i] >= 0) zc.push(m2.t[i - 1] + (60 * -m2.u[i - 1]) / (m2.u[i] - m2.u[i - 1]));
  }
  close((zc[zc.length - 1] - zc[0]) / (zc.length - 1) / 3600, 12.4206012, 1e-6);
  const sn = springNeap({ springPeak: 3.2, neapPeak: 1.6 });
  close(sn.springNeapPeriod / 86400, 14.7652944, 1e-6);
  assert.deepEqual(sn.constituents.map((c) => c.amplitude), [2.4000000000000004, 0.8]);
  close(sn.envelope(0), 3.2, 1e-12);
  close(sn.envelope(sn.springNeapPeriod / 2), 1.6, 1e-12);
  const y = synthesizeTide({ constituents: sn.constituents, duration: 30 * 86400, dt: 300 });
  let mx = 0;
  let neapMax = 0;
  const tNeap = sn.springNeapPeriod / 2;
  for (let i = 0; i < y.u.length; i++) {
    mx = Math.max(mx, Math.abs(y.u[i]));
    if (Math.abs(y.t[i] - tNeap) < 6.3 * 3600) neapMax = Math.max(neapMax, Math.abs(y.u[i]));
  }
  close(mx, 3.2, 1e-3);
  close(neapMax, 1.6, 5e-3);
  const withMean = synthesizeTide({ constituents: [{ name: 'S2', amplitude: 1 }], t: [0, 3 * 3600, 6 * 3600], mean: 0.2 });
  close(withMean.u[0], 1.2, 1e-12);
  close(withMean.u[1], 0.2, 1e-12);
  close(withMean.u[2], -0.8, 1e-12);
});

test('velocity distribution of a sinusoidal tide: mean 2/π, cube-mean 4/(3π), exceedance of U0/2 is 2/3', () => {
  const n = 200000;
  const u = Float64Array.from({ length: n }, (_, i) => 2.5 * Math.cos((2 * Math.PI * (i + 0.5)) / n));
  const d = velocityDistribution(u, { binWidth: 0.25 });
  close(d.mean, (2.5 * 2) / Math.PI, 1e-6);
  close(d.cubeRootMeanCube ** 3, (2.5 ** 3 * 4) / (3 * Math.PI), 1e-6);
  close(d.bins.reduce((s, b) => s + b.prob, 0), 1, 1e-12);
  // |cos θ| > ½ for 2/3 of the cycle
  let idx = 0;
  while (d.speed[idx] > 1.25) idx++;
  close(d.exceedance[idx], 2 / 3, 1e-4);
  const k = kineticPowerDensity(u, 1025);
  close(k.mean, 0.5 * 1025 * 2.5 ** 3 * (4 / (3 * Math.PI)), 1e-6);
  close(k.max, 0.5 * 1025 * 2.5 ** 3, 1e-6);
});

test('vertical profiles: 1/7th-law depth mean is 7/8 of the surface value; rotor averages', () => {
  const h = 40;
  const prof = currentProfile({ uRef: 3, zRef: h, n: 7 });
  let s = 0;
  const N = 100000;
  for (let i = 0; i < N; i++) s += prof(((i + 0.5) / N) * h) / N;
  close(s, (7 / 8) * 3, 1e-6);
  close(currentProfile({ uRef: 3, zRef: h, n: 10 })(h / 2), 3 * 0.5 ** 0.1, 1e-12);
  const lg = currentProfile({ model: 'log', uRef: 2, zRef: 10, z0: 0.01 });
  close(lg(1), (2 * Math.log(100)) / Math.log(1000), 1e-12);
  const uni = rotorAveragedVelocity({ profile: () => 2.2, hubHeight: 15, D: 18 });
  close(uni.mean, 2.2, 1e-12);
  close(uni.powerVelocity, 2.2, 1e-12);
  // linear shear u = c z: disc mean equals the hub value; mean of u³ = c³ (h³ + 3 h R²/4)
  const lin = rotorAveragedVelocity({ profile: (z) => 0.1 * z, hubHeight: 20, D: 20, nSeg: 2000 });
  close(lin.mean, 2, 1e-6);
  close(lin.cubeMean, 1e-3 * (8000 + (3 * 20 * 100) / 4), 1e-5);
  assert.ok(lin.powerVelocity > lin.mean);
});

test('turbine power curve: actuator fallback, cut-in, rated cap and thrust shedding', () => {
  const D = 18;
  const A = (Math.PI * D * D) / 4;
  const c = tidalPowerCurve({ D, ratedPower: 1e6, cutIn: 0.8, uMax: 4, step: 0.1, efficiency: 0.8 });
  const at = (U) => c.u.findIndex((u) => Math.abs(u - U) < 1e-9);
  assert.equal(c.p[at(0.5)], 0);
  // below rated: Cp = 0.8 · 16/27, Ct = 8/9
  close(c.p[at(2)], 0.5 * 1025 * A * 0.8 * (16 / 27) * 8, 1e-12);
  close(c.ct[at(2)], 8 / 9, 1e-12);
  close(c.thrust[at(2)], 0.5 * 1025 * A * (8 / 9) * 4, 1e-12);
  // rated speed: (1e6 / (½ρA·0.474))^(1/3) = 2.53 m/s → first tabulated point at rated is 2.6
  close(Math.cbrt(1e6 / (0.5 * 1025 * A * 0.8 * (16 / 27))), 2.529, 5e-4);
  close(c.ratedSpeed, 2.6, 1e-9);
  assert.equal(c.p[at(3.5)], 1e6);
  assert.ok(c.ct[at(3.5)] < c.ct[at(2)]);
  // injected rotor model
  const inj = tidalPowerCurve({ D, ratedPower: 1e6, rotor: () => ({ Cp: 0.41, Ct: 0.8 }) });
  close(inj.p[at(2)], 0.5 * 1025 * A * 0.41 * 8, 1e-12);
  // adapter for a BEM-style solver called at fixed tip-speed ratio
  const calls = [];
  const rotor = rotorFromSolver({
    solve: (o) => {
      calls.push(o);
      return { power: 0.5 * o.rho * A * 0.4 * o.U ** 3, thrust: 0.5 * o.rho * A * 0.75 * o.U ** 2 };
    },
    rotor: { name: 'r' },
    D,
    tsr: 5,
  });
  const perf = rotor(2);
  close(perf.Cp, 0.4, 1e-12);
  close(perf.Ct, 0.75, 1e-12);
  close(calls[0].omega, (5 * 2) / 9, 1e-12);
  assert.deepEqual(actuatorRotor({ a: 0.25, efficiency: 1 })(), { Cp: 4 * 0.25 * 0.5625, Ct: 0.75 });
});

test('cavitation number at the blade tip by hand', () => {
  // U = 3 m/s, tip speed 12 m/s, D = 18 m, hub 20 m below the surface → tip at 11 m
  // p = 101325 + 1025·9.80665·11 − 1700 = 210194.98 Pa; ½ρW² = ½·1025·(9 + 144) = 78412.5 Pa; σ = 2.68063
  const r = cavitationCheck({ U: 3, tipSpeed: 12, D: 18, hubDepth: 20, cpMin: -2.0 });
  close(r.tipDepth, 11, 1e-12);
  close(r.sigma, 210194.98 / 78412.5, 1e-7);
  assert.equal(r.cavitates, false);
  assert.equal(cavitationCheck({ U: 3, tipSpeed: 12, D: 18, hubDepth: 20, cpMin: -3.0 }).cavitates, true);
  // limiting tip speed: W_max² = 210194.98/(½·1025·2) = 205.068; tip = √(205.068 − 9) = 14.0024
  close(r.maxTipSpeed, Math.sqrt(210194.98 / 1025 - 9), 1e-7);
  const atLimit = cavitationCheck({ U: 3, tipSpeed: r.maxTipSpeed, D: 18, hubDepth: 20, cpMin: -2.0 });
  close(atLimit.sigma, 2.0, 1e-9);
  // shallower submergence lowers the allowable tip speed
  assert.ok(cavitationCheck({ U: 3, tipSpeed: 12, D: 18, hubDepth: 10, cpMin: -2.0 }).maxTipSpeed < r.maxTipSpeed);
});

test('Garrett & Cummins (2007) blocked-channel limit: Cp,max = (16/27)(1 − B)^−2 at α4 = 1/3', () => {
  for (const B of [0, 0.05, 0.1, 0.2, 0.3, 0.5]) {
    const r = channelBlockageLimit({ B });
    close(r.CpMax, 16 / 27 / (1 - B) ** 2, 1e-8, `B = ${B}`);
    close(r.alpha4, 1 / 3, 1e-5);
    // small Froude number tends to the rigid-lid result
    close(channelBlockageLimit({ B, Fr: 0.01 }).CpMax, 16 / 27 / (1 - B) ** 2, 1e-3, `B = ${B}, Fr = 0.01`);
  }
  // unblocked actuator disc: Betz
  const betz = channelBlockage({ B: 0, alpha4: 1 / 3 });
  close(betz.alpha2, 2 / 3, 1e-12);
  close(betz.Ct, 8 / 9, 1e-12);
  close(betz.Cp, 16 / 27, 1e-12);
  // Garrett–Cummins closed forms at B = 0.2, α4 = 1/3: u_t/u = 2/(3(1 + B)) ... checked via their turbine-velocity
  // relation α2 = α4 (β4 + α4)/(β4 + 2α4 − 1)
  const r = channelBlockage({ B: 0.2, alpha4: 0.4 });
  close(r.alpha2, (0.4 * (r.beta4 + 0.4)) / (r.beta4 + 2 * 0.4 - 1), 1e-9);
  assert.ok(r.beta4 > 1 && r.alpha2 > 0.4 && r.alpha2 < 1);
});

test('open-channel blockage (finite Froude number): conservation checks and higher limit than rigid lid', () => {
  const B = 0.25;
  const Fr = 0.25;
  const r = channelBlockage({ B, alpha4: 0.4, Fr });
  // continuity at the pressure-equalisation section
  close(r.h4OverH, (B * r.alpha2) / r.alpha4 + (1 - B * r.alpha2) / r.beta4, 1e-9);
  // bypass Bernoulli
  close(r.h4OverH, 1 - 0.5 * Fr * Fr * (r.beta4 ** 2 - 1), 1e-12);
  // momentum between upstream and section 4
  const lhs = (1 - r.h4OverH ** 2) / (2 * Fr * Fr) - 0.5 * B * r.Ct;
  const rhs = B * r.alpha2 * (r.alpha4 - 1) + (1 - B * r.alpha2) * (r.beta4 - 1);
  close(lhs, rhs, 1e-9);
  // far-downstream depth drop satisfies the overall momentum balance ½(1 − (1 − x)²)/Fr² − ½ B Ct = x/(1 − x)
  const x = r.dhOverH;
  assert.ok(x > 0 && x < 0.05);
  close((1 - (1 - x) ** 2) / (2 * Fr * Fr) - 0.5 * B * r.Ct, x / (1 - x), 1e-10);
  // total head loss exceeds the power extracted (mixing losses): ρ g Q Δh-equivalent ≥ P
  const headLoss = x - (Fr * Fr / 2) * (1 / (1 - x) ** 2 - 1); // (Δ total head)/h
  const powerHead = (0.5 * Fr * Fr * B * r.Cp); // P/(ρ g Q h)
  assert.ok(headLoss > powerHead, `${headLoss} vs ${powerHead}`);
  assert.ok(channelBlockageLimit({ B, Fr }).CpMax > channelBlockageLimit({ B }).CpMax);
  assert.throws(() => channelBlockage({ B: 1.2, alpha4: 0.4 }), RangeError);
});

test('blockage correction to open-water coefficients', () => {
  // no blockage: nothing to correct (Betz point: α2 + Ct/(4α2) = 2/3 + 1/3 = 1)
  const free = blockageCorrection({ Ct: 8 / 9, Cp: 16 / 27, B: 0 });
  close(free.speedRatio, 1, 1e-6);
  close(free.CpFree, 16 / 27, 1e-6);
  // a disc at the blocked optimum, corrected back, cannot exceed the Betz limit
  const lim = channelBlockageLimit({ B: 0.2 });
  const cor = blockageCorrection({ Ct: lim.Ct, Cp: lim.CpMax, B: 0.2 });
  assert.ok(cor.speedRatio < 1);
  assert.ok(cor.CpFree <= 16 / 27 + 1e-9 && cor.CpFree > 0.4);
  assert.ok(cor.CtFree < lim.Ct);
});

test('array wake and bidirectional operation', () => {
  const D = 18;
  const curve = tidalPowerCurve({ D, ratedPower: 1e6, uMax: 5 });
  const t = (x, y) => ({ x, y, D, hubHeight: 15, power: { u: curve.u, v: curve.p }, ct: { u: curve.u, v: curve.ct } });
  const ts = [t(0, 0), t(0, 180)];
  // flood toward the north: the northern turbine is in the wake of the southern one
  const flood = tidalArrayFlow({ turbines: ts, U: 2.2, flowDirDeg: 0 });
  assert.ok(flood.power[1] < flood.power[0]);
  close(flood.U[0], 2.2, 1e-12);
  const ebb = tidalArrayFlow({ turbines: ts, U: 2.2, flowDirDeg: 180 });
  close(ebb.power[0], flood.power[1], 1e-9);
  // staggering the second row by 3D removes most of the loss
  const stag = tidalArrayFlow({ turbines: [t(0, 0), t(54, 180)], U: 2.2, flowDirDeg: 0 });
  assert.ok(stag.power[1] > flood.power[1]);

  const opposed = bidirectionalYield({ floodDirDeg: 40, ebbDirDeg: 220 });
  close(opposed.loss, 0, 1e-12);
  // 20° asymmetry, equal energies: best heading bisects, each phase misaligned by 10° → cos³(10°)
  const asym = bidirectionalYield({ floodDirDeg: 40, ebbDirDeg: 240 });
  close(asym.asymmetryDeg, 20, 1e-12);
  close(asym.headingDeg, 50, 1e-2);
  close(1 - asym.loss, Math.cos((10 * Math.PI) / 180) ** 3, 1e-6);
  // heading aligned with the flood: all the loss is on the ebb
  const fixed = bidirectionalYield({ floodDirDeg: 40, ebbDirDeg: 240, headingDeg: 40, floodEnergy: 3, ebbEnergy: 1 });
  close(fixed.floodFactor, 1, 1e-12);
  close(fixed.ebbFactor, Math.cos((20 * Math.PI) / 180) ** 3, 1e-12);
  close(fixed.loss, (1 - fixed.ebbFactor) / 4, 1e-12);
});

test('annual energy: sinusoidal tide with an uncapped cubic curve has mean power ½ρA Cp U0³ · 4/(3π)', () => {
  const A = (Math.PI * 18 * 18) / 4;
  const k = 0.5 * 1025 * A * 0.4;
  const u = [];
  const p = [];
  for (let i = 0; i <= 600; i++) {
    u.push(i * 0.005);
    p.push(k * (i * 0.005) ** 3);
  }
  const n = 100000;
  const speeds = Float64Array.from({ length: n }, (_, i) => 3 * Math.cos((2 * Math.PI * (i + 0.5)) / n));
  const r = tidalEnergy({ curve: { u, p }, speeds, availability: 0.95 });
  close(r.meanPower, 0.95 * k * 27 * (4 / (3 * Math.PI)), 1e-5);
  close(r.aepWh, r.meanPower * 8760, 1e-12);
  const dist = velocityDistribution(speeds, { binWidth: 0.02 }).bins;
  close(tidalEnergy({ curve: { u, p }, distribution: dist, availability: 0.95 }).meanPower, r.meanPower, 2e-3);
  // rated cap: constant speed above rated
  const capped = tidalEnergy({ curve: { u: [0, 1, 2.5, 4], p: [0, 0, 1e6, 1e6] }, speeds: [3, -3, 3.5] });
  close(capped.capacityFactor, 1, 1e-12);
  close(capped.ratedFraction, 1, 1e-12);
});

test('extreme thrust, wave orbital velocity and support loads by hand', () => {
  const D = 18;
  const A = (Math.PI * D * D) / 4;
  const noWave = extremeThrust({ D, Ct: 0.85, springPeak: 3.5, surgeCurrent: 0.5 });
  close(noWave.thrust, 0.5 * 1025 * A * 0.85 * 16, 1e-12);
  assert.equal(noWave.thrustRange, 0);
  // deep-water wave (T = 6 s, depth 200 m): u = (πH/T) e^(−k d), k = ω²/g
  const depth = 200;
  const hub = 190;
  const k = (2 * Math.PI / 6) ** 2 / 9.80665;
  const w = extremeThrust({ D, Ct: 0.85, springPeak: 3, wave: { H: 4, T: 6, depth }, hubHeight: hub });
  close(w.waveVelocity, ((Math.PI * 4) / 6) * Math.exp(-k * 10), 1e-6);
  close(w.thrustRange, 2 * 1025 * A * 0.85 * 3 * w.waveVelocity, 1e-9);
  close(w.cycles, 600, 1e-12);
  close(w.thrust, 0.5 * 1025 * A * 0.85 * (3 + w.waveVelocity) ** 2, 1e-12);
  // column in a uniform current: drag = ½ρ Cd D U² h acting at h/2
  const s = supportLoads({ thrust: 5e5, hubHeight: 20, columnDiameter: 3, profile: () => 2.5, Cd: 0.9 });
  const drag = 0.5 * 1025 * 0.9 * 3 * 2.5 ** 2 * 20;
  close(s.columnDrag, drag, 1e-12);
  close(s.baseShear, 5e5 + drag, 1e-12);
  close(s.baseMoment, 5e5 * 20 + drag * 10, 1e-12);
  // 1/7th profile: drag ∝ ∫(z/h)^(2/7) dz = 7h/9, moment arm ∫z^(9/7)/∫z^(2/7) = 9h/16
  const prof = currentProfile({ uRef: 2.5, zRef: 20, n: 7 });
  const sp = supportLoads({ thrust: 0, hubHeight: 20, columnDiameter: 3, profile: prof, Cd: 0.9, nSeg: 4000 });
  close(sp.columnDrag, (drag * 7) / 9, 1e-5);
  close(sp.columnDragMoment / sp.columnDrag, (9 * 20) / 16, 1e-5);
  const withWave = supportLoads({ thrust: 0, hubHeight: 20, columnDiameter: 3, profile: () => 2.5, wave: { H: 3, T: 8, depth: 35 } });
  assert.ok(withWave.columnDrag > s.columnDrag && withWave.inertiaShear > 0);
});

test('flux resource and tidal-range energy', () => {
  const f = fluxResource({ meanPowerDensity: 4000, area: 50000, sif: 0.15 });
  close(f.fluxPower, 2e8, 1e-12);
  close(f.extractablePower, 3e7, 1e-12);
  close(f.extractableWh, 3e7 * 8760, 1e-12);
  // potential energy per tide ½ρgAR²: A = 10 km², R = 8 m → 0.5·1025·9.80665·1e7·64 J = 893.5 MWh
  const r = tidalRangeEbb({ basinArea: 1e7, range: 8, turbineArea: 300, sluiceArea: 600, dt: 30 });
  close(r.potentialPerCycleWh, (0.5 * 1025 * 9.80665 * 1e7 * 64) / 3600, 1e-12);
  close(r.potentialPerCycleWh, 893.49e6, 1e-4);
  // ebb generation recovers a fraction of the potential energy and only on the falling tide
  assert.ok(r.efficiency > 0.1 && r.efficiency < 0.6, `efficiency ${r.efficiency}`);
  assert.ok(r.generatingFraction > 0.15 && r.generatingFraction < 0.6);
  for (let i = 0; i < r.t.length; i++) {
    assert.ok(r.basin[i] <= 4 + 1e-9 && r.basin[i] >= -4 - 1e-9);
    if (r.power[i] > 0) assert.ok(r.basin[i] - r.sea[i] >= 1 - 1e-9);
  }
  // energy bookkeeping: generated energy equals η ρ g ∫ Q H dt, and cannot exceed η × potential
  assert.ok(r.energyPerCycleWh < 0.85 * r.potentialPerCycleWh);
  const capped = tidalRangeEbb({ basinArea: 1e7, range: 8, turbineArea: 300, sluiceArea: 600, dt: 30, ratedPower: 100e6 });
  assert.ok(capped.peakPower <= 100e6 * (1 + 1e-12));
  assert.ok(r.peakPower > 100e6);
  close(r.annualWh / r.energyPerCycleWh, (8760 * 3600) / 44714.16, 1e-9);
});
