import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  unbalanceForce, transmissibility, sdofForcedResponse, isolatorSelection, twoDofFoundation, equivalentFootingRadius, footingImpedance,
  separationMargin, vibrationSeverity, harmonicAmplitudes, campbellData,
} from '../../src/engine/structural/machinery.js';
import { solveDense } from '../../src/engine/math/linalg.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);

test('unbalance force from the balance grade e·Ω = G: 1000 kg rotor, G 6.3, 3000 rpm → e = 20.05 µm, F = 1979 N (hand calculation)', () => {
  // Ω = 3000·2π/60 = 314.16 rad/s; e = 6.3e-3/314.16 = 2.005e-5 m; F = m e Ω² = m G Ω = 1000·6.3e-3·314.16 = 1979.2 N
  const r = unbalanceForce({ mass: 1000, gradeMmS: 6.3, speedRpm: 3000 });
  close(r.omega, 314.159265, 1e-8); close(r.eccentricity, 2.0054e-5, 1e-4); close(r.force, 1979.2, 1e-4);
  close(r.unbalance, 1000 * r.eccentricity, 1e-12); close(r.frequencyHz, 50, 1e-12);
  // G 2.5 at 1500 rpm with a service factor of 2: F = 2·1000·2.5e-3·157.08 = 785.4 N
  close(unbalanceForce({ mass: 1000, gradeMmS: 2.5, speedRpm: 1500, serviceFactor: 2 }).force, 785.4, 1e-4);
  // eccentricity given directly: F = m e ω²
  close(unbalanceForce({ mass: 50, eccentricity: 1e-4, omega: 100 }).force, 50, 1e-12);
  // for a fixed grade the force grows linearly with speed
  close(unbalanceForce({ mass: 1000, gradeMmS: 6.3, speedRpm: 6000 }).force / r.force, 2, 1e-12);
  assert.throws(() => unbalanceForce({ mass: 1000, gradeMmS: 6.3 }), RangeError);
});

test('SDOF magnification and transmissibility: D = 1/(2ζ) at resonance, TR = 1 at r = √2 for any damping, TR = 1/(r² − 1) undamped', () => {
  close(transmissibility({ r: 1, zeta: 0.05 }).magnification, 10, 1e-12);
  for (const z of [0, 0.1, 0.5]) close(transmissibility({ r: Math.SQRT2, zeta: z }).transmissibility, 1, 1e-12);
  close(transmissibility({ r: 3 }).transmissibility, 1 / 8, 1e-12);
  close(transmissibility({ r: 3 }).isolationEfficiency, 0.875, 1e-12);
  // r = 3, ζ = 0.1: TR = √(1 + 0.36)/√(64 + 0.36) = 0.14537 — damping reduces isolation above √2
  close(transmissibility({ r: 3, zeta: 0.1 }).transmissibility, Math.sqrt(1.36 / 64.36), 1e-12);
  close(transmissibility({ r: 0, zeta: 0.1 }).magnification, 1, 1e-12);
  close(transmissibility({ r: 1, zeta: 0.05 }).phase, Math.PI / 2, 1e-12);
});

test('SDOF forced response: rotating unbalance X = (m_u e/m) r² D; transmitted force; static deflection', () => {
  // m = 2000 kg on k = 8e6 N/m: ωn = 63.25 rad/s (10.07 Hz); machine at 1500 rpm (157.08 rad/s): r = 2.4836
  const o = { mass: 2000, k: 8e6, zeta: 0.05, speedRpm: 1500, unbalance: 0.02 };
  const r = sdofForcedResponse(o);
  close(r.omegaN, Math.sqrt(4000), 1e-12); close(r.r, 157.0796 / 63.2456, 1e-6);
  const D = 1 / Math.sqrt((1 - r.r ** 2) ** 2 + (0.1 * r.r) ** 2);
  close(r.amplitude, 0.02 / 2000 * r.r ** 2 * D, 1e-12);
  close(r.F0, 0.02 * 157.0796 ** 2, 1e-6);
  close(r.transmittedForce, r.F0 * Math.sqrt(1 + (0.1 * r.r) ** 2) * D, 1e-12);
  close(r.velocityRms, r.amplitude * 157.0796 / Math.SQRT2, 1e-6);
  close(r.staticDeflection, 2000 * 9.80665 / 8e6, 1e-12);
  // well above resonance the amplitude tends to m_u e/m (mass-controlled)
  close(sdofForcedResponse({ ...o, zeta: 0, speedRpm: 60000 }).amplitude, 1e-5, 1e-3);
  // damping given as c: ζ = c/(2 m ωn)
  close(sdofForcedResponse({ mass: 2000, k: 8e6, c: 2 * 2000 * Math.sqrt(4000) * 0.05, speedRpm: 1500, F0: 100 }).zeta, 0.05, 1e-12);
  // static limit
  close(sdofForcedResponse({ mass: 2000, k: 8e6, omega: 1e-3, F0: 100 }).amplitude, 100 / 8e6, 1e-6);
});

test('isolator selection: 90 % isolation at 25 Hz needs fn = 25/√11 = 7.54 Hz and 4.4 mm static deflection', () => {
  const r = isolatorSelection({ efficiency: 0.9, frequencyHz: 25, mass: 500 });
  close(r.r, Math.sqrt(11), 1e-12); close(r.naturalHz, 7.5378, 1e-4);
  close(r.staticDeflection, 9.80665 / (2 * Math.PI * 7.5378) ** 2, 1e-4); close(r.staticDeflection, 4.372e-3, 1e-3);
  close(r.stiffness, 500 * (2 * Math.PI * r.naturalHz) ** 2, 1e-12);
  close(transmissibility({ r: r.r }).isolationEfficiency, 0.9, 1e-12);
});

test('2-DOF machine/foundation: undamped amplitudes against a direct 2×2 solve; natural frequencies; rigid-foundation limit; anti-resonance', () => {
  const p = { m1: 1500, k1: 6e6, m2: 12000, k2: 4e8, F0: 2000 };
  const w = 2 * Math.PI * 25;
  const r = twoDofFoundation({ ...p, omega: w });
  const x = solveDense([[p.k1 - p.m1 * w * w, -p.k1], [-p.k1, p.k1 + p.k2 - p.m2 * w * w]], [p.F0, 0]);
  close(r.x1, Math.abs(x[0]), 1e-10); close(r.x2, Math.abs(x[1]), 1e-10);
  close(r.forceToGround, p.k2 * Math.abs(x[1]), 1e-10); close(r.forceInIsolators, p.k1 * Math.abs(x[0] - x[1]), 1e-10);
  // natural frequencies satisfy det(K − λM) = 0 and λ1 λ2 = k1 k2/(m1 m2)
  for (const wn of r.naturalOmegas) {
    const l = wn * wn;
    const det = (p.k1 - p.m1 * l) * (p.k1 + p.k2 - p.m2 * l) - p.k1 * p.k1;
    assert.ok(Math.abs(det) < 1e-9 * p.k1 * p.k2);
  }
  close((r.naturalOmegas[0] * r.naturalOmegas[1]) ** 2, p.k1 * p.k2 / (p.m1 * p.m2), 1e-10);
  // very heavy and stiff foundation: machine behaves as an SDOF on its isolators
  const rigid = twoDofFoundation({ ...p, m2: 1e12, k2: 1e16, c1: 2 * 0.05 * Math.sqrt(p.k1 * p.m1), omega: w });
  const sd = sdofForcedResponse({ mass: p.m1, k: p.k1, zeta: 0.05, omega: w, F0: p.F0 });
  close(rigid.x1, sd.amplitude, 1e-6); close(rigid.forceInIsolators, sd.transmittedForce, 1e-6);
  // machine motion vanishes when the foundation is tuned to ω² = (k1 + k2)/m2 (it acts as a vibration absorber)
  const anti = twoDofFoundation({ ...p, omega: Math.sqrt((p.k1 + p.k2) / p.m2) });
  assert.ok(anti.x1 < 1e-12 * anti.x2);
  // damping bounds the response at the first natural frequency; power balance: input power = dissipated power
  const d = twoDofFoundation({ ...p, c1: 2e4, c2: 5e5, omega: r.naturalOmegas[0] });
  assert.ok(Number.isFinite(d.x1) && d.x1 > 10 * r.x1);
  const wd = r.naturalOmegas[0];
  const pin = 0.5 * p.F0 * wd * d.x1 * Math.sin(d.phase1);
  const pdis = 0.5 * 2e4 * (wd * d.relative) ** 2 + 0.5 * 5e5 * (wd * d.x2) ** 2;
  close(pin, pdis, 1e-9);
  assert.throws(() => twoDofFoundation({ ...p, omega: r.naturalOmegas[0] }), RangeError);
});

test('rigid footing on an elastic half-space (Richart, Hall & Woods): spring constants and damping ratios by hand', () => {
  // G = 50 MPa, ν = 0.3, ρ = 1800 kg/m³, r = 2 m, m = 80 t
  const r = footingImpedance({ G: 50e6, nu: 0.3, rho: 1800, radius: 2, mass: 80e3, inertiaRocking: 2.0e5, inertiaTorsion: 1.6e5 });
  close(r.kz, 4 * 50e6 * 2 / 0.7, 1e-12); close(r.kz, 5.714e8, 1e-3);
  close(r.kx, 8 * 50e6 * 2 / 1.7, 1e-12); close(r.kRocking, 8 * 50e6 * 8 / 2.1, 1e-12); close(r.kTorsion, 16 * 50e6 * 8 / 3, 1e-12);
  // Bz = 0.7/4·80e3/(1800·8) = 0.9722 → Dz = 0.425/√0.9722 = 0.431
  close(r.Bz, 0.97222, 1e-4); close(r.Dz, 0.4310, 1e-3);
  // Bx = (7 − 2.4)/(32·0.7)·5.5556 = 1.1409 → Dx = 0.288/√1.1409 = 0.2696
  close(r.Bx, 1.14087, 1e-4); close(r.Dx, 0.2696, 1e-3);
  // Bψ = 3·0.7/8·2e5/(1800·32) = 0.9115 → Dψ = 0.15/(1.9115·0.9547) = 0.0822
  close(r.BRocking, 0.91146, 1e-4); close(r.DRocking, 0.0822, 2e-3);
  // Bθ = 1.6e5/(1800·32) = 2.7778 → Dθ = 0.5/6.5556 = 0.0763
  close(r.BTorsion, 2.77778, 1e-5); close(r.DTorsion, 0.07627, 1e-3);
  // vertical natural frequency of the block on the half-space: f = √(kz/m)/2π = 13.45 Hz
  close(Math.sqrt(r.kz / 80e3) / (2 * Math.PI), 13.45, 1e-3);
  // rectangular base 3 × 5 m: equal-area and equal-inertia circles
  const e = equivalentFootingRadius({ B: 3, L: 5 });
  close(Math.PI * e.translation ** 2, 15, 1e-12);
  close(Math.PI * e.rocking ** 4 / 4, 3 * 125 / 12, 1e-12);
  close(Math.PI * e.torsion ** 4 / 2, 3 * 5 * (9 + 25) / 12, 1e-12);
});

test('separation margin, severity zones against supplied limits, harmonic amplitude conversions', () => {
  const s = separationMargin({ excitationHz: 25, naturalHz: 18 });
  close(s.margin, 7 / 18, 1e-12); assert.ok(s.ok); assert.equal(s.tuning, 'low');
  close(s.avoidBandHz[0], 25 / 1.2, 1e-12); close(s.avoidBandHz[1], 25 / 0.8, 1e-12);
  assert.equal(separationMargin({ excitationHz: 25, naturalHz: 27 }).ok, false);
  assert.equal(separationMargin({ excitationHz: 25, naturalHz: 27, required: 0.05 }).ok, true);
  // natural frequencies at the edges of the avoid band have exactly the required margin
  close(separationMargin({ excitationHz: 25, naturalHz: 25 / 0.8 }).margin, 0.2, 1e-12);
  // zone boundaries are inputs (values here are arbitrary test numbers in m/s)
  const limits = { AB: 2e-3, BC: 4e-3, CD: 8e-3 };
  assert.equal(vibrationSeverity({ velocityRms: 1e-3, limits }).zone, 'A');
  assert.equal(vibrationSeverity({ velocityRms: 3e-3, limits }).zone, 'B');
  const c = vibrationSeverity({ velocityRms: 6e-3, limits });
  assert.equal(c.zone, 'C'); assert.equal(c.acceptable, false); close(c.utilisation, 1.5, 1e-12);
  assert.equal(vibrationSeverity({ velocityRms: 9e-3, limits }).zone, 'D');
  assert.throws(() => vibrationSeverity({ velocityRms: 1e-3 }), RangeError);
  // 50 Hz, 20 µm zero-to-peak: v = 2π·50·20e-6 = 6.283 mm/s peak = 4.443 mm/s rms; a = 1.974 m/s²
  const h = harmonicAmplitudes({ frequencyHz: 50, displacement: 20e-6 });
  close(h.velocity, 6.2832e-3, 1e-4); close(h.velocityRms, 4.4429e-3, 1e-4); close(h.acceleration, 1.9739, 1e-4);
  close(harmonicAmplitudes({ frequencyHz: 50, velocity: h.velocity }).displacementPeakToPeak, 40e-6, 1e-12);
  close(harmonicAmplitudes({ frequencyHz: 50, acceleration: h.acceleration }).displacement, 20e-6, 1e-12);
});

test('Campbell data: critical speeds at f_n·60/n for constant modes; speed-dependent whirl pair; operating-range margins', () => {
  const r = campbellData({ modes: [20, 55], orders: [1, 2], speedMinRpm: 0, speedMaxRpm: 4000, operatingRangeRpm: [1400, 1600] });
  const find = (mode, order) => r.crossings.find((c) => c.mode === mode && c.order === order);
  close(find(0, 1).speedRpm, 1200, 1e-9); close(find(0, 2).speedRpm, 600, 1e-9);
  close(find(1, 1).speedRpm, 3300, 1e-9); close(find(1, 2).speedRpm, 1650, 1e-9);
  assert.equal(r.crossings.length, 4);
  assert.equal(find(1, 2).inOperatingRange, false);
  // margins in the operating range: mode 0 vs 1× is closest at 1400 rpm (23.33 Hz): 16.7 %; mode 1 vs 2× at 1600 rpm (53.33 Hz): 3.0 %
  const m = (mode, order) => r.margins.find((c) => c.mode === mode && c.order === order);
  close(m(0, 1).margin, (1400 / 60 - 20) / 20, 1e-9); assert.equal(m(0, 1).ok, false);
  close(m(1, 2).margin, (55 - 2 * 1600 / 60) / 55, 1e-9);
  assert.ok(m(1, 1).ok);
  close(r.orderHz[1][100], 2 * 4000 / 60, 1e-12); assert.equal(r.modeHz[0][50], 20);
  // forward / backward whirl of a rotor with gyroscopic splitting f = f0 ± a·N: crossings with 1× at N = 60 f0/(1 ∓ 60a)
  const f0 = 30; const a = 0.002;
  const w = campbellData({ modes: [(n) => f0 - a * n, (n) => f0 + a * n], orders: [1], speedMinRpm: 0, speedMaxRpm: 3000, operatingRangeRpm: [1500, 1700] });
  close(w.crossings.find((c) => c.mode === 0).speedRpm, 60 * f0 / (1 + 60 * a), 1e-9);
  close(w.crossings.find((c) => c.mode === 1).speedRpm, 60 * f0 / (1 - 60 * a), 1e-9);
  // a crossing inside the operating range gives zero margin
  assert.equal(w.margins.find((c) => c.mode === 0).margin, 0);
  assert.throws(() => campbellData({ modes: [], speedMinRpm: 0, speedMaxRpm: 100 }), RangeError);
});
