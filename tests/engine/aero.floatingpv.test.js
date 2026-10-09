import test from 'node:test';
import assert from 'node:assert/strict';
import {
  floatBuoyancy, islandStability, panelWindLoad, floatWaveLoad, meanDriftForce, currentDrag, combineLoads,
  mooringLoadSharing, mooringLevelRange, anchorCount, hingedRaftResponse, evaporationSaving, coolingYieldGain,
  floatingCableLength, exposureClass, FLOATING_PV_ELECTRICAL_NOTES,
} from '../../src/engine/aero/floatingpv.js';

const G = 9.80665;
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);

test('float buoyancy and freeboard by hand (Archimedes)', () => {
  // float 0.25 m³, 0.25 m high → waterplane 1 m²; dead 60 kg → draft 0.06 m, freeboard 0.19 m, reserve 250/60 − 1
  const r = floatBuoyancy({ floatVolume: 0.25, floatHeight: 0.25, deadMass: 60, liveMass: 90, snowMass: 50, minFreeboard: 0.06 });
  close(r.waterplaneArea, 1, 1e-12);
  close(r.cases[0].draft, 0.06, 1e-12);
  close(r.cases[0].freeboard, 0.19, 1e-12);
  close(r.cases[0].reserveBuoyancy, 250 / 60 - 1, 1e-12);
  close(r.cases[1].draft, 0.15, 1e-12);
  close(r.cases[2].freeboard, 0.05, 1e-12);
  assert.equal(r.cases[1].ok, true);
  assert.equal(r.cases[2].ok, false);
  close(r.cases[2].utilisation, 200 / 250, 1e-12);
  // equilibrium: weight = displaced water weight
  close(1000 * G * r.waterplaneArea * r.cases[1].draft, 150 * G, 1e-12);
});

test('island stability: box-barge metacentric height and heel by hand', () => {
  // 20 m × 6 m raft, full waterplane, 12 t, draft 0.1 m: I = 20·216/12 = 360 m⁴, ∇ = 12 m³, BM = 30 m
  const r = islandStability({ mass: 12000, kg: 0.5, draft: 0.1, length: 20, width: 6, windForce: 6000, windLever: 0.8, personMass: 200, personOffset: 2.5 });
  close(r.waterplaneInertia, 360, 1e-12);
  close(r.BM, 30, 1e-12);
  close(r.GM, 0.05 + 30 - 0.5, 1e-12);
  const M = 6000 * 0.8 + 200 * G * 2.5;
  close(r.heelingMoment, M, 1e-12);
  close(r.heelDeg, (Math.asin(M / (12000 * G * 29.55)) * 180) / Math.PI, 1e-12);
  // two pontoons 1 m wide at ±2.5 m: I = 2 (20·2.5² + 20·1/12) = 253.33 m⁴ (parallel-axis theorem)
  const cat = islandStability({ mass: 12000, kg: 0.5, draft: 0.3, floats: [{ y: 2.5, area: 20, ownInertia: 20 / 12 }, { y: -2.5, area: 20, ownInertia: 20 / 12 }] });
  close(cat.waterplaneInertia, 2 * (125 + 20 / 12), 1e-12);
  assert.equal(islandStability({ mass: 12000, kg: 40, draft: 0.1, length: 20, width: 6 }).stable, false);
});

test('wind load on tilted rows: normal force resolved into drag and lift, sheltered inner rows', () => {
  // q = ½·1.225·30² = 551.25 Pa; F_N = 551.25·1.2·40 = 26460 N; tilt 12°
  const r = panelWindLoad({ windSpeed: 30, panelArea: 40, tiltDeg: 12, cN: 1.2, rows: 5, shelterFactor: 0.4 });
  close(r.dynamicPressure, 551.25, 1e-12);
  close(r.normalForcePerRow, 26460, 1e-12);
  close(r.dragFirstRow, 26460 * Math.sin((12 * Math.PI) / 180), 1e-12);
  close(r.drag, r.dragFirstRow * (1 + 4 * 0.4), 1e-12);
  close(Math.hypot(r.dragFirstRow, r.liftFirstRow), 26460, 1e-12);
  // the pressure coefficient is a mandatory project input
  assert.throws(() => panelWindLoad({ windSpeed: 30, panelArea: 40, tiltDeg: 12 }), RangeError);
});

test('wave, drift and current loads by hand', () => {
  // deep water, H = 1 m, T = 4 s: u = πH/T = 0.7854 m/s, a = ω u = 1.2337 m/s²
  const w = floatWaveLoad({ H: 1, T: 4, submergedVolume: 0.1, projectedArea: 0.1, Cm: 2, Cd: 1 });
  close(w.velocity, Math.PI / 4, 1e-9);
  close(w.acceleration, Math.PI ** 2 / 8, 1e-9);
  close(w.inertia, 1000 * 2 * 0.1 * Math.PI ** 2 / 8, 1e-9);
  close(w.drag, 0.5 * 1000 * 0.1 * (Math.PI / 4) ** 2, 1e-9);
  close(w.wavelength, (G * 16) / (2 * Math.PI), 1e-9);
  // F_I ≥ 2 F_D → the inertia amplitude governs
  assert.ok(w.inertia > 2 * w.drag);
  close(w.max, w.inertia, 1e-12);
  // drag-dominated case: F_D + F_I²/(4 F_D)
  const d = floatWaveLoad({ H: 1, T: 4, submergedVolume: 0.01, projectedArea: 1, Cm: 2, Cd: 1 });
  close(d.max, d.drag + d.inertia ** 2 / (4 * d.drag), 1e-12);
  // direct check of the Morison maximum by scanning the phase
  let mx = 0;
  for (let i = 0; i < 20000; i++) {
    const th = (i / 20000) * 2 * Math.PI;
    mx = Math.max(mx, d.drag * Math.cos(th) * Math.abs(Math.cos(th)) + d.inertia * Math.sin(th));
  }
  close(d.max, mx, 1e-6);
  // mean drift: ρ g R² Hs²/16 per metre: 1000·9.80665·0.25·1/16·50 = 7661.4 N
  close(meanDriftForce({ Hs: 1, width: 50, reflection: 0.5 }), (1000 * G * 0.25 * 50) / 16, 1e-12);
  close(meanDriftForce({ Hs: 1, width: 50, reflection: 1, regular: true }), 0.5 * 1000 * G * 0.25 * 50, 1e-12);
  close(currentDrag({ speed: 0.8, area: 12, Cd: 0.9 }), 0.5 * 1000 * 0.9 * 12 * 0.64, 1e-12);
  const c = combineLoads([{ force: 3000, dirDeg: 0 }, { force: 4000, dirDeg: 90 }]);
  close(c.magnitude, 5000, 1e-12);
  close(c.dirDeg, (Math.atan2(4000, 3000) * 180) / Math.PI, 1e-12);
});

test('mooring load sharing: equilibrium, slack leeward lines and symmetric sharing', () => {
  // four lines N, E, S, W; load toward the north is taken by the southern line alone (others slack or unloaded)
  const lines = [0, 90, 180, 270].map((azimuthDeg) => ({ azimuthDeg, stiffness: 5e4 }));
  const r = mooringLoadSharing({ Fx: 0, Fy: 1e4, lines });
  close(r.tensions[2], 1e4, 1e-9);
  close(r.offsetY, 1e4 / 5e4, 1e-9);
  assert.equal(r.slackLines, 1);
  // force equilibrium Σ T_i n_i + F = 0
  const check = (res, ls, Fx, Fy) => {
    let sx = Fx;
    let sy = Fy;
    ls.forEach((l, i) => {
      sx += res.tensions[i] * Math.sin((l.azimuthDeg * Math.PI) / 180);
      sy += res.tensions[i] * Math.cos((l.azimuthDeg * Math.PI) / 180);
    });
    assert.ok(Math.hypot(sx, sy) < 1e-6 * Math.hypot(Fx, Fy), `residual ${Math.hypot(sx, sy)}`);
  };
  check(r, lines, 0, 1e4);
  // load at 45°: the two windward lines (S and W) share equally, F/√2 each
  const d = mooringLoadSharing({ Fx: 1e4 * Math.SQRT1_2, Fy: 1e4 * Math.SQRT1_2, lines });
  close(d.tensions[2], 1e4 * Math.SQRT1_2, 1e-9);
  close(d.tensions[3], 1e4 * Math.SQRT1_2, 1e-9);
  assert.equal(d.loadedLines, 2);
  check(d, lines, 1e4 * Math.SQRT1_2, 1e4 * Math.SQRT1_2);
  // with pretension all lines stay taut and the offset halves (two springs in parallel per axis)
  const pre = lines.map((l) => ({ ...l, pretension: 2e4 }));
  const p = mooringLoadSharing({ Fx: 0, Fy: 1e4, lines: pre });
  assert.equal(p.slackLines, 0);
  close(p.offsetY, 1e4 / 1e5, 1e-9);
  close(p.tensions[2] - p.tensions[0], 1e4, 1e-9);
  check(p, pre, 0, 1e4);
  assert.throws(() => mooringLoadSharing({ Fx: 1e4, Fy: 0, lines: [{ azimuthDeg: 0, stiffness: 1e4 }, { azimuthDeg: 180, stiffness: 1e4 }] }), RangeError);
});

test('water-level range geometry and anchor count by hand', () => {
  // anchor 30 m away, depth 5 m at low water and 13 m at high water
  const r = mooringLevelRange({ anchorDistance: 30, depthMin: 5, depthMax: 13, load: 2e4, allowableOffset: 2 });
  close(r.length, Math.hypot(30, 13), 1e-12);
  close(r.slackAtLowWater, Math.hypot(30, 13) - Math.hypot(30, 5), 1e-12);
  close(r.excursionAtLowWater, Math.sqrt(30 * 30 + 169 - 25) - 30, 1e-12);
  close(r.elasticStrain, Math.hypot(30, 13) / Math.hypot(30, 5) - 1, 1e-12);
  close(r.requiredStiffness, 1e4, 1e-12);
  // 480 kN design load, γ = 1.5, 100 kN anchors, half of them effective → ceil(720/50) = 15
  const a = anchorCount({ load: 480e3, anchorCapacity: 100e3, loadFactor: 1.5, sharing: 0.5 });
  assert.equal(a.count, 15);
  close(a.utilisation, 720 / 15 / 50, 1e-12);
  assert.ok(a.utilisation <= 1);
});

test('hinged raft on a regular wave: kinematic limits', () => {
  // long waves: module pitch → wave slope, relative rotation → k² a L, no end mismatch
  const long = hingedRaftResponse({ H: 1, wavelength: 400, moduleLength: 2, moduleWidth: 2 });
  const k = (2 * Math.PI) / 400;
  close(long.waveSlope, k * 0.5, 1e-12);
  close(long.modulePitch, k * 0.5, 1e-3);
  close(long.relativeRotation, k * k * 0.5 * 2, 1e-3);
  // end mismatch vanishes as (8/15) a u³, u = kL/2
  close(long.mismatch / ((8 / 15) * 0.5 * k ** 3), 1, 1e-3);
  assert.equal(long.regime, 'contouring');
  // L = λ/2: maximum articulation 2H/L
  const half = hingedRaftResponse({ H: 0.5, wavelength: 8, moduleLength: 4, moduleWidth: 2 });
  close(half.relativeRotation, (2 * 0.5) / 4, 1e-12);
  assert.equal(half.regime, 'articulating');
  // L = λ: every module sees the same wave phase → no relative rotation
  const full = hingedRaftResponse({ H: 0.5, wavelength: 4, moduleLength: 4, moduleWidth: 2 });
  assert.ok(full.relativeRotation < 1e-12);
  // independent check of the free-module fit by numerical least squares at the phase of maximum mismatch
  const L = 4;
  const lam = 10;
  const a = 0.3;
  const kk = (2 * Math.PI) / lam;
  const fit = (xc) => {
    let s0 = 0;
    let s1 = 0;
    const n = 20000;
    for (let i = 0; i < n; i++) {
      const s = -L / 2 + ((i + 0.5) / n) * L;
      const eta = a * Math.cos(kk * (xc + s));
      s0 += eta / n;
      s1 += (eta * s) / n;
    }
    return { heave: s0, slope: (12 * s1) / (L * L) };
  };
  let mis = 0;
  for (let j = 0; j < 400; j++) {
    const xh = (j / 400) * lam;
    const left = fit(xh - L / 2);
    const right = fit(xh + L / 2);
    mis = Math.max(mis, Math.abs(left.heave + (left.slope * L) / 2 - (right.heave - (right.slope * L) / 2)));
  }
  const r = hingedRaftResponse({ H: 2 * a, wavelength: lam, moduleLength: L, moduleWidth: 2 });
  close(r.mismatch, mis, 1e-4);
  close(r.connectorShear, (1000 * G * 2 * L * r.mismatch) / 8, 1e-12);
  assert.equal(hingedRaftResponse({ H: 1, wavelength: 2, moduleLength: 4, moduleWidth: 2 }).regime, 'bridging');
});

test('water saving, cooling gain, cable length and exposure screening', () => {
  // 1.4 m/yr over 2 ha with 70 % reduction = 19 600 m³/yr
  close(evaporationSaving({ evaporationRate: 1.4, coveredArea: 20000, reduction: 0.7 }), 19600, 1e-12);
  // 5 K cooler at −0.37 %/K → +1.85 %
  close(coolingYieldGain({ deltaTempC: -5, gammaPdc: -0.0037 }), 0.0185, 1e-12);
  const c = floatingCableLength({ horizontalDistance: 40, offset: 5, levelRange: 6, slackFactor: 1.1 });
  close(c.length, 1.1 * Math.hypot(45, 6), 1e-12);
  close(c.slack, c.length - 40, 1e-12);
  assert.equal(exposureClass({ Hs: 0.4 }).class, 'inland');
  assert.equal(exposureClass({ Hs: 1 }).class, 'nearshore');
  assert.equal(exposureClass({ Hs: 3.5 }).class, 'offshore');
  assert.equal(exposureClass({ Hs: 2, thresholds: [0.5, 1.5, 2.5], labels: ['a', 'b', 'c', 'd'] }).index, 2);
  assert.ok(FLOATING_PV_ELECTRICAL_NOTES.length >= 3);
});
