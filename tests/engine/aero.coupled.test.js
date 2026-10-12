import test from 'node:test';
import assert from 'node:assert/strict';
import {
  waveNumber, jonswapShape, macCamyFuchsCm, linearIrregularWave, buildSupportStructure, monopileTowerSegments,
  seriesStats, rainflowCount, damageEquivalentLoad, logDecrement, identifyPeaks, simulateCoupled, aeroDampingRatio,
} from '../../src/engine/aero/coupled.js';
import { cpSurface, surfaceLookup } from '../../src/engine/aero/bem.js';
import { controllerFromRotor } from '../../src/engine/aero/control.js';
import { nrel5mw } from '../../src/engine/aero/rotors.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} +- ${tol}`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${tol})`);
const variance = (x) => { let m = 0; for (const v of x) m += v; m /= x.length; let s = 0; for (const v of x) s += (v - m) ** 2; return s / x.length; };

const rotor = nrel5mw();
const surface = cpSurface({ rotor, bem: { tiltDeg: 0, preconeDeg: 0 } });
const cfg = controllerFromRotor({ rotor, surface });
const turbine = { inertia: rotor.drivetrainInertia, hubHeight: 100 };
const A = Math.PI * surface.R ** 2;
// one supplied fore-aft mode: 0.28 Hz, modal mass 4e5 kg, unit hub displacement
const oneMode = (zeta = 0.01) => {
  const m = 4e5; const k = m * (2 * Math.PI * 0.28) ** 2;
  return { zHub: 100, zMudline: -20, zTowerBase: 10, modes: [{ mass: m, stiffness: k, zeta, phiHub: 1, slopeHub: 0.015, phiStations: [] }], stations: [] };
};

test('waveNumber: deep-water and shallow-water limits of the dispersion relation', () => {
  rel(waveNumber(1.0, 5000), 1 / 9.80665, 1e-10);
  rel(waveNumber(0.05, 5), 0.05 / Math.sqrt(9.80665 * 5), 2e-3);
  const k = waveNumber(0.8, 25);
  close(9.80665 * k * Math.tanh(k * 25), 0.64, 1e-12);
});

test('MacCamy-Fuchs: CM -> 2 for slender cylinders; ka = 1 by hand from Bessel tables (|CM| = 1.3716)', () => {
  close(macCamyFuchsCm(1e-4).magnitude, 2, 1e-9);
  // CM rises marginally above 2 before falling: 2 / (1 + (ka)^2/2 ln(ka/2) + ...)
  close(macCamyFuchsCm(0.05).magnitude, 2, 1e-2);
  assert.ok(macCamyFuchsCm(0.05).magnitude >= 2);
  // J1'(1) = J0(1) - J1(1) = 0.7651977 - 0.4400506 = 0.3251471 ; Y1'(1) = Y0(1) - Y1(1) = 0.0882570 + 0.7812128 = 0.8694698
  // |CM| = 4 / (pi sqrt(0.3251471^2 + 0.8694698^2)) = 1.3716 ; lag = atan(J1'/Y1') = 20.50 deg
  const c = macCamyFuchsCm(1);
  close(c.magnitude, 4 / (Math.PI * Math.hypot(0.3251471, 0.8694698)), 1e-6);
  close(c.magnitude, 1.3716, 1e-4);
  close(c.phaseLag, Math.atan(0.3251471 / 0.8694698), 1e-6);
  assert.ok(c.im < 0, 'force lags the acceleration');
  assert.ok(macCamyFuchsCm(2).magnitude < c.magnitude);
});

test('linearIrregularWave: m0 = Hs^2/16, JONSWAP peak at Tp, Airy kinematics decay with depth', () => {
  const w = linearIrregularWave({ Hs: 3, Tp: 9, depth: 30, dt: 0.25, n: 8192, seed: 4, z: [0, -10, -30], diameters: [6, 6, 6] });
  rel(variance(w.eta), 9 / 16, 1e-9);
  rel(w.m0, 9 / 16, 1e-12);
  assert.ok(jonswapShape(1 / 9, 9) > jonswapShape(0.9 / 9, 9) && jonswapShape(1 / 9, 9) > jonswapShape(1.1 / 9, 9));
  // gamma = 1 gives the Pierson-Moskowitz shape: peak enhancement ratio at fp is gamma
  rel(jonswapShape(1 / 9, 9, 3.3) / jonswapShape(1 / 9, 9, 1), 3.3, 1e-12);
  const su = w.u.map((u) => Math.sqrt(variance(u)));
  assert.ok(su[0] > su[1] && su[1] > su[2] && su[2] > 0);
  // acceleration is the time derivative of velocity (central difference at mid-record)
  const i = 4000;
  close(w.acc[1][i], (w.u[1][i + 1] - w.u[1][i - 1]) / 0.5, 0.05 * Math.sqrt(variance(w.acc[1])));
  // diffraction reduces the inertia loading relative to CM = 2
  assert.ok(variance(w.accMF[0]) < 4 * variance(w.acc[0]));
  // seeded
  const w2 = linearIrregularWave({ Hs: 3, Tp: 9, depth: 30, dt: 0.25, n: 8192, seed: 4 });
  assert.equal(w2.eta[77], w.eta[77]);
});

test('buildSupportStructure: uniform cantilever frequencies (beta L = 1.8751, 4.6941) and tip-mass Rayleigh check', () => {
  const D = 5; const t = 0.05; const L = 80; const E = 210e9; const rho = 7850;
  const seg = [{ zBase: 0, zTop: L, dBase: D, dTop: D, tBase: t, tTop: t, density: rho, E }];
  const s = buildSupportStructure({ segments: seg, rnaMass: 0, zHub: L, nModes: 2, nElem: 40 });
  const EI = E * Math.PI / 64 * (D ** 4 - (D - 2 * t) ** 4);
  const m = rho * Math.PI * (D - t) * t;
  rel(s.modes[0].frequency, 1.875104 ** 2 / (2 * Math.PI) * Math.sqrt(EI / (m * L ** 4)), 1e-5);
  rel(s.modes[1].frequency, 4.694091 ** 2 / (2 * Math.PI) * Math.sqrt(EI / (m * L ** 4)), 1e-4);
  close(s.modes[0].phiHub, 1, 1e-12);
  rel(s.mass, m * L, 1e-12);
  // heavy tip mass: f -> sqrt(3 EI / (L^3 (M + 0.2357 m L))) / 2 pi  (Rayleigh, static deflection shape: 33/140 m L)
  const M = 5e5;
  const s2 = buildSupportStructure({ segments: seg, rnaMass: M, zHub: L, nModes: 1 });
  rel(s2.modes[0].frequency, Math.sqrt(3 * EI / (L ** 3 * (M + 33 / 140 * m * L))) / (2 * Math.PI), 2e-3);
  // foundation flexibility lowers the frequency; a very stiff foundation recovers the clamped value
  const soft = buildSupportStructure({ segments: seg, rnaMass: M, zHub: L, nModes: 1, foundation: { Kh: 1e9, Kr: 1e11 } });
  assert.ok(soft.modes[0].frequency < s2.modes[0].frequency);
  const stiff = buildSupportStructure({ segments: seg, rnaMass: M, zHub: L, nModes: 1, foundation: { Kh: 1e12, Kr: 1e14 } });
  rel(stiff.modes[0].frequency, s2.modes[0].frequency, 2e-4);
  assert.ok(stiff.modes[0].frequency <= s2.modes[0].frequency * (1 + 1e-9));
});

test('buildSupportStructure: modal section forces satisfy statics (shear = omega^2 x first mass moment of the mode)', () => {
  const seg = monopileTowerSegments({ tower: { ...rotor.tower, height: 77.6 }, waterDepth: 20, transitionZ: 10, monopileD: 6, monopileT: 0.06 });
  const s = buildSupportStructure({ segments: seg, rnaMass: 350000, zHub: 90, nModes: 2, nElem: 30 });
  assert.equal(s.zMudline, -20); assert.equal(s.zTowerBase, 10); close(s.zTop, 87.6, 1e-12);
  assert.ok(s.modes[0].frequency > 0.2 && s.modes[0].frequency < 0.4, `first fore-aft frequency ${s.modes[0].frequency}`);
  assert.ok(s.modes[1].frequency > 4 * s.modes[0].frequency);
  // wet stations below still water, dry above
  for (const st of s.stations) assert.equal(st.medium, st.z < 0 ? 'water' : 'air');
  // unit hub displacement in mode 1 needs an equivalent tip force k/phi; the mudline moment is of the order force x lever
  const md = s.modes[0];
  const lever = md.momentMudline / md.shearMudline;
  assert.ok(lever > 0.75 * (90 - -20) && lever < 1.05 * (90 - -20), `effective lever ${lever}`);
  assert.ok(md.momentBase < md.momentMudline && md.momentBase > 0);
  // generalised stiffness = omega^2 x generalised mass
  rel(md.stiffness, (2 * Math.PI * md.frequency) ** 2 * md.mass, 1e-12);
});

test('rainflowCount: ASTM E1049-85 worked example (ranges 3, 4, 6, 8, 9 with counts 0.5, 1.5, 0.5, 1.0, 0.5)', () => {
  const c = rainflowCount([-2, 1, -3, 5, -1, 3, -4, 4, -2]);
  const tally = new Map();
  c.ranges.forEach((r, i) => tally.set(r, (tally.get(r) || 0) + c.counts[i]));
  assert.deepEqual([...tally.entries()].sort((a, b) => a[0] - b[0]), [[3, 0.5], [4, 1.5], [6, 0.5], [8, 1.0], [9, 0.5]]);
});

test('damageEquivalentLoad: constant-amplitude signal gives DEL = range (n / nEq)^(1/m)', () => {
  const x = [];
  for (let i = 0; i <= 4000; i++) x.push(7 * Math.sin(2 * Math.PI * i / 40));      // 100 cycles of range 14
  const d = damageEquivalentLoad(x, { m: 4, nEq: 25 });
  close(d.nCycles, 100, 1.01);
  rel(d.del, 14 * Math.pow(100 / 25, 0.25), 2e-3);
  const hook = damageEquivalentLoad(x, { m: 3, nEq: 10, rainflow: () => ({ ranges: [2], counts: [80] }) });
  rel(hook.del, 2 * 2, 1e-12);
});

test('seriesStats, logDecrement and identifyPeaks on analytic signals', () => {
  const st = seriesStats([1, 2, 3, 4, 5]);
  assert.equal(st.mean, 3); assert.equal(st.min, 1); assert.equal(st.max, 5); close(st.std, Math.sqrt(2), 1e-12);
  const z = 0.03; const fd = 0.5; const dt = 0.01; const x = [];
  const wn = 2 * Math.PI * fd / Math.sqrt(1 - z * z);
  for (let i = 0; i < 4000; i++) x.push(2 + Math.exp(-z * wn * i * dt) * Math.cos(2 * Math.PI * fd * i * dt + 0.3));
  const ld = logDecrement(x, dt);
  rel(ld.zeta, z, 0.02); rel(ld.frequency, fd, 2e-3);
  const f = []; const p = [];
  for (let i = 0; i < 200; i++) { f.push(i * 0.01); p.push(1 + (Math.abs(i * 0.01 - 0.31) < 0.006 ? 500 : 0)); }
  const pk = identifyPeaks(f, p, [{ name: 'mode', f: 0.3 }, { name: 'other', f: 1.0 }]);
  close(pk[0].fPeak, 0.31, 1e-9); assert.ok(pk[0].prominence > 100); assert.ok(pk[1].prominence < 2);
});

test('free decay of a supplied mode reproduces its frequency and damping ratio (no aerodynamics)', () => {
  const str = oneMode(0.02);
  const r = simulateCoupled({ surface, turbine, controller: null, structure: str, wind: 0, dt: 0.01, duration: 40,
    initial: { omega: 0, pitchDeg: 90, q: [0.2] }, aeroOverride: { ct: 0 } });
  const ld = logDecrement(r.modal[0], 0.01);
  rel(ld.zeta, 0.02, 0.02);
  rel(ld.frequency, 0.28 * Math.sqrt(1 - 0.0004), 1e-3);
  close(r.towerTopDisp[0], 0.2, 1e-12);
});

test('static thrust gives the static deflection F/k and section moments F x lever (mode-acceleration recovery)', () => {
  const str = oneMode(0.6);
  const ct = 0.8; const V = 10;
  const T = 0.5 * 1.225 * A * ct * V * V;
  const r = simulateCoupled({ surface, turbine, controller: null, structure: str, wind: V, dt: 0.02, duration: 60,
    initial: { omega: 0, pitchDeg: 0, q: [T / str.modes[0].stiffness] }, aeroOverride: { ct } });
  const n = r.t.length - 1;
  rel(r.thrust[n], T, 1e-9);
  rel(r.towerTopDisp[n], T / str.modes[0].stiffness, 1e-6);
  rel(r.mudlineMoment[n], T * 120, 1e-6);
  rel(r.towerBaseMoment[n], T * 90, 1e-6);
  rel(r.mudlineShear[n], T, 1e-6);
  // started from rest instead, the response settles to the same values
  const r0 = simulateCoupled({ surface, turbine, controller: null, structure: str, wind: V, dt: 0.02, duration: 60, initial: { omega: 0, pitchDeg: 0 }, aeroOverride: { ct } });
  rel(r0.towerTopDisp[n], T / str.modes[0].stiffness, 1e-3);
  rel(r0.mudlineMoment[n], T * 120, 1e-3);
});

test('multi-mode beam model: a static tip force is carried to the mudline exactly (statics) with only two modes retained', () => {
  const seg = monopileTowerSegments({ tower: { ...rotor.tower, height: 77.6 }, waterDepth: 20, transitionZ: 10, monopileD: 6, monopileT: 0.06 });
  const str = buildSupportStructure({ segments: seg, rnaMass: 350000, zHub: 90, nModes: 2, nElem: 24, zeta: 0.5 });
  const ct = 0.7; const V = 11;
  const T = 0.5 * 1.225 * A * ct * V * V;
  const r = simulateCoupled({ surface, turbine: { ...turbine, hubHeight: 90 }, controller: null, structure: str, wind: V, dt: 0.01, duration: 60,
    initial: { omega: 0, pitchDeg: 0 }, aeroOverride: { ct } });
  const n = r.t.length - 1;
  rel(r.mudlineMoment[n], T * 110, 1e-4);
  rel(r.towerBaseMoment[n], T * 80, 1e-4);
  rel(r.mudlineShear[n], T, 1e-4);
  // modal truncation: hub deflection from two modes is slightly below sum phi^2/k of all modes, above the first-mode value
  const d1 = T * str.modes[0].phiHub ** 2 / str.modes[0].stiffness;
  const d2 = d1 + T * str.modes[1].phiHub ** 2 / str.modes[1].stiffness;
  rel(r.towerTopDisp[n], d2, 1e-3);
  assert.ok(d2 > d1);
});

test('drivetrain energy balance at steady state below rated: aerodynamic torque = generator torque x gear ratio', () => {
  const r = simulateCoupled({ surface, turbine, controller: cfg, structure: oneMode(0.01), wind: 8, dt: 0.02, duration: 300,
    initial: { omega: 0.9, pitchDeg: 0 }, tTransient: 250 });
  const n = r.t.length - 1;
  rel(r.aeroTorque[n], 97 * r.genTorque[n], 2e-3);
  rel(r.power[n], 0.944 * r.genTorque[n] * r.genSpeed[n], 1e-12);
  rel(r.genSpeed[n], 97 * r.rotorSpeed[n], 1e-12);
  rel(r.rotorSpeed[n] * surface.R / 8, cfg.design.tsrOpt, 0.02);
  // thrust and deflection consistent with the surface at the settled point
  const c = surfaceLookup(surface, r.rotorSpeed[n] * surface.R / 8, r.pitchDeg[n]);
  rel(r.thrust[n], c.ct * 0.5 * 1.225 * A * 64, 2e-3);
  assert.ok(r.stats.towerTopDisp.std < 1e-3 * r.stats.towerTopDisp.mean, 'steady');
  assert.equal(r.mode[n], 'production');
});

test('closed loop above rated holds rated speed and power through a wind step, with the tower in the loop', () => {
  const r = simulateCoupled({ surface, turbine, controller: cfg, structure: oneMode(0.01), wind: (t) => (t < 60 ? 14 : 17), dt: 0.02, duration: 200,
    initial: { omega: rotor.ratedRotorSpeed, pitchDeg: 8.3 }, tTransient: 170 });
  rel(r.stats.rotorSpeed.mean, rotor.ratedRotorSpeed, 3e-3);
  rel(r.stats.power.mean, 5e6, 5e-3);
  assert.ok(r.stats.rotorSpeed.std < 0.005 * rotor.ratedRotorSpeed);
  let mx = 0; for (let i = 3000; i < r.t.length; i++) mx = Math.max(mx, r.rotorSpeed[i]);
  assert.ok(mx < 1.12 * rotor.ratedRotorSpeed && mx > 1.003 * rotor.ratedRotorSpeed);
  // thrust falls with wind speed above rated
  assert.ok(r.stats.thrust.mean < r.thrust[2900]);
});

test('aerodynamic damping emerges from the relative wind: free-decay estimate agrees with phi^2 dT/dV / (2 sqrt(k m))', () => {
  const d = aeroDampingRatio({ surface, turbine, controller: cfg, structure: oneMode(0.01), wind: 8, dt: 0.02, duration: 1, initial: { omega: 0.95, pitchDeg: 0 } },
    { settle: 150 });
  assert.ok(d.zetaDecay > 0.02 && d.zetaDecay < 0.12, `aerodynamic damping ratio ${d.zetaDecay}`);
  rel(d.zetaDecay, d.zetaAnalytic, 0.12);
  rel(d.frequency, 0.28, 0.01);
  assert.ok(d.dTdV > 0);
  // parked/feathered rotor with zero thrust slope: only structural damping remains
  const p = aeroDampingRatio({ surface, turbine, controller: null, structure: oneMode(0.01), wind: 8, dt: 0.02, duration: 1, initial: { omega: 0, pitchDeg: 90 }, aeroOverride: { ct: 0 } },
    { settle: 5 });
  close(p.zetaDecay, 0, 5e-4);
});

test('reduced-order dynamic inflow: same steady state, thrust overshoot on a sudden wind increase', () => {
  const base = { surface, turbine, controller: cfg, structure: oneMode(0.05), wind: (t) => (t < 100 ? 13 : 15), dt: 0.02, duration: 220, initial: { omega: rotor.ratedRotorSpeed, pitchDeg: 6.2 }, tTransient: 200 };
  const a = simulateCoupled(base);
  const b = simulateCoupled({ ...base, dynamicInflow: true });
  rel(b.stats.thrust.mean, a.stats.thrust.mean, 5e-3);
  rel(b.stats.rotorSpeed.mean, a.stats.rotorSpeed.mean, 2e-3);
  // identical before the step; at the step the wake is frozen, so the blades see more of the gust and the thrust overshoots
  rel(b.thrust[4900], a.thrust[4900], 1e-6);
  assert.ok(b.thrust[5000] > 1.005 * a.thrust[5000], `frozen-wake overshoot ${b.thrust[5000]} vs ${a.thrust[5000]}`);
  assert.ok(b.thrust[5000] < 1.2 * a.thrust[5000]);
});

test('wave loading on the monopile: Morison inertia load of a regular-like narrow sea matches the hand estimate; PSD shows the wave peak', () => {
  const seg = [{ zBase: -20, zTop: 10, dBase: 6, dTop: 6, tBase: 0.06, tTop: 0.06 }, { zBase: 10, zTop: 90, dBase: 6, dTop: 4, tBase: 0.03, tTop: 0.02 }];
  const str = buildSupportStructure({ segments: seg, rnaMass: 350000, zHub: 90, nModes: 2, nElem: 30, zeta: 0.02, cd: 0, cm: 2 });
  const r = simulateCoupled({ surface, turbine: { ...turbine, hubHeight: 90 }, controller: null, structure: str, wind: 0, dt: 0.025, duration: 600, tTransient: 60,
    initial: { omega: 0, pitchDeg: 90 }, aeroOverride: { ct: 0 }, waves: { Hs: 1.0, Tp: 12, gamma: 7, depth: 20, seed: 3 } });
  rel(r.stats.waveElevation.std, 0.25, 0.08);
  // quasi-static (wave period 12 s >> structural period): mudline shear ~ sum rho CM A dz * acceleration; compare standard deviations
  // using Airy theory at the peak period: a(z) = (Hs/4 std) * omega^2 cosh(k(z+d))/sinh(kd), integrated over depth = omega^2 /k * std(eta)
  const w = 2 * Math.PI / 12; const k = waveNumber(w, 20);
  const shearStd = 1025 * 2 * Math.PI * 9 * (w * w / k) * r.stats.waveElevation.std;
  rel(r.stats.mudlineShear.std, shearStd, 0.2);
  assert.ok(r.stats.mudlineMoment.std > 5 * r.stats.mudlineShear.std, 'moment lever of several metres');
  const wavePk = r.psd.peaks.find((p) => p.name === 'wave peak');
  assert.ok(wavePk.prominence > 30);
  assert.ok(r.del.mudlineMoment > 0 && r.del.towerBaseMoment > 0 && r.del.mudlineMoment > r.del.towerBaseMoment);
  // MacCamy-Fuchs on short waves (Tp = 3.5 s, ka = 0.99 at the peak): loads scale by about |CM(ka)|/2 (drag-free, linear)
  const short = { surface, turbine: { ...turbine, hubHeight: 90 }, controller: null, structure: str, wind: 0, dt: 0.025, duration: 400, tTransient: 40,
    initial: { omega: 0, pitchDeg: 90 }, aeroOverride: { ct: 0 } };
  const m0 = simulateCoupled({ ...short, waves: { Hs: 0.5, Tp: 3.5, gamma: 5, depth: 20, seed: 3 } });
  const m1 = simulateCoupled({ ...short, waves: { Hs: 0.5, Tp: 3.5, gamma: 5, depth: 20, seed: 3, macCamyFuchs: true } });
  const kp = waveNumber(2 * Math.PI / 3.5, 20);
  const ratio = macCamyFuchsCm(3 * kp).magnitude / 2;
  assert.ok(ratio < 0.95);
  rel(m1.stats.mudlineShear.std / m0.stats.mudlineShear.std, ratio, 0.15);
});

test('floater: surge free decay reproduces 2 pi sqrt((M + A)/K); mean offset and trim under thrust equal F/K and F h/K55', () => {
  const fl = { mass: [8e6, 8e6, 6e9], addedMass: [2e6, 4e6, 2e9], stiffness: [8e4, 3.5e5, 1.5e9], damping: [2e4, 5e4, 5e7], refZ: 0, towerBaseZ: 10, rnaMass: 350000 };
  const r = simulateCoupled({ surface, turbine, controller: null, floater: fl, wind: 0, dt: 0.1, duration: 2000, initial: { omega: 0, pitchDeg: 90, platform: [5, 0, 0] }, aeroOverride: { ct: 0 } });
  const ld = logDecrement(r.platform.surge, 0.1, { maxPeaks: 5 });
  const wn = Math.sqrt(8e4 / 1e7);
  rel(1 / ld.frequency, 2 * Math.PI / (wn * Math.sqrt(1 - (2e4 / (2 * 1e7 * wn)) ** 2)), 2e-3);
  rel(ld.zeta, 2e4 / (2 * 1e7 * wn), 0.02);
  // static thrust
  const ct = 0.8; const V = 10; const T = 0.5 * 1.225 * A * ct * V * V;
  const s = simulateCoupled({ surface, turbine, controller: null, floater: { ...fl, damping: [2e6, 3e6, 3e9] }, wind: V, dt: 0.1, duration: 1500,
    initial: { omega: 0, pitchDeg: 0, platform: [T / 8e4, 0, T * 100 / 1.5e9] }, aeroOverride: { ct } });
  const n = s.t.length - 1;
  rel(s.platform.surge[n], T / 8e4, 1e-6);
  rel(s.platform.pitch[n], T * 100 / 1.5e9, 1e-6);
  close(s.platform.heave[n], 0, 1e-9);
  // tower-base moment: thrust x tower height plus the weight offset of the RNA at the trim angle
  rel(s.towerBaseMoment[n], T * 90 + 350000 * 9.80665 * 90 * Math.sin(T * 100 / 1.5e9), 1e-6);
  // first-order wave forces through the RAO: surge response variance grows from zero
  const wv = simulateCoupled({ surface, turbine, controller: null, floater: { ...fl, forceRao: { f: [0.02, 0.5], surge: { amp: [1e6, 1e6], phaseDeg: [90, 90] } } }, wind: 0, dt: 0.1, duration: 600,
    initial: { omega: 0, pitchDeg: 90 }, aeroOverride: { ct: 0 }, waves: { Hs: 2, Tp: 10, depth: 200, seed: 1 } });
  assert.ok(wv.stats['platform.surge'].std > 0.01);
});

test('input validation: time step too large for the structural mode, missing surface', () => {
  const str = oneMode();
  assert.throws(() => simulateCoupled({ surface, turbine, controller: cfg, structure: str, wind: 8, dt: 2, duration: 10 }), /dt/);
  assert.throws(() => simulateCoupled({ turbine, controller: cfg, wind: 8, dt: 0.02, duration: 10 }), RangeError);
  assert.throws(() => buildSupportStructure({ segments: [], rnaMass: 1, zHub: 1 }), RangeError);
});

test('closed loop above rated in steady wind: rated speed is held without a limit cycle at 12, 13 and 18 m/s (with dynamic inflow)', () => {
  // The run starts 3 % over speed at fine pitch and gets a 1 m/s, 2 s gust at t = 50 s; the last 100 s must be quiet.
  for (const V of [12, 13, 18]) {
    const r = simulateCoupled({ surface, turbine, controller: cfg, structure: oneMode(0.02), dt: 0.02, duration: 300, dynamicInflow: true,
      wind: (t) => V + (t > 50 && t < 52 ? 1 : 0), initial: { omega: rotor.ratedRotorSpeed * 1.03, pitchDeg: 0 } });
    const i0 = Math.round(200 / 0.02);
    const om = seriesStats(r.rotorSpeed, i0); const pt = seriesStats(r.pitchDeg, i0);
    rel(om.mean, rotor.ratedRotorSpeed, 1e-4, `V = ${V}: mean speed`);
    assert.ok((om.max - om.min) / rotor.ratedRotorSpeed < 1e-4, `V = ${V}: speed ripple ${(om.max - om.min) / rotor.ratedRotorSpeed}`);
    assert.ok(pt.max - pt.min < 0.01, `V = ${V}: pitch peak-to-peak ${pt.max - pt.min} deg`);
    assert.ok(pt.mean > 2 && pt.mean < 20);
  }
  // the gains of the plain fit (designed at the vanishing fine-pitch sensitivity) do limit-cycle at 12 m/s
  const stiff = controllerFromRotor({ rotor, surface, pitchDesign: 'fit' });
  const r = simulateCoupled({ surface, turbine, controller: stiff, structure: oneMode(0.02), dt: 0.02, duration: 300, dynamicInflow: true,
    wind: (t) => 12 + (t > 50 && t < 52 ? 1 : 0), initial: { omega: rotor.ratedRotorSpeed * 1.03, pitchDeg: 0 } });
  const pt = seriesStats(r.pitchDeg, Math.round(200 / 0.02));
  assert.ok(pt.max - pt.min > 1, `plain-fit gains: pitch peak-to-peak ${pt.max - pt.min} deg`);
});

test('monopileTowerSegments honours the hub height: the tower top sits at hubHeight − hubOffset and the structure builds', () => {
  // the reference tower record is 87.6 m long (land based); on a transition at +10 m its top would be at 97.6 m, above a 90 m hub
  const base = { tower: rotor.tower, waterDepth: 20, transitionZ: 10, monopileD: 6, monopileT: 0.06 };
  const legacy = monopileTowerSegments(base);
  close(legacy[1].zTop, 10 + rotor.tower.height, 1e-12);
  assert.throws(() => buildSupportStructure({ segments: legacy, rnaMass: 350000, zHub: 90, nModes: 1, nElem: 20 }), /zHub/);
  const seg = monopileTowerSegments({ ...base, hubHeight: 90, hubOffset: 2.4 });
  close(seg[1].zTop, 87.6, 1e-12); close(seg[1].zBase, 10, 1e-12); close(seg[0].zTop, 10, 1e-12); close(seg[0].zBase, -20, 1e-12);
  assert.equal(seg[1].dBase, rotor.tower.baseDiameter); assert.equal(seg[1].dTop, rotor.tower.topDiameter);
  const s = buildSupportStructure({ segments: seg, rnaMass: 350000, zHub: 90, nModes: 2, nElem: 30 });
  close(s.zTop, 87.6, 1e-12); assert.ok(s.modes[0].frequency > 0.2 && s.modes[0].frequency < 0.4);
  // identical to shortening the tower record by hand
  const byHand = buildSupportStructure({ segments: monopileTowerSegments({ ...base, tower: { ...rotor.tower, height: 77.6 } }), rnaMass: 350000, zHub: 90, nModes: 2, nElem: 30 });
  rel(s.modes[0].frequency, byHand.modes[0].frequency, 1e-12);
  // a taller hub lengthens the tower and lowers the frequency; default offset 0 puts the top at the hub
  close(monopileTowerSegments({ ...base, hubHeight: 110 })[1].zTop, 110, 1e-12);
  const tall = buildSupportStructure({ segments: monopileTowerSegments({ ...base, hubHeight: 110, hubOffset: 2.4 }), rnaMass: 350000, zHub: 110, nModes: 1, nElem: 30 });
  assert.ok(tall.modes[0].frequency < s.modes[0].frequency);
  assert.throws(() => monopileTowerSegments({ ...base, hubHeight: 8 }), RangeError);
});

test('dynamic inflow: load statistics do not depend on the time step when the same wind series is used', () => {
  // 14 m/s with 15 % turbulence, one synthesis at 0.1 s interpolated linearly for every step size
  const ws = Float64Array.from({ length: 2400 }, (_, i) => 14 + 1.5 * Math.sin(0.11 * i * 0.1) + 1.0 * Math.sin(0.53 * i * 0.1 + 1) + 0.6 * Math.sin(1.7 * i * 0.1 + 2));
  const wind = (t) => { const x = Math.max(0, t) / 0.1; const i = Math.min(ws.length - 2, Math.floor(x)); return ws[i] + (ws[i + 1] - ws[i]) * (x - i); };
  const std = (dt, dyn) => {
    const r = simulateCoupled({ surface, turbine, controller: cfg, structure: oneMode(0.02), dt, wind, duration: 200, tTransient: 40, dynamicInflow: dyn });
    const i0 = Math.round(40 / dt);
    return { m: seriesStats(r.mudlineMoment, i0).std, t: seriesStats(r.thrust, i0).mean };
  };
  const a = std(0.05, true); const b = std(0.025, true);
  rel(a.m, b.m, 0.015, 'mudline-moment standard deviation'); rel(a.t, b.t, 2e-3, 'mean thrust');
  // dynamic inflow matters (it changes the thrust response to pitch and wind changes) while keeping the mean
  const q = std(0.05, false);
  assert.ok(Math.abs(a.m / q.m - 1) > 0.01, `with ${a.m}, without ${q.m}`); rel(a.t, q.t, 0.05);
});
