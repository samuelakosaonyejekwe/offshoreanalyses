import test from 'node:test';
import assert from 'node:assert/strict';
import {
  region2Gain, piPitchGains, gainCorrection, surfaceOptimum, operatingSchedule, pitchSensitivity, torqueLawConstants,
  generatorTorqueLaw, pitchActuatorStep, controllerFromRotor, createController, yawControllerStep,
  floatingDetuneCheck, simulateClosedLoop,
} from '../../src/engine/aero/control.js';
import { cpSurface, surfaceLookup, bemSolve } from '../../src/engine/aero/bem.js';
import { nrel5mw } from '../../src/engine/aero/rotors.js';

const D2R = Math.PI / 180;
const RPM = Math.PI / 30;
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} +- ${tol}`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${tol})`);

const rotor = nrel5mw();
const surface = cpSurface({ rotor, bem: { tiltDeg: 0, preconeDeg: 0 } });
const cfg = controllerFromRotor({ rotor, surface });

test('region2Gain: k = 0.5 rho pi R^5 Cp/lambda^3 by hand; close to the published NREL 5-MW constant 0.0255764 N m/rpm^2', () => {
  // 0.5*1.225*pi*63^5 = 1.90966e9; * 0.482 / 7.55^3 = 2.1388e6 N m s^2 (LSS); / 97^3 = 2.3434 (HSS); * (pi/30)^2 = 0.025698
  const g = region2Gain({ R: 63, cpMax: 0.482, tsrOpt: 7.55, rho: 1.225, gearboxRatio: 97 });
  rel(g.kLss, 0.5 * 1.225 * Math.PI * 63 ** 5 * 0.482 / 7.55 ** 3, 1e-14);
  rel(g.kLss, 2.1388e6, 1e-4);
  rel(g.kHss, g.kLss / 97 ** 3, 1e-14);
  rel(g.kHssRpm, 0.025698, 1e-4);
  rel(g.kHssRpm, 0.0255764, 0.01, 'published value');
  assert.throws(() => region2Gain({ R: 63, cpMax: 0.7, tsrOpt: 7 }), RangeError);
});

test('piPitchGains: published NREL 5-MW gains Kp = 0.01882681 s, Ki = 0.008068634 from the pole-placement formulas', () => {
  const g = piPitchGains({ inertia: 43784700, omegaRated: 12.1 * RPM, dPdTheta: -25.52e6, gearboxRatio: 97, zeta: 0.7, omegaN: 0.6 });
  // by hand: 2 * 43784700 * 1.26711 * 0.7 * 0.6 / (97 * 25.52e6)
  rel(g.kp, 2 * 43784700 * (12.1 * RPM) * 0.7 * 0.6 / (97 * 25.52e6), 1e-14);
  rel(g.kp, 0.01882681, 1e-3);
  rel(g.ki, 0.008068634, 1e-3);
  // Ki/Kp = omega_n / (2 zeta)
  rel(g.ki / g.kp, 0.6 / 1.4, 1e-12);
  assert.throws(() => piPitchGains({ inertia: 1, omegaRated: 1, dPdTheta: 5 }), RangeError);
});

test('gainCorrection: 1 at fine pitch, 1/2 at thetaK (NREL thetaK = 6.302336 deg)', () => {
  const tk = 6.302336 * D2R;
  assert.equal(gainCorrection(0, tk), 1);
  close(gainCorrection(tk, tk), 0.5, 1e-15);
  close(gainCorrection(3 * tk, tk), 0.25, 1e-15);
});

test('five-region torque law: hand values with the NREL 5-MW constants and continuity at the region boundaries', () => {
  const k = 0.0255764 / (RPM * RPM);
  const tc = torqueLawConstants({ k, ratedGenSpeed: 1173.7 * RPM, ratedPowerMech: 5296610, minGenSpeed: 670 * RPM, region15Span: 0.3, slipPercent: 10 });
  rel(tc.rtGnSp, 121.6805, 1e-6);                          // 0.99 x 1173.7 rpm
  rel(tc.rgn2Sp, 871 * RPM, 1e-12);
  rel(tc.sySp, tc.rtGnSp / 1.1, 1e-12);
  // region 2 at 1000 rpm: 0.0255764 * 1000^2 = 25,576.4 N m
  const r2 = generatorTorqueLaw(tc, 1000 * RPM);
  assert.equal(r2.region, '2'); rel(r2.torque, 25576.4, 1e-9);
  assert.equal(generatorTorqueLaw(tc, 600 * RPM).torque, 0);
  assert.equal(generatorTorqueLaw(tc, 800 * RPM).region, '1.5');
  const eps = 1e-7;
  for (const w of [tc.ctInSp, tc.rgn2Sp, tc.trGnSp, tc.rtGnSp]) {
    close(generatorTorqueLaw(tc, w - eps).torque, generatorTorqueLaw(tc, w + eps).torque, 0.05, `continuity at ${w}`);
  }
  // rated: torque x speed = rated mechanical power; constant power above
  rel(generatorTorqueLaw(tc, tc.rtGnSp).torque * tc.rtGnSp, 5296610, 1e-12);
  rel(generatorTorqueLaw(tc, 130).torque * 130, 5296610, 1e-12);
  rel(generatorTorqueLaw(tc, 1173.7 * RPM).torque, 43093.55, 1e-5);   // published rated generator torque
  // region 3 is forced once the blades have pitched
  assert.equal(generatorTorqueLaw(tc, 1000 * RPM, { pitch: 5 * D2R }).region, '3');
  assert.equal(generatorTorqueLaw(tc, 130, { constantPower: false }).torque, tc.rtTq);
});

test('pitch actuator: second-order response with 4.6 % overshoot (zeta 0.7) and rate limiting', () => {
  const st = { theta: 0, rate: 0 };
  let peak = 0;
  const cmd = 0.5 * D2R;
  for (let i = 0; i < 400; i++) { pitchActuatorStep(st, { cmd, dt: 0.005 }); peak = Math.max(peak, st.theta); }
  close(st.theta, cmd, 1e-4 * cmd);
  close(peak / cmd - 1, Math.exp(-Math.PI * 0.7 / Math.sqrt(1 - 0.49)), 4e-3);
  // large step: the rate limit of 8 deg/s governs -> 20 deg takes at least 2.5 s
  const s2 = { theta: 0, rate: 0 };
  let maxRate = 0; let t20 = null;
  for (let i = 1; i <= 2000; i++) {
    pitchActuatorStep(s2, { cmd: 20 * D2R, dt: 0.005 });
    maxRate = Math.max(maxRate, Math.abs(s2.rate));
    if (t20 === null && s2.theta >= 0.99 * 20 * D2R) t20 = i * 0.005;
  }
  close(maxRate, 8 * D2R, 1e-12);
  assert.ok(t20 > 2.45 && t20 < 3.2, `time to 20 deg ${t20}`);
  // range limit
  const s3 = { theta: 0, rate: 0 };
  for (let i = 0; i < 4000; i++) pitchActuatorStep(s3, { cmd: 2, dt: 0.01 });
  close(s3.theta, Math.PI / 2, 1e-12);
});

test('controllerFromRotor: gains follow the formulas from the rotor\'s own surfaces', () => {
  const d = cfg.design;
  rel(cfg.k, 0.5 * 1.225 * Math.PI * surface.R ** 5 * d.cpMax / d.tsrOpt ** 3 / 97 ** 3, 1e-12);
  rel(cfg.ratedTorque, 5e6 / 0.944 / (12.1 * RPM * 97), 1e-12);
  rel(cfg.ratedTorque, 43093.55, 1e-5);
  assert.ok(d.dPdTheta0 < 0 && cfg.thetaK > 0);
  rel(cfg.kp, 2 * rotor.drivetrainInertia * rotor.ratedRotorSpeed * 0.7 * 0.6 / (97 * -d.dPdTheta0), 1e-12);
  rel(cfg.ki, rotor.drivetrainInertia * rotor.ratedRotorSpeed * 0.36 / (97 * -d.dPdTheta0), 1e-12);
  assert.ok(d.ratedWind > 11 && d.ratedWind < 12);
});

test('pitchSensitivity: dP/dtheta is negative, grows in magnitude with pitch, and matches a direct BEM finite difference', () => {
  const s = pitchSensitivity({ surface, ratedPowerAero: 5e6 / 0.944, ratedRotorSpeed: rotor.ratedRotorSpeed, speeds: [13, 15, 18, 22, 25] });
  for (let i = 0; i < s.U.length; i++) {
    assert.ok(s.dPdTheta[i] < 0);
    if (i > 0) { assert.ok(s.dPdTheta[i] < s.dPdTheta[i - 1]); assert.ok(s.pitchDeg[i] > s.pitchDeg[i - 1]); }
  }
  // independent check at 18 m/s against the BEM itself
  const i = 2; const p = s.pitchDeg[i];
  const P = (pd) => bemSolve({ rotor, U: 18, omega: rotor.ratedRotorSpeed, pitchDeg: pd, tiltDeg: 0, preconeDeg: 0 }).power;
  rel(P(p), 5e6 / 0.944, 0.01, 'operating point holds rated power');
  rel(s.dPdTheta[i], (P(p + 1) - P(p - 1)) / (2 * D2R), 0.03);
  // the linear fit reproduces the points (thetaK = pitch at which the sensitivity has doubled)
  for (let j = 0; j < s.U.length; j++) rel(s.dPdTheta0 * (1 + s.pitchDeg[j] * D2R / s.thetaK), s.dPdTheta[j], 0.12);
});

test('operatingSchedule: regions, rated power held above rated, optimum tip-speed ratio in region 2, peak shaving', () => {
  const args = { surface, ratedPower: 5e6, efficiency: 0.944, ratedRotorSpeed: rotor.ratedRotorSpeed, minRotorSpeed: rotor.minRotorSpeed, gearboxRatio: 97,
    speeds: [2, 4, 8, 11, 13, 18, 25, 27] };
  const s = operatingSchedule(args);
  assert.deepEqual(s.region, ['off', '1.5', '2', '2.5', '3', '3', '3', 'off']);
  rel(s.tsr[2], s.tsrOpt, 1e-9);
  rel(s.cp[2], s.cpMax, 1e-9);
  for (const i of [4, 5, 6]) rel(s.power[i], 5e6, 1e-6);
  rel(s.omega[1], rotor.minRotorSpeed, 1e-12);
  // torque balance at each point: N * Tgen * omega = aerodynamic power
  for (const i of [1, 2, 3, 4, 5]) rel(s.genTorque[i] * 97 * s.omega[i], s.aeroPower[i], 1e-9);
  rel(s.genTorque[5], 43093.55, 1e-4);
  // surface optimum agrees with a direct lookup
  const o = surfaceOptimum(surface);
  assert.ok(o.cpMax >= surfaceLookup(surface, o.tsrOpt, o.pitchOptDeg).cp - 1e-9);
  // peak shaving: thrust clipped, power reduced only where clipped
  const fine = operatingSchedule({ ...args, speeds: [8, 10, 11, 11.4, 12, 14] });
  const limit = 0.85 * Math.max(...fine.thrust);
  const ps = operatingSchedule({ ...args, speeds: [8, 10, 11, 11.4, 12, 14], maxThrust: limit });
  for (let i = 0; i < ps.U.length; i++) {
    assert.ok(ps.thrust[i] <= limit * (1 + 1e-6));
    if (fine.thrust[i] > limit) { assert.ok(ps.power[i] < fine.power[i]); assert.ok(ps.pitchDeg[i] > fine.pitchDeg[i]); assert.match(ps.region[i], /peak-shaved/); } else rel(ps.power[i], fine.power[i], 1e-12);
  }
});

test('closed loop, region 2: steady state tracks the optimum tip-speed ratio and aerodynamic torque equals generator torque', () => {
  const r = simulateClosedLoop({ surface, cfg, inertia: rotor.drivetrainInertia, wind: 8, dt: 0.02, duration: 300, omega0: 0.8, pitch0Deg: cfg.pitchMin / D2R });
  const n = r.t.length - 1;
  rel(r.aeroTorque[n], 97 * r.genTorque[n], 2e-3, 'drivetrain equilibrium');
  rel(r.rotorSpeed[n] * surface.R / 8, cfg.design.tsrOpt, 0.02);
  rel(r.power[n], 0.944 * r.genTorque[n] * r.rotorSpeed[n] * 97, 1e-12);
  // electrical power = eta * Cp,max * 0.5 rho A V^3
  rel(r.power[n], 0.944 * cfg.design.cpMax * 0.5 * 1.225 * Math.PI * surface.R ** 2 * 512, 0.01);
  close(r.pitchDeg[n], cfg.pitchMin / D2R, 1e-6);
});

test('closed loop, region 3: a step in wind above rated is rejected and rated speed and power are held', () => {
  const r = simulateClosedLoop({ surface, cfg, inertia: rotor.drivetrainInertia, wind: (t) => (t < 60 ? 14 : 17), dt: 0.02, duration: 180,
    omega0: rotor.ratedRotorSpeed, pitch0Deg: 8.3 });
  const n = r.t.length - 1;
  let maxDev = 0;
  for (let i = 3000; i <= n; i++) maxDev = Math.max(maxDev, Math.abs(r.rotorSpeed[i] / rotor.ratedRotorSpeed - 1));
  assert.ok(maxDev > 0.005 && maxDev < 0.12, `speed excursion ${maxDev}`);
  rel(r.rotorSpeed[n], rotor.ratedRotorSpeed, 2e-3);
  rel(r.power[n], 5e6, 5e-3);
  // steady pitch equals the schedule value at 17 m/s
  const s = operatingSchedule({ surface, ratedPower: 5e6, efficiency: 0.944, ratedRotorSpeed: rotor.ratedRotorSpeed, gearboxRatio: 97, pitchFineDeg: cfg.pitchMin / D2R, speeds: [14, 17] });
  close(r.pitchDeg[n], s.pitchDeg[1], 0.15);
  close(r.pitchDeg[2999], s.pitchDeg[0], 0.15);
  rel(r.aeroTorque[n], 97 * r.genTorque[n], 5e-3);
});

test('supervisory sequences: emergency stop on grid loss feathers at the pitch rate and removes torque; start-up reaches production', () => {
  const stop = simulateClosedLoop({ surface, cfg, inertia: rotor.drivetrainInertia, wind: 12, dt: 0.02, duration: 60, omega0: rotor.ratedRotorSpeed,
    pitch0Deg: 4, events: (t) => (t >= 10 ? { gridLoss: true } : null) });
  const n = stop.t.length - 1;
  assert.equal(stop.mode[400], 'production');
  assert.equal(stop.mode[n], 'estop');
  assert.equal(stop.genTorque[510], 0);
  close(stop.pitchDeg[n], 90, 0.5);
  assert.ok(stop.rotorSpeed[n] < 0.2 * rotor.ratedRotorSpeed, 'rotor decelerated');
  // pitch travel between 12 s and 17 s is limited by 8 deg/s
  close((stop.pitchDeg[850] - stop.pitchDeg[600]) / 5, 8, 0.1);
  // over-speed stays bounded during the stop
  assert.ok(Math.max(...stop.rotorSpeed) < 1.25 * rotor.ratedRotorSpeed);

  const start = simulateClosedLoop({ surface, cfg, inertia: rotor.drivetrainInertia, wind: 10, dt: 0.02, duration: 400, omega0: 0, mode: 'parked',
    events: (t) => (t >= 5 && t < 5.1 ? { command: 'start' } : null) });
  const m = start.t.length - 1;
  assert.equal(start.mode[100], 'parked');
  assert.equal(start.mode[m], 'production');
  assert.ok(start.power[m] > 2.5e6 && start.power[m] < 4.2e6, `power at 10 m/s ${start.power[m]}`);
  close(start.pitchDeg[m], cfg.pitchMin / D2R, 0.05);
});

test('normal stop command returns the controller to parked', () => {
  const c = createController(cfg, { mode: 'production', pitch: 0.1, genSpeed: cfg.ratedGenSpeed });
  let w = cfg.ratedGenSpeed; let out;
  for (let i = 0; i < 6000; i++) {
    out = c.step({ dt: 0.02, genSpeed: w, pitch: out ? out.pitchCmd : 0.1, command: i === 10 ? 'stop' : undefined });
    if (i > 10) w = Math.max(0, w - 0.05);
  }
  assert.equal(out.mode, 'parked');
  close(out.pitchCmd, Math.PI / 2, 1e-9);
  assert.equal(out.genTorque, 0);
});

test('yaw controller: dead band, delay and constant yaw rate', () => {
  const st = {};
  let t = 0; let started = null; let yawed = 0; let err = 15;
  for (let i = 0; i < 4000; i++) {
    const r = yawControllerStep(st, { errorDeg: err, dt: 0.1, deadbandDeg: 8, delay: 30, rateDegPerS: 0.3 });
    if (r.yawing && started === null) started = t;
    yawed += r.yawRateDegPerS * 0.1; err = 15 - yawed;
    t += 0.1;
  }
  assert.ok(started > 30 && started < 45, `yawing started at ${started} s`);
  assert.ok(Math.abs(err) < 3, `final error ${err}`);
  const st2 = {};
  for (let i = 0; i < 4000; i++) assert.equal(yawControllerStep(st2, { errorDeg: 5, dt: 0.1 }).yawing, false);
});

test('floatingDetuneCheck: 0.6 rad/s fails and 0.2 rad/s passes for a 0.034 Hz platform-pitch mode', () => {
  assert.equal(floatingDetuneCheck({ omegaN: 0.6, platformPitchFreqHz: 0.034 }).ok, false);
  const r = floatingDetuneCheck({ omegaN: 0.2, platformPitchFreqHz: 0.034 });
  assert.equal(r.ok, true);
  close(r.omegaNMax, 2 * Math.PI * 0.034, 1e-12);
});

test('PI loop in region 3 behaves like the designed second-order system: damped speed oscillation near omega_n = 0.6 rad/s', () => {
  // small wind step at high wind where the gain schedule is active
  const r = simulateClosedLoop({ surface, cfg: { ...cfg, filterHz: 5 }, inertia: rotor.drivetrainInertia, wind: (t) => (t < 40 ? 18 : 18.5), dt: 0.01, duration: 80,
    omega0: rotor.ratedRotorSpeed, pitch0Deg: 14.6, actuator: { omegaAct: 30 } });
  // time of the peak speed error after the step ~ rise of a second-order system: between 1 and 6 s
  let iPk = 4000; for (let i = 4000; i < 6000; i++) if (r.rotorSpeed[i] > r.rotorSpeed[iPk]) iPk = i;
  const tPk = r.t[iPk] - 40;
  assert.ok(tPk > 0.8 && tPk < 6, `peak at ${tPk} s`);
  // first undershoot is much smaller than the peak (zeta ~ 0.7 -> ratio exp(-pi zeta/sqrt(1 - zeta^2)) = 0.046)
  let mn = Infinity; for (let i = iPk; i < 8000; i++) mn = Math.min(mn, r.rotorSpeed[i]);
  const pk = r.rotorSpeed[iPk] - rotor.ratedRotorSpeed; const un = rotor.ratedRotorSpeed - mn;
  assert.ok(un / pk < 0.35, `undershoot ratio ${un / pk}`);
  rel(r.rotorSpeed[8000], rotor.ratedRotorSpeed, 1e-3);
});
