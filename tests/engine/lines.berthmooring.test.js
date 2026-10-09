import test from 'node:test';
import assert from 'node:assert/strict';
import * as bm from '../../src/engine/lines/berthmooring.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const DEG = Math.PI / 180;
const vessel = { LBP: 240, beam: 42, draft: 14, AT: 1100, AL: 4200, type: 'tanker', side: 'starboard' };

test('wind and current forces from coefficient curves (hand calculation, OCIMF force definitions)', () => {
  // head wind 30 m/s: Fx = 0.5 * 1.28 * (-0.9) * 900 * 1100 = -570.24 kN
  const w0 = bm.environmentalForce({ kind: 'wind', speed: 30, thetaDeg: 0, vessel });
  assert.ok(rel(w0.Fx, -570240) < 1e-12 && Math.abs(w0.Fy) < 1e-6 && w0.generic === true);
  // wind on the port beam pushes the ship to starboard: Fy = -0.5 * 1.28 * 0.9 * 900 * 4200
  const w90 = bm.environmentalForce({ kind: 'wind', speed: 30, thetaDeg: 90, vessel });
  assert.ok(rel(w90.Fy, -0.5 * 1.28 * 0.9 * 900 * 4200) < 1e-12 && Math.abs(w90.Fx) < 1e-6);
  // user curves over 0..180 with port/starboard symmetry (Cx even, Cy and Cxy odd)
  const curves = { angleDeg: [0, 90, 180], Cx: [-1, 0, 1], Cy: [0, -2, 0], Cxy: [0, -0.1, 0] };
  const a = bm.environmentalForce({ kind: 'current', speed: 1.5, thetaDeg: 45, vessel, curves });
  const b = bm.environmentalForce({ kind: 'current', speed: 1.5, thetaDeg: -45, vessel, curves });
  const q = 0.5 * 1025 * 1.5 * 1.5 * 240 * 14;
  assert.ok(rel(a.Fx, -0.5 * q) < 1e-12 && rel(a.Fy, -1.0 * q) < 1e-12 && rel(a.Mz, -0.05 * q * 240) < 1e-12);
  assert.ok(rel(b.Fx, a.Fx) < 1e-12 && rel(b.Fy, -a.Fy) < 1e-12 && rel(b.Mz, -a.Mz) < 1e-12);
  assert.equal(a.generic, false);
  assert.throws(() => bm.environmentalForce({ kind: 'tide', speed: 1, thetaDeg: 0, vessel }), /kind/);
});

test('two breast lines, off-berth load: T = F/2 and sway = F/(2k) (hand calculation, tension-only)', () => {
  const k = 2e6, F = 6e5;
  const lines = [-60, 60].map((x) => ({ fairlead: [x, -21, 20], bollard: [x, -20, 6], MBL: 1.2e6, stiffness: k, material: 'steel' }));
  // ship origin 23 m off the berth face: starboard side at y = 2 m; fairlead height 20 m above keel, draft 14 m
  const r = bm.berthEquilibrium({ vessel, pose0: { x: 0, y: 23, headingDeg: 0 }, lines, extra: { sway: F } });
  assert.ok(r.converged);
  // geometry is 22 m horizontal, level (z = 0 - 14 + 20 = 6): tension = k * elongation along the line
  for (const ln of r.lines) assert.ok(rel(ln.tension, F / 2) < 1e-8);
  assert.ok(rel(r.excursion.sway, F / (2 * k)) < 1e-6);
  assert.ok(Math.abs(r.excursion.yawDeg) < 1e-8);
  assert.ok(rel(r.lines[0].utilisation, 0.25) < 1e-8 && r.lines[0].limitFraction === 0.55 && r.lines[0].ok);
  // bollard loads equal the line tensions, residual vanishes
  assert.ok(rel(r.bollards[0].load, F / 2) < 1e-8);
  assert.ok(Math.hypot(r.residual[0], r.residual[1]) < 1e-6 * F);
});

const linearFender = (x) => ({ position: [x, 0], height: 1.5, curve: { deflection: [0, 0.75], reaction: [0, 3e6] }, ratedReaction: 3e6 });

test('on-berth load: breast lines go slack, fenders carry it (compression-only), deflection = F/(2 kf)', () => {
  const lines = [-60, 60].map((x) => ({ fairlead: [x, -21, 20], bollard: [x, -20, 6], MBL: 1.2e6, stiffness: 2e6 }));
  const fenders = [linearFender(-40), linearFender(40)];
  const F = 8e5;
  const r = bm.berthEquilibrium({ vessel, pose0: { x: 0, y: 23, headingDeg: 0 }, lines, fenders, extra: { sway: -F } });
  for (const ln of r.lines) assert.equal(ln.tension, 0);
  // fender stiffness 4e6 N/m each: deflection = 8e5/(2*4e6) = 0.1 m; ship side ends 1.4 m from the berth face
  for (const f of r.fenders) assert.ok(rel(f.deflection, 0.1) < 1e-7 && rel(f.reaction, 4e5) < 1e-7);
  assert.ok(rel(r.pose.y, 21 + 1.4) < 1e-8);
  assert.ok(rel(r.checks.maxFenderUc, 4e5 / 3e6) < 1e-7);
  // off-centre load: yaw moment shared by the fender couple, R = F/2 -+ M/(80 m)
  // (two level centreline lines along the berth hold the ship in surge without loading it in sway)
  const M = 8e6;
  const surgeLines = [-200, 200].map((x) => ({ fairlead: [0, 0, 20], bollard: [x, 23, 6], MBL: 1.2e6, stiffness: 2e6 }));
  const r2 = bm.berthEquilibrium({ vessel, pose0: { x: 0, y: 23, headingDeg: 0 }, lines: [...lines, ...surgeLines], fenders, extra: { sway: -F, yaw: M } });
  assert.ok(Math.abs(r2.excursion.surge) < 1e-3);
  // positive yaw moment (bow to port, away from the berth) unloads the forward fender
  assert.ok(rel(r2.fenders[1].reaction, F / 2 - M / 80) < 1e-4 && rel(r2.fenders[0].reaction, F / 2 + M / 80) < 1e-4);
});

function tankerBerth() {
  // fairlead [x, y, z above keel]; bollards on the berth / dolphins (global)
  const mk = (name, type, fl, bl, extra = {}) => ({ name, type, fairlead: fl, bollard: bl, MBL: 1.3e6, EA: 1.1e8, inboardLength: 15, tail: { EA: 1.2e7, L: 11 }, pretension: 1.0e5, material: 'steel', ...extra });
  return {
    vessel,
    pose0: { x: 0, y: 22.9, headingDeg: 0 },
    lines: [
      mk('head 1', 'head', [118, -6, 24], [165, -35, 7]), mk('head 2', 'head', [118, -4, 24], [165, -35, 7], { bollardId: 'MD1' }),
      mk('fwd breast 1', 'breast', [95, -19, 23], [100, -40, 7]), mk('fwd breast 2', 'breast', [90, -20, 23], [100, -40, 7]),
      mk('fwd spring', 'spring', [70, -21, 22], [25, -3, 6]), mk('aft spring', 'spring', [-70, -21, 22], [-25, -3, 6]),
      mk('aft breast 1', 'breast', [-95, -19, 23], [-100, -40, 7]), mk('aft breast 2', 'breast', [-90, -20, 23], [-100, -40, 7]),
      mk('stern 1', 'stern', [-118, -6, 24], [-165, -35, 7]), mk('stern 2', 'stern', [-118, -4, 24], [-165, -35, 7]),
    ],
    fenders: [-60, -20, 20, 60].map((x) => ({ position: [x, 0], height: 2, curve: { deflection: [0, 0.3, 0.7, 1.0], reaction: [0, 1.5e6, 2.0e6, 2.4e6] } })),
    manifold: [5, 0], armEnvelope: { surge: 3, sway: 3 },
  };
}

test('full tanker berth under wind + current: sum of forces and moments = 0 with tension-only lines and fenders', () => {
  const model = tankerBerth();
  const wind = { speed: 60 * bm.KNOT, fromDeg: 250 }; // from the berth side, slightly aft: pushes the ship off the berth
  const current = { speed: 1.5, thetaDeg: 170 };
  const r = bm.berthEquilibrium({ ...model, wind, current });
  assert.ok(r.converged);
  // independent force balance from the reported results (global axes, moments about the displaced ship origin)
  const psi = r.pose.headingDeg * DEG, c = Math.cos(psi), s = Math.sin(psi);
  let Fx = 0, Fy = 0, Mz = 0;
  for (const e of r.environment) { Fx += c * e.Fx - s * e.Fy; Fy += s * e.Fx + c * e.Fy; Mz += e.Mz; }
  r.lines.forEach((ln, i) => {
    const fl = model.lines[i].fairlead;
    const rx = c * fl[0] - s * fl[1], ry = s * fl[0] + c * fl[1];
    const fx = -ln.forceOnBollard[0], fy = -ln.forceOnBollard[1];
    Fx += fx; Fy += fy; Mz += rx * fy - ry * fx;
    assert.ok(ln.tension >= 0);
    // horizontal component follows the vertical angle
    assert.ok(Math.abs(ln.horizontalTension - ln.tension * Math.cos(ln.verticalAngle)) < 1e-6 * (ln.tension + 1));
  });
  r.fenders.forEach((f, i) => {
    assert.ok(f.reaction >= 0);
    if (f.reaction > 0) {
      const px = model.fenders[i].position[0] - r.pose.x, py = 2 - f.deflection - r.pose.y;
      Fy += f.reaction; Mz += px * f.reaction; void py;
    }
  });
  const scale = Math.hypot(...r.environment.map((e) => Math.hypot(e.Fx, e.Fy)));
  assert.ok(Math.abs(Fx) < 1e-6 * scale && Math.abs(Fy) < 1e-6 * scale && Math.abs(Mz) < 1e-6 * scale * 240, `${Fx} ${Fy} ${Mz}`);
  // off-berth wind: breast lines loaded above pretension, ship moved away from the berth
  assert.ok(r.excursion.sway > 0);
  assert.ok(r.lines[2].tension > 1e5);
  assert.ok(r.checks.maxLineUc > 0 && r.checks.armOk !== null);
  // bollard shared by the two forward breast lines carries (nearly) the sum of their tensions
  const shared = r.bollards.find((b) => b.id === '100,-40,7');
  assert.ok(shared.load <= r.lines[2].tension + r.lines[3].tension + 1e-6 && shared.load > 0.98 * (r.lines[2].tension + r.lines[3].tension));
});

test('line stiffness with tail, pretension set-up, tide effect on tension and vertical angle', () => {
  const model = tankerBerth();
  const r0 = bm.berthEquilibrium(model);
  // series stiffness of the first head line: main (geometric + 15 m inboard - 11 m tail) / EA + tail
  const l = Math.hypot(165 - 118, -35 - (22.9 - 6), 7 - (0 - 14 + 24));
  assert.ok(rel(r0.lines[0].stiffness, 1 / ((l + 15 - 11) / 1.1e8 + 11 / 1.2e7)) < 1e-12);
  // the pretensioned lines hold the ship against the fenders: total fender reaction = total line pull
  assert.ok(r0.fenders.some((f) => f.reaction > 0));
  const tot = r0.fenders.reduce((a, f) => a + f.reaction, 0);
  const pull = r0.lines.reduce((a, ln) => a + ln.forceOnBollard[1], 0); // lines pull the berth towards +y
  assert.ok(rel(tot, pull) < 1e-6);
  // rising tide (+2 m) lifts the fairleads: every line gets longer and steeper -> tension up
  const r1 = bm.berthEquilibrium({ ...model, condition: { waterLevel: 2, draft: 14 } });
  r1.lines.forEach((ln, i) => {
    assert.ok(ln.verticalAngle < r0.lines[i].verticalAngle); // bollard further below the fairlead
    assert.ok(ln.tension > r0.lines[i].tension);
  });
  const b = model.lines[2];
  const dz = 7 - (2 - 14 + b.fairlead[2]);
  assert.ok(Math.abs(Math.tan(r1.lines[2].verticalAngle) - dz / Math.hypot(r1.lines[2].length * Math.cos(r1.lines[2].verticalAngle), 0)) < 1e-9);
  // nonlinear load-elongation curve: tension follows MBL * f(strain)
  const curve = { strain: [0, 0.01, 0.02], loadFraction: [0, 0.1, 0.5] };
  const lines = [-60, 60].map((x) => ({ fairlead: [x, -21, 20], bollard: [x, -20, 6], MBL: 1e6, curve }));
  const rc = bm.berthEquilibrium({ vessel, pose0: { x: 0, y: 23, headingDeg: 0 }, lines, extra: { sway: 4e5 } });
  // each line carries 2e5 N = 0.2 MBL -> strain = 0.01 + 0.01 * (0.1/0.4) = 0.0125 on 22 m
  assert.ok(rel(rc.excursion.sway, 0.0125 * 22) < 1e-6);
});

test('OCIMF standard criteria presets and case envelope; MEG line limits', () => {
  const cases = bm.ocimfStandardCases();
  assert.equal(cases.length, 12 * 8);
  assert.ok(rel(cases[0].wind.speed, 30.8667) < 1e-5); // 60 kn
  assert.ok(rel(cases[0].current.speed, 3 * 1852 / 3600) < 1e-12);
  assert.ok(cases.some((c) => rel(c.current.speed, 0.75 * bm.KNOT) < 1e-12 && c.current.thetaDeg === 90));
  assert.ok(cases.some((c) => rel(c.current.speed, 2 * bm.KNOT) < 1e-12 && c.current.thetaDeg === 170));
  assert.deepEqual({ ...bm.MEG_LINE_LIMITS }, { steel: 0.55, synthetic: 0.5, polyamide: 0.45 });
  const env = bm.runBerthCases(tankerBerth(), bm.ocimfStandardCases({ windStepDeg: 90 }));
  assert.equal(env.cases.length, 32);
  assert.ok(env.cases.every((c) => c.converged));
  assert.ok(env.maxLineUc > 0 && env.governing.line !== null && env.maxSway > 0);
});

test('unrestrained vessel and invalid input are reported', () => {
  // a single breast line cannot hold a surge load
  const lines = [{ fairlead: [0, -21, 20], bollard: [0, -20, 6], MBL: 1e6, stiffness: 1e6 }];
  assert.throws(() => bm.berthEquilibrium({ vessel, pose0: { x: 0, y: 23, headingDeg: 0 }, lines, extra: { surge: 1e5 } }), /not restrained|converge/);
  assert.throws(() => bm.berthEquilibrium({ vessel, pose0: { x: 0, y: 23 }, lines: [{ fairlead: [0, 0, 0], bollard: [0, 0, 0], MBL: -1 }] }), /MBL/);
  assert.throws(() => bm.berthEquilibrium({ vessel, pose0: { x: 0, y: 23 }, lines: [{ fairlead: [0, -21, 20], bollard: [0, -20, 6], MBL: 1e6, stiffness: 1e6, material: 'kevlar' }] }), /limit/);
});
