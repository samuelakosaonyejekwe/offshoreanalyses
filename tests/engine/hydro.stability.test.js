import test from 'node:test';
import assert from 'node:assert/strict';
import * as st from '../../src/engine/hydro/stability.js';
import { boxMesh } from '../../src/engine/hydro/hydrostatics.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b}`);
const DEG = Math.PI / 180;
const rho = 1025;
const barge = () => boxMesh({ xMin: -50, xMax: 50, yMin: -10, yMax: 10, zMin: 0, zMax: 10 });
const mass = rho * 100 * 20 * 4;

test('GZ of a box barge equals the wall-sided formula GZ = sin(phi) (GM + ½ BM tan² phi) below bilge emergence', () => {
  const KG = 5; const BM = 400 / 48; const GM = 2 + BM - KG;
  const heels = [0, 2.5, 5, 10, 15, 20];   // bilge emerges at atan(4/10) = 21.8°, deck edge at 31°
  for (const side of [1, -1]) {
    const c = st.gzCurve({ mesh: barge(), mass, cog: [0, 0, KG], rho, heelsDeg: heels, side });
    heels.forEach((h, i) => {
      const phi = h * DEG;
      near(c.gz[i], Math.sin(phi) * (GM + 0.5 * BM * Math.tan(phi) ** 2), 1e-7, `heel ${h}`);
    });
  }
});

test('GZ curve: free trim equals fixed trim for a fore-aft symmetric box; downflooding angle of an opening', () => {
  const heels = [0, 10, 20, 30, 40];
  const a = st.gzCurve({ mesh: barge(), mass, cog: [0, 0, 5], rho, heelsDeg: heels });
  const b = st.gzCurve({ mesh: barge(), mass, cog: [0, 0, 5], rho, heelsDeg: heels, freeTrim: true });
  heels.forEach((h, i) => { near(a.gz[i], b.gz[i], 1e-6); near(b.trimDeg[i], 0, 1e-5); });
  // opening in the starboard side 2 m above the waterline (y = -10, z = 6): in the wall-sided range the waterline
  // pivots about the centreline, so it immerses at atan(2/10) = 11.31°
  const c = st.gzCurve({ mesh: barge(), mass, cog: [0, 0, 5], rho, heelsDeg: Array.from({ length: 41 }, (_, i) => 0.5 * i), openings: [{ point: [0, -10, 6] }] });
  near(c.downfloodDeg, Math.atan(0.2) / DEG, 0.01);
});

test('GZ properties: areas of GZ = GM sin(phi) are GM (cos a - cos b); maximum and vanishing angle of sin(2 phi)', () => {
  const heelDeg = Array.from({ length: 181 }, (_, i) => 0.5 * i);
  const gz = heelDeg.map((h) => 1.5 * Math.sin(h * DEG));
  const p = st.gzProperties({ heelDeg, gz });
  close(p.area0to30, 1.5 * (1 - Math.cos(30 * DEG)), 1e-5);
  close(p.area30to40, 1.5 * (Math.cos(30 * DEG) - Math.cos(40 * DEG)), 1e-5);
  close(p.gm0, 1.5, 1e-4);
  const q = st.gzProperties({ heelDeg: heelDeg.map((h) => h), gz: heelDeg.map((h) => Math.sin(2.5 * h * DEG)) });
  near(q.heelAtGzMaxDeg, 36, 0.01); close(q.gzMax, 1, 1e-6); near(q.vanishingDeg, 72, 0.01);
});

test('IMO 2008 IS Code general criteria: limits 0.055 / 0.09 / 0.03 m·rad, 0.20 m, 25°, 0.15 m', () => {
  const heelDeg = Array.from({ length: 33 }, (_, i) => 2.5 * i);
  const good = st.imoIntactCriteria({ heelDeg, gz: heelDeg.map((h) => 0.6 * Math.sin(2 * h * DEG)) });
  assert.ok(good.pass);
  assert.deepEqual(good.checks.map((c) => c.limit), [0.055, 0.09, 0.03, 0.2, 25, 0.15]);
  // GZ = 0.05 sin(2 phi): GM0 = 0.1 m and all areas too small
  const bad = st.imoIntactCriteria({ heelDeg, gz: heelDeg.map((h) => 0.05 * Math.sin(2 * h * DEG)) });
  assert.ok(!bad.pass);
  assert.ok(!bad.checks[5].ok && !bad.checks[0].ok && bad.checks[4].ok);
  // downflooding at 35° truncates the 0–40° area
  const df = st.imoIntactCriteria({ heelDeg, gz: heelDeg.map((h) => 0.6 * Math.sin(2 * h * DEG)), downfloodDeg: 35 });
  close(df.checks[1].value, 0.3 * (1 - Math.cos(70 * DEG)), 2e-3);
});

test('MODU Code wind criterion: GZ = sin(2 phi), constant lever 0.5 -> intercepts 15° and 75°, ratio 1.4255', () => {
  const heelDeg = Array.from({ length: 181 }, (_, i) => 0.5 * i);
  const gz = heelDeg.map((h) => Math.sin(2 * h * DEG));
  const r = st.moduWindCriterion({ heelDeg, gz, lever0: 0.5 });
  near(r.firstInterceptDeg, 15, 0.01); near(r.secondInterceptDeg, 75, 0.01);
  const aR = (1 - Math.cos(150 * DEG)) / 2; const aH = 0.5 * 75 * DEG;
  close(r.ratio, aR / aH, 1e-3);
  assert.equal(r.required, 1.4); assert.ok(r.pass);
  assert.equal(st.moduWindCriterion({ heelDeg, gz, lever0: 0.5, unitType: 'column-stabilised' }).required, 1.3);
  assert.ok(!st.moduWindCriterion({ heelDeg, gz, lever0: 0.5, downfloodDeg: 40 }).pass); // ratio to 40° = 0.413/0.349 = 1.18
});

test('Wind heeling lever: F = ½ Cs Ch rho V² A with rho = 1.222 kg/m³ (hand calculation)', () => {
  const r = st.windHeelingLever({ V: 51.5, areas: [{ A: 1000, z: 20, Cs: 1.0, Ch: 1.1 }, { A: 200, z: 40, Cs: 0.5, Ch: 1.2 }], zLateralResistance: -5, displacement: 2e7 });
  const q = 0.5 * 1.222 * 51.5 ** 2;
  close(r.force, q * (1100 + 120), 1e-12);
  close(r.lever, (q * (1100 * 25 + 120 * 45)) / (2e7 * 9.80665), 1e-12);
});

test('Damage stability, lost buoyancy: full-breadth midship compartment gives parallel sinkage T·L/(L - mu·l)', () => {
  const comp = { box: { xMin: -10, xMax: 10, yMin: -10, yMax: 10, zMin: 0, zMax: 10 } };
  const d = st.damageStability({ mesh: barge(), mass, cog: [0, 0, 5], rho, compartments: [comp], heelsDeg: [0, 5, 10, 15], marginLine: [[50, 10, 9.924], [-50, -10, 9.924]] });
  near(d.equilibrium.draft, (4 * 100) / 80, 1e-7); near(d.equilibrium.heel, 0, 1e-8); near(d.equilibrium.trim, 0, 1e-8);
  close(d.lostVolume, 20 * 20 * 5, 1e-7);
  near(d.marginLine.minHeight, 4.924, 1e-6); assert.ok(!d.marginLine.submerged);
  // residual stability: intact waterplane inertia of the two end portions only, I = 2·(40·20³/12), V unchanged
  const GMres = 2.5 + (2 * 40 * 8000 / 12) / 8000 - 5;
  close(d.properties.gm0, GMres, 5e-3);
  const half = st.damageStability({ mesh: barge(), mass, cog: [0, 0, 5], rho, compartments: [{ ...comp, permeability: 0.5 }], heelsDeg: [0, 5, 10] });
  near(half.equilibrium.draft, 400 / 90, 1e-7);
});

test('Damage stability: an off-centre wing compartment heels the barge towards the damage', () => {
  const d = st.damageStability({ mesh: barge(), mass, cog: [0, 0, 5], rho, heelsDeg: [0, 5, 10, 15, 20],
    compartments: [{ box: { xMin: 20, xMax: 40, yMin: 5, yMax: 10, zMin: 0, zMax: 10 } }] });
  // damage to port (+y) and forward (+x): port side down = negative heel, bow down = positive trim
  assert.ok(d.equilibrium.heel < -0.005 && d.equilibrium.trim > 0.001 && d.side === -1);
  near(d.curve.gz[0] < 0 ? 1 : 0, 1, 0, 'upright lever is a heeling lever');
  near(d.properties.equilibriumDeg, -d.equilibrium.heelDeg, 0.05);
});

test('Progressive flooding through a bottom opening ends at the lost-buoyancy waterline (T = 5 m), mass conserved', () => {
  const r = st.progressiveFlooding({ mesh: barge(), mass, cog: [0, 0, 5], rho, dt: 5, tMax: 6000,
    compartment: { box: { xMin: -10, xMax: 10, yMin: -10, yMax: 10, zMin: 0, zMax: 10 } }, opening: { point: [0, 0, 0], area: 0.5, Cd: 0.6 } });
  assert.ok(r.equalised && !r.sunk);
  const n = r.t.length - 1;
  near(r.draft[n], 5, 5e-3); close(r.floodMass[n], rho * 20 * 20 * 5, 2e-3);
  // initial inflow: Cd A sqrt(2 g T0)
  close(r.flow[0], 0.6 * 0.5 * Math.sqrt(2 * 9.80665 * 4), 1e-9);
  for (let i = 1; i <= n; i++) assert.ok(r.floodMass[i] >= r.floodMass[i - 1] - 1e-6);
  // hand estimate of the filling time: head h = T - level with dh/dt = -(Q/A_c)(1 - A_c/A_wp), A_c = 400, A_wp = 2000:
  // t = 2 sqrt(h0) A_c / (0.8 Cd A sqrt(2g)) = 2·2·400/(0.8·0.3·4.4287) = 1505 s
  // (explicit Euler converges to it at first order: 1490 s at dt = 5 s, 1502 s at dt = 1 s)
  close(r.timeToEqualise, (2 * 2 * 400) / (0.8 * 0.3 * Math.sqrt(2 * 9.80665)), 0.015);
});

test('Ballast plan: least-squares transfer meets heel and trim targets within tank limits; pump time', () => {
  const hydro = { displacement: mass, waterplaneArea: 2000, lcf: 0, GMt: 5.333, GMl: 205.3, rho };
  const tanks = [[40, 8], [40, -8], [-40, 8], [-40, -8]].map(([x, y], i) => ({ name: `WB${i + 1}`, x, y, capacity: 300, volume: 150 }));
  const p = st.ballastPlan({ tanks, hydro, current: { draft: 4, heelDeg: 1.0, trimDeg: 0.2 }, target: { heelDeg: 0, trimDeg: 0 }, pumpRate: 0.1, nPumps: 2 });
  assert.ok(p.feasible);
  near(p.achieved.heelDeg, 0, 1e-3); near(p.achieved.trimDeg, 0, 1e-3);
  near(p.tanks.reduce((s, t) => s + t.deltaVolume, 0), 0, 0.05);
  for (const t of p.tanks) assert.ok(t.newVolume >= -1e-9 && t.newVolume <= 300 + 1e-9);
  // independent check of the achieved heel from the moved masses
  const heel = 1.0 - p.tanks.reduce((s, t, i) => s + t.deltaMass * tanks[i].y, 0) / (mass * 5.333) / DEG;
  near(heel, p.achieved.heelDeg, 1e-9);
  close(p.pumpTime, p.volumeMoved / 0.2, 1e-12);
  // an unreachable target is reported as infeasible
  assert.ok(!st.ballastPlan({ tanks, hydro, current: { draft: 4, heelDeg: 25, trimDeg: 0 }, target: { heelDeg: 0 } }).feasible);
  // fill/discharge mode reaches a deeper draft: 0.1 m × 2000 m² = 200 m³
  const f = st.ballastPlan({ tanks, hydro, current: { draft: 4, heelDeg: 0, trimDeg: 0 }, target: { draft: 4.1, heelDeg: 0, trimDeg: 0 }, mode: 'fill-discharge' });
  assert.ok(f.feasible); close(f.volumeMoved, 200, 1e-3);
});

test('Inclining test: GM = w d / (Δ tan phi) recovered from pendulum readings', () => {
  const disp = 5e6; const GM = 1.2;
  const shifts = [[10e3, 8], [20e3, 8], [-10e3, 8]].map(([w, d]) => ({ mass: w, distance: d, pendulumLength: 3, deflection: (3 * w * d) / (disp * GM) }));
  const r = st.incliningTest({ displacement: disp, shifts, KM: 7.5 });
  close(r.GM, 1.2, 1e-12); close(r.KG, 6.3, 1e-12); near(r.spread, 0, 1e-12);
});

test('Progressive flooding end states: foundering, capsize and equalisation are told apart; sunk keeps its meaning', () => {
  // (a) the whole hull floods through a bottom opening: draft = 4 m + internal level, so the deck (10 m) goes under
  // when 6 m of water (12 000 m³) is aboard — ship and water then outweigh the 20 000 m³ hull: loss of buoyancy
  const a = st.progressiveFlooding({ mesh: barge(), mass, cog: [0, 0, 5], rho, dt: 20, tMax: 40000,
    compartment: { box: { xMin: -50, xMax: 50, yMin: -10, yMax: 10, zMin: 0, zMax: 10 } }, opening: { point: [0, 0, 0], area: 5 } });
  assert.ok(a.sunk && a.lostBuoyancy && !a.capsized && !a.solverFailed && !a.equalised);
  const na = a.t.length - 1;
  assert.ok(a.draft[na] > 9 && a.draft[na] <= 10 + 1e-6, `last draft ${a.draft[na]}`);
  assert.ok(a.floodMass[na] / rho > 10500 && a.floodMass[na] / rho < 12000 + 1e-6);
  near(a.heelDeg[na], 0, 1e-6);
  // (b) a full-length wing compartment with a high centre of gravity (KG = 9 m, GM = 2 + 8.33 − 9 = 1.33 m intact):
  // buoyancy is ample (under 2000 m³ aboard of 12 000 m³ reserve) but the ship rolls over
  const b = st.progressiveFlooding({ mesh: barge(), mass, cog: [0, 0, 9], rho, dt: 20, tMax: 40000,
    compartment: { box: { xMin: -50, xMax: 50, yMin: 4, yMax: 10, zMin: 0, zMax: 10 } }, opening: { point: [0, 10, 1], area: 3 } });
  assert.ok(b.sunk && b.capsized && !b.lostBuoyancy && !b.solverFailed && !b.equalised);
  assert.ok(Math.abs(b.capsizeHeelDeg) > 40);
  assert.ok(mass + b.floodMass[b.t.length - 1] < 0.6 * rho * 20000);
  // (c) the benign case reports none of the failure causes
  const c = st.progressiveFlooding({ mesh: barge(), mass, cog: [0, 0, 5], rho, dt: 20, tMax: 6000,
    compartment: { box: { xMin: -10, xMax: 10, yMin: -10, yMax: 10, zMin: 0, zMax: 10 } }, opening: { point: [0, 0, 0], area: 0.5, Cd: 0.6 } });
  assert.ok(c.equalised && !c.sunk && !c.lostBuoyancy && !c.capsized && !c.solverFailed);
  assert.equal(c.capsizeHeelDeg, null);
});
