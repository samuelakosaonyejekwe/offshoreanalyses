import test from 'node:test';
import assert from 'node:assert/strict';
import * as C from '../../src/engine/electrical/collection.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const sub = { x: 0, y: 0 };
const line = (n, s = 1000) => Array.from({ length: n }, (_, i) => ({ x: s * (i + 1), y: 0 }));
// own two-entry catalogue: ratings 300 A and 600 A, pure series resistance
const cat = [{ area: 150e-6, I: 300, r: 2e-4, x: 0, c: 0, cost: 100 }, { area: 400e-6, I: 600, r: 1e-4, x: 0, c: 0, cost: 200 }];

test('route length, radial strings and network bookkeeping', () => {
  // 3-4-5 triangle: 500 m × 1.03 + 2 × 20 m of J-tube = 555 m
  close(C.routeLength({ a: sub, b: { x: 300, y: 400 }, endLength: 20 }), 555, 1e-9);
  close(C.routeLength({ a: sub, b: { x: 300, y: 400 }, slackFactor: 1 }), 500, 1e-12);
  assert.deepEqual(C.radialStrings({ nTurbines: 4, strings: [[0, 1], [2, 3]] }), [0, 1, 0, 3]);
  assert.throws(() => C.radialStrings({ nTurbines: 4, strings: [[0, 1], [1, 3]] }), RangeError);
  // chain of three turbines 1 km apart
  const n = C.buildNetwork({ substation: sub, turbines: line(3), parent: [0, 1, 2], slackFactor: 1 });
  assert.deepEqual(Array.from(n.length).slice(1), [1000, 1000, 1000]); assert.deepEqual(Array.from(n.downstream).slice(1), [3, 2, 1]);
  assert.deepEqual(Array.from(n.order), [0, 1, 2, 3]); close(n.totalLength, 3000, 1e-9); assert.equal(n.nodes.length, 4);
  // a branched feeder: turbines 2 and 3 both hang off turbine 1
  const b = C.buildNetwork({ substation: sub, turbines: [{ x: 1000, y: 0 }, { x: 2000, y: 500 }, { x: 2000, y: -500 }], parent: [0, 1, 1], slackFactor: 1 });
  assert.deepEqual(Array.from(b.downstream).slice(1), [3, 1, 1]); rel(b.totalLength, 1000 + 2 * Math.hypot(1000, 500), 1e-12);
  assert.throws(() => C.buildNetwork({ substation: sub, turbines: line(2), parent: [2, 1] }), RangeError); // loop without a path to the substation
});

test('ring layout: closing link between string ends; every ring cable is sized for both strings', () => {
  const turbines = [{ x: 1000, y: 500 }, { x: 2000, y: 500 }, { x: 1000, y: -500 }, { x: 2000, y: -500 }];
  const r = C.ringLayout({ substation: sub, turbines, strings: [[0, 1], [2, 3]], rings: [[0, 1]], slackFactor: 1 });
  assert.equal(r.links.length, 1); assert.deepEqual([r.links[0].a, r.links[0].b].sort(), [2, 4]); close(r.links[0].length, 1000, 1e-9);
  assert.deepEqual(Array.from(r.ringDownstream).slice(1), [4, 4, 4, 4]);
  rel(r.totalLength, r.network.totalLength + 1000, 1e-12);
  assert.deepEqual(Array.from(r.network.downstream).slice(1), [2, 1, 2, 1]);
});

test('Esau–Williams: turbines on a line form one chain when capacity allows and two feeders when it does not', () => {
  const one = C.esauWilliams({ substation: sub, turbines: line(4), maxPerString: 4 });
  assert.deepEqual(one.parent, [0, 1, 2, 3]); close(one.totalDistance, 4000, 1e-9); assert.equal(one.strings.length, 1); assert.equal(one.crossings, 0);
  // capacity 2: the far pair needs its own 3 km gate → 2000 + 3000 + 1000 = 6000 m (collinear cables overlap, so crossings are not penalised here)
  const two = C.esauWilliams({ substation: sub, turbines: line(4), maxPerString: 2, noCrossing: false });
  close(two.totalDistance, 6000, 1e-9); assert.equal(two.strings.length, 2);
  for (const s of two.strings) assert.ok(s.length <= 2);
  // star when every turbine must be alone
  const star = C.esauWilliams({ substation: sub, turbines: line(3), maxPerString: 1 });
  assert.deepEqual(star.parent, [0, 0, 0]); close(star.totalDistance, 6000, 1e-9);
  // two clusters either side of the substation are never joined across it
  const tw = [{ x: 1000, y: 0 }, { x: 2000, y: 0 }, { x: -1000, y: 0 }, { x: -2000, y: 0 }];
  const ew = C.esauWilliams({ substation: sub, turbines: tw, maxPerString: 4 });
  assert.deepEqual(ew.parent, [0, 1, 0, 3]); close(ew.totalDistance, 4000, 1e-9);
  // every capacity limit is respected on a grid; the tree is never longer than the star
  const grid = []; for (let i = 0; i < 4; i++) for (let j = 0; j < 3; j++) grid.push({ x: 1500 + 1000 * i, y: 1000 * (j - 1) });
  const g = C.esauWilliams({ substation: sub, turbines: grid, maxPerString: 4 });
  const net = C.buildNetwork({ substation: sub, turbines: grid, parent: g.parent, slackFactor: 1 });
  assert.ok(Math.max(...Array.from(net.downstream).slice(1)) <= 4);
  assert.ok(g.totalDistance < grid.reduce((s, t) => s + Math.hypot(t.x, t.y), 0)); rel(net.totalLength, g.totalDistance, 1e-9);
  const chain = C.esauWilliams({ substation: sub, turbines: grid, maxPerString: 4, noBranching: true });
  const kids = new Array(13).fill(0); chain.parent.forEach((p) => { if (p > 0) kids[p]++; });
  assert.ok(Math.max(...kids) <= 1, 'no branching: each turbine feeds at most one other');
});

test('cable sizing by turbine count: I = n P/(√3 U cos φ); smallest adequate section; limited number of types', () => {
  // 10 MW at 66 kV, unity power factor: 87.48 A per turbine → 3 turbines 262 A (150 mm²), 4 turbines 350 A (400 mm²)
  const net = C.buildNetwork({ substation: sub, turbines: line(4), parent: [0, 1, 2, 3], slackFactor: 1 });
  const s = C.sizeCables({ network: net, catalogue: cat, turbinePower: 10e6, voltage: 66e3, powerFactor: 1 });
  rel(s.current[4], 10e6 / (Math.sqrt(3) * 66e3), 1e-12); rel(s.current[1], 4 * 87.477, 1e-4);
  assert.deepEqual(Array.from(s.type).slice(1), [1, 0, 0, 0]);
  rel(s.utilisation[1], 349.91 / 600, 1e-4); rel(s.utilisation[2], 262.43 / 300, 1e-4);
  rel(s.totalCost, 1000 * 200 + 3000 * 100, 1e-12); rel(s.totalLength, 4000, 1e-12);
  const q = Object.fromEntries(s.quantities.map((x) => [Math.round(x.area * 1e6), x.length]));
  close(q[150], 3000, 1e-9); close(q[400], 1000, 1e-9);
  assert.deepEqual(C.cableQuantities(s.quantities).sort((a, b) => a.sizeMm2 - b.sizeMm2), [{ sizeMm2: 150, km: 3 }, { sizeMm2: 400, km: 1 }]);
  // one cable type only: everything in 400 mm²
  const one = C.sizeCables({ network: net, catalogue: cat, turbinePower: 10e6, voltage: 66e3, powerFactor: 1, maxTypes: 1 });
  assert.deepEqual(Array.from(one.type).slice(1), [1, 1, 1, 1]); rel(one.totalCost, 4000 * 200, 1e-12);
  // 0.95 power factor raises the current by 1/0.95; too many turbines for the catalogue is an error
  rel(C.sizeCables({ network: net, catalogue: cat, turbinePower: 10e6, voltage: 66e3 }).current[4], 87.477 / 0.95, 1e-4);
  assert.throws(() => C.sizeCables({ network: net, catalogue: cat, turbinePower: 20e6, voltage: 66e3, powerFactor: 1 }), RangeError);
});

test('radial load flow: one resistive branch has the closed form V1 = ½[V0 + √(V0² + 4 R P/3)]; losses are 3 I² R', () => {
  // 5 km of 2e-4 Ω/m → R = 1 Ω; 10 MW at unity power factor into a 33 kV busbar
  const net = C.buildNetwork({ substation: sub, turbines: [{ x: 5000, y: 0 }], parent: [0], slackFactor: 1 });
  const r = C.radialLoadFlow({ network: net, catalogue: cat, type: [0, 0], voltage: 33e3, P: 10e6 });
  const V0 = 33e3 / Math.sqrt(3); const V1 = 0.5 * (V0 + Math.sqrt(V0 * V0 + 4 * 1 * 10e6 / 3)); const I = 10e6 / 3 / V1;
  rel(r.V[1].re, V1, 1e-9); close(r.V[1].im, 0, 1e-6); rel(r.U[1], V1 * Math.sqrt(3), 1e-9); rel(r.I[1], I, 1e-9);
  rel(r.losses, 3 * I * I, 1e-8); rel(r.branchLoss[1], r.losses, 1e-12); rel(r.Psubstation, 10e6 - r.losses, 1e-9);
  rel(r.Pgenerated, 10e6, 1e-12); rel(r.lossFraction, r.losses / 10e6, 1e-9); rel(r.Umax, r.U[1], 1e-12);
  assert.ok(r.iterations < 20);
  // chain of two turbines: the first branch carries both currents — losses about (2² + 1²) I² R per phase-set
  const net2 = C.buildNetwork({ substation: sub, turbines: [{ x: 5000, y: 0 }, { x: 10000, y: 0 }], parent: [0, 1], slackFactor: 1 });
  const r2 = C.radialLoadFlow({ network: net2, catalogue: cat, type: [0, 0, 0], voltage: 33e3, P: 10e6 });
  // (currents evaluated at the nominal voltage; the voltage rise along the string lowers them by 2–3 %)
  rel(r2.losses, 5 * 3 * (10e6 / 3 / V0) ** 2 * 1, 0.06); rel(r2.Psubstation + r2.losses, 20e6, 1e-9);
  // exact: with the solved voltages, I2 = P/(3 V2) and I1 = I2 + P/(3 V1)
  { const i2 = 10e6 / 3 / r2.V[2].re; const i1 = i2 + 10e6 / 3 / r2.V[1].re; rel(r2.losses, 3 * (i1 * i1 + i2 * i2), 1e-8); }
  assert.ok(r2.U[2] > r2.U[1] && r2.U[1] > 33e3);
  // an unloaded cable with capacitance exports its charging power ω C U² to the busbar: 0.25 µF/km × 5 km at 33 kV → 427.6 kvar
  const cc = [{ ...cat[0], c: 0.25e-9 }];
  const q = C.radialLoadFlow({ network: net, catalogue: cc, type: [0, 0], voltage: 33e3, P: 0 });
  rel(Math.abs(q.Qsubstation), 2 * Math.PI * 50 * 0.25e-9 * 5000 * 33e3 ** 2, 2e-3);
  // per-turbine power array and reactive injection
  const pq = C.radialLoadFlow({ network: net2, catalogue: cat, type: [0, 0, 0], voltage: 33e3, P: [10e6, 0], Q: [0, 0] });
  rel(pq.Pgenerated, 10e6, 1e-12); rel(pq.losses, r.losses, 1e-6);
});

test('annual energy loss weights the load-flow losses by the production distribution; energy not delivered by cable faults', () => {
  const net = C.buildNetwork({ substation: sub, turbines: [{ x: 5000, y: 0 }], parent: [0], slackFactor: 1 });
  const lf = { network: net, catalogue: cat, type: [0, 0], voltage: 33e3 };
  const rated = C.radialLoadFlow({ ...lf, P: 10e6 }).losses;
  const half = C.radialLoadFlow({ ...lf, P: 5e6 }).losses;
  // half the year at rated, a quarter at half power, a quarter at standstill
  const a = C.annualEnergyLoss({ ...lf, turbinePower: 10e6, distribution: [{ fraction: 1, probability: 0.5 }, { fraction: 0.5, probability: 0.25 }, { fraction: 0, probability: 0.25 }] });
  rel(a.energyLoss, 8760 * (0.5 * rated + 0.25 * half), 1e-9); rel(a.energyGenerated, 8760 * 10e6 * 0.625, 1e-9);
  rel(a.lossFraction, a.energyLoss / a.energyGenerated, 1e-12); rel(a.lossAtRated, rated, 1e-9);
  // losses scale almost with the square of the load: half power → about a quarter
  rel(half, rated / 4, 0.01);
  // probabilities are normalised
  rel(C.annualEnergyLoss({ ...lf, turbinePower: 10e6, distribution: [{ fraction: 1, probability: 2 }, { fraction: 0, probability: 2 }] }).energyLoss, 8760 * 0.5 * rated, 1e-9);
  // reliability: chain of 3 turbines, 1 km branches, λ = 5e-6 /m/yr (0.005 per km·yr), MTTR 720 h, 10 MW at CF 0.5
  //   E = λ L MTTR P CF (3 + 2 + 1) = 5e-6 × 1000 × 720 × 5e6 × 6 = 108 MWh/yr
  const n3 = C.buildNetwork({ substation: sub, turbines: line(3), parent: [0, 1, 2], slackFactor: 1 });
  const e = C.energyNotDelivered({ network: n3, turbinePower: 10e6, capacityFactor: 0.5, failureRate: 5e-6, mttr: 720 });
  rel(e.energyNotDelivered, 5e-6 * 1000 * 720 * 5e6 * 6, 1e-12); rel(e.expectedFailures, 0.015, 1e-12);
  rel(e.unavailability, e.energyNotDelivered / (3 * 5e6 * 8760), 1e-12);
  // a ring restores supply after the 4 h switching time
  rel(C.energyNotDelivered({ network: n3, turbinePower: 10e6, capacityFactor: 0.5, failureRate: 5e-6, mttr: 720, ring: true }).energyNotDelivered, e.energyNotDelivered * 4 / 720, 1e-12);
});

test('generic array-cable catalogue: ratings rise and resistance falls with section; default cost rule', () => {
  assert.equal(C.ARRAY_CABLE_SIZES.length, 10); close(C.ARRAY_CABLE_SIZES[0], 95e-6, 1e-15);
  const cable = C.genericArrayCable({ voltage: 66e3, area: 400e-6 });
  assert.equal(cable.cores, 3); close(cable.conductor.area, 400e-6, 1e-15);
  const k = C.genericCatalogue({ voltage: 66e3, areas: [150e-6, 400e-6, 800e-6] });
  assert.equal(k.length, 3);
  assert.ok(k[0].I < k[1].I && k[1].I < k[2].I && k[0].r > k[1].r && k[1].r > k[2].r);
  // plausible buried three-core 66 kV copper ratings: a few hundred amperes to about 1 kA
  assert.ok(k[0].I > 250 && k[0].I < 550 && k[2].I > 600 && k[2].I < 1200, `ratings ${k.map((x) => Math.round(x.I))}`);
  // 400 mm² Cu at 90 °C: DC 0.0470 × 1.275 = 0.060 mΩ/m; with skin, proximity, screen and armour losses somewhat higher
  assert.ok(k[1].r > 0.058e-3 && k[1].r < 0.10e-3, `r = ${k[1].r}`);
  rel(k[1].cost, (100 + 0.40 * 400) * 1.15, 1e-9); rel(C.genericCatalogue({ voltage: 33e3, areas: [400e-6] })[0].cost, 260, 1e-9);
  assert.ok(k.every((x) => x.generic === true && x.x > 0 && x.c > 0));
  // aluminium carries less than copper of the same section
  assert.ok(C.genericCatalogue({ voltage: 66e3, areas: [400e-6], material: 'Al' })[0].I < k[1].I);
  rel(C.genericCatalogue({ voltage: 66e3, areas: [400e-6], cost: (a) => 1e6 * a })[0].cost, 400, 1e-9);
});

test('complete design and voltage comparison: 66 kV needs fewer strings and loses less than 33 kV', () => {
  const turbines = []; for (let i = 0; i < 4; i++) for (let j = 0; j < 3; j++) turbines.push({ x: 2000 + 1500 * i, y: 1500 * (j - 1) });
  const d = C.designCollection({ substation: sub, turbines, turbinePower: 10e6, voltage: 66e3, failureRate: 5e-6, mttr: 720, capacityFactor: 0.5,
    distribution: [{ fraction: 1, probability: 0.4 }, { fraction: 0.5, probability: 0.4 }, { fraction: 0, probability: 0.2 }] });
  assert.equal(d.parent.length, 12); assert.equal(d.voltage, 66e3);
  rel(d.loadFlow.Pgenerated, 120e6, 1e-12); assert.ok(d.loadFlow.lossFraction > 0 && d.loadFlow.lossFraction < 0.02);
  assert.ok(Math.max(...Array.from(d.network.downstream).slice(1)) <= d.maxPerString);
  assert.ok(Math.max(...Array.from(d.sizing.utilisation).slice(1)) <= 1 + 1e-9);
  rel(d.cableKm.reduce((s, c) => s + c.km, 0), d.network.totalLength / 1000, 1e-9);
  assert.ok(d.annual.lossFraction < d.loadFlow.lossFraction && d.reliability.energyNotDelivered > 0); rel(d.totalCost, d.sizing.totalCost, 1e-12);
  const cmp = C.compareCollectionVoltages({ substation: sub, turbines, turbinePower: 10e6 });
  assert.equal(cmp.length, 2); assert.equal(cmp[0].voltage, 33e3);
  assert.ok(cmp[1].nStrings <= cmp[0].nStrings && cmp[1].maxPerString >= cmp[0].maxPerString);
  assert.ok(cmp[1].lossFractionAtRated < cmp[0].lossFractionAtRated);
  assert.equal(cmp[0].annualLossFraction, null);
});

test('export system: transformer losses P0 + Pk (S/Sr)², HVDC link losses by hand', () => {
  // 400 MW through two 250 MVA transformers (Pk 600 kW, P0 100 kW each): 2 × (100 + 600 × 0.8²) = 968 kW
  const t = C.exportSystem({ P: 400e6, transformer: { Sr: 250e6, uk: 0.12, Pk: 600e3, P0: 100e3 }, nTransformers: 2 });
  rel(t.transformerLoss, 2 * (100e3 + 600e3 * 0.64), 1e-9); assert.deepEqual(t.options, []);
  // symmetric monopole ±320 kV, 100 km of 0.02 Ω/km per pole, 1 % per converter station, no transformer:
  //   I = 400e6 × 0.99/640e3 = 618.75 A; cable loss I² × 2 × 2 Ω = 1.531 MW; received = (396 − 1.531) × 0.99 = 390.52 MW
  const d = C.exportSystem({ P: 400e6, dc: { configuration: 'symmetricMonopole', U: 320e3, r: 2e-5, length: 100e3, converterLossFraction: 0.01 } });
  const I = 400e6 * 0.99 / 640e3;
  assert.equal(d.options.length, 1); assert.equal(d.options[0].name, 'HVDC'); rel(d.options[0].ImaxCable, I, 1e-12);
  rel(d.options[0].Pdelivered, (396e6 - I * I * 4) * 0.99, 1e-6); rel(d.options[0].losses, 400e6 - d.options[0].Pdelivered, 1e-12);
  rel(d.options[0].efficiency, d.options[0].Pdelivered / 400e6, 1e-12);
});
