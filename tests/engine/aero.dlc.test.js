import test from 'node:test';
import assert from 'node:assert/strict';
import { LOAD_SAFETY_FACTORS, genericSeaState, windWaveBins, runDesignLoadCases } from '../../src/engine/aero/dlc.js';
import { cpSurface } from '../../src/engine/aero/bem.js';
import { controllerFromRotor } from '../../src/engine/aero/control.js';
import { nrel5mw } from '../../src/engine/aero/rotors.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const g = 9.80665;

test('IEC 61400-1 partial safety factors for loads: N 1.35, A 1.1, T 1.5, fatigue 1.0', () => {
  assert.deepEqual(LOAD_SAFETY_FACTORS, { N: 1.35, A: 1.1, T: 1.5, F: 1.0 });
});

test('generic fetch-limited sea state by hand (SPM 1984): 20 m/s at a 100 m hub, 100 km fetch → Hs 2.34 m, Tp 7.06 s', () => {
  // U10 = 20 × (10/100)^0.14 = 14.489 m/s; gF/U² = 4671.4
  //   Hs = (U²/g) × 1.6e-3 × √4671.4 = 21.407 × 0.10936 = 2.341 m;  Tp = (U/g) × 0.2857 × 4671.4^(1/3) = 7.057 s
  const s = genericSeaState({ Vhub: 20, hubHeight: 100 });
  const U = 20 * 0.1 ** 0.14; const x = g * 1e5 / (U * U);
  rel(s.U10, U, 1e-12); rel(s.Hs, U * U / g * 1.6e-3 * Math.sqrt(x), 1e-12); rel(s.Tp, U / g * 0.2857 * Math.cbrt(x), 1e-12);
  close(s.Hs, 2.341, 2e-3); close(s.Tp, 7.057, 5e-3); assert.equal(s.generic, true);
  // unlimited fetch: fully developed sea, g Hs/U² = 0.243 and g Tp/U = 8.134
  const f = genericSeaState({ Vhub: 20, hubHeight: 100, fetch: 1e9 });
  rel(f.Hs, 0.243 * U * U / g, 1e-12); rel(f.Tp, 8.134 * U / g, 1e-12);
  // no shear: the hub speed is the 10 m speed
  rel(genericSeaState({ Vhub: 12, hubHeight: 90, shearExp: 0 }).U10, 12, 1e-12);
  assert.ok(genericSeaState({ Vhub: 25, hubHeight: 100 }).Hs > s.Hs);
});

test('wind–wave bins: Weibull probabilities sum to the operating fraction; user sea-state table is interpolated', () => {
  const o = { weibullA: 10, weibullK: 2, vIn: 4, vOut: 24, hubHeight: 100 };
  const bins = windWaveBins(o);
  assert.equal(bins.length, 10); assert.equal(bins[0].lo, 4); assert.equal(bins[9].hi, 24); assert.equal(bins[0].U, 5);
  // Rayleigh (k = 2): P(4 < V < 24) = e^(−0.16) − e^(−5.76) = 0.84899
  rel(bins.reduce((s, b) => s + b.prob, 0), Math.exp(-0.16) - Math.exp(-5.76), 1e-12);
  rel(bins[0].prob, Math.exp(-0.16) - Math.exp(-0.36), 1e-12);
  assert.ok(bins.every((b) => b.generic));
  rel(bins[3].Hs, genericSeaState({ Vhub: bins[3].U, hubHeight: 100 }).Hs, 1e-12);
  // site table: Hs = 0.1 U, Tp = 4 + 0.2 U at 5 and 25 m/s → linear in between
  const t = windWaveBins({ ...o, binWidth: 4, seaState: [{ U: 5, Hs: 0.5, Tp: 5 }, { U: 25, Hs: 2.5, Tp: 9 }] });
  assert.equal(t.length, 5);
  close(t[1].U, 10, 1e-12); close(t[1].Hs, 1.0, 1e-12); close(t[1].Tp, 6, 1e-12); assert.equal(t[1].generic, false);
  const f = windWaveBins({ ...o, seaState: (U) => ({ Hs: 0.01 * U * U, Tp: 6 }) });
  close(f[2].Hs, 0.01 * 81, 1e-12);
});

// shared model: NREL 5 MW rotor on a single supplied fore–aft mode (0.28 Hz), no waves
const rotor = nrel5mw();
const surface = cpSurface({ rotor, bem: { tiltDeg: 0, preconeDeg: 0 } });
const controller = controllerFromRotor({ rotor, surface });
const mMode = 4e5; const kMode = mMode * (2 * Math.PI * 0.28) ** 2;
const structure = { zHub: 100, zMudline: -20, zTowerBase: 10, modes: [{ mass: mMode, stiffness: kMode, zeta: 0.02, phiHub: 1, slopeHub: 0.015, phiStations: [] }], stations: [] };
const model = { surface, turbine: { inertia: rotor.drivetrainInertia, hubHeight: 100 }, controller, structure, dt: 0.05 };
const site = { hubHeight: 100, rotorDiameter: 2 * surface.R, waterDepth: 0, Vref: 50, Iref: 0.14, ratedWind: 11.4, cutIn: 3, cutOut: 25, weibullA: 10, weibullK: 2 };

test('DLC 6.1 (parked, 50-year wind): characteristic mudline moment brackets the quasi-static thrust × lever; design = 1.35 × characteristic', () => {
  // turbulent extreme wind model: 10-minute mean V50 at hub = Vref = 50 m/s with σ1 = 0.11 Vhub (the 1.4 Vref gust
  // belongs to the steady model). Parked thrust coefficient 0.01 on the swept area:
  //   T = ½ ρ A Ct V² = 0.5 × 1.225 × π 63² × 0.01 × 50² = 190.9 kN; lever hub to mudline 120 m → 22.9 MN·m mean
  // (a 10-minute record: a short window of 11 % turbulence can sit in a lull below the mean for its whole length)
  const r = runDesignLoadCases({ model, site, cases: ['6.1'], parkedCt: 0.01, duration: 600, tTransient: 20, yawErrors61Deg: [0] });
  assert.equal(r.cases.length, 1);
  const c = r.cases[0];
  assert.equal(c.dlc, '6.1'); assert.equal(c.type, 'N'); assert.equal(c.psf, 1.35); assert.equal(c.runs.length, 1);
  const A = Math.PI * surface.R ** 2; const T = 0.5 * 1.225 * A * 0.01 * 50 ** 2;
  assert.ok(c.characteristic.mudlineMoment > 0.95 * T * 120 && c.characteristic.mudlineMoment < 2.5 * T * 120,
    `mudline moment ${c.characteristic.mudlineMoment} vs quasi-static mean ${T * 120}`);
  assert.match(c.runs[0].label, /V = 50\.0 m\/s/);
  assert.ok(c.characteristic.mudlineShear > 0.95 * T && c.characteristic.mudlineShear < 2.5 * T);
  // tower base 10 m above datum: lever 90 m, so the moment is smaller in about that proportion
  assert.ok(c.characteristic.towerBaseMoment < c.characteristic.mudlineMoment && c.characteristic.towerBaseMoment > 0.6 * c.characteristic.mudlineMoment);
  for (const k of ['mudlineMoment', 'towerBaseMoment', 'mudlineShear']) rel(c.design[k], 1.35 * c.characteristic[k], 1e-12);
  close(c.characteristic.maxRotorSpeed, 0, 1e-9); // parked
  assert.equal(r.table.length, 1); rel(r.table[0].mudlineMomentDesign, c.design.mudlineMoment, 1e-12);
  assert.ok(r.notes.some((n) => /seed/.test(n)) && r.notes.some((n) => /without waves|no extreme sea state/i.test(n)));
  assert.equal(r.fatigue, null);
  // a parked rotor with twice the thrust coefficient: about twice the load (same turbulence seed)
  const r2 = runDesignLoadCases({ model, site, cases: ['6.1'], parkedCt: 0.02, duration: 600, tTransient: 20, yawErrors61Deg: [0] });
  rel(r2.cases[0].characteristic.mudlineShear, 2 * c.characteristic.mudlineShear, 0.05);
});

test('DLC 6.2 is an abnormal case (1.1) over yaw errors; unknown cases and a missing parked coefficient are rejected', () => {
  const ct = (yaw) => 0.01 + 0.02 * Math.abs(Math.sin(yaw * Math.PI / 180));
  const r = runDesignLoadCases({ model, site, cases: ['6.2'], parkedCt: ct, duration: 30, tTransient: 10, yawErrors62Deg: [0, 90] });
  const c = r.cases[0];
  assert.equal(c.type, 'A'); assert.equal(c.psf, 1.1); assert.equal(c.runs.length, 2);
  rel(c.design.mudlineMoment, 1.1 * c.characteristic.mudlineMoment, 1e-12);
  // the governing run is the broadside one (three times the thrust coefficient)
  assert.match(c.governing.mudlineMoment, /yaw error 90/);
  assert.ok(c.runs[1].extremes.mudlineMoment > 2 * c.runs[0].extremes.mudlineMoment);
  assert.throws(() => runDesignLoadCases({ model, site, cases: ['6.1'] }), /parkedCt/);
  assert.throws(() => runDesignLoadCases({ model, site, cases: ['9.9'] }), RangeError);
  assert.throws(() => runDesignLoadCases({ model, site: { ...site, Vref: 0 }, cases: ['6.1'], parkedCt: 0.01 }), RangeError);
});

test('DLC 1.3 (power production, extreme turbulence) and DLC 1.2 fatigue bookkeeping', () => {
  const r = runDesignLoadCases({ model, site, cases: ['1.3', '1.2'], duration: 40, tTransient: 20, windSpeeds: [11.4], binWidth: 11, nEq: 1e7, wohlerM: 4 });
  const c = r.cases[0];
  assert.equal(c.dlc, '1.3'); assert.equal(c.type, 'N'); assert.equal(c.runs.length, 1);
  // rated thrust of the 5 MW rotor is of the order of 0.7–0.9 MN: the peak base shear in turbulence lies around it
  assert.ok(c.characteristic.mudlineShear > 4e5 && c.characteristic.mudlineShear < 2.5e6, `shear ${c.characteristic.mudlineShear}`);
  assert.ok(c.characteristic.maxRotorSpeed > 0.9 && c.characteristic.maxRotorSpeed < 1.6, `rotor speed ${c.characteristic.maxRotorSpeed} rad/s`);
  // fatigue: two bins (3–14 and 14–25 m/s); lifetime DEL^m × nEq = Σ (bin 1-Hz DEL)^m × lifetime × probability
  const f = r.fatigue;
  assert.equal(f.bins.length, 2); assert.equal(f.wohlerM, 4); assert.equal(f.psf, 1.0);
  rel(f.lifetimeSeconds, 25 * 365.25 * 86400, 1e-12);
  for (const ch of ['mudlineMoment', 'towerBaseMoment']) {
    const sum = f.bins.reduce((s, b) => s + b.del[ch] ** 4 * f.lifetimeSeconds * b.prob, 0);
    rel(f.del[ch] ** 4 * 1e7, sum, 1e-9, ch); assert.ok(f.del[ch] > 0);
  }
  const row = r.table.find((t) => t.dlc === '1.2');
  rel(row.mudlineMomentDel, f.del.mudlineMoment, 1e-12); assert.equal(row.type, 'F');
  assert.ok(r.notes.some((n) => /GENERIC/.test(n)));
  // a damage hook receives lifetime-scaled cycle counts
  let got = null;
  const h = runDesignLoadCases({ model, site, cases: ['1.2'], duration: 40, tTransient: 20, binWidth: 22, snDamage: (x) => { if (x.channel === 'mudlineMoment') got = x; return 0.25; } });
  assert.equal(h.fatigue.damage.mudlineMoment, 0.25);
  const m4 = got.ranges.reduce((s, rg, i) => s + got.counts[i] * rg ** 4, 0);
  rel(m4, h.fatigue.del.mudlineMoment ** 4 * 1e7, 1e-9);
});
