import test from 'node:test';
import assert from 'node:assert/strict';
import * as b from '../../src/engine/geotech/bearing.js';

const deg = Math.PI / 180;
const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);

test('Prandtl–Reissner factors at phi = 30 deg: Nq = 18.40, Nc = 30.14 (tabulated, e.g. Bowles Table 4-4)', () => {
  close(b.bearingNq(30 * deg), 18.40, 3e-4);
  close(b.bearingNc(30 * deg), 30.14, 3e-4);
  // closed form written out: e^(π tan 30°) tan²(60°) = e^(1.81380) · 3
  close(b.bearingNq(30 * deg), Math.exp(Math.PI / Math.sqrt(3)) * 3, 1e-12);
});

test('N-gamma at 30 deg: Vesic 22.40, Hansen 15.07, Meyerhof 15.67 (tabulated)', () => {
  close(b.bearingNgamma(30 * deg, 'vesic'), 22.40, 3e-4);
  close(b.bearingNgamma(30 * deg, 'hansen'), 15.07, 5e-4);
  close(b.bearingNgamma(30 * deg, 'meyerhof'), 15.67, 5e-4);
});

test('phi = 0: Nc = 2 + pi = 5.14, Nq = 1, Ngamma = 0', () => {
  close(b.bearingNc(0), 2 + Math.PI, 1e-14);
  close(b.bearingNc(0), 5.14, 1e-3);
  close(b.bearingNq(0), 1, 1e-14);
  for (const f of ['vesic', 'hansen', 'meyerhof']) assert.equal(b.bearingNgamma(0, f), 0);
  // continuity of Nc as phi -> 0
  close(b.bearingNc(1e-6), 2 + Math.PI, 1e-5);
});

test('Terzaghi factors: Nc(0) = 5.71; at 30 deg Nq = 22.46, Nc = 37.16; Ngamma approximation within 10 % of 19.13', () => {
  const f0 = b.terzaghiFactors(0);
  close(f0.Nc, 1.5 * Math.PI + 1, 1e-14);
  const f = b.terzaghiFactors(30 * deg);
  close(f.Nq, 22.46, 5e-4);
  close(f.Nc, 37.16, 5e-4);
  assert.ok(Math.abs(f.Ngamma - 19.13) / 19.13 < 0.10);
  // exact Terzaghi expression with a supplied Kp-gamma: Ngamma = 0.5 tan(phi) (Kpg / cos^2(phi) - 1)
  close(b.terzaghiFactors(30 * deg, { Kpgamma: 52 }).Ngamma, 0.5 * Math.tan(30 * deg) * (52 / 0.75 - 1), 1e-12);
});

test('Terzaghi equation by hand: square footing on clay, qult = 1.3 c Nc + q', () => {
  const r = b.terzaghiBearing({ shape: 'square', B: 2, c: 50e3, phi: 0, gamma: 8e3, q: 16e3 });
  close(r.qult, 1.3 * 50e3 * (1.5 * Math.PI + 1) + 16e3, 1e-12);
});

test('Skempton Nc: strip at surface 5.14, square at surface 6.17, deep limits 7.5 and 9.0', () => {
  close(b.skemptonNc({ B: 2 }), 5.14, 1e-3);
  close(b.skemptonNc({ B: 2, L: 2 }), 6.17, 1e-3);
  close(b.skemptonNc({ B: 2, D: 20 }), 7.5, 1e-12);
  close(b.skemptonNc({ B: 2, L: 2, D: 20 }), 9.0, 1e-12);
  close(b.skemptonNc({ B: 2, L: 2, D: 2, form: 'simple' }), 5 * 1.2 * 1.2, 1e-12);
});

test('Meyerhof effective area by hand: 4 x 6 m footing, eB = 0.5 m, eL = 1 m -> 3 x 4 m', () => {
  const r = b.effectiveAreaRect({ B: 4, L: 6, eB: 0.5, eL: 1 });
  assert.deepEqual([r.Beff, r.Leff, r.Aeff], [3, 4, 12]);
});

test('circular footing effective area (API RP 2GEO / DNV): R = 5 m, e = 1 m hand calculation', () => {
  const r = b.effectiveAreaCircle({ R: 5, e: 1 });
  // A' = 2 [25 acos(0.2) - sqrt(24)] = 2 [34.23595 - 4.89898] = 58.6739 m2
  close(r.Aeff, 58.6739, 1e-5);
  // L'/B' = sqrt((R + e)/(R - e)) = sqrt(1.5); B' L' = A'
  close(r.Leff / r.Beff, Math.sqrt(1.5), 1e-12);
  close(r.Leff * r.Beff, r.Aeff, 1e-12);
  close(r.Leff, 8.4771, 1e-4);
  close(r.Beff, 6.9215, 1e-4);
  // no eccentricity: full circle
  close(b.effectiveAreaCircle({ R: 5, e: 0 }).Aeff, Math.PI * 25, 1e-12);
});

test('Vesic general equation, strip on sand, hand calculation', () => {
  // B = 2 m, D = 1 m, gamma = 18 kN/m3, phi = 30 deg: q = 18 kPa, dq = 1 + 2 tan30 (1 - sin30)^2 (D/B) = 1.14434
  const r = b.bearingCapacity({ method: 'vesic', B: 2, D: 1, phi: 30 * deg, gamma: 18e3 });
  const dq = 1 + 2 * Math.tan(30 * deg) * 0.25 * 0.5;
  close(dq, 1.14434, 1e-5);
  const expected = 18e3 * 18.401 * dq + 0.5 * 18e3 * 2 * 22.402;
  close(r.qult, expected, 2e-4);
  close(r.factors.sq, 1, 1e-12);
});

test('Hansen undrained: square surface footing gives 6.17 su; full sliding load halves the strip capacity to (1 + pi/2) su (Green exact)', () => {
  const su = 40e3;
  const sq = b.bearingCapacity({ method: 'hansen', B: 3, L: 3, c: su, gamma: 6e3 });
  close(sq.qult, (2 + Math.PI) * su * 1.2, 1e-12);
  const A = 3;
  const st = b.bearingCapacity({ method: 'hansen', B: 3, c: su, gamma: 6e3, V: 1, HB: A * su });
  close(st.qult, (1 + Math.PI / 2) * su, 1e-12);
});

test('Meyerhof, Hansen and Vesic agree on a vertically loaded surface strip on clay: 5.14 su', () => {
  for (const method of ['meyerhof', 'hansen', 'vesic']) {
    close(b.bearingCapacity({ method, B: 2, c: 30e3, gamma: 7e3 }).qult, (2 + Math.PI) * 30e3, 1e-12, method);
  }
});

test('inclination, base tilt and ground slope reduce capacity monotonically (Hansen and Vesic)', () => {
  for (const method of ['hansen', 'vesic']) {
    const base = { method, B: 3, L: 6, D: 1, phi: 32 * deg, c: 5e3, gamma: 10e3, V: 5e6 };
    const q0 = b.bearingCapacity(base).qult;
    const q1 = b.bearingCapacity({ ...base, HB: 5e5 }).qult;
    const q2 = b.bearingCapacity({ ...base, HB: 1e6 }).qult;
    const q3 = b.bearingCapacity({ ...base, baseTilt: 5 * deg }).qult;
    const q4 = b.bearingCapacity({ ...base, groundSlope: 10 * deg }).qult;
    assert.ok(q0 > q1 && q1 > q2 && q0 > q3 && q0 > q4, method);
  }
});

test('Vesic inclination exponent m by hand: B/L = 0.5, load along B -> m = 2.5/1.5', () => {
  const r = b.bearingCapacity({ method: 'vesic', B: 2, L: 4, phi: 30 * deg, gamma: 10e3, V: 1e6, HB: 2e5 });
  close(r.factors.m, 2.5 / 1.5, 1e-12);
  close(r.factors.iq, 0.8 ** (2.5 / 1.5), 1e-12);
  close(r.factors.ig, 0.8 ** (2.5 / 1.5 + 1), 1e-12);
});

test('API drained capacity equals the Vesic general equation for a central vertical load on the surface', () => {
  const phi = 35 * deg;
  const a = b.apiDrainedCapacity({ Beff: 4, Leff: 8, phi, gamma: 9e3, c: 2e3 });
  const v = b.bearingCapacity({ method: 'vesic', B: 4, L: 8, phi, gamma: 9e3, c: 2e3 });
  close(a.q, v.qult, 1e-12);
  assert.ok(a.Kq > 1 && a.Kgamma < 1 && a.Kc > 1);
});

test('ISO 19901-4 undrained: surface strip 5.14 su B; depth and inclination terms by hand', () => {
  const su = 20e3;
  close(b.isoUndrainedCapacity({ Beff: 5, Leff: Infinity, su0: su }).Qd, (2 + Math.PI) * su * 5, 1e-12);
  // 5 x 10 m, D = 1 m, H = 0.36 A su: ic = 0.5 - 0.5 sqrt(0.64) = 0.1; sc = 0.18 (1 - 0.2) 0.5 = 0.072; dc = 0.3 atan(0.2)
  const r = b.isoUndrainedCapacity({ Beff: 5, Leff: 10, su0: su, D: 1, H: 0.36 * 50 * su });
  close(r.ic, 0.1, 1e-12);
  close(r.sc, 0.072, 1e-12);
  close(r.dc, 0.3 * Math.atan(0.2), 1e-12);
  close(r.Kc, 1 + 0.072 + 0.3 * Math.atan(0.2) - 0.1, 1e-12);
  close(b.isoScv(3), -0.025, 1e-12);
});

test('Davis & Booker form: uniform strength recovers 5.14 su; gradient term F (5.14 su0 + kappa B / 4)', () => {
  close(b.davisBookerStrip({ su0: 10e3, kappa: 0, B: 4 }).qult, (2 + Math.PI) * 10e3, 1e-12);
  close(b.davisBookerStrip({ su0: 10e3, kappa: 2e3, B: 4, F: 1.2 }).qult, 1.2 * ((2 + Math.PI) * 10e3 + 2e3), 1e-12);
});

test('layered soils: punch-through and squeezing hand calculations', () => {
  // Meyerhof & Hanna strip: qb + gamma H^2 (1 + 2D/H) Ks tan(phi) / B - gamma H
  const r = b.punchThroughMeyerhofHanna({ B: 2, H: 1.5, D: 1, gamma1: 10e3, phi1: 35 * deg, Ks: 4, qb: 150e3 });
  close(r.qult, 150e3 + (10e3 * 2.25 / 2) * (1 + 2 / 1.5) * 4 * Math.tan(35 * deg) - 15e3, 1e-12);
  assert.equal(r.governs, 'punch-through');
  // load spread 1:3 below a 10 m spudcan with 6 m of sand: equivalent diameter 14 m
  const ls = b.punchThroughLoadSpread({ B: 10, H: 6, ns: 3, su: 30e3, p0: 50e3 });
  close(ls.Bequiv, 14, 1e-12);
  close(ls.Qv, (Math.PI / 4) * 196 * (1.2 * (2 + Math.PI) * 30e3 + 50e3), 1e-12);
  // squeezing: B = 12, T = 2 -> (5 + 0.33·6) su = 6.98 su > 6 su
  const sq = b.squeezingCapacity({ B: 12, T: 2, su: 20e3 });
  assert.ok(sq.squeezing);
  close(sq.qult, 6.98 * 20e3, 1e-12);
  assert.equal(b.squeezingCapacity({ B: 4, T: 2, su: 20e3 }).squeezing, false);
  const cc = b.punchThroughClayOverClay({ B: 10, H: 2, suTop: 80e3, suBottom: 20e3 });
  close(cc.qPunch, 3 * 0.2 * 80e3 + 6 * (1 + 0.2 * 0.2) * 20e3, 1e-12);
  assert.equal(cc.governs, 'punch-through');
});

test('sliding resistance: undrained su A plus skirt passive; drained V tan(delta) plus Rankine wedge', () => {
  close(b.slidingResistance({ mode: 'undrained', A: 20, suBase: 10e3 }).H, 200e3, 1e-12);
  const u = b.slidingResistance({ mode: 'undrained', A: 20, suBase: 10e3, skirtDepth: 0.5, skirtWidth: 4, suSkirt: 8e3, gamma: 6e3 });
  close(u.Hpassive, (2 * 8e3 + 0.5 * 6e3 * 0.5) * 0.5 * 4, 1e-12);
  assert.equal(u.Hactive, 0);
  const d = b.slidingResistance({ mode: 'drained', A: 20, V: 1e6, delta: 25 * deg, skirtDepth: 1, skirtWidth: 5, phi: 30 * deg, gamma: 10e3 });
  close(d.Hbase, 1e6 * Math.tan(25 * deg), 1e-12);
  close(d.Hpassive - d.Hactive, 0.5 * 10e3 * 5 * (3 - 1 / 3), 1e-9);
});

test('overturning and design-format checks', () => {
  const o = b.overturningCheck({ V: 1000, M: 500, B: 6 });
  close(o.e, 0.5, 1e-12);
  close(o.fsOverturning, 6, 1e-12);
  assert.equal(o.fullContact, true);
  assert.equal(b.overturningCheck({ V: 1000, M: 1500, B: 6 }).fullContact, false);
  const c = b.designCheck({ load: 100, resistance: 300, loadFactor: 1.35, materialFactor: 1.25 });
  close(c.uc, 1.35 * 100 / (300 / 1.25), 1e-12);
});

test('invalid input throws RangeError naming the argument', () => {
  assert.throws(() => b.bearingNq(-0.1), /phi/);
  assert.throws(() => b.bearingCapacity({ B: -1, gamma: 1 }), /B/);
  assert.throws(() => b.effectiveAreaCircle({ R: 1, e: 2 }), RangeError);
  assert.throws(() => b.bearingCapacity({ method: 'x', B: 1, gamma: 1 }), /method/);
});
