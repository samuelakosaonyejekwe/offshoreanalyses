import test from 'node:test';
import assert from 'node:assert/strict';
import * as e from '../../src/engine/geotech/earthpressure.js';

const deg = Math.PI / 180;
const close = (a, x, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - x) <= rel * Math.max(Math.abs(x), 1e-300), `${msg} got ${a}, expected ${x}`);

test('Rankine at 30 deg: Ka = 1/3, Kp = 3', () => {
  const r = e.rankine(30 * deg);
  close(r.Ka, 1 / 3, 1e-12);
  close(r.Kp, 3, 1e-12);
});

test('Coulomb equals Rankine for a smooth vertical wall with level backfill, and for delta = beta on a sloping backfill', () => {
  for (const phi of [20, 30, 40]) {
    close(e.coulombKa({ phi: phi * deg }), e.rankine(phi * deg).Ka, 1e-12);
    close(e.coulombKp({ phi: phi * deg }), e.rankine(phi * deg).Kp, 1e-12);
  }
  close(e.coulombKa({ phi: 30 * deg, delta: 10 * deg, beta: 10 * deg }), e.rankine(30 * deg, 10 * deg).Ka, 1e-12);
  close(e.coulombKa({ phi: 35 * deg, delta: 20 * deg, beta: 20 * deg }), e.rankine(35 * deg, 20 * deg).Ka, 1e-12);
});

test('Coulomb tabulated values: phi = 30, delta = 20 deg -> Ka = 0.297, Kp = 6.105 (e.g. Das, Bowles tables)', () => {
  close(e.coulombKa({ phi: 30 * deg, delta: 20 * deg }), 0.2973, 3e-4);
  close(e.coulombKp({ phi: 30 * deg, delta: 20 * deg }), 6.105, 3e-4);
});

test('at-rest coefficient: Jaky 0.5 at 30 deg; OCR = 4 doubles it (Mayne & Kulhawy)', () => {
  close(e.k0(30 * deg), 0.5, 1e-12);
  close(e.k0(30 * deg, 4), 1.0, 1e-12);
});

test('Mononobe-Okabe reduces to Coulomb at kh = 0 (active and passive, with batter, slope and wall friction)', () => {
  const g = { phi: 34 * deg, delta: 17 * deg, beta: 8 * deg, theta: 5 * deg };
  const mo = e.mononobeOkabe({ ...g, kh: 0 });
  close(mo.KAE, e.coulombKa(g), 1e-12);
  close(mo.KPE, e.coulombKp(g), 1e-12);
});

test('Mononobe-Okabe: seismic increment close to the Seed & Whitman (1970) approximation 0.75 kh; limit of existence', () => {
  const g = { phi: 35 * deg, delta: 17.5 * deg };
  const ka = e.coulombKa(g);
  const kae = e.mononobeOkabe({ ...g, kh: 0.2 }).KAE;
  assert.ok(kae > ka);
  assert.ok(Math.abs((kae - ka) - 0.15) / 0.15 < 0.15, `increment ${kae - ka}`);
  assert.ok(e.mononobeOkabe({ ...g, kh: 0.2 }).KPE < e.coulombKp(g));
  assert.throws(() => e.mononobeOkabe({ phi: 30 * deg, kh: 0.7 }), /kh/);
  // vertical smooth wall, horizontal backfill, psi -> phi: KAE -> cos(0)^2 / (cos^2... ) = 1/cos(psi)·1 limit check
  const lim = e.mononobeOkabe({ phi: 30 * deg, kh: Math.tan(30 * deg) }).KAE;
  close(lim, 1 / (Math.cos(30 * deg) * Math.cos(30 * deg)), 1e-9);
});

test('Westergaard: resultant (7/12) kh gw h^2 equals the integral of the pressure, acting 0.4 h above the base', () => {
  const kh = 0.15, h = 12, gw = 9810;
  const n = 20000;
  let F = 0, M = 0;
  for (let i = 0; i < n; i++) {
    const z = (i + 0.5) * h / n;
    const p = e.westergaard({ kh, h, z, gammaW: gw }).p * h / n;
    F += p; M += p * (h - z);
  }
  const r = e.westergaard({ kh, h, gammaW: gw });
  close(r.force, F, 1e-6);
  close(r.armAboveBase, M / F, 1e-5);
});

test('pressure diagram: dry sand 0.5 Ka gamma H^2 at H/3; surcharge; water table (hand calculation)', () => {
  const dry = e.pressureDiagram({ layers: [{ thickness: 6, gamma: 18e3, phi: 30 * deg }] });
  close(dry.P, 0.5 * (1 / 3) * 18e3 * 36, 1e-12);
  close(dry.armAboveBase, 2, 1e-12);
  const q = e.pressureDiagram({ layers: [{ thickness: 6, gamma: 18e3, phi: 30 * deg }], surcharge: 12e3 });
  close(q.P, 108e3 + 12e3 / 3 * 6, 1e-12);
  close(q.armAboveBase, (108e3 * 2 + 24e3 * 3) / 132e3, 1e-12);
  // water table at 2 m: sigma_v' = 36 kPa at 2 m, 36 + 4 (20 - 9.81) = 76.76 kPa at 6 m; water 0.5·9.81·16
  const w = e.pressureDiagram({ layers: [{ thickness: 6, gamma: 18e3, gammaSat: 20e3, phi: 30 * deg }], waterTableDepth: 2 });
  const sv2 = 36e3, sv6 = 36e3 + 4 * (20e3 - 9810);
  const Pe = 0.5 * (sv2 / 3) * 2 + 0.5 * (sv2 / 3 + sv6 / 3) * 4;
  close(w.Peff, Pe, 1e-12);
  close(w.Pwater, 0.5 * 9810 * 16, 1e-12);
  close(w.P, Pe + 0.5 * 9810 * 16, 1e-12);
});

test('pressure diagram with cohesion: tension crack depth 2c/(gamma sqrt(Ka)); thrust on the cracked section; passive', () => {
  const c = 10e3, gamma = 18e3, H = 8, Ka = 1 / 3;
  const zc = e.tensionCrackDepth({ c, gamma, Ka });
  close(zc, 2 * c / (gamma * Math.sqrt(Ka)), 1e-12);
  const r = e.pressureDiagram({ layers: [{ thickness: H, gamma, phi: 30 * deg, c }] });
  const pBase = Ka * gamma * H - 2 * c * Math.sqrt(Ka);
  close(r.P, 0.5 * pBase * (H - zc), 1e-12);
  close(r.armAboveBase, (H - zc) / 3, 1e-12);
  const noCrack = e.pressureDiagram({ layers: [{ thickness: H, gamma, phi: 30 * deg, c }], tensionCrack: false });
  close(noCrack.P, 0.5 * Ka * gamma * H * H - 2 * c * Math.sqrt(Ka) * H, 1e-12);
  const pas = e.pressureDiagram({ layers: [{ thickness: 3, gamma, phi: 30 * deg, c }], state: 'passive' });
  close(pas.P, 0.5 * 3 * gamma * 9 + 2 * c * Math.sqrt(3) * 3, 1e-12);
});

test('layered pressure diagram: stress jump at an interface between sand and clay layers', () => {
  const r = e.pressureDiagram({ layers: [{ thickness: 3, gamma: 17e3, phi: 32 * deg }, { thickness: 4, gamma: 19e3, phi: 25 * deg }] });
  const k1 = e.rankine(32 * deg).Ka, k2 = e.rankine(25 * deg).Ka;
  close(r.points[1].sigmaH, k1 * 51e3, 1e-12);
  close(r.points[2].sigmaH, k2 * 51e3, 1e-12);
  close(r.P, 0.5 * k1 * 51e3 * 3 + 0.5 * (k2 * 51e3 + k2 * 127e3) * 4, 1e-12);
});

test('strip surcharge covering the whole surface gives sigma_h = q; line load formulas by hand', () => {
  close(e.stripSurchargePressure({ q: 10e3, a: 0, b: 1e9, z: 3 }), 10e3, 1e-6);
  // m = 0.5, n = 0.5: (4Q/(pi H)) m^2 n/(m^2+n^2)^2 = (4Q/(pi H)) 0.125/0.25
  close(e.lineSurchargePressure({ Q: 50e3, x: 3, z: 3, H: 6 }), (4 * 50e3 / (Math.PI * 6)) * 0.5, 1e-12);
  close(e.lineSurchargePressure({ Q: 50e3, x: 1.2, z: 3, H: 6 }), (50e3 / 6) * 0.203 * 0.5 / (0.16 + 0.25) ** 2, 1e-12);
});

test('gravity wall stability by hand: sliding, overturning, eccentricity and base pressures', () => {
  // 4 m wide block, W = 600 kN/m at 2 m; earth thrust 150 kN/m at 2.5 m; uplift 80 kN/m at mid-base
  const r = e.gravityWallStability({ B: 4, vertical: [{ V: 600e3, x: 2 }], horizontal: [{ H: 150e3, y: 2.5 }], U: 80e3, deltaBase: 30 * deg, qult: 800e3 });
  close(r.N, 520e3, 1e-12);
  close(r.fsSliding, 520e3 * Math.tan(30 * deg) / 150e3, 1e-12);
  close(r.fsOverturning, 1200e3 / (375e3 + 160e3), 1e-12);
  const ecc = 2 - (1200e3 - 535e3) / 520e3;
  close(r.e, ecc, 1e-12);
  assert.equal(r.inMiddleThird, false);
  close(r.qMax, 2 * 520e3 / (3 * (2 - ecc)), 1e-12);
  close(r.qMeyerhof, 520e3 / (4 - 2 * ecc), 1e-12);
  close(r.fsBearing, 800e3 / r.qMeyerhof, 1e-12);
  const b = e.bollardPull({ pull: 1000e3, angle: 30 * deg, spreadLength: 20, height: 9 });
  close(b.H, 1000e3 * Math.cos(30 * deg) / 20, 1e-12);
  close(b.M, b.H * 9, 1e-12);
});

test('cantilever sheet pile in dry sand: closed form d0 = H / ((Kp/Ka)^(1/3) - 1) and maximum moment', () => {
  const H = 5, g = 18e3, Ka = 1 / 3, Kp = 3;
  const r = e.sheetPileWall({ H, gamma: g, Ka, Kp });
  const d0 = H / (Math.cbrt(Kp / Ka) - 1);
  close(r.d0, d0, 1e-8);
  close(r.d, 1.2 * d0, 1e-8);
  assert.equal(r.anchorForce, 0);
  // zero shear at x below dredge level: Ka (H + x)^2 = Kp x^2; Mmax = gamma [Ka (H + x)^3 - Kp x^3] / 6
  const x = H / (Math.sqrt(Kp / Ka) - 1);
  close(r.Mmax, g * (Ka * (H + x) ** 3 - Kp * x ** 3) / 6, 1e-5);
  close(r.zMmax, H + x, 2e-3);
});

test('anchored sheet pile (free-earth support): moment and force equilibrium satisfied; factor on Kp deepens the wall; Rowe factor scales the moment', () => {
  const H = 8, g = 18e3, Ka = 1 / 3, Kp = 3, a = 1;
  const r = e.sheetPileWall({ H, gamma: g, Ka, Kp, anchorDepth: a });
  const L = H + r.d0;
  const Pa = 0.5 * Ka * g * L * L;
  const Pp = 0.5 * Kp * g * r.d0 * r.d0;
  // moments about the anchor: Pa (2L/3 - a) = Pp (H + 2 d0/3 - a)
  close(Pa * (2 * L / 3 - a), Pp * (H + 2 * r.d0 / 3 - a), 1e-8);
  close(r.anchorForce, Pa - Pp, 1e-8);
  close(r.Pa, Pa, 1e-9);
  // maximum moment where the shear below the anchor vanishes: 0.5 Ka g z^2 = T  (z < H)
  const z = Math.sqrt(2 * r.anchorForce / (Ka * g));
  assert.ok(z < H);
  close(r.Mmax, r.anchorForce * (z - a) - Ka * g * z ** 3 / 6, 1e-4);
  const safer = e.sheetPileWall({ H, gamma: g, Ka, Kp, anchorDepth: a, fsPassive: 1.5 });
  assert.ok(safer.d0 > r.d0 && safer.anchorForce > r.anchorForce);
  close(e.sheetPileWall({ H, gamma: g, Ka, Kp, anchorDepth: a, roweFactor: 0.7 }).Mmax, 0.7 * r.Mmax, 1e-12);
});

test('sheet pile with water table and surcharge: horizontal equilibrium of the integrated diagrams', () => {
  const o = { H: 6, gamma: 18e3, gammaSub: 10e3, waterDepth: 2, Ka: 0.3, Kp: 3.2, surcharge: 10e3, anchorDepth: 1.5 };
  const r = e.sheetPileWall(o);
  // active thrust by hand: surcharge + dry part + submerged part
  const L = o.H + r.d0;
  const Pa = 0.3 * (10e3 * L + 0.5 * 18e3 * 4 + 36e3 * (L - 2) + 0.5 * 10e3 * (L - 2) ** 2);
  close(r.Pa, Pa, 1e-8);
  close(r.Pp, 0.5 * 3.2 * 10e3 * r.d0 ** 2, 1e-8);
  close(r.anchorForce, r.Pa - r.Pp, 1e-12);
});

test('deadman anchor: 0.5 gamma H^2 (Kp - Ka) per metre; block with shape factor', () => {
  close(e.deadmanAnchor({ gamma: 18e3, H: 2, Ka: 1 / 3, Kp: 3 }).Tult, 0.5 * 18e3 * 4 * (3 - 1 / 3), 1e-12);
  close(e.deadmanAnchor({ gamma: 18e3, H: 2, Ka: 1 / 3, Kp: 3, b: 1.5, shapeFactor: 1.8, safetyFactor: 2 }).Tallow, 0.5 * 18e3 * 4 * (3 - 1 / 3) * 1.5 * 1.8 / 2, 1e-12);
});
