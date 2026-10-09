import test from 'node:test';
import assert from 'node:assert/strict';
import * as a from '../../src/engine/geotech/anchors.js';

const close = (x, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(x - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${x}, expected ${e}`);

test('drag anchor power law UHC = a W^b: 15 t anchor with a = 229, b = 0.92', () => {
  const r = a.dragAnchorUhc({ mass: 15000, a: 229, b: 0.92 });
  close(r.UHC, 229e3 * 15 ** 0.92, 1e-12);
  close(r.efficiency, r.UHC / (15000 * 9.80665), 1e-12);
  assert.ok(Object.keys(a.DRAG_ANCHOR_INDICATIVE).length >= 3);
});

test('Neubecker & Randolph (1995) embedded chain: Ta theta_a^2 / 2 = za Qbar and T/Ta = exp(mu theta_a) by hand', () => {
  // chain d = 0.1 m, En = 2.5, Nc = 7.6, su = 1.5 z kPa, padeye at 10 m, Ta = 1 MN:
  // integral of Q = 2.5·0.1·7.6·1.5e3·10^2/2 = 142.5 kN  ->  theta_a = sqrt(2·142.5/1000) = 0.5339 rad
  const r = a.embeddedChain({ Ta: 1e6, za: 10, d: 0.1, su0: 0, k: 1.5e3, mu: 0.4 });
  close(r.thetaA, 0.5339, 1e-4);
  close(1e6 * r.thetaA ** 2 / 2, 10 * r.Qbar, 1e-12);
  close(r.Tm, 1e6 * Math.exp(0.4 * r.thetaA), 1e-12);
  close(r.Tm, 1.2381e6, 1e-4);
  // non-zero mudline angle
  const m = a.embeddedChain({ Ta: 1e6, za: 10, d: 0.1, su0: 0, k: 1.5e3, mu: 0.4, thetaM: 0.2 });
  close(0.5e6 * (m.thetaA ** 2 - 0.04), 142.5e3, 1e-12);
  // inverse problem
  const inv = a.embeddedChainFromMudline({ Tm: r.Tm, za: 10, d: 0.1, su0: 0, k: 1.5e3, mu: 0.4 });
  close(inv.Ta, 1e6, 1e-8);
  close(inv.thetaA, r.thetaA, 1e-8);
});

test('drag anchor trajectory: ultimate depth z = thetaW^2 Ne Af / (En Nc d) for strength proportional to depth', () => {
  const o = { Af: 10, Ne: 12, thetaW: 0.6, d: 0.1, su0: 0, k: 1500 };
  const r = a.dragAnchorTrajectory(o);
  close(r.zUlt, 0.36 * 12 * 10 / (2.5 * 7.6 * 0.1), 1e-9);
  close(r.TaUlt, 12 * 10 * 1500 * r.zUlt, 1e-9);
  close(r.TmUlt, r.TaUlt * Math.exp(0.4 * 0.6), 1e-12);
  // the chain angle grows as thetaW sqrt(z/zUlt); the trajectory is monotonic
  for (let i = 1; i < r.x.length; i++) assert.ok(r.x[i] > r.x[i - 1] && r.z[i] > r.z[i - 1]);
  close(r.thetaA[r.thetaA.length - 1], 0.6 * Math.sqrt(0.95), 1e-9);
  // a thicker chain cuts less deep
  assert.ok(a.dragAnchorTrajectory({ ...o, d: 0.15 }).zUlt < r.zUlt);
  // drag distance converges with step refinement
  close(a.dragAnchorTrajectory({ ...o, n: 3200 }).dragDistance, a.dragAnchorTrajectory({ ...o, n: 1600 }).dragDistance, 1e-5);
});

test('plate anchor: Nc su A with keying loss; deep strip factor 3 pi + 2 = 11.42 (Rowe & Davis 1982)', () => {
  const r = a.plateAnchorCapacity({ A: 30, B: 3, zInstall: 25, keyingLoss: 1.5, su0: 2e3, k: 1.5e3 });
  close(r.Nc, 11.42, 1e-3);
  close(r.su, 2e3 + 1.5e3 * 23.5, 1e-12);
  close(r.T, (3 * Math.PI + 2) * r.su * 30, 1e-12);
  const sh = a.plateAnchorCapacity({ A: 30, B: 3, zInstall: 3, su0: 10e3, NcOfDepth: (x) => 4 + 3 * x });
  close(sh.Nc, 7, 1e-12);
});

test('anchor pile: clay translation at the optimal padeye and axial pull-out by hand; sand closed forms', () => {
  const D = 2, L = 20;
  const c = a.anchorPileCapacity({ soil: 'clay', D, L, zPadeye: 40 / 3, su0: 0, k: 2e3, NpConst: 9, alpha: 0.6, W: 200e3 });
  close(c.zOptimal, 40 / 3, 1e-9);
  close(c.Hult, 9 * 2e3 * D * L * L / 2, 1e-6);
  close(c.Vult, 0.6 * 0.5 * 2e3 * 400 * Math.PI * D + 200e3, 1e-12);
  // sand: Broms reaction p = 3 Kp gamma z D; load at the mudline: H = c (z0^2 - L^2/2) with z0^3 = L^3/2
  const s = a.anchorPileCapacity({ soil: 'sand', D, L, zPadeye: 0, phi: Math.PI / 6, gamma: 10e3, KtanDelta: 0.4 });
  const cc = 3 * 3 * 10e3 * D;
  const z0 = L / Math.cbrt(2);
  close(s.z0, z0, 1e-4);
  close(s.Hult, cc * (z0 * z0 - L * L / 2), 1e-3);
  close(s.Vult, 0.5 * 0.4 * 10e3 * 400 * Math.PI * D, 1e-12);
  // inclined load lies on the interaction envelope
  const i = a.anchorPileCapacity({ soil: 'clay', D, L, zPadeye: 40 / 3, su0: 0, k: 2e3, NpConst: 9, alpha: 0.6, angle: 0.6 });
  close((i.Tult * Math.cos(0.6) / i.Hult) ** (L / D + 0.5) + (i.Tult * Math.sin(0.6) / i.Vult) ** (L / (3 * D) + 4.5), 1, 1e-9);
});

test('gravity anchor statics: drained Tult = W tan(d) / (cos t + sin t tan(d)); undrained min of sliding and uplift', () => {
  const W = 2e6, d = 0.5, t = 0.3;
  const r = a.gravityAnchorCapacity({ W, angle: t, mode: 'drained', delta: d, T: 5e5 });
  close(r.Tult, W * Math.tan(d) / (Math.cos(t) + Math.sin(t) * Math.tan(d)), 1e-12);
  // at Tult the friction demand equals capacity
  close(r.Tult * Math.cos(t), (W - r.Tult * Math.sin(t)) * Math.tan(d), 1e-12);
  close(r.fs, r.Tult / 5e5, 1e-12);
  const u = a.gravityAnchorCapacity({ W, angle: t, mode: 'undrained', su: 10e3, A: 50, Hkey: 1e5 });
  close(u.Tult, Math.min(6e5 / Math.cos(t), W / Math.sin(t)), 1e-12);
});

test('torpedo anchor penetration: ODE integration matches the energy balance; partially embedded case solved by hand', () => {
  // uniform su, tip resistance only (alpha = 0), stops within one anchor length:
  // 0.5 M v0^2 + W z = Nc su At z  ->  z = 0.5 M v0^2 / (Nc su At - W)
  const o = { mass: 20e3, W: 150e3, D: 1, Lanchor: 15, v0: 10, su0: 100e3, Nc: 12, alpha: 0 };
  const At = Math.PI / 4;
  const z = 0.5 * 20e3 * 100 / (12 * 100e3 * At - 150e3);
  close(a.torpedoPenetrationEnergy(o), z, 1e-8);
  close(a.torpedoPenetration({ ...o, dt: 1e-4 }).zTip, z, 1e-6);
  // realistic case with shaft friction and strength gradient: the two methods agree
  const p = { mass: 80e3, W: 80e3 * 9.81 * 0.87, D: 1, Lanchor: 15, v0: 25, su0: 5e3, k: 1.5e3 };
  const ode = a.torpedoPenetration(p);
  close(ode.zTip, a.torpedoPenetrationEnergy(p), 1e-6);
  assert.ok(ode.zTip > 15);
  // strain-rate enhancement and inertial drag reduce penetration
  assert.ok(a.torpedoPenetration({ ...p, rateModel: 'log', lambda: 0.15 }).zTip < ode.zTip);
  assert.ok(a.torpedoPenetration({ ...p, Cd: 0.3 }).zTip < ode.zTip);
  // a power-law rate model with zero exponent is rate independent
  close(a.torpedoPenetration({ ...p, rateModel: 'power', beta: 0 }).zTip, ode.zTip, 1e-12);
});

test('torpedo anchor pull-out capacity by hand', () => {
  const r = a.torpedoCapacity({ zTip: 30, Lanchor: 15, D: 1, su0: 5e3, k: 1.5e3, alpha: 0.6, Nc: 9, W: 500e3 });
  close(r.friction, 0.6 * (5e3 + 1.5e3 * 22.5) * Math.PI * 15, 1e-12);
  close(r.endBearing, 9 * (5e3 + 1.5e3 * 15) * Math.PI / 4, 1e-12);
  close(r.T, 500e3 + r.friction + r.endBearing, 1e-12);
});

test('design check with set-up and cyclic factors, proof load, and code factor tables', () => {
  const r = a.anchorDesignCheck({ staticCapacity: 8e6, setUpFactor: 1.2, cyclicFactor: 0.9, tension: 5e6, safetyFactor: 1.5, proofLoad: 4e6 });
  close(r.resistance, 8.64e6, 1e-12);
  close(r.fs, 1.728, 1e-12);
  assert.equal(r.ok, true);
  assert.equal(r.proofOk, false);
  assert.equal(a.ANCHOR_SAFETY_FACTORS['API RP 2SK'].dragAnchor.permanent.intact, 1.5);
  const f = a.ANCHOR_SAFETY_FACTORS['DNV-OS-E301'].uls.class2;
  close(a.dnvDesignTension({ Tmean: 2e6, Tdyn: 1e6, ...f }), 1.4 * 2e6 + 2.1 * 1e6, 1e-12);
  assert.throws(() => a.dragAnchorUhc({ mass: -1, a: 1, b: 1 }), /mass/);
});
