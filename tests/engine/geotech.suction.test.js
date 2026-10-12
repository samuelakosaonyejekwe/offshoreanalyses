import test from 'node:test';
import assert from 'node:assert/strict';
import * as s from '../../src/engine/geotech/suction.js';

const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);

test('Houlsby & Byrne (2005a) clay installation: force balance by hand at h = 5 m', () => {
  // Do = 5 m, t = 30 mm, su = 10 + 2 z kPa: integral of su over 5 m = 75 kN/m, su at tip = 20 kPa
  const o = { Do: 5, t: 0.03, h: 5, su: { su0: 10e3, k: 2e3 }, alphaO: 0.5, alphaI: 0.4, Nc: 9, gamma: 6e3, Vprime: 300e3 };
  const r = s.suctionInstallClay(o);
  const out = 0.5 * 75e3 * Math.PI * 5;
  const inn = 0.4 * 75e3 * Math.PI * 4.94;
  const tip = (6e3 * 5 + 9 * 20e3) * Math.PI * 4.97 * 0.03;
  close(r.outside, out, 1e-9);
  close(r.inside, inn, 1e-9);
  close(r.tip, tip, 1e-12);
  // V' + s pi Do^2/4 = R
  close(o.Vprime + r.s * Math.PI * 25 / 4, out + inn + tip, 1e-9);
  // plug heave limit: Nc* su_tip + 4 alpha_i (integral su) / Di
  close(r.sAllow, 9 * 20e3 + 4 * 0.4 * 75e3 / 4.94, 1e-9);
  // with a finite outside stress zone the allowable suction drops
  assert.ok(s.suctionInstallClay({ ...o, m: 1.5 }).sAllow < r.sAllow);
});

test('clay self-weight penetration closes the force balance; layered profile integrates by layers', () => {
  const o = { Do: 4, t: 0.025, su: { su0: 5e3, k: 1.5e3 }, gamma: 5e3, Vprime: 250e3, L: 10 };
  const h = s.selfWeightPenetrationClay(o);
  close(s.suctionInstallClay({ ...o, h }).R, 250e3, 1e-8);
  assert.equal(s.suctionInstallClay({ ...o, h: h * 0.9 }).s, 0);
  const layered = s.suctionInstallClay({ Do: 4, t: 0.025, h: 6, su: [{ zTop: 0, su0: 10e3 }, { zTop: 4, su0: 30e3 }], alphaO: 1, alphaI: 0, Nc: 0 });
  close(layered.outside, (10e3 * 4 + 30e3 * 2) * Math.PI * 4, 1e-6);
  assert.equal(s.strengthProfile([{ zTop: 0, su0: 10e3 }, { zTop: 4, su0: 30e3, k: 1e3 }])(5), 31e3);
});

test('pore-pressure factor a (Houlsby & Byrne 2005b): 0.45 at h/D = 0, 0.09 at depth; permeability ratio by hand', () => {
  close(s.seepageFactorA(0), 0.45, 1e-14);
  close(s.seepageFactorA(100), 0.09, 1e-9);
  const a1 = 0.45 - 0.36 * (1 - Math.exp(-0.5 / 0.48));
  close(s.seepageFactorA(0.5), a1, 1e-14);
  close(s.seepageFactorA(0.5, 3), 3 * a1 / (1 - a1 + 3 * a1), 1e-14);
});

test('Houlsby & Byrne (2005b) sand installation: self-weight terms by hand and force balance with suction at h = 3 m', () => {
  const Do = 6, t = 0.03, h = 3, g = 9e3, Kt = 0.5, Nq = 50, Ng = 45, m = 1.5, V = 500e3;
  const Di = Do - 2 * t, D = Do - t;
  const Zo = Do * (m * m - 1) / (4 * Kt), Zi = Di / (4 * Kt);
  const Fo = Zo * Zo * (Math.exp(h / Zo) - 1 - h / Zo) * Kt * Math.PI * Do;
  const Fi = Zi * Zi * (Math.exp(h / Zi) - 1 - h / Zi) * Kt * Math.PI * Di;
  const Fq = Zi * (Math.exp(h / Zi) - 1) * Nq * Math.PI * D * t;
  const Fg = 0.5 * g * t * Ng * Math.PI * D * t;
  const r = s.suctionInstallSand({ Do, t, h, gamma: g, KtanDeltaO: Kt, KtanDeltaI: Kt, Nq, Ngamma: Ng, m, Vprime: V });
  close(r.Zo, 3.75, 1e-12);
  close(r.R0, g * (Fo + Fi + Fq) + Fg, 1e-10);
  assert.ok(r.R0 > V);
  // balance with seepage: V' + s pi Di^2/4 = (g + a s/h) Fo + (g - (1 - a) s/h)(Fi + Fq) + Fg
  const a = r.a;
  close(a, s.seepageFactorA(0.5), 1e-14);
  const lhs = V + r.s * Math.PI * Di * Di / 4;
  const rhs = (g + a * r.s / h) * Fo + (g - (1 - a) * r.s / h) * (Fi + Fq) + Fg;
  close(lhs, rhs, 1e-10);
  // seepage makes installation possible at far lower suction than without flow (a s/h terms removed)
  const noSeep = (r.R0 - V) / (Math.PI * Di * Di / 4);
  assert.ok(r.s < noSeep);
  // critical suction for piping
  close(r.sCrit, g * h / (1 - a), 1e-12);
  // friction-only sanity: without stress enhancement (h << Z) the outside friction tends to 0.5 g h^2 Ktan(delta) pi Do
  const shallow = s.suctionInstallSand({ Do, t, h: 0.05, gamma: g, KtanDeltaO: Kt, KtanDeltaI: Kt, Nq: 0, Ngamma: 0, m });
  close(shallow.outside, 0.5 * g * 0.0025 * Kt * Math.PI * Do, 1e-2);
});

test('sand self-weight penetration and suction-penetration curve', () => {
  const p = { Do: 6, t: 0.03, gamma: 9e3, KtanDeltaO: 0.5, KtanDeltaI: 0.5, Nq: 50, Ngamma: 45, Vprime: 500e3 };
  const h = s.selfWeightPenetrationSand({ ...p, L: 6 });
  close(s.suctionInstallSand({ ...p, h }).R0, 500e3, 1e-8);
  const c = s.suctionPenetrationCurve({ soil: 'sand', L: 6, n: 30, params: p });
  close(c.selfWeight, h, 1e-12);
  assert.equal(c.h.length, 30);
  assert.ok(c.s[29] > c.s[15] && c.s[0] === 0);
});

test('pump flow (displacement + Darcy seepage) and ring buckling pressure by hand', () => {
  const r = s.pumpFlow({ Di: 5, v: 0.0005, k: 1e-4, s: 30e3, h: 3, a: 0.2, gammaW: 10000 });
  const A = Math.PI * 25 / 4;
  close(r.gradient, 0.8 * 30e3 / (10000 * 3), 1e-12);
  close(r.Q, A * (0.0005 + 1e-4 * 0.8), 1e-12);
  close(s.ringBucklingPressure({ E: 210e9, nu: 0.3, t: 0.03, D: 6 }), 2 * 210e9 / 0.91 * (0.005) ** 3, 1e-12);
});

test('vertical capacity of a caisson in clay: sealed (reverse end bearing Nc = 9), vented and plug pull-out by hand', () => {
  const r = s.caissonVerticalCapacityClay({ Do: 5, t: 0.03, L: 15, su: { su0: 5e3, k: 1.5e3 }, alphaO: 0.65, alphaI: 0.65, W: 600e3, gamma: 6e3 });
  const I = 5e3 * 15 + 0.5 * 1.5e3 * 225;
  close(r.frictionOut, 0.65 * I * Math.PI * 5, 1e-9);
  close(r.endBearing, 9 * 27.5e3 * Math.PI * 25 / 4, 1e-12);
  close(r.reverseBearing, 600e3 + r.frictionOut + r.endBearing, 1e-12);
  close(r.vented, 600e3 + 0.65 * I * Math.PI * (5 + 4.94), 1e-9);
  close(r.plugPullout, 600e3 + r.frictionOut + 6e3 * 15 * Math.PI * 4.94 ** 2 / 4, 1e-12);
  assert.equal(r.sealed, Math.min(r.reverseBearing, r.plugPullout));
});

test('Murff & Hamilton Np: 2 at the mudline, 9 at depth; xi from the strength profile', () => {
  close(s.murffHamiltonNp({ z: 0, D: 2, su0: 10e3, k: 1e3 }), 2, 1e-14);
  close(s.murffHamiltonNp({ z: 500, D: 2, su0: 10e3, k: 1e3 }), 9, 1e-9);
  // lambda = su0/(k D) = 2 -> xi = 0.35
  close(s.murffHamiltonNp({ z: 4, D: 2, su0: 4e3, k: 1e3 }), 9 - 7 * Math.exp(-0.35 * 2), 1e-12);
  // uniform strength: xi = 0.55
  close(s.murffHamiltonNp({ z: 4, D: 2, su0: 4e3 }), 9 - 7 * Math.exp(-1.1), 1e-12);
});

test('rigid pile with uniform reaction: closed-form rotation depth and capacity (Broms-type statics)', () => {
  const p = 1000, L = 10, e = 5;
  const r = s.rigidPileLateralCapacity({ p: () => p, L, zLoad: -e, n: 2000 });
  const z0 = -e + Math.sqrt(e * e + e * L + L * L / 2);
  close(r.z0, z0, 1e-6);
  close(r.H, p * (2 * z0 - L), 1e-5);
  // load at the centroid mobilises the full translational resistance
  close(s.rigidPileLateralCapacity({ p: () => p, L, zLoad: 5 }).H, p * L, 1e-9);
  // load below the centroid: mirror image of a load above it
  close(s.rigidPileLateralCapacity({ p: () => p, L, zLoad: 8, n: 2000 }).H, s.rigidPileLateralCapacity({ p: () => p, L, zLoad: 2, n: 2000 }).H, 1e-6);
});

test('suction anchor horizontal capacity: Np su L D for constant Np; optimal padeye at 2L/3 for strength proportional to depth', () => {
  const r = s.caissonHorizontalCapacityClay({ D: 5, L: 15, su0: 0, k: 1500, NpConst: 10 });
  close(r.Hmax, 10 * (0.5 * 1500 * 15) * 15 * 5, 1e-9);
  close(r.zOptimal, 10, 1e-9);
  close(s.caissonHorizontalCapacityClay({ D: 5, L: 15, su0: 20e3, NpConst: 10 }).zOptimal, 7.5, 1e-9);
  // padeye away from the optimum loses capacity
  assert.ok(s.caissonHorizontalCapacityClay({ D: 5, L: 15, su0: 0, k: 1500, NpConst: 10, zPadeye: 5 }).H < r.Hmax);
  // Murff-Hamilton profile: the average Np lies between the mudline value 2 and the deep value 9
  const mh = s.caissonHorizontalCapacityClay({ D: 5, L: 15, su0: 5e3, k: 1500 });
  assert.ok(mh.NpAverage > 2 && mh.NpAverage < 9);
});

test('inclined-load interaction (Supachawarote et al. 2004): exponents a = L/D + 0.5, b = L/(3D) + 4.5; ultimate load on the envelope', () => {
  const r = s.inclinedInteraction({ Hmax: 10e6, Vmax: 6e6, L: 15, D: 5, angle: 35 * Math.PI / 180 });
  close(r.a, 3.5, 1e-14);
  close(r.b, 5.5, 1e-14);
  const H = r.Tult * Math.cos(35 * Math.PI / 180), V = r.Tult * Math.sin(35 * Math.PI / 180);
  close((H / 10e6) ** 3.5 + (V / 6e6) ** 5.5, 1, 1e-9);
  close(s.inclinedInteraction({ Hmax: 10e6, Vmax: 6e6, L: 15, D: 5, H, V }).f, 1, 1e-9);
  close(s.inclinedInteraction({ Hmax: 10e6, Vmax: 6e6, L: 15, D: 5, angle: 0 }).Tult, 10e6, 1e-9);
  close(s.inclinedInteraction({ Hmax: 10e6, Vmax: 6e6, L: 15, D: 5, angle: Math.PI / 2 }).Tult, 6e6, 1e-9);
});

test('bucket moment capacity: clay closed form with constant Np; sand (Broms) satisfies moment equilibrium', () => {
  const D = 8, L = 6, e = 30, su = 40e3, Np = 9;
  const p = Np * su * D;
  const z0 = -e + Math.sqrt(e * e + e * L + L * L / 2);
  const r = s.bucketMomentCapacity({ soil: 'clay', D, L, e, su0: su, NpConst: Np });
  close(r.H, p * (2 * z0 - L), 2e-4);
  close(r.M, r.H * e, 1e-12);
  // sand: p = c z with c = 3 Kp gamma D; H = c (z0^2 - L^2/2), moment about the load point
  const c = 3 * 3 * 10e3 * D;
  const sd = s.bucketMomentCapacity({ soil: 'sand', D, L, e, phi: Math.PI / 6, gamma: 10e3 });
  close(sd.H, c * (sd.z0 ** 2 - L * L / 2), 2e-4);
  // moment equilibrium about the mudline: H e = c (L^3 - 2 z0^3)/3
  close(sd.H * e, c * (L ** 3 - 2 * sd.z0 ** 3) / 3, 5e-4);
});

test('drained tensile capacity in sand (Houlsby et al. 2005): closed form equals integration of the reduced stress profile; rapid loading capped by cavitation', () => {
  const Do = 4, t = 0.02, h = 4, g = 9e3, Kt = 0.5, m = 1.5;
  const r = s.tensileCapacitySand({ Do, t, h, gamma: g, KtanDeltaO: Kt, KtanDeltaI: Kt, m, W: 100e3, waterDepth: 30 });
  const Zo = Do * (m * m - 1) / (4 * Kt);
  // integrate d(sigma)/dz = gamma - sigma/Z numerically, friction = Ktan(delta) pi Do integral(sigma)
  let sig = 0, integ = 0;
  const n = 200000, dz = h / n;
  for (let i = 0; i < n; i++) {
    const k1 = g - sig / Zo;
    const next = sig + dz * (g - (sig + 0.5 * dz * k1) / Zo);
    integ += 0.5 * (sig + next) * dz;
    sig = next;
  }
  close(r.frictionOut, Kt * Math.PI * Do * integ, 1e-5);
  assert.ok(r.frictionOut < 0.5 * g * h * h * Kt * Math.PI * Do);
  close(r.sCav, 101325 + 10055 * 30 - 2300, 1e-12);
  close(r.rapid, 100e3 + r.frictionOut + r.sCav * Math.PI * (Do - 2 * t) ** 2 / 4, 1e-12);
  close(r.drained, 100e3 + r.frictionOut + r.frictionIn, 1e-12);
});

test('suction in clay: the lid area of the required suction follows from the Houlsby & Byrne force balance (πDo²/4)', () => {
  // Houlsby & Byrne (2005a) Eq. 5: V′ + s·πDi²/4 = h αo su πDo + h αi su πDi + (γ′h − s + Nc su)·π D t
  // (suction acts on the lid over the inner area and also relieves the vertical stress at the skirt tip by s).
  // Moving the s-term of the tip to the left: s·[πDi²/4 + π D t] with D = (Do + Di)/2 the mean diameter;
  // π D t = π (Do² − Di²)/4, so the bracket is exactly πDo²/4 (their Eq. 6) — not the inner area.
  const Do = 5; const t = 0.03; const h = 6; const su = 20e3; const gamma = 6000; const Vp = 2e5;
  const r = s.suctionInstallClay({ Do, t, h, su: { su0: su, k: 0 }, alphaO: 0.5, alphaI: 0.5, Nc: 9, gamma, Vprime: Vp });
  const Di = Do - 2 * t; const D = Do - t;
  const lhs = Vp + r.s * Math.PI * Di * Di / 4;
  const rhs = h * 0.5 * su * Math.PI * Do + h * 0.5 * su * Math.PI * Di + (gamma * h - r.s + 9 * su) * Math.PI * D * t;
  assert.ok(Math.abs(lhs - rhs) <= 1e-12 * rhs, `Eq. 5 residual ${lhs - rhs} N`);
  assert.ok(Math.abs(Math.PI * Di * Di / 4 + Math.PI * D * t - Math.PI * Do * Do / 4) < 1e-12);
  assert.ok(Math.abs(r.s - (r.R - Vp) / (Math.PI * Do * Do / 4)) <= 1e-12 * r.s);
  // dividing by the inner area instead would overstate the suction by Do²/Di² − 1 = 2.4 %
  assert.ok(Math.abs((r.R - Vp) / (Math.PI * Di * Di / 4) / r.s - (Do / Di) ** 2) < 1e-12);
});
