import test from 'node:test';
import assert from 'node:assert/strict';
import {
  G, checkArg, darcyFrictionFactor, mixtureDensity, mixtureViscosity, homogeneous, harmathyVelocity,
  taylorBubbleVelocity, bendiksen, hasanKabir, driftClosure, zuberFindlay, driftFluxHoldup, voidFromSlip,
  slipRatioChisholm, slipRatioZivi, lockhartMartinelli, friedel,
} from '../../src/engine/flow/driftflux.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg ?? ''} got ${a}, expected ${b}`);
const airWater = { rhoL: 1000, rhoG: 1.2, muL: 1e-3, muG: 1.8e-5, sigma: 0.072 };

test('Darcy friction factor: 64/Re laminar and Colebrook–White residual in turbulent flow', () => {
  close(darcyFrictionFactor(1000), 0.064, 1e-14);
  for (const [Re, rr] of [[1e5, 1e-4], [1e7, 1e-3], [5e4, 0]]) {
    const f = darcyFrictionFactor(Re, rr);
    const residual = 1 / Math.sqrt(f) + 2 * Math.log10(rr / 3.7 + 2.51 / (Re * Math.sqrt(f)));
    assert.ok(Math.abs(residual) < 1e-8, `Colebrook residual ${residual}`);
  }
  // Smooth pipe at Re = 1e5: Moody chart value 0.0180.
  close(darcyFrictionFactor(1e5, 0), 0.0180, 0.01);
  assert.equal(darcyFrictionFactor(0), 0);
});

test('homogeneous model: no-slip holdup, density and quality by hand', () => {
  const h = homogeneous({ vsl: 1, vsg: 3, ...airWater });
  close(h.lambdaL, 0.25, 1e-14);
  close(h.rhoNs, 0.25 * 1000 + 0.75 * 1.2, 1e-14);
  close(h.muNs, 0.25 * 1e-3 + 0.75 * 1.8e-5, 1e-14);
  close(h.massFlux, 1003.6, 1e-14);
  close(h.quality, 3.6 / 1003.6, 1e-14);
});

test('single-phase limits: holdup → 1 without gas and → 0 without liquid', () => {
  assert.equal(homogeneous({ vsl: 2, vsg: 0, ...airWater }).holdup, 1);
  assert.equal(homogeneous({ vsl: 0, vsg: 2, ...airWater }).holdup, 0);
  assert.equal(zuberFindlay({ vsl: 2, vsg: 0, C0: 1.2, vd: 0.3 }).holdup, 1);
  assert.equal(mixtureDensity(1, 1000, 1.2), 1000);
  assert.equal(mixtureDensity(0, 1000, 1.2), 1.2);
});

test('no-slip limit: drift flux with C0 = 1, vd = 0 reproduces the homogeneous void fraction', () => {
  const z = zuberFindlay({ vsl: 0.7, vsg: 2.1, C0: 1, vd: 0 });
  close(z.alpha, 0.75, 1e-14);
  close(z.slipRatio, 1, 1e-13);
  const d = driftFluxHoldup({ vsl: 0.7, vsg: 2.1, D: 0.1, ...airWater, regime: 'homogeneous' });
  close(d.holdup, 0.25, 1e-13);
  close(voidFromSlip({ quality: 0.01, rhoL: 1000, rhoG: 10, S: 1 }), (0.01 / 10) / (0.01 / 10 + 0.99 / 1000), 1e-13);
});

test('Zuber–Findlay by hand: vsl = vsg = 1 m/s, C0 = 1.2, vd = 0.245 → α = 1/2.645', () => {
  const z = zuberFindlay({ vsl: 1, vsg: 1, C0: 1.2, vd: 0.245 });
  close(z.alpha, 1 / 2.645, 1e-14);
  close(z.vG, 2.645, 1e-14);
  close(z.vL, 1 / (1 - 1 / 2.645), 1e-14);
  // volumetric identity vm = α vG + (1 − α) vL
  close(z.alpha * z.vG + z.holdup * z.vL, 2, 1e-14);
});

test('Harmathy (1960) bubble rise velocity by hand for air–water', () => {
  // g σ Δρ / ρL² = 9.80665 × 0.072 × 998.8 / 1e6 = 7.0523e-4; fourth root = 0.16296; × 1.53 = 0.2493 m/s
  close(harmathyVelocity(airWater), 0.2493, 1e-3);
  close(harmathyVelocity(airWater), 1.53 * (G * 0.072 * 998.8 / 1e6) ** 0.25, 1e-14);
});

test('Taylor-bubble drift 0.35 √(gD) (Nicklin et al. 1962) by hand for D = 50 mm', () => {
  // √(9.80665 × 0.05) = 0.70024; × 0.35 = 0.24508 m/s
  close(taylorBubbleVelocity({ D: 0.05 }), 0.24508, 1e-4);
  close(taylorBubbleVelocity({ D: 0.05, rhoL: 1000, rhoG: 100 }), 0.35 * Math.sqrt(G * 0.05 * 0.9), 1e-14);
});

test('Bendiksen (1984): horizontal and vertical coefficients and the Froude 3.5 switch', () => {
  const D = 0.05;
  const sgd = Math.sqrt(G * D);
  const h = bendiksen({ vm: 1, D, theta: 0 });
  close(h.C0, 1.05, 1e-14); close(h.vd, 0.54 * sgd, 1e-14);
  const v = bendiksen({ vm: 1, D, theta: Math.PI / 2 });
  close(v.C0, 1.2, 1e-14); close(v.vd, 0.35 * sgd, 1e-13);
  const i = bendiksen({ vm: 1, D, theta: Math.PI / 6 });
  close(i.C0, 1.05 + 0.15 * 0.25, 1e-13); close(i.vd, (0.35 * 0.5 + 0.54 * Math.cos(Math.PI / 6)) * sgd, 1e-13);
  const below = bendiksen({ vm: 3.49 * sgd, D, theta: 0 });
  const above = bendiksen({ vm: 3.51 * sgd, D, theta: 0 });
  assert.equal(below.C0, 1.05);
  assert.equal(above.C0, 1.2);
  assert.equal(above.vd, 0);
  assert.equal(bendiksen({ vm: 1, D, theta: -0.1 }).extrapolated, true);
});

test('regime closures: bubbly, slug, annular and Hasan–Kabir parameter sets', () => {
  const D = 0.1;
  const b = driftClosure({ regime: 'bubbly', D, ...airWater });
  close(b.C0, 1.2, 1e-14); close(b.vd, harmathyVelocity(airWater), 1e-14);
  const s = driftClosure({ regime: 'slug', D, ...airWater });
  close(s.vd, 0.35 * Math.sqrt(G * D * 998.8 / 1000), 1e-13);
  assert.deepEqual(driftClosure({ regime: 'annular', D, ...airWater }), { C0: 1, vd: 0 });
  assert.equal(driftClosure({ regime: 'slug', D, ...airWater, C0: 1.1, vd: 0.2 }).C0, 1.1);
  assert.equal(hasanKabir({ regime: 'churn', D, ...airWater }).C0, 1.15);
  // vertical: inclination factor √(sin 90°)(1 + cos 90°)^1.2 = 1
  close(hasanKabir({ regime: 'slug', D, ...airWater }).vd, taylorBubbleVelocity({ D, rhoL: 1000, rhoG: 1.2 }), 1e-12);
  assert.throws(() => driftClosure({ regime: 'foam', D, ...airWater }), RangeError);
});

test('slip ratio correlations: Chisholm limits, Zivi cube root', () => {
  close(slipRatioChisholm({ quality: 0, rhoL: 1000, rhoG: 10 }), 1, 1e-14);
  close(slipRatioChisholm({ quality: 1, rhoL: 1000, rhoG: 10 }), 10, 1e-14);
  close(slipRatioChisholm({ quality: 0.1, rhoL: 1000, rhoG: 10 }), Math.sqrt(1 + 0.1 * 99), 1e-14);
  close(slipRatioZivi({ rhoL: 1000, rhoG: 8 }), 5, 1e-14);
});

test('mixture viscosity rules reduce to the pure-phase values', () => {
  const mu = { muL: 1e-3, muG: 2e-5 };
  for (const rule of ['linear', 'arrhenius', 'beattie-whalley']) {
    close(mixtureViscosity({ rule, ...mu, holdup: 1 }), 1e-3, 1e-13, rule);
    close(mixtureViscosity({ rule, ...mu, holdup: 0 }), 2e-5, 1e-13, rule);
  }
  for (const rule of ['mcadams', 'cicchitti']) {
    close(mixtureViscosity({ rule, ...mu, quality: 0 }), 1e-3, 1e-13, rule);
    close(mixtureViscosity({ rule, ...mu, quality: 1 }), 2e-5, 1e-13, rule);
  }
  close(mixtureViscosity({ rule: 'mcadams', ...mu, quality: 0.5 }), 1 / (0.5 / 2e-5 + 0.5 / 1e-3), 1e-14);
  assert.throws(() => mixtureViscosity({ rule: 'nope', ...mu, holdup: 0.5 }), RangeError);
});

test('Lockhart–Martinelli / Chisholm: C = 20 turbulent–turbulent, φL²·(dP/dx)L = φG²·(dP/dx)G', () => {
  const r = lockhartMartinelli({ vsl: 1, vsg: 10, D: 0.05, ...airWater });
  assert.equal(r.C, 20);
  close(r.X * r.X, r.gradL / r.gradG, 1e-13);
  close(r.phiL2, 1 + 20 / r.X + 1 / r.X ** 2, 1e-14);
  close(r.phiL2 * r.gradL, r.phiG2 * r.gradG, 1e-12);
  assert.ok(r.gradFriction > r.gradL && r.gradFriction > r.gradG);
  close(r.alphaWallis, (1 + r.X ** 0.8) ** -0.378, 1e-14);
  // viscous–viscous
  assert.equal(lockhartMartinelli({ vsl: 0.01, vsg: 0.1, D: 0.01, ...airWater, muL: 0.1 }).C, 5);
  // Butterworth form of the L–M holdup by hand
  const x = r.X; void x;
  const mL = 1000 * 1; const mG = 1.2 * 10;
  close(r.alphaButterworth, 1 / (1 + 0.28 * (mL / mG) ** 0.64 * (1.2 / 1000) ** 0.36 * (1e-3 / 1.8e-5) ** 0.07), 1e-13);
});

test('Friedel (1979): single-phase limits x = 0 (liquid only) and x = 1 (gas only), and a hand evaluation', () => {
  const base = { massFlux: 500, D: 0.05, ...airWater, rhoG: 10 };
  const liq = friedel({ ...base, quality: 0 });
  close(liq.phiLo2, 1, 1e-13);
  const gas = friedel({ ...base, quality: 1 });
  const fGo = darcyFrictionFactor(500 * 0.05 / 1.8e-5);
  close(gas.gradFriction, fGo * 500 ** 2 / (2 * 0.05 * 10), 1e-12);
  const r = friedel({ ...base, quality: 0.2 });
  const fLo = darcyFrictionFactor(500 * 0.05 / 1e-3);
  const rhoH = 1 / (0.2 / 10 + 0.8 / 1000);
  const E = 0.64 + 0.04 * (1000 * fGo) / (10 * fLo);
  const F = 0.2 ** 0.78 * 0.8 ** 0.224;
  const H = 100 ** 0.91 * 0.018 ** 0.19 * (1 - 0.018) ** 0.7;
  const Fr = 500 ** 2 / (G * 0.05 * rhoH ** 2);
  const We = 500 ** 2 * 0.05 / (0.072 * rhoH);
  close(r.phiLo2, E + 3.24 * F * H / (Fr ** 0.045 * We ** 0.035), 1e-12);
  assert.ok(r.phiLo2 > 1);
});

test('invalid arguments throw RangeError naming the argument', () => {
  assert.throws(() => checkArg('D', -1, { min: 0 }), /D/);
  assert.throws(() => homogeneous({ vsl: -1, vsg: 1, ...airWater }), RangeError);
  assert.throws(() => homogeneous({ vsl: 0, vsg: 0, ...airWater }), RangeError);
  assert.throws(() => zuberFindlay({ vsl: 1, vsg: 1, C0: 0, vd: 0 }), /C0/);
  assert.throws(() => bendiksen({ vm: 1, D: 0 }), /D/);
  assert.throws(() => friedel({ massFlux: 0, quality: 0.5, D: 0.05, ...airWater }), /massFlux/);
});
