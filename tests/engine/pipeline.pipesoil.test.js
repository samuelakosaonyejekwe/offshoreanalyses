import test from 'node:test';
import assert from 'node:assert/strict';
import * as P from '../../src/engine/pipeline/pipesoil.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${a} vs ${b}`);

test('contact geometry: segment area, wedging factor limits 1 and 4/pi', () => {
  const D = 0.5;
  const half = P.contactGeometry({ D, z: 0.25 });
  close(half.B, D, 1e-12);
  close(half.Abm, Math.PI * D * D / 8, 1e-12);
  close(half.zeta, 4 / Math.PI, 1e-12);
  close(P.contactGeometry({ D, z: 1e-9 }).zeta, 1, 1e-6);
  // z = D/4: cos(beta) = 0.5, beta = pi/3 ; B = D sin(60) ; A = D^2/4 (pi/3 - sqrt(3)/4)
  const q = P.contactGeometry({ D, z: 0.125 });
  close(q.beta, Math.PI / 3, 1e-12);
  close(q.B, D * Math.sqrt(3) / 2, 1e-12);
  close(q.Abm, (D * D / 4) * (Math.PI / 3 - Math.sqrt(3) / 4), 1e-12);
  close(q.zeta, 2 * (Math.sqrt(3) / 2) / (Math.PI / 3 + Math.sqrt(3) / 4), 1e-12);
  close(P.contactGeometry({ D, z: 0.35 }).Abm, Math.PI * D * D / 8 + D * 0.1, 1e-12);
  // numerical check of the wedging factor: normal stress ~ cos(theta) distribution on the arc
  let num = 0; let den = 0; const n = 20000;
  for (let i = 0; i < n; i++) { const th = -q.beta + (i + 0.5) * 2 * q.beta / n; num += Math.cos(th); den += Math.cos(th) ** 2; }
  close(q.zeta, num / den, 1e-6);
});

test('undrained vertical bearing: hand calculation', () => {
  // D = 0.5, z = 0.1 (z/D = 0.2), su = 5 kPa uniform, weightless soil:
  // min(6*0.2^0.25 = 4.0124 ; 3.4*sqrt(2) = 4.8083) = 4.0124 -> Qv = 4.0124*0.5*5000 = 10031 N/m
  const r = P.verticalBearingUndrained({ D: 0.5, z: 0.1, su0: 5000 });
  close(r.bearingFactor, 4.0124, 1e-4);
  close(r.Qv, 10031, 1e-4);
  // shallow branch governs for small z/D: z/D = 0.02 -> 3.4*sqrt(0.2) = 1.5205 < 6*0.02^0.25 = 2.2563
  close(P.verticalBearingUndrained({ D: 0.5, z: 0.01, su0: 5000 }).bearingFactor, 3.4 * Math.sqrt(0.2), 1e-12);
  // buoyancy term 1.5 gamma' Abm and strength gradient to the invert
  const b = P.verticalBearingUndrained({ D: 0.5, z: 0.1, su0: 2000, suGrad: 1500, gammaSub: 6000 });
  close(b.suInvert, 2150, 1e-12);
  close(b.Qv, 4.0124 * 0.5 * 2150 + 1.5 * 6000 * P.contactGeometry({ D: 0.5, z: 0.1 }).Abm, 1e-4);
});

test('drained vertical bearing: bearing capacity factors for phi = 30 deg', () => {
  // Nq = e^(pi tan30) tan^2(60) = 6.1337*3 = 18.401 ; Ngamma(Hansen) = 1.5*17.401*0.57735 = 15.070 ; Vesic = 2*19.401*0.57735 = 22.402
  const h = P.verticalBearingDrained({ D: 0.5, z: 0.05, phiDeg: 30, gammaSub: 9000 });
  close(h.Nq, 18.401, 1e-4);
  close(h.Ngamma, 15.070, 1e-4);
  close(P.verticalBearingDrained({ D: 0.5, z: 0.05, phiDeg: 30, gammaSub: 9000, nGamma: 'vesic' }).Ngamma, 22.402, 1e-4);
  const g = P.contactGeometry({ D: 0.5, z: 0.05 });
  const z0 = g.Abm / g.B;
  const dq = 1 + 1.2 * (z0 / g.B) * Math.tan(Math.PI / 6) * 0.25;
  close(h.Qv, 0.5 * 9000 * h.Ngamma * g.B ** 2 + z0 * 9000 * h.Nq * dq * g.B, 1e-12);
});

test('embedment: bearing resistance at the solved depth equals the factored load', () => {
  const a = { D: 0.4, su0: 3000, suGrad: 1500, gammaSub: 6000 };
  const e = P.embedment({ ...a, W: 1500, model: 'undrained', kLay: 2, kDyn: 1.5 });
  close(P.verticalBearingUndrained({ ...a, z: e.zStatic }).Qv, 3000, 1e-9);
  close(e.z, 1.5 * e.zStatic, 1e-12);
  const d = P.embedment({ D: 0.4, W: 1500, model: 'drained', phiDeg: 32, gammaSub: 9000 });
  close(P.verticalBearingDrained({ D: 0.4, z: d.zStatic, phiDeg: 32, gammaSub: 9000 }).Qv, 1500, 1e-9);
  assert.ok(P.embedment({ ...a, W: 3000, model: 'undrained' }).z > P.embedment({ ...a, W: 1500, model: 'undrained' }).z);
  assert.throws(() => P.embedment({ D: 0.4, W: 1e7, model: 'undrained', su0: 100 }), RangeError);
});

test('axial resistance: drained mu W zeta and undrained alpha su arc', () => {
  const dr = P.axialResistance({ model: 'drained', D: 0.5, z: 0.25, W: 2000, mu: 0.5, residualRatio: 0.8 });
  close(dr.breakout, 0.5 * 2000 * 4 / Math.PI, 1e-12);
  close(dr.residual, 0.8 * dr.breakout, 1e-12);
  close(dr.equivalentFriction, 0.5 * 4 / Math.PI, 1e-12);
  // z = D/4: arc = D*pi/3
  close(P.axialResistance({ model: 'undrained', D: 0.5, z: 0.125, alpha: 0.5, su: 4000 }).breakout, 0.5 * 4000 * 0.5 * Math.PI / 3, 1e-12);
});

test('lateral resistance: Coulomb friction plus Rankine passive wedge (hand calculation)', () => {
  // clay: mu W + 0.5 gamma' z^2 + 2 su z = 0.2*2000 + 0.5*6000*0.01 + 2*5000*0.1 = 400 + 30 + 1000 = 1430 N/m
  const c = P.lateralResistance({ soil: 'clay', D: 0.5, z: 0.1, W: 2000, mu: 0.2, su: 5000, gammaSub: 6000 });
  close(c.breakout, 1430, 1e-12);
  close(c.residual, 400, 1e-12);
  // sand phi = 30: Kp = 3 ; 0.6*2000 + 0.5*3*9000*0.01 = 1200 + 135 = 1335 N/m
  const s = P.lateralResistance({ soil: 'sand', D: 0.5, z: 0.1, W: 2000, mu: 0.6, phiDeg: 30, gammaSub: 9000, bermFactorResidual: 0.5 });
  close(s.breakout, 1335, 1e-10);
  close(s.residual, 1200 + 67.5, 1e-10);
  const v = P.lateralResistance({ model: 'verley', soil: 'sand', D: 0.4, z: 0.04, W: 144, mu: 0.6, gammaSub: 9000 });
  close(v.breakout, 0.6 * 144 + 144 * 35 * 0.1 ** 1.25, 1e-10);
});

test('mobilisation curve: tri-linear and odd', () => {
  const m = P.mobilisationCurve({ breakout: 1000, residual: 600, xBreakout: 0.01, xResidual: 0.5 });
  close(m.force(0.005), 500); close(m.force(0.01), 1000); close(m.force(0.255), 800); close(m.force(2), 600);
  close(m.force(-0.005), -500);
  close(m.stiffness, 1e5);
});

test('friction range: lognormal through LE (5 %) and HE (95 %)', () => {
  const f = P.frictionRange({ LE: 0.3, HE: 0.9 });
  close(f.quantile(0.05), 0.3, 1e-9);
  close(f.quantile(0.95), 0.9, 1e-9);
  close(f.BE, Math.sqrt(0.27), 1e-9);                        // median = geometric mean for symmetric quantiles
  close(f.sigma, Math.log(3) / (2 * 1.6448536), 1e-6);       // z(0.95) = 1.6448536
  close(f.mean, f.BE * Math.exp(0.5 * f.sigma ** 2), 1e-12);
  assert.throws(() => P.frictionRange({ LE: 0.9, HE: 0.3 }), RangeError);
});

test('soil stiffness (DNV-RP-F105 form): hand calculation for medium sand', () => {
  // KV = 14.5e6/(1 - 0.35)*(2/3*1.6 + 1/3)*sqrt(0.4) = 22.3077e6*1.4*0.632456 = 19.752e6 N/m2
  // KL = 12.5e6*1.35*1.4*0.632456 = 14.942e6 N/m2
  const k = P.soilStiffness({ D: 0.4, massRatio: 1.6, soilType: 'sand-medium' });
  close(k.KV, 19.752e6, 1e-4);
  close(k.KL, 14.942e6, 1e-4);
  assert.equal(k.KVS, 530e3);
  const c = P.soilStiffness({ D: 0.4, massRatio: 1.6, soilType: 'clay-soft' });
  assert.deepEqual(c.KVSrange, [160e3, 260e3]);
  close(P.soilStiffness({ D: 1, massRatio: 1, CV: 1e6, CL: 2e6, nu: 0 }).KL, 2e6, 1e-12);
  assert.throws(() => P.soilStiffness({ D: 0.4, massRatio: 1.6, soilType: 'peat' }), RangeError);
});

test('uplift resistance (Schaminee) and required cover: hand calculation and round trip', () => {
  // gamma' = 9000, H = 1.0, D = 0.4, f = 0.5: 9000*1*0.4*(1 + 0.5*1/0.4) = 3600*2.25 = 8100 N/m
  close(P.upliftResistanceSand({ D: 0.4, H: 1, gammaSub: 9000, f: 0.5 }), 8100, 1e-12);
  close(P.requiredCover({ soil: 'sand', R: 8100, D: 0.4, gammaSub: 9000, f: 0.5 }), 1, 1e-12);
  close(P.requiredCover({ soil: 'sand', R: 3600, D: 0.4, gammaSub: 9000, f: 0 }), 1, 1e-12);
  // clay: shallow 6000*0.5*0.4 + 2*4000*0.5 = 5200 ; deep 9*4000*0.4 = 14400
  const c = P.upliftResistanceClay({ D: 0.4, H: 0.5, su: 4000, gammaSub: 6000 });
  close(c.R, 5200, 1e-12);
  assert.equal(c.governing, 'shallow');
  close(P.requiredCover({ soil: 'clay', R: 5200, D: 0.4, gammaSub: 6000, su: 4000 }), 0.5, 1e-12);
  assert.equal(P.upliftResistanceClay({ D: 0.4, H: 5, su: 4000, gammaSub: 6000 }).governing, 'deep');
  assert.throws(() => P.requiredCover({ soil: 'clay', R: 20000, D: 0.4, su: 4000 }), RangeError);
  assert.equal(P.requiredCover({ soil: 'sand', R: -5, D: 0.4, gammaSub: 9000, f: 0.5 }), 0);
});

test('cover quantities: trapezoid berm geometry', () => {
  // D = 0.4, H = 0.6 -> h = 1.0 ; crest 1.0, slope 1:3 -> A = (1 + 3)*1 - pi*0.04 = 3.87434 m2 ; base 7 m
  const q = P.coverQuantities({ D: 0.4, H: 0.6, crestWidth: 1, sideSlope: 3, length: 100, rhoSolid: 2650, porosity: 0.4 });
  close(q.area, 4 - Math.PI * 0.04, 1e-12);
  close(q.baseWidth, 7, 1e-12);
  close(q.volume, 100 * q.area, 1e-12);
  close(q.mass, q.volume * 0.6 * 2650, 1e-12);
  close(q.gammaSub, 0.6 * 1625 * 9.80665, 1e-12);
});
