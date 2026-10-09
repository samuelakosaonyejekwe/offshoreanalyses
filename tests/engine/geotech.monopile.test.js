import test from 'node:test';
import assert from 'node:assert/strict';
import {
  rigidPileStiffness, flexiblePileStiffness, flexiblePileStiffnessContinuum, monopileStiffness,
  towerNaturalFrequency, mudlineRotationCheck, cyclicRotationAccumulation, embeddedLengthCriteria,
  apiMinimumWallThickness, monopileSizingScreen,
} from '../../src/engine/geotech/monopile.js';
import { pileHeadStiffness } from '../../src/engine/geotech/pile.js';
import { pyLinear } from '../../src/engine/geotech/py.js';
import { exampleProfile } from '../../src/engine/geotech/soilprofile.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-12), `${a} vs ${b}`);

test('rigid pile on Winkler springs: closed forms (kh D L, kh D L²/2, kh D L³/3; nh L²/2, nh L³/3, nh L⁴/4) and FE check', () => {
  const L = 30, k = 4e7, nh = 3e6;
  const c = rigidPileStiffness({ L, kRef: k });
  close(c.KLL, k * L); close(c.KLR, (-k * L * L) / 2); close(c.KRR, (k * L ** 3) / 3);
  const l = rigidPileStiffness({ L, kRef: nh, n: 1, zRef: 1 });
  close(l.KLL, (nh * L * L) / 2); close(l.KLR, (-nh * L ** 3) / 3); close(l.KRR, (nh * L ** 4) / 4);
  // numerical beam-on-Winkler solution with a practically rigid pile
  const fe = pileHeadStiffness({ L, nElem: 120, EI: 1e16, pyAt: (z) => pyLinear(nh * z), head: { H: 1e6 } }).tangent;
  close(fe.KHH, l.KLL, 2e-3); close(fe.KHM, l.KLR, 2e-3); close(fe.KMM, l.KRR, 2e-3);
});

test('flexible pile: Hetenyi (constant k) and the published linear-modulus coefficients against the FE solution', () => {
  const EI = 2e11, k = 5e7, nh = 4e6;
  const b = Math.pow(k / (4 * EI), 0.25);
  const c = flexiblePileStiffness({ EI, modulus: 'constant', k });
  close(c.KLL, k / b, 1e-12); close(c.KLR, -k / (2 * b * b), 1e-12); close(c.KRR, k / (2 * b ** 3), 1e-12);
  const l = flexiblePileStiffness({ EI, modulus: 'linear', nh });
  close(l.KLL, 1.074 * nh ** 0.6 * EI ** 0.4, 1e-12);
  // long pile (L = 8 T) with k = nh z solved by finite elements: within the rounding of the published
  // coefficients (exact inversion of the Reese–Matlock coefficients gives 1.076, −0.998, 1.497)
  const T = Math.pow(EI / nh, 0.2);
  const fe = pileHeadStiffness({ L: 8 * T, nElem: 400, EI, pyAt: (z) => pyLinear(nh * z), head: { H: 1e6 } }).tangent;
  close(fe.KHH, l.KLL, 0.01); close(fe.KHM, l.KLR, 0.015); close(fe.KMM, l.KRR, 0.02);
  const det = 2.435 * 1.75 - 1.623 ** 2;
  close(fe.KHH / (nh ** 0.6 * EI ** 0.4), 1.75 / det, 3e-3);
  close(fe.KMM / (nh ** 0.2 * EI ** 0.8), 2.435 / det, 3e-3);
});

test('continuum stiffness (EN 1998-5 Table C.1): hand evaluation', () => {
  const D = 2, Es = 50e6, EI = (30e9 * Math.PI * D ** 4) / 64; // Ep = 30 GPa → Ep/Es = 600
  const h = flexiblePileStiffnessContinuum({ D, EI, Es });
  close(h.Ep, 30e9, 1e-12);
  close(h.KLL, 1.08 * D * Es * 600 ** 0.21, 1e-12);
  close(h.KLR, -0.22 * D * D * Es * 600 ** 0.5, 1e-12);
  close(h.KRR, 0.16 * D ** 3 * Es * 600 ** 0.75, 1e-12);
  const p = flexiblePileStiffnessContinuum({ D, EI, Es, profile: 'parabolic' });
  close(p.KLL, 0.79 * D * Es * 600 ** 0.28, 1e-12);
  const l = flexiblePileStiffnessContinuum({ D, EI, Es, profile: 'linear' });
  close(l.KRR, 0.14 * D ** 3 * Es * 600 ** 0.8, 1e-12);
  assert.ok(h.KLL * h.KRR - h.KLR ** 2 > 0);
});

test('tower natural frequency: uniform cantilever f1 = (1.8751²/2π) sqrt(EI/(m L⁴)) and higher modes', () => {
  const L = 80, EI = 2e11, m = 3000;
  const r = towerNaturalFrequency({ segments: [{ L, EI, massPerLength: m }], nModes: 3 });
  const base = Math.sqrt(EI / (m * L ** 4)) / (2 * Math.PI);
  close(r.frequencies[0], 1.875104 ** 2 * base, 1e-5);
  close(r.frequencies[1], 4.694091 ** 2 * base, 1e-4);
  close(r.frequencies[2], 7.854757 ** 2 * base, 1e-3);
  assert.equal(r.modeShapes[0][0], 0);
  close(Math.abs(r.modeShapes[0][r.x.length - 1]), 1, 1e-12);
});

test('tower natural frequency: tip mass (massless beam sqrt(3EI/ML³)), foundation flexibility in series', () => {
  const L = 90, EI = 3e11, M = 4e5;
  const light = { L, EI, massPerLength: 1e-3 }; // practically massless tower
  const f0 = Math.sqrt((3 * EI) / (M * L ** 3)) / (2 * Math.PI);
  close(towerNaturalFrequency({ segments: [light], topMass: M, nModes: 1 }).frequencies[0], f0, 1e-6);
  // with distributed mass: Rayleigh estimate with the 33/140 effective mass (accurate to about 1 %)
  const mt = 2500;
  const fr = Math.sqrt((3 * EI) / ((M + (33 / 140) * mt * L) * L ** 3)) / (2 * Math.PI);
  close(towerNaturalFrequency({ segments: [{ L, EI, massPerLength: mt }], topMass: M, nModes: 1 }).frequencies[0], fr, 0.01);
  // massless tower on a foundation: tip flexibility = L³/3EI + fLL + 2 L fLR + L² fRR (flexibility superposition)
  const K = { KLL: 8e8, KLR: -6e9, KRR: 1.2e11 };
  const det = K.KLL * K.KRR - K.KLR ** 2, fLL = K.KRR / det, fLR = -K.KLR / det, fRR = K.KLL / det;
  const flex = L ** 3 / (3 * EI) + fLL + 2 * L * fLR + L * L * fRR;
  const ff = towerNaturalFrequency({ segments: [light], topMass: M, nModes: 1, ...K }).frequencies[0];
  close(ff, Math.sqrt(1 / (M * flex)) / (2 * Math.PI), 1e-6);
  assert.ok(ff < f0);
  // very stiff foundation recovers the fixed base
  close(towerNaturalFrequency({ segments: [light], topMass: M, nModes: 1, KLL: 1e16, KLR: 0, KRR: 1e18 }).frequencies[0], f0, 1e-6);
  // tapered tube: frequency lies between those of the uniform tubes with the end diameters
  const tube = (Db, Dt) => towerNaturalFrequency({ segments: [{ L, Dbot: Db, Dtop: Dt, t: 0.03 }], topMass: M, nModes: 1 }).frequencies[0];
  assert.ok(tube(6, 4) > tube(4, 4) && tube(6, 4) < tube(6, 6));
});

test('numerical monopile stiffness wraps the pile-head stiffness; embedded-length criteria on an elastic foundation', () => {
  const o = { L: 30, D: 6, t: 0.07, profile: exampleProfile('denseSand'), head: { H: 5e6, M: 150e6 } };
  const s = monopileStiffness(o), ref = pileHeadStiffness(o);
  assert.equal(s.tangent.KLL, ref.tangent.KHH); assert.equal(s.secant.KRR, ref.secant.KMM);
  assert.ok(s.secant.KLR < 0 && s.secant.KLL * s.secant.KRR > s.secant.KLR ** 2);
  // linear Winkler soil, constant k: under H alone the toe deflection of a free-ended beam changes sign as the
  // pile gets longer; Hetenyi classes beams with beta L > pi as long
  const EI = 1e11, k = 5e7, b = Math.pow(k / (4 * EI), 0.25);
  const c = embeddedLengthCriteria({ EI, pyAt: () => k, head: { H: 1e6 }, nElem: 120 }, { Lmin: 0.5 / b, Lmax: 8 / b, nLengths: 31 });
  assert.ok(c.toeY[0] < 0); // short rigid pile kicks back
  assert.ok(c.zeroToeKick !== null && c.verticalTangent !== null && c.critical !== null);
  assert.ok(c.zeroToeKick * b > 2 && c.zeroToeKick * b < 4, `beta L = ${c.zeroToeKick * b}`);
  assert.ok(c.critical * b > 2 && c.critical * b < 5, `beta Lcrit = ${c.critical * b}`);
  assert.ok(c.verticalTangent < c.zeroToeKick);
});

test('serviceability and sizing helpers', () => {
  const ok = mudlineRotationCheck({ thetaTotal: (0.4 * Math.PI) / 180, thetaAccumulated: (0.1 * Math.PI) / 180 });
  close(ok.ucTotal, 0.8, 1e-12); close(ok.ucAccumulated, 0.4, 1e-12); assert.equal(ok.ok, true);
  assert.equal(mudlineRotationCheck({ thetaTotal: 0.01 }).ok, false); // 0.573°
  // LeBlanc et al. (2010): Δθ/θs = Tb Tc N^0.31
  const a = cyclicRotationAccumulation({ thetaStatic: 0.002, N: 1e6, Tb: 0.05, Tc: 1 });
  close(a.ratio, 0.05 * Math.pow(1e6, 0.31), 1e-12); close(a.deltaTheta, 0.002 * a.ratio, 1e-12);
  assert.equal(cyclicRotationAccumulation({ thetaStatic: 0.002, N: 100, Tb: 0.05, Tc: 0 }).deltaTheta, 0);
  // API minimum wall: t = 6.35 + D/100 mm → D = 2000 mm gives 26.35 mm
  close(apiMinimumWallThickness(2), 0.02635, 1e-12);
  const s = monopileSizingScreen({ D: 8, t: 0.08 });
  assert.equal(s.dOverT, 100); assert.equal(s.bucklingCheckRequired, true);
  close(s.tMin, 0.08635, 1e-12); assert.equal(s.wallOk, false);
  close(s.massPerLength, (7850 * Math.PI * (64 - 7.84 ** 2)) / 4, 1e-12);
});
