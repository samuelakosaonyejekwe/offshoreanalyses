import test from 'node:test';
import assert from 'node:assert/strict';
import {
  correctedConeResistance, cptRatios, sbtZone, soilBehaviourType, relativeDensity, frictionAngleDeg,
  undrainedShearStrength, overconsolidationRatio, smallStrainShearModulus, unitWeightEstimate, cptToLayers,
  sptCorrection,
} from '../../src/engine/geotech/cpt.js';
import { createProfile, verticalEffectiveStress } from '../../src/engine/geotech/soilprofile.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-12), `${a} vs ${b}`);

test('corrected cone resistance and basic ratios by hand', () => {
  // qc = 1.0 MPa, u2 = 400 kPa, a = 0.75 → qt = 1.1 MPa
  const qt = correctedConeResistance({ qc: 1e6, u2: 400e3, a: 0.75 });
  close(qt, 1.1e6);
  const r = cptRatios({ qt, fs: 22e3, sigmaV0: 100e3, sigmaV0Eff: 50e3, u2: 400e3, u0: 50e3 });
  close(r.Rf, 2); close(r.qnet, 1e6); close(r.Qt, 20); close(r.Fr, 2.2); close(r.Bq, 0.35);
  assert.throws(() => cptRatios({ qt: 50e3, fs: 1e3, sigmaV0: 100e3, sigmaV0Eff: 50e3 }), RangeError);
});

test('Robertson (2009) Ic: iteration satisfies its defining equations; hand value for the first pass', () => {
  const o = { qt: 1e6, fs: 20e3, sigmaV0: 100e3, sigmaV0Eff: 50e3 };
  // first pass with n = 1: Qtn = 9 × 2 = 18, Fr = 2.2222 %, Ic = sqrt((3.47 − 1.2553)² + (0.3468 + 1.22)²) = 2.713
  const Ic1 = Math.hypot(3.47 - Math.log10(18), Math.log10(20 / 9) + 1.22);
  close(Ic1, 2.713, 5e-4);
  const r = soilBehaviourType(o);
  close(r.Fr, 100 * 20 / 900, 1e-12);
  const n = Math.min(0.381 * r.Ic + 0.05 * 0.5 - 0.15, 1);
  assert.ok(Math.abs(n - r.n) < 1e-5);
  close(r.Qtn, 9 * Math.pow(2, r.n), 1e-12);
  close(r.Ic, Math.hypot(3.47 - Math.log10(r.Qtn), Math.log10(r.Fr) + 1.22), 1e-12);
  // sigma'v < pa: an exponent below 1 reduces the stress normalisation, so Qtn falls and Ic rises slightly
  assert.ok(r.n < 1 && r.Qtn < 18 && r.Ic > Ic1);
  assert.equal(r.zone, 4);
  // dense clean sand: qt = 20 MPa, fs = 100 kPa at sigma'v = 100 kPa → zone 6
  const s = soilBehaviourType({ qt: 20e6, fs: 100e3, sigmaV0: 200e3, sigmaV0Eff: 100e3 });
  assert.equal(s.zone, 6); assert.equal(s.group, 'sand');
  // soft clay: qt = 0.5 MPa, fs = 15 kPa → zone 3 with n = 1
  const c = soilBehaviourType({ qt: 0.5e6, fs: 15e3, sigmaV0: 150e3, sigmaV0Eff: 70e3 });
  assert.equal(c.zone, 3); assert.equal(c.n, 1);
  assert.deepEqual([4, 3.2, 2.7, 2.3, 1.7, 1.0].map((v) => sbtZone(v).zone), [2, 3, 4, 5, 6, 7]);
});

test('relative density and friction angle correlations: published constants', () => {
  // Baldi et al. (1986): qc = 157 sigma'^0.55 exp(2.41 Dr) [kPa] — invert at Dr = 0.6, sigma' = 100 kPa
  const qc = 157 * 100 ** 0.55 * Math.exp(2.41 * 0.6) * 1000;
  close(relativeDensity({ qc, sigmaV0Eff: 100e3 }), 0.6, 1e-12);
  // Jamiolkowski et al. (2003): Dr = 0.268 ln(q/pa / sqrt(sigma'/pa)) − 0.675; q = 20 MPa, sigma' = 100 kPa
  close(relativeDensity({ qc: 20e6, sigmaV0Eff: 100e3, method: 'jamiolkowski' }), 0.268 * Math.log(200) - 0.675, 1e-12);
  // Kulhawy & Mayne (1990): 17.6 + 11 log10(100) = 39.6°
  close(frictionAngleDeg({ qc: 10e6, sigmaV0Eff: 100e3 }), 39.6, 1e-12);
  // Robertson & Campanella (1983): tan phi = (log10(100) + 0.29)/2.68 → 40.5°
  close(frictionAngleDeg({ qc: 10e6, sigmaV0Eff: 100e3, method: 'robertson-campanella' }), (Math.atan(2.29 / 2.68) * 180) / Math.PI, 1e-12);
  assert.throws(() => relativeDensity({ qc, sigmaV0Eff: 100e3, method: 'x' }), RangeError);
});

test('clay strength, OCR, Gmax and unit weight by hand', () => {
  close(undrainedShearStrength({ qt: 1.1e6, sigmaV0: 100e3, Nkt: 14 }), 1e6 / 14);
  close(overconsolidationRatio({ qt: 1.1e6, sigmaV0: 100e3, sigmaV0Eff: 50e3 }), 0.33 * 20);
  // Rix & Stokoe (1991): Gmax = 1634 qc^0.25 sigma'^0.375 [kPa]; qc = 10 MPa, sigma' = 100 kPa → 91.9 MPa
  close(smallStrainShearModulus({ qc: 10e6, sigmaV0Eff: 100e3 }).Gmax, 1634 * 10 * 100 ** 0.375 * 1000, 1e-12);
  assert.ok(Math.abs(smallStrainShearModulus({ qc: 10e6, sigmaV0Eff: 100e3 }).Gmax / 1e6 - 91.9) < 0.1);
  // Robertson (2009): alpha_vs = 10^(0.55·2 + 1.68) = 602.6; Vs = sqrt(602.6 × 9) = 73.6 m/s
  const g = smallStrainShearModulus({ method: 'robertson', qt: 1e6, sigmaV0: 100e3, Ic: 2, rho: 1800 });
  close(g.Vs, Math.sqrt(10 ** 2.78 * 9), 1e-12); close(g.Gmax, 1800 * g.Vs ** 2, 1e-12);
  // Robertson & Cabal (2010): Rf = 1 %, qt = 10 pa → gamma/gamma_w = 0.36 + 1.236
  close(unitWeightEstimate({ qt: 1e6, Rf: 1 }), 9810 * 1.596, 1e-12);
});

test('SPT corrections: N60 and (N1)60 with CN = sqrt(pa/sigma) <= 1.7', () => {
  const r = sptCorrection({ N: 20, energyRatioPct: 72, CR: 0.85, sigmaV0Eff: 50e3 });
  close(r.N60, 20 * 1.2 * 0.85); close(r.CN, Math.SQRT2); close(r.N160, r.N60 * Math.SQRT2);
  close(sptCorrection({ N: 10, sigmaV0Eff: 10e3 }).CN, 1.7);
});

test('layering of a synthetic CPTU sounding: sand over clay recovered and accepted by the soil profile', () => {
  // 0–10 m sand (qc growing 5 → 15 MPa, Rf 0.5 %), 10–30 m clay (qnet = 14 su, su = 40 + 2 (z − 10) kPa, Rf 3 %)
  const z = [], qc = [], fs = [];
  for (let i = 1; i <= 120; i++) {
    const zz = i * 0.25;
    z.push(zz);
    if (zz <= 10) { const q = 5e6 + 1e6 * zz; qc.push(q); fs.push(0.005 * q); } else { const q = 19e3 * zz + 14 * (40e3 + 2e3 * (zz - 10)); qc.push(q); fs.push(0.03 * q); }
  }
  const r = cptToLayers({ z, qc, fs, Nkt: 14, eps50: 0.01 });
  assert.equal(r.layers.length, 2);
  assert.equal(r.layers[0].type, 'sand'); assert.equal(r.layers[1].type, 'clay');
  close(r.layers[0].thickness, 10, 1e-12); close(r.layers[1].top, 10, 1e-12); close(r.layers[1].thickness, 20, 1e-12);
  assert.ok(r.layers[0].phiDeg > 35 && r.layers[0].phiDeg <= 45);
  // su recovered within the accuracy of the estimated total stress (unit weight correlation)
  const l = r.layers[1];
  assert.ok(Math.abs(l.suTop + l.suGrad * 10 - 60e3) < 6e3, `su mid ${l.suTop + l.suGrad * 10}`);
  assert.ok(l.suGrad > 0 && l.OCR >= 1);
  const p = createProfile(r.layers);
  assert.ok(verticalEffectiveStress(p, 30) > 150e3 && verticalEffectiveStress(p, 30) < 350e3);
  // a thin interbed is absorbed by the layer above
  const qc2 = qc.slice(), fs2 = fs.slice();
  for (let i = 60; i < 62; i++) { qc2[i] = 12e6; fs2[i] = 0.005 * 12e6; } // 0.5 m sand seam at 15 m
  assert.equal(cptToLayers({ z, qc: qc2, fs: fs2, Nkt: 14 }).layers.length, 2);
  assert.equal(cptToLayers({ z, qc: qc2, fs: fs2, Nkt: 14, minThickness: 0.25 }).layers.length, 4);
});
