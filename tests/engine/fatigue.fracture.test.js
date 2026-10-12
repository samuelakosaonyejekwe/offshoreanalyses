import test from 'node:test';
import assert from 'node:assert/strict';
import {
  centreCrackK, edgeCrackK, doubleEdgeCrackK, ellipticE, embeddedEllipticalCrackK, newmanRaju, surfaceCrackK,
  parisConstantToSI, kFromNmm32, parisLawBs7910, crackGrowthRate, parisLifeClosedForm, parisCrackSizeClosedForm,
  parisIntegrate, surfaceCrackGrowth, fadOption1, referenceStressSurfaceCrack, referenceStressThroughCrack,
  plasticityInteractionRho, assessSurfaceCrack, assessThroughCrack, criticalCrackSize, limitingStress, ctodToK,
  kToCtod, charpyToK, leakBeforeBreak,
} from '../../src/engine/fatigue/fracture.js';

const MPa = 1e6;
const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('through-crack solutions (Tada): limits and hand values', () => {
  close(centreCrackK({ a: 1e-6, W: 1 }).F, 1, 1e-9);
  // 2a/W = 0.5: F = (1 − 0.025·0.25 + 0.06·0.0625)·√sec(π/4) = 0.9975·1.189207 = 1.18623
  close(centreCrackK({ a: 0.05, W: 0.2 }).F, 1.18623, 1e-5);
  close(centreCrackK({ a: 0.05, W: 0.2, stress: 100 * MPa }).K, 100 * MPa * Math.sqrt(Math.PI * 0.05) * 1.18623, 1e-5);
  // Edge crack: free-surface factor 1.122 for a short crack, in tension and in bending.
  close(edgeCrackK({ a: 1e-7, W: 1, tension: 1 }).Ft, 1.122, 1e-5);
  close(edgeCrackK({ a: 1e-7, W: 1, bending: 1 }).Fb, 1.122, 1e-5);
  // a/W = 0.5 written out: √(tan(π/4)/(π/4))/cos(π/4) = 1.59577;
  // F_t = 1.59577·(0.752 + 1.01 + 0.37·0.29289³) = 2.8266; F_b = 1.59577·(0.923 + 0.199·0.29289⁴) = 1.4752
  close(edgeCrackK({ a: 0.5, W: 1 }).Ft, 2.8266, 1e-4);
  close(edgeCrackK({ a: 0.5, W: 1 }).Fb, 1.4752, 1e-4);
  // Independent check against the Brown & Srawley (ASTM STP 410) polynomial fits, valid for a/W <= 0.6.
  for (const al of [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]) {
    const e = edgeCrackK({ a: al, W: 1 });
    close(e.Ft, 1.122 - 0.231 * al + 10.55 * al ** 2 - 21.72 * al ** 3 + 30.39 * al ** 4, 0.01, `tension a/W ${al}`);
    close(e.Fb, 1.122 - 1.40 * al + 7.33 * al ** 2 - 13.08 * al ** 3 + 14.0 * al ** 4, 0.025, `bending a/W ${al}`);
  }
  // Deep-crack limit in bending: K → 3.975·M/(W − a)^1.5 per unit thickness (σ = 6M/W²).
  const deep = edgeCrackK({ a: 0.98, W: 1, bending: 6 });
  close(deep.K, 3.975 / 0.02 ** 1.5, 0.01);
  close(doubleEdgeCrackK({ a: 1e-7, W: 1 }).F, 1.122, 1e-5);
});

test('embedded crack (Irwin): penny crack K = (2/π)·σ·√(πa); complete elliptic integral values', () => {
  close(ellipticE(0), Math.PI / 2, 1e-14);
  assert.equal(ellipticE(1), 1);
  close(ellipticE(0.5), 1.4674622093394272, 1e-12); // Abramowitz & Stegun Table 17.1 (m = 0.25)
  close(embeddedEllipticalCrackK({ a: 0.01, c: 0.01, stress: 1 }), (2 / Math.PI) * Math.sqrt(Math.PI * 0.01), 1e-12);
  // Long ellipse → through-crack limit at the minor-axis end: K → σ√(πa).
  close(embeddedEllipticalCrackK({ a: 0.01, c: 100, stress: 1 }), Math.sqrt(Math.PI * 0.01), 1e-3);
});

test('Newman–Raju (1981): coefficients written out by hand for a/c = 0.2, a/t = 0.5, deepest point', () => {
  // M1 = 1.13 − 0.018 = 1.112; M2 = −0.54 + 0.89/0.4 = 1.685; M3 = 0.5 − 1/0.85 + 14·0.8^24 = −0.610358
  // F = 1.112 + 1.685·0.25 − 0.610358·0.0625 = 1.495103; Q = 1 + 1.464·0.2^1.65 = 1.102859
  const r = newmanRaju({ a: 0.01, c: 0.05, t: 0.02 });
  close(r.Q, 1.102859, 1e-5);
  close(r.F, 1.495103, 1e-5);
  close(r.Mm, 1.495103 / Math.sqrt(1.102859), 1e-5);
  // Bending: G21 = −1.244, G22 = 0.55 − 1.05·0.2^0.75 + 0.47·0.2^1.5 = 0.278015, H2 = 1 − 0.622 + 0.069504 = 0.447504
  close(r.H, 0.447504, 1e-5);
  // Surface point: H = H1 = 1 − 0.34·0.5 − 0.11·0.2·0.5 = 0.819; g = 1 + 0.1 + 0.35·0.25 = 1.1875; fφ = √0.2
  const s = newmanRaju({ a: 0.01, c: 0.05, t: 0.02, phi: 0 });
  close(s.H, 0.819, 1e-12);
  close(s.F, 1.495103 * 1.1875 * Math.sqrt(0.2), 1e-5);
});

test('Newman–Raju: shallow semicircular crack, finite width and continuity of the a/c > 1 extension', () => {
  // a/c = 1, a/t → 0, deepest point: F = M1 = 1.04, Q = 2.464 → K = 0.6625·σ·√(πa) (penny crack 0.6366 × ~1.04).
  const r = newmanRaju({ a: 1e-6, c: 1e-6, t: 1 });
  close(r.Mm, 1.04 / Math.sqrt(2.464), 1e-6);
  // Finite width: fw = √sec(π·c/(2b)·√(a/t)); c/b = 0.4, a/t = 0.64 → √sec(0.16π) = 1.06834
  const w = newmanRaju({ a: 0.016, c: 0.04, t: 0.025, b: 0.1 });
  close(w.fw, 1 / Math.sqrt(Math.cos(0.16 * Math.PI)), 1e-12);
  // Both branches agree at a/c = 1 (the published coefficient sets differ by < 0.5 % there).
  for (const phi of [0, Math.PI / 4, Math.PI / 2]) {
    const lo = newmanRaju({ a: 0.01, c: 0.01, t: 0.03, phi });
    const hi = newmanRaju({ a: 0.010001, c: 0.01, t: 0.03, phi });
    close(hi.Mm, lo.Mm, 5e-3, `Mm at phi ${phi}`);
    close(hi.Mb, lo.Mb, 5e-3, `Mb at phi ${phi}`);
  }
  const k = surfaceCrackK({ a: 0.01, c: 0.05, t: 0.02, membrane: 100 * MPa, bending: 50 * MPa, Mk: { deepMembrane: 1.2 } });
  const d = newmanRaju({ a: 0.01, c: 0.05, t: 0.02 });
  close(k.Kdeep, (1.2 * d.Mm * 100 * MPa + d.Mb * 50 * MPa) * Math.sqrt(Math.PI * 0.01), 1e-12);
});

test('BS 7910 crack-growth laws: unit conversion and stage A/B transition points of the published two-stage laws', () => {
  // ΔK = 10 MPa√m = 316.23 N/mm^1.5 → 5.21e-13·316.23³ = 1.6476e-5 mm/cycle = 1.6476e-8 m/cycle
  const air = parisLawBs7910();
  close(crackGrowthRate(10e6, air), 1.6476e-8, 1e-4);
  close(air.stages[0].A, parisConstantToSI(5.21e-13, 3), 1e-14);
  close(kFromNmm32(63), 63 / 31.6227766e-6, 1e-8);
  assert.equal(crackGrowthRate(kFromNmm32(60), air), 0); // below the 63 N/mm^1.5 threshold
  close(crackGrowthRate(10e6, parisLawBs7910({ environment: 'marine' })) / crackGrowthRate(10e6, air), 2.3e-12 / 5.21e-13, 1e-12);
  // Published transition ΔK [N/mm^1.5]: the intersection of the two published stages must reproduce it.
  const cases = [
    ['air', 'high', 'mean', 196], ['air', 'high', 'mean+2sd', 144], ['air', 'low', 'mean', 363], ['air', 'low', 'mean+2sd', 315],
    ['marine-free', 'low', 'mean', 1336], ['marine-free', 'low', 'mean+2sd', 993], ['marine-free', 'high', 'mean', 1098],
    ['marine-free', 'high', 'mean+2sd', 748], ['marine-cp850', 'low', 'mean', 462], ['marine-cp850', 'low', 'mean+2sd', 434],
    ['marine-cp850', 'high', 'mean', 323], ['marine-cp850', 'high', 'mean+2sd', 290],
  ];
  for (const [environment, rRatio, curve, kt] of cases) {
    const law = parisLawBs7910({ environment, rRatio, curve });
    close(law.transition, kFromNmm32(kt), 4e-3, `${environment} ${rRatio} ${curve}`);
    // Below the transition stage A governs, above it stage B.
    const lowK = law.transition * 0.8;
    const highK = law.transition * 1.5;
    const noTh = { stages: law.stages };
    close(crackGrowthRate(lowK, noTh), law.stages[0].A * lowK ** law.stages[0].m, 1e-12);
    close(crackGrowthRate(highK, noTh), law.stages[1].A * highK ** law.stages[1].m, 1e-12);
  }
  // Walker correction: R = 0.5, γ = 0.5 → ΔK_eff = ΔK/√0.5.
  const law = { stages: [{ A: 1e-30, m: 3 }] };
  close(crackGrowthRate(1e7, law, { R: 0.5, walkerGamma: 0.5 }), 1e-30 * (1e7 / Math.sqrt(0.5)) ** 3, 1e-12);
});

test('Paris integration against the closed form for a constant geometry factor', () => {
  const A = parisConstantToSI(5.21e-13, 3);
  const args = { a0: 0.5e-3, af: 20e-3, stressRange: 80 * MPa, Y: 1.12, A, m: 3 };
  // Closed form written out: N = (af^-0.5 − a0^-0.5)/(−0.5·A·(Y·Δσ·√π)³)
  const k = A * (1.12 * 80e6 * Math.sqrt(Math.PI)) ** 3;
  const expected = (20e-3 ** -0.5 - 0.5e-3 ** -0.5) / (-0.5 * k);
  close(parisLifeClosedForm(args), expected, 1e-12);
  const law = { stages: [{ A, m: 3 }], threshold: 0 };
  const num = parisIntegrate({ a0: 0.5e-3, af: 20e-3, deltaK: (a) => 1.12 * 80e6 * Math.sqrt(Math.PI * a), law });
  close(num.cycles, expected, 1e-4);
  assert.equal(num.reached, true);
  close(parisCrackSizeClosedForm({ ...args, cycles: expected }), 20e-3, 1e-9);
  // m = 2 (logarithmic form) and m = 4.
  for (const m of [2, 4]) {
    const Am = 1e-11 * (1e-6) ** m;
    const cf = parisLifeClosedForm({ a0: 1e-3, af: 1e-2, stressRange: 50e6, Y: 1, A: Am, m });
    const ni = parisIntegrate({ a0: 1e-3, af: 1e-2, deltaK: (a) => 50e6 * Math.sqrt(Math.PI * a), law: { stages: [{ A: Am, m }] } });
    close(ni.cycles, cf, 1e-4, `m = ${m}`);
  }
  // Below threshold the crack is arrested.
  const arrested = parisIntegrate({ a0: 1e-4, af: 1e-2, deltaK: () => 1e6, law: parisLawBs7910() });
  assert.equal(arrested.arrested, true);
  assert.equal(arrested.cycles, Infinity);
  // Cycle cap.
  const capped = parisIntegrate({ a0: 0.5e-3, af: 20e-3, deltaK: (a) => 1.12 * 80e6 * Math.sqrt(Math.PI * a), law, maxCycles: expected / 2 });
  close(capped.a, parisCrackSizeClosedForm({ ...args, cycles: expected / 2 }), 1e-4);
});

test('surface-crack growth in a and c: single-stage scaling N ∝ Δσ^(−m) and growth to wall break-through', () => {
  const law = { stages: [{ A: parisConstantToSI(5.21e-13, 3), m: 3 }], threshold: 0 };
  const base = { a0: 1e-3, c0: 2e-3, t: 0.025, W: 1, law };
  const r1 = surfaceCrackGrowth({ ...base, membraneRange: 100 * MPa });
  const r2 = surfaceCrackGrowth({ ...base, membraneRange: 200 * MPa });
  assert.equal(r1.status, 'through-thickness');
  close(r1.a, 0.025, 1e-9);
  close(r1.cycles / r2.cycles, 8, 1e-9);
  close(r1.c, r2.c, 1e-9);
  // Under pure tension the aspect ratio tends towards a/c ≈ 0.8–1 before break-through.
  assert.ok(r1.a / r1.c > 0.6 && r1.a / r1.c < 1.1, `a/c = ${r1.a / r1.c}`);
  // Bending keeps the crack shallower and longer than tension does.
  const rb = surfaceCrackGrowth({ ...base, bendingRange: 100 * MPa, throughFraction: 0.8 });
  const rt = surfaceCrackGrowth({ ...base, membraneRange: 100 * MPa, throughFraction: 0.8 });
  assert.ok(rb.a / rb.c < rt.a / rt.c);
  assert.ok(rb.cycles > rt.cycles);
  // Step-size convergence.
  const fine = surfaceCrackGrowth({ ...base, membraneRange: 100 * MPa, stepFraction: 0.001 });
  close(r1.cycles, fine.cycles, 2e-4);
  // Fracture stop and cycle cap.
  const fr = surfaceCrackGrowth({ ...base, membraneRange: 100 * MPa, toughness: 15e6 });
  assert.equal(fr.status, 'fracture');
  assert.ok(fr.a < 0.025);
  const cap = surfaceCrackGrowth({ ...base, membraneRange: 100 * MPa, maxCycles: r1.cycles / 2 });
  assert.equal(cap.status, 'max-cycles');
  close(cap.cycles, r1.cycles / 2, 1e-12);
});

test('BS 7910 Option 1 FAD: f(0) = 1, hand values at Lr = 1 and beyond, cut-off at Lr,max', () => {
  const mat = { yieldStrength: 355 * MPa, tensileStrength: 510 * MPa, E: 207e9 };
  assert.equal(fadOption1(0, mat).Kr, 1);
  // μ = min(0.001·207000/355, 0.6) = 0.583099; f(1) = (1.5)^-0.5·(0.3 + 0.7·e^-0.583099) = 0.563965
  const one = fadOption1(1, mat);
  close(one.mu, 0.583099, 1e-5);
  close(one.Kr, (0.3 + 0.7 * Math.exp(-0.5830986)) / Math.sqrt(1.5), 1e-6);
  close(one.Kr, 0.563965, 1e-4);
  close(one.LrMax, (355 + 510) / 710, 1e-12);
  // Lr = 1.1: N = 0.3·(1 − 355/510) = 0.091176, exponent (N − 1)/(2N) = −4.98387, f = f(1)·1.1^−4.98387 = 0.35072
  close(fadOption1(1.1, mat).Kr, 0.563965 * 1.1 ** -4.98387, 1e-4);
  assert.equal(fadOption1(1.25, mat).Kr, 0);
  // μ capped at 0.6 for low-strength steel.
  assert.equal(fadOption1(1, { yieldStrength: 235 * MPa, tensileStrength: 360 * MPa }).mu, 0.6);
});

test('reference stress and plasticity interaction: hand values', () => {
  // a/t = 0.5, c = t, wide plate: α″ = 0.5/(1 + 1) = 0.25; Pm = 100: σref = √(9·100²·0.5625)/(3·0.5625) = 133.33
  close(referenceStressSurfaceCrack({ a: 0.01, c: 0.02, t: 0.02, Pm: 100 * MPa }), (100 * MPa) / 0.75, 1e-12);
  // With bending Pb = 60: [60 + √(3600 + 9·10000·0.5625)]/(1.6875) = (60 + 232.86)/1.6875 = 173.55
  close(referenceStressSurfaceCrack({ a: 0.01, c: 0.02, t: 0.02, Pm: 100 * MPa, Pb: 60 * MPa }), ((60 + Math.sqrt(3600 + 50625)) / 1.6875) * MPa, 1e-12);
  close(referenceStressThroughCrack({ a: 0.05, W: 0.5, Pm: 100 * MPa }), 125 * MPa, 1e-12);
  // χ = 1 at Lr = 0.5: ρ1 = 0.1 − 0.007 + 0.00003 = 0.09303; tapering to zero at Lr = 1.05.
  close(plasticityInteractionRho({ Kp: 50, Ks: 100, Lr: 0.5 }), 0.09303, 1e-12);
  close(plasticityInteractionRho({ Kp: 50, Ks: 50, Lr: 1.0 }), 4 * 0.09303 * 0.05, 1e-10);
  assert.equal(plasticityInteractionRho({ Kp: 50, Ks: 50, Lr: 1.05 }), 0);
  assert.equal(plasticityInteractionRho({ Kp: 50, Ks: 0, Lr: 0.5 }), 0);
});

test('engineering critical assessment: point, reserve factor, critical size and limiting stress are mutually consistent', () => {
  const base = { a: 4e-3, c: 12e-3, t: 0.03, W: 1, Pm: 150 * MPa, Pb: 30 * MPa, Qm: 200 * MPa, yieldStrength: 355 * MPa, tensileStrength: 510 * MPa, toughness: 120e6 };
  const r = assessSurfaceCrack(base);
  // Recompute the assessment point independently from the building blocks.
  const k = surfaceCrackK({ a: base.a, c: base.c, t: base.t, W: 1, membrane: 1 });
  const kb = surfaceCrackK({ a: base.a, c: base.c, t: base.t, W: 1, bending: 1 });
  const Lr = referenceStressSurfaceCrack(base) / base.yieldStrength;
  close(r.Lr, Lr, 1e-12);
  const pts = ['Kdeep', 'Ksurface'].map((p) => {
    const Kp = k[p] * 150e6 + kb[p] * 30e6;
    const Ks = k[p] * 200e6;
    return (Kp + Ks) / 120e6 + plasticityInteractionRho({ Kp, Ks, Lr });
  });
  close(r.Kr, Math.max(...pts), 1e-12);
  assert.equal(r.acceptable, r.Kr <= r.KrAllowable);
  assert.equal(r.acceptable, true);
  close(r.margin, r.KrAllowable / r.Kr, 1e-12);
  // Scaling the primary stresses by the reserve factor puts the point on the assessment line.
  const lim = limitingStress(base);
  const on = assessSurfaceCrack({ ...base, Pm: lim.Pm, Pb: lim.Pb });
  close(on.Kr, on.KrAllowable, 1e-6);
  assert.ok(r.reserveFactor > 1);
  // Critical depth at fixed aspect ratio: acceptable just below, unacceptable just above.
  const crit = criticalCrackSize({ ...base, aspectRatio: 1 / 3 });
  assert.equal(crit.limitedBy, 'fad');
  assert.ok(crit.a > base.a);
  assert.equal(assessSurfaceCrack({ ...base, a: crit.a * 0.999, c: crit.c * 0.999 }).acceptable, true);
  assert.equal(assessSurfaceCrack({ ...base, a: crit.a * 1.001, c: crit.c * 1.001 }).acceptable, false);
  // Low toughness shrinks the critical size; a primary stress beyond plastic collapse leaves nothing acceptable.
  assert.ok(criticalCrackSize({ ...base, toughness: 20e6, aspectRatio: 1 / 3 }).a < 0.2 * crit.a);
  assert.equal(criticalCrackSize({ ...base, Pm: 450 * MPa, aspectRatio: 1 / 3 }).limitedBy, 'none-acceptable');
});

test('through-crack assessment and leak-before-break', () => {
  const mat = { yieldStrength: 355 * MPa, tensileStrength: 510 * MPa, toughness: 150e6 };
  const r = assessThroughCrack({ a: 0.02, W: 1, Pm: 100 * MPa, ...mat });
  close(r.Kr, (100e6 * centreCrackK({ a: 0.02, W: 1 }).K) / 150e6, 1e-12);
  close(r.Lr, 100 / (1 - 0.04) / 355, 1e-12);
  const lbb = leakBeforeBreak({ c: 0.03, t: 0.02, W: 2, Pm: 120 * MPa, ...mat });
  close(lbb.throughLength, 0.08, 1e-12);
  assert.equal(lbb.leakBeforeBreak, true);
  assert.ok(lbb.margin > 1);
  // The critical through-crack is on the assessment line.
  const atCrit = assessThroughCrack({ a: lbb.criticalThroughLength / 2, W: 2, Pm: 120 * MPa, ...mat });
  close(Math.min(atCrit.Kr / atCrit.KrAllowable, atCrit.Lr / atCrit.LrMax) <= 1 ? 1 : 0, 1, 1e-12);
  close(Math.max(atCrit.Kr / atCrit.KrAllowable, atCrit.Lr / atCrit.LrMax), 1, 1e-4);
  const brittle = leakBeforeBreak({ c: 0.03, t: 0.02, W: 2, Pm: 120 * MPa, ...mat, toughness: 20e6 });
  assert.equal(brittle.leakBeforeBreak, false);
  assert.ok(brittle.margin < 1);
});

test('toughness conversions: CTOD (BS 7910:2019 strength-ratio factor) and Charpy lower bound (Annex J)', () => {
  // m = 1.517·0.8^−0.3188 = 1.62885; K = √(m·σy·δ·E/(1 − ν²))
  const r = ctodToK({ ctod: 0.2e-3, yieldStrength: 400 * MPa, tensileStrength: 500 * MPa, E: 207e9, nu: 0.3 });
  close(r.m, 1.62885, 1e-5);
  close(r.K, Math.sqrt((1.62885 * 400e6 * 0.2e-3 * 207e9) / 0.91), 1e-5);
  close(r.K, 172.2e6, 1e-3);
  close(kToCtod({ K: r.K, yieldStrength: 400 * MPa, tensileStrength: 500 * MPa }), 0.2e-3, 1e-10);
  assert.equal(ctodToK({ ctod: 1e-4, yieldStrength: 400 * MPa }).m, 1.5);
  // Charpy 27 J at 25 mm: (12·√27 − 20)·1 + 20 = 62.354 MPa√m; at 50 mm: 42.354·0.5^0.25 + 20 = 55.615
  close(charpyToK({ charpyEnergy: 27 }), 62.354e6, 1e-5);
  close(charpyToK({ charpyEnergy: 27, thickness: 0.05 }), (42.3538 * 0.5 ** 0.25 + 20) * 1e6, 1e-5);
  close(charpyToK({ charpyEnergy: 100, regime: 'upper-shelf' }), 109e6, 1e-12);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => newmanRaju({ a: 0.03, c: 0.01, t: 0.02 }), /a must not exceed t/);
  assert.throws(() => newmanRaju({ a: 0.01, c: 0.004, t: 0.02 }), /a\/c/);
  assert.throws(() => centreCrackK({ a: 1, W: 1 }), RangeError);
  assert.throws(() => parisLawBs7910({ environment: 'air', curve: 'mean', rRatio: 'medium' }), RangeError);
  assert.throws(() => fadOption1(-1, { yieldStrength: 1, tensileStrength: 2 }), /Lr/);
});

test('growth threshold with an equivalent constant-amplitude range: arrested by default, grows with variableAmplitude', () => {
  // ΔK0 = 63 N/mm^(3/2) = 63 × √1000 × 1e3 Pa·√m = 1.992 MPa·√m. A 1 mm deep, 10 mm long surface flaw under an
  // equivalent range of 20 MPa has ΔK of about 20 × 1.1 × √(π × 0.001) ≈ 1.2 MPa·√m < ΔK0: no growth at all.
  const law = parisLawBs7910({ environment: 'air' });
  close(law.threshold, 63 * Math.sqrt(1000) * 1e3, 1e-12);
  assert.equal(law.variableAmplitude, false);
  const flaw = { a0: 1e-3, c0: 5e-3, t: 0.025, membraneRange: 20 * MPa, maxCycles: 1e7 };
  const ca = surfaceCrackGrowth({ ...flaw, law });
  assert.equal(ca.status, 'arrested'); assert.equal(ca.cycles, Infinity); assert.equal(ca.a, 1e-3);
  // BS 7910 advice for variable-amplitude loading: threshold zero — either on the law or on the integration
  const lawVa = parisLawBs7910({ environment: 'air', variableAmplitude: true });
  assert.equal(lawVa.threshold, 0); assert.equal(lawVa.variableAmplitude, true);
  assert.deepEqual(lawVa.stages, law.stages);
  const va1 = surfaceCrackGrowth({ ...flaw, law: lawVa });
  const va2 = surfaceCrackGrowth({ ...flaw, law, variableAmplitude: true });
  assert.equal(va1.status, 'max-cycles'); assert.ok(va1.a > 1e-3);
  close(va2.a, va1.a, 1e-12); close(va2.c, va1.c, 1e-12);
  assert.equal(law.threshold > 0, true, 'the caller’s law object is not modified');
  // an explicit threshold takes precedence over the flag
  close(parisLawBs7910({ variableAmplitude: true, threshold: 1e6 }).threshold, 1e6, 1e-15);
  // one-parameter integration against the closed form (Y = 1.12, A m = 3): N = (a0^-½ − af^-½)/(½ A Y³ Δσ³ π^1.5)
  const A = law.stages[0].A; const Y = 1.12; const ds = 20 * MPa;
  const dK = (a) => Y * ds * Math.sqrt(Math.PI * a);
  assert.equal(parisIntegrate({ a0: 1e-3, af: 2e-3, deltaK: dK, law }).arrested, true);
  const n = parisIntegrate({ a0: 1e-3, af: 2e-3, deltaK: dK, law, variableAmplitude: true });
  assert.equal(n.arrested, false);
  close(n.cycles, (1e-3 ** -0.5 - 2e-3 ** -0.5) / (0.5 * A * Y ** 3 * ds ** 3 * Math.PI ** 1.5), 1e-4);
  close(n.cycles, parisLifeClosedForm({ a0: 1e-3, af: 2e-3, stressRange: ds, Y, A, m: 3 }), 1e-4);
});
