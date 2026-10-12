import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  ec3ImperfectionFactor, ec3Chi, ec3BucklingCurvesI, ec3Classify, elasticCriticalMoment, ec3CmLinear, ec3MemberCheck,
  aiscCompression, aiscFlexureF2, aiscMemberCheck, deflectionCheck, BOLT_GRADES, BOLT_STRESS_AREA, boltResistance,
  boltGroupElastic, weldCorrelationFactor, filletWeldResistance, weldGroupElastic, padeyeCheck, trunnionCheck,
} from '../../src/engine/structural/steelcheck.js';
import { iSectionByName, tubeSection, weldedISection } from '../../src/engine/structural/sections.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);
const uc = (r, name) => r.checks.find((c) => c.name === name).uc;

/* US customary conversions for the AISC worked examples */
const IN = 0.0254; const KSI = 6.894757293168e6; const KIP = 4448.2216152605;

test('EN 1993-1-1 Fig. 6.4 / ECCS tables: χ at λ̄ = 1.0 for curves a0, a, b, c, d', () => {
  // Φ = 0.5[1 + α(λ̄ − 0.2) + λ̄²]; tabulated χ(1.0): a0 0.7253, a 0.6656, b 0.5970, c 0.5399, d 0.4671
  const table = { a0: 0.7253, a: 0.6656, b: 0.5970, c: 0.5399, d: 0.4671 };
  for (const [curve, chi] of Object.entries(table)) assert.ok(Math.abs(ec3Chi(1.0, curve) - chi) < 6e-5, `curve ${curve}`);
  // tabulated χ(2.0): a 0.2229, c 0.1962
  assert.ok(Math.abs(ec3Chi(2.0, 'a') - 0.2229) < 6e-5);
  assert.ok(Math.abs(ec3Chi(2.0, 'c') - 0.1962) < 6e-5);
  assert.equal(ec3Chi(0.2, 'd'), 1);
  assert.equal(ec3Chi(0.05, 'a'), 1);
  assert.deepEqual(['a0', 'a', 'b', 'c', 'd'].map(ec3ImperfectionFactor), [0.13, 0.21, 0.34, 0.49, 0.76]);
  // every curve tends to the Euler hyperbola from below at large slenderness
  for (const c of ['a0', 'd']) { const x = ec3Chi(6, c); assert.ok(x < 1 / 36 && x > 0.8 / 36); }
  assert.throws(() => ec3Chi(1, 'e'), RangeError);
  assert.throws(() => ec3Chi(-1, 'a'), RangeError);
});

test('EN 1993-1-1 Table 6.2 / 6.4: buckling curve selection for I sections', () => {
  assert.deepEqual(ec3BucklingCurvesI({ h: 0.3, b: 0.15, tf: 0.0107 }), { y: 'a', z: 'b', lt: 'a' }); // IPE 300
  assert.deepEqual(ec3BucklingCurvesI({ h: 0.3, b: 0.3, tf: 0.019 }), { y: 'b', z: 'c', lt: 'a' }); // HEB 300
  assert.deepEqual(ec3BucklingCurvesI({ h: 0.6, b: 0.22, tf: 0.019 }), { y: 'a', z: 'b', lt: 'b' }); // IPE 600, h/b > 2
  assert.deepEqual(ec3BucklingCurvesI({ h: 0.3, b: 0.15, tf: 0.0107, fy: 460e6 }), { y: 'a0', z: 'a0', lt: 'a' });
  assert.deepEqual(ec3BucklingCurvesI({ h: 1.0, b: 0.4, tf: 0.03, welded: true }), { y: 'b', z: 'c', lt: 'd' });
  assert.deepEqual(ec3BucklingCurvesI({ h: 0.6, b: 0.4, tf: 0.05, welded: true }), { y: 'c', z: 'd', lt: 'c' });
});

test('EN 1993-1-1 Table 5.2: classification of IPE 300 (bending class 1, compression class 2 in S235; web class 4 in S355)', () => {
  const s = iSectionByName('IPE 300');
  // c/t flange = (150 − 7.1 − 30)/2/10.7 = 5.28; web = (300 − 21.4 − 30)/7.1 = 35.01
  const b = ec3Classify({ section: s, fy: 235e6 });
  close(b.flangeRatio, 5.276, 1e-3); close(b.webRatio, 35.01, 1e-3);
  assert.equal(b.class, 1);
  const c = ec3Classify({ section: s, fy: 235e6, NEd: Infinity });
  assert.equal(c.webClass, 2); // 33 < 35.01 ≤ 38
  const c355 = ec3Classify({ section: s, fy: 355e6, NEd: Infinity });
  assert.equal(c355.webClass, 4); // 42 ε = 34.2 < 35.01
  // combined: NEd = 300 kN, S235 → α = ½(1 + 300e3/(235e6·0.0071·0.2486)) = 0.8617 → class 1 limit 396/(13α − 1) = 38.8
  const m = ec3Classify({ section: s, fy: 235e6, NEd: 300e3 });
  assert.equal(m.webClass, 1);
  // tubes: d/t limits 50ε², 70ε², 90ε²
  assert.equal(ec3Classify({ section: tubeSection({ D: 0.5, t: 0.0125 }), fy: 235e6 }).class, 1); // 40
  assert.equal(ec3Classify({ section: tubeSection({ D: 0.5, t: 0.0125 }), fy: 355e6 }).class, 2); // 40 > 33.1, ≤ 46.3
  assert.equal(ec3Classify({ section: tubeSection({ D: 1.0, t: 0.0125 }), fy: 355e6 }).class, 4); // 80 > 59.6
});

test('Elastic critical moment: IPE 300, L = 5 m, hand calculation with tabulated properties', () => {
  // Iz = 603.8 cm⁴, It = 20.12 cm⁴, Iw = 125.9e3 cm⁶, E = 210 GPa, G = 80.77 GPa
  // π²EIz/L² = π²·210e9·603.8e-8/25 = 500 577 N
  // Iw/Iz = 125.9e-9/603.8e-8 = 0.020851 m²; G It/(π²EIz/L²) = 80.77e9·20.12e-8/500 577 = 0.032465 m²
  // Mcr = 500 577·√(0.053316) = 115.58 kN·m
  const Mcr = elasticCriticalMoment({ E: 210e9, G: 80.77e9, Iz: 603.8e-8, It: 20.12e-8, Iw: 125.9e-9, L: 5 });
  close(Mcr, 115.58e3, 5e-4);
  close(elasticCriticalMoment({ E: 210e9, G: 80.77e9, Iz: 603.8e-8, It: 20.12e-8, Iw: 125.9e-9, L: 5, C1: 1.127 }), 1.127 * Mcr, 1e-12);
  // warping fixed at both ends (kw = 0.5) raises only the warping term
  const fixedW = elasticCriticalMoment({ E: 210e9, G: 80.77e9, Iz: 603.8e-8, It: 20.12e-8, Iw: 125.9e-9, L: 5, kw: 0.5 });
  close(fixedW, 500577 * Math.sqrt(4 * 0.020851 + 0.032465), 5e-4);
  assert.equal(ec3CmLinear(1), 1); assert.equal(ec3CmLinear(-1), 0.4); close(ec3CmLinear(0), 0.6, 1e-15);
});

test('EN 1993-1-1 §6.3.1: IPE 300 S235 column, Lcr = 3.5 m, hand calculation', () => {
  const s = iSectionByName('IPE 300');
  // section generator against the tabulated profile: A = 53.81 cm², Iy = 8356 cm⁴, Iz = 603.8 cm⁴, Wpl,y = 628.4 cm³
  close(s.A, 53.81e-4, 2e-3); close(s.Iy, 8356e-8, 2e-3); close(s.Iz, 603.8e-8, 2e-3); close(s.Wply, 628.4e-6, 2e-3);
  const fy = 235e6; const E = 210e9; const L = 3.5;
  const r = ec3MemberCheck({ section: s, fy, E, NEd: -400e3, Lcry: L });
  // λ̄z = (L/iz)/(π√(E/fy)) with iz = √(Iz/A); curve b about z: χz from Eq. (6.49)
  const lamZ = L / Math.sqrt(s.Iz / s.A) / (Math.PI * Math.sqrt(E / fy));
  close(lamZ, 1.113, 3e-3, 'λ̄z'); // 3500/33.5/93.9
  const phi = 0.5 * (1 + 0.34 * (lamZ - 0.2) + lamZ * lamZ);
  const chiZ = 1 / (phi + Math.sqrt(phi * phi - lamZ * lamZ));
  close(r.factors.lambdaZ, lamZ, 1e-12); close(r.factors.chiZ, chiZ, 1e-12);
  close(r.resistances.NbzRd, chiZ * s.A * fy, 1e-12);
  close(r.resistances.NbzRd, 665e3, 0.01, 'Nb,z,Rd ≈ 0.526·1264.5 kN');
  close(uc(r, 'Flexural buckling about z'), 400e3 / r.resistances.NbzRd, 1e-12);
  assert.equal(r.governing.name, 'Flexural buckling about z');
  assert.deepEqual(r.curves, { y: 'a', z: 'b', lt: 'a' });
  close(r.resistances.NplRd, s.A * fy, 1e-12);
  // tension member: no buckling checks
  const t = ec3MemberCheck({ section: s, fy, NEd: 400e3, Lcry: L });
  assert.ok(!t.checks.some((c) => c.name.startsWith('Flexural buckling')));
  close(t.governing.uc, 400e3 / (s.A * fy), 1e-12);
});

test('EN 1993-1-1 §6.2.5, §6.2.6, §6.2.8, §6.3.2.2: IPE 300 S235 beam, hand calculation', () => {
  const s = iSectionByName('IPE 300');
  const fy = 235e6;
  const r = ec3MemberCheck({ section: s, fy, MyEd: 80e3, VzEd: 100e3, Lcry: 5, LLT: 5 });
  close(r.resistances.MyRd, 628.4e-6 * fy, 2e-3); // 147.7 kN·m
  // Av = A − 2 b tf + (tw + 2r) tf = 5381 − 3210 + 37.1·10.7 = 2568 mm²; Vpl = 2568·235/√3 = 348.4 kN
  close(r.resistances.VplRd, 348.4e3, 3e-3);
  // LTB: λ̄LT = √(Wpl fy/Mcr), rolled h/b = 2 → curve a
  const Mcr = elasticCriticalMoment({ E: 2.1e11, G: 8.1e10, Iz: s.Iz, It: s.J, Iw: s.Iw, L: 5 });
  const lam = Math.sqrt(s.Wply * fy / Mcr);
  const phi = 0.5 * (1 + 0.21 * (lam - 0.2) + lam * lam);
  const chi = 1 / (phi + Math.sqrt(phi * phi - lam * lam));
  close(r.factors.chiLT, chi, 1e-12);
  close(uc(r, 'Lateral-torsional buckling'), 80e3 / (chi * s.Wply * fy), 1e-12);
  // tabulated It = 20.12 cm⁴ gives λ̄LT = 1.13, χLT = 0.58; the thin-walled J of the generator (no fillets) is slightly conservative
  assert.ok(chi > 0.5 && chi < 0.6, `χLT = ${chi}`);
  // high shear: VEd = 0.8 Vpl → ρ = (2·0.8 − 1)² = 0.36, My,V,Rd = [Wpl − ρ Aw²/(4 tw)] fy with Aw = hw tw
  const V = 0.8 * r.resistances.VplRd;
  const hv = ec3MemberCheck({ section: s, fy, MyEd: 80e3, VzEd: V, Lcry: 5, torsionallyRestrained: true });
  const Aw = (0.3 - 2 * 0.0107) * 0.0071;
  close(hv.resistances.MyRd, (s.Wply - 0.36 * Aw * Aw / (4 * 0.0071)) * fy, 1e-12);
  assert.equal(hv.factors.chiLT, 1);
  // γM0 scales the cross-section resistances only
  const g = ec3MemberCheck({ section: s, fy, MyEd: 80e3, Lcry: 5, gammaM0: 1.1 });
  close(g.resistances.MyRd, r.resistances.MyRd / 1.1, 1e-12);
});

test('EN 1993-1-1 Eq. (6.61)/(6.62) with Annex B Table B.1/B.2 factors: HEB 300 S355 beam-column, hand calculation', () => {
  const s = iSectionByName('HEB 300');
  const fy = 355e6; const E = 2.1e11; const L = 4;
  const N = 1500e3; const My = 200e3; const Mz = 40e3;
  const args = { section: s, fy, E, NEd: -N, MyEd: My, MzEd: Mz, Lcry: L, Cmy: 0.9, Cmz: 0.6, CmLT: 0.9 };
  const r = ec3MemberCheck(args);
  assert.ok(r.classification.class <= 2);
  assert.deepEqual(r.curves, { y: 'b', z: 'c', lt: 'a' });
  const { chiY, chiZ, chiLT, lambdaY, lambdaZ } = r.factors;
  const NRk = s.A * fy;
  const ny = N / (chiY * NRk); const nz = N / (chiZ * NRk);
  const kyy = 0.9 * Math.min(1 + (lambdaY - 0.2) * ny, 1 + 0.8 * ny);
  const kzz = 0.6 * Math.min(1 + (2 * lambdaZ - 0.6) * nz, 1 + 1.4 * nz);
  const kyz = 0.6 * kzz;
  assert.ok(lambdaZ >= 0.4);
  const kzy = Math.max(1 - 0.1 * lambdaZ / (0.9 - 0.25) * nz, 1 - 0.1 / (0.9 - 0.25) * nz);
  close(r.factors.kyy, kyy, 1e-12); close(r.factors.kzz, kzz, 1e-12); close(r.factors.kyz, kyz, 1e-12); close(r.factors.kzy, kzy, 1e-12);
  const my = My / (chiLT * s.Wply * fy); const mz = Mz / (s.Wplz * fy);
  close(uc(r, 'Beam-column interaction (y)'), ny + kyy * my + kyz * mz, 1e-12);
  close(uc(r, 'Beam-column interaction (z)'), nz + kzy * my + kzz * mz, 1e-12);
  // Table B.1 (torsionally restrained): kzy = 0.6 kyy and no LTB reduction
  const tr = ec3MemberCheck({ ...args, torsionallyRestrained: true });
  close(tr.factors.kzy, 0.6 * tr.factors.kyy, 1e-12);
  assert.equal(tr.factors.chiLT, 1);
  assert.ok(uc(r, 'Beam-column interaction (z)') > nz);
  // linear cross-section sum §6.2.1(7)
  close(uc(r, 'Cross-section interaction (linear sum)'), N / NRk + My / (s.Wply * fy) + Mz / (s.Wplz * fy), 1e-12);
  // class 3 welded section uses elastic moduli and the class 3 factors
  const w = weldedISection({ bft: 0.4, tft: 0.016, hw: 0.6, tw: 0.012 });
  const c3 = ec3MemberCheck({ section: w, fy, NEd: -200e3, MyEd: 100e3, Lcry: 6 });
  if (c3.classification.class === 3) {
    close(c3.resistances.MyRd, w.Wely * fy, 1e-12);
    close(c3.factors.kyz, c3.factors.kzz, 1e-12);
  }
  assert.throws(() => ec3MemberCheck({ section: s, fy, Lcry: -1 }), RangeError);
});

test('AISC 360 §E3: W14×132, Lc = 30 ft, Fy = 50 ksi → Fe = 31.2 ksi, Fcr = 25.6 ksi, φcPn = 893 kips (AISC Design Examples E.1A)', () => {
  // A = 38.8 in², ry = 3.76 in → Lc/r = 95.74
  const r = aiscCompression({ A: 38.8 * IN * IN, r: 3.76 * IN, Lc: 360 * IN, fy: 50 * KSI, E: 29000 * KSI });
  close(r.slenderness, 95.74, 1e-3);
  close(r.Fe / KSI, 31.2, 2e-3);
  close(r.Fcr / KSI, 25.6, 2e-3);
  close(r.phiPn / KIP, 893, 2e-3);
  close(r.PnOmega / KIP, 993 / 1.67, 2e-3);
  // elastic branch: Lc/r = 150 > 4.71√(E/Fy) = 113.4 → Fcr = 0.877 Fe = 0.877·π²·29000/150² = 11.16 ksi
  const e = aiscCompression({ A: 1e-2, r: 0.05, Lc: 7.5, fy: 50 * KSI, E: 29000 * KSI });
  close(e.Fcr / KSI, 0.877 * Math.PI ** 2 * 29000 / 22500, 1e-12);
  close(e.Fcr / KSI, 11.16, 1e-3);
  // the two branches meet at Fy/Fe = 2.25: 0.658^2.25 = 0.390 and 0.877/2.25 = 0.390
  close(0.658 ** 2.25, 0.877 / 2.25, 1e-3);
});

/* W18×50 (AISC Manual Table 1-1): Zx 101 in³, Sx 88.9 in³, ry 1.65 in, Iy 40.1 in⁴, J 1.24 in⁴, Cw 3040 in⁶, d 18.0 in, tf 0.570 in, A 14.7 in² */
const W18x50 = {
  type: 'I', h: 18.0 * IN, b: 7.5 * IN, tw: 0.355 * IN, tf: 0.570 * IN, A: 14.7 * IN ** 2, Wply: 101 * IN ** 3, Wely: 88.9 * IN ** 3,
  iz: 1.65 * IN, iy: 7.38 * IN, Iz: 40.1 * IN ** 4, Iy: 800 * IN ** 4, J: 1.24 * IN ** 4, Iw: 3040 * IN ** 6, Wplz: 16.6 * IN ** 3, Welz: 10.7 * IN ** 3,
};

test('AISC 360 §F2: W18×50, Fy = 50 ksi → Lp = 5.83 ft, Lr = 16.9 ft, φbMp = 379 kip·ft (AISC Manual Table 3-2)', () => {
  const fy = 50 * KSI; const E = 29000 * KSI; const FT = 12 * IN; const KFT = KIP * FT;
  const y = aiscFlexureF2({ section: W18x50, fy, Lb: 0, E });
  close(y.Lp / FT, 5.83, 2e-3); close(y.Lr / FT, 16.9, 5e-3); close(y.phiMn / KFT, 379, 2e-3);
  close(y.rts / IN, 1.98, 3e-3);
  assert.equal(y.regime, 'yielding');
  // inelastic LTB at Lb = 11.67 ft (third points of a 35 ft span), Cb = 1:
  // Mn = Mp − (Mp − 0.7 Fy Sx)(Lb − Lp)/(Lr − Lp) = 5050 − (5050 − 3111.5)·(140.0 − 69.94)/(203.4 − 69.94) = 4032 kip·in
  const i = aiscFlexureF2({ section: W18x50, fy, Lb: 11.67 * FT, E });
  assert.equal(i.regime, 'inelastic LTB');
  close(i.Mn / (KIP * IN), 4032, 2e-3);
  // Cb raises the strength but never above Mp
  close(aiscFlexureF2({ section: W18x50, fy, Lb: 11.67 * FT, E, Cb: 1.14 }).Mn, 1.14 * i.Mn, 1e-12);
  close(aiscFlexureF2({ section: W18x50, fy, Lb: 11.67 * FT, E, Cb: 3 }).Mn, y.Mp, 1e-12);
  // continuity at Lr: the elastic formula F2-4 gives Fcr = 0.7 Fy there (definition of Lr; the rounded constants
  // 1.95 and 6.76 of F2-6 against the exact π/√2·√0.078·… leave a 0.1 % step)
  const el = aiscFlexureF2({ section: W18x50, fy, Lb: y.Lr * (1 + 1e-9), E });
  assert.equal(el.regime, 'elastic LTB');
  close(el.Mn, 0.7 * fy * W18x50.Wely, 3e-3);
  // elastic LTB at 35 ft, hand calculation of F2-4: Lb/rts = 420/1.982
  const q = 35 * FT / y.rts; const jc = W18x50.J / (W18x50.Wely * (W18x50.h - W18x50.tf));
  const Fcr = Math.PI ** 2 * E / (q * q) * Math.sqrt(1 + 0.078 * jc * q * q);
  close(aiscFlexureF2({ section: W18x50, fy, Lb: 35 * FT, E }).Mn, Fcr * W18x50.Wely, 1e-12);
  // q = 211.9: π²·29000/q² = 6.374 ksi, √(1 + 0.078·8.016e-4·44 903) = 1.951 → Fcr = 12.44 ksi
  close(Fcr / KSI, 12.44, 3e-3);
});

test('AISC 360 §H1.1: interaction equations H1-1a and H1-1b, LRFD and ASD, hand calculation', () => {
  const fy = 50 * KSI; const E = 29000 * KSI; const FT = 12 * IN;
  const base = { section: W18x50, fy, E, Lcx: 12 * FT, Lcy: 6 * FT, Lb: 6 * FT };
  const c = aiscCompression({ A: W18x50.A, r: W18x50.iz, Lc: 6 * FT, fy, E });
  const cx = aiscCompression({ A: W18x50.A, r: W18x50.iy, Lc: 12 * FT, fy, E });
  assert.ok(c.Pn < cx.Pn); // weak axis governs: 72/1.65 = 43.6 > 144/7.38 = 19.5
  const Pc = 0.9 * c.Pn;
  const f = aiscFlexureF2({ section: W18x50, fy, Lb: 6 * FT, E });
  const Mcx = 0.9 * f.Mn; const Mcy = 0.9 * Math.min(fy * W18x50.Wplz, 1.6 * fy * W18x50.Welz);
  // Pr/Pc ≥ 0.2 → H1-1a
  const a = aiscMemberCheck({ ...base, Pr: -0.4 * Pc, Mrx: 0.3 * Mcx, Mry: 0.1 * Mcy });
  assert.equal(a.equation, 'AISC 360-16 Eq. H1-1a');
  close(a.uc, 0.4 + 8 / 9 * (0.3 + 0.1), 1e-12);
  close(a.Pc, Pc, 1e-12); close(a.Mcx, Mcx, 1e-12); close(a.Mcy, Mcy, 1e-12);
  // Pr/Pc < 0.2 → H1-1b
  const b = aiscMemberCheck({ ...base, Pr: -0.1 * Pc, Mrx: 0.5 * Mcx, Mry: 0.2 * Mcy });
  assert.equal(b.equation, 'AISC 360-16 Eq. H1-1b');
  close(b.uc, 0.05 + 0.7, 1e-12);
  // tension uses yielding on the gross section, φt = 0.90
  const t = aiscMemberCheck({ ...base, Pr: 100 * KIP });
  close(t.Pc, 0.9 * fy * W18x50.A, 1e-12);
  // ASD: Pn/1.67 and Mn/1.67 → for equal nominal demand the ratio LRFD/ASD capacities is 0.9·1.67 = 1.503
  const asd = aiscMemberCheck({ ...base, Pr: -100 * KIP, method: 'ASD' });
  close(asd.Pc, c.Pn / 1.67, 1e-12); close(asd.Mcx, f.Mn / 1.67, 1e-12);
  assert.throws(() => aiscMemberCheck({ ...base, method: 'X' }), RangeError);
});

test('Deflection limit: 20 mm on 6 m against L/250', () => {
  const r = deflectionCheck({ deflection: -0.02, span: 6 });
  close(r.limit, 0.024, 1e-12); close(r.uc, 20 / 24, 1e-12); close(r.actualRatio, 300, 1e-12);
});

test('EN 1993-1-8 Table 3.4: M20 bolts — Fv,Rd 94.1 kN (8.8), 98.0 kN (10.9), Ft,Rd 141.1 kN (8.8); bearing by hand', () => {
  assert.equal(BOLT_GRADES['8.8'].fub, 800e6); assert.equal(BOLT_GRADES['10.9'].fyb, 900e6);
  assert.equal(BOLT_STRESS_AREA.M20, 245e-6);
  const b88 = boltResistance({ d: 0.02, As: 245e-6, grade: '8.8' });
  close(b88.FvRd, 0.6 * 800e6 * 245e-6 / 1.25, 1e-12); close(b88.FvRd, 94.08e3, 1e-9);
  close(b88.FtRd, 141.12e3, 1e-9);
  const b109 = boltResistance({ d: 0.02, As: 245e-6, grade: 10.9 });
  close(b109.FvRd, 98.0e3, 1e-9); assert.equal(b109.alphaV, 0.5);
  // shank in the shear plane, double shear: 2·0.6·1000·π·20²/4/1.25 = 301.6 kN
  close(boltResistance({ d: 0.02, As: 245e-6, grade: '10.9', threadsInShearPlane: false, shearPlanes: 2 }).FvRd, 301.59e3, 1e-4);
  // bearing: t = 10 mm S355 (fu = 490 MPa), d0 = 22, e1 = 40, p1 = 60, e2 = 30, p2 = 60:
  // αd = min(40/66, 60/66 − 0.25) = 0.6061, αb = min(0.6061, 800/490, 1) = 0.6061;
  // k1 = min(2.8·30/22 − 1.7, 1.4·60/22 − 1.7, 2.5) = 2.118 → Fb,Rd = 2.118·0.6061·490·20·10/1.25 = 100.6 kN
  const br = boltResistance({ d: 0.02, As: 245e-6, grade: '8.8', bearing: { t: 0.01, fu: 490e6, d0: 0.022, e1: 0.04, e2: 0.03, p1: 0.06, p2: 0.06 }, FvEd: 50e3, FtEd: 60e3 });
  close(br.alphaB, 40 / 66, 1e-12); close(br.k1, 2.8 * 30 / 22 - 1.7, 1e-12);
  close(br.FbRd, 100.6e3, 1e-3);
  close(br.ucCombined, 50 / 94.08 + 60 / (1.4 * 141.12), 1e-9);
  close(br.ucBearing, 50e3 / br.FbRd, 1e-12);
  assert.throws(() => boltResistance({ d: 0.02, As: 245e-6, grade: '9.9' }), RangeError);
});

test('Bolt group, polar-moment method: 2×2 group under eccentric shear, hand calculation', () => {
  // bolts at (±50, ±75) mm; V = 100 kN down at x = +200 mm from the centroid
  // Ip = 4(50² + 75²) = 32 500 mm²; M = −100·200 = −20 000 kN·mm (clockwise)
  // critical bolt (+50, ±75): Fy = −25 − 20 000·50/32 500 = −55.77 kN, Fx = ±20 000·75/32 500 = 46.15 kN → 72.39 kN
  const bolts = [[-0.05, -0.075], [0.05, -0.075], [0.05, 0.075], [-0.05, 0.075]];
  const r = boltGroupElastic({ bolts, Vy: -100e3, ref: [0.2, 0] });
  close(r.Ip, 32500e-6, 1e-12);
  close(r.maxShear, 72.39e3, 2e-4);
  // equilibrium of the bolt forces with the applied load
  let fx = 0; let fy = 0; let mz = 0;
  r.forces.forEach((f, i) => { fx += f.Fx; fy += f.Fy; mz += bolts[i][0] * f.Fy - bolts[i][1] * f.Fx; });
  assert.ok(Math.abs(fx) < 1e-6); close(fy, -100e3, 1e-12); close(mz, -20e3, 1e-12);
  // out-of-plane: N = 80 kN and Mx = 6 kN·m → 20 ± 6/(4·0.075) = 20 ± 20 kN
  const t = boltGroupElastic({ bolts, N: 80e3, Mx: 6e3 });
  close(t.maxTension, 40e3, 1e-12);
  assert.ok(Math.abs(t.forces[0].tension) < 1e-6);
});

test('EN 1993-1-8 §4.5.3.3: fillet weld S355 (fu = 510 MPa), a = 4.2 mm → fvw,d = 261.7 MPa, 1.10 kN/mm; weld group by hand', () => {
  assert.equal(weldCorrelationFactor('S355'), 0.9); assert.equal(weldCorrelationFactor('S235'), 0.8);
  const w = filletWeldResistance({ a: 0.0042, fu: 510e6, betaW: 0.9 });
  close(w.fvwd, 510e6 / (Math.sqrt(3) * 0.9 * 1.25), 1e-12); close(w.fvwd, 261.7e6, 2e-4); close(w.FwRd, 1.099e6, 1e-3);
  // two parallel vertical welds 200 mm long, 100 mm apart; V = 150 kN along the welds, torsion 10 kN·m
  // L = 400 mm; Ix = 2·200³/12 = 1.3333e6 mm³; Iy = 2·200·50² = 1.0e6 mm³; Ip = 2.3333e6 mm³
  // corner (50, 100): fy = 150/400 + 10 000·50/2.3333e6 = 0.375 + 0.2143 = 0.5893 kN/mm; fx = −10 000·100/2.3333e6 = −0.4286 → 0.7287 kN/mm
  const g = weldGroupElastic({ segments: [{ x1: -0.05, y1: -0.1, x2: -0.05, y2: 0.1 }, { x1: 0.05, y1: -0.1, x2: 0.05, y2: 0.1 }], Vy: 150e3, Mz: 10e3, FwRd: w.FwRd });
  close(g.length, 0.4, 1e-12); close(g.Ix, 1.33333e-3, 1e-5); close(g.Iy, 1.0e-3, 1e-12); close(g.Ip, 2.33333e-3, 1e-5);
  close(g.maxForcePerLength, 0.7287e6, 2e-4);
  close(g.uc, 0.7287 / 1.099, 2e-3);
  // bending about x: f = M y/Ix = 5 000·100/1.3333e6 = 0.375 kN/mm
  close(weldGroupElastic({ segments: [{ x1: -0.05, y1: -0.1, x2: -0.05, y2: 0.1 }, { x1: 0.05, y1: -0.1, x2: 0.05, y2: 0.1 }], Mx: 5e3 }).maxForcePerLength, 0.375e6, 1e-5);
});

test('Padeye: bearing, tear-out, net section, cheek weld and base section by hand', () => {
  // F = 500 kN, main plate t = 40 mm, R = 150 mm, hole 80 mm, pin 76 mm, cheeks 2×15 mm, Rc = 110 mm, weld throat 8 mm
  const F = 500e3; const fy = 355e6;
  const r = padeyeCheck({ F, fy, t: 0.04, R: 0.15, dHole: 0.08, dPin: 0.076, tCheek: 0.015, rCheek: 0.11, aCheekWeld: 0.008, baseWidth: 0.5, holeHeight: 0.2, angleDeg: 30, outOfPlaneFraction: 0.05 });
  // bearing: 500e3/(76·70) = 93.98 MPa
  close(r.stresses.bearing, 93.985e6, 1e-4);
  // tear-out: 2·(150 − 40)·40 + 4·(110 − 40)·15 = 8800 + 4200 = 13 000 mm² → 38.46 MPa
  close(r.stresses.tearOut, 38.4615e6, 1e-5);
  // net tension: (300 − 80)·40 + 2·(220 − 80)·15 = 8800 + 4200 = 13 000 mm² → 38.46 MPa
  close(r.stresses.netTension, 38.4615e6, 1e-5);
  // cheek weld: 500e3·15/70/(2π·110·8) = 19.38 MPa
  close(r.stresses.cheekWeld, 19.378e6, 1e-4);
  // base: axial 500·cos30/(500·40) = 21.65 MPa; in-plane 250e3·200/(40·500²/6) = 30.0 MPa;
  // out-of-plane 25e3·200/(500·40²/6) = 37.5 MPa; shear √(250² + 25²)/20 000 = 12.56 MPa; σeq = √(89.15² + 3·12.56²) = 91.77 MPa
  close(r.stresses.baseAxial, 21.651e6, 1e-4); close(r.stresses.baseInPlaneBending, 30e6, 1e-9);
  close(r.stresses.baseOutOfPlaneBending, 37.5e6, 1e-9); close(r.stresses.baseShear, 12.562e6, 1e-4);
  close(r.stresses.baseVonMises, 91.77e6, 2e-4);
  // utilisations against the default allowables 0.9, 0.4, 0.45, 0.4, 0.6 Fy
  const u = Object.fromEntries(r.checks.map((c) => [c.name, c.uc]));
  close(u['Pin bearing'], 93.985 / 319.5, 1e-4);
  close(u['Tear-out shear'], 38.4615 / 142, 1e-4);
  close(u['Net-section tension at the hole'], 38.4615 / 159.75, 1e-4);
  close(u['Cheek plate weld'], 19.378 / 142, 1e-4);
  close(u['Base section (von Mises)'], 91.77 / 213, 2e-4);
  assert.equal(r.governing.name, 'Base section (von Mises)');
  // user allowables replace the defaults
  const r2 = padeyeCheck({ F, fy, t: 0.04, R: 0.15, dHole: 0.08, dPin: 0.076, allow: { bearing: 0.45 } });
  close(r2.checks[0].uc, F / (0.076 * 0.04) / (0.45 * fy), 1e-12);
  assert.throws(() => padeyeCheck({ F, fy, t: 0.04, R: 0.03, dHole: 0.08, dPin: 0.076 }), RangeError);
  assert.throws(() => padeyeCheck({ F, fy, t: 0.04, R: 0.15, dHole: 0.08, dPin: 0.09 }), RangeError);
  assert.throws(() => padeyeCheck({ F, fy, t: 0.04, R: 0.15, dHole: 0.08, dPin: 0.076, tCheek: 0.01, rCheek: 0.2 }), RangeError);
});

test('Trunnion: tube 508 × 25.4 mm, F = 1 MN at 150 mm, hand calculation', () => {
  // A = π/4 (508² − 457.2²) = 38 510 mm²; W = π/32 (508⁴ − 457.2⁴)/508 = 4.426e6 mm³
  // σb = 1e6·150/4.426e6 = 33.89 MPa; τ = 2e6/38 510 = 51.93 MPa; σeq = √(33.89² + 3·51.93²) = 96.12 MPa
  const r = trunnionCheck({ F: 1e6, arm: 0.15, D: 0.508, t: 0.0254, fy: 355e6 });
  close(r.stresses.bending, 33.89e6, 1e-3); close(r.stresses.shear, 51.93e6, 1e-3); close(r.stresses.vonMises, 96.12e6, 1e-3);
  assert.equal(r.governing.name, 'Shear'); close(r.governing.uc, 51.93 / 142, 1e-3);
});
