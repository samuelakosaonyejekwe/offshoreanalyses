import test from 'node:test';
import assert from 'node:assert/strict';
import * as S from '../../src/engine/pipeline/span.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${a} vs ${b}`);

const Ds = 0.3239; const t = 0.0127; const E = 207e9;
const I = Math.PI / 64 * (Ds ** 4 - (Ds - 2 * t) ** 4);
const EI = E * I;
const base = { D: 0.41, Ds, t, EI, me: 300 };

test('pinned-pinned span without axial force: f1 = (pi/2) sqrt(EI/(m L^4)) (Blevins Table 8-1)', () => {
  const L = 40;
  const r = S.spanResponse({ ...base, L, boundary: 'pinned' });
  const exact = (Math.PI / 2) * Math.sqrt(EI / (300 * L ** 4));
  close(r.f1IL, exact, 1e-3);                  // C1 = 1.57 against pi/2 = 1.5708
  close(r.f1CF, exact, 1e-3);
  close(r.PcrCF, Math.PI ** 2 * EI / L ** 2, 1e-12);       // Euler load
  close(r.LeffCF, L, 1e-12);
});

test('fixed-fixed span: f1 = 3.56 sqrt(EI/(m L^4)), Pcr = 4 pi^2 EI/L^2, static moment q L^2/12', () => {
  const L = 40;
  const r = S.spanResponse({ ...base, L, boundary: 'fixed', q: 1200 });
  // exact clamped-clamped coefficient: (4.73004^2)/(2 pi) = 3.5608
  close(r.f1IL, (4.73004 ** 2 / (2 * Math.PI)) * Math.sqrt(EI / (300 * L ** 4)), 1e-3);
  close(r.PcrIL, 4 * Math.PI ** 2 * EI / L ** 2, 1e-12);
  close(r.momentStatic, 1200 * 1600 / 12, 1e-12);
  close(r.deflection, 1200 * L ** 4 / (384 * EI), 1e-12);
  close(r.stressStatic, r.momentStatic * Ds / (2 * I), 1e-12);
});

test('axial force and sag terms: f ~ sqrt(1 + Seff/Pcr + C3 (delta/D)^2); static response amplification', () => {
  const L = 40;
  const r0 = S.spanResponse({ ...base, L, boundary: 'pinned' });
  const Pcr = r0.PcrCF;
  const rc = S.spanResponse({ ...base, L, boundary: 'pinned', Seff: -0.75 * Pcr });
  close(rc.f1IL, 0.5 * r0.f1IL, 1e-12);
  const rt = S.spanResponse({ ...base, L, boundary: 'pinned', Seff: 3 * Pcr });
  close(rt.f1IL, 2 * r0.f1IL, 1e-12);
  // static: pinned beam 5 q L^4/(384 EI) and q L^2/8, amplified by 1/(1 + Seff/Pcr)
  const rs = S.spanResponse({ ...base, L, boundary: 'pinned', q: 1000, Seff: -0.5 * Pcr });
  close(rs.deflection, 2 * 5 * 1000 * L ** 4 / (384 * EI), 1e-12);
  close(rs.momentStatic, 2 * 1000 * L * L / 8, 1e-12);
  // sag stiffens the cross-flow mode only
  close(rs.f1CF / rs.f1IL, Math.sqrt((0.5 + 0.8 * (rs.deflection / 0.41) ** 2) / 0.5), 1e-12);
  assert.throws(() => S.spanResponse({ ...base, L, boundary: 'pinned', Seff: -1.01 * Pcr }), RangeError);
  // concrete stiffening
  const csf = S.concreteStiffnessFactor({ EIconc: 0.5 * EI, EIsteel: EI, kc: 0.25 });
  close(csf, 0.25 * 0.5 ** 0.75, 1e-12);
  close(S.spanResponse({ ...base, L, boundary: 'pinned', CSF: csf }).f1IL, r0.f1IL * Math.sqrt(1 + csf), 1e-12);
});

test('effective span length: both branches, continuity at beta = 2.7, rigid-shoulder limit', () => {
  // beta = 2.7: 4.73/(-0.066*7.29 + 1.02*2.7 + 0.63) = 4.73/2.90286 = 1.62943 ; 4.73/(0.036*7.29 + 0.61*2.7 + 1) = 4.73/2.90944 = 1.62574
  const K = (b, L = 40) => 10 ** b * EI / L ** 4;
  close(S.effectiveSpanLength({ L: 40, K: K(2.7), EI }).ratio, 1.62943, 1e-4);
  close(S.effectiveSpanLength({ L: 40, K: K(2.6999999), EI }).ratio, 1.62574, 1e-4);
  // beta = 4: 4.73/(-1.056 + 4.08 + 0.63) = 4.73/3.654 = 1.29447
  const e = S.effectiveSpanLength({ L: 40, K: K(4), EI });
  close(e.beta, 4, 1e-12);
  close(e.ratio, 1.29447, 1e-5);
  close(e.Leff, 40 * e.ratio, 1e-12);
  // stiff shoulders approach the fixed-fixed span
  assert.ok(S.effectiveSpanLength({ L: 40, K: K(7.5), EI }).ratio < 1.05);
  const bc = S.boundaryCoefficients('seabed', 1.2);
  close(bc.C4, 14.1 / 1.44, 1e-12);
  close(bc.C5, 1 / (18 * 1.44 - 6), 1e-12);
  assert.deepEqual([S.boundaryCoefficients('pinned').C4, S.boundaryCoefficients('pinned').C6], [4.93, 5 / 384]);
});

test('single span on seabed: softer soil lengthens the effective span and lowers the frequency', () => {
  const a = S.spanResponse({ ...base, L: 40, q: 1200, KV: 2e7, KL: 1.5e7, KVS: 5e5 });
  assert.ok(a.LeffStatic > a.LeffIL && a.LeffIL > a.LeffCF && a.LeffCF > 40);
  const fixed = S.spanResponse({ ...base, L: 40, q: 1200, boundary: 'fixed' });
  assert.ok(a.f1IL < fixed.f1IL);
  // unit stress amplitude A = C4 (1 + CSF) D (Ds - t) E / Leff^2
  close(a.ACF, (14.1 / (a.LeffCF / 40) ** 2) * 0.41 * (Ds - t) * E / a.LeffCF ** 2, 1e-12);
});

test('added mass coefficient and stability parameter', () => {
  close(S.addedMassCoefficient(0), 2.28, 1e-12);
  close(S.addedMassCoefficient(0.4), 0.68 + 1.6 / 3, 1e-12);
  assert.equal(S.addedMassCoefficient(0.8), 1);
  // Ks = 4 pi me zeta/(rho D^2) = 4 pi*300*0.01/(1025*0.1681) = 0.21880
  const k = S.stabilityParameter({ me: 300, zetaT: 0.01, D: 0.41 });
  close(k.Ks, 0.21880, 1e-4);
  close(k.Ksd, k.Ks / 1.15, 1e-12);
});

test('in-line response model: construction points (DNV-RP-F105 Fig. 4-1), hand calculation', () => {
  // Ksd = 0, gamma_on = 1.1: onset 0.9091 ; AY1 = 0.18 ; VR1 = 1.8 + 0.9091 = 2.7091 ; AY2 = 0.13 ; VRend = 4.5 ; VR2 = 4.24
  const p = S.inlineResponse({ VR: 2.7091, Ksd: 0 });
  close(p.VRonset, 1 / 1.1, 1e-12); close(p.VR1, 1.8 + 1 / 1.1, 1e-12); close(p.VR2, 4.24, 1e-12); close(p.VRend, 4.5, 1e-12);
  close(p.AyOverD, 0.18, 1e-4);
  close(S.inlineResponse({ VR: 4.24, Ksd: 0 }).AyOverD, 0.13, 1e-9);
  close(S.inlineResponse({ VR: 4.37, Ksd: 0 }).AyOverD, 0.065, 1e-9);
  close(S.inlineResponse({ VR: (1 / 1.1 + p.VR1) / 2, Ksd: 0 }).AyOverD, 0.09, 1e-9);
  assert.equal(S.inlineResponse({ VR: 0.5, Ksd: 0 }).AyOverD, 0);
  assert.equal(S.inlineResponse({ VR: 5, Ksd: 0 }).AyOverD, 0);
  // Ksd = 1.0: onset (0.6 + 1)/1.1 ; AY1 = max(0.18*(1 - 1/1.2), 0.13*(1 - 1/1.8)) = max(0.03, 0.05778) ; VRend = 3.7
  const q = S.inlineResponse({ VR: 3, Ksd: 1.0 });
  close(q.VRonset, 1.6 / 1.1, 1e-12); close(q.AY1, 0.13 * (1 - 1 / 1.8), 1e-12); close(q.VRend, 3.7, 1e-12);
  close(S.inlineResponse({ VR: 3, Ksd: 2 }).VRonset, 2.2 / 1.1, 1e-12);
  // current-ratio factor
  close(S.inlineResponse({ VR: 3, Ksd: 0, alpha: 0.65 }).psiAlpha, 0.5, 1e-12);
  assert.equal(S.inlineResponse({ VR: 3, Ksd: 0, alpha: 0.4 }).psiAlpha, 0);
  const tr = S.inlineTurbulenceReduction({ Ic: 0.08, thetaRel: Math.PI / 2 });
  close(tr.R2, 1 - 0.05 / 0.17, 1e-12);
  assert.ok(tr.R1 >= 0 && tr.R1 <= 1);
});

test('cross-flow response model: construction points (DNV-RP-F105 Fig. 4-2), hand calculation', () => {
  // current dominated, away from the bed, heavy pipe: onset = 3*(7/6)/1.2 = 2.9167 ; AZ1 = 0.9
  // VR1 = 7 - (7 - 2.9167)/1.15*(1.3 - 0.9) = 5.5797 ; VR2 = 16 - 7/1.3*0.9 = 11.1538
  const p = S.crossflowResponse({ VR: 8 });
  close(p.VRonset, 3.5 / 1.2, 1e-12); close(p.VR1, 5.5797, 1e-4); close(p.VR2, 11.1538, 1e-5);
  close(p.AzOverD, 0.9, 1e-12);
  close(S.crossflowResponse({ VR: 3.5 / 1.2 }).AzOverD, 0.15, 1e-9);
  close(S.crossflowResponse({ VR: (16 + p.VR2) / 2 }).AzOverD, 0.45, 1e-9);
  assert.equal(S.crossflowResponse({ VR: 1.9 }).AzOverD, 0);
  assert.equal(S.crossflowResponse({ VR: 17 }).AzOverD, 0);
  // onset corrections: gap e/D = 0.4 -> (4 + 0.5)/5 = 0.9 ; trench 0.5 -> 1.25 ; mass ratio 1.2 -> 0.9 ; alpha 0.3 -> 1.1
  close(S.crossflowResponse({ VR: 5, gapRatio: 0.4, trenchRatio: 0.5, massRatio: 1.2, alpha: 0.3, KC: 20 }).VRonset, 3 * 0.9 * 1.25 * 0.9 * 1.1 / 1.2, 1e-12);
  close(S.crossflowResponse({ VR: 8, alpha: 0.3, KC: 20 }).AZ1, 0.8, 1e-12);
  close(S.crossflowResponse({ VR: 8, fRatio: 1.9 }).AZ1, 1.1, 1e-12);
  // damping reduction
  close(S.crossflowResponse({ VR: 8, Ksd: 2 }).Rk, 0.7, 1e-12);
  close(S.crossflowResponse({ VR: 8, Ksd: 9 }).Rk, 3.2 / 27, 1e-12);
});

test('screening criteria: hand calculation', () => {
  // IL: required f = 0.6/(0.909*0.41)*(1 - (40/0.41)/250)/1 = 1.6099*0.60976 = 0.98165 Hz ; available 1.2/1.4 = 0.85714
  const r = S.vivScreening({ f1IL: 1.2, f1CF: 1.3, Uc100: 0.6, D: 0.41, L: 40, VRonsetIL: 1 / 1.1, VRonsetCF: 3.5 / 1.2 });
  close(r.ucIL, (0.6 / ((1 / 1.1) * 0.41)) * (1 - 40 / 0.41 / 250) / (1.2 / 1.4), 1e-12);
  close(r.ucIL, 0.98165 / 0.85714, 1e-4);
  assert.equal(r.ilOk, false);
  // CF: required (0.6)/(2.9167*0.41) = 0.50174 Hz ; available 1.3/1.4 = 0.92857
  close(r.ucCF, 0.50174 / 0.92857, 1e-4);
  assert.equal(r.cfOk, true);
  close(S.vivScreening({ f1IL: 1.2, f1CF: 1.3, Uc100: 0.2, Uw1: 0.8, D: 0.41, L: 40, VRonsetIL: 1, VRonsetCF: 3 }).alphaBar, 0.6, 1e-12);
});

test('VIV fatigue: a near-deterministic current reproduces n/N by hand', () => {
  const r = S.spanResponse({ ...base, L: 40, boundary: 'pinned' });
  const Ksd = 0.2; const U0 = 0.5;
  // Weibull with a very large shape is concentrated at the scale value
  const sn = (Sr) => 10 ** (11.61 - 3 * Math.log10(Sr / 1e6));
  const f = S.vivFatigue({ f1IL: r.f1IL, f1CF: r.f1CF, AIL: r.AIL, ACF: r.ACF, D: 0.41, Ksd, current: { scale: U0, shape: 400 }, sn, nBins: 2000 });
  const il = S.inlineResponse({ VR: U0 / (r.f1IL * 0.41), Ksd });
  const SIL = 2 * r.AIL * il.AyOverD * 1.3;
  assert.ok(SIL > 0);
  const year = 365.25 * 24 * 3600;
  close(f.damagePerYearIL, r.f1IL * year / sn(SIL), 0.03);
  close(f.lifeILYears, 1 / f.damagePerYearIL, 1e-12);
  // { logA, m } form is equivalent to the function form
  const g = S.vivFatigue({ f1IL: r.f1IL, f1CF: r.f1CF, AIL: r.AIL, ACF: r.ACF, D: 0.41, Ksd, current: { scale: U0, shape: 400 }, sn: { logA: 11.61, m: 3 }, nBins: 2000 });
  close(g.damagePerYearIL, f.damagePerYearIL, 1e-12);
  // below onset there is no damage
  const none = S.vivFatigue({ f1IL: 5, f1CF: 5, AIL: r.AIL, ACF: r.ACF, D: 0.41, Ksd, current: { scale: 0.2, shape: 2 }, sn });
  assert.equal(none.lifeYears, Infinity);
  assert.throws(() => S.vivFatigue({ f1IL: 1, f1CF: 1, AIL: 1, ACF: 1, D: 0.41, Ksd, current: { scale: 0.2, shape: 2 }, sn, seaStates: [{ prob: 0.5, Uw: 0, Tu: 10 }] }), RangeError);
});

test('maximum span length: static stress and VIV-onset criteria are met with equality', () => {
  const span = { ...base, q: 1200, boundary: 'pinned' };
  // static, pinned: sigma = q L^2/8 * Ds/(2I) = allow -> L = sqrt(16 I allow/(q Ds))
  const a = S.maxSpanLength({ criterion: 'static', span, sigmaAllow: 200e6 });
  close(a.L, Math.sqrt(16 * I * 200e6 / (1200 * Ds)), 1e-9);
  // onset, pinned, no sag effect in-line: f1 = 1.57 sqrt(EI/(m L^4)) = U/(VRonset D)
  const fReq = 0.8 / (0.909 * 0.41);
  const b = S.maxSpanLength({ criterion: 'onset-il', span, U: 0.8, VRonset: 0.909 });
  close(b.L, Math.sqrt(1.57 * Math.sqrt(EI / 300) / fReq), 1e-9);
  close(b.response.f1IL, fReq, 1e-9);
  const c = S.maxSpanLength({ criterion: 'onset-cf', span, U: 0.8, VRonset: 2.9, gammaF: 1.2 });
  assert.ok(c.L > b.L);
  const u = S.spanUls({ response: b.response, D: 0.41, Uextreme: 1.0, sigmaAllow: 300e6 });
  close(u.sigmaTotal, u.sigmaStatic + u.sigmaEnv, 1e-12);
  assert.ok(u.SIL > 0 && u.uc > 0);
});
