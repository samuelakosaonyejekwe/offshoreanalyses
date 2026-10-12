import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  YEAR, concreteProperties, concreteStress, steelStress, rcSection, circularBars, interactionDiagram, momentCapacity, momentCurvature,
  crackingMoment, crackedSection, effectiveModulus, crackedDeflection, shearResistanceNoLinks, shearResistanceLinks, crackWidth,
  chlorideProfile, chlorideInitiationTime, carbonation, corrosionRate, barSectionLoss, residualCapacityHistory, pileCorrosionLoss,
} from '../../src/engine/structural/concrete.js';
import { erfc } from '../../src/engine/math/special.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);
const MPa = 1e6;

test('EN 1992-1-1 Table 3.1: C30/37 fcm 38, fctm 2.9, Ecm 33 GPa; C90/105 fctm 5.0, Ecm 44 GPa, εc2 2.6 ‰, εcu2 2.6 ‰, n 1.4', () => {
  const c30 = concreteProperties({ fck: 30 * MPa });
  close(c30.fcm, 38 * MPa, 1e-12); close(c30.fctm, 2.9 * MPa, 0.005); close(c30.Ecm, 33e9, 0.01); close(c30.fctk005, 2.0 * MPa, 0.02);
  assert.equal(c30.epsC2, 0.002); assert.equal(c30.epsCu2, 0.0035); assert.equal(c30.n, 2);
  const c50 = concreteProperties({ fck: 50 * MPa });
  close(c50.fctm, 4.1 * MPa, 0.01); close(c50.Ecm, 37e9, 0.01);
  const c90 = concreteProperties({ fck: 90 * MPa });
  close(c90.fctm, 5.0 * MPa, 0.01); close(c90.Ecm, 44e9, 0.01); close(c90.epsC2, 2.6e-3, 0.005); close(c90.epsCu2, 2.6e-3, 1e-12); close(c90.n, 1.4, 1e-12);
  const c60 = concreteProperties({ fck: 60 * MPa });
  close(c60.epsC2, 2.3e-3, 0.01); close(c60.epsCu2, 2.9e-3, 0.01); close(c60.n, 1.6, 0.01); close(c60.fctm, 4.4 * MPa, 0.015); // tabulated values are rounded to two figures
  assert.throws(() => concreteProperties({ fck: 5 * MPa }), RangeError);
});

test('material laws: parabola–rectangle and bilinear steel', () => {
  close(concreteStress({ strain: 0.001, fcd: 20 * MPa }), 15 * MPa, 1e-12); // 1 − (1 − ½)² = ¾
  assert.equal(concreteStress({ strain: 0.002, fcd: 20 * MPa }), 20 * MPa);
  assert.equal(concreteStress({ strain: 0.003, fcd: 20 * MPa }), 20 * MPa);
  assert.equal(concreteStress({ strain: -0.001, fcd: 20 * MPa }), 0);
  close(steelStress({ strain: 0.001, fyd: 434.78 * MPa }), 200 * MPa, 1e-12);
  assert.equal(steelStress({ strain: -0.01, fyd: 434.78 * MPa }), -434.78 * MPa);
  // inclined branch: k = 1.08, εuk = 0.05 → at εuk the stress is 1.08 fyd
  close(steelStress({ strain: 0.05, fyd: 400 * MPa, k: 1.08 }), 432 * MPa, 1e-12);
});

const rect = (As, extra = {}) => rcSection({ shape: 'rect', b: 0.3, h: 0.5, bars: [{ depth: 0.45, area: As, diameter: 0.025 }], ...extra });
const mat = { fck: 30 * MPa, fyk: 500 * MPa };

test('singly reinforced rectangular section: MRd vs the stress-block hand calculation (parabola–rectangle: resultant 17/21 fcd b x at 99/238 x)', () => {
  // b = 300, d = 450, As = 1500 mm², fcd = 20 MPa, fyd = 434.8 MPa → Fs = 652.2 kN
  // parabola–rectangle block (exact integrals for εcu2 = 3.5 ‰, εc2 = 2 ‰): Fc = (17/21) fcd b x, acting at (99/238) x from the top
  // x = 652 174/(0.80952·20·300) = 134.27 mm; z = 450 − 0.41597·134.27 = 394.15 mm; MRd = 257.06 kN·m
  const r = momentCapacity({ section: rect(1500e-6), material: mat });
  const Fs = 1500e-6 * 500e6 / 1.15;
  const x = Fs / (17 / 21 * 20e6 * 0.3);
  const M = Fs * (0.45 - 99 / 238 * x);
  close(r.x, x, 2e-4); close(r.MRd, M, 2e-4); close(r.MRd, 257.06e3, 5e-4);
  assert.ok(r.steelYielded); close(r.epsSteel, 0.0035 * (0.45 - x) / x, 1e-3);
  // EN 1992-1-1 §3.1.7(3) rectangular block (λ = 0.8, η = 1): x = Fs/(0.8 fcd b) = 135.9 mm, M = Fs (d − 0.4x) = 258.1 kN·m — within 0.5 %
  close(r.MRd, Fs * (0.45 - 0.4 * Fs / (0.8 * 20e6 * 0.3)), 5e-3);
  // αcc = 0.85 (e.g. UK national annex): fcd = 17 MPa
  const r85 = momentCapacity({ section: rect(1500e-6), material: { ...mat, alphaCC: 0.85 } });
  const x85 = Fs / (17 / 21 * 17e6 * 0.3);
  close(r85.MRd, Fs * (0.45 - 99 / 238 * x85), 2e-4);
  // over-reinforced: steel does not yield (x/d > 0.617 for B500)
  const over = momentCapacity({ section: rect(6000e-6), material: mat });
  assert.equal(over.steelYielded, false); assert.ok(over.xOverD > 0.617);
});

test('N–M interaction: squash load, pure tension, balanced point by hand; momentCapacity lies on the diagram', () => {
  // symmetric 300 × 500 with 2 × 1000 mm² at 50 and 450 mm
  const s = rcSection({ shape: 'rect', b: 0.3, h: 0.5, bars: [{ depth: 0.05, area: 1000e-6 }, { depth: 0.45, area: 1000e-6 }] });
  const d = interactionDiagram({ section: s, material: mat });
  // uniform εc2 = 2 ‰: concrete fcd·Ac = 3000 kN, steel at Es·0.002 = 400 MPa (< fyd) → 800 kN
  close(d.NRdMax, 20e6 * 0.15 + 2000e-6 * 400e6, 1e-9); assert.ok(Math.abs(d.M[0]) < 1e-6 * d.NRdMax);
  // pure tension: −As fyd
  close(d.NRdMin, -2000e-6 * 500e6 / 1.15, 1e-9);
  // balanced strain state: εcu2 at the top, εyd = 2.174 ‰ at d = 450 → x = 450·3.5/5.674 = 277.6 mm
  // Fc = 17/21·20·300·277.6 = 1348.3 kN at 115.5 mm; top steel ε = 3.5(277.6 − 50)/277.6 = 2.87 ‰ → yields (+434.8 kN); bottom −434.8 kN
  // N = 1348.3 kN; M = 1348.3·(0.25 − 0.1155) + 2·434.8·0.2 = 181.4 + 173.9 = 355.3 kN·m
  const x = 0.45 * 3.5 / (3.5 + 500 / 1.15 / 200e3 * 1e3);
  const Fc = 17 / 21 * 20e6 * 0.3 * x;
  const bal = momentCapacity({ section: s, material: mat, N: Fc });
  close(bal.x, x, 2e-3);
  close(bal.MRd, Fc * (0.25 - 99 / 238 * x) + 2 * 1000e-6 * 500e6 / 1.15 * 0.2, 1e-3);
  close(bal.MRd, 355.3e3, 2e-3);
  // the balanced point is (close to) the maximum moment of the diagram
  close(d.balanced.M, bal.MRd, 5e-3);
  // N decreases monotonically along the diagram; capacity at intermediate N interpolates the diagram
  for (let i = 1; i < d.N.length; i++) assert.ok(d.N[i] < d.N[i - 1] + 1e-6);
  const k = 40;
  close(momentCapacity({ section: s, material: mat, N: d.N[k] }).MRd, d.M[k], 1e-8);
  // pure bending of the doubly reinforced section exceeds the singly reinforced one only slightly (compression steel near the neutral axis)
  const m0 = momentCapacity({ section: s, material: mat }).MRd;
  assert.ok(m0 > 1000e-6 * 434.78e6 * 0.4 && m0 < 1000e-6 * 434.78e6 * 0.45);
  assert.throws(() => momentCapacity({ section: s, material: mat, N: 5e6 }), RangeError);
});

test('moment–curvature: initial cracked stiffness, yield and ultimate points, ductility falls with the steel ratio', () => {
  const s = rect(1500e-6);
  const mk = momentCurvature({ section: s, material: mat, nPoints: 80 });
  const ult = momentCapacity({ section: s, material: mat });
  close(mk.Mu, ult.MRd, 1e-12); close(mk.kappaU, 0.0035 / ult.x, 1e-9);
  // first yield: κy = εyd/(d − x) with x from equilibrium; My ≈ Fs·(d − x/3) to within a few % (near-triangular concrete block)
  const i = mk.kappa.findIndex((k) => k > mk.kappaY);
  assert.ok(i > 1);
  const Fs = 1500e-6 * 500e6 / 1.15;
  assert.ok(mk.My > 0.9 * mk.Mu && mk.My < mk.Mu);
  const xy = 0.45 - (500e6 / 1.15 / 200e9) / mk.kappaY;
  assert.ok(Math.abs(mk.My - Fs * (0.45 - xy / 3)) < 0.03 * mk.My, `My = ${mk.My}`);
  // secant stiffness at low moment equals the cracked-section stiffness with the initial concrete modulus of the parabola 2 fcd/εc2 = 20 GPa
  const cr = crackedSection({ section: s, Ec: 2 * 20e6 / 0.002 });
  close(mk.M[1] / mk.kappa[1], 20e9 * cr.inertia, 0.03);
  // moment is non-decreasing along the curve
  for (let j = 1; j < mk.M.length; j++) assert.ok(mk.M[j] >= mk.M[j - 1] - 1e-6 * mk.Mu);
  assert.ok(mk.ductility > 3);
  assert.ok(momentCurvature({ section: rect(3000e-6), material: mat }).ductility < mk.ductility);
  // axial compression raises the moment capacity below the balance point and reduces the ultimate curvature
  const mkN = momentCurvature({ section: s, material: mat, N: 500e3 });
  assert.ok(mkN.Mu > mk.Mu && mkN.kappaU < mk.kappaU);
});

test('T and circular sections: geometry exact; T-beam with the neutral axis in the flange equals the rectangular beam of the flange width', () => {
  const T = rcSection({ shape: 'T', bf: 1.2, hf: 0.15, bw: 0.3, h: 0.6, bars: [{ depth: 0.54, area: 2000e-6 }] });
  // A = 1.2·0.15 + 0.3·0.45 = 0.315 m²; centroid = (0.18·0.075 + 0.135·0.375)/0.315 = 0.20357 m
  close(T.area, 0.315, 1e-12); close(T.centroid, 0.203571, 1e-5);
  // I = 1.2·0.15³/12 + 0.18·0.12857² + 0.3·0.45³/12 + 0.135·0.17143² = 0.0095585 m⁴
  close(T.inertia, 0.0095585, 1e-4);
  const R = rcSection({ shape: 'rect', b: 1.2, h: 0.6, bars: [{ depth: 0.54, area: 2000e-6 }] });
  const mT = momentCapacity({ section: T, material: mat }); const mR = momentCapacity({ section: R, material: mat });
  assert.ok(mT.x < 0.15); close(mT.x, mR.x, 1e-6);
  // same internal couple; moments are both about the respective centroid with N = 0 so they are equal
  close(mT.MRd, mR.MRd, 1e-6);
  const C = rcSection({ shape: 'circle', D: 0.8, bars: circularBars({ D: 0.8, cover: 0.06, nBars: 12, barDiameter: 0.025 }) });
  close(C.area, Math.PI * 0.16, 1e-10); close(C.centroid, 0.4, 1e-10); close(C.inertia, Math.PI * 0.8 ** 4 / 64, 1e-12);
  close(C.As, 12 * Math.PI * 0.025 ** 2 / 4, 1e-12);
  const d = interactionDiagram({ section: C, material: mat });
  close(d.NRdMax, 20e6 * C.area + C.As * 400e6, 1e-9);
  // circular column: moment capacity grows with moderate compression
  assert.ok(momentCapacity({ section: C, material: mat, N: 2e6 }).MRd > momentCapacity({ section: C, material: mat }).MRd);
});

test('cracked and uncracked elastic sections: x = d[√((αρ)² + 2αρ) − αρ], Icr = b x³/3 + αAs(d − x)², Mcr by hand', () => {
  const s = rect(1500e-6);
  const Ec = 33e9; const a = 200 / 33; const rho = 1500e-6 / (0.3 * 0.45);
  const x = 0.45 * (Math.sqrt((a * rho) ** 2 + 2 * a * rho) - a * rho);
  const cr = crackedSection({ section: s, Ec, M: 100e3 });
  close(cr.x, x, 1e-6); close(cr.x, 0.1378, 2e-3);
  const I = 0.3 * x ** 3 / 3 + a * 1500e-6 * (0.45 - x) ** 2;
  close(cr.inertia, I, 1e-5);
  close(cr.sigmaS, a * 100e3 * (0.45 - x) / I, 1e-5); close(cr.sigmaC, 100e3 * x / I, 1e-5);
  // steel stress also from the lever arm z = d − x/3: σs = M/(As z)
  close(cr.sigmaS, 100e3 / (1500e-6 * (0.45 - x / 3)), 1e-5);
  // uncracked transformed: A = 0.15 + (α − 1)As; yu = (0.15·0.25 + 5.0606·0.0015·0.45)/A
  const un = crackingMoment({ section: s, fctm: 2.9e6, Ec });
  const A = 0.15 + (a - 1) * 1500e-6; const yu = (0.15 * 0.25 + (a - 1) * 1500e-6 * 0.45) / A;
  const Iu = 0.3 * 0.5 ** 3 / 12 + 0.15 * (0.25 - yu) ** 2 + (a - 1) * 1500e-6 * (0.45 - yu) ** 2;
  close(un.centroid, yu, 1e-9); close(un.inertia, Iu, 1e-9); close(un.Mcr, 2.9e6 * Iu / (0.5 - yu), 1e-9);
  // the reinforcement raises Mcr above the plain gross-section value fctm b h²/6 = 36.25 kN·m by 10–15 % here
  assert.ok(un.Mcr > 1.10 * 36.25e3 && un.Mcr < 1.15 * 36.25e3);
  // axial compression delays cracking
  assert.ok(crackingMoment({ section: s, fctm: 2.9e6, Ec, N: 300e3 }).Mcr > un.Mcr);
  close(effectiveModulus({ Ecm: 33e9, phi: 2 }), 11e9, 1e-12);
});

test('EN 1992-1-1 §7.4.3 deflection interpolation: ζ = 1 − β(Mcr/M)²', () => {
  const r = crackedDeflection({ M: 100e3, Mcr: 40e3, EIuncracked: 110e6, EIcracked: 50e6, beta: 0.5, L: 6 });
  close(r.zeta, 1 - 0.5 * 0.16, 1e-12);
  const k = 0.92 * 100e3 / 50e6 + 0.08 * 100e3 / 110e6;
  close(r.curvature, k, 1e-12); close(r.deflection, 5 / 48 * 36 * k, 1e-12);
  // uncracked member: elastic deflection 5 M L²/(48 EI)
  const u = crackedDeflection({ M: 30e3, Mcr: 40e3, EIuncracked: 110e6, EIcracked: 50e6, L: 6 });
  assert.equal(u.zeta, 0); close(u.deflection, 5 * 30e3 * 36 / (48 * 110e6), 1e-12);
});

test('EN 1992-1-1 §6.2.2 shear without links: 300 × 450 (d), fck 30, Asl 1500 mm² → VRd,c = 86.9 kN; minimum governs for light reinforcement', () => {
  // k = 1 + √(200/450) = 1.667; ρl = 1500/(300·450) = 0.01111; 0.12·1.667·(100·0.01111·30)^(1/3) = 0.12·1.667·3.218 = 0.6436 MPa
  const r = shearResistanceNoLinks({ bw: 0.3, d: 0.45, fck: 30e6, Asl: 1500e-6 });
  close(r.k, 1.6667, 1e-4); close(r.rhoL, 0.011111, 1e-4); close(r.VRdc, 86.9e3, 2e-3);
  // vmin = 0.035·1.667^1.5·√30 = 0.4125 MPa → 55.7 kN governs when ρl is small
  const m = shearResistanceNoLinks({ bw: 0.3, d: 0.45, fck: 30e6, Asl: 200e-6 });
  assert.ok(m.governedByMinimum); close(m.VRdc, 55.7e3, 3e-3);
  // axial compression 600 kN on 300 × 500: σcp = 4 MPa (= 0.2 fcd limit) → + 0.15·4·0.135 m² = 81 kN
  const n = shearResistanceNoLinks({ bw: 0.3, d: 0.45, fck: 30e6, Asl: 1500e-6, NEd: 600e3, Ac: 0.15 });
  close(n.VRdc - r.VRdc, 0.15 * 4e6 * 0.3 * 0.45, 1e-9);
  // deep member: k limited... d = 1.0 m → k = 1.447; ρl capped at 0.02
  close(shearResistanceNoLinks({ bw: 0.3, d: 1.0, fck: 30e6, Asl: 0.01 }).rhoL, 0.02, 1e-12);
  assert.equal(shearResistanceNoLinks({ bw: 0.3, d: 0.1, fck: 30e6, Asl: 1e-4 }).k, 2);
});

test('EN 1992-1-1 §6.2.3 shear with links: VRd,s and VRd,max by hand; required links with strut angle limited by web crushing', () => {
  // bw 300, d 450 → z = 405 mm; fck 30: ν1 = 0.6(1 − 30/250) = 0.528, fcd = 20 MPa; links H10 two legs at 200: Asw/s = 785.4 mm²/m
  const r = shearResistanceLinks({ bw: 0.3, d: 0.45, fck: 30e6, fywk: 500e6, AswOverS: 785.4e-6 });
  close(r.nu1, 0.528, 1e-12);
  close(r.VRds, 785.4e-6 * 0.405 * 434.78e6 * 2.5, 1e-4); // 345.7 kN
  close(r.VRdmax, 0.3 * 0.405 * 0.528 * 20e6 / (2.5 + 0.4), 1e-12); // 442.4 kN
  close(r.VRd, r.VRds, 1e-12);
  // cot θ = 1: VRd,max = 0.3·0.405·0.528·20/2 = 641.5 kN
  const r1 = shearResistanceLinks({ bw: 0.3, d: 0.45, fck: 30e6, fywk: 500e6, cotTheta: 1 });
  close(r1.VRdmax, 641.52e3, 1e-4);
  // required links: VEd = 300 kN (< VRd,max at 2.5) → Asw/s = 300e3/(0.405·434.78e6·2.5) = 681.5 mm²/m
  close(shearResistanceLinks({ bw: 0.3, d: 0.45, fck: 30e6, fywk: 500e6, VEd: 300e3 }).required.AswOverS, 681.5e-6, 1e-3);
  // VEd = 550 kN: θ = ½ asin(550/641.52) = 29.5°, cot θ = 1.767 → VRd,max(θ) = VEd
  const q = shearResistanceLinks({ bw: 0.3, d: 0.45, fck: 30e6, fywk: 500e6, VEd: 550e3 }).required;
  close(q.cotTheta, 1 / Math.tan(0.5 * Math.asin(550 / 641.52)), 1e-9);
  close(shearResistanceLinks({ bw: 0.3, d: 0.45, fck: 30e6, fywk: 500e6, cotTheta: q.cotTheta }).VRdmax, 550e3, 1e-9);
  assert.ok(q.adequate);
  assert.equal(shearResistanceLinks({ bw: 0.3, d: 0.45, fck: 30e6, fywk: 500e6, VEd: 700e3 }).required.adequate, false);
});

test('EN 1992-1-1 §7.3.4 crack width: hand calculation for a 300 × 500 beam, 3H25, σs = 250 MPa', () => {
  // h − d = 50 mm, x = 138 mm: hc,ef = min(125, 120.7, 250) = 120.7 mm; ρp,eff = 1472.6/(300·120.7) = 0.04068
  // εsm − εcm = [250 − 0.4·2.9/0.04068·(1 + 6.061·0.04068)]/200 000 = (250 − 35.54)/200 000 = 1.0723e-3 (> 0.6·250/200 000 = 0.75e-3)
  // sr,max = 3.4·37.5 + 0.8·0.5·0.425·25/0.04068 = 127.5 + 104.5 = 232.0 mm → wk = 0.249 mm
  const r = crackWidth({ sigmaS: 250e6, As: 1472.6e-6, b: 0.3, h: 0.5, d: 0.45, x: 0.138, cover: 0.0375, barDiameter: 0.025, fctEff: 2.9e6, Ecm: 33e9 });
  close(r.hcEff, 0.12067, 1e-4); close(r.rhoEff, 0.04068, 1e-3); close(r.strainDifference, 1.0723e-3, 1e-3);
  close(r.srMax, 0.2320, 1e-3); close(r.wk, 0.2488e-3, 2e-3);
  assert.equal(r.minimumStrainGoverns, false);
  // low stress: the 0.6 σs/Es floor governs
  const lo = crackWidth({ sigmaS: 60e6, As: 1472.6e-6, b: 0.3, h: 0.5, d: 0.45, x: 0.138, cover: 0.0375, barDiameter: 0.025, fctEff: 2.9e6, Ecm: 33e9 });
  assert.ok(lo.minimumStrainGoverns); close(lo.strainDifference, 0.6 * 60 / 200e3, 1e-12);
  // widely spaced bars: sr,max = 1.3(h − x)
  close(crackWidth({ sigmaS: 250e6, As: 1472.6e-6, b: 0.3, h: 0.5, d: 0.45, x: 0.138, cover: 0.0375, barDiameter: 0.025, fctEff: 2.9e6, Ecm: 33e9, barSpacing: 0.3 }).srMax, 1.3 * 0.362, 1e-12);
});

test('chloride ingress, Fick\'s second law: C/Cs = erfc(x/2√(Dt)) at a point; initiation time closed form with and without ageing', () => {
  // x = 50 mm, D = 1e-12 m²/s, t = 50 years: √(Dt) = √(1e-12·1.57788e9) = 39.72 mm; z = 0.6294; erfc(z) = 0.3734
  const t = 50 * YEAR;
  const c = chlorideProfile({ x: 0.05, t, Cs: 3.0, D0: 1e-12 });
  close(c.C / 3.0, 0.3734, 1e-3);
  close(c.C, 3.0 * erfc(0.05 / (2 * Math.sqrt(1e-12 * t))), 1e-12);
  assert.equal(chlorideProfile({ x: 0, t, Cs: 3.0, Ci: 0.1, D0: 1e-12 }).C, 3.0);
  close(chlorideProfile({ x: 5, t, Cs: 3.0, Ci: 0.1, D0: 1e-12 }).C, 0.1, 1e-9);
  // initiation: the profile evaluated at the returned time gives the critical content at the cover depth
  const ti = chlorideInitiationTime({ cover: 0.05, Cs: 3.0, Ccrit: 0.6, D0: 1e-12 });
  close(chlorideProfile({ x: 0.05, t: ti, Cs: 3.0, D0: 1e-12 }).C, 0.6, 1e-7);
  // erfc(z) = 0.2 → z = 0.9062; t = (0.05/1.8124)²/1e-12 = 7.611e8 s = 24.1 years
  close(ti / YEAR, 24.1, 5e-3);
  // doubling the cover quadruples the time (constant D)
  close(chlorideInitiationTime({ cover: 0.1, Cs: 3.0, Ccrit: 0.6, D0: 1e-12 }) / ti, 4, 1e-9);
  // ageing: Da(t) = D0 (t0/t)^m; consistency with the profile, and longer life than without ageing for t > t0
  const o = { cover: 0.05, Cs: 3.0, Ccrit: 0.6, Ci: 0.05, D0: 8e-12, t0: 28 * 86400, m: 0.4, dx: 0.005 };
  const ta = chlorideInitiationTime(o);
  close(chlorideProfile({ x: 0.05, t: ta, Cs: 3.0, Ci: 0.05, D0: 8e-12, m: 0.4, dx: 0.005 }).C, 0.6, 1e-7);
  assert.ok(ta > chlorideInitiationTime({ ...o, m: 0 }));
  close(chlorideProfile({ x: 0.05, t: ta, Cs: 3.0, D0: 8e-12, m: 0.4 }).Da, 8e-12 * (28 * 86400 / ta) ** 0.4, 1e-12);
  assert.equal(chlorideInitiationTime({ cover: 0.05, Cs: 0.5, Ccrit: 0.6, D0: 1e-12 }), Infinity);
});

test('carbonation √t law and Faraday corrosion rate: 1 µA/cm² = 11.6 µm/year; bar section loss', () => {
  // K = 4 mm/√year, cover 30 mm → 56.25 years; depth after 25 years = 20 mm
  const K = 0.004 / Math.sqrt(YEAR);
  close(K, 4 * 1.780e-7, 1e-3);
  const c = carbonation({ K, t: 25 * YEAR, cover: 0.03 });
  close(c.depth, 0.02, 1e-12); close(c.initiationTime / YEAR, 56.25, 1e-12);
  const rate = corrosionRate({ icorr: 0.01 });
  close(rate * YEAR, 11.6e-6, 5e-3);
  // 20 mm bar, 1 µA/cm², 30 years after initiation at year 20: p = 0.349 mm, d = 19.30 mm, area ratio 0.9314
  const b = barSectionLoss({ diameter: 0.02, rate, t: 50 * YEAR, initiationTime: 20 * YEAR, pittingFactor: 6 });
  close(b.penetration, rate * 30 * YEAR, 1e-12); close(b.areaRatio, ((0.02 - 2 * b.penetration) / 0.02) ** 2, 1e-12);
  close(b.areaRatio, 0.9314, 2e-3); close(b.pitDepth, 6 * b.penetration, 1e-12);
  assert.equal(barSectionLoss({ diameter: 0.02, rate, t: 10 * YEAR, initiationTime: 20 * YEAR }).lossFraction, 0);
});

test('residual capacity with corrosion: under-reinforced beam loses capacity nearly in proportion to the bar area', () => {
  const s = rect(1472.6e-6);
  const rate = corrosionRate({ icorr: 0.02 });
  const h = residualCapacityHistory({ section: s, material: mat, times: [0, 20 * YEAR, 50 * YEAR, 100 * YEAR], initiationTime: 20 * YEAR, rate });
  assert.equal(h[0].capacityRatio, 1); close(h[1].capacityRatio, 1, 1e-12);
  assert.ok(h[2].capacityRatio < 1 && h[3].capacityRatio < h[2].capacityRatio);
  for (const p of h.slice(2)) {
    // M ≈ As fyd (d − 0.416x): the lever arm grows slightly as As falls, so the capacity ratio is a little above the area ratio
    assert.ok(p.capacityRatio > p.areaRatio && p.capacityRatio < p.areaRatio * 1.06, `${p.capacityRatio} vs ${p.areaRatio}`);
    close(p.areaRatio, barSectionLoss({ diameter: 0.025, rate, t: p.t, initiationTime: 20 * YEAR }).areaRatio, 1e-12);
  }
});

test('steel pile corrosion by zone: loss = rate × (life − protection life); remaining tubular section', () => {
  const mm = 1e-3 / YEAR;
  const z = pileCorrosionLoss({
    zones: [{ name: 'splash', rate: 0.3 * mm }, { name: 'immersed', rate: 0.1 * mm, protectionLife: 20 * YEAR }, { name: 'both faces', rate: 0.05 * mm, faces: 2 }],
    t: 50 * YEAR, thickness: 0.02, D: 1.0,
  });
  close(z[0].loss, 0.015, 1e-12); close(z[0].remaining, 0.005, 1e-12);
  close(z[1].loss, 0.003, 1e-12); close(z[2].loss, 0.005, 1e-12);
  // splash zone: outer diameter 0.97 m, wall 5 mm: A = π/4 (0.97² − 0.96²) ; original π/4 (1 − 0.96²)
  close(z[0].areaRatio, (0.97 ** 2 - 0.96 ** 2) / (1 - 0.96 ** 2), 1e-12);
  close(z[0].modulusRatio, ((0.97 ** 4 - 0.96 ** 4) / 0.97) / (1 - 0.96 ** 4), 1e-12);
  // two faces: 2.5 mm from each side
  close(z[2].areaRatio, (0.995 ** 2 - 0.965 ** 2) / (1 - 0.96 ** 2), 1e-12);
  assert.equal(pileCorrosionLoss({ zones: [{ name: 'x', rate: 1 * mm }], t: 100 * YEAR, thickness: 0.02 })[0].remaining, 0);
});
