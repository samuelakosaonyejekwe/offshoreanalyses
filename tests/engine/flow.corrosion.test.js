import test from 'node:test';
import assert from 'node:assert/strict';
import * as C from '../../src/engine/flow/corrosion.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const T60 = 333.15; const bar = 1e5; const F = 96485.33;

test('de Waard–Milliams nomogram equation: 60 °C, 1 bar CO2 → 10^(5.8 − 1710/333.15) = 4.65 mm/yr; rate ∝ p^0.67', () => {
  rel(C.deWaardMilliams({ T: T60, pCO2: bar }), 10 ** (5.8 - 1710 / T60), 1e-12);
  rel(C.deWaardMilliams({ T: T60, pCO2: bar }), 4.647, 1e-3);
  rel(C.deWaardMilliams({ T: T60, pCO2: 2 * bar }) / C.deWaardMilliams({ T: T60, pCO2: bar }), 2 ** 0.67, 1e-12);
  // 20 °C, 0.5 bar: log V = 5.8 − 5.8332 + 0.67 × (−0.30103) = −0.2349 → 0.582 mm/yr
  rel(C.deWaardMilliams({ T: 293.15, pCO2: 0.5 * bar }), 0.582, 2e-3);
  assert.throws(() => C.deWaardMilliams({ T: T60, pCO2: -1 }), RangeError);
});

test('CO2 fugacity coefficient and pH of CO2-saturated water', () => {
  // log φ = (0.0031 − 1.4/333.15) × 100 = −0.11023 → φ = 0.7759
  rel(C.co2FugacityCoefficient({ T: T60, P: 100 * bar }), 10 ** ((0.0031 - 1.4 / T60) * 100), 1e-12);
  rel(C.co2FugacityCoefficient({ T: T60, P: 100 * bar }), 0.7759, 1e-3);
  close(C.co2FugacityCoefficient({ T: T60, P: 0 }), 1, 1e-12);
  // 60 °C, 1 bar: 1995 → 3.82 + 0.2304 = 4.0504; 1991 → 3.71 + 0.2502 = 3.9602; ten times the pressure lowers pH by 0.5
  close(C.phCO2Saturated({ T: T60, pCO2: bar }), 4.0504, 1e-9); close(C.phCO2Saturated({ T: T60, pCO2: bar, version: '1991' }), 3.9602, 1e-9);
  close(C.phCO2Saturated({ T: T60, pCO2: 10 * bar }) - C.phCO2Saturated({ T: T60, pCO2: bar }), -0.5, 1e-12);
});

test('de Waard–Lotz–Milliams correction factors by hand: scale, pH and glycol', () => {
  // 60 °C, 1 bar, ideal gas: below the scaling temperature 2400/6.7 = 358.2 K → no scale factor
  const a = C.deWaardLotzMilliams({ T: T60, pCO2: bar });
  rel(a.nomogram, 4.647, 1e-3); close(a.Fscale, 1, 1e-12); close(a.FpH, 1, 1e-12); close(a.Fglycol, 1, 1e-12);
  rel(a.rate, a.nomogram, 1e-12); rel(a.Tscale, 2400 / 6.7, 1e-12); rel(a.fCO2, bar, 1e-12);
  // 100 °C, 1 bar: log F_scale = 2400/373.15 − 6.7 = −0.2683 → 0.539
  const s = C.deWaardLotzMilliams({ T: 373.15, pCO2: bar });
  rel(s.Fscale, 10 ** (2400 / 373.15 - 6.7), 1e-12); rel(s.Fscale, 0.539, 2e-3); rel(s.rate, s.nomogram * s.Fscale, 1e-12);
  // pH 6 at 60 °C, 1 bar: pH_sat = min(1.36 + 1307/333.15, 5.4) = 5.2832; log F = −0.13 × 0.7168^1.6 = −0.0763 → 0.839
  const p = C.deWaardLotzMilliams({ T: T60, pCO2: bar, pH: 6 });
  close(p.pHsat, 1.36 + 1307 / T60, 1e-9); rel(p.FpH, 10 ** (-0.13 * (6 - p.pHsat) ** 1.6), 1e-12); rel(p.FpH, 0.839, 2e-3);
  // pH 4 (below saturation): log F = 0.32 × 1.2832 → 2.574
  rel(C.deWaardLotzMilliams({ T: T60, pCO2: bar, pH: 4 }).FpH, 10 ** (0.32 * (p.pHsat - 4)), 1e-12);
  // 50 wt % water in glycol: log F = 1.6 (log 50 − 2) = −0.4816 → 0.330
  const g = C.deWaardLotzMilliams({ T: T60, pCO2: bar, waterInGlycol: 50 });
  rel(g.Fglycol, 0.330, 2e-3); rel(g.rate, a.rate * g.Fglycol, 1e-12);
  // total pressure lowers the fugacity; extraFactor multiplies
  const f = C.deWaardLotzMilliams({ T: T60, pCO2: bar, P: 100 * bar, extraFactor: 0.5 });
  rel(f.fCO2, 0.7759 * bar, 1e-3); rel(f.rate, 0.5 * a.rate * 0.7759 ** 0.67, 1e-3);
});

test('de Waard–Lotz–Dugstad resistance model: reaction and mass-transfer rates in series', () => {
  // 60 °C, 1 bar, U = 2 m/s, d = 0.2 m: log Vr = 4.93 − 1119/333.15 = 1.5711 → 37.25; Vm = 2.45 × 2^0.8/0.2^0.2 = 5.886
  const r = C.deWaardLotzDugstad({ T: T60, pCO2: bar, U: 2, d: 0.2 });
  rel(r.Vr, 10 ** (4.93 - 1119 / T60), 1e-12); rel(r.Vm, 2.45 * 2 ** 0.8 / 0.2 ** 0.2, 1e-12);
  rel(r.rate, 1 / (1 / r.Vr + 1 / r.Vm), 1e-12); rel(r.rate, 5.083, 1e-3); close(r.pHCO2, 4.0504, 1e-9);
  // limits: fast flow → reaction controlled; slow flow → mass-transfer controlled
  const fast = C.deWaardLotzDugstad({ T: T60, pCO2: bar, U: 100, d: 0.2 });
  // U = 100 m/s: Vm = 2.45 × 100^0.8/0.2^0.2 = 134.6 → 1/(1/37.25 + 1/134.6) = 29.18, approaching Vr from below
  rel(fast.rate, 29.18, 1e-3); assert.ok(fast.rate < r.Vr && fast.rate > r.rate);
  assert.throws(() => C.deWaardLotzDugstad({ T: T60, pCO2: bar, U: 1e6, d: 0.2 }), RangeError);
  rel(C.deWaardLotzDugstad({ T: T60, pCO2: bar, U: 1e-4, d: 0.2 }).rate, 2.45 * 1e-4 ** 0.8 / 0.2 ** 0.2, 1e-3);
  // one pH unit above saturation lowers the reaction rate by 10^0.34
  rel(C.deWaardLotzDugstad({ T: T60, pCO2: bar, U: 2, d: 0.2, pH: 5.0504 }).Vr, r.Vr / 10 ** 0.34, 1e-9);
});

test('inhibitor availability, corrosion allowance, sour-service and top-of-line screening', () => {
  close(C.inhibitedRate({ uninhibited: 5, inhibited: 0.1, availability: 0.9 }), 0.59, 1e-12);
  close(C.inhibitedRate({ uninhibited: 5, availability: 1 }), 0.1, 1e-12); close(C.inhibitedRate({ uninhibited: 5, availability: 0 }), 5, 1e-12);
  const a = C.corrosionAllowance({ rate: 0.3, life: 25, available: 6 });
  close(a.required, 7.5, 1e-12); assert.equal(a.feasible, true); close(a.lifeAchieved, 20, 1e-12);
  assert.equal(C.corrosionAllowance({ rate: 0.5, life: 25 }).feasible, false); // 12.5 mm > 8 mm
  assert.equal(C.corrosionAllowance({ rate: 0.5, life: 25 }).lifeAchieved, null);
  // ISO 15156 limit 0.3 kPa; regime by the CO2/H2S ratio (500 and 20)
  const s = C.sourServiceScreening({ pH2S: 300, pCO2: 2e5 });
  assert.equal(s.sourService, true); close(s.ratio, 2e5 / 300, 1e-9); assert.equal(s.regime, 'sweet');
  assert.equal(C.sourServiceScreening({ pH2S: 299, pCO2: 2e5 }).sourService, false);
  assert.equal(C.sourServiceScreening({ pH2S: 2000, pCO2: 2e5 }).regime, 'mixed');
  assert.equal(C.sourServiceScreening({ pH2S: 2e4, pCO2: 2e5 }).regime, 'sour');
  assert.equal(C.sourServiceScreening({ pH2S: 0, pCO2: 0 }).regime, 'none');
  const t = C.topOfLineScreening({ condensationRate: 5e-4 });
  assert.equal(t.atRisk, true); close(t.ratio, 2, 1e-12); assert.equal(C.topOfLineScreening({ condensationRate: 1e-4 }).atRisk, false);
});

test('mass transfer (Berger & Hau), oxygen corrosion and Faraday conversion', () => {
  // Re = 1e5, Sc = 500: Sh = 0.0165 × 1e5^0.86 × 500^0.33 = 0.0165 × 19953 × 7.775 = 2560
  const m = C.massTransferBergerHau({ Re: 1e5, Sc: 500, Ddiff: 2e-9, d: 0.1 });
  rel(m.Sh, 0.0165 * 1e5 ** 0.86 * 500 ** 0.33, 1e-12); rel(m.Sh, 2560, 5e-3); rel(m.km, m.Sh * 2e-9 / 0.1, 1e-12);
  // 1 A/m² on iron: 55.845e-3/(2 × 96485 × 7870) = 3.677e-11 m/s = 1.16 mm/yr
  rel(C.currentToRate({ i: 1 }), 1.16, 2e-3);
  rel(C.currentToRate({ i: 1 }), 0.055845 / (2 * F * 7870) * 1000 * 365.25 * 86400, 2e-3);
  // zinc (65.38 g/mol, n = 2, 7140 kg/m³): 1.50 mm/yr per A/m²
  rel(C.currentToRate({ i: 1, M: 0.06538, n: 2, rho: 7140 }), 1.50, 5e-3);
  // 8 ppm O2 (8e-3 kg/m³ = 0.25 mol/m³), km = 5e-5 m/s: i_lim = 4 F × 5e-5 × 0.25 = 4.824 A/m²
  const o = C.oxygenCorrosionRate({ cO2: 8e-3, km: 5e-5 });
  rel(o.iLim, 4 * F * 5e-5 * 8e-3 / 0.0319988, 1e-4); rel(o.iLim, 4.824, 1e-3);
  rel(o.rate, C.currentToRate({ i: o.iLim }), 1e-9);
  // equivalent mole balance: 2 mol Fe per mol O2
  rel(o.rate, 2 * 5e-5 * 0.25 * 0.055845 / 7870 * 1000 * 365.25 * 86400, 3e-3);
});

test('electrode kinetics: Butler–Volmer limits, Tafel slope, mixed potential closed form and diffusion limit', () => {
  // symmetric: i = 2 i0 sinh(F η/(2RT)); small η: i ≈ i0 F η/(RT)
  const fRT = F / (8.314462 * 298.15);
  rel(C.butlerVolmer({ i0: 1e-3, eta: 0.05 }), 2e-3 * Math.sinh(0.5 * fRT * 0.05), 1e-6);
  rel(C.butlerVolmer({ i0: 1e-3, eta: 1e-5 }), 1e-3 * fRT * 1e-5, 1e-6);
  close(C.butlerVolmer({ i0: 1e-3, eta: 0 }), 0, 1e-18); close(C.butlerVolmer({ i0: 1e-3, eta: -0.05 }), -C.butlerVolmer({ i0: 1e-3, eta: 0.05 }), 1e-15);
  // b = 2.303 RT/(αnF) = 0.1183 V/decade at 25 °C; η for a 100-fold current = 2b
  const t = C.tafel({ i: 1e-1, i0: 1e-3 });
  rel(t.b, 0.1183, 1e-3); rel(t.eta, 2 * t.b, 1e-12); assert.equal(C.tafel({}).eta, null);
  rel(C.tafel({ alpha: 0.5, n: 2 }).b, 0.1183 / 2, 1e-3);
  // iron in acid: anodic E0 = −0.44 V, i0 = 1e-4, b = 0.06; cathodic E0 = 0 V, i0 = 1e-3, b = 0.12
  //   E_corr = [−0.44/0.06 + 0 + log10(10)]/(1/0.06 + 1/0.12) = (−7.3333 + 1)/25 = −0.25333 V
  //   i_corr = 1e-4 × 10^((−0.25333 + 0.44)/0.06) = 1e-4 × 10^3.1111 = 0.1292 A/m²
  const an = { E0: -0.44, i0: 1e-4, b: 0.06 }; const ca = { E0: 0, i0: 1e-3, b: 0.12 };
  const m = C.mixedPotential({ anodic: an, cathodic: ca });
  close(m.Ecorr, (-0.44 / 0.06 + 1) / 25, 1e-9); rel(m.icorr, 1e-4 * 10 ** ((m.Ecorr + 0.44) / 0.06), 1e-9); rel(m.icorr, 0.1292, 1e-3);
  rel(m.rate, C.currentToRate({ i: m.icorr }), 1e-12);
  // a diffusion limit far below the activation current controls the corrosion current
  const lim = C.mixedPotential({ anodic: an, cathodic: { ...ca, iLim: 1e-3 } });
  rel(lim.icorr, 1e-3, 5e-3); assert.ok(lim.Ecorr < m.Ecorr);
  // polarisation curve: net current is zero at E_corr and the partial currents cross there
  const pc = C.polarisationCurve({ anodic: an, cathodic: ca, n: 61 });
  assert.equal(pc.E.length, 61); close(pc.Ecorr, m.Ecorr, 1e-12);
  close(pc.E[30], m.Ecorr, 1e-12); rel(pc.iAnodic[30], m.icorr, 1e-9); close(pc.iNet[30], 0, 1e-9);
  assert.ok(pc.iNet[0] < 0 && pc.iNet[60] > 0);
  rel(pc.iAnodic[60] / pc.iAnodic[30], 10 ** (0.3 / 0.06), 1e-9);
});
