import test from 'node:test';
import assert from 'node:assert/strict';
import * as H from '../../src/engine/flow/hydrate.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const psi = 6894.757; const FtoK = (f) => (f - 32) / 1.8 + 273.15;

test('gas-gravity hydrate temperature at 1000 psia: Hammerschmidt 63.7 °F, Towler 64.6 °F and Motiee 63.3 °F (γ = 0.7)', () => {
  const P = 1000 * psi;
  // 8.9 × 1000^0.285 = 8.9 × 7.1614 = 63.74 °F
  close(H.hydrateTemperature({ P, method: 'hammerschmidt' }), FtoK(8.9 * 1000 ** 0.285), 1e-6);
  close(H.hydrateTemperature({ P, method: 'hammerschmidt' }), FtoK(63.74), 0.01);
  // 13.47 × 6.90776 + 34.27 ln 0.7 − 1.675 × 6.90776 ln 0.7 − 20.35 = 93.047 − 12.223 + 4.127 − 20.35 = 64.60 °F
  close(H.hydrateTemperature({ P, gammaG: 0.7, method: 'towler' }), FtoK(64.60), 0.01);
  // −238.2447 + 236.9900 − 48.1729 + 244.6317 − 73.9188 − 57.9685 = 63.32 °F
  close(H.hydrateTemperature({ P, gammaG: 0.7 }), FtoK(63.317), 0.01);
  // heavier gas and higher pressure both raise the hydrate temperature
  assert.ok(H.hydrateTemperature({ P, gammaG: 0.8 }) > H.hydrateTemperature({ P, gammaG: 0.6 }));
  assert.ok(H.hydrateTemperature({ P: 2 * P, gammaG: 0.7 }) > H.hydrateTemperature({ P, gammaG: 0.7 }));
  assert.throws(() => H.hydrateTemperature({ P, method: 'katz' }), RangeError);
});

test('hydrate pressure inverts the temperature correlation; inhibitor shifts the curve', () => {
  for (const method of ['motiee', 'towler', 'hammerschmidt']) {
    const T = H.hydrateTemperature({ P: 8e6, gammaG: 0.7, method });
    rel(H.hydratePressure({ T, gammaG: 0.7, method }), 8e6, 1e-6, method);
  }
  // with 5 K of depression the same pressure is in equilibrium 5 K colder
  const T = H.hydrateTemperature({ P: 8e6, gammaG: 0.7 });
  rel(H.hydratePressure({ T: T - 5, gammaG: 0.7, inhibitorDT: 5 }), 8e6, 1e-6);
  assert.equal(H.hydratePressure({ T: 400, gammaG: 0.7 }), Infinity);
});

test('Hammerschmidt and Nielsen–Bucklin inhibitor relations by hand, with their inverses', () => {
  // 20 wt % methanol: ΔT = 2335 × 20/(32.042 × 80) = 18.22 °F = 10.12 K
  const d = H.hammerschmidtDepression({ w: 20 });
  rel(d, 2335 * 20 / (32.042 * 80) / 1.8, 1e-9); rel(d, 10.12, 1e-3);
  rel(H.hammerschmidtConcentration({ dT: d }), 20, 1e-9);
  // 40 wt % MEG (K_H = 2700, M = 62.068): 2700 × 40/(62.068 × 60)/1.8 = 16.11 K
  rel(H.hammerschmidtDepression({ w: 40, inhibitor: 'MEG' }), 16.11, 1e-3);
  rel(H.hammerschmidtDepression({ w: 40, M: 0.062068, KH: 4000 }), 16.11 * 4000 / 2700, 1e-3);
  assert.equal(H.INHIBITORS.methanol.KH, 2335); close(H.hammerschmidtDepression({ w: 0 }), 0, 1e-15);
  // x = 0.2: −72 ln 0.8 = 16.07 K
  rel(H.nielsenBucklinDepression(0.2), -72 * Math.log(0.8), 1e-12); rel(H.nielsenBucklinDepression(0.2), 16.07, 1e-3);
  rel(H.nielsenBucklinMoleFraction(16.066), 0.2, 1e-4);
  // mole ↔ mass: x = 0.2 methanol → 6.408/(6.408 + 14.412) = 0.3078
  const w = H.moleToMassFraction(0.2, 0.032042);
  rel(w, 0.3078, 1e-3); rel(H.massToMoleFraction(w, 0.032042), 0.2, 1e-12);
  // 10 wt % NaCl: 0.3422 mol ions against 4.996 mol water per 100 g → x_w = 0.9359 → ΔT = 4.77 K
  rel(H.saltDepression({ w: 10 }), 4.77, 2e-3); close(H.saltDepression({ w: 0 }), 0, 1e-15);
});

test('inhibitor injection mass balance', () => {
  // 10.12 K with methanol → 20 wt % in the water phase; 1 kg/s water, pure methanol: 0.2/(1 − 0.2) = 0.25 kg/s
  const dT = H.hammerschmidtDepression({ w: 20 });
  const r = H.inhibitorInjection({ dT, waterRate: 1 });
  rel(r.w, 0.2, 1e-9); rel(r.injectionRate, 0.25, 1e-9); rel(r.inhibitorRate, 0.25, 1e-9); rel(r.aqueousInhibitor, 0.25, 1e-9);
  close(r.vapourLoss, 0, 1e-15);
  // 40 wt % MEG required, injected as 90 wt % lean MEG: ṁ = 0.4 × 1/(0.9 − 0.4) = 0.8 kg/s of solution
  const m = H.inhibitorInjection({ dT: H.hammerschmidtDepression({ w: 40, inhibitor: 'MEG' }), inhibitor: 'MEG', waterRate: 1, purity: 0.9 });
  rel(m.injectionRate, 0.8, 1e-9); rel(m.inhibitorRate, 0.72, 1e-9);
  // check: aqueous phase = 1 + 0.08 water + 0.72 MEG = 1.8 kg/s at 40 %
  rel(m.inhibitorRate / (1 + m.injectionRate), 0.4, 1e-9);
  // vapour loss: 10 Sm³/s × 0.01 kg/Sm³ per unit mass fraction × 0.2 = 0.02 kg/s more methanol, kept in balance
  const v = H.inhibitorInjection({ dT, waterRate: 1, gasRate: 10, kVap: 0.01 });
  rel(v.vapourLoss, 0.02, 1e-9); rel(v.injectionRate, (0.2 * 1 + 0.8 * 0.02) / 0.8, 1e-9);
  // Nielsen–Bucklin route gives the mole fraction first
  const nb = H.inhibitorInjection({ dT: 16.066, method: 'nielsenBucklin', waterRate: 1 });
  rel(nb.x, 0.2, 1e-4); rel(nb.w, 0.3078, 1e-3);
  assert.throws(() => H.inhibitorInjection({ dT: 30, inhibitor: 'MEG', waterRate: 1, purity: 0.3 }), RangeError);
});

test('subcooling along a profile with a custom equilibrium curve: onset by interpolation, maximum and margin', () => {
  // T_eq = 290 K everywhere; the fluid cools linearly from 300 K to 280 K over 10 km → enters at 5 km
  const r = H.subcoolingProfile({ s: [0, 5000, 10000], P: [1e7, 9e6, 8e6], T: [300, 290, 280], curve: () => 290 });
  assert.deepEqual(r.subcooling, [-10, 0, 10]); close(r.onset, 5000, 1e-9); close(r.maxSubcooling, 10, 1e-12);
  close(r.atMax, 10000, 1e-12); assert.equal(r.inHydrateRegion, true); close(r.margin, -10, 1e-12);
  // 4 K of inhibition moves the entry to 7 km
  close(H.subcoolingProfile({ s: [0, 5000, 10000], P: [1e7, 9e6, 8e6], T: [300, 290, 280], curve: () => 290, inhibitorDT: 4 }).onset, 7000, 1e-9);
  const out = H.subcoolingProfile({ s: [0, 1], P: [1e7, 1e7], T: [320, 310], curve: () => 290 });
  assert.equal(out.onset, null); assert.equal(out.inHydrateRegion, false); close(out.margin, 20, 1e-12);
  // default curve: the correlation at each pressure
  const c = H.subcoolingProfile({ s: [0, 1], P: [8e6, 8e6], T: [280, 280], gammaG: 0.7 });
  close(c.Teq[0], H.hydrateTemperature({ P: 8e6, gammaG: 0.7 }), 1e-12);
});

test('kinetics: induction time forms, Englezos and Vysniauskas–Bishnoi rates, transport limits', () => {
  // t = K exp(B/(T ΔT²)): K = 10 s, B = 2.77e4 K³, T = 277 K, ΔT = 5 K → exponent 4 → 546 s
  rel(H.inductionTime({ K: 10, B: 277 * 25 * 4, T: 277, dT: 5 }), 10 * Math.exp(4), 1e-12);
  rel(H.inductionTime({ K: 10, B: 277 ** 3 * Math.log(2) ** 2 * 3, T: 277, S: 2 }), 10 * Math.exp(3), 1e-12);
  assert.equal(H.inductionTime({ K: 10, B: 1e4, T: 277, dT: 0 }), Infinity);
  // 1/K* = 1/kr + 1/kd: kr = 2e-9, kd = 2e-9 → K* = 1e-9; area 100 m², Δf = 1 MPa → 0.1 mol/s
  rel(H.growthRateEnglezos({ area: 100, f: 5e6, feq: 4e6, kr: 2e-9, kd: 2e-9 }), 0.1, 1e-12);
  rel(H.growthRateEnglezos({ area: 100, f: 5e6, feq: 4e6, Kstar: 1e-9 }), 0.1, 1e-12);
  assert.equal(H.growthRateEnglezos({ area: 100, f: 3e6, feq: 4e6, Kstar: 1e-9 }), 0);
  // r = A a exp(−Ea/RT) exp(−a/ΔT^b) P^γ
  const vb = H.growthRateVysniauskasBishnoi({ A: 2, as: 3, Ea: -1e4, T: 280, dT: 4, a: 0.5, b: 2, P: 5e6, gamma: 1.5 });
  rel(vb, 2 * 3 * Math.exp(1e4 / (8.314462 * 280)) * Math.exp(-0.5 / 16) * 5e6 ** 1.5, 1e-6);
  assert.equal(H.growthRateVysniauskasBishnoi({ A: 2, as: 3, Ea: 0, T: 280, dT: 0, a: 0.5, b: 2, P: 5e6, gamma: 1 }), 0);
  // mass transfer 1e-4 × 10 × 50 = 0.05 mol/s; heat transfer 500 × 10 × 5/54000 = 0.463 mol/s → mass transfer limits
  const t = H.transportLimitedGrowth({ area: 10, kL: 1e-4, cBulk: 150, cEq: 100, h: 500, dTsub: 5 });
  rel(t.massTransfer, 0.05, 1e-12); rel(t.heatTransfer, 25000 / 54000, 1e-12); rel(t.rate, 0.05, 1e-12);
  assert.match(t.limiting, /mass/i); rel(t.series, 1 / (1 / 0.05 + 54000 / 25000), 1e-12);
  assert.match(H.transportLimitedGrowth({ area: 10, kL: 1e-4, cBulk: 150, cEq: 100, h: 500, dTsub: 5, intrinsic: 0.01 }).limiting, /intrinsic|kinetic/i);
});

test('hydrate volume and slurry: volume ratio, conversion along a pipe, suspension viscosities, plugging index', () => {
  // (1000/910) × (6 × 18.015 + 16.043)/(6 × 18.015) = 1.0989 × 1.1484 = 1.262
  rel(H.hydrateVolumeRatio(), 1.262, 1e-3);
  // constant subcooling 4 K, k = 1e-4 /(s·K), v = 2 m/s: X(x) = 1 − exp(−k ΔT x/v) → at 10 km 1 − e^(−2) = 0.8647
  const x = Array.from({ length: 201 }, (_, i) => 50 * i);
  const f = H.hydrateFractionAlongPipe({ x, subcooling: x.map(() => 4), v: 2, k: 1e-4, waterFraction: 0.3 });
  rel(f.conversion[200], 1 - Math.exp(-2), 2e-3); close(f.onset, 0, 1e-12);
  const X = f.conversion[200]; const rv = H.hydrateVolumeRatio();
  rel(f.phi[200], X * 0.3 * rv / (1 - 0.3 * X + X * 0.3 * rv), 1e-9);
  // no subcooling before 5 km: formation starts there
  const g = H.hydrateFractionAlongPipe({ x, subcooling: x.map((xx) => (xx < 5000 ? -1 : 4)), v: 2, k: 1e-4, waterFraction: 0.3 });
  close(g.conversion[99], 0, 1e-12); assert.ok(g.conversion[100] < 0.01); assert.ok(g.onset >= 4950 && g.onset <= 5000); assert.ok(g.conversion[200] > 0.5 && g.conversion[200] < X);
  // Krieger–Dougherty: φ = 0.32 = φ_max/2 → 0.5^(−1.6) = 3.031; dilute limit 1 + 2.5 φ (Einstein)
  rel(H.slurryViscosityKriegerDougherty({ phi: 0.32 }), 0.5 ** -1.6, 1e-12);
  rel(H.slurryViscosityKriegerDougherty({ phi: 1e-4 }), 1 + 2.5e-4, 1e-6);
  assert.equal(H.slurryViscosityKriegerDougherty({ phi: 0.7 }), Infinity);
  // Camargo–Palermo without adhesion: no aggregation, Mills' relation: φ = 0.2 → 0.8/(1 − 0.35)² = 1.893
  const cp = H.slurryViscosityCamargoPalermo({ phi: 0.2, dp: 40e-6, muOil: 0.01, shearRate: 100, Fa: 0 });
  close(cp.dRatio, 1, 1e-9); rel(cp.phiEff, 0.2, 1e-9); rel(cp.muRel, 0.8 / (1 - 0.35) ** 2, 1e-9);
  // adhesion makes aggregates: larger effective fraction and viscosity, relieved by shear
  const ag = H.slurryViscosityCamargoPalermo({ phi: 0.2, dp: 40e-6, muOil: 0.01, shearRate: 100, Fa: 1e-6 });
  const hi = H.slurryViscosityCamargoPalermo({ phi: 0.2, dp: 40e-6, muOil: 0.01, shearRate: 1000, Fa: 1e-6 });
  assert.ok(ag.dRatio > 1 && ag.muRel > cp.muRel && hi.muRel < ag.muRel && hi.dRatio < ag.dRatio);
  rel(ag.phiEff, 0.2 * ag.dRatio ** 0.5, 1e-9);
  // plugging index = a b (0.5 + 0.5 c): ΔT 7 K (a = 1), φ 0.2 (b = 0.5), v = 2 m/s (c = 0.5) → 0.375 'moderate'
  const p = H.pluggingRiskIndex({ subcooling: 7, phi: 0.2, velocity: 2 });
  close(p.exposure, 1, 1e-12); close(p.loading, 0.5, 1e-12); close(p.transport, 0.5, 1e-12); close(p.index, 0.375, 1e-12); assert.equal(p.class, 'moderate');
  assert.equal(H.pluggingRiskIndex({ subcooling: -2, phi: 0.2, velocity: 2 }).class, 'none');
  assert.equal(H.pluggingRiskIndex({ subcooling: 7, phi: 0.4, velocity: 0.5 }).class, 'high');
});

test('depressurisation: equilibrium temperature at the reduced pressure, ice limit and dissociation pressure', () => {
  const d = H.depressurisation({ P: 2e6, Tamb: 277.15, gammaG: 0.7 });
  close(d.Teq, H.hydrateTemperature({ P: 2e6, gammaG: 0.7 }), 1e-12); close(d.drivingForce, 277.15 - d.Teq, 1e-12);
  assert.equal(d.dissociates, d.drivingForce > 0); assert.equal(d.iceRisk, d.Teq < 273.15);
  close(H.hydrateTemperature({ P: d.PiceLimit, gammaG: 0.7 }), 273.15, 1e-4);
  close(H.hydrateTemperature({ P: d.Pdissociation, gammaG: 0.7 }), 277.15, 1e-4);
  assert.ok(d.Pdissociation > d.PiceLimit);
  // at high pressure the plug is stable at seabed temperature
  assert.equal(H.depressurisation({ P: 1.5e7, Tamb: 277.15, gammaG: 0.7 }).dissociates, false);
});

test('van der Waals–Platteeuw framework: single guest with equal constant Langmuir constants has a closed-form pressure', () => {
  // sI: ν = (1/23, 3/23). With C_small = C_large = C, Δv = 0 and T = T0:
  //   (4/23) ln(1 + C P) = Δμ0/(R T0)  →  P = [exp(23 Δμ0/(4 R T0)) − 1]/C
  // Δμ0 = 1264 J/mol, C = 1e-5 1/Pa: exponent 3.2002 → P = 23.54/1e-5 = 2.354 MPa; occupancy C P/(1 + C P) = 0.9592
  const R = 8.314462618; const C = 1e-5;
  const r = H.vdwpEquilibriumPressure({ T: 273.15, structure: 'sI', guests: [{ y: 1, C: () => [C, C] }], ref: { dmu0: 1264, dh0: -4858, dv: 0 } });
  const e = Math.exp(23 * 1264 / (4 * R * 273.15));
  rel(r.P, (e - 1) / C, 1e-6); rel(r.P, 2.354e6, 2e-3);
  rel(r.occupancy[0][0], (e - 1) / e, 1e-6); rel(r.occupancy[1][0], (e - 1) / e, 1e-6);
  rel(r.dmuRT, 1264 / (R * 273.15), 1e-6);
  // diluting the guest to half its mole fraction doubles the equilibrium pressure (ideal gas, Δv = 0)
  const half = H.vdwpEquilibriumPressure({ T: 273.15, structure: 'sI', guests: [{ y: 0.5, C: () => [C, C] }, { y: 0.5, C: () => [0, 0] }], ref: { dmu0: 1264, dh0: -4858, dv: 0 } });
  rel(half.P, 2 * r.P, 1e-6);
  // lower water activity (inhibitor) raises the pressure: (4/23) ln(1 + C P) = Δμ0/(R T0) − ln a_w
  const inh = H.vdwpEquilibriumPressure({ T: 273.15, structure: 'sI', guests: [{ y: 1, C: () => [C, C] }], ref: { dmu0: 1264, dh0: -4858, dv: 0 }, waterActivity: 0.9 });
  rel(inh.P, (Math.exp(23 / 4 * (1264 / (R * 273.15) - Math.log(0.9))) - 1) / C, 1e-6);
  // Parrish–Prausnitz form C = (A/T) exp(B/T) reproduces the same constant
  const pp = H.vdwpEquilibriumPressure({ T: 273.15, structure: 'sI', guests: [{ y: 1, A: [C * 273.15, C * 273.15], B: [0, 0] }], ref: { dmu0: 1264, dh0: -4858, dv: 0 } });
  rel(pp.P, r.P, 1e-9);
});
