import test from 'node:test';
import assert from 'node:assert/strict';
import * as W from '../../src/engine/flow/wax.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const R = 8.314462618;

test('solubility curves: table interpolation, exponential fit and zero slope above the WAT', () => {
  const t = W.makeSolubilityCurve({ table: [{ T: 280, C: 0.01 }, { T: 300, C: 0.05 }], WAT: 300 });
  close(t.C(290), 0.03, 1e-12); close(t.dCdT(290), 0.002, 1e-12); assert.equal(t.WAT, 300);
  close(t.C(270), 0.01, 1e-12); close(t.dCdT(270), 0, 1e-15); close(t.C(320), 0.05, 1e-12); close(t.dCdT(310), 0, 1e-15);
  const e = W.makeSolubilityCurve({ exponential: { a: 0.02, b: 0.05, Tref: 290 }, WAT: 305 });
  close(e.C(290), 0.02, 1e-12); rel(e.C(300), 0.02 * Math.exp(0.5), 1e-12); rel(e.dCdT(300), 0.05 * 0.02 * Math.exp(0.5), 1e-9);
  close(e.dCdT(306), 0, 1e-15);
  assert.equal(W.makeSolubilityCurve({ exponential: { a: 0.02, b: 0.05, Tref: 290 } }).WAT, null);
  assert.throws(() => W.makeSolubilityCurve({}), RangeError);
});

test('Won fusion properties of n-C20 (282.55 g/mol): Tf = 310.5 K, ΔHf = 12 511 cal/mol = 52.3 kJ/mol', () => {
  const f = W.wonFusionProperties(0.28255);
  close(f.Tf, 374.5 + 0.02617 * 282.55 - 20172 / 282.55, 1e-9); close(f.Tf, 310.50, 0.01);
  rel(f.dHf, 0.1426 * 282.55 * f.Tf * 4.184, 1e-9); rel(f.dHf, 52345, 1e-3);
  assert.ok(W.wonFusionProperties(0.45).Tf > f.Tf);
});

test('ideal-solution wax equilibrium: single wax in a non-crystallising solvent has a closed-form solubility and WAT', () => {
  // wax: Tf = 320 K, ΔHf = 60 kJ/mol; solvent given a very low melting point so that it never solidifies
  const wax = { M: 0.3, z: 0.1, Tf: 320, dHf: 60000 }; const solv = { M: 0.1, z: 0.9, Tf: 50, dHf: 10000 };
  const T = 300;
  const xsat = Math.exp(-60000 / (R * T) * (1 - T / 320)); // 0.2227: 10 % wax stays dissolved
  const none = W.waxPrecipitationIdeal({ T, components: [wax, solv] });
  close(none.solidFraction, 0, 1e-12); rel(none.K[0], 1 / xsat, 1e-9); rel(none.x[0], 0.1, 1e-12);
  // at 290 K: x_sat = exp(−60000/(R × 290) × (1 − 290/320)) = 0.0972 < 0.1 → S = (z − x_sat)/(1 − x_sat)
  const xs = Math.exp(-60000 / (R * 290) * (1 - 290 / 320));
  const p = W.waxPrecipitationIdeal({ T: 290, components: [wax, solv] });
  rel(p.x[0], xs, 1e-6); rel(p.solidFraction, (0.1 - xs) / (1 - xs), 1e-5); close(p.s[0], 1, 1e-6);
  // colder → more solid
  assert.ok(W.waxPrecipitationIdeal({ T: 270, components: [wax, solv] }).solidFraction > p.solidFraction);
  // WAT: z K(T) = 1 → 1/T = 1/Tf − R ln z/ΔHf → 290.49 K
  const wat = W.waxAppearanceIdeal({ components: [wax, solv] });
  close(wat, 1 / (1 / 320 - R * Math.log(0.1) / 60000), 1e-3);
  assert.ok(wat > 290 && wat < 300);
  // more wax in the feed raises the WAT
  assert.ok(W.waxAppearanceIdeal({ components: [{ ...wax, z: 0.3 }, { ...solv, z: 0.7 }] }) > wat);
});

test('diffusivity correlations by hand: Hayduk–Minhas 2.29e-10 m²/s and Wilke–Chang 1.65e-10 m²/s at 300 K, 5 cP', () => {
  // 13.3e-8 × 300^1.47 × 5^(10.2/430 − 0.791)/430^0.71 = 13.3e-8 × 4379 × 0.2909/74.09 = 2.286e-6 cm²/s
  const hm = W.diffusivityHaydukMinhas({ T: 300, mu: 0.005 });
  rel(hm, 13.3e-8 * 300 ** 1.47 * 5 ** (10.2 / 430 - 0.791) / 430 ** 0.71 * 1e-4, 1e-9); rel(hm, 2.286e-10, 2e-3);
  // 7.4e-8 × √200 × 300/(5 × 430^0.6) = 1.651e-6 cm²/s
  const wc = W.diffusivityWilkeChang({ T: 300, mu: 0.005, MB: 0.2 });
  rel(wc, 7.4e-8 * Math.sqrt(200) * 300 / (5 * 430 ** 0.6) * 1e-4, 1e-9); rel(wc, 1.651e-10, 2e-3);
  // both fall with viscosity and rise with temperature
  assert.ok(W.diffusivityHaydukMinhas({ T: 300, mu: 0.02 }) < hm && W.diffusivityWilkeChang({ T: 330, mu: 0.005, MB: 0.2 }) > wc);
});

test('deposition flux and the Matzain growth rate by hand', () => {
  // ρ D (dC/dT)(dT/dr) = 850 × 2e-10 × 0.002 × 5000 = 1.7e-6 kg/(m²·s); shear dispersion k* C_w γ̇ = 1e-9 × 0.05 × 400
  const f = W.depositionFlux({ rhoOil: 850, Dwo: 2e-10, dCdT: 0.002, dTdr: 5000, kStar: 1e-9, Cw: 0.05, shearRate: 400 });
  rel(f.molecular, 1.7e-6, 1e-12); rel(f.shearDispersion, 2e-8, 1e-12); rel(f.flux, 1.72e-6, 1e-12);
  // Π1 = 15/(1 − 0.4) = 25; N_SR = 850 × 1 × 1e-3/0.005 = 170; Π2 = 0.055 × 170^1.4 = 72.95
  const m = W.matzainGrowthRate({ Dwo: 2e-10, dCdT: 0.002, dTdr: 5000, delta: 1e-3, rhoOil: 850, vOil: 1, muOil: 0.005, oilFraction: 40 });
  close(m.Pi1, 25, 1e-12); close(m.NSR, 170, 1e-9); rel(m.Pi2, 0.055 * 170 ** 1.4, 1e-12); rel(m.Pi2, 72.95, 1e-3);
  rel(m.rate, 25 / (1 + m.Pi2) * 2e-10 * 0.002 * 5000, 1e-12);
  // clean pipe (δ = 0): no shear stripping; annular regime uses √(ρm ρo)
  close(W.matzainGrowthRate({ Dwo: 2e-10, dCdT: 0.002, dTdr: 5000, delta: 0, rhoOil: 850, vOil: 1, muOil: 0.005 }).Pi2, 0, 1e-15);
  rel(W.matzainGrowthRate({ Dwo: 2e-10, dCdT: 0.002, dTdr: 5000, delta: 1e-3, rhoOil: 850, vOil: 1, muOil: 0.005, regime: 'annular', rhoMix: 212.5 }).NSR, 85, 1e-9);
});

test('wax deposition along a line: no deposit above the WAT; growth, insulation feedback and pressure-drop rise below it', () => {
  const sol = W.makeSolubilityCurve({ table: [{ T: 270, C: 0.01 }, { T: 310, C: 0.09 }], WAT: 310 });
  const base = { D: 0.2, L: 20000, n: 40, mDot: 30, Tin: 330, Tamb: 277, Uwall: 10, rhoOil: 850, muOil: 0.005, cpOil: 2000, kOil: 0.13,
    solubility: sol, tEnd: 30 * 86400, dt: 86400 };
  const r = W.waxDepositionLine(base);
  assert.equal(r.x.length, 40); assert.equal(r.thickness.length, 40);
  // clean-pipe arrival temperature: exponential decay with the overall coefficient (inner film in series with Uwall),
  // so it lies between the wall-only estimate and the inlet temperature
  const wallOnly = 277 + 53 * Math.exp(-10 * Math.PI * 0.2 * 20000 / (30 * 2000));
  assert.ok(r.Tarrival[0] > wallOnly && r.Tarrival[0] < 330, `arrival ${r.Tarrival[0]} vs wall-only ${wallOnly}`);
  assert.ok(r.Tarrival[0] - wallOnly < 8);
  // upstream, where the deposit surface is above the WAT, nothing deposits; further down it does
  close(r.thickness[0], 0, 1e-12); assert.ok(Math.max(...r.thickness) > 1e-5);
  const k = r.t.length - 1;
  for (let i = 1; i <= k; i++) { assert.ok(r.thicknessMax[i] >= r.thicknessMax[i - 1] - 1e-15); assert.ok(r.waxVolume[i] >= r.waxVolume[i - 1] - 1e-15); }
  assert.ok(r.dP[k] > r.dP[0], 'the deposit raises the frictional pressure drop');
  assert.ok(r.Tarrival[k] >= r.Tarrival[0] - 1e-9, 'the deposit insulates the line');
  // deposit volume = Σ π (D² − d²)/4 Δx
  let vol = 0; for (const d of r.thickness) vol += Math.PI * (0.2 ** 2 - (0.2 - 2 * d) ** 2) / 4 * 500;
  rel(r.waxVolume[k], vol, 1e-6);
  for (let i = 0; i < 40; i++) assert.ok(r.Tsurface[i] <= r.Tbulk[i] + 1e-9 && r.Tsurface[i] >= 277 - 1e-9);
  // warm surroundings (above the WAT): no deposition at all
  const warm = W.waxDepositionLine({ ...base, Tamb: 315 });
  close(Math.max(...warm.thickness), 0, 1e-15); close(warm.dP[warm.dP.length - 1], warm.dP[0], 1e-9 * warm.dP[0]);
  // a waxier deposit (larger F_w) grows more slowly; the Matzain model runs and deposits too
  assert.ok(Math.max(...W.waxDepositionLine({ ...base, waxFraction: 0.6 }).thickness) < Math.max(...r.thickness));
  assert.ok(Math.max(...W.waxDepositionLine({ ...base, model: 'matzain' }).thickness) > 0);
  // stopping criterion
  const lim = 0.5 * Math.max(...r.thickness);
  const st = W.waxDepositionLine({ ...base, maxThickness: lim });
  assert.ok(st.tStop < base.tEnd && Math.max(...st.thickness) >= lim * (1 - 1e-9));
});

test('pigging interval by interpolation and extrapolation; gel restart pressure 4 τ L/D; gelling screen', () => {
  const h = { t: [0, 10, 20], thicknessMax: [0, 1e-3, 2e-3], dP: [1e5, 1.1e5, 1.3e5] };
  const a = W.piggingInterval({ ...h, maxThickness: 1.5e-3 });
  close(a.interval, 15, 1e-9); assert.equal(a.extrapolated, false); assert.match(a.governing, /thick/i);
  // pressure criterion 1.2 × clean: between 10 and 20 → 15; thickness criterion 3 mm extrapolated to 30
  const b = W.piggingInterval({ ...h, maxThickness: 3e-3, dPratio: 1.2 });
  close(b.byPressure, 15, 1e-9); close(b.byThickness, 30, 1e-9); close(b.interval, 15, 1e-9); assert.match(b.governing, /press/i);
  const c = W.piggingInterval({ ...h, maxThickness: 3e-3 });
  close(c.interval, 30, 1e-9); assert.equal(c.extrapolated, true);
  assert.equal(W.piggingInterval({ t: [0, 10], thicknessMax: [0, 0], maxThickness: 1e-3 }).interval, Infinity);
  // τ_y = 50 Pa over 1 km of 0.2 m pipe: 4 × 50 × 1000/0.2 = 1 MPa
  close(W.gelRestartPressure({ tauY: 50, L: 1000, D: 0.2 }), 1e6, 1e-6);
  const g = W.gellingScreen({ Tmin: 280, pourPoint: 278 });
  close(g.margin, 2, 1e-12); assert.equal(g.gelRisk, true);
  assert.equal(W.gellingScreen({ Tmin: 290, pourPoint: 278 }).gelRisk, false);
  assert.equal(W.gellingScreen({ Tmin: 280, pourPoint: 278, margin: 1 }).gelRisk, false);
});
