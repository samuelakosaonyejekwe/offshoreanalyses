import test from 'node:test';
import assert from 'node:assert/strict';
import {
  makeMixture, zRoots, fugacityCoefficients, properties, wilsonK, rachfordRice, stabilityTest, flashPT,
  bubblePointPressure, dewPointPressure, bubblePointTemperature, dewPointTemperature, upperSaturationPressure,
  vapourPressure, phaseEnvelope, gasViscosityLGE, viscosityLBC, gasThermalConductivity, oilThermalConductivity,
  surfaceTension, waterContentBukacek, STANDARD, standardMolarVolume, standardRateToMass, massToStandardRate,
  mmscfdToSm3s,
} from '../../src/engine/flow/eos.js';
import { makeComponent, R, COMPONENTS } from '../../src/engine/flow/components.js';
import { zFactorReduced } from '../../src/engine/flow/blackoil.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('PR and SRK molar volumes of n-butane at 350 K, 9.4573 bar (Smith, Van Ness & Abbott, Example 3.9)', () => {
  // Textbook constants Tc = 425.1 K, Pc = 37.96 bar, ω = 0.200. Tabulated results (cm³/mol):
  // SRK V^v = 2520, V^l = 127.8; PR V^v = 2486, V^l = 112.6 (experimental 2482 and 115.0)
  const nC4 = makeComponent({ base: 'nC4', Tc: 425.1, Pc: 37.96e5, omega: 0.200 });
  const vol = (eos) => {
    const mix = makeMixture([nC4], { eos });
    const r = zRoots({ mix, T: 350, P: 9.4573e5, z: [1] });
    const f = R * 350 / 9.4573e5 * 1e6;
    return [r.roots[0] * f, r.roots[r.roots.length - 1] * f];
  };
  const pr = vol('PR');
  const srk = vol('SRK');
  assert.ok(Math.abs(pr[0] - 112.6) < 0.1 && Math.abs(pr[1] - 2486) < 1, `PR ${pr}`);
  assert.ok(Math.abs(srk[0] - 127.8) < 0.1 && Math.abs(srk[1] - 2520) < 1, `SRK ${srk}`);
});

test('PR Z-factor of methane at 300 K, 100 bar agrees with the Standing–Katz chart (DAK fit) within 3 %', () => {
  const mix = makeMixture(['C1']);
  const Z = zRoots({ mix, T: 300, P: 1e7, z: [1] }).stable;
  const chart = zFactorReduced({ Tpr: 300 / 190.56, Ppr: 1e7 / 45.99e5 });
  assert.ok(rel(Z, chart) < 0.03, `PR ${Z} vs chart ${chart}`);
  assert.ok(Z > 0.8 && Z < 0.9);
});

test('PR fugacity coefficient of methane equals ∫(Z − 1)/P dP (thermodynamic definition)', () => {
  const mix = makeMixture(['C1']);
  const T = 300;
  const P = 1e7;
  const n = 2000; // Simpson
  let s = 0;
  for (let i = 0; i <= n; i++) {
    const p = P * i / n;
    const g = i === 0 ? null : (zRoots({ mix, T, P: p, z: [1] }).stable - 1) / p;
    // limit of (Z − 1)/P at P → 0 by extrapolation from the first two nodes
    const val = g ?? (2 * (zRoots({ mix, T, P: P / n, z: [1] }).stable - 1) / (P / n)
      - (zRoots({ mix, T, P: 2 * P / n, z: [1] }).stable - 1) / (2 * P / n));
    s += val * (i === 0 || i === n ? 1 : i % 2 ? 4 : 2);
  }
  const lnPhi = s * (P / n) / 3;
  const f = fugacityCoefficients({ mix, T, P, z: [1] });
  assert.ok(Math.abs(f.lnPhi[0] - lnPhi) < 1e-7, `${f.lnPhi[0]} vs ${lnPhi}`);
  assert.ok(f.phi[0] > 0.8 && f.phi[0] < 0.9);
});

test('pure-component vapour pressures by PR: propane at 300 K ≈ 9.98 bar, methane at 150 K ≈ 10.4 bar (NIST)', () => {
  assert.ok(rel(vapourPressure({ component: 'C3', T: 300 }), 9.98e5) < 0.02);
  assert.ok(rel(vapourPressure({ component: 'C1', T: 150 }), 10.4e5) < 0.02);
  assert.ok(rel(vapourPressure({ component: 'C3', T: 300, eos: 'SRK' }), 9.98e5) < 0.03);
  // at the vapour pressure the liquid and vapour roots have equal fugacity
  const mix = makeMixture(['C3']);
  const P = vapourPressure({ component: 'C3', T: 300 });
  const l = fugacityCoefficients({ mix, T: 300, P, z: [1], phase: 'liquid' });
  const v = fugacityCoefficients({ mix, T: 300, P, z: [1], phase: 'vapour' });
  assert.ok(Math.abs(l.lnPhi[0] - v.lnPhi[0]) < 1e-9);
  assert.throws(() => vapourPressure({ component: 'C3', T: 400 }), RangeError);
});

test('PR acentric-factor consistency: Psat(0.7 Tc) = Pc·10^−(1+ω) within 1 %', () => {
  for (const id of ['C1', 'C3', 'nC6', 'nC10', 'CO2']) {
    const c = COMPONENTS[id];
    assert.ok(rel(vapourPressure({ component: id, T: 0.7 * c.Tc }), c.Pc * 10 ** -(1 + c.omega)) < 0.012, id);
  }
});

test('Rachford–Rice: analytic binary and bounds', () => {
  // z = [0.5, 0.5], K = [2, 0.5]: 0.5/(1 + β) − 0.25/(1 − 0.5β) = 0 → β = 0.5
  assert.ok(Math.abs(rachfordRice([0.5, 0.5], [2, 0.5]) - 0.5) < 1e-13);
  assert.equal(rachfordRice([0.5, 0.5], [0.5, 0.2]), 0);
  assert.equal(rachfordRice([0.5, 0.5], [3, 1.2]), 1);
  // Wilson K — hand value for methane at 200 K, 20 bar: (45.99/20)·exp(5.373·1.0115·(1 − 190.56/200))
  const mix = makeMixture(['C1']);
  assert.ok(rel(wilsonK(mix, 200, 20e5)[0], (45.99 / 20) * Math.exp(5.373 * 1.0115 * (1 - 190.56 / 200))) < 1e-12);
});

const oilMix = makeMixture(['C1', 'C3', 'nC5', 'nC10', 'CO2']);
const zOil = [0.6, 0.15, 0.1, 0.1, 0.05];

test('PT flash: component balance and equal fugacities at convergence', () => {
  for (const [T, P] of [[320, 5e6], [280, 2e6], [380, 1.2e7], [330, 1.6e7]]) {
    const f = flashPT({ mix: oilMix, T, P, z: zOil });
    assert.equal(f.phases, 2, `two-phase at ${T} K, ${P} Pa`);
    assert.ok(f.converged);
    const fl = fugacityCoefficients({ mix: oilMix, T, P, z: f.x, phase: 'stable' });
    const fv = fugacityCoefficients({ mix: oilMix, T, P, z: f.y, phase: 'stable' });
    for (let i = 0; i < zOil.length; i++) {
      assert.ok(Math.abs(f.beta * f.y[i] + (1 - f.beta) * f.x[i] - zOil[i]) < 1e-10, 'material balance');
      assert.ok(rel(fl.fugacity[i], fv.fugacity[i]) < 1e-8, 'isofugacity');
    }
    assert.ok(Math.abs(f.x.reduce((a, b) => a + b) - 1) < 1e-12 && Math.abs(f.y.reduce((a, b) => a + b) - 1) < 1e-12);
    assert.ok(f.ZL < f.ZV);
  }
});

test('stability test: single phase far outside the envelope, flash returns the feed', () => {
  assert.ok(stabilityTest({ mix: oilMix, T: 320, P: 5e7, z: zOil }).stable);
  assert.ok(stabilityTest({ mix: oilMix, T: 700, P: 1e5, z: zOil }).stable);
  assert.ok(!stabilityTest({ mix: oilMix, T: 320, P: 5e6, z: zOil }).stable);
  const liq = flashPT({ mix: oilMix, T: 320, P: 5e7, z: zOil });
  assert.equal(liq.phases, 1);
  assert.equal(liq.phase, 'liquid');
  const vap = flashPT({ mix: oilMix, T: 700, P: 1e5, z: zOil });
  assert.equal(vap.phase, 'vapour');
  assert.equal(vap.beta, 1);
});

test('bubble point of a methane–n-butane binary: direct solve vs flash', () => {
  const mix = makeMixture(['C1', 'nC4']);
  const z = [0.3, 0.7];
  const T = 311;
  const bp = bubblePointPressure({ mix, T, z });
  assert.ok(bp.converged);
  assert.ok(Math.abs(bp.y.reduce((a, b) => a + b) - 1) < 1e-9);
  // flash just inside the envelope has a vanishing vapour fraction whose composition is the incipient vapour
  const inside = flashPT({ mix, T, P: bp.P * (1 - 1e-4), z });
  assert.equal(inside.phases, 2);
  assert.ok(inside.beta > 0 && inside.beta < 1e-3, `beta ${inside.beta}`);
  assert.ok(Math.abs(inside.y[0] - bp.y[0]) < 1e-3);
  assert.equal(flashPT({ mix, T, P: bp.P * (1 + 1e-4), z }).phases, 1);
  // bisection on the flash vapour fraction reproduces the direct solve
  let lo = bp.P * 0.9;
  let hi = bp.P * 1.1;
  for (let it = 0; it < 40; it++) {
    const m = 0.5 * (lo + hi);
    if (flashPT({ mix, T, P: m, z }).phases === 2) lo = m; else hi = m;
  }
  assert.ok(rel(0.5 * (lo + hi), bp.P) < 1e-6);
  // dew point: lower boundary, vapour fraction → 1
  const dp = dewPointPressure({ mix, T, z });
  assert.ok(dp.P < bp.P);
  const nearDew = flashPT({ mix, T, P: dp.P * (1 + 1e-4), z });
  assert.ok(nearDew.phases === 2 && nearDew.beta > 0.999);
});

test('saturation temperature solvers invert the saturation pressure solvers', () => {
  const T = 300;
  const bp = bubblePointPressure({ mix: oilMix, T, z: zOil });
  assert.ok(Math.abs(bubblePointTemperature({ mix: oilMix, P: bp.P, z: zOil }).T - T) < 1e-4);
  const dp = dewPointPressure({ mix: oilMix, T: 450, z: zOil });
  assert.ok(Math.abs(dewPointTemperature({ mix: oilMix, P: dp.P, z: zOil }).T - 450) < 1e-4);
  // near-critical bubble point falls back to the stability bisection and agrees with the upper boundary
  const up = upperSaturationPressure({ mix: oilMix, T: 320, z: zOil });
  assert.equal(up.type, 'bubble');
  assert.ok(rel(bubblePointPressure({ mix: oilMix, T: 320, z: zOil }).P, up.P) < 1e-5);
});

test('caloric properties are thermodynamically consistent (Cp = ∂h/∂T, μJT = −(∂h/∂P)/Cp, ideal-gas limits)', () => {
  const mix = makeMixture(['C1', 'C2', 'CO2', 'N2']);
  const z = [0.85, 0.08, 0.04, 0.03];
  const T = 310;
  const P = 8e6;
  const p = properties({ mix, T, P, z });
  const dT = 1e-3;
  const dP = 10;
  const hT = (properties({ mix, T: T + dT, P, z }).h - properties({ mix, T: T - dT, P, z }).h) / (2 * dT);
  const hP = (properties({ mix, T, P: P + dP, z }).h - properties({ mix, T, P: P - dP, z }).h) / (2 * dP);
  const sT = (properties({ mix, T: T + dT, P, z }).s - properties({ mix, T: T - dT, P, z }).s) / (2 * dT);
  assert.ok(rel(hT, p.cp) < 1e-6, `cp ${p.cp} vs dh/dT ${hT}`);
  assert.ok(rel(-hP / p.cp, p.muJT) < 1e-5, `muJT ${p.muJT}`);
  assert.ok(rel(sT, p.cp / T) < 1e-6);
  assert.ok(p.muJT > 2e-6 && p.muJT < 6e-6); // natural gas cools 2–6 K per MPa at these conditions
  // ideal-gas limit: Z → 1, cp − cv → R, w² → γ R T/M
  const ig = properties({ mix: makeMixture(['N2']), T: 300, P: 100, z: [1] });
  assert.ok(Math.abs(ig.Z - 1) < 1e-5);
  assert.ok(rel(ig.cp - ig.cv, R) < 1e-4);
  assert.ok(rel(ig.soundSpeed, Math.sqrt(1.4 * R * 300 / 0.0280134)) < 2e-3); // 353 m/s
  assert.ok(Math.abs(ig.hDep) < 1 && Math.abs(ig.muJT) < 1e-5);
});

test('Péneloux-type volume shift improves liquid densities and leaves phase equilibrium unchanged', () => {
  // saturated liquid propane at 300 K: 489 kg/m³ (NIST); n-decane at 293 K, 1 atm: 730 kg/m³
  for (const [id, T, P, rho] of [['C3', 300, 1.1e6, 489], ['nC10', 293.15, 101325, 730]]) {
    const plain = properties({ mix: makeMixture([id]), T, P, z: [1] }).rho;
    const shifted = properties({ mix: makeMixture([id], { volumeShift: true }), T, P, z: [1] }).rho;
    assert.ok(Math.abs(shifted - rho) < Math.abs(plain - rho), `${id}: ${plain} → ${shifted}`);
    assert.ok(rel(shifted, rho) < 0.025, `${id}: ${shifted}`);
  }
  const a = flashPT({ mix: oilMix, T: 320, P: 5e6, z: zOil });
  const b = flashPT({ mix: makeMixture(oilMix.ids, { volumeShift: true }), T: 320, P: 5e6, z: zOil });
  assert.ok(Math.abs(a.beta - b.beta) < 1e-12);
  // SRK with the Péneloux correlation
  const srk = properties({ mix: makeMixture(['C3'], { eos: 'SRK', volumeShift: true }), T: 300, P: 1.1e6, z: [1] }).rho;
  assert.ok(rel(srk, 489) < 0.04, `SRK shifted ${srk}`);
});

test('phase envelope of a natural gas: closed, with cricondenbar/cricondentherm bounding the critical point', () => {
  const mix = makeMixture(['N2', 'CO2', 'C1', 'C2', 'C3', 'nC4', 'nC5', 'nC6']);
  const z = [0.01, 0.02, 0.85, 0.06, 0.03, 0.015, 0.01, 0.005];
  const env = phaseEnvelope({ mix, z });
  assert.ok(env.complete);
  assert.ok(env.bubble.length > 10 && env.dew.length > 10);
  assert.ok(env.cricondenbar.P >= env.critical.P && env.cricondentherm.T >= env.critical.T);
  // retrograde gas: critical point on the low-temperature side of the cricondenbar
  assert.ok(env.critical.T < env.cricondenbar.T && env.cricondenbar.T < env.cricondentherm.T);
  // traced points are true phase boundaries: two phases just inside, one phase just outside
  const b = env.bubble[Math.floor(env.bubble.length / 2)];
  assert.equal(flashPT({ mix, T: b.T, P: b.P * 0.98, z }).phases, 2);
  assert.equal(flashPT({ mix, T: b.T, P: b.P * 1.02, z }).phases, 1);
  const d = env.dew[3];
  const fd = flashPT({ mix, T: d.T, P: d.P * 1.05, z });
  assert.ok(fd.phases === 2 && fd.beta > 0.97);
  assert.equal(flashPT({ mix, T: d.T, P: d.P * 0.95, z }).phases, 1);
  // outside the cricondentherm and above the cricondenbar the fluid is single phase
  assert.equal(flashPT({ mix, T: env.cricondentherm.T + 3, P: env.cricondentherm.P, z }).phases, 1);
  assert.equal(flashPT({ mix, T: env.cricondenbar.T, P: env.cricondenbar.P * 1.03, z }).phases, 1);
  assert.equal(flashPT({ mix, T: env.cricondenbar.T, P: env.cricondenbar.P * 0.97, z }).phases, 2);
  // approaching the interpolated critical point along the critical isotherm the two phases become identical
  // (K → 1) and the vapour fraction stays intermediate, as it must at a critical point
  const spread = (fac) => {
    const fl = flashPT({ mix, T: env.critical.T, P: env.critical.P * fac, z });
    assert.equal(fl.phases, 2);
    assert.ok(fl.beta > 0.3 && fl.beta < 0.7, `beta ${fl.beta}`);
    return Math.max(...fl.K.map((k) => Math.abs(Math.log(k))));
  };
  const s90 = spread(0.9);
  const s97 = spread(0.97);
  const s99 = spread(0.99);
  assert.ok(s99 < s97 && s97 < s90 && s99 < 0.35 * s90, `${s90} ${s97} ${s99}`);
});

test('phase envelope of a volatile oil closes and agrees with the fixed-temperature solvers', () => {
  const env = phaseEnvelope({ mix: oilMix, z: zOil });
  assert.ok(env.complete);
  const bp = bubblePointPressure({ mix: oilMix, T: 300, z: zOil }).P;
  // interpolate the traced bubble line at 300 K
  let k = 1;
  while (env.bubble[k].T < 300) k++;
  const w = (300 - env.bubble[k - 1].T) / (env.bubble[k].T - env.bubble[k - 1].T);
  const lnP = Math.log(env.bubble[k - 1].P) + w * (Math.log(env.bubble[k].P) - Math.log(env.bubble[k - 1].P));
  assert.ok(rel(Math.exp(lnP), bp) < 0.01);
  assert.throws(() => phaseEnvelope({ mix: makeMixture(['C1']), z: [1] }), RangeError);
});

test('transport properties: LGE, LBC dilute limit and liquid magnitude, Eucken, Cragoe', () => {
  assert.ok(rel(gasViscosityLGE({ T: 300, rho: 0.65, M: 0.016043 }), 1.12e-5) < 0.05);
  const c1 = makeMixture(['C1']);
  // dilute methane at 300 K: 1.12e-5 Pa·s (Stiel–Thodos within 5 %)
  assert.ok(rel(viscosityLBC({ mix: c1, T: 300, z: [1], rhoMolar: 1e-3 }), 1.12e-5) < 0.05);
  // liquid n-decane at 293 K: 0.92 mPa·s measured; LBC with the measured density (730 kg/m³) is within a factor 1.5
  const c10 = makeMixture(['nC10']);
  const mu = viscosityLBC({ mix: c10, T: 293.15, z: [1], rhoMolar: 730 / 0.142285 });
  assert.ok(mu > 0.92e-3 / 1.5 && mu < 0.92e-3 * 1.5, `LBC n-decane ${mu}`);
  // viscosity rises monotonically with density
  assert.ok(viscosityLBC({ mix: c1, T: 300, z: [1], rhoMolar: 8000 }) > viscosityLBC({ mix: c1, T: 300, z: [1], rhoMolar: 4000 }));
  // modified Eucken for methane at 300 K: k = (μ/M)(1.32 Cv + 1.77 R) ≈ 0.034 W/(m·K) (measured 0.034)
  const k = gasThermalConductivity({ mu: 1.12e-5, cv: 35.7 - R, M: 0.016043 });
  assert.ok(rel(k, (1.12e-5 / 0.016043) * (1.32 * (35.7 - R) + 1.77 * R)) < 1e-12);
  assert.ok(rel(k, 0.0343) < 0.08);
  // Cragoe: 0.1172/γ·(1 − 0.00054 t°C) at 20 °C, γ = 0.85
  assert.ok(rel(oilThermalConductivity({ T: 293.15, sg: 0.85 }), 0.1172 / 0.85 * (1 - 0.00054 * 20)) < 2e-3);
});

test('parachor surface tension — hand calculation and vanishing at equal phases', () => {
  // propane at 300 K: ρL = 489 kg/m³, ρV = 21.6 kg/m³ → σ^¼ = 150.3·(489 − 21.6)/44.097e3 = 1.593 → 6.4 mN/m
  // (measured 6.6 mN/m)
  const mix = makeMixture(['C3']);
  const s = surfaceTension({ mix, x: [1], y: [1], rhoMolarL: 489 / 0.044097, rhoMolarV: 21.6 / 0.044097 });
  assert.ok(rel(s, (150.3 * (489 - 21.6) / 44.097e3) ** 4 * 1e-3) < 1e-12);
  assert.ok(rel(s, 6.6e-3) < 0.06);
  assert.equal(surfaceTension({ mix, x: [1], y: [1], rhoMolarL: 5000, rhoMolarV: 5000 }), 0);
  assert.throws(() => surfaceTension({ mix: makeMixture(['H2O']), x: [1], y: [1], rhoMolarL: 5e4, rhoMolarV: 10 }), RangeError);
});

test('Bukacek water content at 100 °F, 1000 psia ≈ 60 lb/MMscf (McKetta–Wehe chart, GPSA Fig. 20-4)', () => {
  const w = waterContentBukacek({ T: (100 + 459.67) / 1.8, P: 1000 * 6894.757293168 });
  assert.ok(rel(w.lbPerMMscf, 60) < 0.05, `${w.lbPerMMscf}`);
  // hand: 47484·0.9503/1000 + 10^(−3083.87/559.67 + 6.69449) = 45.12 + 15.29
  assert.ok(rel(w.lbPerMMscf, 47484 * 0.9503 / 1000 + 10 ** (-3083.87 / 559.67 + 6.69449)) < 2e-3);
  assert.ok(rel(w.kgPerSm3, w.lbPerMMscf * 16.0185e-6) < 1e-5);
});

test('standard-condition conversions', () => {
  assert.ok(rel(standardMolarVolume(), 0.0236445) < 1e-4); // 23.645 m³/kmol at 15 °C, 101.325 kPa
  assert.ok(rel(standardMolarVolume(STANDARD.field) * 453.59237 / 0.028316846592, 379.49) < 1e-4); // scf/lbmol
  assert.ok(rel(standardMolarVolume(STANDARD.normal), 0.022414) < 1e-4);
  const m = standardRateToMass({ qStd: 10, M: 0.018 });
  assert.ok(rel(massToStandardRate({ mass: m.mass, M: 0.018 }), 10) < 1e-12);
  assert.ok(rel(mmscfdToSm3s(1), 0.3277413) < 1e-6);
});

test('input validation', () => {
  assert.throws(() => makeMixture(['C1'], { eos: 'vdW' }), RangeError);
  assert.throws(() => flashPT({ mix: oilMix, T: 300, P: 1e6, z: [0.5, 0.5] }), RangeError);
  assert.throws(() => flashPT({ mix: oilMix, T: 300, P: 1e6, z: [0.5, 0.5, 0.5, 0.5, 0.5] }), RangeError);
  assert.throws(() => zRoots({ mix: oilMix, T: -5, P: 1e6, z: zOil }), RangeError);
});
