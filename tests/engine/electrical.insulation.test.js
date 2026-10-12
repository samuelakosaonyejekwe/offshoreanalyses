import test from 'node:test';
import assert from 'node:assert/strict';
import {
  acStress, stressCheck, dielectricProperties, dcField, dcPolarityReversalStress, thermalLife, electricalLife,
  combinedLife, lossOfLife, weibullBreakdownProbability, weibullEnlargement, IMPULSE_LEVELS, impulseLevels,
} from '../../src/engine/electrical/insulation.js';

const rel = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${a} vs ${b} (tol ${tol})`);

test('coaxial AC stress by hand, 132 kV (U0 = 76.21 kV), r = 20 mm, R = 36 mm: Emax = U0/(r·ln 1.8) = 6.483 kV/mm, Emin = 3.602 kV/mm', () => {
  const U0 = 132e3 / Math.sqrt(3);
  const s = acStress({ U0, ri: 0.020, ro: 0.036, r: 0.028 });
  rel(s.Emax, 6.4828e6, 1e-4); rel(s.Emin, 3.6016e6, 1e-4);
  rel(s.Emax / s.Emin, 1.8, 1e-12);
  rel(s.E, U0 / (0.028 * Math.log(1.8)), 1e-12);
  rel(s.Emean, U0 / 0.016, 1e-12);
  // the field integrates to the voltage: ∫ U0/(r ln(R/ri)) dr = U0
  let integral = 0; const n = 2000;
  for (let i = 0; i < n; i++) { const r = 0.020 + (0.016 * (i + 0.5)) / n; integral += (acStress({ U0, ri: 0.020, ro: 0.036, r }).E * 0.016) / n; }
  rel(integral, U0, 1e-6);
});

test('design-stress check and the thickness that meets a stress limit: t = ri·(exp(U0/(ri·Elim)) − 1)', () => {
  const U0 = 76210;
  const c = stressCheck({ U0, ri: 0.020, ro: 0.036, limitInner: 6e6, limitOuter: 4e6 });
  assert.equal(c.ok, false); rel(c.utilisationInner, 6.4828 / 6, 1e-4);
  const again = acStress({ U0, ri: 0.020, ro: 0.020 + c.thicknessRequired });
  rel(again.Emax, 6e6, 1e-10);
});

test('capacitance, insulation resistance and self-discharge time constant of a coaxial core', () => {
  const d = dielectricProperties({ ri: 0.020, ro: 0.036, epsR: 2.5, resistivity: 1e14, tanDelta: 0.001, U0: 76210, f: 50 });
  rel(d.C, 2.36618e-10, 1e-5);
  rel(d.Rins, (1e14 * Math.log(1.8)) / (2 * Math.PI), 1e-12);
  // R·C of the coaxial geometry equals ρ·ε regardless of dimensions
  rel(d.Rins * d.C, 1e14 * 8.8541878128e-12 * 2.5, 1e-12);
  rel(d.timeConstant, d.Rins * d.C, 1e-12);
  rel(d.Wd, 2 * Math.PI * 50 * d.C * 76210 ** 2 * 0.001, 1e-12);
  rel(d.lossAngle, Math.atan(0.001), 1e-12);
});

test('DC field without temperature drop or field dependence is the capacitive (Laplace) field', () => {
  const f = dcField({ U: 320e3, ri: 0.02, ro: 0.04, alpha: 0.08, beta: 0, deltaT: 0 });
  rel(f.Einner, 320e3 / (0.02 * Math.LN2), 1e-5);
  rel(f.Eouter, 320e3 / (0.04 * Math.LN2), 1e-5);
  assert.equal(f.inverted, false);
});

test('DC field with temperature-dependent conductivity (β = 0) matches the closed form E = δ·U·r^(δ−1)/(ro^δ − ri^δ), δ = α·ΔT/ln(ro/ri)', () => {
  const p = { U: 320e3, ri: 0.02, ro: 0.04, alpha: 0.1, beta: 0, deltaT: 15, n: 801 };
  const f = dcField(p);
  const delta = (0.1 * 15) / Math.LN2; // 2.164 > 1: stress inversion
  const E = (r) => (delta * p.U * r ** (delta - 1)) / (p.ro ** delta - p.ri ** delta);
  rel(f.Einner, E(0.02), 2e-5); rel(f.Eouter, E(0.04), 2e-5);
  rel(f.E[400], E(f.r[400]), 2e-5);
  assert.equal(f.inverted, true);
  rel(f.deltaEoll, delta, 1e-12);
  assert.equal(f.rMax, f.r[800]);
  // δ = 1 exactly gives a uniform field U/(ro − ri)
  const flat = dcField({ ...p, deltaT: Math.LN2 / 0.1 });
  rel(flat.Einner, 320e3 / 0.02, 1e-5); rel(flat.Eouter, 320e3 / 0.02, 1e-5);
});

test('DC field with field-dependent conductivity: voltage integral and current continuity r·σ(T,E)·E = const hold; field is flattened', () => {
  const p = { U: 320e3, ri: 0.02, ro: 0.04, alpha: 0.08, beta: 5e-8, deltaT: 10, n: 801 };
  const f = dcField(p);
  let integral = 0;
  for (let i = 1; i < f.r.length; i++) integral += 0.5 * (f.E[i] + f.E[i - 1]) * (f.r[i] - f.r[i - 1]);
  rel(integral, 320e3, 1e-9);
  const K = (i) => f.r[i] * f.E[i] * Math.exp(-0.08 * 10 * Math.log(f.r[i] / 0.02) / Math.LN2 + 5e-8 * f.E[i]);
  rel(K(0), K(400), 1e-9); rel(K(800), K(400), 1e-9);
  // compared with β = 0 the field dependence reduces the spread between the screens
  const lin = dcField({ ...p, beta: 0 });
  assert.ok(Math.abs(f.Eouter - f.Einner) < Math.abs(lin.Eouter - lin.Einner));
  // no load, with field dependence: the inner stress is below the Laplace value (field grading)
  const cold = dcField({ ...p, deltaT: 0 });
  assert.ok(cold.Einner < cold.Elaplace && cold.Einner > cold.Eouter);
});

test('polarity reversal: capacitive field of the 2U step superposed on the resistive field', () => {
  const U = 320e3; const ri = 0.02; const ro = 0.04;
  const f = dcField({ U, ri, ro, alpha: 0.1, beta: 0, deltaT: 15 });
  const pr = dcPolarityReversalStress({ U, ri, ro, Einner: f.Einner, Eouter: f.Eouter });
  rel(pr.Einner, Math.abs(f.Einner - (2 * U) / (ri * Math.LN2)), 1e-12);
  // with an inverted (loaded) field the conductor-screen stress after reversal exceeds the no-load Laplace stress
  assert.ok(pr.Einner > f.Elaplace);
});

test('thermal life: Montsinger halving interval and Arrhenius by hand (Ea = 1.1 eV, 90 → 100 °C: factor exp(−0.9420) = 0.3898)', () => {
  rel(thermalLife({ thetaC: 100, thetaRefC: 90, lifeRef: 40, halvingIntervalK: 10 }), 20, 1e-12);
  rel(thermalLife({ thetaC: 70, thetaRefC: 90, lifeRef: 40, halvingIntervalK: 10 }), 160, 1e-12);
  // Ea/k = 1.1/8.617333e-5 = 12764.97 K; 1/373.15 − 1/363.15 = −7.3795e-5 1/K
  rel(thermalLife({ thetaC: 100, thetaRefC: 90, lifeRef: 1, activationEnergyEv: 1.1 }), Math.exp(12764.97 * (1 / 373.15 - 1 / 363.15)), 1e-6);
  rel(thermalLife({ thetaC: 100, thetaRefC: 90, lifeRef: 1, activationEnergyEv: 1.1 }), 0.3898, 2e-3);
  assert.throws(() => thermalLife({ thetaC: 100, thetaRefC: 90, lifeRef: 1 }), RangeError);
});

test('inverse power law: 10 % more stress with n = 10 divides the life by 1.1^10 = 2.594; combined model is the product', () => {
  rel(electricalLife({ E: 11, Eref: 10, lifeRef: 40, n: 10 }), 40 / 2.5937425, 1e-7);
  rel(combinedLife({ E: 11, Eref: 10, n: 10, thetaC: 100, thetaRefC: 90, lifeRef: 40, halvingIntervalK: 10 }), 20 / 2.5937425, 1e-7);
});

test('loss of life over a temperature history: constant reference temperature consumes duration/life; a hot spell is weighted', () => {
  const hours = 8760;
  const flat = lossOfLife({ thetaC: new Array(hours).fill(90), dt: 1, thetaRefC: 90, lifeRef: 40 * 8760, halvingIntervalK: 10 });
  rel(flat.lossOfLife, 1 / 40, 1e-12); rel(flat.equivalentAgeingFactor, 1, 1e-12);
  // half the year at 80 °C (rate 0.5) and half at 100 °C (rate 2): mean ageing factor 1.25
  const mix = lossOfLife({ thetaC: [...new Array(4380).fill(80), ...new Array(4380).fill(100)], dt: 1, thetaRefC: 90, lifeRef: 40 * 8760, halvingIntervalK: 10 });
  rel(mix.equivalentAgeingFactor, 1.25, 1e-12);
  rel(mix.lifeAtThisDuty, (40 * 8760) / 1.25, 1e-12);
});

test('Weibull enlargement law: b = 10, length ×1000 → scale stress ×1000^(−1/10) = 0.5012; probabilities consistent', () => {
  const e = weibullEnlargement({ E0: 60e6, b: 10, sizeRatio: 1000, probability: 0.01 });
  rel(e.E632 / 60e6, 0.501187, 1e-6);
  rel(weibullBreakdownProbability({ E: e.E632, E0: 60e6, b: 10, sizeRatio: 1000 }), 1 - Math.exp(-1), 1e-12);
  rel(weibullBreakdownProbability({ E: e.Eprobability, E0: 60e6, b: 10, sizeRatio: 1000 }), 0.01, 1e-10);
  // reference sample at its own scale stress: 63.2 %
  rel(weibullBreakdownProbability({ E: 60e6, E0: 60e6, b: 10 }), 0.6321206, 1e-6);
  // time term: ten times longer with a = 1 is the same as ten times the length
  rel(weibullBreakdownProbability({ E: 30e6, E0: 60e6, b: 10, timeRatio: 10 }), weibullBreakdownProbability({ E: 30e6, E0: 60e6, b: 10, sizeRatio: 10 }), 1e-12);
});

test('impulse withstand levels by Um (IEC 60071-1 / 60502-2 / 60840 / 62067) and arrester co-ordination margin', () => {
  const bil = (um) => impulseLevels({ Um: um }).lightning;
  assert.equal(bil(36e3), 170e3); assert.equal(bil(72.5e3), 325e3); assert.equal(bil(145e3), 650e3);
  assert.equal(bil(245e3), 1050e3); assert.equal(bil(420e3), 1425e3);
  assert.equal(impulseLevels({ Um: 420e3 }).switching, 1050e3);
  // a system voltage between table rows takes the next row up
  assert.equal(impulseLevels({ Um: 66e3 }).Um, 72.5e3);
  const m = impulseLevels({ Um: 145e3, arresterLightningProtectiveLevel: 500e3 });
  rel(m.lightningMargin, 1.3, 1e-12); assert.equal(m.ok, true);
  assert.equal(impulseLevels({ Um: 145e3, arresterLightningProtectiveLevel: 600e3 }).ok, false);
  // the table is editable through the argument
  assert.equal(impulseLevels({ Um: 145e3, table: [{ Um: 145e3, lightning: 550e3, switching: null }] }).lightning, 550e3);
  assert.ok(Object.isFrozen(IMPULSE_LEVELS));
  assert.throws(() => impulseLevels({ Um: 800e3 }), RangeError);
});
