import test from 'node:test';
import assert from 'node:assert/strict';
import {
  R, COMPONENTS, COMPONENT_IDS, getComponent, makeComponent, idealGasCp, idealGasEnthalpyEntropy,
  boilingPointSoreide, keslerLee, riaziDaubert, characterizePlusFraction, defaultKij, defaultKijMatrix,
  mixtureMolarMass,
} from '../../src/engine/flow/components.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('database holds the 21 required components with physically ordered constants', () => {
  const ids = ['N2', 'CO2', 'H2S', 'C1', 'C2', 'C3', 'iC4', 'nC4', 'iC5', 'nC5', 'nC6', 'nC7', 'nC8', 'nC9', 'nC10',
    'H2O', 'methanol', 'MEG', 'H2', 'O2', 'Ar'];
  assert.deepEqual([...COMPONENT_IDS], ids);
  for (const id of ids) {
    const c = COMPONENTS[id];
    assert.ok(c.M > 0 && c.Tc > 0 && c.Pc > 0 && c.Tb > 0 && c.Tb < c.Tc, id);
  }
  // n-alkane series: Tc and ω increase, Pc decreases with carbon number from C2
  const alk = ['C2', 'C3', 'nC4', 'nC5', 'nC6', 'nC7', 'nC8', 'nC9', 'nC10'].map(getComponent);
  for (let i = 1; i < alk.length; i++) {
    assert.ok(alk[i].Tc > alk[i - 1].Tc && alk[i].omega > alk[i - 1].omega && alk[i].Pc < alk[i - 1].Pc);
    // CH2 increment of the molar mass: 14.027 g/mol
    assert.ok(Math.abs(alk[i].M - alk[i - 1].M - 0.014027) < 2e-6);
  }
  assert.throws(() => getComponent('unobtainium'), RangeError);
});

test('acentric factor definition ω = −log10(Psat/Pc at Tr = 0.7) − 1 is met by the boiling-point data (Edmister check)', () => {
  // Edmister estimate ω ≈ (3/7)·θ/(1 − θ)·log10(Pc/atm) − 1 with θ = Tb/Tc reproduces tabulated ω to ~0.02
  for (const id of ['C1', 'C3', 'nC5', 'nC8', 'nC10', 'N2', 'CO2']) {
    if (id === 'CO2') continue; // CO2 sublimes at 1 atm: Tb is not a liquid boiling point
    const c = getComponent(id);
    const th = c.Tb / c.Tc;
    const w = (3 / 7) * th / (1 - th) * Math.log10(c.Pc / 101325) - 1;
    assert.ok(Math.abs(w - c.omega) < 0.025, `${id}: Edmister ${w} vs ${c.omega}`);
  }
});

test('ideal-gas heat capacities at 298.15 K match standard tabulated values', () => {
  // JANAF / TRC values in J/(mol·K)
  const ref = { N2: 29.12, O2: 29.38, H2: 28.84, CO2: 37.13, H2O: 33.59, C1: 35.69, C2: 52.49, C3: 73.60, nC4: 98.49, nC10: 233.1 };
  for (const [id, cp] of Object.entries(ref)) {
    assert.ok(rel(idealGasCp(id, 298.15), cp) < 0.01, `${id}: ${idealGasCp(id, 298.15)}`);
  }
  assert.ok(rel(idealGasCp('Ar', 500), 2.5 * R) < 1e-14); // monatomic gas
  assert.throws(() => idealGasCp('MEG', 300), RangeError);
});

test('ideal-gas enthalpy and entropy integrals are consistent with Cp', () => {
  const T = 420;
  const h = 1e-2;
  const a = idealGasEnthalpyEntropy('C3', T + h);
  const b = idealGasEnthalpyEntropy('C3', T - h);
  assert.ok(rel((a.h - b.h) / (2 * h), idealGasCp('C3', T)) < 1e-7);
  assert.ok(rel((a.s - b.s) / (2 * h), idealGasCp('C3', T) / T) < 1e-7);
  assert.equal(idealGasEnthalpyEntropy('C3', 298.15).h, 0);
});

test('Kesler–Lee and Riazi–Daubert reproduce n-decane critical properties from Tb and SG', () => {
  // n-decane: Tb = 447.3 K, SG = 0.734; Tc = 617.7 K, Pc = 21.1 bar, ω = 0.492
  const kl = keslerLee({ Tb: 447.3, sg: 0.734 });
  assert.ok(rel(kl.Tc, 617.7) < 0.01);
  assert.ok(rel(kl.Pc, 21.1e5) < 0.05);
  assert.ok(Math.abs(kl.omega - 0.4923) < 0.01);
  const rd = riaziDaubert({ Tb: 447.3, sg: 0.734 });
  assert.ok(rel(rd.Tc, 617.7) < 0.01);
  assert.ok(rel(rd.Pc, 21.1e5) < 0.06);
  // the two independent molar-mass correlations agree with each other within 3 %
  assert.ok(rel(kl.M, rd.M) < 0.03);
});

test('Kesler–Lee Tc — hand evaluation for Tb = 900 °R, γ = 0.80', () => {
  // Tc = 341.7 + 811.1·0.8 + (0.4244 + 0.1174·0.8)·900 + (0.4669 − 3.26238·0.8)·1e5/900
  //    = 341.7 + 648.88 + 466.488 − 238.114 = 1218.95 °R
  const kl = keslerLee({ Tb: 500, sg: 0.8 });
  assert.ok(rel(kl.Tc * 1.8, 1218.95) < 1e-5);
  // Watson factor Kw = Tb^(1/3)/γ = 900^(1/3)/0.8 = 12.068
  assert.ok(rel(kl.Kw, 12.0683) < 1e-4);
});

test('Søreide boiling point of n-decane from M and SG is within 1.5 % of 447.3 K', () => {
  assert.ok(rel(boilingPointSoreide(0.142285, 0.734), 447.3) < 0.015);
});

test('plus-fraction characterisation: LBC critical volume by hand, monotonic trends', () => {
  // Vc = 21.573 + 0.015122·200 − 27.656·0.85 + 0.070615·200·0.85 = 13.0943 ft³/lbmol = 8.1745e-4 m³/mol
  const c = characterizePlusFraction({ M: 0.2, sg: 0.85 });
  assert.ok(rel(c.Vc, 13.0943 * 0.028316846592 / 453.59237) < 1e-5);
  assert.ok(c.pseudo);
  assert.ok(rel(c.parachor, 59.3 + 2.34 * 200) < 1e-12);
  const heavy = characterizePlusFraction({ M: 0.35, sg: 0.9 });
  assert.ok(heavy.Tc > c.Tc && heavy.Pc < c.Pc && heavy.omega > c.omega && heavy.Tb > c.Tb);
  // mass-specific ideal-gas Cp equals that of n-decane
  assert.ok(rel(idealGasCp(c, 350) / c.M, idealGasCp('nC10', 350) / COMPONENTS.nC10.M) < 1e-12);
  const rd = characterizePlusFraction({ M: 0.2, sg: 0.85, method: 'riaziDaubert' });
  assert.ok(rel(rd.Tc, c.Tc) < 0.02);
});

test('default binary interaction parameters: symmetric, zero for hydrocarbon pairs', () => {
  assert.equal(defaultKij('N2', 'C1'), 0.025);
  assert.equal(defaultKij('C1', 'N2'), 0.025);
  assert.equal(defaultKij('CO2', 'C1'), 0.105);
  assert.equal(defaultKij('H2S', 'nC10'), 0.050);
  assert.equal(defaultKij('C1', 'nC10'), 0);
  assert.equal(defaultKij('CO2', 'CO2'), 0);
  const plus = characterizePlusFraction({ M: 0.2, sg: 0.85 });
  assert.equal(defaultKij(plus, 'CO2'), 0.115);
  const k = defaultKijMatrix(['N2', 'CO2', 'C1', 'nC5']);
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) assert.equal(k[i][j], k[j][i]);
});

test('component editing and mixture molar mass', () => {
  const c = makeComponent({ base: 'C1', Tc: 191 });
  assert.equal(c.Tc, 191);
  assert.equal(c.M, COMPONENTS.C1.M);
  assert.throws(() => makeComponent({ id: 'x', M: 0.1, Tc: 500, Pc: -1, omega: 0.3 }), RangeError);
  assert.ok(rel(mixtureMolarMass(['C1', 'C2'], [0.9, 0.1]), 0.9 * 0.016043 + 0.1 * 0.03007) < 1e-12);
});
