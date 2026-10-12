import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  hardeningStress, plasticity1d, stressStrainHistory1d, vonMisesStress, plasticity3d, uniaxialResponse3d, rambergOsgoodStrain,
  rambergOsgoodStress, rambergOsgoodExponent, johnsonCook, cowperSymonds, trueFromEngineering, engineeringFromTrue,
  mooneyRivlinUniaxial, neoHookeanUniaxial, manderConfinement, manderConfinedConcrete,
} from '../../src/engine/structural/materials-nl.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b} (abs ${abs})`);
const E = 200e9; const sy0 = 300e6;

test('1-D plasticity: bilinear monotonic curve σ = σy + Et(ε − εy), Et = E H/(E + H); perfect plasticity', () => {
  const H = 2e9;
  const r = plasticity1d({ strain: 0.01, E, sy0, H });
  const Et = E * H / (E + H);
  close(r.stress, sy0 + Et * (0.01 - sy0 / E), 1e-12); close(r.tangent, Et, 1e-12);
  close(r.state.alpha, 0.01 - r.stress / E, 1e-12); close(r.state.epsP, r.state.alpha, 1e-12);
  // consistency: the updated stress lies on the updated yield surface
  close(r.stress, hardeningStress({ alpha: r.state.alpha, sy0, H }).sy, 1e-12);
  const el = plasticity1d({ strain: 0.001, E, sy0, H });
  assert.equal(el.yielded, false); close(el.stress, 200e6, 1e-12); assert.equal(el.tangent, E);
  const pp = plasticity1d({ strain: 0.02, E, sy0 });
  close(pp.stress, sy0, 1e-12); assert.equal(pp.tangent, 0);
  // step-size independence for proportional loading: one step or 100 steps give the same result
  const hist = stressStrainHistory1d({ strains: Array.from({ length: 100 }, (_, i) => 0.01 * (i + 1) / 100), E, sy0, H });
  close(hist.stress[99], r.stress, 1e-12); close(hist.state.alpha, r.state.alpha, 1e-10);
});

test('1-D plasticity, cyclic: isotropic hardening expands the surface, kinematic hardening gives the Bauschinger effect (reverse yield at σmax − 2σy)', () => {
  const H = 5e9;
  const fwd = 0.01;
  const iso = plasticity1d({ strain: fwd, E, sy0, H });
  const kin = plasticity1d({ strain: fwd, E, sy0, Hkin: H });
  close(iso.stress, kin.stress, 1e-12); // same monotonic curve
  // unload elastically: reverse yielding starts when σ = −σmax (isotropic) or σ = σmax − 2σy0 (kinematic)
  const revIso = fwd - 2 * iso.stress / E; const revKin = fwd - 2 * sy0 / E;
  assert.equal(plasticity1d({ strain: revIso * (1 + 1e-9) + 1e-9, E, sy0, H, state: iso.state }).yielded, false);
  assert.equal(plasticity1d({ strain: revIso - 1e-6, E, sy0, H, state: iso.state }).yielded, true);
  const k1 = plasticity1d({ strain: revKin + 1e-9, E, sy0, Hkin: H, state: kin.state });
  assert.equal(k1.yielded, false); close(k1.stress, kin.stress - 2 * sy0, 1e-5);
  assert.equal(plasticity1d({ strain: revKin - 1e-6, E, sy0, Hkin: H, state: kin.state }).yielded, true);
  // closed stabilised loop under kinematic hardening: symmetric strain cycling returns to the same stress, dissipation = loop area
  const cyc = []; const amp = 0.006;
  for (let c = 0; c < 3; c++) for (let i = 0; i < 400; i++) cyc.push(amp * Math.sin(2 * Math.PI * (i + 1) / 400));
  const h = stressStrainHistory1d({ strains: cyc, E, sy0, Hkin: H });
  close(h.stress[99 + 400], h.stress[99 + 800], 1e-9, 'peak stress repeats');
  // loop: parallelogram with plastic strain range Δεp = 2(amp − σmax/E), height 2σy0 → energy per cycle = 2σy0·Δεp
  const sMax = h.stress[99 + 800];
  const loop = 2 * sy0 * 2 * (amp - sMax / E);
  const twoCycles = stressStrainHistory1d({ strains: cyc.slice(0, 800), E, sy0, Hkin: H }).dissipation;
  close(h.dissipation - twoCycles, loop, 2e-3, 'dissipation of the third cycle');
});

test('1-D plasticity with Voce hardening: stress on the saturating yield surface, tangent against finite differences', () => {
  const voce = { sInf: 450e6, delta: 30 };
  const r = plasticity1d({ strain: 0.05, E, sy0, voce });
  close(r.stress, sy0 + 150e6 * (1 - Math.exp(-30 * r.state.alpha)), 1e-12);
  close(r.state.alpha, 0.05 - r.stress / E, 1e-12);
  const d = 1e-8;
  const fd = (plasticity1d({ strain: 0.05 + d, E, sy0, voce }).stress - plasticity1d({ strain: 0.05 - d, E, sy0, voce }).stress) / (2 * d);
  close(r.tangent, fd, 1e-5);
  // saturation
  close(plasticity1d({ strain: 1, E, sy0, voce }).stress, 450e6, 1e-9);
  close(hardeningStress({ alpha: 0, sy0, voce }).slope, 150e6 * 30, 1e-12);
});

test('3-D radial return: stress returned to the yield surface, deviatoric plastic flow, pressure untouched', () => {
  const nu = 0.3; const H = 2e9;
  const strain = [0.004, -0.001, 0.0005, 0.003, -0.002, 0.001];
  const r = plasticity3d({ strain, E, nu, sy0, H });
  assert.ok(r.yielded);
  // consistency condition f = 0 at the end of the step
  close(r.vonMises, hardeningStress({ alpha: r.state.alpha, sy0, H }).sy, 1e-12);
  // plastic strain is deviatoric and the mean stress is elastic: p = K tr ε
  near(r.state.epsP[0] + r.state.epsP[1] + r.state.epsP[2], 0, 1e-18);
  const K = E / (3 * (1 - 2 * nu));
  close((r.stress[0] + r.stress[1] + r.stress[2]) / 3, K * (strain[0] + strain[1] + strain[2]), 1e-12);
  // radial: the final deviatoric stress is parallel to the trial deviatoric stress
  const G = E / 2.6; const tr = strain[0] + strain[1] + strain[2];
  const sTr = strain.map((e, i) => 2 * G * (i < 3 ? e - tr / 3 : e / 2));
  const m = (r.stress[0] + r.stress[1] + r.stress[2]) / 3;
  const s = r.stress.map((v, i) => v - (i < 3 ? m : 0));
  const ratio = s[0] / sTr[0];
  for (let i = 1; i < 6; i++) close(s[i] / sTr[i], ratio, 1e-10);
  assert.ok(ratio > 0 && ratio < 1);
  // equivalent plastic strain: α = √(2/3)‖e_p‖
  let n2 = 0;
  for (let i = 0; i < 6; i++) n2 += (i < 3 ? 1 : 2) * r.state.epsP[i] ** 2;
  close(r.state.alpha, Math.sqrt(2 / 3 * n2), 1e-12);
  // elastic step: Hooke's law, σxy = G γxy
  const e = plasticity3d({ strain: [1e-4, 0, 0, 2e-4, 0, 0], E, nu, sy0 });
  assert.equal(e.yielded, false); close(e.stress[3], G * 2e-4, 1e-12); close(e.stress[0], (K + 4 * G / 3) * 1e-4, 1e-12);
  close(vonMisesStress([100, 0, 0, 0, 0, 0]), 100, 1e-12); close(vonMisesStress([0, 0, 0, 50, 0, 0]), 50 * Math.sqrt(3), 1e-12);
});

test('3-D radial return: consistent tangent equals the finite-difference derivative of the stress update (linear, Voce and kinematic hardening)', () => {
  const nu = 0.3;
  const cases = [{ H: 3e9 }, { voce: { sInf: 420e6, delta: 25 } }, { H: 1e9, Hkin: 4e9 }, {}];
  const base = [0.002, -0.0005, 0.0003, 0.0015, -0.001, 0.0008];
  const strain = [0.005, -0.0012, 0.0002, 0.004, -0.0025, 0.0018];
  for (const c of cases) {
    // a previous plastic step gives a non-trivial state (back stress, plastic strain)
    const st = plasticity3d({ strain: base, E, nu, sy0, ...c }).state;
    const r = plasticity3d({ strain, E, nu, sy0, ...c, state: st });
    assert.ok(r.yielded);
    const d = 1e-9;
    for (let j = 0; j < 6; j++) {
      const ep = strain.slice(); const em = strain.slice(); ep[j] += d; em[j] -= d;
      const sp = plasticity3d({ strain: ep, E, nu, sy0, ...c, state: st, tangent: false }).stress;
      const sm = plasticity3d({ strain: em, E, nu, sy0, ...c, state: st, tangent: false }).stress;
      for (let i = 0; i < 6; i++) near(r.tangent[i][j], (sp[i] - sm[i]) / (2 * d), 2e-5 * E, `C[${i}][${j}]`);
    }
    // major symmetry of the tangent
    for (let i = 0; i < 6; i++) for (let j = 0; j < i; j++) near(r.tangent[i][j], r.tangent[j][i], 1e-9 * E);
  }
});

test('3-D model under uniaxial stress reproduces the 1-D model (monotonic and reversed, isotropic + kinematic); plastic Poisson ratio → ½', () => {
  const p = { E, nu: 0.3, sy0, H: 1.5e9, Hkin: 2.5e9 };
  const strains = [];
  for (let i = 1; i <= 50; i++) strains.push(0.012 * i / 50);
  for (let i = 1; i <= 100; i++) strains.push(0.012 - 0.02 * i / 100);
  const a = uniaxialResponse3d({ ...p, strains });
  const b = stressStrainHistory1d({ ...p, strains });
  for (const i of [5, 20, 49, 60, 100, 149]) close(a.stress[i], b.stress[i], 1e-9, `step ${i}`);
  close(a.alpha[149], b.state.alpha, 1e-8);
  // lateral strain: −ν σ/E − ½ εp
  close(a.lateralStrain[49], -0.3 * a.stress[49] / E - 0.5 * b.epsP[49], 1e-9);
  // Voce hardening as well
  const v = { E, nu: 0.3, sy0, voce: { sInf: 500e6, delta: 40 } };
  close(uniaxialResponse3d({ ...v, strains: [0.002, 0.01, 0.03] }).stress[2], plasticity1d({ ...v, strain: 0.03 }).stress, 1e-9);
  // pure shear: yield at τ = σy/√3
  const G = E / 2.6;
  const sh = plasticity3d({ strain: [0, 0, 0, 0.01, 0, 0], E, nu: 0.3, sy0 });
  close(sh.stress[3], sy0 / Math.sqrt(3), 1e-12);
  assert.equal(plasticity3d({ strain: [0, 0, 0, 0.999 * sy0 / Math.sqrt(3) / G, 0, 0], E, nu: 0.3, sy0 }).yielded, false);
});

test('Ramberg–Osgood: 0.2 % offset at the proof stress, inverse, exponent from two points', () => {
  const o = { E: 200e9, sy: 400e6, n: 10 };
  const y = rambergOsgoodStrain({ ...o, stress: 400e6 });
  close(y.strain, 0.002 + 0.002, 1e-12); close(y.plasticStrain, 0.002, 1e-12);
  // σ = 300 MPa: 0.0015 + 0.002·0.75^10 = 0.0015 + 1.126e-4
  close(rambergOsgoodStrain({ ...o, stress: 300e6 }).strain, 0.0015 + 0.002 * 0.75 ** 10, 1e-12);
  close(rambergOsgoodStrain({ ...o, stress: -300e6 }).strain, -(0.0015 + 0.002 * 0.75 ** 10), 1e-12);
  for (const s of [50e6, 300e6, 400e6, 520e6]) {
    const e = rambergOsgoodStrain({ ...o, stress: s });
    const inv = rambergOsgoodStress({ ...o, strain: e.strain });
    close(inv.stress, s, 1e-10); close(inv.tangent, e.tangent, 1e-8);
    const d = 1e3;
    close(e.tangent, 2 * d / (rambergOsgoodStrain({ ...o, stress: s + d }).strain - rambergOsgoodStrain({ ...o, stress: s - d }).strain), 1e-6);
  }
  close(rambergOsgoodStress({ ...o, strain: -0.004 }).stress, -400e6, 1e-10);
  // n from σy = 400 MPa at 0.2 % and σu = 500 MPa at 8 % plastic strain: ln(40)/ln(1.25) = 16.53
  const n = rambergOsgoodExponent({ sy: 400e6, su: 500e6, plasticStrainU: 0.08 });
  close(n, Math.log(40) / Math.log(1.25), 1e-12); close(n, 16.53, 1e-3);
  close(rambergOsgoodStrain({ E: 200e9, sy: 400e6, n, stress: 500e6 }).plasticStrain, 0.08, 1e-10);
});

test('Johnson–Cook and Cowper–Symonds: factors by hand', () => {
  // generic constants (supplied by the caller): A = 350 MPa, B = 275 MPa, n = 0.36, C = 0.022, m = 1.0, Tm = 1811 K
  const jc = { A: 350e6, B: 275e6, n: 0.36, C: 0.022, m: 1.0, Tmelt: 1811 };
  const q = johnsonCook({ ...jc, epsP: 0.1 });
  close(q.stress, 350e6 + 275e6 * 0.1 ** 0.36, 1e-12);
  const r = johnsonCook({ ...jc, epsP: 0.1, strainRate: 1000, T: 600 });
  close(r.rateFactor, 1 + 0.022 * Math.log(1000), 1e-12);
  close(r.thermalFactor, 1 - (600 - 293.15) / (1811 - 293.15), 1e-12);
  close(r.stress, q.stress * r.rateFactor * r.thermalFactor, 1e-12);
  assert.equal(johnsonCook({ ...jc, epsP: 0.1, strainRate: 0.01 }).rateFactor, 1);
  assert.equal(johnsonCook({ ...jc, epsP: 0.1, T: 2000 }).stress, 0);
  assert.equal(johnsonCook({ ...jc, epsP: 0 }).stress, 350e6);
  close(cowperSymonds({ strainRate: 40.4, D: 40.4, q: 5 }), 2, 1e-12);
  close(cowperSymonds({ strainRate: 100, D: 6500, q: 4 }), 1 + (100 / 6500) ** 0.25, 1e-12);
});

test('true ↔ engineering: σt = σe(1 + εe), εt = ln(1 + εe); Considère point of a power-law material', () => {
  const t = trueFromEngineering({ stress: 500e6, strain: 0.2 });
  close(t.stress, 600e6, 1e-12); close(t.strain, Math.log(1.2), 1e-12);
  const e = engineeringFromTrue(t);
  close(e.stress, 500e6, 1e-12); close(e.strain, 0.2, 1e-12);
  // σt = K εt^n: engineering stress σt e^(−εt) peaks at εt = n (Considère)
  const n = 0.2; const K = 800e6;
  const eng = (et) => engineeringFromTrue({ stress: K * et ** n, strain: et }).stress;
  assert.ok(eng(0.2) > eng(0.19) && eng(0.2) > eng(0.21));
});

test('hyperelastic uniaxial: Neo-Hookean σ = μ(λ² − 1/λ), small-strain modulus 3μ; Mooney–Rivlin 6(C10 + C01); energy consistency', () => {
  const mu = 1e6;
  const r = neoHookeanUniaxial({ stretch: 2, mu });
  close(r.cauchy, mu * (4 - 0.5), 1e-12); close(r.nominal, mu * (2 - 0.25), 1e-12);
  close(neoHookeanUniaxial({ stretch: 1 + 1e-6, mu }).nominal / 1e-6, 3 * mu, 1e-5);
  close(neoHookeanUniaxial({ stretch: 1, mu }).tangent, 3 * mu, 1e-12);
  const mr = { C10: 0.3e6, C01: 0.1e6 };
  close(mooneyRivlinUniaxial({ stretch: 1, ...mr }).tangent, 6 * 0.4e6, 1e-12);
  // λ = 1.5: P = 2(1.5 − 1/2.25)(0.3 + 0.1/1.5) MPa = 2·1.05556·0.36667 = 0.77407 MPa
  close(mooneyRivlinUniaxial({ stretch: 1.5, ...mr }).nominal, 0.774074e6, 1e-5);
  // nominal stress is the derivative of the strain energy with respect to the stretch; tangent is its derivative
  const d = 1e-6;
  for (const l of [0.6, 1.3, 3]) {
    const fd = (mooneyRivlinUniaxial({ stretch: l + d, ...mr }).energy - mooneyRivlinUniaxial({ stretch: l - d, ...mr }).energy) / (2 * d);
    close(mooneyRivlinUniaxial({ stretch: l, ...mr }).nominal, fd, 1e-6);
    const fd2 = (mooneyRivlinUniaxial({ stretch: l + d, ...mr }).nominal - mooneyRivlinUniaxial({ stretch: l - d, ...mr }).nominal) / (2 * d);
    close(mooneyRivlinUniaxial({ stretch: l, ...mr }).tangent, fd2, 1e-6);
  }
  assert.ok(neoHookeanUniaxial({ stretch: 0.5, mu }).cauchy < 0);
});

test('Mander confined concrete: f\'cc = 46.95 MPa and εcc = 0.00765 for f\'co = 30 MPa, f\'l = 3 MPa; peak at εcc; unconfined limit', () => {
  // −1.254 + 2.254 √(1 + 7.94·0.1) − 0.2 = −1.454 + 2.254·1.33940 = 1.5650
  const r = manderConfinedConcrete({ strain: 0.0076503, fco: 30e6, fl: 3e6 });
  close(r.fcc, 46.95e6, 1e-3); close(r.epsCc, 0.00765, 1e-3);
  close(r.stress, r.fcc, 1e-6);
  close(r.Ec, 5000e6 * Math.sqrt(30), 1e-12); close(r.r, r.Ec / (r.Ec - r.fcc / r.epsCc), 1e-12);
  // the curve rises to the peak and falls slowly beyond it; initial slope Ec
  assert.ok(manderConfinedConcrete({ strain: 0.004, fco: 30e6, fl: 3e6 }).stress < r.fcc);
  const post = manderConfinedConcrete({ strain: 0.02, fco: 30e6, fl: 3e6 }).stress;
  assert.ok(post < r.fcc && post > 0.8 * r.fcc);
  close(manderConfinedConcrete({ strain: 1e-7, fco: 30e6, fl: 3e6 }).stress / 1e-7, r.Ec, 1e-3);
  // unconfined: f'cc = f'co at εco
  const u = manderConfinedConcrete({ strain: 0.002, fco: 30e6 });
  close(u.fcc, 30e6, 1e-12); close(u.stress, 30e6, 1e-12); close(u.epsCc, 0.002, 1e-12);
  // confinement by spirals: ds = 500 mm, pitch 75 mm, 12 mm bar (113.1 mm²), fyh = 420 MPa, ρcc = 0.02
  // ρs = 4·113.1/(500·75) = 0.012064; s' = 63 mm; ke = (1 − 63/1000)/(0.98) = 0.9561; f'l = ½·0.9561·0.012064·420 = 2.422 MPa
  const c = manderConfinement({ ds: 0.5, spacing: 0.075, Asp: 113.1e-6, fyh: 420e6, rhoCc: 0.02, type: 'spiral' });
  close(c.rhoS, 0.012064, 1e-4); close(c.ke, 0.9561, 1e-3); close(c.fl, 2.422e6, 1e-3);
  // hoops are less effective than a spiral at the same spacing: ke ratio = (1 − s'/2ds)
  close(manderConfinement({ ds: 0.5, spacing: 0.075, Asp: 113.1e-6, fyh: 420e6, rhoCc: 0.02 }).ke / c.ke, 1 - 0.063 / 1.0, 1e-3);
});
