import test from 'node:test';
import assert from 'node:assert/strict';
import * as M from '../../src/engine/electrical/machines.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const spm = { polePairs: 50, psiM: 10, Ld: 5e-3, Lq: 5e-3 };

test('PMSM torque and current references: zero-d control and MTPA', () => {
  // Te = 1.5 p ψm iq = 1.5 × 50 × 10 × 1000 = 750 kN·m; saliency adds 1.5 p (Ld − Lq) id iq
  rel(M.pmsmTorque({ ...spm, id: 0, iq: 1000 }), 750e3, 1e-12);
  rel(M.pmsmTorque({ ...spm, Ld: 4e-3, Lq: 8e-3, id: -200, iq: 1000 }), 1.5 * 50 * (10 * 1000 + (-4e-3) * (-200) * 1000), 1e-12);
  const z = M.pmsmCurrentReference({ ...spm, torque: -750e3, strategy: 'zeroD' });
  close(z.id, 0, 1e-12); rel(z.iq, -1000, 1e-12); rel(z.Is, 1000, 1e-12);
  // non-salient machine: MTPA is zero-d control
  const m0 = M.pmsmCurrentReference({ ...spm, torque: 750e3 });
  close(m0.id, 0, 1e-6); rel(m0.iq, 1000, 1e-6);
  // salient machine: the reference delivers the torque, obeys the MTPA expression and needs less current than zero-d
  const sal = { polePairs: 50, psiM: 10, Ld: 4e-3, Lq: 8e-3 };
  const r = M.pmsmCurrentReference({ ...sal, torque: 750e3 });
  rel(M.pmsmTorque({ ...sal, id: r.id, iq: r.iq }), 750e3, 1e-6);
  const dL = -4e-3;
  rel(r.id, (-10 + Math.sqrt(100 + 8 * dL * dL * r.Is * r.Is)) / (4 * dL), 1e-6);
  assert.ok(r.id < 0 && r.Is < 1000);
  rel(Math.hypot(r.id, r.iq), r.Is, 1e-9);
});

test('PMSM steady state: dq voltages, power balance and losses by hand (generating, id = 0)', () => {
  // p = 50, ψm = 10 Wb, L = 5 mH, Rs = 10 mΩ, ωm = 1.2 rad/s → ωe = 60 rad/s (9.549 Hz); iq = −1000 A
  const o = { ...spm, Rs: 0.01, speed: 1.2, id: 0, iq: -1000, kh: 100, ke: 2, frictionTorque: 1000, windage: 500 };
  const r = M.pmsmSteadyState(o);
  rel(r.frequency, 60 / (2 * Math.PI), 1e-12);
  // vd = −ωe Lq iq = 300 V; vq = Rs iq + ωe ψm = −10 + 600 = 590 V
  rel(r.vd, 300, 1e-12); rel(r.vq, 590, 1e-12); rel(r.Vs, Math.hypot(300, 590), 1e-12); rel(r.Vll, r.Vs * Math.sqrt(3 / 2), 1e-12);
  rel(r.Is, 1000, 1e-12); rel(r.Irms, 1000 / Math.SQRT2, 1e-12);
  rel(r.torque, -750e3, 1e-12); rel(r.electricalPower, 1.5 * 590 * -1000, 1e-12); // −885 kW
  rel(r.copperLoss, 1.5 * 0.01 * 1e6, 1e-12); // 15 kW
  // terminal power = Te ωm + copper: −900 + 15 = −885 kW
  rel(r.electricalPower, r.torque * 1.2 + r.copperLoss, 1e-12);
  // stator flux √(ψm² + (Lq iq)²) = √(100 + 25) → iron loss (100 × 60 + 2 × 3600) × 1.25 = 16.5 kW
  rel(r.ironLoss, (100 * 60 + 2 * 3600) * 1.25, 1e-9);
  rel(r.mechanicalLoss, 1000 * 1.2 + 500 * 1.2 ** 3, 1e-12);
  // generator: shaft input = |Te ωm| + iron + mechanical; efficiency = electrical out / shaft in
  rel(Math.abs(r.shaftPower), 900e3 + r.ironLoss + r.mechanicalLoss, 1e-9);
  rel(r.efficiency, 885e3 / Math.abs(r.shaftPower), 1e-9);
  rel(Math.abs(r.powerFactor), 885e3 / (1.5 * r.Vs * r.Is), 1e-9);
  // torque demand route gives the same point
  const t = M.pmsmSteadyState({ ...spm, Rs: 0.01, speed: 1.2, torque: -750e3 });
  rel(t.iq, -1000, 1e-6); rel(t.vq, 590, 1e-6);
  // lossless machine: efficiency 1
  rel(M.pmsmSteadyState({ ...spm, Rs: 0, speed: 1.2, id: 0, iq: -1000 }).efficiency, 1, 1e-12);
  const map = M.pmsmEfficiencyMap({ ...spm, Rs: 0.01, speeds: [0.6, 1.2], torques: [375e3, 750e3] });
  assert.equal(map.efficiency.length, 2); assert.equal(map.efficiency[0].length, 2);
  rel(map.efficiency[1][1], 885 / 900, 1e-6); rel(map.power[1][1], 885e3, 1e-6);
  assert.ok(map.efficiency[1][0] < map.efficiency[1][1], 'same copper loss on half the power at half speed');
});

test('induction machine equivalent circuit: Thévenin torque, breakdown point and the Kloss identity', () => {
  // 400 V, 50 Hz, 2 pole pairs; R1 = 0.03, X1 = 0.1, R2 = 0.03, X2 = 0.1, Xm = 4 Ω
  const o = { V: 400, f: 50, polePairs: 2, R1: 0.03, X1: 0.1, R2: 0.03, X2: 0.1, Xm: 4 };
  const ws = 2 * Math.PI * 50 / 2;
  // Thévenin: Vth = V1 jXm/(R1 + j(X1 + Xm)), Zth = jXm (R1 + jX1)/(R1 + j(X1 + Xm))
  const V1 = 400 / Math.sqrt(3); const den = 0.03 ** 2 + 4.1 ** 2;
  const Vth = V1 * 4 / Math.sqrt(den);
  const Rth = (4 * 4 * 0.03) / den; // Re{jXm (R1 + jX1)/(R1 + j(X1 + Xm))} = Xm² R1/den
  const Xth = 4 * (0.03 ** 2 + 0.1 * 4.1) / den;
  const s = 0.02;
  const T = 3 / ws * Vth ** 2 * (0.03 / s) / ((Rth + 0.03 / s) ** 2 + (Xth + 0.1) ** 2);
  const r = M.inductionMachine({ ...o, slip: s });
  rel(r.torque, T, 1e-9); rel(r.torqueKloss, T, 1e-9);
  const sb = 0.03 / Math.hypot(Rth, Xth + 0.1);
  rel(r.breakdownSlip, sb, 1e-9);
  rel(r.breakdownTorque, 3 / ws * Vth ** 2 / (2 * (Rth + Math.hypot(Rth, Xth + 0.1))), 1e-9);
  rel(M.inductionMachine({ ...o, slip: sb }).torque, r.breakdownTorque, 1e-9);
  rel(r.speed, ws * (1 - s), 1e-12);
  // air-gap power = T ωs; mechanical power = (1 − s) × that; efficiency below (1 − s) because of the stator copper loss
  assert.ok(r.efficiency < 1 - s && r.efficiency > 0.9); assert.ok(r.P > 0 && r.Q > 0 && r.powerFactor > 0.7 && r.powerFactor < 1);
  rel(r.P, Math.sqrt(3) * 400 * r.I1 * r.powerFactor, 1e-9);
  // generating at negative slip: torque and stator power reverse, reactive power is still absorbed
  const g = M.inductionMachine({ ...o, slip: -0.02 });
  assert.ok(g.torque < 0 && g.P < 0 && g.Q > 0);
  // simplified Kloss relation: T = Tb at s = sb, T ≈ 2 Tb s/sb for small slip, odd in s
  rel(M.klossTorque({ slip: 0.2, breakdownSlip: 0.2, breakdownTorque: 1000 }), 1000, 1e-12);
  rel(M.klossTorque({ slip: 0.1, breakdownSlip: 0.2, breakdownTorque: 1000 }), 800, 1e-12);
  rel(M.klossTorque({ slip: -0.1, breakdownSlip: 0.2, breakdownTorque: 1000 }), -800, 1e-12);
});

test('doubly-fed generator: mesh equations hold, lossless limit Pr = −s Ps and Pmech = (1 − s) Ps', () => {
  const o = { Vs: 690, Vr: { re: -60, im: 15 }, slip: -0.2, Rs: 1e-9, Xs: 0.05, Rr: 1e-9, Xr: 0.05, Xm: 3 };
  const r = M.dfigSteadyState(o);
  rel(r.rotorPowerRatio, 0.2, 1e-5); rel(r.Pr, 0.2 * r.Ps, 1e-5); rel(r.Pmech, 1.2 * r.Ps, 1e-5); rel(r.Ptotal, r.Ps + r.Pr, 1e-12);
  // check the stator mesh equation Vs = (Rs + jXs) Is + jXm (Is + Ir) with the returned phasors
  const V1 = 690 / Math.sqrt(3);
  const re = -0.05 * r.Is.im - 3 * (r.Is.im + r.Ir.im); const im = 0.05 * r.Is.re + 3 * (r.Is.re + r.Ir.re);
  close(re, V1, 1e-4); close(im, 0, 1e-4);
  // rotor mesh: Vr/s = (Rr/s + jXr) Ir + jXm (Is + Ir)
  close(-0.05 * r.Ir.im - 3 * (r.Is.im + r.Ir.im), -60 / -0.2, 1e-4); close(0.05 * r.Ir.re + 3 * (r.Is.re + r.Ir.re), 15 / -0.2, 1e-4);
  // with resistances the copper loss closes the balance: Pmech = Ps + Pr − copper
  const l = M.dfigSteadyState({ ...o, Rs: 0.005, Rr: 0.005 });
  rel(l.Pmech, l.Ps + l.Pr - l.copperLoss, 1e-12);
  rel(l.copperLoss, 3 * 0.005 * (l.Is.re ** 2 + l.Is.im ** 2 + l.Ir.re ** 2 + l.Ir.im ** 2), 1e-9);
});

test('synchronous machine power–angle relations: round rotor and salient pole', () => {
  // E = 12 kV, V = 11 kV line-to-line, Xd = 10 Ω, δ = 30°: P = E V sin δ/Xd (line values) = 6.6 MW; Q = V (E cos δ − V)/Xd
  const r = M.synchronousMachinePower({ E: 12e3, V: 11e3, delta: Math.PI / 6, Xd: 10 });
  rel(r.P, 12e3 * 11e3 * 0.5 / 10, 1e-12); rel(r.Q, 11e3 * (12e3 * Math.cos(Math.PI / 6) - 11e3) / 10, 1e-12); rel(r.Pmax, 13.2e6, 1e-12);
  // saliency adds the reluctance term V²/2 (1/Xq − 1/Xd) sin 2δ
  const s = M.synchronousMachinePower({ E: 12e3, V: 11e3, delta: Math.PI / 6, Xd: 10, Xq: 6 });
  rel(s.P - r.P, 11e3 ** 2 / 2 * (1 / 6 - 1 / 10) * Math.sin(Math.PI / 3), 1e-9);
  close(M.synchronousMachinePower({ E: 12e3, V: 11e3, delta: 0, Xd: 10 }).P, 0, 1e-6);
});

test('converter losses by the averaged device formulas; DC-link sizing; PQ capability', () => {
  // Î = 1000 A (Irms = 707.1), m = 0.9, cosφ = 1, defaults: per IGBT 1.0 × 1000 × (0.15915 + 0.1125) + 1.2e-3 × 1e6 × (0.125 + 0.095493)
  const r = M.converterLoss({ Irms: 1000 / Math.SQRT2, Vdc: 900 });
  const igbt = 1000 * (1 / (2 * Math.PI) + 0.9 / 8) + 1.2e-3 * 1e6 * (1 / 8 + 0.9 / (3 * Math.PI));
  const diode = 1000 * (1 / (2 * Math.PI) - 0.9 / 8) + 0.8e-3 * 1e6 * (1 / 8 - 0.9 / (3 * Math.PI));
  rel(r.conduction, 6 * (igbt + diode), 1e-9);
  // switching: 6 × 2500 × (0.35 + 0.40 + 0.25)/π × 1 × 1 = 4775 W
  rel(r.switching, 6 * 2500 * 1.0 / Math.PI, 1e-9); rel(r.total, r.conduction + r.switching, 1e-12);
  rel(r.a1 * 1000 + r.a2 * 1e6, r.total, 1e-9);
  // half the DC voltage halves the switching loss; rectifier operation shifts conduction from the IGBTs to the diodes
  rel(M.converterLoss({ Irms: 1000 / Math.SQRT2, Vdc: 450 }).switching, 0.5 * r.switching, 1e-12);
  const rect = M.converterLoss({ Irms: 1000 / Math.SQRT2, Vdc: 900, cosPhi: -1 });
  rel(rect.switching, r.switching, 1e-12); assert.ok(rect.conduction !== r.conduction);
  // two modules in parallel share the current: conduction I² terms halve
  assert.ok(M.converterLoss({ Irms: 1000 / Math.SQRT2, Vdc: 900, parallel: 2 }).conduction < r.conduction);
  // DC link: 5 MVA, 1100 V, τ = 10 ms → C = 2 × 5e6 × 0.01/1100² = 82.6 mF, energy 50 kJ
  const c = M.dcLinkCapacitance({ S: 5e6, Vdc: 1100 });
  rel(c.Cenergy, 2 * 5e6 * 0.01 / 1100 ** 2, 1e-12); rel(c.storedEnergy, 5e4, 1e-9); rel(c.C, c.Cenergy, 1e-12);
  // ripple rule: 10 % of S at 100 Hz with 5 % ripple: C = 5e5/(2π 100 × 1100 × 55) = 13.2 mF
  const cr = M.dcLinkCapacitance({ S: 5e6, Vdc: 1100, tau: 1e-3, rippleRatio: 0.1, rippleFrequency: 100 });
  rel(cr.Cripple, 5e5 / (2 * Math.PI * 100 * 1100 * 55), 1e-12); rel(cr.C, Math.max(cr.Cripple, cr.Cenergy), 1e-12);
  // capability: 690 V, 5000 A → S = 5.976 MVA; at P = 4 MW Q = ±4.44 Mvar by the current limit
  const k = M.converterCapability({ V: 690, Imax: 5000, P: 4e6 });
  rel(k.Smax, Math.sqrt(3) * 690 * 5000, 1e-12); rel(k.Qmax, Math.sqrt(k.Smax ** 2 - 16e12), 1e-12); rel(k.Qmin, -k.Qmax, 1e-12);
  assert.equal(k.feasible, true); assert.match(k.limitedBy, /current/i);
  // voltage limit: X = 0.02 Ω, Vc,max = 720 V: Q ≤ −V²/X + √((V Vc/X)² − P²) = 0.92 Mvar → governs export
  const kv = M.converterCapability({ V: 690, Imax: 5000, P: 4e6, X: 0.02, VcMax: 720 });
  rel(kv.Qmax, -(690 ** 2) / 0.02 + Math.sqrt((690 * 720 / 0.02) ** 2 - 16e12), 1e-9); assert.match(kv.limitedBy, /voltage/i);
  rel(kv.Qmin, -k.Qmax, 1e-9);
  assert.equal(M.converterCapability({ V: 690, Imax: 5000, P: 7e6 }).feasible, false);
});

test('gearbox efficiency, optimal-torque gain and the drive-train efficiency chain', () => {
  // three stages 1 %, 1 %, 1.5 %, no-load 0.5 % of 5 MW; at 2.5 MW: 2.5 × 0.9654 − 0.025 = 2.3885 MW → 95.54 %
  const g = M.gearboxEfficiency({ stageLosses: [0.01, 0.01, 0.015], noLoadLossFraction: 0.005, power: 2.5e6, ratedPower: 5e6 });
  rel(g.outputPower, 2.5e6 * 0.99 * 0.99 * 0.985 - 25e3, 1e-12); rel(g.efficiency, g.outputPower / 2.5e6, 1e-12); rel(g.loss, 2.5e6 - g.outputPower, 1e-12);
  rel(M.gearboxEfficiency({ power: 1e6, ratedPower: 5e6 }).efficiency, 1, 1e-15);
  // k = ½ ρ π R⁵ Cp/λ³: ρ = 1.225, R = 63, Cp = 0.48, λ = 7.5 → 2.1717e6 N·m·s²
  const k = M.optimalTorqueGain({ rho: 1.225, R: 63, CpMax: 0.48, lambdaOpt: 7.5, gearRatio: 97, omega: 1.0 });
  rel(k.k, 0.5 * 1.225 * Math.PI * 63 ** 5 * 0.48 / 7.5 ** 3, 1e-12); rel(k.kGenerator, k.k / 97 ** 3, 1e-12);
  rel(k.torque, k.k, 1e-12); rel(k.power, k.k, 1e-12);
  // the law tracks Cp,max: at ω = λ V/R the power is ½ ρ π R² Cp V³
  const V = 8; const om = 7.5 * V / 63;
  rel(M.optimalTorqueGain({ rho: 1.225, R: 63, CpMax: 0.48, lambdaOpt: 7.5, omega: om }).power, 0.5 * 1.225 * Math.PI * 63 ** 2 * 0.48 * V ** 3, 1e-12);
  assert.equal(M.optimalTorqueGain({ rho: 1.225, R: 63, CpMax: 0.48, lambdaOpt: 7.5 }).torque, null);
  // direct-drive chain without converter or transformer: grid power = generator output of the dq model
  const s = M.operatingSchedule({ rotorSpeed: [1.2], rotorTorque: [750e3], generator: { ...spm, Rs: 0.01, strategy: 'zeroD' }, ratedPower: 900e3 });
  rel(s.shaftPower[0], 900e3, 1e-12); rel(s.generatorPower[0], 885e3, 1e-6); rel(s.gridPower[0], 885e3, 1e-6);
  rel(s.efficiency[0], 885 / 900, 1e-6); rel(s.losses.generator[0], 15e3, 1e-6); close(s.losses.gearbox[0], 0, 1e-9);
  // with a transformer (Pk 9 kW at Sr 1 MVA, P0 1 kW) the output drops by the load and no-load losses
  const t = M.operatingSchedule({ rotorSpeed: [1.2], rotorTorque: [750e3], generator: { ...spm, Rs: 0.01, strategy: 'zeroD' }, transformer: { Sr: 1e6, Pk: 9e3, P0: 1e3 }, ratedPower: 900e3 });
  assert.ok(t.gridPower[0] < s.gridPower[0] && t.losses.transformer[0] > 1e3 && t.losses.transformer[0] < 10e3);
  rel(t.gridPower[0], s.gridPower[0] - t.losses.transformer[0], 1e-9);
});

test('MPPT on a parabolic power curve, linear generator on a resistive load and the hydraulic PTO chain', () => {
  // I(V) = 10 (1 − (V/40)²)... choose P = V I with maximum at dP/dV = 0: I = 10 − V/4 → P = 10V − V²/4, MPP at 20 V, 100 W
  const cur = (V) => Math.max(10 - V / 4, 0);
  const po = M.pvMppt({ current: cur, V0: 5, dV: 0.5, steps: 200 });
  close(po.Vfinal, 20, 1.01); rel(po.Pmean, 100, 5e-3); assert.equal(po.V.length, po.P.length);
  const ic = M.pvMppt({ current: cur, V0: 5, dV: 0.5, steps: 200, method: 'incrementalConductance' });
  close(ic.Vfinal, 20, 0.51); rel(ic.Pfinal, 100, 2e-3);
  assert.ok(ic.Pmean >= po.Pmean - 1e-9, 'incremental conductance settles without the perturb-and-observe oscillation');
  // voltage window respected
  assert.ok(M.pvMppt({ current: cur, V0: 5, dV: 0.5, steps: 100, Vmax: 15 }).Vfinal <= 15 + 1e-12);
  // ke = 400 V·s/m, v = 1.5 m/s, Rg = 2 Ω, RL = 8 Ω: e = 600 V, i = 60 A, F = 1.5 × 400 × 60 = 36 kN, P = 1.5 × 3600 × 8 = 43.2 kW
  const l = M.linearGenerator({ ke: 400, velocity: 1.5, Rg: 2, RL: 8 });
  rel(l.emf, 600, 1e-12); rel(l.current, 60, 1e-12); rel(l.force, 36e3, 1e-12); rel(l.damping, 24e3, 1e-12);
  rel(l.outputPower, 43.2e3, 1e-12); rel(l.copperLoss, 10.8e3, 1e-12); rel(l.efficiency, 0.8, 1e-12); rel(l.meanOutputSinusoidal, 21.6e3, 1e-12);
  // mechanical power F v = output + copper loss
  rel(l.force * 1.5, l.outputPower + l.copperLoss, 1e-12);
  // hydraulic PTO at rated power: product of the rated efficiencies 0.95 × 0.92 × 0.95 × 0.97 = 0.8054
  const h = M.hydraulicPto({ power: 1e6, ratedPower: 1e6 });
  rel(h.efficiency, 0.95 * 0.92 * 0.95 * 0.97, 1e-9); rel(h.outputPower, 1e6 * h.efficiency, 1e-12); rel(h.stages.motor, 0.92, 1e-9);
  // part load: η(x) = x/(x + l0 + l1 x²) with l0 = l1 = ½ (1/η − 1): motor at x = 0.25 → 0.25/(0.25 + 0.043478 + 0.002717) = 0.8440
  const p = M.hydraulicPto({ power: 0.25e6, ratedPower: 1e6 });
  const l0 = 0.5 * (1 / 0.92 - 1);
  rel(p.stages.motor, 0.25 / (0.25 + l0 + l0 * 0.0625), 1e-9); assert.ok(p.efficiency < h.efficiency);
});
