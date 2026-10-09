import test from 'node:test';
import assert from 'node:assert/strict';
import {
  geometricMeanDistance, inductance, coaxialCapacitance, chargingCurrent, criticalLength, surgeImpedance, earthReturn,
  sequenceImpedance, propagation, abcdDistributed, abcdNominalPi, abcdShunt, abcdSeries, abcdCascade, solveTwoPort,
  lineProfile, ferrantiRise, compensatedLink, transferCapability, optimalCompensationSplit, resonanceFrequency,
  hvdcLink, acDcBreakEven, perUnitBase, perUnitChangeBase, threePhasePower, transformer, lossLoadFactor,
  capitalRecoveryFactor, economicConductorSize,
} from '../../src/engine/electrical/circuit.js';
import { cx, cabs, csub, cmul } from '../../src/engine/electrical/complex.js';

const rel = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${a} vs ${b} (tol ${tol})`);
const line = { r: 3e-5, l: 4e-7, c: 2.2e-10, f: 50 };

test('inductance by hand, trefoil s = 100 mm, dc = 40 mm: L = 2e-7·ln(0.1/(0.7788·0.02)) = 0.37189 µH/m; flat GMD = ∛2·s', () => {
  rel(inductance({ dc: 0.04, gmd: geometricMeanDistance({ formation: 'trefoil', s: 0.1 }) }), 3.7189e-7, 1e-4);
  rel(geometricMeanDistance({ formation: 'flat', s: 0.1 }), 0.1259921, 1e-6);
  assert.throws(() => inductance({ dc: 0.04, gmd: 0.01 }), RangeError);
});

test('capacitance and charging current by hand for a 132 kV XLPE cable (εr 2.5, 40/72 mm)', () => {
  // C = 2π·8.8541878e-12·2.5/ln(1.8) = 1.390812e-10/0.587787 = 2.36618e-10 F/m (0.2366 µF/km)
  const C = coaxialCapacitance({ epsR: 2.5, Di: 0.072, dc: 0.040 });
  rel(C, 2.36618e-10, 1e-5);
  // Ic = ω·C·U0 = 314.159·2.36618e-10·76210.2 = 5.6652e-3 A/m = 5.665 A/km; Qc = ω·C·U² = 1.2952 var/m
  const ch = chargingCurrent({ C, U: 132e3, f: 50, length: 1000 });
  rel(ch.Ic, 5.6652, 1e-4);
  rel(ch.Qc, 314.159265 * 2.36618e-10 * 132e3 ** 2 * 1000, 1e-5);
  rel(ch.Qc, Math.sqrt(3) * 132e3 * ch.Ic, 1e-12);
});

test('critical length by hand: 800 A rating / 5.6652 A/km = 141.2 km', () => {
  const C = coaxialCapacitance({ epsR: 2.5, Di: 0.072, dc: 0.040 });
  rel(criticalLength({ C, U: 132e3, f: 50, I: 800 }), 800 / 5.6652e-3, 1e-4);
  rel(criticalLength({ C, U: 132e3, f: 50, I: 800 }), 141.2e3, 1e-3);
});

test('surge impedance and loading: Zc = sqrt(L/C), SIL = U²/Zc, v = 1/sqrt(LC)', () => {
  const s = surgeImpedance({ L: 4e-7, C: 2.5e-10, U: 132e3, f: 50 });
  rel(s.Zc, 40, 1e-12); rel(s.SIL, 132e3 ** 2 / 40, 1e-12); rel(s.velocity, 1e8, 1e-12); rel(s.wavelength, 2e6, 1e-12);
});

test('Carson earth-return terms by hand: Re = π²·50·1e-7 = 4.9348e-5 Ω/m, De = 658.87·sqrt(100/50) = 931.8 m', () => {
  const e = earthReturn({ f: 50, rhoEarth: 100 });
  rel(e.Re, 4.9348e-5, 1e-5); rel(e.De, 931.78, 1e-4);
  // Re equals ω·μ0/8
  rel(e.Re, (2 * Math.PI * 50 * 4e-7 * Math.PI) / 8, 1e-12);
});

test('sequence impedances: open sheaths give Z1 = Rc + jωL; bonded sheaths add Xm²/(Rs + jXm); zero-sequence limits', () => {
  const p = { Rc: 3e-5, Rs: 3e-4, dc: 0.038, ds: 0.08, gmd: 0.1, f: 50, rhoEarth: 0.25 };
  const open = sequenceImpedance({ ...p, sheathsBonded: false });
  rel(open.Z1.re, 3e-5, 1e-12);
  rel(open.Z1.im, 2 * Math.PI * 50 * inductance({ dc: 0.038, gmd: 0.1 }), 1e-12);
  const bonded = sequenceImpedance(p);
  const xm = 2 * Math.PI * 50 * 2e-7 * Math.log(0.1 / 0.04);
  rel(bonded.Z1.re, 3e-5 + (xm * xm * 3e-4) / (9e-8 + xm * xm), 1e-9);
  rel(bonded.Z1.im, open.Z1.im - (xm ** 3) / (9e-8 + xm * xm), 1e-9);
  // bonded sheaths raise the resistance and lower the reactance
  assert.ok(bonded.Z1.re > open.Z1.re && bonded.Z1.im < open.Z1.im);
  // earth-only return: Z0 = Rc + 3Re + j·3·ωμ0/(2π)·ln(De/GMR3c)
  const earth = sequenceImpedance({ ...p, zeroSequenceReturn: 'earth' });
  rel(earth.Z0.re, 3e-5 + 3 * Math.PI ** 2 * 50 * 1e-7, 1e-12);
  // sheath-only return: resistance = Rc + Rs (earth terms cancel), reactance from ln(GMR3s/GMR3c)
  const sh = sequenceImpedance({ ...p, zeroSequenceReturn: 'sheath' });
  rel(sh.Z0.re, 3e-5 + 3e-4, 1e-9);
  const gmrc = Math.exp(-0.25) * 0.019;
  rel(sh.Z0.im, 3 * 2 * Math.PI * 50 * 2e-7 * Math.log(Math.cbrt(0.04 * 0.01) / Math.cbrt(gmrc * 0.01)), 1e-9);
  // parallel return lies between the two in resistance terms
  assert.ok(bonded.Z0.re < sh.Z0.re);
});

test('distributed line tends to the nominal-π model for short lengths: B differs by (γl)²/6, i.e. second order in length', () => {
  const diffB = (length) => {
    const d = abcdDistributed({ ...line, length }); const n = abcdNominalPi({ ...line, length });
    return cabs(csub(d.B, n.B)) / cabs(d.B);
  };
  const { gamma } = propagation(line);
  // B_exact = Z·sinh(γl)/(γl) = Z·(1 + (γl)²/6 + …)
  rel(diffB(1000), (cabs(gamma) * 1000) ** 2 / 6, 1e-3);
  assert.ok(diffB(1000) < 2e-6);
  rel(diffB(2000) / diffB(1000), 4, 1e-3);
  const d = abcdDistributed({ ...line, length: 1000 }); const n = abcdNominalPi({ ...line, length: 1000 });
  for (const k of ['A', 'C', 'D']) assert.ok(cabs(csub(d[k], n[k])) <= 1e-5 * cabs(d[k]), k);
  // at 150 km the two differ visibly but still by less than 1 % in B
  const diff = diffB(150e3);
  assert.ok(diff > 1e-3 && diff < 0.01, `relative difference ${diff}`);
});

test('reciprocity: AD − BC = 1 for the distributed line and for any cascade', () => {
  const d = abcdDistributed({ ...line, length: 80e3 });
  const det = csub(cmul(d.A, d.D), cmul(d.B, d.C));
  rel(det.re, 1, 1e-12); assert.ok(Math.abs(det.im) < 1e-12);
  const c = abcdCascade([abcdShunt(cx(0, -1e-3)), d, abcdSeries(cx(0.5, 3)), abcdShunt(cx(0, -2e-3))]);
  const det2 = csub(cmul(c.A, c.D), cmul(c.B, c.C));
  rel(det2.re, 1, 1e-12); assert.ok(Math.abs(det2.im) < 1e-12);
});

test('no-load Ferranti rise: Vr/Vs = 1/|cosh(γl)|; lossless line 1/cos(ω·l·sqrt(LC))', () => {
  const lossless = { r: 0, l: 4e-7, c: 2.2e-10, f: 50 };
  const f = ferrantiRise({ ...lossless, length: 100e3 });
  rel(f.ratio, 1 / Math.cos(2 * Math.PI * 50 * 100e3 * Math.sqrt(4e-7 * 2.2e-10)), 1e-12);
  const { gamma } = propagation(line);
  const gl = { re: gamma.re * 100e3, im: gamma.im * 100e3 };
  const coshAbs = Math.hypot(Math.cosh(gl.re) * Math.cos(gl.im), Math.sinh(gl.re) * Math.sin(gl.im));
  rel(ferrantiRise({ ...line, length: 100e3 }).ratio, 1 / coshAbs, 1e-12);
  // the same from the two-port with an open receiving end
  const op = solveTwoPort({ abcd: abcdDistributed({ ...line, length: 100e3 }), known: 'receiving', U: 132e3, P: 0, Q: 0 });
  rel(op.Ur / op.Us, 1 / coshAbs, 1e-12);
  assert.ok(f.rise > 0);
});

test('series impedance by hand: 30 MW at 33 kV unity power factor through 1 + j2 Ω → Us = 33.958 kV, losses 826.4 kW', () => {
  // I = 30e6/(√3·33e3) = 524.864 A; Vs = 19052.56 + (1 + j2)·524.864 = 19577.42 + j1049.73 → |Vs| = 19605.5 V
  const op = solveTwoPort({ abcd: abcdSeries(cx(1, 2)), known: 'receiving', U: 33e3, P: 30e6, Q: 0 });
  rel(cabs(op.Ir), 524.864, 1e-5);
  rel(op.Us, 19605.5 * Math.sqrt(3), 1e-5);
  rel(op.losses, 3 * 524.864 ** 2 * 1, 1e-5);
  rel(op.Qs, 3 * 524.864 ** 2 * 2, 1e-5);
  rel(op.efficiency, 30e6 / (30e6 + 826446), 1e-6);
  // the sending-end formulation inverts the receiving-end one
  const back = solveTwoPort({ abcd: abcdSeries(cx(1, 2)), known: 'sending', U: op.Us, P: op.Ps, Q: op.Qs });
  rel(back.Ur, 33e3, 1e-10); rel(back.Pr, 30e6, 1e-10);
});

test('loaded cable: losses = Ps − Pr equal the integral of r·|I(x)|² along the exact current profile; current peaks at an end', () => {
  const length = 90e3;
  const abcd = abcdDistributed({ ...line, length });
  const op = solveTwoPort({ abcd, known: 'receiving', U: 132e3, P: 150e6, Q: 0 });
  const prof = lineProfile({ ...line, length, Vr: op.Vr, Ir: op.Ir, n: 2001 });
  let integral = 0;
  for (let i = 1; i < prof.x.length; i++) integral += 0.5 * (prof.I[i] ** 2 + prof.I[i - 1] ** 2) * (prof.x[i] - prof.x[i - 1]);
  rel(op.losses, 3 * line.r * integral, 1e-5);
  assert.ok(prof.xImax === 0 || prof.xImax === length);
  rel(prof.I[prof.I.length - 1], cabs(op.Is), 1e-10);
  rel(prof.U[prof.U.length - 1], op.Us, 1e-10);
});

test('transfer capability: falls to zero at the critical length when uncompensated; 50/50 reactors restore it; lumped estimate at short length', () => {
  const lossless = { r: 1e-9, l: 4e-7, c: 2.2e-10, f: 50 };
  const U = 132e3; const Irated = 800;
  const lc = criticalLength({ C: 2.2e-10, U, f: 50, I: Irated });
  const cap = transferCapability({ ...lossless, lengths: [20e3, 0.9 * lc, 1.2 * lc], U, Irated });
  // short cable: P ≈ √3·U·sqrt(I² − Ic²) with the whole charging current leaving at one end
  const Ic = chargingCurrent({ C: 2.2e-10, U, f: 50, length: 20e3 }).Ic;
  rel(cap.Psend[0], Math.sqrt(3) * U * Math.sqrt(Irated ** 2 - Ic ** 2), 5e-3);
  assert.ok(cap.Psend[1] > 0 && cap.Psend[1] < 0.6 * cap.Psend[0]);
  assert.equal(cap.Psend[2], 0);
  // 50/50 reactors at 1.2 critical lengths: each end carries 0.6·I of charging current → P ≈ sqrt(1 − 0.36) = 0.8 pu
  const comp = transferCapability({ ...lossless, lengths: [1.2 * lc], U, Irated, compensation: { sending: 0.5, receiving: 0.5 } });
  rel(comp.Psend[0], 0.8 * Math.sqrt(3) * U * Irated, 0.02);
  // mid-point reactor on top improves further
  const mid = transferCapability({ ...lossless, lengths: [1.2 * lc], U, Irated, compensation: { sending: 0.25, mid: 0.5, receiving: 0.25 } });
  assert.ok(mid.Psend[0] > comp.Psend[0]);
});

test('compensated link bookkeeping: reactor ratings are fractions of ω·C·U²·l and the optimal split of full compensation is near one half', () => {
  const lossless = { r: 1e-9, l: 4e-7, c: 2.2e-10, f: 50 };
  const link = compensatedLink({ ...lossless, length: 100e3, U: 132e3, P: 100e6, compensation: { sending: 0.3, receiving: 0.7 } });
  const Qch = 2 * Math.PI * 50 * 2.2e-10 * 132e3 ** 2 * 100e3;
  rel(link.Qcharging, Qch, 1e-12); rel(link.reactors.sending, 0.3 * Qch, 1e-12); rel(link.reactors.receiving, 0.7 * Qch, 1e-12);
  rel(link.Pr, 100e6, 1e-4);
  const opt = optimalCompensationSplit({ ...lossless, length: 150e3, U: 132e3, Irated: 800, steps: 20 });
  assert.ok(Math.abs(opt.sendingShare - 0.5) <= 0.1, `share ${opt.sendingShare}`);
  assert.ok(opt.Preceive >= opt.table[0].Preceive && opt.Preceive >= opt.table[20].Preceive);
});

test('resonance: distributed value tends to 1/(2π·sqrt(Lgrid·C·l)) for a large source inductance and to the quarter wave 1/(4·l·sqrt(LC)) for a stiff source', () => {
  const p = { L: 4e-7, C: 2.2e-10, length: 30e3 };
  const weak = resonanceFrequency({ ...p, Lgrid: 0.5 });
  rel(weak.fLumped, 1 / (2 * Math.PI * Math.sqrt(0.5 * 2.2e-10 * 30e3)), 1e-12);
  rel(weak.fDistributed, weak.fLumped, 0.01);
  const stiff = resonanceFrequency({ ...p, Lgrid: 1e-7 });
  rel(stiff.fQuarterWave, 1 / (4 * 30e3 * Math.sqrt(4e-7 * 2.2e-10)), 1e-12);
  rel(stiff.fDistributed, stiff.fQuarterWave, 1e-3);
  const mid = resonanceFrequency({ ...p, Lgrid: 0.02 });
  assert.ok(mid.fDistributed < mid.fLumped && mid.fDistributed < mid.fQuarterWave);
});

test('HVDC bipole by hand: ±320 kV, 1500 A, 0.01 Ω/km, 200 km → 960 MW, loop 4 Ω, loss 9 MW, drop 6 kV', () => {
  const h = hvdcLink({ configuration: 'bipole', U: 320e3, I: 1500, r: 1e-5, length: 200e3 });
  rel(h.Pdc, 960e6, 1e-12); rel(h.Rloop, 4, 1e-12); rel(h.cableLoss, 9e6, 1e-12); rel(h.voltageDrop, 6000, 1e-12);
  rel(h.efficiency, 951 / 960, 1e-12);
  const c = hvdcLink({ configuration: 'bipole', U: 320e3, I: 1500, r: 1e-5, length: 200e3, converterLossFraction: 0.01 });
  rel(c.Pin, 960e6 / 0.99, 1e-12); rel(c.Pout, 951e6 * 0.99, 1e-12);
  rel(c.Pin - c.Pout, c.converterLoss + c.cableLoss, 1e-12);
  // monopoles
  rel(hvdcLink({ configuration: 'monopoleSeaReturn', U: 320e3, I: 1000, r: 1e-5, length: 100e3, electrodeResistance: 0.2 }).Rloop, 1.2, 1e-12);
  rel(hvdcLink({ configuration: 'monopoleMetallicReturn', U: 320e3, I: 1000, r: 1e-5, length: 100e3, rReturn: 3e-5 }).Rloop, 4, 1e-12);
  rel(hvdcLink({ configuration: 'symmetricMonopole', U: 320e3, I: 1000, r: 1e-5, length: 100e3 }).Pdc, 640e6, 1e-12);
});

test('AC/DC break-even helper: DC overtakes AC in delivered power before the AC critical length; cost-per-watt crossing interpolated', () => {
  const lengths = [20e3, 60e3, 100e3, 140e3, 180e3, 220e3];
  const out = acDcBreakEven({ lengths, ac: { ...line, U: 220e3, Irated: 1000 },
    dc: { configuration: 'symmetricMonopole', U: 160e3, I: 1000, r: 2e-5, converterLossFraction: 0.01 },
    costAc: (L) => 100 + 1.0 * L / 1000, costDc: (L) => 250 + 0.6 * L / 1000 });
  assert.ok(out.Pac[0] > out.Pdc[0]);
  assert.ok(out.breakEvenCapacity > 20e3 && out.breakEvenCapacity < 220e3);
  // hand check of the cost crossing: at 100 km AC costs 200/197.4 MW = 1.013 per MW, DC 310/312.8 MW = 0.991 per MW
  assert.ok(out.breakEvenCost > 60e3 && out.breakEvenCost < 100e3, `cost break-even ${out.breakEvenCost}`);
  assert.equal(out.Pac[3], 0);
  for (let i = 1; i < lengths.length; i++) assert.ok(out.lossFractionDc[i] > out.lossFractionDc[i - 1]);
});

test('per-unit helpers: 132 kV, 100 MVA → Zbase 174.24 Ω, Ibase 437.39 A; change of base', () => {
  const b = perUnitBase({ S: 100e6, U: 132e3 });
  rel(b.Zbase, 174.24, 1e-12); rel(b.Ibase, 437.386, 1e-5);
  rel(perUnitChangeBase({ z: 0.12, Sold: 50e6, Uold: 132e3, Snew: 100e6, Unew: 132e3 }), 0.24, 1e-12);
});

test('three-phase power relations: 800 A at 132 kV, pf 0.95 lagging', () => {
  const p = threePhasePower({ U: 132e3, I: 800, powerFactor: 0.95 });
  rel(p.S, Math.sqrt(3) * 132e3 * 800, 1e-12); rel(p.P, 0.95 * p.S, 1e-12);
  rel(p.Q, p.S * Math.sqrt(1 - 0.9025), 1e-12);
  assert.ok(threePhasePower({ U: 132e3, P: 100e6, powerFactor: 0.95, leading: true }).Q < 0);
});

test('transformer by hand: uk 12 %, Pk 0.5 %, pf 0.9 lag → regulation |1.056762 + j0.105727| − 1 = 6.204 %; efficiency 99.393 %', () => {
  const t = transformer({ Sr: 100e6, uk: 0.12, Pk: 500e3, P0: 50e3, loading: 1, powerFactor: 0.9, U: 132e3 });
  rel(t.er, 0.005, 1e-12); rel(t.ex, 0.119896, 1e-5);
  rel(t.regulation, 0.062038, 1e-4);
  rel(t.efficiency, 90e6 / 90.55e6, 1e-12);
  rel(t.loadingAtMaxEfficiency, Math.sqrt(0.1), 1e-12);
  rel(t.Z.im, 0.119896 * 174.24, 1e-5);
  // efficiency at the optimum loading is the highest
  const at = (x) => transformer({ Sr: 100e6, uk: 0.12, Pk: 500e3, P0: 50e3, loading: x, powerFactor: 0.9 }).efficiency;
  assert.ok(at(Math.sqrt(0.1)) > at(0.25) && at(Math.sqrt(0.1)) > at(0.4));
  rel(t.maxEfficiency, at(Math.sqrt(0.1)), 1e-12);
});

test('loss-load factor from a load-duration curve by hand: loads 1, 0.5, 0.5, 0 → LF 0.5, LLF 0.375; annual loss', () => {
  const l = lossLoadFactor({ load: [1, 0.5, 0.5, 0], lossAtPeak: 100e3 });
  rel(l.loadFactor, 0.5, 1e-12); rel(l.lossLoadFactor, 0.375, 1e-12);
  rel(l.energyLoss, 100e3 * 0.375 * 8760, 1e-12);
  rel(lossLoadFactor({ load: [1, 0.2], weights: [1, 3] }).lossLoadFactor, (1 + 3 * 0.04) / 4, 1e-12);
  // LF² ≤ LLF ≤ LF always
  assert.ok(l.lossLoadFactor >= l.loadFactor ** 2 && l.lossLoadFactor <= l.loadFactor);
});

test("Kelvin's law: at the economic section the annualised section-dependent capital equals the annual cost of losses", () => {
  rel(capitalRecoveryFactor({ rate: 0.08, years: 25 }), 0.0936788, 1e-5);
  rel(capitalRecoveryFactor({ rate: 0, years: 20 }), 0.05, 1e-12);
  const p = { Imax: 700, rho: 2.2e-8, lossLoadFactor: 0.35, energyPrice: 60e-6, costPerArea: 2.5e5, fixedCost: 300, crf: 0.0936788 };
  const e = economicConductorSize(p);
  const lossCost = (3 * 700 ** 2 * 2.2e-8 * 0.35 * 8760 * 60e-6) / e.areaOptimal;
  rel(p.crf * p.costPerArea * e.areaOptimal, lossCost, 1e-10);
  // the continuous optimum is a minimum of the annual cost
  const cost = (S) => p.crf * (p.fixedCost + p.costPerArea * S) + (3 * 700 ** 2 * 2.2e-8 * 0.35 * 8760 * 60e-6) / S;
  assert.ok(cost(e.areaOptimal) < cost(0.9 * e.areaOptimal) && cost(e.areaOptimal) < cost(1.1 * e.areaOptimal));
  rel(e.annualCostOptimal, cost(e.areaOptimal), 1e-12);
  // discrete candidates: the cheapest standard section brackets the continuous optimum
  const sizes = [240, 400, 630, 800, 1000, 1200, 1600].map((a) => ({ area: a * 1e-6, cost: 300 + 2.5e5 * a * 1e-6 }));
  const d = economicConductorSize({ ...p, candidates: sizes });
  assert.ok(Math.abs(Math.log(d.best.area / e.areaOptimal)) < Math.log(1.45));
  assert.equal(d.best.annualCost, Math.min(...d.table.map((r) => r.annualCost)));
});
