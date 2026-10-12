import test from 'node:test';
import assert from 'node:assert/strict';
import {
  sdofResponse, responseSpectrum, dynamicAmplification, ec8Spectrum, EC8_SPECTRUM_PARAMETERS, isoSeismicSpectrum, cqcCoefficient,
  combineModal, combineDirections, responseSpectrumAnalysis, cqcMatrix, modalTimeHistory, directTimeHistory, harmonicResponse,
  spectrumCompatibleAccelerogram, campbellDiagram,
} from '../../src/engine/structural/framedynamics.js';
import { modalAnalysis } from '../../src/engine/structural/framemodal.js';
import { tubeSection } from '../../src/engine/structural/sections.js';

const E = 2.1e11;
const tube = tubeSection({ D: 0.3, t: 0.01 });
const EI = E * tube.Iy;
const near = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

// massless vertical cantilever with a tip mass, free to sway in x only: k = 3EI/L³
const Ls = 8; const Mtip = 4000;
const kS = 3 * EI / Ls ** 3;
const wS = Math.sqrt(kS / Mtip);
const sdofModel = {
  nodes: [[0, 0, 0], [0, 0, Ls]], members: [{ a: 0, b: 1, section: tube, material: { E, rho: 1e-9 } }],
  supports: [{ node: 0, fix: 'fixed' }, { node: 1, fix: [0, 1, 1, 1, 0, 1] }], nodeMasses: [{ node: 1, m: Mtip, inertia: [0, 1e-6, 0] }],
};
// cantilever with distributed mass, in-plane (x–z plane bending, sway in x)
function cantilever(n, L = 12) {
  const nodes = []; const members = []; const supports = [];
  for (let i = 0; i <= n; i++) {
    nodes.push([0, 0, L * i / n]);
    if (i < n) members.push({ a: i, b: i + 1, section: tube });
    supports.push({ node: i, fix: i === 0 ? 'fixed' : [0, 1, 1, 1, 0, 1] });
  }
  return { nodes, members, supports };
}

test('SDOF recurrence is exact for a step load: u(t) = (a/ω²)[1 − e^(−ζωt)(cos ωd t + ζ/√(1−ζ²) sin ωd t)]', () => {
  const omega = 7.3; const zeta = 0.08; const dt = 0.05; const a0 = 2.5;
  const r = sdofResponse({ omega, zeta, dt, a: new Float64Array(200).fill(a0) });
  const wd = omega * Math.sqrt(1 - zeta * zeta);
  for (const i of [1, 17, 60, 199]) {
    const t = i * dt;
    const exact = a0 / omega ** 2 * (1 - Math.exp(-zeta * omega * t) * (Math.cos(wd * t) + zeta / Math.sqrt(1 - zeta * zeta) * Math.sin(wd * t)));
    near(r.u[i], exact, 1e-10, `step ${i}`);
  }
  // ramp load, undamped: u = (a'/ω²)(t − sin ωt/ω)
  const ramp = Float64Array.from({ length: 100 }, (_, i) => 3 * i * dt);
  const rr = sdofResponse({ omega, zeta: 0, dt, a: ramp });
  near(rr.u[99], 3 / omega ** 2 * (99 * dt - Math.sin(omega * 99 * dt) / omega), 1e-10);
});

test('response spectrum of a suddenly applied constant ground acceleration: Sd = 2a/ω² (undamped)', () => {
  const a0 = 1.7; const dt = 0.002;
  const acc = new Float64Array(4000).fill(a0);
  const sp = responseSpectrum({ acc, dt, periods: [0.3, 1.0, 2.5], zeta: 0 });
  [0.3, 1.0, 2.5].forEach((T, i) => {
    const w = 2 * Math.PI / T;
    near(sp.Sd[i], 2 * a0 / w ** 2, 2e-4);
    near(sp.Sa[i], 2 * a0, 2e-4);
  });
});

test('dynamic amplification factor: static, resonance 1/(2ζ), isolation range', () => {
  near(dynamicAmplification({ ratio: 0, zeta: 0.05 }).daf, 1, 1e-12);
  near(dynamicAmplification({ ratio: 1, zeta: 0.05 }).daf, 10, 1e-12);
  near(dynamicAmplification({ ratio: 1, zeta: 0.05 }).phase, Math.PI / 2, 1e-12);
  near(dynamicAmplification({ ratio: Math.SQRT2, zeta: 0.3 }).transmissibility, 1, 1e-12);
  near(dynamicAmplification({ ratio: 3, zeta: 0 }).daf, 1 / 8, 1e-12);
});

test('EN 1998-1 elastic spectra: Table 3.2 / 3.3 values and branch formulas', () => {
  assert.deepEqual(EC8_SPECTRUM_PARAMETERS[1].C, { S: 1.15, TB: 0.20, TC: 0.6, TD: 2.0 });
  assert.deepEqual(EC8_SPECTRUM_PARAMETERS[1].E, { S: 1.4, TB: 0.15, TC: 0.5, TD: 2.0 });
  assert.deepEqual(EC8_SPECTRUM_PARAMETERS[2].D, { S: 1.8, TB: 0.10, TC: 0.30, TD: 1.2 });
  const ag = 2.0;
  const s = ec8Spectrum({ ag, groundType: 'C', type: 1 });
  near(s.Sa(0), ag * 1.15, 1e-12);
  near(s.Sa(0.1), ag * 1.15 * (1 + 0.5 * 1.5), 1e-12);
  near(s.Sa(0.4), 2.5 * ag * 1.15, 1e-12);
  near(s.Sa(1.0), 2.5 * ag * 1.15 * 0.6, 1e-12);
  near(s.Sa(3.0), 2.5 * ag * 1.15 * 0.6 * 2 / 9, 1e-12);
  // 2 % damping: η = √(10/7); 28 % damping hits the floor η = 0.55
  near(ec8Spectrum({ ag, zeta: 0.02 }).eta, Math.sqrt(10 / 7), 1e-12);
  near(ec8Spectrum({ ag, zeta: 0.30 }).eta, 0.55, 1e-12);
  // design spectrum with q = 3 on the plateau and its lower bound β ag
  const d = ec8Spectrum({ ag, groundType: 'B', q: 3 });
  near(d.Sa(0.3), 2.5 * ag * 1.2 / 3, 1e-12);
  near(d.Sa(0), ag * 1.2 * 2 / 3, 1e-12);
  near(d.Sa(10), 0.2 * ag, 1e-12);
});

test('ISO 19901-2 spectrum shape', () => {
  const s = isoSeismicSpectrum({ Sa02: 5, Sa10: 2, Ca: 1.2, Cv: 1.6 });
  near(s.Sa(0), 0.4 * 1.2 * 5, 1e-12);
  near(s.Sa(0.2), 1.2 * 5, 1e-12);
  near(s.Sa(0.4), 1.2 * 5, 1e-12, 'capped at the plateau');
  near(s.Sa(1), 1.6 * 2, 1e-12);
  near(s.Sa(4), 1.6 * 2 / 4, 1e-12);
  near(s.Sa(8), 4 * 1.6 * 2 / 64, 1e-12);
  near(s.cornerPeriod, 1.6 * 2 / 6, 1e-12);
  // 2 % damping: D = ln(50)/ln(20)
  near(isoSeismicSpectrum({ Sa02: 5, Sa10: 2, zeta: 0.02 }).Sa(1), 2 / (Math.log(50) / Math.log(20)), 1e-12);
});

test('CQC coefficient: hand value for r = 0.9, ζ = 5 %; CQC → SRSS for well-separated modes', () => {
  // ρ = 8ζ²(1 + r) r^1.5 / [(1 − r²)² + 4ζ² r (1 + r)²] = 0.032445/0.068590 = 0.4730
  near(cqcCoefficient(10, 9, 0.05), 0.4730, 5e-4);
  near(cqcCoefficient(10, 9, 0.05), cqcCoefficient(9, 10, 0.05), 1e-12, 'symmetric');
  assert.equal(cqcCoefficient(5, 5, 0.05), 1);
  assert.ok(cqcCoefficient(1, 5, 0.05) < 0.01);
  const peaks = [3, -2, 1.5]; const omegas = [2, 11, 40];
  const cqc = combineModal({ peaks, omegas, zeta: 0.05, method: 'CQC' });
  const srss = combineModal({ peaks, method: 'SRSS' });
  near(srss, Math.sqrt(9 + 4 + 2.25), 1e-12);
  near(cqc, srss, 5e-3);
  // identical frequencies: CQC gives the absolute sum of same-sign responses
  near(combineModal({ peaks: [2, 3], omegas: [7, 7], zeta: 0.05 }), 5, 1e-12);
  near(combineDirections([3, 4], 'SRSS'), 5, 1e-12);
  near(combineDirections([10, 4, 2], '100-30-30'), 10 + 0.3 * 6, 1e-12);
});

test('response-spectrum analysis of an SDOF structure is exact: u = Sa/ω², V = M·Sa', () => {
  const modal = modalAnalysis(sdofModel, { nModes: 1, shearDeformation: false });
  near(modal.omegas[0], wS, 1e-6);
  const spec = ec8Spectrum({ ag: 1.5, groundType: 'B' });
  const Sa = spec.Sa(2 * Math.PI / wS);
  for (const combination of ['CQC', 'SRSS']) {
    const r = responseSpectrumAnalysis(modal, { spectrum: spec.Sa, directions: ['x'], combination });
    near(r.displacements[6], Sa / wS ** 2, 1e-5, 'tip displacement');
    near(r.baseShear.x, Mtip * Sa, 1e-5, 'base shear');
    near(r.memberForces[0][4], Mtip * Sa * Ls, 1e-5, 'base moment');
    near(r.massParticipation.x, 1, 1e-6);
  }
});

test('response-spectrum analysis of a cantilever: CQC ≈ SRSS (separated modes); base shear bounded by M·Sa,max', () => {
  const modal = modalAnalysis(cantilever(10), { nModes: 6, shearDeformation: false });
  assert.ok(modal.frequencies[1] / modal.frequencies[0] > 6);
  const spec = ec8Spectrum({ ag: 2, groundType: 'C' });
  const cqc = responseSpectrumAnalysis(modal, { spectrum: spec.Sa, directions: ['x'], combination: 'CQC' });
  const srss = responseSpectrumAnalysis(modal, { spectrum: spec.Sa, directions: ['x'], combination: 'SRSS' });
  near(cqc.baseShear.x, srss.baseShear.x, 1e-2);
  near(cqc.displacements[6 * 10], srss.displacements[6 * 10], 1e-2);
  // hand check of the SRSS base shear from the effective masses
  let s = 0;
  modal.effectiveMass.forEach((me, k) => { s += (me.x * spec.Sa(modal.periods[k])) ** 2; });
  near(srss.baseShear.x, Math.sqrt(s), 1e-9);
  assert.ok(srss.baseShear.x < modal.totalMass.x * 2.5 * 2 * 1.15);
  // two horizontal directions with the same spectrum on an axisymmetric-in-plan restraint: only x is free, so y adds nothing
  const both = responseSpectrumAnalysis(modal, { spectrum: spec.Sa, directions: ['x', 'y'], directional: '100-30-30' });
  near(both.baseShear.x, cqc.baseShear.x, 1e-9);
});

test('modal time history of the SDOF structure: step force u = (F/k)(1 − cos ωt); ground step u = −(a/ω²)(1 − cos ωt)', () => {
  const modal = modalAnalysis(sdofModel, { nModes: 1, shearDeformation: false });
  const F = 5e3; const dt = 0.01; const nSteps = 400;
  const r = modalTimeHistory(modal, { dt, nSteps, zeta: 0, nodalForces: [{ node: 1, dof: 'ux', history: () => F }], outputs: [{ node: 1, dof: 'ux' }] });
  for (const i of [10, 123, 400]) near(r.outputs[0][i], F / kS * (1 - Math.cos(wS * i * dt)), 1e-5);
  near(r.peak[6], 2 * F / kS, 1e-3, 'dynamic load factor 2 for a suddenly applied load');
  const ag = 0.8;
  const g = modalTimeHistory(modal, { dt, nSteps, zeta: 0, baseAcc: { x: () => ag }, outputs: [{ node: 1, dof: 0 }] });
  near(g.outputs[0][250], -ag / wS ** 2 * (1 - Math.cos(wS * 250 * dt)), 1e-5);
});

test('direct Newmark / HHT integration agrees with modal superposition for a cantilever under a step load', () => {
  const model = cantilever(6);
  const modal = modalAnalysis(model, { nModes: 12, shearDeformation: false });
  const T1 = modal.periods[0];
  const dt = T1 / 400; const nSteps = 800; const F = 2e4;
  const forces = [{ node: 6, dof: 'ux', history: () => F }];
  const outputs = [{ node: 6, dof: 'ux' }];
  const ms = modalTimeHistory(modal, { dt, nSteps, zeta: 0, nodalForces: forces, outputs });
  const nm = directTimeHistory(model, { dt, nSteps, nodalForces: forces, outputs, system: modal.system });
  const uStatic = F * 12 ** 3 / (3 * EI);
  let worst = 0;
  for (let i = 0; i <= nSteps; i += 20) worst = Math.max(worst, Math.abs(ms.outputs[0][i] - nm.outputs[0][i]));
  assert.ok(worst < 0.01 * uStatic, `largest difference ${worst / uStatic} of the static deflection`);
  // peak ≈ 2 × static (first mode dominates the tip response)
  assert.ok(ms.peak[36] > 1.95 * uStatic && ms.peak[36] < 2.1 * uStatic);
  // HHT-α damps the high modes but keeps the low-frequency response
  const hht = directTimeHistory(model, { dt, nSteps, alphaHHT: -0.1, nodalForces: forces, outputs, system: modal.system });
  near(hht.outputs[0][200], nm.outputs[0][200], 2e-2);
  // Rayleigh damping makes the response settle on the static deflection
  const damped = directTimeHistory(model, {
    dt: T1 / 50, nSteps: 2000, rayleigh: { alpha: 2 * 0.2 * modal.omegas[0], beta: 0 }, nodalForces: forces, outputs, system: modal.system,
  });
  near(damped.outputs[0][2000], uStatic, 1e-3);
});

test('harmonic response of the SDOF structure: amplitude (F/k)·DAF, resonance F/(2ζk), phase', () => {
  const modal = modalAnalysis(sdofModel, { nModes: 1, shearDeformation: false });
  const fn = wS / (2 * Math.PI); const zeta = 0.03; const F = 1e3;
  const r = harmonicResponse(modal, {
    frequencies: [0, 0.5 * fn, fn, 2 * fn], forces: [{ node: 1, dof: 'ux', amplitude: F }], outputs: [{ node: 1, dof: 'ux' }], zeta,
  });
  near(r.amplitude[0][0], F / kS, 1e-5);
  near(r.amplitude[0][1], F / kS * dynamicAmplification({ ratio: 0.5, zeta }).daf, 1e-5);
  near(r.amplitude[0][2], F / (2 * zeta * kS), 1e-5);
  near(r.phase[0][2], Math.PI / 2, 1e-6);
  near(r.amplitude[0][3], F / kS * dynamicAmplification({ ratio: 2, zeta }).daf, 1e-5);
  assert.ok(r.phase[0][3] > Math.PI / 2);
});

test('spectrum-compatible accelerogram matches the target spectrum and is reproducible for a seed', () => {
  const spec = ec8Spectrum({ ag: 2.5, groundType: 'B' });
  const opts = { spectrum: spec.Sa, duration: 20, dt: 0.01, fMin: 0.4, fMax: 20, nFreq: 60, iterations: 12, seed: 7 };
  const a = spectrumCompatibleAccelerogram(opts);
  assert.ok(a.meanError < 0.05, `mean spectral error ${a.meanError}`);
  assert.ok(a.maxError < 0.20, `max spectral error ${a.maxError}`);
  assert.ok(a.pga > 0.5 * 2.5 * 1.2 && a.pga < 2.5 * 2.5 * 1.2, `pga ${a.pga}`);
  assert.equal(Math.abs(a.acc[0]), 0);
  const b = spectrumCompatibleAccelerogram(opts);
  assert.deepEqual(Array.from(a.acc.subarray(100, 110)), Array.from(b.acc.subarray(100, 110)));
  const c = spectrumCompatibleAccelerogram({ ...opts, seed: 8 });
  assert.notEqual(a.acc[500], c.acc[500]);
  // independent check of the achieved spectrum with the public response-spectrum routine
  const sp = responseSpectrum({ acc: a.acc, dt: 0.01, periods: [0.5], zeta: 0.05 });
  near(sp.Sa[0], spec.Sa(0.5), 0.25);
});

test('Campbell diagram: 1× critical speed 60 f, interference inside the ±10 % band', () => {
  const c = campbellDiagram({ naturalFrequencies: [25, 80], orders: [1, 2], speedMin: 1200, speedMax: 1700, margin: 0.1 });
  const cs = c.criticalSpeeds.find((x) => x.mode === 1 && x.order === 1);
  near(cs.speed, 1500, 1e-12);
  assert.ok(cs.inRange);
  assert.ok(!c.ok);
  assert.ok(c.interferences.some((x) => x.mode === 1 && x.order === 1 && x.separation === 0));
  // 2× line covers 40–56.7 Hz: 25 Hz mode is clear of it by (40 − 25)/25 = 60 %
  assert.ok(!c.interferences.some((x) => x.mode === 1 && x.order === 2));
  const ok = campbellDiagram({ naturalFrequencies: [80], orders: [1, 2], speedMin: 1200, speedMax: 1700, margin: 0.1 });
  assert.ok(ok.ok);
});

test('CQC correlation matrix is computed once and reproduces the pairwise combination', () => {
  const omegas = [5, 5.2, 9, 14, 14.1, 30];
  const peaks = [1.5, -0.8, 0.6, 0.3, 0.25, -0.1];
  const rho = cqcMatrix(omegas, 0.05);
  assert.equal(rho.length, 36);
  for (let i = 0; i < 6; i++) {
    assert.equal(rho[i * 6 + i], 1);
    for (let j = 0; j < 6; j++) {
      assert.equal(rho[i * 6 + j], rho[j * 6 + i]);
      if (i !== j) near(rho[i * 6 + j], cqcCoefficient(omegas[i], omegas[j], 0.05), 1e-12);
    }
  }
  // equal damping: ρ = 8ζ²(1 + r) r^1.5 / [(1 − r²)² + 4ζ² r (1 + r)²]; r = 5.2/5 = 1.04, ζ = 0.05 → 0.8658
  const r = 1.04;
  near(rho[1], 8 * 0.0025 * (1 + r) * r ** 1.5 / ((1 - r * r) ** 2 + 4 * 0.0025 * r * (1 + r) ** 2), 1e-12);
  near(rho[1], 0.8658, 1e-3);
  const direct = combineModal({ peaks, omegas, zeta: 0.05 });
  near(combineModal({ peaks, rho }), direct, 1e-13);
  // hand value of the double sum
  let s = 0;
  for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) s += rho[i * 6 + j] * peaks[i] * peaks[j];
  near(direct, Math.sqrt(s), 1e-13);
  // per-mode damping
  const z = [0.02, 0.05, 0.05, 0.03, 0.03, 0.08];
  near(combineModal({ peaks, rho: cqcMatrix(omegas, z) }), combineModal({ peaks, omegas, zeta: z }), 1e-13);
  assert.throws(() => combineModal({ peaks, rho: new Float64Array(4) }), RangeError);
});
