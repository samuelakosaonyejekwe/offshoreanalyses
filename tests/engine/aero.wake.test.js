import test from 'node:test';
import assert from 'node:assert/strict';
import {
  inductionFromCt, jensenDeficit, frandsenDeficit, wakeGrowthFromTI, gaussianDeficit, yawDeflection,
  meanderedGaussian, addedTurbulenceFrandsen, addedTurbulenceCrespo, effectiveTurbulence, superpose,
  circleOverlapArea, tophatRotorAverage, rotorPoints, gaussianRotorAverage, rotorEquivalentWindSpeed,
  inductionOnAxis, inductionUpstream, deepArrayFrandsen, farmFlow,
} from '../../src/engine/aero/wake.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);

test('Jensen deficit by hand: Ct = 0.8, k = 0.05, D = 100 m, x = 500 m', () => {
  // (1 − √0.2) / (1 + 0.05·500/50)² = 0.5527864 / 2.25 = 0.2456828; Dw = 100 + 2·0.05·500 = 150
  const w = jensenDeficit({ Ct: 0.8, x: 500, D: 100, k: 0.05 });
  close(w.deficit, 0.2456828, 1e-6);
  close(w.wakeDiameter, 150, 1e-12);
});

test('Jensen wake satisfies its mass balance at every distance', () => {
  const Ct = 0.7;
  const D = 120;
  const R = D / 2;
  for (const x of [0, 200, 600, 1500, 4000]) {
    const w = jensenDeficit({ Ct, x, D, k: 0.04 });
    const Rw = w.wakeDiameter / 2;
    // π R² u_r + π (Rw² − R²) U = π Rw² u, with u_r = U √(1 − Ct) and U = 1
    const lhs = R * R * Math.sqrt(1 - Ct) + (Rw * Rw - R * R);
    const rhs = Rw * Rw * (1 - w.deficit);
    close(lhs, rhs, 1e-12, `x = ${x}`);
  }
});

test('Frandsen storage-area wake conserves the momentum deficit', () => {
  const Ct = 0.75;
  const D = 100;
  for (const x of [200, 700, 2000]) {
    const w = frandsenDeficit({ Ct, x, D, alpha: 0.7 });
    const u = 1 - w.deficit;
    // A u (U − u) = Ct A0 U²/2
    close((w.wakeDiameter / D) ** 2 * u * (1 - u), Ct / 2, 1e-12);
  }
  // at x = 0 the area is β A0
  const w0 = frandsenDeficit({ Ct, x: 0, D, alpha: 0.7 });
  close((w0.wakeDiameter / D) ** 2, w0.beta, 1e-12);
});

test('Gaussian wake: ε and k* by hand', () => {
  // Ct = 0.75: √(1 − Ct) = 0.5, β = 0.5·1.5/0.5 = 1.5, ε = 0.2 √1.5 = 0.2449490
  const w = gaussianDeficit({ Ct: 0.75, x: 0, D: 100, kStar: 0.03 });
  close(w.epsilon, 0.244949, 1e-6);
  close(wakeGrowthFromTI(0.1), 0.3837 * 0.1 + 0.003678, 1e-15);
  close(wakeGrowthFromTI(0.1), 0.042048, 1e-12);
  // centre deficit at x = 6D, k* = 0.04: σ/D = 0.24 + 0.244949 = 0.484949; C = 1 − √(1 − 0.75/(8·0.2351755)) = 0.2246...
  const c = gaussianDeficit({ Ct: 0.75, x: 600, D: 100, kStar: 0.04 });
  close(c.sigma, 48.4949, 1e-5);
  close(c.centreDeficit, 1 - Math.sqrt(1 - 0.75 / (8 * 0.484949 ** 2)), 1e-6);
  assert.equal(c.nearWake, false);
});

test('Gaussian wake momentum-deficit integral equals the thrust: ∫ΔU(U∞ − ΔU)dA = ½ Ct A U∞²', () => {
  const D = 126;
  const U = 9;
  for (const [Ct, xD] of [[0.8, 6], [0.8, 10], [0.5, 4], [0.3, 15]]) {
    const x = xD * D;
    const { sigma } = gaussianDeficit({ Ct, x, D, kStar: 0.035 });
    const n = 20000;
    const rMax = 12 * sigma;
    const dr = rMax / n;
    let s = 0;
    for (let i = 0; i < n; i++) {
      const r = (i + 0.5) * dr;
      const dU = U * gaussianDeficit({ Ct, x, D, kStar: 0.035, r }).deficit;
      s += dU * (U - dU) * 2 * Math.PI * r * dr;
    }
    const thrustOverRho = 0.5 * Ct * ((Math.PI * D * D) / 4) * U * U;
    close(s, thrustOverRho, 1e-5, `Ct ${Ct}, x/D ${xD}`);
  }
});

test('Gaussian near wake is flagged and bounded', () => {
  const w = gaussianDeficit({ Ct: 0.95, x: 10, D: 100, kStar: 0.02 });
  assert.equal(w.nearWake, true);
  assert.equal(w.centreDeficit, 1);
});

test('yawed-wake deflection: zero for zero yaw, odd in yaw, continuous at the core end, grows downstream', () => {
  const o = { Ct: 0.8, D: 100, TI: 0.08, ky: 0.03 };
  assert.equal(yawDeflection({ ...o, yawDeg: 0, x: 800 }).deflection, 0);
  const p = yawDeflection({ ...o, yawDeg: 20, x: 800 });
  const m = yawDeflection({ ...o, yawDeg: -20, x: 800 });
  close(p.deflection, -m.deflection, 1e-12);
  const a = yawDeflection({ ...o, yawDeg: 20, x: p.x0 * 0.999999 }).deflection;
  const b = yawDeflection({ ...o, yawDeg: 20, x: p.x0 * 1.000001 }).deflection;
  close(a, b, 1e-4);
  assert.ok(yawDeflection({ ...o, yawDeg: 20, x: 1500 }).deflection > p.deflection);
  // θc0 by hand: 0.3·γ/cos γ (1 − √(1 − Ct cos γ)), γ = 20° = 0.349066 rad
  close(p.theta0, ((0.3 * 0.349066) / Math.cos(0.349066)) * (1 - Math.sqrt(1 - 0.8 * Math.cos(0.349066))), 1e-6);
  // order of magnitude: a 20° yaw deflects the wake by a few tenths of a diameter at 8D
  assert.ok(p.deflection / 100 > 0.2 && p.deflection / 100 < 0.8, `deflection/D = ${p.deflection / 100}`);
});

test('meandering preserves the integrated deficit of a Gaussian wake', () => {
  const r = meanderedGaussian({ centreDeficit: 0.3, sigma: 40, sigmaMeanderY: 30, sigmaMeanderZ: 10 });
  close(r.sigmaY, 50, 1e-12);
  close(r.centreDeficit * r.sigmaY * r.sigmaZ, 0.3 * 40 * 40, 1e-12);
});

test('added turbulence: Frandsen and Crespo–Hernández by hand', () => {
  // Frandsen: 1/(1.5 + 0.8·5/√0.8) = 1/(1.5 + 4.472136) = 0.167444
  close(addedTurbulenceFrandsen({ Ct: 0.8, xOverD: 5 }), 0.167444, 1e-5);
  // Crespo: a = ½(1 − √0.2) = 0.276393; 0.73 · a^0.8325 · 0.1^0.0325 · 7^−0.32
  //   = 0.73 · 0.34283 · 0.92790 · 0.53650 = 0.12458
  const c = addedTurbulenceCrespo({ Ct: 0.8, xOverD: 7, I0: 0.1 });
  close(c.added, 0.12458, 2e-4);
  assert.equal(c.inRange, true);
  assert.equal(addedTurbulenceCrespo({ Ct: 0.8, xOverD: 3, I0: 0.1 }).inRange, false);
});

test('IEC 61400-1 Annex E effective turbulence', () => {
  close(effectiveTurbulence({ I0: 0.1, neighbours: [], m: 4 }).Ieff, 0.1, 1e-12);
  // one neighbour at 5D, Ct = 0.8: I_T = √(0.167444² + 0.1²) = 0.195030
  // I_eff = (0.94·0.1⁴ + 0.06·0.195030⁴)^(1/4) = (9.4e-5 + 8.6806e-5)^(1/4) = 0.115925
  const r = effectiveTurbulence({ I0: 0.1, neighbours: [{ distanceD: 5, Ct: 0.8 }], m: 4 });
  close(r.wakeTI[0], 0.19503, 1e-4);
  close(r.Ieff, 0.115925, 1e-4);
  // larger Wöhler exponent weights the wake condition more
  assert.ok(effectiveTurbulence({ I0: 0.1, neighbours: [{ distanceD: 5, Ct: 0.8 }], m: 10 }).Ieff > r.Ieff);
});

test('wake superposition rules by hand', () => {
  const d = [0.2, 0.1];
  close(superpose(d, 'linear'), 0.3, 1e-12);
  close(superpose(d, 'rss'), Math.sqrt(0.05), 1e-12);
  close(superpose(d, 'max'), 0.2, 1e-12);
  // energy: 1 − (1 − δ)² = 0.36 + 0.19 = 0.55 → δ = 1 − √0.45
  close(superpose(d, 'energy'), 1 - Math.sqrt(0.45), 1e-12);
  assert.equal(superpose([0.8, 0.9], 'linear'), 1);
  assert.throws(() => superpose(d, 'cubic'), RangeError);
});

test('circle overlap: limits and the equal-circle closed form', () => {
  close(circleOverlapArea(1, 1, 1), (2 * Math.PI) / 3 - Math.sqrt(3) / 2, 1e-12);
  close(circleOverlapArea(3, 1, 0.5), Math.PI, 1e-12);
  assert.equal(circleOverlapArea(1, 1, 2.5), 0);
  const r = tophatRotorAverage({ deficit: 0.3, wakeRadius: 100, rotorRadius: 50, offset: 100 });
  close(r.overlapFraction, circleOverlapArea(100, 50, 100) / (Math.PI * 2500), 1e-12);
  close(r.deficit, 0.3 * r.overlapFraction, 1e-12);
});

test('Gaussian rotor average converges to the closed form on the wake axis', () => {
  // disc average of exp(−r²/2σ²) over radius R: (2σ²/R²)(1 − exp(−R²/2σ²))
  const sigma = 60;
  const R = 50;
  const exact = ((2 * sigma * sigma) / (R * R)) * (1 - Math.exp((-R * R) / (2 * sigma * sigma)));
  const coarse = gaussianRotorAverage({ centreDeficit: 1, sigma, dy: 0, rotorRadius: R, points: rotorPoints(2, 6) });
  const fine = gaussianRotorAverage({ centreDeficit: 1, sigma, dy: 0, rotorRadius: R, points: rotorPoints(16, 24) });
  assert.ok(Math.abs(fine - exact) < Math.abs(coarse - exact));
  close(fine, exact, 2e-4);
  close(gaussianRotorAverage({ centreDeficit: 1, sigma, dy: 0, rotorRadius: R }), exact, 5e-3);
});

test('rotor-equivalent wind speed: uniform flow and power-law shear against polar quadrature', () => {
  close(rotorEquivalentWindSpeed({ profile: () => 8, hubHeight: 100, D: 120 }).rews, 8, 1e-12);
  const alpha = 0.2;
  const hub = 100;
  const R = 60;
  const profile = (z) => 10 * (z / hub) ** alpha;
  const r = rotorEquivalentWindSpeed({ profile, hubHeight: hub, D: 2 * R, nSeg: 400 });
  let s = 0;
  const nr = 400;
  const nt = 400;
  for (let i = 0; i < nr; i++) {
    const rr = ((i + 0.5) / nr) * R;
    for (let j = 0; j < nt; j++) {
      const th = ((j + 0.5) / nt) * 2 * Math.PI;
      s += profile(hub + rr * Math.sin(th)) ** 3 * rr * (R / nr) * ((2 * Math.PI) / nt);
    }
  }
  close(r.rews, Math.cbrt(s / (Math.PI * R * R)), 1e-5);
  assert.ok(Math.abs(r.ratio - 1) < 0.02);
});

test('vortex-cylinder induction: limits on the axis and far-field equivalence with the point source', () => {
  const a = 0.3;
  close(inductionOnAxis({ a, x: 0, R: 50 }), 1 - a, 1e-12);
  close(inductionOnAxis({ a, x: -1e7, R: 50 }), 1, 1e-9);
  close(inductionOnAxis({ a, x: 1e7, R: 50 }), 1 - 2 * a, 1e-9);
  // 2.5D upstream: 1 − a(1 − 5/√26) = 1 − 0.01942 a
  close(inductionOnAxis({ a, x: -250, R: 50 }), 1 - a * (1 - 5 / Math.sqrt(26)), 1e-12);
  const far = inductionUpstream({ a, x: -2000, r: 0, R: 50 });
  close(far, inductionOnAxis({ a, x: -2000, R: 50 }) - 1, 2e-3 * Math.abs(far) + 1e-9);
  assert.ok(Math.abs(far - (inductionOnAxis({ a, x: -2000, R: 50 }) - 1)) < 2e-3 * Math.abs(far));
  assert.equal(inductionUpstream({ a, x: 10, r: 0, R: 50 }), 0);
  assert.ok(Math.abs(inductionUpstream({ a, x: -200, r: 300, R: 50 })) < Math.abs(inductionUpstream({ a, x: -200, r: 0, R: 50 })));
});

test('deep-array roughness model: no thrust recovers the surface roughness; z00 by hand', () => {
  const o = { sx: 7, sy: 7, hubHeight: 100, z0: 2e-4, H: 600 };
  const none = deepArrayFrandsen({ ...o, Ct: 0 });
  close(none.z00, 2e-4, 1e-9);
  close(none.speedRatio, 1, 1e-12);
  // ct = π·0.8/(8·49) = 0.0064114; κ/ln(h/z0) = 0.4/13.1224 = 0.0304823
  // z00 = 100 exp(−0.4/√(0.0064114 + 0.00092917)) = 100 exp(−4.66868) = 0.93846
  const r = deepArrayFrandsen({ ...o, Ct: 0.8 });
  close(r.ct, 0.0064114, 1e-5);
  close(r.z00, 0.9385, 2e-3);
  assert.ok(r.speedRatio < 1 && r.speedRatio > 0.5);
  assert.ok(deepArrayFrandsen({ ...o, Ct: 0.8, sx: 10, sy: 10 }).speedRatio > r.speedRatio);
});

const turbine = (x, y) => ({
  x,
  y,
  D: 100,
  hubHeight: 90,
  ct: () => 0.8,
  power: (U) => 0.5 * 1.225 * ((Math.PI * 100 * 100) / 4) * 0.45 * U ** 3,
});

test('two turbines aligned, Jensen: downstream speed by hand', () => {
  // wind from the north (wd = 0): the turbine at y = 0 is 500 m downstream of the one at y = 500
  const r = farmFlow({ turbines: [turbine(0, 0), turbine(0, 500)], U: 10, wdDeg: 0, model: 'jensen', k: 0.05, turbulenceModel: 'none' });
  assert.deepEqual(Array.from(r.order), [1, 0]);
  close(r.U[1], 10, 1e-12);
  close(r.U[0], 10 * (1 - 0.2456828), 1e-6);
  close(r.efficiency, (1 + (1 - 0.2456828) ** 3) / 2, 1e-6);
});

test('two-turbine direction sweep is symmetric about the alignment direction, with the minimum when aligned', () => {
  const ts = [turbine(0, 0), turbine(0, 700)];
  for (const model of ['gaussian', 'jensen']) {
    const p = (wd) => farmFlow({ turbines: ts, U: 9, wdDeg: wd, model, TI: 0.08 }).power[0];
    const aligned = p(0);
    for (const d of [2, 5, 9, 15, 30]) {
      close(p(d), p(-d), 1e-9, `${model} ±${d}`);
      // (the top-hat wake is marginally stronger at small offsets because the streamwise distance shortens)
      assert.ok(p(d) >= aligned * (model === 'jensen' ? 0.999 : 1) - 1e-9);
    }
    // far outside the wake the downstream turbine is undisturbed
    close(p(90), turbine(0, 0).power(9), 1e-9);
    // reversing the wind swaps the roles
    close(farmFlow({ turbines: ts, U: 9, wdDeg: 180, model, TI: 0.08 }).power[1], aligned, 1e-9);
  }
});

test('farm solver: local turbulence rises in the wake, reflection and blockage reduce the yield', () => {
  const ts = [turbine(0, 0), turbine(0, 600), turbine(0, 1200)];
  const base = farmFlow({ turbines: ts, U: 9, wdDeg: 0, TI: 0.06 });
  assert.ok(base.TI[1] > 0.06 && base.TI[0] > 0.06);
  close(base.TI[2], 0.06, 1e-12);
  assert.ok(base.U[0] < base.U[2] && base.U[1] < base.U[2]);
  const refl = farmFlow({ turbines: ts, U: 9, wdDeg: 0, TI: 0.06, reflect: true });
  assert.ok(refl.totalPower < base.totalPower);
  const blk = farmFlow({ turbines: ts, U: 9, wdDeg: 0, TI: 0.06, blockage: true });
  assert.ok(blk.U[2] < 9 && blk.blockage[2] > 0 && blk.blockage[2] < 0.02);
  assert.ok(blk.totalPower < base.totalPower);
  assert.ok(blk.iterations >= 2);
  // linear superposition gives a larger combined deficit than root-sum-square
  const lin = farmFlow({ turbines: ts, U: 9, wdDeg: 0, TI: 0.06, superposition: 'linear' });
  assert.ok(lin.deficit[0] > base.deficit[0]);
  // tabulated curves are accepted
  const tab = farmFlow({ turbines: [{ x: 0, y: 0, D: 100, hubHeight: 90, ct: { u: [3, 25], v: [0.8, 0.8] }, power: { u: [3, 10, 25], v: [0, 2e6, 2e6] } }], U: 6.5, wdDeg: 10 });
  close(tab.totalPower, 1e6, 1e-12);
  assert.throws(() => farmFlow({ turbines: [], U: 9, wdDeg: 0 }), RangeError);
  assert.equal(inductionFromCt(0.75), 0.25);
});
