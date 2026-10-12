import test from 'node:test';
import assert from 'node:assert/strict';
import * as s from '../../src/engine/geotech/siteresponse.js';

const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);

test('uniform undamped layer on a rigid base: |F| = 1/|cos(2 pi f H/Vs)|, resonances at (2n - 1) Vs / 4H (Kramer eq. 7.9)', () => {
  const H = 20, Vs = 200;
  const layers = [{ h: H, rho: 1800, Vs, damping: 0 }];
  const freqs = [0.5, 1, 2, 3.3, 6, 9.1];
  const tf = s.transferFunction({ layers, bedrock: { rigid: true }, freqs });
  freqs.forEach((f, i) => close(tf.amp[i], 1 / Math.abs(Math.cos(2 * Math.PI * f * H / Vs)), 1e-9));
  // locate the first three peaks of a lightly damped layer on a fine frequency grid
  const fine = [];
  for (let f = 0.01; f <= 14; f += 0.005) fine.push(f);
  const d = s.transferFunction({ layers: [{ h: H, rho: 1800, Vs, damping: 0.01 }], bedrock: { rigid: true }, freqs: fine });
  const peaks = [];
  for (let i = 1; i < fine.length - 1; i++) if (d.amp[i] > d.amp[i - 1] && d.amp[i] > d.amp[i + 1]) peaks.push(fine[i]);
  assert.equal(peaks.length, 3);
  [2.5, 7.5, 12.5].forEach((fn, i) => close(peaks[i], fn, 2.5e-3));
});

test('damped layer on a rigid base: peak amplification 2 / ((2n - 1) pi xi) (Kramer eq. 7.13)', () => {
  const xi = 0.02;
  const tf = s.transferFunction({ layers: [{ h: 20, rho: 1800, Vs: 200, damping: xi }], bedrock: { rigid: true }, freqs: [2.5, 7.5] });
  close(tf.amp[0], 2 / (Math.PI * xi), 2e-3);
  close(tf.amp[1], 2 / (3 * Math.PI * xi), 5e-3);
});

test('undamped layer on elastic rock: outcrop amplification at resonance equals the impedance ratio; unity at antiresonance (Kramer eq. 7.20)', () => {
  const tf = s.transferFunction({ layers: [{ h: 20, rho: 1800, Vs: 200 }], bedrock: { rho: 2400, Vs: 1500 }, freqs: [2.5, 5, 1.3] });
  close(tf.amp[0], (2400 * 1500) / (1800 * 200), 1e-9);
  close(tf.amp[1], 1, 1e-9);
  const kH = 2 * Math.PI * 1.3 * 20 / 200;
  const az = (1800 * 200) / (2400 * 1500);
  close(tf.amp[2], 1 / Math.sqrt(Math.cos(kH) ** 2 + (az * Math.sin(kH)) ** 2), 1e-9);
  // 'within' input ignores the rock impedance
  close(s.transferFunction({ layers: [{ h: 20, rho: 1800, Vs: 200 }], bedrock: { rho: 2400, Vs: 1500 }, freqs: [1.3], input: 'within' }).amp[0], 1 / Math.abs(Math.cos(kH)), 1e-9);
});

test('sub-dividing a uniform layer does not change the transfer function; soft-over-stiff two-layer fundamental frequency from the frequency equation', () => {
  const freqs = [0.7, 1.9, 4.2, 8.8];
  const one = s.transferFunction({ layers: [{ h: 30, rho: 1900, Vs: 250, damping: 0.03 }], bedrock: { rho: 2300, Vs: 900, damping: 0.01 }, freqs });
  const three = s.transferFunction({ layers: [{ h: 5, rho: 1900, Vs: 250, damping: 0.03 }, { h: 10, rho: 1900, Vs: 250, damping: 0.03 }, { h: 15, rho: 1900, Vs: 250, damping: 0.03 }], bedrock: { rho: 2300, Vs: 900, damping: 0.01 }, freqs });
  freqs.forEach((f, i) => { close(three.re[i], one.re[i], 1e-9); close(three.im[i], one.im[i], 1e-9); });
  // two undamped layers on a rigid base: tan(w h1/V1) tan(w h2/V2) = rho2 V2 / (rho1 V1)
  const h1 = 10, V1 = 150, r1 = 1700, h2 = 20, V2 = 400, r2 = 2000;
  let lo = 0.1, hi = 3.7;
  const g = (f) => Math.tan(2 * Math.PI * f * h1 / V1) * Math.tan(2 * Math.PI * f * h2 / V2) - (r2 * V2) / (r1 * V1);
  for (let i = 0; i < 80; i++) { const mid = 0.5 * (lo + hi); if (g(mid) < 0) lo = mid; else hi = mid; }
  const f1 = 0.5 * (lo + hi);
  const tf = s.transferFunction({ layers: [{ h: h1, rho: r1, Vs: V1 }, { h: h2, rho: r2, Vs: V2 }], bedrock: { rigid: true }, freqs: [f1 * 0.999, f1, f1 * 1.001] });
  assert.ok(tf.amp[1] > 1e5 && tf.amp[1] > tf.amp[0] && tf.amp[1] > tf.amp[2], `amp ${tf.amp}`);
});

test('linear time-domain run: a harmonic input is scaled by |F(f0)|; strain at mid-layer matches the standing-wave solution', () => {
  const N = 4096, dt = 0.01;
  const f0 = 100 / (N * dt); // integer number of cycles: no leakage
  const acc = new Float64Array(N);
  for (let i = 0; i < N; i++) acc[i] = 0.5 * Math.sin(2 * Math.PI * f0 * i * dt);
  const H = 20, Vs = 200, xi = 0.05;
  const layers = [{ h: H, rho: 1800, Vs, damping: xi }];
  const r = s.equivalentLinearResponse({ layers, bedrock: { rigid: true }, acc, dt });
  const F = s.transferFunction({ layers, bedrock: { rigid: true }, freqs: [f0] }).amp[0];
  close(r.pgaSurface, 0.5 * F, 2e-3);
  close(r.amplification, F, 2e-3);
  assert.equal(r.converged, true);
  // standing wave u = U cos(k z): strain amplitude at mid-height = k U sin(k H/2), U = a_surface / w^2
  // (a trace of damping keeps the transfer function finite at the bin that coincides with the third mode)
  const und = s.equivalentLinearResponse({ layers: [{ h: H, rho: 1800, Vs, damping: 1e-5 }], bedrock: { rigid: true }, acc, dt });
  const w = 2 * Math.PI * f0, k = w / Vs;
  close(und.layers[0].gammaMax, k * (und.pgaSurface / (w * w)) * Math.sin(k * H / 2), 2e-3);
  // PGA at the top of a lower sub-layer follows cos(k z)
  const two = s.equivalentLinearResponse({ layers: [{ h: 8, rho: 1800, Vs, damping: 1e-5 }, { h: 12, rho: 1800, Vs, damping: 1e-5 }], bedrock: { rigid: true }, acc, dt, layerHistories: true });
  close(two.layers[1].pgaTop, two.pgaSurface * Math.abs(Math.cos(k * 8)), 2e-3);
});

test('equivalent-linear iteration: converges, softens the column under strong shaking and stays linear under weak shaking', () => {
  const N = 2048, dt = 0.01;
  const acc = new Float64Array(N);
  // deterministic broadband pulse train
  for (let i = 0; i < N; i++) { const t = i * dt; acc[i] = Math.exp(-((t - 5) ** 2) / 4) * (Math.sin(2 * Math.PI * 1.1 * t) + 0.6 * Math.sin(2 * Math.PI * 2.7 * t + 1) + 0.4 * Math.sin(2 * Math.PI * 5.3 * t + 2)); }
  const mk = () => [0, 1, 2, 3].map((i) => ({ h: 5, rho: 1850, Vs: 180 + 30 * i, curve: { model: 'darendeli', sigmaM: (25 + 50 * i) * 1e3, PI: 15 } }));
  const rock = { rho: 2400, Vs: 800, damping: 0.01 };
  const strong = s.equivalentLinearResponse({ layers: mk(), bedrock: rock, acc: acc.map((a) => 0.5 * a), dt, maxIter: 30 });
  const weak = s.equivalentLinearResponse({ layers: mk(), bedrock: rock, acc: acc.map((a) => 0.001 * a), dt });
  assert.equal(strong.converged, true);
  assert.equal(weak.converged, true);
  for (const l of strong.layers) assert.ok(l.ratio < 0.95 && l.ratio > 0.1 && l.damping > 0.02, `ratio ${l.ratio}, damping ${l.damping}`);
  for (const l of weak.layers) assert.ok(l.ratio > 0.99);
  // softening lowers the fundamental frequency of the column (first peak of the transfer function)
  const firstPeak = (tf) => { for (let i = 1; i < tf.f.length - 1; i++) if (tf.amp[i] > tf.amp[i - 1] && tf.amp[i] > tf.amp[i + 1]) return tf.f[i]; return NaN; };
  assert.ok(firstPeak(strong.tf) < firstPeak(weak.tf));
  // weak-motion fundamental frequency close to Vs,avg / 4H with the travel-time average velocity
  const tt = [0, 1, 2, 3].reduce((a, i) => a + 5 / (180 + 30 * i), 0);
  assert.ok(Math.abs(firstPeak(weak.tf) - 1 / (4 * tt)) / (1 / (4 * tt)) < 0.12);
  // strain-compatible properties are consistent with the curve at the effective strain
  const l0 = strong.layers[0];
  close(l0.gammaEff, 0.65 * l0.gammaMax, 1e-12);
  // linear scaling in the weak regime
  const weak2 = s.equivalentLinearResponse({ layers: mk(), bedrock: rock, acc: acc.map((a) => 0.002 * a), dt });
  close(weak2.pgaSurface / weak.pgaSurface, 2, 5e-3);
});

test('response spectrum (Nigam-Jennings): undamped step input doubles the static response; rigid oscillator returns the PGA', () => {
  const dt = 0.01, n = 2000;
  const acc = new Float64Array(n).fill(2.0);
  acc[0] = 2.0;
  const r = s.responseSpectrum({ acc, dt, periods: [0.5, 1.0], damping: 0 });
  // suddenly applied constant acceleration a: u_max = 2 a / w^2 -> PSA = 2 a
  close(r.PSA[0], 4.0, 1e-9);
  close(r.PSA[1], 4.0, 1e-9);
  const sine = new Float64Array(n);
  for (let i = 0; i < n; i++) sine[i] = 3 * Math.sin(2 * Math.PI * 1.0 * i * dt);
  const stiff = s.responseSpectrum({ acc: sine, dt, periods: [0.01], damping: 0.05 });
  close(stiff.SaAbs[0], 3, 2e-3);
});

test('response spectrum: steady-state resonance of a damped oscillator approaches a0 / (2 xi); time-step independence for a piecewise-linear record; impulse response', () => {
  const dt = 0.005, T = 0.8, xi = 0.05, a0 = 1.5;
  const n = 40000;
  const acc = new Float64Array(n);
  for (let i = 0; i < n; i++) acc[i] = a0 * Math.sin(2 * Math.PI * i * dt / T);
  const r = s.responseSpectrum({ acc, dt, periods: [T], damping: xi });
  // steady-state relative displacement amplitude at resonance: a0 / (2 xi w^2)
  close(r.PSA[0], a0 / (2 * xi), 2e-3);
  // triangular pulse then free vibration: the recurrence is exact for piecewise-linear input, so two time
  // steps give identical spectra when the record is sampled consistently
  const pulse = (h) => { const m = Math.round(2 / h); const p = new Float64Array(m + 1); for (let i = 0; i <= m; i++) { const t = i * h; p[i] = t <= 0.2 ? t / 0.2 : t <= 0.4 ? (0.4 - t) / 0.2 : 0; } return p; };
  const c1 = s.responseSpectrum({ acc: pulse(0.02), dt: 0.02, periods: [1, 2], damping: 0.05, freeVibrationCycles: 3 });
  const c2 = s.responseSpectrum({ acc: pulse(0.001), dt: 0.001, periods: [1, 2], damping: 0.05, freeVibrationCycles: 3 });
  // the coarse run samples the same piecewise-linear record, so only the sampling of the maximum differs
  // (at most 1 - cos(pi dt/T) = 0.2 %)
  for (let i = 0; i < 2; i++) assert.ok(Math.abs(c1.Sd[i] - c2.Sd[i]) / c2.Sd[i] < 3e-3, `T index ${i}: ${c1.Sd[i]} vs ${c2.Sd[i]}`);
  close(c2.PSV[1], Math.PI * c2.Sd[1], 1e-12);
  // undamped oscillator after a short impulse: Sd = I / w with I the impulse area (0.2 here), for T >> 0.4 s
  const long = s.responseSpectrum({ acc: pulse(0.001), dt: 0.001, periods: [40], damping: 0, freeVibrationCycles: 1 });
  close(long.Sd[0], 0.2 / (2 * Math.PI / 40), 1e-3);
});

test('Arias intensity and significant duration of a constant-amplitude sine: Ia = (pi/2g) a0^2 t/2, D5-95 = 0.9 t', () => {
  const dt = 0.001, dur = 20, a0 = 2;
  const n = Math.round(dur / dt) + 1;
  const acc = new Float64Array(n);
  for (let i = 0; i < n; i++) acc[i] = a0 * Math.sin(2 * Math.PI * 2 * i * dt);
  const r = s.ariasIntensity({ acc, dt });
  close(r.arias, (Math.PI / (2 * 9.80665)) * a0 * a0 * dur / 2, 1e-5);
  close(r.significantDuration, 0.9 * dur, 5e-3);
  close(r.pga, a0, 1e-6);
  close(r.husid[n - 1], 1, 1e-12);
});
