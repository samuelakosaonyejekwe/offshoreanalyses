// Coupled vessel–mooring response in the horizontal plane: time-domain slow-drift motion on the nonlinear mooring
// restoring curve under Newman's difference-frequency wave force, with first-order (wave-frequency) motion
// superposed, and the resulting line-tension history.
import { waveComponents } from '../../engine/hydro/timedomain.js';
import { newmanDriftSpectrum, slowDriftResponse } from '../../engine/hydro/drift.js';
import { boxBargeRAO, responseStatistics as raoStats } from '../../engine/hydro/floater.js';
import { jonswap } from '../../engine/metocean/spectra.js';
import { naturalPeriod } from '../../engine/lines/mooringsystem.js';
import { DEG, G, need, run, round, fmt, adv, check, worst, sev, thin, linspace, interp1, systemInputs, buildSystem, fastSystem, steelAdvice, LIMIT_PROPS, CONFIRM } from './_common.js';
import { envInputs, hullInputs, envForce, driftCoeff } from './_vessel.js';

const lin = { key: 'stiffMode', in: ['linear'] };

export default {
  id: 's10-coupled-vessel-mooring',
  suite: 10,
  title: 'Coupled vessel–mooring analysis',
  summary: 'Simulates the slow-drift and wave-frequency motion of the moored vessel in the design sea and returns offsets, six-degree-of-freedom motion statistics and the tension history of the critical line.',
  covers: [3],
  fidelity: 'reduced-order',
  standards: ['DNV-OS-E301 (2018) Ch.2 Sec.2 [2.2]–[2.4]', 'API RP 2SK (3rd ed., 2005) §5.4–5.5', 'DNV-RP-C205 (2019) §7.4 (slowly varying drift forces)', 'ISO 19901-7 (2013) §8'],
  models: ['Newtonian equation of motion along the load direction with nonlinear mooring restoring force', 'Second-order slowly varying drift force by Newman\'s approximation of the quadratic transfer function', 'First-order wave-frequency motion from response amplitude operators (Froude–Krylov surge/sway; heave, pitch and roll of a box hull from a three-dimensional panel solution)', 'JONSWAP random sea with seeded phases', 'Frequency-domain slow-drift response for cross-check', 'Explicit RK4 time integration'],
  limits: `One translational degree of freedom along the environment direction is integrated in time; wave-frequency motions are superposed linearly and the mooring acts quasi-statically (no line dynamics, no current or wave load on lines). The slow-drift equation uses a constant added mass and linearised damping: frequency-dependent radiation (Cummins convolution), full QTFs and yaw/sway coupling are outside this module and should come from a diffraction analysis for design. The wave-frequency heave, pitch and roll in the motion table are those of a rectangular box with the hull's length, beam and draught (panel solution, deep water, zero speed; roll with a gyradius of 0.4 B including added inertia and 5 % of critical added damping): within 5 % of an independent panel code for a box in waves longer than the beam, but not the motions of the actual hull form; a peak period below the wave whose length equals the beam carries a warning. ${LIMIT_PROPS}`,
  publishes: ['mooring.lfSigma', 'mooring.wfSigma', 'mooring.TnSurge'],
  searchTerms: 'coupled vessel mooring slow drift time domain Newman approximation',
  inputs: [
    ...systemInputs(),
    ...hullInputs(),
    { key: 'dispT', label: 'Displacement', unit: 't', type: 'number', default: 50000, min: 500, max: 600000, group: 'Vessel', role: 'geometry', bind: 'vessel.disp', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'CaSurge', label: 'Horizontal added-mass coefficient', unit: '-', type: 'number', default: 0.15, min: 0, max: 2, group: 'Vessel', role: 'load', help: 'Low-frequency added mass / displacement: ≈ 0.05–0.15 surge of a ship, 0.6–1.0 sway; 0.1–0.4 for semi-submersibles.' },
    { key: 'GM', label: 'Transverse metacentric height', unit: 'm', type: 'number', default: 4, min: 0.2, max: 40, group: 'Vessel', role: 'geometry', bind: 'vessel.GM' },
    { key: 'zetaLF', label: 'Low-frequency damping', unit: '% critical', type: 'number', default: 15, min: 1, max: 80, group: 'Vessel', role: 'load', uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'Viscous hull damping, wave-drift damping and mooring-line drag damping together; 10–30 % is common for spread moorings. It controls the slow-drift amplitude: calibrate against model tests or measured motions.' },
    ...envInputs(),
    { key: 'envDirDeg', label: 'Direction the environment travels towards', unit: 'deg', type: 'number', default: 0, min: -180, max: 360, group: 'Environment (design storm)', role: 'bc', help: 'Wind, waves and current collinear. Boundary condition on the vessel: environmental forcing along this direction.' },
    { key: 'stiffMode', label: 'Mooring restoring force', type: 'select', default: 'system', options: [{ value: 'system', label: 'Nonlinear, from the mooring system' }, { value: 'linear', label: 'Linear stiffness (entered)' }], group: 'Mooring restoring', role: 'bc' },
    { key: 'Klin', label: 'Linear mooring stiffness', unit: 'kN/m', type: 'number', default: 200, min: 1, max: 20000, group: 'Mooring restoring', role: 'bc', bind: 'mooring.Kxx', showIf: lin },
    { key: 'startOffset', label: 'Initial offset from the mean position', unit: 'm', type: 'number', default: 0, min: -60, max: 60, group: 'Initial condition', role: 'ic',
      help: '0 starts the vessel at its mean equilibrium under the steady and mean drift forces (a pre-existing environmental response). A non-zero value starts a free-decay on top of the storm.' },
    { key: 'startVel', label: 'Initial velocity', unit: 'm/s', type: 'number', default: 0, min: -3, max: 3, group: 'Initial condition', role: 'ic' },
    { key: 'duration', label: 'Simulated duration', unit: 's', type: 'number', default: 1800, min: 300, max: 10800, group: 'Numerical', role: 'control', help: 'Longer runs give more stable slow-drift statistics; 3 hours (10 800 s) is the design-storm duration but costs proportionally more time.' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 0.5, min: 0.05, max: 2, group: 'Numerical', role: 'control', help: 'Must resolve the wave period (about Tp/20 or better).' },
    { key: 'seed', label: 'Random seed of the sea state', type: 'int', default: 1, min: 1, max: 99999, group: 'Numerical', role: 'control', help: 'Different seeds are different realisations of the same sea; extremes vary between seeds.' },
    { key: 'sfIntact', label: 'Required safety factor on tension', unit: '-', type: 'number', default: 2.0, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: `Quasi-static line model: API RP 2SK intact factor 2.0. ${CONFIRM}` },
    { key: 'offsetLimitPct', label: 'Allowable offset', unit: '% depth', type: 'number', default: 8, min: 0.5, max: 60, group: 'Acceptance criteria', role: 'criteria' },
  ],
  meshTime: { key: 'dt', label: 'Time step', levels: [1.0, 0.5, 0.25], qoi: 'maxOffset', isSize: true },
  compute(inp) {
    const b = buildSystem(inp);
    const th = inp.envDirDeg * DEG, cd = Math.cos(th), sd = Math.sin(th);
    const M = inp.dispT * 1e3 * (1 + inp.CaSurge);                       // t → kg, with added mass
    const steady = envForce({ ...inp, Hs: 0 }, { windDir: inp.envDirDeg, waveDir: inp.envDirDeg, curDir: inp.envDirDeg, headingDeg: 0, rhoW: inp.rhoW });
    const Fs = steady.Fx * cd + steady.Fy * sd;
    // --- sea state: geometric frequency grid (no exact repetition of the record)
    const waves = inp.Hs > 1e-6;
    const nF = 100, wp = 2 * Math.PI / inp.Tp;
    const omega = Array.from({ length: nF }, (_, i) => 0.5 * wp * (7 ** (i / (nF - 1))));
    const S = omega.map((w) => (waves ? jonswap(w, { Hs: inp.Hs, Tp: inp.Tp, gamma: inp.gammaJ }) : 0));
    const comp = waveComponents({ omega, S, seed: inp.seed });
    const width = inp.Bv * Math.abs(cd) + inp.Lv * Math.abs(sd), Lalong = inp.Lv * Math.abs(cd) + inp.Bv * Math.abs(sd);
    const D = driftCoeff(omega, { R2: inp.R2, width, draught: inp.draught, rho: inp.rhoW });
    const FdMean = comp.amp.reduce((a, v, i) => a + v * v * D[i], 0);
    // Froude–Krylov surge/sway RAO of a box of length Lalong and draught T along the wave direction (diffraction neglected)
    const raoX = omega.map((w) => { const k = w * w / G, x = 0.5 * k * Lalong, kT = k * inp.draught; return Math.abs(x > 1e-6 ? Math.sin(x) / x : 1) * (1 - Math.exp(-kT)) / kT / (1 + inp.CaSurge); });
    // --- restoring curve along the environment direction
    let fs = null, calm = null, offs, Rtab, Ttab = null, govLine = 0;
    const Ftot = Fs + FdMean;
    if (inp.stiffMode === 'system') {
      fs = fastSystem(b.system, { maxOffset: Math.max(0.15 * b.depth, 30) });
      calm = fs.solve([0, 0, 0, 0, 0, 0]);
      const mean = fs.solve([Ftot * cd, Ftot * sd, 0, 0, 0, 0], [calm.x, calm.y]);
      const d0 = (mean.x - calm.x) * cd + (mean.y - calm.y) * sd;
      const span = Math.max(3 * d0, 0.1 * b.depth, 20);
      offs = linspace(-0.5 * span, 1.6 * span, 70);
      let cv;
      try { cv = fs.curve(th, offs, [calm.x, calm.y]); } catch (e) {
        offs = linspace(-0.3 * span, Math.max(1.15 * d0, 0.5 * span), 70);
        cv = run(() => fs.curve(th, offs, [calm.x, calm.y]), 'Restoring curve');
      }
      Rtab = cv.restoring; govLine = mean.maxTensionLine;
      Ttab = cv.lineTensions.map((r) => r[govLine]);
    } else {
      const K = inp.Klin * 1e3, d0 = Ftot / K, span = Math.max(3 * d0, 20);
      offs = linspace(-span, 2 * span, 40); Rtab = offs.map((d) => K * d); Ttab = Rtab.map((v) => Math.abs(v));
    }
    const R = (d) => {
      // linear extrapolation beyond the table with the end slope (stiffening is under-estimated there: flagged below)
      const n = offs.length;
      if (d > offs[n - 1]) return Rtab[n - 1] + (Rtab[n - 1] - Rtab[n - 2]) / (offs[n - 1] - offs[n - 2]) * (d - offs[n - 1]);
      if (d < offs[0]) return Rtab[0] + (Rtab[1] - Rtab[0]) / (offs[1] - offs[0]) * (d - offs[0]);
      return interp1(offs, Rtab, d);
    };
    // mean position: R(d) = Ftot by bisection on the monotone table
    need(Ftot <= Rtab[Rtab.length - 1] * 1.5, 'The mean environmental force exceeds the restoring capacity of the mooring within its elastic range. Reduce the environment or strengthen the system.');
    let lo = offs[0], hi = offs[offs.length - 1];
    for (let i = 0; i < 60; i++) { const mid = 0.5 * (lo + hi); if (R(mid) < Ftot) lo = mid; else hi = mid; }
    const dMean = 0.5 * (lo + hi);
    const K0 = (R(dMean + 0.25) - R(dMean - 0.25)) / 0.5;
    need(K0 > 0, 'The mooring has no positive stiffness at the mean position.');
    const Tn = naturalPeriod({ mass: M, stiffness: K0 });
    const c = 2 * (inp.zetaLF / 100) * Math.sqrt(K0 * M);
    // --- time integration (RK4) of M x'' + c x' + R(x) = Fs + Fd(t); forces sampled at half steps
    const dt = inp.dt, N = Math.round(inp.duration / dt);
    need(N >= 100 && N <= 400000, 'Duration / time step must give between 100 and 400 000 steps.');
    const sq = Array.from(D, (v) => Math.sqrt(v));
    const Fd = new Float64Array(2 * N + 1), xw = new Float64Array(N + 1);
    for (let j = 0; j <= 2 * N; j++) {
      const t = 0.5 * j * dt; let re = 0, im = 0, x = 0;
      for (let i = 0; i < nF; i++) { const a = comp.amp[i]; if (a === 0) continue; const ph = omega[i] * t + comp.phase[i], cs = Math.cos(ph); re += a * sq[i] * cs; im += a * sq[i] * Math.sin(ph); x += a * raoX[i] * cs; }
      Fd[j] = re * re + im * im;
      if (j % 2 === 0) xw[j / 2] = x;
    }
    const xl = new Float64Array(N + 1), tt = new Float64Array(N + 1);
    let x = dMean + inp.startOffset, v = inp.startVel;
    const acc = (xx, vv, F) => (Fs + F - c * vv - R(xx)) / M;
    for (let n = 0; n <= N; n++) {
      xl[n] = x; tt[n] = n * dt;
      if (n === N) break;
      const f0 = Fd[2 * n], f1 = Fd[2 * n + 1], f2 = Fd[2 * n + 2];
      const k1v = acc(x, v, f0), k1x = v;
      const k2v = acc(x + 0.5 * dt * k1x, v + 0.5 * dt * k1v, f1), k2x = v + 0.5 * dt * k1v;
      const k3v = acc(x + 0.5 * dt * k2x, v + 0.5 * dt * k2v, f1), k3x = v + 0.5 * dt * k2v;
      const k4v = acc(x + dt * k3x, v + dt * k3v, f2), k4x = v + dt * k3v;
      x += dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x); v += dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v);
      if (!Number.isFinite(x)) throw new Error('The time integration diverged: reduce the time step.');
    }
    // --- statistics after the start-up window
    const i0 = inp.startOffset === 0 && inp.startVel === 0 ? Math.min(Math.round(0.5 * Tn / dt), Math.floor(0.1 * N)) : Math.min(Math.round(1.5 * Tn / dt), Math.floor(0.4 * N));
    let sx = 0, sxx = 0, sw = 0, xmax = -Infinity, xmin = Infinity, Tmax = -Infinity, Tmin = Infinity, sT = 0, sTT = 0, cnt = 0, beyond = false;
    const Tline = new Float64Array(N + 1);
    for (let n = 0; n <= N; n++) {
      const d = xl[n] + xw[n];
      if (d > offs[offs.length - 1]) beyond = true;
      Tline[n] = interp1(offs, Ttab, d);
      if (n < i0) continue;
      cnt++; sx += xl[n]; sxx += xl[n] * xl[n]; sw += xw[n] * xw[n];
      xmax = Math.max(xmax, d); xmin = Math.min(xmin, d);
      Tmax = Math.max(Tmax, Tline[n]); Tmin = Math.min(Tmin, Tline[n]); sT += Tline[n]; sTT += Tline[n] * Tline[n];
    }
    const meanOff = sx / cnt, lfSig = Math.sqrt(Math.max(0, sxx / cnt - meanOff * meanOff)), wfSig = Math.sqrt(sw / cnt);
    const Tmean = sT / cnt, Tstd = Math.sqrt(Math.max(0, sTT / cnt - Tmean * Tmean));
    // --- frequency-domain cross-check of the slow-drift standard deviation (linearised at the mean position)
    let fd = { sigma: 0 };
    if (waves) {
      const mu = linspace(0.2 * 2 * Math.PI / Tn, Math.min(5 * 2 * Math.PI / Tn, 0.6 * wp), 60);
      const og = linspace(omega[0], omega[nF - 1], 240), Sg = og.map((w) => jonswap(w, { Hs: inp.Hs, Tp: inp.Tp, gamma: inp.gammaJ }));
      const SF = newmanDriftSpectrum({ omega: og, S: Sg, D: driftCoeff(og, { R2: inp.R2, width, draught: inp.draught, rho: inp.rhoW }), mu });
      fd = slowDriftResponse({ mass: M, stiffness: K0, damping: c, mu, SF, meanForce: Ftot, duration: inp.duration });
    }
    // --- six-degree-of-freedom wave-frequency statistics from the box-hull RAOs
    const dof = [];
    let hullRange = null;
    if (waves) {
      const st = (amp) => raoStats({ omega, amp, S, duration: 10800 });
      const hull = boxBargeRAO({ L: inp.Lv, B: inp.Bv, T: inp.draught, omega, headingDeg: inp.envDirDeg, GMt: inp.GM, rho: inp.rhoW });
      hullRange = hull.validPeriodMin;
      const sAlong = st(raoX);
      dof.push(['Surge', sAlong.significant * Math.abs(cd), sAlong.mpm * Math.abs(cd), 'm'], ['Sway', sAlong.significant * Math.abs(sd), sAlong.mpm * Math.abs(sd), 'm']);
      const hv = st(hull.heave), pt = st(hull.pitch), rl = hull.roll ? st(hull.roll) : { significant: 0, mpm: 0 };
      dof.push(['Heave', hv.significant, hv.mpm, 'm'], ['Roll', rl.significant / DEG, rl.mpm / DEG, 'deg'], ['Pitch', pt.significant / DEG, pt.mpm / DEG, 'deg'], ['Yaw', 0, 0, 'deg']);
    } else for (const nme of ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']) dof.push([nme, 0, 0, nme.length > 4 || nme === 'Sway' ? 'm' : 'deg']);
    const MBL = b.MBL;
    const lineBased = inp.stiffMode === 'system';
    const checks = [
      check('Maximum offset', xmax, inp.offsetLimitPct / 100 * b.depth, 'm', 'API RP 2SK (3rd ed.) §6.2; project riser limit',
        'Guards against exceeding riser or umbilical excursion limits during the simulated storm.'),
    ];
    if (lineBased) checks.push(check('Maximum tension of the critical line', Tmax * inp.sfIntact / 1e3, MBL / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3 (quasi-static line model)',
      'Guards against breaking the windward line at the largest combined slow-drift and wave-frequency excursion.'));
    const g = worst(checks);
    const [pt, pl, ptot, pT, pF] = thin([tt, xl, Array.from(xl, (vv, i) => vv + xw[i]), Tline, Array.from({ length: N + 1 }, (_, i) => Fd[2 * i])], 900);
    const per = omega.map((w) => 2 * Math.PI / w);
    return {
      kpis: [
        { key: 'Tn', label: 'Surge natural period on the mooring', value: round(Tn, 5), unit: 's', explain: '2π√((M + A)/K) at the mean position. Far above wave periods, so the vessel is excited by slowly varying drift forces, not directly by waves.' },
        { key: 'meanOffset', label: 'Mean offset', value: round(meanOff, 5), unit: 'm', explain: 'Time-average position along the environment direction under wind, current and mean wave drift.' },
        { key: 'lfSigma', label: 'Low-frequency motion, standard deviation', value: round(lfSig, 4), unit: 'm', explain: 'Slow-drift oscillation about the mean position from the time simulation; it scales inversely with the square root of damping and stiffness.' },
        { key: 'lfSigmaFD', label: 'Low-frequency standard deviation, frequency domain', value: round(fd.sigma, 4), unit: 'm', explain: 'Same quantity from the linearised spectral solution with Newman\'s drift-force spectrum. Agreement within sampling scatter confirms the time-domain run.' },
        { key: 'wfSigma', label: 'Wave-frequency motion, standard deviation', value: round(wfSig, 4), unit: 'm', explain: 'First-order motion along the environment direction from the Froude–Krylov RAO; not influenced by the mooring.' },
        { key: 'maxOffset', label: 'Maximum offset in the simulation', value: round(xmax, 5), unit: 'm', limit: round(inp.offsetLimitPct / 100 * b.depth, 5), sense: 'max', explain: 'Largest combined excursion (mean + slow drift + wave frequency) reached in this realisation of the storm.' },
        { key: 'Tmax', label: lineBased ? 'Maximum tension, critical line' : 'Maximum restoring force', value: round(Tmax / 1e3, 5), unit: 'kN', explain: 'Largest quasi-static tension of the windward line during the run (with a linear stiffness: the largest spring force instead).' },
        { key: 'Tmin', label: lineBased ? 'Minimum tension, critical line' : 'Minimum restoring force', value: round(Tmin / 1e3, 5), unit: 'kN', explain: 'Smallest tension of the same line; with the maximum it gives the tension range the line cycles through.' },
        { key: 'Tmean', label: lineBased ? 'Mean tension, critical line' : 'Mean restoring force', value: round(Tmean / 1e3, 5), unit: 'kN', explain: 'Time-average tension, the mean-load part in the DNV partial-factor format.' },
        { key: 'Tstd', label: 'Tension standard deviation', value: round(Tstd / 1e3, 5), unit: 'kN', explain: 'Cyclic tension content (low- plus wave-frequency); the input the fatigue module needs for this sea state.' },
        { key: 'FdriftMean', label: 'Mean wave-drift force', value: round(FdMean / 1e3, 5), unit: 'kN', explain: 'Average of the slowly varying second-order wave force; its fluctuation about this mean is what drives the slow-drift motion.' },
        { key: 'Fsteady', label: 'Steady wind and current force', value: round(Fs / 1e3, 5), unit: 'kN', explain: 'Constant part of the environmental load along the environment direction.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Vessel offset time history', x: { label: 'Time', unit: 's' }, y: { label: 'Offset along environment direction', unit: 'm' },
          series: [{ name: 'Total (LF + WF)', x: pt.map((q) => round(q)), y: ptot.map((q) => round(q)) }, { name: 'Low-frequency', x: pt.map((q) => round(q)), y: pl.map((q) => round(q)) }, { name: 'Allowable', x: [0, inp.duration], y: [inp.offsetLimitPct / 100 * b.depth, inp.offsetLimitPct / 100 * b.depth], dash: true }],
          explain: 'Slow-drift excursions with the faster wave-frequency motion riding on them. Look for groups of large waves pushing the vessel towards the allowable line.' },
        { type: 'line', title: lineBased ? 'Critical-line tension time history' : 'Restoring force time history', x: { label: 'Time', unit: 's' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: lineBased ? b.system.lines[govLine].name : 'Spring force', x: pt.map((q) => round(q)), y: pT.map((q) => round(q / 1e3)) }, ...(lineBased ? [{ name: 'MBL / safety factor', x: [0, inp.duration], y: [MBL / inp.sfIntact / 1e3, MBL / inp.sfIntact / 1e3], dash: true }] : [])],
          explain: 'Quasi-static tension of the windward line as the vessel moves. Peaks are sharper than the offset peaks because the catenary stiffens with offset.' },
        { type: 'line', title: 'Slowly varying wave-drift force', x: { label: 'Time', unit: 's' }, y: { label: 'Force', unit: 'kN' },
          series: [{ name: 'Drift force (Newman)', x: pt.map((q) => round(q)), y: pF.map((q) => round(q / 1e3)) }, { name: 'Mean', x: [0, inp.duration], y: [FdMean / 1e3, FdMean / 1e3], dash: true }],
          explain: 'Second-order wave force following the square of the wave envelope: always positive, large under wave groups. Its energy near the surge natural period is what excites slow drift.' },
        { type: 'line', title: 'Mooring restoring curve', x: { label: 'Offset', unit: 'm' }, y: { label: 'Force', unit: 'kN' },
          series: [{ name: 'Restoring force', x: offs.map((q) => round(q)), y: Rtab.map((q) => round(q / 1e3)) }, { name: 'Mean environmental force', x: [offs[0], offs[offs.length - 1]], y: [Ftot / 1e3, Ftot / 1e3], dash: true }],
          explain: 'Nonlinear force–offset characteristic used in the simulation; the crossing with the dashed line is the mean position and the local slope is the stiffness that sets the natural period.' },
        { type: 'line', title: 'Wave-frequency response amplitude operators', x: { label: 'Wave period', unit: 's' }, y: { label: 'Surge RAO (m/m)', unit: '-' },
          series: [{ name: 'Surge/sway along waves', x: per.map((q) => round(q)), y: raoX.map((q) => round(q)) }],
          explain: 'First-order motion per metre of wave amplitude against wave period (Froude–Krylov estimate). It tends to 1 for long waves, where the hull follows the water particles.' },
      ],
      tables: [
        { title: 'Vessel wave-frequency response in six degrees of freedom (3-hour storm)', columns: [{ key: 'dof', label: 'Motion' }, { key: 'sig', label: 'Significant amplitude' }, { key: 'mpm', label: 'Most probable maximum' }, { key: 'u', label: 'Unit' }],
          rows: dof.map(([d, s, m, u]) => ({ dof: d, sig: round(s, 4), mpm: round(m, 4), u })),
          explain: 'Significant (2σ) and most probable 3-hour maximum amplitude of each first-order motion from box-hull RAOs (panel solution). Yaw is not excited in this collinear, heading-held model and is reported as zero.' },
        { title: 'Offset and tension statistics', columns: [{ key: 'q', label: 'Quantity' }, { key: 'mean', label: 'Mean' }, { key: 'sd', label: 'Std. dev.' }, { key: 'max', label: 'Maximum' }, { key: 'min', label: 'Minimum' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Offset', mean: round(meanOff), sd: round(Math.hypot(lfSig, wfSig)), max: round(xmax), min: round(xmin), u: 'm' }, { q: lineBased ? 'Critical-line tension' : 'Restoring force', mean: round(Tmean / 1e3), sd: round(Tstd / 1e3), max: round(Tmax / 1e3), min: round(Tmin / 1e3), u: 'kN' }],
          explain: 'Summary statistics of the simulated run after the start-up window; feed the standard deviations and periods to the tension-and-offset and fatigue modules.' },
      ],
      notes: [`Natural period ${fmt(Tn)} s with stiffness ${fmt(K0 / 1e3)} kN/m at the mean position and ${inp.zetaLF} % damping; ${fmt(inp.duration / Tn)} slow-drift cycles simulated.`,
        `Time-domain low-frequency σ = ${fmt(lfSig)} m against ${fmt(fd.sigma)} m from the linearised frequency-domain solution.`,
        beyond ? 'The vessel moved beyond the tabulated restoring curve; tension there is capped at the last tabulated value — treat the maximum as a lower bound and reduce the load or strengthen the system.' : 'All excursions stayed within the tabulated restoring curve.'],
      advice: [
        ...(hullRange && inp.Tp < hullRange.heave ? [adv('engineering', 'warn',
          `Peak period ${fmt(inp.Tp)} s is below the range of the box-hull motion model, whose panel solution is carried down to the wave whose length equals the beam (period ${fmt(hullRange.heave)} s). Much of the energy of this sea lies in shorter waves, where the heave, pitch and roll in the motion table use the coefficients of that last period and are not verified.`,
          'Take the wave-frequency motions from the diffraction analysis of the vessel for this sea state. The offset and line tensions of this module do not depend on them.',
          'DNV-RP-C205 (2019) §7.3 (wave-frequency motions by radiation–diffraction analysis); comparison of the box panel model with Capytaine 3.0.0 on barges of 20 × 8 × 2 m, 100 × 20 × 5 m and 80 × 24 × 3 m (cross-solver verification of this application).')] : []),
        adv('engineering', sev(g.uc), `Maximum offset ${fmt(xmax)} m (${fmt(100 * xmax / b.depth)} % of depth): mean ${fmt(meanOff)} m, slow-drift σ ${fmt(lfSig)} m, wave-frequency σ ${fmt(wfSig)} m.${lineBased ? ` Critical-line tension peaks at ${fmt(Tmax / 1e3)} kN (${fmt(100 * Tmax / MBL)} % MBL).` : ''}`,
          g.uc > 1 ? `Reduce the excursion by ${fmt(Math.max(0, xmax - inp.offsetLimitPct / 100 * b.depth))} m: more pretension shortens both the mean offset and the slow drift (σ ∝ 1/√K); confirm the ${inp.zetaLF} % damping, which is the most uncertain input.`
            : `Use σLF = ${fmt(lfSig)} m, σWF = ${fmt(wfSig)} m and Tn = ${fmt(Tn)} s in the tension-and-offset module; repeat with several seeds and take the mean of the maxima for a 3-hour extreme.`,
          'DNV-OS-E301 (2018) Ch.2 Sec.2 [2.2]–[2.4]; API RP 2SK (3rd ed., 2005) §5.4; DNV-RP-C205 (2019) §7.4.'),
        adv('engineering', 'info', 'Hull motions come from the RAOs of a box with the main dimensions of the hull and a constant reflection coefficient, and the lines are quasi-static.',
          'Before design use, replace the RAOs and drift coefficients with diffraction-analysis output for the actual hull, and apply the dynamic amplification from the line-dynamics module to the wave-frequency tension.',
          'DNV-OS-E301 (2018) Ch.2 Sec.2 [2.1] (analysis methods); ISO 19901-7 (2013) §8.6.'),
        steelAdvice(b, lineBased ? Tmax * inp.sfIntact / MBL : 0.5, 'coupled-analysis tension'),
      ],
      publish: { 'mooring.lfSigma': lfSig, 'mooring.wfSigma': wfSig, 'mooring.TnSurge': Tn },
      quantities: { steel_t: round(b.steel_t) },
    };
  },
  verify: [
    { // M = 50 000 t × (1 + 0.1) = 5.5e7 kg, K = 200 kN/m → Tn = 2π·sqrt(5.5e7 / 2e5) = 2π·16.583 = 104.20 s.
      // No waves or current; wind 40 m/s on 1500 m² with Cd = 1 → 1471.2 kN → mean offset 1471.2 / 200 = 7.356 m, no oscillation.
      name: 'Linear mooring, steady load only: natural period and static offset',
      source: 'Single-degree-of-freedom dynamics (hand calculation); DNV-RP-C205 §5.3 wind force',
      inputs: { stiffMode: 'linear', Klin: 200, dispT: 50000, CaSurge: 0.1, Hs: 0, Uc: 0, Uw: 40, Afront: 1500, CdWind: 1.0, envDirDeg: 0 },
      expect: [{ kpi: 'Tn', value: 104.20, tol: 0.001 }, { kpi: 'meanOffset', value: 7.356, tol: 0.002 }, { kpi: 'maxOffset', value: 7.356, tol: 0.002 }, { kpi: 'lfSigma', value: 0, abs: 1e-3 }] },
  ],
};
