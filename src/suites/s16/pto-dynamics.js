// Power-take-off dynamic response: time-domain simulation of a heaving absorber with the Cummins impulse-response
// equation and a linear, Coulomb or hydraulic (piston – accumulator – motor) power take-off, with end stops,
// latching and declutching control; PTO force, pressure and converted power histories.
import { dispersion, waveSpectrum, cylinderHeave, radiationKernel, infiniteAddedMass, simulateCummins } from '../../engine/aero/wec.js';
import { makeRng } from '../../engine/math/rng.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, thin, G, LIMIT_WEC, LIMIT_CODE, CONFIRM, energyCo2Text } from './_common.js';

const hyd = { key: 'ptoType', in: ['hydraulic'] };
const LPM_BAR = 1e-3 / 60 / 1e5;   // 1 L/min/bar in m³/s/Pa

export default {
  id: 's16-pto-dynamics', suite: 16, title: 'Power-take-off dynamic response (time domain)',
  summary: 'Time histories of float motion, PTO force, hydraulic pressure and converted power for a wave-energy converter with a linear, Coulomb or hydraulic power take-off, end stops and latching control.',
  covers: [12], fidelity: 'numerical',
  standards: ['IEC TS 62600-2:2019 (design requirements for marine energy systems)', 'IEC TS 62600-100:2012 (power performance)', 'IEC TS 62600-103:2018 (early-stage scale testing — basis for calibration)', 'Cummins (1962) / Ogilvie (1964) impulse-response formulation', 'Falnes (2002) Ocean Waves and Oscillating Systems'],
  models: ['Cummins impulse-response equation with radiation convolution', 'Hydrodynamic added mass and radiation damping (Ogilvie relations)', 'Power-take-off force models: linear damper/spring, Coulomb, hydraulic', 'Hydraulic PTO dynamics: rectified piston, polytropic gas accumulator, variable-pressure motor', 'Latching control (Budal & Falnes)', 'Declutching control', 'End-stop spring and quadratic viscous drag', 'Random-phase synthesis of irregular wave excitation', 'Fourth-order Runge–Kutta integration'],
  limits: `${LIMIT_WEC} Single degree of freedom (heave); the excitation force is linear (no nonlinear Froude–Krylov correction for the varying waterplane); hydraulic valve losses, oil compressibility and the low-pressure circuit are neglected; the generator is represented by a motor efficiency (see the generator and converter module for the electrical machine). ${LIMIT_CODE}`,
  publishes: ['wec.ptoForceMax', 'wec.meanElectricalPower', 'wec.peakToAverage'],
  searchTerms: 'wave energy hydraulic power take-off Cummins equation latching control time domain',
  inputs: [
    { key: 'radius', label: 'Float radius', unit: 'm', type: 'number', default: 5, min: 0.5, max: 20, group: 'Device', role: 'geometry' },
    { key: 'draft', label: 'Float draft', unit: 'm', type: 'number', default: 6, min: 0.3, max: 40, group: 'Device', role: 'geometry' },
    { key: 'massRatio', label: 'Body mass / displaced mass', unit: '-', type: 'number', default: 1.5, min: 0.5, max: 3, group: 'Device', role: 'ic' },
    { key: 'cdA', label: 'Viscous drag Cd × area', unit: 'm²', type: 'number', default: 40, min: 0, max: 2000, group: 'Device', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Quadratic drag ½ρ·CdA·v|v| on the float (plate area × Cd ≈ 0.5–1).' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 50, min: 10, max: 500, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'waveType', label: 'Incident waves', type: 'select', default: 'irregular', group: 'Sea state', role: 'bc', options: [{ value: 'irregular', label: 'Irregular (JONSWAP, seeded)' }, { value: 'regular', label: 'Regular wave (H = Hs, T = Tp)' }] },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 2.5, min: 0.1, max: 12, group: 'Sea state', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, bind: 'site.Hs1', help: 'Incident-wave boundary condition. Bound to the 1-year sea state when site data are pulled — reduce it to study normal operation.' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 9, min: 3, max: 20, group: 'Sea state', role: 'bc', bind: 'site.Tp1' },
    { key: 'gamma', label: 'Peak-enhancement factor', unit: '-', type: 'number', default: 2, min: 1, max: 7, group: 'Sea state', role: 'bc' },
    { key: 'seed', label: 'Random seed', type: 'int', default: 7, min: 1, max: 99999, group: 'Sea state', role: 'ic', help: 'Fixes the wave phases so that runs are reproducible.' },
    { key: 'x0', label: 'Initial heave displacement', unit: 'm', type: 'number', default: 0, min: -5, max: 5, group: 'Initial state', role: 'ic', help: 'Starting position relative to hydrostatic equilibrium; the float starts at rest.' },
    { key: 'ptoType', label: 'Power take-off', type: 'select', default: 'hydraulic', group: 'Power take-off', role: 'bc', options: [{ value: 'hydraulic', label: 'Hydraulic: piston, accumulator, motor' }, { value: 'linear', label: 'Linear damper (and spring)' }, { value: 'coulomb', label: 'Coulomb (constant force)' }], help: 'Mechanical loading applied by the PTO to the float.' },
    { key: 'Bpto', label: 'PTO damping', unit: 'kN·s/m', type: 'number', default: 500, min: 0, max: 20000, group: 'Power take-off', role: 'bc', showIf: { key: 'ptoType', in: ['linear'] }, design: true, help: 'The frequency-domain module reports the optimal passive damping at the peak period as a starting value.' },
    { key: 'Kpto', label: 'PTO stiffness (negative = reactive)', unit: 'kN/m', type: 'number', default: 0, min: -700, max: 5000, group: 'Power take-off', role: 'bc', showIf: { key: 'ptoType', in: ['linear'] } },
    { key: 'Fc', label: 'Coulomb force', unit: 'kN', type: 'number', default: 300, min: 1, max: 20000, group: 'Power take-off', role: 'bc', showIf: { key: 'ptoType', in: ['coulomb'] }, design: true },
    { key: 'Ap', label: 'Piston area', unit: 'cm²', type: 'number', default: 250, min: 10, max: 5000, group: 'Hydraulic circuit', role: 'geometry', showIf: hyd, design: true, help: 'Effective area of the cylinder(s); force = pressure difference × area.' },
    { key: 'p0', label: 'Accumulator pre-charge pressure', unit: 'bar', type: 'number', default: 80, min: 10, max: 350, group: 'Hydraulic circuit', role: 'ic', showIf: hyd, design: true, help: 'Initial gas pressure in the high-pressure accumulator.' },
    { key: 'V0', label: 'Accumulator gas volume at pre-charge', unit: 'm³', type: 'number', default: 1.5, min: 0.02, max: 50, group: 'Hydraulic circuit', role: 'geometry', showIf: hyd, design: true, help: 'A larger volume smooths the pressure and the electrical output.' },
    { key: 'pLow', label: 'Low-pressure side', unit: 'bar', type: 'number', default: 5, min: 0, max: 50, group: 'Hydraulic circuit', role: 'bc', showIf: hyd },
    { key: 'Cm', label: 'Motor flow coefficient', unit: 'L/min/bar', type: 'number', default: 12, min: 0.1, max: 500, group: 'Hydraulic circuit', role: 'control', showIf: hyd, design: true, help: 'Motor flow per unit pressure difference (displacement × speed control): sets the mean operating pressure.' },
    { key: 'etaM', label: 'Motor and generator efficiency', unit: '-', type: 'number', default: 0.85, min: 0.3, max: 1, group: 'Hydraulic circuit', role: 'operation', uncertain: { dist: 'normal', cov: 0.04 } },
    { key: 'pMax', label: 'Maximum working pressure', unit: 'bar', type: 'number', default: 280, min: 50, max: 420, group: 'Acceptance', role: 'criteria', showIf: hyd, help: 'Relief-valve setting / design pressure of the circuit.' },
    { key: 'stroke', label: 'End-stop stroke (amplitude)', unit: 'm', type: 'number', default: 3.5, min: 0.3, max: 20, group: 'Acceptance', role: 'bc', design: true },
    { key: 'kStop', label: 'End-stop stiffness', unit: 'MN/m', type: 'number', default: 20, min: 0.1, max: 1000, group: 'Acceptance', role: 'bc' },
    { key: 'forceLim', label: 'PTO force rating', unit: 'kN', type: 'number', default: 2500, min: 10, max: 1e5, group: 'Acceptance', role: 'criteria' },
    { key: 'latch', label: 'Latching duration (0 = off)', unit: 's', type: 'number', default: 0, min: 0, max: 5, group: 'Control', role: 'control', help: 'The float is held for this time at each velocity reversal; about (T_wave − T_natural)/2 is the classical setting.' },
    { key: 'declutch', label: 'Declutching velocity threshold (0 = off)', unit: 'm/s', type: 'number', default: 0, min: 0, max: 2, group: 'Control', role: 'control' },
    { key: 'duration', label: 'Simulated duration', unit: 's', type: 'number', default: 400, min: 60, max: 3600, group: 'Solver', role: 'control', help: 'The first 60 s (or a quarter of the run) are discarded as start-up transient. Use 1800 s or more for statistics to report.' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 0.05, min: 0.005, max: 0.25, group: 'Solver', role: 'control' },
  ],
  meshTime: { key: 'dt', label: 'Time step (s)', levels: [0.1, 0.05, 0.025], qoi: 'pAbs', isSize: true },
  compute(inp) {
    const rho = 1025, a = inp.radius, d = inp.draft;
    need(inp.depth > d + 2, 'The draft leaves no water under the float.');
    need(inp.duration / inp.dt <= 2e5, 'Too many time steps: increase the time step or shorten the run (limit 200 000 steps).');
    const disp = rho * Math.PI * a * a * d, mass = disp * inp.massRatio, Kh = rho * G * Math.PI * a * a;
    const geom = { radius: a, draft: d, depth: inp.depth, rho };
    const wp = 2 * Math.PI / inp.Tp, cp = cylinderHeave({ ...geom, omega: wp });
    // radiation kernel from B(ω) up to where it has decayed; A∞ from the Ogilvie relation at the peak frequency
    const wk = linspace(0.02, 8, 800), Bk = wk.map((om) => cylinderHeave({ ...geom, omega: om }).B);
    const tK = 25, kernel = radiationKernel({ omega: wk, B: Bk, dt: Math.min(inp.dt, 0.05), tMax: tK });
    const Ainf = Math.max(0, infiniteAddedMass({ A: cp.A, omega: wp, kernel }));
    // excitation force samples
    const n = Math.round(inp.duration / inp.dt), exc = new Float64Array(n + 2), eta = new Float64Array(n + 2);
    let Jwave;
    if (inp.waveType === 'regular') {
      const amp = inp.Hs / 2;
      for (let i = 0; i < n + 2; i++) { const ramp = Math.min(1, i * inp.dt / (2 * inp.Tp)); eta[i] = ramp * amp * Math.cos(wp * i * inp.dt); exc[i] = cp.Fe * eta[i]; }
      Jwave = 0.5 * rho * G * amp * amp * dispersion(wp, inp.depth).cg;
    } else {
      const nw = 160, om = linspace(0.25 * wp, Math.min(4 * wp, 4), nw), dw = om[1] - om[0];
      const S = waveSpectrum({ Hs: inp.Hs, Tp: inp.Tp, gamma: inp.gamma, omega: om });
      const rng = makeRng(inp.seed);
      // perturb each frequency inside its bin so that the series does not repeat within the run
      const comps = om.map((w0, i) => { const w = w0 + (rng.next() - 0.5) * dw; return { w, a: Math.sqrt(2 * S[i] * dw), ph: 2 * Math.PI * rng.next(), Fe: cylinderHeave({ ...geom, omega: w }).Fe, cg: dispersion(w, inp.depth).cg }; });
      for (let i = 0; i < n + 2; i++) {
        const t = i * inp.dt; let e = 0, f = 0;
        for (const c of comps) { const cs = c.a * Math.cos(c.w * t + c.ph); e += cs; f += c.Fe * cs; }
        eta[i] = e; exc[i] = f;
      }
      Jwave = comps.reduce((s, c) => s + 0.5 * rho * G * c.a * c.a * c.cg, 0);
    }
    const pto = inp.ptoType === 'linear' ? { type: 'linear', B: inp.Bpto * 1e3, K: inp.Kpto * 1e3 }
      : inp.ptoType === 'coulomb' ? { type: 'coulomb', force: inp.Fc * 1e3, vSmooth: 0.02 }
        : { type: 'hydraulic', pistonArea: inp.Ap * 1e-4, p0: inp.p0 * 1e5, V0: inp.V0, gamma: 1.4, pLow: inp.pLow * 1e5, motorCoeff: inp.Cm * LPM_BAR, motorEfficiency: inp.etaM, vSmooth: 0.02 };
    if (inp.ptoType === 'linear') need(Kh + inp.Kpto * 1e3 > 0.05 * Kh, `A PTO stiffness of ${inp.Kpto} kN/m cancels the hydrostatic stiffness (${fmt(Kh / 1e3)} kN/m): the float would have no restoring force. Keep the total stiffness positive.`);
    if (inp.ptoType === 'hydraulic') need(inp.p0 > inp.pLow, 'The accumulator pre-charge must exceed the low-pressure side.');
    const sim = run(() => simulateCummins({ mass, Ainf, stiffness: Kh, kernel, excitation: exc, dt: inp.dt, duration: inp.duration, x0: inp.x0, pto, dragCoeff: 0.5 * rho * inp.cdA,
      endStop: { stroke: inp.stroke, stiffness: inp.kStop * 1e6 }, ...(inp.latch > 0 ? { latching: { duration: inp.latch } } : {}), ...(inp.declutch > 0 ? { declutch: { threshold: inp.declutch } } : {}), kernelLength: tK }), 'Time-domain simulation');
    // statistics after the start-up transient
    const i0 = Math.round(Math.min(60, inp.duration / 4) / inp.dt), N = sim.t.length - i0;
    let eAbs = 0, eGen = 0, maxX = 0, maxV = 0, maxF = 0, pPeak = 0, pMaxBar = 0, pSum = 0, hits = 0, sx = 0;
    const pInst = new Float64Array(sim.t.length);
    for (let i = 0; i < sim.t.length; i++) {
      const pw = -sim.ptoForce[i] * sim.v[i];
      const gen = inp.ptoType === 'hydraulic' ? inp.etaM * pto.motorCoeff * (sim.pressure[i] - pto.pLow) ** 2 : inp.etaM * pw;
      pInst[i] = gen;
      if (i < i0) continue;
      eAbs += pw; eGen += gen; sx += sim.x[i] * sim.x[i];
      maxX = Math.max(maxX, Math.abs(sim.x[i])); maxV = Math.max(maxV, Math.abs(sim.v[i])); maxF = Math.max(maxF, Math.abs(sim.ptoForce[i])); pPeak = Math.max(pPeak, gen);
      if (sim.pressure) { pMaxBar = Math.max(pMaxBar, sim.pressure[i] / 1e5); pSum += sim.pressure[i] / 1e5; }
      if (Math.abs(sim.x[i]) > inp.stroke && Math.abs(sim.x[i - 1]) <= inp.stroke) hits++;
    }
    const pAbs = eAbs / N, pGen = eGen / N;
    need(Number.isFinite(pAbs) && Number.isFinite(maxX), 'The simulation diverged: reduce the time step or the end-stop stiffness.');
    const cwr = pAbs / (Jwave * 2 * a), p2a = pGen > 0 ? pPeak / pGen : 0;
    // frequency-domain reference: optimal passive damper at the peak period (regular wave of the same energy)
    const Mtot = mass + cp.A, Bopt = Math.hypot(cp.B, wp * Mtot - Kh / wp);
    const ampEq = inp.waveType === 'regular' ? inp.Hs / 2 : inp.Hs / (2 * Math.SQRT2);
    const Xopt = cp.Fe * ampEq / Math.hypot(Kh - wp * wp * Mtot, wp * (cp.B + Bopt));
    const pRef = 0.5 * Bopt * wp * wp * Xopt * Xopt;
    const checks = [
      check('Maximum stroke against the end stop', maxX, inp.stroke, 'm', 'IEC TS 62600-2:2019 §7 (PTO end-stop design)', 'Largest heave excursion in the run against the available stroke. Above 1 the float reaches the end stops: impact loads and lost energy.'),
      check('Maximum PTO force against the rating', maxF / 1e3, inp.forceLim, 'kN', 'IEC TS 62600-2:2019; PTO rating', 'Peak force in the PTO against its rating. Exceeding it overloads the cylinder, rod and bearings.'),
      ...(inp.ptoType === 'hydraulic' ? [check('Maximum accumulator pressure against the working pressure', pMaxBar, inp.pMax, 'bar', 'ISO 4413:2010 (hydraulic systems); EN 14359 (gas-loaded accumulators)', 'Highest gas pressure in the run against the design pressure. Above it the relief valve opens and energy is dumped as heat.')] : []),
    ];
    const gov = worst(checks);
    const [tt, xx, ee, ff, pp, pr] = thin([sim.t, sim.x, eta.subarray(0, sim.t.length), sim.ptoForce, pInst, sim.pressure || new Float64Array(sim.t.length)], 900);
    const [xs, fs] = thin([sim.x.subarray(i0), sim.ptoForce.subarray(i0)], 700);
    const hoursYr = 8760 * 0.3;   // order-of-magnitude equivalent hours of this sea state for the quoted annual figure
    return {
      kpis: [
        { key: 'pAbs', label: 'Mean absorbed (mechanical) power', value: round(pAbs / 1e3, 5), unit: 'kW', explain: 'Time-mean of PTO force times float velocity after the start-up transient: the hydrodynamic power taken from the waves.' },
        { key: 'pGen', label: 'Mean converted (electrical) power', value: round(pGen / 1e3, 5), unit: 'kW', explain: 'Mean power after the motor and generator efficiency. For the hydraulic PTO it is the motor flow times pressure difference, smoothed by the accumulator.' },
        { key: 'etaPto', label: 'Conversion efficiency (electrical / absorbed)', value: round(pAbs > 0 ? pGen / pAbs * 100 : 0, 4), unit: '%', explain: 'Electrical over absorbed power in the run. For the hydraulic circuit it differs from the motor efficiency while the accumulator is still charging or discharging.' },
        { key: 'cwr', label: 'Capture-width ratio', value: round(cwr, 5), unit: '-', explain: 'Absorbed power over the wave power crossing the float diameter. Compare with the frequency-domain module: nonlinear PTO and end stops normally reduce it.' },
        { key: 'pRef', label: 'Reference: optimal linear damper at the peak period', value: round(pRef / 1e3, 5), unit: 'kW', explain: 'Absorbed power of an ideal linear damper tuned to the peak period in a regular wave of the same energy (frequency domain). A yardstick for the PTO setting.' },
        { key: 'maxX', label: 'Maximum stroke', value: round(maxX, 4), unit: 'm', limit: inp.stroke, sense: 'max', explain: 'Largest displacement of the float from equilibrium after the transient.' },
        { key: 'maxV', label: 'Maximum velocity', value: round(maxV, 4), unit: 'm/s', explain: 'Peak float velocity: sets the peak flow of a hydraulic PTO or the peak speed of a direct-drive generator.' },
        { key: 'maxF', label: 'Maximum PTO force', value: round(maxF / 1e3, 5), unit: 'kN', limit: inp.forceLim, sense: 'max', explain: 'Peak force applied by the PTO on the float (and reacted by the structure and mooring).' },
        { key: 'p2a', label: 'Peak-to-average electrical power', value: round(p2a, 4), unit: '-', explain: 'Ratio of the highest instantaneous to the mean electrical power. Direct-drive PTOs reach 10 or more; hydraulic accumulators bring it to 1.5–3, which sizes the generator and converter.' },
        { key: 'sigX', label: 'Standard deviation of stroke', value: round(Math.sqrt(sx / N), 4), unit: 'm', explain: 'Root-mean-square heave after the transient; four times this value approximates the significant double amplitude.' },
        { key: 'pMean', label: 'Mean accumulator pressure', value: round(sim.pressure ? pSum / N : 0, 4), unit: 'bar', explain: 'Average gas pressure of the high-pressure accumulator (0 for non-hydraulic PTOs). It sets the PTO force level, like a Coulomb damper.' },
        { key: 'hits', label: 'End-stop contacts', value: hits, unit: '-', explain: 'Number of times the float passed the end-stop stroke in the analysed window. Any contact in an operating sea state calls for a softer setting or more stroke.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Wave elevation and float displacement', x: { label: 'Time', unit: 's' }, y: { label: 'Elevation', unit: 'm' },
          series: [{ name: 'Wave elevation', x: tt.map((v) => round(v, 3)), y: ee.map((v) => round(v, 4)) }, { name: 'Float heave', x: tt.map((v) => round(v, 3)), y: xx.map((v) => round(v, 4)) }, { name: 'End stop', x: [0, inp.duration], y: [inp.stroke, inp.stroke], dash: true }],
          explain: 'Incident wave and float motion. A float that follows the wave absorbs little; a phase lead of the velocity over the wave, or flat tops from latching, signal good energy capture.' },
        { type: 'line', title: 'PTO force', x: { label: 'Time', unit: 's' }, y: { label: 'Force', unit: 'kN' }, series: [{ name: 'PTO force', x: tt.map((v) => round(v, 3)), y: ff.map((v) => round(v / 1e3, 4)) }],
          explain: 'Force of the PTO on the float. A hydraulic or Coulomb PTO gives a square wave that follows the sign of the velocity; a linear damper gives a sinusoid.' },
        { type: 'line', title: 'Converted power', x: { label: 'Time', unit: 's' }, y: { label: 'Power', unit: 'kW' }, series: [{ name: 'Electrical power', x: tt.map((v) => round(v, 3)), y: pp.map((v) => round(v / 1e3, 4)) }, { name: 'Mean', x: [0, inp.duration], y: [pGen / 1e3, pGen / 1e3], dash: true }],
          explain: 'Instantaneous electrical output with its mean. The smoother the trace, the smaller the generator and converter rating needed for the same energy.' },
        ...(sim.pressure ? [{ type: 'line', title: 'Accumulator pressure', x: { label: 'Time', unit: 's' }, y: { label: 'Pressure', unit: 'bar' }, series: [{ name: 'High-pressure accumulator', x: tt.map((v) => round(v, 3)), y: pr.map((v) => round(v / 1e5, 4)) }, { name: 'Working-pressure limit', x: [0, inp.duration], y: [inp.pMax, inp.pMax], dash: true }],
          explain: 'Gas pressure in the accumulator. It rises during wave groups and bleeds through the motor between them; touching the limit line means the relief valve would open.' }] : []),
        { type: 'scatter', title: 'PTO force against stroke (work loop)', x: { label: 'Displacement', unit: 'm' }, y: { label: 'PTO force', unit: 'kN' }, series: [{ name: 'Operating points', x: xs.map((v) => round(v, 4)), y: fs.map((v) => round(v / 1e3, 4)) }],
          explain: 'Force–displacement loop of the PTO: the enclosed area per cycle is the absorbed energy. A rectangle is a Coulomb or hydraulic PTO, an ellipse a linear damper.' },
      ],
      tables: [
        { title: 'Run summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Incident wave power', v: round(Jwave / 1e3, 4), u: 'kW/m' }, { q: 'Mean absorbed power', v: round(pAbs / 1e3, 4), u: 'kW' }, { q: 'Mean electrical power', v: round(pGen / 1e3, 4), u: 'kW' }, { q: 'Peak electrical power', v: round(pPeak / 1e3, 4), u: 'kW' },
            { q: 'Infinite-frequency added mass', v: round(Ainf / 1e3, 4), u: 't' }, { q: 'Radiation damping at the peak period', v: round(cp.B / 1e3, 4), u: 'kN·s/m' }, { q: 'Optimal passive damping at the peak period', v: round(Bopt / 1e3, 4), u: 'kN·s/m' }, { q: 'Analysed window', v: round(N * inp.dt, 4), u: 's' }],
          explain: 'Key inputs and results of the simulation. The optimal passive damping is a starting value for a linear PTO; hydraulic settings are tuned on the mean pressure × piston area instead.' },
      ],
      notes: [
        `${inp.waveType === 'regular' ? `Regular wave H = ${inp.Hs} m, T = ${inp.Tp} s` : `Irregular sea Hs = ${inp.Hs} m, Tp = ${inp.Tp} s (seed ${inp.seed})`}; ${inp.ptoType} PTO${inp.latch > 0 ? `, latching ${inp.latch} s` : ''}${inp.declutch > 0 ? `, declutching below ${inp.declutch} m/s` : ''}; ${fmt(N * inp.dt, 4)} s analysed after the transient.`,
        `Absorbed ${fmt(pAbs / 1e3)} kW (reference optimal linear damper ${fmt(pRef / 1e3)} kW), converted ${fmt(pGen / 1e3)} kW; peak-to-average ${fmt(p2a)}.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}); ${hits} end-stop contact(s).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `The PTO converts ${fmt(pGen / 1e3)} kW on average with a peak force of ${fmt(maxF / 1e3)} kN and a maximum stroke of ${fmt(maxX)} m (${hits} end-stop contact(s)).`,
          gov.uc > 1 ? `'${gov.name}' is exceeded (utilisation ${fmt(gov.uc)}): ${gov === checks[0] ? `raise the PTO force level (higher pre-charge or damping) or provide ${fmt(maxX * 1.1)} m of stroke` : gov === checks[1] ? 'lower the pressure level or piston area, or uprate the cylinder' : `enlarge the accumulator to about ${fmt(inp.V0 * gov.uc * 1.2)} m³ or raise the motor flow coefficient`}.` : (pAbs < 0.6 * pRef ? `Absorbed power is ${fmt(pAbs / pRef * 100)} % of the optimal-damper reference: retune the PTO force level (mean pressure × piston area ≈ ${fmt(Bopt * maxV * 0.6 / 1e3)} kN) or try latching of about ${fmt(Math.max(0, (inp.Tp - 2 * Math.PI * Math.sqrt(Mtot / Kh)) / 2), 2)} s.` : 'The PTO setting is close to the linear optimum and within its limits: confirm with a longer run (1800 s) and several seeds, and calibrate drag and efficiency against bench and basin tests.'),
          'IEC TS 62600-2:2019 (design) with time-domain practice after Cummins (1962) and Falnes (2002); PTO bench testing per IEC TS 62600-103.',
          `A smoother output (peak-to-average ${fmt(p2a)}) allows a smaller generator and converter: less copper and power electronics per kWh. ${energyCo2Text(pGen * hoursYr / 1e6)}`),
        adv('economic', 'info', `Electrical output is ${fmt(pGen / 1e3)} kW in this sea state with a peak-to-average ratio of ${fmt(p2a)}; the electrical plant must be rated for about ${fmt(pPeak / 1e3)} kW unless smoothed.`,
          `Trade accumulator volume against converter rating: doubling the gas volume roughly halves the pressure ripple. If this sea state occurred 30 % of the year the device would yield about ${fmt(pGen * hoursYr / 1e6)} MWh/yr — use the power-matrix module for the annual figure.`,
          'Rating and PTO cost optimisation per IEA-OES LCOE practice: generator and converter cost scale with peak power, revenue with mean power; whole-life costing per ISO 15686-5.'),
      ],
      publish: { 'wec.ptoForceMax': maxF, 'wec.meanElectricalPower': pGen, 'wec.peakToAverage': p2a },
    };
  },
  verify: [
    { name: 'Linear damper in a regular wave: time-domain mean power equals the frequency-domain value ½·B_pto·ω²·|X|²', source: 'Falnes (2002) Ocean Waves and Oscillating Systems, Eq. 3.29–3.41 with the truncated-cylinder coefficients (hand calculation)',
      // a = 5 m, d = 6 m, T = 8 s deep water: k = 0.062902, ka = 0.31451, 2J1(ka)/ka = 0.98769, e^(−kd) = 0.68563
      // Fe = ρgπa²·0.68563·0.98769 = 534.6 kN/m; B = ωk·Fe²/(2ρg²) = 71.6 kN·s/m; M = 483.0 + 170.8 = 653.9 t; K = 789.5 kN/m
      // |Z| = √((K − ω²M)² + (ω(B + 500))²) = √(386.1² + 449.0²) = 592.2 kN/m → |X| = 0.9028 m for ζ = 1 m; P = ½·500e3·0.61685·0.8151 = 125.7 kW
      inputs: { waveType: 'regular', Hs: 2, Tp: 8, depth: 500, massRatio: 1, cdA: 0, ptoType: 'linear', Bpto: 500, Kpto: 0, stroke: 20, latch: 0, declutch: 0, duration: 400, dt: 0.05, x0: 0 },
      expect: [{ kpi: 'pAbs', value: 125.7, tol: 0.02 }, { kpi: 'maxX', value: 0.9028, tol: 0.02 }] },
  ],
};
