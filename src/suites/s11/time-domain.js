// Time-domain simulation of the floater by the Cummins equation: irregular wave excitation, constant or
// convolution (memory) radiation damping, linearised viscous damping, nonlinear station-keeping restoring and a
// numerical free-decay test, with a consistency check against the frequency-domain standard deviations.
import { retardationKernel, waveComponents, waveExcitationForce, simulateCummins, freeDecay, decayAnalysis } from '../../engine/hydro/timedomain.js';
import { hullInputs, hullFrom, motionInputs, motionFrom, seaInputs, seaFrom, spectrum, stats, geomspace, DOF, DOF_UNIT, dofScale,
  G, DEG, LIMIT_HYDRO, LIMIT_CFD, need, run, round, rnd, thin, fmt, adv, check, sev, worst, fuelCo2Text } from './_common.js';

const sub = (A, ids) => ids.map((i) => ids.map((j) => A[i][j]));
const add = (A, B) => A.map((r, i) => r.map((v, j) => v + B[i][j]));

export default {
  id: 's11-time-domain',
  suite: 11,
  title: 'Time-domain coupled dynamic analysis',
  summary: 'Simulates the unit in an irregular sea with the Cummins equation, including nonlinear station-keeping stiffness and a steady force, and reports motion histories, extremes, a free-decay test and agreement with the frequency-domain statistics.',
  covers: [6],
  fidelity: 'numerical',
  standards: ['DNV-RP-C205 (2019) §7.1.5 (time-domain analysis)', 'DNV-OS-E301 (2021) Ch.2 Sec.2', 'ISO 19901-7 (2013) §8', 'ITTC 7.5-02-07-03.1 (decay tests)'],
  models: ['Cummins impulse-response equation of motion', 'Convolution radiation memory from a retardation kernel (Haskind damping, optional)', 'Newton–Euler six-degree-of-freedom rigid body (linearised kinematics)', 'First-order wave excitation by component superposition with seeded phases', 'Stochastically linearised Morison drag damping', 'Nonlinear (hardening) mooring restoring', 'Fixed-step predictor–corrector time integration', 'Free-decay analysis (logarithmic decrement)'],
  limits: `${LIMIT_HYDRO} Small-angle kinematics and constant hydrostatic stiffness; added mass is the constant strip value. The optional memory kernel uses radiation damping from the Haskind relation applied to the slender-member excitation, tapered to zero above the slender-body frequency limit — indicative only. Line dynamics, second-order wave forces and riser reactions are not in this module (see the coupled station-keeping module). The box model simulates heave, roll and pitch as uncoupled oscillators that reproduce its RAOs. ${LIMIT_CFD}`,
  searchTerms: 'Cummins equation time domain floating platform simulation retardation',
  inputs: [
    ...hullInputs({ type: 'semi', withMesh: false }),
    ...seaInputs({ Hs: 8, Tp: 13, ret: '1', heading: 0, withDuration: false }),
    ...motionInputs({ nFreq: false }),
    { key: 'memory', label: 'Radiation damping model', type: 'select', default: 'constant', group: 'Mass distribution and damping', role: 'bc',
      options: [{ value: 'constant', label: 'Constant (user fraction of critical)' }, { value: 'haskind', label: 'Memory kernel from Haskind damping (slender models)' }],
      help: 'The memory option replaces the constant radiation damping by a convolution over the velocity history; it is available for the slender-member models only.' },
    { key: 'Fmean', label: 'Steady force along the wave direction', unit: 'kN', type: 'number', default: 1500, min: 0, max: 50000, group: 'Station keeping', role: 'load', bind: 'vessel.Fmean', help: 'Mean wind, current and wave-drift force applied as a constant load.' },
    { key: 'xHard', label: 'Offset at which the mooring stiffness doubles (0 = linear)', unit: 'm', type: 'number', default: 30, min: 0, max: 500, group: 'Station keeping', role: 'bc', help: 'Cubic hardening of the horizontal restoring force, F = K·x·(1 + (x/x_h)²): catenary and taut moorings stiffen with offset.' },
    { key: 'x0', label: 'Initial surge offset', unit: 'm', type: 'number', default: 0, min: -100, max: 100, group: 'Initial conditions', role: 'ic', help: 'Offset along x at t = 0 (velocities start from zero).' },
    { key: 'z0', label: 'Initial heave offset', unit: 'm', type: 'number', default: 0, min: -10, max: 10, group: 'Initial conditions', role: 'ic', help: 'Heave displacement at t = 0; with Hs = 0 the run is a heave decay test.' },
    { key: 'pitch0', label: 'Initial pitch', unit: 'deg', type: 'number', default: 0, min: -15, max: 15, group: 'Initial conditions', role: 'ic', help: 'Pitch angle at t = 0.' },
    { key: 'duration', label: 'Simulated duration', unit: 's', type: 'number', default: 600, min: 100, max: 10800, group: 'Numerical controls', role: 'control', help: 'Length of the record. Extremes need 3 hours and several seeds; the default is a short screening run.' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 0.25, min: 0.02, max: 1, group: 'Numerical controls', role: 'control', help: 'Fixed step of the integrator; keep it below 1/40 of the shortest natural or wave period of interest.' },
    { key: 'nComp', label: 'Wave components', unit: '-', type: 'int', default: 80, min: 30, max: 600, group: 'Numerical controls', role: 'control', help: 'Regular components representing the spectrum (unequal spacing, so the record does not repeat).' },
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 1, min: 1, max: 99999, group: 'Numerical controls', role: 'control', help: 'Seed of the component phases; change it to obtain another realisation of the same sea.' },
    { key: 'limOffset', label: 'Allowable maximum horizontal offset', unit: 'm', type: 'number', default: 25, min: 1, max: 500, group: 'Acceptance criteria', role: 'criteria', help: 'Riser or umbilical watch-circle limit on the horizontal excursion.' },
    { key: 'limHeaveMax', label: 'Allowable maximum heave', unit: 'm', type: 'number', default: 6, min: 0.2, max: 30, group: 'Acceptance criteria', role: 'criteria', help: 'Stroke or survival limit on heave amplitude.' },
    { key: 'limPitchMax', label: 'Allowable maximum roll / pitch', unit: 'deg', type: 'number', default: 8, min: 0.5, max: 30, group: 'Acceptance criteria', role: 'criteria', help: 'Limit on the largest inclination in the record.' },
  ],
  meshTime: { key: 'dt', label: 'Time step', levels: [0.5, 0.25, 0.125], qoi: 'heaveStd', isSize: true },
  compute(inp) {
    const hull = hullFrom(inp); const mot = motionFrom(inp, hull);
    const sea = seaFrom(inp); const beta = inp.headingDeg;
    const slender = mot.kind === 'slender';
    const Tmin = Math.min(...mot.Tn.filter((v) => v != null));
    need(inp.dt <= Tmin / 12, `Time step ${fmt(inp.dt)} s is too coarse for the shortest natural period (${fmt(Tmin)} s): use ${fmt(Tmin / 40)} s or less.`);
    const om = geomspace(0.12, 2.2, inp.nComp); const S = spectrum(om, sea);
    const r = mot.rao(om, beta, sea);
    const comps = waveComponents({ omega: om, S, seed: inp.seed });
    let ids; let M; let C; let B; let Fre; let Fim; let kernel = null;
    if (slender) {
      need(inp.Kmoor > 0 || hull.type === 'tlp', 'Under-constrained: with zero station-keeping stiffness the unit has no horizontal restoring and drifts away in a time-domain run. Enter the mooring or DP stiffness.');
      ids = [0, 1, 2, 3, 4, 5];
      const m = mot.model; M = add(m.M, m.A); C = add(m.C, m.Kext);
      B = r.viscousDamping.map((row) => row.slice());
      const useMem = inp.memory === 'haskind';
      if (!useMem) for (let d = 0; d < 6; d++) if (C[d][d] > 0) B[d][d] += 2 * inp.zetaRad * Math.sqrt(C[d][d] * M[d][d]);
      if (useMem) {
        // Haskind: B_ii(ω) = ω³/(4πρg³)·∫|X_i(β)|² dβ (deep water), heading integral by 12 headings, drag-free excitation
        const wk = geomspace(0.05, mot.validOmegaMax * 1.3, 60); const nb = 12;
        const acc = Array.from({ length: 6 }, () => new Float64Array(wk.length));
        for (let ib = 0; ib < nb; ib++) { const rr = mot.raoLin(wk, ib * 360 / nb, null); for (let d = 0; d < 6; d++) for (let i = 0; i < wk.length; i++) acc[d][i] += (rr.force.re[d][i] ** 2 + rr.force.im[d][i] ** 2) * (2 * Math.PI / nb); }
        const taper = (w) => { const a = 0.8 * mot.validOmegaMax; const b = 1.3 * mot.validOmegaMax; return w <= a ? 1 : w >= b ? 0 : 0.5 * (1 + Math.cos(Math.PI * (w - a) / (b - a))); };
        const Bw = wk.map((w, i) => Array.from({ length: 6 }, (_, p) => Array.from({ length: 6 }, (_, q) => (p === q ? taper(w) * w ** 3 / (4 * Math.PI * hull.rho * G ** 3) * acc[p][i] : 0))));
        kernel = run(() => retardationKernel({ omega: wk, B: Bw, dt: inp.dt, tMax: 60 }), 'Retardation kernel');
      }
      Fre = r.force.re; Fim = r.force.im;
    } else {
      ids = [2, 3, 4];
      const K = hull.hs.stiffness; const zeta = [inp.zetaRad, inp.zetaRoll ?? 0.06, inp.zetaRad];
      const kk = ids.map((d) => K[d][d]); const mm = ids.map((d, j) => kk[j] * (mot.Tn[d] / (2 * Math.PI)) ** 2); const bb = mm.map((m, j) => 2 * zeta[j] * Math.sqrt(kk[j] * m));
      const diag = (v) => v.map((x, i) => v.map((_, j) => (i === j ? x : 0)));
      M = diag(mm); C = diag(kk); B = diag(bb);
      // excitation that reproduces the box-model RAO: F = (C − ω²M + iωB)·ξ
      Fre = ids.map((d, j) => Float64Array.from(om, (w, i) => (kk[j] - w * w * mm[j]) * r.re[d][i] - w * bb[j] * r.im[d][i]));
      Fim = ids.map((d, j) => Float64Array.from(om, (w, i) => (kk[j] - w * w * mm[j]) * r.im[d][i] + w * bb[j] * r.re[d][i]));
    }
    const n = ids.length;
    const wave = inp.Hs > 0 ? waveExcitationForce({ components: comps, Fre, Fim }) : null;
    const Kh = slender ? inp.Kmoor * 1e3 : 0; const cb = Math.cos(beta * DEG); const sb = Math.sin(beta * DEG);
    const Fm = inp.Fmean * 1e3; const ramp = 2 * inp.Tp;
    // the integrator evaluates the force several times at each grid time: tabulate the wave force once per step
    const cache = new Map(); const zero = new Array(n).fill(0);
    const waveAt = (t) => { if (!wave) return zero; const k = Math.round(t / inp.dt); if (Math.abs(t - k * inp.dt) > 1e-9) return wave(t); let f = cache.get(k); if (!f) { f = wave(t); cache.set(k, f); } return f; };
    const force = (t, x) => {
      const f = waveAt(t);
      const g = Math.min(1, t / ramp);
      const out = new Array(n);
      for (let i = 0; i < n; i++) out[i] = g * f[i];
      if (slender) {
        out[0] += g * Fm * cb; out[1] += g * Fm * sb;
        if (inp.xHard > 0 && Kh > 0) { const rr = (x[0] * x[0] + x[1] * x[1]) / (inp.xHard * inp.xHard); out[0] -= Kh * x[0] * rr; out[1] -= Kh * x[1] * rr; }
      }
      return out;
    };
    const x0 = new Array(n).fill(0);
    if (slender) { x0[0] = inp.x0; x0[2] = inp.z0; x0[4] = inp.pitch0 * DEG; } else { x0[0] = inp.z0; x0[2] = inp.pitch0 * DEG; }
    const nSteps = Math.round(inp.duration / inp.dt);
    need(nSteps <= 120000, 'The run is too long for an interactive simulation: reduce the duration or increase the time step (limit 120 000 steps).');
    const sim = run(() => simulateCummins({ M, Blin: B, C, kernel: kernel || undefined, force, x0, dt: inp.dt, nSteps }), 'Time integration');
    const i0 = Math.floor(0.2 * sim.t.length);
    const st = sim.x.map((x) => { let s = 0; let s2 = 0; let mx = 0; let c = 0; for (let i = i0; i < x.length; i++) { s += x[i]; c++; } const mean = s / c; for (let i = i0; i < x.length; i++) { s2 += (x[i] - mean) ** 2; mx = Math.max(mx, Math.abs(x[i])); } return { mean, std: Math.sqrt(s2 / c), max: mx }; });
    for (const s of st) need(Number.isFinite(s.std) && s.max < 1e6, 'The simulation diverged: reduce the time step (the integrator needs about 40 steps per shortest period) or check the stiffness and damping.');
    const idx = (d) => ids.indexOf(d);
    const g = (d) => (idx(d) >= 0 ? st[idx(d)] : { mean: 0, std: 0, max: 0 });
    const fd = (d) => stats(om, r.amp[d], S).sigma;
    const ratio = fd(2) > 0 ? g(2).std / fd(2) : 1;
    const offMax = (() => { if (!slender) return 0; let m = 0; const X = sim.x[0]; const Y = sim.x[1]; for (let i = i0; i < X.length; i++) m = Math.max(m, Math.hypot(X[i], Y[i])); return m; })();
    const offStatic = Kh > 0 ? (() => { // F = K·x·(1 + (x/xh)²), solved by Newton
      let x = Fm / Kh; if (inp.xHard > 0) for (let i = 0; i < 50; i++) { const f = Kh * x * (1 + (x / inp.xHard) ** 2) - Fm; x -= f / (Kh * (1 + 3 * (x / inp.xHard) ** 2)); } return x; })() : 0;
    // numerical heave decay test on the same inertia, stiffness and linear damping
    const jh = idx(2); const mh = M[jh][jh]; const kh = C[jh][jh]; const Tn3 = 2 * Math.PI * Math.sqrt(mh / kh);
    const dec = freeDecay({ mass: mh, stiffness: kh, b1: slender && inp.memory === 'haskind' ? 2 * inp.zetaRad * Math.sqrt(kh * mh) : B[jh][jh], x0: 1, dt: Tn3 / 400, duration: 8 * Tn3 });
    const da = decayAnalysis({ t: dec.t, x: dec.x, inertia: mh });
    const rotMax = Math.max(g(3).max, g(4).max) / DEG;
    const checks = [
      check('Maximum horizontal offset', offMax, inp.limOffset, 'm', 'ISO 19901-7 (2013) §8; riser watch-circle limit (project)', 'Guards risers, umbilicals and gangways, which tolerate only a limited horizontal excursion before they overstress or must be disconnected.'),
      check('Maximum heave', g(2).max, inp.limHeaveMax, 'm', 'Riser tensioner / compensator stroke (project criterion)', 'Guards tensioner and compensator stroke-out, which transfers the full motion to the riser or load.'),
      check('Maximum roll / pitch', rotMax, inp.limPitchMax, 'deg', 'DNV-OS-C301 / topsides equipment limits (project criterion)', 'Guards process equipment and sea-fastening, which are designed for a maximum inclination.'),
    ];
    const gov = worst(checks);
    const show = slender ? [0, 2, 4] : [2, 3, 4];
    const hist = thin([sim.t, ...show.map((d) => sim.x[idx(d)])], 900);
    return {
      kpis: [
        { key: 'heaveStd', label: 'Heave standard deviation', value: round(g(2).std, 4), unit: 'm', explain: 'Standard deviation of heave after the start-up transient; four times this is the significant double amplitude.' },
        { key: 'heaveMax', label: 'Maximum heave', value: round(g(2).max, 4), unit: 'm', limit: inp.limHeaveMax, sense: 'max', explain: 'Largest absolute heave in the record. One realisation only: design extremes need several seeds or a fitted extreme-value distribution.' },
        { key: 'surgeMean', label: 'Mean surge offset', value: round(g(0).mean, 4), unit: 'm', explain: 'Average surge after start-up: the static offset under the steady force on the (hardening) mooring.' },
        { key: 'surgeStd', label: 'Surge standard deviation', value: round(g(0).std, 4), unit: 'm', explain: 'Standard deviation of wave-frequency surge about the mean offset.' },
        { key: 'offsetMax', label: 'Maximum horizontal offset', value: round(offMax, 4), unit: 'm', limit: inp.limOffset, sense: 'max', explain: 'Largest horizontal distance from the reference position: mean offset plus dynamic excursion.' },
        { key: 'pitchStd', label: 'Pitch standard deviation', value: round(g(4).std / DEG, 4), unit: 'deg', explain: 'Standard deviation of pitch after start-up.' },
        { key: 'pitchMax', label: 'Maximum pitch', value: round(g(4).max / DEG, 4), unit: 'deg', explain: 'Largest absolute pitch in the record.' },
        { key: 'rollMax', label: 'Maximum roll', value: Math.abs(g(3).max / DEG) < 1e-6 ? 0 : round(g(3).max / DEG, 4), unit: 'deg', explain: 'Largest absolute roll in the record (zero in head or following seas). A roll below 1e-6 deg is round-off of a symmetric case, far below anything measurable, and is reported as 0.' },
        { key: 'fdRatio', label: 'Heave σ: time domain / frequency domain', value: round(ratio, 4), unit: '-', explain: 'Standard deviation from the simulation divided by the spectral value from the same RAO. Values near 1.0 confirm the time integration; a short record scatters by ±10 %.' },
        { key: 'offsetStatic', label: 'Static offset from the steady force', value: round(offStatic, 4), unit: 'm', explain: 'Offset at which the hardening mooring balances the steady force, solved directly: the simulated mean surge should match it.' },
        { key: 'decayPeriod', label: 'Heave period from the decay test', value: round(da.period, 5), unit: 's', explain: 'Damped heave period measured from a simulated free decay; compare with the natural period 2π√((M + A)/C).' },
        { key: 'decayZeta', label: 'Heave damping ratio from the decay test', value: round(da.dampingRatio, 5), unit: '-', explain: 'Equivalent linear damping ratio from the logarithmic decrement of the simulated decay: the number to compare with a model-test decay.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Motion histories', x: { label: 'Time', unit: 's' }, y: { label: 'Motion', unit: 'm or deg' },
          series: show.map((d, j) => ({ name: `${DOF[d]} (${DOF_UNIT[d]})`, x: rnd(hist[0]), y: rnd(hist[j + 1].map((v) => v * dofScale(d))) })),
          explain: 'Simulated motions against time. The first two peak periods are a force ramp; look for the slow surge build-up to the mean offset and for wave groups driving the largest heave and pitch.' },
        { type: 'line', title: 'Heave free decay (numerical test)', x: { label: 'Time', unit: 's' }, y: { label: 'Heave', unit: 'm' },
          series: [(() => { const [t, x] = thin([dec.t, dec.x], 400); return { name: 'Heave', x: rnd(t), y: rnd(x) }; })()],
          explain: 'Heave released from 1 m in calm water with the model damping. The spacing of the peaks is the natural period and their decay rate the damping: the pattern a tank decay test should reproduce.' },
        ...(kernel ? [{ type: 'line', title: 'Heave retardation function', x: { label: 'Lag', unit: 's' }, y: { label: 'K33', unit: 'MN/m' },
          series: [{ name: 'K33(t)', x: rnd(Array.from(kernel.t).filter((_, i) => i % Math.ceil(kernel.t.length / 300) === 0)), y: rnd(kernel.K.filter((_, i) => i % Math.ceil(kernel.t.length / 300) === 0).map((k) => k[2 * 6 + 2] / 1e6)) }],
          explain: 'Impulse-response (memory) function of heave radiation: how a past heave velocity still pushes on the hull through the waves it made. It must have decayed to near zero by the end of the lag window.' }] : []),
      ],
      tables: [
        { title: 'Time-domain statistics', columns: [{ key: 'dof', label: 'DOF' }, { key: 'mean', label: 'Mean' }, { key: 'std', label: 'Std dev' }, { key: 'max', label: 'Max |x|' }, { key: 'fd', label: 'Std dev (frequency domain)' }, { key: 'u', label: 'Unit' }],
          rows: ids.map((d, j) => ({ dof: DOF[d], mean: round(st[j].mean * dofScale(d), 4), std: round(st[j].std * dofScale(d), 4), max: round(st[j].max * dofScale(d), 4), fd: round(fd(d) * dofScale(d), 4), u: DOF_UNIT[d] })),
          explain: 'Statistics of each simulated degree of freedom after discarding the first 20 % of the record, beside the frequency-domain standard deviation for the same sea.' },
      ],
      notes: [
        `${slender ? 'Six-degree-of-freedom' : 'Heave–roll–pitch (box model)'} simulation of ${fmt(inp.duration)} s at dt = ${fmt(inp.dt)} s (${nSteps} steps, ${inp.nComp} wave components, seed ${inp.seed}); radiation damping: ${kernel ? 'memory kernel (Haskind)' : 'constant'}.`,
        `Heave σ = ${fmt(g(2).std)} m against ${fmt(fd(2))} m from the spectrum (ratio ${fmt(ratio)}); mean surge ${fmt(g(0).mean)} m against the static solution ${fmt(offStatic)} m.`,
        `Decay test: period ${fmt(da.period)} s, damping ${fmt(100 * da.dampingRatio)} % of critical. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit} in this ${fmt(inp.duration / 60)}-minute realisation.`,
          gov.uc > 1 ? 'Reduce the excursion with more station-keeping stiffness or pretension (offset), a longer natural period (heave, pitch) or added damping; then confirm with 3-hour runs.'
            : `A single ${fmt(inp.duration / 60)}-minute record under-samples the extreme: run at least 10 seeds of 3 hours and fit the maxima (Gumbel) before using them as design values; refine dt until heave σ changes by less than 2 %.`,
          'DNV-RP-C205 (2019) §7.1.5; DNV-OS-E301 (2021) Ch.2 Sec.2 (time-domain analysis, number of realisations); ISO 19901-7 (2013).'),
        adv('economic', 'info', `Time-domain analysis is the expensive tier: this screening run uses ${nSteps} steps; a design matrix of 10 seeds × 3 h × load cases is roughly ${fmt(10 * 10800 / inp.duration)} times that per case.`,
          'Use the frequency-domain modules to rank load cases and keep time-domain runs for the few governing ones where nonlinear mooring or large motions matter.',
          'Analysis-tier selection per DNV-OS-E301 (2021) and API RP 2SK (3rd ed.); value-of-information reasoning for model tests.',
          `Demonstrating a smaller offset can remove a mooring upgrade; each avoided day of anchor-handling vessel time saves about 30 t of fuel: ${fuelCo2Text(30)}`),
      ],
    };
  },
  verify: [
    { // Spar D = 30 m, T = 90 m: Tn = 2π√((90 + 4500/706.858)/9.80665) = 19.696 s (hand, see the 6-DOF module); damping set to 4 %
      // of critical → log decrement 2πζ/√(1−ζ²) and measured ratio 0.04; damped period 19.696/√(1−0.0016) = 19.712 s.
      name: 'Heave free decay of a spar: period and damping ratio recovered',
      source: 'Single-degree-of-freedom free vibration (logarithmic decrement); Faltinsen (1990) ch. 3 for the natural period',
      inputs: { hullType: 'spar', Dspar: 30, Hspar: 105, draft: 90, KG: 38, Hs: 0, zetaRad: 0.04, Fmean: 0, duration: 200 },
      expect: [{ kpi: 'decayPeriod', value: 19.712, tol: 0.005 }, { kpi: 'decayZeta', value: 0.04, tol: 0.03 }, { kpi: 'heaveStd', value: 0, abs: 1e-9 }] },
    { // Calm water, linear mooring 250 kN/m, steady 1500 kN: mean offset = F/K = 6.0 m.
      name: 'Static offset under a steady force on a linear mooring: x = F/K',
      source: 'Hooke\'s law; force balance of the time integration',
      inputs: { Hs: 0, Fmean: 1500, Kmoor: 250, xHard: 0, duration: 1500, zetaRad: 0.3 },
      expect: [{ kpi: 'offsetStatic', value: 6.0, tol: 1e-6 }, { kpi: 'surgeMean', value: 6.0, tol: 0.03 }] },
  ],
};
