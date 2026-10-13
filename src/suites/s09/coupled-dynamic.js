// Coupled dynamic riser response in an irregular sea: vessel motion, riser inertia, drag and seabed contact.
// Default method is hybrid (frequency-domain transfer functions on the beam model, then time-history synthesis);
// an opt-in nonlinear 3-D lumped-mass time-domain run is available for regular waves.
import { buildLumpedLine, simulateLine } from '../../engine/lines/lumpedmass.js';
import { makeRng } from '../../engine/math/rng.js';
import {
  DEG, need, run, fmt, round, thin, scaled, adv, check, kpi, plot, meshInput, finish, handWeight, currentInputs, currentFn, steelTonnes, steelCo2Text, CONFIRM,
} from './_common.js';
import { scrInputs, criteriaInputs, buildScr, currentOf, stressOf, tdpNode } from './_scr.js';
import { vesselInputs, hangoffRao, vesselRangeAdvice, LIMIT_VESSEL } from './_vessel.js';
import { dynInputs, stressAmp, sweep, seaGrid, responseStats, interpAt, mpm } from './_dyn.js';

const lm = { key: 'method', in: ['lumped'] }, hy = { key: 'method', in: ['hybrid'] };

export default {
  id: 's09-coupled-dynamic', suite: 9, title: 'Coupled dynamic riser response',
  summary: 'Time histories and statistics of top tension, hang-off reactions, riser motion and touchdown stress for a catenary riser coupled to host-vessel motion in an irregular sea with current.',
  covers: [2], fidelity: 'reduced-order',
  standards: ['API RP 2RD (1998, R2013) §5.5', 'DNV-ST-F201 (2018) Sec. 4, App. A', 'DNV-RP-C205 (2019) §3.3, §6.2', 'DNV-RP-F205 (2019) coupled analysis'],
  models: ['Hybrid frequency-domain / time-domain response synthesis', 'Rod dynamics with consistent mass and added mass', 'Linearised Morison drag (hybrid) and full quadratic drag (lumped mass)', 'Vessel response amplitude operators', 'JONSWAP spectrum and random-phase wave synthesis', 'Lumped-mass cable dynamics with seabed contact and Airy wave kinematics (option)'],
  limits: 'Hybrid method: linear about the static shape (frozen seabed contact, equivalent-linear drag, in-plane motion, vessel motion one-way coupled to the riser, no direct wave load on the riser); response phases are not resolved, so combined quantities are statistical. Lumped-mass option: 3-D nonlinear cable dynamics with seabed contact, quadratic drag and Airy kinematics, regular waves only and explicit integration (slow). Fully coupled vessel–mooring–riser simulation with a diffraction database is outside this module.' + LIMIT_VESSEL,
  publishes: ['riser.dynTensionSigma', 'riser.dynTensionMpm'],
  searchTerms: 'coupled dynamic riser analysis time domain irregular sea',
  inputs: [
    { key: 'method', label: 'Dynamic method', type: 'select', default: 'hybrid', group: 'Dynamic model', role: 'control',
      options: [{ value: 'hybrid', label: 'Hybrid frequency/time domain, irregular sea (fast)' }, { value: 'lumped', label: 'Nonlinear lumped-mass time domain, regular wave (slow)' }],
      help: 'The lumped-mass run integrates explicitly at millisecond steps and can take tens of seconds.' },
    ...scrInputs(),
    ...currentInputs(),
    ...vesselInputs(),
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 10, min: 0.1, max: 22, group: 'Waves', role: 'bc', bind: 'site.Hs100', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'Tp', label: 'Spectral peak period', unit: 's', type: 'number', default: 14, min: 3, max: 30, group: 'Waves', role: 'bc', bind: 'site.Tp100' },
    { key: 'gamma', label: 'JONSWAP peak factor', unit: '-', type: 'number', default: 2.5, min: 1, max: 7, group: 'Waves', role: 'bc' },
    { key: 'duration', label: 'Simulated duration', unit: 's', type: 'number', default: 600, min: 60, max: 1800, group: 'Initial state and time stepping', role: 'ic', showIf: hy,
      help: 'The run starts from the converged static equilibrium with zero dynamic velocity; statistics for 3 h are extrapolated from the spectra.' },
    { key: 'dt', label: 'Output time step', unit: 's', type: 'number', default: 0.5, min: 0.1, max: 2, group: 'Initial state and time stepping', role: 'control', showIf: hy },
    { key: 'seed', label: 'Random-phase seed', type: 'int', default: 1, min: 1, max: 99999, group: 'Initial state and time stepping', role: 'ic', showIf: hy, help: 'Fixes the wave realisation so runs are repeatable.' },
    ...dynInputs(),
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.2, min: 1, max: 5, group: 'Dynamic model', role: 'material' },
    { key: 'nSeg', label: 'Lumped-mass segments', type: 'int', default: 24, min: 8, max: 120, group: 'Lumped-mass option', role: 'control', showIf: lm },
    { key: 'nPeriods', label: 'Wave periods simulated', type: 'int', default: 4, min: 2, max: 20, group: 'Lumped-mass option', role: 'control', showIf: lm },
    { key: 'curHeading', label: 'Current heading from the riser plane', unit: 'deg', type: 'number', default: 0, min: 0, max: 180, group: 'Lumped-mass option', role: 'bc', showIf: lm, help: 'The lumped-mass model is three-dimensional: 90° gives a cross current.' },
    ...criteriaInputs().map((c) => (c.key === 'Cf' ? { ...c, default: 1.2 } : c)),
    meshInput(160),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [80, 160, 320, 640], qoi: 'sigT' },
  meshTime: { key: 'dt', label: 'Output time step', levels: [1, 0.5, 0.25], qoi: 'TtopMaxSim', isSize: true },

  async compute(inp, ctx) {
    // engineering units → SI in buildScr; wave heights m, periods s
    const { pp, base, pressure } = buildScr(inp, ctx, { nElem: inp.nElem, current: currentOf(inp, 1) });
    const n = base.s.length - 1, kT = tdpNode(base), Ts = base.topTension;
    const sw = sweep(base, pp, inp, inp.Hs / 2);
    const grid = seaGrid(sw, inp, 160);
    const col = (f) => sw.rows.map(f);
    const hT = col((r) => r.resp.tensionAmp[0] / r.a);
    const sAll = sw.rows.map((r) => stressAmp(r.resp, pp, inp.scf, r.mv));
    // most stressed touchdown-zone node at the period nearest the spectral peak
    const iPk = sw.rows.reduce((b, r, i) => (Math.abs(r.T - inp.Tp) < Math.abs(sw.rows[b].T - inp.Tp) ? i : b), 0);
    const sT = base.s[kT]; let kS = kT, best = 0;
    for (let k = 0; k <= n; k++) if (Math.abs(base.s[k] - sT) < 0.25 * sT && sAll[iPk][k] > best) { best = sAll[iPk][k]; kS = k; }
    const hS = sAll.map((a, i) => a[kS] / sw.rows[i].a);
    const hTd = col((r) => r.resp.tensionAmp[kT] / r.a);
    const hZ = col((r) => r.heave), hX = col((r) => r.surge);
    const stT = responseStats(sw, grid, hT), stS = responseStats(sw, grid, hS), stZ = responseStats(sw, grid, hZ), stX = responseStats(sw, grid, hX), stTd = responseStats(sw, grid, hTd);
    const stV = responseStats(sw, grid, col((r) => r.heave * r.omega)), stA = responseStats(sw, grid, col((r) => r.heave * r.omega ** 2));
    // standard deviation of motion and tension along the riser
    const sigN = new Float64Array(n + 1), sigTe = new Float64Array(n + 1), sigVn = new Float64Array(n + 1), sigAn = new Float64Array(n + 1);
    const nf = sw.rows.length, q = new Float64Array(nf), qv = new Float64Array(nf), qa = new Float64Array(nf), qt = new Float64Array(nf);
    for (let k = 0; k <= n; k++) {
      for (let i = 0; i < nf; i++) { const r = sw.rows[i], v = r.resp.normalAmp[k] / r.a; q[i] = v; qv[i] = v * r.omega; qa[i] = v * r.omega * r.omega; qt[i] = r.resp.tensionAmp[k] / r.a; }
      sigN[k] = responseStats(sw, grid, q).sigma; sigVn[k] = responseStats(sw, grid, qv).sigma; sigAn[k] = responseStats(sw, grid, qa).sigma; sigTe[k] = responseStats(sw, grid, qt).sigma;
    }
    // time-history synthesis: random-phase components; response phase relative to the wave is not resolved
    const rng = makeRng(inp.seed), M = grid.omega.length;
    const amp = grid.S.map((S) => Math.sqrt(2 * S * grid.dw)), ph = grid.omega.map(() => 2 * Math.PI * rng.next());
    const gT = grid.omega.map((w) => interpAt(sw.omega, hT, w)), gS = grid.omega.map((w) => interpAt(sw.omega, hS, w)), gZ = grid.omega.map((w) => interpAt(sw.omega, hZ, w)), gX = grid.omega.map((w) => interpAt(sw.omega, hX, w));
    const nt = Math.min(3600, Math.floor(inp.duration / inp.dt) + 1);
    const t = new Array(nt), eta = new Array(nt), Tt = new Array(nt), St = new Array(nt), Zt = new Array(nt), Xt = new Array(nt);
    const ramp = Math.min(30, 0.1 * inp.duration);   // smooth start from the static state (zero initial velocity)
    for (let i = 0; i < nt; i++) {
      const ti = i * inp.dt, r = ti < ramp ? 0.5 * (1 - Math.cos(Math.PI * ti / ramp)) : 1;
      let e = 0, a = 0, b = 0, z = 0, x = 0;
      for (let j = 0; j < M; j++) { const c = amp[j] * Math.cos(grid.omega[j] * ti + ph[j]); e += c; a += c * gT[j]; b += c * gS[j]; z += c * gZ[j]; x += amp[j] * gX[j] * Math.sin(grid.omega[j] * ti + ph[j]); }
      t[i] = ti; eta[i] = r * e; Tt[i] = Ts + r * a; St[i] = r * b; Zt[i] = r * z; Xt[i] = r * x;
    }
    const TmaxSim = Math.max(...Tt), TminSim = Math.min(...Tt);
    const mpmT = mpm(stT.sigma, stT.Tz), mpmS = mpm(stS.sigma, stS.Tz), mpmTd = mpm(stTd.sigma, stTd.Tz);
    const st = stressOf(base, pp, pressure);
    const fy = pp.mat.fy, allow = inp.Cf * inp.CaAllow * fy;
    const vmX = st.vonMises[kS] + mpmS / inp.scf;
    const th0 = base.hangoffAngleDeg * DEG;

    const kpis = [
      kpi('TtopStatic', 'Static top tension', Ts / 1e3, 'kN', 'Top effective tension at the converged static equilibrium (weight, buoyancy, contents, current): the initial condition of the dynamic run.'),
      kpi('sigT', 'Dynamic top-tension standard deviation', stT.sigma / 1e3, 'kN', 'Standard deviation of the fluctuating top tension in the sea state. It measures how strongly vessel motion is transmitted into axial load.'),
      kpi('TtopMaxSim', 'Maximum top tension in the simulated record', TmaxSim / 1e3, 'kN', 'Largest top tension in the synthesised time history. It depends on the realisation (seed) and the duration; the 3-hour most probable maximum is the design value.'),
      kpi('TtopMinSim', 'Minimum top tension in the simulated record', TminSim / 1e3, 'kN', 'Smallest top tension in the record. A value near zero would mean the riser unloads at the porch.'),
      kpi('mpmT', 'Most probable maximum top tension in 3 h', (Ts + mpmT) / 1e3, 'kN', 'Static tension plus the most probable largest dynamic tension in a three-hour storm (Rayleigh peaks). Use it as the extreme hang-off tension for this sea state.'),
      kpi('daf', 'Amplification (3-h maximum / static)', (Ts + mpmT) / Ts, '-', 'Ratio of the 3-hour extreme top tension to the static value; compare with the dynamic factor assumed in the connection design.'),
      kpi('TzT', 'Mean zero-crossing period of top tension', stT.Tz, 's', 'Average period of tension cycles; with the duration it gives the number of load cycles for fatigue.'),
      kpi('minTeTdp', 'Minimum effective tension at touchdown (3 h)', (base.Te[kT] - mpmTd) / 1e3, 'kN', 'Static tension at the touchdown point minus its 3-hour most probable dynamic amplitude. Negative means compression at the sag bend during the storm.', 0, 'min'),
      kpi('sigStressTdz', 'Touchdown hot-spot stress standard deviation', stS.sigma / 1e6, 'MPa', 'Standard deviation of the cyclic hot-spot stress (with SCF) at the most loaded touchdown-zone node; the basis for fatigue in this sea state.'),
      kpi('vmExtreme', 'Extreme von Mises stress at the touchdown zone', vmX / 1e6, 'MPa', 'Static von Mises stress at the critical touchdown node plus the 3-hour most probable dynamic stress (without SCF), against the extreme-case allowable.', allow / 1e6, 'max'),
      kpi('sigHeave', 'Vertical hang-off displacement, standard deviation', stZ.sigma, 'm', 'Standard deviation of vertical porch motion (heave plus lever arm × pitch) in the sea state — the prescribed top boundary motion.'),
      kpi('sigSurge', 'Horizontal hang-off displacement, standard deviation', stX.sigma, 'm', 'Standard deviation of in-plane horizontal porch motion in the sea state.'),
      kpi('sigVel', 'Vertical hang-off velocity, standard deviation', stV.sigma, 'm/s', 'Standard deviation of the vertical velocity imposed on the riser top; high values relative to the pipe sinking speed lead to tension loss.'),
      kpi('sigAcc', 'Vertical hang-off acceleration, standard deviation', stA.sigma, 'm/s²', 'Standard deviation of vertical acceleration at the porch, comparable with measured vessel motion records.'),
      kpi('sigDispMax', 'Largest riser displacement standard deviation', Math.max(...sigN), 'm', 'Largest standard deviation of motion normal to the riser axis anywhere along its length; it shows how far the imposed motion penetrates down the riser.'),
      kpi('sigRiserVel', 'Largest riser normal velocity, standard deviation', Math.max(...sigVn), 'm/s', 'Largest standard deviation of the riser velocity normal to its axis; it controls the drag damping and the relative-velocity loading.'),
      kpi('sigRiserAcc', 'Largest riser normal acceleration, standard deviation', Math.max(...sigAn), 'm/s²', 'Largest standard deviation of the riser acceleration normal to its axis, which scales the inertia and added-mass forces.'),
    ];
    const checks = [
      check('Extreme combined stress at touchdown (von Mises)', vmX / 1e6, allow / 1e6, 'MPa', 'API RP 2RD (1998) §5.2.3, Table 2 (extreme case)', 'Guards against yielding at the touchdown zone under static stress plus the 3-hour extreme dynamic stress.'),
      check('Dynamic tension at touchdown vs static tension', mpmTd / 1e3, Math.max(base.Te[kT], 1) / 1e3, 'kN', 'DNV-ST-F201 (2018) Sec. 5 — effective compression', 'Compares the extreme dynamic tension amplitude at touchdown with the static tension. Above 1.0 the sag bend goes into compression during the storm and may buckle.'),
    ];

    const dec = (a) => thin([t, a], 1200);
    const [tt, e1] = dec(eta), T1 = dec(Tt)[1], S1 = dec(St)[1], Z1 = dec(Zt)[1], X1 = dec(Xt)[1];
    const [sx, sN, sTe] = thin([base.s, sigN, sigTe], 300);
    const plots = [
      plot('Top tension history', 'Time', 's', 'Top effective tension', 'kN', [{ name: 'Top tension', x: tt, y: scaled(T1, 1e-3) }, { name: 'Static', x: [0, tt[tt.length - 1]], y: [Ts / 1e3, Ts / 1e3], dash: true }], 'Top tension against time for the simulated realisation, starting from static equilibrium. Look for the largest excursions and whether the trace approaches zero.'),
      plot('Hang-off reaction histories', 'Time', 's', 'Reaction', 'kN', [{ name: 'Vertical reaction', x: tt, y: T1.map((v) => v * Math.cos(th0) / 1e3) }, { name: 'Horizontal reaction', x: tt, y: T1.map((v) => v * Math.sin(th0) / 1e3) }], 'Vertical and horizontal components of the hang-off force (top tension resolved at the mean riser angle). These are the dynamic interface loads on the hull.'),
      plot('Wave elevation and hang-off motion', 'Time', 's', 'Elevation / displacement', 'm', [{ name: 'Wave elevation', x: tt, y: e1 }, { name: 'Vertical hang-off motion', x: tt, y: Z1 }, { name: 'Horizontal hang-off motion', x: tt, y: X1 }], 'The wave realisation and the resulting porch motions that are imposed on the riser top.'),
      plot('Touchdown hot-spot stress history', 'Time', 's', 'Dynamic stress', 'MPa', [{ name: 'Hot-spot stress', x: tt, y: scaled(S1, 1e-6) }], 'Cyclic hot-spot stress at the most loaded touchdown-zone node. Rainflow counting of this trace feeds time-domain fatigue.'),
      plot('Dynamic response along the riser', 'Arc length from hang-off', 'm', 'Standard deviation', 'm | kN', [{ name: 'Normal displacement σ (m)', x: sx, y: sN }, { name: 'Effective tension σ (kN ÷ 100)', x: sx, y: scaled(sTe, 1e-5) }], 'Standard deviation of normal displacement and of effective tension along the riser. Displacement decays with depth through drag; tension variation is carried almost uniformly down to the touchdown.'),
      plot('Transfer functions', 'Wave period', 's', 'Response per metre of wave amplitude', 'kN/m | MPa/m', [{ name: 'Top tension (kN/m)', x: col((r) => r.T), y: hT.map((v) => v / 1e3) }, { name: 'Touchdown stress (MPa/m ×10)', x: col((r) => r.T), y: hS.map((v) => v / 1e5) }], 'Transfer functions of top tension and touchdown stress used to build the time histories. Their peaks show the wave periods the riser responds to.'),
    ];
    const tables = [{ title: 'Response statistics', columns: [{ key: 'q', label: 'Quantity' }, { key: 'sigma', label: 'Standard deviation' }, { key: 'Tz', label: 'Zero-crossing period', unit: 's' }, { key: 'mpm', label: '3-h most probable maximum' }, { key: 'unit', label: 'Unit' }],
      rows: [['Top tension (dynamic)', stT, 1e-3, 'kN'], ['Tension at touchdown (dynamic)', stTd, 1e-3, 'kN'], ['Touchdown hot-spot stress', stS, 1e-6, 'MPa'], ['Vertical hang-off motion', stZ, 1, 'm'], ['Horizontal hang-off motion', stX, 1, 'm'], ['Vertical hang-off velocity', stV, 1, 'm/s'], ['Vertical hang-off acceleration', stA, 1, 'm/s²']]
        .map(([q, s, f, unit]) => ({ q, sigma: round(s.sigma * f, 4), Tz: round(s.Tz, 3), mpm: round(mpm(s.sigma, s.Tz) * f, 4), unit })),
      explain: 'Spectral statistics of each response in the sea state: standard deviation, mean zero-crossing period and the most probable maximum in three hours.' }];
    const notes = [`Hybrid method: ${inp.nFreq} transfer-function periods (${sw.source}), ${M} wave components, seed ${inp.seed}, ${nt} time steps of ${inp.dt} s after a ${fmt(ramp)} s ramp from static equilibrium.`,
      `Sea state Hs = ${inp.Hs} m, Tp = ${inp.Tp} s; static configuration includes current ${inp.Us} m/s at the surface.`];

    if (inp.method === 'lumped') {
      // opt-in nonlinear check: regular wave of height 1.86·Hs (most probable maximum in ~1000 waves) at Tp
      ctx.progress?.(0.5, 'Lumped-mass time-domain run');
      const H = 1.86 * inp.Hs, om = 2 * Math.PI / inp.Tp, r1 = hangoffRao(inp, [om]);
      const cat = base.catenary, L = cat.suspendedLength + inp.groundedLength, cf = currentFn(inp, inp.depth);
      const zs = [], us = []; for (let i = 0; i <= 10; i++) { zs.push(-inp.depth + inp.depth * i / 10); us.push(cf(inp.depth * i / 10)); }
      const model = run(() => buildLumpedLine({ anchor: [cat.horizontalSpan + inp.groundedLength, 0, -inp.depth], fairlead: [0, 0, -inp.hangoffDepth],
        sections: [{ L, EA: pp.section.EA, mass: pp.section.m, d: pp.section.D, nSeg: inp.nSeg, Cdn: inp.Cdn, Cdt: inp.Cdt, Can: inp.Ca, zeta: 0.05 }],
        seabed: { z: -inp.depth, mu: inp.mu }, rho: pp.rhoW, current: { z: zs, speed: us, direction: inp.curHeading * DEG }, wave: { H, T: inp.Tp, direction: 0, depth: inp.depth } }), 'Lumped-mass model');
      const dur = inp.nPeriods * inp.Tp;
      const sim = run(() => simulateLine(model, { motion: { type: 'harmonic', amplitude: [r1.surge[0] * H / 2, 0, r1.heave[0] * H / 2], period: inp.Tp, phase: [Math.PI / 2, 0, 0] }, duration: dur, dt: inp.Tp / 40, statsFrom: dur - inp.Tp, integrator: 'semi-implicit' }), 'Lumped-mass simulation');
      need(Number.isFinite(sim.stats.fairleadMax), 'The lumped-mass simulation became unstable; increase the number of segments or reduce the wave height.');
      kpis.push(kpi('lmTmax', 'Lumped-mass maximum top tension (regular wave)', sim.stats.fairleadMax / 1e3, 'kN', `Largest top tension over the last wave period of the nonlinear run with a regular wave of ${fmt(H)} m at Tp. It includes quadratic drag, seabed lift-off and direct wave loading that the hybrid method linearises.`),
        kpi('lmTmin', 'Lumped-mass minimum top tension (regular wave)', sim.stats.fairleadMin / 1e3, 'kN', 'Smallest top tension over the last wave period of the nonlinear run; compare with the hybrid 3-hour minimum.'));
      plots.push(plot('Lumped-mass top tension (regular wave)', 'Time', 's', 'Top tension', 'kN', [{ name: 'Nonlinear time domain', x: Array.from(sim.t), y: scaled(sim.fairleadTension, 1e-3) }], 'Top tension from the nonlinear lumped-mass run in the equivalent regular wave. A non-sinusoidal trace or sharp spikes indicate slack–snap behaviour that the linearised method cannot show.'));
      notes.push(`Lumped-mass check: ${inp.nSeg} segments, ${sim.steps} internal steps of ${fmt(sim.dtInternal * 1e3)} ms, regular wave H = ${fmt(H)} m, T = ${inp.Tp} s, current heading ${inp.curHeading}°.`);
    }

    const ucMax = Math.max(...checks.map((c) => c.uc));
    const cyclesYr = stT.Tz > 0 ? 31557600 / stT.Tz : 0;
    const advice = [
      ...vesselRangeAdvice(inp, inp.Tp),
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.85 ? 'warn' : 'ok',
        `In Hs = ${inp.Hs} m the top tension varies with σ = ${fmt(stT.sigma / 1e3)} kN; 3-hour extreme ${fmt((Ts + mpmT) / 1e3)} kN (${fmt((Ts + mpmT) / Ts)} × static). Touchdown stress σ = ${fmt(stS.sigma / 1e6)} MPa; minimum touchdown tension ${fmt((base.Te[kT] - mpmTd) / 1e3)} kN.`,
        checks[1].uc > 1 ? 'Compression at touchdown is predicted in this sea state. Run the nonlinear lumped-mass option (or a full time-domain FE model) to confirm, and raise the hang-off angle or reduce porch motion if it persists.'
          : checks[0].uc > 1 ? `Reduce the dynamic stress at the touchdown zone or increase the wall there to about ${fmt(pp.t * 1e3 * vmX / (0.9 * allow))} mm.`
            : 'Repeat with several seeds and with the governing 100-year and 1-year sea states; use at least the 3-hour most probable maximum, not the maximum of one short record, as the design value.',
        `DNV-ST-F201 (2018) App. A (dynamic analysis methods) and API RP 2RD §5.5. ${CONFIRM}`),
      adv('economic', 'info',
        `Each year in this sea state would apply about ${fmt(cyclesYr)} tension cycles. A fully nonlinear time-domain campaign over a fatigue scatter diagram takes orders of magnitude more computing than this hybrid screening.`,
        'Use the hybrid method to rank sea states and configurations, then spend nonlinear time-domain runs only on the few cases that govern extremes or show compression; that keeps the analysis budget on the cases that change the design.',
        'Value-of-information reasoning for analysis scope; estimate maturity per AACE RP 18R-97.',
        `If dynamics force a thicker riser, each extra millimetre of wall over ${fmt(base.s[n])} m is ${fmt(steelTonnes(pp.D, 0.001, base.s[n], pp.mat.rho))} t of steel. ${steelCo2Text(steelTonnes(pp.D, 0.001, base.s[n], pp.mat.rho))}`),
    ];
    return finish({
      kpis, checks, plots, tables, notes, advice,
      publish: { 'riser.dynTensionSigma': stT.sigma, 'riser.dynTensionMpm': Ts + mpmT },
      quantities: { steel_t: steelTonnes(pp.D, pp.t, base.s[n], pp.mat.rho) },
    });
  },

  verify: [
    {
      // Calm-sea limit: with a negligible sea state (Hs = 0.1 m at Tp = 4 s) and no current the response collapses to
      // static equilibrium, T = w·h/(1 − sinθ) for the inextensible catenary (h = 985 m, θ = 12°).
      name: 'Calm-sea limit: top tension → static catenary w·h/(1 − sin θ)', source: 'Limiting case; Irvine, Cable Structures (1981) §1.2',
      inputs: { Hs: 0.1, Tp: 4, Us: 0, Ub: 0, topType: 'pinned', pInt: 0, nElem: 320 },
      expect: [{ kpi: 'TtopStatic', value: handWeight() * 985 / (1 - Math.sin(12 * DEG)) / 1e3, tol: 0.005 }, { kpi: 'mpmT', value: handWeight() * 985 / (1 - Math.sin(12 * DEG)) / 1e3, tol: 0.01 }],
    },
  ],
};
