// Coupled dynamic response in the time domain: 3-D lumped-mass line driven by the host's hang-off motion (surge,
// sway, heave from RAOs) in a regular Airy wave with optional steady current from any direction.
import { buildLumpedLine, relaxLine, simulateLine, detectSnapLoads } from '../../engine/lines/lumpedmass.js';
import { buoyancyModules } from '../../engine/xsection/ancillaries.js';
import { lineInputs, configInputs, lineFrom, adv, check, worst, sev, fmt, round, run, need, thin, DEG, LIMIT_RIGID, STD_UMB, CO2_STEEL } from './_common.js';

export default {
  id: 's12-coupled-dynamics', suite: 12, title: 'Coupled dynamic response (time domain, 3-D)',
  summary: 'Tension, motion, velocity and acceleration histories of the umbilical under combined wave loading and host motion, with dynamic amplification and slack screening.',
  covers: [2], fidelity: 'numerical',
  standards: [`${STD_UMB} §7.3 (dynamic global analysis)`, 'DNV-RP-C205 (2021) §3 and §6 (Airy kinematics, Morison load)', 'DNV-RP-N103 (2021) §4.7 (slack / snap criterion)'],
  models: ['Dynamic cable finite elements (3-D lumped mass; Hall & Goupee 2015)', "Morison's equation with relative velocity", 'Airy wave kinematics', 'Vessel-motion response by RAO', 'Seabed spring–damper contact with Coulomb friction', 'Three-dimensional cable equilibrium'],
  limits: 'Regular wave and single-harmonic host motion; the host is not affected by the umbilical (one-way coupling, adequate for a single umbilical on a floater). Bending enters only as a discrete penalty, so curvature is not reported here (see the vessel-motion module); torsion and stick–slip hysteresis are not modelled. Irregular-sea extreme statistics need a 3-hour simulation in a dedicated global-analysis tool.',
  publishes: ['umbilical.dynTensionMax', 'umbilical.daf'],
  searchTerms: 'umbilical dynamic analysis lumped mass time domain',
  inputs: [
    ...lineInputs(),
    ...configInputs().filter((i) => i.key !== 'kn'),
    { key: 'kn', label: 'Seabed vertical stiffness', unit: 'kPa', type: 'number', default: 200, min: 5, max: 1e5, group: 'Seabed', role: 'soil', help: 'Force per metre of umbilical per metre of penetration.' },
    { key: 'H', label: 'Wave height (regular)', unit: 'm', type: 'number', default: 6, min: 0, max: 30, group: 'Wave and current', role: 'bc', bind: 'site.Hs1', uncertain: { dist: 'lognormal', cov: 0.1 },
      help: 'Regular design wave. For a sea state use Hmax ≈ 1.86·Hs with an associated period.' },
    { key: 'Tw', label: 'Wave period', unit: 's', type: 'number', default: 10, min: 3, max: 25, group: 'Wave and current', role: 'bc', bind: 'site.Tp1' },
    { key: 'waveDir', label: 'Wave direction from the umbilical plane', unit: 'deg', type: 'number', default: 0, min: -180, max: 180, group: 'Wave and current', role: 'bc', help: '0 = travelling from hang-off towards the termination.' },
    { key: 'Uc', label: 'Steady current speed (uniform)', unit: 'm/s', type: 'number', default: 0, min: 0, max: 2.5, group: 'Wave and current', role: 'bc',
      help: 'Opt-in: a cross current lengthens the static relaxation considerably (several seconds per run). Zero by default.' },
    { key: 'curDir', label: 'Current direction from the umbilical plane', unit: 'deg', type: 'number', default: 90, min: -180, max: 180, group: 'Wave and current', role: 'bc' },
    { key: 'raoSurge', label: 'Hang-off surge RAO', unit: 'm/m', type: 'number', default: 0.6, min: 0, max: 5, group: 'Host motion at the hang-off', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 },
      help: 'Motion amplitude per unit wave amplitude at the hang-off (including pitch × lever arm), at the wave period.' },
    { key: 'raoSway', label: 'Hang-off sway RAO', unit: 'm/m', type: 'number', default: 0.3, min: 0, max: 5, group: 'Host motion at the hang-off', role: 'bc' },
    { key: 'raoHeave', label: 'Hang-off heave RAO', unit: 'm/m', type: 'number', default: 0.9, min: 0, max: 5, group: 'Host motion at the hang-off', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'phHeave', label: 'Heave phase lead over surge', unit: 'deg', type: 'number', default: 90, min: -180, max: 180, group: 'Host motion at the hang-off', role: 'bc' },
    { key: 'monitor', label: 'Monitored point (fraction of length from hang-off)', unit: '-', type: 'number', default: 0.55, min: 0.02, max: 0.98, group: 'Output', role: 'control', help: 'Node whose displacement, velocity and acceleration histories are plotted (0.55 ≈ buoyancy section of the default wave).' },
    { key: 'Tallow', label: 'Allowable tension', unit: 'kN', type: 'number', default: 350, min: 5, max: 5000, group: 'Acceptance criteria', role: 'criteria', bind: 'umbilical.allowTension' },
    { key: 'slackFrac', label: 'Slack threshold (fraction of local static tension)', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.5, group: 'Acceptance criteria', role: 'criteria', help: 'Minimum dynamic tension must stay above this fraction of the static value (slack-line criterion in the style of DNV-RP-N103).' },
    { key: 'zeta', label: 'Axial structural damping ratio', unit: '-', type: 'number', default: 0.05, min: 0, max: 0.3, group: 'Solver', role: 'ic', help: 'Damping of each segment\'s own axial mode; removes numerical ringing without affecting wave-frequency response.' },
    { key: 'duration', label: 'Simulated time', unit: 'wave periods', type: 'number', default: 4, min: 3, max: 40, group: 'Solver', role: 'control', help: 'The first half is start-up (ramp); statistics use the second half. Longer runs are opt-in.' },
    { key: 'nSeg', label: 'Line segments', type: 'int', default: 40, min: 12, max: 400, group: 'Solver', role: 'control' },
    { key: 'cfl', label: 'Time step / critical step', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 0.9, group: 'Solver', role: 'control', help: 'Explicit integration sub-step as a fraction of the axial-wave stability limit.' },
  ],
  mesh: { key: 'nSeg', label: 'Line segments', levels: [20, 40, 80], qoi: 'Tmax' },
  meshTime: { key: 'cfl', label: 'Time step / critical step', levels: [0.125, 0.25, 0.5], qoi: 'Tmax', isSize: true },
  compute(inp, ctx) {
    const line = lineFrom(inp);
    need(line.w > 0, 'The umbilical is buoyant in sea water: it cannot hang as a catenary.');
    const h = inp.depth - inp.hangoffDepth, L = inp.length;
    need(h > 10, 'Hang-off depth must leave at least 10 m between the hang-off and the seabed.');
    need(L > Math.hypot(inp.span, h) * 1.002, `Under-constrained geometry: the length (${L} m) is shorter than the straight distance hang-off to termination (${fmt(Math.hypot(inp.span, h), 5)} m).`);
    const f = inp.nSeg / L;
    const mk = (len, mass, d) => ({ L: len, EA: line.EA, mass, d, nSeg: Math.max(2, Math.round(len * f)), Cdn: inp.Cdn, Cdt: inp.Cdt, Can: inp.Ca, zeta: inp.zeta, EI: line.EI });
    let sections;
    if (inp.config === 'lazy-wave') {
      const lower = L - inp.upperLength - inp.buoyLength;
      need(lower > 5, 'Over-constrained lengths: upper + buoyancy section leave no lower catenary within the total length.');
      const b = run(() => buoyancyModules({ outerDiameter: inp.modOD, innerDiameter: line.D, moduleLength: inp.modL, foamDensity: inp.foamRho, rhoW: line.rhoW }), 'Buoyancy module');
      const AeB = line.Ae + b.foamVolume / inp.modPitch;       // smeared displaced area → volume-equivalent diameter
      sections = [mk(lower, line.m, line.D), mk(inp.buoyLength, line.m + inp.foamRho * b.foamVolume / inp.modPitch, Math.sqrt(4 * AeB / Math.PI)), mk(inp.upperLength, line.m, line.D)];
    } else sections = [mk(L, line.m, line.D)];
    const a = inp.H / 2, per = inp.Tw, om = 2 * Math.PI / per;
    const model = run(() => buildLumpedLine({
      anchor: [inp.span, 0, -inp.depth], fairlead: [0, 0, -inp.hangoffDepth], sections, rho: line.rhoW,
      seabed: { z: -inp.depth, k: inp.kn * 1e3 / line.D, mu: inp.muAx },
      current: inp.Uc > 0 ? { speed: inp.Uc, direction: inp.curDir * DEG } : undefined,
      wave: inp.H > 0 ? { H: inp.H, T: per, direction: inp.waveDir * DEG, depth: inp.depth } : undefined, waterLevel: 0,
    }), 'The line model could not be built (check length, span and buoyancy)');
    const st = run(() => relaxLine(model, { tol: 1e-4, maxIter: inp.Uc > 0 ? 400000 : 120000 }), 'Static relaxation failed');
    need(Number.isFinite(st.fairleadTension), 'The static state is not finite: the configuration is not feasible. Check that length, span, depth and buoyancy describe a shape that can hang between the two ends.');
    const N = model.N, nn = N + 1;
    const sStatic = Float64Array.from(model.s), r0 = Float64Array.from(model.r), Tst = Float64Array.from(st.segmentTension);
    // node index from the hang-off fraction (nodes are numbered from the anchor)
    const iMon = Math.min(N - 1, Math.max(1, Math.round((1 - inp.monitor) * N)));
    const amp = [inp.raoSurge * a, inp.raoSway * a, inp.raoHeave * a];
    const motion = { type: 'harmonic', amplitude: amp, period: per, phase: [0, 0, inp.phHeave * DEG], ramp: per };
    const dur = inp.duration * per, nOut = Math.min(600, Math.max(120, Math.round(dur / (per / 40)))), dt = dur / nOut, from = dur / 2;
    const t = [], Tf = [], Ta = [], dx = [], dy = [], dz = [], vel = [], acc = [];
    const env = { Tmax: new Float64Array(N).fill(-Infinity), Tmin: new Float64Array(N).fill(Infinity), yAbs: new Float64Array(nn), rAbs: new Float64Array(nn) };
    let tdMin = Infinity, tdMax = -Infinity, vPrev = [0, 0, 0];
    // advance in output-interval chunks so that node kinematics can be sampled between chunks
    for (let k = 1; k <= nOut; k++) {
      const res = run(() => simulateLine(model, { motion, duration: dt, dt, integrator: 'semi-implicit', cfl: inp.cfl, quasiStatic: false, statsFrom: 0 }), 'Time integration failed');
      const tk = k * dt, j = 3 * iMon, v = [model.v[j], model.v[j + 1], model.v[j + 2]];
      const TfK = res.fairleadTension[res.fairleadTension.length - 1];
      if (!Number.isFinite(TfK)) throw new Error('The time integration became unstable. Reduce the time-step fraction or increase the axial damping ratio.');
      t.push(tk); Tf.push(TfK / 1e3); Ta.push(res.anchorTension[res.anchorTension.length - 1] / 1e3);
      dx.push(model.r[j] - r0[j]); dy.push(model.r[j + 1] - r0[j + 1]); dz.push(model.r[j + 2] - r0[j + 2]);
      vel.push(Math.hypot(...v)); acc.push(Math.hypot(v[0] - vPrev[0], v[1] - vPrev[1], v[2] - vPrev[2]) / dt); vPrev = v;
      if (tk >= from) {
        for (let s = 0; s < N; s++) { env.Tmax[s] = Math.max(env.Tmax[s], res.envelopes.Tmax[s]); env.Tmin[s] = Math.min(env.Tmin[s], res.envelopes.Tmin[s]); }
        for (let i = 0; i < nn; i++) {
          env.yAbs[i] = Math.max(env.yAbs[i], Math.abs(res.envelopes.yMin[i] - r0[3 * i + 1]), Math.abs(res.envelopes.yMax[i] - r0[3 * i + 1]));
          env.rAbs[i] = Math.max(env.rAbs[i], Math.hypot(model.r[3 * i] - r0[3 * i], model.r[3 * i + 1] - r0[3 * i + 1], model.r[3 * i + 2] - r0[3 * i + 2]));
        }
        tdMin = Math.min(tdMin, res.stats.touchdownMin); tdMax = Math.max(tdMax, res.stats.touchdownMax);
      }
    }
    const i0 = t.findIndex((v) => v >= from);
    const win = (arr) => arr.slice(i0);
    const Tstat = st.fairleadTension / 1e3, Tmax = Math.max(...win(Tf)), Tmin = Math.min(...win(Tf));
    const TaMax = Math.max(...win(Ta));
    // lowest tension anywhere relative to its static value: slack screening per segment
    let slackWorst = Infinity, sSlack = 0, TminLine = Infinity;
    for (let s = 0; s < N; s++) { const ratio = env.Tmin[s] / Math.max(Tst[s], 1); if (ratio < slackWorst) { slackWorst = ratio; sSlack = L - 0.5 * (sStatic[s] + sStatic[s + 1]); } TminLine = Math.min(TminLine, env.Tmin[s]); }
    const snap = detectSnapLoads({ t, T: Tf.map((v) => v * 1e3), Tref: st.fairleadTension, slackFraction: inp.slackFrac });
    const daf = Tmax / Tstat, dynAmp = 0.5 * (Tmax - Tmin);
    const vMax = Math.max(...win(vel)), aMax = Math.max(...win(acc)), yMax = Math.max(...env.yAbs), rMax = Math.max(...env.rAbs);
    const checks = [
      check('Maximum dynamic tension at the hang-off', Tmax, inp.Tallow, 'kN', `${STD_UMB} §7.3 (tension within the capacity curve)`, 'Guards against overstressing tubes and armour at the peak of the wave cycle; exceeding it yields the most loaded element.'),
      check('Slack (minimum tension ratio)', inp.slackFrac, Math.max(slackWorst, 1e-6), '-', 'DNV-RP-N103 (2021) §4.7 slack criterion (threshold is a project input)', 'Demand is the required minimum ratio of dynamic to static tension, capacity the lowest ratio found along the line. Below it the line goes slack and re-tensions with a snap load and a compression kink.'),
    ];
    const gov = worst(checks);
    const sFromTop = Array.from({ length: N }, (_, s) => L - 0.5 * (sStatic[s] + sStatic[s + 1]));
    const [th, TfH, TaH, dxh, dyh, dzh, vh, ah] = thin([t, Tf, Ta, dx, dy, dz, vel, acc], 400);
    const nodes = [], lines = [], vals = [], deformed = [];
    for (let i = 0; i < nn; i++) {
      nodes.push([r0[3 * i], r0[3 * i + 1], r0[3 * i + 2]]); deformed.push([model.r[3 * i], model.r[3 * i + 1], model.r[3 * i + 2]]);
      vals.push((env.Tmax[Math.min(i, N - 1)]) / 1e3); if (i) lines.push([i - 1, i]);
    }
    return {
      kpis: [
        { key: 'Tstatic', label: 'Static hang-off tension', value: round(Tstat, 6), unit: 'kN', explain: 'Top tension of the relaxed static shape (with current if specified): the reference for the dynamic amplification.' },
        { key: 'Tmax', label: 'Maximum dynamic hang-off tension', value: round(Tmax, 6), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Peak top tension in the steady part of the simulation; it must stay inside the capacity curve of the cross-section.' },
        { key: 'Tmin', label: 'Minimum dynamic hang-off tension', value: round(Tmin), unit: 'kN', explain: 'Trough of the top tension cycle. A value near zero warns of slack and snap loading.' },
        { key: 'daf', label: 'Dynamic amplification Tmax / Tstatic', value: round(daf), unit: '-', explain: 'Ratio of peak to static top tension. Values well above 1.3 show that drag locks the line and hang-off heave is taken by axial stretch.' },
        { key: 'dynAmp', label: 'Dynamic tension amplitude', value: round(dynAmp), unit: 'kN', explain: 'Half the tension range per wave cycle: the input to tension fatigue of tubes and armour at the hang-off.' },
        { key: 'TtermMax', label: 'Maximum tension at the termination', value: round(TaMax), unit: 'kN', explain: 'Peak load arriving at the subsea termination; a lazy wave should keep it close to the static value.' },
        { key: 'slackRatio', label: 'Lowest dynamic / static tension ratio', value: round(slackWorst, 3), unit: '-', limit: inp.slackFrac, sense: 'min', explain: 'Smallest ratio of minimum dynamic to static tension along the line. Below the threshold the segment goes slack; negative means compression.' },
        { key: 'tdTravel', label: 'Touchdown travel per cycle', value: round(Number.isFinite(tdMax - tdMin) ? tdMax - tdMin : 0, 3), unit: 'm', explain: 'Length of line that lifts off and lands each wave cycle: the zone of sheath abrasion and trenching.' },
        { key: 'yMax', label: 'Maximum out-of-plane excursion', value: round(yMax, 3), unit: 'm', explain: 'Largest sideways displacement from the static plane; it consumes clearance to neighbouring risers and moorings.' },
        { key: 'rMax', label: 'Maximum node displacement from static', value: round(rMax, 3), unit: 'm', explain: 'Largest total dynamic displacement of any node, normally in the upper catenary or buoyant section.' },
        { key: 'vMax', label: 'Peak velocity at the monitored point', value: round(vMax, 3), unit: 'm/s', explain: 'Structural velocity at the monitored node; it sets the relative-velocity drag and the Keulegan–Carpenter regime there.' },
        { key: 'aMax', label: 'Peak acceleration at the monitored point', value: round(aMax, 3), unit: 'm/s²', explain: 'Structural acceleration at the monitored node; multiplied by clamp mass it gives the inertia load on buoyancy-module clamps and ancillaries.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Tension histories', x: { label: 'Time', unit: 's' }, y: { label: 'Effective tension', unit: 'kN' },
          series: [{ name: 'Hang-off', x: th, y: TfH }, { name: 'Termination', x: th, y: TaH }, { name: 'Static hang-off', x: [0, dur], y: [Tstat, Tstat], dash: true }],
          explain: 'Top and termination tension through the run. The first half is the start-up ramp; a steady periodic response in the second half means the statistics are representative. Sharp spikes after near-zero troughs are snap loads.' },
        { type: 'line', title: 'Displacement history at the monitored point', x: { label: 'Time', unit: 's' }, y: { label: 'Displacement from static', unit: 'm' },
          series: [{ name: 'In-plane horizontal', x: th, y: dxh }, { name: 'Out-of-plane', x: th, y: dyh }, { name: 'Vertical', x: th, y: dzh }],
          explain: 'Three components of motion of the monitored node relative to its static position. Compare the amplitudes with the hang-off motion to see how much of it reaches that point.' },
        { type: 'line', title: 'Velocity and acceleration at the monitored point', x: { label: 'Time', unit: 's' }, y: { label: 'Magnitude', unit: 'm/s, m/s²' },
          series: [{ name: 'Velocity [m/s]', x: th, y: vh }, { name: 'Acceleration [m/s²]', x: th, y: ah }],
          explain: 'Magnitudes of structural velocity and acceleration at the monitored node. Peaks give the inertia loads on clamps and the relative-velocity drag.' },
        { type: 'line', title: 'Tension envelope along the line', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Effective tension', unit: 'kN' },
          series: [{ name: 'Maximum', x: sFromTop, y: Array.from(env.Tmax, (v) => v / 1e3) }, { name: 'Static', x: sFromTop, y: Array.from(Tst, (v) => v / 1e3), dash: true }, { name: 'Minimum', x: sFromTop, y: Array.from(env.Tmin, (v) => v / 1e3) }],
          explain: 'Range of tension at every segment during the steady part. Where the minimum curve approaches zero the line is close to slack; the widest band marks the tension-fatigue hot spot.' },
        { type: 'model3d', title: 'Static shape and final dynamic position', nodes, lines, tris: [], nodeValues: vals, valueLabel: 'Maximum tension', valueUnit: 'kN', deformed,
          explain: 'Three-dimensional static configuration (z = 0 at still water) coloured by maximum dynamic tension, with the position at the end of the run as the deformed shape. Out-of-plane displacement shows the effect of sway and oblique waves or current.' },
      ],
      tables: [{ title: 'Dynamic response summary', explain: 'Static, maximum and minimum tension at the two ends and the slack screening result. Use the hang-off row for the capacity check and the termination row for the connector and foundation loads.',
        columns: [{ key: 'loc', label: 'Location' }, { key: 'st', label: 'Static', unit: 'kN' }, { key: 'mx', label: 'Maximum', unit: 'kN' }, { key: 'mn', label: 'Minimum', unit: 'kN' }, { key: 'daf', label: 'Max / static', unit: '-' }],
        rows: [{ loc: 'Hang-off', st: round(Tstat), mx: round(Tmax), mn: round(Tmin), daf: round(daf) }, { loc: 'Termination', st: round(st.anchorTension / 1e3), mx: round(TaMax), mn: round(Math.min(...win(Ta))), daf: round(TaMax / Math.max(st.anchorTension / 1e3, 1e-9)) }] }],
      notes: [
        `Hang-off motion amplitudes: surge ${fmt(amp[0])} m, sway ${fmt(amp[1])} m, heave ${fmt(amp[2])} m at ${per} s. Statistics from the last ${fmt(dur - from)} s of ${fmt(dur)} s.`,
        `Lowest dynamic / static tension ratio ${fmt(slackWorst)} at ${fmt(sSlack, 4)} m from the hang-off${snap.snap ? '; a snap (slack → re-tension) event was detected at the hang-off' : ''}.`,
        st.converged ? 'Static relaxation converged.' : `Static relaxation stopped at residual ${fmt(st.residual)}: the steady current had not fully deflected the line, so the early part of the run still drifts. Increase the simulated time.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Dynamic amplification is ${fmt(daf)} (peak ${fmt(Tmax)} kN against static ${fmt(Tstat)} kN); ${gov.name.toLowerCase()} governs at utilisation ${fmt(gov.uc)}.`,
          gov.uc > 1 ? 'Reduce the heave transmitted to the line: move the hang-off closer to the vessel motion centre, add buoyancy to form or enlarge the wave, or restrict the operating sea state; then repeat with an irregular 3-hour sea in a full global analysis.' : 'Confirm the extreme with an irregular sea (3-hour, several seeds) and directional cases; use the tension amplitude and cycle count in the fatigue module.',
          `${STD_UMB} §7.3; DNV-RP-C205 (2021) §6; DNV-RP-N103 (2021) §4.7`),
        adv('economic', 'info', `Tension range ${fmt(2 * dynAmp)} kN per ${per} s cycle is about ${fmt(31557600 / per / 1e6)} million cycles a year if this sea state persisted.`,
          'Fatigue damage rises with the 3rd–4th power of this range: a configuration change that halves the range extends fatigue life roughly 8–16 times and can remove a mid-life replacement campaign (umbilical supply plus a lay-vessel mobilisation).',
          'Whole-life costing ISO 15686-5; risk-based inspection planning DNV-RP-G101',
          `Avoiding one mid-life replacement of a ${fmt(L)} m dynamic section saves its full embodied material (≈ ${fmt(L * 31 / 1e3)} t steel ≈ ${fmt(L * 31 / 1e3 * CO2_STEEL)} t CO₂e at the screening factor ${CO2_STEEL} t/t) and the vessel fuel of the campaign.`),
      ],
      publish: { 'umbilical.dynTensionMax': Tmax * 1e3, 'umbilical.daf': daf },
      quantities: { cable_km: L / 1e3 },
    };
  },
  verify: [
    // No wave, no host motion: the dynamic model must stay at the static catenary. Inextensible catenary with seabed
    // contact (Irvine §1.2), h = 340 m, L = 600 m, X = 330 m → a = 32.00 m, Ttop = w(a + h) = 235.11 × 372.0 = 87.46 kN.
    { name: 'Still water, no host motion → static catenary top tension, DAF = 1',
      source: 'Limiting case; Irvine, Cable Structures (1981) §1.2 (Ttop = w·(a + h), a = 32.00 m)',
      inputs: { config: 'catenary', length: 600, span: 330, H: 0, raoSurge: 0, raoSway: 0, raoHeave: 0, EI: 200, muAx: 0 },
      expect: [{ kpi: 'Tstatic', value: 87.46, tol: 0.015 }, { kpi: 'daf', value: 1, tol: 0.01 }] },
  ],
};
