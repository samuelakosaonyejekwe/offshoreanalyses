// Extreme storm and vessel-motion response of a dynamic cable: 100-year current with near, mean and far platform
// offsets on the nonlinear rod, then the short-term extreme of tension and curvature in the storm sea state from
// frequency-domain transfer functions (hang-off motion through the platform RAOs, drag linearised) integrated over
// a JONSWAP spectrum, with a regular design-wave cross-check.
import { rodHarmonicResponse } from '../../engine/lines/rod.js';
import { riserOffsetCases } from '../../engine/lines/riser.js';
import { jonswap, jonswapGamma, tabulatedMoments, spectralParameters } from '../../engine/metocean/spectra.js';
import { lineInputs, configInputs, currentInputs, solveConfig, rodCurrent, adv, check, worst, sev, fmt, round, run, need, thin, linspace, DEG, LIMIT_2D, STD_CABLE, CO2_COPPER, vesselCo2Text } from './_common.js';

// Indicative hang-off motion RAOs of a semi-submersible floating-wind platform (heave resonance near 20 s).
const DEFAULT_RAO = [
  { T: 4, heave: 0.02, surge: 0.05 }, { T: 6, heave: 0.08, surge: 0.15 }, { T: 8, heave: 0.2, surge: 0.3 }, { T: 10, heave: 0.35, surge: 0.45 }, { T: 12, heave: 0.5, surge: 0.58 },
  { T: 14, heave: 0.62, surge: 0.68 }, { T: 16, heave: 0.75, surge: 0.75 }, { T: 18, heave: 0.9, surge: 0.8 }, { T: 20, heave: 1.1, surge: 0.85 }, { T: 25, heave: 1.0, surge: 0.9 },
];

function interp(rows, key, T) {
  if (T <= rows[0].T) return rows[0][key];
  for (let i = 1; i < rows.length; i++) if (T <= rows[i].T) return rows[i - 1][key] + (rows[i][key] - rows[i - 1][key]) * (T - rows[i - 1].T) / (rows[i].T - rows[i - 1].T);
  return rows[rows.length - 1][key];
}
const trapz = (x, y) => { let s = 0; for (let i = 1; i < x.length; i++) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]); return s; };

export default {
  id: 's13-extreme-storm-response', suite: 13, title: 'Extreme storm and vessel-motion response',
  summary: 'Extreme tension, curvature and compression of a dynamic cable in the design storm: extreme current and platform offsets with the most probable maximum wave-frequency response from the platform motions.',
  covers: [15], fidelity: 'reduced-order',
  standards: [`${STD_CABLE} (extreme load cases of dynamic cables)`, 'CIGRE TB 862 (2022) §4 (dynamic cable design loads)', 'DNV-RP-C205 (2021) §3.5 (JONSWAP spectrum, short-term extremes) and §6.2 (Morison drag, linearisation)', 'DNV-RP-F205 (2019) (coupled analysis; vessel motion as input to slender structures)', 'API RP 2RD / API RP 17B §5 (near, far and cross offset cases)'],
  models: ['Elastic rod with geometric nonlinearity for the offset static shapes', 'Vessel motion analysis: hang-off motion from platform RAOs', 'JONSWAP wave spectrum', 'Frequency-domain dynamic tension and curvature transfer functions about each static shape', 'Morison drag with harmonic linearisation at the significant motion amplitude', 'Wave–current interaction through the static current-loaded shape', 'Rayleigh short-term extreme: most probable maximum σ·√(2·ln N)', 'Regular design-wave cross-check with Hmax'],
  limits: `${LIMIT_2D} The wave-frequency response is linear about each static shape with drag linearised at the significant motion level — not a full stochastic linearisation — and extremes are Rayleigh (narrow-band, Gaussian). Wave kinematics acting directly on the cable, slow-drift motion beyond the entered offsets, snatch after slack, out-of-plane seas and the non-Gaussian tension tail need a nonlinear irregular-wave time-domain analysis, which remains the design method for final verification. The default RAOs are indicative.`,
  publishes: ['cable.stormTension', 'cable.stormMinBendRadius', 'cable.stormTensionSig'],
  searchTerms: 'dynamic cable extreme response 100-year storm floating wind most probable maximum tension',
  inputs: [
    ...lineInputs(), ...configInputs(), ...currentInputs({ Us: 1.2, Ub: 0.4, bind: 'site.Uc100', group: 'Storm current' }),
    { key: 'Hs', label: 'Significant wave height of the design storm', unit: 'm', type: 'number', default: 10, min: 0, max: 20, group: 'Storm sea state', role: 'bc', bind: 'site.Hs100', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 14, min: 4, max: 25, group: 'Storm sea state', role: 'bc', bind: 'site.Tp100', uncertain: { dist: 'normal', cov: 0.08 }, help: 'Check the range of periods on the 100-year contour: a shorter period with the same height is often worse for a cable.' },
    { key: 'gamma', label: 'JONSWAP peak-enhancement factor', unit: '-', type: 'number', default: 0, min: 0, max: 7, group: 'Storm sea state', role: 'bc', help: '0 = from Hs and Tp per DNV-RP-C205 §3.5.5.' },
    { key: 'duration', label: 'Storm duration', unit: 'h', type: 'number', default: 3, min: 0.5, max: 12, group: 'Storm sea state', role: 'bc' },
    { key: 'offset', label: 'Platform offset in the storm (near and far)', unit: 'm', type: 'number', default: 18, min: 0, max: 400, group: 'Platform motion', role: 'bc', bind: 'vessel.offsetMax', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Mean plus slow-drift excursion of the hang-off from the mooring analysis (intact or one line broken, as the case requires).' },
    { key: 'rao', label: 'Hang-off motion RAOs', type: 'table', group: 'Platform motion', role: 'bc', columns: [{ key: 'T', label: 'Period', unit: 's', type: 'number' }, { key: 'heave', label: 'Heave', unit: 'm/m', type: 'number' }, { key: 'surge', label: 'Surge', unit: 'm/m', type: 'number' }],
      default: DEFAULT_RAO.map((r) => ({ ...r })), help: 'Motion of the hang-off point per metre of wave amplitude (platform heave and surge combined with pitch × lever arm), ascending periods. Values outside the table are held at the end values.' },
    { key: 'phase', label: 'Phase of heave relative to surge', unit: 'deg', type: 'number', default: 90, min: -180, max: 180, group: 'Platform motion', role: 'bc' },
    { key: 'zeta', label: 'Structural damping ratio', unit: '-', type: 'number', default: 0.03, min: 0, max: 0.3, group: 'Platform motion', role: 'material' },
    { key: 'v0', label: 'Initial velocity of the cable', unit: 'm/s', type: 'number', default: 0, min: 0, max: 0, group: 'Platform motion', role: 'ic', help: 'Stationary response: the cable starts from its static equilibrium in the storm current; start-up transients are not computed.' },
    { key: 'Tallow', label: 'Allowable tension', unit: 'kN', type: 'number', default: 250, min: 2, max: 5000, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.allowTension' },
    { key: 'mbrAllow', label: 'Allowable minimum bend radius (dynamic)', unit: 'm', type: 'number', default: 3.2, min: 0.3, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.mbrDynamic' },
    { key: 'nFreq', label: 'Frequencies for the transfer functions', type: 'int', default: 8, min: 3, max: 40, group: 'Solver', role: 'control', help: 'Number of wave frequencies solved between 0.45·Tp and 1.8·Tp.' },
    { key: 'nElem', label: 'Rod elements', type: 'int', default: 80, min: 30, max: 600, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Rod elements', levels: [40, 80, 160], qoi: 'Tmax' },
  meshTime: { key: 'nFreq', label: 'Frequencies', levels: [4, 8, 16], qoi: 'tensionSig' },
  compute(inp) {
    const rao = (inp.rao || []).filter((r) => r.T > 0).map((r) => ({ T: Number(r.T), heave: Math.max(0, Number(r.heave) || 0), surge: Math.max(0, Number(r.surge) || 0) })).sort((a, b) => a.T - b.T);
    need(rao.length >= 1, 'The RAO table needs at least one row with a positive period.');
    const { res, line } = solveConfig(inp, { nElem: inp.nElem, current: rodCurrent(inp, inp.rhoW) });
    const cases = inp.offset > 0 ? run(() => riserOffsetCases(res, [-inp.offset, inp.offset]), 'Storm offset cases did not solve (the line is too short for the far offset or piles up in the near offset)') : [];
    const shapes = [{ label: 'Near', r: cases[0]?.result }, { label: 'Mean', r: res }, { label: 'Far', r: cases[1]?.result }].filter((c) => c.r);
    const wave = inp.config === 'lazy-wave';
    // sea state on a fine grid: wave and hang-off motion statistics (no rod solves needed)
    const wp = 2 * Math.PI / inp.Tp, storm = inp.Hs > 0;
    const gam = inp.gamma > 0 ? inp.gamma : (storm ? jonswapGamma(inp.Hs, inp.Tp) : 1);
    const wf = linspace(0.3 * wp, 5 * wp, 400), Sf = wf.map((w) => (storm ? jonswap(w, { Hs: inp.Hs, Tp: inp.Tp, gamma: gam }) : 0));
    const seaPar = storm ? spectralParameters(tabulatedMoments(wf, Sf)) : { Tz: inp.Tp / 1.3, Hm0: 0 };
    const dur = inp.duration * 3600, Nw = Math.max(dur / seaPar.Tz, 2), Hmax = inp.Hs * Math.sqrt(0.5 * Math.log(Nw));
    const hv = wf.map((w) => interp(rao, 'heave', 2 * Math.PI / w)), sg = wf.map((w) => interp(rao, 'surge', 2 * Math.PI / w));
    const heaveSig = Math.sqrt(trapz(wf, Sf.map((s, i) => s * hv[i] ** 2))), surgeSig = Math.sqrt(trapz(wf, Sf.map((s, i) => s * sg[i] ** 2)));
    const mpm = (m0, m2) => { const tz = m2 > 0 ? 2 * Math.PI * Math.sqrt(m0 / m2) : inp.Tp; return { tz, c: Math.sqrt(2 * Math.log(Math.max(dur / tz, 2))) }; };
    const heaveM2 = trapz(wf, Sf.map((s, i) => s * hv[i] ** 2 * wf[i] ** 2));
    const heaveMpm = heaveSig * mpm(heaveSig ** 2, heaveM2).c;
    // transfer functions about each static shape, drag linearised at the significant wave amplitude
    const aLin = inp.Hs / 2, wg = linspace(2 * Math.PI / (1.8 * inp.Tp), 2 * Math.PI / Math.max(0.45 * inp.Tp, 2.5), inp.nFreq);
    const Sg = wg.map((w) => (storm ? jonswap(w, { Hs: inp.Hs, Tp: inp.Tp, gamma: gam }) : 0));
    const moving = storm && rao.some((r) => r.heave + r.surge > 0);
    const dynOf = (r) => {
      const n = r.s.length, m0T = new Float64Array(n), m0K = new Float64Array(n), HT = [];
      let m2Top = 0;
      if (moving) {
        const H = wg.map((w) => {
          const T = 2 * Math.PI / w, xa = interp(rao, 'surge', T) * aLin, za = interp(rao, 'heave', T) * aLin;
          if (!(xa + za > 0)) return null;
          const h = run(() => rodHarmonicResponse(r, { omega: w, end: 'A', xAmp: xa, zAmp: za, phaseZ: inp.phase * DEG, rho: line.rhoW, zeta: inp.zeta }), `Harmonic response at ${fmt(T)} s`);
          need(h.converged !== false, `The drag linearisation did not converge at ${fmt(T)} s: increase the damping or reduce the wave height.`);
          return h;
        });
        // trapezoidal integration of |H|²·S over the solved frequencies, node by node
        for (let j = 0; j < wg.length; j++) {
          const wt = (j === 0 ? wg[1] - wg[0] : j === wg.length - 1 ? wg[j] - wg[j - 1] : wg[j + 1] - wg[j - 1]) / 2;
          const h = H[j]; HT.push(h ? h.tensionAmp[0] / aLin : 0);
          if (!h) continue;
          for (let i = 0; i < n; i++) { m0T[i] += (h.tensionAmp[i] / aLin) ** 2 * Sg[j] * wt; m0K[i] += (h.curvatureAmp[i] / aLin) ** 2 * Sg[j] * wt; }
          m2Top += (h.tensionAmp[0] / aLin) ** 2 * Sg[j] * wt * wg[j] ** 2;
        }
      } else for (let j = 0; j < wg.length; j++) HT.push(0);
      const st = mpm(m0T[0], m2Top), c = moving ? st.c : 0;
      let kTot = 0, sK = 0, slack = 0, tMin = Infinity;
      const tEnvMax = new Float64Array(n), tEnvMin = new Float64Array(n), kEnv = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        const dT = c * Math.sqrt(m0T[i]), dK = c * Math.sqrt(m0K[i]);
        tEnvMax[i] = r.Te[i] + dT; tEnvMin[i] = r.Te[i] - dT; kEnv[i] = Math.abs(r.curvature[i]) + dK;
        if (kEnv[i] > kTot) { kTot = kEnv[i]; sK = r.s[i]; }
        slack = Math.max(slack, dT / Math.max(r.Te[i], 1)); tMin = Math.min(tMin, tEnvMin[i]);
      }
      return { sigTop: Math.sqrt(m0T[0]), c, tz: st.tz, Tmax: r.topTension + c * Math.sqrt(m0T[0]), kTot, sK, slack, tMin, tEnvMax, tEnvMin, kEnv, HT };
    };
    const rows = shapes.map((c) => ({ ...c, d: dynOf(c.r) }));
    const far = rows[rows.length - 1], gT = rows.reduce((a, q) => (q.d.Tmax > a.d.Tmax ? q : a), rows[0]), gK = rows.reduce((a, q) => (q.d.kTot > a.d.kTot ? q : a), rows[0]), gS = rows.reduce((a, q) => (q.d.slack > a.d.slack ? q : a), rows[0]);
    // regular design wave Hmax at Tp on the far shape (deterministic cross-check)
    let TmaxWave = far.r.topTension;
    if (moving) {
      const xa = interp(rao, 'surge', inp.Tp) * Hmax / 2, za = interp(rao, 'heave', inp.Tp) * Hmax / 2;
      if (xa + za > 0) { const h = run(() => rodHarmonicResponse(far.r, { omega: wp, end: 'A', xAmp: xa, zAmp: za, phaseZ: inp.phase * DEG, rho: line.rhoW, zeta: inp.zeta }), 'Design-wave response'); TmaxWave += h.tensionAmp[0]; }
    }
    const sagMin = wave ? Math.min(...rows.map((q) => (q.r.sag ? q.r.sag.z : Infinity))) : null;
    if (wave && !(sagMin > 0 && Number.isFinite(sagMin))) throw new Error('No lazy wave forms in at least one storm position: the buoyancy section rests on the seabed or the wave is pulled flat. Increase buoyancy, lengthen the cable or reduce the offset.');
    const checks = [
      check('Extreme tension at the hang-off', gT.d.Tmax / 1e3, inp.Tallow, 'kN', `${STD_CABLE}; CIGRE TB 862 §4 (extreme load case)`, 'Static tension in the governing offset position plus the most probable maximum dynamic tension of the storm against the allowable tension. Exceeding it yields armour wires or overstrains the conductors in the storm peak.'),
      check('Extreme curvature', gK.d.kTot, 1 / inp.mbrAllow, '1/m', `${STD_CABLE} (MBR); CIGRE TB 862`, 'Static curvature plus the most probable maximum dynamic curvature at the worst point against the allowable. Tighter bending in the storm kinks the cable at the touchdown or the sag bend.'),
      check('Effective compression in the storm', gS.d.slack, 1, '-', 'DNV-RP-0360 (no compression in dynamic cables); project criterion', 'Largest ratio of extreme dynamic tension to local static tension. At 1.0 the cable goes slack in the largest cycles; beyond it it loops and snatches, which the linear analysis cannot follow.'),
    ];
    const gov = worst(checks);
    const sh = (c) => { const [xx, zz] = thin([c.r.x, c.r.z], 200); return { name: c.label, x: xx, y: zz, dash: c.label !== 'Mean' }; };
    const [sT, tSt, tMx, tMn] = thin([gT.r.s, Array.from(gT.r.Te, (v) => v / 1e3), Array.from(gT.d.tEnvMax, (v) => v / 1e3), Array.from(gT.d.tEnvMin, (v) => v / 1e3)], 300);
    const [sK, kSt, kEx] = thin([gK.r.s, Array.from(gK.r.curvature, Math.abs), Array.from(gK.d.kEnv)], 300);
    const periods = wg.map((w) => 2 * Math.PI / w), Smax = Math.max(...Sg, 1e-12), HTmax = Math.max(...far.d.HT, 1e-12);
    const angles = rows.map((q) => q.r.hangoffAngleDeg);
    return {
      kpis: [
        { key: 'Tstatic', label: 'Static top tension, far offset with storm current', value: round(far.r.topTension / 1e3, 6), unit: 'kN', explain: 'Hang-off tension of the static shape in the far position under the extreme current: the mean about which the wave-frequency tension varies.' },
        { key: 'tensionSig', label: 'Standard deviation of dynamic top tension', value: round(gT.d.sigTop / 1e3, 5), unit: 'kN', explain: 'Root of the area under the tension response spectrum in the governing position. Multiply by the peak factor for the storm extreme.' },
        { key: 'peakFactor', label: 'Peak factor √(2·ln N)', value: round(gT.d.c), unit: '-', explain: 'Most probable maximum over standard deviation for the number of response cycles in the storm (Rayleigh); typically 3.6–3.9 for three hours.' },
        { key: 'Tmax', label: 'Extreme top tension', value: round(gT.d.Tmax / 1e3, 6), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Static tension plus the most probable maximum dynamic tension in the governing offset position.' },
        { key: 'TmaxWave', label: 'Top tension in the regular design wave', value: round(TmaxWave / 1e3, 5), unit: 'kN', explain: 'Cross-check with a regular wave of the most probable maximum height at the peak period on the far shape. Large differences from the spectral value point to strong drag nonlinearity or a period-sensitive response.' },
        { key: 'Tmin', label: 'Lowest effective tension in the storm', value: round(Math.min(...rows.map((q) => q.d.tMin)) / 1e3, 4), unit: 'kN', explain: 'Smallest static minus extreme dynamic tension anywhere along the cable; a negative value is effective compression.' },
        { key: 'slackRatio', label: 'Extreme dynamic / static tension (worst point)', value: round(gS.d.slack), unit: '-', limit: 1, sense: 'max', explain: 'Above 1 the cable goes slack in the largest cycles.' },
        { key: 'mbrStorm', label: 'Minimum bend radius in the storm', value: round(1 / Math.max(gK.d.kTot, 1e-9)), unit: 'm', limit: inp.mbrAllow, sense: 'min', explain: 'Reciprocal of the largest static plus extreme dynamic curvature over all positions.' },
        { key: 'sMbr', label: 'Location of the tightest bend', value: round(gK.d.sK), unit: 'm', explain: 'Arc length from the hang-off of the governing curvature — hang-off region, sag or hog bend, or touchdown.' },
        { key: 'heaveSig', label: 'Standard deviation of hang-off heave', value: round(heaveSig, 5), unit: 'm', explain: 'Vertical motion of the hang-off in the storm from the heave RAO and the wave spectrum; four times this value is the significant double amplitude.' },
        { key: 'heaveMax', label: 'Most probable maximum hang-off heave', value: round(heaveMpm), unit: 'm', explain: 'Largest heave amplitude expected once in the storm duration.' },
        { key: 'surgeSig', label: 'Standard deviation of hang-off surge', value: round(surgeSig, 5), unit: 'm', explain: 'Wave-frequency horizontal motion of the hang-off, additional to the static offset.' },
        { key: 'Hmax', label: 'Most probable maximum wave height', value: round(Hmax), unit: 'm', explain: 'Hs·√(½·ln N) for the number of waves in the storm; used for the regular design-wave cross-check.' },
        { key: 'angleRange', label: 'Static hang-off angle range near → far', value: round(Math.max(...angles) - Math.min(...angles)), unit: 'deg', explain: 'Change of the departure angle between the storm offset positions: with the dynamic angle it sizes the bend stiffener.' },
        ...(wave ? [{ key: 'sagClearance', label: 'Lowest sag-bend clearance in the storm positions', value: round(sagMin), unit: 'm', explain: 'Static height of the sag bend above the seabed in its lowest offset position; deduct the dynamic motion and any buoyancy loss.' }] : []),
      ],
      checks,
      plots: [
        { type: 'line', title: 'Storm positions: near, mean and far', x: { label: 'Horizontal distance from the mean hang-off', unit: 'm' }, y: { label: 'Elevation above seabed', unit: 'm' }, series: rows.map(sh),
          explain: 'Static shapes under the extreme current at the storm offsets. The far shape stretches the cable (tension), the near shape compresses the wave and tightens the bends.' },
        { type: 'line', title: `Tension envelope (${gT.label.toLowerCase()} position)`, x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Effective tension', unit: 'kN' },
          series: [{ name: 'Static', x: sT, y: tSt }, { name: 'Storm maximum', x: sT, y: tMx, dash: true }, { name: 'Storm minimum', x: sT, y: tMn, dash: true }, { name: 'Allowable', x: [sT[0], sT[sT.length - 1]], y: [inp.Tallow, inp.Tallow], dash: true }],
          explain: 'Static tension with the most probable extreme envelope. The upper envelope must stay below the allowable and the lower one above zero along the whole length.' },
        { type: 'line', title: `Curvature envelope (${gK.label.toLowerCase()} position)`, x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Curvature', unit: '1/m' },
          series: [{ name: 'Static |κ|', x: sK, y: kSt }, { name: 'Storm extreme', x: sK, y: kEx, dash: true }, { name: 'Allowable (1/MBR)', x: [sK[0], sK[sK.length - 1]], y: [1 / inp.mbrAllow, 1 / inp.mbrAllow], dash: true }],
          explain: 'Static curvature and its storm extreme. Peaks at the hang-off, the sag and hog bends and the touchdown must stay below the dashed allowable.' },
        { type: 'line', title: 'Top-tension transfer function and wave spectrum (far position)', x: { label: 'Wave period', unit: 's' }, y: { label: 'Normalised to the maximum', unit: '-' },
          series: [{ name: 'Tension per metre of wave amplitude', x: periods, y: far.d.HT.map((v) => v / HTmax) }, { name: 'Wave spectrum', x: periods, y: Sg.map((v) => v / Smax) }],
          explain: 'Where the tension transfer function overlaps the wave spectrum the response is generated. A transfer function rising towards short periods means drag-locked behaviour: check shorter peak periods on the contour.' },
      ],
      tables: [
        { title: 'Storm cases', explain: 'One row per offset position: static top tension, standard deviation and extreme of the dynamic tension, slack ratio, static and storm bend radius and the zero-crossing period of the tension response.',
          columns: [{ key: 'case', label: 'Position' }, { key: 'Ts', label: 'Static top tension', unit: 'kN' }, { key: 'sig', label: 'σ tension', unit: 'kN' }, { key: 'Tmax', label: 'Extreme tension', unit: 'kN' }, { key: 'slack', label: 'Slack ratio', unit: '-' }, { key: 'mbrS', label: 'Static min radius', unit: 'm' }, { key: 'mbr', label: 'Storm min radius', unit: 'm' }, { key: 'angle', label: 'Hang-off angle', unit: 'deg' }, { key: 'tz', label: 'Response period', unit: 's' }],
          rows: rows.map((q) => ({ case: q.label, Ts: round(q.r.topTension / 1e3), sig: round(q.d.sigTop / 1e3), Tmax: round(q.d.Tmax / 1e3), slack: round(q.d.slack), mbrS: round(q.r.minBendRadius), mbr: round(1 / Math.max(q.d.kTot, 1e-9)), angle: round(q.r.hangoffAngleDeg), tz: round(q.d.tz) })) },
        { title: 'Transfer functions (far position)', explain: 'Hang-off motion amplitudes used for the drag linearisation and the top-tension amplitude per metre of wave amplitude at each solved period.',
          columns: [{ key: 'T', label: 'Period', unit: 's' }, { key: 'heave', label: 'Heave RAO', unit: 'm/m' }, { key: 'surge', label: 'Surge RAO', unit: 'm/m' }, { key: 'H', label: 'Tension RAO', unit: 'kN/m' }, { key: 'S', label: 'Wave spectrum', unit: 'm²·s/rad' }],
          rows: periods.map((T, j) => ({ T: round(T), heave: round(interp(rao, 'heave', T)), surge: round(interp(rao, 'surge', T)), H: round(far.d.HT[j] / 1e3), S: round(Sg[j]) })) },
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}), ${gT.label.toLowerCase()} position for tension and ${gK.label.toLowerCase()} position for curvature.`,
        `Sea state Hs ${inp.Hs} m, Tp ${inp.Tp} s, γ ${fmt(gam)}, Tz ${fmt(seaPar.Tz)} s; ${fmt(Nw)} waves in ${inp.duration} h, Hmax ${fmt(Hmax)} m.`,
        `Hang-off motion: heave σ ${fmt(heaveSig)} m (maximum ${fmt(heaveMpm)} m), surge σ ${fmt(surgeSig)} m, on top of a static offset of ±${inp.offset} m.`,
        `Extreme top tension ${fmt(gT.d.Tmax / 1e3)} kN from the spectrum against ${fmt(TmaxWave / 1e3)} kN in the regular design wave.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}: extreme top tension ${fmt(gT.d.Tmax / 1e3)} kN of ${inp.Tallow} kN, storm bend radius ${fmt(1 / Math.max(gK.d.kTot, 1e-9))} m against ${inp.mbrAllow} m at ${fmt(gK.d.sK)} m, slack ratio ${fmt(gS.d.slack)}.`,
          gov.uc > 1 ? (gov.name.startsWith('Effective') ? 'Remove the compression: more buoyancy or a longer buoyancy section to decouple the lower cable, or a longer upper catenary so that heave is taken up in geometry; then confirm in the time domain.' : gov.name.startsWith('Extreme curvature') ? `Open the governing bend at ${fmt(gK.d.sK)} m: a bend stiffener at the hang-off, a longer cable or re-tuned buoyancy for the near position; the mooring offset limit may need tightening.` : 'Lengthen the cable or move the anchor point to relieve the far position, reduce the offset through the mooring design, or use a cable with stronger armour.')
            : 'Confirm the governing case with an irregular-wave time-domain analysis over the 100-year contour (several peak periods and seeds), including the one-line-broken offset, and carry the extreme hang-off loads to the bend-stiffener and hang-off design.',
          `${STD_CABLE}; CIGRE TB 862 (2022) §4; DNV-RP-C205 (2021) §3.5; DNV-RP-F205 (2019)`),
        adv('economic', 'info', `The storm case uses ${fmt(gov.uc * 100)} % of the governing capacity; the far offset alone accounts for ${fmt(far.r.topTension / Math.max(gT.d.Tmax, 1) * 100)} % of the extreme tension and the wave-frequency part for the rest.`,
          gov.uc > 1 ? `A configuration change is the cheaper fix: about 10 m more cable or a few more buoyancy modules cost far less than a stronger cross-section over the whole ${fmt(inp.length)} m dynamic length, and a storm-damaged dynamic cable means a replacement campaign of the order of 10–20 vessel days plus lost production.`
            : `There is reserve of ${fmt((1 - gov.uc) * 100)} % in the storm case: before spending it on a shorter cable or less buoyancy, check that fatigue — usually the governing limit state of a dynamic cable — does not consume it. A storm failure means a replacement campaign of the order of 10–20 vessel days plus lost production.`,
          'Whole-life costing ISO 15686-5; risk-based contingency AACE RP 44R-08',
          `Replacing the dynamic section repeats its manufacture (about ${fmt(17 * inp.length / 1e3)} t of copper ≈ ${fmt(17 * inp.length / 1e3 * CO2_COPPER)} t CO₂e for a 3 × 630 mm² cable, screening ${CO2_COPPER} t CO₂e/t) and a vessel campaign: ${vesselCo2Text(15)}`),
      ],
      publish: { 'cable.stormTension': gT.d.Tmax, 'cable.stormMinBendRadius': 1 / Math.max(gK.d.kTot, 1e-9), 'cable.stormTensionSig': gT.d.sigTop },
      quantities: {},
    };
  },
  verify: [
    // Flat heave RAO of 0.6 m/m: the heave spectrum is 0.6²·S(ω), so σ_heave = 0.6·√m0 = 0.6·Hs/4 = 0.6 × 10/4 = 1.5 m
    // (m0 = Hs²/16 by definition of the significant wave height; the JONSWAP normalisation is within 0.3 % for γ ≤ 5).
    { name: 'Hang-off heave in an irregular sea with a flat RAO: σ = RAO·Hs/4',
      source: 'DNV-RP-C205 §3.5.3 (Hm0 = 4·√m0); linear response spectrum S_R = |RAO|²·S (hand calculation)',
      inputs: { rao: [{ T: 2, heave: 0.6, surge: 0 }, { T: 60, heave: 0.6, surge: 0 }], Hs: 10, Tp: 14, nFreq: 4, nElem: 40 },
      expect: [{ kpi: 'heaveSig', value: 1.5, tol: 0.01 }] },
    // Limiting case: no waves, no current, no offset → the extreme tension is the static catenary value.
    // Irvine (1981) §1.2: h = 138 m, L = 300 m, X = 200 m → a ≈ 21.5 m; w = 415.7 N/m → Ttop = w·(a + h) = 66.3 kN.
    { name: 'Calm limit: extreme tension equals the static catenary tension w·(a + h)',
      source: 'Limiting-case physics; Irvine, Cable Structures (1981) §1.2 (hand solution a ≈ 21.5 m)',
      inputs: { config: 'catenary', length: 300, span: 200, Us: 0, Ub: 0, Hs: 0, offset: 0, EI: 200, muAx: 0, nElem: 160, mbrAllow: 0.5, Tallow: 5000 },
      expect: [{ kpi: 'tensionSig', value: 0, abs: 1e-9 }, { kpi: 'Tmax', value: 66.3, tol: 0.015 }] },
  ],
};
