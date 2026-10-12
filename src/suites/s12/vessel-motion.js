// Vessel-motion-induced dynamic response in the frequency domain: hang-off surge and heave from the host RAOs
// applied to the static shape, with linearised drag; tension and bending amplification over the quasi-static response.
import { rodHarmonicResponse, rodModes } from '../../engine/lines/rod.js';
import { riserOffsetCases } from '../../engine/lines/riser.js';
import { lineInputs, configInputs, solveConfig, adv, check, worst, sev, fmt, round, run, need, thin, scaled, DEG, YEAR_S, LIMIT_2D, STD_UMB } from './_common.js';

const interp = (xs, ys, x) => { if (x <= xs[0]) return ys[0]; for (let i = 1; i < xs.length; i++) if (x <= xs[i]) return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[ys.length - 1]; };

export default {
  id: 's12-vessel-motion-response', suite: 12, title: 'Vessel-motion-induced dynamic response',
  summary: 'Dynamic tension and curvature along the umbilical caused by host surge and heave at the hang-off, across wave periods, with amplification over the quasi-static response.',
  covers: [9], fidelity: 'numerical',
  standards: [`${STD_UMB} §7.3 (dynamic analysis, host motions)`, 'DNV-RP-C205 (2021) §6.2 (drag linearisation)', 'API RP 2RD (1998) §5.5 (quasi-static and frequency-domain methods)'],
  models: ['Vessel-motion response by RAO at the hang-off', 'Dynamic cable finite elements linearised about the static shape (frequency domain)', 'Added-mass and harmonically linearised drag formulations', 'Modal analysis of the tensioned rod', 'Quasi-static offset method'],
  limits: `${LIMIT_2D} Response is linear about the static shape with frozen seabed contact; direct wave loading on the line, steady current in the drag linearisation and stick–slip bending hysteresis are not included. Sway and roll (out of plane) need the 3-D time-domain module.`,
  publishes: ['umbilical.tensionRange', 'umbilical.curvatureRange', 'umbilical.dafVesselMotion'],
  searchTerms: 'umbilical vessel motion RAO dynamic tension curvature frequency domain',
  inputs: [
    ...lineInputs(), ...configInputs(),
    { key: 'H', label: 'Wave height (regular)', unit: 'm', type: 'number', default: 6, min: 0.1, max: 30, group: 'Wave and host motion', role: 'bc', bind: 'site.Hs1', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'Td', label: 'Design wave period', unit: 's', type: 'number', default: 10, min: 3, max: 300, group: 'Wave and host motion', role: 'bc', bind: 'site.Tp1' },
    { key: 'rao', label: 'Hang-off motion RAOs', type: 'table', group: 'Wave and host motion', role: 'bc',
      columns: [{ key: 'T', label: 'Period', unit: 's', type: 'number' }, { key: 'surge', label: 'Surge RAO', unit: 'm/m', type: 'number' }, { key: 'heave', label: 'Heave RAO', unit: 'm/m', type: 'number' }],
      default: [{ T: 5, surge: 0.05, heave: 0.03 }, { T: 7, surge: 0.2, heave: 0.15 }, { T: 9, surge: 0.45, heave: 0.5 }, { T: 11, surge: 0.65, heave: 0.85 }, { T: 13, surge: 0.8, heave: 1.0 }, { T: 16, surge: 0.9, heave: 1.05 }, { T: 20, surge: 0.95, heave: 1.0 }],
      help: 'Motion amplitude at the hang-off per unit wave amplitude (include pitch × lever arm in heave and surge), from the host motion analysis. Ascending periods.' },
    { key: 'phase', label: 'Heave phase lead over surge', unit: 'deg', type: 'number', default: 90, min: -180, max: 180, group: 'Wave and host motion', role: 'bc' },
    { key: 'zeta', label: 'Structural damping ratio', unit: '-', type: 'number', default: 0.02, min: 0, max: 0.3, group: 'Wave and host motion', role: 'ic', help: 'Stiffness-proportional damping at the wave frequency; bending hysteresis of an umbilical gives a few per cent.' },
    { key: 'Tallow', label: 'Allowable tension', unit: 'kN', type: 'number', default: 350, min: 5, max: 5000, group: 'Acceptance criteria', role: 'criteria', bind: 'umbilical.allowTension' },
    { key: 'mbrAllow', label: 'Allowable minimum bend radius (dynamic)', unit: 'm', type: 'number', default: 11, min: 0.5, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'umbilical.mbrDynamic' },
    { key: 'nElem', label: 'Rod elements', type: 'int', default: 120, min: 30, max: 1200, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Rod elements', levels: [60, 120, 240], qoi: 'tensionAmp' },
  compute(inp) {
    const rows = (inp.rao || []).filter((r) => r.T > 0).sort((a, b) => a.T - b.T);
    need(rows.length >= 2, 'Give at least two RAO rows (period, surge, heave) so that the design period can be interpolated.');
    for (const r of rows) need(r.surge >= 0 && r.heave >= 0, 'RAO amplitudes cannot be negative (use the phase input for the sign).');
    const { res, line } = solveConfig(inp, { nElem: inp.nElem });
    const a = inp.H / 2, n = res.s.length - 1, ph = inp.phase * DEG;
    const Ts = rows.map((r) => r.T), raoS = (T) => interp(Ts, rows.map((r) => r.surge), T), raoH = (T) => interp(Ts, rows.map((r) => r.heave), T);
    // the rod x-axis points from the hang-off to the termination: positive surge = towards the termination
    const harm = (T) => {
      const r = run(() => rodHarmonicResponse(res, { omega: 2 * Math.PI / T, end: 'A', xAmp: raoS(T) * a, zAmp: raoH(T) * a, phaseZ: ph, rho: line.rhoW, zeta: inp.zeta }), `Harmonic response at ${T} s`);
      let kMax = 0, sK = 0, slack = 0, sSlack = 0, kTot = 0, sTot = 0;
      for (let i = 0; i <= n; i++) {
        if (r.curvatureAmp[i] > kMax) { kMax = r.curvatureAmp[i]; sK = res.s[i]; }
        const tot = Math.abs(res.curvature[i]) + r.curvatureAmp[i]; if (tot > kTot) { kTot = tot; sTot = res.s[i]; }
        const ratio = r.tensionAmp[i] / Math.max(res.Te[i], 1); if (ratio > slack) { slack = ratio; sSlack = res.s[i]; }
      }
      return { r, T, tTop: r.tensionAmp[0], kMax, sK, slack, sSlack, kTot, sTot };
    };
    // quasi-static sensitivities of the top tension and curvature to unit hang-off heave and surge (one-sided re-solves)
    const dz = Math.max(0.2, raoH(inp.Td) * a), dx = Math.max(0.2, raoS(inp.Td) * a);
    const up = run(() => riserOffsetCases(res, [0], { heave: dz }), 'Quasi-static heave case')[0].result;
    const sg = run(() => riserOffsetCases(res, [-dx]), 'Quasi-static surge case')[0].result;   // negative offset = towards the termination = +x
    const kz = (up.topTension - res.topTension) / dz, kx = (sg.topTension - res.topTension) / dx;
    const qsT = (T) => { const x = raoS(T) * a, z = raoH(T) * a; return Math.sqrt(Math.max(0, (kx * x) ** 2 + (kz * z) ** 2 + 2 * kx * x * kz * z * Math.cos(ph))); };
    let qsK = 0;
    { const x = raoS(inp.Td) * a / dx, z = raoH(inp.Td) * a / dz;
      for (let i = 0; i <= n; i++) { const cx = (sg.curvature[i] - res.curvature[i]) * x, cz = (up.curvature[i] - res.curvature[i]) * z; qsK = Math.max(qsK, Math.sqrt(Math.max(0, cx * cx + cz * cz + 2 * cx * cz * Math.cos(ph)))); } }
    const d = harm(inp.Td);
    need(d.r.converged !== false, 'The drag linearisation did not converge at the design period.');
    const sweep = rows.map((r) => (Math.abs(r.T - inp.Td) < 1e-9 ? d : harm(r.T)));
    const modes = run(() => rodModes(res, { nModes: 4 }), 'Modal analysis');
    const Tn = modes.modes.map((m) => m.periodS);
    const qs = qsT(inp.Td), dafT = qs > 1 ? d.tTop / qs : 1, dafK = qsK > 1e-9 ? d.kMax / qsK : 1;
    const Tmax = res.topTension + d.tTop, Tmin = res.topTension - d.tTop;
    const checks = [
      check('Maximum dynamic tension at the hang-off', Tmax / 1e3, inp.Tallow, 'kN', `${STD_UMB} §7.3`, 'Static plus dynamic top tension at the design wave against the allowable tension of the cross-section.'),
      check('Effective compression (dynamic amplitude / static tension)', d.slack, 1, '-', `${STD_UMB} §7.3 (no compression in dynamic sections); project criterion`, 'Largest ratio of dynamic tension amplitude to local static tension along the line. At 1.0 the tension reaches zero each cycle; beyond it the line goes into compression, buckles locally and snatches.'),
      check('Dynamic bend radius', d.kTot, 1 / inp.mbrAllow, '1/m', `${STD_UMB} §7.3 (dynamic MBR)`, 'Static curvature plus dynamic curvature amplitude at the worst point against the allowable dynamic curvature. Exceeding it overstresses tubes each wave cycle.'),
    ];
    const gov = worst(checks);
    const [s, tA, kA, kS, nA] = thin([res.s, scaled(d.r.tensionAmp, 1e-3), d.r.curvatureAmp, Array.from(res.curvature, Math.abs), d.r.normalAmp], 300);
    const cyclesYr = YEAR_S / inp.Td;
    return {
      kpis: [
        { key: 'Tstatic', label: 'Static top tension', value: round(res.topTension / 1e3), unit: 'kN', explain: 'Hang-off tension of the static shape about which the motion is applied.' },
        { key: 'tensionAmp', label: 'Dynamic top-tension amplitude', value: round(d.tTop / 1e3, 6), unit: 'kN', explain: 'Amplitude of the hang-off tension variation at the design wave. Twice this value is the tension range per cycle used in fatigue.' },
        { key: 'Tmax', label: 'Maximum top tension', value: round(Tmax / 1e3), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Static tension plus the dynamic amplitude.' },
        { key: 'Tmin', label: 'Minimum top tension', value: round(Tmin / 1e3), unit: 'kN', explain: 'Static tension minus the dynamic amplitude; it should stay comfortably positive.' },
        { key: 'daf', label: 'Dynamic amplification Tmax / Tstatic', value: round(Tmax / res.topTension, 4), unit: '-', explain: 'Maximum over static top tension at the design wave: the factor a static design tension must be multiplied by.' },
        { key: 'dafT', label: 'Tension amplification over quasi-static', value: round(dafT, 3), unit: '-', explain: 'Dynamic tension amplitude divided by the amplitude a slow (quasi-static) motion of the same size would cause. Above 1 the line cannot follow the hang-off: drag and inertia add tension. Compliant wave configurations show very large ratios because their quasi-static tension change is almost nil — the quasi-static method is then not usable for them.' },
        { key: 'curvAmp', label: 'Maximum dynamic curvature amplitude', value: round(d.kMax, 3), unit: '1/m', explain: `Largest curvature variation along the line (at ${fmt(d.sK, 4)} m from the hang-off). Twice this value is the curvature range per cycle for bending fatigue.` },
        { key: 'dafK', label: 'Bending amplification over quasi-static', value: round(dafK, 3), unit: '-', explain: 'Dynamic curvature amplitude over the quasi-static curvature change for the same hang-off motion. Large values show a travelling bending wave rather than a slow change of shape.' },
        { key: 'minRadiusDyn', label: 'Minimum dynamic bend radius', value: round(1 / d.kTot), unit: 'm', limit: inp.mbrAllow, sense: 'min', explain: `Inverse of static plus dynamic curvature at the worst point (${fmt(d.sTot, 4)} m from the hang-off).` },
        { key: 'slackRatio', label: 'Largest dynamic / static tension ratio', value: round(d.slack, 3), unit: '-', limit: 1, sense: 'max', explain: `Tension amplitude over local static tension, worst at ${fmt(d.sSlack, 4)} m from the hang-off (usually the sag bend or touchdown). At 1.0 the line goes slack.` },
        { key: 'T1', label: 'First natural period (in plane)', value: round(Tn[0]), unit: 's', explain: 'Longest in-plane natural period of the suspended line. Wave periods near a natural period amplify motion; long catenaries have many modes inside the wave band.' },
        { key: 'motionMax', label: 'Largest lateral motion amplitude along the line', value: round(Math.max(...d.r.normalAmp), 3), unit: 'm', explain: 'Largest displacement amplitude normal to the line axis; it shows how far the hang-off motion propagates down the line.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Top-tension amplitude vs wave period', x: { label: 'Wave period', unit: 's' }, y: { label: 'Tension amplitude', unit: 'kN' },
          series: [{ name: 'Dynamic', x: sweep.map((q) => q.T), y: sweep.map((q) => q.tTop / 1e3) }, { name: 'Quasi-static', x: sweep.map((q) => q.T), y: sweep.map((q) => qsT(q.T) / 1e3), dash: true }, { name: 'Design period', x: [inp.Td], y: [d.tTop / 1e3] }],
          explain: 'Dynamic and quasi-static tension amplitude for the same wave height at each RAO period. The gap between the curves is the dynamic amplification; a peak marks the period at which the host RAO and the line dynamics combine worst.' },
        { type: 'line', title: 'Maximum curvature amplitude vs wave period', x: { label: 'Wave period', unit: 's' }, y: { label: 'Curvature amplitude', unit: '1/m' },
          series: [{ name: 'Dynamic', x: sweep.map((q) => q.T), y: sweep.map((q) => q.kMax) }],
          explain: 'Largest dynamic curvature anywhere along the line at each period. Short periods excite bending waves near the hang-off and touchdown even when the motion RAO is small.' },
        { type: 'line', title: 'Dynamic tension amplitude along the line', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Tension amplitude', unit: 'kN' }, series: [{ name: `T = ${inp.Td} s`, x: s, y: tA }],
          explain: 'How the dynamic tension decays (or not) down the line. A lazy wave should cut it sharply beyond the buoyant section, protecting the touchdown and termination.' },
        { type: 'line', title: 'Static and dynamic curvature along the line', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Curvature', unit: '1/m' },
          series: [{ name: 'Static', x: s, y: kS }, { name: 'Dynamic amplitude', x: s, y: kA }, { name: 'Allowable (1/MBR)', x: [s[0], s[s.length - 1]], y: [1 / inp.mbrAllow, 1 / inp.mbrAllow], dash: true }],
          explain: 'Static curvature with the dynamic amplitude superimposed. The sum must stay below the dashed line; the location of the largest dynamic amplitude is the bending-fatigue hot spot.' },
        { type: 'line', title: 'Lateral motion amplitude along the line', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Normal displacement amplitude', unit: 'm' }, series: [{ name: `T = ${inp.Td} s`, x: s, y: nA }],
          explain: 'Displacement amplitude normal to the axis. It shows where the motion is absorbed by drag and where the buoyant section or the seabed pins the line.' },
      ],
      tables: [{ title: 'Response by wave period', explain: 'One row per RAO period for the same wave height: hang-off motion amplitudes, dynamic and quasi-static top-tension amplitude, their ratio, the largest curvature amplitude and the slack ratio.',
        columns: [{ key: 'T', label: 'Period', unit: 's' }, { key: 'x', label: 'Surge amp.', unit: 'm' }, { key: 'z', label: 'Heave amp.', unit: 'm' }, { key: 'dT', label: 'Tension amp.', unit: 'kN' }, { key: 'qs', label: 'Quasi-static', unit: 'kN' }, { key: 'daf', label: 'Amplification', unit: '-' }, { key: 'k', label: 'Curvature amp.', unit: '1/m' }, { key: 'slack', label: 'Slack ratio', unit: '-' }],
        rows: sweep.map((q) => ({ T: q.T, x: round(raoS(q.T) * a), z: round(raoH(q.T) * a), dT: round(q.tTop / 1e3), qs: round(qsT(q.T) / 1e3), daf: qsT(q.T) > 1 ? round(q.tTop / qsT(q.T), 3) : null, k: round(q.kMax, 3), slack: round(q.slack, 3) })) }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `At ${inp.Td} s the hang-off moves ±${fmt(raoS(inp.Td) * a)} m in surge and ±${fmt(raoH(inp.Td) * a)} m in heave; top tension varies by ±${fmt(d.tTop / 1e3)} kN, ${fmt(dafT)} × the quasi-static value.`,
        `In-plane natural periods: ${Tn.map((t) => fmt(t)).join(', ')} s.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Tension amplification ${fmt(dafT)}, bending amplification ${fmt(dafK)}; curvature range ${fmt(2 * d.kMax)} 1/m at ${fmt(d.sK, 4)} m from the hang-off.`,
          gov.uc > 1 ? 'Reduce the motion delivered to the line (hang-off nearer the motion centre, different heading) or make the configuration more compliant (more buoyancy, longer line); then confirm in the 3-D time-domain module.' : 'Feed the tension range and curvature range per sea state into the fatigue module, and size the bend stiffener for the hang-off angle range; verify the extreme with an irregular-wave time-domain run.',
          `${STD_UMB} §7.3; DNV-RP-C205 (2021) §6.2`),
        adv('economic', 'info', `Each cycle at the design wave carries a tension range of ${fmt(2 * d.tTop / 1e3)} kN and a curvature range of ${fmt(2 * d.kMax)} 1/m; a continuous year at this period would give ${fmt(cyclesYr / 1e6)} million cycles.`,
          'Because fatigue life varies with the 3rd–5th power of these ranges, a hang-off relocation or configuration change that cuts them by 20 % typically doubles the life of the dynamic section and defers its replacement campaign.',
          'Whole-life costing ISO 15686-5; value-of-design-change assessment at FEED (AACE RP 18R-97 Class 3)',
          'Longer fatigue life of the dynamic section avoids re-manufacturing its steel tubes and armour and a replacement lay campaign with its vessel fuel; reducing motion at source is the lowest-material route to that life.'),
      ],
      publish: { 'umbilical.tensionRange': 2 * d.tTop, 'umbilical.curvatureRange': 2 * d.kMax, 'umbilical.dafVesselMotion': dafT },
    };
  },
  verify: [
    { name: 'Very slow hang-off heave → dynamic response equals the quasi-static response (amplification → 1)',
      source: 'Limiting case: as ω → 0 inertia and drag vanish and (K − ω²M + iωC)U → K·U, the tangent-stiffness (quasi-static) solution',
      inputs: { config: 'catenary', length: 600, span: 330, Td: 250, H: 1, phase: 0, rao: [{ T: 200, surge: 0, heave: 1 }, { T: 300, surge: 0, heave: 1 }] },
      expect: [{ kpi: 'dafT', value: 1, tol: 0.06 }] },
  ],
};
