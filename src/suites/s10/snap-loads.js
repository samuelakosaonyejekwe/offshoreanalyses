// Transient snap loading of a mooring line: slack → re-tension events in a lumped-mass simulation with large
// fairlead motion, compared with an energy-based snap estimate.
import { detectSnapLoads } from '../../engine/lines/lumpedmass.js';
import { dynInputs, runLine, LIMIT_LM } from './_lumped.js';
import { round, fmt, adv, check, sev, thin, steelCo2Text, LIMIT_PROPS, CONFIRM } from './_common.js';

export default {
  id: 's10-snap-loads',
  suite: 10,
  title: 'Transient snap loading',
  summary: 'Detects whether a line goes slack and snaps taut under large fairlead motion and quantifies the snap peak against the static tension and the break load.',
  covers: [12],
  fidelity: 'numerical',
  standards: ['DNV-OS-E301 (2018) Ch.2 Sec.2 [2.7]', 'DNV-RP-N103 (2017) §4.7 (snap forces)', 'API RP 2SK (3rd ed., 2005) §7.3'],
  models: ['Lumped-mass line with tension-only segments', 'Slack / re-tension event detection', 'Morison drag and added mass', 'Seabed contact with regularised Coulomb friction', 'Work–energy snap estimate F = v·√(K·M)'],
  limits: `${LIMIT_LM} Segments carry no compression, so slack is represented, but the re-tensioning peak is sensitive to segmentation, axial damping and time step: treat the peak as an estimate and confirm it with the refinement study. ${LIMIT_PROPS}`,
  searchTerms: 'mooring line snap load slack re-tension dynamic',
  inputs: [
    ...dynInputs({ X: 880, ax: 4, az: 6, period: 7, nCycles: 4 }),
    { key: 'slackFraction', label: 'Slack threshold (fraction of static tension)', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.5, group: 'Acceptance criteria', role: 'criteria',
      help: 'The line is called slack when its fairlead tension falls below this share of the static tension (10 % is customary).' },
    { key: 'sfSnap', label: 'Required safety factor on the snap peak', unit: '-', type: 'number', default: 1.67, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria',
      help: `Applied to the largest tension in the simulation. ${CONFIRM}` },
  ],
  mesh: { key: 'nSeg', label: 'Line segments', levels: [10, 20, 40], qoi: 'Tmax' },
  meshTime: { key: 'cfl', label: 'Time step / critical step', levels: [0.8, 0.4, 0.2], qoi: 'Tmax', isSize: true },
  compute(inp, ctx) {
    const { p, cat, model, res, period } = runLine(inp, ctx);
    const st = res.stats;
    const Ts = res.static.fairleadTension;
    const snap = detectSnapLoads({ t: res.t, T: res.fairleadTension, Tref: Ts, slackFraction: inp.slackFraction, spikeFactor: 1.0 });
    const events = snap.events.filter((e) => Number.isFinite(e.peak));
    const worstEv = events.reduce((a, e) => (!a || e.peak > a.peak ? e : a), null);
    // Energy estimate (DNV-RP-N103 §4.7.3 form): the fairlead velocity along the line is arrested by the axial spring EA/L
    // acting on the suspended line mass moving with a linear velocity profile (effective mass = one third).
    const om = 2 * Math.PI / period;
    const vAx = om * Math.hypot(inp.ax * Math.cos(cat.angle), inp.az * Math.sin(cat.angle));
    const K = p.EA / inp.L, Meff = p.mass * cat.suspendedLength / 3;
    const Fsnap = vAx * Math.sqrt(K * Meff);
    const uc = st.fairleadMax * inp.sfSnap / p.MBL;
    const slackUc = Math.max(0, inp.slackFraction * Ts - st.fairleadMin) / Math.max(inp.slackFraction * Ts, 1);
    const [t, Tf, Tq] = thin([res.t, res.fairleadTension, res.quasiStaticTension], 800);
    const dynInc = st.fairleadMax - Ts, qsInc = (st.quasiStaticMax ?? Ts) - Ts;
    const daf = qsInc > 1e-6 * Ts ? dynInc / qsInc : 1;
    const sMid = Array.from({ length: model.N }, (_, i) => 0.5 * (model.s[i] + model.s[i + 1]));
    const mass = p.mass * inp.L / 1e3;
    const checks = [
      check('Snap / maximum tension against break load', st.fairleadMax * inp.sfSnap / 1e3, p.MBL / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3; DNV-OS-E301 Ch.2 Sec.2 [4]',
        'Guards against the re-tensioning shock breaking the line or a connector; snap peaks are short and can exceed the wave-frequency maximum several times.'),
      { name: 'No slack line', demand: round(Math.max(0, inp.slackFraction * Ts - st.fairleadMin) / 1e3, 4), capacity: round(inp.slackFraction * Ts / 1e3, 4), unit: 'kN', uc: slackUc,
        ref: 'DNV-RP-N103 (2017) §4.7.2 slack criterion; DNV-OS-E301 Ch.2 Sec.2',
        explain: 'Guards against the slack condition itself: a utilisation above zero means tension fell below the threshold, and at 1.0 it reached zero. Slack chain can knot, hook the fairlead or load connectors out of plane.' },
    ];
    const g = checks[0].uc >= checks[1].uc ? checks[0] : checks[1];
    return {
      kpis: [
        { key: 'Tstatic', label: 'Static fairlead tension', value: round(Ts / 1e3, 5), unit: 'kN', explain: 'Tension at rest, the reference for the slack threshold and the snap ratio.' },
        { key: 'Tmax', label: 'Maximum tension (snap peak)', value: round(st.fairleadMax / 1e3, 5), unit: 'kN', limit: round(p.MBL / inp.sfSnap / 1e3, 5), sense: 'max', explain: 'Largest fairlead tension after start-up; when slack events occur this is the re-tensioning (snap) peak.' },
        { key: 'Tmin', label: 'Minimum tension', value: round(st.fairleadMin / 1e3, 5), unit: 'kN', limit: round(inp.slackFraction * Ts / 1e3, 5), sense: 'min', explain: 'Smallest fairlead tension; below the slack threshold the line is considered slack.' },
        { key: 'snapRatio', label: 'Peak tension / static tension', value: round(st.fairleadMax / Ts, 4), unit: '-', explain: 'Transient amplification of the static tension. Values of 2–5 are typical once a chain line goes slack in a short-period motion.' },
        { key: 'nSlack', label: 'Slack events', value: snap.events.length, unit: '-', explain: 'Number of times the tension dropped below the slack threshold during the simulation (normally once per motion cycle when it happens at all).' },
        { key: 'slackDuration', label: 'Longest slack duration', value: round(events.reduce((a, e) => Math.max(a, e.slackDuration || 0), 0), 4), unit: 's', explain: 'Longest time the line stayed below the slack threshold; longer slack means more free fall and a harder re-tensioning.' },
        { key: 'daf', label: 'Dynamic amplification over quasi-static', value: round(daf, 4), unit: '-', explain: 'Dynamic tension increase divided by the quasi-static increase for the same motion; very large values confirm the response is transient, not catenary-like.' },
        { key: 'snapEstimate', label: 'Energy-based snap force estimate', value: round(Fsnap / 1e3, 5), unit: 'kN', explain: 'Screening estimate v·√(K·M) of the tension needed to arrest the line at the fairlead velocity; compare its order of magnitude with the simulated dynamic increment.' },
        { key: 'util', label: 'Tmax / MBL', value: round(st.fairleadMax / p.MBL, 4), unit: '-', explain: 'Share of the break load used at the snap peak.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Fairlead tension with slack threshold', x: { label: 'Time', unit: 's' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Dynamic', x: t.map((v) => round(v)), y: Tf.map((v) => round(v / 1e3)) }, { name: 'Quasi-static', x: t.map((v) => round(v)), y: Tq.map((v) => round(v / 1e3)), dash: true },
            { name: 'Slack threshold', x: [0, res.t[res.t.length - 1]], y: [snap.threshold / 1e3, snap.threshold / 1e3], dash: true }, { name: 'Allowable', x: [0, res.t[res.t.length - 1]], y: [p.MBL / inp.sfSnap / 1e3, p.MBL / inp.sfSnap / 1e3], dash: true }],
          explain: 'Tension history against the slack threshold and the allowable. A snap shows as a trough touching the threshold followed by a sharp spike well above the quasi-static curve.' },
        { type: 'line', title: 'Tension envelope along the line', x: { label: 'Arc length from anchor', unit: 'm' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Maximum', x: sMid.map((v) => round(v)), y: Array.from(res.envelopes.Tmax, (v) => round(v / 1e3)) }, { name: 'Minimum', x: sMid.map((v) => round(v)), y: Array.from(res.envelopes.Tmin, (v) => round(v / 1e3)) }],
          explain: 'Where along the line the slack (minimum near zero) and the snap peak occur. Slack usually starts near the touchdown; the peak travels up to the fairlead.' },
        { type: 'line', title: 'Energy diagnostics', x: { label: 'Time', unit: 's' }, y: { label: 'Energy', unit: 'kJ' },
          series: ['kinetic', 'strain'].map((k) => { const [tt, e] = thin([res.t, res.energy[k]], 400); return { name: `${k[0].toUpperCase()}${k.slice(1)} energy`, x: tt.map((v) => round(v)), y: e.map((v) => round(v / 1e3)) }; }),
          explain: 'Kinetic energy built up while the line is slack is converted into strain energy at re-tensioning: a strain-energy spike right after a kinetic-energy peak is the signature of a snap.' },
      ],
      tables: [
        { title: 'Slack and snap events', columns: [{ key: 'tS', label: 'Slack starts', unit: 's' }, { key: 'dur', label: 'Slack duration', unit: 's' }, { key: 'tP', label: 'Peak at', unit: 's' }, { key: 'peak', label: 'Re-tension peak', unit: 'kN' }, { key: 'ratio', label: 'Peak / static' }],
          rows: events.slice(0, 60).map((e) => ({ tS: round(e.tSlack), dur: round(e.slackDuration), tP: round(e.tPeak), peak: round(e.peak / 1e3), ratio: round(e.ratio, 3) })),
          explain: 'One row per slack event with the peak that followed it. An empty table means the line stayed above the slack threshold throughout.' },
      ],
      notes: [snap.events.length ? `${snap.events.length} slack event(s); largest re-tension peak ${fmt((worstEv ? worstEv.peak : st.fairleadMax) / 1e3)} kN = ${fmt((worstEv ? worstEv.peak : st.fairleadMax) / Ts)} × static.` : `No slack: minimum tension ${fmt(st.fairleadMin / 1e3)} kN stays above the ${fmt(snap.threshold / 1e3)} kN threshold.`,
        `Energy estimate ${fmt(Fsnap / 1e3)} kN from axial fairlead velocity ${fmt(vAx)} m/s, line stiffness ${fmt(K / 1e3)} kN/m and effective mass ${fmt(Meff / 1e3)} t.`,
        'Seabed friction is regularised Coulomb; snap peaks depend on segmentation and damping — check the refinement study.'],
      advice: [
        adv('engineering', sev(g.uc), snap.events.length ? `The line goes slack ${snap.events.length} time(s) and snaps to ${fmt(st.fairleadMax / 1e3)} kN (${fmt(st.fairleadMax / Ts)} × static, ${fmt(100 * st.fairleadMax / p.MBL)} % MBL).` : `No slack at this motion: minimum tension is ${fmt(100 * st.fairleadMin / Ts)} % of static, maximum ${fmt(100 * st.fairleadMax / p.MBL)} % MBL.`,
          snap.events.length ? `Remove the slack rather than design for the snap: raise pretension so the static tension exceeds the dynamic amplitude (≈ ${fmt((Ts - st.fairleadMin) / 1e3)} kN here), add a clump weight near touchdown, or insert a compliant polyester length to lower the axial stiffness (snap force ∝ √(EA/L)).`
            : `Margin to slack is ${fmt((st.fairleadMin - snap.threshold) / 1e3)} kN; repeat for the leeward line at its lowest mean tension, where slack appears first.`,
          'DNV-RP-N103 (2017) §4.7 (slack criterion and snap force); DNV-OS-E301 (2018) Ch.2 Sec.2 [2.7].'),
        adv('economic', snap.events.length ? 'warn' : 'info', snap.events.length ? `Repeated snap loading consumes fatigue life roughly with the cube of the tension range: a range of ${fmt((st.fairleadMax - st.fairleadMin) / 1e3)} kN is ${fmt(((st.fairleadMax - st.fairleadMin) / Math.max(2 * (Ts - st.fairleadMin), 1)) ** 3)} × as damaging per cycle as the same motion without the snap.` : `No snap-related fatigue penalty; the line holds ${fmt(mass)} t and has ${fmt(100 * (1 - uc))} % strength reserve.`,
          snap.events.length ? 'Budget for shorter inspection intervals and earlier top-chain replacement on lines that snap, or spend once on pretension or compliance to eliminate it.' : 'Keep the present pretension; lowering it to save tensioning time would bring the slack threshold closer.',
          'Risk-based inspection planning per DNV-RP-G101 / API RP 2I (in-service inspection of mooring hardware); whole-life costing per ISO 15686-5.',
          `${steelCo2Text(mass)} Early replacement of a snapped or worn top section repeats part of that footprint and adds vessel emissions.`),
      ],
      quantities: p.family === 'fibre' ? { polymer_t: round(mass) } : { steel_t: round(mass) },
    };
  },
  verify: [
    { // Same taut vertical line as the line-dynamics benchmark, at rest: top tension 1645.9 kN, no slack events.
      name: 'No motion: peak equals the static tension and there is no slack',
      source: 'Hand calculation (Hooke\'s law plus self-weight); limiting case of zero excitation',
      inputs: { lineType: 'custom', wCustom: 1500, EAcustom: 800, MBLcustom: 9000, massCustom: 176, L: 279.5, X: 1, h: 280, ax: 0, az: 0, nCycles: 2, period: 5 },
      expect: [{ kpi: 'Tmax', value: 1645.9, tol: 0.01 }, { kpi: 'nSlack', value: 0, abs: 0 }, { kpi: 'snapRatio', value: 1, tol: 0.005 }] },
    { // Energy estimate for the taut vertical line with 0.5 m heave at 5 s: v = 2π/5·0.5 = 0.6283 m/s (line is vertical),
      // K = EA/L = 800e6/279.5 = 2.8623e6 N/m, M = 176·279.5/3 = 16 397 kg → F = v·sqrt(K·M) = 0.6283·216 640 = 136.1 kN.
      name: 'Energy-based snap estimate F = v·√(K·M)',
      source: 'DNV-RP-N103 (2017) §4.7.3 formula evaluated by hand',
      inputs: { lineType: 'custom', wCustom: 1500, EAcustom: 800, MBLcustom: 9000, massCustom: 176, L: 279.5, X: 1, h: 280, ax: 0, az: 0.5, nCycles: 2, period: 5 },
      expect: [{ kpi: 'snapEstimate', value: 136.1, tol: 0.01 }] },
  ],
};
