// Dynamic response of a dynamic power cable to wave-frequency motion of its hang-off: frequency-domain rod solution
// about the static shape with harmonically linearised Morison drag, natural periods, and a period sweep.
import { rodHarmonicResponse, rodModes } from '../../engine/lines/rod.js';
import { quasiStaticHeave } from '../../engine/lines/riser.js';
import { lineInputs, configInputs, currentInputs, solveConfig, rodCurrent, adv, check, worst, sev, fmt, round, run, need, thin, scaled, linspace, DEG, LIMIT_2D, STD_CABLE, CO2_COPPER } from './_common.js';

export default {
  id: 's13-dynamic-response', suite: 13, title: 'Dynamic cable response',
  summary: 'Dynamic tension, curvature and motion along a dynamic cable for a design wave acting through the platform motion at the hang-off.',
  covers: [2], fidelity: 'numerical',
  standards: [`${STD_CABLE} (dynamic analysis of cables at floating units)`, 'CIGRE TB 862 (2022) §4 (dynamic cable design loads)', 'DNV-RP-C205 (2021) §6.2 (Morison load, drag linearisation)', 'DNV-RP-F205 (2019) (coupled analysis — hang-off motion as input)'],
  models: ['Elastic rod / geometrically nonlinear finite elements (tangent stiffness about the static shape)', 'Consistent mass with added mass', 'Morison hydrodynamic loading with harmonic drag linearisation', 'Vessel motion transfer (RAO) at the hang-off', 'Quasi-static comparison (dynamic amplification)', 'Subspace-iteration eigenmodes'],
  limits: `${LIMIT_2D} Linearised response about the static shape for one regular wave: contact status is frozen, wave kinematics on the cable itself and snatch loads after slack are not modelled; a nonlinear time-domain analysis with irregular waves is needed for final design.`,
  publishes: ['cable.dynTensionAmp', 'cable.dynCurvatureAmp'],
  searchTerms: 'dynamic cable response floating wind hang-off motion tension curvature',
  inputs: [
    ...lineInputs(), ...configInputs(), ...currentInputs({ Us: 0.4, Ub: 0.1, bind: 'site.Uc1' }),
    { key: 'H', label: 'Design wave height', unit: 'm', type: 'number', default: 6, min: 0, max: 35, group: 'Wave and platform motion', role: 'bc', bind: 'site.Hs1', uncertain: { dist: 'lognormal', cov: 0.1 }, help: 'Regular-wave height; the motion amplitude at the hang-off is RAO × H/2.' },
    { key: 'T', label: 'Wave period', unit: 's', type: 'number', default: 10, min: 3, max: 30, group: 'Wave and platform motion', role: 'bc', bind: 'site.Tp1' },
    { key: 'raoHeave', label: 'Hang-off heave RAO', unit: 'm/m', type: 'number', default: 0.6, min: 0, max: 3, group: 'Wave and platform motion', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Vertical motion of the hang-off per metre of wave amplitude at the wave period (platform heave plus pitch × lever arm).' },
    { key: 'raoSurge', label: 'Hang-off surge RAO', unit: 'm/m', type: 'number', default: 0.5, min: 0, max: 3, group: 'Wave and platform motion', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'phase', label: 'Phase of heave relative to surge', unit: 'deg', type: 'number', default: 90, min: -180, max: 180, group: 'Wave and platform motion', role: 'bc' },
    { key: 'zeta', label: 'Structural damping ratio', unit: '-', type: 'number', default: 0.03, min: 0, max: 0.3, group: 'Wave and platform motion', role: 'material', help: 'Stick–slip hysteresis of the armour gives a few per cent; measured in bending tests.' },
    { key: 'v0', label: 'Initial velocity of the cable', unit: 'm/s', type: 'number', default: 0, min: 0, max: 0, group: 'Wave and platform motion', role: 'ic', help: 'Steady-state harmonic solution: the cable starts at rest in its static equilibrium and transients are not computed.' },
    { key: 'Tallow', label: 'Allowable tension', unit: 'kN', type: 'number', default: 250, min: 2, max: 5000, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.allowTension' },
    { key: 'mbrAllow', label: 'Allowable minimum bend radius (dynamic)', unit: 'm', type: 'number', default: 3.2, min: 0.3, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.mbrDynamic' },
    { key: 'nSweep', label: 'Periods in the sweep', type: 'int', default: 5, min: 0, max: 15, group: 'Solver', role: 'control', help: 'Additional wave periods solved for the response curve (0 = design period only).' },
    { key: 'nElem', label: 'Rod elements', type: 'int', default: 60, min: 30, max: 600, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Rod elements', levels: [40, 80, 160], qoi: 'Tmax' },
  compute(inp) {
    const { res, line } = solveConfig(inp, { nElem: inp.nElem, current: rodCurrent(inp, inp.rhoW) });
    const a = inp.H / 2, ph = inp.phase * DEG, n = res.s.length - 1;
    const harm = (T) => {
      if (!(a > 0) || !(inp.raoHeave + inp.raoSurge > 0)) return { T, tTop: 0, kMax: 0, sK: 0, slack: 0, kTot: Math.max(...Array.from(res.curvature, Math.abs)), nMax: 0, r: null };
      const r = run(() => rodHarmonicResponse(res, { omega: 2 * Math.PI / T, end: 'A', xAmp: inp.raoSurge * a, zAmp: inp.raoHeave * a, phaseZ: ph, rho: line.rhoW, zeta: inp.zeta }), `Harmonic response at ${T} s`);
      let kMax = 0, sK = 0, slack = 0, kTot = 0, nMax = 0;
      for (let i = 0; i <= n; i++) {
        if (r.curvatureAmp[i] > kMax) { kMax = r.curvatureAmp[i]; sK = res.s[i]; }
        kTot = Math.max(kTot, Math.abs(res.curvature[i]) + r.curvatureAmp[i]);
        slack = Math.max(slack, r.tensionAmp[i] / Math.max(res.Te[i], 1));
        nMax = Math.max(nMax, r.normalAmp[i]);
      }
      return { r, T, tTop: r.tensionAmp[0], kMax, sK, slack, kTot, nMax };
    };
    const d = harm(inp.T);
    need(!d.r || d.r.converged !== false, 'The drag linearisation did not converge at the design period: reduce the wave height or increase the damping.');
    const qs = a * inp.raoHeave > 0 ? run(() => quasiStaticHeave(res, { heaveAmp: a * inp.raoHeave }), 'Quasi-static heave case') : { tensionAmplitude: 0 };
    const daf = qs.tensionAmplitude > 1 ? d.tTop / qs.tensionAmplitude : 1;
    const modes = run(() => rodModes(res, { nModes: 4 }), 'Modal analysis');
    const Tn = modes.modes.map((m) => m.periodS);
    const sweepT = inp.nSweep > 0 ? linspace(Math.max(3, 0.5 * inp.T), 2 * inp.T, inp.nSweep) : [];
    const sweep = sweepT.map(harm);
    const Tmax = res.topTension + d.tTop, Tmin = res.topTension - d.tTop;
    const checks = [
      check('Maximum dynamic tension at the hang-off', Tmax / 1e3, inp.Tallow, 'kN', `${STD_CABLE}; CIGRE TB 862 §4`, 'Static plus dynamic top tension in the design wave against the allowable tension. Exceeding it overstrains armour wires and conductors every wave cycle.'),
      check('Effective compression (dynamic amplitude / static tension)', d.slack, 1, '-', 'DNV-RP-0360 (no compression in dynamic sections); project criterion', 'Largest ratio of dynamic tension amplitude to local static tension along the cable. At 1.0 the tension reaches zero each cycle; beyond it the cable goes slack, loops and snatches.'),
      check('Dynamic bend radius', d.kTot, 1 / inp.mbrAllow, '1/m', `${STD_CABLE} (dynamic MBR); CIGRE TB 862`, 'Static curvature plus dynamic curvature amplitude at the worst point against the allowable dynamic curvature. Exceeding it works the lead sheath and armour beyond their qualified bending range.'),
    ];
    const gov = worst(checks);
    const tA = d.r ? scaled(d.r.tensionAmp, 1e-3) : Array.from(res.s, () => 0), kA = d.r ? d.r.curvatureAmp : Array.from(res.s, () => 0), nA = d.r ? d.r.normalAmp : Array.from(res.s, () => 0);
    const [s, tAm, kAm, kS, nAm, te] = thin([res.s, tA, kA, Array.from(res.curvature, Math.abs), nA, scaled(res.Te, 1e-3)], 300);
    const cyclesYr = 31557600 / inp.T;
    return {
      kpis: [
        { key: 'Tstatic', label: 'Static top tension', value: round(res.topTension / 1e3, 6), unit: 'kN', explain: 'Hang-off tension of the static shape about which the motion is applied (initial equilibrium).' },
        { key: 'tensionAmp', label: 'Dynamic top-tension amplitude', value: round(d.tTop / 1e3, 5), unit: 'kN', explain: 'Amplitude of the hang-off tension variation in the design wave. Twice this value is the tension range per cycle that enters fatigue.' },
        { key: 'Tmax', label: 'Maximum top tension', value: round(Tmax / 1e3, 6), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Static tension plus the dynamic amplitude; compare with the allowable tension of the cable.' },
        { key: 'Tmin', label: 'Minimum top tension', value: round(Tmin / 1e3, 5), unit: 'kN', explain: 'Static tension minus the dynamic amplitude; it must stay comfortably positive.' },
        { key: 'daf', label: 'Dynamic amplification of top tension', value: round(daf), unit: '-', explain: 'Dynamic tension amplitude divided by the quasi-static amplitude for the same hang-off heave. Values well above 1 mean drag and inertia lock the cable so that motion is taken as axial stretch.' },
        { key: 'curvAmp', label: 'Largest dynamic curvature amplitude', value: round(d.kMax, 4), unit: '1/m', explain: 'Peak curvature variation along the cable (outside the bend stiffener). Twice this value is the curvature range that drives armour and sheath fatigue.' },
        { key: 'sCurv', label: 'Location of the largest curvature amplitude', value: round(d.sK), unit: 'm', explain: 'Arc length from the hang-off at which the curvature varies most — normally the hang-off region, the sag or hog bend or the touchdown.' },
        { key: 'mbrDyn', label: 'Minimum dynamic bend radius', value: round(1 / Math.max(d.kTot, 1e-9)), unit: 'm', limit: inp.mbrAllow, sense: 'min', explain: 'Reciprocal of the largest static plus dynamic curvature; it must stay above the allowable dynamic bend radius.' },
        { key: 'slackRatio', label: 'Dynamic / static tension ratio (worst point)', value: round(d.slack), unit: '-', limit: 1, sense: 'max', explain: 'Above 1 the local effective tension becomes compressive during the cycle.' },
        { key: 'motionAmp', label: 'Largest transverse motion amplitude', value: round(d.nMax), unit: 'm', explain: 'Peak movement of the cable normal to its axis; use it for clearance to the seabed, mooring lines and neighbouring cables.' },
        { key: 'Tn1', label: 'First natural period', value: round(Tn[0] ?? 0), unit: 's', explain: 'Longest natural period of the suspended cable in its plane. Wave periods close to a natural period are amplified until drag limits the motion.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Tension envelope along the cable', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Effective tension', unit: 'kN' },
          series: [{ name: 'Static', x: s, y: te }, { name: 'Maximum', x: s, y: te.map((v, i) => v + tAm[i]), dash: true }, { name: 'Minimum', x: s, y: te.map((v, i) => v - tAm[i]), dash: true }],
          explain: 'Static tension with the dynamic envelope. Where the lower envelope crosses zero the cable goes into compression each cycle.' },
        { type: 'line', title: 'Curvature: static and dynamic amplitude', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Curvature', unit: '1/m' },
          series: [{ name: 'Static |κ|', x: s, y: kS }, { name: 'Dynamic amplitude', x: s, y: kAm }, { name: 'Allowable (1/MBR)', x: [s[0], s[s.length - 1]], y: [1 / inp.mbrAllow, 1 / inp.mbrAllow], dash: true }],
          explain: 'The dynamic amplitude curve shows where fatigue concentrates; the sum of the two curves must stay below the dashed line.' },
        { type: 'line', title: 'Transverse motion amplitude', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Amplitude', unit: 'm' }, series: [{ name: 'Normal displacement amplitude', x: s, y: nAm }],
          explain: 'Dynamic cable movement normal to its axis. A lazy wave shows the motion dying out before the touchdown; a catenary carries it to the seabed.' },
        ...(sweep.length ? [{ type: 'line', title: 'Response against wave period', x: { label: 'Wave period', unit: 's' }, y: { label: 'Top-tension amplitude', unit: 'kN' },
          series: [{ name: 'Tension amplitude (same motion amplitude)', x: sweepT, y: sweep.map((q) => q.tTop / 1e3) }, { name: 'Design period', x: [inp.T, inp.T], y: [0, Math.max(d.tTop, ...sweep.map((q) => q.tTop)) / 1e3], dash: true }],
          explain: 'Dynamic top tension for the same hang-off amplitudes at other periods. A rise towards short periods is drag locking; a local peak marks a natural period.' }] : []),
      ],
      tables: [
        { title: 'Natural periods (in-plane)', explain: 'Lowest in-plane modes of the suspended cable about its static shape. Compare them with the wave-period range of the site.',
          columns: [{ key: 'mode', label: 'Mode' }, { key: 'T', label: 'Period', unit: 's' }, { key: 'f', label: 'Frequency', unit: 'Hz' }], rows: modes.modes.map((m, i) => ({ mode: i + 1, T: round(m.periodS), f: round(m.fHz) })) },
        ...(sweep.length ? [{ title: 'Period sweep', explain: 'One row per wave period with the same hang-off amplitudes: tension amplitude, largest curvature amplitude and the dynamic-to-static tension ratio.',
          columns: [{ key: 'T', label: 'Period', unit: 's' }, { key: 'dT', label: 'Tension amplitude', unit: 'kN' }, { key: 'k', label: 'Curvature amplitude', unit: '1/m' }, { key: 'slack', label: 'Slack ratio', unit: '-' }],
          rows: sweep.map((q) => ({ T: round(q.T), dT: round(q.tTop / 1e3), k: round(q.kMax), slack: round(q.slack) })) }] : []),
      ],
      notes: [
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Hang-off motion: ${fmt(a * inp.raoHeave)} m heave and ${fmt(a * inp.raoSurge)} m surge amplitude at ${inp.T} s; quasi-static tension amplitude ${fmt(qs.tensionAmplitude / 1e3)} kN, dynamic ${fmt(d.tTop / 1e3)} kN.`,
        `Natural periods ${Tn.slice(0, 3).map((t) => fmt(t)).join(', ')} s.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}: top tension ${fmt(Tmin / 1e3)}–${fmt(Tmax / 1e3)} kN, curvature amplitude ${fmt(d.kMax)} 1/m at ${fmt(d.sK)} m, amplification ${fmt(daf)}.`,
          gov.uc > 1 ? (d.slack > 1 ? 'Remove the compression: add a lazy wave or more buoyancy to decouple the lower cable, or lengthen the upper catenary so that heave is absorbed in geometry.' : 'Reduce curvature at the governing location (bend stiffener at the hang-off, longer buoyancy section at the hog bend) and re-run.')
            : 'Feed the tension and curvature ranges into the fatigue module, and confirm with an irregular-wave time-domain analysis covering the scatter diagram before design freeze.',
          `${STD_CABLE}; CIGRE TB 862 (2022) §4; DNV-RP-C205 (2021) §6`),
        adv('economic', 'info', `About ${fmt(cyclesYr / 1e6)} million wave cycles a year act on the cable; the curvature range of ${fmt(2 * d.kMax)} 1/m at ${fmt(d.sK)} m sets the fatigue life and therefore the replacement interval of the dynamic section.`,
          'Spend on configuration (buoyancy, bend stiffener) before spending on a heavier cross-section: halving the curvature range typically multiplies armour fatigue life by 8–30 (S–N slope 3–5), which defers a dynamic-cable change-out and its vessel campaign.',
          'Whole-life costing ISO 15686-5; risk-based maintenance planning DNV-RP-G101 principles',
          `A dynamic cable replacement repeats the manufacture of the whole dynamic length (about ${fmt(17 * inp.length / 1e3)} t of copper ≈ ${fmt(17 * inp.length / 1e3 * CO2_COPPER)} t CO₂e for a 3 × 630 mm² cable, screening factor ${CO2_COPPER} t CO₂e/t) plus a lay-vessel campaign: fatigue life is a sustainability driver.`),
      ],
      publish: { 'cable.dynTensionAmp': d.tTop, 'cable.dynCurvatureAmp': d.kMax },
      quantities: {},
    };
  },
  verify: [
    // Limiting case: no wave → no dynamic component, and the maximum top tension equals the static catenary value.
    // Irvine (1981) §1.2: h = 138 m, L = 300 m, X = 200 m → a ≈ 21.5 m; w = 415.7 N/m → Ttop = w·(a + h) = 66.3 kN.
    { name: 'Zero wave height: dynamic amplitude vanishes and Tmax = static catenary tension w·(a + h)',
      source: 'Limiting-case physics; Irvine, Cable Structures (1981) §1.2 (hand solution a ≈ 21.5 m)',
      inputs: { config: 'catenary', length: 300, span: 200, Us: 0, Ub: 0, H: 0, EI: 200, muAx: 0, nElem: 160, nSweep: 0, mbrAllow: 0.5, Tallow: 5000 },
      expect: [{ kpi: 'tensionAmp', value: 0, abs: 1e-6 }, { kpi: 'Tmax', value: 66.3, tol: 0.015 }] },
  ],
};
