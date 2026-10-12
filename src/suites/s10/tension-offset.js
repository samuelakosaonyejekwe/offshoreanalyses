// Maximum line tension and vessel offset by the API RP 2SK combination of mean, low-frequency and
// wave-frequency components, with API safety-factor and DNV partial-factor checks.
import { responseStatistics, combineApi2sk, quasiDynamicTension, checkApi2sk, checkDnvE301, checkOffset } from '../../engine/lines/mooringsystem.js';
import { DEG, run, round, fmt, adv, check, worst, sev, linspace, systemInputs, loadInputs, criteriaInputs, meanState, lineTable, steelAdvice,
  LIMIT_QS, LIMIT_PROPS } from './_common.js';

export default {
  id: 's10-tension-offset',
  suite: 10,
  title: 'Line tension and vessel offset',
  summary: 'Combines mean, low-frequency and wave-frequency motion into the design offset and the maximum tension of every line, and checks them against break-load and offset limits.',
  covers: [4],
  fidelity: 'reduced-order',
  standards: ['API RP 2SK (3rd ed., 2005) §5.4, §6.2–6.3', 'DNV-OS-E301 (2018) Ch.2 Sec.2 [4]', 'ISO 19901-7 (2013) §8–10'],
  models: ['Quasi-static catenary system', 'Narrow-band Gaussian response statistics (2σ significant, σ√(2 ln N) maximum)', 'API RP 2SK mean + LF + WF combination rule', 'Dynamic amplification factor on the wave-frequency tension increment', 'DNV partial-factor design equation'],
  limits: `${LIMIT_QS} Low- and wave-frequency motion standard deviations are inputs (from the coupled-analysis module, model tests or a diffraction analysis); line dynamics enter only through the amplification factor on the wave-frequency increment. ${LIMIT_PROPS}`,
  publishes: ['mooring.Tmax', 'mooring.anchorH', 'mooring.anchorV', 'mooring.dynRatio', 'vessel.offsetMax'],
  searchTerms: 'mooring line tension offset low frequency wave frequency combination API RP 2SK',
  inputs: [
    ...systemInputs(),
    ...loadInputs(),
    { key: 'lfSigma', label: 'Low-frequency motion, standard deviation', unit: 'm', type: 'number', default: 2.0, min: 0, max: 40, group: 'Vessel motions', role: 'bc', bind: 'mooring.lfSigma',
      uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Slow-drift surge/sway along the load direction. Take it from the coupled vessel–mooring module or from model tests.' },
    { key: 'lfPeriod', label: 'Low-frequency period', unit: 's', type: 'number', default: 120, min: 20, max: 1000, group: 'Vessel motions', role: 'bc', bind: 'mooring.TnSurge',
      help: 'Surge natural period on the mooring; sets the number of slow-drift cycles in the storm.' },
    { key: 'wfSigma', label: 'Wave-frequency fairlead motion, standard deviation', unit: 'm', type: 'number', default: 0.9, min: 0, max: 15, group: 'Vessel motions', role: 'bc', bind: 'mooring.wfSigma',
      uncertain: { dist: 'normal', cov: 0.1 }, help: 'Horizontal first-order motion of the fairlead along the load direction (from RAOs and the wave spectrum).' },
    { key: 'wfPeriod', label: 'Wave-frequency zero-crossing period', unit: 's', type: 'number', default: 10, min: 3, max: 30, group: 'Vessel motions', role: 'bc' },
    { key: 'stormHours', label: 'Storm duration', unit: 'h', type: 'number', default: 3, min: 0.5, max: 12, group: 'Vessel motions', role: 'environment',
      help: 'Stationary sea-state duration used for the expected maxima (3 h is customary).' },
    { key: 'daf', label: 'Dynamic amplification of WF tension', unit: '-', type: 'number', default: 1.2, min: 0.5, max: 5, group: 'Vessel motions', role: 'load', bind: 'mooring.daf',
      uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Ratio of dynamic to quasi-static wave-frequency tension amplitude, from the line-dynamics module. 1.0 is a pure quasi-static analysis (then use the quasi-static safety factor 2.0).' },
    ...criteriaInputs(),
  ],
  compute(inp) {
    const { b, sys, fs, Fext, calm, eq } = meanState(inp);
    const T = inp.stormHours * 3600;
    const lf = inp.lfSigma > 0 ? responseStatistics({ sigma: inp.lfSigma, N: Math.max(T / inp.lfPeriod, 2) }) : { sig: 0, max: 0 };
    const wf = inp.wfSigma > 0 ? responseStatistics({ sigma: inp.wfSigma, N: Math.max(T / inp.wfPeriod, 2) }) : { sig: 0, max: 0 };
    const F = Math.hypot(Fext[0], Fext[1]);
    // dynamic excursion acts along the mean offset (i.e. the load) direction
    const dir = F > 0 && eq.offset > 1e-6 ? Math.atan2(eq.y - calm.y, eq.x - calm.x) : inp.loadDirDeg * DEG;
    const comb = combineApi2sk({ mean: eq.offset, lfSig: lf.sig, lfMax: lf.max, wfSig: wf.sig, wfMax: wf.max });
    const qd = run(() => quasiDynamicTension(sys, { meanPose: eq.pose, direction: dir, lf, wf, daf: inp.daf }),
      'Dynamic excursion (a windward line reaches its elastic limit at the extreme offset — reduce the load or lengthen the lines)');
    const dExt = comb.max - eq.offset;
    const ext = fs.forces(eq.x + dExt * Math.cos(dir), eq.y + dExt * Math.sin(dir));
    const gov = qd.lines.find((l) => l.index === qd.maxTensionLine);
    const api = checkApi2sk({ tension: qd.maxTension, MBL: b.MBL, condition: 'intact', analysis: 'dynamic', requiredSF: inp.sfIntact });
    const dnv = checkDnvE301({ Tmean: gov.Tmean, Tdyn: gov.Tdyn, MBL: b.MBL, limitState: 'ULS', consequenceClass: 1, strengthFactor: inp.strengthFactor, gammaMean: inp.gammaMean, gammaDyn: inp.gammaDyn });
    const off = checkOffset({ offset: comb.max, depth: b.depth, limitFraction: inp.offsetLimitPct / 100 });
    const aH = Math.max(...ext.lines.map((l) => l.anchorH)), aV = Math.max(...ext.lines.map((l) => l.anchorV));
    const grMin = Math.min(...ext.lines.map((l) => l.groundedLength));
    const checks = [
      check('Line tension, API safety factor (intact)', qd.maxTension * api.required / 1e3, b.MBL / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3 tension limits',
        'Guards against breaking the most loaded line in the design storm: the maximum tension times the required safety factor must not exceed the catalogue break load.'),
      check('Line tension, DNV partial factors (ULS)', dnv.Td / 1e3, dnv.Sc / 1e3, 'kN', 'DNV-OS-E301 Ch.2 Sec.2 [4.2]: Sc − Tmean·γmean − Tdyn·γdyn ≥ 0',
        'Same failure mode in partial-factor format: mean and dynamic tension carry separate load factors and the strength is the characteristic 0.95·MBL.'),
      check('Maximum offset', comb.max, inp.offsetLimitPct / 100 * b.depth, 'm', 'API RP 2SK (3rd ed.) §6.2; project riser limit',
        'Guards against exceeding riser, umbilical or gangway stroke; beyond it the riser system, not the mooring, fails first.'),
    ];
    const g = worst(checks);
    const sfAch = b.MBL / qd.maxTension;
    const tensOffs = linspace(0, Math.max(comb.max * 1.15, 1), 30);
    const cv = fs.curve(dir, tensOffs, [calm.x, calm.y]);
    const iG = qd.maxTensionLine;
    const advice = [
      adv('engineering', sev(g.uc), `Maximum offset ${fmt(comb.max)} m (${fmt(100 * off.offsetFraction)} % of depth, governed by ${comb.governing}); maximum tension ${fmt(qd.maxTension / 1e3)} kN in ${sys.lines[iG].name}, safety factor ${fmt(sfAch)} against ${inp.sfIntact} required. Governing check: ${g.name} (${fmt(g.uc)}).`,
        g.uc > 1 ? (g.name.startsWith('Maximum offset') ? `Reduce the offset by about ${fmt(comb.max - inp.offsetLimitPct / 100 * b.depth)} m: raise pretension (move anchors out by a few metres), add lines, or add slow-drift damping; the low-frequency part contributes ${fmt(Math.max(lf.max, lf.sig))} m.`
          : `Raise the break load by about ${fmt(100 * (g.uc - 1))} % (diameter +${fmt(100 * (Math.sqrt(g.uc) - 1))} %, or one grade up), or soften the windward lines with more length so the offset produces less tension.`)
          : `Reserve is ${fmt(100 * (1 - g.uc))} % on the governing check. Confirm the motion standard deviations with the coupled analysis and the amplification factor ${inp.daf} with the line-dynamics module.`,
        'API RP 2SK (3rd ed., 2005) §5.4.4 combination and §6.3 limits; DNV-OS-E301 (2018) Ch.2 Sec.2 [4]. The DNV factors were revised in the 2021 edition — confirm the governing edition.'),
      adv('engineering', aV > 0 ? 'warn' : 'info', `At the extreme offset the windward anchor sees ${fmt(aH / 1e3)} kN horizontal and ${fmt(aV / 1e3)} kN vertical; ${fmt(grMin)} m of line stays grounded.`,
        aV > 0 ? 'The anchor is lifted: use a suction, pile or plate anchor rated for uplift, or add grounded chain until the vertical load vanishes.' : 'Use these anchor loads in the anchor holding-capacity module.',
        'API RP 2SK (3rd ed.) §6.5 anchor requirements; DNV-RP-E301 (2017) for fluke anchors.'),
      steelAdvice(b, Math.max(checks[0].uc, checks[1].uc), 'line-tension'),
    ];
    return {
      kpis: [
        { key: 'meanOffset', label: 'Mean offset', value: round(eq.offset, 5), unit: 'm', explain: 'Static offset under the mean wind, current and wave-drift force; the base around which the vessel oscillates.' },
        { key: 'maxOffset', label: 'Maximum offset', value: round(comb.max, 5), unit: 'm', limit: round(inp.offsetLimitPct / 100 * b.depth, 5), sense: 'max',
          explain: 'Mean offset plus the larger of (LF maximum + WF significant) and (LF significant + WF maximum). This is the excursion risers and umbilicals must survive.' },
        { key: 'maxOffsetPct', label: 'Maximum offset / depth', value: round(100 * off.offsetFraction, 4), unit: '%', explain: 'The design excursion as a share of water depth, the form in which riser limits are quoted.' },
        { key: 'Tmean', label: 'Mean tension, governing line', value: round(gov.Tmean / 1e3, 5), unit: 'kN', explain: 'Tension of the most loaded line at the mean offset; the part multiplied by γmean in the DNV format.' },
        { key: 'Tdyn', label: 'Dynamic tension, governing line', value: round(gov.Tdyn / 1e3, 5), unit: 'kN', explain: 'Increase from the mean to the maximum tension due to low- and wave-frequency motion (amplified for line dynamics); multiplied by γdyn.' },
        { key: 'Tmax', label: 'Maximum line tension', value: round(qd.maxTension / 1e3, 5), unit: 'kN', explain: 'Largest fairlead tension in the storm from the governing motion combination. Compare with the break load of the weakest segment.' },
        { key: 'sfApi', label: 'Safety factor MBL / Tmax', value: round(sfAch, 4), unit: '-', limit: inp.sfIntact, sense: 'min',
          explain: 'Break load over maximum tension. API RP 2SK asks for 1.67 intact with dynamic analysis (tension ≤ 60 % MBL).' },
        { key: 'ucDnv', label: 'DNV ULS utilisation', value: round(dnv.uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Factored tension over characteristic strength; above 1.0 the line does not meet the ultimate limit state.' },
        { key: 'lfMax', label: 'Low-frequency maximum amplitude', value: round(lf.max, 4), unit: 'm', explain: 'Expected largest slow-drift excursion in the storm, σ·√(2 ln N) with N the number of slow-drift cycles.' },
        { key: 'wfMax', label: 'Wave-frequency maximum amplitude', value: round(wf.max, 4), unit: 'm', explain: 'Expected largest first-order fairlead excursion in the storm, σ·√(2 ln N) with N the number of wave cycles.' },
        { key: 'anchorH', label: 'Anchor horizontal load at maximum offset', value: round(aH / 1e3, 5), unit: 'kN', explain: 'Horizontal design load for the windward anchor after seabed friction on the grounded line.' },
        { key: 'anchorV', label: 'Anchor vertical load at maximum offset', value: round(aV / 1e3, 5), unit: 'kN', explain: 'Uplift on the windward anchor at the extreme excursion; must be zero for a drag-embedment anchor.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Tension per line: mean, maximum and allowable', categories: qd.lines.map((l) => sys.lines[l.index].name),
          series: [{ name: 'Mean (kN)', y: qd.lines.map((l) => round(l.Tmean / 1e3)) }, { name: 'Maximum (kN)', y: qd.lines.map((l) => round(l.Tmax / 1e3)) }, { name: 'Allowable MBL/SF (kN)', y: qd.lines.map(() => round(b.MBL / inp.sfIntact / 1e3)) }],
          explain: 'Each line\'s mean and storm-maximum tension against the allowable. The windward lines carry the peak; leeward lines far below the mean show how unevenly the load is shared.' },
        { type: 'line', title: 'Governing-line tension vs offset', x: { label: 'Offset along load direction', unit: 'm' }, y: { label: 'Fairlead tension', unit: 'kN' },
          series: [{ name: sys.lines[iG].name, x: tensOffs.map((v) => round(v)), y: cv.lineTensions.map((r) => round(r[iG] / 1e3)) },
            { name: 'Mean offset', x: [eq.offset, eq.offset], y: [0, gov.Tmean / 1e3], dash: true }, { name: 'Maximum offset', x: [comb.max, comb.max], y: [0, qd.maxTension / 1e3], dash: true },
            { name: 'Allowable', x: [0, tensOffs[tensOffs.length - 1]], y: [b.MBL / inp.sfIntact / 1e3, b.MBL / inp.sfIntact / 1e3], dash: true }],
          explain: 'Quasi-static tension of the most loaded line as the vessel moves downwind. The steepening curve shows why the last metres of excursion cost the most tension; the crossing with the allowable is the offset the line can tolerate.' },
        { type: 'bar', title: 'Offset components', categories: ['Mean', 'LF significant', 'LF maximum', 'WF significant', 'WF maximum', 'Combined maximum'],
          series: [{ name: 'Offset (m)', y: [eq.offset, lf.sig, lf.max, wf.sig, wf.max, comb.max].map((v) => round(v)) }],
          explain: 'What the design offset is made of. A dominant low-frequency bar points to slow-drift damping or stiffness as the lever; a dominant mean bar points to pretension.' },
      ],
      tables: [
        { title: 'Line tensions by motion combination', columns: [{ key: 'line', label: 'Line' }, { key: 'Tmean', label: 'Mean', unit: 'kN' }, { key: 'c1', label: 'Mean + LFmax + WFsig', unit: 'kN' }, { key: 'c2', label: 'Mean + LFsig + WFmax', unit: 'kN' }, { key: 'Tmax', label: 'Maximum', unit: 'kN' }, { key: 'sf', label: 'MBL / Tmax' }, { key: 'gov', label: 'Governing combination' }],
          rows: qd.lines.map((l) => ({ line: sys.lines[l.index].name, Tmean: round(l.Tmean / 1e3), c1: round(l.case1 / 1e3), c2: round(l.case2 / 1e3), Tmax: round(l.Tmax / 1e3), sf: round(b.MBL / l.Tmax, 3), gov: l.governing })),
          explain: 'For every line the tension at the mean position and in the two API motion combinations. The safety-factor column identifies the critical lines to report.' },
        lineTable(ext.lines, b.MBL, 'Quasi-static line and anchor loads at the maximum offset'),
      ],
      notes: [`${comb.governing} governs the offset; ${sys.lines[iG].name} governs tension with ${gov.governing}.`,
        `Storm of ${inp.stormHours} h gives ${fmt(T / inp.lfPeriod)} low-frequency and ${fmt(T / inp.wfPeriod)} wave-frequency cycles.`,
        'Safety factors and partial factors are overridable defaults; confirm them against the governing edition of the standard.'],
      advice,
      publish: { 'mooring.Tmax': qd.maxTension, 'mooring.anchorH': aH, 'mooring.anchorV': aV, 'mooring.dynRatio': gov.Tmean > 0 ? gov.Tdyn / gov.Tmean : 0, 'vessel.offsetMax': comb.max },
      quantities: { steel_t: round(b.steel_t) },
    };
  },
  verify: [
    { // No mean load. σLF = 2 m, 3 h / 120 s = 90 cycles → LFmax = 2·sqrt(2 ln 90) = 6.000 m, LFsig = 4 m;
      // σWF = 1 m, 10800 / 10.8 = 1000 cycles → WFmax = sqrt(2 ln 1000) = 3.717 m, WFsig = 2 m.
      // Combination: max(6.000 + 2, 4 + 3.717) = 8.000 m.
      name: 'API RP 2SK offset combination with zero mean load',
      source: 'API RP 2SK (3rd ed.) §5.4.4; hand calculation with σ√(2 ln N)',
      inputs: { Fmean: 0, lfSigma: 2, lfPeriod: 120, wfSigma: 1, wfPeriod: 10.8, stormHours: 3 },
      expect: [{ kpi: 'maxOffset', value: 8.0, tol: 0.001 }, { kpi: 'wfMax', value: 3.717, tol: 0.001 }, { kpi: 'meanOffset', value: 0, abs: 1e-3 }] },
  ],
};
