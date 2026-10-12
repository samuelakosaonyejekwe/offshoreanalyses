// Static equilibrium of a moored floater under a mean environmental load: catenary shapes, pretensions,
// anchor loads, restoring curve and linearised stiffness (API RP 2SK §5.3 mean-offset analysis).
import { equilibrium, restoringCurve, stiffnessMatrix } from '../../engine/lines/mooringsystem.js';
import { DEG, need, run, round, fmt, adv, check, worst, sev, linspace, systemInputs, loadInputs, fextFrom, buildSystem, planPlot, elevationPlot, lineTable,
  steelAdvice, vesselCo2Text, LIMIT_QS, LIMIT_PROPS, CONFIRM } from './_common.js';

export default {
  id: 's10-static-equilibrium',
  suite: 10,
  title: 'Static mooring equilibrium',
  summary: 'Finds the vessel position where the mooring balances the mean environmental load and reports line shapes, pretensions, anchor loads and the restoring stiffness.',
  covers: [1],
  fidelity: 'analytical',
  standards: ['API RP 2SK (3rd ed., 2005) §5.3', 'DNV-OS-E301 (2018) Ch.2 Sec.2', 'ISO 19901-7 (2013) §8'],
  models: ['Elastic catenary (Irvine 1981; Jonkman 2007)', 'Multi-segment line with clump weights and buoys', 'Force and moment equilibrium by damped Newton iteration', 'Seabed friction on the grounded length', 'Linearised stiffness by central differences'],
  limits: `${LIMIT_QS} ${LIMIT_PROPS} Dynamic tension, low-frequency motion and line bending stiffness are outside this module; three-dimensional rod and Cosserat formulations are not used.`,
  publishes: ['mooring.pretension', 'mooring.Kxx', 'mooring.anchorH', 'mooring.anchorV', 'mooring.MBL', 'mooring.lineLength', 'quantities.steel_t'],
  searchTerms: 'quasi-static mooring analysis catenary equilibrium restoring force',
  inputs: [
    ...systemInputs(),
    ...loadInputs(),
    { key: 'x0', label: 'Initial vessel position x', unit: 'm', type: 'number', default: 0, min: -500, max: 500, group: 'Initial condition', role: 'ic',
      help: 'Starting position of the Newton iteration and the reference from which the offset is measured (the as-installed, calm-water position).' },
    { key: 'y0', label: 'Initial vessel position y', unit: 'm', type: 'number', default: 0, min: -500, max: 500, group: 'Initial condition', role: 'ic' },
    { key: 'yaw0Deg', label: 'Initial vessel heading', unit: 'deg', type: 'number', default: 0, min: -180, max: 180, group: 'Initial condition', role: 'ic' },
    { key: 'freeYaw', label: 'Yaw free to respond', type: 'bool', default: true, group: 'Initial condition', role: 'bc',
      help: 'On: surge, sway and yaw are solved. Off: heading is held (e.g. by thrusters) and only surge and sway are solved.' },
    { key: 'sfStatic', label: 'Required safety factor on mean tension', unit: '-', type: 'number', default: 2.0, min: 1, max: 5, group: 'Acceptance criteria', role: 'criteria',
      help: `API RP 2SK quasi-static intact factor 2.0 is used as a first screen on the mean tension only. ${CONFIRM}` },
    { key: 'offsetLimitPct', label: 'Allowable mean offset', unit: '% depth', type: 'number', default: 5, min: 0.5, max: 60, group: 'Acceptance criteria', role: 'criteria',
      help: 'Project limit on the mean (static) offset; the dynamic excursion is added in the tension-and-offset module.' },
    { key: 'upliftAllowDeg', label: 'Allowable anchor uplift angle', unit: 'deg', type: 'number', default: 0, min: 0, max: 60, group: 'Acceptance criteria', role: 'criteria',
      help: 'Zero for drag-embedment anchors; suction, pile and plate anchors accept uplift (set their design angle).' },
  ],
  compute(inp) {
    const b = buildSystem(inp);
    const sys = b.system;
    const Fext = fextFrom(inp);                       // kN → N
    const pose0 = [inp.x0, inp.y0, 0, 0, 0, inp.yaw0Deg * DEG];
    const dofs = inp.freeYaw ? [0, 1, 5] : [0, 1];
    const calm = run(() => equilibrium(sys, { dofs, pose0 }), 'Calm-water equilibrium');
    const eq = run(() => equilibrium(sys, { Fext, dofs, pose0: calm.pose }), 'Loaded equilibrium');
    const K = run(() => stiffnessMatrix(sys, { pose: calm.pose }), 'Stiffness');
    const live = eq.lines;
    const pre = calm.lines.map((l) => l.T);
    const preMean = pre.reduce((a, v) => a + v, 0) / pre.length;
    const gov = live[eq.maxTensionLine];
    const aH = Math.max(...live.map((l) => l.anchorH)), aV = Math.max(...live.map((l) => l.anchorV));
    const upMax = Math.max(...live.map((l) => l.upliftAngle)), grMin = Math.min(...live.map((l) => l.groundedLength));
    const F = Math.hypot(Fext[0], Fext[1]);
    const dirLoad = inp.loadDirDeg * DEG;
    // restoring curve along the load direction out to twice the solved offset (at least 4 % of depth)
    const dMax = Math.max(2 * eq.offset, 0.04 * b.depth);
    let curve;
    const offs = linspace(0, dMax, 25);
    try { curve = restoringCurve(sys, { direction: dirLoad, offsets: offs, pose0: calm.pose }); } catch (e) {
      // the far end of the sweep may exceed the taut limit: fall back to the solved range only
      curve = run(() => restoringCurve(sys, { direction: dirLoad, offsets: linspace(0, Math.max(eq.offset, 0.01 * b.depth), 25), pose0: calm.pose }), 'Restoring curve');
    }
    const residual = Math.hypot(eq.F[0] + Fext[0], eq.F[1] + Fext[1]);
    const util = eq.maxTension / b.MBL;
    const checks = [
      check('Mean tension against break load (quasi-static screen)', eq.maxTension * inp.sfStatic / 1e3, b.MBL / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3, quasi-static intact factor',
        'Guards against line breakage. Only the mean load is included here, so a pass is necessary but not sufficient: dynamic tension is added in the tension-and-offset module.'),
      check('Mean offset', eq.offset, inp.offsetLimitPct / 100 * b.depth, 'm', 'Project riser / umbilical offset limit (API RP 2SK §6.2)',
        'Guards against over-bending or over-tensioning risers and umbilicals; beyond the limit the riser stroke or bend radius is exceeded.'),
      { name: 'Anchor uplift angle', demand: round(upMax / DEG, 4), capacity: inp.upliftAllowDeg, unit: 'deg', uc: inp.upliftAllowDeg > 0 ? upMax / DEG / inp.upliftAllowDeg : (upMax > 1e-6 ? 99 : 0), ref: 'API RP 2SK (3rd ed.) §6.5',
        explain: 'A drag anchor loses holding power if the line lifts it; a non-zero angle with zero allowable means the grounded length has been used up and the anchor is pulled upwards.' },
    ];
    const g = worst(checks);
    const quant = { steel_t: round(b.steel_t), vessel_days: { mooring_hookup: round(0.75 * b.nLines + 2, 3) } };
    if (b.polymer_t > 0) quant.polymer_t = round(b.polymer_t);
    const advice = [];
    advice.push(adv('engineering', sev(g.uc), `Mean offset ${fmt(eq.offset)} m (${fmt(100 * eq.offset / b.depth)} % of depth) with ${gov.name} at ${fmt(eq.maxTension / 1e3)} kN = ${fmt(100 * util)} % of MBL; governing check: ${g.name} (utilisation ${fmt(g.uc)}).`,
      g.uc > 1 ? (g.name.startsWith('Mean offset') ? `Stiffen the system: raising pretension by moving anchors out or shortening lines increases Kxx (now ${fmt(K[0][0] / 1e3)} kN/m); roughly ${fmt(F / 1e3 / (inp.offsetLimitPct / 100 * b.depth))} kN/m secant stiffness is needed for this load.`
        : g.name.startsWith('Anchor') ? `Add grounded chain or move the anchors out: the windward line needs about ${fmt(Math.max(50, 0.1 * b.Ltot))} m more grounded length to keep the anchor load horizontal.`
          : `Increase the weakest segment's diameter by about ${fmt(100 * (Math.sqrt(g.uc) - 1))} % (strength ∝ d²) or move to a higher grade.`)
        : `Carry this mean state into the dynamic tension-and-offset check; smallest grounded length is ${fmt(grMin)} m, so ${grMin > 50 ? 'the anchors see horizontal load only.' : 'there is little reserve before anchor uplift.'}`,
      'API RP 2SK (3rd ed., 2005) §5.3 and §6; DNV-OS-E301 (2018) Ch.2 Sec.2 [2].'));
    advice.push(adv('engineering', 'info', `Force balance closes to ${fmt(residual)} N on a ${fmt(F / 1e3)} kN load; the lines carry no current or wave drag in this solver.`,
      'For deep water or strong current, add the line drag to the mean force (or check with the line-dynamics module with current) before relying on the offset.',
      'ISO 19901-7 (2013) §8.3: line drag to be considered where significant.'));
    advice.push(steelAdvice(b, Math.max(checks[0].uc, 0.01), 'mean-tension'));
    advice.push(adv('economic', 'info', `Hook-up of ${b.nLines} lines is taken as ${fmt(0.75 * b.nLines + 2)} anchor-handling vessel days (screening: 0.75 day per line plus 2 days mobilisation).`,
      'Fewer, stronger lines shorten the offshore campaign; check the redundancy module before reducing the line count.',
      'Marine-operations planning per DNV-ST-N001 (2021) §2 (operation duration and weather windows); AACE RP 18R-97 Class 4.',
      vesselCo2Text(0.75 * b.nLines + 2)));
    return {
      kpis: [
        { key: 'pretension', label: 'Calm-water pretension (mean of lines)', value: round(preMean / 1e3, 5), unit: 'kN',
          explain: 'Fairlead tension with no environmental load. Typically 10–20 % of MBL for catenary systems; higher pretension stiffens the mooring but uses up strength and fatigue margin.' },
        { key: 'offset', label: 'Mean offset under load', value: round(eq.offset, 5), unit: 'm', limit: round(inp.offsetLimitPct / 100 * b.depth, 5), sense: 'max',
          explain: 'Horizontal movement of the vessel from its calm-water position to the point where line forces balance the mean environmental force.' },
        { key: 'offsetPct', label: 'Mean offset / water depth', value: round(100 * eq.offset / b.depth, 4), unit: '%',
          explain: 'Offset expressed as a share of depth, the usual way riser limits are written. Rigid-riser systems normally need single digits.' },
        { key: 'Tmax', label: 'Largest mean line tension', value: round(eq.maxTension / 1e3, 5), unit: 'kN',
          explain: 'Fairlead tension of the most loaded (windward) line under the mean load only. Dynamic tension comes on top.' },
        { key: 'util', label: 'Largest mean tension / MBL', value: round(util, 4), unit: '-', limit: round(1 / inp.sfStatic, 4), sense: 'max',
          explain: 'Share of the weakest segment\'s breaking load used by the mean tension. The limit is the reciprocal of the required safety factor.' },
        { key: 'Kxx', label: 'Surge stiffness at calm-water position', value: round(K[0][0] / 1e3, 5), unit: 'kN/m',
          explain: 'Restoring force per metre of surge offset. It sets the slow-drift natural period (with vessel mass) and how far a given load pushes the vessel.' },
        { key: 'Kyy', label: 'Sway stiffness at calm-water position', value: round(K[1][1] / 1e3, 5), unit: 'kN/m',
          explain: 'Restoring force per metre of sway offset; equal to the surge value for a symmetric spread.' },
        { key: 'Kyaw', label: 'Yaw stiffness', value: round(K[5][5] / 1e6, 5), unit: 'MN·m/rad',
          explain: 'Restoring moment per radian of heading change; near zero for a turret, which is what lets the hull weathervane.' },
        { key: 'anchorH', label: 'Largest anchor horizontal load', value: round(aH / 1e3, 5), unit: 'kN',
          explain: 'Horizontal force delivered to the windward anchor after seabed friction on the grounded length; the demand for the anchor holding-capacity check.' },
        { key: 'anchorV', label: 'Largest anchor vertical load', value: round(aV / 1e3, 5), unit: 'kN',
          explain: 'Uplift at the anchor. Zero while line still rests on the seabed; any value above zero rules out a conventional drag anchor.' },
        { key: 'grounded', label: 'Smallest grounded length', value: round(grMin, 5), unit: 'm',
          explain: 'Line length still on the seabed for the most loaded line. This is the reserve against anchor uplift.' },
        { key: 'residual', label: 'Force-balance residual', value: round(residual, 3), unit: 'N',
          explain: 'Unbalanced horizontal force left by the solver; it should be negligible against the applied load (verification of equilibrium).' },
      ],
      checks,
      plots: [
        planPlot(sys, eq.pose),
        elevationPlot(sys, eq.pose, [eq.maxTensionLine, live.reduce((m, l, i) => (l.T < live[m].T ? i : m), 0)], 'Most and least loaded line (elevation)'),
        { type: 'line', title: 'Restoring force and line tension vs offset', x: { label: 'Offset along load direction', unit: 'm' }, y: { label: 'Force', unit: 'kN' },
          series: [{ name: 'Restoring force', x: curve.offset.map((v) => round(v)), y: curve.restoring.map((v) => round(v / 1e3)) },
            { name: 'Largest line tension', x: curve.offset.map((v) => round(v)), y: curve.maxTension.map((v) => round(v / 1e3)) },
            { name: 'Applied mean force', x: [0, curve.offset[curve.offset.length - 1]], y: [F / 1e3, F / 1e3], dash: true }],
          explain: 'The mooring force–offset characteristic. The vessel sits where the restoring curve crosses the applied force; a curve that steepens sharply means the lines are going taut and tension will climb fast with any further offset.' },
        { type: 'bar', title: 'Line tensions: calm water vs loaded', categories: live.map((l) => l.name), series: [{ name: 'Pretension (kN)', y: pre.map((v) => round(v / 1e3)) }, { name: 'Under mean load (kN)', y: live.map((l) => round(l.T / 1e3)) }],
          explain: 'How the mean load redistributes tension: windward lines pick up load, leeward lines slacken. A leeward line falling close to zero is a snap-load candidate.' },
      ],
      tables: [
        lineTable(live, b.MBL, 'Line results under mean load'),
        { title: 'Line make-up used', columns: [{ key: 'name', label: 'Segment' }, { key: 'type', label: 'Type' }, { key: 'L', label: 'Length', unit: 'm' }, { key: 'w', label: 'Submerged weight', unit: 'N/m' }, { key: 'EA', label: 'EA', unit: 'MN' }, { key: 'MBL', label: 'MBL', unit: 'kN' }],
          rows: b.segs.map((s) => ({ name: s.name, type: s.type, L: s.L, w: round(s.w), EA: round(s.EA / 1e6), MBL: round(s.MBL / 1e3) })),
          explain: 'Properties actually used for each segment, from the estimators or your overrides. Check them against the manufacturer certificates before design use.' },
        { title: 'Horizontal stiffness matrix at the calm-water position', columns: [{ key: 'row', label: '' }, { key: 'x', label: 'surge', unit: 'kN/m | kN/rad' }, { key: 'y', label: 'sway', unit: 'kN/m | kN/rad' }, { key: 'yaw', label: 'yaw', unit: 'kN | kN·m/rad' }],
          rows: [0, 1, 5].map((i) => ({ row: ['surge', 'sway', '', '', '', 'yaw'][i], x: round(K[i][0] / 1e3), y: round(K[i][1] / 1e3), yaw: round(K[i][5] / 1e3) })),
          explain: 'Linearised restoring matrix for surge, sway and yaw. Off-diagonal terms show coupling (non-zero for asymmetric layouts); use it as the mooring stiffness in vessel-motion analyses.' },
      ],
      notes: [`${gov.name} governs with ${fmt(eq.maxTension / 1e3)} kN; the mean load of ${fmt(F / 1e3)} kN moves the vessel ${fmt(eq.offset)} m.`,
        `Weakest segment MBL ${fmt(b.MBL / 1e3)} kN; total line length ${fmt(b.Ltot, 5)} m per line.`,
        'Line properties are indicative estimators unless overridden in the make-up table.'],
      advice,
      publish: { 'mooring.pretension': preMean, 'mooring.Kxx': K[0][0], 'mooring.anchorH': aH, 'mooring.anchorV': aV, 'mooring.MBL': b.MBL, 'mooring.lineLength': b.Ltot, 'quantities.steel_t': b.steel_t },
      quantities: quant,
    };
  },
  verify: [
    { name: 'Zero load on a symmetric spread: no offset',
      source: 'Limiting case: symmetric pretensions cancel, so equilibrium is at the reference position (statics).',
      inputs: { Fmean: 0 }, expect: [{ kpi: 'offset', value: 0, abs: 1e-3 }] },
    { // Inextensible catenary with seabed: w = 1500 N/m, fairlead 280 m above seabed, H = 800 kN → a = H/w = 533.33 m;
      // suspended s = sqrt(h² + 2ha) = 614.06 m; x = a·acosh(1 + h/a) = 525.04 m; with L = 1000 m the grounded length is
      // 385.94 m, so span X = 910.98 m and fairlead tension T = H + w·h = 800 + 420 = 1220 kN (no friction).
      name: 'Pretension of an inextensible chain catenary, T = H + w·h',
      source: 'Faltinsen (1990) Sea Loads on Ships and Offshore Structures, Eqs. 8.14–8.20',
      inputs: { Fmean: 0, depth: 300, fairleadDraft: 20, fairleadRadius: 40, anchorRadius: 950.98, CB: 0,
        segments: [{ name: 'chain', type: 'custom', d: 0, grade: '', L: 1000, w: 1500, EA: 1e7, MBL: 9000 }] },
      expect: [{ kpi: 'pretension', value: 1220, tol: 0.002 }, { kpi: 'grounded', value: 385.94, tol: 0.003 }] },
  ],
};
