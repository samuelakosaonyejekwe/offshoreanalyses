// Mooring and station-keeping of a floating wind turbine (or other floating renewable device): quasi-static
// catenary spread mooring — pretension, restoring curve and stiffness, mean and maximum offset, line tensions in
// the intact and one-line-broken conditions, partial-factor tension checks and anchor holding capacity.
import { spreadLayout, equilibrium, lineForces, restoringCurve, stiffnessMatrix, brokenLineAnalysis, checkDnvE301, naturalPeriod } from '../../engine/lines/mooringsystem.js';
import { chainProperties, toSegment } from '../../engine/lines/lineprops.js';
import { dragAnchorUhc, anchorDesignCheck, DRAG_ANCHOR_INDICATIVE } from '../../engine/geotech/anchors.js';
import { caissonHorizontalCapacityClay } from '../../engine/geotech/suction.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, DEG, LIMIT_CODE, CONFIRM, steelCo2Text, vesselCo2Text } from './_common.js';

const drag = { key: 'anchorType', in: ['drag'] }, suction = { key: 'anchorType', in: ['suction'] };
const SOILS = Object.keys(DRAG_ANCHOR_INDICATIVE);

export default {
  id: 's16-mooring', suite: 16, title: 'Mooring and station-keeping (catenary spread)',
  summary: 'Pretension, stiffness, offsets and line tensions of a catenary spread mooring for a floating wind platform, intact and with one line broken, with DNV partial-factor checks and anchor holding capacity.',
  covers: [8], fidelity: 'analytical',
  standards: ['DNV-ST-0119 (2021) Floating wind turbine structures §8 (station keeping)', 'DNV-OS-E301 (2018) Position mooring Ch.2 Sec.2 (design equation and partial factors)', 'IEC TS 61400-3-2:2019 §14 (station-keeping systems)', 'API RP 2SK (2005, reaff. 2015) §5–7 (quasi-static analysis, damaged condition)', 'DNV-RP-E301 / DNV-RP-E303 (anchors)'],
  models: ['Elastic catenary mooring equations with seabed friction', 'Quasi-static equilibrium of the moored floater (surge, sway, yaw)', 'Linearised mooring stiffness matrix', 'One-line-broken equilibrium with work–energy transient overshoot', 'DNV partial-factor tension format', 'Drag-embedment anchor holding power law', 'Suction-anchor lateral capacity in clay (Murff & Hamilton)', 'Single-degree-of-freedom surge natural period'],
  limits: `Quasi-static analysis: line dynamics (drag and inertia of the line, snap loads), low-frequency platform response and vortex-induced motion are represented only through the user-entered dynamic offset. Chain properties are indicative estimators; anchor capacities are indicative chart-type values to be replaced by the anchor supplier's design. Fatigue of the lines is not covered here. ${LIMIT_CODE}`,
  publishes: ['mooring.Tmax', 'mooring.Kxx', 'mooring.anchorH', 'mooring.anchorV', 'mooring.pretension', 'vessel.offsetMax', 'mooring.MBL', 'mooring.TmaxDamaged', 'mooring.TnSurge'],
  searchTerms: 'floating offshore wind catenary mooring quasi-static DNV-ST-0119 anchor holding capacity',
  inputs: [
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 200, min: 40, max: 2000, group: 'Site', role: 'environment', bind: 'site.depth', help: 'Catenary chain moorings suit about 60–400 m; shallower sites need heavy chain, clump weights or synthetic ropes.' },
    { key: 'nLines', label: 'Mooring lines', type: 'int', default: 6, min: 3, max: 12, group: 'Mooring layout', role: 'geometry', design: true, help: 'Equally spaced in azimuth. Three lines have no redundancy: loss of one lets the platform drift on the remaining pair.' },
    { key: 'anchorRadius', label: 'Anchor radius', unit: 'm', type: 'number', default: 900, min: 80, max: 6000, group: 'Mooring layout', role: 'bc', design: true, help: 'Horizontal distance of the anchors from the platform centre (fixed seabed end of each line).' },
    { key: 'fairleadRadius', label: 'Fairlead radius', unit: 'm', type: 'number', default: 40.9, min: 2, max: 100, group: 'Mooring layout', role: 'geometry' },
    { key: 'fairleadDepth', label: 'Fairlead depth below still water', unit: 'm', type: 'number', default: 14, min: 0, max: 150, group: 'Mooring layout', role: 'geometry' },
    { key: 'lineLength', label: 'Unstretched line length', unit: 'm', type: 'number', default: 905, min: 60, max: 7000, group: 'Mooring line', role: 'ic', design: true, help: 'Sets the pretension: the initial line tension must be consistent with the calm-water platform position.' },
    { key: 'chain', label: 'Chain grade', type: 'select', source: 'materials:chain', default: 'chain-r4', group: 'Mooring line', role: 'material' },
    { key: 'dChain', label: 'Chain nominal diameter', unit: 'mm', type: 'number', default: 95, min: 40, max: 220, group: 'Mooring line', role: 'geometry', design: true, help: 'Studless chain; mass ≈ 19.9·d² t/m (d in m).' },
    { key: 'corr', label: 'Corrosion and wear allowance on diameter', unit: 'mm', type: 'number', default: 8, min: 0, max: 40, group: 'Mooring line', role: 'material', help: 'Design life × rate (0.2–0.8 mm/yr by zone, DNV-OS-E301): the break load is taken on the reduced diameter.' },
    { key: 'cb', label: 'Seabed friction coefficient', unit: '-', type: 'number', default: 0.6, min: 0, max: 1.2, group: 'Mooring line', role: 'soil', help: 'Chain on sand about 0.7–1.0, on clay 0.5–0.7 (lower bound for anchor load).' },
    { key: 'thrust', label: 'Rotor thrust (mean, at rated)', unit: 'kN', type: 'number', default: 750, min: 0, max: 10000, group: 'Mean environmental load', role: 'bc', bind: 'owt.thrustRated', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'Fenv', label: 'Mean wave drift, current and wind force on the platform', unit: 'kN', type: 'number', default: 450, min: 0, max: 20000, group: 'Mean environmental load', role: 'bc', bind: 'vessel.Fmean', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Steady force on hull and tower in the design sea state, co-linear with the thrust.' },
    { key: 'dirDeg', label: 'Load direction relative to line 1', unit: 'deg', type: 'number', default: 0, min: 0, max: 180, group: 'Mean environmental load', role: 'bc', help: '0° pulls directly against line 1 (highest line tension); 180/n degrees loads two lines equally.' },
    { key: 'dynOffset', label: 'Dynamic offset amplitude (wave- and low-frequency motion)', unit: 'm', type: 'number', default: 6, min: 0, max: 60, group: 'Mean environmental load', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Most probable maximum excursion beyond the mean offset in the design storm, from the platform motion analysis.' },
    { key: 'daf', label: 'Dynamic amplification of the dynamic tension', unit: '-', type: 'number', default: 1.3, min: 1, max: 3, group: 'Mean environmental load', role: 'bc', help: 'Ratio of dynamic to quasi-static line tension increment (line drag and inertia); 1.1–1.6 for catenary chain in moderate depth.' },
    { key: 'pMass', label: 'Platform mass (displacement)', unit: 't', type: 'number', default: 14070, min: 200, max: 100000, group: 'Platform', role: 'geometry', bind: 'vessel.disp' },
    { key: 'aSurge', label: 'Surge added mass', unit: 't', type: 'number', default: 8500, min: 0, max: 100000, group: 'Platform', role: 'bc' },
    { key: 'anchorType', label: 'Anchor type', type: 'select', default: 'drag', group: 'Anchors', role: 'bc', options: [{ value: 'drag', label: 'Drag-embedment anchor' }, { value: 'suction', label: 'Suction anchor in clay' }] },
    { key: 'anchorMass', label: 'Drag-anchor mass', unit: 't', type: 'number', default: 18, min: 1, max: 80, group: 'Anchors', role: 'geometry', showIf: drag, design: true },
    { key: 'anchorSoil', label: 'Anchor / soil class', type: 'select', default: SOILS[1], group: 'Anchors', role: 'soil', showIf: drag, options: SOILS.map((s) => ({ value: s, label: s })), help: 'Indicative holding-power coefficients of high-efficiency fluke anchors — use the supplier chart for design.' },
    { key: 'sucD', label: 'Suction-anchor diameter', unit: 'm', type: 'number', default: 5, min: 1, max: 12, group: 'Anchors', role: 'geometry', showIf: suction, design: true },
    { key: 'sucL', label: 'Suction-anchor embedded length', unit: 'm', type: 'number', default: 15, min: 2, max: 40, group: 'Anchors', role: 'geometry', showIf: suction, design: true },
    { key: 'su0', label: 'Undrained shear strength at the mudline', unit: 'kPa', type: 'number', default: 5, min: 0, max: 200, group: 'Anchors', role: 'soil', showIf: suction, bind: 'site.su0' },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 1.8, min: 0, max: 20, group: 'Anchors', role: 'soil', showIf: suction, bind: 'site.suGrad' },
    { key: 'cc', label: 'Consequence class', type: 'select', default: '1', group: 'Acceptance', role: 'criteria', options: [{ value: '1', label: 'Class 1 (unmanned, no unacceptable consequence of failure)' }, { value: '2', label: 'Class 2' }], help: `DNV-ST-0119 normally assigns class 1 to unmanned floating wind units with redundant mooring. ${CONFIRM}` },
    { key: 'anchorSF', label: 'Required anchor safety factor (intact)', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance', role: 'criteria', help: `API RP 2SK: 1.5 intact / 1.0 damaged for drag anchors (permanent moorings); DNV uses partial factors. ${CONFIRM}` },
    { key: 'offsetLimit', label: 'Allowable platform offset / water depth', unit: '%', type: 'number', default: 15, min: 2, max: 60, group: 'Acceptance', role: 'criteria', help: 'Set by the dynamic power cable (typically 10–30 % of depth, larger in shallow water).' },
  ],
  compute(inp, ctx) {
    const cm = ctx.material(inp.chain);
    need(cm && cm.grade, `Chain grade '${inp.chain}' is not in the materials library.`);
    const d = inp.dChain / 1e3, dNet = (inp.dChain - inp.corr) / 1e3;
    need(dNet > 0.02, 'The corrosion allowance leaves less than 20 mm of chain diameter.');
    const ch = chainProperties({ d, grade: cm.grade }), MBL = chainProperties({ d: dNet, grade: cm.grade }).MBL;
    const h = inp.depth - inp.fairleadDepth, span = inp.anchorRadius - inp.fairleadRadius;
    need(span > 10 && h > 5, 'The anchors must lie outside the fairleads and the fairleads above the seabed.');
    need(inp.lineLength > Math.hypot(span, h) * 1.0005, `Under-constrained geometry: the line (${inp.lineLength} m) is shorter than the straight distance from anchor to fairlead (${fmt(Math.hypot(span, h), 5)} m).`);
    need(inp.lineLength < span + h, `Over-length: a ${inp.lineLength} m line between an anchor ${fmt(span, 4)} m away and a fairlead ${fmt(h, 4)} m above the seabed hangs vertically with slack on the seabed (limit ${fmt(span + h, 5)} m). Shorten the line or increase the anchor radius.`);
    const sys = spreadLayout({ nLines: inp.nLines, anchorRadius: inp.anchorRadius, fairleadRadius: inp.fairleadRadius, fairleadZ: -inp.fairleadDepth, depth: inp.depth, heading: Math.PI, segments: [toSegment(ch, inp.lineLength, 'chain')], MBL, CB: inp.cb });
    const e0 = run(() => equilibrium(sys, {}), 'Calm-water equilibrium');
    const L0 = e0.lines[0], pre = L0.T;
    const F = (inp.thrust + inp.Fenv) * 1e3, a = inp.dirDeg * DEG;
    // line 1 runs towards −x: a load in +x at 0° pulls directly against it
    const Fext = [F * Math.cos(a), F * Math.sin(a), 0, 0, 0, 0];
    const eq = run(() => equilibrium(sys, { Fext }), 'Equilibrium under the mean load (the mooring cannot balance this force — increase chain size or pretension)');
    need(eq.converged, 'The mean-load equilibrium did not converge: the mooring is too soft or too slack for the applied force.');
    const mean = eq.offset, xMax = mean + inp.dynOffset;
    const offs = linspace(0, Math.max(xMax * 1.25, 10), 26);
    const rc = run(() => restoringCurve(sys, { direction: a, offsets: [...offs, mean, xMax] }), 'Restoring curve (a line is pulled taut within the plotted offsets — shorten the offset range or lengthen the lines)');
    const nO = offs.length, Tmean = rc.maxTension[nO], TqsMax = rc.maxTension[nO + 1];
    const Tdyn = (TqsMax - Tmean) * inp.daf, Tmax = Tmean + Tdyn;
    const K = run(() => stiffnessMatrix(sys, {}), 'Mooring stiffness');
    const Kxx = K[0][0], Tn = naturalPeriod({ mass: inp.pMass * 1e3, addedMass: inp.aSurge * 1e3, stiffness: Kxx });
    const uls = checkDnvE301({ Tmean, Tdyn, MBL, limitState: 'ULS', consequenceClass: Number(inp.cc), analysis: 'dynamic' });
    // one line broken: the most loaded line is removed
    let als = null, dmg = null;
    try {
      const b = brokenLineAnalysis(sys, { Fext, lines: [eq.maxTensionLine], overshootReduction: 0.7 });
      dmg = b.damaged[0];
      if (dmg && dmg.equilibrium.converged) { const Tm = dmg.equilibrium.maxTension; als = { ...checkDnvE301({ Tmean: Tm, Tdyn: Math.max(dmg.maxTransientTension - Tm, Tdyn), MBL, limitState: 'ALS', consequenceClass: Number(inp.cc), analysis: 'dynamic' }), offset: Math.max(dmg.transientOffset, dmg.equilibrium.offset), Tmax: Math.max(dmg.maxTransientTension, Tm + Tdyn) }; }
    } catch { als = null; }
    // anchor load at maximum intact tension: tension at the anchor of the most loaded line at the maximum offset
    const eMax = run(() => lineForces(sys, [xMax * Math.cos(a), xMax * Math.sin(a), 0, 0, 0, 0]), 'Line state at the maximum offset');
    const worstLine = eMax.lines.reduce((p, l) => (l.T > p.T ? l : p), eMax.lines[0]);
    const anchorH = worstLine.anchorH + (Tmax - TqsMax), anchorV = worstLine.anchorV, grounded = worstLine.groundedLength;
    let anc, anchorSteel;
    if (inp.anchorType === 'drag') {
      const c = DRAG_ANCHOR_INDICATIVE[inp.anchorSoil];
      anc = { cap: dragAnchorUhc({ mass: inp.anchorMass * 1e3, a: c.a, b: c.b }).UHC, name: `${inp.anchorMass} t drag anchor` };
      anchorSteel = inp.anchorMass;
    } else {
      need(inp.su0 + inp.suGrad * inp.sucL > 1, 'The clay strength is too low for a suction anchor: enter the strength profile.');
      const s = run(() => caissonHorizontalCapacityClay({ D: inp.sucD, L: inp.sucL, su0: inp.su0 * 1e3, k: inp.suGrad * 1e3 }), 'Suction-anchor capacity');
      anc = { cap: s.Hmax, name: `${inp.sucD} m × ${inp.sucL} m suction anchor`, zOpt: s.zOptimal };
      anchorSteel = Math.PI * inp.sucD * inp.sucL * 0.035 * 7.85 * 1.25;   // 35 mm skirt with 25 % for lid and stiffening
    }
    const ad = anchorDesignCheck({ staticCapacity: anc.cap, tension: Math.max(anchorH, 1), safetyFactor: inp.anchorSF });
    const offLim = inp.offsetLimit / 100 * inp.depth;
    const checks = [
      check('Line tension, intact (ULS)', uls.Td / 1e6, uls.Sc / 1e6, 'MN', `DNV-OS-E301 Ch.2 Sec.2 (γmean = ${uls.gammaMean}, γdyn = ${uls.gammaDyn}, Sc = 0.95·MBL); DNV-ST-0119 §8`, 'Factored mean plus dynamic tension of the most loaded line against the characteristic strength of the corroded chain. Exceeding it means the line can part in the design storm.'),
      ...(als ? [check('Line tension with one line broken (ALS)', als.Td / 1e6, als.Sc / 1e6, 'MN', `DNV-OS-E301 Ch.2 Sec.2 (ALS, γmean = ${als.gammaMean}, γdyn = ${als.gammaDyn}); API RP 2SK §7 (transient)`, 'Tension in the remaining lines after the most loaded line fails, including the transient overshoot. If this fails, one line failure cascades into loss of station.')] : []),
      check('Maximum platform offset, intact', xMax, offLim, 'm', 'DNV-ST-0119 §8 (offset limits from the power cable); project criterion', 'Mean plus dynamic excursion against the offset the dynamic cable and neighbouring assets can accept. Beyond it the cable is over-tensioned or over-bent at the hang-off or touchdown.'),
      check('Anchor holding capacity', ad.uc * anc.cap / 1e6, anc.cap / 1e6, 'MN', 'API RP 2SK §7.4 / DNV-RP-E301 (drag anchors), DNV-RP-E303 (suction anchors)', 'Horizontal line load at the anchor times the required safety factor against the anchor holding capacity. A dragging anchor lets the platform walk off station and unloads the other lines unevenly.'),
      ...(inp.anchorType === 'drag' ? [check('Vertical load at the drag anchor (no uplift)', anchorV / 1e3 + 1e-9, Math.max(0.02 * anchorH / 1e3, 1), 'kN', 'API RP 2SK §7.4.2 (drag anchors: zero uplift in the intact condition)', 'Vertical component at the anchor at maximum tension against a nominal tolerance. Drag anchors lose capacity quickly with uplift: the line must stay on the seabed in front of the anchor.')] : []),
    ];
    const gov = worst(checks);
    // plan view and most-loaded line profile
    const nodes = [], lines = [], vals = [];
    eq.lines.forEach((l, i) => { const A = sys.lines[i].anchor, f = l.fairlead || sys.lines[i].fairlead; nodes.push([round(A[0], 2), round(A[1], 2), -inp.depth], [round(f[0] + eq.pose[0], 2), round(f[1] + eq.pose[1], 2), -inp.fairleadDepth]); lines.push([2 * i, 2 * i + 1]); vals.push(round(l.T / 1e3, 4), round(l.T / 1e3, 4)); });
    const aCat = worstLine.H / ch.w, xs = linspace(0, Math.acosh(1 + h / aCat) * aCat, 40);
    const chainT = ch.mass * inp.lineLength * inp.nLines / 1e3, steelT = chainT + anchorSteel * inp.nLines;
    const dT = MBL / 1e6 * (((inp.dChain + 5 - inp.corr) / (inp.dChain - inp.corr)) ** 2 - 1);
    return {
      kpis: [
        { key: 'pretension', label: 'Fairlead pretension (calm water)', value: round(pre / 1e3, 5), unit: 'kN', explain: 'Line tension at the fairlead with no environmental load: the initial condition of the mooring. Higher pretension gives a stiffer system and smaller offsets but higher mean line and anchor loads.' },
        { key: 'Hpre', label: 'Horizontal pretension component', value: round(L0.H / 1e3, 5), unit: 'kN', explain: 'Horizontal tension of each line in calm water; the lines balance each other so the platform sits at the centre.' },
        { key: 'catenaryCheck', label: 'Catenary identity (T − H) / (w·h)', value: round((pre - L0.H) / (ch.w * h), 5), unit: '-', explain: 'For a chain with a horizontal tangent on the seabed, fairlead tension minus horizontal tension equals submerged weight per metre times fairlead height. A value of 1 (slightly below for an elastic line) confirms the catenary solution.' },
        { key: 'Kxx', label: 'Surge stiffness of the mooring', value: round(Kxx / 1e3, 5), unit: 'kN/m', explain: 'Restoring force per metre of offset at the calm-water position. With mass and added mass it sets the surge natural period.' },
        { key: 'Tn', label: 'Surge natural period', value: round(Tn, 4), unit: 's', explain: 'Period of the platform on its mooring. It should be well above wave periods (typically 60–150 s for catenary systems) so that waves do not excite it directly.' },
        { key: 'offsetMean', label: 'Mean offset under thrust and environment', value: round(mean, 4), unit: 'm', explain: 'Static excursion at which the mooring restoring force balances the mean load.' },
        { key: 'restoringAtMean', label: 'Restoring force at the mean offset', value: round(rc.restoring[nO] / 1e3, 5), unit: 'kN', explain: 'Mooring force opposing the offset at the equilibrium position; it must equal the applied mean load (equilibrium check).' },
        { key: 'offsetMax', label: 'Maximum offset, intact', value: round(xMax, 4), unit: 'm', limit: round(offLim, 4), sense: 'max', explain: 'Mean offset plus the dynamic excursion: the design offset for the power cable.' },
        { key: 'Tmax', label: 'Maximum line tension, intact', value: round(Tmax / 1e3, 5), unit: 'kN', explain: 'Characteristic tension of the most loaded line at the maximum offset, with the dynamic part amplified.' },
        { key: 'ucUls', label: 'Line utilisation, intact (ULS)', value: round(uls.uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Factored tension over characteristic strength (0.95·MBL of the corroded chain).' },
        { key: 'ucAls', label: 'Line utilisation, one line broken (ALS)', value: round(als ? als.uc : 0, 4), unit: '-', limit: 1, sense: 'max', explain: 'Utilisation of the most loaded remaining line after a line failure (0 if the damaged system finds no equilibrium, which the notes then state).' },
        { key: 'offsetDamaged', label: 'Offset with one line broken (transient)', value: round(als ? als.offset : 0, 4), unit: 'm', explain: 'Largest excursion after the most loaded line fails, including overshoot. With three lines it can be hundreds of metres: the cable must then be designed to part or the system given redundancy.' },
        { key: 'anchorH', label: 'Horizontal load at the anchor', value: round(anchorH / 1e3, 5), unit: 'kN', explain: 'Line tension at the anchor of the most loaded line at maximum tension, reduced by seabed friction on the grounded chain.' },
        { key: 'anchorFs', label: 'Anchor factor of safety', value: round(Math.min(ad.fs, 999), 4), unit: '-', limit: inp.anchorSF, sense: 'min', explain: 'Anchor holding capacity over the horizontal anchor load.' },
        { key: 'grounded', label: 'Grounded chain length at maximum tension', value: round(grounded, 4), unit: 'm', explain: 'Chain still lying on the seabed in front of the anchor of the most loaded line. Zero means the anchor sees uplift.' },
        { key: 'steel', label: 'Chain and anchor steel', value: round(steelT, 5), unit: 't', explain: 'Mass of all mooring chain and anchors: the main quantity of the station-keeping system.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Restoring force and line tension versus offset', x: { label: 'Offset in the load direction', unit: 'm' }, y: { label: 'Force (kN)', unit: 'kN' },
          series: [{ name: 'Restoring force', x: offs.map((v) => round(v, 3)), y: rc.restoring.slice(0, nO).map((v) => round(v / 1e3, 4)) }, { name: 'Highest line tension', x: offs.map((v) => round(v, 3)), y: rc.maxTension.slice(0, nO).map((v) => round(v / 1e3, 4)) }, { name: 'Mean load', x: [0, offs[nO - 1]], y: [F / 1e3, F / 1e3], dash: true }, { name: 'Maximum offset', x: [xMax, xMax], y: [0, TqsMax / 1e3], dash: true }],
          explain: 'Force–offset characteristic of the mooring with the tension of the most loaded line. The system stiffens as the windward chain lifts off the seabed; the mean load line gives the mean offset where it crosses the restoring curve.' },
        { type: 'line', title: 'Most loaded line at maximum offset', x: { label: 'Horizontal distance from touchdown', unit: 'm' }, y: { label: 'Height above seabed', unit: 'm' }, series: [{ name: 'Suspended chain', x: xs.map((v) => round(v, 3)), y: xs.map((x) => round(aCat * (Math.cosh(x / aCat) - 1), 3)) }, { name: 'Grounded chain', x: [round(-grounded, 3), 0], y: [0, 0] }],
          explain: 'Catenary of the most loaded line from its touchdown point to the fairlead, with the chain still on the seabed towards the anchor. A short grounded length means the anchor is close to seeing uplift.' },
        { type: 'model3d', title: 'Mooring layout under the mean load', nodes, lines, tris: [], nodeValues: vals, valueLabel: 'Fairlead tension', valueUnit: 'kN',
          explain: 'Anchors and fairleads of the lines with the platform at its mean offset, coloured by line tension. The windward line carries the load; leeward lines slacken.' },
        { type: 'bar', title: 'Line tensions', categories: eq.lines.map((l, i) => `Line ${i + 1}`), series: [{ name: 'Calm water (kN)', y: e0.lines.map((l) => round(l.T / 1e3, 4)) }, { name: 'Mean load (kN)', y: eq.lines.map((l) => round(l.T / 1e3, 4)) }, { name: 'Maximum offset (kN)', y: eMax.lines.map((l) => round(l.T / 1e3, 4)) }], y: { label: 'Fairlead tension', unit: 'kN' },
          explain: 'Fairlead tension of every line in calm water, at the mean offset and at the maximum offset (quasi-static). The spread between lines shows how unevenly the system shares load in this direction.' },
      ],
      tables: [
        { title: 'Lines at the mean offset', columns: [{ key: 'i', label: 'Line' }, { key: 'T', label: 'Fairlead tension', unit: 'kN' }, { key: 'H', label: 'Horizontal', unit: 'kN' }, { key: 'V', label: 'Vertical', unit: 'kN' }, { key: 'aH', label: 'Anchor horizontal', unit: 'kN' }, { key: 'g', label: 'Grounded length', unit: 'm' }, { key: 'u', label: 'T / MBL', unit: '-' }],
          rows: eq.lines.map((l, i) => ({ i: i + 1, T: round(l.T / 1e3, 5), H: round(l.H / 1e3, 5), V: round(l.V / 1e3, 5), aH: round(l.anchorH / 1e3, 5), g: round(l.groundedLength, 4), u: round(l.T / MBL, 4) })),
          explain: 'Tension components, anchor load and grounded length of each line at the mean position. T/MBL above about 0.3 in the mean condition leaves little room for dynamics.' },
        { title: 'Mooring summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Chain', v: `${inp.dChain} mm ${cm.grade} studless`, u: '' }, { q: 'Mass in air / submerged weight', v: `${fmt(ch.mass)} / ${fmt(ch.w)}`, u: 'kg/m / N/m' }, { q: 'MBL new / corroded', v: `${fmt(ch.MBL / 1e6, 4)} / ${fmt(MBL / 1e6, 4)}`, u: 'MN' }, { q: 'Axial stiffness EA', v: fmt(ch.EA / 1e6, 4), u: 'MN' }, { q: 'Mean / dynamic / maximum tension', v: `${fmt(Tmean / 1e3, 4)} / ${fmt(Tdyn / 1e3, 4)} / ${fmt(Tmax / 1e3, 4)}`, u: 'kN' },
            { q: 'Anchor', v: anc.name, u: '' }, { q: 'Anchor capacity', v: fmt(anc.cap / 1e6, 4), u: 'MN' }, ...(anc.zOpt ? [{ q: 'Optimal padeye depth', v: fmt(anc.zOpt), u: 'm' }] : []), { q: 'Damaged condition', v: als ? `offset ${fmt(als.offset)} m, tension ${fmt(als.Tmax / 1e3)} kN` : 'no equilibrium found with the most loaded line removed', u: '' }],
          explain: 'Line make-up, strength and the tension components behind the checks, with the anchor selected.' },
      ],
      notes: [
        `${inp.nLines} × ${inp.lineLength} m of ${inp.dChain} mm ${cm.grade} chain to anchors at ${inp.anchorRadius} m in ${inp.depth} m of water; pretension ${fmt(pre / 1e3)} kN (${fmt(pre / MBL * 100)} % of corroded MBL).`,
        `Mean load ${fmt(F / 1e3)} kN at ${inp.dirDeg}° gives ${fmt(mean)} m mean offset; with ${inp.dynOffset} m of dynamic motion the maximum tension is ${fmt(Tmax / 1e3)} kN.`,
        als ? `One line broken: transient offset ${fmt(als.offset)} m, tension ${fmt(als.Tmax / 1e3)} kN.` : 'One line broken: the remaining lines cannot hold the platform against the mean load in this direction — the system has no redundancy (acceptable only with consequence class and cable design agreed for drift-off).',
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Intact line utilisation ${fmt(uls.uc)}, ${als ? `damaged ${fmt(als.uc)}` : 'no damaged equilibrium'}; maximum offset ${fmt(xMax)} m against ${fmt(offLim)} m; anchor factor of safety ${fmt(ad.fs)} against ${inp.anchorSF}.`,
          gov.uc > 1 ? (gov.name.startsWith('Line tension') ? `Increase the chain to about ${fmt(inp.dChain * Math.sqrt(gov.uc) + 2, 3)} mm (break load scales with d²) or reduce pretension by lengthening the lines.` : gov.name.startsWith('Maximum platform') ? `Stiffen the mooring: shorten the lines by about ${fmt(Math.max(1, (xMax - offLim) * 0.5), 2)} m or use heavier chain near the touchdown; otherwise raise the cable offset allowance.` : gov.name.startsWith('Anchor') ? `Increase the anchor capacity to ${fmt(anchorH * inp.anchorSF / 1e6)} MN (${inp.anchorType === 'drag' ? `about ${fmt(inp.anchorMass * (ad.uc) ** (1 / 0.92), 3)} t` : 'longer or larger caisson'}).` : 'Lengthen the lines so that chain remains grounded at the anchor at maximum tension.') : 'The mooring meets the tension, offset and anchor criteria in the quasi-static analysis: confirm with a coupled dynamic analysis (line dynamics, low-frequency motion, all directions) and add the fatigue check.',
          'DNV-ST-0119 (2021) §8 with DNV-OS-E301 partial factors; API RP 2SK quasi-static method for preliminary design.',
          `Mooring steel: ${fmt(steelT)} t. ${steelCo2Text(steelT)} Chain on the seabed sweeps an area of roughly ${fmt(inp.nLines * grounded * 2 * xMax / 1e4)} ha as the platform moves — taut or semi-taut synthetic lines cut both steel and footprint.`),
        adv('economic', 'info', `The system uses ${fmt(chainT)} t of chain and ${inp.nLines} anchors per platform; one 5 mm step up in chain diameter adds ${fmt(chainProperties({ d: d + 0.005, grade: cm.grade }).mass * inp.lineLength * inp.nLines / 1e3 - chainT)} t and about ${fmt(dT)} MN of break load.`,
          `Trade chain size against line length and anchor radius: a longer, lighter line lowers tension but enlarges the lease footprint. Anchor installation is about one anchor-handling vessel day per anchor (${inp.nLines} days per platform) plus hook-up.`,
          'Class 4 quantity-based estimate (AACE RP 18R-97): chain tonnage, anchors and anchor-handler days are the station-keeping cost drivers; whole-life costing for inspection and chain replacement per ISO 15686-5.',
          vesselCo2Text(inp.nLines + 1, 40)),
      ],
      publish: { 'mooring.Tmax': Tmax, 'mooring.Kxx': Kxx, 'mooring.anchorH': anchorH, 'mooring.anchorV': anchorV, 'mooring.pretension': pre, 'vessel.offsetMax': xMax, 'mooring.MBL': MBL, 'mooring.TmaxDamaged': als ? als.Tmax : Tmax, 'mooring.TnSurge': Tn },
      quantities: { steel_t: round(steelT, 4), vessel_days: { anchor_handling: inp.nLines + 1 } },
    };
  },
  verify: [
    { name: 'Hanging chain with a horizontal tangent on the seabed: T_fairlead − H = w·h', source: 'Catenary statics (Faltinsen 1990, Sea Loads on Ships and Offshore Structures, Eq. 8.20)',
      // 95 mm chain: w = 19.9e3·0.095²·g·(1 − 1025/7850) = 1531 N/m; w·h = 1531 × (200 − 14) = 284.8 kN; the elastic line gives about 0.2 % less than the inextensible identity
      inputs: {}, expect: [{ kpi: 'catenaryCheck', value: 1.0, tol: 0.005 }] },
    { name: 'Static equilibrium: restoring force at the mean offset equals the applied mean load', source: 'Force balance of the moored floater (API RP 2SK §5.3)',
      inputs: { thrust: 800, Fenv: 400, dirDeg: 0 }, expect: [{ kpi: 'restoringAtMean', value: 1200, tol: 0.005 }] },
  ],
};
