// Riser installation and deployment: J-lay catenary, staged tensions (empty, flooded, abandonment and recovery,
// transfer to the host), dynamic tension bounds from lay-vessel motion, limiting sea state and vessel days.
import { jLay, abandonmentRecovery, layDynamicTension, limitingSeaState, layCriteriaF101 } from '../../engine/lines/pipelay.js';
import { catenaryRiser } from '../../engine/lines/riser.js';
import { harmonic } from './_dyn.js';
import {
  G, DEG, need, run, fmt, round, thin, scaled, pipeInputs, waterInputs, hydroInputs, pipeFrom, steelTonnes, adv, check, kpi, plot, meshInput, finish,
  handWeight, steelCo2Text, LIMIT_2D, CONFIRM, MGO_T_PER_DAY, CO2_MGO,
} from './_common.js';

export default {
  id: 's09-installation', suite: 9, title: 'Riser installation and deployment',
  summary: 'J-lay deployment tensions and sag-bend loading, staged installation cases, abandonment and recovery, dynamic tension from lay-vessel motion, limiting sea state and installation vessel days.',
  covers: [11], fidelity: 'numerical',
  standards: ['DNV-ST-F101 (2021) Sec. 10 installation', 'DNV-ST-N001 (2021) marine operations', 'DNV-RP-N103 (2019) §5 lowering', 'API RP 2RD (1998, R2013) §7 installation'],
  models: ['Deployment catenary with bending stiffness (rod model) and seabed springs', 'Closed-form catenary for staged cases', 'Abandonment and recovery wire tension', 'Frequency-domain response of the lay catenary to tower-top motion (linearised drag)', 'Quasi-dynamic tension bounds (inertia lower bound, drag-locked elastic upper bound)', 'Simplified laying criteria (sag-bend equivalent stress)', 'Weibull sea-state operability'],
  limits: `${LIMIT_2D} Dynamic tension is a linearised single-harmonic response bracketed by two closed-form bounds, not a time-domain simulation; lay-tower roller geometry, reel-lay plastic strain, buoyancy-module installation sequence and pull-in winch dynamics are represented only by their static stages. Weather downtime uses a Weibull Hs distribution without persistence (no weather-window statistics).`,
  publishes: ['install.vesselDays', 'install.operabilityFrac', 'riser.layTension'],
  searchTerms: 'steel catenary riser J-lay installation tension abandonment recovery limiting sea state',
  inputs: [
    ...pipeInputs(),
    ...waterInputs(),
    ...hydroInputs(),
    { key: 'riserLength', label: 'Riser length to install', unit: 'm', type: 'number', default: 1620, min: 100, max: 20000, group: 'Installation sequence', role: 'geometry', bind: 'riser.length' },
    { key: 'towerAngle', label: 'J-lay tower angle from vertical', unit: 'deg', type: 'number', default: 8, min: 1, max: 30, group: 'Installation support', role: 'bc', design: true,
      help: 'Sets the lay catenary: a steeper tower gives lower tension and a tighter sag bend.' },
    { key: 'towerHeight', label: 'Tower exit height above still water', unit: 'm', type: 'number', default: 10, min: 0, max: 60, group: 'Installation support', role: 'bc' },
    { key: 'kn', label: 'Seabed vertical stiffness', unit: 'kPa', type: 'number', default: 1000, min: 5, max: 1e5, group: 'Installation support', role: 'soil' },
    { key: 'tensionerCap', label: 'Lay tensioner / hang-off clamp capacity', unit: 'kN', type: 'number', default: 2500, min: 100, max: 30000, group: 'Installation support', role: 'criteria' },
    { key: 'arCap', label: 'Abandonment and recovery winch capacity', unit: 'kN', type: 'number', default: 2500, min: 100, max: 30000, group: 'Installation support', role: 'criteria' },
    { key: 'wWire', label: 'A&R wire submerged weight', unit: 'N/m', type: 'number', default: 300, min: 20, max: 2000, group: 'Installation support', role: 'load' },
    { key: 'headWeight', label: 'Pull head / end fitting submerged weight', unit: 'kN', type: 'number', default: 30, min: 0, max: 1000, group: 'Installation support', role: 'load' },
    { key: 'hangoffAngle', label: 'In-place hang-off angle (transfer stage)', unit: 'deg', type: 'number', default: 12, min: 3, max: 30, group: 'Installation sequence', role: 'bc', bind: 'riser.hangoffAngleDeg' },
    { key: 'hangoffDepth', label: 'Host hang-off depth below still water', unit: 'm', type: 'number', default: 15, min: 0, max: 200, group: 'Installation sequence', role: 'bc' },
    { key: 'heaveRao', label: 'Tower-top heave per metre of Hs', unit: 'm/m', type: 'number', default: 0.25, min: 0, max: 1.5, group: 'Lay vessel motion', role: 'bc', uncertain: { dist: 'normal', cov: 0.2 },
      help: 'Single-amplitude vertical motion of the tower top divided by significant wave height (including pitch lever arm), from the lay-vessel motion analysis.' },
    { key: 'surgeRao', label: 'Tower-top surge per metre of Hs', unit: 'm/m', type: 'number', default: 0.1, min: 0, max: 1.5, group: 'Lay vessel motion', role: 'bc' },
    { key: 'Tz', label: 'Wave period for dynamic tension', unit: 's', type: 'number', default: 8, min: 3, max: 20, group: 'Lay vessel motion', role: 'bc' },
    { key: 'HsOp', label: 'Planned limiting sea state', unit: 'm', type: 'number', default: 2.5, min: 0.5, max: 8, group: 'Lay vessel motion', role: 'operation', bind: 'site.Hs1' },
    { key: 'weibA', label: 'Hs Weibull scale (installation season)', unit: 'm', type: 'number', default: 1.8, min: 0.3, max: 6, group: 'Weather and schedule', role: 'environment' },
    { key: 'weibK', label: 'Hs Weibull shape', unit: '-', type: 'number', default: 1.5, min: 0.7, max: 3, group: 'Weather and schedule', role: 'environment' },
    { key: 'layRate', label: 'Lay rate', unit: 'm/day', type: 'number', default: 1500, min: 100, max: 8000, group: 'Weather and schedule', role: 'operation', help: 'J-lay of welded double or quad joints: roughly 1–3 km/day.' },
    { key: 'fixedDays', label: 'Fixed operations (initiation, A&R, transfer, hand-over)', unit: 'days', type: 'number', default: 5, min: 0, max: 60, group: 'Weather and schedule', role: 'operation' },
    { key: 'dayRate', label: 'Indicative vessel day rate', unit: 'USD', type: 'number', default: 350000, min: 0, max: 2e6, group: 'Weather and schedule', role: 'cost', help: 'Indicative unit rate for the lay spread; replace with the contract rate.' },
    { key: 'sagFactor', label: 'Sag-bend equivalent stress factor', unit: '-', type: 'number', default: 0.87, min: 0.5, max: 1, group: 'Acceptance criteria', role: 'criteria', help: `Allowable equivalent stress / yield in the sag bend (simplified laying criterion). ${CONFIRM}` },
    meshInput(160),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [80, 160, 320, 640], qoi: 'layTension' },

  compute(inp, ctx) {
    // engineering units → SI: mm → m, kN → N, kPa → Pa
    const empty = pipeFrom(inp, ctx, { rhoInt: 0 }), flooded = pipeFrom(inp, ctx, { rhoInt: empty.rhoW }), oper = pipeFrom(inp, ctx);
    const w = empty.section.w;
    need(w > 0, `The empty pipe is buoyant (apparent weight ${Math.round(w)} N/m): it cannot be laid empty without added weight. Lay flooded or increase the wall.`);
    const h = inp.depth + inp.towerHeight;
    const lay = run(() => jLay({ pipe: empty.section, D: empty.D, depth: h, towerAngleDeg: inp.towerAngle, seabed: { kn: inp.kn * 1e3 }, nElem: inp.nElem }), 'The J-lay configuration did not solve');
    const rod = lay.rod;
    need(rod.converged, 'The J-lay solution did not converge; change the tower angle.');
    const pe = empty.rhoW * G * inp.depth;
    const crit = layCriteriaF101({ overbendStrain: 0, strainLimit: 1, D: empty.D, t: empty.t, fy: empty.mat.fy, sagbendTe: lay.bottomTension, sagbendMoment: lay.sagbendMoment, pe, stressFactor: inp.sagFactor });
    // dynamic bounds at the planned limiting sea state and the allowable sea state
    const alpha = 90 - inp.towerAngle, mAx = empty.section.m;
    const dyn = (Hs) => layDynamicTension({ staticTension: lay.topTension, w, EA: empty.section.EA, m: mAx, suspendedLength: lay.suspendedLength, topAngleDeg: alpha, heaveAmp: inp.heaveRao * Hs, surgeAmp: inp.surgeRao * Hs, period: inp.Tz });
    const d0 = dyn(inp.HsOp);
    const cap = inp.tensionerCap * 1e3;
    need(cap > lay.topTension, `Static lay tension ${fmt(lay.topTension / 1e3)} kN already exceeds the tensioner capacity ${inp.tensionerCap} kN. Use a larger spread, a steeper tower or lay with less weight.`);
    // resolved dynamic tension: harmonic response of the lay catenary to tower-top motion (drag linearised per Hs)
    const om = 2 * Math.PI / inp.Tz, nr = rod.s.length;
    const resolved = (Hs) => {
      const r = harmonic(rod, { omega: om, xAmp: inp.surgeRao * Hs, zAmp: inp.heaveRao * Hs, rho: empty.rhoW, zeta: 0.005 });
      let minTe = Infinity; for (let k = 0; k < nr; k++) minTe = Math.min(minTe, rod.Te[k] - r.tensionAmp[k]);
      return { amp: r.tensionAmp[0], minTe };
    };
    const hsGrid = [0.5, 1.5, 2.5, 3.5, 4.5, 6], rs = hsGrid.map(resolved), rOp = hsGrid.includes(inp.HsOp) ? rs[hsGrid.indexOf(inp.HsOp)] : resolved(inp.HsOp);
    // allowable Hs: first crossing of maximum tension = capacity or minimum effective tension = 0
    const margin = (i) => Math.min(cap - (lay.topTension + rs[i].amp), rs[i].minTe);
    let hsAllow = margin(0) <= 0 ? 0 : hsGrid[hsGrid.length - 1];
    for (let i = 1; i < hsGrid.length && margin(0) > 0; i++) if (margin(i) <= 0) { hsAllow = hsGrid[i - 1] + margin(i - 1) * (hsGrid[i] - hsGrid[i - 1]) / (margin(i - 1) - margin(i)); break; }
    const lim = limitingSeaState({ staticTension: lay.topTension, maxTension: cap, minTension: 0 });
    const operab = 1 - Math.exp(-((Math.max(hsAllow, 1e-6) / inp.weibA) ** inp.weibK));
    // staged static cases from the closed-form catenary
    const cE = catenaryRiser({ h, w, hangoffAngleDeg: inp.towerAngle });
    const stage = (name, ww, hh, ang, note) => { const c = catenaryRiser({ h: hh, w: ww, hangoffAngleDeg: ang }); return { name, T: c.topTension, H: c.H, S: c.suspendedLength, R: c.a, note }; };
    const stages = [
      stage('Normal lay, empty pipe', w, h, inp.towerAngle, 'J-lay tower'),
      stage('Lay with flooded pipe (wet-buckle contingency)', flooded.section.w, h, inp.towerAngle, 'J-lay tower'),
      stage('Transfer to host hang-off, empty', w, inp.depth - inp.hangoffDepth, inp.hangoffAngle, 'pull-in to porch'),
      stage('In place, operating contents', oper.section.w, inp.depth - inp.hangoffDepth, inp.hangoffAngle, 'host porch'),
      stage('In place, flooded (hydrotest)', flooded.section.w, inp.depth - inp.hangoffDepth, inp.hangoffAngle, 'host porch'),
    ];
    const heads = []; for (let i = 0; i <= 10; i++) heads.push(inp.depth * i / 10);
    const ar = abandonmentRecovery({ H: cE.H, wPipe: w, wWire: inp.wWire, depth: h, headHeights: heads, headWeight: inp.headWeight * 1e3 });
    const arMax = Math.max(...ar.map((r) => r.vesselTension));
    const arFl = abandonmentRecovery({ H: stages[1].H, wPipe: flooded.section.w, wWire: inp.wWire, depth: h, headHeights: heads, headWeight: inp.headWeight * 1e3 });
    const arFlMax = Math.max(...arFl.map((r) => r.vesselTension));
    const need_ = inp.riserLength / inp.layRate + inp.fixedDays, days = need_ / Math.max(operab, 0.05);
    const steel = steelTonnes(empty.D, empty.t, inp.riserLength, empty.mat.rho);
    const Tdyn = lay.topTension + rOp.amp;

    const kpis = [
      kpi('layTension', 'Static lay tension at the tower', lay.topTension / 1e3, 'kN', 'Tension the lay tower holds during normal lay with the pipe empty. It must stay inside the tensioner or clamp capacity with margin for dynamics.'),
      kpi('layTensionDyn', 'Maximum lay tension at the planned sea state', Tdyn / 1e3, 'kN', 'Static tension plus the dynamic amplitude from the frequency-domain response of the lay catenary to tower-top motion at the planned limiting Hs (drag linearised).', inp.tensionerCap, 'max'),
      kpi('dafResolved', 'Dynamic amplification (resolved response)', Tdyn / lay.topTension, '-', `Maximum over static top tension at the planned sea state from the frequency-domain solution; the limiting value set by the tensioner capacity is ${fmt(lim.limitingDaf)}.`),
      kpi('dafUpper', 'Dynamic amplification, upper bound', d0.dafUpper, '-', 'Upper bound on (static + dynamic)/static top tension at the planned sea state: the axial motion of the tower top is taken by elastic stretch of the suspended pipe.'),
      kpi('dafLower', 'Dynamic amplification, lower bound', d0.dafLower, '-', 'Lower bound: the catenary follows the vessel quasi-statically and only weight change and axial inertia add tension. The true value lies between the two bounds.'),
      kpi('minTensionDyn', 'Minimum effective tension along the pipe (planned sea state)', rOp.minTe / 1e3, 'kN', 'Lowest static minus dynamic effective tension anywhere in the suspended pipe at the planned sea state. Negative means compression in the sag bend during lay — a buckling risk.', 0, 'min'),
      kpi('bottomTension', 'Bottom (touchdown) tension during lay', lay.bottomTension / 1e3, 'kN', 'Horizontal tension left in the pipe at touchdown. It becomes residual lay tension in the flowline and controls the sag-bend radius.'),
      kpi('sagStrain', 'Sag-bend bending strain', lay.sagbendStrain * 100, '%', 'Outer-fibre bending strain at the tightest point of the lay catenary; installation normally keeps it elastic.'),
      kpi('sagUc', 'Sag-bend equivalent stress utilisation', crit.ucSagbend, '-', 'Equivalent stress from axial force, bending and external pressure in the empty pipe over the allowable fraction of yield. Above 1.0 the sag bend is over-stressed.', 1, 'max'),
      kpi('layBack', 'Lay-back distance', lay.layBack, 'm', 'Horizontal distance from the tower to the touchdown point; it positions the vessel relative to the target and the host.'),
      kpi('suspended', 'Suspended length during lay', lay.suspendedLength, 'm', 'Length of pipe hanging between the tower and the seabed during normal lay.'),
      kpi('arTension', 'Maximum A&R wire tension (empty pipe)', arMax / 1e3, 'kN', 'Largest vessel-side wire tension while the pipe head is lowered to or recovered from the seabed with the pipe empty.', inp.arCap, 'max'),
      kpi('arTensionFlooded', 'Maximum A&R wire tension (flooded pipe)', arFlMax / 1e3, 'kN', 'Same with the pipe flooded: the contingency case that usually sizes the A&R system.'),
      kpi('transferTension', 'Tension at transfer to the host', stages[2].T / 1e3, 'kN', 'Top tension when the empty riser is handed over to the host porch at the in-place hang-off angle; sizes the pull-in winch and rigging.'),
      kpi('hsAllow', 'Allowable significant wave height', hsAllow, 'm', 'Largest Hs at which the resolved dynamic tension stays below the tensioner capacity and the pipe stays in tension everywhere. It is the operational limit before any alpha factor for forecast uncertainty.'),
      kpi('operability', 'Weather operability', operab * 100, '%', 'Fraction of time Hs is below the allowable limit for the installation season (Weibull distribution, no persistence).'),
      kpi('vesselDays', 'Installation vessel days', days, 'days', 'Lay time plus fixed operations divided by operability: the expected duration of the lay spread on site.'),
    ];
    const checks = [
      check('Lay tension vs tensioner capacity (planned sea state)', Tdyn / 1e3, inp.tensionerCap, 'kN', 'DNV-ST-F101 (2021) Sec. 10; DNV-ST-N001 (2021) equipment capacity', 'Guards against losing the pipe: if the tension demand exceeds the tensioner or clamp capacity the pipe can slip through the tower.'),
      check('Sag-bend equivalent stress', crit.sagbendStress / 1e6, crit.sagbendAllowable / 1e6, 'MPa', 'DNV-ST-F101 (2021) Sec. 10 simplified laying criteria', 'Guards against yielding or buckling of the empty pipe in the sag bend under bending and external pressure; a buckle here floods the pipe (wet buckle).'),
      check('A&R winch capacity (flooded contingency)', arFlMax / 1e3, inp.arCap, 'kN', 'DNV-ST-N001 (2021) contingency and recovery capability', 'Checks that a flooded pipe can still be held and recovered after a wet buckle. If not, the pipe would have to be laid down and recovered by another method.'),
      check('Planned sea state vs allowable', inp.HsOp, Math.max(hsAllow, 1e-6), 'm', 'DNV-ST-N001 (2021) operational limiting criteria', 'Compares the sea state the operation is planned for with the limit set by dynamic tension. Planning above the allowable sea state leaves no margin for motions.'),
    ];

    const [x, z] = thin([rod.x, rod.z], 250), [s1, te, mo] = thin([rod.s, rod.Te, rod.moment], 300);
    const plots = [
      plot('Deployment configuration', 'Horizontal distance from the tower', 'm', 'Height above seabed', 'm', [{ name: 'Lay catenary (empty)', x, y: z }, { name: 'Sea surface', x: [0, x[x.length - 1]], y: [inp.depth, inp.depth], dash: true }], 'Shape of the pipe from the J-lay tower to the seabed during normal lay. The sag bend just above touchdown is where bending is highest.'),
      plot('Tension and moment during lay', 'Arc length from the tower', 'm', 'Tension | moment', 'kN | kN·m', [{ name: 'Effective tension (kN)', x: s1, y: scaled(te, 1e-3) }, { name: 'Bending moment (kN·m)', x: s1, y: scaled(mo, 1e-3) }], 'Effective tension falls from the tower to the touchdown; the bending moment peaks in the sag bend. Both are needed for the sag-bend stress check.'),
      plot('A&R wire tension', 'Pipe head height above seabed', 'm', 'Vessel wire tension', 'kN', [{ name: 'Empty pipe', x: heads, y: ar.map((r) => r.vesselTension / 1e3) }, { name: 'Flooded pipe', x: heads, y: arFl.map((r) => r.vesselTension / 1e3) }, { name: 'Winch capacity', x: [0, inp.depth], y: [inp.arCap, inp.arCap], dash: true }], 'Wire tension at the vessel as the pipe head is lowered or raised. The flooded curve against the capacity line shows whether recovery after a wet buckle is possible.'),
      plot('Dynamic tension bounds against sea state', 'Significant wave height', 'm', 'Maximum top tension', 'kN', [{ name: 'Resolved (frequency domain)', x: hsGrid, y: rs.map((q) => (lay.topTension + q.amp) / 1e3) }, { name: 'Elastic upper bound', x: hsGrid, y: hsGrid.map((v) => (lay.topTension + dyn(v).upper) / 1e3), dash: true }, { name: 'Quasi-static lower bound', x: hsGrid, y: hsGrid.map((v) => (lay.topTension + dyn(v).lower) / 1e3), dash: true }, { name: 'Tensioner capacity', x: [hsGrid[0], hsGrid[hsGrid.length - 1]], y: [inp.tensionerCap, inp.tensionerCap], dash: true }], 'Maximum lay tension against Hs: the resolved response between its closed-form bounds. The allowable sea state is where the resolved curve meets the capacity line (or where minimum tension reaches zero).'),
    ];
    plots.push({ type: 'bar', title: 'Top tension by installation stage', categories: stages.map((s) => s.name), x: { label: 'Stage', unit: '' }, y: { label: 'Top tension', unit: 'kN' }, series: [{ name: 'Static top tension', y: stages.map((s) => s.T / 1e3) }],
      explain: 'Static top tension at each stage of the installation sequence. The flooded cases show the contingency loads the equipment must survive.' });
    const tables = [
      { title: 'Installation stages', columns: [{ key: 'stage', label: 'Stage' }, { key: 'at', label: 'Held by' }, { key: 'T', label: 'Top tension', unit: 'kN' }, { key: 'H', label: 'Bottom tension', unit: 'kN' }, { key: 'S', label: 'Suspended length', unit: 'm' }, { key: 'R', label: 'Sag-bend radius', unit: 'm' }],
        rows: stages.map((s) => ({ stage: s.name, at: s.note, T: round(s.T / 1e3, 5), H: round(s.H / 1e3, 4), S: round(s.S, 5), R: round(s.R, 4) })),
        explain: 'Critical stages of the deployment sequence with the static tension at the support, the residual bottom tension, the suspended length and the catenary sag-bend radius H/w.' },
      { title: 'Abandonment and recovery', columns: [{ key: 'z', label: 'Head height', unit: 'm' }, { key: 'Th', label: 'Head tension', unit: 'kN' }, { key: 'Tv', label: 'Vessel tension (empty)', unit: 'kN' }, { key: 'Tf', label: 'Vessel tension (flooded)', unit: 'kN' }, { key: 'L', label: 'Pipe suspended', unit: 'm' }],
        rows: ar.map((r, i) => ({ z: round(r.headHeight, 5), Th: round(r.headTension / 1e3, 4), Tv: round(r.vesselTension / 1e3, 5), Tf: round(arFl[i].vesselTension / 1e3, 5), L: round(r.pipeSuspendedLength, 5) })),
        explain: 'Wire and head tensions as the pipe head is lowered to the seabed (abandonment) or recovered; use the maximum for winch and wire selection.' },
    ];
    const ucMax = Math.max(...checks.map((c) => c.uc)), gov = checks.find((c) => c.uc === ucMax);
    const fuel = days * MGO_T_PER_DAY;
    const wd = days - need_;
    const advice = [
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.85 ? 'warn' : 'ok',
        `${gov.name} governs at ${fmt(ucMax)}. Static lay tension ${fmt(lay.topTension / 1e3)} kN, up to ${fmt(Tdyn / 1e3)} kN at Hs = ${inp.HsOp} m; allowable sea state ${fmt(hsAllow)} m; flooded A&R tension ${fmt(arFlMax / 1e3)} kN.`,
        ucMax > 1 ? (gov.name.startsWith('Sag') ? `Reduce sag-bend loading by laying at a larger tower angle than ${inp.towerAngle}° (higher bottom tension) or with a thicker wall.`
          : gov.name.startsWith('A&R') ? `Provide an A&R system rated for at least ${fmt(arFlMax / 1e3)} kN, or define a lay-down procedure that avoids lifting a flooded string.`
            : gov.name.startsWith('Planned') ? `Plan the operation for Hs ≤ ${fmt(hsAllow)} m, or reduce tower-top motion (heading control, larger vessel).`
              : `Use a spread with at least ${fmt(Tdyn / 1e3 * 1.1)} kN tension capacity or reduce the tower angle.`)
          : 'Confirm the dynamic amplification with a time-domain installation analysis for the selected vessel, and apply the forecast alpha factor to the allowable sea state before fixing the operational limit.',
        `DNV-ST-F101 (2021) Sec. 10 and DNV-ST-N001 (2021) operational limits. ${CONFIRM}`),
      adv('economic', operab < 0.6 ? 'warn' : 'info',
        `Installing ${fmt(inp.riserLength)} m takes ${fmt(need_)} working days and ${fmt(days)} vessel days at ${fmt(operab * 100)} % operability, of which ${fmt(wd)} days are waiting on weather. At the indicative rate of ${fmt(inp.dayRate)} USD/day that is about ${fmt(days * inp.dayRate / 1e6)} million USD (indicative unit rate × days, not a quote).`,
        `Raising the allowable Hs is worth more than raising the lay rate here: each 0.25 m of extra limit changes operability by about ${fmt(100 * (Math.exp(-((Math.max(hsAllow, 0.01) / inp.weibA) ** inp.weibK)) - Math.exp(-(((hsAllow + 0.25) / inp.weibA) ** inp.weibK))))} percentage points. Check a higher-capacity tensioner or a better-heading vessel against the day-rate difference.`,
        'AACE RP 18R-97 class 4 estimate from vessel days × indicative day rate; weather risk allowance per DNV-ST-N001 planning of weather-restricted operations.',
        `Lay spread fuel ≈ ${fmt(fuel)} t marine gas oil ≈ ${fmt(fuel * CO2_MGO)} t CO₂ (${MGO_T_PER_DAY} t/day and ${CO2_MGO} t CO₂ per t fuel, screening factors); waiting on weather accounts for ${fmt(wd * MGO_T_PER_DAY * CO2_MGO)} t of it. Pipe steel installed: ${steelCo2Text(steel)}`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`Empty apparent weight ${fmt(w)} N/m, flooded ${fmt(flooded.section.w)} N/m, operating ${fmt(oper.section.w)} N/m; tower exit ${fmt(h)} m above the seabed.`, `${gov.name} governs (${fmt(ucMax)}). Dynamic tension is for a single harmonic of period ${inp.Tz} s; closed-form bounds DAF ${fmt(d0.dafLower)}–${fmt(d0.dafUpper)}.`],
      advice,
      publish: { 'install.vesselDays': days, 'install.operabilityFrac': operab, 'riser.layTension': lay.topTension },
      quantities: { steel_t: steel, pipe_km: inp.riserLength / 1e3, vessel_days: { riser_lay: days } },
    });
  },

  verify: [
    {
      // J-lay catenary with the tower aligned to the natural top angle: T = w·h/(1 − sinθ),
      // h = 1000 + 10 = 1010 m, θ = 8°, w = empty apparent weight from the hand formula (≈ 666 N/m) → T ≈ 781 kN.
      name: 'J-lay top tension T = w·h/(1 − sin θ)', source: 'Irvine, Cable Structures (1981) §1.2; DNV-ST-F101 Sec. 10 catenary lay',
      inputs: { nElem: 320 },
      expect: [{ kpi: 'layTension', value: handWeight({ rhoInt: 0 }) * 1010 / (1 - Math.sin(8 * DEG)) / 1e3, tol: 0.005 },
        // transfer stage: empty riser on the host porch, h = 985 m, θ = 12° (closed form)
        { kpi: 'transferTension', value: handWeight({ rhoInt: 0 }) * 985 / (1 - Math.sin(12 * DEG)) / 1e3, tol: 0.002 }],
    },
  ],
};
