// Cable laying and installation: catenary lay from a chute, design (test) tension, chute sidewall pressure and bend
// radius, tensioner grip, slack control, vessel-motion dynamic tension and limiting sea state, landfall / J-tube
// pull-in, and the installation duration.
import { cableLay, layDynamicTension, limitingSeaState, shorePullIn } from '../../engine/lines/pipelay.js';
import { cableLayTension, cableHandlingLimits, telecomSinkAngle, telecomSlack } from '../../engine/xsection/cablemech.js';
import { lineInputs, lineFrom, polylineLength, adv, check, worst, sev, fmt, round, run, need, linspace, DEG, G, STD_CABLE, CONFIRM, MGO_T_PER_DAY, CO2_MGO, vesselCo2Text } from './_common.js';

const DEFAULT_PULL = [
  { type: 'straight', length: 60, w: 620, mu: 0.25, slope: 0, angle: 0, radius: 0 },
  { type: 'bend', length: 0, w: 0, mu: 0.25, slope: 0, angle: 30, radius: 6 },
  { type: 'straight', length: 400, w: 620, mu: 0.25, slope: 2, angle: 0, radius: 0 },
  { type: 'bend', length: 0, w: 0, mu: 0.25, slope: 0, angle: 15, radius: 8 },
  { type: 'straight', length: 150, w: 620, mu: 0.3, slope: 5, angle: 0, radius: 0 },
];

export default {
  id: 's13-laying-installation', suite: 13, title: 'Cable laying and installation',
  summary: 'Lay tension, departure angle, layback, chute and touchdown bend radii, sidewall pressure, tensioner grip, slack, weather limit, pull-in tension and installation duration.',
  covers: [6], fidelity: 'analytical',
  standards: [`${STD_CABLE} §4 (cable installation)`, 'CIGRE TB 623 (2015) §5 (tensile-bending test tension, sidewall pressure)', 'CIGRE TB 883 (2022) Installation of submarine power cables', 'DNV-ST-N001 (2021) (marine operations, weather-restricted operations)', 'Zajac (1957) Bell Syst. Tech. J. 36 (cable laying kinematics)'],
  models: ['Catenary laying (nonlinear catenary equations)', 'Bend radius and sidewall bearing pressure SWP = T/R', 'Capstan friction and cable-pulling equations', 'Tensioner friction grip', 'Vessel-motion dynamic cable tension (quasi-static plus inertia)', 'Steady-state sinking and slack control (Zajac)'],
  limits: 'Two-dimensional steady lay without current on the suspended span or bending stiffness (the touchdown module adds the boundary layer); the dynamic tension is the quasi-static plus inertia estimate, not a time-domain lay simulation. Handling limits are overridable defaults to be confirmed with the cable manufacturer; burial-tool interaction is in the burial module.',
  publishes: ['cable.layTension', 'cable.layDays', 'install.vesselDays'],
  searchTerms: 'submarine power cable laying installation tension chute tensioner',
  inputs: [
    ...lineInputs({ hydro: false }),
    { key: 'route', label: 'Imported route (optional)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Route', role: 'geometry', help: 'CSV / GeoJSON / DXF route polyline; its length replaces the route length below. Leave empty to use the form value.' },
    { key: 'L', label: 'Route length', unit: 'km', type: 'number', default: 25, min: 0.05, max: 800, group: 'Route', role: 'geometry', bind: 'cable.length' },
    { key: 'depth', label: 'Water depth at the lay position', unit: 'm', type: 'number', default: 40, min: 2, max: 3000, group: 'Route', role: 'environment', bind: 'site.depth' },
    { key: 'slope', label: 'Seabed slope in the lay direction', unit: 'deg', type: 'number', default: 0, min: -20, max: 20, group: 'Route', role: 'environment', help: 'Positive = seabed descending ahead of the vessel.' },
    { key: 'chuteH', label: 'Chute height above the water line', unit: 'm', type: 'number', default: 6, min: 0, max: 30, group: 'Lay vessel', role: 'geometry' },
    { key: 'chuteR', label: 'Chute / lay-wheel radius', unit: 'm', type: 'number', default: 4, min: 0.5, max: 15, group: 'Lay vessel', role: 'geometry', design: true },
    { key: 'H', label: 'Bottom (touchdown) tension', unit: 'kN', type: 'number', default: 12, min: 0.2, max: 2000, group: 'Lay vessel', role: 'bc', design: true, help: 'Prescribed lay tension at the seabed, controlled by the tensioner and layback.' },
    { key: 'laySpeed', label: 'Lay speed', unit: 'm/h', type: 'number', default: 400, min: 10, max: 5000, group: 'Lay vessel', role: 'operation', design: true },
    { key: 'slack', label: 'Target bottom slack', unit: '%', type: 'number', default: 1, min: -2, max: 10, group: 'Lay vessel', role: 'operation', help: 'Extra cable length per seabed length; 0–2 % for power cables before burial.' },
    { key: 'Cdn', label: 'Normal drag coefficient', unit: '-', type: 'number', default: 1.2, min: 0.3, max: 3, group: 'Lay vessel', role: 'load' },
    { key: 'nTracks', label: 'Tensioner tracks', type: 'int', default: 4, min: 2, max: 4, group: 'Tensioner', role: 'geometry' },
    { key: 'squeeze', label: 'Squeeze per track', unit: 'kN/m', type: 'number', default: 20, min: 0.5, max: 300, group: 'Tensioner', role: 'bc', design: true },
    { key: 'padMu', label: 'Pad friction coefficient', unit: '-', type: 'number', default: 0.15, min: 0.03, max: 0.6, group: 'Tensioner', role: 'material', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Wet polyurethane pad on a polyethylene serving: 0.1–0.2.' },
    { key: 'tensL', label: 'Tensioner contact length', unit: 'm', type: 'number', default: 6, min: 0.5, max: 30, group: 'Tensioner', role: 'geometry' },
    { key: 'Hs', label: 'Installation sea state Hs', unit: 'm', type: 'number', default: 2.0, min: 0, max: 8, group: 'Vessel motion', role: 'bc', bind: 'site.Hs1' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 8, min: 3, max: 20, group: 'Vessel motion', role: 'bc' },
    { key: 'rao', label: 'Chute heave RAO', unit: 'm/m', type: 'number', default: 0.8, min: 0, max: 3, group: 'Vessel motion', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Vertical motion of the chute per metre of wave amplitude (heave plus pitch × lever arm).' },
    { key: 'pull', label: 'Pull-in route (tail → winch)', type: 'table', group: 'Landfall / J-tube pull-in', role: 'bc',
      columns: [{ key: 'type', label: 'Type (straight / bend)', type: 'text' }, { key: 'length', label: 'Length', unit: 'm', type: 'number' }, { key: 'w', label: 'Weight in the local medium', unit: 'N/m', type: 'number' }, { key: 'mu', label: 'Friction', type: 'number' }, { key: 'slope', label: 'Slope up', unit: 'deg', type: 'number' }, { key: 'angle', label: 'Bend angle', unit: 'deg', type: 'number' }, { key: 'radius', label: 'Bend radius', unit: 'm', type: 'number' }],
      default: DEFAULT_PULL.map((r) => ({ ...r })), help: 'Sequence of straights and bends from the cable tail to the winch (beach pull, HDD duct or J-tube). An empty table skips the pull-in check.' },
    { key: 'tail', label: 'Back tension at the tail', unit: 'kN', type: 'number', default: 2, min: 0, max: 500, group: 'Landfall / J-tube pull-in', role: 'bc' },
    { key: 'Tallow', label: 'Maximum installation tension', unit: 'kN', type: 'number', default: 250, min: 2, max: 5000, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.allowTension' },
    { key: 'mbrFactor', label: 'Minimum bend radius under tension / diameter', unit: '-', type: 'number', default: 20, min: 5, max: 40, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'swpLimit', label: 'Allowable sidewall pressure', unit: 'kN/m', type: 'number', default: 50, min: 2, max: 500, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'squeezeLimit', label: 'Allowable squeeze per track', unit: 'kN/m', type: 'number', default: 40, min: 1, max: 500, group: 'Acceptance criteria', role: 'criteria', help: `From the crush test. ${CONFIRM}` },
    { key: 'dynFactor', label: 'Dynamic factor on w·d (test tension)', unit: '-', type: 'number', default: 1.3, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria', help: `CIGRE TB 623 tensile-bending test tension T = 1.3·w·d + H for depths to 500 m. ${CONFIRM}` },
    { key: 'gripSF', label: 'Safety factor on tensioner grip', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'operability', label: 'Weather operability', unit: '-', type: 'number', default: 0.75, min: 0.1, max: 1, group: 'Schedule', role: 'operation', bind: 'install.operabilityFrac', help: 'Fraction of time the limiting sea state is not exceeded.' },
    { key: 'fixedDays', label: 'Fixed operations (mobilisation, pull-ins, jointing, testing)', unit: 'days', type: 'number', default: 8, min: 0, max: 200, group: 'Schedule', role: 'operation' },
    { key: 'dayRate', label: 'Indicative lay-spread day rate', unit: 'USD/day', type: 'number', default: 180000, min: 0, max: 2e6, group: 'Schedule', role: 'cost', help: 'Indicative unit rate for the economic comparison only.' },
  ],
  compute(inp) {
    const line = lineFrom(inp), w = line.w, D = line.D;
    need(w > 0, `The cable is buoyant (apparent weight ${fmt(w)} N/m) and cannot be laid as a catenary.`);
    const L = (polylineLength(inp.route) ?? inp.L * 1e3), H = inp.H * 1e3, h = inp.depth + inp.chuteH;
    const mbr = inp.mbrFactor * D, V = inp.laySpeed / 3600;
    const sink = run(() => telecomSinkAngle({ weightWater: w, outerDiameter: D, shipSpeed: V, Cd: inp.Cdn, rho: line.rhoW, depth: inp.depth }), 'Sinking angle');
    const slk = run(() => telecomSlack({ sinkAngle: Math.min(sink.sinkAngle, 1.55), seabedSlope: inp.slope * DEG, shipSpeed: V, slack: inp.slack / 100 }), 'Slack control');
    const lay = run(() => cableLay({ w, h, H, chuteRadius: inp.chuteR, mbr, swpLimit: inp.swpLimit * 1e3, payoutSpeed: slk.payoutSpeed, vesselSpeed: V }), 'Catenary lay');
    // the span between chute and water line hangs with its weight in air: add the difference to the top tension
    const Ttop = lay.topTension + (line.m * G - w) * inp.chuteH;
    const test = cableLayTension({ submergedWeight: w, depth: inp.depth, bottomTension: H, dynamicFactor: inp.dynFactor });
    // vessel motion: most probable maximum heave amplitude in a 3 h sea state ≈ 1.86 × significant amplitude
    const heaveOf = (Hs) => 1.86 * inp.rao * Hs / 2;
    const dyn = (Hs) => layDynamicTension({ staticTension: Ttop, w, EA: line.EA, m: line.m, suspendedLength: lay.suspendedLength, topAngleDeg: lay.departureAngleDeg, heaveAmp: Math.max(heaveOf(Hs), 1e-6), period: inp.Tp });
    const d0 = dyn(inp.Hs), dT = inp.Hs > 0 ? d0.lower : 0;
    const Tmax = Ttop + dT, Tmin = H - dT;
    const hsGrid = linspace(0.25, 6, 24);
    // a static tension already above the allowable leaves no workable sea state
    const overStatic = Ttop >= inp.Tallow * 1e3;
    const lim = overStatic ? { allowableHs: 0 } : limitingSeaState({ staticTension: Ttop, maxTension: inp.Tallow * 1e3, minTension: 0, hsGrid, dafOfHs: (Hs) => dyn(Hs).dafLower });
    // limiting Hs also from the touchdown going slack (dynamic amplitude reaching the bottom tension)
    let hsSlack = hsGrid[hsGrid.length - 1];
    for (const Hs of hsGrid) if (dyn(Hs).lower >= H) { hsSlack = Hs; break; }
    const hsLimit = overStatic ? 0 : Math.min(Number.isFinite(lim.allowableHs) && lim.allowableHs > 0 ? lim.allowableHs : hsGrid[hsGrid.length - 1], hsSlack);
    const hand = cableHandlingLimits({ outerDiameter: D, installFactor: inp.mbrFactor, tension: Tmax, bendRadius: inp.chuteR, sidewallLimit: inp.swpLimit * 1e3, squeeze: inp.squeeze * 1e3, squeezeLimit: inp.squeezeLimit * 1e3, nTracks: inp.nTracks, mu: inp.padMu, safetyFactor: inp.gripSF });
    // pull-in
    const segs = (inp.pull || []).filter((r) => r && r.type);
    let pullT = null, pullRows = [];
    if (segs.length) {
      const es = segs.map((r, i) => {
        const type = String(r.type).trim().toLowerCase();
        need(type === 'straight' || type === 'bend', `Pull-in row ${i + 1}: type must be 'straight' or 'bend'.`);
        need(r.mu >= 0, `Pull-in row ${i + 1}: friction cannot be negative.`);
        if (type === 'straight') { need(r.length > 0 && r.w >= 0, `Pull-in row ${i + 1}: a straight needs a positive length and a non-negative weight.`); return { type, length: r.length, w: r.w, mu: r.mu, slopeDeg: r.slope || 0 }; }
        need(r.angle > 0 && r.radius > 0, `Pull-in row ${i + 1}: a bend needs a positive angle and radius.`);
        return { type, mu: r.mu, angleDeg: r.angle };
      });
      const p = run(() => shorePullIn({ tailTension: inp.tail * 1e3, segments: es }), 'Pull-in tension');
      pullT = p.winchTension;
      pullRows = segs.map((r, i) => ({ seg: i + 1, type: r.type, T: round(p.tensions[i] / 1e3), swp: String(r.type).trim().toLowerCase() === 'bend' ? round(p.tensions[i] / r.radius / 1e3) : null, mbrOk: String(r.type).trim().toLowerCase() === 'bend' ? (r.radius >= mbr ? 'yes' : 'no') : '' }));
    }
    const pullSwp = Math.max(0, ...pullRows.map((r) => r.swp ?? 0));
    const minPullR = Math.min(Infinity, ...segs.filter((r) => String(r.type).trim().toLowerCase() === 'bend').map((r) => r.radius));
    const checks = [
      check('Top tension incl. vessel motion', Tmax / 1e3, inp.Tallow, 'kN', `${STD_CABLE} §4; CIGRE TB 623 §5.2`, 'Static lay tension at the chute plus the dynamic amplitude in the installation sea state against the maximum installation tension. Overtension stretches conductors and armour.'),
      check('Chute bend radius', mbr, inp.chuteR, 'm', `${STD_CABLE} (MBR under tension); CIGRE TB 623`, 'Demand is the minimum bend radius under tension (multiple of the diameter), capacity the chute radius.'),
      check('Sidewall pressure on the chute', hand.sidewallPressure / 1e3, inp.swpLimit, 'kN/m', 'CIGRE TB 623 §5.4; SWP = T/R', 'Bearing load per metre where the tensioned cable wraps the chute. Excess pressure flattens the cores and indents the insulation.'),
      check('Bend radius at touchdown', mbr, lay.tdpBendRadius, 'm', `${STD_CABLE} §4 (lay-back and bottom tension control)`, 'Demand is the minimum bend radius, capacity the catenary radius H/w at the seabed. Too little bottom tension or a lay-back that is too short kinks the cable at touchdown.'),
      check('Touchdown tension stays positive', dT / 1e3, Math.max(H, 1) / 1e3, 'kN', 'DNV-RP-0360 §4 (no compression or loop formation at touchdown)', 'Dynamic tension amplitude against the bottom tension. When the amplitude exceeds the bottom tension the touchdown goes slack each cycle and the cable can form loops.'),
      check('Tensioner squeeze', inp.squeeze, inp.squeezeLimit, 'kN/m', 'CIGRE TB 623 §5.4 (crush test); manufacturer limit', 'Radial pad load per track per metre against the tested crush capacity of the cable.'),
      check('Tensioner grip length', hand.gripLength, inp.tensL, 'm', 'Friction grip μ·n·q·L ≥ SF·T', 'Contact length needed to hold the top tension by pad friction with the safety factor, against the tensioner length. Too short and the cable slips through the tensioner.'),
      ...(pullT != null ? [
        check('Pull-in winch tension', pullT / 1e3, inp.Tallow, 'kN', 'Cable-pulling equations (Rifenburg 1953); CIGRE TB 883', 'Tension at the pulling head after friction on straights and capstan build-up round bends, against the maximum installation tension.'),
        check('Pull-in sidewall pressure at bends', pullSwp, inp.swpLimit, 'kN/m', 'CIGRE TB 623 §5.4; SWP = T/R', 'Largest bearing load at a bend of the pull-in route (rollers, bellmouth, J-tube bend).'),
        ...(Number.isFinite(minPullR) ? [check('Pull-in bend radius', mbr, minPullR, 'm', `${STD_CABLE} (MBR under tension)`, 'Minimum bend radius under tension against the tightest bend of the pull-in route.')] : []),
      ] : []),
    ];
    const gov = worst(checks);
    const layDays = L / (inp.laySpeed * 24) / inp.operability, totalDays = layDays + inp.fixedDays;
    const Hgrid = linspace(Math.max(0.3 * H, w * mbr * 0.5), 3 * H, 30);
    const fuel = totalDays * MGO_T_PER_DAY;
    return {
      kpis: [
        { key: 'topTension', label: 'Static top tension at the chute', value: round(Ttop / 1e3, 5), unit: 'kN', explain: 'Bottom tension plus the weight of the suspended cable (submerged below the water line, in air above it).' },
        { key: 'Tmax', label: 'Top tension incl. vessel motion', value: round(Tmax / 1e3), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Static tension plus the lower-bound dynamic amplitude for the most probable maximum chute heave in the installation sea state.' },
        { key: 'testTension', label: 'Tensile-bending test tension', value: round(test.designTension / 1e3, 5), unit: 'kN', explain: 'Factor × w·d + H: the tension the cable should have been type-tested to for this laying depth.' },
        { key: 'departure', label: 'Departure angle from horizontal', value: round(lay.departureAngleDeg, 5), unit: 'deg', explain: 'Angle of the cable leaving the chute; the chute must cover it with margin for vessel pitch.' },
        { key: 'layback', label: 'Layback', value: round(lay.layBack, 5), unit: 'm', explain: 'Horizontal distance from the chute to the touchdown point — the quantity monitored by ROV or USBL to control bottom tension.' },
        { key: 'suspended', label: 'Suspended length', value: round(lay.suspendedLength), unit: 'm', explain: 'Cable length between chute and touchdown.' },
        { key: 'tdpRadius', label: 'Bend radius at touchdown', value: round(lay.tdpBendRadius, 5), unit: 'm', limit: round(mbr), sense: 'min', explain: 'Catenary radius H/w at the seabed; it must stay above the minimum bend radius.' },
        { key: 'swp', label: 'Sidewall pressure on the chute', value: round(hand.sidewallPressure / 1e3), unit: 'kN/m', limit: inp.swpLimit, sense: 'max', explain: 'Top tension divided by chute radius.' },
        { key: 'gripLength', label: 'Tensioner grip length required', value: round(hand.gripLength), unit: 'm', limit: inp.tensL, sense: 'max', explain: 'Pad contact length needed to hold the top tension with the grip safety factor.' },
        { key: 'sinkAngle', label: 'Steady sinking angle at lay speed', value: round(sink.sinkAngleDeg), unit: 'deg', explain: 'Straight-line angle at which weight balances normal drag at the lay speed (Zajac). Heavy power cables at low speed sink almost vertically, so the catenary governs.' },
        { key: 'payoutRatio', label: 'Pay-out / vessel speed', value: round(slk.speedRatio, 5), unit: '-', explain: 'Speed ratio that lays the target slack on the given seabed slope; below the zero-slack ratio the bottom tension rises.' },
        { key: 'hsLimit', label: 'Limiting sea state', value: round(hsLimit), unit: 'm', explain: 'Largest Hs for which the top tension stays below the limit and the touchdown stays in tension (lower-bound dynamics). It sets the weather criterion of the lay.' },
        ...(pullT != null ? [{ key: 'pullTension', label: 'Pull-in winch tension', value: round(pullT / 1e3), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Tension at the winch end of the pull-in route.' }] : []),
        { key: 'layDays', label: 'Lay duration incl. weather', value: round(layDays), unit: 'days', explain: 'Route length over lay speed, divided by the weather operability.' },
        { key: 'totalDays', label: 'Installation campaign', value: round(totalDays), unit: 'days', explain: 'Lay duration plus fixed operations; the vessel-day quantity for the estimate.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Lay catenary', x: { label: 'Distance from touchdown', unit: 'm' }, y: { label: 'Height above seabed', unit: 'm' },
          series: (() => { const a = H / w, xs = linspace(0, lay.layBack, 50); return [{ name: 'Cable', x: xs, y: xs.map((x) => a * (Math.cosh(x / a) - 1)) }, { name: 'Water line', x: [0, lay.layBack], y: [inp.depth, inp.depth], dash: true }]; })(),
          explain: 'Shape of the suspended cable from the touchdown (origin) to the chute. The slope at the top is the departure angle; the curvature at the origin is the touchdown bend.' },
        { type: 'line', title: 'Effect of bottom tension', x: { label: 'Bottom tension', unit: 'kN' }, y: { label: 'Value', unit: 'm or kN' },
          series: [{ name: 'Touchdown radius (m)', x: Hgrid.map((v) => v / 1e3), y: Hgrid.map((v) => v / w) }, { name: 'Top tension (kN)', x: Hgrid.map((v) => v / 1e3), y: Hgrid.map((v) => (v + w * h) / 1e3) }, { name: 'Minimum bend radius (m)', x: [Hgrid[0] / 1e3, Hgrid[29] / 1e3], y: [mbr, mbr], dash: true }],
          explain: 'Raising the bottom tension opens the touchdown radius linearly but adds the same amount to the top tension. The operating window is where the radius is above the dashed limit and the top tension below the allowable.' },
        { type: 'line', title: 'Dynamic tension against sea state', x: { label: 'Hs', unit: 'm' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Maximum top tension (lower bound)', x: hsGrid, y: hsGrid.map((v) => (Ttop + dyn(v).lower) / 1e3) }, { name: 'Minimum touchdown tension', x: hsGrid, y: hsGrid.map((v) => (H - dyn(v).lower) / 1e3) }, { name: 'Allowable', x: [hsGrid[0], hsGrid[23]], y: [inp.Tallow, inp.Tallow], dash: true }],
          explain: 'How the tension envelope opens with wave height. The limiting sea state is where the upper curve meets the allowable or the lower curve reaches zero.' },
      ],
      tables: [
        ...(pullRows.length ? [{ title: 'Pull-in tension build-up', explain: 'Tension after each segment from the tail to the winch, with sidewall pressure and bend-radius compliance at the bends.',
          columns: [{ key: 'seg', label: 'Segment' }, { key: 'type', label: 'Type' }, { key: 'T', label: 'Tension after', unit: 'kN' }, { key: 'swp', label: 'Sidewall pressure', unit: 'kN/m' }, { key: 'mbrOk', label: 'Radius ≥ MBR' }], rows: pullRows }] : []),
        { title: 'Installation compliance summary', explain: 'Each row is one handling limit with its demand, capacity and utilisation, for the installation compliance report.',
          columns: [{ key: 'name', label: 'Check' }, { key: 'demand', label: 'Demand' }, { key: 'capacity', label: 'Capacity' }, { key: 'unit', label: 'Unit' }, { key: 'uc', label: 'Utilisation' }], rows: checks.map((c) => ({ name: c.name, demand: c.demand, capacity: c.capacity, unit: c.unit, uc: round(c.uc, 3) })) },
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Lay at ${inp.laySpeed} m/h with ${fmt(lay.layBack)} m layback; dynamic tension amplitude ${fmt(dT / 1e3)} kN at Hs ${inp.Hs} m (quasi-static plus inertia; the drag-locked elastic bound is not credible for a slack lay catenary and is not reported).`,
        `Campaign: ${fmt(layDays)} days laying ${fmt(L / 1e3)} km at ${fmt(inp.operability * 100)} % operability + ${inp.fixedDays} fixed days.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Top tension ${fmt(Tmax / 1e3)} kN of ${inp.Tallow} kN, touchdown radius ${fmt(lay.tdpBendRadius)} m against ${fmt(mbr)} m, limiting Hs ${fmt(hsLimit)} m.`,
          gov.uc > 1 ? (gov.name.includes('touchdown') || gov.name.includes('Touchdown') ? `Raise the bottom tension to at least ${fmt(Math.max(w * mbr, dT) / 1e3 * 1.1)} kN or lower the weather limit.` : gov.name.includes('grip') ? `Use a longer tensioner (${fmt(hand.gripLength)} m needed), more tracks or a higher squeeze within the crush limit.` : gov.name.includes('Chute') || gov.name.includes('Sidewall') ? `Increase the chute radius to at least ${fmt(Math.max(mbr, hand.minRadiusForSidewall ?? 0))} m.` : 'Reduce tension (shorter pull lengths, rollers, lubricant) or re-route the pull-in.')
            : `Set the weather criterion at Hs ≤ ${fmt(hsLimit)} m (with the alpha factor of DNV-ST-N001 on the forecast) and monitor layback and touchdown by ROV; confirm with a time-domain lay analysis for the actual vessel RAOs.`,
          `${STD_CABLE} §4; CIGRE TB 883 (2022); DNV-ST-N001 (2021)`),
        adv('economic', 'info', `The installation needs about ${fmt(totalDays)} vessel days (${fmt(layDays)} laying, ${inp.fixedDays} fixed) — indicatively ${fmt(totalDays * inp.dayRate / 1e6)} M USD at the entered day rate of ${fmt(inp.dayRate, 6)} USD/day.`,
          `Weather operability is the lever: raising it from ${fmt(inp.operability * 100)} % to ${fmt(Math.min(100, inp.operability * 100 + 10))} % saves ${fmt(layDays - L / (inp.laySpeed * 24) / Math.min(1, inp.operability + 0.1))} days; a higher limiting sea state (larger chute, higher bottom tension) buys that.`,
          'AACE RP 18R-97 (Class 4 estimate from quantities); weather-window planning per DNV-ST-N001',
          `${vesselCo2Text(totalDays)} Fuel ≈ ${fmt(fuel)} t; each lay day saved avoids about ${fmt(MGO_T_PER_DAY * CO2_MGO)} t CO₂.`),
      ],
      publish: { 'cable.layTension': Tmax, 'cable.layDays': totalDays, 'install.vesselDays': totalDays },
      quantities: { cable_km: round(L / 1e3, 6), vessel_days: { cableLay: round(totalDays) } },
    };
  },
  verify: [
    // Catenary lay (Irvine 1981 §1.2): w = 415.73 N/m, d = 60 m, chute at the water line, H = 15 kN → a = H/w = 36.08 m;
    // T = H + w·d = 15000 + 24944 = 39.94 kN; layback = a·acosh(1 + d/a) = 36.08 × acosh(2.663) = 36.08 × 1.636 = 59.03 m;
    // departure angle = acos(H/T) = acos(0.3755) = 67.94°; CIGRE test tension 1.3·w·d + H = 47.43 kN.
    { name: 'Catenary lay: T = H + w·d, layback a·acosh(1 + d/a), departure angle acos(H/T), test tension 1.3·w·d + H',
      source: 'Irvine, Cable Structures (1981) §1.2; CIGRE TB 623 (2015) §5.2 test-tension formula; hand calculation',
      inputs: { depth: 60, chuteH: 0, H: 15, Hs: 0, pull: [] },
      expect: [{ kpi: 'topTension', value: 39.94, tol: 0.002 }, { kpi: 'layback', value: 59.03, tol: 0.003 }, { kpi: 'departure', value: 67.94, tol: 0.002 }, { kpi: 'tdpRadius', value: 36.08, tol: 0.002 }, { kpi: 'testTension', value: 47.43, tol: 0.002 }] },
  ],
};
