// Installation and lay configuration: catenary lay over a chute, overboarding curvature and sidewall pressure,
// tensioner grip, dynamic tension in the lay sea state, the allowable deployment envelope (bottom tension / lay-back
// window and limiting sea state), landing of the termination head and the lay duration.
import { cableLay, layDynamicTension, limitingSeaState } from '../../engine/lines/pipelay.js';
import { cableHandlingLimits } from '../../engine/xsection/cablemech.js';
import { terminationLandingLoad } from '../../engine/xsection/ancillaries.js';
import { heaveResponse } from '../../engine/lines/lowering.js';
import { lineInputs, lineFrom, adv, check, worst, sev, fmt, round, run, need, linspace, G, CO2_MGO, MGO_T_PER_DAY, STD_UMB, vesselCo2Text } from './_common.js';

export default {
  id: 's12-installation-lay', suite: 12, title: 'Installation and lay configuration',
  summary: 'Lay tension, overboarding and touchdown curvature, tensioner grip, the allowable lay window and sea state, termination landing and lay duration.',
  covers: [8], fidelity: 'analytical',
  standards: [`${STD_UMB} §15 (load-out and installation)`, 'DNV-RP-0360 (2016) §4 (cable installation: MBR, sidewall pressure)', 'DNV-RP-N103 (2021) §4.3 and §5 (wave amplitude 0.9·Hs, lowering)', 'DNV-ST-N001 (2021) (weather-restricted operations)'],
  models: ['Catenary lay geometry (Irvine 1981)', 'Minimum bend-radius relationships at chute and touchdown', 'Sidewall pressure T/R and capstan / tensioner grip', 'Quasi-static plus inertia dynamic lay tension', 'Lowering-line landing impact (energy method)'],
  limits: 'Catenary without bending stiffness or current, flat seabed, single-harmonic vessel motion: the dynamic tension is the quasi-static-plus-inertia estimate, with the drag-locked elastic value reported as an upper bound only. A time-domain lay analysis with vessel RAOs is needed to set final weather limits. Handling limits are overridable defaults.',
  publishes: ['umbilical.layTension', 'umbilical.layHs', 'install.vesselDays'],
  searchTerms: 'umbilical installation lay analysis chute tensioner weather window',
  inputs: [
    ...lineInputs().filter((i) => !['Cdn', 'Cdt', 'Ca', 'EI'].includes(i.key)),
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 350, min: 5, max: 3500, group: 'Lay configuration', role: 'environment', bind: 'site.depth' },
    { key: 'chuteH', label: 'Chute departure height above still water', unit: 'm', type: 'number', default: 8, min: 0, max: 60, group: 'Lay configuration', role: 'bc', help: 'Installation vessel geometry: elevation of the chute / lay wheel exit.' },
    { key: 'H', label: 'Bottom (touchdown) tension', unit: 'kN', type: 'number', default: 20, min: 0.05, max: 2000, group: 'Lay configuration', role: 'ic', design: true, uncertain: { dist: 'lognormal', cov: 0.2 },
      help: 'The lay-control variable: set through tensioner tension or lay-back. Starting stage of the deployment.' },
    { key: 'chuteR', label: 'Chute / lay-wheel radius', unit: 'm', type: 'number', default: 9, min: 0.5, max: 30, group: 'Lay configuration', role: 'bc', design: true },
    { key: 'vVessel', label: 'Vessel lay speed', unit: 'm/h', type: 'number', default: 400, min: 10, max: 3000, group: 'Lay configuration', role: 'operation', design: true },
    { key: 'slackPct', label: 'Pay-out slack', unit: '%', type: 'number', default: 1, min: -2, max: 10, group: 'Lay configuration', role: 'operation', help: 'Product length paid out over vessel distance, minus one.' },
    { key: 'nTracks', label: 'Tensioner tracks', type: 'select', default: 4, options: [{ value: 2, label: '2-track' }, { value: 3, label: '3-track' }, { value: 4, label: '4-track' }], group: 'Tensioner', role: 'operation' },
    { key: 'squeeze', label: 'Squeeze load per track', unit: 'kN/m', type: 'number', default: 40, min: 1, max: 500, group: 'Tensioner', role: 'bc', design: true },
    { key: 'squeezeLimit', label: 'Allowable squeeze per track', unit: 'kN/m', type: 'number', default: 60, min: 1, max: 1000, group: 'Tensioner', role: 'criteria', bind: 'umbilical.allowSqueeze', help: 'From the crush analysis or crush test of the cross-section.' },
    { key: 'padMu', label: 'Pad / sheath friction coefficient', unit: '-', type: 'number', default: 0.15, min: 0.03, max: 0.6, group: 'Tensioner', role: 'operation', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Wet polyurethane pad on polyethylene sheath is low (0.1–0.2); verify by a grip test.' },
    { key: 'tensL', label: 'Tensioner contact length', unit: 'm', type: 'number', default: 9, min: 0.5, max: 30, group: 'Tensioner', role: 'operation' },
    { key: 'gripSF', label: 'Grip safety factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Tensioner', role: 'criteria' },
    { key: 'Hs', label: 'Significant wave height during lay', unit: 'm', type: 'number', default: 2, min: 0, max: 8, group: 'Sea state and vessel motion', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'Tp', label: 'Motion period at the chute', unit: 's', type: 'number', default: 8, min: 3, max: 25, group: 'Sea state and vessel motion', role: 'bc', bind: 'site.Tp1' },
    { key: 'raoHeave', label: 'Chute heave RAO (incl. pitch lever)', unit: 'm/m', type: 'number', default: 0.8, min: 0, max: 5, group: 'Sea state and vessel motion', role: 'bc', uncertain: { dist: 'normal', cov: 0.2 } },
    { key: 'raoSurge', label: 'Chute surge RAO', unit: 'm/m', type: 'number', default: 0.3, min: 0, max: 5, group: 'Sea state and vessel motion', role: 'bc' },
    { key: 'ampFactor', label: 'Design wave amplitude / Hs', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1.1, group: 'Sea state and vessel motion', role: 'criteria', help: 'ζa = 0.9·Hs in DNV-RP-N103 §4.3.2 (most probable largest amplitude in a short operation).' },
    { key: 'mbrInstall', label: 'Installation minimum bend radius', unit: 'm', type: 'number', default: 8.2, min: 0.3, max: 60, group: 'Mechanical limits', role: 'criteria', bind: 'umbilical.mbrInstall' },
    { key: 'Tallow', label: 'Allowable installation tension', unit: 'kN', type: 'number', default: 300, min: 5, max: 5000, group: 'Mechanical limits', role: 'criteria', bind: 'umbilical.allowTension' },
    { key: 'swpLimit', label: 'Allowable sidewall pressure on the chute', unit: 'kN/m', type: 'number', default: 25, min: 1, max: 500, group: 'Mechanical limits', role: 'criteria', help: 'Line load T/R the lay-up tolerates on a curved support (crush analysis or manufacturer data).' },
    { key: 'utaMass', label: 'Termination head mass', unit: 't', type: 'number', default: 6, min: 0.1, max: 200, group: 'Termination landing', role: 'load' },
    { key: 'utaVol', label: 'Termination head displaced volume', unit: 'm³', type: 'number', default: 1.2, min: 0.01, max: 200, group: 'Termination landing', role: 'load' },
    { key: 'wireEA', label: 'Lowering wire axial stiffness', unit: 'MN', type: 'number', default: 60, min: 1, max: 2000, group: 'Termination landing', role: 'material', help: 'A&R or crane wire used to lower the termination head.' },
    { key: 'wireM', label: 'Lowering wire mass', unit: 'kg/m', type: 'number', default: 6, min: 0.1, max: 200, group: 'Termination landing', role: 'material' },
    { key: 'utaCdA', label: 'Termination head vertical drag area Cd·A', unit: 'm²', type: 'number', default: 5, min: 0, max: 500, group: 'Termination landing', role: 'load' },
    { key: 'vLand', label: 'Landing velocity', unit: 'm/s', type: 'number', default: 0.4, min: 0, max: 3, group: 'Termination landing', role: 'ic', help: 'Lowering speed plus heave-induced velocity at touchdown of the termination head.' },
    { key: 'kLand', label: 'Seabed / mudmat contact stiffness', unit: 'kN/m', type: 'number', default: 1500, min: 10, max: 1e6, group: 'Termination landing', role: 'soil' },
    { key: 'landCap', label: 'Allowable landing load', unit: 'kN', type: 'number', default: 250, min: 1, max: 1e5, group: 'Termination landing', role: 'criteria' },
    { key: 'routeKm', label: 'Route length', unit: 'km', type: 'number', default: 5, min: 0.05, max: 300, group: 'Schedule', role: 'operation', bind: 'cable.length' },
    { key: 'fixedDays', label: 'Mobilisation, load-out, first- and second-end operations', unit: 'days', type: 'number', default: 6, min: 0, max: 90, group: 'Schedule', role: 'operation' },
    { key: 'operability', label: 'Weather operability for the lay limit', unit: '-', type: 'number', default: 0.8, min: 0.05, max: 1, group: 'Schedule', role: 'operation', bind: 'install.operabilityFrac', help: 'Fraction of time the sea state is below the lay limit (from the metocean persistence analysis).' },
    { key: 'dayRate', label: 'Indicative lay-spread day rate', unit: 'USD', type: 'number', default: 180000, min: 0, max: 2e6, group: 'Schedule', role: 'cost', help: 'Indicative unit rate for the economic comparison only — replace with the contract rate.' },
  ],
  compute(inp) {
    const line = lineFrom({ ...inp, EI: 1, Cdn: 1.2, Cdt: 0.01, Ca: 1 });
    need(line.w > 0, 'The umbilical is buoyant in sea water: it cannot be laid as a catenary without ballast.');
    const h = inp.depth + inp.chuteH, w = line.w;
    // the part above water hangs with its weight in air: equivalent uniform weight over the whole height
    const wEq = (w * inp.depth + line.m * G * inp.chuteH) / h;
    const lay = run(() => cableLay({ w: wEq, h, H: inp.H * 1e3, chuteRadius: inp.chuteR, mbr: inp.mbrInstall, swpLimit: inp.swpLimit * 1e3, payoutSpeed: 1 + inp.slackPct / 100, vesselSpeed: 1 }), 'Lay catenary');
    const Tt = lay.topTension, per = inp.Tp;
    const dynAt = (hs) => layDynamicTension({ staticTension: Tt, w: wEq, EA: line.EA, m: line.m + 0 * line.ma, suspendedLength: lay.suspendedLength, topAngleDeg: lay.departureAngleDeg, heaveAmp: inp.raoHeave * inp.ampFactor * hs, surgeAmp: inp.raoSurge * inp.ampFactor * hs, period: per });
    const dyn = dynAt(inp.Hs);
    const Tmax = Tt + dyn.lower, Tmin = Tt - dyn.lower;
    // lowest admissible top tension: touchdown radius a = H/w must not fall below the MBR → Ttop ≥ w(h + MBR)
    const TminAllow = wEq * (h + inp.mbrInstall), TmaxAllow = Math.min(inp.Tallow * 1e3, inp.swpLimit * 1e3 * inp.chuteR);
    const hsGrid = linspace(0, 8, 81);
    const lim = TmaxAllow > Tt && TminAllow < Tt ? limitingSeaState({ staticTension: Tt, maxTension: TmaxAllow, minTension: TminAllow, hsGrid, dafOfHs: (hs) => dynAt(hs).dafLower }) : { limitingDaf: 1, allowableHs: 0 };
    // static deployment envelope in terms of bottom tension
    const Hmin = wEq * inp.mbrInstall, Hmax = TmaxAllow - wEq * h;
    const envOk = Hmax > Hmin;
    const lbMin = cableLay({ w: wEq, h, H: Hmin }).layBack, lbMax = envOk ? cableLay({ w: wEq, h, H: Hmax }).layBack : lbMin;
    const hl = run(() => cableHandlingLimits({ outerDiameter: line.D, tension: Tmax, bendRadius: inp.chuteR, sidewallLimit: inp.swpLimit * 1e3, squeeze: inp.squeeze * 1e3, squeezeLimit: inp.squeezeLimit * 1e3, nTracks: Number(inp.nTracks), mu: inp.padMu, safetyFactor: inp.gripSF }), 'Tensioner');
    const utaW = (inp.utaMass * 1e3 - line.rhoW * inp.utaVol) * G;
    need(utaW > 0, 'The termination head is buoyant: check its mass and displaced volume.');
    const land = terminationLandingLoad({ mass: inp.utaMass * 1e3, addedMass: line.rhoW * inp.utaVol, submergedWeight: utaW, velocity: inp.vLand, stiffness: inp.kLand * 1e3, capacity: inp.landCap * 1e3 });
    // lowering-line dynamics just above the seabed: SDOF head on the elastic wire, linearised drag (DNV-RP-N103 §5.3)
    const low = run(() => heaveResponse({ mass: inp.utaMass * 1e3, addedMass: line.rhoW * inp.utaVol, EA: inp.wireEA * 1e6, length: Math.max(inp.depth, 5), cableMass: inp.wireM,
      heaveAmp: inp.raoHeave * inp.ampFactor * inp.Hs, period: per, Cd: 1, area: inp.utaCdA, rho: line.rhoW, staticTension: utaW }), 'Lowering-line response');
    const layDays = inp.routeKm * 1e3 / inp.vVessel / 24, days = inp.fixedDays + layDays / inp.operability;
    const checks = [
      check('Overboarding radius (chute)', inp.mbrInstall, inp.chuteR, 'm', `${STD_UMB} §15; DNV-RP-0360 §4`, 'Demand is the installation MBR, capacity the chute radius. A chute tighter than the MBR overbends the tubes under lay tension.'),
      check('Touchdown bend radius', inp.mbrInstall, lay.tdpBendRadius, 'm', `${STD_UMB} §15; catenary a = H/w`, 'Demand is the installation MBR, capacity the catenary radius at touchdown (H/w). Too little bottom tension lets the umbilical kink or loop on the seabed.'),
      check('Sidewall pressure on the chute (dynamic)', Tmax / inp.chuteR / 1e3, inp.swpLimit, 'kN/m', 'DNV-RP-0360 (2016) §4 (SWP = T/R)', 'Line load pressing the umbilical onto the chute at the dynamic top tension. Above the limit the lay-up is crushed against the curved support.'),
      check('Maximum dynamic top tension', Tmax / 1e3, inp.Tallow, 'kN', `${STD_UMB} §15 (installation tension within the capacity curve)`, 'Static top tension plus the dynamic amplitude in the lay sea state against the allowable installation tension.'),
      check('Minimum dynamic top tension (touchdown radius kept)', TminAllow / 1e3, Math.max(Tmin, 1) / 1e3, 'kN', 'Catenary relation Ttop = w(h + a) with a ≥ MBR', 'Demand is the top tension at which the touchdown radius equals the MBR, capacity the lowest dynamic top tension. If the trough falls below it the touchdown zone overbends or goes into compression each cycle.'),
      check('Tensioner squeeze', inp.squeeze, inp.squeezeLimit, 'kN/m', `${STD_UMB} §7 (crush load) — limit from the crush analysis`, 'Radial squeeze per track against the allowable crush load of the lay-up. Excess squeeze ovalises tubes and damages cables.'),
      check('Tensioner grip length', hl.gripLength, inp.tensL, 'm', 'Capstan / Coulomb grip: L = SF·T/(μ·n·q)', 'Contact length needed to hold the dynamic top tension by friction with the safety factor, against the tensioner length. Too short and the umbilical slips through the tensioner.'),
      check('Lowering wire: no slack', low.dynamicForce / 1e3, 0.9 * utaW / 1e3, 'kN', 'DNV-RP-N103 (2021) §4.4.3 / §5.3 (dynamic force ≤ 0.9 × static to avoid slack)', 'Dynamic force in the lowering wire from vessel heave against 90 % of the submerged weight of the termination head. Beyond it the wire goes slack and snatches, risking the head and its umbilical tail.'),
      check('Termination landing load', land.peakLoad / 1e3, inp.landCap, 'kN', 'DNV-RP-N103 (2021) §5–6 (landing on seabed, energy method)', 'Peak contact force when the termination head lands at the given velocity against its structural / foundation allowable.'),
    ];
    const gov = worst(checks);
    const Hs = linspace(Math.max(0.2, Hmin / 1e3 * 0.3), Math.max(inp.H * 2, (envOk ? Hmax : Hmin) / 1e3 * 1.1), 40);
    const curve = Hs.map((v) => cableLay({ w: wEq, h, H: v * 1e3 }));
    const dafs = hsGrid.filter((_, i) => i % 2 === 0).map((hs) => dynAt(hs).dafLower);
    const fuelT = days * MGO_T_PER_DAY;
    return {
      kpis: [
        { key: 'topTension', label: 'Static top tension at the chute', value: round(Tt / 1e3), unit: 'kN', explain: 'Catenary tension H + w·h at the departure point: what the tensioner must hold in calm water.' },
        { key: 'Tmax', label: 'Maximum dynamic top tension', value: round(Tmax / 1e3), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Static tension plus the quasi-static and inertia amplitude in the lay sea state.' },
        { key: 'daf', label: 'Dynamic amplification (estimate)', value: round(dyn.dafLower), unit: '-', explain: 'Ratio of maximum to static top tension when the catenary follows the chute motion. The drag-locked elastic bound is reported separately and is far higher for a stiff umbilical.' },
        { key: 'dafUpper', label: 'Drag-locked upper bound of the amplification', value: round(dyn.dafUpper), unit: '-', explain: 'Bound that assumes drag prevents any change of shape within a cycle, so the motion along the line is taken by elastic stretch. Real response lies between the two; resolve it with a time-domain lay analysis if the bound matters.' },
        { key: 'departure', label: 'Departure angle from horizontal', value: round(lay.departureAngleDeg), unit: 'deg', explain: 'Angle of the umbilical leaving the chute; steep angles (above about 80°) mean little bottom tension and a tight touchdown bend.' },
        { key: 'layBack', label: 'Lay-back', value: round(lay.layBack), unit: 'm', explain: 'Horizontal distance from the chute to the touchdown point — the quantity the ROV monitors to control bottom tension.' },
        { key: 'suspended', label: 'Suspended length', value: round(lay.suspendedLength), unit: 'm', explain: 'Length of umbilical hanging between chute and seabed.' },
        { key: 'tdpRadius', label: 'Touchdown bend radius', value: round(lay.tdpBendRadius), unit: 'm', limit: inp.mbrInstall, sense: 'min', explain: 'Catenary radius H/w at the seabed. It must stay above the installation MBR.' },
        { key: 'swp', label: 'Sidewall pressure on the chute (static)', value: round(lay.sidewallPressure / 1e3), unit: 'kN/m', limit: inp.swpLimit, sense: 'max', explain: 'Line load T/R on the chute in calm water.' },
        { key: 'Hmin', label: 'Lay window: minimum bottom tension', value: round(Hmin / 1e3), unit: 'kN', explain: 'Bottom tension at which the touchdown radius equals the MBR: the lower edge of the allowable deployment envelope.' },
        { key: 'Hmax', label: 'Lay window: maximum bottom tension', value: round(Math.max(Hmax, 0) / 1e3), unit: 'kN', explain: 'Bottom tension at which the top tension reaches the allowable tension or the chute sidewall limit: the upper edge of the envelope.' },
        { key: 'layBackWindow', label: 'Allowable lay-back window', value: round(Math.max(lbMax - lbMin, 0)), unit: 'm', explain: `Range of lay-back (${fmt(lbMin)} to ${fmt(lbMax)} m) inside the static envelope. A narrow window demands tight vessel positioning and touchdown monitoring.` },
        { key: 'hsLimit', label: 'Limiting sea state for lay', value: round(lim.allowableHs), unit: 'm', explain: 'Largest Hs at which both the maximum and the minimum dynamic top tension stay inside the envelope at the chosen bottom tension. It is the operational limit before the alpha factor of DNV-ST-N001 is applied.' },
        { key: 'gripLength', label: 'Required tensioner grip length', value: round(hl.gripLength), unit: 'm', limit: inp.tensL, sense: 'max', explain: 'Contact length needed to hold the dynamic tension by pad friction with the safety factor.' },
        { key: 'lowerDaf', label: 'Lowering-wire dynamic amplification', value: round(low.daf, 3), unit: '-', explain: `Maximum over static wire tension while the termination head hangs near the seabed (wire natural period ${fmt(low.naturalPeriod)} s against a motion period of ${per} s). Close to resonance the value rises sharply.` },
        { key: 'landing', label: 'Termination landing load', value: round(land.peakLoad / 1e3), unit: 'kN', limit: inp.landCap, sense: 'max', explain: 'Peak seabed contact force of the termination head: submerged weight plus the force that absorbs its kinetic energy.' },
        { key: 'vesselDays', label: 'Lay-spread duration', value: round(days), unit: 'days', explain: 'Fixed operations plus lay time divided by weather operability.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Deployment envelope: top tension vs bottom tension', x: { label: 'Bottom tension H', unit: 'kN' }, y: { label: 'Top tension', unit: 'kN' },
          series: [{ name: 'Static top tension', x: Hs, y: curve.map((c) => c.topTension / 1e3) }, { name: 'Upper limit (tension / sidewall)', x: [Hs[0], Hs[Hs.length - 1]], y: [TmaxAllow / 1e3, TmaxAllow / 1e3], dash: true }, { name: 'Lower limit (touchdown MBR)', x: [Hmin / 1e3, Hmin / 1e3], y: [0, TmaxAllow / 1e3], dash: true },
            { name: 'Operating point ± dynamic', x: [inp.H, inp.H, inp.H], y: [Tmin / 1e3, Tt / 1e3, Tmax / 1e3] }],
          explain: 'The lay is acceptable between the vertical line (touchdown radius = MBR) and the horizontal line (allowable top tension or sidewall pressure). The operating point with its dynamic range must sit inside.' },
        { type: 'line', title: 'Touchdown radius and lay-back vs bottom tension', x: { label: 'Bottom tension H', unit: 'kN' }, y: { label: 'Length', unit: 'm' },
          series: [{ name: 'Touchdown radius', x: Hs, y: curve.map((c) => c.tdpBendRadius) }, { name: 'Lay-back', x: Hs, y: curve.map((c) => c.layBack) }, { name: 'Installation MBR', x: [Hs[0], Hs[Hs.length - 1]], y: [inp.mbrInstall, inp.mbrInstall], dash: true }],
          explain: 'Both grow with bottom tension. Read the lay-back the ROV should see for a target bottom tension, and the lay-back below which the touchdown radius violates the MBR.' },
        { type: 'line', title: 'Dynamic amplification vs sea state', x: { label: 'Hs', unit: 'm' }, y: { label: 'Tmax / Tstatic', unit: '-' },
          series: [{ name: 'Estimated DAF', x: hsGrid.filter((_, i) => i % 2 === 0), y: dafs }, { name: 'Limiting DAF', x: [0, 8], y: [lim.limitingDaf, lim.limitingDaf], dash: true }],
          explain: 'Amplification grows linearly with wave height. The crossing with the limiting value — set by whichever of the maximum- or minimum-tension limits is closer — is the limiting sea state for the lay.' },
      ],
      tables: [{ title: 'Lay summary and deployment sequence', explain: 'Stage-by-stage quantities for the installation procedure: what is controlled, the calculated value and the limit it is held against.',
        columns: [{ key: 'stage', label: 'Stage / item' }, { key: 'value', label: 'Value' }, { key: 'limit', label: 'Limit' }, { key: 'unit', label: 'Unit' }],
        rows: [
          { stage: 'Load-out / tensioner: squeeze per track', value: inp.squeeze, limit: inp.squeezeLimit, unit: 'kN/m' },
          { stage: 'Tensioner: grip length', value: round(hl.gripLength), limit: inp.tensL, unit: 'm' },
          { stage: 'Overboarding: chute radius', value: inp.chuteR, limit: `≥ ${inp.mbrInstall}`, unit: 'm' },
          { stage: 'Overboarding: sidewall pressure (dynamic)', value: round(Tmax / inp.chuteR / 1e3), limit: inp.swpLimit, unit: 'kN/m' },
          { stage: 'Normal lay: top tension static / max', value: `${fmt(Tt / 1e3)} / ${fmt(Tmax / 1e3)}`, limit: inp.Tallow, unit: 'kN' },
          { stage: 'Normal lay: bottom tension window', value: inp.H, limit: `${fmt(Hmin / 1e3)} – ${fmt(Math.max(Hmax, 0) / 1e3)}`, unit: 'kN' },
          { stage: 'Normal lay: lay-back window', value: round(lay.layBack), limit: `${fmt(lbMin)} – ${fmt(lbMax)}`, unit: 'm' },
          { stage: 'Normal lay: limiting Hs', value: inp.Hs, limit: round(lim.allowableHs), unit: 'm' },
          { stage: 'Second end: lowering-wire dynamic force', value: round(low.dynamicForce / 1e3), limit: round(0.9 * utaW / 1e3), unit: 'kN' },
          { stage: 'Second end: termination landing load', value: round(land.peakLoad / 1e3), limit: inp.landCap, unit: 'kN' },
        ] }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Allowable bottom-tension window ${fmt(Hmin / 1e3)}–${fmt(Math.max(Hmax, 0) / 1e3)} kN (lay-back ${fmt(lbMin)}–${fmt(lbMax)} m); limiting sea state Hs ≈ ${fmt(lim.allowableHs)} m at H = ${inp.H} kN.`,
        `Pay-out slack ${inp.slackPct} % adds ${fmt(inp.routeKm * 10 * inp.slackPct)} m of product over the route.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}; static top tension ${fmt(Tt / 1e3)} kN, dynamic ${fmt(Tmin / 1e3)}–${fmt(Tmax / 1e3)} kN in Hs ${inp.Hs} m, limiting Hs ${fmt(lim.allowableHs)} m.`,
          gov.uc > 1 ? (gov.name.includes('Minimum dynamic') || gov.name.includes('Touchdown') ? `Raise the bottom tension to about ${fmt((Hmin + dyn.lower) / 1e3 * 1.05)} kN so that the trough of the cycle keeps the touchdown radius above the MBR, or lower the lay sea-state limit.` : gov.name.includes('grip') ? `Raise the squeeze towards its limit or use a longer / additional tensioner: ${fmt(hl.gripLength)} m of contact is needed.` : 'Adjust the governing item (larger chute, lower tension or lower sea state) and re-run.')
            : `Write the bottom-tension window and lay-back window into the lay procedure, with touchdown monitoring by ROV; confirm the limiting sea state with a time-domain lay analysis using the vessel RAOs.`,
          `${STD_UMB} §15; DNV-RP-0360 (2016) §4; DNV-ST-N001 (2021) weather-restricted operations`),
        adv('economic', 'info', `The campaign takes about ${fmt(days)} lay-spread days (${fmt(layDays)} days of laying at ${inp.vVessel} m/h, operability ${inp.operability}); at the indicative rate of ${inp.dayRate} USD/day that is about ${fmt(days * inp.dayRate / 1e6)} million USD.`,
          `Weather downtime is ${fmt(layDays / inp.operability - layDays)} days. Raising the limiting sea state (now ${fmt(lim.allowableHs)} m) by increasing bottom tension within the window, or laying in the calmer season, cuts waiting on weather directly.`,
          'DNV-ST-N001 (2021) operational limits and alpha factors; AACE RP 18R-97 (indicative unit rate, Class 4)',
          `${vesselCo2Text(days)} Shortening the campaign by one day saves about ${fmt(MGO_T_PER_DAY * CO2_MGO)} t CO₂.`),
      ],
      publish: { 'umbilical.layTension': Tmax, 'umbilical.layHs': lim.allowableHs, 'install.vesselDays': days },
      quantities: { cable_km: inp.routeKm * (1 + inp.slackPct / 100), vessel_days: { umbilicalLay: days } },
    };
  },
  verify: [
    // Catenary (Irvine §1.2) with chute at the water line: w = (40.9 − 1025·π·0.145²/4)·9.80665 = 235.11 N/m, h = 350 m,
    // H = 6 kN → a = H/w = 25.520 m; Ttop = H + w·h = 6000 + 82 288 = 88.29 kN; S = √(h(h + 2a)) = √(350 × 401.04) = 374.65 m;
    // lay-back = a·asinh(S/a) = 25.520 × asinh(14.681) = 25.520 × 3.3805 = 86.27 m; cos α = H/T → α = 86.10°.
    { name: 'Lay catenary: Ttop = H + w·h, a = H/w, lay-back = a·asinh(S/a)',
      source: 'Irvine, Cable Structures (1981) §1.2; hand calculation',
      inputs: { chuteH: 0, H: 6, Hs: 0 },
      expect: [{ kpi: 'topTension', value: 88.29, tol: 0.002 }, { kpi: 'tdpRadius', value: 25.52, tol: 0.002 }, { kpi: 'layBack', value: 86.27, tol: 0.003 }, { kpi: 'departure', value: 86.10, tol: 0.002 }, { kpi: 'daf', value: 1, tol: 1e-9 }] },
  ],
};
