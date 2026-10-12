// Subsea lowering and landing: static hook load against depth, axial resonance of the wire–object system, dynamic
// amplification under crane-tip heave, slack-wire check, heave compensation, current offset and landing velocity.
import { staticTension, axialResonance, heaveResponse, heaveCompensator, pendulumPeriod, currentOffset, landingVelocity } from '../../engine/lines/lowering.js';
import { G, need, run, round, fmt, kpi, check, adv, line, tbl, worst, linspace, vesselCo2Text, CONFIRM, N103, LIMIT_N103 } from './_common.js';

export default {
  id: 's15-subsea-lowering', suite: 15, title: 'Subsea lowering and landing',
  summary: 'Follows a subsea structure down the water column on the crane wire: static and dynamic wire tension, the depth at which heave resonance occurs, slack-wire risk, horizontal offset in current and the landing velocity on the seabed.',
  covers: [9], fidelity: 'analytical',
  standards: [`${N103} §5 (deepwater lowering operations) and §6 (landing on seabed)`, 'DNV-ST-N001 (2021) Sec. 16 (lifting: wire and crane capacity)', 'ISO 13628-1 / API RP 17A (subsea structure installation tolerances)'],
  models: ['Cable tension relationships (static hook load with wire weight)', 'Nonlinear lowering-line dynamics reduced to a single degree of freedom with linearised quadratic drag', 'Axial resonance of an elastic wire with end mass (exact and lumped)', 'Vessel heave coupling through prescribed crane-tip motion', 'Passive heave-compensator transmissibility', 'Pendulum motion', 'Steady current offset', 'Landing (touchdown) velocity'],
  limits: `${LIMIT_N103} Single harmonic crane-tip motion and a one-degree-of-freedom vertical model: no wave forces on the object (it is assumed below the wave zone), no horizontal dynamics, no wire bending or sheave friction. The crane-tip motion is an input (from the vessel-motion module); a time-domain coupled lowering analysis is needed where the resonance depth lies within the lowering range.`,
  publishes: ['install.loweringTensionMax', 'install.landingVelocity'],
  searchTerms: 'deepwater lowering subsea installation wire resonance heave compensation landing velocity DNV-RP-N103',
  inputs: [
    { key: 'mass', label: 'Object mass in air', unit: 't', type: 'number', default: 180, min: 0.5, max: 20000, group: 'Object', role: 'load', bind: 'subsea.liftMass', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'vol', label: 'Displaced volume', unit: 'm³', type: 'number', default: 28, min: 0, max: 20000, group: 'Object', role: 'geometry' },
    { key: 'A33', label: 'Heave added mass', unit: 't', type: 'number', default: 320, min: 0, max: 1e5, group: 'Object', role: 'load', uncertain: { dist: 'lognormal', cov: 0.25 }, help: `Mudmats and suction cans trap large water masses: often one to three times the mass in air. ${CONFIRM}` },
    { key: 'Ap', label: 'Projected horizontal area', unit: 'm²', type: 'number', default: 120, min: 0.1, max: 5000, group: 'Object', role: 'geometry' },
    { key: 'Cd', label: 'Vertical drag coefficient', unit: '-', type: 'number', default: 2.5, min: 0.5, max: 8, group: 'Object', role: 'load', uncertain: { dist: 'normal', cov: 0.2 }, help: `Not less than 2.5 in oscillatory flow unless documented. ${CONFIRM}` },
    { key: 'AdH', label: 'Horizontal drag area Cd·A of the object', unit: 'm²', type: 'number', default: 60, min: 0, max: 5000, group: 'Object', role: 'load' },
    { key: 'depth', label: 'Water depth at the landing site', unit: 'm', type: 'number', default: 800, min: 20, max: 4000, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'Uc', label: 'Current (uniform over depth)', unit: 'm/s', type: 'number', default: 0.4, min: 0, max: 3, group: 'Site', role: 'environment', bind: 'site.Uc1' },
    { key: 'dWire', label: 'Lowering wire diameter', unit: 'mm', type: 'number', default: 100, min: 10, max: 180, group: 'Lowering system', role: 'geometry', design: true },
    { key: 'mat', label: 'Wire construction', type: 'select', source: 'materials:wire', default: 'wire-sixstrand', group: 'Lowering system', role: 'material' },
    { key: 'sfWire', label: 'Wire safety factor on MBL', unit: '-', type: 'number', default: 3.0, min: 1.5, max: 6, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'craneCap', label: 'Crane / winch dynamic capacity at this configuration', unit: 't', type: 'number', default: 400, min: 1, max: 20000, group: 'Lowering system', role: 'bc' },
    { key: 'kTip', label: 'Crane tip / boom stiffness', unit: 'kN/m', type: 'number', default: 20000, min: 100, max: 1e6, group: 'Lowering system', role: 'bc', help: 'In series with the wire; governs the system stiffness at shallow depth.' },
    { key: 'heave', label: 'Crane-tip vertical motion single amplitude', unit: 'm', type: 'number', default: 0.8, min: 0, max: 8, group: 'Crane-tip motion', role: 'bc', bind: 'install.craneTipHeave', design: true, help: 'Boundary condition: prescribed motion of the wire upper end in the limiting sea state (from the vessel-motion module).' },
    { key: 'Tct', label: 'Crane-tip motion period', unit: 's', type: 'number', default: 8, min: 3, max: 25, group: 'Crane-tip motion', role: 'bc', bind: 'install.craneTipPeriod' },
    { key: 'ahc', label: 'Heave compensation', type: 'select', default: 'active', group: 'Lowering system', role: 'operation', options: [{ value: 'none', label: 'None' }, { value: 'passive', label: 'Passive compensator (spring–damper)' }, { value: 'active', label: 'Active heave compensation (residual motion)' }] },
    { key: 'kPhc', label: 'Passive compensator stiffness', unit: 'kN/m', type: 'number', default: 150, min: 1, max: 1e5, group: 'Lowering system', role: 'bc', showIf: { key: 'ahc', in: ['passive'] } },
    { key: 'zetaPhc', label: 'Passive compensator damping ratio', unit: '-', type: 'number', default: 0.2, min: 0, max: 2, group: 'Lowering system', role: 'bc', showIf: { key: 'ahc', in: ['passive'] } },
    { key: 'ahcRes', label: 'Residual motion with active compensation', unit: '%', type: 'number', default: 15, min: 1, max: 100, group: 'Lowering system', role: 'bc', showIf: { key: 'ahc', in: ['active'] }, help: 'Fraction of the crane-tip motion that the active system leaves at the hook; 10–20 % is typical in its rated sea state.' },
    { key: 'vLower', label: 'Lowering speed at landing', unit: 'm/s', type: 'number', default: 0.1, min: 0, max: 2, group: 'Landing', role: 'ic', help: 'Initial condition of the touchdown: winch pay-out speed in the final approach.' },
    { key: 'vLandMax', label: 'Allowable landing velocity', unit: 'm/s', type: 'number', default: 0.5, min: 0.05, max: 3, group: 'Criteria', role: 'criteria', help: 'From the structure and foundation design (mudmat bearing, skirt penetration, connector tolerances).' },
    { key: 'offMax', label: 'Allowable horizontal offset at landing', unit: 'm', type: 'number', default: 15, min: 0.5, max: 500, group: 'Criteria', role: 'criteria', help: 'What the vessel can take out by repositioning before guide-post or target-box tolerances apply.' },
    { key: 'hoist', label: 'Hoist length in air (crane tip to object CoG)', unit: 'm', type: 'number', default: 35, min: 2, max: 200, group: 'Lowering system', role: 'geometry' },
    { key: 'nDepth', label: 'Depth stations', unit: '-', type: 'int', default: 40, min: 5, max: 400, group: 'Numerics', role: 'control' },
  ],
  compute(inp, ctx) {
    // t → kg, mm → m, kN/m → N/m.
    const wire = ctx.material(inp.mat); need(wire, `Wire '${inp.mat}' is not in the materials library.`);
    const rho = 1025, d = inp.dWire / 1e3;
    const M = inp.mass * 1e3, A33 = inp.A33 * 1e3, mL = wire.massCoeff * d * d, EA = wire.eaCoeff * d * d, MBL = wire.mblCoeff * d * d, wSub = mL * G * 0.87;
    const depths = linspace(20, inp.depth, inp.nDepth);
    const st = run(() => staticTension({ mass: M, volume: inp.vol, cableWeight: wSub, depths, rho }), 'Static tension');
    const W = st.objectWeightInWater;
    need(W > 0, 'The object is buoyant in water (displaced volume too large for its mass): it cannot be lowered on a wire without added weight.');
    // Crane-tip motion seen by the wire top after any compensation.
    let eta = inp.heave, comp = 1;
    if (inp.ahc === 'active') comp = inp.ahcRes / 100;
    else if (inp.ahc === 'passive') { const k = inp.kPhc * 1e3, mm = M + A33; comp = run(() => heaveCompensator({ mass: mm, k, c: inp.zetaPhc * 2 * Math.sqrt(k * mm), period: inp.Tct }), 'Heave compensator').transmissibility; }
    eta *= Math.min(comp, 3);
    const kc = inp.kTip * 1e3;
    const rowsOut = [], Tn = [], dafs = [], tmax = [], tmin = [];
    let TmaxAll = 0, dafMax = 0, zRes = NaN, slackAny = false, last = null, zDafMax = 0;
    for (let i = 0; i < depths.length; i++) {
      const L = depths[i];
      const r = run(() => heaveResponse({ mass: M, addedMass: A33, EA, length: L, cableMass: mL, kc, heaveAmp: eta, period: inp.Tct, Cd: inp.Cd, area: inp.Ap, rho, staticTension: W }), 'Heave response');
      const top = st.rows[i].tension + r.dynamicForce;
      Tn.push(r.naturalPeriod); dafs.push(r.daf); tmax.push(top / 1e3); tmin.push(r.minTension / 1e3);
      if (top > TmaxAll) TmaxAll = top;
      if (r.daf > dafMax) { dafMax = r.daf; zDafMax = L; }
      if (r.slack) slackAny = true;
      if (i > 0 && Number.isNaN(zRes) && (Tn[i - 1] - inp.Tct) * (Tn[i] - inp.Tct) <= 0) zRes = depths[i - 1] + (depths[i] - depths[i - 1]) * (inp.Tct - Tn[i - 1]) / (Tn[i] - Tn[i - 1]);
      rowsOut.push({ z: round(L), Ts: round(st.rows[i].tension / 1e3), Tn: round(r.naturalPeriod, 3), x: round(r.motionAmp, 3), daf: round(r.daf, 3), Tmax: round(top / 1e3), Tmin: round(r.minTension / 1e3) });
      last = r;
    }
    const res = run(() => axialResonance({ mass: M, addedMass: A33, EA, length: inp.depth, cableMass: mL }), 'Axial resonance');
    const off = run(() => currentOffset({ length: inp.depth, objectWeight: W, cableWeight: wSub, cableDiameter: d, objectDragArea: inp.AdH, current: inp.Uc, rho }), 'Current offset');
    const land = landingVelocity({ loweringSpeed: inp.vLower, objectMotionAmp: last.motionAmp, period: inp.Tct });
    const Tpend = pendulumPeriod({ length: inp.hoist });
    const minT = Math.min(...tmin) * 1e3;
    const checks = [
      check('Wire: maximum dynamic tension against MBL / safety factor', TmaxAll / 1e3, MBL / inp.sfWire / 1e3, 'kN', `${N103} §5.3; DNV-ST-N001 Sec. 16 (wire safety factor)`, 'Largest tension at the crane tip over the whole lowering, including the weight of the paid-out wire. Exceedance risks parting the wire with the load in mid-water.'),
      check('Crane / winch dynamic capacity', TmaxAll / G / 1e3, inp.craneCap, 't', 'DNV-ST-N001 Sec. 16 (crane capacity incl. deployed wire weight)', 'In deep water the wire weight uses up a growing share of the crane capacity; the check is at the worst depth.'),
      check('Slack wire: dynamic force against static object weight', (1 - minT / W) * W / 1e3, 0.9 * W / 1e3, 'kN', `${N103} §4.4.3 / §5.3 (dynamic force not above 90 % of the static submerged weight)`, 'If the dynamic force unloads the wire it goes slack and then snatches, with loads several times the static value.'),
      check('Landing velocity', land.velocity, inp.vLandMax, 'm/s', `${N103} §6.2 (landing on seabed: impact velocity)`, 'Vertical velocity at touchdown: lowering speed plus the object heave velocity. High values overload mudmats, skirts and connectors on impact.'),
      check('Horizontal offset in current at full depth', off.offset, inp.offMax, 'm', `${N103} §5.2.2 (horizontal offset due to current)`, 'Steady drift of the object down-current of the crane tip. It must be within what vessel repositioning can correct before landing.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const resInRange = Number.isFinite(zRes);
    const dReq = Math.sqrt(TmaxAll * inp.sfWire / wire.mblCoeff) * 1e3;
    return {
      kpis: [
        kpi('W', 'Object weight in water', W / 1e3, 'kN', 'Mass less displaced water: the static load on the wire at the object.'),
        kpi('Tstat', 'Static top tension at full depth', st.rows[st.rows.length - 1].tension / 1e3, 'kN', 'Object weight in water plus the submerged weight of the paid-out wire.'),
        kpi('Tmax', 'Maximum dynamic top tension', TmaxAll / 1e3, 'kN', 'Largest static-plus-dynamic tension at the crane tip over the lowering range.', MBL / inp.sfWire / 1e3, 'max'),
        kpi('dafMax', 'Maximum dynamic amplification at the object', dafMax, '-', '(static + dynamic)/static at the object. Peaks where the natural period of the wire–object system meets the crane-tip period.'),
        kpi('zDaf', 'Depth of maximum amplification', zDafMax, 'm', 'Lowered length at which the dynamic force is largest; pass through it without stopping.'),
        kpi('zRes', 'Resonance depth for the crane-tip period', resInRange ? zRes : 0, 'm', 'Lowered length at which the axial natural period equals the crane-tip motion period. Zero when resonance is not reached within the water depth.'),
        kpi('T0', 'Axial natural period at full depth (lumped)', res.periodLumped, 's', '2π·√((M + A33 + mL/3)/K) with the wire stiffness EA/L. Compare with the wave-period band of the operation.'),
        kpi('T0exact', 'Axial natural period at full depth (elastic rod)', res.periodExact, 's', 'Exact solution for a wire of distributed mass with an end mass; the lumped value should agree closely unless the wire mass dominates.'),
        kpi('minT', 'Minimum wire tension at the object', minT / 1e3, 'kN', 'Static object weight less the dynamic force at the worst depth; it must stay well above zero.', 0.1 * W / 1e3, 'min'),
        kpi('vLand', 'Landing velocity', land.velocity, 'm/s', 'Lowering speed plus object heave velocity amplitude at the seabed.', inp.vLandMax, 'max'),
        kpi('offset', 'Horizontal offset in current', off.offset, 'm', 'Offset of the object from the vertical below the crane tip at full depth.', inp.offMax, 'max'),
        kpi('comp', 'Motion passed by the heave compensation', comp, '-', 'Ratio of wire-top motion to crane-tip motion; 1.0 without compensation. A passive system amplifies near its own natural period.'),
        kpi('Tpend', 'Pendulum period in air', Tpend, 's', 'Swing period of the object on the hoist before it enters the water.'),
        kpi('dReq', 'Required wire diameter', dReq, 'mm', 'Diameter at which the maximum dynamic tension equals MBL over the safety factor for this construction.'),
      ],
      checks,
      plots: [
        line('Wire tension against lowered depth', 'Tension at the crane tip', 'kN', 'Depth', 'm', [{ name: 'Static', x: st.rows.map((r) => round(r.tension / 1e3)), y: depths.map((z) => -round(z)) }, { name: 'Static + dynamic', x: tmax.map((v) => round(v)), y: depths.map((z) => -round(z)) }, { name: 'MBL / safety factor', x: [MBL / inp.sfWire / 1e3, MBL / inp.sfWire / 1e3], y: [0, -inp.depth], dash: true }],
          'Top tension grows with depth as wire is paid out; the dynamic envelope bulges at the resonance depth. The envelope must stay left of the wire limit.'),
        line('Natural period and amplification against depth', 'Depth', 'm', 'Period (s) · DAF (-)', '', [{ name: 'Axial natural period (s)', x: depths.map((z) => round(z)), y: Tn.map((v) => round(v, 4)) }, { name: 'Crane-tip period (s)', x: [depths[0], inp.depth], y: [inp.Tct, inp.Tct], dash: true }, { name: 'DAF at the object', x: depths.map((z) => round(z)), y: dafs.map((v) => round(v, 4)) }],
          'The natural period lengthens with depth. Where it crosses the crane-tip period the amplification peaks: this is the depth band to pass quickly or to compensate.'),
      ],
      tables: [tbl('Lowering stations', [['z', 'Lowered length', 'm'], ['Ts', 'Static top tension', 'kN'], ['Tn', 'Natural period', 's'], ['x', 'Object motion amplitude', 'm'], ['daf', 'DAF at object'], ['Tmax', 'Max top tension', 'kN'], ['Tmin', 'Min tension at object', 'kN']], rowsOut,
        'One row per depth station: static tension, the natural period of the wire–object system and the dynamic response to the crane-tip motion. Use it to mark hold points and no-stop zones in the procedure.')],
      notes: [resInRange ? `Resonance with the ${inp.Tct} s crane-tip motion occurs near ${fmt(zRes)} m lowered length.` : `The axial natural period stays ${Tn[Tn.length - 1] < inp.Tct ? 'below' : 'above'} the crane-tip period over the whole depth: no resonance crossing.`, slackAny ? 'Slack wire is predicted at one or more depths.' : 'No slack wire predicted.', `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 || resInRange ? 'warn' : 'ok', `Lowering to ${inp.depth} m is within limits: maximum tension ${fmt(TmaxAll / 1e3)} kN, peak DAF ${fmt(dafMax)} at ${fmt(zDafMax)} m, landing at ${fmt(land.velocity)} m/s with ${fmt(off.offset)} m offset.`,
          resInRange ? `Do not hold the load between about ${fmt(0.8 * zRes)} and ${fmt(1.2 * zRes)} m; check the range of periods around ${inp.Tct} s and consider heave compensation through that band.` : 'Check the full range of crane-tip periods for the operating sea states, since a longer period moves the resonance into the lowering range.', `${N103} §5.3`)
          : adv('engineering', 'fail', `"${gov.name}" is exceeded (${fmt(gov.uc)}): maximum tension ${fmt(TmaxAll / 1e3)} kN, minimum ${fmt(minT / 1e3)} kN, landing velocity ${fmt(land.velocity)} m/s, offset ${fmt(off.offset)} m.`,
            gov.name.startsWith('Wire') ? `Use a wire of at least ${fmt(Math.ceil(dReq / 2) * 2)} mm, or reduce the crane-tip motion with heave compensation or a lower limiting sea state.` : gov.name.startsWith('Slack') || gov.name.startsWith('Landing') ? `Limit the crane-tip motion to about ${fmt(inp.heave / gov.uc)} m (lower sea state or heave compensation), or reduce added mass with perforated mudmats.` : gov.name.startsWith('Horizontal') ? 'Wait for a weaker current, add clump weight, or plan the vessel offset up-current.' : 'Use a crane or winch with more capacity at depth, or a lighter (fibre) lowering line.',
            `${N103} §5–6`),
        adv('economic', 'info', `Lowering through ${inp.depth} m at typical winch speeds takes one to two hours; the limiting crane-tip motion (${inp.heave} m) sets the workable sea state.`,
          inp.ahc === 'none' ? `Heave compensation would cut the dynamic force by the factor shown for the compensator and widen the weather window: compare its hire against the waiting-on-weather days it removes (weather-window module).` : `Compensation passes ${fmt(comp * 100)} % of the crane-tip motion; confirm that this performance is guaranteed at the natural period of ${fmt(res.periodLumped)} s before relying on a higher limiting sea state.`,
          'Construction vessel day rate × waiting on weather (DNV-ST-N001 Sec. 2); value of the wider window against equipment hire',
          `Each day of waiting on weather for this operation costs about ${vesselCo2Text(1, 35)}`),
      ],
      publish: { 'install.loweringTensionMax': TmaxAll, 'install.landingVelocity': land.velocity },
      quantities: { vessel_days: { subseaConstruction: 0.5 } },
    };
  },
  verify: [
    { name: 'Static top tension: F = (M − ρ·V)·g + w·L', source: 'DNV-RP-N103 §5.2.1, hand calculation: (180 000 − 1025·28)·9.80665 = 1483.7 kN; wire 90 mm six-strand: 4000·0.09² = 32.4 kg/m, submerged 0.87·32.4·9.80665 = 276.43 N/m; at 1200 m: 1483.7 + 331.7 = 1815.5 kN',
      inputs: { depth: 1200, dWire: 90, ahc: 'none' }, expect: [{ kpi: 'W', value: (180000 - 1025 * 28) * 9.80665 / 1e3, tol: 1e-9 }, { kpi: 'Tstat', value: ((180000 - 1025 * 28) * 9.80665 + 0.87 * 32.4 * 9.80665 * 1200) / 1e3, tol: 1e-9 }] },
    { name: 'Axial natural period: T0 = 2π·√((M + A33 + m·L/3)·L/EA)', source: 'DNV-RP-N103 §5.3.3 lumped form: M + A33 = 500 t, m·L/3 = 12.96 t, EA = 4.5e10·0.09² = 364.5 MN, L = 1200 m → K = 303.75 kN/m → T0 = 2π·√(512 960/303 750) = 8.165 s',
      inputs: { depth: 1200, dWire: 90, ahc: 'none' }, expect: [{ kpi: 'T0', value: 2 * Math.PI * Math.sqrt((500e3 + 32.4 * 1200 / 3) * 1200 / (4.5e10 * 0.0081)), tol: 1e-9 }, { kpi: 'T0exact', value: 2 * Math.PI * Math.sqrt((500e3 + 32.4 * 1200 / 3) * 1200 / (4.5e10 * 0.0081)), tol: 0.002 }] },
    { name: 'Landing velocity with no crane-tip motion equals the lowering speed', source: 'Limiting case: zero heave → object at rest relative to the wire, v = v_lower = 0.1 m/s',
      inputs: { heave: 0 }, expect: [{ kpi: 'vLand', value: 0.1, tol: 1e-9 }, { kpi: 'dafMax', value: 1, tol: 1e-9 }] },
  ],
};
