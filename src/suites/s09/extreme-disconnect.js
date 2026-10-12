// Extreme storm and accidental disconnection of a drilling riser: operating window against flex-joint limits,
// drift-off and watch circles, emergency disconnect recoil, and storm hang-off of the disconnected string.
import { drillingRiserOffsets, driftOff } from '../../engine/lines/riser.js';
import { axialResonance, heaveResponse } from '../../engine/lines/lowering.js';
import {
  G, DEG, need, run, fmt, round, pipeInputs, waterInputs, hydroInputs, currentInputs, currentFn, pipeFrom, steelTonnes, adv, check, kpi, plot, meshInput, finish,
  steelCo2Text, vesselCo2Text, CONFIRM,
} from './_common.js';

export default {
  id: 's09-extreme-disconnect', suite: 9, title: 'Extreme storm and emergency disconnect',
  summary: 'Safe operating window of a drilling riser: flex-joint angles against vessel offset, drift-off time and watch circles, recoil after emergency disconnect, and tension in storm hang-off.',
  covers: [13], fidelity: 'reduced-order',
  standards: ['API RP 16Q (2017) §3–4 operating limits, emergency disconnect, hang-off', 'ISO 13624-1 (2009) marine drilling riser', 'DNV-RP-N103 (2019) §5.3 axial dynamics', 'IMCA M 220 / DP operations guidance (watch circles)'],
  models: ['Tensioned beam with flex joints and current drag', 'Closed-form drift-off under constant force and quadratic resistance', 'Watch-circle (yellow / red) offsets', 'Rigid-body recoil after disconnect with overpull release', 'Axial mass–spring response of a hung-off riser with linearised drag', 'Tensioner stiffness'],
  limits: 'Quasi-static in-plane riser response per offset; drift-off is one degree of freedom with constant force and resistance coefficient; recoil is a rigid-body estimate without tensioner valve control, mud column discharge or elastic wave propagation; storm hang-off is a single-degree-of-freedom axial model. Wellhead and conductor bending capacity, and a transient nonlinear contact–impact simulation of the disconnect, are outside the module. Operating-limit angles are overridable defaults.',
  publishes: ['riser.watchCircleRed', 'riser.watchCircleYellow', 'riser.offsetLimit'],
  searchTerms: 'drilling riser emergency disconnect drift-off watch circle recoil hang-off',
  inputs: [
    ...pipeInputs({ D: 533.4, t: 22.2, mat: 'api5l-x80', rhoInt: 1300, coat: 0 }),
    ...waterInputs(),
    ...hydroInputs({ Cdn: 1.0 }),
    { key: 'length', label: 'Riser length, lower flex joint to tensioner ring', unit: 'm', type: 'number', default: 980, min: 50, max: 3600, group: 'Configuration', role: 'geometry' },
    { key: 'buoyD', label: 'Buoyancy module outer diameter', unit: 'mm', type: 'number', default: 1250, min: 560, max: 1700, group: 'Buoyancy system', role: 'geometry' },
    { key: 'buoyRho', label: 'Module density', unit: 'kg/m³', type: 'number', default: 500, min: 300, max: 800, group: 'Buoyancy system', role: 'material' },
    { key: 'buoyFill', label: 'Fraction of the riser fitted with buoyancy', unit: '-', type: 'number', default: 0.55, min: 0, max: 1, group: 'Buoyancy system', role: 'geometry', design: true },
    { key: 'auxMass', label: 'Auxiliary lines and couplings mass', unit: 'kg/m', type: 'number', default: 150, min: 0, max: 600, group: 'Buoyancy system', role: 'load', help: 'Choke, kill, booster and hydraulic lines plus flanges, averaged per metre (steel, submerged weight ≈ 0.87 × dry).' },
    { key: 'tensionFactor', label: 'Top tension factor', unit: '-', type: 'number', default: 2.0, min: 1.05, max: 4, group: 'Top end', role: 'bc', design: true, help: 'Top tension over apparent weight of the riser string alone. It must be high enough to leave overpull at the LMRP connector after carrying the LMRP on disconnect.' },
    { key: 'tensionerK', label: 'Tensioner stiffness', unit: 'kN/m', type: 'number', default: 0, min: 0, max: 3000, group: 'Top end', role: 'bc' },
    { key: 'kLower', label: 'Lower flex-joint stiffness', unit: 'kN·m/deg', type: 'number', default: 100, min: 0, max: 5000, group: 'Bottom end', role: 'bc' },
    { key: 'lmrpMass', label: 'LMRP mass', unit: 't', type: 'number', default: 120, min: 10, max: 500, group: 'Bottom end', role: 'load' },
    { key: 'lmrpVol', label: 'LMRP displaced volume', unit: 'm³', type: 'number', default: 15, min: 1, max: 100, group: 'Bottom end', role: 'load' },
    { key: 'lmrpArea', label: 'LMRP vertical drag area', unit: 'm²', type: 'number', default: 25, min: 1, max: 200, group: 'Bottom end', role: 'load' },
    ...currentInputs({ Us: 1.0, Ub: 0.1 }),
    { key: 'offsetMaxPct', label: 'Largest offset in the sweep', unit: '% depth', type: 'number', default: 8, min: 2, max: 15, group: 'Vessel and drift-off', role: 'bc' },
    { key: 'vesselMass', label: 'Vessel mass incl. added mass', unit: 't', type: 'number', default: 60000, min: 5000, max: 400000, group: 'Vessel and drift-off', role: 'bc', bind: 'vessel.disp' },
    { key: 'Fenv', label: 'Net environmental force after loss of position', unit: 'kN', type: 'number', default: 1500, min: 50, max: 20000, group: 'Vessel and drift-off', role: 'bc', bind: 'vessel.Fmean', uncertain: { dist: 'lognormal', cov: 0.25 } },
    { key: 'cDrag', label: 'Hull resistance coefficient', unit: 'kN·s²/m²', type: 'number', default: 600, min: 20, max: 10000, group: 'Vessel and drift-off', role: 'bc', help: 'F = c·v²: ½ρ·Cd·A of the hull in the drift direction (beam-on drillship ≈ 300–1500).' },
    { key: 'edsTime', label: 'Emergency disconnect sequence time', unit: 's', type: 'number', default: 45, min: 10, max: 180, group: 'Vessel and drift-off', role: 'operation' },
    { key: 'reactTime', label: 'Reaction time between yellow and red', unit: 's', type: 'number', default: 60, min: 0, max: 300, group: 'Vessel and drift-off', role: 'operation' },
    { key: 'clearance', label: 'LMRP lift needed to clear the BOP', unit: 'm', type: 'number', default: 3, min: 0.5, max: 15, group: 'Disconnect and hang-off', role: 'operation' },
    { key: 'heaveAmp', label: 'Storm heave amplitude at the hang-off', unit: 'm', type: 'number', default: 4, min: 0, max: 12, group: 'Disconnect and hang-off', role: 'bc', bind: 'vessel.heaveSig', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'heaveT', label: 'Storm heave period', unit: 's', type: 'number', default: 11, min: 4, max: 25, group: 'Disconnect and hang-off', role: 'bc' },
    { key: 'hangCap', label: 'Hang-off (spider / gimbal) capacity', unit: 'kN', type: 'number', default: 9000, min: 500, max: 40000, group: 'Disconnect and hang-off', role: 'criteria' },
    { key: 'meanLimit', label: 'Lower flex-joint mean angle limit (drilling)', unit: 'deg', type: 'number', default: 2, min: 0.5, max: 6, group: 'Acceptance criteria', role: 'criteria', help: `Commonly 2° mean while drilling. ${CONFIRM}` },
    { key: 'maxLimit', label: 'Lower flex-joint maximum angle limit', unit: 'deg', type: 'number', default: 4, min: 1, max: 10, group: 'Acceptance criteria', role: 'criteria', help: `Commonly 4° maximum connected, non-drilling; the physical joint limit is typically 9–10°. ${CONFIRM}` },
    { key: 'physLimit', label: 'Physical limit angle (point of disconnect)', unit: 'deg', type: 'number', default: 9, min: 2, max: 15, group: 'Acceptance criteria', role: 'criteria', help: 'Angle at which equipment damage starts (flex-joint stop, wellhead moment or stroke-out). It defines the offset by which disconnect must be complete.' },
    meshInput(60, '30 / 60 / 120'),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [30, 60, 120, 240], qoi: 'angleRef' },

  compute(inp, ctx) {
    // engineering units → SI: mm → m, t → kg, kN → N, kN·m/deg → N·m/rad
    const outer = inp.D / 1e3;
    need(inp.buoyD / 1e3 > outer, 'Buoyancy module diameter must exceed the riser pipe diameter.');
    const pp = pipeFrom(inp, ctx, { buoyancy: inp.buoyFill > 0 ? { Dout: inp.buoyD / 1e3, rho: inp.buoyRho, fill: inp.buoyFill } : undefined });
    const sec = { ...pp.section };
    // auxiliary lines: mass added, submerged weight reduced by their displaced water (steel: 1 − ρw/ρs)
    sec.m += inp.auxMass; sec.w += inp.auxMass * G * (1 - pp.rhoW / pp.mat.rho);
    need(sec.w > 0, `The riser string is net buoyant (${Math.round(sec.w)} N/m): reduce the buoyancy coverage (now ${inp.buoyFill}).`);
    const L = inp.length; need(L <= inp.depth + 40, `Riser length ${L} m exceeds the water depth plus deck elevation.`);
    const W = sec.w * L, cf = currentFn(inp, inp.depth);
    const offMax = inp.offsetMaxPct / 100 * inp.depth, offsets = [];
    for (let i = 0; i <= 8; i++) offsets.push(offMax * i / 8);
    const kRad = inp.kLower > 0 ? inp.kLower * 1e3 / DEG : undefined;
    const cfg = { length: L, pipe: sec, tensionFactor: inp.tensionFactor, kFlexLower: kRad, tensionerStiffness: inp.tensionerK * 1e3, nElem: inp.nElem, current: { rho: pp.rhoW, profile: (z) => cf(z) } };
    const sw = run(() => drillingRiserOffsets({ ...cfg, offsets, meanLimitDeg: inp.meanLimit, maxLimitDeg: inp.maxLimit }), 'The drilling-riser offset sweep did not solve');
    const ang = sw.rows.map((r) => Math.abs(r.lowerAngleDeg));
    const cross = (lim) => { for (let k = 1; k < ang.length; k++) if ((ang[k - 1] - lim) * (ang[k] - lim) <= 0 && ang[k] !== ang[k - 1]) return offsets[k - 1] + (lim - ang[k - 1]) * (offsets[k] - offsets[k - 1]) / (ang[k] - ang[k - 1]); return ang[0] > lim ? 0 : NaN; };
    const offMean = cross(inp.meanLimit), offMaxA = cross(inp.maxLimit), offPhys = cross(inp.physLimit);
    // beyond the sweep: extrapolate linearly on the last interval (angle is nearly linear in offset)
    const slope = (ang[8] - ang[7]) / (offsets[8] - offsets[7]);
    const ext = (lim, v) => (Number.isFinite(v) ? v : slope > 0 ? offsets[8] + (lim - ang[8]) / slope : offsets[8]);
    const oMean = ext(inp.meanLimit, offMean), oMax = ext(inp.maxLimit, offMaxA), oPhys = ext(inp.physLimit, offPhys);
    need(oPhys > 0, 'The physical limit angle is exceeded at zero offset: the current alone over-rotates the lower flex joint. Increase top tension.');
    const dr = run(() => driftOff({ mass: inp.vesselMass * 1e3, force: inp.Fenv * 1e3, dragCoefficient: inp.cDrag * 1e3, limitOffset: oPhys, edsTime: inp.edsTime, reactionTime: inp.reactTime, nPoints: 60 }), 'Drift-off calculation failed');
    const tRed = dr.timeToLimit - inp.edsTime, tYel = tRed - inp.reactTime;
    // recoil after disconnect: tensioner pull minus string weight accelerates riser + LMRP (rigid body, no tensioner control)
    const Ttop = inp.tensionFactor * W;
    const mL = inp.lmrpMass * 1e3, wL = (mL - pp.rhoW * inp.lmrpVol) * G;
    need(wL > 0, `The LMRP is buoyant with ${inp.lmrpMass} t mass and ${inp.lmrpVol} m³ displaced volume; check the two values.`);
    const mString = sec.m * L + mL, over = Ttop - W - wL;
    const a0 = over / (mString + pp.rhoW * inp.lmrpVol);   // LMRP added mass ≈ displaced mass
    const tClear = a0 > 0 ? Math.sqrt(2 * inp.clearance / a0) : Infinity, vClear = a0 > 0 ? a0 * tClear : 0;
    const stretch = (Ttop - 0.5 * W) * L / sec.EA;   // elastic stretch stored in the string at mean tension
    // storm hang-off: riser + LMRP hanging from the spider, axial mass–spring with linearised drag on the LMRP
    const hangStatic = W + wL;
    const res = axialResonance({ mass: mL, addedMass: pp.rhoW * inp.lmrpVol, EA: sec.EA, length: L, cableMass: sec.m });
    const hr = inp.heaveAmp > 0 ? run(() => heaveResponse({ mass: mL, addedMass: pp.rhoW * inp.lmrpVol, EA: sec.EA, length: L, cableMass: sec.m, heaveAmp: inp.heaveAmp, period: inp.heaveT, Cd: 1.2, area: inp.lmrpArea, rho: pp.rhoW, staticTension: wL }), 'Hang-off dynamic response failed') : { dynamicForce: 0, motionRatio: 1 };
    // top dynamic force: inertia of the whole string following the heave (rigid) scaled by the end-mass amplification
    const om = 2 * Math.PI / inp.heaveT, Fdyn = (sec.m * L + mL + pp.rhoW * inp.lmrpVol) * om * om * inp.heaveAmp * Math.max(1, hr.motionRatio);
    const hangMax = hangStatic + Fdyn, hangMin = hangStatic - Fdyn;
    const angRef = ang[4];

    const kpis = [
      kpi('topTension', 'Top tension', Ttop / 1e3, 'kN', 'Tensioner setting: tension factor times the apparent weight of the riser string in mud. It must also lift the LMRP clear on disconnect.'),
      kpi('overpull', 'Overpull at the LMRP connector', over / 1e3, 'kN', 'Tension left at the bottom of the riser after carrying the string and the LMRP. Positive overpull is needed for a clean lift-off at emergency disconnect.', 0, 'min'),
      kpi('angleZero', 'Lower flex-joint angle at zero offset', ang[0], 'deg', 'Angle caused by current alone with the vessel over the well. It uses up part of the operating limit before any offset.'),
      kpi('angleRef', 'Lower flex-joint angle at half the swept offset', angRef, 'deg', `Angle at an offset of ${fmt(offsets[4])} m, used to follow mesh convergence.`),
      kpi('offsetDrill', 'Offset at the drilling (mean angle) limit', oMean, 'm', `Vessel offset at which the lower flex-joint angle reaches ${inp.meanLimit}°. Drilling should stop outside this radius.`),
      kpi('offsetMax', 'Offset at the maximum connected angle', oMax, 'm', `Offset at which the angle reaches ${inp.maxLimit}°: the edge of the connected, non-drilling envelope.`),
      kpi('offsetPhys', 'Offset at the physical limit (point of disconnect)', oPhys, 'm', `Offset at which the angle reaches ${inp.physLimit}°. The emergency disconnect must be complete before the vessel gets here.`),
      kpi('driftTime', 'Time to reach the physical limit after loss of position', dr.timeToLimit, 's', 'Time for the vessel to drift from the well centre to the point of disconnect under the net environmental force. It must exceed the disconnect sequence time with margin.', inp.edsTime + inp.reactTime, 'min'),
      kpi('redCircle', 'Red watch circle (initiate disconnect)', dr.redOffset, 'm', 'Offset at which the emergency disconnect must be started so that it completes before the physical limit is reached.'),
      kpi('yellowCircle', 'Yellow watch circle (prepare to disconnect)', dr.yellowOffset, 'm', 'Offset reached a reaction time before the red circle: drilling stops and the crew prepares to disconnect.'),
      kpi('vTerminal', 'Terminal drift velocity', dr.terminalVelocity, 'm/s', 'Steady drift speed √(F/c) the vessel approaches. It bounds how fast the offset grows late in a drift-off.'),
      kpi('recoilAcc', 'Initial recoil acceleration', a0, 'm/s²', 'Upward acceleration of the riser and LMRP at the instant of disconnect if the tensioners keep pulling. High values need an anti-recoil system to avoid slamming the telescopic joint.'),
      kpi('clearTime', 'Time for the LMRP to clear the BOP', Number.isFinite(tClear) ? tClear : 1e6, 's', 'Time to lift the LMRP by the required clearance at the initial recoil acceleration. Compare with the vessel heave period: it should clear before the next downward heave.'),
      kpi('recoilVel', 'Recoil velocity at clearance', vClear, 'm/s', 'Upward velocity of the string when the LMRP has cleared; the anti-recoil system must arrest it within the remaining tensioner stroke.'),
      kpi('stretch', 'Elastic stretch released on disconnect', stretch, 'm', 'Axial stretch of the tensioned string that is released when the bottom lets go: the first, fast part of the recoil.'),
      kpi('hangMax', 'Maximum hang-off load in the storm', hangMax / 1e3, 'kN', 'Static weight of the disconnected riser and LMRP plus the inertia force from storm heave (hard hang-off, riser follows the vessel).', inp.hangCap, 'max'),
      kpi('hangMin', 'Minimum hang-off load in the storm', hangMin / 1e3, 'kN', 'Static weight minus the dynamic force. A negative value means the string goes into compression each cycle and the hard hang-off must be replaced by a soft (tensioner) hang-off.', 0, 'min'),
      kpi('axialPeriod', 'Axial natural period of the hung-off string', res.periodExact, 's', 'Natural period of the riser stretching with the LMRP as an end mass. If it approaches the heave period the hang-off loads are amplified by resonance.'),
    ];
    const checks = [
      check('Disconnect time available in a drift-off', inp.edsTime + inp.reactTime, Math.max(dr.timeToLimit, 1e-6), 's', 'API RP 16Q (2017) §4 emergency disconnect; DP watch-circle practice (IMCA M 220)', 'Compares the time needed to react and run the disconnect sequence with the time the drift-off leaves. If the sequence is longer, the riser, wellhead or BOP is damaged before the LMRP is released.'),
      check('Lower flex-joint angle while drilling, vessel on location', ang[0], inp.meanLimit, 'deg', 'API RP 16Q (2017) §3 operating limits (mean angle)', 'Checks the mean flex-joint angle from current alone against the drilling limit; beyond it the drill string wears the flex joint and wellhead.'),
      check('Storm hang-off load', hangMax / 1e3, inp.hangCap, 'kN', 'API RP 16Q (2017) §4 hang-off; equipment rating', 'Guards the spider, gimbal and upper riser joints against the peak dynamic load of the hung-off string in the storm.'),
      check('Compression in the hung-off string', Fdyn / 1e3, hangStatic / 1e3, 'kN', 'API RP 16Q (2017) §4 hang-off (avoid dynamic compression)', 'Compares the dynamic force with the static weight. Above 1.0 the riser unloads completely on the downstroke and then snatches, which can part or buckle the string.'),
    ];

    const plots = [
      plot('Lower flex-joint angle against vessel offset', 'Vessel offset', 'm', 'Lower flex-joint angle', 'deg', [{ name: 'Lower flex-joint angle', x: offsets, y: ang }, { name: 'Upper angle', x: offsets, y: sw.rows.map((r) => Math.abs(r.upperAngleDeg)) }, { name: 'Drilling limit', x: [0, offMax], y: [inp.meanLimit, inp.meanLimit], dash: true }, { name: 'Maximum connected', x: [0, offMax], y: [inp.maxLimit, inp.maxLimit], dash: true }],
        'Flex-joint angle as the vessel moves off location. The offsets where the curve crosses the dashed limits define the drilling and connected envelopes.'),
      plot('Drift-off', 'Time after loss of position', 's', 'Vessel offset', 'm', [{ name: 'Offset', x: dr.t, y: dr.x }, { name: 'Red circle', x: [0, dr.timeToLimit], y: [dr.redOffset, dr.redOffset], dash: true }, { name: 'Yellow circle', x: [0, dr.timeToLimit], y: [dr.yellowOffset, dr.yellowOffset], dash: true }],
        'Vessel offset against time after losing station keeping, up to the point of disconnect. The yellow and red lines are the alarm offsets that leave time for the disconnect sequence.'),
      plot('Drift velocity', 'Time after loss of position', 's', 'Velocity', 'm/s', [{ name: 'Drift velocity', x: dr.t, y: dr.v }], 'Vessel speed during the drift-off, approaching the terminal velocity. A fast early rise means little warning time.'),
      plot('Top tension and stroke against offset', 'Vessel offset', 'm', 'Top tension | stroke', 'kN | m×100', [{ name: 'Top tension (kN)', x: offsets, y: sw.rows.map((r) => r.topTension / 1e3) }, { name: 'Stroke-out (m × 100)', x: offsets, y: sw.rows.map((r) => r.stroke * 100) }], 'Tensioner load and stroke as the riser leans over. Stroke-out of the telescopic joint is the other physical limit besides the flex-joint angle.'),
      { type: 'polar', title: 'Watch circles', series: [['Drilling limit', oMean], ['Yellow', dr.yellowOffset], ['Red', dr.redOffset], ['Point of disconnect', oPhys]].map(([name, r]) => ({ name, thetaDeg: Array.from({ length: 37 }, (_, i) => i * 10), r: Array.from({ length: 37 }, () => r) })),
        explain: 'Plan view of the offsets around the well centre: drilling limit, yellow and red alarm circles and the point of disconnect. Their spacing is the time margin the crew has.' },
    ];
    const tables = [
      { title: 'Offset sweep', columns: [{ key: 'off', label: 'Offset', unit: 'm' }, { key: 'pct', label: 'Offset', unit: '% depth' }, { key: 'low', label: 'Lower angle', unit: 'deg' }, { key: 'up', label: 'Upper angle', unit: 'deg' }, { key: 'T', label: 'Top tension', unit: 'kN' }, { key: 'stroke', label: 'Stroke', unit: 'm' }, { key: 'uc', label: 'Angle / maximum limit', unit: '-' }],
        rows: sw.rows.map((r) => ({ off: round(r.offset, 4), pct: round(100 * r.offset / inp.depth, 3), low: round(Math.abs(r.lowerAngleDeg), 3), up: round(Math.abs(r.upperAngleDeg), 3), T: round(r.topTension / 1e3, 5), stroke: round(r.stroke, 3), uc: round(r.ucMax, 3) })),
        explain: 'Riser response at each vessel offset: flex-joint angles, tension and telescopic stroke, with the utilisation of the maximum connected angle.' },
      { title: 'Safe operating limits', columns: [{ key: 'item', label: 'Limit' }, { key: 'value', label: 'Offset', unit: 'm' }, { key: 'pct', label: 'Offset', unit: '% depth' }, { key: 't', label: 'Time after loss of position', unit: 's' }],
        rows: [['Drilling limit (mean angle)', oMean, NaN], ['Yellow watch circle', dr.yellowOffset, Math.max(0, tYel)], ['Red watch circle', dr.redOffset, Math.max(0, tRed)], ['Maximum connected angle', oMax, NaN], ['Point of disconnect', oPhys, dr.timeToLimit]].map(([item, v, t]) => ({ item, value: round(v, 4), pct: round(100 * v / inp.depth, 3), t: Number.isFinite(t) ? round(t, 4) : '–' })),
        explain: 'Summary of the operating envelope: each limit as an offset, as a percentage of water depth and, for the alarm circles, the time after loss of position at which it is reached.' },
    ];
    const ucMax = Math.max(...checks.map((c) => c.uc)), gov = checks.find((c) => c.uc === ucMax);
    const tfNeed = (W + wL) * 1.05 / W;
    const advice = [
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.8 ? 'warn' : 'ok',
        `${gov.name} governs at ${fmt(ucMax)}. Point of disconnect ${fmt(oPhys)} m (${fmt(100 * oPhys / inp.depth)} % of depth), reached ${fmt(dr.timeToLimit)} s after loss of position; red circle ${fmt(dr.redOffset)} m, yellow ${fmt(dr.yellowOffset)} m; overpull ${fmt(over / 1e3)} kN.`,
        ucMax > 1 ? (gov.name.startsWith('Disconnect') ? `Shorten the sequence (now ${inp.edsTime} s) or widen the envelope with more top tension; alternatively restrict operations to weather in which the net force is below about ${fmt(inp.Fenv * (dr.timeToLimit / (inp.edsTime + inp.reactTime)) ** 2)} kN.`
          : gov.name.startsWith('Lower') ? `Raise the tension factor above ${inp.tensionFactor} or suspend drilling in this current; the angle falls roughly in proportion to 1/tension.`
            : 'Use a soft hang-off on the tensioners, or pull the riser before the storm; hard hang-off is not acceptable for this heave.')
          : over <= 0 ? `Overpull is not positive: raise the tension factor to at least ${fmt(tfNeed)} so the LMRP lifts off on disconnect.`
            : 'Set the DP alarm circles to the red and yellow offsets above, and confirm recoil with the tensioner anti-recoil system model before operations.',
        `API RP 16Q (2017) §3–4. Limits ${inp.meanLimit}° / ${inp.maxLimit}° are overridable defaults. ${CONFIRM}`),
      adv('economic', 'info',
        `The drilling envelope is ${fmt(oMean)} m radius. Every disconnect costs the trip to recover or re-run ${fmt(L)} m of riser — typically 2 to 4 rig days — plus the well downtime; an uncontrolled drift-off beyond the point of disconnect risks the wellhead and BOP.`,
        'Use the watch circles to set weather-specific operating criteria: restricting operations in the few conditions where drift-off time is short is cheaper than a single damaged wellhead, and a wider envelope (more tension) is paid for in tensioner capacity.',
        'Risk-based operating criteria (API RP 16Q, IMCA M 220) and non-productive-time accounting in the drilling AFE.',
        `One avoided riser trip of about 3 rig days saves ${vesselCo2Text(3)} The riser string itself: ${steelCo2Text(steelTonnes(pp.D, pp.t, L, pp.mat.rho))}`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`Riser string apparent weight ${fmt(sec.w)} N/m in ${inp.rhoInt} kg/m³ mud (${fmt(W / 1e3)} kN), LMRP ${fmt(wL / 1e3)} kN submerged; buoyancy on ${fmt(inp.buoyFill * 100)} % of the length.`,
        `${Number.isFinite(offPhys) ? 'Limits found inside the swept offsets.' : 'The physical limit lies beyond the swept offsets and was extrapolated linearly.'} ${gov.name} governs.`],
      advice,
      publish: { 'riser.watchCircleRed': dr.redOffset, 'riser.watchCircleYellow': dr.yellowOffset, 'riser.offsetLimit': oPhys },
      quantities: { steel_t: steelTonnes(pp.D, pp.t, L, pp.mat.rho), vessel_days: { riser_trip: 3 } },
    });
  },

  verify: [
    {
      // Drift-off under constant force with quadratic resistance: terminal velocity v∞ = sqrt(F/c)
      // = sqrt(1500e3 / 600e3) = 1.5811 m/s (first principles, M·ẍ = F − c·ẋ²).
      name: 'Terminal drift velocity v∞ = √(F/c)', source: 'First principles (Newton’s second law with quadratic resistance)',
      inputs: {},
      expect: [{ kpi: 'vTerminal', value: Math.sqrt(1500 / 600), tol: 1e-6 }],
    },
    {
      // Axial period of the hung-off string with a small end mass. Bare rod fixed at the top: T = 4L/c = 4L·√(m/EA),
      // m = steel + mud = π/4·(0.5334² − 0.489²)·7850 + π/4·0.489²·1300 = 279.8 + 244.2 kg/m,
      // EA = 207e9·π/4·(0.5334² − 0.489²) = 7.379e9 N → 4·980·√(524.0/7.379e9) = 1.0446 s.
      // An end mass M' = 10 t + 1.025 t added mass is μ = 2.15 % of the string mass; (ωL/c)·tan(ωL/c) = 1/μ gives
      // ωL/c ≈ (π/2)/(1 + μ), i.e. the period is longer by the factor (1 + μ) ≈ 1.0215 → 1.067 s.
      name: 'Hung-off string axial period: rod with a small end mass, T = (4L/c)·(1 + μ)', source: 'Timoshenko, Vibration Problems in Engineering — longitudinal vibration of a rod with end mass; c = √(EA/m)',
      inputs: { lmrpMass: 10, lmrpVol: 1, buoyFill: 0, auxMass: 0 },
      expect: [{ kpi: 'axialPeriod', value: 4 * 980 * Math.sqrt((Math.PI / 4 * (0.5334 ** 2 - 0.489 ** 2) * 7850 + Math.PI / 4 * 0.489 ** 2 * 1300) / (207e9 * Math.PI / 4 * (0.5334 ** 2 - 0.489 ** 2))) * 1.0215, tol: 0.01 }],
    },
  ],
};
