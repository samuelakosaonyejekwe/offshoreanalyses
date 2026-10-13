// Coupled mooring–hull–riser station keeping: mean offset under wind, current and wave drift on the quasi-static
// catenary mooring with riser restoring, low-frequency slow-drift response, wave-frequency motion from the hull
// RAOs, the API combination of offsets and line tensions, and the riser reaction at the extreme offset.
import { chainProperties } from '../../engine/lines/lineprops.js';
import { spreadLayout, equilibrium, stiffnessMatrix, restoringCurve, quasiDynamicTension, lineForces } from '../../engine/lines/mooringsystem.js';
import { meanDriftForce, newmanDriftSpectrum, slowDriftResponse } from '../../engine/hydro/drift.js';
import { hullInputs, hullFrom, motionInputs, motionFrom, motionRangeAdvice, seaInputs, seaFrom, spectrum, omegaGrid, stats, windCurrentInputs, windCurrentForce, driftCoefficient,
  G, DEG, LIMIT_HYDRO, need, run, round, rnd, fmt, adv, check, sev, worst, linspace, geomspace, steelCo2Text } from './_common.js';
import { boxFroudeKrylov } from './wave-loads.js';

export default {
  id: 's11-coupled-mooring-riser',
  suite: 11,
  title: 'Mooring–hull–riser coupled station keeping',
  summary: 'Mean offset, slow-drift and wave-frequency excursions of the moored unit with its risers, the resulting extreme line tension against the break load, and the riser reaction and top angle at the extreme offset.',
  covers: [10],
  fidelity: 'reduced-order',
  standards: ['API RP 2SK (3rd ed., 2005; reaffirmed 2015) §5–7', 'DNV-OS-E301 (2021) Ch.2 Sec.2', 'ISO 19901-7 (2013) §8–10', 'API RP 2RD (2nd ed., 2013) / API STD 2RD (riser interface)', 'DNV-RP-C205 (2019) §7.4'],
  models: ['Elastic catenary mooring lines, quasi-static system equilibrium', 'Linearised mooring stiffness matrix by finite differences', 'Riser restoring as a tensioned-string spring T/L with static hang-off pull', 'Mean wave drift and Newman-approximation slow-drift force spectrum', 'Single-degree-of-freedom low-frequency response with linear damping', 'Wave-frequency motion from the hull RAOs', 'API RP 2SK combination of mean, low-frequency and wave-frequency responses'],
  limits: `Frequency-domain, quasi-static coupling: line dynamics enter only through a user amplification factor, risers act as linear springs with a static pull, and the low-frequency damping is a user fraction of critical. Mean drift uses a screening coefficient, not quadratic transfer functions. ${LIMIT_HYDRO} A fully coupled time-domain analysis with line and riser dynamics is outside the module.`,
  publishes: ['vessel.offsetMax', 'vessel.Fmean'],
  searchTerms: 'coupled mooring riser floater offset slow drift API RP 2SK',
  inputs: [
    ...hullInputs({ type: 'semi', withMesh: false }),
    ...seaInputs({ Hs: 12, Tp: 15, ret: '100', heading: 0 }),
    ...windCurrentInputs({ Uw: 36, Uc: 1.0, ret: '100' }),
    ...motionInputs().filter((i) => !['Kmoor'].includes(i.key)),
    { key: 'nLines', label: 'Number of mooring lines', unit: '-', type: 'int', default: 12, min: 3, max: 24, group: 'Mooring system', role: 'bc', design: true, help: 'Lines equally spaced in azimuth (spread mooring).' },
    { key: 'chainD', label: 'Chain diameter', unit: 'mm', type: 'number', default: 127, min: 40, max: 220, group: 'Mooring system', role: 'geometry', design: true, help: 'Nominal diameter of the studless chain; mass, stiffness and break load follow from indicative estimators.' },
    { key: 'grade', label: 'Chain grade', type: 'select', default: 'R4', options: ['R3', 'R3S', 'R4', 'R4S', 'R5'].map((g) => ({ value: g, label: g })), group: 'Mooring system', role: 'material', help: 'Offshore mooring chain grade (IACS UR W22).' },
    { key: 'lineL', label: 'Line length', unit: 'm', type: 'number', default: 1700, min: 100, max: 8000, group: 'Mooring system', role: 'geometry', design: true, help: 'Unstretched length of each line from anchor to fairlead.' },
    { key: 'anchorR', label: 'Anchor radius', unit: 'm', type: 'number', default: 1655, min: 100, max: 8000, group: 'Mooring system', role: 'bc', design: true, help: 'Horizontal distance from the platform centre to the anchors; with the line length it sets the pretension.' },
    { key: 'fairR', label: 'Fairlead radius', unit: 'm', type: 'number', default: 50, min: 1, max: 200, group: 'Mooring system', role: 'geometry', help: 'Horizontal distance from the platform centre to the fairleads.' },
    { key: 'fairZ', label: 'Fairlead depth below waterline', unit: 'm', type: 'number', default: 15, min: -20, max: 200, group: 'Mooring system', role: 'geometry', help: 'Positive below the still waterline.' },
    { key: 'nRisers', label: 'Number of risers', unit: '-', type: 'int', default: 8, min: 0, max: 60, group: 'Riser interface', role: 'bc', help: 'Risers and umbilicals hung off the unit.' },
    { key: 'riserT', label: 'Top tension per riser', unit: 'kN', type: 'number', default: 1800, min: 0, max: 30000, group: 'Riser interface', role: 'bc', bind: 'riser.topTension', help: 'Effective tension at the hang-off of each riser.' },
    { key: 'riserH0', label: 'Net static horizontal riser pull', unit: 'kN', type: 'number', default: 250, min: 0, max: 20000, group: 'Riser interface', role: 'bc', bind: 'riser.reactionH', help: 'Resultant horizontal force of all risers on the unit at zero offset (zero for a symmetric arrangement).' },
    { key: 'riserDirDeg', label: 'Direction of the static riser pull', unit: 'deg', type: 'number', default: 180, min: 0, max: 360, group: 'Riser interface', role: 'bc', help: 'Direction towards which the risers pull the unit (same convention as the wave heading).' },
    { key: 'riserAngle0', label: 'Riser top angle at zero offset', unit: 'deg', type: 'number', default: 6, min: 0, max: 30, group: 'Riser interface', role: 'ic', bind: 'riser.hangoffAngleDeg', help: 'Static hang-off angle from vertical.' },
    { key: 'lfDamping', label: 'Low-frequency damping (fraction of critical)', unit: '-', type: 'number', default: 0.15, min: 0.02, max: 0.6, group: 'Dynamics', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Surge/sway damping of the slow-drift mode from current drag on hull, mooring lines and risers and from wave-drift damping (typically 0.1–0.3).' },
    { key: 'daf', label: 'Dynamic amplification of wave-frequency line tension', unit: '-', type: 'number', default: 1.3, min: 1, max: 4, group: 'Dynamics', role: 'bc', help: 'Ratio of dynamic to quasi-static wave-frequency tension, from a line-dynamics analysis (1.0 = quasi-static).' },
    { key: 'sfLine', label: 'Required safety factor on line tension', unit: '-', type: 'number', default: 1.67, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria', help: 'API RP 2SK intact condition: 1.67 for dynamic analysis, 2.0 for quasi-static.' },
    { key: 'offLimPct', label: 'Allowable offset (% of water depth)', unit: '%', type: 'number', default: 12, min: 1, max: 30, group: 'Acceptance criteria', role: 'criteria', help: 'Riser-driven excursion limit in the survival condition: 5–10 % of depth for rigid risers, 10–15 % for flexibles.' },
    { key: 'riserAngleMax', label: 'Allowable riser top angle', unit: 'deg', type: 'number', default: 15, min: 2, max: 40, group: 'Acceptance criteria', role: 'criteria', help: 'Flex-joint or bend-stiffener limit on the top angle from vertical.' },
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [50, 100, 200], qoi: 'wfSig' },
  compute(inp) {
    const hull = hullFrom(inp); const depth = inp.depth; const rho = hull.rho;
    const ch = run(() => chainProperties({ d: inp.chainD / 1e3, grade: inp.grade, type: 'studless', rhoWater: rho }), 'Chain properties');
    const hF = depth - inp.fairZ; const span = inp.anchorR - inp.fairR;
    need(span > 0, 'Anchor radius must exceed the fairlead radius.');
    need(inp.lineL > hF, `Under-length line: ${fmt(inp.lineL)} m cannot reach from the seabed to the fairlead ${fmt(hF)} m above it.`);
    need(inp.lineL * (1 + ch.MBL / ch.EA) > Math.hypot(span, hF), `Over-constrained: anchors ${fmt(Math.hypot(span, hF), 5)} m from the fairleads cannot be reached by a ${fmt(inp.lineL)} m line even at breaking strain. Reduce the anchor radius or lengthen the lines.`);
    need(inp.lineL < span + hF, `Slack system: the line (${fmt(inp.lineL)} m) is longer than the seabed run plus the fairlead height (${fmt(span + hF)} m), so it carries no pretension. Shorten the lines or move the anchors out.`);
    const system = run(() => spreadLayout({ nLines: inp.nLines, anchorRadius: inp.anchorR, fairleadRadius: inp.fairR, fairleadZ: -inp.fairZ, depth, heading: 0.5 * 2 * Math.PI / inp.nLines, segments: [{ L: inp.lineL, w: ch.w, EA: ch.EA, name: 'chain' }], MBL: ch.MBL, CB: 0 }), 'Mooring layout');
    const pose0 = [0, 0, 0, 0, 0, 0];
    const f0 = run(() => lineForces(system, pose0), 'Mooring pretension');
    const pret = f0.lines.reduce((a, l) => a + l.T, 0) / f0.lines.length;
    // risers: tensioned-string lateral stiffness n·T/L and the static pull
    const Kr = inp.nRisers * inp.riserT * 1e3 / depth;
    const rd = inp.riserDirDeg * DEG; const H0 = inp.riserH0 * 1e3;
    const beta = inp.headingDeg; const cb = Math.cos(beta * DEG); const sb = Math.sin(beta * DEG);
    const sea = seaFrom(inp); const om = omegaGrid(inp.nFreq); const S = spectrum(om, sea); const dur = inp.durationH * 3600;
    const wc = windCurrentForce(inp, hull, beta);
    const Dw = driftCoefficient(inp, hull, om, beta);
    const Fdrift = meanDriftForce({ omega: om, S, D: Dw });
    const Fwind = Math.hypot(wc.wind.Fx, wc.wind.Fy); const Fcur = Math.hypot(wc.cur.Fx, wc.cur.Fy);
    const Fenv = [wc.wind.Fx + wc.cur.Fx + Fdrift * cb, wc.wind.Fy + wc.cur.Fy + Fdrift * sb];
    const Fmean = Math.hypot(Fenv[0], Fenv[1]);
    const hydroK = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]); hydroK[0][0] = Kr; hydroK[1][1] = Kr;
    const Fext = [Fenv[0] + H0 * Math.cos(rd), Fenv[1] + H0 * Math.sin(rd), 0, 0, 0, 0];
    const eq = run(() => equilibrium(system, { Fext, dofs: [0, 1, 5], pose0, hydroK }), 'Mean offset equilibrium');
    need(eq.converged, 'The mean equilibrium did not converge: the environmental force probably exceeds what the mooring can restore. Stiffen the mooring or reduce the load.');
    const dirOff = eq.offset > 1e-6 ? eq.offsetDirection : beta * DEG; const dx = Math.cos(dirOff); const dy = Math.sin(dirOff);
    const Km = run(() => stiffnessMatrix(system, { pose: eq.pose }), 'Mooring stiffness');
    const Kmoor = Km[0][0] * dx * dx + 2 * Km[0][1] * dx * dy + Km[1][1] * dy * dy;
    const Ktot = Kmoor + Kr;
    need(Ktot > 0, 'The station-keeping system has no restoring stiffness at the mean offset.');
    const mot = motionFrom({ ...inp, Kmoor: Ktot / 1e3 }, hull);
    const r = mot.rao(om, beta, sea);
    // wave-frequency excursion along the offset direction
    let wfAmp; let mEff;
    if (mot.kind === 'slender') {
      wfAmp = om.map((_, i) => Math.hypot(dx * r.re[0][i] + dy * r.re[1][i], dx * r.im[0][i] + dy * r.im[1][i]));
      mEff = (mot.model.M[0][0] + mot.model.A[0][0]) * dx * dx + (mot.model.M[1][1] + mot.model.A[1][1]) * dy * dy;
    } else {
      // box hull: mass-controlled horizontal response to the Froude–Krylov force (added mass 10 % surge, 80 % sway — screening)
      mEff = hull.mass * (1 + 0.1 * dx * dx + 0.8 * dy * dy);
      wfAmp = om.map((w) => { const f = boxFroudeKrylov({ L: hull.Lc, B: hull.Bc, T: hull.T, omega: w, headingDeg: beta, rho }); return Math.hypot(f[0] * dx / (1.1 * hull.mass), f[1] * dy / (1.8 * hull.mass)) / (w * w); });
    }
    const wf = stats(om, wfAmp, S, dur);
    const mu = geomspace(0.002, 0.25, 120);
    const SF = newmanDriftSpectrum({ omega: om, S, D: Dw, mu });
    const lf = run(() => slowDriftResponse({ mass: mEff, stiffness: Ktot, dampingRatio: inp.lfDamping, mu, SF, meanForce: 0, duration: dur }), 'Slow-drift response');
    const lfMax = lf.sigma * Math.sqrt(2 * Math.log(Math.max(lf.nCycles, 2)));
    const lfSig = 2 * lf.sigma;
    const dyn = Math.max(lfMax + wf.significant, lfSig + wf.mpm);
    const offMax = eq.offset + dyn;
    const qd = run(() => quasiDynamicTension(system, { meanPose: eq.pose, direction: dirOff, lf: { sig: lfSig, max: lfMax }, wf: { sig: wf.significant, max: wf.mpm }, daf: inp.daf }), 'Line tensions at the extreme offset');
    const Tmax = qd.maxTension; const Tmean = eq.maxTension;
    const riserR = H0 + Kr * offMax;
    const angle = inp.riserAngle0 + Math.atan(offMax / depth) / DEG;
    const offLim = inp.offLimPct / 100 * depth;
    const checks = [
      check('Maximum line tension (intact)', Tmax * inp.sfLine / 1e3, ch.MBL / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3 Table 1; ISO 19901-7 (2013) §10.2', 'Guards against breaking the most loaded line in the design storm; a line failure shifts load to its neighbours and can cascade.'),
      check('Maximum offset against the riser limit', offMax, offLim, 'm', 'API RP 2SK §6.2; riser watch circle per API STD 2RD (project limit)', 'Guards the risers: beyond the allowable excursion they overstress at the touchdown or the hang-off, or must be disconnected.'),
      check('Riser top angle at the extreme offset', angle, inp.riserAngleMax, 'deg', 'API STD 2RD (2013) §5 (flex-joint / stress-joint rotation limits; project value)', 'Guards the riser hang-off hardware, which accepts only a limited rotation before the flex joint bottoms out or the bend stiffener overstresses.'),
    ];
    const gov = worst(checks);
    const offs = linspace(0, Math.max(1.3 * offMax, 0.03 * depth), 14);
    const rc = run(() => restoringCurve(system, { direction: dirOff, offsets: offs }), 'Restoring curve');
    const Tn = 2 * Math.PI * Math.sqrt(mEff / Ktot);
    const chainT = ch.mass * inp.lineL * inp.nLines / 1e3;
    // chain size that would just meet the tension factor: MBL ∝ d² to first order
    const dReq = inp.chainD * Math.sqrt(Math.max(checks[0].uc, 1e-6));
    const lineRows = qd.lines.map((l, i) => ({ line: i + 1, pre: round(f0.lines[i].T / 1e3, 5), mean: round(l.Tmean / 1e3, 5), max: round(l.Tmax / 1e3, 5), sf: round(ch.MBL / Math.max(l.Tmax, 1), 4), gov: l.governing }));
    return {
      kpis: [
        { key: 'Fmean', label: 'Mean environmental force', value: round(Fmean / 1e3, 5), unit: 'kN', explain: 'Steady wind, current and wave-drift force on the unit, acting along the environment direction.' },
        { key: 'offsetMean', label: 'Mean offset', value: round(eq.offset, 4), unit: 'm', explain: 'Static excursion where mooring and riser restoring balance the mean force and the static riser pull.' },
        { key: 'lfSigma', label: 'Low-frequency offset, standard deviation', value: round(lf.sigma, 4), unit: 'm', explain: 'Slow-drift surge about the mean at the natural period of the moored system; inversely proportional to the square root of damping.' },
        { key: 'wfSig', label: 'Wave-frequency offset, significant amplitude', value: round(wf.significant, 4), unit: 'm', explain: 'First-order horizontal motion along the offset direction, from the hull RAOs.' },
        { key: 'offsetMax', label: 'Maximum offset', value: round(offMax, 4), unit: 'm', limit: round(offLim, 4), sense: 'max', explain: 'Mean plus the larger of (LF maximum + WF significant) and (LF significant + WF maximum), per API RP 2SK.' },
        { key: 'offsetPct', label: 'Maximum offset / water depth', value: round(100 * offMax / depth, 4), unit: '%', explain: 'The usual way of quoting excursion for riser design.' },
        { key: 'pretension', label: 'Line pretension', value: round(pret / 1e3, 5), unit: 'kN', explain: 'Fairlead tension of each line at zero offset. Typically 10–20 % of the break load for a chain catenary.' },
        { key: 'Tmean', label: 'Most loaded line, mean tension', value: round(Tmean / 1e3, 5), unit: 'kN', explain: 'Largest fairlead tension at the mean offset.' },
        { key: 'Tmax', label: 'Most loaded line, maximum tension', value: round(Tmax / 1e3, 5), unit: 'kN', limit: round(ch.MBL / inp.sfLine / 1e3, 5), sense: 'max', explain: 'Largest fairlead tension at the extreme offset with the dynamic amplification on the wave-frequency part.' },
        { key: 'sfLine', label: 'Line safety factor (MBL / Tmax)', value: round(ch.MBL / Tmax, 4), unit: '-', limit: inp.sfLine, sense: 'min', explain: 'Break load over maximum tension; API RP 2SK requires 1.67 intact with dynamic analysis.' },
        { key: 'Kmoor', label: 'Mooring stiffness at the mean offset', value: round(Kmoor / 1e3, 5), unit: 'kN/m', explain: 'Tangent horizontal stiffness of the mooring along the offset direction; it rises as the leeward lines lift off.' },
        { key: 'Kriser', label: 'Riser restoring stiffness', value: round(Kr / 1e3, 5), unit: 'kN/m', explain: 'Lateral stiffness of the tensioned risers, n·T/depth. In deep water with many risers it is a real share of the total.' },
        { key: 'TnSurge', label: 'Natural period of the moored system', value: round(Tn, 4), unit: 's', explain: 'Surge/sway period on the combined stiffness with added mass. It should lie well above wave periods (60–300 s).' },
        { key: 'riserReaction', label: 'Riser horizontal reaction at the extreme offset', value: round(riserR / 1e3, 5), unit: 'kN', explain: 'Static pull plus restoring force of the risers at the extreme offset: the load the hang-off structure must take.' },
        { key: 'riserAngle', label: 'Riser top angle at the extreme offset', value: round(angle, 4), unit: 'deg', limit: inp.riserAngleMax, sense: 'max', explain: 'Hang-off angle from vertical after adding the rotation caused by the offset.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Restoring force and line tension against offset', x: { label: 'Offset', unit: 'm' }, y: { label: 'Force', unit: 'kN' },
          series: [{ name: 'Mooring restoring', x: rnd(offs), y: rnd(rc.restoring.map((v) => Math.abs(v) / 1e3)) }, { name: 'Mooring + risers', x: rnd(offs), y: rnd(rc.restoring.map((v, i) => (Math.abs(v) + Kr * offs[i]) / 1e3)) },
            { name: 'Most loaded line tension', x: rnd(offs), y: rnd(rc.maxTension.map((v) => v / 1e3)), dash: true }, { name: 'Mean force', x: [0, round(offs[offs.length - 1])], y: [round(Fmean / 1e3), round(Fmean / 1e3)], dash: true }],
          explain: 'Force–offset curve of the station-keeping system along the offset direction. The mean offset is where the restoring curve crosses the mean force; the hardening of the curve and the line tension beyond it show how close the system is to its limit.' },
        { type: 'bar', title: 'Offset budget', categories: ['Mean', 'Low-frequency max', 'Wave-frequency max', 'Total (API combination)', 'Allowable'], y: { label: 'Offset', unit: 'm' },
          series: [{ name: 'Offset', y: [eq.offset, lfMax, wf.mpm, offMax, offLim].map((v) => round(v, 4)) }],
          explain: 'Contributions to the extreme excursion. A large mean calls for more stiffness or less windage; a large low-frequency part calls for damping or stiffness; the wave-frequency part cannot be changed by the mooring.' },
        { type: 'bar', title: 'Line tensions', categories: lineRows.map((l) => `L${l.line}`), y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Pretension', y: lineRows.map((l) => l.pre) }, { name: 'Mean', y: lineRows.map((l) => l.mean) }, { name: 'Maximum', y: lineRows.map((l) => l.max) }],
          explain: 'Fairlead tension of every line at rest, at the mean offset and at the extreme. Windward lines pick up the load while leeward lines slacken.' },
      ],
      tables: [
        { title: 'Line results', columns: [{ key: 'line', label: 'Line' }, { key: 'pre', label: 'Pretension', unit: 'kN' }, { key: 'mean', label: 'Mean', unit: 'kN' }, { key: 'max', label: 'Maximum', unit: 'kN' }, { key: 'sf', label: 'MBL / max' }, { key: 'gov', label: 'Governing combination' }], rows: lineRows,
          explain: 'One row per mooring line with the tension at each stage and the safety factor on the break load.' },
        { title: 'Force and stiffness components', columns: [{ key: 'c', label: 'Component' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Wind force', Fwind / 1e3, 'kN'], ['Current force', Fcur / 1e3, 'kN'], ['Mean wave drift', Fdrift / 1e3, 'kN'], ['Static riser pull', H0 / 1e3, 'kN'], ['Mooring stiffness', Kmoor / 1e3, 'kN/m'], ['Riser stiffness', Kr / 1e3, 'kN/m'], ['Effective mass incl. added mass', mEff / 1e3, 't'],
            ['Low-frequency damping', lf.damping / 1e3, 'kN·s/m'], ['Chain break load', ch.MBL / 1e3, 'kN'], ['Chain submerged weight', ch.w, 'N/m'], ['Chain axial stiffness EA', ch.EA / 1e6, 'MN']].map(([c, v, u]) => ({ c, v: round(v, 5), u })),
          explain: 'The separate contributions that the coupled solution combines — use it to compare each component with its own stand-alone analysis.' },
      ],
      notes: [
        `Mean force ${fmt(Fmean / 1e3)} kN (wind ${fmt(Fwind / 1e3)}, current ${fmt(Fcur / 1e3)}, drift ${fmt(Fdrift / 1e3)} kN) gives a mean offset of ${fmt(eq.offset)} m; extreme ${fmt(offMax)} m = ${fmt(100 * offMax / depth)} % of depth.`,
        `Stiffness: mooring ${fmt(Kmoor / 1e3)} kN/m + risers ${fmt(Kr / 1e3)} kN/m (${fmt(100 * Kr / Ktot)} % of the total); natural period ${fmt(Tn)} s.`,
        `Chain ${fmt(inp.chainD)} mm ${inp.grade}: break load ${fmt(ch.MBL / 1e3)} kN (indicative estimator), pretension ${fmt(100 * pret / ch.MBL)} % MBL. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        mot.kind === 'box' ? 'Wave-frequency surge and sway of the box hull are mass-controlled estimates from the Froude–Krylov force.' : `Wave-frequency motion from the slender-member RAOs with the solved stiffness.`,
      ],
      advice: [
        ...motionRangeAdvice(mot, inp.Tp, inp.headingDeg),
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}: Tmax ${fmt(Tmax / 1e3)} kN (factor ${fmt(ch.MBL / Tmax)} on MBL), offset ${fmt(offMax)} m against ${fmt(offLim)} m.`,
          gov.uc > 1 ? (gov === checks[0] ? `Increase the chain to about ${fmt(Math.ceil(dReq))} mm or add lines; also check whether a lower pretension (larger anchor radius margin) reduces the peak tension without breaching the offset limit.` : `Stiffen the system: shorten the lines or move the anchors out to raise pretension (now ${fmt(100 * pret / ch.MBL)} % MBL), add lines on the weather side, or use taut polyester in deep water.`)
            : `Margins: tension ${fmt(100 * (1 - checks[0].uc))} %, offset ${fmt(100 * (1 - checks[1].uc))} %. Repeat for the one-line-broken condition (factor 1.25) and for 45° and 90° environments, and confirm the low-frequency damping (${fmt(inp.lfDamping)}), which controls the slow-drift part.`,
          'API RP 2SK (3rd ed.) §5–7; DNV-OS-E301 (2021) Ch.2 Sec.2; ISO 19901-7 (2013).'),
        adv('economic', 'info', `The mooring holds ${fmt(chainT)} t of chain. Tension scales with offset stiffness and chain weight; the present factor is ${fmt(ch.MBL / Tmax)} against ${inp.sfLine} required.`,
          checks[0].uc < 0.8 ? `The chain could come down to about ${fmt(Math.ceil(dReq))} mm on tension alone, saving roughly ${fmt(chainT * (1 - (dReq / inp.chainD) ** 2))} t — check fatigue, corrosion allowance and the damaged case before taking it.` : 'Keep the chain size; spend effort on documenting damping and line dynamics, which move the result more cheaply than steel.',
          'Mooring optimisation practice of API RP 2SK; whole-life costing per ISO 15686-5 (chain replacement interval); AACE RP 18R-97.',
          `${steelCo2Text(chainT)} Each avoided line also saves an anchor and about a day of anchor-handling vessel time.`),
      ],
      publish: { 'vessel.offsetMax': offMax, 'vessel.Fmean': Fmean },
      quantities: { steel_t: round(chainT, 4) },
    };
  },
  verify: [
    { // No wind, current, waves or riser pull: the symmetric spread mooring keeps the unit at the origin and every line carries
      // the same pretension; with risers removed the riser stiffness is zero.
      name: 'Component check: zero environment gives zero offset and equal pretensions',
      source: 'Symmetry and static equilibrium (limiting case of the coupled solution)',
      inputs: { Hs: 0, Uw: 0, Uc: 0, riserH0: 0, nRisers: 0 },
      expect: [{ kpi: 'offsetMean', value: 0, abs: 1e-3 }, { kpi: 'offsetMax', value: 0, abs: 1e-3 }, { kpi: 'Kriser', value: 0, abs: 1e-9 }, { kpi: 'Fmean', value: 0, abs: 1e-6 }] },
    { // Riser spring by hand: 8 risers × 1800 kN / 350 m = 41.14 kN/m. Current only on a barge 100 × 20 × 5 m, heading 0:
      // F = ½·1025·1.0·(20·5)·1.5² = 115.31 kN.
      name: 'Riser tensioned-string stiffness n·T/L and current drag by hand',
      source: 'String statics; drag equation F = ½ρCdAU² (DNV-RP-C205 §6) — hand calculation',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, Hs: 0, Uw: 0, Uc: 1.5, Ccur: 1, riserH0: 0, nRisers: 8, riserT: 1800, depth: 350, fairR: 30 },
      expect: [{ kpi: 'Kriser', value: 41.143, tol: 1e-3 }, { kpi: 'Fmean', value: 115.31, tol: 1e-3 }] },
  ],
};
