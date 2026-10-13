// Installation vessel motion: box-barge RAOs in heave, pitch and roll (panel solution or closed form), response spectra in a JONSWAP sea,
// motion of the crane tip (or any deck point), and the limiting sea state for motion criteria.
import { boxBargeRAO, waveSpectrum, responseStatistics } from '../../engine/hydro/floater.js';
import { G, DEG, need, run, round, fmt, kpi, check, adv, line, tbl, worst, linspace, thin, vesselCo2Text, seaInputs, CONFIRM } from './_common.js';

const HEADINGS = [{ value: 180, label: 'Head sea (180°)' }, { value: 135, label: 'Bow quartering (135°)' }, { value: 90, label: 'Beam sea (90°)' }];

export default {
  id: 's15-vessel-motion', suite: 15, title: 'Installation vessel motion',
  summary: 'Computes heave, pitch and roll response of a barge-shaped installation vessel in the operating sea state, the vertical motion, velocity and acceleration at the crane tip, and the sea state at which motion limits are reached.',
  covers: [11], fidelity: 'reduced-order',
  standards: ['DNV-RP-C205 (2019) §3.5 (JONSWAP spectrum), §7 (wave-induced motions and short-term statistics)', 'DNV-RP-N103 (2017) §4.3.3 (crane-tip motion), §9 (operational criteria)', 'DNV-ST-N001 (2021) Sec. 11 (motion criteria for transport), Sec. 2 (operational limits)', 'Newman (1977) Marine Hydrodynamics ch. 6 (radiation and diffraction of a floating body)', 'Jensen, Mansour & Olsen (2004) Estimation of ship motions using closed-form expressions, Ocean Eng. 31 (optional closed-form model)'],
  models: ['Response amplitude operators of a box hull from a three-dimensional panel solution with the free-surface Green function (frequency-dependent added mass, wave damping, Haskind excitation)', 'Coupled sway–roll with viscous roll damping as a fraction of critical', 'Optional closed-form heave and pitch (Jensen, Mansour & Olsen) with one-degree-of-freedom roll', 'JONSWAP wave spectrum', 'Linear spectral response and Rayleigh extremes', 'Rigid-body point motion (crane tip)', 'Hydrostatic restoring and metacentric stability'],
  limits: 'RAOs of a rectangular box with the length, beam, draft and displacement of the vessel in deep water at zero speed, long-crested waves: panel solution of the box (pitch gyradius 0.25 L, roll with the gyradius entered, wave-radiation damping and a user fraction of critical), within 5 % of an independent panel code for a box in waves longer than the beam; a peak period below that wave carries a warning. The closed-form option is within 5 % in heave only for wave periods above twice 2π√(2T/g) (about 15 % low at 1.5 times, 30–40 % low near it). Heave, pitch and roll contributions at a point are combined as the square root of the sum of squares — conservative checks should also consider the plain sum. Ship-shaped hulls, semi-submersibles and moored or DP response need a hydrodynamic database from a diffraction code, which this module does not produce.',
  publishes: ['install.craneTipHeave', 'install.craneTipVelocity', 'install.craneTipPeriod', 'install.craneTipHeavePerHs', 'install.HsOperation', 'install.transportRollDeg', 'install.transportRollT', 'install.transportPitchDeg', 'install.transportHeaveG'],
  searchTerms: 'installation vessel motion RAO crane tip heave roll barge limiting sea state',
  inputs: [
    { key: 'L', label: 'Vessel length', unit: 'm', type: 'number', default: 160, min: 30, max: 400, group: 'Vessel', role: 'geometry' },
    { key: 'B', label: 'Beam', unit: 'm', type: 'number', default: 42, min: 8, max: 120, group: 'Vessel', role: 'geometry' },
    { key: 'T', label: 'Draft', unit: 'm', type: 'number', default: 8, min: 1, max: 30, group: 'Vessel', role: 'ic', help: 'Operating draft: the initial floating condition for the motion analysis.' },
    { key: 'Cb', label: 'Block coefficient', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1, group: 'Vessel', role: 'geometry' },
    { key: 'GM', label: 'Transverse metacentric height', unit: 'm', type: 'number', default: 6, min: 0.2, max: 60, group: 'Vessel', role: 'ic', bind: 'install.bargeGM', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Sets the roll natural period; from the stability module for the operating condition.' },
    { key: 'kxxB', label: 'Roll radius of gyration / beam (incl. added inertia)', unit: '-', type: 'number', default: 0.4, min: 0.25, max: 0.6, group: 'Vessel', role: 'geometry', uncertain: { dist: 'normal', cov: 0.08 } },
    { key: 'zetaRoll', label: 'Roll damping ratio', unit: '-', type: 'number', default: 0.08, min: 0.01, max: 0.4, group: 'Vessel', role: 'load', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Damping in addition to the wave-radiation damping that the model computes: about 0.03–0.10 of critical for barges with sharp bilges or bilge keels; governs resonant roll.' },
    ...seaInputs({ Hs: 1.5, Tp: 8 }),
    { key: 'gamma', label: 'JONSWAP peak enhancement', unit: '-', type: 'number', default: 3.3, min: 1, max: 7, group: 'Sea state', role: 'environment' },
    { key: 'heading', label: 'Wave heading', type: 'select', default: 135, options: HEADINGS, group: 'Sea state', role: 'bc', help: 'Direction the waves travel relative to the vessel: boundary condition of the response.' },
    { key: 'hours', label: 'Exposure duration for the extreme', unit: 'h', type: 'number', default: 3, min: 0.1, max: 72, group: 'Sea state', role: 'operation' },
    { key: 'xP', label: 'Point x from midship (forward positive)', unit: 'm', type: 'number', default: -70, min: -250, max: 250, group: 'Crane tip / point of interest', role: 'geometry' },
    { key: 'yP', label: 'Point y from the centreline', unit: 'm', type: 'number', default: 45, min: -150, max: 150, group: 'Crane tip / point of interest', role: 'geometry', help: 'An over-the-side lift puts the crane tip well outboard, where roll contributes most.' },
    { key: 'vLim', label: 'Limit on crane-tip vertical velocity (significant amplitude)', unit: 'm/s', type: 'number', default: 0.6, min: 0.05, max: 5, group: 'Criteria', role: 'criteria', help: `Typical lift-off / set-down criterion; take it from the lift dynamic analysis. ${CONFIRM}` },
    { key: 'rollLim', label: 'Limit on roll (significant amplitude)', unit: 'deg', type: 'number', default: 2, min: 0.1, max: 30, group: 'Criteria', role: 'criteria', help: `Crane operating limit or cargo criterion. ${CONFIRM}` },
    { key: 'purpose', label: 'Publish results as', type: 'select', default: 'operation', group: 'Criteria', role: 'control', options: [{ value: 'operation', label: 'Operating sea state (crane-tip motion for lifting modules)' }, { value: 'transport', label: 'Transport design sea state (motion criteria for seafastening)' }], help: 'Chooses which downstream inputs are fed from this run.' },
    { key: 'raoMethod', label: 'Motion model of the box hull', type: 'select', default: 'panel', group: 'Numerics', role: 'control', options: [{ value: 'panel', label: 'Panel solution (frequency-dependent added mass and damping)' }, { value: 'closed-form', label: 'Closed form (Jensen, Mansour & Olsen 2004; one-degree-of-freedom roll)' }], help: 'The panel solution is the accurate one for wave periods near and below the natural periods; the closed form is kept for hand checks and under-predicts heave below twice 2π√(2T/g).' },
    { key: 'nOmega', label: 'Frequency points', unit: '-', type: 'int', default: 200, min: 40, max: 2000, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nOmega', label: 'Frequency points', levels: [50, 100, 200, 400], qoi: 'tipSig' },
  compute(inp) {
    const head = Number(inp.heading), dur = inp.hours * 3600;
    const om = linspace(0.2, 2.2, inp.nOmega);
    const rao = (h) => run(() => boxBargeRAO({ L: inp.L, B: inp.B, T: inp.T, omega: om, headingDeg: h, Cb: inp.Cb, GMt: inp.GM, kxx: inp.kxxB * inp.B, rollDampingRatio: inp.zetaRoll, method: inp.raoMethod }), 'Vessel RAO');
    const r = rao(head);
    const roll = r.roll || new Float64Array(om.length);
    // Vertical motion RAO at the point: heave − x·pitch + y·roll, phases unknown → root-sum-square of the three contributions.
    const tipAmp = om.map((o, i) => Math.hypot(r.heave[i], inp.xP * r.pitch[i], inp.yP * roll[i]));
    const stats = (amp, Hs, Tp) => responseStatistics({ omega: om, amp, S: waveSpectrum(om, { Hs, Tp, gamma: inp.gamma }), duration: dur });
    const velAmp = tipAmp.map((a, i) => a * om[i]), accAmp = tipAmp.map((a, i) => a * om[i] * om[i]);
    const at = (Hs, Tp) => ({ heave: stats(r.heave, Hs, Tp), pitch: stats(r.pitch, Hs, Tp), roll: stats(roll, Hs, Tp), tip: stats(tipAmp, Hs, Tp), vel: stats(velAmp, Hs, Tp), acc: stats(accAmp, Hs, Tp) });
    need(inp.Hs > 0, 'Enter a significant wave height above zero.');
    const s = at(inp.Hs, inp.Tp);
    const range = r.validPeriodMin;
    const rollSig = s.roll.significant / DEG, pitchSig = s.pitch.significant / DEG, rollMpm = s.roll.mpm / DEG, pitchMpm = s.pitch.mpm / DEG;
    // Linear response: limiting Hs scales with the ratio of limit to response per metre of Hs.
    const tps = linspace(4, 16, 25);
    const limAt = (Tp) => { const q = at(1, Tp); return Math.min(inp.vLim / Math.max(q.vel.significant, 1e-9), inp.rollLim * DEG / Math.max(q.roll.significant, 1e-9)); };
    const lims = tps.map(limAt), hsLim = limAt(inp.Tp), hsLimMin = Math.min(...lims);
    const checks = [
      check('Crane-tip vertical velocity', s.vel.significant, inp.vLim, 'm/s', 'DNV-RP-N103 (2017) §4.3.3 and §9 (operational criteria on crane-tip motion)', 'Velocity of the crane tip sets lift-off snatch, slamming velocity and landing velocity. Above the limit the lift dynamic factors assumed in design are exceeded.'),
      check('Roll', rollSig, inp.rollLim, 'deg', 'Crane operating manual / DNV-ST-N001 Sec. 16 (operational limits)', 'Roll swings the suspended load sideways and side-loads the crane boom; cranes have a certified limit.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const Troll = r.rollNaturalPeriod ?? 0, Theave = r.heaveNaturalPeriod;
    const nearRes = Troll > 0 && Math.abs(inp.Tp - Troll) / Troll < 0.2;
    const [op, ah, ap, ar, atp] = thin([Array.from(om).map((o) => 2 * Math.PI / o), Array.from(r.heave), Array.from(r.pitch).map((v) => v / DEG), Array.from(roll).map((v) => v / DEG), tipAmp], 200);
    const pub = { 'install.craneTipHeave': s.tip.significant, 'install.craneTipVelocity': s.vel.significant, 'install.craneTipPeriod': s.tip.Tz, 'install.craneTipHeavePerHs': s.tip.significant / inp.Hs, 'install.HsOperation': inp.Hs };
    const pubT = { 'install.transportRollDeg': rollMpm, 'install.transportRollT': Troll > 0 ? Troll : inp.Tp, 'install.transportPitchDeg': pitchMpm, 'install.transportHeaveG': s.heave.mpm * (2 * Math.PI / Math.max(s.heave.Tz, 1)) ** 2 / G };
    return {
      kpis: [
        kpi('tipSig', 'Crane-tip vertical motion, significant amplitude', s.tip.significant, 'm', 'Two standard deviations of the vertical motion at the point of interest; the characteristic single amplitude for lifting calculations.'),
        kpi('tipMpm', 'Crane-tip vertical motion, most probable maximum', s.tip.mpm, 'm', 'Expected largest single amplitude in the exposure duration (Rayleigh maxima).'),
        kpi('tipVel', 'Crane-tip vertical velocity, significant amplitude', s.vel.significant, 'm/s', 'Feeds lift-off DAF, slamming velocity and landing velocity.', inp.vLim, 'max'),
        kpi('tipAcc', 'Crane-tip vertical acceleration, significant amplitude', s.acc.significant, 'm/s²', 'Feeds the inertia term of the splash-zone force.'),
        kpi('tipTz', 'Mean period of the crane-tip motion', s.tip.Tz, 's', 'Zero-crossing period of the vertical motion; compare with the natural period of the hoisting system and of the lowered object.'),
        kpi('heaveSig', 'Heave, significant amplitude', s.heave.significant, 'm', 'Vertical motion of the centre of the vessel.'),
        kpi('rollSig', 'Roll, significant amplitude', rollSig, 'deg', 'Zero in a pure head sea. Resonant when the wave period is near the roll natural period.', inp.rollLim, 'max'),
        kpi('pitchSig', 'Pitch, significant amplitude', pitchSig, 'deg', 'Largest in head seas with wavelengths near the vessel length.'),
        kpi('rollMpm', 'Roll, most probable maximum', rollMpm, 'deg', 'Expected extreme roll amplitude in the exposure duration; the value used as a transport motion criterion when this run is the transport design sea state.'),
        kpi('Troll', 'Roll natural period', Troll, 's', '2π·kxx/√(g·GM). Keep it away from the dominant wave periods of the operation.'),
        kpi('Theave', 'Heave natural period', Theave, 's', 'From the hydrostatic stiffness and mass with added mass of the box hull.'),
        kpi('hsLim', 'Limiting Hs at the given period and heading', hsLim, 'm', 'Significant wave height at which the first motion criterion is reached (linear scaling). Design limit; apply the alpha-factor for the forecast limit.'),
        kpi('hsLimMin', 'Limiting Hs over peak periods 4–16 s', hsLimMin, 'm', 'Lowest limit over the period range at this heading.'),
      ],
      checks,
      plots: [
        line('Response amplitude operators', 'Wave period', 's', 'Heave, point (m/m) · pitch, roll (deg/m)', '', [{ name: 'Heave (m/m)', x: op, y: ah }, { name: 'Pitch (deg/m)', x: op, y: ap }, { name: 'Roll (deg/m)', x: op, y: ar }, { name: 'Vertical at the point (m/m)', x: op, y: atp }],
          'Motion per metre of wave amplitude against wave period. The roll peak marks the roll natural period; the point curve shows how much pitch and roll add at the crane tip.'),
        line('Limiting sea state against peak period', 'Peak period Tp', 's', 'Limiting Hs', 'm', [{ name: 'Limiting Hs', x: tps.map((t) => round(t)), y: lims.map((v) => round(Math.min(v, 15))) }, { name: 'Operating point', x: [inp.Tp], y: [inp.Hs] }],
          'Sea states below the curve satisfy both motion criteria at this heading. The dip is where the wave period meets the roll or heave resonance: avoid planning the operation there.'),
      ],
      tables: [tbl('Short-term response statistics', [['q', 'Response'], ['sig', 'Significant amplitude'], ['mpm', 'Most probable maximum'], ['tz', 'Mean period', 's'], ['u', 'Unit']], [
        { q: 'Heave', sig: round(s.heave.significant), mpm: round(s.heave.mpm), tz: round(s.heave.Tz), u: 'm' },
        { q: 'Pitch', sig: round(pitchSig), mpm: round(pitchMpm), tz: round(s.pitch.Tz), u: 'deg' },
        { q: 'Roll', sig: round(rollSig), mpm: round(rollMpm), tz: round(s.roll.Tz || 0), u: 'deg' },
        { q: 'Vertical motion at the point', sig: round(s.tip.significant), mpm: round(s.tip.mpm), tz: round(s.tip.Tz), u: 'm' },
        { q: 'Vertical velocity at the point', sig: round(s.vel.significant), mpm: round(s.vel.mpm), tz: round(s.vel.Tz), u: 'm/s' },
        { q: 'Vertical acceleration at the point', sig: round(s.acc.significant), mpm: round(s.acc.mpm), tz: round(s.acc.Tz), u: 'm/s²' },
      ], 'Significant (two standard deviations) and expected extreme single amplitudes of each response in the sea state, for use as inputs to the lifting, lowering and transport modules.')],
      notes: [`Heading ${head}°, Hs ${inp.Hs} m, Tp ${inp.Tp} s: crane-tip motion ${fmt(s.tip.significant)} m, velocity ${fmt(s.vel.significant)} m/s (significant amplitudes).`, nearRes ? `The peak period is within 20 % of the roll natural period (${fmt(Troll)} s): resonant roll.` : `Roll natural period ${fmt(Troll)} s, heave natural period ${fmt(Theave)} s.`, `"${gov.name}" governs at ${fmt(gov.uc)}; limiting Hs ${fmt(hsLim)} m.`],
      advice: [
        ...(range && r.method === 'panel' && inp.Tp < Math.max(range.heave, range.roll || 0) ? [adv('engineering', 'warn',
          `Peak period ${fmt(inp.Tp)} s is below the range of the box panel model, which is carried down to the wave whose length equals the beam (period ${fmt(Math.max(range.heave, range.roll || 0))} s). Much of the energy of this sea lies in shorter waves, where the model uses the coefficients of that last period: the motions there are small but not verified.`,
          `Use the vessel's diffraction RAOs for sea states with Tp below ${fmt(Math.max(range.heave, range.roll || 0))} s before setting an operational limit.`,
          'DNV-RP-C205 (2019) §7.3 (wave-frequency motions by radiation–diffraction analysis); comparison of the box panel model with Capytaine 3.0.0 on barges of 20 × 8 × 2 m, 100 × 20 × 5 m and 80 × 24 × 3 m (cross-solver verification of this application).')] : []),
        ...(range && r.method !== 'panel' && inp.Tp < range.heave ? [adv('engineering', 'warn',
          `Peak period ${fmt(inp.Tp)} s is below twice the heave natural period of the closed-form model (${fmt(range.heave / 2)} s, limit ${fmt(range.heave)} s). There its heave is low against a panel solution, about 15 % at 1.5 times the natural period and 30–40 % near it, and pitch 20–60 % low near it: crane-tip motion and the limiting sea state are on the unsafe side.`,
          `Switch the motion model to the panel solution, or use the vessel's diffraction RAOs for sea states with Tp below ${fmt(range.heave)} s.`,
          'DNV-RP-C205 (2019) §7.3 (range of strip and closed-form motion estimates); comparison of the closed-form box model with Capytaine 3.0.0 on a 20 × 8 × 2 m barge (cross-solver verification of this application).')] : []),
        ...(range && r.method !== 'panel' && range.roll && inp.Tp < range.roll && Math.abs(Math.sin(head * DEG)) > 0.1 ? [adv('engineering', 'warn',
          `Peak period ${fmt(inp.Tp)} s is below 1.3 times the roll natural period (${fmt(range.roll / 1.3)} s, limit ${fmt(range.roll)} s). There the one-degree-of-freedom roll was within 15 % of a panel solution instead of 5 % (compared down to 0.75 times the natural period), and near resonance it follows the roll gyradius and the added damping that were entered.`,
          'Confirm the roll gyradius and damping from the stability booklet and a roll-decay record, and vary both before relying on the resonant roll.',
          'DNV-RP-C205 (2019) §7.3 (range of strip and closed-form motion estimates); comparison of the closed-form box model with Capytaine 3.0.0 on a 20 × 8 × 2 m barge (cross-solver verification of this application).')] : []),
        ok ? adv('engineering', gov.uc > 0.9 || nearRes ? 'warn' : 'ok', `Motions are within the criteria: crane-tip velocity ${fmt(s.vel.significant)} m/s against ${inp.vLim} m/s, roll ${fmt(rollSig)}° against ${inp.rollLim}°. The design limit at this period and heading is Hs ${fmt(hsLim)} m.`,
          `Check the other headings and the period range (lowest limit ${fmt(hsLimMin)} m over 4–16 s), and replace these box-hull RAOs with the diffraction-analysis RAOs of the actual hull before fixing the operational limits.`, 'DNV-RP-N103 (2017) §4.3.3, §9; DNV-RP-C205 (2019) §7')
          : adv('engineering', 'fail', `"${gov.name}" is exceeded (${fmt(gov.uc)}): crane-tip velocity ${fmt(s.vel.significant)} m/s, roll ${fmt(rollSig)}°. The limit at this period and heading is Hs ${fmt(hsLim)} m.`,
            `Work below Hs ${fmt(hsLim)} m, turn the vessel towards a head sea to remove roll, bring the lift inboard (y = ${inp.yP} m contributes ${fmt(inp.yP * s.roll.significant)} m), or change GM so the roll period of ${fmt(Troll)} s moves away from Tp.`, 'DNV-RP-N103 (2017) §4.3.3, §9; DNV-RP-C205 (2019) §7'),
        adv('economic', 'info', `The motion-limited sea state is Hs ${fmt(hsLim)} m at Tp ${inp.Tp} s and ${fmt(hsLimMin)} m in the worst period band.`,
          'Operability at the site follows directly from this limit: run the weather-window module with it. Heading control (DP or mooring spread) and a longer roll period are usually cheaper ways to gain workable days than a larger vessel.',
          'Vessel selection by operability: DNV-ST-N001 Sec. 2 (operational limits and alpha-factor); day rate × expected waiting on weather',
          `Each day of weather standby avoided saves about ${vesselCo2Text(1, 40)}`),
      ],
      publish: inp.purpose === 'transport' ? pubT : pub,
      quantities: {},
    };
  },
  verify: [
    { name: 'Roll natural period: T = 2π·kxx/√(g·GM)', source: 'One-degree-of-freedom roll (Biran & López-Pulido ch. 6): kxx = 0.4·42 = 16.8 m, GM = 6 m → 2π·16.8/√(9.80665·6) = 13.762 s',
      inputs: {}, expect: [{ kpi: 'Troll', value: 2 * Math.PI * 16.8 / Math.sqrt(9.80665 * 6), tol: 1e-6 }] },
    { name: 'Closed-form model: heave natural period of a box, T = 2π·√(2T/g)', source: 'Jensen, Mansour & Olsen (2004) heave equation (sectional added mass equal to the displaced mass): T = 8 m → 2π·√(16/9.80665) = 8.026 s',
      inputs: { raoMethod: 'closed-form' }, expect: [{ kpi: 'Theave', value: 2 * Math.PI * Math.sqrt(16 / 9.80665), tol: 1e-6 }] },
    { name: 'Long-wave limit: the vessel follows the wave, heave significant amplitude → Hs/2', source: 'Hydrostatic limit of a small body in long waves: RAO → 1, so the significant heave amplitude equals 2·√m0 = Hs/2; 30 m barge in Tp = 20 s swell, head sea (within 3 %)',
      inputs: { L: 30, B: 10, T: 2, GM: 3, Hs: 2, Tp: 20, gamma: 5, heading: 180, xP: 0, yP: 0, vLim: 5, rollLim: 30 }, expect: [{ kpi: 'heaveSig', value: 1.0, tol: 0.03 }] },
  ],
};
