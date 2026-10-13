// Six-degree-of-freedom rigid-body motions: natural periods, regular-wave response, irregular-sea statistics for
// every degree of freedom and for a point on the unit, and the variation with wave heading.
import { pointMotion } from '../../engine/hydro/floater.js';
import { hullInputs, hullFrom, motionInputs, motionFrom, motionRangeAdvice, seaInputs, seaFrom, spectrum, omegaGrid, stats, pointInputs, DOF, DOF_UNIT, dofScale,
  G, LIMIT_HYDRO, LIMIT_CFD, round, rnd, fmt, adv, check, sev, worst, fuelCo2Text } from './_common.js';

const at = (xs, ys, x) => { for (let i = 1; i < xs.length; i++) if (x <= xs[i]) return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[ys.length - 1]; };

export default {
  id: 's11-motions-6dof',
  suite: 11,
  title: 'Six-degree-of-freedom vessel motions',
  summary: 'Natural periods and surge, sway, heave, roll, pitch and yaw response of the unit in a regular wave and in an irregular sea, with motions, velocities and accelerations at a chosen point and their variation with heading.',
  covers: [3],
  fidelity: 'reduced-order',
  standards: ['DNV-RP-C205 (2019) §7', 'DNV-RP-C103 (2015) §3', 'ISO 19904-1 (2019) §9', 'NORDFORSK (1987) seakeeping criteria', 'CAP 437 (helideck motion limits)'],
  models: ['Newton–Euler six-degree-of-freedom equations of motion (linearised, frequency domain)', 'Strip added mass and hydrostatic restoring', 'Wave excitation by long-wave Froude–Krylov and Morison inertia', 'Stochastic linearisation of viscous drag', 'Rigid-body point kinematics', 'Rayleigh short-term extremes'],
  limits: `${LIMIT_HYDRO} Responses are linear and Gaussian; the box model returns heave, roll and pitch only, so point motions of ship-shaped hulls omit surge, sway and yaw. Second-order slow-drift motions are treated in the coupled station-keeping module. ${LIMIT_CFD}`,
  publishes: ['vessel.heaveSig', 'vessel.rollSigDeg', 'vessel.Tn3', 'vessel.Tn4'],
  searchTerms: 'floating platform six degree of freedom motions natural periods',
  inputs: [
    ...hullInputs({ type: 'semi', withMesh: false }),
    ...seaInputs({ Hs: 6, Tp: 11, ret: '1', heading: 30 }),
    { key: 'Hreg', label: 'Regular wave height', unit: 'm', type: 'number', default: 5, min: 0, max: 40, group: 'Regular wave', role: 'environment', help: 'Height of the regular wave for the harmonic response table.' },
    { key: 'Treg', label: 'Regular wave period', unit: 's', type: 'number', default: 10, min: 3, max: 30, group: 'Regular wave', role: 'environment', help: 'Period of the regular wave.' },
    ...motionInputs(),
    ...pointInputs(),
    { key: 'limHeave', label: 'Allowable significant heave amplitude', unit: 'm', type: 'number', default: 2, min: 0.1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: 'Operating limit on the significant single amplitude of heave (drilling, crane or gangway operations).' },
    { key: 'limAngle', label: 'Allowable significant roll / pitch amplitude', unit: 'deg', type: 'number', default: 3, min: 0.2, max: 20, group: 'Acceptance criteria', role: 'criteria', help: 'Operating limit on significant single-amplitude roll and pitch (CAP 437 helideck limits are of this order).' },
    { key: 'limAcc', label: 'Allowable significant vertical acceleration at the point', unit: 'm/s²', type: 'number', default: 1.0, min: 0.1, max: 10, group: 'Acceptance criteria', role: 'criteria', help: 'NORDFORSK-type limit; 0.1 g significant amplitude is a common comfort / equipment value.' },
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [50, 100, 200], qoi: 'heaveSig' },
  compute(inp) {
    const hull = hullFrom(inp); const mot = motionFrom(inp, hull);
    const sea = seaFrom(inp); const om = omegaGrid(inp.nFreq); const S = spectrum(om, sea); const dur = inp.durationH * 3600;
    const r = mot.rao(om, inp.headingDeg, sea);
    const st = r.amp.map((a) => stats(om, a, S, dur));
    const pt = [inp.px, inp.py, inp.pz - hull.T];
    const pm = pointMotion({ rao: r, point: pt });
    const pd = pm.amp.map((a) => stats(om, a, S, dur)); const pv = pm.velAmp.map((a) => stats(om, a, S, dur)); const pa = pm.accAmp.map((a) => stats(om, a, S, dur));
    // regular wave: harmonic amplitudes at the wave frequency (drag damping linearised for the sea state is kept)
    const wr = 2 * Math.PI / inp.Treg;
    const reg = r.amp.map((a) => at(om, a, wr) * inp.Hreg / 2);
    // heading sweep with the linearised drag damping of the design sea
    const heads = Array.from({ length: 13 }, (_, i) => i * 30);
    const sweep = heads.map((h) => { const rr = mot.raoLin(om, h, r.viscousDamping); return [2, 3, 4].map((d) => stats(om, rr.amp[d], S, dur).significant * dofScale(d)); });
    const worstHead = (j) => heads[sweep.reduce((b, v, i) => (v[j] > sweep[b][j] + 1e-12 ? i : b), 0)] % 360;
    const rot = Math.max(st[3].significant, st[4].significant) * dofScale(3);
    const checks = [
      check('Significant heave amplitude', st[2].significant, inp.limHeave, 'm', 'Operational limiting criterion (DNV-ST-N001 (2021) §2.6; project value)',
        'Guards the operation, not the structure: above the limit the heave compensator, crane or gangway runs out of stroke and work must stop.'),
      check('Significant roll / pitch amplitude', rot, inp.limAngle, 'deg', 'CAP 437 helideck motion limits; NORDFORSK (1987)',
        'Guards helicopter operations, personnel safety and process separation efficiency, which all degrade with platform inclination.'),
      check('Significant vertical acceleration at the point', pa[2].significant, inp.limAcc, 'm/s²', 'NORDFORSK (1987) seakeeping criteria; ISO 19901-1 operability practice',
        'Guards equipment foundations and people at the point of interest: inertia loads on cargo and fatigue of supports scale directly with this acceleration.'),
    ];
    const gov = worst(checks);
    const Tn = mot.Tn; const T = om.map((w) => round(2 * Math.PI / w, 4));
    const Tz = st[2].Tz || inp.Tp / 1.3;
    const fuelDay = 12; // screening support-vessel standby burn, t/day
    return {
      kpis: [
        ...[0, 1, 2, 3, 4, 5].filter((d) => Tn[d] != null).map((d) => ({ key: `Tn${d + 1}`, label: `${DOF[d]} natural period`, value: round(Tn[d], 5), unit: 's',
          explain: d < 2 || d === 5 ? `Natural period of ${DOF[d].toLowerCase()} on the station-keeping stiffness. It lies far above wave periods, so the mode is excited by slowly varying wind and drift forces, not by first-order waves.`
            : `Undamped natural period in ${DOF[d].toLowerCase()}. Keep it outside the 5–20 s band of energetic waves (above about 20 s for a semi-submersible or spar, below about 4 s for a TLP).` })),
        ...[0, 1, 2, 3, 4, 5].map((d) => ({ key: `${DOF[d].toLowerCase()}Sig`, label: `${DOF[d]} significant amplitude`, value: round(st[d].significant * dofScale(d), 4), unit: DOF_UNIT[d],
          ...(d === 2 ? { limit: inp.limHeave, sense: 'max' } : d === 3 || d === 4 ? { limit: inp.limAngle, sense: 'max' } : {}),
          explain: `Twice the standard deviation of wave-frequency ${DOF[d].toLowerCase()} in the sea state: the average of the highest third of single amplitudes.` })),
        { key: 'heaveMpm', label: 'Most probable maximum heave', value: round(st[2].mpm, 4), unit: 'm', explain: 'Largest heave amplitude expected in the sea-state duration, σ·√(2 ln N) for Rayleigh-distributed maxima.' },
        { key: 'rollMpm', label: 'Most probable maximum roll', value: round(st[3].mpm * dofScale(3), 4), unit: 'deg', explain: 'Largest roll amplitude expected in the sea-state duration.' },
        { key: 'pitchMpm', label: 'Most probable maximum pitch', value: round(st[4].mpm * dofScale(4), 4), unit: 'deg', explain: 'Largest pitch amplitude expected in the sea-state duration.' },
        { key: 'ptDispSig', label: 'Vertical motion at the point, significant', value: round(pd[2].significant, 4), unit: 'm', explain: 'Significant vertical displacement amplitude of the point: heave combined in phase with the roll and pitch levers.' },
        { key: 'ptVelSig', label: 'Vertical velocity at the point, significant', value: round(pv[2].significant, 4), unit: 'm/s', explain: 'Significant vertical velocity amplitude at the point; the quantity that limits load landing and helicopter touchdown.' },
        { key: 'ptAccSig', label: 'Vertical acceleration at the point, significant', value: round(pa[2].significant, 4), unit: 'm/s²', limit: inp.limAcc, sense: 'max', explain: 'Significant vertical acceleration amplitude at the point; multiply by equipment mass for the inertia load.' },
        { key: 'ptAccHSig', label: 'Horizontal acceleration at the point, significant', value: round(Math.hypot(pa[0].significant, pa[1].significant), 4), unit: 'm/s²', explain: 'Significant horizontal acceleration amplitude from surge, sway and the rotation levers. The gravity component g·sin(angle) is not included and must be added for sea-fastening loads.' },
        { key: 'heaveReg', label: 'Heave amplitude in the regular wave', value: round(reg[2], 4), unit: 'm', explain: 'Harmonic heave amplitude in the regular wave of the given height and period: RAO times wave amplitude.' },
        { key: 'pitchReg', label: 'Pitch amplitude in the regular wave', value: round(reg[4] * dofScale(4), 4), unit: 'deg', explain: 'Harmonic pitch amplitude in the regular wave.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Significant motion amplitudes', categories: DOF.map((n, d) => `${n} (${DOF_UNIT[d]})`), y: { label: 'Significant amplitude', unit: 'm or deg' },
          series: [{ name: 'Significant', y: st.map((s, d) => round(s.significant * dofScale(d), 4)) }, { name: 'Most probable maximum', y: st.map((s, d) => round(s.mpm * dofScale(d), 4)) }],
          explain: 'All six degrees of freedom side by side for this heading. The ratio of maximum to significant is about 1.9 for a 3-hour sea state; the dominant bars show which motion limits operations.' },
        { type: 'line', title: 'Response spectra of heave, roll and pitch', x: { label: 'Wave period', unit: 's' }, y: { label: 'Normalised spectral density', unit: '-' },
          series: [2, 3, 4].map((d) => { const y = r.amp[d].map((a, i) => a * a * S[i]); const m = Math.max(...y, 1e-300); return { name: DOF[d], x: T, y: rnd(y.map((v) => v / m)) }; }).concat([{ name: 'Wave spectrum', x: T, y: rnd(Array.from(S, (v) => v / Math.max(...S, 1e-300))), dash: true }]),
          explain: 'Each response spectrum scaled to a peak of one, with the wave spectrum dashed. A response peak away from the wave peak is resonance at a natural period; a peak on the wave peak is direct wave forcing.' },
        { type: 'polar', title: 'Significant heave, roll and pitch against wave heading', series: [{ name: 'Heave (m)', thetaDeg: heads, r: sweep.map((v) => round(v[0], 4)) }, { name: 'Roll (deg)', thetaDeg: heads, r: sweep.map((v) => round(v[1], 4)) }, { name: 'Pitch (deg)', thetaDeg: heads, r: sweep.map((v) => round(v[2], 4)) }],
          explain: 'Heading rosette: radius is the significant amplitude for waves travelling towards that direction. Use it to choose the unit orientation or the weather-vaning heading with the least motion.' },
      ],
      tables: [
        { title: 'Motion statistics by degree of freedom', columns: [{ key: 'dof', label: 'DOF' }, { key: 'tn', label: 'Natural period', unit: 's' }, { key: 'reg', label: 'Regular-wave amplitude' }, { key: 'sig', label: 'Significant' }, { key: 'mpm', label: 'Most probable max' }, { key: 'tz', label: 'Tz', unit: 's' }, { key: 'u', label: 'Unit' }],
          rows: DOF.map((n, d) => ({ dof: n, tn: Tn[d] != null ? round(Tn[d], 4) : '–', reg: round(reg[d] * dofScale(d), 4), sig: round(st[d].significant * dofScale(d), 4), mpm: round(st[d].mpm * dofScale(d), 4), tz: round(st[d].Tz, 4), u: DOF_UNIT[d] })),
          explain: 'One row per degree of freedom about the waterline origin: natural period, harmonic amplitude in the regular wave, and irregular-sea statistics as single amplitudes.' },
        { title: 'Motion at the point of interest', columns: [{ key: 'dir', label: 'Direction' }, { key: 'd', label: 'Displacement sig.', unit: 'm' }, { key: 'v', label: 'Velocity sig.', unit: 'm/s' }, { key: 'a', label: 'Acceleration sig.', unit: 'm/s²' }, { key: 'am', label: 'Acceleration max', unit: 'm/s²' }],
          rows: ['x', 'y', 'z'].map((n, i) => ({ dir: n, d: round(pd[i].significant, 4), v: round(pv[i].significant, 4), a: round(pa[i].significant, 4), am: round(pa[i].mpm, 4) })),
          explain: 'Displacement, velocity and acceleration of the chosen point in the three directions (significant single amplitudes, plus the most probable maximum acceleration).' },
      ],
      notes: [
        `${mot.kind === 'slender' ? `Slender-member model, valid for wave periods above ${fmt(2 * Math.PI / mot.validOmegaMax)} s` : 'Box-barge model (heave, roll, pitch only)'}; drag linearisation ${r.converged ? `converged in ${r.iterations} iterations` : 'did NOT converge — treat the damping as approximate'}.`,
        `Natural periods: heave ${fmt(Tn[2])} s, roll ${fmt(Tn[3])} s, pitch ${fmt(Tn[4])} s; sea state Tp = ${fmt(inp.Tp)} s.`,
        `Worst headings: heave at ${worstHead(0)}°, roll at ${worstHead(1)}°, pitch at ${worstHead(2)}°. Governing criterion: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        ...motionRangeAdvice(mot, Math.min(inp.Tp, inp.Treg), inp.headingDeg),
        adv('engineering', sev(gov.uc), `${gov.name} is ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit} (utilisation ${fmt(gov.uc)}) in Hs ${fmt(inp.Hs)} m, Tp ${fmt(inp.Tp)} s.`,
          gov.uc > 1 ? `Reduce the response at source: move the natural period away from ${fmt(inp.Tp)} s (heave ${fmt(Tn[2])} s, roll ${fmt(Tn[3])} s), turn the unit towards the best heading in the rosette, or lower the operating sea-state limit to about Hs ${fmt(inp.Hs / gov.uc)} m.`
            : `The limit is reached at about Hs ${fmt(inp.Hs / Math.max(gov.uc, 1e-6))} m for this period. Confirm the damping fraction with a decay test or diffraction analysis: resonant response scales with 1/damping.`,
          'DNV-RP-C205 (2019) §7 (floater motions); DNV-ST-N001 (2021) for operational limits; CAP 437 for helideck motion.'),
        adv('economic', 'info', `Motion limits translate into weather downtime: every hour above the limit idles the spread. With a zero-crossing period of ${fmt(Tz)} s the unit sees about ${fmt(3600 / Tz)} motion cycles per hour.`,
          'Take the limiting Hs from this module into the frequency-domain operability module to convert it to downtime hours per year before investing in motion-compensation equipment.',
          'Operability and weather-window practice of DNV-ST-N001 (2021) §2.6; whole-life costing per ISO 15686-5.',
          `Each day a supply or support vessel waits on weather burns about ${fuelDay} t of fuel: ${fuelCo2Text(fuelDay)} Better motion behaviour cuts those standby days.`),
      ],
      publish: { 'vessel.heaveSig': st[2].significant, 'vessel.rollSigDeg': st[3].significant * dofScale(3), 'vessel.Tn3': Tn[2] ?? 0, 'vessel.Tn4': Tn[3] ?? 0 },
    };
  },
  verify: [
    { // Spar D = 30 m, T = 90 m: m = ρ·A·T, heave added mass = (2/3)·π·ρ·r³ (hemisphere of water below the keel), C33 = ρ·g·A.
      // Added mass over ρ·A = (2/3)·r = 10 m, so Tn = 2π·sqrt((T + 2r/3)/g) = 2π·sqrt(100/9.80665) = 2π·sqrt(10.1972) = 20.064 s.
      name: 'Spar heave natural period T = 2π √((m + a33)/(ρ g Aw))',
      source: 'Faltinsen (1990) ch. 3; hemisphere end added mass (2/3)πρr³ of a deep column — hand calculation',
      inputs: { hullType: 'spar', Dspar: 30, Hspar: 105, draft: 90, KG: 38 },
      expect: [{ kpi: 'Tn3', value: 20.064, tol: 0.003 }] },
    { // A 30 × 10 × 1 m barge in a 20 s wave (λ = 624 m ≫ L) follows the surface: heave amplitude = H/2 = 1.0 m.
      name: 'Small barge in a long regular wave follows the water surface',
      source: 'Long-wave limit of linear seakeeping (heave RAO → 1), Newman (1977) Marine Hydrodynamics ch. 6',
      inputs: { hullType: 'barge', L: 30, B: 10, Dhull: 4, draft: 1, KG: 1.5, Hreg: 2, Treg: 20, headingDeg: 180, Hs: 2, Tp: 12 },
      expect: [{ kpi: 'heaveReg', value: 1.0, tol: 0.02 }] },
    { name: 'Flat calm gives no motion', source: 'Limiting case Hs = 0 of the spectral response', inputs: { Hs: 0 },
      expect: [{ kpi: 'heaveSig', value: 0, abs: 1e-9 }, { kpi: 'ptAccSig', value: 0, abs: 1e-9 }] },
  ],
};
