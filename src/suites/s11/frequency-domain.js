// Frequency-domain spectral analysis: response spectra and statistics in the design sea, motion-induced loads,
// limiting sea states and operability over a wave scatter diagram.
import { pointMotion, operability, limitingHs } from '../../engine/hydro/floater.js';
import { hullInputs, hullFrom, motionInputs, motionFrom, motionRangeAdvice, seaInputs, seaFrom, spectrum, omegaGrid, stats, pointInputs, dofScale,
  G, LIMIT_HYDRO, LIMIT_CFD, need, round, rnd, fmt, adv, checkMin, check, sev, worst, fuelCo2Text } from './_common.js';

export const SCATTER = [
  { Hs: 0.5, Tp: 6, prob: 0.08 }, { Hs: 1, Tp: 7, prob: 0.14 }, { Hs: 1.5, Tp: 8, prob: 0.15 }, { Hs: 2, Tp: 8, prob: 0.12 }, { Hs: 2, Tp: 10, prob: 0.08 }, { Hs: 2.5, Tp: 9, prob: 0.09 },
  { Hs: 3, Tp: 10, prob: 0.08 }, { Hs: 3.5, Tp: 10, prob: 0.06 }, { Hs: 4, Tp: 11, prob: 0.055 }, { Hs: 4.5, Tp: 12, prob: 0.04 }, { Hs: 5, Tp: 12, prob: 0.03 }, { Hs: 6, Tp: 13, prob: 0.03 },
  { Hs: 7, Tp: 13, prob: 0.018 }, { Hs: 8, Tp: 14, prob: 0.012 }, { Hs: 9, Tp: 15, prob: 0.008 }, { Hs: 11, Tp: 16, prob: 0.004 }, { Hs: 13, Tp: 17, prob: 0.002 }, { Hs: 15, Tp: 18, prob: 0.001 },
];
export const scatterInput = (group = 'Wave scatter diagram', rows = SCATTER) => ({ key: 'scatter', label: 'Sea states (Hs, Tp, probability)', type: 'table', group, role: 'environment',
  columns: [{ key: 'Hs', label: 'Hs', unit: 'm', type: 'number' }, { key: 'Tp', label: 'Tp', unit: 's', type: 'number' }, { key: 'prob', label: 'Probability', unit: '-', type: 'number' }],
  default: rows, help: 'Long-term distribution of sea states at the site, one row per Hs–Tp cell; probabilities are normalised automatically.' });
export function scatterFrom(inp) {
  const rows = (inp.scatter || []).filter((r) => r.Hs > 0 && r.Tp > 0 && r.prob > 0);
  need(rows.length > 0, 'The scatter diagram is empty: add at least one sea state with positive Hs, Tp and probability.');
  const tot = rows.reduce((a, r) => a + r.prob, 0);
  return rows.map((r) => ({ Hs: r.Hs, Tp: r.Tp, prob: r.prob / tot }));
}

export default {
  id: 's11-frequency-domain',
  suite: 11,
  title: 'Frequency-domain motion and load analysis',
  summary: 'Spectral response of the unit in the design sea, motion-induced inertia loads on equipment, limiting sea state against period, and operability and downtime over the site scatter diagram.',
  covers: [5],
  fidelity: 'reduced-order',
  standards: ['DNV-RP-C205 (2019) §3.5, §7', 'DNV-ST-N001 (2021) §2.6 (operational limiting criteria)', 'ISO 19901-1 (2015) (metocean, scatter diagrams)', 'NORDFORSK (1987)'],
  models: ['Wave spectral response S_R = |RAO|²·S (JONSWAP)', 'Spectral moments, bandwidth and Rayleigh extremes', 'Linear transfer functions for motion-induced inertia loads', 'Limiting significant wave height by linear scaling', 'Operability over a scatter diagram'],
  limits: `${LIMIT_HYDRO} Long-crested seas from one heading; the drag damping linearised for the design sea is reused for every cell of the scatter diagram (conservative below the design sea, slightly optimistic above). Responses are linear and Gaussian. ${LIMIT_CFD}`,
  publishes: ['install.operabilityFrac'],
  searchTerms: 'spectral response operability scatter diagram floating platform',
  inputs: [
    ...hullInputs({ type: 'semi', withMesh: false }),
    ...seaInputs({ Hs: 6, Tp: 12, ret: '1', heading: 30 }),
    scatterInput(),
    ...motionInputs(),
    ...pointInputs(),
    { key: 'equipMass', label: 'Equipment mass at the point', unit: 't', type: 'number', default: 400, min: 0, max: 50000, group: 'Point of interest', role: 'load', help: 'Mass of the module, crane load or flare boom whose motion-induced inertia load is wanted.' },
    { key: 'limHeave', label: 'Limit: significant heave amplitude', unit: 'm', type: 'number', default: 1.5, min: 0.1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: 'Operation stops above this significant single-amplitude heave.' },
    { key: 'limAngle', label: 'Limit: significant roll / pitch amplitude', unit: 'deg', type: 'number', default: 2, min: 0.2, max: 20, group: 'Acceptance criteria', role: 'criteria', help: 'Operation stops above this significant single-amplitude rotation.' },
    { key: 'limAcc', label: 'Limit: significant vertical acceleration at the point', unit: 'm/s²', type: 'number', default: 0.7, min: 0.05, max: 10, group: 'Acceptance criteria', role: 'criteria', help: 'Operation stops above this significant acceleration amplitude.' },
    { key: 'opReq', label: 'Required operability', unit: '%', type: 'number', default: 90, min: 10, max: 100, group: 'Acceptance criteria', role: 'criteria', help: 'Share of time the operation must be workable for the project economics.' },
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [50, 100, 200], qoi: 'heaveSig' },
  compute(inp) {
    const hull = hullFrom(inp); const mot = motionFrom(inp, hull);
    const sea = seaFrom(inp); const om = omegaGrid(inp.nFreq); const S = spectrum(om, sea); const dur = inp.durationH * 3600;
    const r = mot.rao(om, inp.headingDeg, sea);
    const pm = pointMotion({ rao: r, point: [inp.px, inp.py, inp.pz - hull.T] });
    const rotAmp = Float64Array.from(om, (_, i) => Math.max(r.amp[3][i], r.amp[4][i]));
    const h = stats(om, r.amp[2], S, dur); const ro = stats(om, r.amp[3], S, dur); const pi = stats(om, r.amp[4], S, dur); const ac = stats(om, pm.accAmp[2], S, dur);
    const m4 = (() => { let s = 0; for (let i = 1; i < om.length; i++) s += 0.5 * (r.amp[2][i] ** 2 * S[i] * om[i] ** 4 + r.amp[2][i - 1] ** 2 * S[i - 1] * om[i - 1] ** 4) * (om[i] - om[i - 1]); return s; })();
    const eps = h.m0 > 0 && m4 > 0 ? Math.sqrt(Math.max(0, 1 - h.m2 * h.m2 / (h.m0 * m4))) : 0;
    const Fz = inp.equipMass * 1e3 * ac.mpm; // inertia force amplitude on the equipment, N
    const sc = scatterFrom(inp);
    const lim = { heave: inp.limHeave, rotation: inp.limAngle, acceleration: inp.limAcc };
    const evaluate = (Hs, Tp) => { const Sx = spectrum(om, { Hs, Tp, gamma: sea.gamma }); return { heave: stats(om, r.amp[2], Sx).significant, rotation: Math.max(stats(om, r.amp[3], Sx).significant, stats(om, r.amp[4], Sx).significant) * dofScale(3), acceleration: stats(om, pm.accAmp[2], Sx).significant }; };
    const op = operability({ seaStates: sc, evaluate, limits: lim });
    const Tps = Array.from({ length: 15 }, (_, i) => 5 + i);
    const lh = [['Heave', r.amp[2], inp.limHeave], ['Roll / pitch', rotAmp.map((v) => v * dofScale(3)), inp.limAngle], ['Point acceleration', pm.accAmp[2], inp.limAcc]]
      .map(([name, amp, limit]) => ({ name, ...limitingHs({ omega: om, amp, limit, Tp: Tps, gamma: sea.gamma, statistic: 'significant' }) }));
    const env = Tps.map((_, i) => Math.min(...lh.map((l) => l.HsLimit[i])));
    const HsLimAtTp = (() => { const i = Math.max(0, Math.min(Tps.length - 1, Math.round(inp.Tp - 5))); return env[i]; })();
    const checks = [
      checkMin('Operability over the scatter diagram', 100 * op.operability, inp.opReq, '%', 'DNV-ST-N001 (2021) §2.6; project availability target',
        'Guards the schedule and revenue rather than the structure: below the target the operation waits on weather for more of the year than the project can afford.'),
      check('Significant heave amplitude in the design sea', h.significant, inp.limHeave, 'm', 'Operational limiting criterion (project value)', 'Whether the operation can continue in the stated design sea; exceeding it means the design sea is beyond the working limit.'),
      check('Significant roll / pitch in the design sea', Math.max(ro.significant, pi.significant) * dofScale(3), inp.limAngle, 'deg', 'CAP 437; NORDFORSK (1987)', 'Inclination limit for helicopter, crane and process operations in the design sea.'),
      check('Significant vertical acceleration at the point in the design sea', ac.significant, inp.limAcc, 'm/s²', 'NORDFORSK (1987)', 'Acceleration limit for people and equipment at the point of interest in the design sea.'),
    ];
    const gov = worst(checks);
    const T = om.map((w) => round(2 * Math.PI / w, 4));
    const downH = op.downtime * 8766;
    const govCount = Object.entries(op.governingCounts || {}).sort((a, b) => b[1] - a[1])[0];
    const hsVals = [...new Set(sc.map((s) => s.Hs))].sort((a, b) => a - b); const tpVals = [...new Set(sc.map((s) => s.Tp))].sort((a, b) => a - b);
    const z = hsVals.map((H) => tpVals.map((Tp) => { const c = op.perState.find((s) => s.Hs === H && s.Tp === Tp); return c ? round(c.utilisation, 3) : 0; }));
    return {
      kpis: [
        { key: 'operability', label: 'Operability', value: round(100 * op.operability, 4), unit: '%', limit: inp.opReq, sense: 'min', explain: 'Share of the year in which every motion criterion is met, summed over the scatter diagram.' },
        { key: 'downtime', label: 'Weather downtime', value: round(downH, 4), unit: 'h/yr', explain: 'Hours per year in which at least one criterion is exceeded: the time the operation waits on weather.' },
        { key: 'HsLimit', label: 'Limiting Hs at the design Tp', value: round(HsLimAtTp, 4), unit: 'm', explain: 'Largest significant wave height at this peak period for which all three criteria are met; the working limit to write into the operating procedure.' },
        { key: 'heaveSig', label: 'Significant heave amplitude (design sea)', value: round(h.significant, 4), unit: 'm', limit: inp.limHeave, sense: 'max', explain: 'Twice the standard deviation of heave in the design sea state.' },
        { key: 'rollSig', label: 'Significant roll amplitude (design sea)', value: round(ro.significant * dofScale(3), 4), unit: 'deg', explain: 'Twice the standard deviation of roll in the design sea.' },
        { key: 'pitchSig', label: 'Significant pitch amplitude (design sea)', value: round(pi.significant * dofScale(4), 4), unit: 'deg', explain: 'Twice the standard deviation of pitch in the design sea.' },
        { key: 'accSig', label: 'Significant vertical acceleration at the point', value: round(ac.significant, 4), unit: 'm/s²', limit: inp.limAcc, sense: 'max', explain: 'Twice the standard deviation of vertical acceleration at the point in the design sea.' },
        { key: 'heaveMpm', label: 'Most probable maximum heave', value: round(h.mpm, 4), unit: 'm', explain: 'Largest heave amplitude expected during the sea-state duration.' },
        { key: 'TzHeave', label: 'Heave zero-crossing period', value: round(h.Tz, 4), unit: 's', explain: 'Mean period of the heave response, 2π√(m0/m2); it sets the number of cycles for extremes and fatigue.' },
        { key: 'bandwidth', label: 'Heave spectral bandwidth ε', value: round(eps, 4), unit: '-', explain: '0 for a pure sine, towards 1 for broad-band response. Rayleigh extremes are accurate for narrow-band response (ε below about 0.6).' },
        { key: 'inertiaForce', label: 'Equipment vertical inertia force (maximum)', value: round(Fz / 1e3, 4), unit: 'kN', explain: 'Equipment mass times the most probable maximum vertical acceleration: the dynamic load on its supports in addition to weight.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Wave and heave response spectra', x: { label: 'Frequency', unit: 'rad/s' }, y: { label: 'Spectral density', unit: 'm²·s/rad' },
          series: [{ name: 'Wave spectrum', x: rnd(om), y: rnd(S) }, { name: 'Heave response spectrum', x: rnd(om), y: rnd(r.amp[2].map((a, i) => a * a * S[i])) }],
          explain: 'The heave spectrum is the wave spectrum multiplied by the squared RAO. Where it rises above the wave spectrum the unit amplifies the waves (resonance); the area under each curve is the variance.' },
        { type: 'line', title: 'Limiting significant wave height', x: { label: 'Peak period', unit: 's' }, y: { label: 'Limiting Hs', unit: 'm' },
          series: [...lh.map((l) => ({ name: l.name, x: Tps, y: rnd(l.HsLimit.map((v) => Math.min(v, 30))) })), { name: 'Design sea', x: [inp.Tp], y: [inp.Hs] }],
          explain: 'Workability envelope: below every curve the operation is allowed. The lowest curve at each period governs; a dip marks the period band where the unit responds most.' },
        { type: 'heatmap', title: 'Criterion utilisation over the scatter diagram', xs: tpVals, ys: hsVals, z, zLabel: 'Utilisation', zUnit: '-', x: { label: 'Tp', unit: 's' }, y: { label: 'Hs', unit: 'm' },
          explain: 'Largest response-to-limit ratio in each sea state of the scatter diagram (0 where no sea state is listed). Cells above 1.0 are downtime.' },
      ],
      tables: [
        { title: 'Sea-state results', columns: [{ key: 'Hs', label: 'Hs', unit: 'm' }, { key: 'Tp', label: 'Tp', unit: 's' }, { key: 'p', label: 'Probability', unit: '%' }, { key: 'u', label: 'Utilisation' }, { key: 'g', label: 'Governing' }, { key: 'ok', label: 'Workable' }],
          rows: op.perState.map((s) => ({ Hs: s.Hs, Tp: s.Tp, p: round(100 * s.prob, 4), u: round(s.utilisation, 4), g: s.governing || '–', ok: s.ok ? 'yes' : 'no' })),
          explain: 'One row per sea state: its probability, the highest utilisation among the criteria, the criterion responsible and whether the operation can proceed.' },
      ],
      notes: [
        `Operability ${fmt(100 * op.operability)} % (${fmt(downH)} h/yr downtime)${govCount ? `; most downtime is caused by the ${govCount[0]} criterion` : ''}.`,
        `Design sea Hs ${fmt(inp.Hs)} m / Tp ${fmt(inp.Tp)} s: heave ${fmt(h.significant)} m, roll ${fmt(ro.significant * dofScale(3))}°, pitch ${fmt(pi.significant * dofScale(3))}°, acceleration ${fmt(ac.significant)} m/s² (significant amplitudes).`,
        'The scatter evaluation reuses the drag damping of the design sea; set the design sea near the limiting sea state for the best accuracy.',
      ],
      advice: [
        ...motionRangeAdvice(mot, inp.Tp, inp.headingDeg),
        adv('engineering', sev(gov.uc), `${gov.name}: utilisation ${fmt(gov.uc)}. Limiting Hs at Tp ${fmt(inp.Tp)} s is ${fmt(HsLimAtTp)} m.`,
          gov.uc > 1 ? `Raise workability by relaxing the governing limit with better equipment (heave compensation, motion-compensated gangway), by choosing a heading from the 6-DOF module rosette, or by moving the natural period out of the ${fmt(Tps[env.indexOf(Math.min(...env))])} s band where the envelope is lowest.`
            : `Write the limiting-Hs envelope into the operating procedure as a function of Tp rather than a single Hs: the limit varies from ${fmt(Math.min(...env))} m to ${fmt(Math.min(Math.max(...env), 30))} m over 5–19 s.`,
          'DNV-ST-N001 (2021) §2.6 (operational criteria, alpha factors apply to forecast limits); DNV-RP-C205 (2019) §7.'),
        adv('economic', op.operability * 100 < inp.opReq ? 'warn' : 'info', `Downtime of ${fmt(downH)} h/yr is ${fmt(downH / 24)} spread-days per year waiting on weather.`,
          `Price the downtime with the day rate of the operation and compare with the cost of the mitigation: one percentage point of operability is worth ${fmt(87.66)} h/yr.`,
          'Weather-downtime estimating per DNV-ST-N001 and AACE RP 18R-97 allowances; whole-life costing per ISO 15686-5.',
          `Standby at about 12 t of fuel per day gives ${fuelCo2Text(12 * downH / 24)} per year for one support vessel waiting on weather.`),
      ],
      publish: { 'install.operabilityFrac': op.operability },
      quantities: { vessel_days: { weather_standby: round(downH / 24, 4) } },
    };
  },
  verify: [
    { // 20 × 8 × 0.8 m barge, Hs = 4 m, Tp = 20 s: the wave energy lies at λ ≫ L, heave RAO ≈ 1 → σ_heave = σ_wave = Hs/4,
      // significant amplitude = 2σ = Hs/2 = 2.0 m.
      name: 'Small body in long waves: significant heave amplitude = Hs/2',
      source: 'Definition Hs = 4·√m0 with unit RAO (DNV-RP-C205 §3.5); long-wave limit',
      inputs: { hullType: 'barge', L: 20, B: 8, Dhull: 4, draft: 0.8, KG: 1.5, Hs: 4, Tp: 20, gamma: 1, headingDeg: 180, nFreq: 200 },
      expect: [{ kpi: 'heaveSig', value: 2.0, tol: 0.03 }] },
    { name: 'All sea states workable when the limits are very large', source: 'Limiting case of the operability sum (Σp = 1)', inputs: { limHeave: 20, limAngle: 20, limAcc: 10 },
      expect: [{ kpi: 'operability', value: 100, tol: 1e-9 }, { kpi: 'downtime', value: 0, abs: 1e-9 }] },
  ],
};
