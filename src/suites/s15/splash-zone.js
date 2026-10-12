// Dynamic lift through the splash zone by the DNV-RP-N103 Simplified Method: slamming, varying buoyancy, mass
// and drag forces at successive submergence stages, slack-sling and capacity criteria, snap load and the
// limiting sea state against zero-crossing period.
import { craneTipMotion, waveKinematicsN103, splashZoneN103, snapLoadN103 } from '../../engine/lines/lowering.js';
import { G, need, run, round, fmt, kpi, check, adv, line, bar, tbl, worst, linspace, vesselCo2Text, CONFIRM, N103, LIMIT_N103 } from './_common.js';

export default {
  id: 's15-splash-zone', suite: 15, title: 'Dynamic lift through the splash zone',
  summary: 'Evaluates the hydrodynamic force on an object lowered through the wave zone stage by stage, checks slack slings and hook-load capacity, and finds the limiting significant wave height for each wave period.',
  covers: [10], fidelity: 'screening',
  standards: [`${N103} §4 (lifting through wave zone — Simplified Method), §4.7 (snap forces)`, 'DNV-ST-N001 (2021) Sec. 16 (dynamic hook load, rigging design)', 'DNV-RP-C205 (2019) §8 (slamming), App. D (added mass)'],
  models: ['Morison-type wave loading in characteristic form', 'Slamming impact force ½·ρ·Cs·As·vs²', 'Time-varying buoyancy', 'Added-mass and inertia force', 'Quadratic drag', 'Vessel heave coupling through the crane-tip motion', 'Slack-sling criterion and snap load'],
  limits: `${LIMIT_N103} The Simplified Method combines characteristic amplitudes rather than simulating the passage: it assumes a small object compared with the wavelength, vertical motion only, and wave kinematics unaffected by the vessel. It is intentionally conservative; a time-domain lowering analysis with measured or computed coefficients (model tests, CFD — outside this module) usually raises the limiting sea state.`,
  publishes: ['install.splashHsLimit', 'install.splashDaf'],
  searchTerms: 'splash zone lifting simplified method slamming slack sling limiting sea state DNV-RP-N103',
  inputs: [
    { key: 'mass', label: 'Object mass in air', unit: 't', type: 'number', default: 180, min: 0.5, max: 20000, group: 'Object', role: 'load', bind: 'subsea.liftMass', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'vol', label: 'Displaced volume when fully submerged', unit: 'm³', type: 'number', default: 28, min: 0, max: 20000, group: 'Object', role: 'geometry' },
    { key: 'A33', label: 'Heave added mass when fully submerged', unit: 't', type: 'number', default: 320, min: 0, max: 1e5, group: 'Object', role: 'load', uncertain: { dist: 'lognormal', cov: 0.25 }, help: `From DNV-RP-N103 §4.6 / DNV-RP-C205 App. D for the object shape, with perforation reduction where applicable. ${CONFIRM}` },
    { key: 'hObj', label: 'Object height', unit: 'm', type: 'number', default: 6, min: 0.2, max: 80, group: 'Object', role: 'geometry' },
    { key: 'As', label: 'Slamming area (bottom, exposed to water entry)', unit: 'm²', type: 'number', default: 40, min: 0, max: 5000, group: 'Object', role: 'geometry', design: true, help: 'Horizontal projected area that hits the water at first entry. Perforated or sloped bottoms reduce it.' },
    { key: 'Cs', label: 'Slamming coefficient', unit: '-', type: 'number', default: 5, min: 1, max: 12, group: 'Coefficients', role: 'criteria', uncertain: { dist: 'normal', cov: 0.2 }, help: `Not less than 3.0 for smooth circular cylinders and 5.0 for other shapes. ${CONFIRM}` },
    { key: 'Ap', label: 'Vertical drag area (horizontal projection)', unit: 'm²', type: 'number', default: 120, min: 0.1, max: 5000, group: 'Object', role: 'geometry' },
    { key: 'CD', label: 'Drag coefficient in oscillatory flow', unit: '-', type: 'number', default: 2.5, min: 0.5, max: 8, group: 'Coefficients', role: 'criteria', uncertain: { dist: 'normal', cov: 0.2 }, help: `Not less than 2.5 unless documented by tests or CFD. ${CONFIRM}` },
    { key: 'Aw', label: 'Waterplane area while partly submerged', unit: 'm²', type: 'number', default: 12, min: 0, max: 5000, group: 'Object', role: 'geometry', help: 'Mean water-line area of the members piercing the surface: drives the varying buoyancy force.' },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 1.5, min: 0.1, max: 8, group: 'Sea state', role: 'environment', bind: 'install.HsOperation', design: true, uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'Tz', label: 'Zero-crossing period', unit: 's', type: 'number', default: 6, min: 2.5, max: 14, group: 'Sea state', role: 'environment', help: 'The method is stated for 8.9·√(Hs/g) ≤ Tz ≤ 13 s; the full range of periods at the site must be covered.' },
    { key: 'kAmp', label: 'Characteristic wave amplitude / Hs', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1.2, group: 'Coefficients', role: 'criteria', help: `0.9 for a splash-zone passage of up to 30 minutes, 1.0 for longer. ${CONFIRM}` },
    { key: 'ctPerHs', label: 'Crane-tip vertical motion amplitude per metre of Hs', unit: 'm/m', type: 'number', default: 0.45, min: 0, max: 3, group: 'Crane-tip motion', role: 'bc', bind: 'install.craneTipHeavePerHs', uncertain: { dist: 'normal', cov: 0.2 }, help: 'Boundary condition: characteristic single amplitude of the crane tip divided by Hs, from the vessel-motion module at the lift heading.' },
    { key: 'Tct', label: 'Crane-tip motion period', unit: 's', type: 'number', default: 7.5, min: 3, max: 25, group: 'Crane-tip motion', role: 'bc', bind: 'install.craneTipPeriod' },
    { key: 'vc', label: 'Lowering speed through the splash zone', unit: 'm/s', type: 'number', default: 0.3, min: 0, max: 2, group: 'Operation', role: 'ic', design: true, help: `Initial condition of water entry. The method assumes about 0.5 m/s if unknown. ${CONFIRM}` },
    { key: 'dafRig', label: 'DAF the rigging and lift points are designed for', unit: '-', type: 'number', default: 1.6, min: 1, max: 4, group: 'Capacity', role: 'criteria', bind: 'install.daf', help: 'The converted dynamic amplification from this analysis must not exceed the factor used in the lift design.' },
    { key: 'craneCap', label: 'Crane dynamic capacity at the lift radius', unit: 't', type: 'number', default: 400, min: 1, max: 20000, group: 'Capacity', role: 'bc' },
    { key: 'Krig', label: 'Hoisting system stiffness (for snap load)', unit: 'kN/m', type: 'number', default: 15000, min: 50, max: 1e6, group: 'Capacity', role: 'bc' },
  ],
  compute(inp) {
    // t → kg, kN/m → N/m.
    const rho = 1025, M = inp.mass * 1e3, A33 = inp.A33 * 1e3, Wair = M * G;
    need(M - rho * inp.vol > 0, 'The object is buoyant when submerged: it cannot be lowered through the splash zone on slings alone.');
    // Submergence stages: fraction of volume and added mass active, distance of the item centre below still water.
    const stages = [
      { name: '1 Bottom at the surface (slamming)', fV: 0, fA: 0, slam: 1, wp: 0, drag: 0, d: 0 },
      { name: '2 Half submerged', fV: 0.5, fA: 0.5, slam: 0, wp: 1, drag: 0.5, d: 0 },
      { name: '3 Just submerged', fV: 1, fA: 1, slam: 0, wp: 0, drag: 1, d: inp.hObj / 2 },
      { name: '4 One height below the surface', fV: 1, fA: 1, slam: 0, wp: 0, drag: 1, d: 1.5 * inp.hObj },
    ];
    const evalAt = (Hs, Tz) => {
      const ct = craneTipMotion({ heaveAmp: inp.ctPerHs * Hs, T3: inp.Tct });
      return stages.map((s) => {
        const wk = waveKinematicsN103({ Hs, Tz, d: s.d, amplitudeFactor: inp.kAmp });
        const r = splashZoneN103({ mass: M, volume: inp.vol * s.fV, addedMass: A33 * s.fA, slamArea: inp.As * s.slam, Cs: inp.Cs, dragArea: inp.Ap * s.drag, CD: inp.CD, waterplaneArea: inp.Aw * s.wp, vc: inp.vc,
          etaCt: ct.eta, vCt: ct.velocity, aCt: ct.acceleration, zetaA: wk.zetaA, vW: wk.velocity, aW: wk.acceleration, rho });
        return { ...r, stage: s.name, inRange: wk.tzInRange };
      });
    };
    const res = run(() => evalAt(inp.Hs, inp.Tz), 'Simplified Method');
    const slackUc = Math.max(...res.map((r) => r.slackUc)), dafConv = Math.max(...res.map((r) => r.dafConv)), Ftot = Math.max(...res.map((r) => r.Ftotal));
    const iSlack = res.findIndex((r) => r.slackUc === slackUc), iCap = res.findIndex((r) => r.Ftotal === Ftot);
    const ucAt = (Hs, Tz) => { const r = evalAt(Hs, Tz); return Math.max(...r.map((q) => q.slackUc), ...r.map((q) => q.dafConv / inp.dafRig), ...r.map((q) => q.Ftotal / (inp.craneCap * 1e3 * G))); };
    // Limiting Hs per zero-crossing period: bisection on the governing utilisation (monotonic in Hs).
    const hsLimit = (Tz) => { if (ucAt(0.05, Tz) > 1) return 0; let lo = 0.05, hi = 8; if (ucAt(hi, Tz) <= 1) return hi; for (let i = 0; i < 40; i++) { const m = 0.5 * (lo + hi); if (ucAt(m, Tz) <= 1) lo = m; else hi = m; } return lo; };
    const tzs = linspace(4, 12, 17), lims = tzs.map((t) => run(() => hsLimit(t), 'Limiting sea state'));
    const hsLim = run(() => hsLimit(inp.Tz), 'Limiting sea state'), hsLimMin = Math.min(...lims);
    // Snap load if a sling goes slack: free-fall velocity from weight in water against drag.
    const wk0 = waveKinematicsN103({ Hs: inp.Hs, Tz: inp.Tz, d: 0, amplitudeFactor: inp.kAmp }), ct0 = craneTipMotion({ heaveAmp: inp.ctPerHs * inp.Hs, T3: inp.Tct });
    const vr = inp.vc + Math.hypot(ct0.velocity, wk0.velocity);
    const vff = Math.sqrt(2 * (M - rho * inp.vol) * G / (rho * inp.CD * inp.Ap));
    const snap = vr > 0 ? snapLoadN103({ mass: M, addedMass: A33, K: inp.Krig * 1e3, vff, vr }) : { Fsnap: 0, vSnap: 0 };
    const checks = [
      check(`Slack sling (stage ${iSlack + 1})`, res[iSlack].Fhyd / 1e3, 0.9 * res[iSlack].staticMin / 1e3, 'kN', `${N103} §4.4.3 (Fhyd ≤ 0.9·Fstatic-min)`, 'If the hydrodynamic force exceeds 90 % of the static weight the slings can go slack and then snatch: the snap load may be several times the static load.'),
      check(`Converted DAF against the rigging design DAF (stage ${iCap + 1})`, dafConv, inp.dafRig, '-', `${N103} §4.4.2–4.4.4 (DAFconv = Ftotal/(M·g)); DNV-ST-N001 Sec. 16`, 'The total force in the wave zone, expressed as a factor on the weight in air, must not exceed the dynamic factor the slings, shackles and padeyes were designed for.'),
      check('Total force against crane dynamic capacity', Ftot / G / 1e3, inp.craneCap, 't', 'DNV-ST-N001 Sec. 16 (crane capacity); crane load chart for subsea lifts', 'Largest characteristic hook load during the passage against the crane rating at this radius.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const inRange = res[0].inRange;
    return {
      kpis: [
        kpi('Fhyd', 'Maximum hydrodynamic force', Math.max(...res.map((r) => r.Fhyd)) / 1e3, 'kN', 'Characteristic hydrodynamic force √((FD + Fslam)² + (FM − Fρ)²), largest over the stages.'),
        kpi('Fslam', 'Slamming force at first water entry', res[0].Fslam / 1e3, 'kN', '½·ρ·Cs·As·vs² with the slamming velocity from lowering speed, crane-tip velocity and wave particle velocity.'),
        kpi('vs', 'Slamming impact velocity', res[0].slamVelocity, 'm/s', 'vc + √(vct² + vw²): relative velocity between the object bottom and the water surface.'),
        kpi('slackUc', 'Slack-sling utilisation', slackUc, '-', 'Fhyd over 0.9 × the minimum static sling load. Above 1.0 slack slings and snap loads must be expected.', 1, 'max'),
        kpi('dafConv', 'Converted dynamic amplification factor', dafConv, '-', '(static + hydrodynamic force)/(weight in air), largest over the stages. Compare with the DAF used for rigging design.', inp.dafRig, 'max'),
        kpi('Ftot', 'Maximum total hook load', Ftot / 1e3, 'kN', 'Maximum static load plus hydrodynamic force.', inp.craneCap * G, 'max'),
        kpi('hsLim', 'Limiting Hs at the given period', hsLim, 'm', 'Largest significant wave height that satisfies slack, DAF and crane criteria at this Tz. This is the design limit; the forecast limit is lower by the alpha-factor.'),
        kpi('hsLimMin', 'Limiting Hs over all periods 4–12 s', hsLimMin, 'm', 'Lowest limit over the period range: the value to use when the period cannot be restricted in the procedure.'),
        kpi('Fsnap', 'Snap load if a sling goes slack', snap.Fsnap / 1e3, 'kN', 'vsnap·√(K·(M + A33)): the load when a slack sling re-tightens. Shown for information; it applies only if the slack criterion fails.'),
        kpi('etaCt', 'Crane-tip motion amplitude at this Hs', ct0.eta, 'm', 'Characteristic single amplitude used for the crane-tip velocity and acceleration.'),
      ],
      checks,
      plots: [
        bar('Force components by submergence stage', stages.map((s) => s.name), 'Force', 'kN', [{ name: 'Slamming', y: res.map((r) => round(r.Fslam / 1e3)) }, { name: 'Drag', y: res.map((r) => round(r.FD / 1e3)) }, { name: 'Mass (inertia)', y: res.map((r) => round(r.FM / 1e3)) }, { name: 'Varying buoyancy', y: res.map((r) => round(r.Frho / 1e3)) }, { name: 'Hydrodynamic total', y: res.map((r) => round(r.Fhyd / 1e3)) }, { name: '0.9 × static minimum', y: res.map((r) => round(0.9 * r.staticMin / 1e3)) }],
          'Which force governs at each stage: slamming at first entry, then inertia and drag as added mass builds up. The hydrodynamic total must stay below the 0.9 × static bar at every stage.'),
        line('Limiting sea state against wave period', 'Zero-crossing period Tz', 's', 'Limiting Hs', 'm', [{ name: 'Limiting Hs', x: tzs.map((t) => round(t)), y: lims.map((v) => round(v)) }, { name: 'Operating point', x: [inp.Tz], y: [inp.Hs] }],
          'The operation is acceptable below the curve. Short periods give high particle velocities and accelerations, so the curve usually rises with period; plan the procedure for the lowest part of the curve that can occur.'),
      ],
      tables: [tbl('Simplified Method results by stage', [['s', 'Stage'], ['Fst', 'Static load', 'kN'], ['Fslam', 'Slam', 'kN'], ['FD', 'Drag', 'kN'], ['FM', 'Mass', 'kN'], ['Fr', 'Varying buoyancy', 'kN'], ['Fh', 'Fhyd', 'kN'], ['uc', 'Slack utilisation'], ['daf', 'DAFconv']],
        res.map((r) => ({ s: r.stage, Fst: round(r.staticMin / 1e3), Fslam: round(r.Fslam / 1e3), FD: round(r.FD / 1e3), FM: round(r.FM / 1e3), Fr: round(r.Frho / 1e3), Fh: round(r.Fhyd / 1e3), uc: round(r.slackUc, 3), daf: round(r.dafConv, 3) })),
        'One row per submergence level with the force components of DNV-RP-N103 §4.3 and the two acceptance measures. The worst row sets the limiting sea state.')],
      notes: [`Stage ${iSlack + 1} governs the slack criterion (${fmt(slackUc)}); stage ${iCap + 1} governs the hook load (DAFconv ${fmt(dafConv)}).`, inRange ? 'Tz is within the stated validity range of the Simplified Method.' : 'Tz is outside the stated validity range 8.9·√(Hs/g) ≤ Tz ≤ 13 s: treat the result with caution.', `Limiting Hs ${fmt(hsLim)} m at Tz ${inp.Tz} s; ${fmt(hsLimMin)} m over 4–12 s.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 || !inRange ? 'warn' : 'ok', `The splash-zone passage is acceptable at Hs ${inp.Hs} m, Tz ${inp.Tz} s: slack utilisation ${fmt(slackUc)}, converted DAF ${fmt(dafConv)}; the design limit at this period is Hs ${fmt(hsLim)} m.`,
          `Write the limiting sea state into the procedure as Hs ${fmt(hsLimMin)} m for all periods (or the curve by period), apply the alpha-factor for the forecast limit, and confirm Cs, CD and added mass against the governing edition or model tests.`, `${N103} §4.3–4.4`)
          : adv('engineering', 'fail', `"${gov.name}" is exceeded (${fmt(gov.uc)}) at Hs ${inp.Hs} m, Tz ${inp.Tz} s; the limit at this period is Hs ${fmt(hsLim)} m. A slack sling would re-tighten with a snap load of about ${fmt(snap.Fsnap / 1e3)} kN.`,
            `Reduce the sea state to Hs ${fmt(hsLim)} m, or cut the hydrodynamic force: perforate or slope the bottom to reduce the slamming area (currently ${inp.As} m²) and added mass, change the lowering speed, or move the lift to a position with less crane-tip motion. A time-domain analysis usually recovers margin.`, `${N103} §4.3–4.7`),
        adv('economic', 'info', `The operation is limited to Hs ${fmt(hsLimMin)} m (design value) over all periods; every 0.5 m of extra limit typically adds several percent of workable time at an exposed site.`,
          'Use the weather-window module with this limit to price the waiting on weather, and compare it with the cost of measures that raise the limit (perforated mudmats, a time-domain analysis, a vessel with better motions or a moonpool).',
          'Operability-driven vessel time: DNV-ST-N001 Sec. 2 (operational limits, alpha-factor); day rate × expected waiting days',
          `A day of waiting on weather with the construction vessel on station burns about ${vesselCo2Text(1, 35)}`),
      ],
      publish: { 'install.splashHsLimit': hsLimMin, 'install.splashDaf': dafConv },
      quantities: {},
    };
  },
  verify: [
    { name: 'Slamming at first entry: Fslam = ½·ρ·Cs·As·(vc + √(vct² + vw²))²', source: 'DNV-RP-N103 §4.3.5 hand calculation: Hs 1.5 m, Tz 6 s → ζa = 1.35 m, vw = 1.35·2π/6 = 1.4137 m/s; crane tip 0.675 m at 7.5 s → vct = 0.5655 m/s; vs = 0.3 + √(0.5655² + 1.4137²) = 1.8226 m/s; Fslam = ½·1025·5·40·1.8226² = 340.5 kN',
      inputs: {}, expect: [{ kpi: 'vs', value: 0.3 + Math.hypot(2 * Math.PI * 0.675 / 7.5, 1.35 * 2 * Math.PI / 6), tol: 1e-9 }, { kpi: 'Fslam', value: 0.5 * 1025 * 5 * 40 * (0.3 + Math.hypot(2 * Math.PI * 0.675 / 7.5, 1.35 * 2 * Math.PI / 6)) ** 2 / 1e3, tol: 1e-9 }] },
    { name: 'Calm-water limit: no waves and no crane-tip motion leave only the entry slam at the lowering speed', source: 'Limiting case: Hs → 0 gives vs = vc and Fslam = ½·ρ·Cs·As·vc² = ½·1025·5·40·0.09 = 9.225 kN (Hs = 0.1 m with the amplitude factor and crane-tip motion set to zero)',
      inputs: { Hs: 0.1, kAmp: 0.5, ctPerHs: 0, Tz: 13 }, expect: [{ kpi: 'Fslam', value: 0.5 * 1025 * 5 * 40 * (0.3 + 0.05 * 2 * Math.PI / 13) ** 2 / 1e3, tol: 1e-9 }] },
  ],
};
