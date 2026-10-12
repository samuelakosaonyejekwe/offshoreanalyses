// Dynamic response of a single mooring line to prescribed fairlead motion: lumped-mass time-domain model with
// Morison loading, seabed contact and friction; tension statistics and dynamic amplification.
import { dynInputs, runLine, LIMIT_LM } from './_lumped.js';
import { round, fmt, adv, check, sev, thin, steelCo2Text, LIMIT_PROPS, CONFIRM } from './_common.js';

export default {
  id: 's10-line-dynamics',
  suite: 10,
  title: 'Dynamic mooring line response',
  summary: 'Simulates one line driven at its fairlead and reports the dynamic tension history, its extremes and the amplification over the quasi-static tension.',
  covers: [2],
  fidelity: 'numerical',
  standards: ['DNV-OS-E301 (2018) Ch.2 Sec.2 [2.6]–[2.7]', 'API RP 2SK (3rd ed., 2005) §5.4, §7.3', 'DNV-RP-C205 (2019) §6 (Morison loading)'],
  models: ['Lumped-mass line model (Hall & Goupee 2015)', 'Morison equation with nonlinear relative-velocity drag', 'Added mass and axial material damping', 'Seabed contact spring with regularised Coulomb friction', 'Dynamic relaxation for the static start state', 'Explicit RK4 sub-stepped on the axial wave speed'],
  limits: `${LIMIT_LM} ${LIMIT_PROPS}`,
  publishes: ['mooring.daf'],
  searchTerms: 'mooring line dynamics lumped mass dynamic tension amplification',
  inputs: [
    ...dynInputs(),
    { key: 'sfDyn', label: 'Required safety factor on dynamic tension', unit: '-', type: 'number', default: 1.67, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria',
      help: `API RP 2SK intact condition with dynamic analysis. ${CONFIRM}` },
  ],
  mesh: { key: 'nSeg', label: 'Line segments', levels: [10, 20, 40], qoi: 'Tmax' },
  meshTime: { key: 'cfl', label: 'Time step / critical step', levels: [0.8, 0.4, 0.2], qoi: 'Tmax', isSize: true },
  compute(inp, ctx) {
    const { p, cat, model, relax, res } = runLine(inp, ctx);
    const st = res.stats;
    const Ts = res.static.fairleadTension;
    const dynInc = st.fairleadMax - Ts;
    const qsInc = st.quasiStaticMax != null ? st.quasiStaticMax - Ts : 0;
    const daf = qsInc > 1e-6 * Ts ? dynInc / qsInc : 1;
    const [t, Tf, Ta, Tq] = thin([res.t, res.fairleadTension, res.anchorTension, res.quasiStaticTension], 600);
    const sMid = Array.from({ length: model.N }, (_, i) => 0.5 * (model.s[i] + model.s[i + 1]));
    const uc = st.fairleadMax * inp.sfDyn / p.MBL;
    const checks = [
      check('Maximum dynamic tension against break load', st.fairleadMax * inp.sfDyn / 1e3, p.MBL / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3 (dynamic analysis)',
        'Guards against breaking the line at the tension peak of the motion cycle; dynamic peaks can far exceed what the static catenary predicts.'),
      { name: 'Line stays taut (no slack)', demand: round(Math.max(0, 0.1 * Ts - st.fairleadMin) / 1e3, 4), capacity: round(0.1 * Ts / 1e3, 4), unit: 'kN', uc: Math.max(0, 0.1 * Ts - st.fairleadMin) / (0.1 * Ts),
        ref: 'DNV-OS-E301 Ch.2 Sec.2 (slack-line criterion: minimum tension above 10 % of static)',
        explain: 'Guards against the line going slack and snapping taut again. A utilisation above zero means the minimum tension dropped below 10 % of the static value: run the snap-load module.' },
    ];
    const mass = p.mass * inp.L / 1e3;
    return {
      kpis: [
        { key: 'Tstatic', label: 'Static fairlead tension', value: round(Ts / 1e3, 5), unit: 'kN', explain: 'Tension of the lumped-mass line at rest; it should match the analytic catenary and is the reference for amplification.' },
        { key: 'Tcat', label: 'Analytic catenary tension', value: round(cat.T / 1e3, 5), unit: 'kN', explain: 'Closed-form elastic catenary tension for the same geometry (no current); the difference from the static lumped-mass value measures discretisation error.' },
        { key: 'Tmax', label: 'Maximum dynamic tension', value: round(st.fairleadMax / 1e3, 5), unit: 'kN', limit: round(p.MBL / inp.sfDyn / 1e3, 5), sense: 'max', explain: 'Largest fairlead tension after start-up. This is the value to hold against the break load with the dynamic safety factor.' },
        { key: 'Tmin', label: 'Minimum dynamic tension', value: round(st.fairleadMin / 1e3, 5), unit: 'kN', explain: 'Smallest fairlead tension in the cycle. Close to zero means the line slackens and may snap.' },
        { key: 'Tmean', label: 'Mean tension', value: round(st.fairleadMean / 1e3, 5), unit: 'kN', explain: 'Time-average fairlead tension; drag makes it differ slightly from the static value.' },
        { key: 'Trange', label: 'Cyclic tension range', value: round((st.fairleadMax - st.fairleadMin) / 1e3, 5), unit: 'kN', explain: 'Peak-to-trough tension per motion cycle; this range, not the maximum, drives fatigue.' },
        { key: 'daf', label: 'Dynamic amplification factor', value: round(daf, 4), unit: '-', explain: 'Dynamic tension increase over the quasi-static increase for the same fairlead motion. Near 1 the catenary follows the motion; well above 1 drag "locks" the line and it stretches instead of changing shape.' },
        { key: 'TanchorMax', label: 'Maximum anchor tension', value: round(st.anchorMax / 1e3, 5), unit: 'kN', explain: 'Largest tension reaching the anchor, after friction on the grounded length; the dynamic anchor demand.' },
        { key: 'tdTravel', label: 'Touchdown travel', value: round(st.touchdownMax - st.touchdownMin, 4), unit: 'm', explain: 'Distance the touchdown point moves along the line each cycle; this length of chain is repeatedly lifted and laid, with wear and trenching.' },
        { key: 'util', label: 'Tmax / MBL', value: round(st.fairleadMax / p.MBL, 4), unit: '-', explain: 'Share of the break load used at the dynamic peak (60 % is the API intact limit with dynamic analysis).' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Tension time history', x: { label: 'Time', unit: 's' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Fairlead (dynamic)', x: t.map((v) => round(v)), y: Tf.map((v) => round(v / 1e3)) }, { name: 'Fairlead (quasi-static)', x: t.map((v) => round(v)), y: Tq.map((v) => round(v / 1e3)), dash: true }, { name: 'Anchor (dynamic)', x: t.map((v) => round(v)), y: Ta.map((v) => round(v / 1e3)) }],
          explain: 'Dynamic fairlead and anchor tension against the quasi-static prediction for the same motion. Look for peaks rising above the dashed curve (amplification), phase shift (drag) and troughs approaching zero (slack).' },
        { type: 'line', title: 'Tension envelope along the line', x: { label: 'Arc length from anchor', unit: 'm' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Maximum', x: sMid.map((v) => round(v)), y: Array.from(res.envelopes.Tmax, (v) => round(v / 1e3)) }, { name: 'Minimum', x: sMid.map((v) => round(v)), y: Array.from(res.envelopes.Tmin, (v) => round(v / 1e3)) }, { name: 'Static', x: sMid.map((v) => round(v)), y: Array.from(res.static.segmentTension, (v) => round(v / 1e3)), dash: true }],
          explain: 'Largest and smallest tension reached in each segment. The band is widest where the line is cycled hardest; the grounded part near the anchor is shielded by friction.' },
        { type: 'line', title: 'Line shape envelope', x: { label: 'Horizontal distance from anchor', unit: 'm' }, y: { label: 'Height above seabed', unit: 'm' },
          series: [{ name: 'Highest position', x: Array.from(res.envelopes.xMax, (v) => round(v)), y: Array.from(res.envelopes.zMax, (v) => round(v + inp.h)) }, { name: 'Lowest position', x: Array.from(res.envelopes.xMin, (v) => round(v)), y: Array.from(res.envelopes.zMin, (v) => round(v + inp.h)) }],
          explain: 'Envelope of node positions over the simulation. The gap between the curves near the seabed marks the touchdown zone that is repeatedly lifted off and re-laid.' },
        { type: 'line', title: 'Energy diagnostics', x: { label: 'Time', unit: 's' }, y: { label: 'Energy', unit: 'kJ' },
          series: ['kinetic', 'strain', 'fairleadWork'].map((k) => { const [tt, e] = thin([res.t, res.energy[k]], 300); return { name: k === 'fairleadWork' ? 'Work done by fairlead' : `${k[0].toUpperCase()}${k.slice(1)} energy`, x: tt.map((v) => round(v)), y: e.map((v) => round(v / 1e3)) }; }),
          explain: 'Kinetic and elastic strain energy of the line and the cumulative work put in at the fairlead. Work rising steadily while the stored energies stay bounded is the energy dissipated by drag and friction — a consistency check on the integration.' },
      ],
      tables: [
        { title: 'Tension statistics', columns: [{ key: 'q', label: 'Quantity' }, { key: 'fair', label: 'Fairlead', unit: 'kN' }, { key: 'anch', label: 'Anchor', unit: 'kN' }],
          rows: [{ q: 'Static', fair: round(Ts / 1e3), anch: round(res.static.anchorTension / 1e3) }, { q: 'Maximum', fair: round(st.fairleadMax / 1e3), anch: round(st.anchorMax / 1e3) }, { q: 'Minimum', fair: round(st.fairleadMin / 1e3), anch: round(st.anchorMin / 1e3) },
            { q: 'Mean', fair: round(st.fairleadMean / 1e3), anch: round(0.5 * (st.anchorMax + st.anchorMin) / 1e3) }, { q: 'Cyclic range', fair: round((st.fairleadMax - st.fairleadMin) / 1e3), anch: round((st.anchorMax - st.anchorMin) / 1e3) }],
          explain: 'Maximum, minimum, mean and cyclic tension at both ends after the start-up ramp. The anchor column mid-value is the mid-range, not a time average.' },
      ],
      notes: [`Static lumped-mass tension ${fmt(Ts / 1e3)} kN vs analytic catenary ${fmt(cat.T / 1e3)} kN (${fmt(100 * (Ts / cat.T - 1))} % difference${inp.Uc ? ', including current' : ''}).`,
        `${res.steps} internal steps of ${fmt(res.dtInternal * 1e3)} ms; static relaxation converged in ${relax.iterations} iterations.`,
        'Drag and added-mass coefficients are indicative defaults from DNV-OS-E301; chain properties come from estimators unless Custom is chosen.'],
      advice: [
        adv('engineering', sev(Math.max(uc, checks[1].uc)), `Dynamic maximum ${fmt(st.fairleadMax / 1e3)} kN (${fmt(100 * st.fairleadMax / p.MBL)} % MBL), minimum ${fmt(st.fairleadMin / 1e3)} kN, amplification factor ${fmt(daf)} over quasi-static.`,
          uc > 1 ? `Increase the break load by ${fmt(100 * (uc - 1))} % or reduce the dynamic increment: a lighter or more compliant upper section (wire or polyester insert) cuts drag-induced tension.`
            : daf > 1.5 ? `Use DAF = ${fmt(daf)} in the tension-and-offset module instead of a quasi-static analysis; quasi-static tension would under-predict the wave-frequency increment by ${fmt(100 * (1 - 1 / daf))} %.`
              : 'Quasi-static analysis is adequate for this line and motion; carry the amplification factor forward and re-check at the shortest governing wave period.',
          'DNV-OS-E301 (2018) Ch.2 Sec.2 [2.7] (dynamic line analysis required when quasi-static is non-conservative); API RP 2SK (3rd ed.) §7.3.'),
        adv('engineering', 'info', `Touchdown point travels ${fmt(st.touchdownMax - st.touchdownMin)} m per cycle; segmentation ${inp.nSeg} segments.`,
          'Run the discretisation study (segments 10/20/40 and time-step ratio) before quoting the maximum to better than a few per cent; put the thickest chain or a wear allowance in the touchdown zone.',
          'DNV-OS-E301 (2018) Ch.2 Sec.2 [2.7.4] (model discretisation to be verified by convergence).'),
        adv('economic', 'info', `The line holds ${fmt(mass)} t of ${p.family === 'fibre' ? 'rope' : 'steel'}; ${fmt(100 * (1 - uc))} % strength reserve at the dynamic peak with safety factor ${inp.sfDyn}.`,
          uc < 0.8 ? `A diameter about ${fmt(100 * (1 - Math.sqrt(Math.max(uc, 0.05) / 0.95)))} % smaller would still pass this load case, saving about ${fmt(mass * (1 - Math.max(uc, 0.05) / 0.95))} t per line; check fatigue and corrosion allowance first.` : 'Little to trim on strength; cost reduction would have to come from line length or a cheaper grade.',
          'Whole-life costing per ISO 15686-5; quantity-based Class 4 estimate per AACE RP 18R-97.',
          p.family === 'fibre' ? 'Fibre rope carries roughly a tenth of the mass of chain for the same strength; its embodied carbon per tonne is higher but the total is usually lower.' : steelCo2Text(mass * Math.max(0, 1 - Math.max(uc, 0.05) / 0.95) || mass)),
      ],
      publish: { 'mooring.daf': daf },
      quantities: p.family === 'fibre' ? { polymer_t: round(mass) } : { steel_t: round(mass) },
    };
  },
  verify: [
    { // Taut near-vertical elastic line, no motion: anchor 280 m below the fairlead and 1 m aside, L = 279.5 m, EA = 800 MN,
      // w = 1500 N/m. Straight length = sqrt(280² + 1²) = 280.0018 m → stretch 0.5018 m → mean tension EA·ΔL/L =
      // 800e6·0.5018/279.5 = 1436.3 kN; tension varies linearly with the weight, so top = mean + w·L/2 = 1436.3 + 209.6 = 1645.9 kN.
      name: 'Static top tension of a taut vertical elastic line',
      source: 'Hand calculation: Hooke\'s law plus self-weight (Irvine 1981, Cable Structures, §1)',
      inputs: { lineType: 'custom', wCustom: 1500, EAcustom: 800, MBLcustom: 9000, massCustom: 176, L: 279.5, X: 1, h: 280, ax: 0, az: 0, nCycles: 2, period: 5 },
      expect: [{ kpi: 'Tmax', value: 1645.9, tol: 0.01 }, { kpi: 'Tstatic', value: 1645.9, tol: 0.01 }] },
    { // No fairlead motion on the default catenary: the dynamic model must stay at its static state, so the cyclic range is zero.
      name: 'No motion: zero cyclic tension range',
      source: 'Limiting case (static equilibrium is a fixed point of the equations of motion)',
      inputs: { ax: 0, az: 0, nCycles: 2 }, expect: [{ kpi: 'Trange', value: 0, abs: 0.5 }] },
  ],
};
