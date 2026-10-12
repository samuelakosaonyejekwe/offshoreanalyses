// Interlayer friction and contact: contact pressures with initial gaps, friction per wire, slip-onset curvature,
// slip per bending cycle, friction stress, hysteresis damping, anti-wear tape wear and sheath creep into gaps.
import { axisymmetricResponse, bendingProperties } from '../../engine/xsection/flexible.js';
import { wireFrictionStress, momentCurvatureLoop } from '../../engine/xsection/helix.js';
import { wireSlipPerCycle, archardWear, sheathCreepIntoGaps } from '../../engine/xsection/flexfatigue.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, muInput, kpi, chk, adv, plot, bar, tbl, col, run, need, fmt, round, linspace, sev, worst, LIMIT_AXI, GENERIC, STD_17J, STD_17B, CONFIRM, DEG, CO2_MGO, MGO_T_PER_DAY } from './_common.js';

export default {
  id: 's14-interlayer-contact',
  suite: 14,
  title: 'Interlayer friction and contact',
  summary: 'Solves the contact pressures between layers and derives friction, armour slip, hysteresis, tape wear and sheath creep from them.',
  covers: [11],
  fidelity: 'analytical',
  standards: ['API RP 17B (5th ed.) §8 (wear, creep, hysteresis)', 'API Spec 17J (4th ed.) §5.3.1 (sheath thickness reduction by creep)'],
  models: ['Unilateral interlayer contact with initial gaps (active-set)', 'Coulomb frictional slip–stick model of helical wires', 'Slip amplitude R²cos²α·Δκ/sinα (Féret & Bournazel 1987)', 'Archard wear', 'Bilinear moment–curvature hysteresis and equivalent damping', 'Volume-balance creep of the sheath into armour gaps'],
  limits: LIMIT_AXI + ' Contact is frictionless in the axisymmetric solve and uniform round the circumference; friction enters only the bending response. Wear coefficient, friction coefficient and creep intrusion depth are test data — the defaults are indicative. ' + GENERIC,
  searchTerms: 'flexible pipe interlayer friction contact pressure wear anti-wear tape',
  inputs: [
    ...stackInputs(),
    ...loadInputs(),
    { key: 'gap0', label: 'Initial radial gap outside the outer tensile armour', unit: 'mm', type: 'number', default: 0, min: 0, max: 5, group: 'Interlayer contact', role: 'ic',
      help: 'Manufacturing gap between the outer tensile armour and the layer outside it. The gap must close before that interface carries pressure.' },
    muInput(),
    { key: 'dKappa', label: 'Cyclic curvature range', unit: '1/m', type: 'number', default: 0.04, min: 0, max: 1, group: 'Global loading', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Peak-to-peak curvature of the dominant bending cycle at the section.' },
    { key: 'cyclesPerYear', label: 'Bending cycles per year', type: 'number', default: 4e6, min: 0, max: 1e8, group: 'Global loading', role: 'bc', help: 'About 4–5 million wave cycles per year for a riser.' },
    { key: 'years', label: 'Service period', unit: 'years', type: 'number', default: 25, min: 0.1, max: 60, group: 'Global loading', role: 'operation', bind: 'life.designYears' },
    { key: 'wearK', label: 'Wear coefficient of the tape / wire pair', unit: 'm²/N', type: 'number', default: 2e-17, min: 0, max: 1e-12, group: 'Interlayer contact', role: 'material',
      uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Dimensional Archard coefficient from wear tests on the actual material pair; reported values span about 1e-18 to 1e-13 m²/N.' },
    { key: 'gapWidth', label: 'Pressure-armour gap width', unit: 'mm', type: 'number', default: 1.5, min: 0, max: 10, group: 'Sheath creep', role: 'geometry', help: 'Gap between adjacent pressure-armour profiles into which the pressure sheath creeps.' },
    { key: 'intrusion', label: 'Creep intrusion depth', unit: 'mm', type: 'number', default: 2, min: 0, max: 15, group: 'Sheath creep', role: 'material', help: 'Depth to which the sheath material extrudes into the gap at design pressure and temperature, from creep tests.' },
    { key: 'creepAllow', label: 'Allowable sheath thickness reduction', unit: '%', type: 'number', default: 30, min: 5, max: 50, group: 'Acceptance criteria', role: 'criteria', help: 'API 17J limits thinning of the pressure sheath by creep to 30 %. ' + CONFIRM },
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    need(S.idx.ta.length > 0, 'The stack has no tensile armour layer.');
    const N = S.stack.layers.length;
    const gaps = new Array(N - 1).fill(0);
    const last = S.idx.ta[S.idx.ta.length - 1];
    if (last < N - 1) gaps[last] = inp.gap0 / 1e3;
    const loads = loadsFrom(inp);
    const r = run(() => axisymmetricResponse(S.stack, { ...loads, initialGaps: gaps }), 'Axisymmetric contact solution');
    const bend = run(() => bendingProperties(S.stack, { response: r, mu: inp.mu }), 'Bending properties');
    const amp = inp.dKappa / 2, nCyc = inp.cyclesPerYear * inp.years;
    const rows = S.idx.ta.map((i, k) => {
      const L = S.stack.layers[i], b = bend.layers[k], n = L.n ?? 1, A = L.nA / n;
      const slips = amp > (b.kappaCr ?? 0);
      const slip = wireSlipPerCycle({ R: L.R, alpha: L.alpha, kappaRange: inp.dKappa });
      const inner = S.stack.layers[i - 1];
      const tape = inner && (inner.type === 'antiWear' || inner.type === 'tape') ? inner : null;
      const q = Math.max(r.layers[i].contactIn, r.layers[i].contactOut, 0) / L.fill; // pressure on the wire footprint
      const wear = archardWear({ wearCoefficient: inp.wearK, contactPressure: q, slidingPerCycle: slips ? slip.slidingPerCycle : 0, cycles: nCyc, thickness: tape ? tape.t : (L.wireThickness ?? L.t) });
      const st = wireFrictionStress({ E: L.E, A, R: L.R, alpha: L.alpha, kappa: amp, psi: 0, frictionPerLength: b.frictionPerLength ?? 0 });
      return { L, b, A, slip, slips, tape, q, wear, st };
    });
    const out = rows[rows.length - 1];
    const kY = bend.frictionMoment / Math.max(bend.EIstick - bend.EIslip, 1e-9);
    const loop = amp > kY * 1.0001 ? momentCurvatureLoop({ EIstick: bend.EIstick, EIslip: bend.EIslip, frictionMoment: bend.frictionMoment, kappaAmp: amp, nPoints: 21 }) : null;
    const maxContact = Math.max(0, ...r.interfaces.map((f) => f.contactPressure));
    const open = r.interfaces.filter((f) => !f.closed).length;
    const checks = rows.map((y) => chk(`${y.tape ? y.tape.name : y.L.name}: wear depth over ${fmt(inp.years)} years`, y.wear.wearDepth * 1e3, (y.tape ? y.tape.t : (y.L.wireThickness ?? y.L.t)) * 1e3, 'mm', 'Archard (1953) wear law; API RP 17B §8 (wear of anti-wear layers)',
      y.tape ? 'Guards against the anti-wear tape wearing through, after which steel slides on steel, friction rises and the armour wires wear and fret.' : 'Guards against loss of wire section by wear where there is no anti-wear tape; a thinner wire carries a higher stress.'));
    let creep = null;
    if (S.idx.pa >= 0 && S.idx.sheath >= 0) {
      const P = S.stack.layers[S.idx.pa], Sh = S.stack.layers[S.idx.sheath];
      const pitch = 2 * Math.PI * P.R / Math.tan(P.alpha);
      creep = run(() => sheathCreepIntoGaps({ gapWidth: Math.min(inp.gapWidth / 1e3, pitch), intrusionDepth: inp.intrusion / 1e3, pitch, sheathThickness: Sh.t, allowableReduction: inp.creepAllow / 100 }), 'Sheath creep');
      checks.push(chk('Pressure sheath thinning by creep into armour gaps', creep.reduction * 100, inp.creepAllow, '%', `${STD_17J} §5.3.1 (thickness reduction ≤ ${inp.creepAllow} %)`,
        'Guards against the pressure sheath thinning as it extrudes into the gaps of the pressure armour; excessive thinning leads to sheath rupture and loss of containment.'));
    }
    const gov = worst(checks);
    const names = S.stack.layers.map((L) => L.name);
    const psis = linspace(0, Math.PI, 37);
    const fr = psis.map((p) => wireFrictionStress({ E: out.L.E, A: out.A, R: out.L.R, alpha: out.L.alpha, kappa: amp, psi: p, frictionPerLength: out.b.frictionPerLength ?? 0 }).sigma);
    const slipMax = Math.max(...rows.map((y) => y.slip.slipRange));
    const wearMax = rows.reduce((a, c) => (c.wear.wearDepth > a.wear.wearDepth ? c : a), rows[0]);
    const yearsToWear = wearMax.wear.wearDepth > 0 ? inp.years * (wearMax.tape ? wearMax.tape.t : wearMax.L.t) / wearMax.wear.wearDepth : Infinity;
    return {
      kpis: [
        kpi('maxContact', 'Largest interface contact pressure', maxContact / 1e6, 'MPa', 'Highest layer-average pressure between adjacent layers under the applied tension and pressures.'),
        kpi('openInterfaces', 'Open interfaces', open, '-', 'Number of interfaces carrying no contact pressure because a gap has not closed or the layers have separated.'),
        kpi('frictionPerLength', 'Friction per unit wire length, outer armour', (out.b.frictionPerLength ?? 0) / 1e3, 'kN/m', 'Friction force available along one wire of the outer armour: μ × contact pressure on both faces × tributary width.'),
        kpi('kappaCr', 'Slip-onset curvature', bend.kappaCrMin ?? 0, '1/m', 'Curvature at which the first armour layer starts to slide. Cycles smaller than this stay in stick and cause no wear.'),
        kpi('slipRange', 'Wire slip per cycle (range)', slipMax * 1e3, 'mm', 'Peak-to-peak relative sliding between a fully slipped wire and the layer beneath, at the neutral axis of bending: R²·cos²α·Δκ/sinα.'),
        kpi('frictionStress', 'Friction stress at full slip', Math.max(...rows.map((y) => y.b.sigmaFrictionMax ?? 0)) / 1e6, 'MPa', 'Largest axial stress friction can build up in a wire; twice this is the friction contribution to the fatigue stress range.'),
        kpi('frictionMoment', 'Friction moment', bend.frictionMoment / 1e3, 'kN·m', 'Bending moment carried by friction at full slip — the height of the hysteresis loop.'),
        kpi('damping', 'Equivalent hysteretic damping', loop ? loop.dampingRatio : 0, '-', 'Equivalent viscous damping ratio of the bending cycle from friction; zero if the cycle stays in stick.'),
        kpi('energyPerCycle', 'Energy dissipated per cycle', loop ? loop.energyPerCycle : 0, 'J/m', 'Friction work per bending cycle and per metre of pipe; it ends up as heat and wear.'),
        kpi('wearDepth', 'Largest wear depth', wearMax.wear.wearDepth * 1e3, 'mm', `Archard wear depth over the service period at ${wearMax.tape ? wearMax.tape.name : wearMax.L.name}.`),
        kpi('wearUc', 'Wear utilisation', Math.max(...rows.map((y, k) => checks[k].uc)), '-', 'Wear depth over the available thickness. At 1.0 the anti-wear tape is worn through.', 1, 'max'),
        kpi('creepUc', 'Sheath creep utilisation', creep ? creep.utilisation : 0, '-', 'Thinning of the pressure sheath by creep into the armour gaps over the allowable reduction.', 1, 'max'),
      ],
      checks,
      plots: [
        bar('Interface contact pressure', r.interfaces.map((_, j) => `${names[j]} / ${names[j + 1]}`), 'Contact pressure', 'MPa', [{ name: 'Contact pressure', y: r.interfaces.map((f) => round(f.contactPressure / 1e6)) }],
          'Pressure at every interface from the bore outwards. Bore pressure pushes the inner interfaces together; armour tension squeezes inwards. A zero bar is an open interface.'),
        plot('line', 'Friction stress along the wire, outer armour', 'Angle from extrados', 'deg', 'Friction stress', 'MPa', [{ name: out.st.state, x: psis.map((p) => round(p / DEG)), y: fr.map((v) => round(v / 1e6)) }],
          'Axial stress built up by friction in a wire between the extrados (0°) and intrados (180°) at the curvature amplitude. A cosine shape means stick; a straight line means the wire has slipped.'),
        plot('line', 'Moment–curvature hysteresis loop', 'Curvature', '1/m', 'Moment', 'kN·m', loop
          ? [{ name: 'Steady cycle', x: loop.kappa.map((k) => round(k)), y: loop.moment.map((v) => round(v / 1e3)) }]
          : [{ name: 'Stick (no slip)', x: [-round(amp), round(amp)], y: [-round(bend.EIstick * amp / 1e3), round(bend.EIstick * amp / 1e3)] }],
        'Bending moment against curvature for the cyclic amplitude. The enclosed area is the friction energy per cycle; a single line means the cycle is too small to slip.'),
      ],
      tables: [
        tbl('Interfaces', [col('iface', 'Interface'), col('p', 'Contact pressure', 'MPa'), col('gap', 'Gap', 'mm'), col('state', 'State')],
          r.interfaces.map((f, j) => ({ iface: `${names[j]} / ${names[j + 1]}`, p: round(f.contactPressure / 1e6), gap: round(f.gap * 1e3), state: f.closed ? 'in contact' : 'open' })),
          'Every layer interface with its contact pressure and gap. Open interfaces transmit nothing and provide no friction.'),
        tbl('Tensile armour friction and slip', [col('layer', 'Layer'), col('f', 'Friction per wire length', 'kN/m'), col('kcr', 'Slip-onset curvature', '1/m'), col('kfs', 'Full-slip curvature', '1/m'), col('sf', 'Max friction stress', 'MPa'), col('slip', 'Slip range', 'mm'), col('state', 'State at amplitude'), col('wear', 'Wear depth', 'mm')],
          rows.map((y) => ({ layer: y.L.name, f: round((y.b.frictionPerLength ?? 0) / 1e3), kcr: round(y.b.kappaCr ?? 0), kfs: round(y.b.kappaFullSlip ?? 0), sf: round((y.b.sigmaFrictionMax ?? 0) / 1e6), slip: round(y.slip.slipRange * 1e3), state: y.st.state, wear: round(y.wear.wearDepth * 1e3) })),
          'One row per tensile armour layer: available friction, the curvatures at which slip starts and is complete, and the resulting slip and wear.'),
      ],
      notes: [`Outer armour is in ${out.st.state} at the curvature amplitude ${fmt(amp)} 1/m (slip starts at ${fmt(bend.kappaCrMin ?? 0)} 1/m).`, `${gov.name} governs (${fmt(gov.uc)}).`, open ? `${open} interface(s) are open.` : 'All interfaces are in contact.', GENERIC],
      advice: [
        adv('engineering', sev(gov.uc), `Largest contact pressure ${fmt(maxContact / 1e6)} MPa; wire slip ${fmt(slipMax * 1e3)} mm per cycle; wear ${fmt(wearMax.wear.wearDepth * 1e3)} mm in ${fmt(inp.years)} years (utilisation ${fmt(gov.uc)}).`,
          gov.uc > 1 ? `The governing layer lasts about ${fmt(yearsToWear)} years: specify a thicker or more wear-resistant tape, or reduce the curvature range by ${fmt((1 - 1 / gov.uc) * 100)} %.` : 'Calibrate the friction coefficient against a full-scale bending-hysteresis test and the wear coefficient against a tape wear test at the calculated contact pressure.',
          `${STD_17B} §8; ${STD_17J} §5.3.1; Féret & Bournazel (1987); Archard (1953)`),
        adv('economic', 'info', `Friction dissipates ${fmt(loop ? loop.energyPerCycle : 0)} J per metre and cycle and gives ${fmt((loop ? loop.dampingRatio : 0) * 100)} % equivalent damping.`,
          Number.isFinite(yearsToWear) && yearsToWear < 2 * inp.years ? `Tape wear-out is predicted at about ${fmt(yearsToWear)} years: plan the first annulus inspection before year ${fmt(yearsToWear / 2)} rather than at the end of the design life.` : 'Wear does not drive the inspection interval; set it from the fatigue and annulus-integrity results instead.',
          'Risk-based inspection planning per DNV-RP-G101 / API RP 580',
          `Wear debris and a worn tape cannot be repaired in place — the consequence is early replacement of the pipe and a campaign of roughly 10 vessel days ≈ ${fmt(10 * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t). A more durable tape adds a few kilograms of polymer per kilometre.`),
      ],
    };
  },
  verify: [
    {
      // Slip range of a fully slipped wire: R²·cos²α·Δκ/sinα. Outer armour of the default stack: R = 76.2 + 7 + 8 + 8 + 1.5 +
      // 4 + 1.5 + 2 = 108.2 mm, α = 35°, Δκ = 0.02 1/m → 0.011707·0.67101·0.02/0.57358 = 2.739e-4 m = 0.2739 mm.
      name: 'Wire slip range of a fully slipped armour wire',
      source: 'Féret & Bournazel (1987) J. Offshore Mech. Arct. Eng. 109; hand calculation',
      inputs: { dKappa: 0.02, mu: 0.01 },
      expect: [{ kpi: 'slipRange', value: 0.2739, tol: 0.002 }],
    },
    {
      // Separation: with no load at all, an interface given an initial gap must stay open and carry no pressure,
      // while every other interface remains closed at zero pressure.
      name: 'Unloaded pipe with an initial gap: exactly one open interface, zero contact pressure',
      source: 'Unilateral contact condition (gap ≥ 0, pressure ≥ 0, gap × pressure = 0)',
      inputs: { Te: 0, pInt: 0, depth: 0, pAnn: 0, gap0: 0.5 },
      expect: [{ kpi: 'openInterfaces', value: 1, abs: 0 }, { kpi: 'maxContact', value: 0, abs: 1e-9 }],
    },
  ],
};
