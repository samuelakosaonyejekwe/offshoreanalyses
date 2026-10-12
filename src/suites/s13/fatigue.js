// Fatigue of a dynamic cable under cyclic tension and curvature: armour-wire stress ranges (stick–slip model),
// strain-based fatigue of the metallic sheath, conductor strand stress ranges, Miner summation over a long-term
// load spectrum and an optional rainflow-counted tension history.
import { armourStressRange, snCycles, genericArmourSN, minerDamage } from '../../engine/xsection/flexfatigue.js';
import { helixContactPressure, criticalCurvature, wireStrain } from '../../engine/xsection/helix.js';
import { rainflowFourPoint, rangeHistogram } from '../../engine/fatigue/rainflow.js';
import { cableInputs, buildMech, adv, check, worst, sev, fmt, round, run, need, logspace, STD_CABLE, LIMIT_GENERIC, CONFIRM, vesselCo2Text } from './_common.js';

const DEFAULT_SPECTRUM = [
  { name: 'Calm (Hs < 1 m)', meanT: 55, dT: 2, dK: 0.0012, cycles: 2.2e6 },
  { name: 'Moderate (1–2 m)', meanT: 56, dT: 5, dK: 0.003, cycles: 1.8e6 },
  { name: 'Fresh (2–3 m)', meanT: 57, dT: 9, dK: 0.005, cycles: 7.0e5 },
  { name: 'Rough (3–5 m)', meanT: 58, dT: 14, dK: 0.0085, cycles: 2.0e5 },
  { name: 'Storm (5–8 m)', meanT: 60, dT: 22, dK: 0.013, cycles: 2.0e4 },
  { name: 'Extreme (> 8 m)', meanT: 62, dT: 32, dK: 0.019, cycles: 6.0e2 },
];

export default {
  id: 's13-fatigue', suite: 13, title: 'Cable fatigue and cyclic loading',
  summary: 'Miner fatigue damage and life of armour wires, metallic sheath and conductors of a dynamic cable for a long-term spectrum of tension and curvature ranges.',
  covers: [4], fidelity: 'analytical',
  standards: ['CIGRE TB 862 (2022) §6 (fatigue analysis and testing of dynamic cables)', `${STD_CABLE}`, 'DNV-RP-C203 (2021) §2.2 (Palmgren–Miner), §2.4 (S–N format)', 'ASTM E1049-85 (rainflow counting)'],
  models: ['S–N curves and Palmgren–Miner damage', 'Armour-wire fatigue: stick–slip friction stress and local wire bending (helical-tendon model)', 'Strain-based fatigue of the metallic sheath (ε–N power law)', 'Conductor strand stress range', 'Goodman mean-stress correction', 'Rainflow cycle counting of a tension history'],
  limits: `${LIMIT_GENERIC} Curvature ranges are taken as given (outside the bend stiffener, no local stiffener model); fretting between strands, corrosion fatigue beyond the S–N knock-down and insulation fatigue are not modelled. The default S–N and ε–N parameters are indicative placeholders — fatigue design needs test curves of the actual wires, sheath alloy and conductor (CIGRE TB 862).`,
  publishes: ['life.fatigueYears', 'cable.fatigueLife'],
  searchTerms: 'dynamic submarine cable fatigue armour wire lead sheath Miner',
  inputs: [
    ...cableInputs({ mech: true }),
    { key: 'spectrum', label: 'Long-term load spectrum at the critical section', type: 'table', group: 'Cyclic loading', role: 'load',
      columns: [{ key: 'name', label: 'Sea-state class', type: 'text' }, { key: 'meanT', label: 'Mean tension', unit: 'kN', type: 'number' }, { key: 'dT', label: 'Tension range', unit: 'kN', type: 'number' }, { key: 'dK', label: 'Curvature range', unit: '1/m', type: 'number' }, { key: 'cycles', label: 'Cycles per year', type: 'number' }],
      default: DEFAULT_SPECTRUM.map((r) => ({ ...r })), help: 'Ranges (peak to trough) per class from the dynamic analysis over the wave scatter diagram, at the hot spot (hang-off, hog bend or touchdown).' },
    { key: 'history', label: 'Tension time history (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Time', xUnit: 's', yLabel: 'Tension', yUnit: 'kN', group: 'Cyclic loading', role: 'load',
      help: 'Measured or simulated tension record. Its rainflow cycles, scaled to a year, are added to the armour and conductor damage (no curvature).' },
    { key: 'D0', label: 'Damage already accumulated', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'Cyclic loading', role: 'ic', help: 'Initial condition: Miner sum consumed so far (installation, earlier service, from inspection or monitoring).' },
    { key: 'ends', label: 'Torsional restraint', type: 'select', default: 'fixed', group: 'Cyclic loading', role: 'bc', options: [{ value: 'fixed', label: 'Ends restrained against rotation' }, { value: 'free', label: 'Torque-free' }] },
    { key: 'mu', label: 'Friction coefficient armour / bedding', unit: '-', type: 'number', default: 0.2, min: 0.02, max: 0.6, group: 'Cyclic loading', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'snModel', label: 'Armour S–N curve', type: 'select', default: 'custom', group: 'Fatigue resistance', role: 'material', options: [{ value: 'custom', label: 'Reference range at 2·10⁶ cycles and slope (seawater)' }, { value: 'generic-air', label: 'Generic textbook estimate from tensile strength (air, with endurance limit)' }] },
    { key: 'snRef', label: 'Armour stress range at 2·10⁶ cycles', unit: 'MPa', type: 'number', default: 180, min: 20, max: 800, group: 'Fatigue resistance', role: 'material', uncertain: { dist: 'lognormal', cov: 0.15 }, showIf: { key: 'snModel', in: ['custom'] }, help: 'Indicative for galvanised armour wire in sea water; replace with the wire test curve.' },
    { key: 'snM', label: 'Armour S–N slope m', unit: '-', type: 'number', default: 4, min: 2.5, max: 12, group: 'Fatigue resistance', role: 'material', showIf: { key: 'snModel', in: ['custom'] } },
    { key: 'epsRef', label: 'Sheath strain range at 10⁶ cycles', unit: '%', type: 'number', default: 0.1, min: 0.005, max: 2, group: 'Fatigue resistance', role: 'material', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'ε–N power law N = 10⁶·(ε_ref/Δε)^m of the sheath alloy. Indicative for lead alloy; copper and aluminium sheaths are far better — use test data.' },
    { key: 'epsM', label: 'Sheath ε–N slope m', unit: '-', type: 'number', default: 3.5, min: 1.5, max: 10, group: 'Fatigue resistance', role: 'material' },
    { key: 'condRef', label: 'Conductor stress range at 2·10⁶ cycles', unit: 'MPa', type: 'number', default: 60, min: 5, max: 300, group: 'Fatigue resistance', role: 'material', help: 'Indicative for annealed copper strand including fretting; aluminium is lower. Use conductor fatigue test data.' },
    { key: 'condM', label: 'Conductor S–N slope m', unit: '-', type: 'number', default: 5, min: 2.5, max: 12, group: 'Fatigue resistance', role: 'material' },
    { key: 'strandD', label: 'Conductor strand diameter', unit: 'mm', type: 'number', default: 3, min: 0.5, max: 6, group: 'Fatigue resistance', role: 'geometry' },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.0, min: 1, max: 3, group: 'Fatigue resistance', role: 'criteria', help: 'E.g. at the hang-off termination or wire welds.' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 10, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: `10 for non-inspectable dynamic cable components. ${CONFIRM}` },
    { key: 'life', label: 'Design life', unit: 'years', type: 'number', default: 25, min: 1, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
  ],
  compute(inp, ctx) {
    const b = buildMech(inp, ctx), cab = b.cab, arm = b.armourLayers;
    const rows = (inp.spectrum || []).filter((r) => r.cycles > 0);
    need(rows.length > 0, 'The load spectrum is empty: add at least one class with cycles per year.');
    for (const r of rows) need(r.dT >= 0 && r.dK >= 0 && r.meanT >= 0, `Spectrum row '${r.name}': mean tension, tension range and curvature range cannot be negative.`);
    const K11 = cab.EAfixed, epsPerN = inp.ends === 'fixed' ? 1 / K11 : 1 / cab.EA;
    const twistPerN = inp.ends === 'fixed' ? 0 : -cab.coupling / (cab.EAfixed * cab.GJfixed - cab.coupling ** 2);
    const curveA = inp.snModel === 'generic-air' ? genericArmourSN({ uts: (inp.armUts ?? 500) * 1e6 }) : { m: inp.snM, a: 2e6 * (inp.snRef * 1e6) ** inp.snM };
    const curveC = { m: inp.condM, a: 2e6 * (inp.condRef * 1e6) ** inp.condM };
    // armour: governing layer is evaluated bin by bin (outer layer has the largest radius, inner the largest contact pressure)
    const wireStressPerN = (L) => L.E * wireStrain({ alpha: L.alpha, R: L.R, lay: L.lay, epsAxial: epsPerN, twist: twistPerN }).strain;
    const armBin = (r) => {
      let best = { range: 0, mean: 0, slipped: false, layer: '–' };
      arm.forEach((L, i) => {
        const d = Math.sqrt(4 * L.wireArea / Math.PI), sPerN = wireStressPerN(L);
        const qSelf = helixContactPressure({ n: L.nWires, wireForce: Math.max(sPerN * r.meanT * 1e3, 0) * L.wireArea, R: L.R, alpha: L.alpha }).pressure;
        const qOut = arm.slice(i + 1).reduce((a, Lo) => a + helixContactPressure({ n: Lo.nWires, wireForce: Math.max(wireStressPerN(Lo) * r.meanT * 1e3, 0) * Lo.wireArea, R: Lo.R, alpha: Lo.alpha }).pressure, 0);
        const f = criticalCurvature({ E: L.E, A: L.wireArea, R: L.R, alpha: L.alpha, n: L.nWires, mu: inp.mu, qOuter: qOut, qInner: qOut + qSelf }).frictionPerLength;
        const sr = r.dK > 0 ? armourStressRange({ E: L.E, A: L.wireArea, R: L.R, alpha: L.alpha, thickness: d, width: d, frictionPerLength: f, kappaRange: r.dK }) : { governing: 0, slipped: false };
        const range = Math.abs(sPerN) * r.dT * 1e3 + sr.governing;
        if (range >= best.range) best = { range, mean: sPerN * r.meanT * 1e3, slipped: sr.slipped, layer: L.name };
      });
      return best;
    };
    const sheath = cab.layers.find((L) => L.type === 'sheath'), rSheath = sheath ? sheath.ro : null;
    const cond = cab.layers.find((L) => L.type === 'conductor'), cosC = cond.layAngle ? Math.cos(cond.layAngle) : 1;
    const bins = rows.map((r) => {
      const a = arm.length ? armBin(r) : { range: 0, mean: 0, slipped: false, layer: '–' };
      const epsS = rSheath ? r.dK * rSheath : 0;
      const sC = cond.E * (epsPerN * r.dT * 1e3 * cosC * cosC + r.dK * inp.strandD / 2e3);
      return { r, a, epsS, sC, Ns: epsS > 0 ? 1e6 * (inp.epsRef / 100 / epsS) ** inp.epsM : Infinity };
    });
    // optional tension history → rainflow → cycles per year (tension only)
    const hist = inp.history, extraA = [], extraC = [];
    let histCycles = 0;
    if (hist?.y?.length > 3) {
      need(hist.x.length === hist.y.length, 'The tension history needs the same number of time and tension values.');
      const dur = hist.x[hist.x.length - 1] - hist.x[0];
      need(dur > 0, 'The tension history must span a positive duration.');
      const rf = rainflowFourPoint(Array.from(hist.y)), hg = rangeHistogram(rf.cycles ?? rf, { bins: 24 }), scale = 31557600 / dur;
      for (let i = 0; i < hg.counts.length; i++) if (hg.counts[i] > 0) {
        const dT = hg.centres[i] * 1e3, n = hg.counts[i] * scale; histCycles += n;
        if (arm.length) extraA.push({ stressRange: Math.max(...arm.map((L) => Math.abs(wireStressPerN(L)))) * dT, cycles: n });
        extraC.push({ stressRange: cond.E * epsPerN * dT * cosC * cosC, cycles: n });
      }
    }
    const uts = (inp.armUts ?? 500) * 1e6;
    const mdA = arm.length ? run(() => minerDamage({ bins: [...bins.map((q) => ({ stressRange: q.a.range, cycles: q.r.cycles, meanStress: Math.max(q.a.mean, 0) })), ...extraA], curve: curveA, uts, meanStressMethod: 'goodman', scf: inp.scf, dff: inp.dff, exposureYears: 1 }), 'Armour fatigue') : { damagePerYear: 0, bins: bins.map(() => ({ damage: 0, cyclesToFailure: Infinity })) };
    const mdC = run(() => minerDamage({ bins: [...bins.map((q) => ({ stressRange: q.sC, cycles: q.r.cycles })), ...extraC], curve: curveC, scf: inp.scf, dff: inp.dff, exposureYears: 1 }), 'Conductor fatigue');
    const dS = bins.reduce((a, q) => a + q.r.cycles / q.Ns, 0);
    const comp = [{ name: 'Armour wires', d: mdA.damagePerYear }, { name: 'Metallic sheath', d: dS }, { name: 'Conductor strands', d: mdC.damagePerYear }];
    const lifeOf = (d) => (d > 0 ? Math.min((1 - inp.D0) / d, 1e6) : 1e6);
    const govC = comp.reduce((a, c) => (c.d > a.d ? c : a), comp[0]);
    const lifeChar = lifeOf(govC.d), lifeDes = lifeChar / inp.dff;
    const checks = comp.map((c) => check(`${c.name}: Miner sum × DFF over the design life`, inp.D0 + c.d * inp.life * inp.dff, 1, '-', 'DNV-RP-C203 §2.2 (Miner); CIGRE TB 862 §6 (design fatigue factor)',
      c.name.startsWith('Armour') ? 'Accumulated armour-wire damage with the design fatigue factor. Above 1.0 wire breaks are expected within the design life, starting at the hot spot, reducing tensile capacity and torque balance.'
        : c.name.startsWith('Metallic') ? 'Accumulated strain-fatigue damage of the radial water barrier. A cracked sheath lets water reach the insulation, leading to water treeing and electrical breakdown.'
          : 'Accumulated damage of the conductor strands (bending and fretting). Broken strands raise local resistance and temperature and end in an open-circuit fault.'));
    const gov = worst(checks);
    const sAxis = logspace(10e6, 800e6, 40), nA = sAxis.map((s) => snCycles({ stressRange: s, curve: curveA }).cycles).map((n) => (Number.isFinite(n) ? n : 1e14));
    const maxRange = Math.max(...bins.map((q) => q.a.range));
    return {
      kpis: [
        { key: 'lifeDesign', label: 'Fatigue life with design factor', value: round(lifeDes), unit: 'years', limit: inp.life, sense: 'min', explain: 'Years until the governing component reaches a Miner sum of 1/DFF, allowing for damage already accumulated. It must exceed the design life.' },
        { key: 'lifeChar', label: 'Characteristic fatigue life (no factor)', value: round(lifeChar), unit: 'years', explain: 'Life at a Miner sum of 1.0 on the characteristic curves; the ratio to the design life is the available fatigue factor.' },
        { key: 'damageArmour', label: 'Armour damage per year', value: round(mdA.damagePerYear, 4), unit: '1/yr', explain: 'Unfactored Miner damage rate of the most loaded armour layer from tension and curvature cycling.' },
        { key: 'damageSheath', label: 'Sheath damage per year', value: round(dS, 4), unit: '1/yr', explain: 'Unfactored strain-fatigue damage rate of the metallic sheath; lead alloys usually govern dynamic cable design.' },
        { key: 'damageCond', label: 'Conductor damage per year', value: round(mdC.damagePerYear, 4), unit: '1/yr', explain: 'Unfactored damage rate of the conductor strands.' },
        { key: 'lifeArmour', label: 'Armour life (no factor)', value: round(lifeOf(mdA.damagePerYear)), unit: 'years', explain: 'Characteristic armour-wire life; 10⁶ means no damage is predicted.' },
        { key: 'lifeSheath', label: 'Sheath life (no factor)', value: round(lifeOf(dS)), unit: 'years', explain: 'Characteristic life of the metallic sheath from the ε–N power law.' },
        { key: 'lifeCond', label: 'Conductor life (no factor)', value: round(lifeOf(mdC.damagePerYear)), unit: 'years', explain: 'Characteristic life of the conductor strands.' },
        { key: 'maxArmourRange', label: 'Largest armour stress range', value: round(maxRange / 1e6), unit: 'MPa', explain: 'Stress range of the most severe class: axial range from tension plus friction and local bending range from curvature.' },
        { key: 'maxSheathStrain', label: 'Largest sheath strain range', value: round(Math.max(...bins.map((q) => q.epsS)) * 100, 4), unit: '%', explain: 'Curvature range × sheath radius of one core in the most severe class.' },
        { key: 'cycles', label: 'Cycles per year', value: round(rows.reduce((a, r) => a + r.cycles, 0) + histCycles, 5), unit: '-', explain: 'Total number of load cycles a year in the spectrum (plus the rainflow-counted history).' },
        { key: 'uc', label: 'Fatigue utilisation', value: round(gov.uc), unit: '-', limit: 1, sense: 'max', explain: 'Largest factored Miner sum over the design life; below 1.0 the fatigue criterion is met.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Armour S–N curve and spectrum classes', x: { label: 'Cycles', unit: '-' }, y: { label: 'Stress range', unit: 'MPa' },
          series: [{ name: 'S–N curve', x: nA.map((n) => Math.log10(n)), y: sAxis.map((s) => s / 1e6) }, { name: 'Spectrum classes (log₁₀ cycles in design life)', x: bins.map((q) => Math.log10(q.r.cycles * inp.life)), y: bins.map((q) => q.a.range / 1e6) }],
          explain: 'Horizontal axis is log₁₀ of cycles. Classes plotted close to or right of the curve consume life; classes far to the left are harmless. The few severe classes and the many mild ones can contribute equally.' },
        { type: 'bar', title: 'Damage per year by sea-state class', categories: rows.map((r) => r.name),
          series: [{ name: 'Armour', y: bins.map((q, i) => mdA.bins[i].damage) }, { name: 'Sheath', y: bins.map((q) => q.r.cycles / q.Ns) }, { name: 'Conductor', y: bins.map((q, i) => mdC.bins[i].damage) }],
          explain: 'Which sea states consume the fatigue life of each component. If moderate seas dominate, the remedy is the configuration (curvature in everyday conditions) rather than storm survival.' },
        { type: 'line', title: 'Accumulated factored damage over time', x: { label: 'Time', unit: 'years' }, y: { label: 'Miner sum × DFF', unit: '-' },
          series: [...comp.map((c) => ({ name: c.name, x: [0, inp.life], y: [inp.D0, inp.D0 + c.d * inp.dff * inp.life] })), { name: 'Limit', x: [0, inp.life], y: [1, 1], dash: true }],
          explain: 'Linear damage growth of each component including the design fatigue factor. A line crossing 1.0 before the end of the design life fails the criterion at that year.' },
      ],
      tables: [{ title: 'Fatigue by spectrum class', explain: 'One row per class: armour stress range (and whether the wires slip), sheath strain range, conductor stress range and the unfactored damage each contributes per year.',
        columns: [{ key: 'name', label: 'Class' }, { key: 'cycles', label: 'Cycles/yr' }, { key: 'sA', label: 'Armour range', unit: 'MPa' }, { key: 'slip', label: 'Wires slip' }, { key: 'eS', label: 'Sheath strain range', unit: '%' }, { key: 'sC', label: 'Conductor range', unit: 'MPa' }, { key: 'dA', label: 'Armour damage', unit: '1/yr' }, { key: 'dS', label: 'Sheath damage', unit: '1/yr' }, { key: 'dC', label: 'Conductor damage', unit: '1/yr' }],
        rows: bins.map((q, i) => ({ name: q.r.name, cycles: q.r.cycles, sA: round(q.a.range / 1e6), slip: q.a.slipped ? 'yes' : 'no', eS: round(q.epsS * 100, 3), sC: round(q.sC / 1e6), dA: round(mdA.bins[i].damage, 3), dS: round(q.r.cycles / q.Ns, 3), dC: round(mdC.bins[i].damage, 3) })) }],
      notes: [
        `${govC.name} govern with ${fmt(govC.d)} damage per year (unfactored): characteristic life ${fmt(lifeChar)} years, ${fmt(lifeDes)} years with DFF ${inp.dff}.`,
        inp.snModel === 'generic-air' ? 'The generic air S–N curve has an endurance limit and is not valid in sea water: use it only for comparison.' : `Armour S–N: ${inp.snRef} MPa at 2·10⁶ cycles, slope ${inp.snM} (indicative).`,
        histCycles > 0 ? `The tension history adds ${fmt(histCycles)} rainflow cycles per year.` : 'No tension history supplied: damage comes from the spectrum table only.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} = ${fmt(gov.uc)}; design fatigue life ${fmt(lifeDes)} years against ${inp.life} years required.`,
          gov.uc > 1 ? `Cut the governing range by ${fmt((1 - gov.uc ** (-1 / (govC.name.startsWith('Metallic') ? inp.epsM : govC.name.startsWith('Armour') ? inp.snM : inp.condM))) * 100)} % — through a longer bend stiffener, more buoyancy in the wave, or ${govC.name.startsWith('Metallic') ? 'a copper or corrugated sheath instead of lead' : 'smaller strands / wires'} — and re-run.`
            : 'Confirm the S–N and ε–N curves with full-scale flex-fatigue tests of the cable (CIGRE TB 862) and repeat the analysis with the scatter-diagram ranges of the final configuration.',
          'CIGRE TB 862 (2022) §6; DNV-RP-C203 (2021) §2.2; DNV-RP-0360'),
        adv('economic', gov.uc > 1 ? 'warn' : 'info', `At the present loading the dynamic section would need replacing after about ${fmt(Math.min(lifeDes, 999))} years on the design basis (${fmt(Math.min(lifeChar, 9999))} years characteristic).`,
          lifeDes < inp.life ? `Plan ${Math.ceil(inp.life / Math.max(lifeDes, 0.1)) - 1} change-out(s) of the dynamic cable within the ${inp.life}-year life, or invest in the configuration now: each change-out costs a spare dynamic length, a cable-lay spread for roughly 5–10 days and the lost production of the connected turbine(s).` : `The margin (${fmt(lifeDes / inp.life)} × design life) can be traded against a lighter armour or shorter buoyancy section; alternatively keep it and reduce the inspection frequency of the dynamic section.`,
          'Whole-life costing ISO 15686-5; risk-based inspection planning DNV-RP-G101 principles',
          `A change-out campaign of about 7 vessel days: ${vesselCo2Text(7)} A longer-lived dynamic section avoids both the vessel emissions and the remanufacture of the cable.`),
      ],
      publish: { 'life.fatigueYears': Math.min(lifeDes, 1e6), 'cable.fatigueLife': Math.min(lifeDes, 1e6) },
      quantities: {},
    };
  },
  verify: [
    // Strain-based sheath fatigue, two-block loading. Single core: sheath outer radius = 14.927 + 1 + 9 + 1 + 1.8 = 27.727 mm.
    // Block 1: Δκ = 0.01 1/m → Δε = 2.7727e-4; N1 = 1e6·(1e-3/2.7727e-4)^3.5 = 8.909e7; n1 = 1e6 → d1 = 0.011225.
    // Block 2: Δκ = 0.02 1/m → Δε = 5.5454e-4; N2 = 1e6·(1.8033)^3.5 = 7.875e6; n2 = 1e5 → d2 = 0.012698.
    // Miner sum per year = 0.023923 → life = 41.80 years.
    { name: 'Two-block Miner sum on an ε–N power law (sheath strain = Δκ·r)',
      source: 'Palmgren–Miner rule, DNV-RP-C203 Eq. 2.2.1; hand calculation',
      inputs: { cores: '1', spectrum: [{ name: 'A', meanT: 20, dT: 0, dK: 0.01, cycles: 1e6 }, { name: 'B', meanT: 20, dT: 0, dK: 0.02, cycles: 1e5 }] },
      expect: [{ kpi: 'damageSheath', value: 0.023923, tol: 0.005 }, { kpi: 'lifeSheath', value: 41.80, tol: 0.005 }] },
  ],
};
