// Nonlinear pushover: the 100-year environmental load set of one heading is scaled, with gravity held, until the
// jacket collapses. Reserve strength ratio, failure sequence, redundancy and an implied annual failure probability.
import { need, run, fmt, round, adv, fsCheck, check, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, PUSHOVER_CAVEAT, CONFIRM, BASIS_COST } from './_common.js';
import { stdInputs, buildJacket, stormCase, pushCases, pushCurve, nodeValuesFromMembers, plot3d, steelMass, groupSteel, foundNote, testFrame } from './_model.js';
import { pushoverAnalysis } from '../../engine/structural/pushover.js';
import { form } from '../../engine/reliability/form.js';

export default {
  id: 's01-pushover', suite: 1, title: 'Pushover and reserve strength ratio',
  summary: 'Ultimate lateral capacity of the jacket for one storm heading: pushover curve, sequence of member failures, reserve strength ratio, redundancy and the annual probability of collapse it implies.',
  covers: [11], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §12.6 (ultimate strength), §A.9.9 and §24.9 (RSR acceptance)', 'API RP 2SIM (2014) §8 (ultimate strength assessment)', 'API RP 2A-WSD 22nd ed. §C17.7', 'ISO 2394:2015 (reliability)'],
  models: ['Nonlinear pushover analysis (event-to-event)', 'Plastic hinges with tube N–M interaction', 'Brace buckling with constant post-buckling force', 'P–Δ geometric stiffness', 'Reserve strength ratio RSR = collapse load / 100-year load', 'First-order reliability method (FORM) with lognormal resistance and Gumbel annual load'],
  limits: `${PUSHOVER_CAVEAT} Joints, piles and soil are not failure elements here (linear pile-head springs): check joint and foundation capacity against the collapse load separately. Load pattern is fixed in shape (no wave-in-deck growth with height). ${MORISON_CAVEAT} ${CODE_CAVEAT}`,
  searchTerms: 'jacket pushover reserve strength ratio RSR ultimate strength',
  inputs: [
    ...stdInputs({ crit: { gG: 1.0, gE: 1.0 }, mesh: { nSub: 1, nPhase: 16 } }),
    { key: 'postBuck', label: 'Post-buckling residual force ratio', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Design criteria', role: 'criteria', uncertain: { dist: 'uniform', min: 0.2, max: 0.6 }, help: 'Compression kept by a buckled brace as a fraction of its buckling capacity. A simple post-buckling model: 0.3–0.5 is usual for D/t 30–50; confirm critical members with shell FE or test data.' },
    { key: 'pDelta', label: 'Include P–Δ', type: 'bool', default: true, group: 'Design criteria', role: 'control' },
    { key: 'rsrReq', label: 'Required reserve strength ratio', unit: '-', type: 'number', default: 1.85, min: 1, max: 4, group: 'Design criteria', role: 'criteria', help: `1.85 for manned (L1) new structures in ISO 19902 A.9.9; 1.6 / 1.5 / 0.8 by exposure category in API RP 2SIM assessments. ${CONFIRM}` },
    { key: 'covLoad', label: 'CoV of the annual maximum base shear', unit: '-', type: 'number', default: 0.35, min: 0.1, max: 1, group: 'Reliability', role: 'environment', help: 'About 0.3–0.4 for the North Sea, 0.6–0.8 for hurricane areas. Gumbel distribution fitted so that the 100-year load equals the reference base shear.' },
    { key: 'covR', label: 'CoV of the collapse resistance', unit: '-', type: 'number', default: 0.1, min: 0.02, max: 0.4, group: 'Reliability', role: 'material' },
    { key: 'biasR', label: 'Bias of the collapse resistance (mean / computed)', unit: '-', type: 'number', default: 1.1, min: 0.7, max: 1.5, group: 'Reliability', role: 'material', help: 'Mean yield over nominal and model bias of the frame analysis.' },
    { key: 'pfTarget', label: 'Target annual probability of collapse', unit: '-', type: 'number', default: 3e-5, min: 1e-7, max: 1e-2, group: 'Reliability', role: 'criteria', help: `Order of 3×10⁻⁵ per year for manned installations in ISO 19902 A.9.9.3. ${CONFIRM}` },
  ],
  mesh: { key: 'nPhase', label: 'Wave phase steps', levels: [8, 16, 32], qoi: 'baseShear' },
  compute(inp, ctx) {
    // SI conversion inside the shared model. Action factors default to 1.0: the reference is the unfactored 100-year load.
    const J = buildJacket(inp, ctx), st = stormCase(J, inp, { name: 'Storm' });
    need(st.baseShear > 1 || st.otm > 1, 'The environmental load pattern is zero: give a wave, current, wind or interface load to push the structure with.');
    pushCases(J, inp, st);
    const po = run(() => pushoverAnalysis(J.model, { pattern: 'Env', base: 'Gravity', pDelta: inp.pDelta, postBucklingRatio: inp.postBuck }), 'Pushover analysis');
    need(po.status !== 'collapse-under-base-load', 'The structure collapses under gravity loads alone: check member sizes and topsides weight.');
    need(Number.isFinite(po.ultimateLoadFactor) && po.ultimateLoadFactor > 0, 'The pushover did not reach a collapse load: the environmental load pattern is zero or the structure has no failure mode. Give a wave, current or interface load.');
    const rsr = po.ultimateLoadFactor, first = Number.isFinite(po.firstEventLoadFactor) ? po.firstEventLoadFactor : rsr, Vref = Math.abs(po.baseShearRef), Vult = Math.abs(po.ultimateBaseShear), c = pushCurve(po);
    const oOf = (e) => J.elems[e]?.o ?? e, ev = po.events.map((e) => ({ ...e, o: oOf(e.member) }));
    const failed = new Array(J.orig.length).fill(0);
    po.memberStates.forEach((s, e) => { if (s.axial !== 'elastic' || s.hingeA || s.hingeB) failed[oOf(e)] = 1; });
    const nFail = failed.reduce((s, v) => s + v, 0), byG = {}; ev.forEach((e) => { const g = J.orig[e.o].label; byG[g] = (byG[g] || 0) + 1; });
    // annual reliability: R lognormal (bias·RSR·V100), S Gumbel with its 99 % fractile at V100: F⁻¹(0.99) = μ(1 + 3.137·CoV)
    const muS = 1 / (1 + 3.137 * inp.covLoad), rel = run(() => form((x) => x.R - x.S, { R: { dist: 'lognormal', mean: inp.biasR * rsr, cov: inp.covR }, S: { dist: 'gumbel', mean: muS, std: muS * inp.covLoad } }), 'Reliability analysis (FORM)');
    const dUlt = Math.max(...po.curve.displacement.map(Math.abs)), steel = steelMass(J, inp), g1 = ev[0] ? J.orig[ev[0].o] : null, ok = rsr >= inp.rsrReq && rel.pf <= inp.pfTarget;
    const gSteel = g1 ? groupSteel(J, g1.group) : steel;
    const checks = [
      fsCheck('Reserve strength ratio', rsr, inp.rsrReq, '-', 'ISO 19902:2020 §A.9.9.3 / §24.9; API RP 2SIM §8.5', 'Collapse load of the structure over the 100-year environmental load. Below the required ratio the platform does not have the system reserve the standards assume beyond first member failure.'),
      fsCheck('First member failure above the 100-year load', first, 1, 'load factor', 'ISO 19902:2020 §12.6; design-level check', 'Load factor at which the first brace buckles, member yields or hinge forms. Below 1.0 a component fails inside the design storm, even if the system survives.'),
      check('Annual probability of collapse', rel.pf, inp.pfTarget, '1/yr', 'FORM; ISO 2394:2015; ISO 19902:2020 §A.9.9.3', 'Probability per year that the annual maximum storm load exceeds the collapse resistance, from the RSR and the stated load and resistance statistics. It turns the RSR into a risk that can be compared between regions with different load variability.'),
    ];
    return {
      kpis: [
        { key: 'rsr', label: 'Reserve strength ratio RSR', unit: '-', value: rsr, limit: inp.rsrReq, sense: 'min', explain: 'Ultimate environmental load the jacket can carry divided by the reference (100-year) load of this heading, with gravity held. The standards ask for about 1.5–1.85 on new manned platforms.' },
        { key: 'lfFirst', label: 'Load factor at first member failure', unit: '-', value: first, limit: 1, sense: 'min', explain: 'Multiple of the reference load at the first buckling, yielding or hinge event. The gap up to the RSR is the system reserve provided by redundancy.' },
        { key: 'redundancy', label: 'Redundancy factor (ultimate / first failure)', unit: '-', value: rsr / first, explain: 'Close to 1.0 the structure fails when its first member fails (K-braced or weakly redundant frames); X-braced jackets typically reach 1.2–1.5.' },
        { key: 'baseShear', label: 'Reference base shear', unit: 'MN', value: Vref / 1e6, explain: 'Horizontal environmental load of the reference case for this heading (wave, current and wind at the given action factor).' },
        { key: 'ultBaseShear', label: 'Ultimate base shear', unit: 'MN', value: Vult / 1e6, explain: 'Base shear at collapse: the lateral capacity of the structure for this load pattern.' },
        { key: 'dispUlt', label: 'Control displacement at collapse', unit: 'mm', value: dUlt * 1e3, explain: 'Displacement of the control node when the mechanism forms. A large value relative to the elastic displacement means ductile behaviour with warning before collapse.' },
        { key: 'nFailed', label: 'Members failed at collapse', unit: '-', value: nFail, explain: 'Number of members that have buckled, yielded or formed hinges when the mechanism develops.' },
        { key: 'beta', label: 'Annual reliability index β', unit: '-', value: rel.beta, explain: 'Reliability index of collapse in any one year. About 4.0 corresponds to 3×10⁻⁵ per year.' },
        { key: 'pfAnnual', label: 'Annual probability of collapse', unit: '1/yr', value: rel.pf, limit: inp.pfTarget, sense: 'max', explain: 'Probability that the yearly maximum load exceeds the collapse resistance. It depends strongly on the load CoV, so the same RSR is safer in the North Sea than in a hurricane region.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Pushover curve', x: { label: 'Control displacement', unit: 'mm' }, y: { label: 'Base shear', unit: 'MN' }, series: [{ name: 'Pushover', x: c.d, y: c.v }, { name: '100-year load', x: [0, Math.max(...c.d)], y: [round(Vref / 1e6), round(Vref / 1e6)], dash: true }, { name: `Required (RSR ${inp.rsrReq})`, x: [0, Math.max(...c.d)], y: [round(Vref * inp.rsrReq / 1e6), round(Vref * inp.rsrReq / 1e6)], dash: true }, { name: 'Member failures', x: ev.slice(0, 60).map((e) => round(Math.abs(e.displacement) * 1e3)), y: ev.slice(0, 60).map((e) => round(e.loadFactor * Vref / 1e6)) }],
          explain: 'Base shear against displacement as the environmental load is scaled. The knees are member failures; the peak is the ultimate capacity. It must clear the dashed required line, and a long plateau after the first knee shows a ductile, redundant frame.' },
        plot3d(J, { title: 'Failed members and displaced shape at collapse', values: nodeValuesFromMembers(J, failed), label: 'Member failed (1) or intact (0)', unit: '-', u: po.displacements, explain: 'Nodes of members that have buckled, yielded or hinged at the collapse load, on the displaced shape. The failed members outline the collapse mechanism: usually one bay of braces followed by portal hinges in the legs.' }),
        { type: 'bar', title: 'Failure events by member group', categories: Object.keys(byG).length ? Object.keys(byG) : ['none'], series: [{ name: 'Events', y: Object.keys(byG).length ? Object.values(byG) : [0] }], y: { label: 'Events', unit: '-' }, explain: 'Number of plastic events in each member group up to collapse. Failures concentrated in braces indicate a brace-dominated mechanism; leg hinges indicate portal action.' },
      ],
      tables: [
        { title: 'Failure sequence', columns: [{ key: 'k', label: 'No.' }, { key: 'lf', label: 'Load factor' }, { key: 'V', label: 'Base shear', unit: 'MN' }, { key: 'm', label: 'Member' }, { key: 'g', label: 'Group' }, { key: 'type', label: 'Event' }, { key: 'd', label: 'Displacement', unit: 'mm' }],
          rows: ev.slice(0, 60).map((e, k) => ({ k: k + 1, lf: round(e.loadFactor), V: round(e.loadFactor * Vref / 1e6), m: e.o, g: J.orig[e.o].label, type: e.type, d: round(Math.abs(e.displacement) * 1e3) })), explain: 'Order in which members fail as the load increases, with the load factor of each event. The first rows identify the members to strengthen if the RSR is short.' },
        { title: 'Reliability summary', columns: [{ key: 'c', label: 'Quantity' }, { key: 'v', label: 'Value' }], rows: [{ c: 'Mean resistance / 100-year load', v: round(inp.biasR * rsr) }, { c: 'Mean annual maximum load / 100-year load', v: round(muS) }, { c: 'Design point: resistance', v: round(rel.designPoint.R) }, { c: 'Design point: load', v: round(rel.designPoint.S) }, { c: 'Importance of load uncertainty (α²)', v: round(rel.importance.S) }, { c: 'Importance of resistance uncertainty (α²)', v: round(rel.importance.R) }],
          explain: 'Inputs and design point of the FORM analysis in units of the 100-year load. The importance factors show that the storm load uncertainty, not steel strength, controls the collapse risk.' },
      ],
      notes: [`Heading ${fmt(st.beta)}°: collapse at ${fmt(rsr)} × the reference load (${fmt(Vult / 1e6)} MN against ${fmt(Vref / 1e6)} MN), status "${po.status}", ${ev.length} event(s), ${nFail} member(s) failed.`,
        g1 ? `First failure: ${ev[0].type} of member ${ev[0].o} (${g1.label.toLowerCase()}) at load factor ${fmt(first)}; redundancy factor ${fmt(rsr / first)}.` : 'No member failed before the mechanism formed.',
        `Annual reliability index ${fmt(rel.beta)} (Pf ${rel.pf.toExponential(2)} per year) for a load CoV of ${fmt(inp.covLoad)}. Wave by ${st.wf.label}.`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `Reserve strength is adequate for this heading: RSR ${fmt(rsr)} ≥ ${fmt(inp.rsrReq)}, annual collapse probability ${rel.pf.toExponential(1)}.`, `Repeat for the other headings (the lowest RSR governs) and check joints and piles against the collapse base shear of ${fmt(Vult / 1e6)} MN, since they are not failure elements here. ${rsr > 1.5 * inp.rsrReq ? 'The reserve is large: bracing or leg wall could be optimised if fatigue and installation allow.' : ''}`, 'ISO 19902:2020 §12.6, §A.9.9; API RP 2SIM (2014) §8', `System reserve replaces steel: an X-braced frame reaches the required RSR with lighter members than a non-redundant one. ${g1 ? `The ${g1.label.toLowerCase()} that fail first hold ${steelCo2(gSteel)}` : `The jacket holds ${steelCo2(steel)}`}`)
          : adv('engineering', 'fail', `Reserve strength is not adequate: RSR ${fmt(rsr)} against ${fmt(inp.rsrReq)}, annual collapse probability ${rel.pf.toExponential(1)} against ${inp.pfTarget.toExponential(1)}.`, g1 ? `Strengthen the members that fail first (${g1.label.toLowerCase()}, member ${ev[0].o}): capacity must rise by about ${fmt(Math.max(inp.rsrReq / rsr - 1, 0.05) * 100)} %; for an existing platform consider load reduction (marine growth removal, deck raising) before grouting or clamps.` : 'Increase leg and brace capacity or reduce the environmental load.', 'ISO 19902:2020 §A.9.9; API RP 2SIM (2014) §8', `Strengthening the first-failing group by that ratio is about ${steelCo2(gSteel * Math.max(inp.rsrReq / rsr - 1, 0.05))}`),
        adv('economic', 'info', `Each 0.1 of RSR is worth ${fmt(0.1 * Vref / 1e6)} MN of collapse capacity; the redundancy factor is ${fmt(rsr / first)}.`, 'For life extension, a demonstrated RSR above the requirement is the cheapest way to accept local over-utilisation: analysis hours instead of offshore clamps. Spend on reducing load uncertainty (site-specific metocean) before adding steel.', `${BASIS_COST}; risk-based assessment per API RP 2SIM / ISO 19902 §24`),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Cantilever tube pushed at the tip: collapse when the base section becomes fully plastic, H = Mp/L.
    // Zp = (D³ − (D − 2t)³)/6 = (1.5³ − 1.41³)/6 = 0.095296 m³, Mp = 355e6·0.095296 = 33.830 MN·m, L = 95 m,
    // reference load 100 kN → RSR = 33.830e6/(1e5·95) = 3.561.
    { name: 'Plastic collapse of a cantilever: λ = Mp/(H·L)', source: 'Plastic hinge theory (Neal, The Plastic Methods of Structural Analysis, 1977, §2)',
      inputs: { geom: testFrame('column'), H: 0, Uc: 0, Uw: 0, windA: 0, swCont: 0, Wtop: 0, FxTop: 100, found: 'fixed', heading: 0, nSub: 1 },
      expect: [{ kpi: 'rsr', value: 3.561, tol: 0.005 }, { kpi: 'baseShear', value: 0.1, tol: 1e-6 }] },
  ],
};
