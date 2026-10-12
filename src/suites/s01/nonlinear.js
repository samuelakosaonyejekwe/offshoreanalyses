// Nonlinear analysis of the jacket at and above the design load: geometric nonlinearity (second-order P–Δ / P–δ on
// the geometric stiffness) and material nonlinearity (plastic hinges and brace buckling, event to event).
import { need, run, fmt, round, adv, check, fsCheck, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, PUSHOVER_CAVEAT, BASIS_MEMBER, BASIS_COST } from './_common.js';
import { stdInputs, buildJacket, stormCase, loadCase, solve, memberRows, byGroup, displacements, nodeValuesFromMembers, plot3d, memberCols, memberTableRows, groupChecks, linearBuckling, steelMass, groupSteel, foundNote, testFrame } from './_model.js';
import { pDeltaAnalysis } from '../../engine/structural/framebuckling.js';
import { pushoverAnalysis } from '../../engine/structural/pushover.js';

export default {
  id: 's01-nonlinear', suite: 1, title: 'Nonlinear structural analysis',
  summary: 'Second-order (P–Δ) amplification of displacements and member utilisation, elastic critical load factor, and the elastic-plastic response up to a target multiple of the design storm load.',
  covers: [3], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §12.3 (nonlinear analysis), §12.6', 'EN 1993-1-1:2005 §5.2 (second-order effects, αcr)', 'API RP 2A-WSD 22nd ed. §C17.7'],
  models: ['Geometric nonlinear analysis (geometric stiffness iteration, P–Δ and P–δ)', 'Material nonlinear analysis with plastic hinges (elastic–perfectly plastic)', 'Elastoplastic tube interaction surface cos(πN/2Np) − M/Mp', 'Column buckling with constant post-buckling force', 'Eigenvalue buckling', 'Event-to-event incremental solution'],
  limits: `Moderate displacements: the geometry is not updated (no corotational large-displacement formulation), there is no strain hardening, no local shell buckling and no joint failure. ${PUSHOVER_CAVEAT} ${MORISON_CAVEAT} ${CODE_CAVEAT}`,
  searchTerms: 'jacket second-order analysis plastic hinge geometric nonlinearity',
  inputs: [
    ...stdInputs({ mesh: { nSub: 2, nPhase: 16 } }),
    { key: 'lfTarget', label: 'Target load factor on the factored environmental load', unit: '-', type: 'number', default: 2, min: 1, max: 10, group: 'Solver', role: 'control', help: 'The elastic-plastic trace is followed up to this multiple of the factored storm load with gravity held constant (use the pushover module for the full collapse load).' },
    { key: 'postBuck', label: 'Post-buckling residual force ratio', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Solver', role: 'criteria', help: 'Compression force kept by a buckled member as a fraction of its buckling capacity: a simple post-buckling model, typically 0.3–0.5 for jacket braces.' },
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'amp' },
  compute(inp, ctx) {
    // SI conversion inside the shared model. Gravity is held at its factored value; the environmental set is scaled.
    const J = buildJacket(inp, ctx), st = stormCase(J, inp, { name: 'Storm' });
    const G0 = loadCase('Gravity', J, inp, { betaDeg: st.beta, fG: inp.gG, fE: 0 }), E0 = loadCase('Env', J, inp, { betaDeg: st.beta, fG: 0, fE: inp.gE, env: st.env, sw: false, W: 0 });
    const lin = solve(J, [st.lc, G0, E0]).Storm, bk = linearBuckling(J, 'Storm', 1), aCr = bk.factors[0];
    need(aCr > 1.02, `The factored load set is at or above the elastic critical load of the frame (αcr = ${fmt(aCr)}): a second-order solution does not exist. Stiffen the slender members (see the stability module for the buckling mode) or reduce the load.`);
    const pd = run(() => pDeltaAnalysis(J.model, { loadCase: 'Storm' }), 'Second-order analysis (the design load may exceed the elastic critical load)');
    need(pd.converged, 'The second-order analysis did not converge: the design load is at or above the elastic critical load of the frame.');
    const rowsL = memberRows(J, lin, inp, { name: 'first order' }), rowsN = memberRows(J, pd.result, inp, { name: 'second order' });
    const dL = displacements(J, lin), dN = displacements(J, pd.result);
    const ucL = Math.max(...rowsL.map((r) => r.uc)), ucN = Math.max(...rowsN.map((r) => r.uc)), groups = byGroup(rowsN);
    ctx.progress?.(0.5, 'Elastic-plastic trace');
    // The hinge model works per member (capacities with K and the buckling length), so the plastic trace runs on one
    // element per member with the same loads; the refined mesh is kept for the second-order elastic results.
    const J1 = J.nSub > 1 ? buildJacket(inp, ctx, { nSub: 1 }) : J, st1 = J1 === J ? st : stormCase(J1, inp, { name: 'Storm' });
    J1.model.loadCases = [loadCase('Gravity', J1, inp, { betaDeg: st.beta, fG: inp.gG, fE: 0 }), loadCase('Env', J1, inp, { betaDeg: st.beta, fG: 0, fE: inp.gE, env: st1.env, sw: false, W: 0 })];
    const po = run(() => pushoverAnalysis(J1.model, { pattern: 'Env', base: 'Gravity', pDelta: true, postBucklingRatio: inp.postBuck, maxLoadFactor: inp.lfTarget, maxStep: inp.lfTarget / 3 }), 'Elastic-plastic analysis');
    const first = Number.isFinite(po.firstEventLoadFactor) ? po.firstEventLoadFactor : inp.lfTarget, nEv = po.events.length, collapsed = po.status.startsWith('collapse');
    // second-order elastic curve: a few load levels of the environmental set
    const lfs = [...new Set([1, Math.min(inp.lfTarget, 0.8 * aCr)])].filter((v) => v >= 1), el2 = lfs.map((lf) => {
      if (lf === 1) return dN.top;
      const lc = loadCase(`L${lf}`, J, inp, { betaDeg: st.beta, fG: inp.gG, fE: inp.gE * lf, env: st.env });
      J.model.loadCases = [lc];
      try { return displacements(J, pDeltaAnalysis(J.model, { loadCase: lc.name }).result).top; } catch { return null; }
    });
    const k1 = dL.top - displacements(J, solve(J, [G0]).Gravity).top;   // first-order deck displacement per unit load factor
    const steel = steelMass(J, inp), gg = groups.reduce((a, g) => (g.uc > a.uc ? g : a)), gSteel = groupSteel(J, gg.group), amp = pd.amplification, ok = ucN <= 1 && first >= 1 && !collapsed;
    const checks = [
      ...groupChecks(groups, inp, 'second order'),
      fsCheck('No plasticity or buckling at the factored design load', first, 1, 'load factor', 'ISO 19902:2020 §12.3; elastic design basis', 'Load factor on the factored storm load at which the first plastic hinge, tension yield or brace buckling occurs. Below 1.0 the structure leaves the elastic range within the design load and a linear analysis is no longer valid.'),
      fsCheck('Elastic critical load factor', aCr, 1.0, 'load factor', 'EN 1993-1-1:2005 §5.2.1 (second-order effects significant below αcr = 10)', 'Multiple of the factored load set at which the perfect elastic frame buckles. It must be well above 1; below 10 second-order effects should be included, as they are here.'),
    ];
    return {
      kpis: [
        { key: 'amp', label: 'Second-order amplification of displacement', unit: '-', value: amp, explain: 'Largest second-order displacement over the first-order value. Up to about 1.1 a linear analysis with code moment amplifiers is adequate; above 1.3 the frame is sensitive to sway and needs a nonlinear design basis.' },
        { key: 'alphaCr', label: 'Elastic critical load factor αcr', unit: '-', value: aCr, explain: 'Factor on the whole factored load set at elastic buckling of the frame. The approximate sway amplifier is 1/(1 − 1/αcr).' },
        { key: 'ucLin', label: 'Highest utilisation — first order', unit: '-', value: ucL, limit: 1, sense: 'max', explain: 'Member utilisation from the linear analysis, for comparison.' },
        { key: 'ucNL', label: 'Highest utilisation — second order', unit: '-', value: ucN, limit: 1, sense: 'max', explain: 'Member utilisation with second-order forces. The difference from the first-order value is the cost of frame slenderness.' },
        { key: 'dispLin', label: 'Deck displacement — first order', unit: 'mm', value: dL.top * 1e3, explain: 'Lateral displacement of the leg tops from the linear analysis under the factored storm load.' },
        { key: 'dispNL', label: 'Deck displacement — second order', unit: 'mm', value: dN.top * 1e3, explain: 'Lateral displacement of the leg tops including the geometric stiffness of the member axial forces.' },
        { key: 'lfFirst', label: 'Load factor at first plastic event', unit: '-', value: first, limit: 1, sense: 'min', explain: `Multiple of the factored environmental load at the first hinge, yield or buckling event (capped at the target ${fmt(inp.lfTarget)} when none occurs). The margin above 1.0 is the elastic reserve.` },
        { key: 'nEvents', label: 'Plastic events up to the target load', unit: '-', value: nEv, explain: 'Number of hinges formed and members yielded or buckled before the target load factor. Zero means the structure is still elastic at the target.' },
        { key: 'lfReached', label: 'Load factor reached', unit: '-', value: po.ultimateLoadFactor, explain: 'Highest load factor carried in the elastic-plastic trace: the target, or the collapse load if a mechanism forms before it.' },
        { key: 'iterations', label: 'Second-order iterations', unit: '-', value: pd.iterations, explain: 'Iterations on the member axial forces until the geometric stiffness stops changing; more than about 10 indicates a load close to the critical load.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Load factor against deck displacement', x: { label: 'Control displacement', unit: 'mm' }, y: { label: 'Load factor on factored environmental load', unit: '-' },
          series: [{ name: 'Elastic-plastic with P–Δ', x: po.curve.displacement.map((v) => round(Math.abs(v) * 1e3)), y: po.curve.loadFactor.map((v) => round(v)) },
            { name: 'Second-order elastic (deck)', x: el2.filter((v) => v != null).map((v) => round(v * 1e3)), y: lfs.filter((v, i) => el2[i] != null) },
            { name: 'First-order elastic (deck)', x: [round((dL.top - k1) * 1e3), round((dL.top - k1 + k1 * inp.lfTarget) * 1e3)], y: [0, inp.lfTarget], dash: true }],
          explain: 'Response as the environmental load is scaled with gravity held. The curves separate where geometric (bending away from the straight line) and material (kinks at plastic events) nonlinearity start to matter; a horizontal branch is collapse.' },
        plot3d(J, { title: 'Second-order utilisation and deformed shape', values: nodeValuesFromMembers(J, rowsN.map((r) => r.uc)), label: 'Utilisation', unit: '-', u: pd.result.displacements, explain: 'Member utilisation with second-order forces on the second-order deflected shape at the factored design load.' }),
        { type: 'bar', title: 'First-order and second-order utilisation by group', categories: groups.map((g) => g.label), series: [{ name: 'First order', y: byGroup(rowsL).map((g) => round(g.uc)) }, { name: 'Second order', y: groups.map((g) => round(g.uc)) }], y: { label: 'Utilisation', unit: '-' },
          explain: 'Effect of second-order forces on each member group. Legs normally gain most because they carry the gravity load that drives P–Δ.' },
      ],
      tables: [
        { title: 'Plastic events up to the target load factor', columns: [{ key: 'k', label: 'No.' }, { key: 'lf', label: 'Load factor' }, { key: 'm', label: 'Member' }, { key: 'g', label: 'Group' }, { key: 'type', label: 'Event' }, { key: 'd', label: 'Control displacement', unit: 'mm' }],
          rows: po.events.slice(0, 40).map((e, k) => { const o = J1.memberMap.findIndex((l) => l.includes(e.member)); return { k: k + 1, lf: round(e.loadFactor), m: o, g: J.orig[o]?.label ?? '', type: e.type, d: round(Math.abs(e.displacement) * 1e3) }; }),
          explain: 'Sequence of hinge formation, tension yield and brace buckling as the load grows. An empty table means the frame stays elastic up to the target load factor.' },
        { title: 'Most utilised members (second order)', columns: memberCols, rows: memberTableRows(rowsN, 12), explain: 'Members with the highest code utilisation using second-order section forces.' },
      ],
      notes: [`Second-order effects raise displacements by ${fmt((amp - 1) * 100)} % and move the highest utilisation from ${fmt(ucL)} to ${fmt(ucN)}; αcr = ${fmt(aCr)} (1/(1 − 1/αcr) = ${fmt(1 / (1 - 1 / aCr))} bounds the amplification from above when the first buckling mode is a single brace rather than frame sway).`,
        nEv ? `First plastic event at ${fmt(first)} × the factored environmental load (${po.events[0].type}); ${nEv} event(s) up to ${fmt(po.ultimateLoadFactor)}${collapsed ? ', where a collapse mechanism forms' : ''}.` : `No plastic event up to ${fmt(inp.lfTarget)} × the factored environmental load: the structure is elastic over the whole range.`,
        `Wave by ${st.wf.label}, heading ${fmt(st.beta)}°.`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `The frame is stable and elastic at the design load: amplification ${fmt(amp)}, αcr ${fmt(aCr)}, second-order utilisation ${fmt(ucN)}, first plastic event at load factor ${fmt(first)}.`, amp < 1.1 ? 'Second-order effects are small, so the linear in-place analysis with the code beam-column amplifiers is a sound design basis. Run the pushover module for the collapse load and reserve strength ratio.' : 'Second-order effects are significant: design the legs with second-order forces or stiffen the foundation and bracing, and run the pushover module for the collapse load.', 'ISO 19902:2020 §12.3; EN 1993-1-1:2005 §5.2', `Elastic reserve of ${fmt((first - 1) * 100)} % before any member yields means no permanent set in the design storm and no post-storm repair steel. The jacket holds ${steelCo2(steel)}`)
          : adv('engineering', 'fail', `Nonlinear check failed: second-order utilisation ${fmt(ucN)}, first plastic event at load factor ${fmt(first)}${collapsed ? `, collapse at ${fmt(po.ultimateLoadFactor)}` : ''}.`, `Strengthen the ${gg.label.toLowerCase()} (about ${fmt(Math.max(gg.uc / 0.9 - 1, 0.1) * 100)} % more wall) or stiffen the foundation; then repeat with the second-order forces.`, BASIS_MEMBER, `Strengthening is roughly ${steelCo2(gSteel * Math.max(gg.uc / 0.9 - 1, 0.1))}`),
        adv('economic', 'info', `Second-order forces change the governing member utilisation by ${fmt((ucN / ucL - 1) * 100)} % and displacements by ${fmt((amp - 1) * 100)} % on a ${fmt(steel)} t jacket.`, 'Where the amplification is below about 1.1, do not pay for nonlinear analysis in every design cycle: use linear analysis for sizing and reserve the nonlinear run for the final check and for reassessment, where proving reserve strength avoids offshore strengthening.', `${BASIS_COST}; ISO 19902:2020 §12.3`),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Cantilever beam-column, axial P at the tip with lateral H: tip deflection = (HL³/3EI)·3(tan u − u)/u³, u = L√(P/EI).
    // L = 95 m, EI = 2.1e11·0.054485 = 1.14418e10 N·m², P = 1.0 MN → u = 95·√(1e6/1.14418e10) = 0.88813,
    // tan u = 1.22989, amplification = 3(1.22989 − 0.88813)/0.88813³ = 1.4636.
    { name: 'Cantilever beam-column: amplification 3(tan u − u)/u³', source: 'Timoshenko & Gere, Theory of Elastic Stability (1961) §1.11',
      inputs: { geom: testFrame('column'), H: 0, Uc: 0, Uw: 0, windA: 0, swCont: 0, Wtop: 1, FxTop: 10, found: 'fixed', gG: 1, gE: 1, nSub: 8, lfTarget: 1 },
      expect: [{ kpi: 'amp', value: 1.4636, tol: 0.005 }, { kpi: 'alphaCr', value: 3.128, tol: 0.005 }] },
  ],
};
