// Progressive collapse and redundancy: members are removed one at a time (severed brace, failed joint, missing leg
// section), the damaged frame is re-analysed with the wave load on the remaining members, and the worst case is
// pushed to collapse to give the damaged strength ratio.
import { need, run, fmt, round, adv, fsCheck, check, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, PUSHOVER_CAVEAT, CONFIRM, BASIS_COST } from './_common.js';
import { stdInputs, buildJacket, stormCase, envLoader, loadCase, solve, memberRows, displacements, pushCases, pushCurve, nodeValuesFromMembers, plot3d, steelMass, foundNote, testFrame } from './_model.js';
import { pushoverAnalysis } from '../../engine/structural/pushover.js';

const GROUPS = [{ value: 'brace', label: 'Diagonal braces' }, { value: 'horizontal', label: 'Horizontal braces' }, { value: 'leg', label: 'Legs' }, { value: 'all', label: 'All members' }];

export default {
  id: 's01-progressive-collapse', suite: 1, title: 'Progressive collapse and redundancy',
  summary: 'Member-removal study: utilisation and displacement of the damaged jacket for each removed member, and the residual collapse strength and damaged strength ratio for the most critical removal.',
  covers: [15], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §10 (accidental situations), §12.6, §A.9.9 (damaged condition)', 'API RP 2SIM (2014) §8.5', 'NORSOK N-004 (2022) §A.3 (accidental limit state, post-damage check)', 'DNV-RP-C204 (2019) §2'],
  models: ['Alternate load path method (notional member removal)', 'Direct stiffness re-analysis of the damaged frame', 'Nonlinear pushover of intact and damaged structure', 'Damaged strength ratio and redundancy factor', 'Morison loading on the remaining members'],
  limits: `Static redistribution: the dynamic overshoot of a sudden failure is represented only by the amplification factor on the loads. Removal is complete (no residual strength of the damaged member) and joints and piles are not failure elements. ${PUSHOVER_CAVEAT} ${MORISON_CAVEAT} ${CODE_CAVEAT}`,
  searchTerms: 'jacket redundancy member removal damaged strength ratio progressive collapse',
  inputs: [
    ...stdInputs({ crit: { gG: 1.0, gE: 1.0 }, mesh: { nSub: 1, nPhase: 16 } }),
    { key: 'removed', label: 'Members to remove together (optional)', type: 'text', default: '', group: 'Damage scenario', role: 'ic', help: 'Comma-separated member numbers of one damage scenario, as listed in the member tables. Leave empty to screen single-member removals automatically.' },
    { key: 'group', label: 'Members screened', type: 'select', default: 'brace', group: 'Damage scenario', role: 'ic', options: GROUPS, help: 'Each of the most utilised members of this group is removed in turn.' },
    { key: 'nCand', label: 'Number of single removals screened', type: 'int', default: 6, min: 1, max: 40, group: 'Damage scenario', role: 'control', help: 'The most utilised members of the group in the intact structure are tried first.' },
    { key: 'daf', label: 'Dynamic amplification for sudden removal', unit: '-', type: 'number', default: 1.0, min: 1, max: 2, group: 'Damage scenario', role: 'criteria', help: '1.0 for pre-existing damage found by inspection; up to 2.0 for sudden loss under load (upper bound of a static analysis).' },
    { key: 'postBuck', label: 'Post-buckling residual force ratio', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Design criteria', role: 'criteria' },
    { key: 'rsrDam', label: 'Required load factor of the damaged structure', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 3, group: 'Design criteria', role: 'criteria', help: `The damaged structure must carry the specified environmental load set (normally a reduced return period with action factors of 1.0) with a collapse factor of at least this value. ${CONFIRM}` },
  ],
  mesh: { key: 'nPhase', label: 'Wave phase steps', levels: [8, 16, 32], qoi: 'baseShear' },
  compute(inp, ctx) {
    // SI conversion inside the shared model. Loads carry the given action factors times the removal amplification.
    need(inp.nSub === 1 || inp.nSub == null, 'Member-removal analysis uses one element per member: set elements per member to 1.');
    const fi = { ...inp, gG: inp.gG * inp.daf, gE: inp.gE * inp.daf };
    const J = buildJacket(inp, ctx), st = stormCase(J, fi, { name: 'Intact' }), tF = st.sc.t[st.sc.iF];
    need(st.baseShear > 1 || st.otm > 1, 'The environmental load pattern is zero: give a wave, current, wind or interface load for the damaged-condition check.');
    const res0 = solve(J, [st.lc]).Intact, rows0 = memberRows(J, res0, inp, { name: 'intact' }), uc0 = Math.max(...rows0.map((r) => r.uc)), d0 = displacements(J, res0).top;
    const damaged = (rm) => {
      // rebuild without the removed members and reload: they no longer attract wave load
      try {
        const Jd = buildJacket(inp, ctx, { remove: rm }), env = envLoader(Jd, fi, st.wf, st.beta).at(tF, true);
        const r = solve(Jd, [loadCase('Damaged', Jd, fi, { betaDeg: st.beta, fG: fi.gG, fE: fi.gE, env })]).Damaged, rw = memberRows(Jd, r, inp, { name: 'damaged' }), g = rw.reduce((a, q) => (q.uc > a.uc ? q : a));
        // map the row index of the damaged model back to the intact numbering
        const keep = J.orig.map((m, i) => i).filter((i) => !rm.includes(i));
        return { Jd, r, rw, uc: g.uc, gov: keep[g.i], ref: g.ref, disp: displacements(Jd, r).top, keep, ok: true };
      } catch (e) { return { uc: 99, gov: -1, ref: 'mechanism', disp: Infinity, ok: false, msg: e.message }; }
    };
    let scen;
    const explicit = String(inp.removed || '').split(/[,;\s]+/).filter(Boolean).map(Number);
    if (explicit.length) {
      need(explicit.every((i) => Number.isInteger(i) && i >= 0 && i < J.orig.length), `Member numbers to remove must be whole numbers between 0 and ${J.orig.length - 1}.`);
      scen = [{ rm: explicit, name: `Members ${explicit.join(', ')}` }];
    } else {
      const pool = rows0.filter((r) => inp.group === 'all' || r.group === inp.group);
      need(pool.length, `The model has no members in the group "${inp.group}": choose another group or list the members to remove.`);
      scen = pool.sort((a, b) => b.uc - a.uc).slice(0, inp.nCand).map((r) => ({ rm: [r.i], name: `Member ${r.i} (${r.label})` }));
    }
    scen.forEach((s) => Object.assign(s, damaged(s.rm)));
    // most critical removal: highest utilisation of the remaining members, weighted by the loss of lateral stiffness
    const score = (s) => (s.ok ? s.uc * (d0 > 0 ? s.disp / d0 : 1) : Infinity), worst = scen.reduce((a, s) => (score(s) > score(a) ? s : a));
    ctx.progress?.(0.5, 'Pushover of the intact and damaged structure');
    // collapse strength, intact and for the worst removal, for the same environmental pattern
    // the geometric stiffness is updated at every member failure (no sub-stepping between events) to keep two pushovers affordable
    const opts = { pattern: 'Env', base: 'Gravity', pDelta: true, postBucklingRatio: inp.postBuck, maxStep: 1e9 };
    pushCases(J, fi, st);
    const pi = run(() => pushoverAnalysis(J.model, opts), 'Pushover of the intact structure');
    need(pi.status !== 'collapse-under-base-load' && pi.ultimateLoadFactor > 0 && Number.isFinite(pi.ultimateLoadFactor), 'The intact structure has no collapse load to compare with: it fails under gravity alone or the environmental load pattern is zero. Check member sizes, topsides weight and the wave, current or interface loads.');
    let pdm = { ultimateLoadFactor: 0, events: [], curve: { loadFactor: [0], displacement: [0], baseShear: [0] }, status: 'collapse' };
    if (worst.ok) {
      const env = envLoader(worst.Jd, fi, st.wf, st.beta).at(tF, true);
      pushCases(worst.Jd, fi, { beta: st.beta, env });
      try { pdm = pushoverAnalysis(worst.Jd.model, opts); } catch (e) { if (!e.mechanism) throw new Error(`Pushover of the damaged structure: ${e.message}`); }
    }
    const lfI = pi.ultimateLoadFactor, lfD = Number.isFinite(pdm.ultimateLoadFactor) ? pdm.ultimateLoadFactor : 0, dsr = lfI > 0 ? lfD / lfI : 0, Vref = Math.abs(pi.baseShearRef);
    const first = Number.isFinite(pi.firstEventLoadFactor) ? pi.firstEventLoadFactor : lfI, nSurv = scen.filter((s) => s.uc <= 1).length, ci = pushCurve(pi), cd = pushCurve(pdm);
    const steel = steelMass(J, inp), ok = lfD >= inp.rsrDam, wm = J.orig[worst.rm[0]];
    const vals = worst.ok ? (() => { const full = new Array(J.orig.length).fill(0); worst.rw.forEach((r, k) => { full[worst.keep[k]] = r.uc; }); return full; })() : rows0.map((r) => r.uc);
    const nv = nodeValuesFromMembers(J, vals);
    return {
      kpis: [
        { key: 'dsr', label: 'Damaged strength ratio (worst removal)', unit: '-', value: dsr, explain: 'Collapse load of the structure without the removed member(s) over the intact collapse load. Above about 0.8 the frame is robust to this damage; below 0.5 the member is a key element.' },
        { key: 'rsrIntact', label: 'Collapse load factor — intact', unit: '-', value: lfI, explain: 'Multiple of the specified environmental load set that the intact structure carries (the RSR when the load set is the unfactored 100-year condition).' },
        { key: 'rsrDamaged', label: 'Collapse load factor — damaged', unit: '-', value: lfD, limit: inp.rsrDam, sense: 'min', explain: 'Multiple of the specified load set carried after the worst removal. It must stay above the required factor for the platform to survive until repair.' },
        { key: 'redundancy', label: 'Redundancy factor (intact ultimate / first failure)', unit: '-', value: lfI / first, explain: 'Reserve between the first member failure and collapse of the intact frame. Values near 1.0 mean little alternative load path.' },
        { key: 'ucIntact', label: 'Highest utilisation — intact', unit: '-', value: uc0, explain: 'Linear code utilisation of the intact structure for the same load set, for reference.' },
        { key: 'ucDamaged', label: 'Highest utilisation after the worst removal', unit: '-', value: Math.min(worst.uc, 99), explain: 'Linear code utilisation of the remaining members. Above 1.0 further members fail and the pushover decides whether the collapse progresses.' },
        { key: 'nSurvive', label: 'Removals with all members within code limits', unit: '-', value: nSurv, explain: `Number of the ${scen.length} removal case(s) in which no remaining member exceeds a utilisation of 1.0.` },
        { key: 'baseShear', label: 'Reference base shear', unit: 'MN', value: Vref / 1e6, explain: 'Environmental base shear of the load set applied to intact and damaged structures (with action and amplification factors).' },
        { key: 'dispRatio', label: 'Deck displacement, damaged / intact', unit: '-', value: Number.isFinite(worst.disp) && d0 > 0 ? worst.disp / d0 : 99, explain: 'Loss of lateral stiffness caused by the worst removal; a ratio well above 1 shows the removed member was a main load path.' },
      ],
      checks: [
        fsCheck('Damaged structure carries the specified load', lfD, inp.rsrDam, 'load factor', 'ISO 19902:2020 §10.1.3, §A.9.9; NORSOK N-004 §A.3 (post-damage resistance)', 'Collapse factor of the structure without the removed member(s). Below the requirement a single member failure can grow into global collapse in the specified storm.'),
        ...scen.slice(0, 12).map((s) => check(`Remaining members after removing ${s.name.toLowerCase()}`, Math.min(s.uc, 99), 1, 'utilisation', `${inp.code} member check on the damaged frame (${s.ref})`, 'Highest code utilisation of the remaining members after this removal. Above 1.0 a second member is overloaded and the failure may progress; the pushover result decides.')),
      ],
      plots: [
        { type: 'line', title: 'Pushover: intact and after the worst removal', x: { label: 'Control displacement', unit: 'mm' }, y: { label: 'Load factor', unit: '-' }, series: [{ name: 'Intact', x: ci.d, y: ci.lf }, { name: `Without ${worst.name.toLowerCase()}`, x: cd.d, y: cd.lf }, { name: 'Required', x: [0, Math.max(...ci.d)], y: [inp.rsrDam, inp.rsrDam], dash: true }],
          explain: 'Load factor against displacement for the intact structure and for the worst damage case. The drop of the peak is the strength lost; the damaged curve must still clear the dashed required level.' },
        { type: 'bar', title: 'Highest utilisation after each removal', categories: scen.map((s) => (s.rm.length === 1 ? `${s.rm[0]}` : s.rm.join('+'))), series: [{ name: 'Damaged', y: scen.map((s) => round(Math.min(s.uc, 5))) }, { name: 'Intact', y: scen.map(() => round(uc0)) }], y: { label: 'Utilisation', unit: '-' },
          explain: 'Governing utilisation of the remaining members for each removed member (bars capped at 5). Members whose removal pushes the bar far above the intact level are the critical load paths to inspect first.' },
        plot3d(J, { title: `Utilisation after removing ${worst.name.toLowerCase()}`, values: nv, label: 'Utilisation', unit: '-', u: worst.ok && worst.Jd.model.nodes.length === J.model.nodes.length ? worst.r.displacements : res0.displacements, explain: 'Utilisation of the remaining members after the worst removal, on the displaced shape. The removed member shows zero; neighbours that pick up its load show the alternative load path.' }),
      ],
      tables: [
        { title: 'Member-removal cases', columns: [{ key: 'n', label: 'Removed' }, { key: 'uc', label: 'Highest utilisation' }, { key: 'gov', label: 'Governing member' }, { key: 'ref', label: 'Check' }, { key: 'd', label: 'Deck displacement', unit: 'mm' }, { key: 'r', label: 'Displacement ratio' }, { key: 's', label: 'Within limits' }],
          rows: scen.map((s) => ({ n: s.name, uc: round(Math.min(s.uc, 99)), gov: s.gov, ref: s.ref, d: Number.isFinite(s.disp) ? round(s.disp * 1e3) : '–', r: Number.isFinite(s.disp) && d0 > 0 ? round(s.disp / d0) : '–', s: s.uc <= 1 ? 'yes' : 'no' })), explain: 'One row per damage case: which remaining member governs, its utilisation and the change in deck displacement. Use it to rank members by consequence of failure for the inspection plan.' },
        { title: 'Failure sequence of the damaged structure', columns: [{ key: 'k', label: 'No.' }, { key: 'lf', label: 'Load factor' }, { key: 'm', label: 'Member (intact numbering)' }, { key: 'type', label: 'Event' }],
          rows: pdm.events.slice(0, 30).map((e, k) => ({ k: k + 1, lf: round(e.loadFactor), m: worst.keep ? worst.keep[worst.Jd.elems[e.member]?.o ?? e.member] : e.member, type: e.type })), explain: 'Order of member failures after the worst removal as the load grows: how the collapse would progress.' },
      ],
      notes: [`Worst case: ${worst.name} removed — utilisation ${fmt(Math.min(worst.uc, 99))} (intact ${fmt(uc0)}), collapse factor ${fmt(lfD)} against ${fmt(lfI)} intact, damaged strength ratio ${fmt(dsr)}.`,
        `${nSurv} of ${scen.length} removal case(s) leave every remaining member within code limits${inp.daf > 1 ? `; loads include a sudden-removal amplification of ${fmt(inp.daf)}` : ''}.`,
        `Heading ${fmt(st.beta)}°, wave by ${st.wf.label}; removed members attract no wave load.`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `The frame tolerates the loss of ${worst.name.toLowerCase()}: collapse factor ${fmt(lfD)} ≥ ${fmt(inp.rsrDam)}, damaged strength ratio ${fmt(dsr)}.`, `Rank inspection by consequence: the members with the largest utilisation jump in the table are the key elements. Repeat for other headings and for leg damage, and check the joints of the members that pick up load (${worst.gov >= 0 ? `member ${worst.gov}` : 'see table'}).`, 'ISO 19902:2020 §10 and §A.9.9; API RP 2SIM (2014) §8.5', `Redundancy is what allows a damaged member to wait for a planned repair campaign instead of an emergency mobilisation. ${wm ? `A replacement of the critical member is about ${steelCo2(wm.sec.A * wm.L * J.steel.rho / 1e3)}` : ''}`)
          : adv('engineering', 'fail', `Loss of ${worst.name.toLowerCase()} leads to collapse below the specified load: factor ${fmt(lfD)} against ${fmt(inp.rsrDam)} (damaged strength ratio ${fmt(dsr)}).`, 'Treat the member as a key element: protect it from impact, raise its inspection frequency, or add an alternative load path (X instead of K or single-diagonal bracing, heavier horizontals).', 'ISO 19902:2020 §10; NORSOK N-004 §A.3', `Adding an alternative load path in one bay is of the order of ${steelCo2(steel * 0.03)}`),
        adv('economic', 'info', `Damaged strength ratio ${fmt(dsr)}; redundancy factor of the intact frame ${fmt(lfI / first)}.`, 'Use the removal ranking to set inspection scope: flooded-member detection on the high-consequence braces every campaign, the low-consequence ones on a longer interval. Each avoided diver or ROV day is a direct saving, and a robust frame defers repair to a planned vessel visit.', `${BASIS_COST}; risk-based inspection per API RP 2SIM / DNV-RP-G101`),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Fixed-base portal (two columns 95 m, beam 20 m, same tube): sway mechanism needs four hinges, λ·H·L = 4Mp;
    // with one column removed the other is a cantilever, λ·H·L = Mp. Damaged strength ratio = 1/4.
    // Mp = 355e6·0.095296 = 33.83 MN·m, H = 200 kN: intact λ ≈ 4·33.83e6/(2e5·95) = 7.12 (slightly less through
    // axial force in the columns and P–Δ, which ends the analysis just before the fourth hinge), damaged λ = 1.78.
    { name: 'Portal frame: column removal reduces the sway collapse load from 4Mp/L to Mp/L', source: 'Plastic collapse of portal frames (Neal, The Plastic Methods of Structural Analysis, 1977, §3)',
      inputs: { geom: testFrame('portal'), H: 0, Uc: 0, Uw: 0, windA: 0, swCont: 0, Wtop: 0, FxTop: 200, found: 'fixed', heading: 0, removed: '1' },
      expect: [{ kpi: 'rsrDamaged', value: 1.7805, tol: 0.01 }, { kpi: 'dsr', value: 0.25, tol: 0.05 }] },
  ],
};
