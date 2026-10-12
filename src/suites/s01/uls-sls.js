// Ultimate and serviceability limit states: the permanent, variable and environmental actions are solved once and
// combined with the partial action factors of a load-combination table; ULS combinations are checked for member
// strength and pile reactions, SLS combinations for deck sway, deck tilt and brace deflection.
import { need, fmt, round, adv, check, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, CONFIRM, BASIS_MEMBER, BASIS_COST } from './_common.js';
import { geomInputs, topsidesInputs, envInputs, foundationInputs, criteriaInputs, meshInputs, buildJacket, designWave, loadCase, solve, memberRows, envelope, byGroup, pileForces, displacements, nodeValuesFromMembers, plot3d, memberCols, memberTableRows, groupChecks, steelMass, groupSteel, foundNote, testFrame } from './_model.js';
import { deflectionCheck } from '../../engine/structural/steelcheck.js';

const crit = criteriaInputs().filter((i) => i.key !== 'gG' && i.key !== 'gE');

export default {
  id: 's01-uls-sls', suite: 1, title: 'Ultimate and serviceability limit states',
  summary: 'Applies the partial-factor load combinations (gravity-dominated, operating storm, extreme storm, minimum gravity) for the ultimate limit state and unfactored combinations for serviceability: member utilisation, pile reactions, deck sway, tilt and brace deflection.',
  covers: [5], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §7 (limit states), §9.10 (partial action factors), §13', 'ISO 19900:2019 §8', 'DNV-OS-C101 (2021) Ch.2 Sec.1–2 (ULS a/b, SLS)', 'NORSOK N-001 (2020) §7.2'],
  models: ['Limit state design with partial action and resistance factors', 'Linear superposition of permanent, variable and environmental load cases', 'Direct stiffness method', 'Morison equation', 'ISO 19902 tubular member checks', 'Deflection limits (span ratio)'],
  limits: `Linear elastic combinations of one extreme and one operating regular wave on a single heading; accidental and fatigue limit states are in their own modules. Serviceability limits (sway, tilt, deflection) are project criteria, not code values. ${MORISON_CAVEAT} ${CODE_CAVEAT}`,
  searchTerms: 'limit state design jacket partial factors ULS SLS load combinations',
  inputs: [
    ...geomInputs(), ...topsidesInputs(),
    { key: 'permFrac', label: 'Permanent share of the topsides weight', unit: '-', type: 'number', default: 0.7, min: 0, max: 1, group: 'Topsides and gravity loads', role: 'bc', help: 'Structure, equipment and piping are permanent (G); stored liquids, consumables, live and drilling loads are variable (Q).' },
    ...envInputs({ group: 'Environment — extreme (100-year)' }).slice(0, 5),
    { key: 'Hop', label: 'Operating wave height', unit: 'm', type: 'number', default: 10, min: 0, max: 40, group: 'Environment — operating (1-year)', role: 'bc', help: 'Maximum individual wave of the 1-year sea state.' },
    { key: 'Top', label: 'Operating wave period', unit: 's', type: 'number', default: 10.5, min: 3, max: 25, group: 'Environment — operating (1-year)', role: 'bc', bind: 'site.Tp1' },
    { key: 'UcOp', label: 'Operating surface current', unit: 'm/s', type: 'number', default: 0.6, min: 0, max: 4, group: 'Environment — operating (1-year)', role: 'bc', bind: 'site.Uc1' },
    { key: 'UwOp', label: 'Operating 1-hour wind at 10 m', unit: 'm/s', type: 'number', default: 22, min: 0, max: 80, group: 'Environment — operating (1-year)', role: 'bc', bind: 'site.Uw1' },
    ...envInputs().slice(5),
    { key: 'combos', label: 'Load combinations', type: 'table', group: 'Limit states', role: 'criteria',
      columns: [{ key: 'name', label: 'Combination', type: 'text' }, { key: 'ls', label: 'Limit state (ULS / SLS)', type: 'text' }, { key: 'fG', label: 'γ permanent', type: 'number' }, { key: 'fQ', label: 'γ variable', type: 'number' }, { key: 'fE', label: 'γ environmental', type: 'number' }, { key: 'env', label: 'Environment (extreme / operating / none)', type: 'text' }],
      default: [{ name: 'Gravity dominated', ls: 'ULS', fG: 1.3, fQ: 1.5, fE: 0, env: 'none' }, { name: 'Operating storm', ls: 'ULS', fG: 1.3, fQ: 1.5, fE: 1.215, env: 'operating' }, { name: 'Extreme storm', ls: 'ULS', fG: 1.1, fQ: 1.1, fE: 1.35, env: 'extreme' }, { name: 'Extreme storm, minimum gravity', ls: 'ULS', fG: 0.9, fQ: 0.8, fE: 1.35, env: 'extreme' }, { name: 'Serviceability, operating', ls: 'SLS', fG: 1, fQ: 1, fE: 1, env: 'operating' }, { name: 'Serviceability, extreme', ls: 'SLS', fG: 1, fQ: 1, fE: 1, env: 'extreme' }],
      help: `Partial action factors as recalled from ISO 19902 §9.10 (1.3 G + 1.5 Q; 1.3 G + 1.5 Q + 0.9·1.35 E operating; 1.1 G + 1.1 Q + 1.35 E extreme; 0.9 G + 0.8 Q + 1.35 E when gravity relieves). ${CONFIRM}` },
    ...foundationInputs(), ...crit,
    { key: 'swayRatio', label: 'Deck sway limit: height above mudline / ratio', unit: '-', type: 'number', default: 250, min: 50, max: 1000, group: 'Serviceability criteria', role: 'criteria', help: 'Project criterion for lateral displacement at the jacket top in serviceability combinations (riser, conductor and bridge tolerances).' },
    { key: 'tiltLim', label: 'Deck tilt limit', unit: 'deg', type: 'number', default: 0.5, min: 0.05, max: 3, group: 'Serviceability criteria', role: 'criteria', help: 'Project criterion for process equipment and drilling (often 0.25°–0.5°).' },
    { key: 'deflRatio', label: 'Brace deflection limit: span / ratio', unit: '-', type: 'number', default: 250, min: 100, max: 1000, group: 'Serviceability criteria', role: 'criteria' },
    ...meshInputs(),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'deckSway' },
  compute(inp, ctx) {
    // SI conversion inside the shared model. Base cases are solved unfactored; combinations are factored sums.
    const combos = (inp.combos || []).map((c, i) => ({ name: String(c.name || `Combination ${i + 1}`), ls: String(c.ls || 'ULS').trim().toUpperCase(), fG: Number(c.fG), fQ: Number(c.fQ), fE: Number(c.fE), env: String(c.env || 'none').trim().toLowerCase() }));
    need(combos.length, 'Give at least one load combination.');
    need(combos.every((c) => ['ULS', 'SLS'].includes(c.ls) && ['extreme', 'operating', 'none'].includes(c.env) && [c.fG, c.fQ, c.fE].every((v) => Number.isFinite(v) && v >= 0)), 'Each combination needs a limit state (ULS or SLS), non-negative factors and an environment (extreme, operating or none).');
    need(new Set(combos.map((c) => c.name)).size === combos.length, 'Combination names must be unique.');
    need(combos.some((c) => c.ls === 'ULS') && combos.some((c) => c.ls === 'SLS'), 'Give at least one ULS and one SLS combination.');
    const J = buildJacket(inp, ctx), beta = inp.heading ?? 0, Wt = inp.Wtop * 1e6;
    const iop = { ...inp, H: inp.Hop, T: inp.Top, Uc: inp.UcOp, Uw: inp.UwOp };
    const dX = designWave(J, inp, beta), dO = designWave(J, iop, beta);
    const base = [loadCase('G', J, { ...inp, FxTop: 0, MyTop: 0, windA: 0 }, { betaDeg: beta, fG: 1, fE: 0, W: Wt * inp.permFrac }),
      { name: 'Q', nodal: loadCase('Q', J, { ...inp, FxTop: 0, MyTop: 0, windA: 0 }, { betaDeg: beta, fG: 1, fE: 0, sw: false, W: Wt * (1 - inp.permFrac) }).nodal, member: [] },
      loadCase('Eext', J, inp, { betaDeg: beta, fG: 0, fE: 1, sw: false, W: 0, env: dX.ld.at(dX.sc.t[dX.sc.iF], true) }),
      loadCase('Eop', J, iop, { betaDeg: beta, fG: 0, fE: 1, sw: false, W: 0, env: dO.ld.at(dO.sc.t[dO.sc.iF], true) })];
    J.model.combinations = combos.map((c) => ({ name: c.name, factors: { G: c.fG, Q: c.fQ, ...(c.env === 'extreme' ? { Eext: c.fE } : c.env === 'operating' ? { Eop: c.fE } : {}) } }));
    const res = solve(J, base), P = J.model.nodes, hTop = J.zTop - J.zMud;
    for (const c of combos) need(res[c.name] && Number.isFinite(res[c.name].displacements[0]), `Combination "${c.name}" could not be solved.`);
    // brace deflection relative to the chord line between member ends (from the nodes along each member)
    const sag = (r) => { let best = { ratio: 0, defl: 0, i: -1, L: 1 }; J.orig.forEach((m, i) => { if (m.group === 'leg') return; const ids = [J.model.members[J.memberMap[i][0]].a, ...J.memberMap[i].map((e) => J.model.members[e].b)], u = r.displacements, a = ids[0], b = ids[ids.length - 1]; for (let k = 1; k < ids.length - 1; k++) { const w = k / (ids.length - 1), n = ids[k], d = Math.hypot(...[0, 1, 2].map((q) => u[6 * n + q] - ((1 - w) * u[6 * a + q] + w * u[6 * b + q]))); if (d / m.L > best.ratio) best = { ratio: d / m.L, defl: d, i, L: m.L }; } }); return best; };
    const out = combos.map((c) => {
      const r = res[c.name], rows = memberRows(J, r, inp, { name: c.name, hydro: true }), g = rows.reduce((a, q) => (q.uc > a.uc ? q : a)), pf = pileForces(J, r), dp = displacements(J, r), u = r.displacements;
      // tilt of the plane through the leg tops: largest differential vertical displacement over the horizontal distance
      let tilt = 0; for (const a of J.topNodes) for (const b of J.topNodes) { const d = Math.hypot(P[a][0] - P[b][0], P[a][1] - P[b][1]); if (a < b && d > 0.1) tilt = Math.max(tilt, Math.abs(u[6 * a + 2] - u[6 * b + 2]) / d); }
      return { ...c, r, rows, g, pf, dp, tilt: tilt || dp.rotMax * 0, sag: sag(r) };
    });
    const uls = out.filter((c) => c.ls === 'ULS'), sls = out.filter((c) => c.ls === 'SLS'), envU = envelope(uls.map((c) => c.rows)), groups = byGroup(envU);
    const gU = uls.reduce((a, c) => (c.g.uc > a.g.uc ? c : a)), gS = sls.reduce((a, c) => (c.dp.top > a.dp.top ? c : a)), gT = sls.reduce((a, c) => (c.tilt > a.tilt ? c : a)), gD = sls.reduce((a, c) => (c.sag.ratio > a.sag.ratio ? c : a));
    const swayLim = hTop / inp.swayRatio, tiltDeg = gT.tilt * 180 / Math.PI, pileC = Math.max(...uls.map((c) => c.pf.axialMax)), pileT = Math.max(...uls.map((c) => c.pf.tensionMax)), pileS = Math.max(...uls.map((c) => c.pf.shearMax));
    const dc = gD.sag.i >= 0 ? deflectionCheck({ deflection: gD.sag.defl, span: gD.sag.L, ratio: inp.deflRatio }) : { uc: 0, limit: 1, actualRatio: Infinity };
    const checks = [
      ...groupChecks(groups, inp, 'ULS envelope'),
      check(`SLS deck sway (${gS.name})`, gS.dp.top * 1e3, swayLim * 1e3, 'mm', `Project criterion: height / ${inp.swayRatio}; ISO 19902:2020 §7.4 (serviceability)`, 'Lateral displacement of the jacket top under unfactored loads. Excess sway overstresses conductors, risers and bridge bearings and is felt by personnel.'),
      check(`SLS deck tilt (${gT.name})`, tiltDeg, inp.tiltLim, 'deg', 'Project criterion for process and drilling equipment', 'Rotation of the plane through the leg tops from differential leg shortening and pile-head settlement. Tilt disturbs drainage, rotating equipment and drilling operations.'),
      check(`SLS brace deflection (${gD.name})`, gD.sag.defl * 1e3, dc.limit * 1e3, 'mm', `Span / ${inp.deflRatio} (project criterion; EN 1990 A1.4 practice)`, 'Deflection of the most flexible brace relative to its ends under unfactored loads. Large values indicate a slender member prone to vibration and to damage of attachments and anodes.'),
    ];
    const steel = steelMass(J, inp), gg = groups.reduce((a, g) => (g.uc > a.uc ? g : a)), gSteel = groupSteel(J, gg.group), okU = gU.g.uc <= 1, okS = checks.slice(groups.length).every((c) => c.uc <= 1), ok = okU && okS;
    return {
      kpis: [
        { key: 'ucULS', label: 'Highest ULS member utilisation', unit: '-', value: gU.g.uc, limit: 1, sense: 'max', explain: 'Largest code utilisation over all ULS combinations. It must not exceed 1.0; the combination that produces it is the critical load case.' },
        { key: 'deckSway', label: 'SLS deck sway', unit: 'mm', value: gS.dp.top * 1e3, limit: swayLim * 1e3, sense: 'max', explain: 'Largest lateral displacement of the leg tops in the serviceability combinations, including foundation flexibility.' },
        { key: 'deckTilt', label: 'SLS deck tilt', unit: 'deg', value: tiltDeg, limit: inp.tiltLim, sense: 'max', explain: 'Largest inclination between any two leg tops in the serviceability combinations.' },
        { key: 'braceDefl', label: 'SLS brace deflection as span ratio (L / value)', unit: '-', value: gD.sag.ratio > 0 ? Math.min(1 / gD.sag.ratio, 1e6) : 1e6, limit: inp.deflRatio, sense: 'min', explain: 'Span divided by the mid-span deflection of the most flexible brace; larger is stiffer. Needs at least two elements per member to be evaluated.' },
        { key: 'pileComp', label: 'Largest factored pile-head compression', unit: 'MN', value: pileC / 1e6, explain: 'Largest ULS leg reaction in compression: the design action for the pile axial resistance check in a partial-factor format.' },
        { key: 'pileTens', label: 'Largest factored pile-head tension', unit: 'MN', value: pileT / 1e6, explain: 'Largest ULS uplift on a pile. It comes from the minimum-gravity combination, which is why that combination exists.' },
        { key: 'pileShear', label: 'Largest factored pile-head shear', unit: 'MN', value: pileS / 1e6, explain: 'Largest ULS horizontal reaction of one leg.' },
        { key: 'sumHmax', label: 'Largest factored base shear', unit: 'MN', value: Math.max(...uls.map((c) => c.pf.sumH)) / 1e6, explain: 'Largest resultant horizontal reaction over the ULS combinations: environmental base shear times its action factor.' },
        { key: 'nFail', label: 'ULS combinations with a member above 1.0', unit: '-', value: uls.filter((c) => c.g.uc > 1).length, limit: 0, sense: 'max', explain: 'Number of ultimate limit state combinations in which at least one member fails its code check.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Highest member utilisation by combination', categories: out.map((c) => `${c.ls}: ${c.name}`), series: [{ name: 'Utilisation', y: out.map((c) => round(c.g.uc)) }], y: { label: 'Utilisation', unit: '-' }, explain: 'Governing member utilisation of each combination (SLS shown for information with the same resistance factors). The tallest ULS bar is the critical load case.' },
        { type: 'bar', title: 'Serviceability utilisation', categories: sls.map((c) => c.name), series: [{ name: 'Sway', y: sls.map((c) => round(c.dp.top / swayLim)) }, { name: 'Tilt', y: sls.map((c) => round(c.tilt * 180 / Math.PI / inp.tiltLim)) }, { name: 'Brace deflection', y: sls.map((c) => round(c.sag.ratio * inp.deflRatio)) }], y: { label: 'Value / limit', unit: '-' }, explain: 'Each serviceability response divided by its limit for every SLS combination; bars above 1.0 exceed the project criterion.' },
        plot3d(J, { title: `ULS utilisation envelope and deformed shape of "${gS.name}"`, values: nodeValuesFromMembers(J, envU.map((r) => r.uc)), label: 'ULS utilisation', unit: '-', u: gS.r.displacements, explain: 'Envelope of ULS member utilisation at the nodes, drawn on the deflected shape of the governing serviceability combination.' }),
      ],
      tables: [
        { title: 'Results by combination', columns: [{ key: 'ls', label: 'Limit state' }, { key: 'name', label: 'Combination' }, { key: 'f', label: 'Factors G / Q / E' }, { key: 'uc', label: 'Member utilisation' }, { key: 'mem', label: 'Governing member' }, { key: 'V', label: 'Base shear', unit: 'MN' }, { key: 'pc', label: 'Pile compression', unit: 'MN' }, { key: 'pt', label: 'Pile tension', unit: 'MN' }, { key: 'sway', label: 'Deck sway', unit: 'mm' }, { key: 'tilt', label: 'Tilt', unit: 'deg' }],
          rows: out.map((c) => ({ ls: c.ls, name: c.name, f: `${c.fG} / ${c.fQ} / ${c.fE} (${c.env})`, uc: round(c.g.uc), mem: `${c.g.label} ${c.g.i}`, V: round(c.pf.sumH / 1e6), pc: round(c.pf.axialMax / 1e6), pt: round(c.pf.tensionMax / 1e6), sway: round(c.dp.top * 1e3), tilt: round(c.tilt * 180 / Math.PI, 3) })), explain: 'One row per combination with its factors and governing responses. ULS rows are judged on utilisation and pile reactions, SLS rows on sway, tilt and deflection.' },
        { title: 'Most utilised members (ULS envelope)', columns: memberCols, rows: memberTableRows(envU), explain: 'Members with the highest utilisation over all ULS combinations and the combination that governs each.' },
        { title: 'Base load cases', columns: [{ key: 'c', label: 'Case' }, { key: 'V', label: 'Horizontal resultant', unit: 'MN' }, { key: 'Fz', label: 'Vertical reaction', unit: 'MN' }], rows: base.map((b) => { const pf = pileForces(J, res[b.name]); return { c: { G: 'Permanent (jacket + permanent topsides)', Q: 'Variable topsides', Eext: 'Extreme environment', Eop: 'Operating environment' }[b.name], V: round(pf.sumH / 1e6), Fz: round(pf.sumV / 1e6) }; }), explain: 'Unfactored reactions of the four base cases that the combinations are built from.' },
      ],
      notes: [`Critical ULS combination: "${gU.name}" — utilisation ${fmt(gU.g.uc)} in ${gU.g.label.toLowerCase()} member ${gU.g.i} (${gU.g.ref}).`,
        `SLS: sway ${fmt(gS.dp.top * 1e3)} mm against ${fmt(swayLim * 1e3)} mm, tilt ${fmt(tiltDeg)}° against ${fmt(inp.tiltLim)}°, brace deflection L/${gD.sag.ratio > 0 ? fmt(1 / gD.sag.ratio) : '∞'} against L/${fmt(inp.deflRatio)}.`,
        `Extreme wave by ${dX.wf.label}, operating wave by ${dO.wf.label}, heading ${fmt(beta)}°. ${pileT > 0 ? `Pile uplift of ${fmt(pileT / 1e6)} MN occurs in the minimum-gravity combination.` : 'No pile goes into tension.'}`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `Both limit states are met: ULS utilisation ${fmt(gU.g.uc)} ("${gU.name}"), SLS sway ${fmt(gS.dp.top * 1e3)} mm of ${fmt(swayLim * 1e3)} mm.`, 'Check the other headings, then take the factored pile-head reactions to the foundation design. If serviceability has much more margin than strength, the frame is stiffness-rich and bracing could be lightened; if the reverse, foundation stiffness is the lever.', 'ISO 19902:2020 §7 and §9.10; ISO 19900:2019 §8', `The ${gg.label.toLowerCase()} govern ULS at ${fmt(gg.uc)}; sizing them to 0.9 changes steel by about ${steelCo2(gSteel * (1 - gg.uc / 0.9))}`)
          : adv('engineering', 'fail', `${okU ? 'Serviceability' : 'Ultimate strength'} is not met: ${okU ? checks.filter((c) => c.uc > 1).map((c) => `${c.name} at ${fmt(c.uc)}`).join('; ') : `utilisation ${fmt(gU.g.uc)} in "${gU.name}"`}.`, okU ? 'Stiffen the foundation or widen the base for sway and tilt (pile-head springs dominate deck movement); add a node or diameter to the flexible brace.' : `Increase the wall of the ${gg.label.toLowerCase()} by about ${fmt((gg.uc / 0.9 - 1) * 100)} % or reduce the governing action.`, BASIS_MEMBER, `The change is of the order of ${steelCo2(gSteel * Math.max(gg.uc / 0.9 - 1, 0.05))}`),
        adv('economic', 'info', `"${gU.name}" governs strength and "${gS.name}" serviceability for ${fmt(steel)} t of jacket steel.`, 'Know which combination pays for the steel: when the gravity-dominated combination governs, topsides weight control is worth more than metocean refinement; when the extreme storm governs, the reverse. Serviceability failures are cheapest to cure in the foundation, not in jacket steel.', `${BASIS_COST}; ISO 19901-5:2021 (weight management)`),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Cantilever tube (L = 95 m, EI = 1.14418e10 N·m²) with a 10 kN tip load: SLS sway = PL³/3EI + PL/(GA/2) = 249.89 mm (factor 1.0);
    // the ULS combination with γE = 1.5 gives a factored base shear of 15 kN.
    { name: 'Factored combination of a tip load: sway PL³/3EI at γ = 1.0, base shear γE·P', source: 'Euler–Bernoulli cantilever (Roark Table 8.1 case 1a) and linear superposition',
      inputs: { geom: testFrame('column'), H: 0, Hop: 0, Uc: 0, UcOp: 0, Uw: 0, UwOp: 0, windA: 0, swCont: 0, Wtop: 0, FxTop: 10, found: 'fixed', nSub: 2, combos: [{ name: 'U', ls: 'ULS', fG: 1.1, fQ: 1.1, fE: 1.5, env: 'extreme' }, { name: 'S', ls: 'SLS', fG: 1, fQ: 1, fE: 1, env: 'extreme' }] },
      expect: [{ kpi: 'deckSway', value: 249.89, tol: 0.002 }, { kpi: 'sumHmax', value: 0.015, tol: 1e-6 }] },
  ],
};
