// In-place analysis of the jacket under operating and extreme environmental conditions: wave + current + wind by
// Morison on every wetted member, stepped through the wave phase for each heading, with topsides and gravity loads.
import { need, fmt, round, adv, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, BASIS_MEMBER, BASIS_COST } from './_common.js';
import { geomInputs, envInputs, topsidesInputs, foundationInputs, criteriaInputs, meshInputs, headingInput, headingList, buildJacket, conditionCases, loadCase, solve,
  memberRows, envelope, byGroup, pileForces, displacements, nodeValuesFromMembers, plot3d, memberCols, memberTableRows } from './_model.js';

export default {
  id: 's01-inplace', suite: 1, title: 'In-place analysis — operating and extreme storm',
  summary: 'Base shear, overturning moment, member utilisation, displacements and pile-head reactions of the jacket for the operating and the extreme wave, current and wind condition over several headings.',
  covers: [1], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §9 (actions), §12 (in-place), §13 (members)', 'ISO 19901-1:2015 (metocean)', 'API RP 2A-WSD 22nd ed. §5.3, §6', 'DNV-RP-C205 (2019) §6.2 (Morison)'],
  models: ['Direct stiffness method (3-D Euler–Bernoulli / Timoshenko beam elements)', 'Morison equation with drag and inertia', 'Stokes 5th-order, stream-function and Airy (Wheeler-stretched) wave kinematics', 'Power-law current profile with blockage', 'Wind force on the topsides area', 'Linear pile-head springs (6×6 with lateral–rocking coupling)', 'ISO 19902 tubular beam-column interaction'],
  limits: `Quasi-static linear elastic frame with rigid joints; one regular design wave per condition stepped through its phase (no irregular sea, no dynamic amplification — see the dynamic module). ${MORISON_CAVEAT} ${CODE_CAVEAT} Joint cans, conductors and appurtenances are represented only through the self-weight contingency and the hydrodynamic coefficients.`,
  publishes: ['jacket.baseShear', 'jacket.otm', 'jacket.pileAxialMax', 'jacket.pileShearMax', 'quantities.steel_t'],
  searchTerms: 'jacket in-place analysis extreme storm base shear overturning moment',
  inputs: [
    ...geomInputs(), ...topsidesInputs(),
    ...envInputs({ heading: false, group: 'Environment — extreme (100-year)' }).slice(0, 5), headingInput(),
    { key: 'Hop', label: 'Operating wave height', unit: 'm', type: 'number', default: 10, min: 0, max: 40, group: 'Environment — operating (1-year)', role: 'bc', help: 'Maximum individual wave of the 1-year sea state (about 1.86 Hs for a 3-hour storm).' },
    { key: 'Top', label: 'Operating wave period', unit: 's', type: 'number', default: 10.5, min: 3, max: 25, group: 'Environment — operating (1-year)', role: 'bc', bind: 'site.Tp1' },
    { key: 'UcOp', label: 'Operating surface current', unit: 'm/s', type: 'number', default: 0.6, min: 0, max: 4, group: 'Environment — operating (1-year)', role: 'bc', bind: 'site.Uc1' },
    { key: 'UwOp', label: 'Operating 1-hour wind at 10 m', unit: 'm/s', type: 'number', default: 22, min: 0, max: 80, group: 'Environment — operating (1-year)', role: 'bc', bind: 'site.Uw1' },
    ...envInputs({ heading: false }).slice(5),
    ...foundationInputs(), ...criteriaInputs(),
    { key: 'gGop', label: 'Operating: action factor on gravity loads', unit: '-', type: 'number', default: 1.3, min: 0.8, max: 2, group: 'Design criteria', role: 'criteria', help: 'ISO 19902 operating combination 1.3 G + 1.5 Q + 0.9·1.35 E, with G and Q lumped here. Confirm against the governing edition.' },
    { key: 'gEop', label: 'Operating: action factor on environmental loads', unit: '-', type: 'number', default: 1.215, min: 0.5, max: 2.5, group: 'Design criteria', role: 'criteria' },
    ...meshInputs(),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxDisp' },
  meshTime: { key: 'nPhase', label: 'Wave phase steps', levels: [8, 16, 32], qoi: 'baseShear' },
  compute(inp, ctx) {
    // SI conversion: mm → m, MN → N, kN → N inside the shared model.
    const J = buildJacket(inp, ctx), heads = headingList(inp.nDir);
    need(inp.H >= inp.Hop, 'The extreme wave height must not be smaller than the operating wave height.');
    const ext = conditionCases(J, inp, { tag: 'Extreme', H: inp.H, T: inp.T, Uc: inp.Uc, Uw: inp.Uw, fG: inp.gG, fE: inp.gE, headings: heads });
    const op = conditionCases(J, inp, { tag: 'Operating', H: inp.Hop, T: inp.Top, Uc: inp.UcOp, Uw: inp.UwOp, fG: inp.gGop, fE: inp.gEop, headings: heads });
    // Unfactored extreme cases give the characteristic pile-head loads handed to the foundation design.
    const chr = ext.dirs.map((d) => loadCase(`char ${d.beta}`, J, inp, { betaDeg: d.beta, env: d.ld.at(d.sc.t[d.sc.iM], true) }));
    ctx.progress?.(0.4, 'Solving load cases');
    const res = solve(J, [...ext.cases, ...op.cases, ...chr]);
    const table = [], all = { Extreme: [], Operating: [] };
    let gov = null;
    for (const [tag, cond] of [['Extreme', ext], ['Operating', op]]) for (const d of cond.dirs) {
      const rows = envelope(d.names.map((n) => memberRows(J, res[n], inp, { name: n })));
      const pf = d.names.map((n) => pileForces(J, res[n])), dp = d.names.map((n) => displacements(J, res[n]));
      const g = rows.reduce((a, r) => (r.uc > a.uc ? r : a));
      const row = { tag, beta: d.beta, baseShear: d.baseShear, otm: d.otm, uc: g.uc, member: g, axial: Math.max(...pf.map((p) => p.axialMax)), tension: Math.max(...pf.map((p) => p.tensionMax)), shear: Math.max(...pf.map((p) => p.shearMax)), top: Math.max(...dp.map((x) => x.top)), max: Math.max(...dp.map((x) => x.max)), rows, name: g.case };
      table.push(row); all[tag].push(rows);
      if (!gov || row.uc > gov.uc) gov = row;
    }
    const envExt = envelope(all.Extreme), envOp = envelope(all.Operating), envAll = envelope([envExt, envOp]);
    const tExt = table.filter((r) => r.tag === 'Extreme'), tOp = table.filter((r) => r.tag === 'Operating');
    const mx = (list, k) => Math.max(...list.map((r) => r[k]));
    const pfChr = chr.map((c) => pileForces(J, res[c.name]));
    const pileAx = Math.max(...pfChr.map((p) => p.axialMax)), pileTen = Math.max(...pfChr.map((p) => p.tensionMax)), pileSh = Math.max(...pfChr.map((p) => p.shearMax));
    const baseShear = mx(tExt, 'baseShear'), otm = mx(tExt, 'otm'), maxUc = mx(tExt, 'uc'), maxUcOp = mx(tOp, 'uc');
    const govRes = res[gov.name], dGov = displacements(J, govRes), groups = byGroup(envAll), gExt = byGroup(envExt), gOp = byGroup(envOp);
    const steel = J.steelT * (inp.swCont > 0 ? inp.swCont : 1), gg = groups.reduce((a, g) => (g.uc > a.uc ? g : a));
    const gSteel = J.orig.filter((m) => m.group === gg.group).reduce((s, m) => s + m.sec.A * m.L * J.steel.rho, 0) / 1e3;
    const ok = maxUc <= 1 && maxUcOp <= 1;
    const checks = [
      ...gExt.map((g) => ({ name: `Member strength, extreme — ${g.label}`, demand: round(g.uc), capacity: 1, unit: 'interaction', uc: g.uc, ref: `${inp.code} member check (${g.ref}); ${g.case}`, explain: 'Governing tension, compression, bending and beam-column interaction of this member group in the 100-year storm. Above 1.0 the member yields or buckles before the factored load is reached.' })),
      ...gOp.map((g) => ({ name: `Member strength, operating — ${g.label}`, demand: round(g.uc), capacity: 1, unit: 'interaction', uc: g.uc, ref: `${inp.code} member check (${g.ref}); ${g.case}`, explain: 'Same check under the operating storm with the higher gravity-load factor: it governs members that mainly carry deck weight.' })),
    ];
    const ph = ext.dirs[0].sc;
    return {
      kpis: [
        { key: 'baseShear', label: 'Extreme base shear (unfactored)', unit: 'MN', value: baseShear / 1e6, explain: 'Largest total horizontal wave, current and wind force over the headings in the 100-year condition. It sizes the piles laterally and is the reference load for the reserve strength ratio.' },
        { key: 'otm', label: 'Extreme overturning moment at mudline (unfactored)', unit: 'MN·m', value: otm / 1e6, explain: 'Largest moment of the environmental loads about the mudline centre. Divided by the leg spacing it gives the push–pull on the piles.' },
        { key: 'baseShearOp', label: 'Operating base shear (unfactored)', unit: 'MN', value: mx(tOp, 'baseShear') / 1e6, explain: 'Total horizontal environmental force in the operating (1-year) condition.' },
        { key: 'maxUc', label: 'Highest member utilisation — extreme', unit: '-', value: maxUc, limit: 1, sense: 'max', explain: 'Largest code utilisation of any member under the factored extreme combination. Above 1.0 the member is overstressed.' },
        { key: 'maxUcOp', label: 'Highest member utilisation — operating', unit: '-', value: maxUcOp, limit: 1, sense: 'max', explain: 'Largest code utilisation under the factored operating combination.' },
        { key: 'pileAxialMax', label: 'Largest pile-head compression (characteristic)', unit: 'MN', value: pileAx / 1e6, explain: 'Largest unfactored axial reaction of a leg at the mudline in the extreme storm: the demand on pile axial capacity.' },
        { key: 'pileTensionMax', label: 'Largest pile-head tension (characteristic)', unit: 'MN', value: pileTen / 1e6, explain: 'Largest unfactored uplift on a pile; zero means all piles stay in compression. Tension needs shaft friction only and usually sets penetration.' },
        { key: 'pileShearMax', label: 'Largest pile-head shear (characteristic)', unit: 'MN', value: pileSh / 1e6, explain: 'Largest unfactored horizontal reaction of one leg: the head load for the lateral pile analysis.' },
        { key: 'maxDisp', label: 'Largest nodal displacement (governing case)', unit: 'mm', value: dGov.max * 1e3, explain: 'Largest translation of any node in the governing factored case, including foundation flexibility.' },
        { key: 'deckDisp', label: 'Deck-level lateral displacement (governing case)', unit: 'mm', value: dGov.top * 1e3, explain: 'Horizontal movement of the leg tops; it drives conductor, riser and bridge interface design.' },
        { key: 'sumV', label: 'Total vertical reaction (governing case)', unit: 'MN', value: pileForces(J, govRes).sumV / 1e6, explain: 'Sum of vertical pile-head reactions: must equal factored topsides weight plus buoyant jacket weight (equilibrium check).' },
        { key: 'steelT', label: 'Jacket steel mass', unit: 't', value: steel, explain: 'Modelled tubular steel times the self-weight contingency; piles and topsides are not included.' },
      ],
      checks,
      plots: [
        plot3d(J, { title: 'Member utilisation (envelope) and deformed shape of the governing case', values: nodeValuesFromMembers(J, envAll.map((r) => r.uc)), label: 'Utilisation', unit: '-', u: govRes.displacements, explain: `Envelope of code utilisation over all cases, shown at the nodes, with the deflected shape for ${gov.name}. Look for a hot bay or a leg that works much harder than the others.` }),
        { type: 'line', title: 'Environmental force through the wave cycle (extreme)', x: { label: 'Time in wave period', unit: 's' }, y: { label: 'Hydrodynamic base shear', unit: 'MN' },
          series: ext.dirs.map((d) => ({ name: `${d.beta}°`, x: d.sc.t.map((v) => round(v)), y: d.sc.Fh.map((v) => round(v / 1e6)) })),
          explain: 'Wave and current force on the jacket as the crest passes (crest at the centre at t = 0). The peak just before or at the crest is the design base shear; a flat peak means the phase step is fine enough.' },
        { type: 'polar', title: 'Base shear by heading', series: [{ name: 'Extreme (MN)', thetaDeg: tExt.map((r) => r.beta), r: tExt.map((r) => round(r.baseShear / 1e6)) }, { name: 'Operating (MN)', thetaDeg: tOp.map((r) => r.beta), r: tOp.map((r) => round(r.baseShear / 1e6)) }],
          explain: 'Total base shear against the direction the wave travels. Diagonal headings usually load a single corner pile hardest even when the base shear is similar.' },
        { type: 'bar', title: 'Utilisation by member group', categories: groups.map((g) => g.label), series: [{ name: 'Extreme', y: gExt.map((g) => round(g.uc)) }, { name: 'Operating', y: gOp.map((g) => round(g.uc)) }], y: { label: 'Utilisation', unit: '-' },
          explain: 'Highest utilisation in each member group for the two conditions. An uneven picture shows where steel can be moved from one group to another.' },
      ],
      tables: [
        { title: 'Load cases by condition and heading', columns: [{ key: 'tag', label: 'Condition' }, { key: 'beta', label: 'Heading', unit: 'deg' }, { key: 'bs', label: 'Base shear', unit: 'MN' }, { key: 'otm', label: 'OTM', unit: 'MN·m' }, { key: 'uc', label: 'Max utilisation' }, { key: 'mem', label: 'Governing member' }, { key: 'ax', label: 'Pile compression (factored)', unit: 'MN' }, { key: 'ten', label: 'Pile tension (factored)', unit: 'MN' }, { key: 'sh', label: 'Pile shear (factored)', unit: 'MN' }, { key: 'top', label: 'Deck displacement', unit: 'mm' }],
          rows: table.map((r) => ({ tag: r.tag, beta: r.beta, bs: round(r.baseShear / 1e6), otm: round(r.otm / 1e6), uc: round(r.uc), mem: `${r.member.label} ${r.member.i}`, ax: round(r.axial / 1e6), ten: round(r.tension / 1e6), sh: round(r.shear / 1e6), top: round(r.top * 1e3) })),
          explain: 'One row per condition and heading: unfactored base shear and overturning moment, and the factored member and pile results. The critical load case is the row with the highest utilisation.' },
        { title: 'Most utilised members (envelope)', columns: memberCols, rows: memberTableRows(envAll), explain: 'The fifteen most utilised members with their extreme section forces and the case that governs each: the list to take into joint and local checks.' },
        { title: 'Pile-head reactions, governing case', columns: [{ key: 'node', label: 'Leg node' }, { key: 'axial', label: 'Axial (compression +)', unit: 'MN' }, { key: 'shear', label: 'Shear', unit: 'MN' }, { key: 'moment', label: 'Moment', unit: 'MN·m' }],
          rows: pileForces(J, govRes).piles.map((p) => ({ node: p.node, axial: round(p.axial / 1e6), shear: round(p.shear / 1e6), moment: round(p.moment / 1e6) })), explain: `Factored reactions at each leg for ${gov.name}. Negative axial values are uplift on that pile.` },
      ],
      notes: [`Governing: ${gov.name}, ${gov.member.label.toLowerCase()} member ${gov.member.i} (${gov.member.ref}), utilisation ${fmt(gov.uc)}.`, `Extreme wave by ${ext.wf.label}; operating wave by ${op.wf.label}. Largest leg diameter over wavelength D/λ = ${fmt(Math.max(...J.orig.map((m) => m.D)) / ext.wf.L)} (Morison needs < 0.2).`,
        `Foundation: ${J.K ? `pile-head springs (${J.K.source}) Kh ${fmt(J.K.Kh / 1e6)} MN/m, Kr ${fmt(J.K.Kr / 1e9)} GN·m/rad, Kv ${fmt(J.K.Kv / 1e6)} MN/m` : `${J.found} supports at the mudline`}.`,
        pileTen > 0 ? `A pile goes into tension (${fmt(pileTen / 1e6)} MN): check tensile shaft capacity in the pile module.` : 'All piles remain in compression in the extreme storm.'],
      advice: [
        ok ? adv('engineering', 'ok', `All members pass: highest utilisation ${fmt(Math.max(maxUc, maxUcOp))} in the ${gg.label.toLowerCase()} (${gov.name}); extreme base shear ${fmt(baseShear / 1e6)} MN, overturning moment ${fmt(otm / 1e6)} MN·m.`,
          gg.uc < 0.6 ? `The ${gg.label.toLowerCase()} have reserve: wall thickness could come down by about ${fmt((1 - gg.uc / 0.85) * 100)} % before reaching 0.85. Run joints, fatigue and the installation cases first — they often govern over in-place strength.` : 'Utilisation is in an efficient range. Confirm tubular joints, fatigue at the splash-zone braces and the pile capacities with the loads published here.',
          BASIS_MEMBER, gg.uc < 0.6 ? `Trimming the ${gg.label.toLowerCase()} to 0.85 utilisation would save roughly ${steelCo2(gSteel * (1 - gg.uc / 0.85))}` : `The jacket holds ${steelCo2(steel)}`)
          : adv('engineering', 'fail', `Member strength is exceeded: utilisation ${fmt(Math.max(maxUc, maxUcOp))} in ${gov.member.label.toLowerCase()} member ${gov.member.i} (${gov.member.ref}) for ${gov.name}.`,
            `Increase the wall thickness of the ${gg.label.toLowerCase()} by about ${fmt((gg.uc / 0.9 - 1) * 100)} % (to reach 0.9), or shorten the buckling length with an extra bay; a larger diameter also raises wave load, so thickness or grade is the first lever.`,
            BASIS_MEMBER, `Strengthening the ${gg.label.toLowerCase()} adds about ${steelCo2(gSteel * Math.max(gg.uc / 0.9 - 1, 0.05))} A higher grade reaches the same utilisation with less mass where buckling does not govern.`),
        adv('economic', ok ? 'info' : 'warn', `The jacket carries ${fmt(steel)} t of primary steel for ${fmt(baseShear / 1e6)} MN of base shear; the ${gg.label.toLowerCase()} (${fmt(gSteel)} t) set the governing utilisation.`,
          'Fabricated jacket steel is paid, transported and lifted per tonne: keep the lift weight under the capacity step of the intended crane vessel and standardise on a few can sizes. Drag-dominated loading scales with member diameter, so removing redundant plan bracing and appurtenances near the surface cuts both steel and wave load.',
          `${BASIS_COST}; lift-weight control per DNV-ST-N001 (2021)`),
      ],
      publish: { 'jacket.baseShear': baseShear, 'jacket.otm': otm, 'jacket.pileAxialMax': pileAx, 'jacket.pileShearMax': pileSh, 'quantities.steel_t': steel },
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    { name: 'Global equilibrium without environment: ΣRz = W, base shear = topsides interface load', source: 'Statics: 40 MN deck weight on four legs with self-weight off gives ΣRz = 40 MN; a 1000 kN interface load is the only horizontal force, so the base shear is 1.0 MN and the overturning moment 1.0 MN × 95 m = 95 MN·m',
      inputs: { H: 0, Hop: 0, Uc: 0, UcOp: 0, Uw: 0, UwOp: 0, windA: 0, swCont: 0, found: 'fixed', Wtop: 40, FxTop: 1000, gG: 1, gE: 1, gGop: 1, gEop: 1, nDir: '1' },
      expect: [{ kpi: 'sumV', value: 40, tol: 1e-6 }, { kpi: 'baseShear', value: 1.0, tol: 1e-6 }, { kpi: 'otm', value: 95, tol: 1e-6 }] },
  ],
};
