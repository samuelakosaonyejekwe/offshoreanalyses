// Reserve and residual strength of the topsides frame: event-to-event plastic collapse (pushover) analysis under
// gravity or lateral load, reserve strength ratio and redundancy, and the damaged condition with one member
// removed (dropped-object, blast, fire or fatigue loss) to show the nonlinear redistribution of load.
import { damagedStrength, progressiveCollapse } from '../../engine/structural/pushover.js';
import { deckInputs, buildDeck, loadCase, solveDeck, memberForces, memberChecks, equipmentOf, need, run, fmt, round, thin, adv, fsCheck, steelCo2Text, CONFIRM, DECK_LIMITS } from './_common.js';

const GROUPS = ['Brace', 'Column', 'Deck leg', 'Girder', 'Beam'];

export default {
  id: 's02-reserve-strength', suite: 2, title: 'Reserve strength and damaged-condition redundancy',
  summary: 'Plastic collapse load factor of the deck frame, reserve strength ratio and redundancy, and the residual strength after losing a member to an accidental event or deterioration.',
  covers: [10, 14], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §7.10, §12.6 (robustness, non-linear analysis, reserve strength ratio)', 'ISO 19901-3:2014 §8 (accidental situations: post-damage resistance)', 'NORSOK N-004 (2022) Annex A.2 (accidental limit state, two-step check)', 'NORSOK N-006 (2015) §7 (damaged-condition assessment)', 'DNV-RP-C204 (2019) §2 (residual strength after accidental actions)'],
  models: ['Plastic-hinge formulation (event-to-event)', 'Plastic collapse analysis', 'Progressive failure simulation (member removal)', 'Nonlinear load redistribution', 'Residual strength assessment', 'Second-order (P–Δ) geometric stiffness', 'ISO 19902 column curve with post-buckling residual capacity'],
  limits: `${DECK_LIMITS} Members are elastic–perfectly-plastic with concentrated hinges at their ends and an axial limit from the ISO 19902 column curve with a fixed post-buckling residual; hinges within the span, joint failure, strain hardening, local buckling of hinges and the dynamics of sudden member loss (other than through a load factor) are not modelled. Explicit nonlinear finite-element collapse analysis is outside this module.`,
  searchTerms: 'pushover reserve strength ratio redundancy damaged condition residual strength offshore deck frame',
  inputs: [
    ...deckInputs(),
    { key: 'pattern', label: 'Load increased to collapse', type: 'select', default: 'gravity', group: 'Loading', role: 'bc', options: [{ value: 'gravity', label: 'Gravity (dead, equipment and live together)' }, { value: 'lateral', label: 'Lateral load along x, with gravity held' }] },
    { key: 'Hx', label: 'Reference lateral load along x', unit: 'kN', type: 'number', default: 900, min: 1, max: 5e4, group: 'Loading', role: 'bc', bind: 'topsides.windFx', showIf: { key: 'pattern', in: ['lateral'] }, help: 'Design lateral load at load factor 1.0, shared between levels in proportion to elevation.' },
    { key: 'gG', label: 'Factor on permanent actions at the design level', unit: '-', type: 'number', default: 1.3, min: 0.8, max: 1.6, group: 'Loading', role: 'criteria', help: `The collapse load factor multiplies this design level. For the damaged condition all factors are normally 1.0. ${CONFIRM}` },
    { key: 'gQ', label: 'Factor on variable actions at the design level', unit: '-', type: 'number', default: 1.5, min: 0, max: 1.8, group: 'Loading', role: 'criteria', help: CONFIRM },
    { key: 'damage', label: 'Member group damaged', type: 'select', default: 'Brace', group: 'Damaged condition', role: 'ic', options: GROUPS.map((v) => ({ value: v, label: v })), help: 'The most utilised member of this group (from the elastic analysis) is removed, representing loss to impact, blast, fire or a through-thickness crack.' },
    { key: 'pDelta', label: 'Include P–Δ effects', type: 'bool', default: false, group: 'Solver', role: 'control', help: 'Adds the geometric stiffness of the current axial forces at every event; needed for sway collapse of unbraced frames, and slower.' },
    { key: 'postBuck', label: 'Post-buckling residual compression capacity', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Solver', role: 'material', help: 'Fraction of the buckling load a strut still carries after it buckles (0.3–0.5 for typical tubular braces).' },
    { key: 'screen', label: 'Screen every member of the group (slower)', type: 'bool', default: false, group: 'Solver', role: 'control', help: 'Removes each member of the damaged group in turn and reports the worst case; several seconds for large frames.' },
    { key: 'rsrReq', label: 'Required collapse load factor, intact', unit: '-', type: 'number', default: 1.3, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: `Margin on the factored design loads at plastic collapse. ISO 19902 quotes RSR targets of about 1.85 for unfactored extreme environmental load on L1 jackets; for factored loads a value somewhat above 1.0 shows robustness. ${CONFIRM}` },
    { key: 'dmgReq', label: 'Required collapse load factor, damaged', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 3, group: 'Acceptance criteria', role: 'criteria', help: `The damaged structure must carry the loads present after the event; with the factored design level entered above, 1.0 is conservative. ${CONFIRM}` },
  ],
  compute(inp, ctx) {
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp), zs = deck.levelsZ, zsum = zs.reduce((s, z) => s + z, 0), lat = inp.pattern === 'lateral';
    const grav = loadCase(deck, inp, { name: 'grav', gDead: inp.gG, gEquip: inp.gG, gLive: inp.gQ, eqs }), nodal = [];
    deck.levelNodes.forEach((ns, k) => { for (const n of ns) nodal.push({ node: n, F: [(inp.Hx * 1e3 * zs[k]) / zsum / ns.length, 0, 0, 0, 0, 0] }); });
    // Elastic analysis picks the member to remove: the most utilised of the chosen group under the design loads.
    const el = solveDeck(deck, lat ? [{ ...grav, name: 'design', nodal: [...grav.nodal, ...nodal] }] : [{ ...grav, name: 'design' }], 1);
    // The chosen group may not exist in an imported frame: fall back to the most utilised member of any group.
    const rows = memberChecks(deck, memberForces(deck, el.ref, el.sol.results.design)), own = rows.filter((r) => r.group === inp.damage), inGroup = own.length ? own : rows;
    const victim = inGroup.reduce((a, c) => (c.uc > a.uc ? c : a)), ucEl = Math.max(...rows.map((r) => r.uc));
    const Wdes = Math.abs(el.sol.results.design.reactions.reduce((sum, r) => sum + r.R[2], 0));   // factored gravity on the frame [N]
    // Generic imported sections have no plastic moduli: their capacities are set from the elastic modulus at an
    // extreme fibre of 1.4 radii of gyration (conservative for bending), squash load A·fy.
    const fy = deck.steel.fy, members = deck.model.members.map((m, i) => { const o = deck.orig[i]; if (o.kind !== 'other' || o.sec.Wply) return m; const cy = 1.4 * Math.sqrt(o.sec.Iy / o.sec.A), cz = 1.4 * Math.sqrt(o.sec.Iz / o.sec.A); return { ...m, capacity: { Np: o.sec.A * fy, Mpy: (fy * o.sec.Iy) / cy, Mpz: (fy * o.sec.Iz) / cz } }; });
    const model = { ...deck.model, members, loadCases: [grav, { name: 'lat', nodal }] };
    const opts = lat ? { pattern: 'lat', base: 'grav', pDelta: inp.pDelta, postBucklingRatio: inp.postBuck } : { pattern: 'grav', pDelta: inp.pDelta, postBucklingRatio: inp.postBuck };
    const ds = run(() => damagedStrength(model, { remove: [victim.i], ...opts }), 'Collapse analysis (check that the frame is stable under the held loads)');
    const A = ds.intact, D = ds.damaged, lamI = A.ultimateLoadFactor, lamD = D.ultimateLoadFactor;
    need(Number.isFinite(lamI) && lamI > 0, 'The intact frame has no positive collapse load factor: it is unstable under the held gravity loads.');
    let scr = null;
    if (inp.screen) scr = run(() => progressiveCollapse(model, { members: inGroup.map((r) => r.i), daf: inp.dmgReq, ...opts }), 'Member-removal screening');
    const lamDmin = scr ? Math.min(lamD, scr.critical.ultimateLoadFactor) : lamD;
    const checks = [
      fsCheck('Intact frame: collapse load factor', lamI, inp.rsrReq, 'ISO 19902:2020 §7.10, §12.6 (reserve strength by non-linear analysis)', `Multiplier on the ${lat ? 'design lateral load (with factored gravity held)' : 'factored gravity loads'} at which a plastic mechanism forms. Below the requirement the frame has too little reserve beyond first member failure.`),
      fsCheck(`Damaged frame (${inp.damage.toLowerCase()} ${victim.i + 1} removed${scr ? ', worst of the screened group' : ''}): collapse load factor`, Math.max(lamDmin, 1e-6), inp.dmgReq, 'NORSOK N-004 Annex A.2 (ALS step 2: resistance in the damaged condition); ISO 19901-3 §8', 'Collapse load factor after the member is lost. Below the requirement the loss of one member leads to progressive collapse: the member is non-redundant and must be protected, duplicated or designed to survive the event.'),
    ];
    const ok = checks.every((c) => c.uc <= 1), gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a));
    // Member state plot of the intact frame at collapse: 0 elastic, 1 plastic hinge, 2 axial limit (buckled or yielded).
    const P = deck.model.nodes, st = A.memberStates, val = new Array(P.length).fill(0);
    deck.model.members.forEach((m, i) => { const s = st[i]; if (!s) return; const v = s.axial && s.axial !== 'elastic' ? 2 : 0; val[m.a] = Math.max(val[m.a], v, s.hingeA ? 1 : 0); val[m.b] = Math.max(val[m.b], v, s.hingeB ? 1 : 0); });
    let umax = 0, size = 0; for (let i = 0; i < P.length; i++) { umax = Math.max(umax, Math.hypot(A.displacements[6 * i], A.displacements[6 * i + 1], A.displacements[6 * i + 2])); size = Math.max(size, ...P[i].map(Math.abs)); }
    const sc = umax > 0 ? 0.06 * size / umax : 0;
    const curve = (r) => { const [x, y] = thin([r.curve.displacement.map((v) => Math.abs(v) * 1e3), r.curve.loadFactor], 300); return { x: x.map((v) => round(v)), y: y.map((v) => round(v)) }; };
    const ci = curve(A), cd = curve(D), nHinge = st.filter((s) => s && (s.hingeA || s.hingeB)).length, nAx = st.filter((s) => s && s.axial && s.axial !== 'elastic').length;
    const gT = deck.orig.filter((o) => o.group === inp.damage).reduce((s, o) => s + o.sec.A * o.L * deck.steel.rho, 0) / 1e3;
    return {
      kpis: [
        { key: 'lamUlt', label: 'Collapse load factor, intact', unit: '-', value: lamI, limit: inp.rsrReq, sense: 'min', explain: `Multiplier on the ${lat ? 'design lateral load' : 'factored gravity loads'} at plastic collapse of the intact frame (the reserve strength ratio against the design level entered).` },
        { key: 'lamFirst', label: 'Load factor at first member failure', unit: '-', value: A.firstEventLoadFactor, explain: 'Load factor at which the first hinge forms or the first strut buckles: about the inverse of the elastic utilisation.' },
        { key: 'redundancy', label: 'Redundancy factor', unit: '-', value: ds.redundancyFactor, explain: 'Collapse load over the load at first member failure. 1.0 means the first failure is the collapse (no redistribution); braced, statically indeterminate frames reach 1.2–2.' },
        { key: 'lamDamaged', label: 'Collapse load factor, damaged', unit: '-', value: lamDmin, limit: inp.dmgReq, sense: 'min', explain: 'Collapse load factor with the selected member removed: the residual strength after the accidental event.' },
        { key: 'dsr', label: 'Damaged strength ratio', unit: '-', value: lamDmin / lamI, explain: 'Residual over intact collapse load. Close to 1.0 the member is redundant; a low value marks a member whose loss the frame cannot absorb.' },
        { key: 'ultLoad', label: lat ? 'Lateral load at collapse' : 'Total gravity load at collapse', unit: 'kN', value: lamI * (lat ? inp.Hx * 1e3 : Wdes) / 1e3, explain: lat ? 'Collapse load factor times the reference lateral load: the ultimate lateral resistance of the frame.' : 'Collapse load factor times the factored gravity load on the frame.' },
        { key: 'ucElastic', label: 'Highest elastic member utilisation', unit: '-', value: round(ucEl, 5), explain: 'Code utilisation of the most loaded member at the design level, for comparison with the load factor at first failure.' },
        { key: 'events', label: 'Member failures before collapse', unit: '-', value: A.events.length, explain: 'Number of hinge formations and strut failures between first yield and the collapse mechanism of the intact frame.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Collapse (pushover) curves', x: { label: 'Control displacement', unit: 'mm' }, y: { label: 'Load factor', unit: '-' }, series: [{ name: 'Intact', ...ci }, { name: `${inp.damage} ${victim.i + 1} removed`, ...cd }, { name: 'Required, intact', x: [0, Math.max(...ci.x, 1)], y: [inp.rsrReq, inp.rsrReq], dash: true }, { name: 'Design level', x: [0, Math.max(...ci.x, 1)], y: [1, 1], dash: true }],
          explain: 'Load factor against displacement as members yield and buckle one after another. Each kink is a member failure; the plateau is the collapse load. The gap between the curves is the strength lost with the damaged member.' },
        { type: 'model3d', title: 'Member states at collapse (intact frame)', nodes: P.map((p) => [round(p[0], 6), round(p[1], 6), round(p[2], 6)]), lines: deck.model.members.map((m) => [m.a, m.b]), nodeValues: val, valueLabel: 'State (0 elastic, 1 plastic hinge, 2 strut failed)', valueUnit: '-', deformed: P.map((p, i) => [round(p[0] + sc * A.displacements[6 * i], 6), round(p[1] + sc * A.displacements[6 * i + 1], 6), round(p[2] + sc * A.displacements[6 * i + 2], 6)]),
          explain: `Where the frame has yielded when the mechanism forms: ${nHinge} member(s) with plastic hinges and ${nAx} at their axial limit. The pattern shows the collapse mechanism and which members to strengthen for more reserve. Deformed shape exaggerated ${fmt(sc, 2)} times.` },
      ],
      tables: [
        { title: 'Sequence of member failures (intact frame)', columns: [{ key: 'n', label: 'Event' }, { key: 'lam', label: 'Load factor' }, { key: 'member', label: 'Member' }, { key: 'group', label: 'Group' }, { key: 'type', label: 'Failure' }, { key: 'd', label: 'Displacement', unit: 'mm' }],
          rows: A.events.slice(0, 40).map((e, i) => ({ n: i + 1, lam: round(e.loadFactor), member: e.member + 1, group: deck.orig[e.member]?.group ?? '-', type: e.type, d: round(Math.abs(e.displacement) * 1e3) })), explain: 'Order in which hinges form and struts fail as the load rises (first 40 events). The first rows identify the weakest members; a long list before collapse means good redistribution.' },
        ...(scr ? [{ title: `Member-removal screening — ${inp.damage} group`, columns: [{ key: 'member', label: 'Member removed' }, { key: 'lam', label: 'Collapse load factor' }, { key: 'uc', label: 'Required / achieved' }, { key: 'ok', label: 'Survives' }], rows: scr.cases.map((c) => ({ member: c.member + 1, lam: round(c.ultimateLoadFactor), uc: round(c.uc), ok: c.survives ? 'yes' : 'no' })), explain: 'Each member of the group removed in turn with the collapse load factor of what remains. Members marked "no" are critical: their loss leads to collapse at less than the required load.' }] : []),
      ],
      notes: [`Intact collapse at load factor ${fmt(lamI)} (${A.status}); first failure at ${fmt(A.firstEventLoadFactor)} in member ${A.events[0] ? A.events[0].member + 1 : '-'}.`, `With ${own.length ? inp.damage.toLowerCase() : `member (no '${inp.damage}' group in this frame: most utilised ${victim.group.toLowerCase()})`} ${victim.i + 1} removed (elastic utilisation ${fmt(victim.uc)}) the frame collapses at ${fmt(lamD)}: damaged strength ratio ${fmt(lamD / lamI)}.`,
        inp.pDelta ? 'P–Δ effects are included.' : 'P–Δ effects are not included: switch them on for unbraced frames collapsing in sway.'],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The frame reaches ${fmt(lamI)} times the design load before collapse (redundancy ${fmt(ds.redundancyFactor)}) and ${fmt(lamDmin)} with ${inp.damage.toLowerCase()} ${victim.i + 1} lost (${fmt(lamDmin / lamI * 100)} % of the intact strength).`,
          'Use the failure sequence to place inspection effort and passive protection on the first members to fail, and confirm joints and connections can develop the member capacities assumed (joint failure would cut the reserve shown here).', 'ISO 19902:2020 §7.10, §12.6; NORSOK N-004 (2022) Annex A.2',
          `Redundancy is cheap robustness: the ${inp.damage.toLowerCase()} group weighs ${fmt(gT)} t (${steelCo2Text(gT)}) and gives a reserve that would otherwise need heavier primary members throughout.`)
          : adv('engineering', 'fail', `${gov.name} is ${fmt(gov.capacity)} against ${fmt(gov.demand)} required: ${lamI < inp.rsrReq ? 'the intact frame has too little reserve' : `losing ${inp.damage.toLowerCase()} ${victim.i + 1} leaves ${fmt(lamDmin / lamI * 100)} % of the intact strength`}.`,
            lamI < inp.rsrReq ? `Strengthen the members that fail first (member ${A.events[0] ? A.events[0].member + 1 : '-'} at load factor ${fmt(A.firstEventLoadFactor)}) by about ${fmt((inp.rsrReq / lamI - 1) * 100, 2)} %, or add a braced bay so that load can redistribute.` : `Provide an alternative load path around ${inp.damage.toLowerCase()} ${victim.i + 1} (a second brace or column in the same line), protect it against the event (impact barrier, blast or fire protection), or design it to survive the event itself.`,
            'NORSOK N-004 (2022) Annex A.2; ISO 19901-3:2014 §8; ISO 19902:2020 §7.10', `An added brace or column to restore redundancy is of the order of ${steelCo2Text(Math.max(gT / Math.max(inGroup.length, 1), 0.5))} — far less than strengthening every member to an elastic accidental design.`),
        adv('economic', 'info', `A redundancy factor of ${fmt(ds.redundancyFactor)} and a damaged strength ratio of ${fmt(lamDmin / lamI)} decide how urgent a repair is after damage: with reserve, production can continue with restrictions until a planned intervention.`,
          'Record the reserve and damaged-condition factors in the structural integrity file: they justify longer inspection intervals for redundant members and focus spending on the few critical ones.', 'Risk-based inspection per ISO 19901-9:2019 and DNV-RP-G101; whole-life costing per ISO 15686-5'),
      ],
    };
  },
  verify: [{
    // One-bay deck on four pinned legs CHS 914 × 32 (h = 4 m) pushed sideways: a sway mechanism with a plastic hinge
    // at the head of each leg (the girders are stronger). Wpl = (D³ − d³)/6 = (0.914³ − 0.850³)/6 = 0.024905 m³,
    // Mp = 355 MPa · Wpl = 8841 kN·m, H_ult = 4 Mp / h = 8841 kN. Axial forces are below 2 % of the squash load.
    name: 'Sway collapse of a portal on pinned legs: H = Σ Mp / h', source: 'Plastic analysis of frames (Neal, The Plastic Methods of Structural Analysis, ch. 3; hand calculation)',
    inputs: { pattern: 'lateral', Hx: 1000, gG: 1, gQ: 0, live: 0, tPlate: 0, secSteel: 0, nx: 1, ny: 1, nLevels: 1, bracing: 'none', columns: 'corners', support: 'pinned', damage: 'Beam', pDelta: false, equip: [] },
    expect: [{ kpi: 'ultLoad', value: 8841, tol: 0.03 }],
  }],
};
