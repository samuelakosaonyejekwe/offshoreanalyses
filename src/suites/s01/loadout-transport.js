// Load-out, transportation and seafastening: the jacket lies on one face (or stands on its leg bases) on the barge. Barge roll, pitch
// and heave give inertia loads through an acceleration field plus the gravity component of the heel angle;
// load-out is checked with one support settled out of level. Seafastening forces, weld length and tipping margin follow.
import { G, DEG, need, fmt, round, adv, check, steelCo2, vesselCo2, CODE_CAVEAT, CONFIRM, BASIS_MEMBER, BASIS_COST } from './_common.js';
import { geomInputs, criteriaInputs, meshInputs, buildJacket, solve, memberRows, envelope, byGroup, displacements, nodeValuesFromMembers, plot3d, memberCols, memberTableRows, groupChecks, testFrame } from './_model.js';
import { filletWeldResistance } from '../../engine/structural/steelcheck.js';

const crit = criteriaInputs({ gG: 1.3 }).filter((i) => i.key !== 'gE');

export default {
  id: 's01-loadout-transport', suite: 1, title: 'Load-out, transportation and seafastening',
  summary: 'Member utilisation, support reactions and seafastening forces of the jacket on its transport barge under roll, pitch and heave, the tipping margin, and the racking stresses of a load-out with one support out of level.',
  covers: [16], fidelity: 'numerical',
  standards: ['DNV-ST-N001 (2021) §5 (load-out), §11 (sea voyages, default motion criteria), §11.9 (seafastening)', 'ISO 19901-6:2009 §12 (transportation), §11 (load-out)', 'ISO 19902:2020 §8 and §22', 'API RP 2A-WSD 22nd ed. §5.4.2', 'EN 1993-1-8:2005 §4.5.3.3 (fillet welds)'],
  models: ['Rigid-body barge motion as an acceleration field (roll, pitch, heave) with the gravity component of the heel angle', 'Newton\'s second law (d\'Alembert inertia loads on distributed member mass)', 'Direct stiffness method with prescribed support displacement', 'Coulomb friction for the skidding force', 'Fillet weld resistance (simplified method)', 'Tubular member code check'],
  limits: `The jacket rests on point supports (the framing nodes of one face, or the leg bases) on a rigid barge: barge flexibility, the distributed support of launch rails, wave slam on overhanging members, wind on the cargo, fatigue during the tow and barge stability are not modelled, and motions are single-amplitude design values, not a motion-response analysis. An imported model is analysed in the position it is imported in. Motion criteria and factors are to be confirmed against the governing edition and the tow route. ${CODE_CAVEAT}`,
  searchTerms: 'jacket transportation seafastening barge motions load-out DNV-ST-N001',
  inputs: [
    ...geomInputs().filter((i) => !['depth'].includes(i.key)),
    { key: 'depth', label: 'Jacket height below the top reference (water depth in place)', unit: 'm', type: 'number', default: 80, min: 5, max: 400, group: 'Jacket geometry', role: 'geometry', bind: 'site.depth' },
    { key: 'pos', label: 'Transport position', type: 'select', default: 'side', group: 'Jacket geometry', role: 'ic', options: [{ value: 'side', label: 'On its side, on one face (launch or lift-and-upend jackets)' }, { value: 'upright', label: 'Upright on the leg bases (short jackets)' }], help: 'On its side the jacket rests on the two legs of one face at every framing level, with its axis along the barge. An imported model is analysed as imported, supported at its lowest nodes.' },
    { key: 'swCont', label: 'Weight contingency', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Weight', role: 'load', help: 'Factor on modelled steel for joint cans, anodes, mudmats and appurtenances (weight report class).' },
    { key: 'roll', label: 'Roll single amplitude', unit: 'deg', type: 'number', default: 20, min: 0, max: 40, group: 'Barge motions (prescribed)', role: 'bc', help: `Default criteria of DNV-ST-N001 for large barges in unrestricted tows: 20° roll and 12.5° pitch in 10 s with ±0.2 g heave. Route-specific motion analysis usually gives less. ${CONFIRM}` },
    { key: 'pitch', label: 'Pitch single amplitude', unit: 'deg', type: 'number', default: 12.5, min: 0, max: 30, group: 'Barge motions (prescribed)', role: 'bc' },
    { key: 'Tmot', label: 'Roll and pitch period', unit: 's', type: 'number', default: 10, min: 3, max: 30, group: 'Barge motions (prescribed)', role: 'bc' },
    { key: 'heave', label: 'Heave acceleration', unit: 'g', type: 'number', default: 0.2, min: 0, max: 1, group: 'Barge motions (prescribed)', role: 'bc' },
    { key: 'hRc', label: 'Supports above the centre of rotation', unit: 'm', type: 'number', default: 5, min: -10, max: 30, group: 'Barge motions (prescribed)', role: 'geometry', help: 'Height of the support level (barge deck and grillage) above the roll centre, normally taken at the barge waterline.' },
    { key: 'sfType', label: 'Support restraint on the barge', type: 'select', default: 'pinned', group: 'Seafastening (boundary condition)', role: 'bc', options: [{ value: 'pinned', label: 'Pinned (grillage bearing + seafastening braces)' }, { value: 'fixed', label: 'Fixed (welded into a stiff support frame)' }] },
    { key: 'weldA', label: 'Seafastening fillet weld throat', unit: 'mm', type: 'number', default: 10, min: 3, max: 30, group: 'Seafastening (boundary condition)', role: 'geometry' },
    { key: 'weldL', label: 'Weld length available per support', unit: 'mm', type: 'number', default: 3000, min: 100, max: 20000, group: 'Seafastening (boundary condition)', role: 'geometry', design: true, help: 'Total effective fillet length of the seafastening braces and tie-down plates at one support.' },
    { key: 'settle', label: 'Load-out support level tolerance', unit: 'mm', type: 'number', default: 25, min: 0, max: 300, group: 'Load-out', role: 'bc', help: 'One support is taken this much lower than the others while the jacket moves from quay to barge (ballasting and skidway tolerance).' },
    { key: 'mu', label: 'Skidding friction coefficient (break-out)', unit: '-', type: 'number', default: 0.2, min: 0.01, max: 0.6, group: 'Load-out', role: 'criteria', help: `About 0.2 for steel on teflon or greased timber at break-out, 0.1 moving (DNV-ST-N001 load-out friction table). ${CONFIRM}` },
    { key: 'towDays', label: 'Tow duration', unit: 'days', type: 'number', default: 6, min: 0.5, max: 90, group: 'Load-out', role: 'operation' },
    ...crit, ...meshInputs({ nSub: 2, phase: false }),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxDisp' },
  compute(inp, ctx) {
    // SI conversion: deg → rad, g → m/s², mm → m. The model is in air: no buoyancy, no water in the legs.
    const J = buildJacket({ ...inp, found: inp.sfType }, ctx, { found: inp.sfType, layDown: inp.pos === 'side' });
    J.model.members.forEach((m) => { m.flooded = false; });
    need(J.mudNodes.length >= 3 || inp.sfType === 'fixed', 'A pinned cargo needs at least three supports that are not in one line: choose the fixed restraint or import a model with more support nodes.');
    const P = J.model.nodes, zB = J.zMud, org = [J.centre[0], J.centre[1], zB - inp.hRc], w2 = (2 * Math.PI / inp.Tmot) ** 2, f = inp.gG * inp.swCont;
    const zmid = (m) => 0.5 * (J.base.nodes[m.a][2] + J.base.nodes[m.b][2]);
    const W = J.steelT * 1e3 * inp.swCont * G, zc = J.orig.reduce((s, m) => s + m.sec.A * m.L * zmid(m), 0) / J.orig.reduce((s, m) => s + m.sec.A * m.L, 0) - zB;
    // motion cases: at the extreme angle the angular acceleration is −ω²θ, so gravity component and inertia act together
    const mk = (name, th, axis, sgn, hv) => { const t = th * DEG * sgn, a = [0, 0, 0], g = [0, 0, -(G * Math.cos(t) + hv * inp.heave * G)]; if (axis === 'roll') { a[0] = -w2 * t; g[1] = -G * Math.sin(t); } else { a[1] = -w2 * t; g[0] = G * Math.sin(t); } return { name, nodal: [], member: [], selfWeight: { g, factor: f }, acceleration: { alpha: a.map((v) => v * f), origin: org } }; };
    const cases = [];
    for (const [axis, th] of [['roll', inp.roll], ['pitch', inp.pitch]]) for (const s of [1, -1]) for (const hv of [1, -1]) cases.push(mk(`${axis} ${s > 0 ? '+' : '−'} / heave ${hv > 0 ? 'down' : 'up'}`, th, axis, s, hv));
    const still = { name: 'Still water', nodal: [], member: [], selfWeight: { factor: f } };
    const lo = { name: 'Load-out, one support low', nodal: [], member: [], selfWeight: { factor: f }, settlements: inp.settle > 0 ? [{ node: J.mudNodes[0], d: [0, 0, -inp.settle / 1e3, 0, 0, 0] }] : [] };
    const res = solve(J, [...cases, still, lo], 'Transport analysis');
    const rowsT = cases.map((c) => memberRows(J, res[c.name], inp, { name: c.name })), envT = envelope(rowsT), rowsL = memberRows(J, res[lo.name], inp, { name: lo.name }), groups = byGroup(envT), gL = rowsL.reduce((a, r) => (r.uc > a.uc ? r : a)), gT = envT.reduce((a, r) => (r.uc > a.uc ? r : a));
    // support reactions per case
    const sup = cases.map((c) => { const r = res[c.name].reactions; return { name: c.name, up: Math.max(0, ...r.map((x) => -x.R[2])), dn: Math.max(...r.map((x) => x.R[2])), sh: Math.max(...r.map((x) => Math.hypot(x.R[0], x.R[1]))), H: Math.hypot(r.reduce((s, x) => s + x.R[0], 0), r.reduce((s, x) => s + x.R[1], 0)), V: r.reduce((s, x) => s + x.R[2], 0) }; });
    const up = Math.max(...sup.map((s) => s.up)), sh = Math.max(...sup.map((s) => s.sh)), dn = Math.max(...sup.map((s) => s.dn)), Hmax = Math.max(...sup.map((s) => s.H)), gS = sup.reduce((a, s) => (s.sh > a.sh ? s : a));
    // tipping margin of the unfastened cargo about the line of supports (rigid body, unfactored)
    const bx = Math.max(...J.mudNodes.map((n) => P[n][0])) - Math.min(...J.mudNodes.map((n) => P[n][0])), by = Math.max(...J.mudNodes.map((n) => P[n][1])) - Math.min(...J.mudNodes.map((n) => P[n][1]));
    const tip = (th, b) => { const t = th * DEG, lat = G * Math.sin(t) + w2 * t * (inp.hRc + zc), ver = G * Math.cos(t) - inp.heave * G; return lat > 0 ? ver * (b / 2) / (lat * zc) : 999; };
    const tipR = tip(inp.roll, by), tipP = tip(inp.pitch, bx), tipMin = Math.min(tipR, tipP), aTop = (G * Math.sin(inp.roll * DEG) + w2 * inp.roll * DEG * (inp.hRc + J.zTop - zB)) / G;
    const fw = filletWeldResistance({ a: inp.weldA / 1e3, fu: J.steel.fu, betaW: 0.9 }), Fsf = Math.hypot(sh, up), Lreq = Fsf / fw.FwRd, dp = Math.max(...cases.map((c) => displacements(J, res[c.name]).max)), pull = inp.mu * W;
    const checks = [
      ...groupChecks(groups, inp, 'transport envelope'),
      check('Member strength during load-out (one support low)', gL.uc, 1, 'interaction', `${inp.code} member check with a support ${fmt(inp.settle)} mm out of level; DNV-ST-N001 §5`, 'Racking of the frame when the supports are not in one plane during the move onto the barge. A stiff, pinned-base jacket attracts large brace forces from small level differences.'),
      check('Seafastening weld at the most loaded support', Lreq * 1e3, inp.weldL, 'mm', `EN 1993-1-8 §4.5.3.3 (fvw,d = fu/(√3·βw·γM2), βw = 0.9, γM2 = 1.25); DNV-ST-N001 §11.9`, 'Fillet weld length needed to carry the resultant of shear and uplift at one support against the length available. If the seafastening fails the cargo slides or tips in the design roll.'),
    ];
    const steel = J.steelT * inp.swCont, ok = checks.every((c) => c.uc <= 1), bad = checks.filter((c) => c.uc > 1), sfT = Lreq * inp.weldA / 1e3 * 0.03 * 7850 / 1e3 * J.mudNodes.length;
    const gg = groups.reduce((a, g) => (g.uc > a.uc ? g : a));
    return {
      kpis: [
        { key: 'weight', label: 'Transport weight', unit: 't', value: steel, explain: 'Modelled jacket steel times the weight contingency: the cargo weight for barge selection, load-out and seafastening.' },
        { key: 'cogZ', label: 'Centre of gravity above the supports', unit: 'm', value: zc, explain: 'Height of the centre of gravity above the barge deck. Together with the support spacing it decides tipping and the uplift on the seafastening.' },
        { key: 'ucTransport', label: 'Highest member utilisation in transport', unit: '-', value: gT.uc, limit: 1, sense: 'max', explain: 'Largest code utilisation over the roll and pitch cases with heave. Members near the supports and long horizontal braces usually govern.' },
        { key: 'ucLoadout', label: 'Highest member utilisation in load-out', unit: '-', value: gL.uc, limit: 1, sense: 'max', explain: 'Largest utilisation under self-weight with one support out of level.' },
        { key: 'maxUplift', label: 'Largest uplift at a support', unit: 'MN', value: up / 1e6, explain: 'Largest factored tension on a support: the tie-down load of the seafastening. Zero means the cargo stays seated.' },
        { key: 'maxShearSupport', label: 'Largest horizontal reaction at a support', unit: 'MN', value: sh / 1e6, explain: 'Largest factored horizontal force at one support: the design load of the seafastening braces and stoppers.' },
        { key: 'maxBearing', label: 'Largest bearing reaction at a support', unit: 'MN', value: dn / 1e6, explain: 'Largest factored compression at one support: the load on the grillage and on the barge frames below.' },
        { key: 'baseShear', label: 'Total horizontal force on the cargo', unit: 'MN', value: Hmax / 1e6, explain: 'Resultant factored horizontal load from heel and angular acceleration in the governing motion case.' },
        { key: 'tipFoS', label: 'Tipping factor without seafastening', unit: '-', value: Math.min(tipMin, 999), explain: 'Unfactored restoring over overturning moment of the cargo as a rigid body at the design roll or pitch. Below 1.0 the cargo would lift off without tie-downs: the uplift reaction then is a structural seafastening load, which the weld check covers.' },
        { key: 'aLatTop', label: 'Lateral acceleration at the jacket top in roll', unit: 'g', value: aTop, explain: 'Heel component plus angular acceleration at the top of the cargo. Values near 0.7–0.9 g are typical of the default criteria for tall cargo.' },
        { key: 'weldReq', label: 'Seafastening weld length required per leg', unit: 'mm', value: Lreq * 1e3, limit: inp.weldL, sense: 'max', explain: 'Fillet length of the given throat needed at the most loaded support.' },
        { key: 'pullForce', label: 'Skidding force at break-out', unit: 'MN', value: pull / 1e6, explain: 'Friction coefficient times the transport weight: the winch or strand-jack capacity needed to start the load-out.' },
        { key: 'maxDisp', label: 'Largest displacement in transport', unit: 'mm', value: dp * 1e3, explain: 'Largest nodal displacement relative to the barge deck over the motion cases.' },
      ],
      checks,
      plots: [
        plot3d(J, { title: `Transport utilisation envelope and deformed shape (${gT.case})`, values: nodeValuesFromMembers(J, envT.map((r) => r.uc)), label: 'Utilisation', unit: '-', u: res[gT.case].displacements, explain: 'Envelope of member utilisation over all motion cases at the nodes, on the deflected shape of the governing case. Loads act sideways on a structure designed for vertical loads, so the lower bays work hardest.' }),
        { type: 'bar', title: 'Support reactions by motion case', categories: sup.map((s) => s.name), series: [{ name: 'Uplift (MN)', y: sup.map((s) => round(s.up / 1e6)) }, { name: 'Shear (MN)', y: sup.map((s) => round(s.sh / 1e6)) }, { name: 'Bearing (MN)', y: sup.map((s) => round(s.dn / 1e6)) }], y: { label: 'Largest reaction at one support', unit: 'MN' }, explain: 'Largest uplift, horizontal and bearing reaction at a single support for each motion case. Roll normally governs; heave up with roll gives the largest uplift.' },
        { type: 'bar', title: 'Member utilisation: transport and load-out', categories: groups.map((g) => g.label), series: [{ name: 'Transport', y: groups.map((g) => round(g.uc)) }, { name: 'Load-out', y: groups.map((g) => round(Math.max(...rowsL.filter((r) => r.group === g.group).map((r) => r.uc)))) }], y: { label: 'Utilisation', unit: '-' }, explain: 'Highest utilisation by member group in the two temporary conditions. Temporary conditions often govern the lower horizontal bracing.' },
      ],
      tables: [
        { title: 'Motion cases and reactions', columns: [{ key: 'name', label: 'Case' }, { key: 'H', label: 'Horizontal resultant', unit: 'MN' }, { key: 'V', label: 'Vertical resultant', unit: 'MN' }, { key: 'up', label: 'Uplift', unit: 'MN' }, { key: 'sh', label: 'Shear', unit: 'MN' }, { key: 'dn', label: 'Bearing', unit: 'MN' }, { key: 'uc', label: 'Member utilisation' }],
          rows: sup.map((s, i) => ({ name: s.name, H: round(s.H / 1e6), V: round(s.V / 1e6), up: round(s.up / 1e6), sh: round(s.sh / 1e6), dn: round(s.dn / 1e6), uc: round(Math.max(...rowsT[i].map((r) => r.uc))) })), explain: 'Factored resultants and the largest single-support reactions for each combination of roll or pitch with heave. These are the design loads of grillage, seafastening and barge deck.' },
        { title: 'Most utilised members in transport', columns: memberCols, rows: memberTableRows(envT, 12), explain: 'Members with the highest utilisation over the motion cases and the case that governs each.' },
        { title: 'Tipping and weight', columns: [{ key: 'c', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }], rows: [{ c: 'Support spacing along x / y', v: `${fmt(bx)} / ${fmt(by)}`, u: 'm' }, { c: 'Tipping factor in roll', v: round(Math.min(tipR, 999)), u: '-' }, { c: 'Tipping factor in pitch', v: round(Math.min(tipP, 999)), u: '-' }, { c: 'Angular acceleration in roll', v: round(w2 * inp.roll * DEG), u: 'rad/s²' }, { c: 'Weld design resistance', v: round(fw.FwRd / 1e6), u: 'MN/m' }, { c: 'Load-out utilisation, governing member', v: `${round(gL.uc)} (member ${gL.i})`, u: '-' }], explain: 'Geometry of the support base, rigid-body tipping margins and the seafastening weld resistance used.' },
      ],
      notes: [`Transport weight ${fmt(steel)} t, centre of gravity ${fmt(zc)} m above the deck; governing case "${gT.case}" with utilisation ${fmt(gT.uc)} in ${gT.label.toLowerCase()} member ${gT.i}.`,
        `Seafastening: ${fmt(sh / 1e6)} MN shear and ${fmt(up / 1e6)} MN uplift at one support (case "${gS.name}") need ${fmt(Lreq * 1e3)} mm of ${fmt(inp.weldA)} mm fillet; tipping factor ${fmt(Math.min(tipMin, 999))} without tie-downs.`,
        `Load-out with one support ${fmt(inp.settle)} mm low: utilisation ${fmt(gL.uc)}; break-out skidding force ${fmt(pull / 1e6)} MN. ${J.mudNodes.length} ${inp.sfType} supports, jacket ${J.imported ? 'as imported' : inp.pos === 'side' ? 'on its side' : 'upright'}.`],
      advice: [
        ok ? adv('engineering', 'ok', `Transport and load-out are acceptable: utilisation ${fmt(gT.uc)} (transport) and ${fmt(gL.uc)} (load-out), seafastening weld ${fmt(Lreq * 1e3)} mm of ${fmt(inp.weldL)} mm available.`, tipMin < 1 ? `The cargo would tip without tie-downs (factor ${fmt(tipMin)}): detail the uplift connection for ${fmt(up / 1e6)} MN per support and check the barge deck for it. A route-specific motion analysis would normally cut the roll loads.` : 'Check the barge deck and frames for the bearing reactions, and replace the default motions by a motion analysis for the actual barge and route if weight is to be saved.', 'DNV-ST-N001 (2021) §11; ISO 19901-6:2009 §12', `Seafastening steel is scrapped after the tow: about ${fmt(sfT)} t of weld metal and several tens of tonnes of braces. Design it for reuse or minimum cutting; the tow itself is about ${vesselCo2(inp.towDays)}`)
          : adv('engineering', 'fail', `Temporary-condition criteria are not met: ${bad.map((c) => `${c.name} (${fmt(c.uc)})`).join('; ')}.`, bad[0].name.startsWith('Seafastening') ? `Provide at least ${fmt(Lreq * 1e3)} mm of ${fmt(inp.weldA)} mm fillet per leg or a larger throat.` : bad[0].name.includes('load-out') ? `Tighten the level tolerance (utilisation scales with it: about ${fmt(inp.settle / gL.uc * 0.9)} mm for 0.9) or use hydraulic load equalisation on the skid shoes.` : (gg.uc >= 50 ? `The ${gg.label.toLowerCase()} buckle: their compression exceeds the elastic critical load. Add temporary transport bracing across the long panels, support more nodes on the barge, or use route-specific motions.` : `Add temporary transport bracing or use route-specific motions; the ${gg.label.toLowerCase()} need about ${fmt((gg.uc / 0.9 - 1) * 100)} % more capacity as it stands.`), BASIS_MEMBER + '; DNV-ST-N001 (2021) §11', `Permanent steel added for a temporary condition is carried for the whole life: temporary bracing of about ${steelCo2(J.steelT * 0.02)} is usually the better choice.`),
        adv('economic', 'info', `A ${fmt(steel)} t cargo with a ${fmt(inp.towDays)}-day tow; default motion criteria give ${fmt(aTop)} g at the jacket top.`, 'Barge day rate and weather waiting dominate transport cost. Route-specific motions and a restricted (weather-routed) tow shrink seafastening and transport bracing, at the price of waiting time: compare the two with the weather-window statistics of the route.', `${BASIS_COST}; marine operation planning per DNV-ST-N001 §2`),
      ],
      quantities: { steel_t: round(steel), vessel_days: { tow: inp.towDays } },
    };
  },
  verify: [
    // Vertical tube (D 1.5 m, t 45 mm, L 20 m, m = 1614.71·20 = 32 294 kg) fixed to the deck 5 m above the roll centre;
    // roll 20° in 10 s, no heave, factors 1.0. Lateral force = m·[g·sin20° + (2π/10)²·0.34907·(5 + 10)]
    // = 32 294·(3.3541 + 2.0671) = 175.07 kN; the top sees (3.3541 + 0.13781·25)/9.80665 = 0.6933 g.
    { name: 'Roll inertia of a uniform column: F = m·(g·sinθ + ω²θ·h_cg)', source: 'Rigid-body kinematics of simple harmonic roll (DNV-ST-N001 §11 simplified motion loads; Noble Denton 0030/ND §7)',
      inputs: { geom: testFrame('column', { H: 20 }), swCont: 1, gG: 1, roll: 20, pitch: 0, heave: 0, Tmot: 10, hRc: 5, sfType: 'fixed', settle: 0, nSub: 4 },
      expect: [{ kpi: 'baseShear', value: 0.17507, tol: 0.002 }, { kpi: 'aLatTop', value: 0.6933, tol: 0.002 }] },
  ],
};
