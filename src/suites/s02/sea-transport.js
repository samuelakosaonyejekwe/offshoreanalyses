// Sea transport of the topsides on a barge or heavy-transport vessel: quasi-static accelerations from roll, pitch
// and heave at the cargo centre of gravity, inertia load cases on the deck frame standing on its transport
// supports, support reactions, uplift, seafastening forces and weld lengths, and member utilisation.
import { filletWeldResistance } from '../../engine/structural/steelcheck.js';
import { deckInputs, solverInput, buildDeck, loadCase, weights, solveDeck, memberForces, memberChecks, groupSummary, framePlot, nodeDisp, nodeValuesFromMembers, equipmentOf, need, fmt, round, linspace, adv, check, steelCo2Text, vesselCo2Text, CONFIRM, DECK_LIMITS, G, DEG } from './_common.js';

export default {
  id: 's02-sea-transport', suite: 2, title: 'Sea transport and seafastening',
  summary: 'Transport accelerations, inertia forces, support reactions and uplift of the topsides on a barge, the seafastening forces and weld lengths needed at each support, and member utilisation in the transport condition.',
  covers: [13], fidelity: 'numerical',
  standards: ['DNV-ST-N001 (2021) §11 (sea voyages: motion criteria, loads, seafastening and grillage design)', 'ISO 19901-6:2009 §12–§13 (transportation, seafastening)', 'Noble Denton 0030/ND (2016) §7 (default motion criteria)', 'EN 1993-1-8:2005 §4.5.3.3 (fillet welds, simplified method)', 'ISO 19902:2020 §13 / EN 1993-1-1 §6 (member checks)'],
  models: ['Transportation acceleration model (harmonic roll, pitch and heave about the motion centre)', 'Seafastening force equilibrium', 'Quasi-static inertia loading (d’Alembert)', 'Matrix stiffness method (3-D direct stiffness)', 'Fillet weld resistance (simplified method)'],
  limits: `${DECK_LIMITS} Accelerations are the harmonic single-amplitude values at the cargo centre of gravity applied uniformly to all mass, with roll and pitch taken separately and each combined with heave; phase relations, barge flexibility, wave slamming on overhangs, wind heel and friction are not modelled. Motion amplitudes are inputs: default criteria or a vessel motion analysis supply them.`,
  publishes: ['install.seafasteningForce', 'install.transportUplift'],
  searchTerms: 'sea transport barge motion criteria seafastening grillage topsides transportation accelerations',
  inputs: [
    ...deckInputs({ support: 'fixed', legs: 'all' }),
    { key: 'roll', label: 'Roll single amplitude', unit: 'deg', type: 'number', default: 20, min: 0, max: 40, group: 'Vessel motions', role: 'bc', bind: 'vessel.rollSigDeg', uncertain: { dist: 'normal', cov: 0.15 }, help: `Default criteria for large barges in unrestricted waters: 20° in 10 s roll, 10° (12.5°) in 10 s pitch, 0.2 g heave. A motion analysis for the route normally gives less. ${CONFIRM}` },
    { key: 'Troll', label: 'Roll period', unit: 's', type: 'number', default: 10, min: 3, max: 40, group: 'Vessel motions', role: 'bc', bind: 'vessel.Tn4' },
    { key: 'pitch', label: 'Pitch single amplitude', unit: 'deg', type: 'number', default: 10, min: 0, max: 25, group: 'Vessel motions', role: 'bc' },
    { key: 'Tpitch', label: 'Pitch period', unit: 's', type: 'number', default: 10, min: 3, max: 40, group: 'Vessel motions', role: 'bc' },
    { key: 'heave', label: 'Heave acceleration', unit: 'g', type: 'number', default: 0.2, min: 0, max: 1, group: 'Vessel motions', role: 'bc' },
    { key: 'zRoll', label: 'Transport supports above the motion centre', unit: 'm', type: 'number', default: 6, min: -10, max: 60, group: 'Stowage', role: 'geometry', help: 'Height of the leg bases (top of grillage) above the roll and pitch axis, usually taken at the waterline.' },
    { key: 'xOff', label: 'Deck centre from midship along the vessel', unit: 'm', type: 'number', default: 15, min: -150, max: 150, group: 'Stowage', role: 'geometry' },
    { key: 'yOff', label: 'Deck centre from the vessel centreline', unit: 'm', type: 'number', default: 0, min: -40, max: 40, group: 'Stowage', role: 'geometry' },
    { key: 'orient', label: 'Deck x axis points', type: 'select', default: 'along', group: 'Stowage', role: 'geometry', options: [{ value: 'along', label: 'Along the vessel (roll loads the deck in y)' }, { value: 'across', label: 'Across the vessel (roll loads the deck in x)' }] },
    { key: 'Fwind', label: 'Wind force on the cargo (transverse)', unit: 'kN', type: 'number', default: 300, min: 0, max: 2e4, group: 'Stowage', role: 'environment', help: 'Acts with the roll cases at the centre of the deck levels.' },
    { key: 'gF', label: 'Load factor on the transport loads', unit: '-', type: 'number', default: 1.3, min: 1, max: 1.7, group: 'Acceptance criteria', role: 'criteria', help: `DNV-ST-N001 LRFD: 1.3 on the combined weight and inertia in the ULS-b type combination. ${CONFIRM}` },
    { key: 'aWeld', label: 'Seafastening fillet weld throat', unit: 'mm', type: 'number', default: 10, min: 3, max: 30, group: 'Seafastening', role: 'geometry', design: true },
    { key: 'Lweld', label: 'Weld length available at each support', unit: 'mm', type: 'number', default: 3000, min: 100, max: 30000, group: 'Seafastening', role: 'geometry', design: true, help: 'Total effective length of fillet weld connecting the seafastening stoppers or braces to the barge deck at one support.' },
    { key: 'betaW', label: 'Weld correlation factor βw', unit: '-', type: 'number', default: 0.9, min: 0.8, max: 1, group: 'Seafastening', role: 'material', help: 'EN 1993-1-8 Table 4.1: 0.8 S235, 0.85 S275, 0.9 S355, 1.0 S420/S460 (the weaker of cargo and barge steel governs).' },
    { key: 'gM2', label: 'γM2 for welds', unit: '-', type: 'number', default: 1.25, min: 1, max: 1.5, group: 'Seafastening', role: 'criteria', help: CONFIRM },
    { key: 'Rdeck', label: 'Barge deck capacity under one support', unit: 'kN', type: 'number', default: 9000, min: 100, max: 2e5, group: 'Seafastening', role: 'criteria', help: 'Allowable load on the barge frame or bulkhead under a grillage support, from the vessel operator.' },
    { key: 'days', label: 'Transit duration', unit: 'days', type: 'number', default: 8, min: 0.5, max: 120, group: 'Stowage', role: 'operation' },
    solverInput(1),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxDisp' },
  compute(inp, ctx) {
    // SI: deg → rad, g → fraction of g in the load direction vector, kN → N, mm → m.
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp), w = weights(deck, inp, { eqs }), M = w.total / G;
    const th = inp.roll * DEG, ph = inp.pitch * DEG, wr = 2 * Math.PI / inp.Troll, wp = 2 * Math.PI / inp.Tpitch, zc = inp.zRoll + w.cog[2];
    // Accelerations at the cargo centre of gravity in vessel-fixed axes [m/s²]: tangential rotation term plus the
    // gravity component along the inclined deck; vertical term from the lever arm to the axis.
    const aT = wr * wr * th * zc + G * Math.sin(th), aTz = wr * wr * th * Math.abs(inp.yOff + (inp.orient === 'along' ? w.cog[1] - deck.y0 - deck.Ly / 2 : w.cog[0] - deck.x0 - deck.Lx / 2));
    const aL = wp * wp * ph * zc + G * Math.sin(ph), aLz = wp * wp * ph * Math.abs(inp.xOff);
    const tr = inp.orient === 'along' ? 1 : 0, zs = deck.levelsZ;   // index of the deck axis loaded by roll
    const windOn = (sgn) => { const e = []; deck.levelNodes.forEach((ns) => { for (const n of ns) { const F = [0, 0, 0, 0, 0, 0]; F[tr] = sgn * inp.gF * inp.Fwind * 1e3 / zs.length / ns.length; e.push({ node: n, F }); } }); return e; };
    const defs = [];
    for (const s of [1, -1]) for (const hv of [1, -1]) {
      const dr = [0, 0, -(Math.cos(th) * (1 + hv * inp.heave) + hv * aTz / G)]; dr[tr] = s * aT / G; defs.push({ name: `Roll ${s > 0 ? '+' : '−'}, heave ${hv > 0 ? 'down' : 'up'}`, dir: dr, extra: windOn(s), roll: true });
      const dp = [0, 0, -(Math.cos(ph) * (1 + hv * inp.heave) + hv * aLz / G)]; dp[1 - tr] = s * aL / G; defs.push({ name: `Pitch ${s > 0 ? '+' : '−'}, heave ${hv > 0 ? 'down' : 'up'}`, dir: dp, extra: [], roll: false });
    }
    const cases = defs.map((d, i) => loadCase(deck, inp, { name: `c${i}`, gDead: inp.gF, gEquip: inp.gF, gLive: 0, dir: d.dir, extra: d.extra, eqs }));
    cases.push(loadCase(deck, inp, { name: 'still', gDead: 1, gEquip: 1, gLive: 0, eqs }));
    const { ref, sol } = solveDeck(deck, cases, inp.nSub), nN = ref.model.nodes.length;
    const res = defs.map((d, i) => {
      const r = sol.results[`c${i}`], rows = memberChecks(deck, memberForces(deck, ref, r)), legs = r.reactions.filter((s) => deck.supNodes.includes(s.node));
      return { ...d, r, rows, gov: rows.reduce((a, c) => (c.uc > a.uc ? c : a)), legs, H: legs.map((l) => Math.hypot(l.R[0], l.R[1])), Ht: Math.hypot(legs.reduce((s, l) => s + l.R[0], 0), legs.reduce((s, l) => s + l.R[1], 0)), V: legs.map((l) => l.R[2]) };
    });
    const still = sol.results.still.reactions.filter((s) => deck.supNodes.includes(s.node)).map((l) => l.R[2]);
    const Hmax = Math.max(...res.map((c) => Math.max(...c.H))), Vmax = Math.max(...res.map((c) => Math.max(...c.V))), Vmin = Math.min(...res.map((c) => Math.min(...c.V)));
    const gm = res.reduce((a, c) => (c.gov.uc > a.gov.uc ? c : a)), gH = res.reduce((a, c) => (Math.max(...c.H) > Math.max(...a.H) ? c : a));
    const Hroll = Math.max(...res.filter((c) => c.roll).map((c) => c.Ht)), Hpitch = Math.max(...res.filter((c) => !c.roll).map((c) => c.Ht));
    // Seafastening weld at the most loaded support: shear, and uplift where the support lifts, on the same weld group.
    const fw = filletWeldResistance({ a: inp.aWeld / 1e3, fu: deck.steel.fu, betaW: inp.betaW, gammaM2: inp.gM2 }), uplift = Math.max(-Vmin, 0);
    const Fsf = Math.hypot(Hmax, uplift), Lreq = Fsf / fw.FwRd, groups = groupSummary(gm.rows);
    const checks = [
      check('Seafastening weld at the most loaded support', Lreq * 1e3, inp.Lweld, 'mm of weld', `EN 1993-1-8 §4.5.3.3, throat ${inp.aWeld} mm, βw = ${inp.betaW}, γM2 = ${inp.gM2}; DNV-ST-N001 §11.9 (friction not relied on)`, `Fillet weld length needed to carry the horizontal force of ${fmt(Hmax / 1e3)} kN${uplift > 0 ? ` with ${fmt(uplift / 1e3)} kN of uplift` : ''} at the most loaded support against the length available. Too little weld lets the cargo shift, which changes the stability of the barge.`),
      check(`Member strength in transport (${gm.gov.group.toLowerCase()} ${gm.gov.i + 1}, ${gm.name})`, gm.gov.uc, 1, 'interaction', gm.gov.ref, 'Governing member utilisation under the factored weight and inertia forces. Transport loads act sideways on a structure designed for gravity: legs, braces and their joints often need temporary transport bracing.'),
      check('Barge deck load under a support', Vmax / 1e3, inp.Rdeck, 'kN', 'DNV-ST-N001 §11.8 (grillage and vessel strength)', 'Largest factored vertical reaction passed through the grillage into the barge against the allowable load at a frame or bulkhead. Exceeding it needs load-spreading grillage beams or re-positioning over stronger structure.'),
    ];
    if (Vmin < 0) checks.push(check('Support uplift', uplift / 1e3, Math.max(...still) / 1e3, 'kN', 'DNV-ST-N001 §11.9 (uplift to be held by seafastening)', 'A support lifts in a roll or pitch case: the seafastening must hold the cargo down (the weld check above includes it; shown here against the still-water reaction for scale).'));
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = gov.uc <= 1;
    const u = gm.r.displacements, maxDisp = Math.max(...nodeDisp(u, nN));
    const sfT = deck.supNodes.length * 1.5;   // nominal seafastening steel per support [t], for the quantity estimate
    const rolls = linspace(0, Math.max(inp.roll * 1.5, 5), 16);
    return {
      kpis: [
        { key: 'ayRoll', label: 'Transverse acceleration in roll', unit: 'g', value: aT / G, explain: 'Acceleration along the deck at the cargo centre of gravity at maximum roll: rotation term ω²θz plus the gravity component g·sinθ.' },
        { key: 'axPitch', label: 'Longitudinal acceleration in pitch', unit: 'g', value: aL / G, explain: 'Same for pitch, along the vessel.' },
        { key: 'Hroll', label: 'Total transverse force in roll (factored)', unit: 'kN', value: Hroll / 1e3, explain: 'Sum of the horizontal support reactions in the worst roll case: cargo mass times the transverse acceleration, plus wind, times the load factor.' },
        { key: 'Hpitch', label: 'Total longitudinal force in pitch (factored)', unit: 'kN', value: Hpitch / 1e3, explain: 'Sum of the horizontal support reactions in the worst pitch case.' },
        { key: 'Hsupport', label: 'Largest horizontal force at one support', unit: 'kN', value: Hmax / 1e3, explain: `Design shear for the seafastening of the most loaded support (${gH.name}).` },
        { key: 'Vmax', label: 'Largest vertical support reaction', unit: 'kN', value: Vmax / 1e3, limit: inp.Rdeck, sense: 'max', explain: 'Largest factored compression on a transport support: the load into the grillage and barge deck.' },
        { key: 'Vmin', label: 'Smallest vertical support reaction', unit: 'kN', value: Vmin / 1e3, explain: 'Least loaded support over all cases; negative is uplift to be held by the seafastening.' },
        { key: 'LweldReq', label: 'Seafastening weld length required per support', unit: 'mm', value: Lreq * 1e3, limit: inp.Lweld, sense: 'max', explain: `Fillet weld of ${inp.aWeld} mm throat needed at the most loaded support (resistance ${fmt(fw.FwRd / 1e3)} kN per metre).` },
        { key: 'maxUc', label: 'Highest member utilisation in transport', unit: '-', value: round(gm.gov.uc, 5), limit: 1, sense: 'max', explain: 'Governing member utilisation over the eight roll and pitch cases.' },
        { key: 'maxDisp', label: 'Largest displacement in the governing case', unit: 'mm', value: maxDisp * 1e3, explain: 'Largest nodal displacement in the case that governs the members: sway of the upper decks relative to the grillage.' },
        { key: 'weight', label: 'Transport weight', unit: 't', value: M / 1e3, explain: 'Structure and equipment without live load.' },
      ],
      checks,
      plots: [
        framePlot({ title: `Utilisation and deformed shape — ${gm.name}`, ref, u, values: nodeValuesFromMembers(ref, gm.rows.map((x) => x.uc)), valueLabel: 'Utilisation', valueUnit: '-', explain: 'Frame under the governing transport case. The whole deck leans on its supports: look at the legs and braces between the grillage and the lowest deck.' }),
        { type: 'bar', title: 'Support reactions by transport case', categories: res.map((c) => c.name), series: [{ name: 'Largest vertical', y: res.map((c) => round(Math.max(...c.V) / 1e3)) }, { name: 'Smallest vertical', y: res.map((c) => round(Math.min(...c.V) / 1e3)) }, { name: 'Largest horizontal', y: res.map((c) => round(Math.max(...c.H) / 1e3)) }], y: { label: 'Reaction at one support', unit: 'kN' }, explain: 'Extreme support reactions for each roll and pitch case. Roll normally governs; a negative smallest vertical reaction is uplift.' },
        { type: 'line', title: 'Transverse acceleration against roll amplitude', x: { label: 'Roll single amplitude', unit: 'deg' }, y: { label: 'Transverse acceleration', unit: 'g' }, series: [{ name: 'Total', x: rolls.map((v) => round(v)), y: rolls.map((v) => round((wr * wr * v * DEG * zc + G * Math.sin(v * DEG)) / G)) }, { name: 'Gravity component g·sinθ', x: rolls.map((v) => round(v)), y: rolls.map((v) => round(Math.sin(v * DEG))), dash: true }, { name: 'Design roll', x: [inp.roll], y: [round(aT / G)], mode: 'markers' }],
          explain: 'How the transverse acceleration at the cargo centre of gravity grows with roll amplitude at the given roll period. The gap between the curves is the rotation term, which grows with height above the roll axis: a lower stowage or a motion analysis with smaller roll reduces every seafastening force.' },
        { type: 'bar', title: 'Member utilisation in transport by group', categories: groups.map((g) => g.group), series: [{ name: 'Utilisation', y: groups.map((g) => round(g.uc)) }], y: { label: 'Utilisation', unit: '-' }, explain: 'Governing utilisation of each member group in the governing transport case.' },
      ],
      tables: [{ title: 'Transport load cases', columns: [{ key: 'name', label: 'Case' }, { key: 'ah', label: 'Horizontal acceleration', unit: 'g' }, { key: 'av', label: 'Vertical acceleration', unit: 'g' }, { key: 'H', label: 'Total horizontal force', unit: 'kN' }, { key: 'Vmax', label: 'Max support', unit: 'kN' }, { key: 'Vmin', label: 'Min support', unit: 'kN' }, { key: 'uc', label: 'Max UC' }],
        rows: res.map((c) => ({ name: c.name, ah: round(Math.hypot(c.dir[0], c.dir[1])), av: round(-c.dir[2]), H: round(c.Ht / 1e3), Vmax: round(Math.max(...c.V) / 1e3), Vmin: round(Math.min(...c.V) / 1e3), uc: round(c.gov.uc) })),
        explain: 'Each roll and pitch case with its accelerations in units of g (before the load factor), the total horizontal force and the extreme support reactions (factored). These are the seafastening and grillage design loads.' }],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}. Centre of gravity ${fmt(zc)} m above the motion centre: roll gives ${fmt(aT / G)} g transverse, pitch ${fmt(aL / G)} g longitudinal.`, `Seafastening force at the most loaded support ${fmt(Fsf / 1e3)} kN → ${fmt(Lreq * 1e3)} mm of ${inp.aWeld} mm fillet weld.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The transport condition is acceptable: member utilisation ${fmt(gm.gov.uc)}, seafastening weld ${fmt(Lreq * 1e3)} mm needed of ${inp.Lweld} mm, support reactions ${fmt(Vmin / 1e3)} to ${fmt(Vmax / 1e3)} kN.`,
          `Design each seafastening for ${fmt(Hmax / 1e3)} kN horizontal${uplift > 0 ? ` and ${fmt(uplift / 1e3)} kN uplift` : ''}, align it with a barge frame or bulkhead, and re-run with route-specific motions if a motion analysis is available (they are usually lower than the default criteria).`, 'DNV-ST-N001 (2021) §11; ISO 19901-6:2009 §12–§13',
          `Transit of ${inp.days} days: ${vesselCo2Text(inp.days)} Seafastening and grillage steel (about ${fmt(sfT)} t here) is temporary: specify it for cutting and re-use or recycling.`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): horizontal force ${fmt(Hmax / 1e3)} kN at one support, member utilisation ${fmt(gm.gov.uc)}, support reactions ${fmt(Vmin / 1e3)} to ${fmt(Vmax / 1e3)} kN.`,
            gov.name.startsWith('Seafastening') ? `Provide ${fmt(Math.ceil(Lreq * 10) * 100, 3)} mm of ${inp.aWeld} mm fillet weld per support, or a ${fmt(Math.ceil(inp.aWeld * gov.uc), 2)} mm throat on the present length, or add seafastening braces to more supports.` : gov.name.startsWith('Member') ? `Add temporary transport braces between the grillage and the lowest deck, or strengthen the ${gm.gov.group.toLowerCase()} group by about ${fmt((gm.gov.uc - 1) * 100, 2)} %; lowering the roll criterion to ${fmt(inp.roll / gm.gov.uc, 2)}° by a route-specific motion analysis is the alternative.` : gov.name.startsWith('Barge') ? 'Spread the support over two barge frames with a longer grillage beam, or move the cargo over a bulkhead.' : `Hold the lifting support down for ${fmt(uplift / 1e3)} kN with tie-down plates or braces.`,
            'DNV-ST-N001 (2021) §11.6–11.9; EN 1993-1-8:2005 §4.5.3.3', `Temporary bracing of about ${steelCo2Text(Math.max(gm.gov.uc - 1, 0.03) * deck.frameSteelT * 0.3)} is used once: design it to be removed and recycled, and prefer a calmer route or season where the schedule allows.`),
        adv('economic', 'info', `Default motion criteria give ${fmt(aT / G)} g transverse on a ${fmt(M / 1e3)} t cargo; seafastening is designed for ${fmt(Hroll / 1e3)} kN in roll and the barge is occupied for ${inp.days} days of transit.`,
          'A route- and season-specific motion analysis typically cuts the design accelerations well below the default criteria and with them the seafastening, grillage and transport bracing; weigh its cost against that steel and the welding and cutting time on the barge.', 'DNV-ST-N001 (2021) §11.3 (motion analysis versus default criteria); AACE RP 18R-97'),
      ],
      publish: { 'install.seafasteningForce': Hmax, 'install.transportUplift': uplift },
      quantities: { steel_t: sfT, vessel_days: { transportBarge: inp.days } },
    };
  },
  verify: [{
    // 2059.63 t deck (hand take-off of the in-place module) with its centre of gravity 3.989 m above the supports,
    // supports 10 m above the roll axis: z = 13.989 m. Roll 20° (0.34907 rad) in 10 s, no heave, wind or load factor:
    //   a = ω²θz + g sinθ = 0.39478·0.34907·13.989 + 9.80665·0.34202 = 1.9278 + 3.3541 = 5.2819 m/s² = 0.5386 g
    //   H = M a = 2 059 630·5.2819 = 10 879 kN.
    name: 'Roll inertia force on the cargo: H = M (ω²θz + g sinθ)', source: 'DNV-ST-N001 §11 quasi-static motion loads (hand calculation)',
    inputs: { roll: 20, Troll: 10, heave: 0, zRoll: 10, yOff: 0, xOff: 0, Fwind: 0, gF: 1, live: 0, tPlate: 0, secSteel: 0, nx: 1, ny: 1, nLevels: 1, bracing: 'none', columns: 'corners', support: 'fixed', equip: [{ name: 'Block', W: 2000, x: 15, y: 10, level: 1, a: 0, b: 0, zc: 0 }] },
    expect: [{ kpi: 'ayRoll', value: 0.5386, tol: 2e-3 }, { kpi: 'Hroll', value: 10879, tol: 2e-3 }],
  }],
};
